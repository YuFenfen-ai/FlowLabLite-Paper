from __future__ import annotations

import argparse
import csv
import json
import platform
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence

from validation.ghia_study import GHIA_REFERENCE, interpolate, parse_centerline, write_summary


_CENTERLINE_PREFIX = "# FlowLabLite centerline velocity"
_STEP_PATTERN = re.compile(r"(?:^|\|)\s*step=(\d+)(?:\s*\||$)")


def split_checkpoint_stream(stream: str) -> dict[int, str]:
    """Split concatenated centerline CSV blocks and return them in emitted order."""
    blocks: dict[int, str] = {}
    current: list[str] = []

    def finish() -> None:
        if not current:
            return
        header = current[0]
        match = _STEP_PATTERN.search(header)
        if match is None or len(current) < 2 or current[1] != "component,index,coord,value":
            raise ValueError("checkpoint stream contains a malformed centerline block")
        step = int(match.group(1))
        if step in blocks:
            raise ValueError(f"checkpoint stream contains duplicate step {step}")
        blocks[step] = "\n".join(current).rstrip() + "\n"

    for raw_line in stream.splitlines():
        if raw_line.startswith(_CENTERLINE_PREFIX):
            finish()
            current = [raw_line]
        elif current:
            if raw_line.strip() or len(current) == 1:
                current.append(raw_line)
    finish()
    if not blocks:
        raise ValueError("solver produced no centerline checkpoint blocks")
    if any(right <= left for left, right in zip(blocks, tuple(blocks)[1:])):
        raise ValueError("solver checkpoint blocks are not strictly increasing")
    return blocks


def build_solver_command(
    moon: Path, grid: int, checkpoints: Sequence[int], reynolds: int, dt: float = 0.001
) -> list[str]:
    return [
        str(moon), "run", "cmd/main", "--release", "--target", "wasm", "--",
        "--format", "centerline", "--solver", "chorin", "--grid", str(grid),
        "--checkpoints", ",".join(str(step) for step in checkpoints), "--re", str(reynolds),
        "--dt", format(dt, "g"),
    ]


def _write_pointwise(selected_path: Path, output_path: Path, reynolds: int, dt: float) -> None:
    data = parse_centerline(selected_path, expected_re=reynolds, expected_dt=dt)
    reference = GHIA_REFERENCE[reynolds]
    u_predicted = interpolate(data.u_coord, data.u_value, reference["u_coord"])
    v_predicted = interpolate(data.v_coord, data.v_value, reference["v_coord"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("re", "component", "coord_over_l", "reference", "flowlablite", "signed_error", "abs_error"))
        for component, coordinates, reference_values, predicted in (
            ("u", reference["u_coord"], reference["u_value"], u_predicted),
            ("v", reference["v_coord"], reference["v_value"], v_predicted),
        ):
            for coordinate, expected, actual in zip(coordinates, reference_values, predicted, strict=True):
                error = actual - expected
                writer.writerow((reynolds, component, coordinate / 2.0, expected, actual, error, abs(error)))


def run_case(
    repo: Path,
    moon: Path,
    output_root: Path,
    grid: int,
    checkpoints: Sequence[int],
    reynolds: int,
    *,
    dt: float,
    resume: bool,
) -> dict[str, object]:
    case_dir = output_root / f"re{reynolds}" / f"g{grid}"
    case_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_paths = {step: case_dir / f"centerline_s{step}.csv" for step in checkpoints}
    can_resume = resume and all(path.is_file() for path in checkpoint_paths.values())
    command = build_solver_command(moon, grid, checkpoints, reynolds, dt)
    if can_resume:
        for path in checkpoint_paths.values():
            data = parse_centerline(path, expected_re=reynolds, expected_dt=dt)
            if data.metadata["grid"] != f"{grid}x{grid}":
                raise ValueError(f"resume file has the wrong grid: {path}")
        elapsed_seconds = None
    else:
        started = time.perf_counter()
        completed = subprocess.run(command, cwd=repo, capture_output=True, text=True, encoding="utf-8")
        elapsed_seconds = time.perf_counter() - started
        (case_dir / "solver_stderr.log").write_text(completed.stderr, encoding="utf-8")
        (case_dir / "checkpoint_stream.txt").write_text(completed.stdout, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"solver failed for Re={reynolds}; see {case_dir / 'solver_stderr.log'}")
        blocks = split_checkpoint_stream(completed.stdout)
        if tuple(blocks) != tuple(checkpoints):
            raise ValueError(f"expected checkpoints {tuple(checkpoints)}, received {tuple(blocks)}")
        for step, text in blocks.items():
            path = checkpoint_paths[step]
            path.write_text(text, encoding="utf-8")
            data = parse_centerline(path, expected_re=reynolds, expected_dt=dt)
            if data.metadata["grid"] != f"{grid}x{grid}":
                raise ValueError(f"solver emitted the wrong grid for step {step}")

    summary_csv = case_dir / "time_convergence.csv"
    summary_json = case_dir / "time_convergence.json"
    summary = write_summary(tuple(checkpoint_paths.values()), summary_csv, summary_json, expected_dt=dt)
    _write_pointwise(checkpoint_paths[checkpoints[-1]], case_dir / "pointwise_errors.csv", reynolds, dt)
    manifest = {
        "re": reynolds,
        "grid": grid,
        "checkpoints": list(checkpoints),
        "command": command,
        "elapsed_seconds": elapsed_seconds,
        "resumed": can_resume,
        "python": sys.version,
        "platform": platform.platform(),
        "summary_json": str(summary_json.resolve()),
        "time_converged": summary["time_converged"],
        "selected_time_step": summary["selected_time_step"],
    }
    (case_dir / "run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run reproducible 129x129 Ghia validation cases.")
    parser.add_argument("--moon", type=Path, required=True)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--output-root", type=Path, default=Path("validation/results/ghia"))
    parser.add_argument("--grid", type=int, default=129)
    parser.add_argument("--checkpoints", type=int, nargs="+", default=(5000, 10000, 20000))
    parser.add_argument("--re", type=int, nargs="+", default=(100, 400, 1000))
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--resume", action="store_true")
    arguments = parser.parse_args(argv)
    checkpoints = tuple(arguments.checkpoints)
    if any(step <= 0 for step in checkpoints) or any(
        right <= left for left, right in zip(checkpoints, checkpoints[1:])
    ):
        parser.error("--checkpoints must be strictly increasing positive integers")
    if arguments.dt <= 0.0:
        parser.error("--dt must be positive")
    unsupported = set(arguments.re) - set(GHIA_REFERENCE)
    if unsupported:
        parser.error(f"no Ghia data configured for Re={sorted(unsupported)}")
    manifests = [
        run_case(
            arguments.repo.resolve(), arguments.moon.resolve(), arguments.output_root.resolve(),
            arguments.grid, checkpoints, reynolds, dt=arguments.dt, resume=arguments.resume,
        )
        for reynolds in arguments.re
    ]
    print(json.dumps(manifests, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
