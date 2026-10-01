from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from pathlib import Path
from typing import Sequence

from validation.ghia_study import GHIA_REFERENCE
from validation.projection_tvd_study import (
    parse_projection_centerline,
    write_projection_summary,
)
from validation.run_ghia_study import split_checkpoint_stream


def build_solver_command(
    moon: Path,
    target: str,
    grid: int,
    checkpoints: Sequence[int],
    reynolds: int,
    scheme: str,
    timestep: float,
    pressure_tolerance: float,
    pressure_maximum_iterations: int,
    pressure_omega: float,
) -> list[str]:
    return [
        str(moon),
        "run",
        "cmd/main",
        "--release",
        "--target",
        target,
        "--",
        "--format",
        "centerline",
        "--solver",
        "projection-tvd",
        "--grid",
        str(grid),
        "--checkpoints",
        ",".join(str(step) for step in checkpoints),
        "--re",
        str(reynolds),
        "--scheme",
        scheme,
        "--dt",
        str(timestep),
        "--pressure-tol",
        str(pressure_tolerance),
        "--pressure-max-iter",
        str(pressure_maximum_iterations),
        "--pressure-omega",
        str(pressure_omega),
    ]


def run_case(
    repo: Path,
    moon: Path,
    target: str,
    output_root: Path,
    checkpoints: Sequence[int],
    reynolds: int,
    scheme: str,
    timestep: float,
    pressure_tolerance: float,
    pressure_maximum_iterations: int,
    pressure_omega: float,
    *,
    resume: bool,
) -> dict[str, object]:
    case_dir = output_root / f"re{reynolds}" / scheme
    case_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_paths = {step: case_dir / f"centerline_s{step}.csv" for step in checkpoints}
    command = build_solver_command(
        moon,
        target,
        129,
        checkpoints,
        reynolds,
        scheme,
        timestep,
        pressure_tolerance,
        pressure_maximum_iterations,
        pressure_omega,
    )
    can_resume = resume and all(path.is_file() for path in checkpoint_paths.values())
    elapsed_seconds: float | None = None
    if can_resume:
        for path in checkpoint_paths.values():
            parse_projection_centerline(path, expected_re=reynolds)
    else:
        started = time.perf_counter()
        completed = subprocess.run(
            command,
            cwd=repo,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        elapsed_seconds = time.perf_counter() - started
        (case_dir / "solver_stdout.log").write_text(completed.stdout, encoding="utf-8")
        (case_dir / "solver_stderr.log").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(
                f"projection-tvd failed for Re={reynolds}; see {case_dir / 'solver_stderr.log'}"
            )
        blocks = split_checkpoint_stream(completed.stdout)
        if tuple(blocks) != tuple(checkpoints):
            raise ValueError(f"expected checkpoints {tuple(checkpoints)}, received {tuple(blocks)}")
        for step, text in blocks.items():
            path = checkpoint_paths[step]
            path.write_text(text, encoding="utf-8")
            parse_projection_centerline(path, expected_re=reynolds)

    summary_csv = case_dir / "time_convergence.csv"
    summary_json = case_dir / "time_convergence.json"
    summary = write_projection_summary(
        tuple(checkpoint_paths.values()),
        summary_csv,
        summary_json,
        case_dir / "pointwise_errors.csv",
    )
    manifest: dict[str, object] = {
        "re": reynolds,
        "grid": 129,
        "target": target,
        "scheme": scheme,
        "checkpoints": list(checkpoints),
        "dt": timestep,
        "pressure_tol": pressure_tolerance,
        "pressure_max_iter": pressure_maximum_iterations,
        "pressure_omega": pressure_omega,
        "command": command,
        "elapsed_seconds": elapsed_seconds,
        "resumed": can_resume,
        "python": sys.version,
        "platform": platform.platform(),
        "strict_validation_pass": summary["strict_validation_pass"],
        "summary_json": str(summary_json.resolve()),
    }
    (case_dir / "run_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run strict 129x129 Ghia validation for the fifth projection-TVD solver."
    )
    parser.add_argument("--moon", type=Path, required=True)
    parser.add_argument("--target", choices=("wasm", "wasm-gc", "js"), default="wasm")
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--output-root", type=Path, default=Path("validation/results/projection_tvd")
    )
    parser.add_argument("--checkpoints", type=int, nargs="+", default=(500, 1000, 2000))
    parser.add_argument("--re", type=int, nargs="+", default=(400, 1000))
    parser.add_argument("--scheme", choices=("upwind", "tvd-vanleer"), default="tvd-vanleer")
    parser.add_argument("--dt", type=float, default=0.001)
    parser.add_argument("--pressure-tol", type=float, default=1e-6)
    parser.add_argument("--pressure-max-iter", type=int, default=12000)
    parser.add_argument("--pressure-omega", type=float, default=1.95)
    parser.add_argument("--resume", action="store_true")
    arguments = parser.parse_args(argv)
    checkpoints = tuple(arguments.checkpoints)
    if any(step <= 0 for step in checkpoints) or any(
        right <= left for left, right in zip(checkpoints, checkpoints[1:])
    ):
        parser.error("--checkpoints must be strictly increasing positive integers")
    unsupported = set(arguments.re) - set(GHIA_REFERENCE)
    if unsupported:
        parser.error(f"no Ghia data configured for Re={sorted(unsupported)}")
    if arguments.dt <= 0 or arguments.pressure_tol <= 0:
        parser.error("--dt and --pressure-tol must be positive")
    if arguments.pressure_max_iter <= 0 or not 0 < arguments.pressure_omega < 2:
        parser.error("pressure iteration settings are invalid")
    manifests = [
        run_case(
            arguments.repo.resolve(),
            arguments.moon.resolve(),
            arguments.target,
            arguments.output_root.resolve(),
            checkpoints,
            reynolds,
            arguments.scheme,
            arguments.dt,
            arguments.pressure_tol,
            arguments.pressure_max_iter,
            arguments.pressure_omega,
            resume=arguments.resume,
        )
        for reynolds in arguments.re
    ]
    print(json.dumps(manifests, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
