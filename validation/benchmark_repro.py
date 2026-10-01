from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import statistics
import subprocess
import time
from pathlib import Path
from typing import Sequence


WORKLOADS: tuple[tuple[str, int], ...] = (
    ("chorin", 500),
    ("simple", 100),
    ("pcg", 500),
    ("mac", 500),
)


def summarize_samples(samples: Sequence[float]) -> dict[str, float | int]:
    if len(samples) < 2:
        raise ValueError("at least two measured samples are required")
    mean = statistics.fmean(samples)
    sample_std = statistics.stdev(samples)
    return {
        "n": len(samples),
        "mean_ms": mean,
        "sample_std_ms": sample_std,
        "median_ms": statistics.median(samples),
        "min_ms": min(samples),
        "max_ms": max(samples),
        "cv_percent": 100.0 * sample_std / mean,
    }


def build_runtime_command(
    moonrun: Path, artifact: Path, solver: str, steps: int
) -> list[str]:
    return [
        str(moonrun), str(artifact), "--", "--format", "centerline",
        "--solver", solver, "--grid", "41", "--steps", str(steps),
    ]


def build_prebuild_command(moon: Path) -> list[str]:
    """Build the release WASM artifact without launching the solver."""
    return [str(moon), "build", "cmd/main", "--release", "--target", "wasm"]


def _run_once(command: Sequence[str], cwd: Path) -> float:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=cwd,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    if completed.returncode != 0:
        raise RuntimeError(f"benchmark command failed with exit code {completed.returncode}")
    return elapsed_ms


def run_benchmark(
    repo: Path,
    moon: Path,
    moonrun: Path,
    artifact: Path,
    output_dir: Path,
    *,
    warmups: int,
    repeats: int,
) -> dict[str, object]:
    if warmups < 1 or repeats < 2:
        raise ValueError("benchmark requires at least one warm-up and two measured repetitions")
    output_dir.mkdir(parents=True, exist_ok=True)
    build_command = build_prebuild_command(moon)
    build = subprocess.run(build_command, cwd=repo, capture_output=True, text=True, encoding="utf-8")
    (output_dir / "prebuild_stderr.log").write_text(build.stderr, encoding="utf-8")
    if build.returncode != 0 or not artifact.is_file():
        raise RuntimeError("release WASM prebuild failed")

    commands = {
        solver: build_runtime_command(moonrun, artifact, solver, steps)
        for solver, steps in WORKLOADS
    }
    for _ in range(warmups):
        for solver, _steps in WORKLOADS:
            _run_once(commands[solver], repo)

    raw_rows: list[dict[str, object]] = []
    solvers = [solver for solver, _steps in WORKLOADS]
    step_by_solver = dict(WORKLOADS)
    for round_index in range(repeats):
        offset = round_index % len(solvers)
        round_order = solvers[offset:] + solvers[:offset]
        for order_index, solver in enumerate(round_order):
            elapsed_ms = _run_once(commands[solver], repo)
            raw_rows.append({
                "round": round_index + 1,
                "order": order_index + 1,
                "solver": solver,
                "grid": "41x41",
                "steps_or_iterations": step_by_solver[solver],
                "elapsed_ms": elapsed_ms,
            })

    raw_csv = output_dir / "benchmark_samples.csv"
    with raw_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(raw_rows[0]))
        writer.writeheader()
        writer.writerows(raw_rows)
    summary = {
        solver: summarize_samples(tuple(
            float(row["elapsed_ms"]) for row in raw_rows if row["solver"] == solver
        ))
        for solver, _steps in WORKLOADS
    }
    moon_version = subprocess.run(
        [str(moon), "version"], cwd=repo, capture_output=True, text=True, encoding="utf-8"
    ).stdout.strip()
    payload: dict[str, object] = {
        "method": {
            "timed_scope": "prebuilt release WASM executed by moonrun; process startup, solver run, and centerline serialization; compilation excluded",
            "warmups_discarded_per_solver": warmups,
            "measured_repetitions_per_solver": repeats,
            "ordering": "rotating solver order per measured round",
            "grid": "41x41",
            "workloads": dict(WORKLOADS),
        },
        "environment": {
            "platform": platform.platform(),
            "processor": platform.processor(),
            "logical_cpu_count": os.cpu_count(),
            "python": platform.python_version(),
            "moon": moon_version,
            "moonrun": str(moonrun.resolve()),
            "artifact": str(artifact.resolve()),
        },
        "summary": summary,
        "raw_csv": str(raw_csv.resolve()),
    }
    (output_dir / "benchmark_summary.json").write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reproducible FlowLabLite release-WASM benchmark.")
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--moon", type=Path, required=True)
    parser.add_argument("--moonrun", type=Path, required=True)
    parser.add_argument(
        "--artifact",
        type=Path,
        default=Path("_build/wasm/release/build/cmd/main/main.wasm"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("validation/results/performance"))
    parser.add_argument("--warmups", type=int, default=3)
    parser.add_argument("--repeats", type=int, default=15)
    arguments = parser.parse_args(argv)
    payload = run_benchmark(
        arguments.repo.resolve(),
        arguments.moon.resolve(),
        arguments.moonrun.resolve(),
        (arguments.repo / arguments.artifact).resolve() if not arguments.artifact.is_absolute() else arguments.artifact,
        arguments.output_dir.resolve(),
        warmups=arguments.warmups,
        repeats=arguments.repeats,
    )
    print(json.dumps(payload, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
