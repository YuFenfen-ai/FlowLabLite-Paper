from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


_SOLVERS = (
    ("chorin", "Chorin-GS"),
    ("simple", "SIMPLE"),
    ("pcg", "Chorin-PCG"),
    ("mac", "MAC"),
)
_FIELDS = (
    "solver",
    "grid",
    "steps_or_iterations",
    "n",
    "mean_ms",
    "sample_std_ms",
    "median_ms",
    "min_ms",
    "max_ms",
    "cv_percent",
)
_NUMERIC_FIELDS = (
    "mean_ms",
    "sample_std_ms",
    "median_ms",
    "min_ms",
    "max_ms",
    "cv_percent",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate_and_rows(payload: dict[str, object]) -> tuple[dict[str, object], list[dict[str, object]]]:
    method = payload.get("method")
    summary = payload.get("summary")
    if not isinstance(method, dict) or not isinstance(summary, dict):
        raise ValueError("benchmark summary requires method and summary objects")
    warmups = int(method.get("warmups_discarded_per_solver", 0))
    repeats = int(method.get("measured_repetitions_per_solver", 0))
    if warmups < 1:
        raise ValueError("at least one discarded warm-up is required")
    if repeats < 2:
        raise ValueError("at least two measured repetitions are required")
    if method.get("grid") != "41x41":
        raise ValueError("performance report requires the 41x41 workload")
    workloads = method.get("workloads")
    expected_workloads = {"chorin": 500, "simple": 100, "pcg": 500, "mac": 500}
    if workloads != expected_workloads:
        raise ValueError("performance workloads do not match the thesis configuration")

    rows: list[dict[str, object]] = []
    for key, label in _SOLVERS:
        metrics = summary.get(key)
        if not isinstance(metrics, dict):
            raise ValueError(f"benchmark summary is missing {key}")
        if int(metrics.get("n", -1)) != repeats:
            raise ValueError(f"{key} sample count does not match measured repetitions")
        numeric = {name: float(metrics.get(name, math.nan)) for name in _NUMERIC_FIELDS}
        if not all(math.isfinite(value) and value >= 0.0 for value in numeric.values()):
            raise ValueError(f"{key} contains an invalid performance statistic")
        if not numeric["min_ms"] <= numeric["median_ms"] <= numeric["max_ms"]:
            raise ValueError(f"{key} min/median/max are inconsistent")
        rows.append(
            {
                "solver": label,
                "grid": method["grid"],
                "steps_or_iterations": expected_workloads[key],
                "n": repeats,
                **numeric,
            }
        )
    return method, rows


def render_performance_report(
    summary_path: str | Path,
    output_csv: str | Path,
    output_png: str | Path,
    output_pdf: str | Path,
    output_manifest: str | Path,
) -> dict[str, object]:
    source = Path(summary_path).resolve()
    payload = json.loads(source.read_text(encoding="utf-8"))
    method, rows = _validate_and_rows(payload)

    csv_path = Path(output_csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot

    figure, axis = pyplot.subplots(figsize=(7.2, 4.0), constrained_layout=True)
    try:
        labels = [str(row["solver"]) for row in rows]
        means = [float(row["mean_ms"]) / 1000.0 for row in rows]
        errors = [float(row["sample_std_ms"]) / 1000.0 for row in rows]
        ymax_seconds = 1.25 * max(mean + error for mean, error in zip(means, errors, strict=True))
        bars = axis.bar(
            labels,
            means,
            yerr=errors,
            capsize=4,
            color=("#2563eb", "#16a34a", "#f59e0b", "#dc2626"),
            edgecolor="0.25",
            linewidth=0.6,
        )
        axis.set_ylabel("Wall-clock time (s)")
        axis.set_title("Release WASM solver runtime (mean ± sample SD, n=15)")
        axis.set_ylim(0.0, ymax_seconds)
        axis.grid(axis="y", color="0.88", linewidth=0.8)
        axis.set_axisbelow(True)
        for bar, row, error in zip(bars, rows, errors, strict=True):
            axis.text(
                bar.get_x() + bar.get_width() / 2.0,
                bar.get_height() + error + 0.015 * ymax_seconds,
                f"median={float(row['median_ms']) / 1000.0:.2f}s\nCV={float(row['cv_percent']):.1f}%",
                ha="center",
                va="bottom",
                fontsize=8,
            )
        axis.text(
            0.0,
            -0.22,
            "41×41; Chorin/PCG/MAC: 500 steps; SIMPLE: 100 iterations; "
            "3 warm-ups discarded; rotating order; compilation excluded",
            transform=axis.transAxes,
            fontsize=8,
        )
        png_path, pdf_path = Path(output_png), Path(output_pdf)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(png_path, dpi=220)
        figure.savefig(pdf_path)
    finally:
        pyplot.close(figure)

    manifest: dict[str, object] = {
        "source_summary": str(source),
        "source_sha256": _sha256(source),
        "method": method,
        "rows": rows,
        "plot": {"ymax_seconds": ymax_seconds},
        "csv": str(csv_path.resolve()),
        "png": str(Path(output_png).resolve()),
        "pdf": str(Path(output_pdf).resolve()),
    }
    manifest_path = Path(output_manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest
