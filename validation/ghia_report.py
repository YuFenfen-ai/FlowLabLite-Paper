from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Sequence

from validation.ghia_study import GHIA_REFERENCE, _load_selected_case, compute_errors


_REPORT_FIELDS = (
    "re", "grid", "steps", "time", "time_converged", "div_norm",
    "u_rms", "u_relative_l2_percent", "u_max_abs", "u_max_abs_percent", "u_pass_5pct",
    "v_rms", "v_relative_l2_percent", "v_max_abs", "v_max_abs_percent", "v_pass_5pct",
    "source_csv", "source_sha256",
)


def render_report(
    summary_paths: Sequence[str | Path],
    output_csv: str | Path,
    output_png: str | Path,
    output_pdf: str | Path,
    output_manifest: str | Path,
) -> dict[str, object]:
    """Create the three-Re quantitative table and 3x2 u/v centerline figure."""
    if len(summary_paths) != 3:
        raise ValueError("the thesis report requires exactly Re=100, 400, and 1000 summaries")
    loaded: list[tuple[int, object, Path, str, str]] = []
    rows: list[dict[str, object]] = []
    for summary_path in summary_paths:
        summary_file = Path(summary_path).resolve()
        summary = json.loads(summary_file.read_text(encoding="utf-8"))
        data, source_path, source_hash, status, _selected_step = _load_selected_case(summary_file)
        reynolds = int(round(float(data.metadata["Re"])))
        if data.metadata["grid"] != "129x129":
            raise ValueError("the thesis quantitative table requires 129x129 results")
        reference = GHIA_REFERENCE[reynolds]
        u_metrics = compute_errors(data.u_coord, data.u_value, reference["u_coord"], reference["u_value"])
        v_metrics = compute_errors(data.v_coord, data.v_value, reference["v_coord"], reference["v_value"])
        row = {
            "re": reynolds,
            "grid": data.metadata["grid"],
            "steps": data.metadata["step"],
            "time": float(data.metadata["step"]) * float(data.metadata["dt"]),
            "time_converged": status == "converged",
            "div_norm": data.metadata["div_norm"],
            "u_rms": u_metrics.rms,
            "u_relative_l2_percent": 100.0 * u_metrics.relative_l2,
            "u_max_abs": u_metrics.max_abs,
            "u_max_abs_percent": 100.0 * u_metrics.max_abs,
            "u_pass_5pct": u_metrics.max_abs <= 0.05,
            "v_rms": v_metrics.rms,
            "v_relative_l2_percent": 100.0 * v_metrics.relative_l2,
            "v_max_abs": v_metrics.max_abs,
            "v_max_abs_percent": 100.0 * v_metrics.max_abs,
            "v_pass_5pct": v_metrics.max_abs <= 0.05,
            "source_csv": str(source_path),
            "source_sha256": source_hash,
        }
        rows.append(row)
        loaded.append((reynolds, data, source_path, source_hash, status))
    loaded.sort(key=lambda item: item[0])
    rows.sort(key=lambda item: int(item["re"]))
    if [row["re"] for row in rows] != [100, 400, 1000]:
        raise ValueError("summaries must cover Re=100, 400, and 1000 exactly once")

    csv_path = Path(output_csv)
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_REPORT_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot

    figure, axes = pyplot.subplots(2, 3, figsize=(11.0, 7.17))
    figure.subplots_adjust(left=0.07, right=0.98, bottom=0.10, top=0.86, wspace=0.34, hspace=0.42)
    try:
        for column_index, (reynolds, data, _source, _hash, status) in enumerate(loaded):
            reference = GHIA_REFERENCE[reynolds]
            u_axis, v_axis = axes[0, column_index], axes[1, column_index]
            u_axis.plot(data.u_value, tuple(value / 2.0 for value in data.u_coord), color="#2563eb", linewidth=1.6, label="FlowLabLite")
            u_axis.plot(reference["u_value"], tuple(value / 2.0 for value in reference["u_coord"]), "s", color="#dc2626", markersize=4, label="Ghia (1982)")
            v_axis.plot(tuple(value / 2.0 for value in data.v_coord), data.v_value, color="#16a34a", linewidth=1.6, label="FlowLabLite")
            v_axis.plot(tuple(value / 2.0 for value in reference["v_coord"]), reference["v_value"], "s", color="#dc2626", markersize=4, label="Ghia (1982)")
            u_axis.axvline(0.0, color="0.5", linewidth=0.7)
            v_axis.axhline(0.0, color="0.5", linewidth=0.7)
            u_axis.set(xlabel="u/U_lid", ylabel="y/L", title=f"Re={reynolds}: vertical centerline u")
            v_axis.set(xlabel="x/L", ylabel="v/U_lid", title=f"Re={reynolds}: horizontal centerline v")
            u_axis.set_ylim(0.0, 1.0)
            v_axis.set_xlim(0.0, 1.0)
            for axis in (u_axis, v_axis):
                axis.grid(True, color="0.88", linewidth=0.8)
                axis.legend(fontsize=8)
            row = rows[column_index]
            annotation = (
                f"{row['steps']} steps; rel-L2 u/v="
                f"{row['u_relative_l2_percent']:.2f}%/{row['v_relative_l2_percent']:.2f}%\n"
                f"max u/v={row['u_max_abs_percent']:.2f}%/{row['v_max_abs_percent']:.2f}%; {status}"
            )
            v_axis.text(0.02, 0.03, annotation, transform=v_axis.transAxes, fontsize=6.5, va="bottom", bbox={"facecolor": "white", "alpha": 0.72, "edgecolor": "none"})
        figure.suptitle(
            "FlowLabLite vs Ghia centerline velocities\n"
            "Chorin-GS; 129x129 uniform grid; L=2; U_lid=1; rho=1; dt=0.001; GS50",
            fontsize=12,
        )
        png_path, pdf_path = Path(output_png), Path(output_pdf)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(png_path, dpi=220)
        figure.savefig(pdf_path)
    finally:
        pyplot.close(figure)

    manifest: dict[str, object] = {
        "settings": {
            "solver": "Chorin-GS",
            "grid": "129x129",
            "domain": "2x2",
            "u_lid": 1.0,
            "rho": 1.0,
            "dt": 0.001,
            "pressure_iterations_per_step": 50,
            "reference": "Ghia, Ghia and Shin (1982), Tables I and II",
            "max_error_pass_rule": "100*max_abs/U_lid <= 5%",
        },
        "rows": rows,
        "combined_csv": str(csv_path.resolve()),
        "png": str(Path(output_png).resolve()),
        "pdf": str(Path(output_pdf).resolve()),
    }
    manifest_path = Path(output_manifest)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest
