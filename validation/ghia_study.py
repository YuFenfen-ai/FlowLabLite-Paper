from __future__ import annotations

import bisect
import csv
import hashlib
import json
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


@dataclass(frozen=True)
class CenterlineData:
    metadata: dict[str, str | int | float]
    u_coord: tuple[float, ...]
    u_value: tuple[float, ...]
    v_coord: tuple[float, ...]
    v_value: tuple[float, ...]


@dataclass(frozen=True)
class ErrorMetrics:
    rms: float
    relative_l2: float
    max_abs: float


_FLOAT_METADATA = frozenset(("Re", "h", "dt", "div_norm"))
_INTEGER_METADATA = frozenset(("step", "pressure_iters"))
_ENDPOINT_TOLERANCE = 1e-10
_SETTINGS_TOLERANCE = 1e-12
_EXPECTED_DT = 0.001
_EXPECTED_PRESSURE_ITERS = 50
_EXPECTED_SOLVER = "chorin"


def parse_centerline(path: str | Path, *, expected_re: int | None = None) -> CenterlineData:
    """Read and validate one self-describing FlowLabLite centerline CSV."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    metadata_line = next((line for line in lines if line.startswith("#")), None)
    if metadata_line is None:
        raise ValueError("centerline CSV is missing its metadata comment")
    metadata = _parse_metadata(metadata_line)
    header_index = lines.index(metadata_line) + 1
    while header_index < len(lines) and not lines[header_index].strip():
        header_index += 1
    if header_index == len(lines):
        raise ValueError("centerline CSV is missing its data header")

    components: dict[str, list[tuple[float, float]]] = {"u_vcl": [], "v_hcl": []}
    for row in csv.DictReader(lines[header_index:]):
        component = row.get("component")
        if component not in components:
            raise ValueError(f"unknown centerline component: {component!r}")
        try:
            coordinate = float(row["coord"])
            value = float(row["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("centerline rows must contain numeric coord and value fields") from exc
        if not math.isfinite(coordinate) or not math.isfinite(value):
            raise ValueError("centerline coordinates and values must be finite")
        components[component].append((coordinate, value))

    u_coord, u_value = _validate_component("u_vcl", components["u_vcl"], (0.0, 1.0))
    v_coord, v_value = _validate_component("v_hcl", components["v_hcl"], (0.0, 0.0))
    data = CenterlineData(metadata, u_coord, u_value, v_coord, v_value)
    _validate_study_settings(data, expected_re=expected_re)
    return data


def interpolate(
    source_coord: Sequence[float], source_value: Sequence[float], target_coord: Sequence[float]
) -> tuple[float, ...]:
    """Linearly interpolate samples, rejecting targets outside the sampled range."""
    _validate_profile(source_coord, source_value, "source")
    interpolated: list[float] = []
    lower, upper = source_coord[0], source_coord[-1]
    for target in target_coord:
        if not math.isfinite(target):
            raise ValueError("target coordinates must be finite")
        if target < lower or target > upper:
            raise ValueError(f"target coordinate {target} is outside [{lower}, {upper}]")
        right = bisect.bisect_left(source_coord, target)
        if right < len(source_coord) and source_coord[right] == target:
            interpolated.append(source_value[right])
            continue
        left = right - 1
        fraction = (target - source_coord[left]) / (source_coord[right] - source_coord[left])
        interpolated.append(source_value[left] * (1.0 - fraction) + source_value[right] * fraction)
    return tuple(interpolated)


def compute_errors(
    solver_coord: Sequence[float],
    solver_value: Sequence[float],
    reference_coord: Sequence[float],
    reference_value: Sequence[float],
) -> ErrorMetrics:
    """Calculate RMS, relative L2, and maximum absolute profile error."""
    _validate_profile(reference_coord, reference_value, "reference")
    predicted = interpolate(solver_coord, solver_value, reference_coord)
    errors = [actual - expected for actual, expected in zip(predicted, reference_value, strict=True)]
    squared_error = sum(error * error for error in errors)
    denominator = math.sqrt(sum(value * value for value in reference_value))
    if denominator == 0.0:
        raise ValueError("relative L2 is undefined for an all-zero reference profile")
    return ErrorMetrics(
        math.sqrt(squared_error / len(errors)),
        math.sqrt(squared_error) / denominator,
        max(abs(error) for error in errors),
    )


def profile_change(
    older_coord: Sequence[float],
    older_value: Sequence[float],
    newer_coord: Sequence[float],
    newer_value: Sequence[float],
) -> float:
    """Return the relative L2 change after sampling the newer profile on older points."""
    _validate_profile(older_coord, older_value, "older")
    newer_at_older = interpolate(newer_coord, newer_value, older_coord)
    numerator = math.sqrt(
        sum((newer - older) ** 2 for older, newer in zip(older_value, newer_at_older, strict=True))
    )
    denominator = math.sqrt(sum(value * value for value in newer_at_older))
    if denominator == 0.0:
        raise ValueError("profile change is undefined for an all-zero newer profile")
    return numerator / denominator


def _parse_metadata(line: str) -> dict[str, str | int | float]:
    metadata: dict[str, str | int | float] = {}
    for segment in line[1:].split("|"):
        if "=" not in segment:
            continue
        key, raw_value = (part.strip() for part in segment.split("=", 1))
        if key in _FLOAT_METADATA:
            value = float(raw_value)
            if not math.isfinite(value):
                raise ValueError(f"metadata {key} must be finite")
            metadata[key] = value
        elif key in _INTEGER_METADATA:
            metadata[key] = int(raw_value)
        else:
            metadata[key] = raw_value
    return metadata


def _validate_component(
    name: str, samples: list[tuple[float, float]], expected_endpoints: tuple[float, float]
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if not samples:
        raise ValueError(f"centerline CSV is missing {name} samples")
    coordinates, values = zip(*samples, strict=True)
    _validate_profile(coordinates, values, name)
    if (
        abs(values[0] - expected_endpoints[0]) > _ENDPOINT_TOLERANCE
        or abs(values[-1] - expected_endpoints[1]) > _ENDPOINT_TOLERANCE
    ):
        raise ValueError(f"{name} does not satisfy its boundary endpoint values")
    return tuple(coordinates), tuple(values)


def _validate_profile(coord: Sequence[float], value: Sequence[float], label: str) -> None:
    if len(coord) != len(value) or len(coord) < 2:
        raise ValueError(f"{label} profile must have at least two equally sized samples")
    if not all(math.isfinite(item) for item in (*coord, *value)):
        raise ValueError(f"{label} profile contains a non-finite value")
    if any(right <= left for left, right in zip(coord, coord[1:])):
        raise ValueError(f"{label} coordinates must be strictly increasing")


def _validate_study_settings(data: CenterlineData, *, expected_re: int | None) -> None:
    solver = _required_metadata(data, "solver", str)
    reynolds = _required_metadata(data, "Re", float)
    grid = _required_metadata(data, "grid", str)
    spacing = _required_metadata(data, "h", float)
    timestep = _required_metadata(data, "dt", float)
    pressure_iters = _required_metadata(data, "pressure_iters", int)
    step = _required_metadata(data, "step", int)
    divergence = _required_metadata(data, "div_norm", float)
    if solver != _EXPECTED_SOLVER:
        raise ValueError(f"validation requires solver={_EXPECTED_SOLVER}")
    if reynolds <= 0.0 or (expected_re is not None and not math.isclose(reynolds, expected_re)):
        raise ValueError(f"centerline Re={reynolds} does not match expected Re={expected_re}")
    if not math.isclose(timestep, _EXPECTED_DT, rel_tol=0.0, abs_tol=_SETTINGS_TOLERANCE):
        raise ValueError(f"validation requires dt={_EXPECTED_DT}")
    if pressure_iters != _EXPECTED_PRESSURE_ITERS:
        raise ValueError(f"validation requires pressure_iters={_EXPECTED_PRESSURE_ITERS}")
    if step <= 0 or divergence < 0.0:
        raise ValueError("step must be positive and div_norm must be non-negative")
    try:
        nx_text, ny_text = grid.lower().split("x", 1)
        nx, ny = int(nx_text), int(ny_text)
    except (ValueError, AttributeError) as exc:
        raise ValueError("grid metadata must use NxN notation") from exc
    if nx != ny or nx < 3:
        raise ValueError("validation requires a square grid with at least three nodes")
    expected_h = 2.0 / (nx - 1)
    if not math.isclose(spacing, expected_h, rel_tol=0.0, abs_tol=_SETTINGS_TOLERANCE):
        raise ValueError("grid metadata and h are inconsistent with the 0..2 domain")
    if len(data.u_coord) != nx or len(data.v_coord) != nx:
        raise ValueError("centerline sample counts do not match grid metadata")
    for label, coordinates in (("u", data.u_coord), ("v", data.v_coord)):
        if not (
            math.isclose(coordinates[0], 0.0, abs_tol=_ENDPOINT_TOLERANCE)
            and math.isclose(coordinates[-1], 2.0, abs_tol=_ENDPOINT_TOLERANCE)
        ):
            raise ValueError(f"{label} centerline must span the physical domain [0,2]")


_U_COORD = tuple(2.0 * value for value in (
    0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000,
    0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000,
))
_V_COORD = tuple(2.0 * value for value in (
    0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000,
    0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000,
))


GHIA_REFERENCE: dict[int, dict[str, tuple[float, ...]]] = {
    100: {
        "u_coord": _U_COORD,
        "u_value": (0.0, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
                    -0.21090, -0.20581, -0.13641, 0.00332, 0.23151, 0.68717, 0.73722,
                    0.78871, 0.84123, 1.0),
        "v_coord": _V_COORD,
        "v_value": (0.0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077, 0.17507,
                    0.17527, 0.05454, -0.24533, -0.22445, -0.16914, -0.10313,
                    -0.08864, -0.07391, -0.05906, 0.0),
    },
    400: {
        "u_coord": _U_COORD,
        "u_value": (0.0, -0.08186, -0.09266, -0.10338, -0.14612, -0.24299, -0.32726,
                    -0.17119, -0.11477, 0.02135, 0.16256, 0.29093, 0.55892, 0.61756,
                    0.68439, 0.75837, 1.0),
        "v_coord": _V_COORD,
        "v_value": (0.0, 0.18360, 0.19713, 0.20920, 0.22965, 0.28124, 0.30203,
                    0.30174, 0.05188, -0.38598, -0.44993, -0.23827, -0.22847,
                    -0.19254, -0.15663, -0.12146, 0.0),
    },
    1000: {
        "u_coord": _U_COORD,
        "u_value": (0.0, -0.18109, -0.20196, -0.22220, -0.29730, -0.38289, -0.27805,
                    -0.10648, -0.06080, 0.05702, 0.18719, 0.33304, 0.46604, 0.51117,
                    0.57492, 0.65928, 1.0),
        "v_coord": _V_COORD,
        "v_value": (0.0, 0.27485, 0.29012, 0.30353, 0.32627, 0.37095, 0.33075,
                    0.32235, 0.02526, -0.31966, -0.42665, -0.51550, -0.39188,
                    -0.33714, -0.27669, -0.21388, 0.0),
    },
}

GHIA_RE100_U_COORD = GHIA_REFERENCE[100]["u_coord"]
GHIA_RE100_U_VALUE = GHIA_REFERENCE[100]["u_value"]
GHIA_RE100_V_COORD = GHIA_REFERENCE[100]["v_coord"]
GHIA_RE100_V_VALUE = GHIA_REFERENCE[100]["v_value"]

_SUMMARY_FIELDS = (
    "re", "grid", "steps", "time", "div_norm", "u_rms", "u_relative_l2",
    "u_max_abs", "u_max_abs_percent", "u_max_pass_5pct", "v_rms", "v_relative_l2",
    "v_max_abs", "v_max_abs_percent", "v_max_pass_5pct", "u_min", "u_min_y_over_l",
    "u_change", "v_change", "time_converged",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _is_time_converged(
    u_change: float,
    v_change: float,
    previous_divergence: float,
    divergence: float,
) -> bool:
    roundoff = max(1e-15, abs(previous_divergence) * 1e-12)
    return u_change < 0.005 and v_change < 0.005 and divergence <= previous_divergence + roundoff


def _settings_signature(data: CenterlineData) -> tuple[object, ...]:
    return (
        data.metadata["solver"], data.metadata["Re"], data.metadata["grid"],
        data.metadata["h"], data.metadata["dt"], data.metadata["pressure_iters"],
        data.u_coord[0], data.u_coord[-1], data.v_coord[0], data.v_coord[-1],
    )


def write_summary(
    cases: Sequence[str | Path | CenterlineData], output_csv: str | Path, output_json: str | Path
) -> dict[str, object]:
    """Write per-case Ghia metrics and strict adjacent-checkpoint convergence evidence."""
    if not cases:
        raise ValueError("at least one centerline case is required")
    csv_path = Path(output_csv).resolve()
    json_path = Path(output_json).resolve()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)

    normalized: list[tuple[CenterlineData, Path | None, str | None]] = []
    expected_re: int | None = None
    for case in cases:
        if isinstance(case, CenterlineData):
            data, source_path, source_hash = case, None, None
        else:
            source_path = Path(case).resolve()
            source_hash = _sha256(source_path)
            data = parse_centerline(source_path)
        case_re = int(round(float(data.metadata["Re"])))
        if case_re not in GHIA_REFERENCE:
            raise ValueError(f"no Ghia reference data are configured for Re={case_re}")
        if expected_re is None:
            expected_re = case_re
        elif case_re != expected_re:
            raise ValueError("all time-convergence checkpoints must use the same Reynolds number")
        _validate_study_settings(data, expected_re=expected_re)
        normalized.append((data, source_path, source_hash))
    assert expected_re is not None

    reference = GHIA_REFERENCE[expected_re]
    rows: list[dict[str, object]] = []
    first_signature: tuple[object, ...] | None = None
    previous: CenterlineData | None = None
    previous_steps: int | None = None
    previous_divergence: float | None = None
    for data, source_path, source_hash in normalized:
        signature = _settings_signature(data)
        if first_signature is None:
            first_signature = signature
        elif signature != first_signature:
            raise ValueError("time-convergence checkpoints must use identical solver settings and grid")
        steps = int(data.metadata["step"])
        if previous_steps is not None and steps <= previous_steps:
            raise ValueError("time-convergence checkpoints must have strictly increasing step counts")
        divergence = float(data.metadata["div_norm"])
        u_metrics = compute_errors(data.u_coord, data.u_value, reference["u_coord"], reference["u_value"])
        v_metrics = compute_errors(data.v_coord, data.v_value, reference["v_coord"], reference["v_value"])
        u_change = None if previous is None else profile_change(
            previous.u_coord, previous.u_value, data.u_coord, data.u_value
        )
        v_change = None if previous is None else profile_change(
            previous.v_coord, previous.v_value, data.v_coord, data.v_value
        )
        time_converged = bool(
            u_change is not None and v_change is not None and previous_divergence is not None
            and _is_time_converged(u_change, v_change, previous_divergence, divergence)
        )
        minimum_index = min(range(len(data.u_value)), key=data.u_value.__getitem__)
        row: dict[str, object] = {
            "re": expected_re,
            "grid": data.metadata["grid"],
            "steps": steps,
            "time": steps * float(data.metadata["dt"]),
            "div_norm": divergence,
            "u_rms": u_metrics.rms,
            "u_relative_l2": u_metrics.relative_l2,
            "u_max_abs": u_metrics.max_abs,
            "u_max_abs_percent": 100.0 * u_metrics.max_abs,
            "u_max_pass_5pct": u_metrics.max_abs <= 0.05,
            "v_rms": v_metrics.rms,
            "v_relative_l2": v_metrics.relative_l2,
            "v_max_abs": v_metrics.max_abs,
            "v_max_abs_percent": 100.0 * v_metrics.max_abs,
            "v_max_pass_5pct": v_metrics.max_abs <= 0.05,
            "u_min": data.u_value[minimum_index],
            "u_min_y_over_l": data.u_coord[minimum_index] / 2.0,
            "u_change": u_change,
            "v_change": v_change,
            "time_converged": time_converged,
        }
        if source_path is not None and source_hash is not None:
            row["source_csv"] = os.path.relpath(source_path, json_path.parent)
            row["source_sha256"] = source_hash
        rows.append(row)
        previous, previous_steps, previous_divergence = data, steps, divergence

    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_SUMMARY_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    selected_case = rows[-1]
    selected_time_step = selected_case["steps"] if selected_case["time_converged"] else None
    payload: dict[str, object] = {
        "re": expected_re,
        "reference": "Ghia, Ghia and Shin (1982), Tables I and II, 129x129",
        "cases": rows,
        "selected_time_step": selected_time_step,
        "selected_steps": selected_time_step,
        "time_converged": selected_case["time_converged"],
        "selected_case": selected_case,
    }
    json_path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return payload


def _load_re20_case(path: str | Path) -> tuple[tuple[float, ...], tuple[float, ...], str]:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    config = payload.get("config", {})
    required = {"re": 20, "nx": 41, "ny": 41, "nt": 500}
    for key, expected in required.items():
        if int(config.get(key, -1)) != expected:
            raise ValueError(f"Re=20 JSON requires config.{key}={expected}")
    if not math.isclose(float(config.get("dt", math.nan)), 0.001, abs_tol=_SETTINGS_TOLERANCE):
        raise ValueError("Re=20 JSON requires config.dt=0.001")
    nx = int(config["nx"])
    samples = [point for point in payload.get("grid", ()) if int(point["j"]) == nx // 2]
    if not samples:
        raise ValueError("Re=20 results contain no vertical-centerline samples")
    samples.sort(key=lambda point: int(point["i"]))
    coordinates = tuple(float(point["y"]) for point in samples)
    values = tuple(float(point["u"]) for point in samples)
    _validate_profile(coordinates, values, "Re=20")
    label = f"FlowLabLite Re={config['re']} ({config['nx']}\u00d7{config['ny']}, {config['nt']} steps)"
    return coordinates, values, label


def load_re20_centerline(path: str | Path) -> tuple[tuple[float, ...], tuple[float, ...]]:
    coordinates, values, _ = _load_re20_case(path)
    return coordinates, values


def _load_selected_case(path: str | Path) -> tuple[CenterlineData, Path, str, str, int | None]:
    candidate = Path(path).resolve()
    if candidate.suffix.lower() != ".json":
        source_hash = _sha256(candidate)
        return parse_centerline(candidate), candidate, source_hash, "not assessed", None
    summary = json.loads(candidate.read_text(encoding="utf-8"))
    selected_case = summary.get("selected_case")
    if not isinstance(selected_case, dict):
        raise ValueError("summary is missing selected_case")
    top_converged = bool(summary.get("time_converged"))
    row_converged = bool(selected_case.get("time_converged"))
    selected_time_step = summary.get("selected_time_step")
    row_steps = selected_case.get("steps")
    if top_converged != row_converged:
        raise ValueError("summary convergence flags are inconsistent")
    if top_converged and selected_time_step != row_steps:
        raise ValueError("summary selected time step is inconsistent with selected_case")
    if not top_converged and selected_time_step is not None:
        raise ValueError("a non-converged summary cannot select a time step")
    source_path = Path(str(selected_case.get("source_csv", "")))
    if not source_path.is_absolute():
        source_path = (candidate.parent / source_path).resolve()
    expected_hash = selected_case.get("source_sha256")
    actual_hash = _sha256(source_path)
    if expected_hash != actual_hash:
        raise ValueError("selected centerline checksum does not match the summary")
    expected_re = int(summary.get("re"))
    data = parse_centerline(source_path, expected_re=expected_re)
    if data.metadata["grid"] != selected_case.get("grid") or data.metadata["step"] != row_steps:
        raise ValueError("selected centerline metadata does not match the summary")
    status = "converged" if top_converged else "not converged"
    return data, source_path, actual_hash, status, selected_time_step


def plot_comparison(
    re20_json: str | Path,
    selected_re100_csv: str | Path,
    output_png: str | Path,
    output_pdf: str | Path,
) -> dict[str, object]:
    """Render Re=20 context plus one same-Re FlowLabLite/Ghia quantitative comparison."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as pyplot

    re20_coord, re20_value, re20_label = _load_re20_case(re20_json)
    selected, source_path, source_hash, convergence_status, selected_time_step = _load_selected_case(
        selected_re100_csv
    )
    reynolds = int(round(float(selected.metadata["Re"])))
    reference = GHIA_REFERENCE[reynolds]
    u_metrics = compute_errors(
        selected.u_coord, selected.u_value, reference["u_coord"], reference["u_value"]
    )
    labels = (
        re20_label,
        f"FlowLabLite Re={reynolds} ({selected.metadata['grid']}, {selected.metadata['step']} steps)",
        f"Ghia Re={reynolds}",
    )
    figure, axis = pyplot.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    try:
        axis.plot(re20_value, tuple(coord / 2.0 for coord in re20_coord), "o-", color="#2563eb", label=labels[0])
        axis.plot(selected.u_value, tuple(coord / 2.0 for coord in selected.u_coord), "^-", color="#16a34a", label=labels[1])
        axis.plot(reference["u_value"], tuple(coord / 2.0 for coord in reference["u_coord"]), "s--", color="#dc2626", label=labels[2])
        axis.axvline(0.0, color="0.45", linewidth=0.8)
        axis.set(xlabel="u/U_lid", ylabel="y/L", ylim=(0.0, 1.0))
        axis.grid(True, color="0.88", linewidth=0.8)
        axis.legend(loc="lower right")
        settings = "\n".join((
            f"grid: {selected.metadata['grid']}",
            f"steps: {selected.metadata['step']}",
            f"dt={selected.metadata['dt']}, nu={2.0 / reynolds:g}",
            f"pressure: GS{selected.metadata['pressure_iters']}",
            f"u relative L2: {u_metrics.relative_l2:.4g}",
            f"u max abs: {u_metrics.max_abs:.4g}",
            f"time convergence: {convergence_status}",
        ))
        axis.text(0.02, 0.98, settings, transform=axis.transAxes, va="top", fontsize=8,
                  bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.9, "edgecolor": "0.65"})
        png_path, pdf_path = Path(output_png), Path(output_pdf)
        png_path.parent.mkdir(parents=True, exist_ok=True)
        pdf_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(png_path, dpi=200)
        figure.savefig(pdf_path)
    finally:
        pyplot.close(figure)
    return {
        "labels": labels,
        "source_csv": str(source_path),
        "source_sha256": source_hash,
        "selected_time_step": selected_time_step,
        "re": reynolds,
        "u_y_over_l": tuple(coord / 2.0 for coord in selected.u_coord),
    }


def _required_metadata(data: CenterlineData, key: str, expected_type: type[object]) -> object:
    value = data.metadata.get(key)
    if not isinstance(value, expected_type) or isinstance(value, bool):
        raise ValueError(f"centerline metadata {key!r} is missing or has the wrong type")
    return value


def main(argv: Sequence[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Analyze FlowLabLite Ghia centerline studies.")
    commands = parser.add_subparsers(dest="command", required=True)
    summarize = commands.add_parser("summarize", help="write Ghia and convergence metrics")
    summarize.add_argument("--inputs", nargs="+", required=True)
    summarize.add_argument("--csv", required=True)
    summarize.add_argument("--json", required=True)
    plot = commands.add_parser("plot", help="render a same-Re quantitative comparison")
    plot.add_argument("--re20", required=True)
    plot.add_argument("--re100", required=True, help="selected centerline CSV or summary JSON")
    plot.add_argument("--png", required=True)
    plot.add_argument("--pdf", required=True)
    arguments = parser.parse_args(argv)
    if arguments.command == "summarize":
        manifest = write_summary(arguments.inputs, arguments.csv, arguments.json)
    else:
        manifest = plot_comparison(arguments.re20, arguments.re100, arguments.png, arguments.pdf)
    print(json.dumps(manifest, ensure_ascii=False))
    return 0
