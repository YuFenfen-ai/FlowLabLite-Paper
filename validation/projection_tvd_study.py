from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Sequence

from validation.ghia_study import (
    CenterlineData,
    GHIA_REFERENCE,
    compute_errors,
    interpolate,
    profile_change,
)


_FLOAT_FIELDS = frozenset(
    (
        "Re",
        "h",
        "dt",
        "pressure_tol",
        "pressure_omega",
        "pressure_residual",
        "div_norm",
        "cfl_max",
    )
)
_INTEGER_FIELDS = frozenset(("step", "pressure_max_iter", "pressure_iters"))
_REQUIRED_FIELDS = frozenset(
    (
        "solver",
        "scheme",
        "Re",
        "step",
        "grid",
        "h",
        "dt",
        "pressure_tol",
        "pressure_max_iter",
        "pressure_omega",
        "pressure_iters",
        "pressure_residual",
        "pressure_converged",
        "div_norm",
        "cfl_max",
        "finite",
    )
)
_ENDPOINT_TOLERANCE = 1e-10
_SETTINGS_TOLERANCE = 1e-12
_GRID = 129
_MAX_ABS_LIMIT = 0.05
_PROFILE_CHANGE_LIMIT = 0.005


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _parse_metadata(line: str) -> dict[str, str | int | float]:
    metadata: dict[str, str | int | float] = {}
    for segment in line[1:].split("|"):
        if "=" not in segment:
            continue
        key, raw_value = (part.strip() for part in segment.split("=", 1))
        try:
            if key in _FLOAT_FIELDS:
                value = float(raw_value)
                if not math.isfinite(value):
                    raise ValueError(f"metadata {key} must be finite")
                metadata[key] = value
            elif key in _INTEGER_FIELDS:
                metadata[key] = int(raw_value)
            else:
                metadata[key] = raw_value
        except ValueError as exc:
            raise ValueError(f"metadata {key} has an invalid value") from exc
    missing = sorted(_REQUIRED_FIELDS - metadata.keys())
    if missing:
        raise ValueError(f"centerline metadata are missing: {', '.join(missing)}")
    return metadata


def _validate_profile(
    label: str,
    samples: list[tuple[int, float, float]],
    expected_endpoints: tuple[float, float],
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    if len(samples) != _GRID:
        raise ValueError(f"{label} must contain exactly {_GRID} samples")
    indices, coordinates, values = zip(*samples, strict=True)
    if tuple(indices) != tuple(range(_GRID)):
        raise ValueError(f"{label} indices must be exactly 0..{_GRID - 1}")
    if not all(math.isfinite(item) for item in (*coordinates, *values)):
        raise ValueError(f"{label} contains a non-finite value")
    if any(right <= left for left, right in zip(coordinates, coordinates[1:])):
        raise ValueError(f"{label} coordinates must be strictly increasing")
    if not (
        math.isclose(coordinates[0], 0.0, abs_tol=_ENDPOINT_TOLERANCE)
        and math.isclose(coordinates[-1], 2.0, abs_tol=_ENDPOINT_TOLERANCE)
    ):
        raise ValueError(f"{label} must span the physical domain [0,2]")
    if not (
        math.isclose(values[0], expected_endpoints[0], abs_tol=_ENDPOINT_TOLERANCE)
        and math.isclose(values[-1], expected_endpoints[1], abs_tol=_ENDPOINT_TOLERANCE)
    ):
        raise ValueError(f"{label} does not satisfy its boundary endpoint values")
    return tuple(coordinates), tuple(values)


def _require_number(metadata: dict[str, object], key: str, expected_type: type) -> object:
    value = metadata.get(key)
    if not isinstance(value, expected_type) or isinstance(value, bool):
        raise ValueError(f"metadata {key} has the wrong type")
    return value


def _validate_settings(
    data: CenterlineData, *, expected_re: int | None = None
) -> None:
    metadata = data.metadata
    solver = _require_number(metadata, "solver", str)
    scheme = _require_number(metadata, "scheme", str)
    reynolds = float(_require_number(metadata, "Re", float))
    grid = _require_number(metadata, "grid", str)
    spacing = float(_require_number(metadata, "h", float))
    timestep = float(_require_number(metadata, "dt", float))
    tolerance = float(_require_number(metadata, "pressure_tol", float))
    maximum_iterations = int(_require_number(metadata, "pressure_max_iter", int))
    omega = float(_require_number(metadata, "pressure_omega", float))
    iterations = int(_require_number(metadata, "pressure_iters", int))
    residual = float(_require_number(metadata, "pressure_residual", float))
    step = int(_require_number(metadata, "step", int))
    divergence = float(_require_number(metadata, "div_norm", float))
    cfl = float(_require_number(metadata, "cfl_max", float))

    if solver != "projection-tvd":
        raise ValueError("validation requires solver=projection-tvd")
    if scheme not in ("upwind", "tvd-vanleer"):
        raise ValueError("validation requires scheme=upwind or tvd-vanleer")
    if reynolds <= 0.0 or (
        expected_re is not None and not math.isclose(reynolds, expected_re)
    ):
        raise ValueError(f"centerline Re={reynolds} does not match expected Re={expected_re}")
    if grid != f"{_GRID}x{_GRID}":
        raise ValueError(f"strict Ghia validation requires {_GRID}x{_GRID}")
    expected_h = 2.0 / (_GRID - 1)
    if not math.isclose(spacing, expected_h, rel_tol=0.0, abs_tol=_SETTINGS_TOLERANCE):
        raise ValueError("grid metadata and h are inconsistent with the [0,2] domain")
    if timestep <= 0.0 or tolerance <= 0.0:
        raise ValueError("dt and pressure_tol must be positive")
    if maximum_iterations <= 0 or not 0 < omega < 2:
        raise ValueError("pressure iteration settings are invalid")
    if not 0 <= iterations <= maximum_iterations:
        raise ValueError("pressure_iters is outside the configured iteration limit")
    if residual < 0.0 or residual > tolerance * (1.0 + 1e-12):
        raise ValueError("pressure residual is inconsistent with a converged solve")
    if step <= 0 or divergence < 0.0 or cfl < 0.0:
        raise ValueError("step, divergence, or CFL metadata are invalid")
    if metadata["finite"] != "true":
        raise ValueError("strict validation rejects non-finite solver states")
    if metadata["pressure_converged"] != "true":
        raise ValueError("strict validation rejects unconverged pressure solves")


def parse_projection_centerline(
    path: str | Path, *, expected_re: int | None = None
) -> CenterlineData:
    """Read one complete, finite and pressure-converged 129x129 solver block."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    metadata_index = next(
        (index for index, line in enumerate(lines) if line.startswith("# FlowLabLite centerline velocity")),
        None,
    )
    if metadata_index is None:
        raise ValueError("centerline CSV is missing its metadata comment")
    metadata = _parse_metadata(lines[metadata_index])
    header_index = metadata_index + 1
    while header_index < len(lines) and not lines[header_index].strip():
        header_index += 1
    if header_index == len(lines):
        raise ValueError("centerline CSV is missing its data header")

    components: dict[str, list[tuple[int, float, float]]] = {"u_vcl": [], "v_hcl": []}
    for row in csv.DictReader(lines[header_index:]):
        component = row.get("component")
        if component not in components:
            raise ValueError(f"unknown centerline component: {component!r}")
        try:
            index = int(row["index"])
            coordinate = float(row["coord"])
            value = float(row["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("centerline rows require numeric index, coord and value") from exc
        components[component].append((index, coordinate, value))

    u_coord, u_value = _validate_profile("u_vcl", components["u_vcl"], (0.0, 1.0))
    v_coord, v_value = _validate_profile("v_hcl", components["v_hcl"], (0.0, 0.0))
    data = CenterlineData(metadata, u_coord, u_value, v_coord, v_value)
    _validate_settings(data, expected_re=expected_re)
    return data


def strict_max_abs_pass(max_abs: float) -> bool:
    """Apply the thesis requirement literally: strict inequality, not <=."""
    return math.isfinite(max_abs) and max_abs < _MAX_ABS_LIMIT


def pointwise_errors(data: CenterlineData, reynolds: int) -> list[dict[str, object]]:
    if reynolds not in GHIA_REFERENCE:
        raise ValueError(f"no Ghia reference data are configured for Re={reynolds}")
    _validate_settings(data, expected_re=reynolds)
    reference = GHIA_REFERENCE[reynolds]
    rows: list[dict[str, object]] = []
    for component, solver_coord, solver_value in (
        ("u", data.u_coord, data.u_value),
        ("v", data.v_coord, data.v_value),
    ):
        reference_coord = reference[f"{component}_coord"]
        reference_value = reference[f"{component}_value"]
        predicted = interpolate(solver_coord, solver_value, reference_coord)
        for coordinate, expected, actual in zip(
            reference_coord, reference_value, predicted, strict=True
        ):
            signed_error = actual - expected
            rows.append(
                {
                    "re": reynolds,
                    "component": component,
                    "coord_over_l": coordinate / 2.0,
                    "reference": expected,
                    "flowlablite": actual,
                    "signed_error": signed_error,
                    "abs_error": abs(signed_error),
                }
            )
    return rows


def _time_converged(
    u_change: float, v_change: float, previous_divergence: float, divergence: float
) -> bool:
    roundoff = max(1e-15, abs(previous_divergence) * 1e-12)
    return (
        u_change < _PROFILE_CHANGE_LIMIT
        and v_change < _PROFILE_CHANGE_LIMIT
        and divergence <= previous_divergence + roundoff
    )


def _settings_signature(data: CenterlineData) -> tuple[object, ...]:
    keys = (
        "solver",
        "scheme",
        "Re",
        "grid",
        "h",
        "dt",
        "pressure_tol",
        "pressure_max_iter",
        "pressure_omega",
    )
    return tuple(data.metadata[key] for key in keys)


def _worst(rows: list[dict[str, object]], component: str) -> dict[str, object]:
    candidates = [row for row in rows if row["component"] == component]
    return dict(max(candidates, key=lambda row: float(row["abs_error"])))


def write_projection_summary(
    cases: Sequence[str | Path | CenterlineData],
    output_csv: str | Path,
    output_json: str | Path,
    pointwise_csv: str | Path,
) -> dict[str, object]:
    """Write strict 17-point Ghia metrics, convergence evidence, and provenance."""
    if not cases:
        raise ValueError("at least one centerline case is required")
    csv_path = Path(output_csv).resolve()
    json_path = Path(output_json).resolve()
    pointwise_path = Path(pointwise_csv).resolve()
    for path in (csv_path, json_path, pointwise_path):
        path.parent.mkdir(parents=True, exist_ok=True)

    normalized: list[tuple[CenterlineData, Path | None, str | None]] = []
    expected_re: int | None = None
    for case in cases:
        if isinstance(case, CenterlineData):
            data, source, source_hash = case, None, None
        else:
            source = Path(case).resolve()
            source_hash = _sha256(source)
            data = parse_projection_centerline(source)
        reynolds = int(round(float(data.metadata["Re"])))
        if reynolds not in GHIA_REFERENCE:
            raise ValueError(f"no Ghia reference data are configured for Re={reynolds}")
        if expected_re is None:
            expected_re = reynolds
        elif expected_re != reynolds:
            raise ValueError("all checkpoints must use the same Reynolds number")
        _validate_settings(data, expected_re=expected_re)
        normalized.append((data, source, source_hash))
    assert expected_re is not None

    reference = GHIA_REFERENCE[expected_re]
    rows: list[dict[str, object]] = []
    signature: tuple[object, ...] | None = None
    previous: CenterlineData | None = None
    previous_step: int | None = None
    previous_divergence: float | None = None
    final_pointwise: list[dict[str, object]] = []
    for data, source, source_hash in normalized:
        current_signature = _settings_signature(data)
        if signature is None:
            signature = current_signature
        elif signature != current_signature:
            raise ValueError("checkpoints must use identical numerical settings")
        step = int(data.metadata["step"])
        if previous_step is not None and step <= previous_step:
            raise ValueError("checkpoint steps must be strictly increasing")
        divergence = float(data.metadata["div_norm"])
        u_metrics = compute_errors(
            data.u_coord, data.u_value, reference["u_coord"], reference["u_value"]
        )
        v_metrics = compute_errors(
            data.v_coord, data.v_value, reference["v_coord"], reference["v_value"]
        )
        u_change = None if previous is None else profile_change(
            previous.u_coord, previous.u_value, data.u_coord, data.u_value
        )
        v_change = None if previous is None else profile_change(
            previous.v_coord, previous.v_value, data.v_coord, data.v_value
        )
        time_converged = bool(
            previous_divergence is not None
            and u_change is not None
            and v_change is not None
            and _time_converged(u_change, v_change, previous_divergence, divergence)
        )
        case_pointwise = pointwise_errors(data, expected_re)
        row: dict[str, object] = {
            "re": expected_re,
            "scheme": data.metadata["scheme"],
            "grid": data.metadata["grid"],
            "steps": step,
            "time": step * float(data.metadata["dt"]),
            "dt": data.metadata["dt"],
            "pressure_tol": data.metadata["pressure_tol"],
            "pressure_max_iter": data.metadata["pressure_max_iter"],
            "pressure_omega": data.metadata["pressure_omega"],
            "pressure_iters": data.metadata["pressure_iters"],
            "pressure_residual": data.metadata["pressure_residual"],
            "div_norm": divergence,
            "cfl_max": data.metadata["cfl_max"],
            "u_rms": u_metrics.rms,
            "u_relative_l2": u_metrics.relative_l2,
            "u_max_abs": u_metrics.max_abs,
            "u_max_pass_strict_5pct": strict_max_abs_pass(u_metrics.max_abs),
            "v_rms": v_metrics.rms,
            "v_relative_l2": v_metrics.relative_l2,
            "v_max_abs": v_metrics.max_abs,
            "v_max_pass_strict_5pct": strict_max_abs_pass(v_metrics.max_abs),
            "worst_u": _worst(case_pointwise, "u"),
            "worst_v": _worst(case_pointwise, "v"),
            "u_change": u_change,
            "v_change": v_change,
            "time_converged": time_converged,
        }
        row["ghia_max_abs_pass"] = bool(
            row["u_max_pass_strict_5pct"] and row["v_max_pass_strict_5pct"]
        )
        row["strict_validation_pass"] = bool(
            time_converged and row["ghia_max_abs_pass"]
        )
        if source is not None and source_hash is not None:
            row["source_csv"] = os.path.relpath(source, json_path.parent)
            row["source_sha256"] = source_hash
        rows.append(row)
        previous = data
        previous_step = step
        previous_divergence = divergence
        final_pointwise = case_pointwise

    csv_fields = (
        "re",
        "scheme",
        "grid",
        "steps",
        "time",
        "dt",
        "pressure_tol",
        "pressure_max_iter",
        "pressure_omega",
        "pressure_iters",
        "pressure_residual",
        "div_norm",
        "cfl_max",
        "u_rms",
        "u_relative_l2",
        "u_max_abs",
        "u_max_pass_strict_5pct",
        "v_rms",
        "v_relative_l2",
        "v_max_abs",
        "v_max_pass_strict_5pct",
        "u_change",
        "v_change",
        "time_converged",
        "ghia_max_abs_pass",
        "strict_validation_pass",
    )
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=csv_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    with pointwise_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(final_pointwise[0]))
        writer.writeheader()
        writer.writerows(final_pointwise)

    selected = rows[-1]
    payload: dict[str, object] = {
        "re": expected_re,
        "reference": "Ghia, Ghia and Shin (1982), Tables I and II, 129x129",
        "reference_point_count_per_component": 17,
        "acceptance_rule": "u_max_abs < 0.05 and v_max_abs < 0.05",
        "time_convergence_rule": "adjacent u/v relative L2 change < 0.005 and div_norm non-increasing",
        "cases": rows,
        "selected_case": selected,
        "ghia_max_abs_pass": selected["ghia_max_abs_pass"],
        "time_converged": selected["time_converged"],
        "strict_validation_pass": selected["strict_validation_pass"],
    }
    json_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return payload
