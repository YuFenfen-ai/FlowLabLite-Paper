from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from validation.projection_tvd_study import (
    pointwise_errors,
    parse_projection_centerline,
    strict_max_abs_pass,
    write_projection_summary,
)


def projection_centerline_text(
    *,
    step: int = 100,
    solver: str = "projection-tvd",
    scheme: str = "tvd-vanleer",
    reynolds: int = 400,
    grid: int = 129,
    finite: str = "true",
    pressure_converged: str = "true",
    shift: float = 0.0,
    nonfinite: bool = False,
    omit: str | None = None,
) -> str:
    h = 2.0 / (grid - 1)
    metadata = {
        "solver": solver,
        "scheme": scheme,
        "Re": str(reynolds),
        "step": str(step),
        "grid": f"{grid}x{grid}",
        "h": str(h),
        "dt": "0.001",
        "pressure_tol": "1e-6",
        "pressure_max_iter": "5000",
        "pressure_omega": "1.9",
        "pressure_iters": "712",
        "pressure_residual": "9.9e-7",
        "pressure_converged": pressure_converged,
        "div_norm": "1.2e-5",
        "cfl_max": "0.064",
        "finite": finite,
    }
    if omit is not None:
        metadata.pop(omit)
    header = "# FlowLabLite centerline velocity | " + " | ".join(
        f"{key}={value}" for key, value in metadata.items()
    )
    lines = [header, "component,index,coord,value"]
    for index in range(grid):
        coordinate = index * h
        value = coordinate / 2.0
        if index == 0:
            value = 0.0
        elif index == grid - 1:
            value = 1.0
        else:
            value += shift * coordinate * (2.0 - coordinate)
        if nonfinite and index == 10:
            value = float("nan")
        lines.append(f"u_vcl,{index},{coordinate},{value}")
    for index in range(grid):
        coordinate = index * h
        value = 0.0 if index in (0, grid - 1) else shift * coordinate * (2.0 - coordinate)
        lines.append(f"v_hcl,{index},{coordinate},{value}")
    return "\n".join(lines) + "\n"


class ProjectionTvdStudyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_case(self, name: str, **kwargs: object) -> Path:
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(projection_centerline_text(**kwargs), encoding="utf-8")
        return path

    def test_parse_accepts_only_complete_converged_129_grid_output(self) -> None:
        data = parse_projection_centerline(self.write_case("valid.csv"), expected_re=400)
        self.assertEqual(data.metadata["solver"], "projection-tvd")
        self.assertEqual(data.metadata["scheme"], "tvd-vanleer")
        self.assertEqual(len(data.u_coord), 129)
        self.assertEqual(len(data.v_coord), 129)

    def test_parse_rejects_invalid_identity_grid_status_and_data(self) -> None:
        cases = (
            ("solver.csv", {"solver": "chorin"}),
            ("scheme.csv", {"scheme": "central"}),
            ("grid.csv", {"grid": 41}),
            ("finite.csv", {"finite": "false"}),
            ("pressure.csv", {"pressure_converged": "false"}),
            ("missing.csv", {"omit": "pressure_residual"}),
            ("nan.csv", {"nonfinite": True}),
        )
        for name, kwargs in cases:
            with self.subTest(name=name), self.assertRaises(ValueError):
                parse_projection_centerline(self.write_case(name, **kwargs), expected_re=400)

    def test_five_percent_gate_is_strict(self) -> None:
        self.assertTrue(strict_max_abs_pass(0.049999999))
        self.assertFalse(strict_max_abs_pass(0.05))
        self.assertFalse(strict_max_abs_pass(0.050000001))

    def test_pointwise_errors_keep_all_17_ghia_points_per_component(self) -> None:
        data = parse_projection_centerline(self.write_case("points.csv"), expected_re=400)
        rows = pointwise_errors(data, 400)
        self.assertEqual(len(rows), 34)
        self.assertEqual(sum(row["component"] == "u" for row in rows), 17)
        self.assertEqual(sum(row["component"] == "v" for row in rows), 17)

    def test_summary_records_strict_gate_worst_points_and_provenance(self) -> None:
        first = self.write_case("cases/s100.csv", step=100)
        second = self.write_case("cases/s200.csv", step=200, shift=1e-7)
        output_csv = self.root / "summary/time_convergence.csv"
        output_json = self.root / "summary/time_convergence.json"
        pointwise_csv = self.root / "summary/pointwise_errors.csv"
        summary = write_projection_summary(
            (first, second), output_csv, output_json, pointwise_csv
        )
        self.assertEqual(summary["re"], 400)
        self.assertEqual(summary["reference_point_count_per_component"], 17)
        self.assertIn("worst_u", summary["selected_case"])
        self.assertIn("worst_v", summary["selected_case"])
        self.assertIn("ghia_max_abs_pass", summary["selected_case"])
        self.assertEqual(
            summary["selected_case"]["strict_validation_pass"],
            summary["selected_case"]["ghia_max_abs_pass"]
            and summary["selected_case"]["time_converged"],
        )
        self.assertIn("source_sha256", summary["selected_case"])
        self.assertEqual(len(pointwise_csv.read_text(encoding="utf-8").splitlines()), 35)
        persisted = json.loads(output_json.read_text(encoding="utf-8"))
        self.assertEqual(persisted["acceptance_rule"], "u_max_abs < 0.05 and v_max_abs < 0.05")
        self.assertEqual(
            persisted["ghia_max_abs_pass"],
            persisted["selected_case"]["ghia_max_abs_pass"],
        )


if __name__ == "__main__":
    unittest.main()
