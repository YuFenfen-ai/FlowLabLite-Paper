from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from validation import re100_study
from validation.re100_study import compute_errors, interpolate, parse_centerline, profile_change


def centerline_text(*, step=100, grid="41x41", solver="chorin", re=100, dt=0.001,
                    pressure_iters=50, div_norm=0.01, shift=0.0, start=0.0,
                    include_v=True, nonmonotone=False, nonfinite=False, bad_endpoint=False) -> str:
    coordinates = [start + 2.0 * index / 40.0 for index in range(41)]
    if nonmonotone:
        coordinates[20] = coordinates[19]
    lines = [
        f"# FlowLabLite centerline velocity | solver={solver} | Re={re} | step={step} | grid={grid} | h=0.05 | dt={dt} | pressure_iters={pressure_iters} | div_norm={div_norm}",
        "component,index,coord,value",
    ]
    for index, coordinate in enumerate(coordinates):
        value = 0.0 if index == 0 else 1.0 if index == 40 else -0.1 + 0.02 * coordinate + shift * coordinate * (2.0 - coordinate)
        if nonfinite and index == 10:
            value = float("nan")
        if bad_endpoint and index == 40:
            value = 0.9
        lines.append(f"u_vcl,{index},{coordinate},{value}")
    if include_v:
        for index, coordinate in enumerate(coordinates):
            value = 0.0 if index in (0, 40) else 0.05 + shift * coordinate * (2.0 - coordinate)
            lines.append(f"v_hcl,{index},{coordinate},{value}")
    return "\n".join(lines) + "\n"


class Re100StudyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary_directory.name)
        self.fixture = self.write_case("centerline.csv")

    def tearDown(self) -> None:
        self.temporary_directory.cleanup()

    def write_case(self, name: str, **kwargs: object) -> Path:
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(centerline_text(**kwargs), encoding="utf-8")
        return path

    def test_parse_centerline_metadata_and_components(self) -> None:
        data = parse_centerline(self.fixture)
        self.assertEqual(data.metadata["Re"], 100.0)
        self.assertEqual(data.metadata["grid"], "41x41")
        self.assertEqual(data.u_coord[0], 0.0)
        self.assertEqual(data.u_coord[-1], 2.0)
        self.assertEqual(data.v_value[-1], 0.0)

    def test_parse_centerline_rejects_malformed_profiles_and_settings(self) -> None:
        for name, kwargs in (
            ("missing.csv", {"include_v": False}),
            ("monotone.csv", {"nonmonotone": True}),
            ("finite.csv", {"nonfinite": True}),
            ("endpoint.csv", {"bad_endpoint": True}),
            ("domain.csv", {"start": 0.1}),
            ("solver.csv", {"solver": "simple"}),
            ("dt.csv", {"dt": 0.002}),
            ("pressure.csv", {"pressure_iters": 40}),
        ):
            with self.subTest(name=name), self.assertRaises(ValueError):
                parse_centerline(self.write_case(name, **kwargs))
        with self.assertRaises(ValueError):
            parse_centerline(self.write_case("re.csv", re=20), expected_re=100)

    def test_compute_errors_validates_reference_and_interpolation_range(self) -> None:
        with self.assertRaises(ValueError):
            interpolate((0.0, 1.0), (0.0, 1.0), (-0.01,))
        with self.assertRaises(ValueError):
            compute_errors((0.0, 1.0), (0.0, 1.0), (0.0, float("nan")), (0.0, 1.0))

    def test_compute_errors_exact_match_is_zero(self) -> None:
        result = compute_errors((0.0, 1.0), (0.0, 1.0), (0.0, 1.0), (0.0, 1.0))
        self.assertEqual((result.rms, result.relative_l2, result.max_abs), (0.0, 0.0, 0.0))

    def test_profile_change_uses_relative_l2(self) -> None:
        change = profile_change((0.0, 1.0), (1.0, 2.0), (0.0, 1.0), (1.0, 2.02))
        self.assertGreater(change, 0.0)
        self.assertLess(change, 0.02)

    def test_time_convergence_requires_strict_half_percent_and_no_divergence_increase(self) -> None:
        self.assertFalse(re100_study._is_time_converged(0.005, 0.004, 0.01, 0.01))
        self.assertFalse(re100_study._is_time_converged(0.004, 0.004, 0.01, 0.010001))
        self.assertTrue(re100_study._is_time_converged(0.004, 0.004, 0.01, 0.01 + 1e-15))

    def test_summary_rejects_cross_grid_and_out_of_order_cases(self) -> None:
        csv_path, json_path = self.root / "summary.csv", self.root / "summary.json"
        with self.assertRaises(ValueError):
            re100_study.write_summary((self.write_case("a.csv", step=100), self.write_case("b.csv", step=200, grid="81x81")), csv_path, json_path)
        with self.assertRaises(ValueError):
            re100_study.write_summary((self.write_case("c.csv", step=200), self.write_case("d.csv", step=100)), csv_path, json_path)

    def test_summary_persists_only_a_strictly_converged_time_selection_with_relative_provenance(self) -> None:
        first = self.write_case("cases/first.csv", step=100)
        second = self.write_case("cases/second.csv", step=200, shift=1e-6)
        output_csv, output_json = self.root / "summary/summary.csv", self.root / "summary/summary.json"
        summary = re100_study.write_summary((first, second), output_csv, output_json)
        self.assertEqual(summary["selected_time_step"], 200)
        self.assertEqual(summary["selected_case"]["steps"], 200)
        self.assertFalse(Path(summary["selected_case"]["source_csv"]).is_absolute())
        self.assertIn("source_sha256", summary["selected_case"])

    def test_plot_validates_summary_checksum_and_prints_cli_manifest(self) -> None:
        first = self.write_case("cases/first.csv", step=100)
        second = self.write_case("cases/second.csv", step=200, shift=1e-6)
        summary_json = self.root / "summary.json"
        re100_study.write_summary((first, second), self.root / "summary.csv", summary_json)
        second.write_text(centerline_text(step=200, shift=2e-6), encoding="utf-8")
        with self.assertRaises(ValueError):
            re100_study.plot_comparison(Path(__file__).parents[1] / "data" / "sample_results.json", summary_json, self.root / "bad.png", self.root / "bad.pdf")
        re100_study.write_summary((first, second), self.root / "summary.csv", summary_json)
        result = subprocess.run([sys.executable, "-m", "validation.re100_study", "plot", "--re20", str(Path(__file__).parents[1] / "data" / "sample_results.json"), "--re100", str(summary_json), "--png", str(self.root / "ok.png"), "--pdf", str(self.root / "ok.pdf")], capture_output=True, text=True, check=True)
        self.assertEqual(json.loads(result.stdout)["selected_time_step"], 200)

    def test_plot_comparison_writes_png_and_pdf(self) -> None:
        output_png, output_pdf = self.root / "comparison.png", self.root / "comparison.pdf"
        manifest = re100_study.plot_comparison(Path(__file__).parents[1] / "data" / "sample_results.json", self.fixture, output_png, output_pdf)
        self.assertGreater(output_png.stat().st_size, 0)
        self.assertGreater(output_pdf.stat().st_size, 0)
        self.assertEqual(manifest["labels"][0], "FlowLabLite Re=20 (41\u00d741, 500 steps)")


    def test_ghia_reference_contains_three_129_grid_cases(self) -> None:
        self.assertEqual(set(re100_study.GHIA_REFERENCE), {100, 400, 1000})
        for reynolds, reference in re100_study.GHIA_REFERENCE.items():
            with self.subTest(reynolds=reynolds):
                self.assertEqual(len(reference["u_coord"]), 17)
                self.assertEqual(len(reference["u_value"]), 17)
                self.assertEqual(len(reference["v_coord"]), 17)
                self.assertEqual(len(reference["v_value"]), 17)

if __name__ == "__main__":
    unittest.main()