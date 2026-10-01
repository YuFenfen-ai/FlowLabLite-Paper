from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from validation.performance_report import render_performance_report


class PerformanceReportTests(unittest.TestCase):
    def test_report_validates_method_and_writes_table_plot_manifest(self) -> None:
        summary = {
            "method": {
                "timed_scope": "prebuilt release WASM executed by moonrun; compilation excluded",
                "warmups_discarded_per_solver": 3,
                "measured_repetitions_per_solver": 15,
                "ordering": "rotating solver order per measured round",
                "grid": "41x41",
                "workloads": {"chorin": 500, "simple": 100, "pcg": 500, "mac": 500},
            },
            "summary": {
                solver: {
                    "n": 15,
                    "mean_ms": 1000.0 + index * 200.0,
                    "sample_std_ms": 25.0 + index,
                    "median_ms": 995.0 + index * 200.0,
                    "min_ms": 950.0 + index * 200.0,
                    "max_ms": 1050.0 + index * 200.0,
                    "cv_percent": 2.5,
                }
                for index, solver in enumerate(("chorin", "simple", "pcg", "mac"))
            },
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "summary.json"
            source.write_text(json.dumps(summary), encoding="utf-8")
            manifest = render_performance_report(
                source,
                root / "table.csv",
                root / "plot.png",
                root / "plot.pdf",
                root / "manifest.json",
            )
            self.assertTrue((root / "plot.png").is_file())
            self.assertTrue((root / "plot.pdf").is_file())
            self.assertTrue((root / "manifest.json").is_file())
            with (root / "table.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([row["solver"] for row in rows], ["Chorin-GS", "SIMPLE", "Chorin-PCG", "MAC"])
            self.assertEqual(manifest["method"]["measured_repetitions_per_solver"], 15)
            self.assertIn("plot", manifest)
            highest_error_bar_seconds = max(
                (float(row["mean_ms"]) + float(row["sample_std_ms"])) / 1000.0
                for row in rows
            )
            self.assertGreaterEqual(
                float(manifest["plot"]["ymax_seconds"]),
                1.2 * highest_error_bar_seconds,
            )

    def test_report_rejects_missing_warmup_or_wrong_sample_count(self) -> None:
        payload = {
            "method": {
                "warmups_discarded_per_solver": 0,
                "measured_repetitions_per_solver": 15,
                "grid": "41x41",
                "workloads": {"chorin": 500, "simple": 100, "pcg": 500, "mac": 500},
            },
            "summary": {},
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "summary.json"
            source.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "warm-up"):
                render_performance_report(
                    source,
                    root / "table.csv",
                    root / "plot.png",
                    root / "plot.pdf",
                    root / "manifest.json",
                )


if __name__ == "__main__":
    unittest.main()
