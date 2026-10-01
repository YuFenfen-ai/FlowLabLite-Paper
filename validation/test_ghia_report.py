from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from validation.ghia_report import render_report
from validation.ghia_study import write_summary


def synthetic_centerline(reynolds: int, step: int = 20000, grid: int = 129) -> str:
    h = 2.0 / (grid - 1)
    lines = [
        f"# FlowLabLite centerline velocity | solver=chorin | Re={reynolds} | step={step} | grid={grid}x{grid} | h={h} | dt=0.001 | pressure_iters=50 | div_norm=0.01",
        "component,index,coord,value",
    ]
    for index in range(grid):
        coord = index * h
        value = 0.0 if index == 0 else 1.0 if index == grid - 1 else -0.1 * coord * (2.0 - coord)
        lines.append(f"u_vcl,{index},{coord},{value}")
    for index in range(grid):
        coord = index * h
        value = 0.0 if index in (0, grid - 1) else 0.05 * (1.0 - coord) * coord * (2.0 - coord)
        lines.append(f"v_hcl,{index},{coord},{value}")
    return "\n".join(lines) + "\n"


class GhiaReportTests(unittest.TestCase):
    def test_report_writes_three_case_table_plot_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            summaries: list[Path] = []
            for reynolds in (100, 400, 1000):
                case = root / f"re{reynolds}" / "centerline.csv"
                case.parent.mkdir(parents=True)
                case.write_text(synthetic_centerline(reynolds), encoding="utf-8")
                summary = case.parent / "summary.json"
                write_summary((case,), case.parent / "summary.csv", summary)
                summaries.append(summary)
            manifest = render_report(
                summaries,
                root / "combined.csv",
                root / "comparison.png",
                root / "comparison.pdf",
                root / "manifest.json",
            )
            self.assertEqual([row["re"] for row in manifest["rows"]], [100, 400, 1000])
            self.assertGreater((root / "comparison.png").stat().st_size, 0)
            self.assertGreater((root / "comparison.pdf").stat().st_size, 0)
            saved = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["settings"]["grid"], "129x129")


if __name__ == "__main__":
    unittest.main()
