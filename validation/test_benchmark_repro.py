from __future__ import annotations

import unittest
from pathlib import Path

from validation.benchmark_repro import build_runtime_command, summarize_samples


class BenchmarkReproTests(unittest.TestCase):
    def test_summary_uses_sample_standard_deviation_and_reports_cv(self) -> None:
        summary = summarize_samples((1.0, 2.0, 3.0))
        self.assertEqual(summary["n"], 3)
        self.assertEqual(summary["mean_ms"], 2.0)
        self.assertEqual(summary["sample_std_ms"], 1.0)
        self.assertEqual(summary["median_ms"], 2.0)
        self.assertEqual(summary["min_ms"], 1.0)
        self.assertEqual(summary["max_ms"], 3.0)
        self.assertEqual(summary["cv_percent"], 50.0)

    def test_runtime_command_uses_prebuilt_wasm_and_one_solver(self) -> None:
        command = build_runtime_command(Path("moonrun.exe"), Path("main.wasm"), "simple", 100)
        self.assertEqual(command[:3], ["moonrun.exe", "main.wasm", "--"])
        self.assertEqual(command[command.index("--solver") + 1], "simple")
        self.assertEqual(command[command.index("--steps") + 1], "100")
        self.assertEqual(command[command.index("--format") + 1], "centerline")


if __name__ == "__main__":
    unittest.main()
