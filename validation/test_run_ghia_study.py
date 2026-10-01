from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from validation.run_ghia_study import build_solver_command, split_checkpoint_stream
from validation.test_re100_study import centerline_text


class GhiaStudyRunnerTests(unittest.TestCase):
    def test_split_checkpoint_stream_returns_one_parseable_block_per_step(self) -> None:
        stream = centerline_text(step=2) + "\n" + centerline_text(step=5, shift=1e-6)
        blocks = split_checkpoint_stream(stream)
        self.assertEqual(tuple(blocks), (2, 5))
        with tempfile.TemporaryDirectory() as directory:
            for step, text in blocks.items():
                path = Path(directory) / f"s{step}.csv"
                path.write_text(text, encoding="utf-8")
                self.assertIn(f"step={step}", path.read_text(encoding="utf-8").splitlines()[0])

    def test_split_checkpoint_stream_rejects_duplicate_or_non_csv_output(self) -> None:
        with self.assertRaises(ValueError):
            split_checkpoint_stream(centerline_text(step=2) + centerline_text(step=2))
        with self.assertRaises(ValueError):
            split_checkpoint_stream("not a centerline stream")

    def test_solver_command_is_release_wasm_chorin_with_cumulative_checkpoints(self) -> None:
        command = build_solver_command(Path("moon.exe"), 129, (5000, 10000, 20000), 400)
        self.assertIn("--release", command)
        self.assertEqual(command[command.index("--target") + 1], "wasm")
        self.assertEqual(command[command.index("--solver") + 1], "chorin")
        self.assertEqual(command[command.index("--grid") + 1], "129")
        self.assertEqual(command[command.index("--checkpoints") + 1], "5000,10000,20000")
        self.assertEqual(command[command.index("--re") + 1], "400")

    def test_solver_command_passes_runtime_dt(self) -> None:
        command = build_solver_command(Path("moon.exe"), 257, (5000,), 100, 0.0005)
        self.assertEqual(command[command.index("--dt") + 1], "0.0005")


if __name__ == "__main__":
    unittest.main()
