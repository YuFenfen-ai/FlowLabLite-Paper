from __future__ import annotations

import unittest
from pathlib import Path

from validation.run_projection_tvd_study import build_solver_command
from validation.run_ghia_study import split_checkpoint_stream
from validation.test_projection_tvd_study import projection_centerline_text


class ProjectionTvdRunnerTests(unittest.TestCase):
    def test_command_preserves_all_numerical_settings_explicitly(self) -> None:
        command = build_solver_command(
            Path("moon.exe"),
            "wasm-gc",
            129,
            (500, 1000, 2000),
            1000,
            "tvd-vanleer",
            0.0005,
            1e-6,
            12000,
            1.95,
        )
        expected = {
            "--target": "wasm-gc",
            "--format": "centerline",
            "--solver": "projection-tvd",
            "--grid": "129",
            "--checkpoints": "500,1000,2000",
            "--re": "1000",
            "--scheme": "tvd-vanleer",
            "--dt": "0.0005",
            "--pressure-tol": "1e-06",
            "--pressure-max-iter": "12000",
            "--pressure-omega": "1.95",
        }
        self.assertIn("--release", command)
        for option, value in expected.items():
            with self.subTest(option=option):
                self.assertEqual(command[command.index(option) + 1], value)

    def test_existing_stream_splitter_accepts_projection_blocks(self) -> None:
        stream = projection_centerline_text(step=10) + projection_centerline_text(step=20)
        blocks = split_checkpoint_stream(stream)
        self.assertEqual(tuple(blocks), (10, 20))


if __name__ == "__main__":
    unittest.main()
