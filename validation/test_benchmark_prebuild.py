from __future__ import annotations

import unittest
from pathlib import Path

from validation.benchmark_repro import build_prebuild_command


class BenchmarkPrebuildTests(unittest.TestCase):
    def test_prebuild_uses_moon_build_and_never_moon_run(self) -> None:
        command = build_prebuild_command(Path("moon.exe"))
        self.assertEqual(command[:3], ["moon.exe", "build", "cmd/main"])
        self.assertIn("--release", command)
        self.assertEqual(command[command.index("--target") + 1], "wasm")
        self.assertNotIn("run", command)


if __name__ == "__main__":
    unittest.main()
