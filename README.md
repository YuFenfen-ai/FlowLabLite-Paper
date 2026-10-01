# FlowLabLite Thesis Final Submission (2D)

This is the frozen two-dimensional source snapshot used by the thesis quantitative validation.

## Version boundary

- Git commit: `4d0162b76e85b53413336789a87a165fc25becaf`
- Snapshot base: aa8c89c (2026-08-11); compatibility-fix commit: 4d0162b
- Package date: 2026-08-17
- Included: tracked 2D source, MoonBit tests, HTML pages with embedded JavaScript, validation scripts/results, benchmark sources, documentation, configuration, and sample data.
- Excluded: the complete `cmd/main3d/` package and `build_wasm_3d.sh`.
- No MoonBit solver source was edited while creating this package.

The historical commit stored `README.md` with an erroneous symbolic-link mode and described the removed experimental 3D package. This submission README replaces it so the package description matches the 2D-only boundary.

## Thesis line-count definition

The thesis figure **5843 lines** means physical lines in the 14 non-test MoonBit files under `cmd/main`. It is the `Lines` column reported by scc 3.7.0, not the `Code` column.

| Scope | Files | Lines | Code | Comments | Blanks |
|---|---:|---:|---:|---:|---:|
| 2D non-test MoonBit | 14 | 5843 | 3955 | 1357 | 531 |
| 2D MoonBit tests | 3 | 2852 | 1979 | 672 | 201 |
| All 2D MoonBit | 17 | 8695 | 5934 | 2029 | 732 |
| HTML with embedded JavaScript | 2 | 3928 | 3514 | 68 | 346 |
| All Shell scripts in the package | 4 | 458 | 338 | 80 | 40 |

The four Shell files are `build_wasm.sh`, `run_local.sh`, `scan_re.sh`, and `bench/run_bench.sh`. Running scc on the entire package also counts tests, validation programs, benchmarks, and documents, so that total is intentionally different from 5843.

Reproduce the core statistic from the package root:

```powershell
scc.exe cmd\main --include-ext mbt -M ".*wbtest\.mbt$" --no-cocomo --no-size
```

## Verification

```bash
moon test
```

Expected for this frozen 2D snapshot: 156 passed, 0 failed.

See `SUBMISSION_MANIFEST.txt`, `SCC_NONTEST.csv`, `SCC_TEST.csv`, and `SHA256SUMS.txt` for reproducible evidence.
