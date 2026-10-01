# Browser launch guide

Reproduces the browser demo and the paper's quantitative runs.

## 1. Build the WASM bundle

A plain `moon build` does **not** export the API functions on moon
0.1.20260309 (the version used for the baseline build); re-link with the
explicit export list instead:

```bash
bash build_wasm.sh release
```

The script prints the export count when finished. Validation-era runs in
the paper were executed with moon 0.1.20260618.

## 2. Serve the repository root over HTTP

```bash
python -m http.server 8000
```

## 3. Open the visualization page

```
http://localhost:8000/cmd/main/main.html
```

(Chrome / Edge / Brave 115+; the page uses wasm-gc with stringref.)
The page shows velocity/pressure heatmaps with a streamline overlay for
the four solvers (Chorin / SIMPLE / Chorin-PCG / MAC) from the built-in
41 x 41 demonstration run.

## 4. Reproduce the paper's Re = 100 baseline (Section 4.3)

The quantitative runs are driven from the command line:

```bash
moon run cmd/main --release --target wasm -- \
  --format centerline --solver chorin --grid 129 --re 100 \
  --checkpoints 5000,10000,20000,40000
```

Baseline settings: grid = 129 x 129, dt = 0.001, 40,000 steps,
50 Gauss-Seidel pressure iterations per step.

## 5. Reproduce the grid-convergence and natural-convection studies

```bash
# Section 4.4 grid convergence (see data/paper/sec4_6_natural_convection/commands.txt
# and validation/run_ghia_study.py for the exact commands)
python -m validation.run_ghia_study --grid 65 --re 100 --checkpoints 5000 10000 20000 40000

# Section 4.6 Boussinesq natural convection (De Vahl Davis benchmark):
moon run cmd/main --release --target wasm -- --format centerline --solver boussinesq \
  --bous-nu 0.071 --bous-alpha 0.1 --bous-beta 0.8875 \
  --checkpoints 4000,8000,12000,16000,20000,24000,28000,32000
```

Full command records and outputs: `data/paper/` (curated) and
`validation/results/` (raw run archives).
