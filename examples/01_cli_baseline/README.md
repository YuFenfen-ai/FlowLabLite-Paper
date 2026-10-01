# Example 1 — CLI Baseline Run

This example shows how to run all four FlowLabLite solvers from the command line
and inspect the results without a browser.

## What you will get

- A `results.json` file containing full-field data (1681 points) for all four solvers
- Console output with timing and field statistics
- A local viewer to explore the results interactively

## Steps

### 1. From the project root, run the simulation

```bash
bash run_local.sh
```

Expected console output:
```
[run_local] Building and running FlowLabLite locally...
[run_local] Target: wasm (Wasmtime runtime)
[run_local] Simulation completed. Extracting JSON...
[run_local] Valid JSON saved to: results.json
[run_local]   Chorin grid  : 1681 points
[run_local]   SIMPLE grid  : 1681 points
```

### 2. Inspect key statistics with Node.js

```bash
node -e "
const r = require('./results.json');
const s = r.statistics;
console.log('=== Chorin solver (Re=20, 500 steps) ===');
console.log('  max u:          ', s.max_u.toFixed(6));
console.log('  max v:          ', s.max_v.toFixed(6));
console.log('  centre u:       ', s.center_u.toFixed(6), '  (negative = backflow vortex confirmed)');
console.log('  divergence norm:', s.divergence_norm.toExponential(3));
"
```

Expected output:
```
=== Chorin solver (Re=20, 500 steps) ===
  max u:           1.000000
  max v:           0.253067
  centre u:        -0.101540   (negative = backflow vortex confirmed)
  divergence norm: 1.437e-2
```

### 3. Compare all four solvers

```bash
node -e "
const r = require('./results.json');
console.log('Solver         | centre_u      | max_v         | div_norm');
console.log('---------------|---------------|---------------|----------');
const cs = r.statistics;
console.log('Chorin         |', cs.center_u.toFixed(6), '  |', cs.max_v.toFixed(6), '  |', cs.divergence_norm.toExponential(2));
const ss = r.simple_statistics;
console.log('SIMPLE         |', ss.center_u.toFixed(6), '  |', ss.max_v.toFixed(6), '  |', ss.divergence_norm.toExponential(2));
const ps = r.pcg_statistics;
console.log('Chorin-PCG     |', ps.center_u.toFixed(6), '  |', ps.max_v.toFixed(6), '  |', ps.divergence_norm.toExponential(2));
const ms = r.mac_statistics;
console.log('MAC            |', ms.center_u.toFixed(6), '  |', ms.max_v.toFixed(6), '  |', ms.divergence_norm.toExponential(2));
"
```

### 4. Open the local viewer

Open `cmd/main/local_viewer.html` in a browser and drag `results.json` onto it.
Click each solver tab to compare flow fields side by side.

## What the results mean

| Metric | Physical meaning | Expected range |
|---|---|---|
| `max_u = 1.0` | Lid velocity boundary condition maintained | = 1.000 |
| `centre_u < 0` | Clockwise primary vortex exists | − 0.05 to − 0.15 |
| `max_v > 0` | Flow turns the corner at the lid | 0.1 to 0.4 |
| `divergence_norm` | Mass conservation error | < 0.05 (GS), < 1e-4 (PCG/MAC) |

## Saving a custom run

```bash
bash run_local.sh re20_500steps.json    # custom filename
```
