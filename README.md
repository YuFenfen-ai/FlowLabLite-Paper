# FlowLabLite

FlowLabLite is a lightweight, AI-native CFD solver for 2D lid-driven cavity flow,
written in [MoonBit](https://www.moonbitlang.com/).
It compiles to WebAssembly (wasm-gc) and ships an interactive browser-based
visualization page (`cmd/main/main.html`) — including velocity heatmaps, pressure
heatmaps, and **streamline** overlay — with no runtime dependencies.

> **Online Demo**: serve the repo over HTTP and open `cmd/main/main.html` in Chrome 115+.  
> See [INSTALL.md](INSTALL.md) for one-command local serving instructions.

You can also view local data by opening `cmd/main/local_viewer.html` in a browser.

![Application Screenshot](./images/local_viewer.png) 

Four independent solvers are provided:

| Solver | Algorithm | Grid | Use case |
|---|---|---|---|
| **Chorin** | Projection method, explicit time-marching | 41×41 nodes | Transient simulation, 500 time steps |
| **SIMPLE** | Pressure-correction iteration, under-relaxation | 41×41 nodes | Steady-state iterative solve |
| **Chorin-PCG** | Projection method + PCG pressure solve | 41×41 nodes (collocated) | Transient, faster pressure convergence |
| **MAC** | Harlow-Welch staggered grid, PCG pressure | 40×40 cells | Transient, divergence-free by construction |

---


## Code Structure

```
FlowLabLite/
├── cmd/main/
│   ├── main.mbt              # Solver control, all solver paths + WASM API
│   ├── main_bench.mbt        # Micro-benchmark (fib baseline)
│   ├── main_wbtest.mbt       # White-box tests T1–T58 (original suite)
│   ├── main_ext_wbtest.mbt   # White-box tests T59–T154 (extended suite)
│   ├── solver_projection_tvd.mbt  # Projection-TVD extension (Section 6)
│   ├── main.html             # Browser visualization (heatmaps + streamline overlay)
│   ├── local_viewer.html     # Browser viewer for locally computed JSON (4-solver tabs)
│   ├── moon.pkg.json         # Package config: imports, link/exports (121 functions)
│   └── ...                   # Other solver and I/O modules
├── lib/
│   └── moon.pkg.json         # Library package placeholder
├── docs/                     # Documentation; launch guide: docs/launch.md
├── examples/                 # Example CLI and AI workflow
├── validation/               # Validation harness + raw run archives (Sections 4.3-4.6)
├── data/                     # Paper data: data/paper/ (curated) + 17-point Ghia CSVs
├── build_wasm.sh             # WASM build script (explicit 121-function export list)
├── run_local.sh              # Run solver locally + extract JSON results
├── TEST_RESULTS.log          # Recorded moon test output (176 passed / 0 failed)
├── LICENSE                   # Apache-2.0
└── ...```

---


## Test Suite

**176 tests, all passing** (`moon test cmd/main --target wasm`, recorded in `TEST_RESULTS.log`)

- `cmd/main/main_wbtest.mbt` — T1–T58 (core/physics/preconditioner)
- `cmd/main/main_ext_wbtest.mbt` — T59–T154 (extended suite)
- `cmd/main/validation_cli_wbtest.mbt` — 3 (CLI filtering, checkpoint and time-step constraints)
- `cmd/main/solver_projection_tvd_wbtest.mbt` — 19 (Projection-TVD extension, Section 6)

See [docs/test_report.md](docs/test_report.md) for the recorded per-test report and
[data/paper/](data/paper/) for the paper's quantitative data.

---

## Runtime & Build Environment

| Tool/Env         | Minimum Version | Install/Usage |
|------------------|-----------------|--------------|
| moon (MoonBit)   | 0.1.20260309    | [Official website](https://www.moonbitlang.com/download/) |
| moonc (compiler) | 0.8.3           | bundled      |
| node             | 18+             | [Node.js official](https://nodejs.org/) (export check in build_wasm.sh) |
| python           | 3.x             | (for python http.server) |
| Chrome/Edge/Brave| 115+            | Browser support for wasm-gc stringref |
| Any HTTP server  |                 | Recommended: python |

### Starting Local Server

1. **Python**:
	```bash
	python -m http.server 8080
	# Open http://localhost:8080/cmd/main/main.html
	```


---

## docs/ Directory Reference

| File | Description |
|---|---|
| [arch.md](docs/arch.md) | System architecture layers, modules, data flow |
| [api_reference.md](docs/api_reference.md) | All WASM exported function signatures |
| [dev_guide.md](docs/dev_guide.md) | Development environment, AI-assisted workflow |
| [ghia_validation.md](docs/ghia_validation.md) | Ghia (1982) numerical benchmark validation |
| [test_report.md](docs/test_report.md) | Two-dimensional test-suite report |
| [launch.md](docs/launch.md) | Browser build/launch guide + paper reproduction commands |
| [preconditioner_theory.md](docs/preconditioner_theory.md) | Preconditioner mathematical theory |
| [preconditioner_plan.md](docs/preconditioner_plan.md) | Preconditioner implementation plan |
| [flow.md](docs/flow.md) | Execution flow diagrams (Mermaid) |
| [validation_report.md](docs/validation_report.md) | Cross-version consistency validation |
| [design_pcg_solver.md](docs/design_pcg_solver.md) | PCG solver design notes |
| [cfd_terminology.md](docs/cfd_terminology.md) | CFD terminology glossary |
| [SKILL.md](docs/SKILL.md) | Common bugs and architecture decisions |
| [solver_monitor.md](docs/solver_monitor.md) | Monitor quantities definition |
| [theory.md](docs/theory.md) | Theoretical background and derivations |

---

### Key source file: `cmd/main/main.mbt`

**Chorin solver**

| Section | Description |
|---|---|
| `init_simulation()` | Zeroes Chorin state (`g_u`, `g_v`, `g_p`), resets step counter |
| `run_n_steps(n)` | Advances Chorin by *n* time steps |
| `cavity_flow_array()` | Core Navier-Stokes loop (explicit projection, finite-difference) |
| `pressure_poisson_array()` | Gauss-Seidel pressure solve (`nit` iterations) |
| `get_u/v/p_at(i,j)` | Point access to Chorin fields |
| `get_divergence_norm()` | Mean \|div u\| over interior cells |

**SIMPLE solver**

| Section | Description |
|---|---|
| `init_simple()` | Zeroes SIMPLE state (`g_u_s`, `g_v_s`, `g_p_s`) |
| `run_simple_n_iter(n)` | Runs *n* SIMPLE iterations with under-relaxation |
| `simple_one_iter()` | One SIMPLE sweep: predictor → pressure correction → update |
| `get_u/v/p_simple_at(i,j)` | Point access to SIMPLE fields |

**Chorin-PCG solver**

| Section | Description |
|---|---|
| `init_chorin_pcg()` | Zeroes PCG state (`g_u_pcg`, `g_v_pcg`, `g_p_pcg`) |
| `run_chorin_pcg_n_steps(n)` | Advances PCG solver by *n* time steps |
| `pressure_poisson_pcg()` | PCG pressure solve with Jacobi preconditioner |
| `laplacian_apply()` | Matrix-free 5-point Laplacian (collocated grid) |
| `get_u/v/p_pcg_at(i,j)` | Point access to PCG fields |

**MAC staggered-grid solver**

| Section | Description |
|---|---|
| `init_mac()` | Zeroes MAC state (`g_u_mac`, `g_v_mac`, `g_p_mac`) |
| `run_mac_n_steps(n)` | Advances MAC solver by *n* time steps |
| `cavity_flow_mac()` | MAC main loop: predictor → divergence → PCG → correction |
| `mac_u_predictor()` / `mac_v_predictor()` | Advection-diffusion step with ghost-cell BCs |
| `pressure_poisson_pcg_mac()` | PCG pressure solve on `mac_nc × mac_nc` grid |
| `mac_correct_velocity()` | Velocity correction: u -= dt/ρ/dx · ∂p/∂x |
| `apply_pressure_bcs_mac()` | Dirichlet p=0 at top lid; Neumann on other walls |
| `get_u/v/p_mac_at(i,j)` | Point access to MAC fields |
| `get_mac_divergence_norm()` | Mean \|div u\| over PCG-interior cells (i,j = 1..nc-2) |

---

## Prerequisites

| Tool | Minimum version | Install |
|---|---|---|
| `moon` (MoonBit build tool) | 0.1.20260309 | https://www.moonbitlang.com/download/ |
| `moonc` (MoonBit compiler) | 0.8.3 | bundled with moon |
| `node` (for export verification only) | 18+ | https://nodejs.org/ |
| Chrome / Edge / Brave | 115+ | required to run wasm-gc with stringref |
| Any HTTP server | — | see "Serving the page" below |

---

## Quick Start — From Source to Browser

### 1. Install dependencies

```bash
moon version
moonc -v
node --version   # 18+ required
```

### 2. Clone and enter the project

```bash
git clone https://github.com/YuFenfen-ai/FlowLabLite.git
cd FlowLabLite
```

### 3. Run the tests

```bash
moon test cmd/main --target wasm    # runs all 176 tests
```

Expected output:
```
Total tests: 176, passed: 176, failed: 0
```

### 4. Build the WASM binary

```bash
bash build_wasm.sh release    # optimised build
bash build_wasm.sh            # debug build
```

The script prints the export list (121 functions + `_start` + `memory` = 123 symbols).

### 5. Serve over HTTP

```bash
python -m http.server 8080
```

### 6. Open in Chrome 115+

```
http://localhost:8080/cmd/main/main.html
```

---

## How to Build

### Step 1 – Compile

```bash
moon build --target wasm-gc          # debug build
moon build --target wasm-gc --release  # release build
```

> **Known issue in moon 0.1.20260309:** `moon build` alone does not export
> the API functions. Use `build_wasm.sh` to re-link with full exports.

### Step 2 – Re-link with exported functions (required)

```bash
bash build_wasm.sh          # debug WASM
bash build_wasm.sh release  # release WASM
```

### Step 3 – Run the tests

```bash
moon test --target wasm
```

---

## How to Use – Local Run + Browser Viewer

### Step 1 – Configure the run (optional)

Edit the constants at the top of `cmd/main/main.mbt`:

| Constant | Default | Description |
|---|---|---|
| `timing_enabled` | `true` | Master timing switch |
| `run_chorin` | `true` | Run Chorin solver |
| `local_nt` | `nt` (500) | Chorin time steps |
| `run_simple` | `true` | Run SIMPLE solver |
| `local_simple_n` | `100` | SIMPLE iterations |
| `run_pcg` | `true` | Run Chorin-PCG solver |
| `local_pcg_nt` | `nt` (500) | PCG time steps |
| `run_mac` | `true` | Run MAC solver |
| `local_mac_nt` | `nt` (500) | MAC time steps |

### Step 2 – Run the simulation locally

```bash
bash run_local.sh                  # → results.json (all four solvers)
bash run_local.sh my_results.json  # → custom filename
```

### Step 3 – Open the local viewer

Open `cmd/main/local_viewer.html` in a browser.

### Step 4 – Load results

Drag & drop `data/sample_results.json` onto the page. Solver tabs **Chorin / SIMPLE / PCG / MAC**
appear when the JSON contains multiple solvers. Clicking a tab switches all plots
and statistics to that solver. For MAC, the grid size shown is 40×40 (cell centres).

---

## Solver Algorithms

### Chorin Projection Method (1968)

Explicit fractional-step time-marching:
1. Solve momentum equations explicitly for intermediate velocity u*
2. Solve pressure Poisson equation: ∇²p = ρ/dt · div(u*)  (Gauss-Seidel, `nit` iterations)
3. Correct velocity: u = u* − (dt/ρ) · ∇p
4. Apply boundary conditions

### SIMPLE (Patankar & Spalding 1972)

Steady-state pressure-correction iteration with under-relaxation:
1. Momentum predictor (explicit, under-relaxation α_u = 0.7)
2. Build pressure-correction source: b = −∇·u*
3. Solve ∇²p' = b (Gauss-Seidel, 50 inner iterations)
4. Update: p ← p* + α_p · p'  (α_p = 0.3); u/v corrected

### Chorin-PCG

Same fractional-step structure as Chorin but the pressure Poisson solve uses
**Preconditioned Conjugate Gradient** (PCG) with Jacobi preconditioner instead
of Gauss-Seidel:
- Matrix-free 5-point Laplacian (`laplacian_apply`)
- Jacobi preconditioner: M⁻¹r = r / (2/dx² + 2/dy²)
- Tolerance: 1×10⁻⁶, max iterations: 200
- Boundary conditions: Dirichlet p=0 at top, Neumann elsewhere

### MAC Staggered Grid (Harlow & Welch 1965)

Staggered arrangement with PCG pressure solve:
- **p** at cell centres — mac_nc × mac_nc (40×40)
- **u** at x-face centres — mac_nc × (mac_nc+1)
- **v** at y-face centres — (mac_nc+1) × mac_nc
- Ghost-cell BCs: top lid u = 2·U_lid − u[ny-1][j] (linear interpolation), no-slip on other walls
- Divergence-free guarantee: PCG-interior cells (i,j = 1..38) have exact ∇·u = 0 by construction

---

## WASM Exports (121 functions)

**Chorin solver (27):**
`init_simulation`, `run_all_steps`, `run_n_steps`, `get_nx`, `get_ny`, `get_nt`, `get_nit`, `get_re`, `get_dx`, `get_dy`, `get_dt`, `get_rho`, `get_nu`, `get_step_count`, `get_u_at`, `get_velocity_magnitude_at`, `get_v_at`, `get_p_at`, `get_u_center`, `get_v_center`, `get_p_center`, `get_divergence_norm`, `get_max_velocity_magnitude`, `get_max_u`, `get_max_v`, `get_max_p`, `get_min_p`

**SIMPLE solver (9):**
`init_simple`, `run_simple_n_iter`, `get_simple_step_count`, `get_simple_residual`, `get_u_simple_at`, `get_v_simple_at`, `get_p_simple_at`, `get_max_u_simple`, `get_simple_divergence_norm`

**Chorin-PCG collocated solver (13):**
`init_chorin_pcg`, `run_chorin_pcg_n_steps`, `get_pcg_step_count`, `get_pcg_last_iters`, `get_u_pcg_at`, `get_v_pcg_at`, `get_p_pcg_at`, `get_velocity_magnitude_pcg_at`, `get_max_u_pcg`, `get_max_v_pcg`, `get_max_p_pcg`, `get_min_p_pcg`, `get_pcg_divergence_norm`

**MAC staggered solver (14):**
`init_mac`, `run_mac_n_steps`, `get_mac_step_count`, `get_mac_last_iters`, `get_mac_nc`, `get_u_mac_at`, `get_v_mac_at`, `get_p_mac_at`, `get_velocity_magnitude_mac_at`, `get_max_u_mac`, `get_max_v_mac`, `get_max_p_mac`, `get_min_p_mac`, `get_mac_divergence_norm`

**RK3 solver (7):**
`init_rk3`, `run_rk3_n_steps`, `get_rk3_step_count`, `get_u_rk3_at`, `get_v_rk3_at`, `get_p_rk3_at`, `get_rk3_divergence_norm`

**Finite-volume verification (7):**
`init_fvm`, `run_fvm_n_iter`, `get_fvm_step_count`, `get_u_fvm_at`, `get_v_fvm_at`, `get_p_fvm_at`, `get_fvm_divergence_norm`

**Scalar transport (7):**
`init_scalar`, `set_scalar_alpha`, `get_scalar_alpha`, `run_scalar_n_steps`, `get_scalar_step_count`, `get_scalar_at`, `get_scalar_mean`

**Boussinesq natural convection (14):**
`init_boussinesq`, `set_bous_nu`, `set_bous_alpha`, `set_bous_beta`, `get_bous_nu`, `get_bous_alpha`, `get_bous_beta`, `run_boussinesq_n_steps`, `get_bous_step_count`, `get_u_bous_at`, `get_v_bous_at`, `get_p_bous_at`, `get_phi_bous_at`, `get_bous_nusselt`

**Channel flow (10):**
`init_channel`, `set_channel_fx`, `set_channel_nu`, `get_channel_fx`, `run_channel_n_steps`, `get_channel_step_count`, `get_u_channel_at`, `get_v_channel_at`, `get_p_channel_at`, `get_channel_analytic_u`

**Grid utilities & curvilinear coordinates (13):**
`set_grid_size`, `get_tanh_y_at`, `get_tanh_dy_at`, `get_tanh_mean_dy`, `set_curv_c`, `get_curv_c`, `get_curv_x_at`, `get_curv_y_at`, `get_curv_jacobian_at`, `get_curv_nonortho_at`, `get_curv_g11_at`, `get_curv_g12_at`, `get_curv_g22_at`

---

## Test Report

```bash
moon test cmd/main --target wasm      # 176 tests on the wasm target
```

All 176 tests pass in the recorded run (see `TEST_RESULTS.log`):
- `cmd/main/main_wbtest.mbt` — T1–T58 (core solvers + preconditioners)
- `cmd/main/main_ext_wbtest.mbt` — T59–T154 (extended suite)
- `cmd/main/validation_cli_wbtest.mbt` — 3 (CLI constraints, Sections 4.3-4.4)
- `cmd/main/solver_projection_tvd_wbtest.mbt` — 19 (Projection-TVD extension, Section 6)

The per-test listing of the original suite is available in
[docs/test_report.md](docs/test_report.md); the paper's quantitative
benchmarks are in [data/paper/](data/paper/).

### Physical validation

**Chorin solver (Re = 20, 500 time steps):**
- Centre u-velocity ≈ −0.06 (backflow confirms clockwise primary vortex)
- Max u-velocity ≥ 1.0 (lid velocity maintained)
- Divergence norm decreases to O(10⁻²)

**SIMPLE solver (Re = 20):**
- Negative u at centre after 200 iterations confirms vortex
- Residual ‖div u‖ non-negative and decreasing
- Under-relaxation keeps max \|u\| bounded within [0, 2]

**Chorin-PCG solver (Re = 20):**
- PCG pressure convergence < 1×10⁻⁶ in typically 50–130 iterations per step
- Interior divergence norm < 1×10⁻⁴ after 50 steps
- Vortex structure (negative u at centre) confirmed after 200 steps

**MAC staggered solver (Re = 20):**
- PCG-interior cells (i,j = 1..38) guaranteed divergence-free by construction
- Divergence norm < 1×10⁻⁴ after 50 steps
- Vortex confirmed after 200 steps
- PCG iteration count decreases as flow reaches quasi-steady state (~120 → ~50 iters/step)

---

## Simulation Parameters

| Parameter | Value | Description |
|---|---|---|
| `nx`, `ny` | 41 × 41 | Grid nodes in x and y (Chorin / SIMPLE / PCG) |
| `mac_nc` | 40 | MAC cell count per direction (= nx − 1) |
| `nt` | 500 | Chorin / PCG / MAC time steps |
| `nit` | 50 | Pressure-Poisson inner iterations (Chorin / SIMPLE) |
| `dx`, `dy` | 0.05 | Grid spacing (domain = 2 × 2) |
| `dt` | 0.001 | Time step |
| `rho` | 1.0 | Fluid density |
| `nu` | 0.1 | Kinematic viscosity |
| Re | 20 | Reynolds number (u_lid × L / nu = 1 × 2 / 0.1) |
| `simple_alpha_p` | 0.3 | SIMPLE pressure under-relaxation factor |
| `simple_alpha_u` | 0.7 | SIMPLE velocity under-relaxation factor |
| `pcg_tol` | 1×10⁻⁶ | PCG convergence tolerance |
| `pcg_max_iter` | 200 | PCG maximum iterations per step |

---

## License

Apache-2.0 — see [LICENSE](LICENSE).

## How to cite

If you use FlowLabLite in your work, please cite the thesis that documents
this software:

> 余芬芬. 基于MoonBit的CFD演示原型AI原生构建方法[D]. 北京: 北京航空航天大学, 2026.

```bibtex
@misc{yu2026flowlablite,
  author       = {Yu, Fenfen},
  title        = {An AI-native construction method for a MoonBit-based CFD
                  demonstration prototype},
  howpublished = {Bachelor's thesis, School of Continuing Education,
                  Beihang University, Beijing},
  year         = {2026},
  note         = {In Chinese}
}
```

Developed by Fenfen Yu (余芬芬), School of Continuing Education, Beihang
University (北京航空航天大学), in collaboration with Ezhou Hi-Modeling
Technology Co., Ltd. (鄂州海慕科技有限公司).

---

## Paper artifacts and release structure

This repository is the submission artifact for the manuscript
"Theory-Code-Test Co-Design for AI-Assisted Scientific Software"
(IEEE Computing in Science & Engineering draft). Three tagged states:

| Tag | Contents | Non-test lines | Tests | Exports |
|---|---|---|---|---|
| `2d-baseline` | Four-solver 2D artifact (Sections 3, 4.1-4.2) | 5,843 | 156 | 121 |
| `2d-projection-tvd-v1` | + Projection-TVD extension (Section 6) | 6,920 | 175 | 121 |
| `main` (this tree) | + validation harness & paper data (Sections 4.3-4.6) | 6,982 | 176 | 121 |

- `data/paper/` — curated paper data package. Its README maps every
  table and figure in Sections 4.3-4.6 to the exact data file, and
  `verify_paper_numbers.py` re-checks every published number against
  the raw data (74/74 checks; see `verification_report.txt`).
- `data/ghia_re100_u.csv`, `data/ghia_re100_v.csv` — the 17-point
  Ghia (1982) Re = 100 comparison with reference values and errors
  (Table 3, Re = 100 row).
- `TEST_RESULTS.log` — full recorded output of
  `moon test cmd/main --target wasm` on this tree
  (176 passed / 0 failed).
- MoonBit toolchain: the baseline artifact was built with moon
  0.1.20260309; the validation-era runs with moon 0.1.20260618.
  Browser build and launch steps: `docs/launch.md`.
- Zenodo snapshot DOI: to be added upon archival.
