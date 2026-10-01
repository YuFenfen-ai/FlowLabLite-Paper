# FlowLabLite — System Architecture

**Version**: HEAD (`ba1505e`)  
**Author**: Fenfen Yu (余芬芬)  
**Date**: 2026-04-17

---

## 1. Three-Layer Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│  Layer 3 — User Layer (Browser)                                 │
│                                                                 │
│  ┌──────────────────────┐   ┌──────────────────────────────┐    │ 
│  │  • Parameter display │   │  • Drag-drop JSON loader     │    │
│  │  • Run / Reset btns  │   │  • 4-solver tab selector     │    │
│  │  • Velocity heatmap  │   │  • Heatmaps + statistics     │    │
│  │  • Pressure heatmap  │   │                              │    │
│  │  • Streamline overlay│   │                              │    │
│  └──────────┬───────────┘   └──────────────────────────────┘    │
│             │  WebAssembly.instantiate() + JS getter calls      │
└─────────────┼───────────────────────────────────────────────────┘
              │
┌─────────────┼───────────────────────────────────────────────────┐
│  Layer 2 — Application Layer (MoonBit → WASM-GC)                │
│                                                                 │
│  ┌──────────▼──────────────────────────────────────────────┐    │
│  │  WASM Bridge (moon.pkg.json exports, build_wasm.sh)     │    │
│  │  121 exported functions — all solvers, monitors, getters│    │
│  └──────────┬──────────────────────────────────────────────┘    │
│             │                                                   │
│  ┌──────────▼──────────────────────────────────────────────┐    │
│  │  Solver Main Control (cmd/main/main.mbt)                │    │
│  │                                                         │    │
│  │  ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌────────┐   │    │
│  │  │  Chorin  │  │  SIMPLE  │  │Chorin-PCG│  │  MAC   │   │    │
│  │  │ 500 steps│  │ 100 iter │  │ 500 steps│  │500 step│   │    │
│  │  │ GS press │  │ GS press │  │PCG press │  │PCG+stag│   │    │
│  │  └────┬─────┘  └────┬─────┘  └────┬─────┘  └───┬────┘   │    │
│  │       └─────────────┴──────────────┴─────────────┘      │    │
│  │                         │                               │    │
│  │              ┌──────────▼──────────┐                    │    │
│  │              │  Shared Modules     │                    │    │
│  │              │  • laplacian_apply  │                    │    │
│  │              │  • apply_pressure_bcs│                   │    │
│  │              │  • PCG solver       │                    │    │
│  │              │  • DILU/DIC/GAMG    │                    │    │
│  │              │    preconditioners  │                    │    │
│  │              └─────────────────────┘                    │    │
│  └─────────────────────────────────────────────────────────┘    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
              │
┌─────────────┼───────────────────────────────────────────────────┐
│  Layer 1 — Foundation Layer (MoonBit stdlib)                    │
│                                                                 │
│  Array[T], Double arithmetic, Int, Bool — no external deps      │
└─────────────────────────────────────────────────────────────────┘
```

---

## 2. Core Module Responsibilities

### Module 1 — Frontend Interaction (`cmd/main/main.html`)

| Responsibility | Implementation |
|---|---|
| Trigger computation | JS calls `init_*()` then `run_*_n_steps(n)` via WASM |
| Display parameters | Reads `get_nx/ny/re/nu/dt/dx/dy` from WASM on load |
| Velocity heatmap | Reads `get_velocity_magnitude_at(i,j)` for all grid points, maps to jet colormap |
| Pressure heatmap | Reads `get_p_at(i,j)`, maps to blue→red diverging colormap |
| **Streamline overlay** | RK4 integration of WASM velocity field; drawn on Canvas as curves |
| Statistics panel | Reads `get_max_u`, `get_divergence_norm`, etc. |

**Key design**: the page is stateless — every "Run" press calls `init_*()` to zero state
before advancing, so results are always reproducible.

---

### Module 2 — WASM Bridge (`cmd/main/moon.pkg.json` + `build_wasm.sh`)

| Responsibility | Implementation |
|---|---|
| Export function list | `"exports": [...]` in `moon.pkg.json` (121 functions) |
| Re-linking with exports | `build_wasm.sh` invokes `moonc link-core -exported_functions ...` |
| Import object | JS `Proxy`-based import satisfying `spectest/env/moonbit:ffi` |

**Known limitation**: `moon build` alone (v0.1.20260309) does not propagate
the export list to the linker. `build_wasm.sh` is the required build path.

---

### Module 3 — Solver Main Control (`cmd/main/main.mbt`, lines 1–200)

Four independent solvers, each with its own global state arrays and WASM API:

| Solver | Global state | Init function | Run function |
|---|---|---|---|
| Chorin | `g_u`, `g_v`, `g_p`, `g_steps` | `init_simulation()` | `run_n_steps(n)` |
| SIMPLE | `g_u_s`, `g_v_s`, `g_p_s`, `g_simple_steps` | `init_simple()` | `run_simple_n_iter(n)` |
| Chorin-PCG | `g_u_pcg`, `g_v_pcg`, `g_p_pcg`, `g_pcg_steps` | `init_chorin_pcg()` | `run_chorin_pcg_n_steps(n)` |
| MAC | `g_u_mac`, `g_v_mac`, `g_p_mac`, `g_mac_steps` | `init_mac()` | `run_mac_n_steps(n)` |

The four state arrays are completely independent; running one solver never modifies another.

---

### Module 4 — Physics Model (`cmd/main/main.mbt`, `cavity_flow_array`, `simple_one_iter`, `cavity_flow_mac`)

Implements the Navier-Stokes equations in discrete form:

| Solver | Convection scheme | Diffusion scheme | Pressure solve |
|---|---|---|---|
| Chorin | 1st-order upwind (explicit) | Central difference (explicit) | Gauss-Seidel (50 iter) |
| SIMPLE | 1st-order upwind (explicit) | Central difference | Gauss-Seidel (50 iter) |
| Chorin-PCG | 1st-order upwind | Central difference | PCG + Jacobi precond. |
| MAC | 1st-order upwind (staggered) | Central difference (staggered) | PCG + Jacobi precond. |

Boundary conditions:
- Velocity: side walls no-slip → lid U = 1 (order matters — see `CLAUDE.md`)
- Pressure: Neumann `∂p/∂n = 0` on walls, Dirichlet `p = 0` at top lid

---

### Module 5 — Linear Algebra Solver (`cmd/main/main.mbt`, `pressure_poisson_pcg*`)

The PCG framework with four preconditioner options:

```
pressure_poisson_pcg_prec(p, dx, dy, b, precond_type)
  precond_type = 0 → Jacobi (diagonal scaling)
  precond_type = 1 → DILU  (modified incomplete LU)
  precond_type = 2 → DIC   (IC(0) incomplete Cholesky)
  precond_type = 3 → GAMG  (2-level geometric multigrid)
```

The GAMG outer loop uses **stationary Richardson iteration** (`p -= z`)
rather than standard PCG — see `docs/preconditioner_theory.md` for the
mathematical derivation.

---

### Module 6 — Grid and Field Management (`cmd/main/main.mbt`, utilities)

| Function | Purpose |
|---|---|
| `create_zeros_2d(rows, cols)` | Allocate 2D array with independent rows |
| `copy_2d_array(src)` | Deep copy (each row individually) |
| `generate_mesh_grid()` | Return (X, Y) coordinate grids |
| Global `let` arrays | 41×41 (or 40×40 MAC) Double arrays |

Array convention: `a[i][j]` where `i` = row (y-direction, 0 = bottom),
`j` = column (x-direction, 0 = left).

---

### Module 7 — Validation and Output (`cmd/main/main.mbt`, `run_local.sh`)

| Component | Purpose |
|---|---|
| `get_divergence_norm()` | Mean |div u| over interior cells — mass conservation check |
| `get_max_u/v/p/min_p` | Field statistics for convergence monitoring |
| JSON output (via `moon run`) | Full field data serialised for `local_viewer.html` |
| `run_local.sh` | Orchestrates run + JSON extraction |
| `docs/ghia_validation.md` | Centreline data vs Ghia 1982 benchmark |

---

## 3. Data Flow

### 3.1 Browser real-time visualisation

```
User opens main.html
  → JS loads WASM binary (fetch + WebAssembly.instantiate)
  → JS reads config: get_nx/ny/re/nu/dt from WASM
  → User clicks "Run N Steps"
  → JS calls init_*()/run_*_n_steps(n) on WASM instance
  → JS reads get_*_at(i,j) for all grid points  [O(nx×ny) calls]
  → JS maps values to RGB via colormap
  → JS draws Canvas 2D heatmap
  → JS integrates streamlines via RK4 on velocity field
  → JS draws streamline curves on Canvas overlay
```

### 3.2 Local JSON workflow

```
bash run_local.sh
  → moon build --target wasm (compile)
  → wasmtime run main.wasm (execute all 4 solvers)
  → stdout contains JSON block
  → shell script extracts JSON → results.json

User opens local_viewer.html
  → Drag-drop results.json
  → JS parses JSON, populates 4-solver tabs
  → Renders heatmaps and statistics
```

### 3.3 Test execution

```
moon test --target wasm
  → Discovers main_wbtest.mbt (T1–T58), main_ext_wbtest.mbt (T59–T154),
  validation_cli_wbtest.mbt and solver_projection_tvd_wbtest.mbt
  → Compiles to wasm-gc and runs all test functions
  → Reports 176/176 passed
```

---

## 4. Key Constraints and Design Decisions

| Decision | Rationale | Reference |
|---|---|---|
| Global mutable state via `Array[T]` | MoonBit disallows `let mut` at global scope | `CLAUDE.md` §7 |
| `build_wasm.sh` for exports | moon 0.1.20260309 linker bug | `CLAUDE.md` Known Issues |
| BC order: walls → lid | Lid BC must be applied last to preserve corner values | `CLAUDE.md` §1 |
| GAMG uses Richardson, not PCG | CG-based coarse solver is non-linear → breaks PCG orthogonality | `docs/preconditioner_theory.md` |
| SND Laplacian (diagonal < 0) | Standard finite-difference convention; p\* < 0 for b > 0 | `CLAUDE.md` §1 |
| Separate test files | T1–T58 (functional), T59–T154 (engineering), CLI + extension tests | `docs/test_report.md` |

---

## 5. Technology Stack

| Component | Technology | Version |
|---|---|---|
| Solver language | MoonBit | 0.8.3 |
| Build tool | moon | 0.1.20260309 |
| Compilation target | wasm-gc (WebAssembly GC) | — |
| Browser runtime | Chrome/Edge/Brave | 115+ (wasm-gc stringref) |
| Test runtime | Wasmtime (via moon test) | bundled |
| Local run runtime | Wasmtime | bundled |
| Visualisation | Canvas 2D API | browser native |
| No external JS dependencies | — | — |
