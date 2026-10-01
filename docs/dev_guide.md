# FlowLabLite — Developer Guide

**Author**: Fenfen Yu (余芬芬)  
**Date**: 2026-04-17

---

## 1. Development Environment Setup

### 1.1 Required Tools

| Tool | Version | Install |
|---|---|---|
| `moon` (MoonBit build tool) | ≥ 0.1.20260309 | https://www.moonbitlang.com/download/ |
| `moonc` (compiler, bundled) | ≥ 0.8.3 | bundled with moon |
| `wasmtime` (test runtime) | any | bundled with moon |
| `node` (export verification) | ≥ 18 | https://nodejs.org/ |
| Chrome / Edge / Brave | ≥ 115 | wasm-gc stringref support |

### 1.2 Verify Installation

```bash
moon version      # e.g. 0.1.20260309
moonc -v          # e.g. v0.8.3+...
node --version    # e.g. v20.x
```

### 1.3 Clone and Build

```bash
git clone https://github.com/YuFenfen-ai/FlowLabLite.git
cd FlowLabLite
moon test --target wasm          # 165 tests should all pass
bash build_wasm.sh release       # build the production WASM
```

---

## 2. Project Layout Walkthrough

```
cmd/main/main.mbt          ← all solver code (~3000+ lines)
cmd/main/main_wbtest.mbt   ← tests T1–T58 (core solvers)
cmd/main/main_ext_wbtest.mbt ← tests T59–T152 (extended features/3D/advanced)
cmd/main3d/main3d_wbtest.mbt ← tests 3D solvers (13 basic tests)
cmd/main/moon.pkg.json     ← WASM exports list (121 functions)
docs/arch.md               ← architecture overview
docs/api_reference.md      ← WASM function signatures
```

Code is organised in layers within `main.mbt`:
1. Global constants (`nx`, `ny`, `dt`, `nu`, …)
2. Global state arrays per solver (`g_u`, `g_u_s`, `g_u_pcg`, `g_u_mac`, …)
3. Utility functions (`create_zeros_2d`, `copy_2d_array`)
4. Chorin solver (GS pressure)
5. SIMPLE solver
6. Chorin-PCG solver + DILU/DIC/GAMG preconditioners
7. MAC staggered solver
8. RK3, FVM, Scalar, Boussinesq, Channel solvers
9. WASM exports

---

## 3. Running Tests

```bash
# Run all 165 tests on the wasm target (required — native target differs)
moon test --target wasm

# Run a single named test (prefix match)
moon test --target wasm --filter dilu_zero_rhs

# See verbose output
moon test --target wasm --verbose
```

Test files are auto-discovered: any `*_wbtest.mbt` in the same package is included.

### Test categories (T59–T83 in `main_ext_wbtest.mbt`)

| Category | Tests | What it checks |
|---|---|---|
| Unit — array utilities | T59–T60 | Deep-copy independence |
| Unit — Laplacian | T61–T63 | Polynomial exactness, harmonic zero |
| Unit — BC | T64–T65 | Idempotency |
| Unit — RHS | T66–T67 | Source-term formula |
| Unit — GAMG parts | T68–T70 | Prolongate, restrict /4, smoother descent |
| Integration — PCG residual | T71–T74 | All 4 preconditioners converge |
| Integration — cross-solver | T75 | < 1% difference between all preconditioners |
| System — physics | T76–T79 | SND sign, determinism, additivity |
| Regression | T80–T83 | Numerical stability properties |

---

## 4. How to Add a New Solver Module

Suppose you want to add a new **Runge-Kutta 4 time integration** variant.

### Step 1 — Add global state

```moonbit
// After the last existing solver's state (search "[Layer RK4]"):
let g_u_rk4 : Array[Array[Double]] = create_zeros_2d(ny, nx)
let g_v_rk4 : Array[Array[Double]] = create_zeros_2d(ny, nx)
let g_p_rk4 : Array[Array[Double]] = create_zeros_2d(ny, nx)
let g_rk4_steps : Array[Int] = [0]
```

### Step 2 — Implement the solver function

```moonbit
fn cavity_flow_rk4(u, v, p) -> (Array[Array[Double]], Array[Array[Double]], Array[Array[Double]]) {
  // RK4 stages: k1, k2, k3, k4
  // ... your implementation
}
```

### Step 3 — Add WASM lifecycle functions

```moonbit
pub fn init_rk4() -> Unit {
  g_rk4_steps[0] = 0
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      g_u_rk4[i][j] = 0.0; g_v_rk4[i][j] = 0.0; g_p_rk4[i][j] = 0.0
    }
  }
}

pub fn run_rk4_n_steps(n : Int) -> Unit {
  for step = 0; step < n; step = step + 1 {
    let (u_new, v_new, p_new) = cavity_flow_rk4(g_u_rk4, g_v_rk4, g_p_rk4)
    for i = 0; i < ny; i = i + 1 {
      for j = 0; j < nx; j = j + 1 {
        g_u_rk4[i][j] = u_new[i][j]
        g_v_rk4[i][j] = v_new[i][j]
        g_p_rk4[i][j] = p_new[i][j]
      }
    }
  }
  g_rk4_steps[0] = g_rk4_steps[0] + n
}

pub fn get_u_rk4_at(i : Int, j : Int) -> Double {
  if i < 0 || i >= ny || j < 0 || j >= nx { 0.0 } else { g_u_rk4[i][j] }
}
```

### Step 4 — Register WASM exports

In `cmd/main/moon.pkg.json`, add to the `exports` array:
```json
"init_rk4",
"run_rk4_n_steps",
"get_u_rk4_at",
"get_v_rk4_at",
"get_p_rk4_at"
```

Also add the same names to the `EXPORTS` variable in `build_wasm.sh`.

### Step 5 — Write tests first (TDD)

```moonbit
test "rk4_init_resets_state" {
  init_rk4()
  assert_eq(get_rk4_step_count(), 0)
  assert_eq(get_u_rk4_at(20, 20), 0.0)
}

test "rk4_boundary_conditions" {
  init_rk4()
  run_rk4_n_steps(10)
  // lid u = 1
  for j = 0; j < nx; j = j + 1 {
    assert_true((get_u_rk4_at(ny - 1, j) - 1.0).abs() < 1.0e-10)
  }
}
```

Run: `moon test --target wasm --filter rk4`

---

## 5. AI-Assisted Development Workflow (Reproducible Example)

This section documents the exact AI-assisted workflow used to implement the
`build_modified_diag` function (shared by DILU and DIC preconditioners).
The workflow is **fully reproducible** — anyone can follow these steps.

### 5.1 Context: what was needed

The DILU preconditioner requires a "modified diagonal" array where each
interior node's diagonal coefficient is reduced by contributions from its
upstream neighbours:

```
d[i][j] = a_diag - (inv_dx²)² / d[i][j-1] - (inv_dy²)² / d[i-1][j]
```
with fallback: if `d[i][j] ≤ 1e-15`, reset to `a_diag` (avoid division by zero).

### 5.2 Prompt used (verbatim)

```
Context:
- Language: MoonBit (similar syntax to Rust/Swift, uses fn, let, Array[T])
- Grid: nx=41, ny=41, domain 2×2, dx=dy=0.05, uniform
- The discrete 5-point Laplacian has diagonal coefficient a_diag = 2/dx² + 2/dy² = 1600
- inv_dx2 = 1/dx² = 400.0, inv_dy2 = 1/dy² = 400.0

Task: Implement the MoonBit function build_modified_diag that computes the
modified-diagonal array for the DILU/DIC preconditioner. The function should:
1. Take inv_dx2 (Double), inv_dy2 (Double), a_diag (Double) as inputs
2. Return Array[Array[Double]] of size ny × nx, all zeros on boundary rows/cols
3. For interior nodes (i in 1..ny-2, j in 1..nx-2), compute:
   d[i][j] = a_diag - inv_dx2*inv_dx2/d[i][j-1] - inv_dy2*inv_dy2/d[i-1][j]
4. If d[i][j] <= 1e-15, set it back to a_diag (numerical safety)
5. Do NOT use let mut (MoonBit global scope restriction)
6. Use row-major scan order (i outer, j inner)

Generate the complete MoonBit function.
```

### 5.3 AI-generated code (lightly edited)

The following was generated by Claude Sonnet 4.6 in response to the above prompt:

```moonbit
fn build_modified_diag(
  inv_dx2 : Double,
  inv_dy2 : Double,
  a_diag  : Double,
) -> Array[Array[Double]] {
  let d = create_zeros_2d(ny, nx)
  // Boundary rows/cols remain 0 (ghost values, not used in apply step)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let contrib_x = inv_dx2 * inv_dx2 / d[i][j - 1]   // left neighbour
      let contrib_y = inv_dy2 * inv_dy2 / d[i - 1][j]   // bottom neighbour
      // Guard against zero (boundary cells or first row/col of interior)
      let cx = if d[i][j - 1] > 1.0e-15 { contrib_x } else { 0.0 }
      let cy = if d[i - 1][j] > 1.0e-15 { contrib_y } else { 0.0 }
      let dval = a_diag - cx - cy
      d[i][j] = if dval > 1.0e-15 { dval } else { a_diag }
    }
  }
  d
}
```

### 5.4 Validation (unit test written immediately after)

```moonbit
test "dilu_modified_diag_positive" {
  // T46: all modified diagonal values must be strictly positive
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let d = build_modified_diag(inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(d[i][j] > 0.0)
    }
  }
}
```

Test result: **PASS** (T46 in the 83-test suite).

### 5.5 Lessons from this workflow

| Principle | Application |
|---|---|
| **Precise context in prompt** | Specified language syntax constraints (`let mut` restriction), exact formula, numerical safety requirement |
| **Test immediately** | Wrote T46 before integrating the function into the PCG solver |
| **Minimal edits** | Only edited variable names for consistency with existing code; logic was correct first try |
| **AI as "senior autocomplete"** | AI generated the loop structure and guard logic; domain knowledge (the formula itself) came from the mathematical derivation |

---

## 6. AI-Assisted Code Template: QUICK Convection Scheme

This section shows how AI was used to generate a code template for the
**QUICK (Quadratic Upstream Interpolation for Convective Kinematics)** scheme
(Leonard, 1979), a 3rd-order upwind method for the convective terms.

### 6.1 Why QUICK?

Current solvers use 1st-order upwind differencing for convective terms.
QUICK reduces the truncation error from O(Δx) to O(Δx³), improving accuracy
at higher Re without changing the grid.

The QUICK face-value interpolation (1D, flow in +x direction):

```
φ_e = (3/8)·φ_E + (6/8)·φ_P − (1/8)·φ_W    when u_e > 0
φ_e = (3/8)·φ_P + (6/8)·φ_E − (1/8)·φ_EE   when u_e < 0
```

### 6.2 Prompt used

```
Context:
- Language: MoonBit
- Grid: 41×41, uniform dx = 0.05
- Implementing QUICK scheme for the x-convective term of u-momentum: ∂(uu)/∂x
- phi[i][j] is the cell value
- Face velocity u_e = 0.5*(u[i][j] + u[i][j+1]) for east face between j and j+1
- Handle j=1 (no west-west neighbor) with fallback to 1st-order upwind

Generate a MoonBit function:
  fn quick_convection_x(phi: Array[Array[Double]], u: Array[Array[Double]])
    -> Array[Array[Double]]
that returns the QUICK convective flux divergence ∂(u·φ)/∂x at each interior node.
Only fill interior nodes (i in 1..ny-2, j in 1..nx-2).
```

### 6.3 AI-generated QUICK template

```moonbit
/// QUICK scheme: convective flux divergence ∂(u·φ)/∂x at interior nodes.
/// Uses 3rd-order upstream quadratic interpolation (Leonard 1979).
/// Falls back to 1st-order upwind at j=1 (no west-west neighbour).
fn quick_convection_x(
  phi : Array[Array[Double]],
  u   : Array[Array[Double]],
) -> Array[Array[Double]] {
  let flux_div = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // ── East face (between j and j+1) ──────────────────────────────────
      let ue = 0.5 * (u[i][j] + u[i][j + 1])
      let phi_e = if ue > 0.0 {
        if j >= 2 {
          (3.0 / 8.0) * phi[i][j + 1] + (6.0 / 8.0) * phi[i][j] - (1.0 / 8.0) * phi[i][j - 1]
        } else {
          // Fallback: 1st-order upwind (no west-west node)
          phi[i][j]
        }
      } else {
        if j + 2 < nx {
          (3.0 / 8.0) * phi[i][j] + (6.0 / 8.0) * phi[i][j + 1] - (1.0 / 8.0) * phi[i][j + 2]
        } else {
          phi[i][j + 1]
        }
      }

      // ── West face (between j-1 and j) ──────────────────────────────────
      let uw = 0.5 * (u[i][j - 1] + u[i][j])
      let phi_w = if uw > 0.0 {
        if j >= 3 {
          (3.0 / 8.0) * phi[i][j] + (6.0 / 8.0) * phi[i][j - 1] - (1.0 / 8.0) * phi[i][j - 2]
        } else {
          phi[i][j - 1]
        }
      } else {
        if j + 1 < nx {
          (3.0 / 8.0) * phi[i][j - 1] + (6.0 / 8.0) * phi[i][j] - (1.0 / 8.0) * phi[i][j + 1]
        } else {
          phi[i][j]
        }
      }

      // Flux divergence: (F_e - F_w) / dx
      flux_div[i][j] = (ue * phi_e - uw * phi_w) / dx
    }
  }
  flux_div
}
```

### 6.4 Verification test (AI-generated, then validated)

For a uniform flow u = 1.0 and φ = φ₀ (constant), the convective flux
divergence must be zero (constant advection, no change in φ).

```moonbit
test "quick_convection_uniform_zero_divergence" {
  // Constant phi with uniform rightward flow → zero divergence
  let phi = create_zeros_2d(ny, nx)
  let u   = create_zeros_2d(ny, nx)
  let phi0 = 2.5
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      phi[i][j] = phi0
      u[i][j]   = 1.0
    }
  }
  let div = quick_convection_x(phi, u)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(div[i][j].abs() < 1.0e-12)
    }
  }
}
```

### 6.5 Integration path

To use `quick_convection_x` in the u-momentum equation, replace the
current 1st-order upwind term in `cavity_flow_array`:

```moonbit
// Current (1st-order upwind):
let un = u[i][j]
let adv_u_x = if un > 0.0 { un * (un - u[i][j-1]) / dx }
              else         { un * (u[i][j+1] - un) / dx }

// QUICK replacement:
let quick_div = quick_convection_x(u, u)   // φ = u, advected by u
let adv_u_x = quick_div[i][j]
```

A similar `quick_convection_y` handles the ∂(vu)/∂y term.

---

## 7. Coding Conventions

### 7.1 MoonBit-specific rules

| Rule | Rationale |
|---|---|
| No `let mut` at global scope | Language limitation — use `Array[Int]` or `Array[Double]` wrapper |
| No leading `+` on continuation lines | Parsed as unary `+`, not line continuation |
| No `Double.sin()` | Use polynomial approximations in test RHS |
| Array allocation: always loop for 2D | `Array::make(n, Array::make(m, 0.0))` shares rows — allocate per row |

### 7.2 Naming conventions

```
g_*           Global state array  (e.g. g_u, g_u_pcg)
g_*_steps     Step counter array  (single element [0])
*_at(i,j)     Point accessor with OOB guard
*_n_steps(n)  Run n time steps
*_n_iter(n)   Run n SIMPLE iterations
```

### 7.3 Commit message format

```
type(scope): short description

Types: feat, fix, docs, test, refactor, chore, perf
Scopes: solver, viewer, wasm, pkg, readme, test
```
