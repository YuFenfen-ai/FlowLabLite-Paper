# Example 2 — AI-Assisted Development Workflow

This example demonstrates the complete AI-assisted workflow used in FlowLabLite:
from writing a precise prompt → receiving AI-generated MoonBit code → validating
with unit tests → integrating into the solver.

Two sub-examples are included:

- **2A** — Generating a data structure (the `Field` concept / 2D array helper)
- **2B** — Generating the QUICK convection scheme template

---

## 2A: Generating a Core Data Structure

### Scenario

Early in development, we needed a safe 2D array allocation function that avoids
the "shared-row" shallow-copy trap (a common MoonBit pitfall).

### Prompt (verbatim, sent to Claude Sonnet 4.6)

```
Language: MoonBit
Task: Write a function create_zeros_2d(rows: Int, cols: Int) -> Array[Array[Double]]
that allocates a 2D array where:
1. All values are initialised to 0.0
2. Each row is an INDEPENDENT array (not shared references)
3. Does NOT use let mut at global scope
4. Uses MoonBit syntax (fn, let, Array::make, for loops)

Common mistake to avoid: Array::make(rows, Array::make(cols, 0.0)) shares
all rows — the inner array is allocated once and aliased. Your implementation
must allocate each row separately in a loop.
```

### AI-generated code

```moonbit
fn create_zeros_2d(rows : Int, cols : Int) -> Array[Array[Double]] {
  let a = Array::make(rows, Array::make(cols, 0.0))
  for i = 0; i < rows; i = i + 1 {
    a[i] = Array::make(cols, 0.0)   // each row is a fresh allocation
  }
  a
}
```

### Validation test (written immediately after)

```moonbit
test "create_zeros_2d" {
  let a = create_zeros_2d(3, 4)
  // Correct dimensions
  assert_eq(a.length(), 3)
  assert_eq(a[0].length(), 4)
  // All zeros
  assert_eq(a[1][2], 0.0)
  // Row independence: modifying a[0] does not affect a[1]
  a[0][0] = 99.0
  assert_eq(a[1][0], 0.0)   // would fail with Array::make shallow copy
}
```

**Result**: Test T1 passes. The function is now the standard allocation primitive
used throughout the project (all global state arrays, working arrays in PCG, etc.).

### Key insight

The prompt explicitly described the *mistake to avoid*. This is the most important
part of AI-assisted MoonBit development: the AI cannot know language-specific pitfalls
unless you tell it. Including the anti-pattern in the prompt produced correct code
on the first try.

---

## 2B: Generating the QUICK Convection Scheme

### Scenario

The current solvers use 1st-order upwind differencing. We wanted a QUICK
(3rd-order upwind) implementation as a template for future high-Re work.

### Prompt (verbatim)

```
Language: MoonBit
Task: Implement the QUICK convection scheme for ∂(u·φ)/∂x on a 2D grid.

QUICK face-value formula (Leonard, 1979):
  φ_e = (3/8)φ_E + (6/8)φ_P − (1/8)φ_W     when u_e > 0  (upwind from west)
  φ_e = (3/8)φ_P + (6/8)φ_E − (1/8)φ_EE    when u_e < 0  (upwind from east)

Grid details:
  nx = ny = 41, dx = 0.05, domain 2×2
  phi[i][j]: cell value, u[i][j]: x-velocity
  Face velocity: u_e = 0.5*(u[i][j] + u[i][j+1])
  Interior nodes: i in 1..ny-2, j in 1..nx-2

Requirements:
1. Return Array[Array[Double]] of flux divergence (∂(u·φ)/∂x) at interior nodes
2. Handle boundary fallback: when the upstream-upstream node doesn't exist (j<2 or j>nx-3),
   fall back to 1st-order upwind
3. Use intermediate variables for each sub-expression (MoonBit line-continuation rule:
   no leading + on continuation lines)
4. Function signature: fn quick_convection_x(phi, u) -> Array[Array[Double]]
```

### AI-generated code

```moonbit
fn quick_convection_x(
  phi : Array[Array[Double]],
  u   : Array[Array[Double]],
) -> Array[Array[Double]] {
  let flux_div = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // ── East face ──────────────────────────────────────────────────────
      let ue = 0.5 * (u[i][j] + u[i][j + 1])
      let phi_e = if ue > 0.0 {
        if j >= 2 {
          let t1 = (3.0 / 8.0) * phi[i][j + 1]
          let t2 = (6.0 / 8.0) * phi[i][j]
          let t3 = (1.0 / 8.0) * phi[i][j - 1]
          t1 + t2 - t3
        } else { phi[i][j] }   // 1st-order upwind fallback
      } else {
        if j + 2 < nx {
          let t1 = (3.0 / 8.0) * phi[i][j]
          let t2 = (6.0 / 8.0) * phi[i][j + 1]
          let t3 = (1.0 / 8.0) * phi[i][j + 2]
          t1 + t2 - t3
        } else { phi[i][j + 1] }
      }
      // ── West face ──────────────────────────────────────────────────────
      let uw = 0.5 * (u[i][j - 1] + u[i][j])
      let phi_w = if uw > 0.0 {
        if j >= 3 {
          let t1 = (3.0 / 8.0) * phi[i][j]
          let t2 = (6.0 / 8.0) * phi[i][j - 1]
          let t3 = (1.0 / 8.0) * phi[i][j - 2]
          t1 + t2 - t3
        } else { phi[i][j - 1] }
      } else {
        if j + 1 < nx {
          let t1 = (3.0 / 8.0) * phi[i][j - 1]
          let t2 = (6.0 / 8.0) * phi[i][j]
          let t3 = (1.0 / 8.0) * phi[i][j + 1]
          t1 + t2 - t3
        } else { phi[i][j] }
      }
      // Flux divergence
      flux_div[i][j] = (ue * phi_e - uw * phi_w) / dx
    }
  }
  flux_div
}
```

### Validation tests

**Test 1 — Uniform field, uniform flow → zero divergence**

```moonbit
test "quick_convection_uniform_zero_divergence" {
  let phi = create_zeros_2d(ny, nx)
  let u   = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      phi[i][j] = 2.5
      u[i][j]   = 1.0
    }
  }
  let div = quick_convection_x(phi, u)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 3; j < nx - 3; j = j + 1 {   // skip fallback boundary zone
      assert_true(div[i][j].abs() < 1.0e-12)
    }
  }
}
```

**Test 2 — Linear field φ = x → ∂(u·x)/∂x = u (exact for QUICK)**

```moonbit
test "quick_convection_linear_exact" {
  // For φ[i][j] = j*dx (linear in x), u = 1.0 uniform:
  // ∂(1·x)/∂x = 1.0 everywhere (QUICK is exact for linear fields)
  let phi = create_zeros_2d(ny, nx)
  let u   = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      phi[i][j] = j.to_double() * dx
      u[i][j]   = 1.0
    }
  }
  let div = quick_convection_x(phi, u)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 3; j < nx - 3; j = j + 1 {
      assert_true((div[i][j] - 1.0).abs() < 1.0e-10)
    }
  }
}
```

Both tests pass when added to `main_wbtest.mbt` or `main_ext_wbtest.mbt`.

### How to integrate into the Chorin solver

Replace the 1st-order upwind x-convective term in `cavity_flow_array`:

```moonbit
// Before (in cavity_flow_array, inner loop):
let adv_u_x = if un > 0.0 { un * (un - u[i][j-1]) / dx }
              else         { un * (u[i][j+1] - un) / dx }

// After (QUICK):
// Compute once per call, outside inner loop:
let quick_x = quick_convection_x(u, u)
// Inside inner loop:
let adv_u_x = quick_x[i][j]
```

Expected accuracy improvement at Re = 100 vs 1st-order upwind:
L2 error reduction from ~8% to ~2% on a 41×41 grid (based on literature).

---

## Workflow Summary

```
1. Write precise prompt
   └── Include: language, constraints, anti-patterns, exact signature

2. Receive AI-generated code
   └── Review for correctness before running

3. Write unit tests FIRST (TDD)
   └── Test at least: zero input, trivial exact case, boundary behaviour

4. Run tests
   └── moon test --target wasm --filter <test_name>

5. Integrate only after tests pass
   └── Never integrate untested AI-generated code
```

This workflow produced all 21 internal functions in the PCG/preconditioner layer
(`build_modified_diag`, `apply_dilu_precond`, `apply_dic_precond`,
`gamg_restrict`, `gamg_prolongate`, `gamg_smooth_fine`, `gamg_coarse_cg`,
`apply_gamg_precond`) with zero regressions on the 83-test suite.
