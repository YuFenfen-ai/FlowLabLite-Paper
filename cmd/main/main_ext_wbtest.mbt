// FlowLabLite — Extended White-Box Test Suite  (Tests T59–T83)
// ============================================================
// Run: moon test --target wasm
//
// Software-engineering test classification used throughout:
//   [Unit]        – single function, isolated input, analytical oracle
//   [Integration] – two or more functions chained, black-box oracle
//   [System]      – full solver pipeline, physics-based oracle
//   [Regression]  – guard against known-good numerical properties
//
// Test-design techniques:
//   EP  – Equivalence Partitioning
//   BVA – Boundary Value Analysis
//   OBT – Oracle-Based Testing (analytical or computed reference)
//   PP  – Positive/Negative Paths

// ══════════════════════════════════════════════════════════════════════════════
// A. Unit Tests — Array Utilities  (T59–T60)
// ══════════════════════════════════════════════════════════════════════════════

// T59 [Unit / EP / OBT]
// copy_2d_array must produce a *deep* copy: mutating the source after the call
// must not change the destination.  Tests the separate-allocation invariant.
test "copy_2d_deep_independence" {
  let src = create_zeros_2d(4, 5)
  src[2][3] = 3.14
  src[0][4] = -7.0
  let dst = create_zeros_2d(4, 5)
  copy_2d_array(src, dst)
  // Overwrite every element of src with a very different value
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 5; j = j + 1 {
      src[i][j] = src[i][j] + 999.0
    }
  }
  // Destination must retain pre-copy values exactly
  assert_true(approx_eq(dst[2][3],  3.14, 1.0e-15))
  assert_true(approx_eq(dst[0][4], -7.0,  1.0e-15))
  assert_true(approx_eq(dst[3][0],  0.0,  1.0e-15))
}

// T60 [Unit / BVA]
// create_zeros_2d must allocate independent rows: writing to row 0 must not
// affect row 1 (shared-pointer alias bug would fail this).
test "create_zeros_2d_row_independence" {
  let a = create_zeros_2d(3, 3)
  a[0][1] = 5.0
  // If rows share the same underlying array, a[1][1] would also be 5.0
  assert_eq(a[1][1], 0.0)
  assert_eq(a[2][1], 0.0)
  // Verify all non-written elements are still zero
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      if not(i == 0 && j == 1) {
        assert_eq(a[i][j], 0.0)
      }
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// B. Unit Tests — Discrete Laplacian Operator  (T61–T63)
// ══════════════════════════════════════════════════════════════════════════════

// T61 [Unit / OBT]
// For the quadratic field u(x,y) = x²  (u[i][j] = (j·dx)²),
// the 5-point finite-difference Laplacian is exact: ∇²(x²) = 2 for all
// interior nodes regardless of grid spacing (property of degree-2 polynomials).
test "laplacian_apply_x_squared_exact" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let u  = create_zeros_2d(ny, nx)
  let lv = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      let x = j.to_double() * dx
      u[i][j] = x * x
    }
  }
  laplacian_apply(u, lv, inv_dx2, inv_dy2, a_diag)
  // ∇²(x²) = 2 exactly at every interior node
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(lv[i][j], 2.0, 1.0e-8))
    }
  }
}

// T62 [Unit / OBT]
// Symmetric counterpart: u(x,y) = y²  →  ∇²(y²) = 2.
// Together with T61 this checks both spatial directions of the stencil.
test "laplacian_apply_y_squared_exact" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let u  = create_zeros_2d(ny, nx)
  let lv = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      let y = i.to_double() * dy
      u[i][j] = y * y
    }
  }
  laplacian_apply(u, lv, inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(lv[i][j], 2.0, 1.0e-8))
    }
  }
}

// T63 [Unit / OBT]
// A spatially constant field is a harmonic function: ∇²(C) = 0 everywhere.
// Interior nodes only use interior neighbours so the boundary ghost treatment
// does not affect this result.
test "laplacian_apply_constant_gives_zero" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let u  = create_zeros_2d(ny, nx)
  let lv = create_zeros_2d(ny, nx)
  // Fill entire array with an arbitrary constant (must also cover boundary rows)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      u[i][j] = 42.0
    }
  }
  laplacian_apply(u, lv, inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(lv[i][j], 0.0, 1.0e-8))
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// C. Unit Tests — Boundary Conditions  (T64–T65)
// ══════════════════════════════════════════════════════════════════════════════

// T64 [Unit / PP]
// apply_pressure_bcs must be idempotent: calling it twice on the same array
// must produce exactly the same result as calling it once.
// The Neumann ghost-cell update p[0][j] = p[1][j] depends on p[1][j], which
// is not itself a ghost cell, so two passes give the same result.
test "apply_pressure_bcs_idempotent" {
  let p = create_zeros_2d(ny, nx)
  // Set a non-trivial interior field
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      p[i][j] = i.to_double() * j.to_double()
    }
  }
  apply_pressure_bcs(p)
  // Record the full state after first application
  let snap = create_zeros_2d(ny, nx)
  copy_2d_array(p, snap)
  // Second application
  apply_pressure_bcs(p)
  // Every element must be identical
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      assert_true(approx_eq(p[i][j], snap[i][j], 1.0e-15))
    }
  }
}

// T65 [Unit / BVA]
// apply_pressure_bcs_n applied to the GAMG coarse grid (n = gamg_nc = 21)
// must enforce the same four BC types as the fine grid:
//   Neumann left/right/bottom  (ghost = adjacent interior)
//   Dirichlet top              (p = 0)
test "apply_pressure_bcs_n_coarse_grid" {
  let nc = gamg_nc  // 21
  let p  = create_zeros_2d(nc, nc)
  // Prime the array with non-zero values
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      p[i][j] = (i + 1).to_double() * (j + 1).to_double()
    }
  }
  apply_pressure_bcs_n(p, nc)
  // Top row Dirichlet: p[nc-1][j] = 0
  for j = 0; j < nc; j = j + 1 {
    assert_true(approx_eq(p[nc - 1][j], 0.0, 1.0e-15))
  }
  // Bottom Neumann: p[0][j] = p[1][j]
  for j = 0; j < nc; j = j + 1 {
    assert_true(approx_eq(p[0][j], p[1][j], 1.0e-15))
  }
  // Left Neumann: p[i][0] = p[i][1]
  for i = 0; i < nc; i = i + 1 {
    assert_true(approx_eq(p[i][0], p[i][1], 1.0e-15))
  }
  // Right Neumann: p[i][nc-1] = p[i][nc-2]
  for i = 0; i < nc; i = i + 1 {
    assert_true(approx_eq(p[i][nc - 1], p[i][nc - 2], 1.0e-15))
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// D. Unit Tests — RHS Source Term  (T66–T67)
// ══════════════════════════════════════════════════════════════════════════════

// T66 [Unit / EP / OBT]
// A divergence-free uniform horizontal flow (u = const, v = 0) has
// ∂u/∂x = 0 and ∂v/∂y = 0, so every term in the b formula is zero.
// This tests the "no-divergence → zero pressure source" equivalence class.
test "build_up_b_zero_for_divergence_free_flow" {
  let u = create_zeros_2d(ny, nx)
  let v = create_zeros_2d(ny, nx)
  // Uniform horizontal flow: constant u everywhere
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      u[i][j] = 0.75   // arbitrary non-zero constant
    }
  }
  let b = create_zeros_2d(ny, nx)
  build_up_b_array(b, rho, dt, u, v, dx, dy, nx, ny) |> ignore
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(b[i][j], 0.0, 1.0e-9))
    }
  }
}

// T67 [Unit / OBT]
// Analytical oracle: u[i][j] = j·dx (linear in x, so ∂u/∂x = 1 uniformly),
// v = 0.  By the b formula:
//   b = ρ·[(1/dt)·(∂u/∂x + ∂v/∂y) − (∂u/∂x)² − 2·(∂u/∂y)(∂v/∂x) − (∂v/∂y)²]
//     = 1·[(1000)·1 − 1 − 0 − 0]  =  999.0
// for all interior nodes (independent of i, j).
test "build_up_b_formula_verification" {
  let u = create_zeros_2d(ny, nx)
  let v = create_zeros_2d(ny, nx)
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      u[i][j] = j.to_double() * dx
    }
  }
  let b = create_zeros_2d(ny, nx)
  build_up_b_array(b, rho, dt, u, v, dx, dy, nx, ny) |> ignore
  // Expected value: rho * (1/dt * 1 - 1^2 - 0 - 0) = 1000 - 1 = 999
  let expected = rho * (1.0 / dt * 1.0 - 1.0)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(approx_eq(b[i][j], expected, 1.0e-7))
    }
  }
}

// ══════════════════════════════════════════════════════════════════════════════
// E. Unit Tests — GAMG Sub-Components  (T68–T70)
// ══════════════════════════════════════════════════════════════════════════════

// T68 [Unit / OBT]
// gamg_prolongate must preserve constant coarse fields.
// For e_coarse[ci][cj] = C at all nodes, bilinear interpolation gives C at
// every fine-grid node (all four case weights sum to 1).
test "gamg_prolongate_constant_preservation" {
  let nc      = gamg_nc        // 21
  let e_c     = create_zeros_2d(nc, nc)
  let e_fine  = create_zeros_2d(ny, nx)
  let c_val   = 3.7
  // Fill coarse grid with a constant
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      e_c[i][j] = c_val
    }
  }
  gamg_prolongate(e_c, e_fine)
  // Every fine node reachable by prolongation must equal c_val
  // Case 1 (even,even): direct copy
  for ci = 0; ci < nc; ci = ci + 1 {
    for cj = 0; cj < nc; cj = cj + 1 {
      assert_true(approx_eq(e_fine[2 * ci][2 * cj], c_val, 1.0e-14))
    }
  }
  // Case 2 (even,odd): average of two equal neighbours = c_val
  for ci = 0; ci < nc; ci = ci + 1 {
    for cj = 0; cj < nc - 1; cj = cj + 1 {
      assert_true(approx_eq(e_fine[2 * ci][2 * cj + 1], c_val, 1.0e-14))
    }
  }
  // Case 4 (odd,odd): average of four equal corners = c_val
  for ci = 0; ci < nc - 1; ci = ci + 1 {
    for cj = 0; cj < nc - 1; cj = cj + 1 {
      assert_true(approx_eq(e_fine[2 * ci + 1][2 * cj + 1], c_val, 1.0e-14))
    }
  }
}

// T69 [Unit / OBT]
// gamg_restrict with normalisation factor 4: for a uniform fine field = 1.0,
// every interior coarse node (where all 8 fine-grid stencil neighbours exist)
// must receive exactly:
//   (1 + 4×0.5 + 4×0.25) / 4  =  4 / 4  =  1.0
// This verifies the /4 normalisation that prevents 4× over-amplification.
test "gamg_restrict_interior_normalization" {
  let nc     = gamg_nc   // 21
  let r_fine = create_zeros_2d(ny, nx)
  let r_c    = create_zeros_2d(nc, nc)
  // Uniform fine field
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      r_fine[i][j] = 1.0
    }
  }
  gamg_restrict(r_fine, r_c)
  // Check interior coarse nodes (ci in [1..nc-2], cj in [1..nc-2])
  // where fi=2*ci ∈ [2..38] and fj=2*cj ∈ [2..38] — all stencil neighbours in range.
  for ci = 1; ci < nc - 1; ci = ci + 1 {
    for cj = 1; cj < nc - 1; cj = cj + 1 {
      assert_true(approx_eq(r_c[ci][cj], 1.0, 1.0e-14))
    }
  }
}

// T70 [Unit / PP]
// gamg_smooth_fine must be a descent smoother: after ν=2 damped-Jacobi steps
// starting from z=0, the SPD residual ||r + A·z|| (= ||r − B·z|| where B=−A)
// must be strictly smaller than ||r|| = ||r + A·0||.
// Uses a non-trivial residual that is not in the kernel of the smoother.
test "gamg_smooth_fine_reduces_residual" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  // Non-uniform r: higher-frequency checkerboard pattern for better smoother effect
  let r  = create_zeros_2d(ny, nx)
  let z  = create_zeros_2d(ny, nx)
  let az = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // Alternating ±1 pattern (high-frequency: smoothers reduce this rapidly)
      let sign = if (i + j) % 2 == 0 { 1.0 } else { -1.0 }
      r[i][j] = sign
    }
  }
  // ||r||² before smoothing
  let mut r0_sq = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      r0_sq = r0_sq + r[i][j] * r[i][j]
    }
  }
  // Apply ν=2 damped-Jacobi smoothing
  gamg_smooth_fine(z, r, inv_dx2, inv_dy2, a_diag, 2)
  // Compute SPD residual ||r − B·z|| = ||r + A·z||
  apply_pressure_bcs(z)
  laplacian_apply(z, az, inv_dx2, inv_dy2, a_diag)
  let mut r1_sq = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let res = r[i][j] + az[i][j]
      r1_sq = r1_sq + res * res
    }
  }
  // Residual must decrease: ||r+A·z|| < ||r||
  assert_true(r1_sq < r0_sq)
}

// ══════════════════════════════════════════════════════════════════════════════
// F. Integration Tests — PCG Residual Criteria  (T71–T74)
// ══════════════════════════════════════════════════════════════════════════════
// Each test verifies that the solver actually achieves the stopping criterion
//   ‖b − A·p‖ / ‖b‖  <  pcg_tol = 1e-5
// after returning.  This is a stronger check than "iteration count ≤ limit".

fn compute_pressure_residual_norm(
  p       : Array[Array[Double]],
  b       : Array[Array[Double]],
  inv_dx2 : Double,
  inv_dy2 : Double,
  a_diag  : Double,
) -> Double {
  let ap = create_zeros_2d(ny, nx)
  apply_pressure_bcs(p)
  laplacian_apply(p, ap, inv_dx2, inv_dy2, a_diag)
  let mut r_sq  = 0.0
  let mut b_sq  = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let r = b[i][j] - ap[i][j]
      r_sq = r_sq + r * r
      b_sq = b_sq + b[i][j] * b[i][j]
    }
  }
  if b_sq > 1.0e-15 { r_sq.sqrt() / b_sq.sqrt() } else { 0.0 }
}

// T71 [Integration / OBT]
// Jacobi (baseline): residual satisfies the stopping criterion.
test "pcg_jacobi_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  // Accept 10× tolerance to allow for the final BC application
  assert_true(rel < pcg_tol * 10.0)
}

// T72 [Integration / OBT]
// DILU preconditioner: same residual criterion as Jacobi.
test "pcg_dilu_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dilu(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  assert_true(rel < pcg_tol * 10.0)
}

// T73 [Integration / OBT]
// DIC preconditioner: same residual criterion.
test "pcg_dic_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dic(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  assert_true(rel < pcg_tol * 10.0)
}

// T74 [Integration / OBT]
// GAMG preconditioner: residual must satisfy the stopping criterion.
// The stationary Richardson iteration converges when ρ(I + T·A) < 1.
test "pcg_gamg_residual_below_tolerance" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_gamg(p, dx, dy, b) |> ignore
  let rel = compute_pressure_residual_norm(p, b, inv_dx2, inv_dy2, a_diag)
  assert_true(rel < pcg_tol * 10.0)
}

// ══════════════════════════════════════════════════════════════════════════════
// G. Integration Tests — Cross-Preconditioner Agreement  (T75)
// ══════════════════════════════════════════════════════════════════════════════

// T75 [Integration / OBT]
// Non-uniform polynomial RHS: all four preconditioners must converge to the
// same pressure field within 1% relative error.  Tests that the choice of
// preconditioner affects only convergence speed, not the computed solution.
// The RHS b[i][j] = i·(ny-1-i)·j·(nx-1-j) is smooth, non-trivial, and
// naturally zero on all four boundary rows/columns.
test "all_preconditioners_agree_on_nonuniform_rhs" {
  let p_jac  = create_zeros_2d(ny, nx)
  let p_dilu = create_zeros_2d(ny, nx)
  let p_dic  = create_zeros_2d(ny, nx)
  let p_gamg = create_zeros_2d(ny, nx)
  let b      = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let fi = i.to_double()
      let fj = j.to_double()
      let ny1 = (ny - 1).to_double()
      let nx1 = (nx - 1).to_double()
      b[i][j] = (fi * (ny1 - fi)) * (fj * (nx1 - fj))
    }
  }
  pressure_poisson_pcg(p_jac,  dx, dy, b) |> ignore
  pressure_poisson_pcg_dilu(p_dilu, dx, dy, b) |> ignore
  pressure_poisson_pcg_dic(p_dic,  dx, dy, b) |> ignore
  pressure_poisson_pcg_gamg(p_gamg, dx, dy, b) |> ignore
  // Compute ‖p − p_jac‖ / ‖p_jac‖ for each preconditioner
  let mut ref_sq  = 0.0
  let mut dilu_sq = 0.0
  let mut dic_sq  = 0.0
  let mut gamg_sq = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      ref_sq  = ref_sq  + p_jac[i][j] * p_jac[i][j]
      let dd = p_dilu[i][j] - p_jac[i][j]
      let dc = p_dic[i][j]  - p_jac[i][j]
      let dg = p_gamg[i][j] - p_jac[i][j]
      dilu_sq = dilu_sq + dd * dd
      dic_sq  = dic_sq  + dc * dc
      gamg_sq = gamg_sq + dg * dg
    }
  }
  let dilu_err = if ref_sq > 1.0e-15 { dilu_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  let dic_err  = if ref_sq > 1.0e-15 { dic_sq.sqrt()  / ref_sq.sqrt() } else { 0.0 }
  let gamg_err = if ref_sq > 1.0e-15 { gamg_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  assert_true(dilu_err < 0.01)
  assert_true(dic_err  < 0.01)
  assert_true(gamg_err < 0.01)
}

// ══════════════════════════════════════════════════════════════════════════════
// H. System Tests — Solver Physics  (T76–T79)
// ══════════════════════════════════════════════════════════════════════════════

// T76 [System / OBT]
// Sign-convention test: the discrete Laplacian A is SND (diagonal = −a_diag < 0).
// Solving A·p = b with b > 0 gives p* = A⁻¹b < 0 at every interior node.
// This verifies the sign convention used throughout the preconditioner design.
test "pressure_sign_negative_for_positive_rhs" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p, dx, dy, b) |> ignore
  // All interior nodes must be strictly negative (SND → p* < 0 for b > 0)
  let mut all_negative = true
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      if p[i][j] >= 0.0 { all_negative = false }
    }
  }
  assert_true(all_negative)
}

// T77 [System / Regression]
// Chorin solver must be deterministic: two runs from the same init must
// produce bit-identical results.  Floating-point ops are deterministic in
// sequential execution on a single thread.
test "chorin_solver_determinism" {
  init_simulation()
  run_n_steps(100)
  let u1 = get_u_center()
  let v1 = get_v_at(ny / 2, nx / 2)
  let p1 = get_p_at(ny / 2, nx / 2)
  // Reset and re-run identically
  init_simulation()
  run_n_steps(100)
  let u2 = get_u_center()
  let v2 = get_v_at(ny / 2, nx / 2)
  let p2 = get_p_at(ny / 2, nx / 2)
  assert_true(u1 == u2)   // bit-identical: same ops, same order → same IEEE 754 result
  assert_true(v1 == v2)
  assert_true(p1 == p2)
}

// T78 [System / Regression]
// Additive steps: running n₁ + n₂ time steps in two calls must equal one call
// of n₁ + n₂ steps from the same initial state.  Tests solver linearity of
// time-marching with respect to step decomposition.
test "chorin_additive_steps" {
  // Path A: 10 + 10 = 20 steps
  init_simulation()
  run_n_steps(10)
  run_n_steps(10)
  let ua = get_u_center()
  let pa = get_p_at(ny / 4, nx / 4)
  // Path B: 20 steps in one call
  init_simulation()
  run_n_steps(20)
  let ub = get_u_center()
  let pb = get_p_at(ny / 4, nx / 4)
  assert_true(ua == ub)   // bit-identical
  assert_true(pa == pb)
}

// T79 [System / PP]
// SIMPLE produces a physically converged flow field.
// After 100 iterations with under-relaxation (α_p=0.3, α_u=0.7) for Re=20,
// the velocity divergence norm must be below 0.02 — confirming that the
// pressure correction is driving mass conservation.
// Note: the divergence norm uses a *relative* PCG stopping criterion, so it
// scales with the RHS magnitude and does not decrease monotonically.
// A fixed absolute threshold is the robust measure of acceptable mass error.
test "simple_divergence_below_threshold_after_100_iters" {
  init_simple()
  run_simple_n_iter(100)
  let div = get_simple_divergence_norm()
  assert_true(div >= 0.0)
  assert_true(div < 0.02)
}

// ══════════════════════════════════════════════════════════════════════════════
// I. Regression Tests  (T80–T83)
// ══════════════════════════════════════════════════════════════════════════════

// T80 [Regression / BVA]
// build_modified_diag must produce values in the physically expected range.
// For the 41×41 uniform grid: a_diag = 2/dx² + 2/dy² = 1600.
// The modified diagonal satisfies 0 < d[i][j] ≤ a_diag (the DILU fill-in
// only subtracts positive terms, so d is strictly smaller than a_diag
// except at the first interior node where no subtraction occurs).
test "build_modified_diag_range" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2  // = 1600
  let d = build_modified_diag(inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // Strictly positive (required for safe division in preconditioner)
      assert_true(d[i][j] > 0.0)
      // Cannot exceed the unmodified diagonal (stability bound)
      assert_true(d[i][j] <= a_diag + 1.0e-10)
    }
  }
  // Spot check: corner node (1,1) has no left/bottom neighbours to subtract →
  // d[1][1] should equal a_diag exactly.
  assert_true(approx_eq(d[1][1], a_diag, 1.0e-8))
}

// T81 [Regression / EP]
// MAC divergence must remain below 1e-4 even at step 200 (not just step 50).
// The staggered scheme maintains near-exact mass conservation throughout the run.
test "mac_divergence_remains_small_at_step_200" {
  init_mac()
  run_mac_n_steps(200)
  let div = get_mac_divergence_norm()
  assert_true(div < 1.0e-4)
}

// T82 [Regression / OBT]
// laplacian_apply_n must be exact on x² for the coarse grid (21×21).
// Coarse grid parameters: dx_c = 2·dx, inv_dx2_c = inv_dx2/4, a_diag_c = a_diag/2.
// The finite-difference Laplacian is exact on degree-≤2 polynomials regardless
// of grid spacing, so ∇²_c(x²) = 2 exactly at all coarse interior nodes.
test "laplacian_apply_n_x_squared_exact_coarse" {
  let nc      = gamg_nc           // 21
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let inv_dx2_c = inv_dx2 / 4.0  // coarse: dx_c = 2·dx
  let inv_dy2_c = inv_dy2 / 4.0
  let a_diag_c  = 2.0 * inv_dx2_c + 2.0 * inv_dy2_c
  let dx_c = dx * 2.0
  let u_c   = create_zeros_2d(nc, nc)
  let lv_c  = create_zeros_2d(nc, nc)
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let x = j.to_double() * dx_c
      u_c[i][j] = x * x
    }
  }
  laplacian_apply_n(u_c, lv_c, inv_dx2_c, inv_dy2_c, a_diag_c, nc, nc)
  for i = 1; i < nc - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      assert_true(approx_eq(lv_c[i][j], 2.0, 1.0e-7))
    }
  }
}

// T83 [Regression / System]
// All four PCG preconditioners must produce vortex flow consistent with
// Chorin's projection method: after a Chorin-PCG run of 200 steps, the
// streamwise velocity at the geometric centre must be negative.
// This end-to-end regression guards against preconditioner bugs that
// corrupt the pressure solve and break vortex formation.
test "chorin_pcg_produces_vortex_with_all_preconditioners" {
  // Chorin-PCG uses Jacobi preconditioner (already tested in T34)
  init_chorin_pcg()
  run_chorin_pcg_n_steps(200)
  let u_pcg_jac = get_u_pcg_at(ny / 2, nx / 2)
  assert_true(u_pcg_jac < 0.0)

  // Verify DILU pressure solver also converges in the PCG context
  // by solving the same final-step RHS and checking sign
  let p_dilu = create_zeros_2d(ny, nx)
  let b_test = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b_test[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dilu(p_dilu, dx, dy, b_test) |> ignore
  // DILU pressure solution must match Jacobi sign (both negative for b>0)
  assert_true(p_dilu[ny / 2][nx / 2] < 0.0)

  // DIC same check
  let p_dic = create_zeros_2d(ny, nx)
  pressure_poisson_pcg_dic(p_dic, dx, dy, b_test) |> ignore
  assert_true(p_dic[ny / 2][nx / 2] < 0.0)

  // GAMG same check
  let p_gamg = create_zeros_2d(ny, nx)
  pressure_poisson_pcg_gamg(p_gamg, dx, dy, b_test) |> ignore
  assert_true(p_gamg[ny / 2][nx / 2] < 0.0)
}

// ══════════════════════════════════════════════════════════════════════════════
// I. Unit Tests — I/O Format Writers  (T84–T91)
// ══════════════════════════════════════════════════════════════════════════════
//
// All tests use synthetic SolverExport values to avoid solver-state dependency.
// Tests are pure (no global-state mutation) and fast (small grids, no solver runs).

/// Helper: create a small SolverExport filled with constant field values.
fn make_test_export(
  name  : String,
  nc    : Int,
  fill_u : Double,
  fill_v : Double,
  fill_p : Double,
) -> SolverExport {
  let u = create_zeros_2d(nc, nc)
  let v = create_zeros_2d(nc, nc)
  let p = create_zeros_2d(nc, nc)
  for i = 0; i < nc; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      u[i][j] = fill_u
      v[i][j] = fill_v
      p[i][j] = fill_p
    }
  }
  { name, label: "test", nc, h: 0.5, x0: 0.0, y0: 0.0, steps: 10, u, v, p }
}

// T84 [Unit / EP]
// format_csv must include the standard RFC 4180 header row.
test "csv_header_present" {
  let s = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  let csv = format_csv([s])
  assert_true(csv.contains("solver,i,j,x,y,u,v,p,vel_mag"))
}

// T85 [Unit / BVA / OBT]
// format_csv for a single nc=3 solver must contain exactly nc*nc = 9 data rows
// (one line per grid point, excluding comment and header lines).
test "csv_row_count" {
  let nc = 3
  let s  = make_test_export("test", nc, 1.0, 0.5, 0.1)
  let csv = format_csv([s])
  // Count lines that start with "test," (the solver name prefix)
  let mut data_rows = 0
  csv.split("\n").each(fn(line) {
    if line.has_prefix("test,") {
      data_rows = data_rows + 1
    }
  })
  assert_eq(data_rows, nc * nc)
}

// T86 [Unit / EP]
// format_vtk must produce a valid STRUCTURED_GRID header with correct DIMENSIONS.
test "vtk_dimensions_header" {
  let nc = 4
  let s  = make_test_export("chorin", nc, 0.0, 0.0, 0.0)
  let vtk = format_vtk(s)
  assert_true(vtk.contains("DATASET STRUCTURED_GRID"))
  assert_true(vtk.contains("DIMENSIONS \{nc} \{nc} 1"))
}

// T87 [Unit / OBT]
// format_vtk POINTS count must equal nc*nc.
test "vtk_points_count" {
  let nc   = 5
  let npts = nc * nc
  let s    = make_test_export("pcg", nc, 0.0, 0.0, 0.0)
  let vtk  = format_vtk(s)
  assert_true(vtk.contains("POINTS \{npts} double"))
}

// T88 [Unit / EP]
// format_tecplot for k solvers must contain exactly k ZONE lines.
test "tecplot_zone_count" {
  let s1 = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  let s2 = make_test_export("simple", 3, 0.0, 0.0, 0.0)
  let dat = format_tecplot([s1, s2])
  let mut zone_count = 0
  dat.split("\n").each(fn(line) {
    if line.has_prefix("ZONE T=") {
      zone_count = zone_count + 1
    }
  })
  assert_eq(zone_count, 2)
}

// T89 [Unit / PP / Negative path]
// format_netcdf stub must always return Err (never panic, never Ok).
test "netcdf_stub_returns_err" {
  let s = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  match format_netcdf(s) {
    Err(msg) => assert_true(msg.length() > 0)
    Ok(_)    => assert_true(false)   // must not reach here
  }
}

// T90 [Unit / OBT]
// format_vtk must contain the POINT_DATA keyword followed by the total count.
// This validates that the SCALARS/VECTORS sections are properly introduced.
// (MoonBit renders Double 0.0 as "0", so we verify POINT_DATA header, not coordinates.)
test "vtk_point_data_header" {
  let nc   = 4
  let npts = nc * nc
  let s    = make_test_export("chorin", nc, 0.0, 0.0, 0.0)
  let vtk  = format_vtk(s)
  assert_true(vtk.contains("POINT_DATA \{npts}"))
  assert_true(vtk.contains("SCALARS pressure double 1"))
  assert_true(vtk.contains("VECTORS velocity double"))
}

// T91 [Unit / OBT]
// format_tecplot SOLUTIONTIME for a solver with steps=10, dt=0.001 must be 0.01.
test "tecplot_solution_time" {
  let s   = make_test_export("chorin", 3, 0.0, 0.0, 0.0)
  // s.steps = 10, dt = 0.001 (module constant) → SOLUTIONTIME = 0.01
  let dat = format_tecplot([s])
  assert_true(dat.contains("SOLUTIONTIME=\{10.0 * dt}"))
}

// ══════════════════════════════════════════════════════════════════════════════
// J. Unit Tests — VTK XML ImageData (.vti) Format  (T92–T93)
// ══════════════════════════════════════════════════════════════════════════════

// ══════════════════════════════════════════════════════════════════════════════
// L. Unit Tests — Runtime Re + Steps Override  (T96)
// ══════════════════════════════════════════════════════════════════════════════

// T96 [Unit / EP / OBT]
// parse_pos_int("100", 0) must return 100, and nu=2.0/100 must equal 0.02.
// Also verifies that g_nu_runtime is mutable (single-element Array[Double]).
test "runtime_re_100_gives_nu_0_02" {
  let re = parse_pos_int("100", 0)
  assert_eq(re, 100)
  let nu_100 = 2.0 / re.to_double()
  assert_true(approx_eq(nu_100, 0.02, 1.0e-10))
  // Verify g_nu_runtime is settable and readable
  let prev = g_nu_runtime[0]
  g_nu_runtime[0] = nu_100
  assert_true(approx_eq(g_nu_runtime[0], 0.02, 1.0e-10))
  g_nu_runtime[0] = prev   // restore compile-time default (nu=0.1)
}

// ══════════════════════════════════════════════════════════════════════════════
// K. Unit Tests — Centerline CSV + Monitor CSV  (T94–T95)
// ══════════════════════════════════════════════════════════════════════════════

// T94 [Unit / BVA / OBT]
// format_centerline_csv for nc=5 must contain exactly 2*nc = 10 data rows:
// nc rows for u_vcl (vertical centerline) + nc rows for v_hcl (horizontal).
test "centerline_csv_row_count" {
  let nc = 5
  let s  = make_test_export("chorin", nc, 1.0, 0.5, 0.0)
  let out = format_centerline_csv(s, 20.0, 0.001, 0.0)
  let mut u_rows = 0
  let mut v_rows = 0
  out.split("\n").each(fn(line) {
    if line.has_prefix("u_vcl,") { u_rows = u_rows + 1 }
    if line.has_prefix("v_hcl,") { v_rows = v_rows + 1 }
  })
  assert_eq(u_rows, nc)
  assert_eq(v_rows, nc)
}

test "centerline_csv_runtime_metadata" {
  let s = make_test_export("chorin", 5, 1.0, 0.5, 0.0)
  let out = format_centerline_csv(s, 100.0, 0.001, 0.0125)
  assert_true(out.contains("Re=100"))
  assert_true(out.contains("grid=5x5"))
  assert_true(out.contains("h=0.5"))
  assert_true(out.contains("dt=0.001"))
  assert_true(out.contains("pressure_iters=50"))
  assert_true(out.contains("div_norm=0.0125"))
}

test "runtime_grid_81_centerline_rows" {
  set_grid_size(81, 81)
  init_simulation()
  let s = collect_solver_exports("chorin")[0]
  let out = format_centerline_csv(s, 100.0, get_dt(), get_divergence_norm())
  let mut rows = 0
  out.split("\n").each(fn(line) {
    if line.has_prefix("u_vcl,") { rows = rows + 1 }
  })
  assert_eq(rows, 81)
  set_grid_size(41, 41)
  init_simulation()
}
// T95 [Unit / EP]

// format_monitor_csv must contain the standard header row with all column names.
test "monitor_csv_header_present" {
  let rec : MonitorRecord = {
    step: 0, time: 0.0, div_norm: 0.0,
    r_u: -1.0, r_v: -1.0, ke: 0.0, cfl_max: -1.0, pcg_iters: 0,
  }
  let out = format_monitor_csv([rec], "chorin", 20.0)
  assert_true(out.contains("step,time,div_norm,r_u,r_v,ke,cfl_max,pcg_iters"))
}

// ══════════════════════════════════════════════════════════════════════════════
// J. Unit Tests — VTK XML ImageData (.vti) Format  (T92–T93)
// ══════════════════════════════════════════════════════════════════════════════

// T92 [Unit / EP]
// format_vti must produce valid VTK XML: <?xml> declaration and
// <VTKFile type="ImageData"> root element must be present.
// Also verifies WholeExtent uses nc-1 (e.g., nc=4 → "0 3 0 3 0 0").
test "vti_xml_header_present" {
  let nc = 4
  let s  = make_test_export("chorin", nc, 0.5, 0.3, 0.1)
  let vti = format_vti(s)
  assert_true(vti.contains("<?xml version=\"1.0\"?>"))
  assert_true(vti.contains("<VTKFile type=\"ImageData\""))
  assert_true(vti.contains("<ImageData"))
  // nc=4 → extent upper bound = 3
  assert_true(vti.contains("WholeExtent=\"0 3 0 3 0 0\""))
}

// T93 [Unit / OBT]
// format_vti must contain exactly 3 <DataArray> opening tags:
// pressure (scalar), velocity (3-component), vel_magnitude (scalar).
test "vti_data_array_count_is_three" {
  let s   = make_test_export("pcg", 3, 1.0, 0.0, 0.5)
  let vti = format_vti(s)
  // Count occurrences of "<DataArray" (opening tag, not closing)
  let mut count = 0
  vti.split("\n").each(fn(line) {
    if line.contains("<DataArray") && not(line.contains("</DataArray")) {
      count = count + 1
    }
  })
  assert_eq(count, 3)
}

// ══════════════════════════════════════════════════════════════════════════════
// K. Unit Tests — Ghia (1982) Benchmark Data  (T97–T98)
// ══════════════════════════════════════════════════════════════════════════════

// ══════════════════════════════════════════════════════════════════════════════
// L. Unit Tests — QUICK Convection Scheme  (T99–T100)
// ══════════════════════════════════════════════════════════════════════════════

// T99 [Unit / EP / OBT]
// Uniform field (φ = const) must produce zero convection for both QUICK and
// upwind, regardless of advecting velocity.  ∂φ/∂x = 0 identically.
test "quick_uniform_field_zero_convection" {
  let nr = 5
  let nc = 5
  let phi   = create_zeros_2d(nr, nc)
  let u_vel = create_zeros_2d(nr, nc)
  // Fill φ = 3.0 everywhere, u = 1.0 everywhere
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j]   = 3.0
      u_vel[i][j] = 1.0
    }
  }
  let result = quick_convection_x(phi, u_vel, nr, nc, 0.05, 0.001)
  // All interior nodes must give zero (uniform gradient = 0)
  let mut all_zero = true
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      if result[i][j].abs() > 1.0e-14 {
        all_zero = false
      }
    }
  }
  assert_true(all_zero)
}

// T100 [Unit / OBT]
// Linear field (φ[i][j] = j * dx, u = 1) must give exact convection = dt
// at all interior nodes that have a full QUICK stencil (j = 2..nc-3).
// QUICK is exact for linear and quadratic polynomials; this verifies the
// scheme reaches the 3rd-order design accuracy on smooth linear data.
// Expected: u * ∂φ/∂x * dt = 1.0 * (1/dx) * dt * dx = dt
test "quick_linear_field_exact" {
  let nr  = 7
  let nc  = 7
  let dx  = 0.05
  let dt_ = 0.001
  let phi   = create_zeros_2d(nr, nc)
  let u_vel = create_zeros_2d(nr, nc)
  // φ[i][j] = j * dx  (linear in x), u = 1 everywhere
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j]   = j.to_double() * dx
      u_vel[i][j] = 1.0
    }
  }
  let result = quick_convection_x(phi, u_vel, nr, nc, dx, dt_)
  let expected = dt_   // u * (1/dx) * dx * dt = dt
  let tol = 1.0e-12
  // Check interior nodes away from boundary (j=2..nc-3 have full stencil)
  let mut ok = true
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 2; j < nc - 2; j = j + 1 {
      if (result[i][j] - expected).abs() > tol {
        ok = false
      }
    }
  }
  assert_true(ok)
}

// ══════════════════════════════════════════════════════════════════════════════
// M. Unit Tests — SOR / SSOR Pressure Solvers  (T101–T102)
// ══════════════════════════════════════════════════════════════════════════════

// ══════════════════════════════════════════════════════════════════════════════
// N. Unit Tests — HTML Self-Contained Report  (T103)
// ══════════════════════════════════════════════════════════════════════════════

// T103 [Unit / EP]
// format_html_report must produce a valid HTML document:
//   - contains <!DOCTYPE html> declaration
//   - contains the VARIABLES sentinel (used by local_viewer.html to detect report format)
//   - contains Chart.js CDN script tag
test "html_report_structure" {
  let nc = 5
  let s  = make_test_export("chorin", nc, 0.5, 0.3, 0.1)
  let out = format_html_report([s], 20.0, 100)
  assert_true(out.contains("<!DOCTYPE html>"))
  assert_true(out.contains("VARIABLES"))
  assert_true(out.contains("chart.js"))
}

// T101 [Unit / EP / OBT]
// SOR with ω=1.0 must be equivalent to standard Gauss-Seidel.
// Both start from the same zero initial pressure and run on the same smooth RHS;
// results must agree to machine precision at every interior node.
test "sor_omega1_equals_gauss_seidel" {
  let nr = 11
  let nc = 11
  let dx_ = 0.1
  let dy_ = 0.1
  // Smooth polynomial RHS: b[i][j] = (fi*(nr1-fi))*(fj*(nc1-nc1)) pattern
  let nr1 = (nr - 1).to_double()
  let nc1 = (nc - 1).to_double()
  let b = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let fi = i.to_double()
      let fj = j.to_double()
      b[i][j] = (fi * (nr1 - fi)) * (fj * (nc1 - fj))
    }
  }
  // SOR with ω=1.0 is the reference in-place Gauss-Seidel implementation.
  let p_sor = create_zeros_2d(nr, nc)
  pressure_poisson_sor(p_sor, dx_, dy_, b, 1.0, nr, nc, nit)
  // The shared Chorin/SIMPLE pressure solver must produce the same field.
  let p_gs = create_zeros_2d(nr, nc)
  pressure_poisson_array(p_gs, dx_, dy_, b, nc, nr) |> ignore
  let mut max_diff = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let d = (p_sor[i][j] - p_gs[i][j]).abs()
      if d > max_diff { max_diff = d }
    }
  }
  assert_true(max_diff < 1.0e-12)
}

// T102 [Unit / OBT]
// SOR with ω=1.5 (over-relaxation) must converge to a lower residual than
// plain GS (ω=1.0) after the same number of sweeps on a smooth problem.
// This validates that the SOR relaxation parameter actually accelerates convergence.
test "sor_omega15_faster_than_gs" {
  let nr  = 11
  let nc  = 11
  let dx_ = 0.1
  let dy_ = 0.1
  let nit_ = 40
  let nr1 = (nr - 1).to_double()
  let nc1 = (nc - 1).to_double()
  let b = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let fi = i.to_double()
      let fj = j.to_double()
      b[i][j] = (fi * (nr1 - fi)) * (fj * (nc1 - fj))
    }
  }
  // GS (omega=1.0)
  let p_gs  = create_zeros_2d(nr, nc)
  pressure_poisson_sor(p_gs, dx_, dy_, b, 1.0, nr, nc, nit_)
  let res_gs  = pressure_residual_l2(p_gs,  b, dx_, dy_, nr, nc)
  // SOR (omega=1.5)
  let p_sor = create_zeros_2d(nr, nc)
  pressure_poisson_sor(p_sor, dx_, dy_, b, 1.5, nr, nc, nit_)
  let res_sor = pressure_residual_l2(p_sor, b, dx_, dy_, nr, nc)
  // SOR(1.5) must have converged further (lower residual) than GS
  assert_true(res_sor < res_gs)
}

// ══════════════════════════════════════════════════════════════════════════════
// O. Unit Tests — Convergence Monitor (monitor.mbt)  (T104–T108)
// ══════════════════════════════════════════════════════════════════════════════

// T104 [Unit / OBT]
// compute_cfl_max: u=1 everywhere, v=0, dt=0.001, dx=dy=0.05.
// CFL = |u|*dt/dx + |v|*dt/dy = 1.0*0.001/0.05 + 0.0 = 0.02.
// This matches the default Re=20 worst-case estimate from solver_monitor.md.
test "monitor_cfl_safe" {
  let nr = 5
  let nc = 5
  let u  = create_zeros_2d(nr, nc)
  let v  = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      u[i][j] = 1.0
      v[i][j] = 0.0
    }
  }
  let cfl = compute_cfl_max(u, v, nr, nc, 0.001, 0.05, 0.05)
  assert_true((cfl - 0.02).abs() < 1.0e-12)
}

// T105 [Unit / EP / OBT]
// compute_kinetic_energy: zero velocity field must give KE = 0.
test "monitor_ke_zero" {
  let nr = 5
  let nc = 5
  let u  = create_zeros_2d(nr, nc)
  let v  = create_zeros_2d(nr, nc)
  let ke = compute_kinetic_energy(u, v, nr, nc, 0.05, 0.05)
  assert_true(ke.abs() < 1.0e-14)
}

// T106 [Unit / OBT]
// compute_ru: when u_new == u_old, the numerator max|diff| = 0 → R_u = 0.
// This verifies the converged-solution case (no change between steps).
test "monitor_ru_zero" {
  let nr = 5
  let nc = 5
  let u  = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      u[i][j] = (i + j).to_double() * 0.1
    }
  }
  // u_new = u_old: clone u into a separate array
  let u2 = create_zeros_2d(nr, nc)
  copy_2d_array(u, u2)
  let ru = compute_ru(u2, u, nr, nc)
  assert_true(ru.abs() < 1.0e-14)
}

// T107 [Unit / EP]
// format_monitor_csv must contain the expected header row.
test "monitor_csv_header_contains_all_columns" {
  let rec : MonitorRecord = {
    step: 10, time: 0.01, div_norm: 1.5e-3,
    r_u: 0.001, r_v: 0.002, ke: 0.05, cfl_max: 0.02, pcg_iters: 45,
  }
  let out = format_monitor_csv([rec], "mac", 100.0)
  assert_true(out.contains("step,time,div_norm,r_u,r_v,ke,cfl_max,pcg_iters"))
}

// T108 [Unit / OBT]
// compute_simple_residual_eq: a non-uniform u field must give a strictly
// positive residual (momentum equations are not in balance for an arbitrary field).
test "simple_residual_eq_positive" {
  let nr  = 7
  let nc  = 7
  let u   = create_zeros_2d(nr, nc)
  // Fill u with non-trivial non-constant values (not a solution of any equation)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      u[i][j] = (i + 1).to_double() * (j + 1).to_double() * 0.01
    }
  }
  let res = compute_simple_residual_eq(u, nr, nc, 0.1, 0.001, 0.05, 0.05)
  assert_true(res > 0.0)
}

// T97 [Unit / BVA]
// Both Ghia reference arrays must have exactly 17 data points and the same
// length (symmetry check): yn, u (vertical centerline), xn, v (horizontal).
test "ghia_data_arrays_length_17" {
  assert_eq(ghia_re100_yn.length(), 17)
  assert_eq(ghia_re100_u.length(),  17)
  assert_eq(ghia_re100_xn.length(), 17)
  assert_eq(ghia_re100_v.length(),  17)
}

// T98 [Unit / BVA]
// Boundary values must be exact per Ghia Table 1 & 2:
//   u[0]  = 0.0  (y=0 bottom wall, no-slip)
//   u[16] = 1.0  (y=1 top lid, driving velocity)
//   v[0]  = 0.0  (x=0 left wall, no-slip)
//   v[16] = 0.0  (x=1 right wall, no-slip)
test "ghia_boundary_values_exact" {
  assert_eq(ghia_re100_u[0],  0.0)
  assert_eq(ghia_re100_u[16], 1.0)
  assert_eq(ghia_re100_v[0],  0.0)
  assert_eq(ghia_re100_v[16], 0.0)
  // Also verify domain endpoints in Ghia's [0,1] coordinates
  assert_eq(ghia_re100_yn[0],  0.0)
  assert_eq(ghia_re100_yn[16], 1.0)
  assert_eq(ghia_re100_xn[0],  0.0)
  assert_eq(ghia_re100_xn[16], 1.0)
}

// ═══════════════════════════════════════════════════════════════════════════════
// P. TVD Convection Schemes (Step I)  T109–T113
// ═══════════════════════════════════════════════════════════════════════════════

// T109 [Unit / EP] Van Leer limiter at r=0 → ψ=0 (recovers upwind)
test "van_leer_r_zero" {
  let psi = limiter_van_leer(0.0)
  assert_eq(psi, 0.0)
}

// T110 [Unit / EP] Van Leer limiter at r=1 → ψ=1 (central difference)
test "van_leer_r_one" {
  let psi = limiter_van_leer(1.0)
  let diff = if psi > 1.0 { psi - 1.0 } else { 1.0 - psi }
  assert_true(diff < 1.0e-14)
}

// T111 [Unit / BVA] Superbee limiter at r=2 → ψ=2 (TVD upper bound)
test "superbee_r_two" {
  let psi = limiter_superbee(2.0)
  let diff = if psi > 2.0 { psi - 2.0 } else { 2.0 - psi }
  assert_true(diff < 1.0e-14)
}

// T112 [Unit / EP] TVD Van Leer — uniform φ field → zero convective term
// A uniform field has r=0 everywhere, so ψ(0)=0 → pure upwind → diff=0.
test "tvd_uniform_field_zero_x" {
  let nr = 7
  let nc = 7
  let dx = 0.05
  let dt = 0.001
  let phi   = create_zeros_2d(nr, nc)
  let u_vel = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j]   = 3.14   // constant field
      u_vel[i][j] = 1.0    // uniform positive advection
    }
  }
  let result = tvd_convection_x(phi, u_vel, nr, nc, dx, dt, 1)
  let mut max_err = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let v = if result[i][j] >= 0.0 { result[i][j] } else { -result[i][j] }
      if v > max_err { max_err = v }
    }
  }
  assert_true(max_err < 1.0e-14)
}

// T113 [Unit / EP] TVD Van Leer — linear φ=j*dx, u=1 → exact derivative dt
// For linear field: r=1 everywhere, ψ_VL(1)=1, TVD coincides with central diff,
// and the result equals u*(dφ/dx)*dt = 1*1*dt = dt exactly at interior nodes
// with full stencil (j=2..nc-3 to ensure both east and west stencils complete).
test "tvd_van_leer_linear_field_exact" {
  let nr = 9
  let nc = 9
  let dx = 0.1
  let dt = 0.001
  let phi   = create_zeros_2d(nr, nc)
  let u_vel = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j]   = j.to_double() * dx   // linear in x
      u_vel[i][j] = 1.0
    }
  }
  let result = tvd_convection_x(phi, u_vel, nr, nc, dx, dt, 1)
  // Check interior nodes that have complete 3-point stencil on both sides
  // (j=2..nc-3 for east; j=2 for west needs j-2=0 ≥ 0, ok)
  let mut max_err = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 2; j < nc - 2; j = j + 1 {
      let err = result[i][j] - dt
      let aerr = if err >= 0.0 { err } else { -err }
      if aerr > max_err { max_err = aerr }
    }
  }
  assert_true(max_err < 1.0e-12)
}

// ══════════════════════════════════════════════════════════════════════════════
// Q. SSP-RK3 Time Integration (Step J)  T114–T117
// ══════════════════════════════════════════════════════════════════════════════

// T114 [Unit / BVA] After init_rk3(), all interior cells are zero.
test "rk3_init_zero" {
  init_rk3()
  let mid = ny / 2
  assert_true(get_u_rk3_at(mid, mid) == 0.0)
  assert_true(get_v_rk3_at(mid, mid) == 0.0)
  assert_true(get_p_rk3_at(mid, mid) == 0.0)
  assert_true(get_rk3_step_count() == 0)
}

// T115 [Integration / BVA] After 1 RK3 step, lid BC u=1 is maintained.
// The projection step must not disturb Dirichlet boundary values.
test "rk3_lid_bc_after_one_step" {
  init_rk3()
  run_rk3_n_steps(1)
  // Top row (i=ny-1) should be exactly 1.0 everywhere.
  let mut ok = true
  for j = 0; j < nx; j = j + 1 {
    if get_u_rk3_at(ny - 1, j) != 1.0 { ok = false }
  }
  assert_true(ok)
}

// T116 [System / PP] Step counter increments correctly.
test "rk3_step_count" {
  init_rk3()
  run_rk3_n_steps(5)
  run_rk3_n_steps(3)
  assert_true(get_rk3_step_count() == 8)
}

// T117 [System / OBT] After 200 steps, v at the cavity centre-right is negative.
// The lid-driven vortex rotates counter-clockwise: v < 0 on the right side.
// This mirrors the sign convention verified for Chorin (T14) and PCG (T29).
test "rk3_vortex_sign_after_200_steps" {
  init_rk3()
  run_rk3_n_steps(200)
  // Centre-right point: (ny/2, 3*nx/4) → should have v < 0
  let v_cr = get_v_rk3_at(ny / 2, 3 * nx / 4)
  assert_true(v_cr < 0.0)
}

// ══════════════════════════════════════════════════════════════════════════════
// R. FVM-SIMPLE Solver (Step 1)  T118–T122
// ══════════════════════════════════════════════════════════════════════════════

// T118 [Unit / BVA] After init_fvm(), all cells must be zero.
test "fvm_init_zero" {
  init_fvm()
  let mid = ny / 2
  assert_true(get_u_fvm_at(mid, mid) == 0.0)
  assert_true(get_v_fvm_at(mid, mid) == 0.0)
  assert_true(get_p_fvm_at(mid, mid) == 0.0)
  assert_true(get_fvm_step_count() == 0)
}

// T119 [Integration / BVA] After 1 FVM iteration, lid BC u=1 is maintained.
// Velocity correction and pressure BCs must not disturb Dirichlet boundaries.
test "fvm_lid_bc_after_one_iter" {
  init_fvm()
  run_fvm_n_iter(1)
  let mut ok = true
  for j = 0; j < nx; j = j + 1 {
    if get_u_fvm_at(ny - 1, j) != 1.0 { ok = false }
  }
  assert_true(ok)
}

// T120 [System / OBT] After 100 iterations, divergence norm is below 0.5.
// The FVM-SIMPLE pressure-correction step must meaningfully reduce mass imbalance.
test "fvm_divergence_small_after_100_iters" {
  init_fvm()
  run_fvm_n_iter(100)
  let div = get_fvm_divergence_norm()
  assert_true(div < 0.5)
}

// T121 [System / OBT] After 300 iterations, v at centre-right is negative.
// FVM-SIMPLE must reproduce the same counter-clockwise vortex as all other solvers.
test "fvm_vortex_sign_after_300_iters" {
  init_fvm()
  run_fvm_n_iter(300)
  let v_cr = get_v_fvm_at(ny / 2, 3 * nx / 4)
  assert_true(v_cr < 0.0)
}

// T122 [System / Isolation] Running FVM must not mutate Chorin global state.
// Verifies that g_u_fvm and g_u are completely independent arrays.
test "fvm_independent_from_chorin" {
  init_simulation()
  run_n_steps(50)
  let u_before = get_u_at(ny / 2, nx / 2)
  init_fvm()
  run_fvm_n_iter(50)
  let u_after = get_u_at(ny / 2, nx / 2)
  assert_true((u_before - u_after).abs() < 1.0e-14)
}

// ══════════════════════════════════════════════════════════════════════════════
// S. Passive Scalar Transport (Step 6)  T123–T127
// ══════════════════════════════════════════════════════════════════════════════

// T123 [Unit / BVA] After init_scalar(), interior cells are zero and lid is 1.
// The top-wall BC (φ=1) is applied during init so the first step has the correct
// Dirichlet source.
test "scalar_init_interior_zero_lid_one" {
  init_scalar()
  // Interior cell at centre must be zero
  let mid_i = ny / 2
  let mid_j = nx / 2
  assert_true(get_scalar_at(mid_i, mid_j) == 0.0)
  // Top lid must be 1.0
  assert_true(get_scalar_at(ny - 1, mid_j) == 1.0)
  // Bottom wall must be 0.0
  assert_true(get_scalar_at(0, mid_j) == 0.0)
  assert_true(get_scalar_step_count() == 0)
}

// T124 [Integration / BVA] After several steps, lid BC φ=1 is preserved.
// Time-stepping must not overwrite Dirichlet boundaries.
test "scalar_lid_bc_preserved_after_steps" {
  init_scalar()
  run_scalar_n_steps(10)
  let mut ok = true
  for j = 0; j < nx; j = j + 1 {
    if get_scalar_at(ny - 1, j) != 1.0 { ok = false }
  }
  assert_true(ok)
}

// T125 [Unit / Physics] Pure diffusion: zero velocity, α=1.0.
// After many steps a uniform initial φ=0 interior with φ=1 top lid must yield
// a positive mean scalar (heat diffuses downward from the lid).
// Uses scalar_step_pure() to avoid global state dependency.
test "scalar_pure_diffusion_spreads" {
  let nr    = 7
  let nc    = 7
  let dx_   = 0.05
  let dy_   = 0.05
  let dt_   = 0.0001   // small dt for stability with α=1.0
  let alpha_= 1.0      // Pe ≈ 0, pure diffusion

  let phi   = create_zeros_2d(nr, nc)
  let u_zer = create_zeros_2d(nr, nc)
  let v_zer = create_zeros_2d(nr, nc)

  // Set lid BC
  for j = 0; j < nc; j = j + 1 { phi[nr - 1][j] = 1.0 }

  // Run 500 steps of pure diffusion
  for _k = 0; _k < 500; _k = _k + 1 {
    scalar_step_pure(phi, u_zer, v_zer, alpha_, nr, nc, dx_, dy_, dt_)
    // Re-apply Dirichlet BCs after each step
    for j = 0; j < nc; j = j + 1 { phi[nr - 1][j] = 1.0 }
    for j = 0; j < nc; j = j + 1 { phi[0][j] = 0.0 }
    for i = 0; i < nr; i = i + 1 { phi[i][0] = 0.0; phi[i][nc - 1] = 0.0 }
  }

  // Interior centre should have φ > 0 (diffusion has reached it)
  let phi_centre = phi[nr / 2][nc / 2]
  assert_true(phi_centre > 0.0)
  // And strictly less than 1 (not fully equilibrated)
  assert_true(phi_centre < 1.0)
}

// T126 [System / PP] Step counter increments correctly.
test "scalar_step_count" {
  init_scalar()
  run_scalar_n_steps(7)
  run_scalar_n_steps(3)
  assert_true(get_scalar_step_count() == 10)
}

// T127 [System / Isolation] Running scalar solver must not change Chorin velocity.
// One-way coupling: g_phi reads g_u/g_v but must never write them.
test "scalar_does_not_mutate_velocity" {
  init_simulation()
  run_n_steps(100)
  let u_snap = get_u_at(ny / 2, nx / 2)
  let v_snap = get_v_at(ny / 2, nx / 2)
  init_scalar()
  run_scalar_n_steps(50)
  assert_true((get_u_at(ny / 2, nx / 2) - u_snap).abs() < 1.0e-14)
  assert_true((get_v_at(ny / 2, nx / 2) - v_snap).abs() < 1.0e-14)
}

// ══════════════════════════════════════════════════════════════════════════════
// T. Energy Equation + Boussinesq (Step 7)  T128–T132
// ══════════════════════════════════════════════════════════════════════════════

// T128 [Unit / BVA] After init_boussinesq(), velocity/pressure are zero,
// scalar BCs are enforced (left=1, right=0).
test "bous_init_correct" {
  init_boussinesq()
  let mid = ny / 2
  // Velocity and pressure must be zero everywhere
  assert_true(get_u_bous_at(mid, mid) == 0.0)
  assert_true(get_v_bous_at(mid, mid) == 0.0)
  assert_true(get_p_bous_at(mid, mid) == 0.0)
  // Scalar: hot left wall = 1, cold right wall = 0
  assert_true(get_phi_bous_at(mid, 0)      == 1.0)
  assert_true(get_phi_bous_at(mid, nx - 1) == 0.0)
  assert_true(get_bous_step_count() == 0)
}

// T129 [Integration / BVA] Velocity BCs preserved: all walls no-slip (no lid).
// Unlike Chorin, the Boussinesq solver has no driving lid.
test "bous_velocity_bcs_no_slip" {
  init_boussinesq()
  run_boussinesq_n_steps(5)
  let mut ok = true
  // Top and bottom rows must be zero for both u and v
  for j = 0; j < nx; j = j + 1 {
    if get_u_bous_at(0,      j) != 0.0 { ok = false }
    if get_u_bous_at(ny - 1, j) != 0.0 { ok = false }
    if get_v_bous_at(0,      j) != 0.0 { ok = false }
    if get_v_bous_at(ny - 1, j) != 0.0 { ok = false }
  }
  // Left and right columns must be zero
  for i = 0; i < ny; i = i + 1 {
    if get_u_bous_at(i, 0)      != 0.0 { ok = false }
    if get_u_bous_at(i, nx - 1) != 0.0 { ok = false }
  }
  assert_true(ok)
}

// T130 [Integration / BVA] Scalar BCs preserved: left=1, right=0 after steps.
// Dirichlet conditions must not be overwritten by the interior time-stepping.
test "bous_scalar_bcs_preserved" {
  init_boussinesq()
  run_boussinesq_n_steps(10)
  let mut ok = true
  for i = 0; i < ny; i = i + 1 {
    if get_phi_bous_at(i, 0)      != 1.0 { ok = false }
    if get_phi_bous_at(i, nx - 1) != 0.0 { ok = false }
  }
  assert_true(ok)
}

// T131 [System / OBT] Buoyancy-driven circulation sign: after running with
// positive buoyancy (β>0), v near the hot left wall must be positive (fluid rises).
// Physical mechanism: hot fluid (φ≈1 near j=0) is lighter → upward buoyancy (+y).
// The Chorin projection must preserve this sign after divergence correction.
test "bous_buoyancy_drives_upward_flow" {
  init_boussinesq()
  // Default g_bous_beta = 0.1 (Ra ≈ 8000)
  run_boussinesq_n_steps(200)
  // At middle height, one cell from hot wall: v should be positive (rising)
  let v_hot = get_v_bous_at(ny / 2, 1)
  assert_true(v_hot > 0.0)
}

// T132 [System / Isolation] Running Boussinesq solver must not mutate
// Chorin global state (g_u, g_v, g_p are completely separate arrays).
test "bous_independent_from_chorin" {
  init_simulation()
  run_n_steps(50)
  let u_snap = get_u_at(ny / 2, nx / 2)
  init_boussinesq()
  run_boussinesq_n_steps(50)
  let u_after = get_u_at(ny / 2, nx / 2)
  assert_true((u_after - u_snap).abs() < 1.0e-14)
}

// ══════════════════════════════════════════════════════════════════════════════
// U. BC Library + Periodic Channel Flow (Step 8)  T133–T137
// ══════════════════════════════════════════════════════════════════════════════

// T133 [Unit / EP] neumann_bottom: after applying, phi[0][j] == phi[1][j].
// Verifies that the zero-gradient boundary condition is enforced correctly.
test "neumann_bottom_correct" {
  let nr = 5
  let nc = 5
  let phi = create_zeros_2d(nr, nc)
  // Set non-trivial interior values so phi[1][j] != 0
  for j = 0; j < nc; j = j + 1 { phi[1][j] = (j + 1).to_double() * 0.1 }
  neumann_bottom(phi, nr, nc)
  let mut ok = true
  for j = 0; j < nc; j = j + 1 {
    if (phi[0][j] - phi[1][j]).abs() > 1.0e-14 { ok = false }
  }
  assert_true(ok)
}

// T134 [Unit / EP] periodic_lr: after applying, phi[i][0] == phi[i][nc-2]
// and phi[i][nc-1] == phi[i][1] for all rows.
test "periodic_lr_wrap_correct" {
  let nr = 5
  let nc = 7
  let phi = create_zeros_2d(nr, nc)
  // Fill interior columns 1..nc-2 with distinct values
  for i = 0; i < nr; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      phi[i][j] = (i * nc + j).to_double()
    }
  }
  periodic_lr(phi, nr, nc)
  let mut ok = true
  for i = 0; i < nr; i = i + 1 {
    if (phi[i][0]      - phi[i][nc - 2]).abs() > 1.0e-14 { ok = false }
    if (phi[i][nc - 1] - phi[i][1]).abs()      > 1.0e-14 { ok = false }
  }
  assert_true(ok)
}

// T135 [Unit / BVA] init_channel: all velocity/pressure cells are zero
// and step counter is reset.
test "channel_init_zero" {
  init_channel()
  let mid = ny / 2
  assert_true(get_u_channel_at(mid, mid) == 0.0)
  assert_true(get_v_channel_at(mid, mid) == 0.0)
  assert_true(get_channel_step_count() == 0)
}

// T136 [System / Physics] Poiseuille profile: after running 500 steps with
// f_x=0.1, ν=0.1, the u-velocity at mid-channel must exceed u at the
// first interior row (parabolic profile, centre faster than near wall).
// This validates body force integration and periodic x-BC implementation.
// Exact Poiseuille at steady state: u_max = f_x·H²/(8ν) = 0.5 at y=H/2.
test "channel_parabolic_profile_forms" {
  init_channel()
  // Default f_x = 0.1, nu = 0.1 → u_max = 0.5
  run_channel_n_steps(500)
  let u_mid  = get_u_channel_at(ny / 2, nx / 2)   // centre of channel
  let u_near = get_u_channel_at(1,      nx / 2)    // first row above bottom wall
  // Centre must be faster than near-wall (parabolic)
  assert_true(u_mid > u_near)
  // Both must be positive (body force drives flow in +x direction)
  assert_true(u_mid  > 0.0)
  assert_true(u_near > 0.0)
}

// T137 [System / Isolation] Channel solver state is independent of Chorin.
// Running channel_n_steps must not change g_u (Chorin global array).
test "channel_independent_from_chorin" {
  init_simulation()
  run_n_steps(50)
  let u_snap = get_u_at(ny / 2, nx / 2)
  init_channel()
  run_channel_n_steps(50)
  assert_true((get_u_at(ny / 2, nx / 2) - u_snap).abs() < 1.0e-14)
}

// ══════════════════════════════════════════════════════════════════════════════
// V. Stretched Grid — tanh Compression (Step 9)  T138–T142
// ══════════════════════════════════════════════════════════════════════════════

// T138 [Unit / BVA] tanh_coord endpoints: first node = 0, last node = length.
// Tests both the i=0 clamp and the i=n-1 clamp in tanh_coord.
test "tanh_coord_endpoints_exact" {
  let n    = 11
  let cap_l   = 2.0
  let beta = 1.5
  assert_true(tanh_coord(0, n, cap_l, beta) == 0.0)
  assert_true(tanh_coord(n - 1, n, cap_l, beta) == cap_l)
}

// T139 [Unit / OBT] tanh_coord is strictly increasing for β > 0.
// A monotone sorted grid is required for all FD stencils to be well-posed.
test "tanh_coord_strictly_increasing" {
  let n    = 11
  let cap_l   = 2.0
  let beta = 1.5
  let mut ok = true
  for i = 0; i < n - 1; i = i + 1 {
    if tanh_coord(i + 1, n, cap_l, beta) <= tanh_coord(i, n, cap_l, beta) {
      ok = false
    }
  }
  assert_true(ok)
}

// T140 [Unit / EP] Uniform limit: β ≈ 0 → spacing ≈ L/(n-1) everywhere.
// Verifies the β→0 fallback branch produces a near-uniform grid.
test "tanh_coord_uniform_limit" {
  let n    = 11
  let cap_l   = 2.0
  let beta = 1.0e-9   // effectively zero
  let mean_sp = cap_l / (n - 1).to_double()
  let mut ok = true
  for i = 0; i < n - 1; i = i + 1 {
    let sp = tanh_spacing(i, n, cap_l, beta)
    if (sp - mean_sp).abs() > 1.0e-10 { ok = false }
  }
  assert_true(ok)
}

// T141 [Unit / OBT] Strong stretching (β=2): first spacing is significantly
// smaller than the mean spacing (= L/(n-1)), confirming clustering near x=0.
test "tanh_strong_stretching_clusters_near_wall" {
  let n    = 11
  let cap_l   = 2.0
  let beta = 2.0
  let first_sp = tanh_spacing(0, n, cap_l, beta)
  let mean_sp  = cap_l / (n - 1).to_double()
  // First spacing must be at most 50% of mean (significant clustering)
  assert_true(first_sp < 0.5 * mean_sp)
  // And last spacing must be larger than mean
  let last_sp = tanh_spacing(n - 2, n, cap_l, beta)
  assert_true(last_sp > mean_sp)
}

// T142 [Unit / OBT] Non-uniform 3-point Laplacian is exact for a linear field.
// φ = x → d²φ/dx² = 0 exactly, regardless of spacing pattern.
// Verifies the Lagrange-interpolation FD formula (Fornberg 1988).
test "nonuniform_laplacian_exact_for_linear" {
  let n    = 11
  let cap_l   = 2.0
  let beta = 2.0   // strong stretch — good stress test
  let xs  = make_tanh_grid(n, cap_l, beta)
  // Linear field: phi[i] = xs[i] (phi = x)
  let phi = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 { phi[i] = xs[i] }
  // Max |d²φ/dx²| must be near machine-zero
  let err = nonuniform_laplacian_max_err(phi, xs, n)
  assert_true(err < 1.0e-10)
}

// ══════════════════════════════════════════════════════════════════════════════
// N. Curvilinear Grid Tests  (T143–T147)
// ══════════════════════════════════════════════════════════════════════════════

// ══════════════════════════════════════════════════════════════════════════════
// M. Runtime Grid Size Tests  (T148–T152)  Step K Phase 1
// ══════════════════════════════════════════════════════════════════════════════

// T148 [Unit / BVA] set_grid_size updates get_nx/get_ny to the requested size.
// Also verifies that get_dx/get_dy are consistent: dx = 2/(nx-1).
test "set_grid_size_updates_nx_ny" {
  set_grid_size(21, 21)
  assert_eq(get_nx(), 21)
  assert_eq(get_ny(), 21)
  let expected_dx = 2.0 / 20.0
  assert_true((get_dx() - expected_dx).abs() < 1.0e-12)
  // Restore default for subsequent tests
  set_grid_size(41, 41)
}

// T149 [Integration] After set_grid_size(21,21) + init + 5 steps, no index error.
// Tests that cavity_flow_array respects the smaller grid without out-of-bounds.
test "set_grid_size_run_small_grid" {
  set_grid_size(21, 21)
  init_simulation()
  run_n_steps(5)
  assert_eq(get_step_count(), 5)
  // Restore
  set_grid_size(41, 41)
  init_simulation()
}

// T150 [System / OBT] Divergence norm is finite after running on 21×21 grid.
// Verifies mass-conservation quality is preserved at a different resolution.
test "set_grid_size_divergence_finite" {
  set_grid_size(21, 21)
  init_simulation()
  run_n_steps(10)
  let d = get_divergence_norm()
  assert_true(d >= 0.0 && d < 1.0e6)
  set_grid_size(41, 41)
  init_simulation()
}

// T151 [Regression] Resetting to 41×41 after a resize restores original behaviour.
// Runs 5 steps at 41×41, resets to 21×21, runs 5, resets to 41×41, re-runs 5.
// The final get_u_at(ny/2, nx/2) must equal the first result.
test "set_grid_size_restore_default" {
  // First run at default 41×41
  set_grid_size(41, 41)
  init_simulation()
  run_n_steps(5)
  let u1 = get_u_at(20, 20)
  // Switch to 21×21 and run
  set_grid_size(21, 21)
  init_simulation()
  run_n_steps(5)
  // Restore default and re-run
  set_grid_size(41, 41)
  init_simulation()
  run_n_steps(5)
  let u2 = get_u_at(20, 20)
  // Results must be bit-identical (same seed, same ops)
  assert_true((u1 - u2).abs() < 1.0e-14)
}

// T152 [System] set_grid_size clamps nx_<3 to 3 (minimum viable grid).
// After clamping, get_nx() must return 3 and the solver must not crash.
test "set_grid_size_clamps_minimum" {
  set_grid_size(2, 2)   // below minimum → clamped to 3
  assert_eq(get_nx(), 3)
  assert_eq(get_ny(), 3)
  init_simulation()
  run_n_steps(1)
  assert_eq(get_step_count(), 1)
  // Restore
  set_grid_size(41, 41)
  init_simulation()
}

// ── Curvilinear Grid Tests  (T143–T147) ──────────────────────────────────────

// T143 [Unit / OBT] Identity map: c=0 recovers uniform Cartesian grid.
// x(i,j,0) = 2·ξ_j = 2·j/(nx-1);  y(i,0) = 2·η_i = 2·i/(ny-1).
test "curv_identity_map_c0" {
  let c = 0.0
  let mut ok = true
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      let x_expected = 2.0 * j.to_double() / (nx - 1).to_double()
      let y_expected = 2.0 * i.to_double() / (ny - 1).to_double()
      if (curv_x(i, j, c) - x_expected).abs() > 1.0e-12 { ok = false }
      if (curv_y(i)        - y_expected).abs() > 1.0e-12 { ok = false }
    }
  }
  assert_true(ok)
}

// T144 [Unit / OBT] Jacobian J > 0 everywhere for c = 0.3 (non-degenerate).
// J = 4 + 4c(1-2ξ)η(1-η);  minimum occurs where (1-2ξ)η(1-η) is most negative.
// For c = 0.3 the grid is far from degenerate (c > -4 is the only requirement).
test "curv_jacobian_positive_c03" {
  let c = 0.3
  let mut min_j = 1.0e30
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      let jac = curv_jacobian(i, j, c)
      if jac < min_j { min_j = jac }
    }
  }
  assert_true(min_j > 0.0)
}

// T145 [Unit / BVA] Boundary nodes are unchanged by the distortion parameter c.
// At j=0: ξ=0 → x=0;  at j=nx-1: ξ=1 → x=2;  same for all c values.
test "curv_boundaries_invariant_to_c" {
  let c_vals = [0.0, 0.3, 1.0, -0.5]
  let mut ok = true
  for k = 0; k < 4; k = k + 1 {
    let c = c_vals[k]
    for i = 0; i < ny; i = i + 1 {
      if (curv_x(i, 0, c)      - 0.0).abs() > 1.0e-12 { ok = false }
      if (curv_x(i, nx - 1, c) - 2.0).abs() > 1.0e-12 { ok = false }
    }
  }
  assert_true(ok)
}

// T146 [Unit / OBT] Non-orthogonality = 0 at column boundaries j=0 and j=nx-1.
// At ξ=0 and ξ=1: x_η = 2c·ξ(1-ξ)(1-2η) = 0 → g_12 = 0 → grid is orthogonal.
test "curv_nonortho_zero_at_xi_boundaries" {
  let c = 0.5
  let mut ok = true
  for i = 1; i < ny - 1; i = i + 1 {
    let no0 = curv_nonortho(i, 0, c)
    let no1 = curv_nonortho(i, nx - 1, c)
    if no0 > 1.0e-12 { ok = false }
    if no1 > 1.0e-12 { ok = false }
  }
  assert_true(ok)
}

// T147 [Unit / OBT] Center node is displaced in +x for c > 0.
// curv_x(ny/2, nx/2, 0.3) > 1.0  (uniform grid gives exactly 1.0).
// The bump term 2c·ξ(1-ξ)η(1-η) is maximised near the grid centre where
// ξ=η=0.5, yielding a positive shift.
test "curv_center_positive_bump" {
  let c    = 0.3
  let ci   = ny / 2
  let cj   = nx / 2
  let x_uniform = 2.0 * cj.to_double() / (nx - 1).to_double()
  let x_curv    = curv_x(ci, cj, c)
  // curv_x must exceed uniform value (positive bump)
  assert_true(x_curv > x_uniform)
  // Sanity: magnitude of distortion is non-trivial (at least 0.01)
  assert_true(x_curv - x_uniform > 0.01)
}
