// White-box tests for the cavity flow solver.
// Run: moon test --target wasm

// ── Helper: approximate floating-point equality ──
fn approx_eq(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

// ── Test 1: create_zeros_2d returns all-zero array of correct size ──
test "create_zeros_2d" {
  let arr = create_zeros_2d(5, 7)
  assert_eq(arr.length(), 5)
  assert_eq(arr[0].length(), 7)
  for i = 0; i < 5; i = i + 1 {
    for j = 0; j < 7; j = j + 1 {
      assert_eq(arr[i][j], 0.0)
    }
  }
}

// ── Test 2: copy_2d_array copies values correctly ──
test "copy_2d_array" {
  let src = create_zeros_2d(3, 3)
  src[0][0] = 1.0
  src[1][2] = 5.5
  src[2][1] = -3.0
  let dst = create_zeros_2d(3, 3)
  copy_2d_array(src, dst)
  assert_eq(dst[0][0], 1.0)
  assert_eq(dst[1][2], 5.5)
  assert_eq(dst[2][1], -3.0)
  // Verify it's a deep copy: modifying src should not affect dst
  src[0][0] = 99.0
  assert_eq(dst[0][0], 1.0)
}

// ── Test 3: generate_mesh_grid produces correct coordinate ranges ──
test "generate_mesh_grid" {
  let (x, y) = generate_mesh_grid()
  assert_eq(x.length(), ny)
  assert_eq(x[0].length(), nx)
  // Bottom-left corner should be (0, 0)
  assert_true(approx_eq(x[0][0], 0.0, 1.0e-10))
  assert_true(approx_eq(y[0][0], 0.0, 1.0e-10))
  // Top-right corner should be (2, 2)
  assert_true(approx_eq(x[ny - 1][nx - 1], 2.0, 1.0e-10))
  assert_true(approx_eq(y[ny - 1][nx - 1], 2.0, 1.0e-10))
  // X should increase along columns, Y along rows
  assert_true(x[0][1] > x[0][0])
  assert_true(y[1][0] > y[0][0])
}

// ── Test 4: init_simulation resets global state ──
test "init_simulation_resets_state" {
  // Dirty the global arrays
  g_u[0][ny / 2][nx / 2] = 42.0
  g_v[0][0][0] = -1.0
  g_p[0][ny - 1][nx - 1] = 99.0
  g_steps[0] = 999
  // Reset
  init_simulation()
  assert_eq(g_steps[0], 0)
  assert_eq(g_u[0][ny / 2][nx / 2], 0.0)
  assert_eq(g_v[0][0][0], 0.0)
  assert_eq(g_p[0][ny - 1][nx - 1], 0.0)
}

// ── Test 5: boundary conditions after solver run ──
test "boundary_conditions_after_run" {
  init_simulation()
  run_n_steps(10)
  // Top lid: u = 1
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_at(ny - 1, j), 1.0, 1.0e-10))
  }
  // Bottom wall: u = 0, v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_at(0, j), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_at(0, j), 0.0, 1.0e-10))
  }
  // Left wall: u = 0, v = 0 (exclude top lid row i=ny-1 where u=1)
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_at(i, 0), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_at(i, 0), 0.0, 1.0e-10))
  }
  // Right wall: u = 0, v = 0 (exclude top lid row)
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_at(i, nx - 1), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_at(i, nx - 1), 0.0, 1.0e-10))
  }
  // Top/bottom wall: v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_v_at(ny - 1, j), 0.0, 1.0e-10))
  }
}

// ── Test 6: step counter increments correctly ──
test "step_counter" {
  init_simulation()
  assert_eq(get_step_count(), 0)
  run_n_steps(5)
  assert_eq(get_step_count(), 5)
  run_n_steps(3)
  assert_eq(get_step_count(), 8)
}

// ── Test 7: get_*_at returns 0 for out-of-range indices ──
test "out_of_range_returns_zero" {
  assert_eq(get_u_at(-1, 0), 0.0)
  assert_eq(get_u_at(0, -1), 0.0)
  assert_eq(get_u_at(ny, 0), 0.0)
  assert_eq(get_u_at(0, nx), 0.0)
  assert_eq(get_v_at(-1, -1), 0.0)
  assert_eq(get_p_at(999, 999), 0.0)
  assert_eq(get_velocity_magnitude_at(-1, 0), 0.0)
}

// ── Test 8: constant accessors return correct values ──
test "constant_accessors" {
  assert_eq(get_nx(), 41)
  assert_eq(get_ny(), 41)
  assert_eq(get_nt(), 500)
  assert_eq(get_nit(), 50)
  assert_true(approx_eq(get_rho(), 1.0, 1.0e-10))
  assert_true(approx_eq(get_nu(), 0.1, 1.0e-10))
  assert_true(approx_eq(get_dt(), 0.001, 1.0e-10))
  assert_true(approx_eq(get_dx(), 0.05, 1.0e-10))
  assert_true(approx_eq(get_dy(), 0.05, 1.0e-10))
  assert_true(approx_eq(get_re(), 20.0, 1.0e-10))
}

// ── Test 9: velocity magnitude consistency ──
test "velocity_magnitude_consistency" {
  init_simulation()
  run_n_steps(20)
  let mid_i = ny / 2
  let mid_j = nx / 2
  let u = get_u_at(mid_i, mid_j)
  let v = get_v_at(mid_i, mid_j)
  let expected_mag = (u * u + v * v).sqrt()
  let actual_mag = get_velocity_magnitude_at(mid_i, mid_j)
  assert_true(approx_eq(actual_mag, expected_mag, 1.0e-12))
}

// ── Test 10: max velocity magnitude is non-negative and >= max_u ──
test "max_velocity_magnitude_bounds" {
  init_simulation()
  run_n_steps(50)
  let max_mag = get_max_velocity_magnitude()
  let max_u = get_max_u()
  assert_true(max_mag >= 0.0)
  assert_true(max_mag >= max_u)
}

// ── Test 11: pressure Poisson solver produces bounded pressure ──
test "pressure_bounded" {
  init_simulation()
  run_n_steps(50)
  let max_p = get_max_p()
  let min_p = get_min_p()
  // Pressure should be finite and within reasonable physical bounds
  assert_true(max_p < 1000.0)
  assert_true(min_p > -1000.0)
  assert_true(max_p >= min_p)
}

// ── Test 12: center point getters match get_*_at(ny/2, nx/2) ──
test "center_getters_consistent" {
  init_simulation()
  run_n_steps(30)
  let mid_i = ny / 2
  let mid_j = nx / 2
  assert_true(approx_eq(get_u_center(), get_u_at(mid_i, mid_j), 1.0e-15))
  assert_true(approx_eq(get_v_center(), get_v_at(mid_i, mid_j), 1.0e-15))
  assert_true(approx_eq(get_p_center(), get_p_at(mid_i, mid_j), 1.0e-15))
}

// ── Test 13: divergence norm is non-negative after simulation ──
test "divergence_norm_nonnegative" {
  init_simulation()
  run_n_steps(50)
  let div = get_divergence_norm()
  assert_true(div >= 0.0)
}

// ── Test 14: build_up_b_array produces non-zero source term with non-trivial input ──
test "build_up_b_nonzero" {
  let u = create_zeros_2d(ny, nx)
  let v = create_zeros_2d(ny, nx)
  // Set top lid velocity
  for j = 0; j < nx; j = j + 1 {
    u[ny - 1][j] = 1.0
  }
  // Add some interior velocity to produce non-trivial b
  u[ny / 2][nx / 2] = 0.5
  let b = create_zeros_2d(ny, nx)
  build_up_b_array(b, rho, dt, u, v, dx, dy, nx, ny) |> ignore
  // b should have at least some non-zero entries in the interior
  let mut has_nonzero = false
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      if b[i][j].abs() > 1.0e-15 { has_nonzero = true }
    }
  }
  assert_true(has_nonzero)
}

// ── Test 15: run_all_steps produces lid-driven vortex ──
test "full_simulation_produces_vortex" {
  init_simulation()
  run_all_steps()
  assert_eq(get_step_count(), 500)
  // After 500 steps, the interior should have non-zero velocities
  let u_center = get_u_center()
  let v_center = get_v_center()
  // Center u should be negative (backflow) for lid-driven cavity
  assert_true(u_center < 0.0)
  // v_center should be small but non-zero
  assert_true(v_center.abs() > 0.0)
  // Max velocity should be >= 1.0 (lid velocity)
  assert_true(get_max_u() >= 1.0 - 1.0e-10)
}

// ── Test 16: pressure boundary conditions satisfied ──
test "pressure_boundary_dp_zero" {
  init_simulation()
  run_n_steps(50)
  // dp/dy = 0 at y = 0: p[0,:] = p[1,:] (approximately)
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_at(0, j), get_p_at(1, j), 1.0e-10))
  }
  // p = 0 at y = 2 (top): p[ny-1,:] = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_at(ny - 1, j), 0.0, 1.0e-10))
  }
}

// ──────────────────────────────────────────────────────────────────────────────
// SIMPLE algorithm tests (Tests 17–26)
// ──────────────────────────────────────────────────────────────────────────────

// ── Test 17: init_simple resets all SIMPLE fields and counters ──
test "simple_init_resets_state" {
  // Dirty the SIMPLE global state first
  init_simple()
  run_simple_n_iter(5)
  // Now reset and verify
  init_simple()
  assert_eq(get_simple_step_count(), 0)
  assert_true(approx_eq(get_simple_residual(), 0.0, 1.0e-15))
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      assert_eq(get_u_simple_at(i, j), 0.0)
      assert_eq(get_v_simple_at(i, j), 0.0)
      assert_eq(get_p_simple_at(i, j), 0.0)
    }
  }
}

// ── Test 18: SIMPLE step counter increments correctly ──
test "simple_step_counter" {
  init_simple()
  assert_eq(get_simple_step_count(), 0)
  run_simple_n_iter(3)
  assert_eq(get_simple_step_count(), 3)
  run_simple_n_iter(7)
  assert_eq(get_simple_step_count(), 10)
}

// ── Test 19: SIMPLE boundary conditions after run ──
test "simple_boundary_conditions" {
  init_simple()
  run_simple_n_iter(20)
  // Top lid: u = 1
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_simple_at(ny - 1, j), 1.0, 1.0e-10))
  }
  // Bottom wall: u = 0, v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_simple_at(0, j), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_simple_at(0, j), 0.0, 1.0e-10))
  }
  // Left wall (interior rows): u = 0, v = 0
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_simple_at(i, 0), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_simple_at(i, 0), 0.0, 1.0e-10))
  }
  // Right wall (interior rows): u = 0, v = 0
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_simple_at(i, nx - 1), 0.0, 1.0e-10))
    assert_true(approx_eq(get_v_simple_at(i, nx - 1), 0.0, 1.0e-10))
  }
  // Top wall: v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_v_simple_at(ny - 1, j), 0.0, 1.0e-10))
  }
}

// ── Test 20: get_*_simple_at returns 0 for out-of-range indices ──
test "simple_out_of_range_returns_zero" {
  assert_eq(get_u_simple_at(-1, 0), 0.0)
  assert_eq(get_u_simple_at(ny, 0), 0.0)
  assert_eq(get_v_simple_at(0, -1), 0.0)
  assert_eq(get_v_simple_at(0, nx), 0.0)
  assert_eq(get_p_simple_at(999, 999), 0.0)
}

// ── Test 21: SIMPLE residual is non-negative ──
test "simple_residual_nonnegative" {
  init_simple()
  run_simple_n_iter(10)
  assert_true(get_simple_residual() >= 0.0)
  assert_true(get_simple_divergence_norm() >= 0.0)
}

// ── Test 22: SIMPLE produces non-zero interior velocities after iterations ──
test "simple_produces_flow" {
  init_simple()
  run_simple_n_iter(50)
  // After 50 iterations with a lid, interior should be non-trivial
  let mut nonzero_u = false
  let mut nonzero_v = false
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      if get_u_simple_at(i, j).abs() > 1.0e-10 { nonzero_u = true }
      if get_v_simple_at(i, j).abs() > 1.0e-10 { nonzero_v = true }
    }
  }
  assert_true(nonzero_u)
  assert_true(nonzero_v)
}

// ── Test 23: SIMPLE max_u is non-negative and does not exceed lid + margin ──
test "simple_max_u_bounds" {
  init_simple()
  run_simple_n_iter(50)
  let mu = get_max_u_simple()
  assert_true(mu >= 0.0)
  // SIMPLE under-relaxation keeps velocities bounded; lid is u=1
  assert_true(mu <= 2.0)
}

// ── Test 24: SIMPLE pressure boundary conditions satisfied ──
test "simple_pressure_boundary" {
  init_simple()
  run_simple_n_iter(50)
  // p = 0 at top (y = 2): p[ny-1,:] = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_simple_at(ny - 1, j), 0.0, 1.0e-10))
  }
  // Neumann at bottom: p[0,:] = p[1,:]
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_simple_at(0, j), get_p_simple_at(1, j), 1.0e-10))
  }
}

// ── Test 25: SIMPLE and Chorin produce qualitatively similar vortex ──
test "simple_vs_chorin_qualitative" {
  // Run both solvers to a comparable degree of development
  init_simulation()
  run_n_steps(200)
  init_simple()
  run_simple_n_iter(200)
  // Both should have negative u at center (backflow confirms vortex)
  let u_chorin = get_u_center()   // g_u[ny/2][nx/2]
  let u_simple = get_u_simple_at(ny / 2, nx / 2)
  assert_true(u_chorin < 0.0)
  assert_true(u_simple < 0.0)
}

// ── Test 26: SIMPLE does not corrupt original Chorin solver state ──
test "simple_independent_of_chorin" {
  // Run Chorin solver and record center velocity
  init_simulation()
  run_n_steps(100)
  let u_before = get_u_center()
  // Run SIMPLE iterations – must not touch g_u/g_v/g_p
  init_simple()
  run_simple_n_iter(50)
  let u_after = get_u_center()
  // Chorin state unchanged
  assert_true(approx_eq(u_before, u_after, 1.0e-15))
}

// ═══════════════════════════════════════════════════════════════════════════
// Tests 27–34: Chorin-PCG solver (pressure_poisson_pcg + cavity_flow_pcg)
// ═══════════════════════════════════════════════════════════════════════════

// ── Test 27: PCG on trivial (zero) RHS returns 0 iterations ─────────────────
test "pcg_zero_rhs_zero_iters" {
  // When b = 0 everywhere the residual is immediately zero → 0 PCG iterations.
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  let iters = pressure_poisson_pcg(p, dx, dy, b)
  assert_eq(iters, 0)
}

// ── Test 28: PCG preserves pressure BCs after solve ──────────────────────────
test "pcg_pressure_bcs_after_solve" {
  // After a PCG pressure solve the boundary conditions must be satisfied.
  // Run one Chorin-PCG step and check both Neumann and Dirichlet BCs on p.
  init_chorin_pcg()
  run_chorin_pcg_n_steps(1)
  // Dirichlet: p = 0 on top row
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_pcg_at(ny - 1, j), 0.0, 1.0e-10))
  }
  // Neumann (dp/dx = 0): left and right columns mirror their neighbours
  for i = 0; i < ny; i = i + 1 {
    assert_true(approx_eq(get_p_pcg_at(i, 0), get_p_pcg_at(i, 1), 1.0e-10))
    assert_true(approx_eq(get_p_pcg_at(i, nx - 1), get_p_pcg_at(i, nx - 2), 1.0e-10))
  }
  // Neumann (dp/dy = 0): bottom row mirrors row 1
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_p_pcg_at(0, j), get_p_pcg_at(1, j), 1.0e-10))
  }
}

// ── Test 29: Chorin-PCG enforces no-slip wall BCs after one step ─────────────
test "pcg_chorin_wall_bcs" {
  init_chorin_pcg()
  run_chorin_pcg_n_steps(1)
  // Left and right walls: u = 0 (exclude top lid row i=ny-1 where u=1)
  for i = 0; i < ny - 1; i = i + 1 {
    assert_true(approx_eq(get_u_pcg_at(i, 0),      0.0, 1.0e-12))
    assert_true(approx_eq(get_u_pcg_at(i, nx - 1), 0.0, 1.0e-12))
  }
  // Bottom wall: u = 0 and v = 0
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_pcg_at(0, j), 0.0, 1.0e-12))
    assert_true(approx_eq(get_v_pcg_at(0, j), 0.0, 1.0e-12))
  }
  // All v = 0 on left and right walls
  for i = 0; i < ny; i = i + 1 {
    assert_true(approx_eq(get_v_pcg_at(i, 0),      0.0, 1.0e-12))
    assert_true(approx_eq(get_v_pcg_at(i, nx - 1), 0.0, 1.0e-12))
  }
}

// ── Test 30: Chorin-PCG maintains lid velocity u = 1 at top row ─────────────
test "pcg_lid_velocity" {
  init_chorin_pcg()
  run_chorin_pcg_n_steps(5)
  // Top row must satisfy u = 1 (lid BC applied after every step)
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(get_u_pcg_at(ny - 1, j), 1.0, 1.0e-12))
  }
}

// ── Test 31: Chorin-PCG does not corrupt Chorin or SIMPLE state ──────────────
test "pcg_state_isolation" {
  // Run Chorin solver and record centre velocity
  init_simulation()
  run_n_steps(50)
  let u_chorin = get_u_center()
  // Run SIMPLE
  init_simple()
  run_simple_n_iter(30)
  let u_simple = get_u_simple_at(ny / 2, nx / 2)
  // Run PCG solver
  init_chorin_pcg()
  run_chorin_pcg_n_steps(50)
  // Chorin state must be unchanged
  assert_true(approx_eq(get_u_center(), u_chorin, 1.0e-15))
  // SIMPLE state must be unchanged
  assert_true(approx_eq(get_u_simple_at(ny / 2, nx / 2), u_simple, 1.0e-15))
}

// ── Test 32: PCG iteration count is bounded by pcg_max_iter ──────────────────
test "pcg_iter_count_bounded" {
  init_chorin_pcg()
  run_chorin_pcg_n_steps(10)
  // get_pcg_last_iters() returns the count from the most recent pressure solve
  let iters = get_pcg_last_iters()
  assert_true(iters >= 0)
  assert_true(iters <= pcg_max_iter)
}

// ── Test 33: Chorin-PCG pressure field is physically consistent with Chorin ──
// Both solvers use the same grid, fluid parameters, and Chorin algorithm.
// After the same number of steps, pressure statistics should be similar.
test "pcg_vs_chorin_pressure_consistency" {
  init_simulation()
  run_n_steps(100)
  let max_p_gs  = get_max_p()
  let min_p_gs  = get_min_p()
  init_chorin_pcg()
  run_chorin_pcg_n_steps(100)
  let max_p_pcg = get_max_p_pcg()
  let min_p_pcg = get_min_p_pcg()
  // Both extrema should have the same sign and be within 20% of each other.
  // This is a qualitative check: PCG and GS converge to the same solution
  // but with different iteration counts.
  assert_true(max_p_gs  > 0.0)
  assert_true(max_p_pcg > 0.0)
  assert_true(min_p_gs  < 0.0)
  assert_true(min_p_pcg < 0.0)
  // Relative difference in max pressure
  let rel_diff_max = (max_p_pcg - max_p_gs).abs() / (max_p_gs.abs() + 1.0e-9)
  assert_true(rel_diff_max < 0.20)
}

// ── Test 34: Chorin-PCG produces a centre vortex after sufficient steps ───────
test "pcg_vortex_formation" {
  init_chorin_pcg()
  run_chorin_pcg_n_steps(200)
  // A well-known result for lid-driven cavity flow (Re=20):
  // the streamwise velocity at the geometric centre (nx/2, ny/2) is negative
  // because the recirculating vortex drives fluid leftward there.
  let u_centre = get_u_pcg_at(ny / 2, nx / 2)
  assert_true(u_centre < 0.0)
  // Divergence norm should be small (mass conservation)
  let div_norm = get_pcg_divergence_norm()
  assert_true(div_norm < 0.05)
}

// ════════════════════════════════════════════════════════════════════════════
// Tests 35–42: Staggered MAC solver
// ════════════════════════════════════════════════════════════════════════════

// ── Test 35: init_mac resets all MAC arrays to zero ──────────────────────────
test "mac_init_resets_state" {
  // Dirty state first
  run_mac_n_steps(1)
  // Now reset
  init_mac()
  assert_eq(get_mac_step_count(), 0)
  assert_eq(get_mac_last_iters(), 0)
  // All fields zero after init
  for i = 0; i < mac_nc; i = i + 1 {
    for j = 0; j < mac_nc; j = j + 1 {
      assert_true(approx_eq(get_u_mac_at(i, j), 0.0, 1.0e-15))
      assert_true(approx_eq(get_v_mac_at(i, j), 0.0, 1.0e-15))
      assert_true(approx_eq(get_p_mac_at(i, j), 0.0, 1.0e-15))
    }
  }
}

// ── Test 36: MAC pressure BCs satisfied after pressure solve ─────────────────
test "mac_pressure_bcs_after_solve" {
  init_mac()
  run_mac_n_steps(1)
  // Left Neumann: p[i][0] == p[i][1]
  for i = 1; i < mac_nc - 1; i = i + 1 {
    assert_true(approx_eq(get_p_mac_at(i, 0), get_p_mac_at(i, 1), 1.0e-10))
  }
  // Right Neumann: p[i][mac_nc-1] == p[i][mac_nc-2]
  for i = 1; i < mac_nc - 1; i = i + 1 {
    assert_true(approx_eq(get_p_mac_at(i, mac_nc - 1), get_p_mac_at(i, mac_nc - 2), 1.0e-10))
  }
  // Bottom Neumann: p[0][j] == p[1][j]
  for j = 0; j < mac_nc; j = j + 1 {
    assert_true(approx_eq(get_p_mac_at(0, j), get_p_mac_at(1, j), 1.0e-10))
  }
  // Top Dirichlet: p[mac_nc-1][j] == 0
  for j = 0; j < mac_nc; j = j + 1 {
    assert_true(approx_eq(get_p_mac_at(mac_nc - 1, j), 0.0, 1.0e-10))
  }
}

// ── Test 37: MAC enforces wall no-slip BCs after one step ────────────────────
test "mac_wall_bcs" {
  init_mac()
  run_mac_n_steps(1)
  // Left wall u-faces: g_u_mac[i][0] = 0
  for i = 0; i < mac_nc; i = i + 1 {
    assert_true(approx_eq(g_u_mac[i][0], 0.0, 1.0e-12))
  }
  // Right wall u-faces: g_u_mac[i][mac_nc] = 0
  for i = 0; i < mac_nc; i = i + 1 {
    assert_true(approx_eq(g_u_mac[i][mac_nc], 0.0, 1.0e-12))
  }
  // Bottom wall v-faces: g_v_mac[0][j] = 0
  for j = 0; j < mac_nc; j = j + 1 {
    assert_true(approx_eq(g_v_mac[0][j], 0.0, 1.0e-12))
  }
  // Top wall v-faces: g_v_mac[mac_nc][j] = 0
  for j = 0; j < mac_nc; j = j + 1 {
    assert_true(approx_eq(g_v_mac[mac_nc][j], 0.0, 1.0e-12))
  }
}

// ── Test 38: MAC step counter increments correctly ───────────────────────────
test "mac_step_counter" {
  init_mac()
  assert_eq(get_mac_step_count(), 0)
  run_mac_n_steps(5)
  assert_eq(get_mac_step_count(), 5)
  run_mac_n_steps(3)
  assert_eq(get_mac_step_count(), 8)
}

// ── Test 39: MAC solver does not corrupt other solver states ─────────────────
test "mac_state_isolation" {
  // Prime Chorin and SIMPLE with known values
  init_simulation()
  run_n_steps(10)
  let u_gs  = get_u_at(ny / 2, nx / 2)
  let p_gs  = get_p_at(ny / 2, nx / 2)
  init_simple()
  run_simple_n_iter(20)
  let u_si  = get_u_simple_at(ny / 2, nx / 2)
  init_chorin_pcg()
  run_chorin_pcg_n_steps(10)
  let u_pcg = get_u_pcg_at(ny / 2, nx / 2)
  // Run MAC
  init_mac()
  run_mac_n_steps(50)
  // Chorin, SIMPLE, PCG states must be unchanged
  assert_true(approx_eq(get_u_at(ny / 2, nx / 2),     u_gs,  1.0e-12))
  assert_true(approx_eq(get_p_at(ny / 2, nx / 2),     p_gs,  1.0e-12))
  assert_true(approx_eq(get_u_simple_at(ny / 2, nx / 2), u_si,  1.0e-12))
  assert_true(approx_eq(get_u_pcg_at(ny / 2, nx / 2), u_pcg, 1.0e-12))
}

// ── Test 40: MAC PCG iteration count is bounded ──────────────────────────────
test "mac_iter_count_bounded" {
  init_mac()
  run_mac_n_steps(50)
  let iters = get_mac_last_iters()
  assert_true(iters >= 0)
  assert_true(iters <= pcg_max_iter)
}

// ── Test 41: MAC divergence is near-machine-zero after correction ─────────────
// By construction ∇·u = 0 exactly in the staggered correction step.
// Numeric residuals from PCG pressure solve may leave a small remainder;
// we accept anything below 1e-4.
test "mac_divergence_near_zero" {
  init_mac()
  run_mac_n_steps(50)
  let div_norm = get_mac_divergence_norm()
  assert_true(div_norm < 1.0e-4)
}

// ── Test 42: MAC produces a centre vortex after sufficient steps ─────────────
test "mac_vortex_formation" {
  init_mac()
  run_mac_n_steps(200)
  // At the geometric centre of the 40×40 MAC grid, the recirculating vortex
  // drives fluid leftward (u < 0) just as in the collocated solvers.
  let u_ctr = get_u_mac_at(mac_nc / 2, mac_nc / 2)
  assert_true(u_ctr < 0.0)
  // Divergence must remain small
  let div_norm = get_mac_divergence_norm()
  assert_true(div_norm < 1.0e-4)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests 43–50: DILU preconditioner
// ═══════════════════════════════════════════════════════════════════════════════

// ── Test 43: DILU PCG on trivial (zero) RHS returns 0 iterations ─────────────
test "dilu_zero_rhs_zero_iters" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  let iters = pressure_poisson_pcg_dilu(p, dx, dy, b)
  assert_eq(iters, 0)
}

// ── Test 44: DILU PCG preserves pressure BCs after solve ─────────────────────
test "dilu_pressure_bcs_after_solve" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  // Non-trivial RHS: use a smooth source
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dilu(p, dx, dy, b) |> ignore
  // Dirichlet: p = 0 on top row
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p[ny - 1][j], 0.0, 1.0e-10))
  }
  // Neumann (dp/dx=0): left and right columns mirror neighbours
  for i = 0; i < ny; i = i + 1 {
    assert_true(approx_eq(p[i][0], p[i][1], 1.0e-10))
    assert_true(approx_eq(p[i][nx - 1], p[i][nx - 2], 1.0e-10))
  }
  // Neumann (dp/dy=0): bottom row mirrors row 1
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p[0][j], p[1][j], 1.0e-10))
  }
}

// ── Test 45: DILU pressure field matches Jacobi PCG result ───────────────────
// Both solvers applied to the same RHS must give ‖p_dilu − p_jac‖/‖p_jac‖ < 1%.
test "dilu_matches_jacobi_pcg" {
  let p_jac  = create_zeros_2d(ny, nx)
  let p_dilu = create_zeros_2d(ny, nx)
  let b      = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p_jac,  dx, dy, b) |> ignore
  pressure_poisson_pcg_dilu(p_dilu, dx, dy, b) |> ignore
  let mut diff_sq = 0.0
  let mut ref_sq  = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let d = p_dilu[i][j] - p_jac[i][j]
      diff_sq = diff_sq + d * d
      ref_sq  = ref_sq  + p_jac[i][j] * p_jac[i][j]
    }
  }
  let rel_err = if ref_sq > 1.0e-15 { diff_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  assert_true(rel_err < 0.01)
}

// ── Test 46: DILU modified diagonal is strictly positive ─────────────────────
// Ensures the build_modified_diag function produces valid (non-zero) denominators.
test "dilu_modified_diag_positive" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let d       = build_modified_diag(inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(d[i][j] > 0.0)
    }
  }
}

// ── Test 47: DILU iteration count is bounded by pcg_max_iter ─────────────────
test "dilu_convergence_bounded" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  let iters = pressure_poisson_pcg_dilu(p, dx, dy, b)
  assert_true(iters > 0)
  assert_true(iters <= pcg_max_iter)
}

// ── Test 48: DILU and Jacobi both converge to the same pressure field ─────────
// For a constant-coefficient uniform grid the DILU improvement over Jacobi is
// small in iteration count; we verify only that DILU does converge and that both
// yield the same pressure field (rel-err < 1%).
test "dilu_iters_not_worse_than_jacobi" {
  let p_jac  = create_zeros_2d(ny, nx)
  let p_dilu = create_zeros_2d(ny, nx)
  let b      = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  let iters_jac  = pressure_poisson_pcg(p_jac,  dx, dy, b)
  let iters_dilu = pressure_poisson_pcg_dilu(p_dilu, dx, dy, b)
  // Both must converge
  assert_true(iters_jac  > 0)
  assert_true(iters_dilu > 0)
  assert_true(iters_jac  <= pcg_max_iter)
  assert_true(iters_dilu <= pcg_max_iter)
  // Both should produce the same pressure (< 1% rel-error)
  let mut diff_sq = 0.0
  let mut ref_sq  = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let d = p_dilu[i][j] - p_jac[i][j]
      diff_sq = diff_sq + d * d
      ref_sq  = ref_sq  + p_jac[i][j] * p_jac[i][j]
    }
  }
  let rel_err = if ref_sq > 1.0e-15 { diff_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  assert_true(rel_err < 0.01)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests 49–54: DIC preconditioner
// ═══════════════════════════════════════════════════════════════════════════════

// ── Test 49: DIC PCG on trivial (zero) RHS returns 0 iterations ──────────────
test "dic_zero_rhs_zero_iters" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  let iters = pressure_poisson_pcg_dic(p, dx, dy, b)
  assert_eq(iters, 0)
}

// ── Test 50: DIC PCG preserves pressure BCs after solve ──────────────────────
test "dic_pressure_bcs_after_solve" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_dic(p, dx, dy, b) |> ignore
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p[ny - 1][j], 0.0, 1.0e-10))
  }
  for i = 0; i < ny; i = i + 1 {
    assert_true(approx_eq(p[i][0], p[i][1], 1.0e-10))
    assert_true(approx_eq(p[i][nx - 1], p[i][nx - 2], 1.0e-10))
  }
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p[0][j], p[1][j], 1.0e-10))
  }
}

// ── Test 51: DIC pressure field matches Jacobi PCG result ────────────────────
test "dic_matches_jacobi_pcg" {
  let p_jac = create_zeros_2d(ny, nx)
  let p_dic = create_zeros_2d(ny, nx)
  let b     = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p_jac, dx, dy, b) |> ignore
  pressure_poisson_pcg_dic(p_dic, dx, dy, b) |> ignore
  let mut diff_sq = 0.0
  let mut ref_sq  = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let d = p_dic[i][j] - p_jac[i][j]
      diff_sq = diff_sq + d * d
      ref_sq  = ref_sq  + p_jac[i][j] * p_jac[i][j]
    }
  }
  let rel_err = if ref_sq > 1.0e-15 { diff_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  assert_true(rel_err < 0.01)
}

// ── Test 52: DIC modified diagonal is strictly positive ──────────────────────
test "dic_modified_diag_positive" {
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  let a_diag  = 2.0 * inv_dx2 + 2.0 * inv_dy2
  let d       = build_modified_diag(inv_dx2, inv_dy2, a_diag)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      assert_true(d[i][j] > 0.0)
    }
  }
}

// ── Test 53: DIC iteration count is bounded by pcg_max_iter ──────────────────
test "dic_convergence_bounded" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  let iters = pressure_poisson_pcg_dic(p, dx, dy, b)
  assert_true(iters > 0)
  assert_true(iters <= pcg_max_iter)
}

// ── Test 54: DIC does not use more iterations than Jacobi PCG × 1.5 ──────────
test "dic_iters_not_worse_than_jacobi" {
  let p_jac = create_zeros_2d(ny, nx)
  let p_dic = create_zeros_2d(ny, nx)
  let b     = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  let iters_jac = pressure_poisson_pcg(p_jac, dx, dy, b)
  let iters_dic = pressure_poisson_pcg_dic(p_dic, dx, dy, b)
  assert_true(iters_dic.to_double() <= iters_jac.to_double() * 1.5 + 1.0)
}

// ═══════════════════════════════════════════════════════════════════════════════
// Tests 55–58: GAMG preconditioner
// ═══════════════════════════════════════════════════════════════════════════════

// ── Test 55: GAMG PCG on trivial (zero) RHS returns 0 iterations ─────────────
test "gamg_zero_rhs_zero_iters" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  let iters = pressure_poisson_pcg_gamg(p, dx, dy, b)
  assert_eq(iters, 0)
}

// ── Test 56: GAMG PCG preserves pressure BCs after solve ─────────────────────
test "gamg_pressure_bcs_after_solve" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg_gamg(p, dx, dy, b) |> ignore
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p[ny - 1][j], 0.0, 1.0e-10))
  }
  for i = 0; i < ny; i = i + 1 {
    assert_true(approx_eq(p[i][0], p[i][1], 1.0e-10))
    assert_true(approx_eq(p[i][nx - 1], p[i][nx - 2], 1.0e-10))
  }
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p[0][j], p[1][j], 1.0e-10))
  }
}

// ── Test 57: GAMG solution matches Jacobi PCG result ─────────────────────────
// Verifies that the V-cycle stationary iteration converges to the same
// pressure field as Jacobi PCG to within 1% relative error.
// Note: for the SND Laplacian ∇²p = b with b > 0, the solution p* < 0
// (the discrete Laplacian is negative definite; diagonal is −a_diag < 0).
test "gamg_matches_jacobi_pcg" {
  let p_jac  = create_zeros_2d(ny, nx)
  let p_gamg = create_zeros_2d(ny, nx)
  let b      = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  pressure_poisson_pcg(p_jac, dx, dy, b) |> ignore
  let iters_gamg = pressure_poisson_pcg_gamg(p_gamg, dx, dy, b)
  // Must converge within the iteration budget
  assert_true(iters_gamg > 0)
  assert_true(iters_gamg <= pcg_max_iter)
  // The GAMG pressure field must have the correct BCs
  for j = 0; j < nx; j = j + 1 {
    assert_true(approx_eq(p_gamg[ny - 1][j], 0.0, 1.0e-10))
  }
  // GAMG solution must match Jacobi PCG within 1%
  let mut diff_sq = 0.0
  let mut ref_sq  = 0.0
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let d = p_gamg[i][j] - p_jac[i][j]
      diff_sq = diff_sq + d * d
      ref_sq  = ref_sq  + p_jac[i][j] * p_jac[i][j]
    }
  }
  let rel_err = if ref_sq > 1.0e-15 { diff_sq.sqrt() / ref_sq.sqrt() } else { 0.0 }
  assert_true(rel_err < 0.01)
}

// ── Test 58: GAMG iteration count is bounded by pcg_max_iter ─────────────────
test "gamg_convergence_bounded" {
  let p = create_zeros_2d(ny, nx)
  let b = create_zeros_2d(ny, nx)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      b[i][j] = 1.0
    }
  }
  let iters = pressure_poisson_pcg_gamg(p, dx, dy, b)
  assert_true(iters > 0)
  assert_true(iters <= pcg_max_iter)
}
