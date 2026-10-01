// Fifth-solver white-box tests. Production code lives in
// solver_projection_tvd.mbt and is intentionally not exported to WASM.

test "projection_upwind_negative_x_uses_positive_index" {
  let phi = create_zeros_2d(4, 4)
  phi[1][1] = 10.0
  phi[1][2] = 20.0
  let face = projection_face_state_x(
    phi,
    1,
    1,
    -1.0,
    projection_scheme_upwind,
    4,
    4,
  )
  assert_eq(face, 20.0)
}

test "projection_upwind_negative_y_uses_positive_index" {
  let phi = create_zeros_2d(4, 4)
  phi[1][1] = 10.0
  phi[2][1] = 20.0
  let face = projection_face_state_y(
    phi,
    1,
    1,
    -1.0,
    projection_scheme_upwind,
    4,
    4,
  )
  assert_eq(face, 20.0)
}

test "projection_convection_uniform_field_is_zero" {
  let nr = 6
  let nc = 6
  let phi = create_zeros_2d(nr, nc)
  let adv = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j] = 3.25
      adv[i][j] = 0.7
    }
  }
  let cx_up = projection_convection_x(
    phi, adv, 1.0, projection_scheme_upwind, nr, nc,
  )
  let cy_up = projection_convection_y(
    phi, adv, 1.0, projection_scheme_upwind, nr, nc,
  )
  let cx_tvd = projection_convection_x(
    phi, adv, 1.0, projection_scheme_tvd_van_leer, nr, nc,
  )
  let cy_tvd = projection_convection_y(
    phi, adv, 1.0, projection_scheme_tvd_van_leer, nr, nc,
  )
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      assert_true(approx_eq(cx_up[i][j], 0.0, 1.0e-12))
      assert_true(approx_eq(cy_up[i][j], 0.0, 1.0e-12))
      assert_true(approx_eq(cx_tvd[i][j], 0.0, 1.0e-12))
      assert_true(approx_eq(cy_tvd[i][j], 0.0, 1.0e-12))
    }
  }
}

test "projection_tvd_linear_x_matches_velocity_sign" {
  let nr = 7
  let nc = 7
  let phi = create_zeros_2d(nr, nc)
  let adv_pos = create_zeros_2d(nr, nc)
  let adv_neg = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j] = j.to_double()
      adv_pos[i][j] = 1.0
      adv_neg[i][j] = -1.0
    }
  }
  let pos = projection_convection_x(
    phi, adv_pos, 1.0, projection_scheme_tvd_van_leer, nr, nc,
  )
  let neg = projection_convection_x(
    phi, adv_neg, 1.0, projection_scheme_tvd_van_leer, nr, nc,
  )
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 2; j < nc - 2; j = j + 1 {
      assert_true(approx_eq(pos[i][j], 1.0, 1.0e-12))
      assert_true(approx_eq(neg[i][j], -1.0, 1.0e-12))
    }
  }
}

test "projection_tvd_linear_y_matches_velocity_sign" {
  let nr = 7
  let nc = 7
  let phi = create_zeros_2d(nr, nc)
  let adv_pos = create_zeros_2d(nr, nc)
  let adv_neg = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      phi[i][j] = i.to_double()
      adv_pos[i][j] = 1.0
      adv_neg[i][j] = -1.0
    }
  }
  let pos = projection_convection_y(
    phi, adv_pos, 1.0, projection_scheme_tvd_van_leer, nr, nc,
  )
  let neg = projection_convection_y(
    phi, adv_neg, 1.0, projection_scheme_tvd_van_leer, nr, nc,
  )
  for i = 2; i < nr - 2; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      assert_true(approx_eq(pos[i][j], 1.0, 1.0e-12))
      assert_true(approx_eq(neg[i][j], -1.0, 1.0e-12))
    }
  }
}

test "projection_discrete_div_grad_matches_laplacian" {
  let nr = 8
  let nc = 9
  let dx_ = 0.25
  let dy_ = 0.20
  let phi = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let x = j.to_double() * dx_
      let y = i.to_double() * dy_
      phi[i][j] = x * x + 0.5 * y * y + 0.1 * x * y
    }
  }
  let gx = projection_gradient_x(phi, dx_, nr, nc)
  let gy = projection_gradient_y(phi, dy_, nr, nc)
  let div_grad = projection_divergence(gx, gy, dx_, dy_, nr, nc)
  let lap = projection_laplacian(phi, dx_, dy_, nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      assert_true(approx_eq(div_grad[i][j], lap[i][j], 1.0e-12))
    }
  }
}

test "projection_remove_interior_mean_is_compatible" {
  let nr = 7
  let nc = 6
  let rhs = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      rhs[i][j] = 2.0 * i.to_double() - 3.0 * j.to_double() + 7.0
    }
  }
  let removed = projection_remove_interior_mean(rhs, nr, nc)
  let mut sum = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      sum = sum + rhs[i][j]
    }
  }
  assert_true(removed.abs() > 1.0e-12)
  assert_true(sum.abs() < 1.0e-14)
}

test "projection_velocity_bc_applies_walls_then_lid" {
  let nr = 6
  let nc = 7
  let u = create_zeros_2d(nr, nc)
  let v = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      u[i][j] = 3.0
      v[i][j] = -2.0
    }
  }
  projection_apply_velocity_bc(u, v, nr, nc)
  for j = 0; j < nc; j = j + 1 {
    assert_eq(u[0][j], 0.0)
    assert_eq(v[0][j], 0.0)
    assert_eq(u[nr - 1][j], 1.0)
    assert_eq(v[nr - 1][j], 0.0)
  }
  for i = 0; i < nr - 1; i = i + 1 {
    assert_eq(u[i][0], 0.0)
    assert_eq(v[i][0], 0.0)
    assert_eq(u[i][nc - 1], 0.0)
    assert_eq(v[i][nc - 1], 0.0)
  }
}

test "projection_neumann_bc_copies_adjacent_values" {
  let nr = 6
  let nc = 7
  let phi = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      phi[i][j] = (10 * i + j).to_double()
    }
  }
  projection_apply_neumann_bc(phi, nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    assert_eq(phi[i][0], phi[i][1])
    assert_eq(phi[i][nc - 1], phi[i][nc - 2])
  }
  for j = 1; j < nc - 1; j = j + 1 {
    assert_eq(phi[0][j], phi[1][j])
    assert_eq(phi[nr - 1][j], phi[nr - 2][j])
  }
}

test "projection_pressure_sor_reaches_manufactured_residual" {
  let nr = 9
  let nc = 10
  let dx_ = 0.20
  let dy_ = 0.25
  let exact = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let x = j.to_double() * dx_
      let y = i.to_double() * dy_
      exact[i][j] = x * x * y + 0.3 * x * y * y
    }
  }
  projection_apply_neumann_bc(exact, nr, nc)
  projection_remove_interior_mean(exact, nr, nc) |> ignore
  projection_apply_neumann_bc(exact, nr, nc)
  let rhs = projection_laplacian(exact, dx_, dy_, nr, nc)
  projection_remove_interior_mean(rhs, nr, nc) |> ignore
  let phi = create_zeros_2d(nr, nc)
  let result = projection_pressure_sor(
    phi, rhs, dx_, dy_, nr, nc, 1.0e-8, 20000, 1.7,
  )
  assert_true(result.converged)
  assert_true(result.residual <= 1.0e-8)
  assert_true(result.iterations > 0)
  assert_true(result.iterations < 20000)
  let mut mean = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      mean = mean + phi[i][j]
    }
  }
  mean = mean / ((nr - 2) * (nc - 2)).to_double()
  assert_true(mean.abs() < 1.0e-12)
}

test "projection_pressure_settings_validation" {
  assert_true(projection_pressure_settings_valid(1.0e-6, 1000, 1.7))
  assert_false(projection_pressure_settings_valid(0.0, 1000, 1.7))
  assert_false(projection_pressure_settings_valid(-1.0e-6, 1000, 1.7))
  assert_false(projection_pressure_settings_valid(1.0e-6, 0, 1.7))
  assert_false(projection_pressure_settings_valid(1.0e-6, 1000, 0.0))
  assert_false(projection_pressure_settings_valid(1.0e-6, 1000, 2.0))
}

test "projection_correction_reduces_compatible_divergence_100x" {
  let nr = 9
  let nc = 9
  let dx_ = 0.25
  let dy_ = 0.25
  let psi = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let x = j.to_double() * dx_
      let y = i.to_double() * dy_
      psi[i][j] = x * x * y + 0.2 * x * y * y
    }
  }
  projection_apply_neumann_bc(psi, nr, nc)
  projection_remove_interior_mean(psi, nr, nc) |> ignore
  projection_apply_neumann_bc(psi, nr, nc)
  let u = projection_gradient_x(psi, dx_, nr, nc)
  let v = projection_gradient_y(psi, dy_, nr, nc)
  projection_apply_velocity_bc(u, v, nr, nc)
  let p = create_zeros_2d(nr, nc)
  let before = projection_divergence_l2(u, v, dx_, dy_, nr, nc)
  let config : ProjectionTvdConfig = {
    scheme: projection_scheme_upwind,
    dt: 0.01,
    pressure_tol: 1.0e-10,
    pressure_max_iter: 20000,
    pressure_omega: 1.7,
  }
  let diagnostics = projection_correct_predictor(
    u, v, p, 1.0, dx_, dy_, nr, nc, config,
  )
  assert_true(diagnostics.pressure_converged)
  assert_true(diagnostics.finite)
  assert_true(before > 1.0e-8)
  assert_true(diagnostics.divergence_norm < before * 0.01)
}

test "projection_correction_reduces_arbitrary_wall_bounded_divergence_100x" {
  let nr = 9
  let nc = 9
  let dx_ = 0.25
  let dy_ = 0.25
  let u = create_zeros_2d(nr, nc)
  let v = create_zeros_2d(nr, nc)
  let p = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let x = j.to_double() * dx_
      let y = i.to_double() * dy_
      u[i][j] = x * (2.0 - x) * (0.2 + y)
      v[i][j] = y * (2.0 - y) * (0.1 - 0.3 * x)
    }
  }
  projection_apply_velocity_bc(u, v, nr, nc)
  let before = projection_divergence_l2(u, v, dx_, dy_, nr, nc)
  let config : ProjectionTvdConfig = {
    scheme: projection_scheme_upwind,
    dt: 0.01,
    pressure_tol: 1.0e-10,
    pressure_max_iter: 20000,
    pressure_omega: 1.7,
  }
  let diagnostics = projection_correct_predictor(
    u, v, p, 1.0, dx_, dy_, nr, nc, config,
  )
  assert_true(diagnostics.pressure_converged)
  assert_true(before > 1.0e-8)
  assert_true(diagnostics.divergence_norm < before * 0.01)
}

test "projection_tvd_step_keeps_boundary_and_fields_finite" {
  let nr = 11
  let nc = 11
  let dx_ = 0.20
  let dy_ = 0.20
  let u = create_zeros_2d(nr, nc)
  let v = create_zeros_2d(nr, nc)
  let p = create_zeros_2d(nr, nc)
  projection_apply_velocity_bc(u, v, nr, nc)
  let config : ProjectionTvdConfig = {
    scheme: projection_scheme_tvd_van_leer,
    dt: 0.001,
    pressure_tol: 1.0e-8,
    pressure_max_iter: 5000,
    pressure_omega: 1.7,
  }
  let diagnostics = projection_tvd_step(
    u, v, p, 1.0, 0.02, dx_, dy_, nr, nc, config,
  )
  assert_true(diagnostics.pressure_converged)
  assert_true(diagnostics.finite)
  for j = 0; j < nc; j = j + 1 {
    assert_eq(u[0][j], 0.0)
    assert_eq(v[0][j], 0.0)
    assert_eq(u[nr - 1][j], 1.0)
    assert_eq(v[nr - 1][j], 0.0)
  }
  for i = 0; i < nr - 1; i = i + 1 {
    assert_eq(u[i][0], 0.0)
    assert_eq(v[i][0], 0.0)
    assert_eq(u[i][nc - 1], 0.0)
    assert_eq(v[i][nc - 1], 0.0)
  }
}

test "projection_diagnostics_reject_nonfinite_derived_quantities" {
  let nr = 5
  let nc = 5
  let u = create_zeros_2d(nr, nc)
  let v = create_zeros_2d(nr, nc)
  let p = create_zeros_2d(nr, nc)
  u[2][2] = 1.0e240
  let pressure_result : ProjectionPressureResult = {
    iterations: 1,
    residual: 1.0e-7,
    converged: true,
  }
  let diagnostics = projection_compute_diagnostics(
    u, v, p, pressure_result, 0.01, 0.5, 0.5, nr, nc,
  )
  assert_false(diagnostics.finite)
}

test "projection_tvd_runtime_state_is_isolated_and_collectable" {
  let config : ProjectionTvdConfig = {
    scheme: projection_scheme_upwind,
    dt: 0.001,
    pressure_tol: 1.0e-6,
    pressure_max_iter: 5000,
    pressure_omega: 1.7,
  }
  init_projection_tvd(13, config, 100.0)
  assert_eq(g_projection_steps[0], 0)
  assert_eq(g_projection_u[0].length(), 13)
  assert_eq(g_projection_u[0][0].length(), 13)
  assert_eq(g_projection_phi[0].length(), 13)
  run_projection_tvd_steps(2)
  assert_eq(g_projection_steps[0], 2)
  let diagnostics = get_projection_tvd_diagnostics()
  assert_true(diagnostics.finite)
  assert_true(diagnostics.pressure_converged)
  let solver_data = collect_projection_tvd_export()
  assert_eq(solver_data.name, "projection-tvd")
  assert_eq(solver_data.nc, 13)
  assert_eq(solver_data.steps, 2)
  let mut phi_norm = 0.0
  for i = 1; i < 12; i = i + 1 {
    for j = 1; j < 12; j = j + 1 {
      phi_norm = phi_norm + g_projection_phi[0][i][j].abs()
    }
  }
  assert_true(phi_norm > 0.0)
}

test "projection_cli_parses_positive_doubles" {
  assert_true(approx_eq(parse_positive_double("0.001", -1.0), 0.001, 1.0e-15))
  assert_true(approx_eq(parse_positive_double("1e-6", -1.0), 1.0e-6, 1.0e-18))
  assert_true(approx_eq(parse_positive_double("1.7", -1.0), 1.7, 1.0e-15))
  assert_eq(parse_positive_double("0", -1.0), -1.0)
  assert_eq(parse_positive_double("-0.1", -1.0), -1.0)
  assert_eq(parse_positive_double("bad", -1.0), -1.0)
}

test "projection_solver_is_explicit_not_part_of_all" {
  assert_true(cli_solver_selected("all", "chorin"))
  assert_true(cli_solver_selected("all", "simple"))
  assert_true(cli_solver_selected("all", "pcg"))
  assert_true(cli_solver_selected("all", "mac"))
  assert_false(cli_solver_selected("all", "projection-tvd"))
  assert_true(cli_solver_selected("projection-tvd", "projection-tvd"))
}

test "projection_centerline_formatter_is_self_describing" {
  let config : ProjectionTvdConfig = {
    scheme: projection_scheme_upwind,
    dt: 0.001,
    pressure_tol: 1.0e-6,
    pressure_max_iter: 1000,
    pressure_omega: 1.7,
  }
  init_projection_tvd(9, config, 100.0)
  run_projection_tvd_steps(1)
  let text = format_projection_tvd_centerline_csv(
    collect_projection_tvd_export(),
    100.0,
    config,
    get_projection_tvd_diagnostics(),
  )
  assert_true(text.contains("solver=projection-tvd"))
  assert_true(text.contains("scheme=upwind"))
  assert_true(text.contains("grid=9x9"))
  assert_true(text.contains("dt=0.001"))
  assert_true(text.contains("pressure_tol=0.000001"))
  assert_true(text.contains("pressure_max_iter=1000"))
  assert_true(text.contains("pressure_omega=1.7"))
  assert_true(text.contains("pressure_iters="))
  assert_true(text.contains("pressure_residual="))
  assert_true(text.contains("pressure_converged=true"))
  assert_true(text.contains("div_norm="))
  assert_true(text.contains("cfl_max="))
  assert_true(text.contains("finite=true"))
}
