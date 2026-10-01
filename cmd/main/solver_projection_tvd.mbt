// FlowLabLite fifth solver: sign-aware convection and incremental projection.
// This file is CLI-only and deliberately adds no WASM exports.

let projection_scheme_upwind = 1
let projection_scheme_tvd_van_leer = 2

fn projection_van_leer(r : Double) -> Double {
  (r + r.abs()) / (1.0 + r.abs())
}

// Return the transported state at the x-face between columns `left` and
// `left + 1`. A negative face velocity is upstream on the positive-index side.
fn projection_face_state_x(
  phi : Array[Array[Double]],
  row : Int,
  left : Int,
  face_velocity : Double,
  scheme : Int,
  nr : Int,
  nc : Int,
) -> Double {
  let _nr = nr
  if row < 0 || row >= phi.length() || left < 0 || left + 1 >= nc {
    return 0.0
  }
  let upstream_col = if face_velocity >= 0.0 { left } else { left + 1 }
  let downstream_col = if face_velocity >= 0.0 { left + 1 } else { left }
  let next_upstream_col = if face_velocity >= 0.0 { left - 1 } else { left + 2 }
  let upwind = phi[row][upstream_col]
  if scheme != projection_scheme_tvd_van_leer ||
    next_upstream_col < 0 || next_upstream_col >= nc {
    return upwind
  }
  let downstream_delta = phi[row][downstream_col] - upwind
  if downstream_delta.abs() < 1.0e-30 { return upwind }
  let r = (upwind - phi[row][next_upstream_col]) / downstream_delta
  upwind + 0.5 * projection_van_leer(r) * downstream_delta
}

// Return the transported state at the y-face between rows `lower` and
// `lower + 1`. A negative face velocity is upstream on the positive-index side.
fn projection_face_state_y(
  phi : Array[Array[Double]],
  lower : Int,
  col : Int,
  face_velocity : Double,
  scheme : Int,
  nr : Int,
  nc : Int,
) -> Double {
  let _nc = nc
  if lower < 0 || lower + 1 >= nr || col < 0 || col >= phi[0].length() {
    return 0.0
  }
  let upstream_row = if face_velocity >= 0.0 { lower } else { lower + 1 }
  let downstream_row = if face_velocity >= 0.0 { lower + 1 } else { lower }
  let next_upstream_row = if face_velocity >= 0.0 { lower - 1 } else { lower + 2 }
  let upwind = phi[upstream_row][col]
  if scheme != projection_scheme_tvd_van_leer ||
    next_upstream_row < 0 || next_upstream_row >= nr {
    return upwind
  }
  let downstream_delta = phi[downstream_row][col] - upwind
  if downstream_delta.abs() < 1.0e-30 { return upwind }
  let r = (upwind - phi[next_upstream_row][col]) / downstream_delta
  upwind + 0.5 * projection_van_leer(r) * downstream_delta
}

fn projection_convection_x(
  phi : Array[Array[Double]],
  adv : Array[Array[Double]],
  dx_ : Double,
  scheme : Int,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let velocity_right = 0.5 * (adv[i][j] + adv[i][j + 1])
      let velocity_left = 0.5 * (adv[i][j - 1] + adv[i][j])
      let state_right = projection_face_state_x(
        phi, i, j, velocity_right, scheme, nr, nc,
      )
      let state_left = projection_face_state_x(
        phi, i, j - 1, velocity_left, scheme, nr, nc,
      )
      result[i][j] = (
        velocity_right * state_right - velocity_left * state_left
      ) / dx_
    }
  }
  result
}

fn projection_convection_y(
  phi : Array[Array[Double]],
  adv : Array[Array[Double]],
  dy_ : Double,
  scheme : Int,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let velocity_top = 0.5 * (adv[i][j] + adv[i + 1][j])
      let velocity_bottom = 0.5 * (adv[i - 1][j] + adv[i][j])
      let state_top = projection_face_state_y(
        phi, i, j, velocity_top, scheme, nr, nc,
      )
      let state_bottom = projection_face_state_y(
        phi, i - 1, j, velocity_bottom, scheme, nr, nc,
      )
      result[i][j] = (
        velocity_top * state_top - velocity_bottom * state_bottom
      ) / dy_
    }
  }
  result
}

// Conservative transport driven by the already projected normal face fluxes.
// Keeping these fluxes as state avoids rebuilding a non-solenoidal advection
// field from the collocated momentum values at the next time step.
fn projection_convection_x_from_faces(
  phi : Array[Array[Double]],
  u_face : Array[Array[Double]],
  dx_ : Double,
  scheme : Int,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let velocity_right = u_face[i][j + 1]
      let velocity_left = u_face[i][j]
      let state_right = projection_face_state_x(
        phi, i, j, velocity_right, scheme, nr, nc,
      )
      let state_left = projection_face_state_x(
        phi, i, j - 1, velocity_left, scheme, nr, nc,
      )
      result[i][j] = (
        velocity_right * state_right - velocity_left * state_left
      ) / dx_
    }
  }
  result
}

fn projection_convection_y_from_faces(
  phi : Array[Array[Double]],
  v_face : Array[Array[Double]],
  dy_ : Double,
  scheme : Int,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let velocity_top = v_face[i + 1][j]
      let velocity_bottom = v_face[i][j]
      let state_top = projection_face_state_y(
        phi, i, j, velocity_top, scheme, nr, nc,
      )
      let state_bottom = projection_face_state_y(
        phi, i - 1, j, velocity_bottom, scheme, nr, nc,
      )
      result[i][j] = (
        velocity_top * state_top - velocity_bottom * state_bottom
      ) / dy_
    }
  }
  result
}

// Forward pressure gradient. Together with the backward divergence below,
// D(G(phi)) is the five-point Laplacian on every interior node.
fn projection_gradient_x(
  phi : Array[Array[Double]],
  dx_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc - 1; j = j + 1 {
      result[i][j] = (phi[i][j + 1] - phi[i][j]) / dx_
    }
  }
  result
}

fn projection_gradient_y(
  phi : Array[Array[Double]],
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 0; i < nr - 1; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      result[i][j] = (phi[i + 1][j] - phi[i][j]) / dy_
    }
  }
  result
}

fn projection_divergence(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = (u[i][j] - u[i][j - 1]) / dx_ +
        (v[i][j] - v[i - 1][j]) / dy_
    }
  }
  result
}

fn projection_laplacian(
  phi : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  let inv_dx2 = 1.0 / (dx_ * dx_)
  let inv_dy2 = 1.0 / (dy_ * dy_)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]) * inv_dx2 +
        (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j]) * inv_dy2
    }
  }
  result
}

fn projection_remove_interior_mean(
  field : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Double {
  let count = (nr - 2) * (nc - 2)
  if count <= 0 { return 0.0 }
  let mut sum = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      sum = sum + field[i][j]
    }
  }
  let mean = sum / count.to_double()
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      field[i][j] = field[i][j] - mean
    }
  }
  mean
}

fn projection_apply_velocity_bc(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Unit {
  for i = 0; i < nr; i = i + 1 {
    u[i][0] = 0.0
    u[i][nc - 1] = 0.0
    v[i][0] = 0.0
    v[i][nc - 1] = 0.0
  }
  for j = 0; j < nc; j = j + 1 {
    u[0][j] = 0.0
    u[nr - 1][j] = 0.0
    v[0][j] = 0.0
    v[nr - 1][j] = 0.0
  }
  // Apply the moving lid last, including both top corners.
  for j = 0; j < nc; j = j + 1 {
    u[nr - 1][j] = 1.0
  }
}

// Build normal velocities on finite-volume faces. Face index 1 and the last
// index are physical wall faces and remain exactly zero; only internal faces
// are interpolated from the collocated predictor.
fn projection_build_face_velocities(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let u_face = create_zeros_2d(nr, nc)
  let v_face = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 2; j < nc - 1; j = j + 1 {
      u_face[i][j] = 0.5 * (u[i][j - 1] + u[i][j])
    }
  }
  for i = 2; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      v_face[i][j] = 0.5 * (v[i - 1][j] + v[i][j])
    }
  }
  (u_face, v_face)
}

fn projection_face_divergence(
  u_face : Array[Array[Double]],
  v_face : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = (u_face[i][j + 1] - u_face[i][j]) / dx_ +
        (v_face[i + 1][j] - v_face[i][j]) / dy_
    }
  }
  result
}

fn projection_face_divergence_l2(
  u_face : Array[Array[Double]],
  v_face : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Double {
  let divergence = projection_face_divergence(
    u_face, v_face, dx_, dy_, nr, nc,
  )
  let count = (nr - 2) * (nc - 2)
  if count <= 0 { return 0.0 }
  let mut sum = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      sum = sum + divergence[i][j] * divergence[i][j]
    }
  }
  (sum / count.to_double()).sqrt()
}

fn projection_reconstruct_cell_velocities(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  u_face : Array[Array[Double]],
  v_face : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Unit {
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      u[i][j] = 0.5 * (u_face[i][j] + u_face[i][j + 1])
      v[i][j] = 0.5 * (v_face[i][j] + v_face[i + 1][j])
    }
  }
  projection_apply_velocity_bc(u, v, nr, nc)
}

fn projection_cell_gradient_x(
  phi : Array[Array[Double]],
  dx_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = (phi[i][j + 1] - phi[i][j - 1]) / (2.0 * dx_)
    }
  }
  result
}

fn projection_cell_gradient_y(
  phi : Array[Array[Double]],
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = (phi[i + 1][j] - phi[i - 1][j]) / (2.0 * dy_)
    }
  }
  result
}

fn projection_apply_neumann_bc(
  phi : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Unit {
  for i = 1; i < nr - 1; i = i + 1 {
    phi[i][0] = phi[i][1]
    phi[i][nc - 1] = phi[i][nc - 2]
  }
  for j = 1; j < nc - 1; j = j + 1 {
    phi[0][j] = phi[1][j]
    phi[nr - 1][j] = phi[nr - 2][j]
  }
  phi[0][0] = phi[1][1]
  phi[0][nc - 1] = phi[1][nc - 2]
  phi[nr - 1][0] = phi[nr - 2][1]
  phi[nr - 1][nc - 1] = phi[nr - 2][nc - 2]
}

struct ProjectionPressureResult {
  iterations : Int
  residual : Double
  converged : Bool
}

fn projection_pressure_settings_valid(
  tol : Double,
  max_iter : Int,
  omega : Double,
) -> Bool {
  tol > 0.0 && max_iter > 0 && omega > 0.0 && omega < 2.0
}

fn projection_pressure_relative_residual(
  phi : Array[Array[Double]],
  rhs : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Double {
  let lap = projection_laplacian(phi, dx_, dy_, nr, nc)
  let mut residual_sum = 0.0
  let mut rhs_sum = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let delta = lap[i][j] - rhs[i][j]
      residual_sum = residual_sum + delta * delta
      rhs_sum = rhs_sum + rhs[i][j] * rhs[i][j]
    }
  }
  let rhs_norm = rhs_sum.sqrt()
  let denominator = if rhs_norm > 1.0e-30 { rhs_norm } else { 1.0e-30 }
  residual_sum.sqrt() / denominator
}

fn projection_pressure_sor(
  phi : Array[Array[Double]],
  rhs : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
  tol : Double,
  max_iter : Int,
  omega : Double,
) -> ProjectionPressureResult {
  if !projection_pressure_settings_valid(tol, max_iter, omega) {
    return { iterations: 0, residual: 1.0e300, converged: false }
  }
  let dx2 = dx_ * dx_
  let dy2 = dy_ * dy_
  let denominator = 2.0 * (dx2 + dy2)
  let mut residual = projection_pressure_relative_residual(
    phi, rhs, dx_, dy_, nr, nc,
  )
  if residual <= tol {
    return { iterations: 0, residual, converged: true }
  }
  for iteration = 1; iteration <= max_iter; iteration = iteration + 1 {
    for i = 1; i < nr - 1; i = i + 1 {
      for j = 1; j < nc - 1; j = j + 1 {
        let gs_value = (
          (phi[i][j + 1] + phi[i][j - 1]) * dy2 +
          (phi[i + 1][j] + phi[i - 1][j]) * dx2 -
          rhs[i][j] * dx2 * dy2
        ) / denominator
        phi[i][j] = (1.0 - omega) * phi[i][j] + omega * gs_value
      }
    }
    projection_apply_neumann_bc(phi, nr, nc)
    projection_remove_interior_mean(phi, nr, nc) |> ignore
    projection_apply_neumann_bc(phi, nr, nc)
    residual = projection_pressure_relative_residual(
      phi, rhs, dx_, dy_, nr, nc,
    )
    if residual != residual || residual.abs() > 1.0e300 {
      return { iterations: iteration, residual, converged: false }
    }
    if residual <= tol {
      return { iterations: iteration, residual, converged: true }
    }
  }
  { iterations: max_iter, residual, converged: false }
}

struct ProjectionTvdConfig {
  scheme : Int
  dt : Double
  pressure_tol : Double
  pressure_max_iter : Int
  pressure_omega : Double
}

struct ProjectionTvdDiagnostics {
  pressure_iterations : Int
  pressure_residual : Double
  pressure_converged : Bool
  divergence_norm : Double
  cfl_max : Double
  finite : Bool
}

fn projection_divergence_l2(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> Double {
  let divergence = projection_divergence(u, v, dx_, dy_, nr, nc)
  let count = (nr - 2) * (nc - 2)
  if count <= 0 { return 0.0 }
  let mut sum = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      sum = sum + divergence[i][j] * divergence[i][j]
    }
  }
  (sum / count.to_double()).sqrt()
}

fn projection_fields_finite(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Bool {
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let u_value = u[i][j]
      let v_value = v[i][j]
      let p_value = p[i][j]
      if u_value != u_value || v_value != v_value || p_value != p_value ||
        u_value.abs() > 1.0e300 || v_value.abs() > 1.0e300 ||
        p_value.abs() > 1.0e300 {
        return false
      }
    }
  }
  true
}

fn projection_compute_diagnostics(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  pressure_result : ProjectionPressureResult,
  dt_ : Double,
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
) -> ProjectionTvdDiagnostics {
  let (u_face, v_face) = projection_build_face_velocities(u, v, nr, nc)
  projection_compute_diagnostics_with_divergence(
    u,
    v,
    p,
    pressure_result,
    dt_,
    dx_,
    dy_,
    nr,
    nc,
    projection_face_divergence_l2(u_face, v_face, dx_, dy_, nr, nc),
  )
}

fn projection_compute_diagnostics_with_divergence(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  pressure_result : ProjectionPressureResult,
  dt_ : Double,
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
  divergence_norm : Double,
) -> ProjectionTvdDiagnostics {
  let mut cfl_max = 0.0
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let local_cfl = dt_ * (u[i][j].abs() / dx_ + v[i][j].abs() / dy_)
      if local_cfl > cfl_max { cfl_max = local_cfl }
    }
  }
  let derived_quantities_finite =
    divergence_norm == divergence_norm && divergence_norm.abs() < 1.0e300 &&
    cfl_max == cfl_max && cfl_max.abs() < 1.0e300 &&
    pressure_result.residual == pressure_result.residual &&
    pressure_result.residual.abs() < 1.0e300
  {
    pressure_iterations: pressure_result.iterations,
    pressure_residual: pressure_result.residual,
    pressure_converged: pressure_result.converged,
    divergence_norm,
    cfl_max,
    finite: projection_fields_finite(u, v, p, nr, nc) &&
      derived_quantities_finite,
  }
}

fn projection_correct_predictor_with_faces(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  projected_u_face : Array[Array[Double]],
  projected_v_face : Array[Array[Double]],
  pressure_correction : Array[Array[Double]],
  rho_ : Double,
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
  config : ProjectionTvdConfig,
) -> ProjectionTvdDiagnostics {
  projection_apply_velocity_bc(u, v, nr, nc)
  let (u_face, v_face) = projection_build_face_velocities(u, v, nr, nc)
  let divergence = projection_face_divergence(
    u_face, v_face, dx_, dy_, nr, nc,
  )
  let predictor_divergence = projection_face_divergence_l2(
    u_face, v_face, dx_, dy_, nr, nc,
  )
  let rhs = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      rhs[i][j] = (rho_ / config.dt) * divergence[i][j]
    }
  }
  projection_remove_interior_mean(rhs, nr, nc) |> ignore
  let phi = pressure_correction
  projection_apply_neumann_bc(phi, nr, nc)
  projection_remove_interior_mean(phi, nr, nc) |> ignore
  projection_apply_neumann_bc(phi, nr, nc)
  let pressure_result = projection_pressure_sor(
    phi,
    rhs,
    dx_,
    dy_,
    nr,
    nc,
    config.pressure_tol,
    config.pressure_max_iter,
    config.pressure_omega,
  )
  if !pressure_result.converged {
    return projection_compute_diagnostics_with_divergence(
      u, v, p, pressure_result, config.dt, dx_, dy_, nr, nc,
      predictor_divergence,
    )
  }
  let correction_scale = config.dt / rho_
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 2; j < nc - 1; j = j + 1 {
      u_face[i][j] = u_face[i][j] - correction_scale *
        (phi[i][j] - phi[i][j - 1]) / dx_
    }
  }
  for i = 2; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      v_face[i][j] = v_face[i][j] - correction_scale *
        (phi[i][j] - phi[i - 1][j]) / dy_
    }
  }
  let corrected_divergence = projection_face_divergence_l2(
    u_face, v_face, dx_, dy_, nr, nc,
  )
  copy_2d_array(u_face, projected_u_face)
  copy_2d_array(v_face, projected_v_face)
  let gradient_x = projection_cell_gradient_x(phi, dx_, nr, nc)
  let gradient_y = projection_cell_gradient_y(phi, dy_, nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      u[i][j] = u[i][j] - correction_scale * gradient_x[i][j]
      v[i][j] = v[i][j] - correction_scale * gradient_y[i][j]
      p[i][j] = p[i][j] + phi[i][j]
    }
  }
  projection_apply_velocity_bc(u, v, nr, nc)
  projection_remove_interior_mean(p, nr, nc) |> ignore
  projection_apply_neumann_bc(p, nr, nc)
  projection_compute_diagnostics_with_divergence(
    u, v, p, pressure_result, config.dt, dx_, dy_, nr, nc,
    corrected_divergence,
  )
}

fn projection_correct_predictor(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  rho_ : Double,
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
  config : ProjectionTvdConfig,
) -> ProjectionTvdDiagnostics {
  let projected_u_face = create_zeros_2d(nr, nc)
  let projected_v_face = create_zeros_2d(nr, nc)
  let pressure_correction = create_zeros_2d(nr, nc)
  projection_correct_predictor_with_faces(
    u, v, p, projected_u_face, projected_v_face, pressure_correction,
    rho_, dx_, dy_, nr, nc, config,
  )
}

fn projection_tvd_step_with_faces(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  projected_u_face : Array[Array[Double]],
  projected_v_face : Array[Array[Double]],
  pressure_correction : Array[Array[Double]],
  rho_ : Double,
  nu_ : Double,
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
  config : ProjectionTvdConfig,
) -> ProjectionTvdDiagnostics {
  let previous_u = create_zeros_2d(nr, nc)
  let previous_v = create_zeros_2d(nr, nc)
  let previous_p = create_zeros_2d(nr, nc)
  copy_2d_array(u, previous_u)
  copy_2d_array(v, previous_v)
  copy_2d_array(p, previous_p)
  let convection_ux = projection_convection_x_from_faces(
    previous_u, projected_u_face, dx_, config.scheme, nr, nc,
  )
  let convection_uy = projection_convection_y_from_faces(
    previous_u, projected_v_face, dy_, config.scheme, nr, nc,
  )
  let convection_vx = projection_convection_x_from_faces(
    previous_v, projected_u_face, dx_, config.scheme, nr, nc,
  )
  let convection_vy = projection_convection_y_from_faces(
    previous_v, projected_v_face, dy_, config.scheme, nr, nc,
  )
  let pressure_gradient_x = projection_cell_gradient_x(previous_p, dx_, nr, nc)
  let pressure_gradient_y = projection_cell_gradient_y(previous_p, dy_, nr, nc)
  let laplacian_u = projection_laplacian(previous_u, dx_, dy_, nr, nc)
  let laplacian_v = projection_laplacian(previous_v, dx_, dy_, nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      u[i][j] = previous_u[i][j] - config.dt * (
        convection_ux[i][j] + convection_uy[i][j]
      ) - (config.dt / rho_) * pressure_gradient_x[i][j] +
        config.dt * nu_ * laplacian_u[i][j]
      v[i][j] = previous_v[i][j] - config.dt * (
        convection_vx[i][j] + convection_vy[i][j]
      ) - (config.dt / rho_) * pressure_gradient_y[i][j] +
        config.dt * nu_ * laplacian_v[i][j]
    }
  }
  projection_apply_velocity_bc(u, v, nr, nc)
  projection_correct_predictor_with_faces(
    u, v, p, projected_u_face, projected_v_face, pressure_correction,
    rho_, dx_, dy_, nr, nc, config,
  )
}

fn projection_tvd_step(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
  rho_ : Double,
  nu_ : Double,
  dx_ : Double,
  dy_ : Double,
  nr : Int,
  nc : Int,
  config : ProjectionTvdConfig,
) -> ProjectionTvdDiagnostics {
  let (projected_u_face, projected_v_face) =
    projection_build_face_velocities(u, v, nr, nc)
  let pressure_correction = create_zeros_2d(nr, nc)
  projection_tvd_step_with_faces(
    u, v, p, projected_u_face, projected_v_face, pressure_correction,
    rho_, nu_, dx_, dy_, nr, nc, config,
  )
}

let g_projection_u : Array[Array[Array[Double]]] = [create_zeros_2d(41, 41)]
let g_projection_v : Array[Array[Array[Double]]] = [create_zeros_2d(41, 41)]
let g_projection_p : Array[Array[Array[Double]]] = [create_zeros_2d(41, 41)]
let g_projection_u_face : Array[Array[Array[Double]]] = [create_zeros_2d(41, 41)]
let g_projection_v_face : Array[Array[Array[Double]]] = [create_zeros_2d(41, 41)]
let g_projection_phi : Array[Array[Array[Double]]] = [create_zeros_2d(41, 41)]
let g_projection_nc : Array[Int] = [41]
let g_projection_dx : Array[Double] = [0.05]
let g_projection_re : Array[Double] = [100.0]
let g_projection_nu : Array[Double] = [0.02]
let g_projection_steps : Array[Int] = [0]
let g_projection_config : Array[ProjectionTvdConfig] = [{
  scheme: projection_scheme_tvd_van_leer,
  dt: 0.001,
  pressure_tol: 1.0e-6,
  pressure_max_iter: 1000,
  pressure_omega: 1.7,
}]
let g_projection_last_diagnostics : Array[ProjectionTvdDiagnostics] = [{
  pressure_iterations: 0,
  pressure_residual: 0.0,
  pressure_converged: true,
  divergence_norm: 0.0,
  cfl_max: 0.0,
  finite: true,
}]

fn init_projection_tvd(
  requested_nc : Int,
  config : ProjectionTvdConfig,
  re_ : Double,
) -> Unit {
  let nc_ = if requested_nc >= 3 { requested_nc } else { 3 }
  g_projection_nc[0] = nc_
  g_projection_dx[0] = 2.0 / (nc_ - 1).to_double()
  g_projection_re[0] = re_
  g_projection_nu[0] = 2.0 / re_
  g_projection_steps[0] = 0
  g_projection_config[0] = config
  g_projection_u[0] = create_zeros_2d(nc_, nc_)
  g_projection_v[0] = create_zeros_2d(nc_, nc_)
  g_projection_p[0] = create_zeros_2d(nc_, nc_)
  g_projection_phi[0] = create_zeros_2d(nc_, nc_)
  projection_apply_velocity_bc(
    g_projection_u[0], g_projection_v[0], nc_, nc_,
  )
  let (initial_u_face, initial_v_face) = projection_build_face_velocities(
    g_projection_u[0], g_projection_v[0], nc_, nc_,
  )
  g_projection_u_face[0] = initial_u_face
  g_projection_v_face[0] = initial_v_face
  g_projection_last_diagnostics[0] = {
    pressure_iterations: 0,
    pressure_residual: 0.0,
    pressure_converged: true,
    divergence_norm: projection_face_divergence_l2(
      g_projection_u_face[0],
      g_projection_v_face[0],
      g_projection_dx[0],
      g_projection_dx[0],
      nc_,
      nc_,
    ),
    cfl_max: config.dt / g_projection_dx[0],
    finite: true,
  }
}

fn run_projection_tvd_steps(requested_steps : Int) -> Unit {
  let nc_ = g_projection_nc[0]
  let dx_ = g_projection_dx[0]
  for step_index = 0; step_index < requested_steps; step_index = step_index + 1 {
    let diagnostics = projection_tvd_step_with_faces(
      g_projection_u[0],
      g_projection_v[0],
      g_projection_p[0],
      g_projection_u_face[0],
      g_projection_v_face[0],
      g_projection_phi[0],
      1.0,
      g_projection_nu[0],
      dx_,
      dx_,
      nc_,
      nc_,
      g_projection_config[0],
    )
    g_projection_last_diagnostics[0] = diagnostics
    if !diagnostics.finite || !diagnostics.pressure_converged {
      return
    }
    g_projection_steps[0] = g_projection_steps[0] + 1
  }
}

fn get_projection_tvd_diagnostics() -> ProjectionTvdDiagnostics {
  g_projection_last_diagnostics[0]
}

fn collect_projection_tvd_export() -> SolverExport {
  {
    name: "projection-tvd",
    label: "Incremental Projection (Sign-Aware Upwind/TVD)",
    nc: g_projection_nc[0],
    h: g_projection_dx[0],
    x0: 0.0,
    y0: 0.0,
    steps: g_projection_steps[0],
    u: g_projection_u[0],
    v: g_projection_v[0],
    p: g_projection_p[0],
  }
}
