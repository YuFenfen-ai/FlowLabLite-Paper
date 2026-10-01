// ── FVM-SIMPLE Solver (Step 1) ────────────────────────────────────────────────
//
// Implements the SIMPLE algorithm using finite-volume face fluxes (Patankar 1980).
//
// Distinguishing features vs. the FDM SIMPLE in main.mbt (simple_one_iter):
//   1. Momentum predictor: implicit Gauss-Seidel using FVM coefficients
//      (upwind convection + central diffusion, face-flux based).
//      FDM SIMPLE uses an explicit Euler update.
//   2. Face mass fluxes explicitly computed:
//      F_e = ρ u_e Δy,  D_e = ρν Δy/Δx   (upwind scheme)
//   3. Pressure/velocity correction step reuses existing infrastructure
//      (build_up_b_array + pressure_poisson_array) for consistency.
//
// For a uniform Cartesian grid, FVM and FDM converge to the same steady state.
// The FVM formulation is more physically rigorous and generalises to
// non-uniform / unstructured grids.
//
// Reference: Patankar, S.V. (1980). Numerical Heat Transfer and Fluid Flow.
//            §5.3 (SIMPLE algorithm), §5.2 (upwind scheme coefficients).

let fvm_alpha_u : Double = 0.7   // velocity under-relaxation (same as FDM SIMPLE)
let fvm_alpha_p : Double = 0.3   // pressure under-relaxation

// One SIMPLE iteration using FVM momentum predictor.
fn fvm_one_iter(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  p : Array[Array[Double]],
) -> Unit {
  let nu_  = g_nu_runtime[0]
  // Diffusion coefficient for one face: D = ρν(Δy/Δx) = ρν  (Δx = Δy)
  let diff = rho * nu_

  let u_star = create_zeros_2d(ny, nx)
  let v_star = create_zeros_2d(ny, nx)
  let b_corr = create_zeros_2d(ny, nx)
  let p_corr = create_zeros_2d(ny, nx)

  // ── (a) u-momentum predictor ────────────────────────────────────────────
  // One Gauss-Seidel sweep of: a_P u*_P = Σ a_nb u_nb + b_P
  // a_E = D_e + max(-F_e, 0),  etc.  (upwind scheme)
  // b_P = face-interpolated pressure source: −Δy(p_e − p_w)/2
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      // Face-centre velocities (linear interpolation)
      let ue  = (u[i][j] + u[i][j + 1]) / 2.0
      let uw  = (u[i][j] + u[i][j - 1]) / 2.0
      let vn  = (v[i][j] + v[i + 1][j]) / 2.0
      let vs_ = (v[i][j] + v[i - 1][j]) / 2.0
      // Face mass fluxes: F = ρ u_face A_face (A = Δy for E/W, Δx for N/S)
      let fe  = rho * ue  * dy
      let fw  = rho * uw  * dy
      let fn_ = rho * vn  * dx
      let fs_ = rho * vs_ * dx
      // Neighbour coefficients: a_nb = D + max(±F, 0)
      let ae  = diff + (if fe  < 0.0 { -fe  } else { 0.0 })
      let aw  = diff + (if fw  > 0.0 {  fw  } else { 0.0 })
      let an  = diff + (if fn_ < 0.0 { -fn_ } else { 0.0 })
      let as_ = diff + (if fs_ > 0.0 {  fs_ } else { 0.0 })
      // Central coefficient (neglect continuity residual for stability)
      let ap  = ae + aw + an + as_
      // Pressure-gradient source at node P (central-difference face pressures)
      let bp  = -0.5 * dy * (p[i][j + 1] - p[i][j - 1])
      // Gauss-Seidel update (reads from u, writes to u_star)
      let rhs   = ae * u[i][j + 1] + aw * u[i][j - 1] + an * u[i + 1][j] + as_ * u[i - 1][j] + bp
      let ap_s  = if ap > 1.0e-14 { ap } else { 4.0 * diff }
      let u_new = rhs / ap_s
      u_star[i][j] = fvm_alpha_u * u_new + (1.0 - fvm_alpha_u) * u[i][j]
    }
  }

  // ── (b) v-momentum predictor ────────────────────────────────────────────
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let ue  = (u[i][j] + u[i][j + 1]) / 2.0
      let uw  = (u[i][j] + u[i][j - 1]) / 2.0
      let vn  = (v[i][j] + v[i + 1][j]) / 2.0
      let vs_ = (v[i][j] + v[i - 1][j]) / 2.0
      let fe  = rho * ue  * dy
      let fw  = rho * uw  * dy
      let fn_ = rho * vn  * dx
      let fs_ = rho * vs_ * dx
      let ae  = diff + (if fe  < 0.0 { -fe  } else { 0.0 })
      let aw  = diff + (if fw  > 0.0 {  fw  } else { 0.0 })
      let an  = diff + (if fn_ < 0.0 { -fn_ } else { 0.0 })
      let as_ = diff + (if fs_ > 0.0 {  fs_ } else { 0.0 })
      let ap  = ae + aw + an + as_
      let bp  = -0.5 * dx * (p[i + 1][j] - p[i - 1][j])
      let rhs   = ae * v[i][j + 1] + aw * v[i][j - 1] + an * v[i + 1][j] + as_ * v[i - 1][j] + bp
      let ap_s  = if ap > 1.0e-14 { ap } else { 4.0 * diff }
      let v_new = rhs / ap_s
      v_star[i][j] = fvm_alpha_u * v_new + (1.0 - fvm_alpha_u) * v[i][j]
    }
  }

  // BCs for predicted velocities — walls first, lid last
  for j = 0; j < nx; j = j + 1 {
    u_star[0][j]      = 0.0
    v_star[0][j]      = 0.0
    v_star[ny - 1][j] = 0.0
  }
  for i = 0; i < ny; i = i + 1 {
    u_star[i][0]      = 0.0
    u_star[i][nx - 1] = 0.0
    v_star[i][0]      = 0.0
    v_star[i][nx - 1] = 0.0
  }
  for j = 0; j < nx; j = j + 1 { u_star[ny - 1][j] = 1.0 }

  // ── (c) Pressure-correction source ──────────────────────────────────────
  // b = −(ρ/Δt) ∇·u*  (reuses Chorin/FDM infrastructure, same Δt constant)
  build_up_b_array(b_corr, rho, dt, u_star, v_star, dx, dy, nx, ny) |> ignore

  // ── (d) Pressure-correction Gauss-Seidel ────────────────────────────────
  pressure_poisson_array(p_corr, dx, dy, b_corr, nx, ny) |> ignore

  // ── (e) Velocity and pressure update ────────────────────────────────────
  // Velocity correction consistent with build_up_b_array scaling:
  //   u = u* − (Δt / 2ρΔx) ∇p'
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      u[i][j] = u_star[i][j] - dt / (2.0 * rho * dx) * (p_corr[i][j + 1] - p_corr[i][j - 1])
      v[i][j] = v_star[i][j] - dt / (2.0 * rho * dy) * (p_corr[i + 1][j] - p_corr[i - 1][j])
      p[i][j] = p[i][j] + fvm_alpha_p * p_corr[i][j]
    }
  }

  // Pressure BCs
  for i = 0; i < ny; i = i + 1 {
    p[i][0]      = p[i][1]
    p[i][nx - 1] = p[i][nx - 2]
  }
  for j = 0; j < nx; j = j + 1 {
    p[0][j]      = p[1][j]
    p[ny - 1][j] = 0.0
  }

  // Velocity BCs — walls first, lid last
  for j = 0; j < nx; j = j + 1 {
    u[0][j]      = 0.0
    v[0][j]      = 0.0
    v[ny - 1][j] = 0.0
  }
  for i = 0; i < ny; i = i + 1 {
    u[i][0]      = 0.0
    u[i][nx - 1] = 0.0
    v[i][0]      = 0.0
    v[i][nx - 1] = 0.0
  }
  for j = 0; j < nx; j = j + 1 { u[ny - 1][j] = 1.0 }
}

// ── Global FVM state ──────────────────────────────────────────────────────────
let g_u_fvm : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_v_fvm : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_p_fvm : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_fvm_iters : Array[Int] = [0]

// ── WASM API – FVM solver ─────────────────────────────────────────────────────

pub fn init_fvm() -> Unit {
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      g_u_fvm[i][j] = 0.0
      g_v_fvm[i][j] = 0.0
      g_p_fvm[i][j] = 0.0
    }
  }
  g_fvm_iters[0] = 0
}

pub fn run_fvm_n_iter(n : Int) -> Unit {
  for _k = 0; _k < n; _k = _k + 1 {
    fvm_one_iter(g_u_fvm, g_v_fvm, g_p_fvm)
  }
  g_fvm_iters[0] = g_fvm_iters[0] + n
}

pub fn get_fvm_step_count() -> Int { g_fvm_iters[0] }

pub fn get_u_fvm_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_u_fvm[i][j] } else { 0.0 }
}

pub fn get_v_fvm_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_v_fvm[i][j] } else { 0.0 }
}

pub fn get_p_fvm_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_p_fvm[i][j] } else { 0.0 }
}

pub fn get_fvm_divergence_norm() -> Double {
  let mut sum = 0.0
  let cnt = (ny - 2) * (nx - 2)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let du = (g_u_fvm[i][j + 1] - g_u_fvm[i][j - 1]) / (2.0 * dx)
      let dv = (g_v_fvm[i + 1][j] - g_v_fvm[i - 1][j]) / (2.0 * dy)
      let d  = du + dv
      sum    = sum + d * d
    }
  }
  if cnt > 0 { (sum / cnt.to_double()).sqrt() } else { 0.0 }
}
