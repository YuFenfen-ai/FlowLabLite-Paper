// ── Energy Equation + Boussinesq Approximation (Step 7) ──────────────────────
//
// Solves the differentially heated square cavity (De Vahl Davis 1983 benchmark)
// using the Boussinesq approximation.  Two-way coupling:
//   velocity → scalar  (convection of φ)
//   scalar   → velocity (buoyancy body force on v)
//
// Governing equations (dimensional, Boussinesq):
//   ∂u/∂t + u·∇u = -∇p/ρ + ν·∇²u
//   ∂v/∂t + u·∇v = -∇p/ρ + ν·∇²v + β_g·φ     ← buoyancy (+y = upward)
//   ∂φ/∂t + u·∇φ = α·∇²φ                        ← scalar (temperature)
//   ∇·u = 0
//
// Time scheme: Chorin fractional-step (explicit Euler + pressure projection).
// Convection: upwind; Diffusion: central difference.
//
// Boundary conditions:
//   Left  wall (j=0):    φ = 1.0  (hot),  u = v = 0
//   Right wall (j=nx-1): φ = 0.0  (cold), u = v = 0
//   Top   wall (i=ny-1): ∂φ/∂n = 0  (insulated, Neumann), u = v = 0
//   Bottom wall (i=0):   ∂φ/∂n = 0  (insulated, Neumann), u = v = 0
//   (No lid — all walls no-slip)
//
// Buoyancy convention: f_v = β_g · φ  (reference temperature = 0, cold wall).
//   Hot regions (φ > 0) experience upward buoyancy (+y direction).
//   Drives clockwise circulation: rising near hot left wall, sinking near cold right wall.
//
// Default parameters → Ra = β_g · ΔT · L³ / (ν · α) = 0.1 · 1 · 8 / (10⁻⁴) ≈ 8000
//   L = 2 (domain), ΔT = 1, ν = α = 0.01, β_g = 0.1
//
// Reference: De Vahl Davis, G. (1983). Natural convection of air in a square cavity:
//            a benchmark numerical solution. Int. J. Numer. Meth. Fluids 3(3), 249-264.

let g_bous_nu    : Array[Double] = [0.01]   // kinematic viscosity  (= Pr · α in non-dim)
let g_bous_alpha : Array[Double] = [0.01]   // thermal diffusivity
let g_bous_beta  : Array[Double] = [0.1]    // effective buoyancy: g · β_thermal

// ── Global state ──────────────────────────────────────────────────────────────

let g_u_bous : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_v_bous : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_p_bous : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_phi_bous : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_bous_steps : Array[Int] = [0]

// ── Boundary conditions ───────────────────────────────────────────────────────

fn apply_bous_velocity_bcs(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
) -> Unit {
  // All four walls: no-slip (no lid unlike Chorin)
  for j = 0; j < nx; j = j + 1 {
    u[0][j]      = 0.0;  v[0][j]      = 0.0
    u[ny - 1][j] = 0.0;  v[ny - 1][j] = 0.0
  }
  for i = 0; i < ny; i = i + 1 {
    u[i][0]      = 0.0;  v[i][0]      = 0.0
    u[i][nx - 1] = 0.0;  v[i][nx - 1] = 0.0
  }
}

// Apply scalar BCs: Neumann at top/bottom first, Dirichlet at left/right last
// (so corners belong to the hot/cold side walls, not the insulated top/bottom).
fn apply_bous_scalar_bcs(phi : Array[Array[Double]]) -> Unit {
  // Top/bottom: zero-gradient (insulated)
  for j = 1; j < nx - 1; j = j + 1 {
    phi[0][j]      = phi[1][j]
    phi[ny - 1][j] = phi[ny - 2][j]
  }
  // Left/right: Dirichlet (overwrites corners → hot/cold wins)
  for i = 0; i < ny; i = i + 1 {
    phi[i][0]      = 1.0   // hot wall
    phi[i][nx - 1] = 0.0   // cold wall
  }
}

// ── One coupled Boussinesq step ───────────────────────────────────────────────
//
// Sequence: u* predictor → v* predictor (+ buoyancy) → BCs → pressure solve →
//           velocity correction → pressure BCs → velocity BCs → scalar step → scalar BCs.
fn bous_one_step(
  u      : Array[Array[Double]],
  v      : Array[Array[Double]],
  p      : Array[Array[Double]],
  phi    : Array[Array[Double]],
  nu_    : Double,
  alpha_ : Double,
  beta_  : Double,
) -> Unit {
  let u_star = create_zeros_2d(ny, nx)
  let v_star = create_zeros_2d(ny, nx)
  let b_corr = create_zeros_2d(ny, nx)
  let p_corr = create_zeros_2d(ny, nx)

  let inv_dx  = 1.0 / dx
  let inv_dy  = 1.0 / dy
  let coef_xx = nu_ * dt / (dx * dx)
  let coef_yy = nu_ * dt / (dy * dy)

  // ── u-momentum predictor ─────────────────────────────────────────────────
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let uij = u[i][j]
      let vij = v[i][j]
      let conv_x = if uij >= 0.0 {
        uij * dt * inv_dx * (u[i][j] - u[i][j - 1])
      } else {
        uij * dt * inv_dx * (u[i][j + 1] - u[i][j])
      }
      let conv_y = if vij >= 0.0 {
        vij * dt * inv_dy * (u[i][j] - u[i - 1][j])
      } else {
        vij * dt * inv_dy * (u[i + 1][j] - u[i][j])
      }
      let visc_x = coef_xx * (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1])
      let visc_y = coef_yy * (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j])
      u_star[i][j] = uij + visc_x + visc_y - conv_x - conv_y
    }
  }

  // ── v-momentum predictor (+ buoyancy source) ─────────────────────────────
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let uij = u[i][j]
      let vij = v[i][j]
      let conv_x = if uij >= 0.0 {
        uij * dt * inv_dx * (v[i][j] - v[i][j - 1])
      } else {
        uij * dt * inv_dx * (v[i][j + 1] - v[i][j])
      }
      let conv_y = if vij >= 0.0 {
        vij * dt * inv_dy * (v[i][j] - v[i - 1][j])
      } else {
        vij * dt * inv_dy * (v[i + 1][j] - v[i][j])
      }
      let visc_x = coef_xx * (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1])
      let visc_y = coef_yy * (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j])
      // Boussinesq buoyancy: β_g · φ  (reference temperature = 0)
      let buoyancy = beta_ * phi[i][j]
      v_star[i][j] = vij + visc_x + visc_y - conv_x - conv_y + dt * buoyancy
    }
  }

  // BCs for intermediate velocities
  apply_bous_velocity_bcs(u_star, v_star)

  // ── Pressure correction ───────────────────────────────────────────────────
  build_up_b_array(b_corr, rho, dt, u_star, v_star, dx, dy, nx, ny) |> ignore
  pressure_poisson_array(p_corr, dx, dy, b_corr, nx, ny) |> ignore

  // ── Velocity and pressure update ─────────────────────────────────────────
  let corr_u = dt / (2.0 * rho * dx)
  let corr_v = dt / (2.0 * rho * dy)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      u[i][j] = u_star[i][j] - corr_u * (p_corr[i][j + 1] - p_corr[i][j - 1])
      v[i][j] = v_star[i][j] - corr_v * (p_corr[i + 1][j] - p_corr[i - 1][j])
      p[i][j] = p[i][j] + p_corr[i][j]
    }
  }

  // Pressure BCs (same reference as Chorin: p[ny-1][j]=0)
  for i = 0; i < ny; i = i + 1 {
    p[i][0]      = p[i][1]
    p[i][nx - 1] = p[i][nx - 2]
  }
  for j = 0; j < nx; j = j + 1 {
    p[0][j]      = p[1][j]
    p[ny - 1][j] = 0.0
  }

  // Velocity BCs (all no-slip)
  apply_bous_velocity_bcs(u, v)

  // ── Scalar transport step ─────────────────────────────────────────────────
  let phi_new = create_zeros_2d(ny, nx)
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let uij = u[i][j]
      let vij = v[i][j]
      let sc_cx = if uij >= 0.0 {
        uij * (phi[i][j] - phi[i][j - 1]) * inv_dx
      } else {
        uij * (phi[i][j + 1] - phi[i][j]) * inv_dx
      }
      let sc_cy = if vij >= 0.0 {
        vij * (phi[i][j] - phi[i - 1][j]) * inv_dy
      } else {
        vij * (phi[i + 1][j] - phi[i][j]) * inv_dy
      }
      let sc_dx = alpha_ * (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]) * inv_dx2
      let sc_dy = alpha_ * (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j]) * inv_dy2
      phi_new[i][j] = phi[i][j] + dt * (sc_dx + sc_dy - sc_cx - sc_cy)
    }
  }
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      phi[i][j] = phi_new[i][j]
    }
  }

  apply_bous_scalar_bcs(phi)
}

// ── WASM API ──────────────────────────────────────────────────────────────────

pub fn init_boussinesq() -> Unit {
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      g_u_bous[i][j]   = 0.0
      g_v_bous[i][j]   = 0.0
      g_p_bous[i][j]   = 0.0
      g_phi_bous[i][j] = 0.0
    }
  }
  // Set initial scalar BCs so first step has correct Dirichlet sources
  apply_bous_scalar_bcs(g_phi_bous)
  g_bous_steps[0] = 0
}

pub fn set_bous_nu(val : Double)    -> Unit { g_bous_nu[0]    = val }
pub fn set_bous_alpha(val : Double) -> Unit { g_bous_alpha[0] = val }
pub fn set_bous_beta(val : Double)  -> Unit { g_bous_beta[0]  = val }

pub fn get_bous_nu()    -> Double { g_bous_nu[0]    }
pub fn get_bous_alpha() -> Double { g_bous_alpha[0] }
pub fn get_bous_beta()  -> Double { g_bous_beta[0]  }

pub fn run_boussinesq_n_steps(n : Int) -> Unit {
  let nu_    = g_bous_nu[0]
  let alpha_ = g_bous_alpha[0]
  let beta_  = g_bous_beta[0]
  for _k = 0; _k < n; _k = _k + 1 {
    bous_one_step(g_u_bous, g_v_bous, g_p_bous, g_phi_bous, nu_, alpha_, beta_)
  }
  g_bous_steps[0] = g_bous_steps[0] + n
}

pub fn get_bous_step_count() -> Int { g_bous_steps[0] }

pub fn get_u_bous_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_u_bous[i][j] } else { 0.0 }
}

pub fn get_v_bous_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_v_bous[i][j] } else { 0.0 }
}

pub fn get_p_bous_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_p_bous[i][j] } else { 0.0 }
}

pub fn get_phi_bous_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_phi_bous[i][j] } else { 0.0 }
}

// Nusselt number at the hot left wall (j=0): Nu = -L · ∂φ/∂x|_{x=0} / ΔT
// For unit ΔT and L=2: Nu = -2 · (φ[i][1]-φ[i][0])/dx  averaged over interior rows.
pub fn get_bous_nusselt() -> Double {
  let mut sum = 0.0
  let cnt = ny - 2
  for i = 1; i < ny - 1; i = i + 1 {
    // One-sided difference: ∂φ/∂x ≈ (φ[i][1] - φ[i][0]) / dx  (both > 0 since φ[i][0]=1)
    let dphi_dx = (g_phi_bous[i][1] - g_phi_bous[i][0]) / dx
    // Heat flux: -α · ∂φ/∂x  → Nusselt = -(L/ΔT) · mean(∂φ/∂x)
    sum = sum + dphi_dx
  }
  if cnt > 0 { -2.0 * sum / cnt.to_double() } else { 0.0 }
}
