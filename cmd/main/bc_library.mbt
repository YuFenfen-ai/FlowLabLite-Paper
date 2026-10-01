// ── Boundary Condition Library (Step 8) ──────────────────────────────────────
//
// Pure BC functions (operate on arbitrary nr×nc arrays) and a periodic channel
// flow solver driven by a constant body force (Poiseuille validation problem).
//
// ── BC Types implemented ──────────────────────────────────────────────────────
//
// 1. Neumann (zero-gradient, ∂φ/∂n = 0):
//    First-order one-sided difference: boundary = first interior node.
//    Used for outlets, insulated walls, symmetry planes.
//
// 2. Periodic (left-right wrap):
//    φ[i][0] = φ[i][nc-2]   (left ghost ← right interior)
//    φ[i][nc-1] = φ[i][1]   (right ghost ← left interior)
//    Enables channel flow without pressure-driven inlet/outlet.
//
// 3. Symmetric (left or right plane):
//    Normal velocity = 0 (impermeable), tangential velocity Neumann.
//    e.g. symmetry at j=0: u[i][0] = u[i][1], v[i][0] = 0.
//
// ── Channel flow solver ───────────────────────────────────────────────────────
//
// Plane Poiseuille flow in a 2×2 domain driven by body force f_x (m/s²):
//   ∂u/∂t + u·∇u = ν∇²u + f_x     (x-momentum, body force replaces −∂p/∂x)
//   ∂v/∂t + u·∇v = ν∇²v            (y-momentum)
//   ∇·u = 0
//
// Boundary conditions:
//   Top/Bottom walls: u = v = 0  (no-slip)
//   Left/Right:       periodic   (wrap-around in x)
//   Pressure:         Neumann left/right, Neumann bottom, p[ny-1][j]=0 (reference)
//
// Analytical steady-state (Poiseuille):
//   u(y) = f_x / (2ν) · y · (H − y),   H = 2
//   u_max = f_x · H² / (8ν) = f_x · 4 / (8ν)  at y = H/2 = 1
//
// With f_x = 0.1, ν = 0.1: u_max = 0.5,  u_mean = u_max · 2/3 ≈ 0.333
//
// Reference: Batchelor, G.K. (1967). Introduction to Fluid Dynamics. §4.2.

// ── Pure Neumann BC functions ─────────────────────────────────────────────────

// ∂φ/∂y = 0 at bottom boundary (i = 0)
pub fn neumann_bottom(phi : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  for j = 0; j < nc; j = j + 1 { phi[0][j] = phi[1][j] }
  let _ = nr
}

// ∂φ/∂y = 0 at top boundary (i = nr-1)
pub fn neumann_top(phi : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  for j = 0; j < nc; j = j + 1 { phi[nr - 1][j] = phi[nr - 2][j] }
}

// ∂φ/∂x = 0 at left boundary (j = 0)
pub fn neumann_left(phi : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  for i = 0; i < nr; i = i + 1 { phi[i][0] = phi[i][1] }
  let _ = nc
}

// ∂φ/∂x = 0 at right boundary (j = nc-1)
pub fn neumann_right(phi : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  for i = 0; i < nr; i = i + 1 { phi[i][nc - 1] = phi[i][nc - 2] }
}

// ── Pure periodic BC functions ────────────────────────────────────────────────

// Periodic in x (left-right wrap-around).
// Node 0 mirrors node nc-2; node nc-1 mirrors node 1.
pub fn periodic_lr(phi : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  for i = 0; i < nr; i = i + 1 {
    phi[i][0]      = phi[i][nc - 2]   // left ghost ← right interior
    phi[i][nc - 1] = phi[i][1]        // right ghost ← left interior
  }
}

// Periodic in y (top-bottom wrap-around).
pub fn periodic_tb(phi : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  let _ = nc
  for j = 0; j < nc; j = j + 1 {
    phi[0][j]      = phi[nr - 2][j]
    phi[nr - 1][j] = phi[1][j]
  }
}

// ── Symmetric BC functions ────────────────────────────────────────────────────

// Symmetry plane at left wall (j = 0):
//   Normal component v[i][0] = 0 (no flow through symmetry plane)
//   Tangential component u: Neumann (∂u/∂x = 0)
pub fn symmetric_left(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Unit {
  for i = 0; i < nr; i = i + 1 {
    u[i][0] = u[i][1]   // Neumann for tangential
    v[i][0] = 0.0       // zero normal
  }
  let _ = nc
}

// Symmetry plane at right wall (j = nc-1)
pub fn symmetric_right(
  u : Array[Array[Double]],
  v : Array[Array[Double]],
  nr : Int,
  nc : Int,
) -> Unit {
  for i = 0; i < nr; i = i + 1 {
    u[i][nc - 1] = u[i][nc - 2]
    v[i][nc - 1] = 0.0
  }
  let _ = nr
}

// ── Channel flow global state ─────────────────────────────────────────────────

let g_u_ch : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_v_ch : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_p_ch : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_ch_steps : Array[Int]    = [0]
let g_ch_fx    : Array[Double] = [0.1]   // body force (m/s²)
let g_ch_nu    : Array[Double] = [0.1]   // kinematic viscosity

// ── One channel time step ─────────────────────────────────────────────────────
//
// Chorin projection + body force f_x in u-momentum + periodic x-BCs.
fn channel_one_step(
  u   : Array[Array[Double]],
  v   : Array[Array[Double]],
  p   : Array[Array[Double]],
  f_x : Double,
  nu_ : Double,
) -> Unit {
  let u_star = create_zeros_2d(ny, nx)
  let v_star = create_zeros_2d(ny, nx)
  let b_corr = create_zeros_2d(ny, nx)
  let p_corr = create_zeros_2d(ny, nx)
  let coef_xx = nu_ * dt / (dx * dx)
  let coef_yy = nu_ * dt / (dy * dy)
  let inv_dx  = 1.0 / dx
  let inv_dy  = 1.0 / dy

  // ── u-momentum predictor (+ body force) ──────────────────────────────────
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
      u_star[i][j] = uij + visc_x + visc_y - conv_x - conv_y + dt * f_x
    }
  }

  // ── v-momentum predictor (no body force in y) ─────────────────────────────
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
      v_star[i][j] = vij + visc_x + visc_y - conv_x - conv_y
    }
  }

  // Top/bottom no-slip BCs for intermediate velocity
  for j = 0; j < nx; j = j + 1 {
    u_star[0][j]      = 0.0;  v_star[0][j]      = 0.0
    u_star[ny - 1][j] = 0.0;  v_star[ny - 1][j] = 0.0
  }
  // Periodic x-BCs for intermediate velocity
  periodic_lr(u_star, ny, nx)
  periodic_lr(v_star, ny, nx)

  // ── Pressure correction (Neumann all sides except top reference) ──────────
  build_up_b_array(b_corr, rho, dt, u_star, v_star, dx, dy, nx, ny) |> ignore
  pressure_poisson_array(p_corr, dx, dy, b_corr, nx, ny) |> ignore

  // ── Velocity and pressure update ──────────────────────────────────────────
  let corr_u = dt / (2.0 * rho * dx)
  let corr_v = dt / (2.0 * rho * dy)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      u[i][j] = u_star[i][j] - corr_u * (p_corr[i][j + 1] - p_corr[i][j - 1])
      v[i][j] = v_star[i][j] - corr_v * (p_corr[i + 1][j] - p_corr[i - 1][j])
      p[i][j] = p[i][j] + p_corr[i][j]
    }
  }

  // Pressure BCs (Neumann left/right/bottom, p=0 at top as reference)
  for i = 0; i < ny; i = i + 1 {
    p[i][0]      = p[i][1]
    p[i][nx - 1] = p[i][nx - 2]
  }
  for j = 0; j < nx; j = j + 1 {
    p[0][j]      = p[1][j]
    p[ny - 1][j] = 0.0
  }

  // Velocity BCs: no-slip top/bottom, periodic x
  for j = 0; j < nx; j = j + 1 {
    u[0][j]      = 0.0;  v[0][j]      = 0.0
    u[ny - 1][j] = 0.0;  v[ny - 1][j] = 0.0
  }
  periodic_lr(u, ny, nx)
  periodic_lr(v, ny, nx)
}

// ── WASM API – Channel solver ─────────────────────────────────────────────────

pub fn init_channel() -> Unit {
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      g_u_ch[i][j] = 0.0
      g_v_ch[i][j] = 0.0
      g_p_ch[i][j] = 0.0
    }
  }
  g_ch_steps[0] = 0
}

pub fn set_channel_fx(f_x : Double) -> Unit { g_ch_fx[0] = f_x }
pub fn set_channel_nu(nu_ : Double) -> Unit { g_ch_nu[0] = nu_ }
pub fn get_channel_fx() -> Double { g_ch_fx[0] }

pub fn run_channel_n_steps(n : Int) -> Unit {
  let f_x = g_ch_fx[0]
  let nu_ = g_ch_nu[0]
  for _k = 0; _k < n; _k = _k + 1 {
    channel_one_step(g_u_ch, g_v_ch, g_p_ch, f_x, nu_)
  }
  g_ch_steps[0] = g_ch_steps[0] + n
}

pub fn get_channel_step_count() -> Int { g_ch_steps[0] }

pub fn get_u_channel_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_u_ch[i][j] } else { 0.0 }
}

pub fn get_v_channel_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_v_ch[i][j] } else { 0.0 }
}

pub fn get_p_channel_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_p_ch[i][j] } else { 0.0 }
}

// Analytical Poiseuille u at row i: u_analytic = f_x/(2ν) · y_i · (H - y_i)
// y_i = i * dy,  H = (ny-1) * dy = 2.0
pub fn get_channel_analytic_u(i : Int) -> Double {
  let f_x = g_ch_fx[0]
  let nu_ = g_ch_nu[0]
  let yi  = i.to_double() * dy
  let cap_h  = (ny - 1).to_double() * dy
  f_x / (2.0 * nu_) * yi * (cap_h - yi)
}
