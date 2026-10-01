// ── Passive Scalar Transport (Step 6) ────────────────────────────────────────
//
// Solves the convection-diffusion equation for a passive scalar φ (e.g. temperature
// or concentration) one-way coupled to the Chorin velocity field (g_u, g_v):
//
//   ∂φ/∂t + u ∂φ/∂x + v ∂φ/∂y = α (∂²φ/∂x² + ∂²φ/∂y²)
//
// Discretisation:
//   - Convection: first-order upwind (sign of u/v selects donor cell)
//   - Diffusion:  central difference (second-order)
//   - Time:       explicit Euler (same Δt as Chorin solver)
//
// Boundary conditions (lid-driven thermal cavity analogy):
//   - Top lid   (i = ny-1): φ = 1.0  (hot driving surface)
//   - Bottom    (i = 0):    φ = 0.0  (cold wall)
//   - Left/Right (j = 0/nx-1): φ = 0.0 (cold walls)
//
// One-way coupling: velocity drives scalar, scalar does not affect velocity.
// This allows independent initialisation and step-counting.
//
// Péclet number: Pe = U·L / α = 1·2 / α  (U=1 lid velocity, L=2 domain size)
// Default α = 0.01 → Pe = 200  (convection-dominated, thin thermal boundary layers)
//
// Reference: Patankar (1980) §5.2 (upwind scheme), §3.3 (scalar transport).

// Runtime-settable diffusivity (single-element array, default Pe≈200)
let g_scalar_alpha : Array[Double] = [0.01]

// Global scalar field and step counter
let g_phi : Array[Array[Double]] = {
  let a = Array::make(ny, Array::make(nx, 0.0))
  for i = 0; i < ny; i = i + 1 { a[i] = Array::make(nx, 0.0) }
  a
}

let g_scalar_steps : Array[Int] = [0]

// ── Boundary conditions for φ ─────────────────────────────────────────────────
fn apply_scalar_bcs(phi : Array[Array[Double]]) -> Unit {
  // Bottom and side walls: φ = 0 (cold)
  for j = 0; j < nx; j = j + 1 { phi[0][j] = 0.0 }
  for i = 0; i < ny; i = i + 1 {
    phi[i][0]      = 0.0
    phi[i][nx - 1] = 0.0
  }
  // Top lid: φ = 1 (hot) — set last so corners belong to lid
  for j = 0; j < nx; j = j + 1 { phi[ny - 1][j] = 1.0 }
}

// ── One explicit-Euler time step for the scalar field ─────────────────────────
//
// u_vel, v_vel: advecting velocity arrays (typically g_u from Chorin solver)
// alpha_:       scalar diffusivity
fn scalar_one_step(
  phi   : Array[Array[Double]],
  u_vel : Array[Array[Double]],
  v_vel : Array[Array[Double]],
  alpha_ : Double,
) -> Unit {
  let phi_new = create_zeros_2d(ny, nx)
  let inv_dx  = 1.0 / dx
  let inv_dy  = 1.0 / dy
  let inv_dx2 = 1.0 / (dx * dx)
  let inv_dy2 = 1.0 / (dy * dy)

  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      let u_ij = u_vel[i][j]
      let v_ij = v_vel[i][j]

      // Upwind convection in x: donor-cell selection
      let conv_x = if u_ij >= 0.0 {
        u_ij * (phi[i][j] - phi[i][j - 1]) * inv_dx
      } else {
        u_ij * (phi[i][j + 1] - phi[i][j]) * inv_dx
      }

      // Upwind convection in y
      let conv_y = if v_ij >= 0.0 {
        v_ij * (phi[i][j] - phi[i - 1][j]) * inv_dy
      } else {
        v_ij * (phi[i + 1][j] - phi[i][j]) * inv_dy
      }

      // Central-difference diffusion
      let diff_x = alpha_ * (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]) * inv_dx2
      let diff_y = alpha_ * (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j]) * inv_dy2

      phi_new[i][j] = phi[i][j] + dt * (diff_x + diff_y - conv_x - conv_y)
    }
  }

  // Copy updated interior back
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      phi[i][j] = phi_new[i][j]
    }
  }

  apply_scalar_bcs(phi)
}

// ── WASM API ──────────────────────────────────────────────────────────────────

pub fn init_scalar() -> Unit {
  for i = 0; i < ny; i = i + 1 {
    for j = 0; j < nx; j = j + 1 {
      g_phi[i][j] = 0.0
    }
  }
  apply_scalar_bcs(g_phi)
  g_scalar_steps[0] = 0
}

pub fn set_scalar_alpha(alpha_ : Double) -> Unit {
  g_scalar_alpha[0] = alpha_
}

pub fn get_scalar_alpha() -> Double { g_scalar_alpha[0] }

// Advance scalar n steps using the current Chorin velocity field (g_u, g_v).
pub fn run_scalar_n_steps(n : Int) -> Unit {
  let alpha_ = g_scalar_alpha[0]
  for _k = 0; _k < n; _k = _k + 1 {
    scalar_one_step(g_phi, g_u[0], g_v[0], alpha_)
  }
  g_scalar_steps[0] = g_scalar_steps[0] + n
}

pub fn get_scalar_step_count() -> Int { g_scalar_steps[0] }

pub fn get_scalar_at(i : Int, j : Int) -> Double {
  if i >= 0 && i < ny && j >= 0 && j < nx { g_phi[i][j] } else { 0.0 }
}

// Mean scalar value over the interior domain (area average).
pub fn get_scalar_mean() -> Double {
  let mut s = 0.0
  let cnt = (ny - 2) * (nx - 2)
  for i = 1; i < ny - 1; i = i + 1 {
    for j = 1; j < nx - 1; j = j + 1 {
      s = s + g_phi[i][j]
    }
  }
  if cnt > 0 { s / cnt.to_double() } else { 0.0 }
}

// ── Pure helper: scalar one-step on supplied arrays (for unit testing) ─────────
//
// Identical to scalar_one_step but operates on arbitrary nr×nc arrays.
// Used by T125 (pure diffusion verification) without touching global state.
pub fn scalar_step_pure(
  phi   : Array[Array[Double]],
  u_vel : Array[Array[Double]],
  v_vel : Array[Array[Double]],
  alpha_ : Double,
  nr    : Int,
  nc    : Int,
  dx_   : Double,
  dy_   : Double,
  dt_   : Double,
) -> Unit {
  let phi_new = create_zeros_2d(nr, nc)
  let inv_dx  = 1.0 / dx_
  let inv_dy  = 1.0 / dy_
  let inv_dx2 = 1.0 / (dx_ * dx_)
  let inv_dy2 = 1.0 / (dy_ * dy_)

  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let u_ij = u_vel[i][j]
      let v_ij = v_vel[i][j]
      let conv_x = if u_ij >= 0.0 {
        u_ij * (phi[i][j] - phi[i][j - 1]) * inv_dx
      } else {
        u_ij * (phi[i][j + 1] - phi[i][j]) * inv_dx
      }
      let conv_y = if v_ij >= 0.0 {
        v_ij * (phi[i][j] - phi[i - 1][j]) * inv_dy
      } else {
        v_ij * (phi[i + 1][j] - phi[i][j]) * inv_dy
      }
      let diff_x = alpha_ * (phi[i][j + 1] - 2.0 * phi[i][j] + phi[i][j - 1]) * inv_dx2
      let diff_y = alpha_ * (phi[i + 1][j] - 2.0 * phi[i][j] + phi[i - 1][j]) * inv_dy2
      phi_new[i][j] = phi[i][j] + dt_ * (diff_x + diff_y - conv_x - conv_y)
    }
  }
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      phi[i][j] = phi_new[i][j]
    }
  }
}
