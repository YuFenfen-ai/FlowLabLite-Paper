// FlowLabLite — Convection Schemes
//
// Provides alternative convection discretisations for the momentum equations.
// Currently implemented:
//   - First-order upwind  (baseline, used by default in main.mbt)
//   - QUICK               (Quadratic Upstream Interpolation for Convective
//                          Kinematics, Leonard 1979)
//
// Reference: Leonard, B.P. (1979). "A stable and accurate convective
//   modelling procedure based on quadratic upstream interpolation."
//   Comput. Methods Appl. Mech. Eng. 19, 59–98.
//
// Interface contract
// ------------------
// Both functions return a 2D array of size [n_rows][n_cols] containing the
// convective contribution  u * (∂φ/∂x) * dt  or  v * (∂φ/∂y) * dt  at every
// node (i, j).  The caller subtracts this array from the velocity field, exactly
// as the existing conv_u_x / conv_u_y terms are used in main.mbt.
//
// Interior nodes (i = 1..nr-2, j = 1..nc-2) use QUICK interpolation.
// Boundary rows/columns are clamped to 0.0 (BCs are applied elsewhere).
// Nodes adjacent to boundaries that lack the full QUICK stencil fall back to
// first-order upwind to avoid out-of-bounds access.

// ─────────────────────────────────────────────────────────────────────────────
// QUICK helper — interpolate φ at east face (j + ½) using upwind quadratic
// ─────────────────────────────────────────────────────────────────────────────

// East face (j + ½) for a positive-u flow (upstream = j-1, j, j+1):
//   φ_e = (6 φ_j + 3 φ_{j+1} - φ_{j-1}) / 8
// East face (j + ½) for a negative-u flow (upstream = j+2, j+1, j):
//   φ_e = (6 φ_{j+1} + 3 φ_j - φ_{j+2}) / 8
fn quick_east(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nc   : Int,
  u_ij : Double,
) -> Double {
  if u_ij >= 0.0 {
    // Positive flow: upstream stencil needs j-1
    if j == 0 {
      return phi[i][j]   // boundary: degenerate to node value
    }
    (6.0 * phi[i][j] + 3.0 * phi[i][j + 1] - phi[i][j - 1]) / 8.0
  } else {
    // Negative flow: upstream stencil needs j+2
    if j + 2 >= nc {
      return phi[i][j + 1]   // boundary: degenerate to next node value
    }
    (6.0 * phi[i][j + 1] + 3.0 * phi[i][j] - phi[i][j + 2]) / 8.0
  }
}

// West face (j − ½):
//   Positive u: upstream stencil is j-2, j-1, j  → (6 φ_{j-1} + 3 φ_j - φ_{j-2}) / 8
//   Negative u: upstream stencil is j+1, j, j-1  → (6 φ_j + 3 φ_{j-1} - φ_{j+1}) / 8
fn quick_west(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nc   : Int,
  u_ij : Double,
) -> Double {
  if u_ij >= 0.0 {
    if j - 2 < 0 {
      return phi[i][j - 1]   // boundary: degenerate
    }
    (6.0 * phi[i][j - 1] + 3.0 * phi[i][j] - phi[i][j - 2]) / 8.0
  } else {
    if j + 1 >= nc {
      return phi[i][j]   // boundary: degenerate
    }
    (6.0 * phi[i][j] + 3.0 * phi[i][j - 1] - phi[i][j + 1]) / 8.0
  }
}

// North face (i + ½) for y-direction:
//   Positive v: upstream stencil i-1, i, i+1 → (6 φ_i + 3 φ_{i+1} - φ_{i-1}) / 8
//   Negative v: upstream stencil i+2, i+1, i → (6 φ_{i+1} + 3 φ_i - φ_{i+2}) / 8
fn quick_north(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nr   : Int,
  v_ij : Double,
) -> Double {
  if v_ij >= 0.0 {
    if i == 0 {
      return phi[i][j]
    }
    (6.0 * phi[i][j] + 3.0 * phi[i + 1][j] - phi[i - 1][j]) / 8.0
  } else {
    if i + 2 >= nr {
      return phi[i + 1][j]
    }
    (6.0 * phi[i + 1][j] + 3.0 * phi[i][j] - phi[i + 2][j]) / 8.0
  }
}

// South face (i − ½):
fn quick_south(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nr   : Int,
  v_ij : Double,
) -> Double {
  if v_ij >= 0.0 {
    if i - 2 < 0 {
      return phi[i - 1][j]
    }
    (6.0 * phi[i - 1][j] + 3.0 * phi[i][j] - phi[i - 2][j]) / 8.0
  } else {
    if i + 1 >= nr {
      return phi[i][j]
    }
    (6.0 * phi[i][j] + 3.0 * phi[i - 1][j] - phi[i + 1][j]) / 8.0
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Public API
// ─────────────────────────────────────────────────────────────────────────────

/// QUICK convective term in the x-direction: u * ∂φ/∂x * dt.
///
/// Returns a [nr][nc] array.  Interior nodes (1..nr-2, 1..nc-2) use the QUICK
/// quadratic upstream interpolation; nodes that lack a full 3-point stencil
/// (next to walls) fall back to first-order upwind automatically via the
/// quick_east/west helpers.  Boundary rows/columns return 0.0.
///
/// Parameters
///   phi   – field being advected (e.g. u-velocity)
///   u_vel – advecting velocity in x-direction (same layout as phi)
///   nr    – number of rows
///   nc    – number of columns
///   dx    – x-cell size
///   dt    – time-step
fn quick_convection_x(
  phi   : Array[Array[Double]],
  u_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dx    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let u_ij  = u_vel[i][j]
      let phi_e = quick_east(phi, i, j, nc, u_ij)
      let phi_w = quick_west(phi, i, j, nc, u_ij)
      result[i][j] = u_ij * (phi_e - phi_w) / dx * dt
    }
  }
  result
}

/// QUICK convective term in the y-direction: v * ∂φ/∂y * dt.
///
/// Same conventions as quick_convection_x but operates on rows (i-direction).
fn quick_convection_y(
  phi   : Array[Array[Double]],
  v_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dy    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let v_ij  = v_vel[i][j]
      let phi_n = quick_north(phi, i, j, nr, v_ij)
      let phi_s = quick_south(phi, i, j, nr, v_ij)
      result[i][j] = v_ij * (phi_n - phi_s) / dy * dt
    }
  }
  result
}

/// First-order upwind convective term in x: u * ∂φ/∂x * dt  (backward difference).
/// Provided for comparison with QUICK; used when --scheme upwind (default).
fn upwind_convection_x(
  phi   : Array[Array[Double]],
  u_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dx    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = u_vel[i][j] * dt / dx * (phi[i][j] - phi[i][j - 1])
    }
  }
  result
}

// ─────────────────────────────────────────────────────────────────────────────
// SOR / SSOR pressure Poisson solvers
// ─────────────────────────────────────────────────────────────────────────────
//
// These are standalone alternatives to the Gauss-Seidel inner-loop in
// pressure_poisson_array.  They accept generic nr/nc so they can be unit-tested
// on small grids without touching the global solver state.
//
// Boundary conditions follow the lid-driven cavity convention:
//   p[:, nc-1] = p[:, nc-2]   right wall (dp/dx = 0)
//   p[0, :]    = p[1, :]      bottom wall (dp/dy = 0)
//   p[:, 0]    = p[:, 1]      left wall  (dp/dx = 0)
//   p[nr-1, :] = 0.0          top lid (Dirichlet p=0)
//
// Reference: Ferziger & Perić, "Computational Methods for Fluid Dynamics",
//   §3.7 SOR; Golub & Van Loan, "Matrix Computations" §10.1 SSOR.

/// Apply cavity pressure boundary conditions in-place.
fn apply_pressure_bc_sor(p : Array[Array[Double]], nr : Int, nc : Int) -> Unit {
  // Right wall: dp/dx = 0
  for i = 0; i < nr; i = i + 1 { p[i][nc - 1] = p[i][nc - 2] }
  // Bottom wall: dp/dy = 0
  for j = 0; j < nc; j = j + 1 { p[0][j] = p[1][j] }
  // Left wall: dp/dx = 0
  for i = 0; i < nr; i = i + 1 { p[i][0] = p[i][1] }
  // Top lid: p = 0 (Dirichlet)
  for j = 0; j < nc; j = j + 1 { p[nr - 1][j] = 0.0 }
}

/// Successive Over-Relaxation (SOR) pressure Poisson solver.
///
/// Solves ∇²p = b in-place using `nit_` sweeps of the SOR iteration with
/// relaxation factor ω:
///   p_GS  = (dy²(p_e + p_w) + dx²(p_n + p_s) − dx²dy²·b) / (2(dx²+dy²))
///   p_new = (1−ω)·p_old + ω·p_GS
///
/// ω = 1.0  →  reduces exactly to Gauss-Seidel.
/// ω ∈ (1,2) →  over-relaxation; optimal ω ≈ 2/(1+π/n) for Poisson on uniform grid.
/// ω > 2    →  diverges.
///
/// Parameters:  p – pressure field (modified in-place)
///              dx, dy – cell sizes; b – RHS; omega – relaxation factor
///              nr, nc – grid dimensions; nit_ – number of sweeps
fn pressure_poisson_sor(
  p     : Array[Array[Double]],
  dx    : Double,
  dy    : Double,
  b     : Array[Array[Double]],
  omega : Double,
  nr    : Int,
  nc    : Int,
  nit_  : Int,
) -> Unit {
  let dx2      = dx * dx
  let dy2      = dy * dy
  let inv_denom = 1.0 / (2.0 * (dx2 + dy2))
  for _q = 0; _q < nit_; _q = _q + 1 {
    for i = 1; i < nr - 1; i = i + 1 {
      for j = 1; j < nc - 1; j = j + 1 {
        let laplacian = (p[i][j + 1] + p[i][j - 1]) * dy2
        let laplacian_y = (p[i + 1][j] + p[i - 1][j]) * dx2
        let p_gs = (laplacian + laplacian_y - dx2 * dy2 * b[i][j]) * inv_denom
        p[i][j] = (1.0 - omega) * p[i][j] + omega * p_gs
      }
    }
    apply_pressure_bc_sor(p, nr, nc)
  }
}

/// Symmetric SOR (SSOR) pressure Poisson solver.
///
/// Each SSOR sweep consists of:
///   1. Forward  SOR sweep (i=1..nr-2, j=1..nc-2 ascending)
///   2. Backward SOR sweep (i=nr-2..1, j=nc-2..1 descending)
/// Both half-sweeps use the same ω.  SSOR is symmetric and can reduce
/// effective spectral radius compared to forward-only SOR.
fn pressure_poisson_ssor(
  p     : Array[Array[Double]],
  dx    : Double,
  dy    : Double,
  b     : Array[Array[Double]],
  omega : Double,
  nr    : Int,
  nc    : Int,
  nit_  : Int,
) -> Unit {
  let dx2       = dx * dx
  let dy2       = dy * dy
  let inv_denom = 1.0 / (2.0 * (dx2 + dy2))
  for _q = 0; _q < nit_; _q = _q + 1 {
    // Forward sweep
    for i = 1; i < nr - 1; i = i + 1 {
      for j = 1; j < nc - 1; j = j + 1 {
        let laplacian = (p[i][j + 1] + p[i][j - 1]) * dy2
        let laplacian_y = (p[i + 1][j] + p[i - 1][j]) * dx2
        let p_gs = (laplacian + laplacian_y - dx2 * dy2 * b[i][j]) * inv_denom
        p[i][j] = (1.0 - omega) * p[i][j] + omega * p_gs
      }
    }
    apply_pressure_bc_sor(p, nr, nc)
    // Backward sweep
    let mut i = nr - 2
    while i >= 1 {
      let mut j = nc - 2
      while j >= 1 {
        let laplacian = (p[i][j + 1] + p[i][j - 1]) * dy2
        let laplacian_y = (p[i + 1][j] + p[i - 1][j]) * dx2
        let p_gs = (laplacian + laplacian_y - dx2 * dy2 * b[i][j]) * inv_denom
        p[i][j] = (1.0 - omega) * p[i][j] + omega * p_gs
        j = j - 1
      }
      i = i - 1
    }
    apply_pressure_bc_sor(p, nr, nc)
  }
}

/// Compute the L2 residual of ∇²p − b on a uniform grid.
/// Used to compare convergence quality between GS, SOR, and SSOR.
fn pressure_residual_l2(
  p  : Array[Array[Double]],
  b  : Array[Array[Double]],
  dx : Double,
  dy : Double,
  nr : Int,
  nc : Int,
) -> Double {
  let dx2 = dx * dx
  let dy2 = dy * dy
  let mut sum_sq = 0.0
  let mut cnt    = 0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let lap = (p[i][j + 1] - 2.0 * p[i][j] + p[i][j - 1]) / dx2
      let lap_y = (p[i + 1][j] - 2.0 * p[i][j] + p[i - 1][j]) / dy2
      let res = lap + lap_y - b[i][j]
      sum_sq = sum_sq + res * res
      cnt    = cnt + 1
    }
  }
  if cnt == 0 { return 0.0 }
  (sum_sq / cnt.to_double()).sqrt()
}

/// First-order upwind convective term in y: v * ∂φ/∂y * dt  (backward difference).
fn upwind_convection_y(
  phi   : Array[Array[Double]],
  v_vel : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
  dy    : Double,
  dt    : Double,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      result[i][j] = v_vel[i][j] * dt / dy * (phi[i][j] - phi[i - 1][j])
    }
  }
  result
}

// ─────────────────────────────────────────────────────────────────────────────
// TVD (Total Variation Diminishing) convection schemes
//
// Reference: Sweby P.K. (1984). "High resolution schemes using flux limiters
//   for hyperbolic conservation laws." SIAM J. Numer. Anal. 21, 995-1011.
//
// Limiter ID convention used in tvd_convection_x / tvd_convection_y:
//   1 = Van Leer  (smooth, C¹, TVD-optimal for smooth flows)
//   2 = Superbee  (aggressive, minimal dissipation near shocks)
//
// Face interpolation formula (u > 0, east face j+½ as example):
//   r      = (φ_C − φ_U) / (φ_D − φ_C)          [smoothness ratio]
//   φ_e    = φ_C + ½ · ψ(r) · (φ_C − φ_U)
// where φ_U = φ[i][j-1], φ_C = φ[i][j], φ_D = φ[i][j+1].
// ─────────────────────────────────────────────────────────────────────────────

fn tvd_dabs(x : Double) -> Double {
  if x >= 0.0 { x } else { -x }
}

fn tvd_dmax(a : Double, b : Double) -> Double {
  if a >= b { a } else { b }
}

fn tvd_dmin(a : Double, b : Double) -> Double {
  if a <= b { a } else { b }
}

/// Van Leer flux limiter: ψ(r) = (r + |r|) / (1 + |r|)
/// Properties: ψ(0)=0 (upwind), ψ(1)=1 (central), ψ(∞)→2; smooth (C¹).
fn limiter_van_leer(r : Double) -> Double {
  let ar = tvd_dabs(r)
  (r + ar) / (1.0 + ar)
}

/// Superbee flux limiter: ψ(r) = max(0, min(2r, 1), min(r, 2))
/// Properties: most compressive TVD limiter; minimal numerical diffusion.
fn limiter_superbee(r : Double) -> Double {
  let a = tvd_dmax(0.0, tvd_dmin(2.0 * r, 1.0))
  let b = tvd_dmax(0.0, tvd_dmin(r, 2.0))
  tvd_dmax(a, b)
}

// Apply limiter by id; returns ψ(r).
fn tvd_apply_limiter(r : Double, limiter_id : Int) -> Double {
  if limiter_id == 2 {
    limiter_superbee(r)
  } else {
    limiter_van_leer(r)
  }
}

// Smoothness ratio r for TVD; returns 0.0 when denominator ≈ 0 (uniform region).
fn tvd_ratio(diff_up : Double, diff_dn : Double) -> Double {
  if tvd_dabs(diff_dn) < 1.0e-12 { 0.0 } else { diff_up / diff_dn }
}

// TVD east-face interpolation φ_{j+½} for x-direction.
// Falls back to upwind when the 3-point stencil is incomplete at boundaries.
fn tvd_east_x(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nc   : Int,
  u_ij : Double,
  lid  : Int,
) -> Double {
  if u_ij >= 0.0 {
    // positive flow: upwind=j-1, center=j, downstream=j+1
    if j + 1 >= nc { return phi[i][j] }        // boundary: pure upwind
    let diff_up = phi[i][j]     - phi[i][j - 1]
    let diff_dn = phi[i][j + 1] - phi[i][j]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i][j] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  } else {
    // negative flow: upwind=j+1, center=j+1, upstream-upstream=j+2
    if j + 2 >= nc { return phi[i][j + 1] }    // boundary: pure upwind
    let diff_up = phi[i][j + 1] - phi[i][j + 2]
    let diff_dn = phi[i][j]     - phi[i][j + 1]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i][j + 1] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  }
}

// TVD west-face interpolation φ_{j-½} for x-direction.
fn tvd_west_x(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nc   : Int,
  u_ij : Double,
  lid  : Int,
) -> Double {
  if u_ij >= 0.0 {
    // positive flow: upwind=j-2, center=j-1, downstream=j
    if j - 2 < 0 { return phi[i][j - 1] }      // boundary: pure upwind
    let diff_up = phi[i][j - 1] - phi[i][j - 2]
    let diff_dn = phi[i][j]     - phi[i][j - 1]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i][j - 1] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  } else {
    // negative flow: upwind=j, center=j, upstream-upstream=j+1
    if j + 1 >= nc { return phi[i][j] }         // boundary: pure upwind
    let diff_up = phi[i][j]     - phi[i][j + 1]
    let diff_dn = phi[i][j - 1] - phi[i][j]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i][j] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  }
}

// TVD north-face interpolation φ_{i+½} for y-direction.
fn tvd_north_y(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nr   : Int,
  v_ij : Double,
  lid  : Int,
) -> Double {
  if v_ij >= 0.0 {
    if i + 1 >= nr { return phi[i][j] }
    let diff_up = phi[i][j]     - phi[i - 1][j]
    let diff_dn = phi[i + 1][j] - phi[i][j]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i][j] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  } else {
    if i + 2 >= nr { return phi[i + 1][j] }
    let diff_up = phi[i + 1][j] - phi[i + 2][j]
    let diff_dn = phi[i][j]     - phi[i + 1][j]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i + 1][j] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  }
}

// TVD south-face interpolation φ_{i-½} for y-direction.
fn tvd_south_y(
  phi  : Array[Array[Double]],
  i    : Int,
  j    : Int,
  nr   : Int,
  v_ij : Double,
  lid  : Int,
) -> Double {
  if v_ij >= 0.0 {
    if i - 2 < 0 { return phi[i - 1][j] }
    let diff_up = phi[i - 1][j] - phi[i - 2][j]
    let diff_dn = phi[i][j]     - phi[i - 1][j]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i - 1][j] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  } else {
    if i + 1 >= nr { return phi[i][j] }
    let diff_up = phi[i][j]     - phi[i + 1][j]
    let diff_dn = phi[i - 1][j] - phi[i][j]
    let r       = tvd_ratio(diff_up, diff_dn)
    phi[i][j] + 0.5 * tvd_apply_limiter(r, lid) * diff_up
  }
}

/// TVD convective term in x: u · ∂φ/∂x · dt.
///
/// limiter_id: 1 = Van Leer, 2 = Superbee.
/// Interior nodes with full stencil use TVD face interpolation.
/// Boundary-adjacent nodes lacking the 3-point stencil fall back to upwind.
fn tvd_convection_x(
  phi        : Array[Array[Double]],
  u_vel      : Array[Array[Double]],
  nr         : Int,
  nc         : Int,
  dx         : Double,
  dt         : Double,
  limiter_id : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let u_ij  = u_vel[i][j]
      let phi_e = tvd_east_x(phi, i, j, nc, u_ij, limiter_id)
      let phi_w = tvd_west_x(phi, i, j, nc, u_ij, limiter_id)
      result[i][j] = u_ij * (phi_e - phi_w) / dx * dt
    }
  }
  result
}

/// TVD convective term in y: v · ∂φ/∂y · dt.
///
/// limiter_id: 1 = Van Leer, 2 = Superbee.
fn tvd_convection_y(
  phi        : Array[Array[Double]],
  v_vel      : Array[Array[Double]],
  nr         : Int,
  nc         : Int,
  dy         : Double,
  dt         : Double,
  limiter_id : Int,
) -> Array[Array[Double]] {
  let result = create_zeros_2d(nr, nc)
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let v_ij  = v_vel[i][j]
      let phi_n = tvd_north_y(phi, i, j, nr, v_ij, limiter_id)
      let phi_s = tvd_south_y(phi, i, j, nr, v_ij, limiter_id)
      result[i][j] = v_ij * (phi_n - phi_s) / dy * dt
    }
  }
  result
}
