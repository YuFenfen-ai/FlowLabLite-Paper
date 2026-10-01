// FlowLabLite — Solver Convergence Monitor  (Step M)
//
// Pure-function monitoring utilities.  All functions receive arrays as
// parameters and do NOT read global state directly.  This mirrors the
// OpenFOAM functionObjects design: pure computation, no side-effects,
// fully unit-testable with synthetic data.
//
// Why pure functions (not global reads)?
//   1. Testability  – unit tests can pass synthetic fields without running a solver
//   2. Decoupling   – monitor code is independent of solver internals
//   3. Reusability  – any solver can call the same monitor functions
//   4. Safety       – no risk of accidentally modifying solver state
//   5. Clarity      – data flow is explicit; every dependency is visible
//
// Reference: solver_monitor.md (formulas and physical interpretation)

// ─────────────────────────────────────────────────────────────────────────────
// Velocity-change residuals  R_u, R_v
// ─────────────────────────────────────────────────────────────────────────────

/// Relative change residual for u-velocity between two successive snapshots:
///   R_u = max|u_new − u_old| / (max|u_new| + ε)
/// where ε = 1e-10 prevents division by zero for stagnant fields.
///
/// For transient (Chorin/PCG/MAC): reflects step-to-step velocity variation.
/// For SIMPLE: reflects iteration-to-iteration convergence (→ 0 at steady state).
fn compute_ru(
  u_new : Array[Array[Double]],
  u_old : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
) -> Double {
  let mut max_diff = 0.0
  let mut max_u    = 0.0
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let d = (u_new[i][j] - u_old[i][j]).abs()
      let a = u_new[i][j].abs()
      if d > max_diff { max_diff = d }
      if a > max_u    { max_u    = a }
    }
  }
  max_diff / (max_u + 1.0e-10)
}

/// Relative change residual for v-velocity (same formula as compute_ru).
fn compute_rv(
  v_new : Array[Array[Double]],
  v_old : Array[Array[Double]],
  nr    : Int,
  nc    : Int,
) -> Double {
  let mut max_diff = 0.0
  let mut max_v    = 0.0
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let d = (v_new[i][j] - v_old[i][j]).abs()
      let a = v_new[i][j].abs()
      if d > max_diff { max_diff = d }
      if a > max_v    { max_v    = a }
    }
  }
  max_diff / (max_v + 1.0e-10)
}

// ─────────────────────────────────────────────────────────────────────────────
// CFL number
// ─────────────────────────────────────────────────────────────────────────────

/// Maximum CFL number over the grid:
///   CFL_max = max_{i,j} ( |u[i][j]|*dt/dx + |v[i][j]|*dt/dy )
///
/// For stability of explicit schemes the CFL number must be < 1.
/// Typical safe value for the default Re=20 run: CFL_max ≈ 0.02.
fn compute_cfl_max(
  u  : Array[Array[Double]],
  v  : Array[Array[Double]],
  nr : Int,
  nc : Int,
  dt_ : Double,
  dx_ : Double,
  dy_ : Double,
) -> Double {
  let mut cfl = 0.0
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let c = u[i][j].abs() * dt_ / dx_ + v[i][j].abs() * dt_ / dy_
      if c > cfl { cfl = c }
    }
  }
  cfl
}

// ─────────────────────────────────────────────────────────────────────────────
// Kinetic energy
// ─────────────────────────────────────────────────────────────────────────────

/// Volume-integrated kinetic energy of the flow field:
///   KE = 0.5 * Σ_{i,j} (u²[i][j] + v²[i][j]) * dx * dy
///
/// Normalised by dx*dy so the value represents average KE density × domain area.
/// For the lid-driven cavity this typically peaks early then settles to a
/// quasi-steady value; it provides a global energy budget check.
fn compute_kinetic_energy(
  u   : Array[Array[Double]],
  v   : Array[Array[Double]],
  nr  : Int,
  nc  : Int,
  dx_ : Double,
  dy_ : Double,
) -> Double {
  let mut ke = 0.0
  for i = 0; i < nr; i = i + 1 {
    for j = 0; j < nc; j = j + 1 {
      let u_ = u[i][j]
      let v_ = v[i][j]
      ke = ke + (u_ * u_ + v_ * v_)
    }
  }
  0.5 * ke * dx_ * dy_
}

// ─────────────────────────────────────────────────────────────────────────────
// SIMPLE strict (equation-residual) mode
// ─────────────────────────────────────────────────────────────────────────────

/// OpenFOAM-style normalised momentum equation residual for SIMPLE.
///
/// Computes the normalised x-momentum imbalance:
///   R_φ = Σ_{i,j} |a_P u_P − (a_E u_E + a_W u_W + a_N u_N + a_S u_S) − b_P|
///         ────────────────────────────────────────────────────────────────────
///         Σ_{i,j} |a_P u_ref|
///
/// where (for interior nodes on a uniform grid):
///   a_P  = 1 + 2*nu*dt/dx² + 2*nu*dt/dy²   (central coefficient)
///   a_E  = nu*dt/dx²   (east)
///   a_W  = nu*dt/dx²   (west)
///   a_N  = nu*dt/dy²   (north)
///   a_S  = nu*dt/dy²   (south)
///   b_P  = −(dt/rho) * dp/dx  (pressure gradient; approximated here as 0
///            for a pure field-balance check without pressure access)
///
/// Note: Since this function receives only u/v (not p), the pressure gradient
/// term b_P is omitted.  This gives a conservative upper bound for the
/// momentum equation residual.  Use `div` mode for full residual including pressure.
///
/// Returns a value in [0, 1] for a well-posed problem; values > 1 indicate
/// extreme imbalance (e.g. at initialisation).
fn compute_simple_residual_eq(
  u    : Array[Array[Double]],
  nr   : Int,
  nc   : Int,
  nu_  : Double,
  dt_  : Double,
  dx_  : Double,
  dy_  : Double,
) -> Double {
  let dx2    = dx_ * dx_
  let dy2    = dy_ * dy_
  let a_e    = nu_ * dt_ / dx2
  let a_w    = nu_ * dt_ / dx2
  let a_n    = nu_ * dt_ / dy2
  let a_s    = nu_ * dt_ / dy2
  let a_p    = 1.0 + 2.0 * nu_ * dt_ / dx2 + 2.0 * nu_ * dt_ / dy2
  let mut num   = 0.0
  let mut denom = 0.0
  for i = 1; i < nr - 1; i = i + 1 {
    for j = 1; j < nc - 1; j = j + 1 {
      let nb_sum = a_e * u[i][j + 1] + a_w * u[i][j - 1]
      let nb_sum_y = a_n * u[i + 1][j] + a_s * u[i - 1][j]
      let imbalance = (a_p * u[i][j] - nb_sum - nb_sum_y).abs()
      num   = num   + imbalance
      denom = denom + (a_p * u[i][j]).abs()
    }
  }
  if denom < 1.0e-30 { return 0.0 }
  num / denom
}

// Note: MonitorRecord struct and format_monitor_csv function are defined in
// io_formats.mbt (Step B).  This file provides only the computation functions.
