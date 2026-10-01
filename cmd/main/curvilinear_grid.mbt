// ── Curvilinear Grid Basics (Step 10) ────────────────────────────────────────
//
// Implements a body-fitted curvilinear coordinate framework for structured grids.
//
// ── Coordinate Transformation ─────────────────────────────────────────────────
//
// Computational domain: (ξ, η) ∈ [0,1]×[0,1], uniform spacing.
//   ξ_j = j/(nx−1)   (x-direction)
//   η_i = i/(ny−1)   (y-direction)
//
// Physical domain: (x,y) ∈ [0,2]×[0,2], mapped by an algebraic "bump" formula:
//
//   x(ξ, η) = 2ξ + 2c · ξ(1−ξ) · η(1−η)
//   y(ξ, η) = 2η
//
// Properties:
//   - All four boundaries are identical to the uniform grid (boundary-conforming).
//   - c = 0: identity map (uniform grid recovered exactly).
//   - c > 0: interior nodes shift in +x direction (positive bump).
//   - c > −4: grid is non-degenerate (Jacobian J > 0 everywhere).
//
// Choice of algebraic (polynomial) mapping avoids transcendental functions
// (sin/cos are unavailable in MoonBit stdlib — see SKILL.md §2.2).
//
// ── Jacobian Matrix and Metric Coefficients ────────────────────────────────────
//
// Covariant base vectors:
//   ∂x/∂ξ = x_ξ = 2 + 2c(1−2ξ)η(1−η)
//   ∂x/∂η = x_η = 2c·ξ(1−ξ)(1−2η)
//   ∂y/∂ξ = y_ξ = 0           (y depends only on η)
//   ∂y/∂η = y_η = 2
//
// Jacobian determinant: J = x_ξ·y_η − x_η·y_ξ = 2·x_ξ
//
// Inverse metrics (from Cramer's rule):
//   ξ_x = y_η/J = 2/J = 1/x_ξ
//   ξ_y = −x_η/J
//   η_x = −y_ξ/J = 0
//   η_y = x_ξ/J = 1          (since J = 2·x_ξ and y_η=2)
//
// Covariant metric tensor: g_αβ = (∂x/∂ξ^α)·(∂x/∂ξ^β)
//   g_11 = x_ξ² + y_ξ² = x_ξ²
//   g_12 = x_ξ·x_η + y_ξ·y_η = x_ξ·x_η
//   g_22 = x_η² + y_η² = x_η² + 4
//
// Non-orthogonality measure: cos θ = g_12 / sqrt(g_11·g_22)
//   = 0 when x_ξ·x_η = 0 (orthogonal grid lines)
//
// ── FVM control-volume framework (sketch) ─────────────────────────────────────
//
// For a curvilinear FVM cell centred at (i,j):
//   Face areas:   A_e ≈ J·Δη·e_ξ̂,  A_n ≈ J·Δξ·e_η̂   (contravariant face normals)
//   Cell volume:  V_ij = J·Δξ·Δη
//   The contravariant velocity: U = ξ_x·u + ξ_y·v  (flux through ξ=const face)
//
// Full curvilinear FVM solver is deferred to future work (requires solving the
// transformed N-S equations with metric source terms). This file provides the
// grid generation and metric infrastructure needed for that solver.
//
// Reference: Thompson, J.F., Thames, F.C., Mastin, C.W. (1974). Automatic
//            numerical generation of body-fitted curvilinear coordinate system
//            for field containing any number of arbitrary two-dimensional bodies.
//            J. Comput. Phys. 15, 299–319.

// Global distortion parameter (c ∈ (−4, ∞) for non-degenerate grid)
let g_curv_c : Array[Double] = [0.3]

// ── Pure coordinate and metric functions ──────────────────────────────────────

// Physical x-coordinate at grid node (i = row/y-index, j = col/x-index).
pub fn curv_x(i : Int, j : Int, c : Double) -> Double {
  let xi_j  = j.to_double() / (nx - 1).to_double()
  let eta_i = i.to_double() / (ny - 1).to_double()
  2.0 * xi_j + 2.0 * c * xi_j * (1.0 - xi_j) * eta_i * (1.0 - eta_i)
}

// Physical y-coordinate (no distortion in y: y = 2·η_i).
pub fn curv_y(i : Int) -> Double {
  i.to_double() / (ny - 1).to_double() * 2.0
}

// Analytical Jacobian J = ∂x/∂ξ · ∂y/∂η = 2·(2 + 2c(1−2ξ_j)·η_i(1−η_i)) / 2
// Simplified: J = 2·x_ξ = 4 + 4c·(1−2ξ_j)·η_i·(1−η_i)
pub fn curv_jacobian(i : Int, j : Int, c : Double) -> Double {
  let xi_j  = j.to_double() / (nx - 1).to_double()
  let eta_i = i.to_double() / (ny - 1).to_double()
  4.0 + 4.0 * c * (1.0 - 2.0 * xi_j) * eta_i * (1.0 - eta_i)
}

// Covariant metric components at node (i, j).
pub fn curv_g11(i : Int, j : Int, c : Double) -> Double {
  let xi_j  = j.to_double() / (nx - 1).to_double()
  let eta_i = i.to_double() / (ny - 1).to_double()
  let x_xi  = 2.0 + 2.0 * c * (1.0 - 2.0 * xi_j) * eta_i * (1.0 - eta_i)
  x_xi * x_xi   // g_11 = x_ξ² + y_ξ² = x_ξ² (since y_ξ = 0)
}

pub fn curv_g12(i : Int, j : Int, c : Double) -> Double {
  let xi_j  = j.to_double() / (nx - 1).to_double()
  let eta_i = i.to_double() / (ny - 1).to_double()
  let x_xi  = 2.0 + 2.0 * c * (1.0 - 2.0 * xi_j) * eta_i * (1.0 - eta_i)
  let x_eta = 2.0 * c * xi_j * (1.0 - xi_j) * (1.0 - 2.0 * eta_i)
  x_xi * x_eta  // g_12 = x_ξ·x_η (off-diagonal = 0 when orthogonal)
}

pub fn curv_g22(i : Int, j : Int, c : Double) -> Double {
  let xi_j  = j.to_double() / (nx - 1).to_double()
  let eta_i = i.to_double() / (ny - 1).to_double()
  let x_eta = 2.0 * c * xi_j * (1.0 - xi_j) * (1.0 - 2.0 * eta_i)
  x_eta * x_eta + 4.0   // g_22 = x_η² + y_η² = x_η² + 4
}

// Non-orthogonality cosine: cos θ = |g_12| / sqrt(g_11 · g_22)
// = 0 for orthogonal grid, = 1 for fully skewed (degenerate).
pub fn curv_nonortho(i : Int, j : Int, c : Double) -> Double {
  let g11 = curv_g11(i, j, c)
  let g12 = curv_g12(i, j, c)
  let g22 = curv_g22(i, j, c)
  let denom = (g11 * g22).sqrt()
  if denom < 1.0e-14 { 0.0 } else { (g12 / denom).abs() }
}

// Inverse metric ξ_x = y_η / J = 2 / J  (reciprocal of x_ξ).
pub fn curv_xi_x(i : Int, j : Int, c : Double) -> Double {
  2.0 / curv_jacobian(i, j, c)
}

// Inverse metric ξ_y = −x_η / J.
pub fn curv_xi_y(i : Int, j : Int, c : Double) -> Double {
  let xi_j  = j.to_double() / (nx - 1).to_double()
  let eta_i = i.to_double() / (ny - 1).to_double()
  let x_eta = 2.0 * c * xi_j * (1.0 - xi_j) * (1.0 - 2.0 * eta_i)
  -x_eta / curv_jacobian(i, j, c)
}

// ── WASM API ──────────────────────────────────────────────────────────────────

pub fn set_curv_c(c : Double) -> Unit { g_curv_c[0] = c }
pub fn get_curv_c() -> Double { g_curv_c[0] }

pub fn get_curv_x_at(i : Int, j : Int) -> Double { curv_x(i, j, g_curv_c[0]) }
pub fn get_curv_y_at(i : Int) -> Double { curv_y(i) }
pub fn get_curv_jacobian_at(i : Int, j : Int) -> Double { curv_jacobian(i, j, g_curv_c[0]) }
pub fn get_curv_nonortho_at(i : Int, j : Int) -> Double { curv_nonortho(i, j, g_curv_c[0]) }
pub fn get_curv_g11_at(i : Int, j : Int) -> Double { curv_g11(i, j, g_curv_c[0]) }
pub fn get_curv_g12_at(i : Int, j : Int) -> Double { curv_g12(i, j, g_curv_c[0]) }
pub fn get_curv_g22_at(i : Int, j : Int) -> Double { curv_g22(i, j, g_curv_c[0]) }
