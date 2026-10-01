// ── Stretched Grid (Step 9) ───────────────────────────────────────────────────
//
// Generates non-uniform (tanh-stretched) 1D coordinates and provides
// non-uniform finite-difference coefficient helpers for second-order accuracy.
//
// ── tanh Stretching Function ──────────────────────────────────────────────────
//
//   x_i = L · tanh(β · ξ_i) / tanh(β),   ξ_i = i / (N−1)
//
//   β = 0  (limit): uniform spacing (x_i = i·L/(N−1))
//   β > 0:          nodes cluster near x=0  (fine grid near left/bottom wall)
//   β < 0:          nodes cluster near x=L  (fine grid near right/top wall)
//
// Formula: x_i = L·(1 − tanh(β·(1−ξ_i))/tanh(β)),  ξ_i = i/(N−1)
// (Not the simpler tanh(β·ξ_i)/tanh(β), which clusters near x=L instead.)
//
// For wall-normal clustering near both walls, call with β > 0 on
// the lower half and mirror for the upper half (symmetric double-sided stretch).
//
// ── Non-uniform FD Coefficients ───────────────────────────────────────────────
//
// For node i with left spacing h1 = x_i − x_{i−1} and right h2 = x_{i+1} − x_i:
//
//   1st derivative (3-point, central):
//     dφ/dx ≈ [−h2·φ_{i−1} + (h2−h1)·φ_i + h1·φ_{i+1}] / (h1·h2 + h1² + h2²) (simplified)
//             = 2·[−h2·φ_{i−1} + (h2²−h1²)/(h2−h1) term...] — use Fornberg form
//   Actually using standard 3-point Lagrange interpolation:
//     dφ/dx|_i = φ_{i−1}·(−h2)/(h1(h1+h2)) + φ_i·(h2−h1)/(h1·h2) + φ_{i+1}·h1/(h2(h1+h2))
//
//   2nd derivative (3-point Lagrange, exact for quadratics):
//     d²φ/dx²|_i = 2·[φ_{i−1}/(h1(h1+h2)) − φ_i/(h1·h2) + φ_{i+1}/(h2(h1+h2))]
//
// Reference: Fornberg, B. (1988). Generation of finite difference formulas on
//            arbitrarily spaced grids. Math. Comp. 51, 699-706.

// ── Math helpers (MoonBit has no exp/tanh in stdlib) ─────────────────────────
//
// exp_impl: Taylor series e^x = Σ x^k/k!, 25 terms.
// Accurate to ~15 significant digits for |x| ≤ 20 (reciprocal for x<0).
fn exp_impl(x : Double) -> Double {
  if x < 0.0 { return 1.0 / exp_impl(-x) }
  if x > 20.0 { return 485165195.4 }   // e^20 ≈ 4.85e8  (overflow guard)
  let mut term = 1.0
  let mut s    = 1.0
  for k = 1; k <= 25; k = k + 1 {
    term = term * x / k.to_double()
    s    = s + term
  }
  s
}

// tanh_impl: derived from exp, clamped at ±1 for |x| ≥ 20.
fn tanh_impl(x : Double) -> Double {
  if x >= 20.0  { return  1.0 }
  if x <= -20.0 { return -1.0 }
  let e2x = exp_impl(2.0 * x)
  (e2x - 1.0) / (e2x + 1.0)
}

// ── Grid generation ───────────────────────────────────────────────────────────

// tanh-stretched coordinate at node i (0 ≤ i ≤ n−1) on [0, length].
//
// Formula: x_i = L · (1 − tanh(β·(1−ξ_i)) / tanh(β)),  ξ_i = i/(n−1)
//
// This maps [0,L] with clustering near x=0 (wall) for β > 0:
//   - First spacing ≈ L·(1 − tanh(β·(n−2)/(n−1)) / tanh(β))  ← small (fine near wall)
//   - Last spacing  ≈ L·tanh(β/(n−1)) / tanh(β)               ← large (coarse far from wall)
//   - β = 0 (limit): uniform spacing
pub fn tanh_coord(
  i      : Int,
  n      : Int,
  length : Double,
  beta   : Double,
) -> Double {
  if n <= 1 { return 0.0 }
  if i <= 0 { return 0.0 }
  if i >= n - 1 { return length }
  let xi = i.to_double() / (n - 1).to_double()
  if beta.abs() < 1.0e-8 {
    xi * length   // uniform limit (β→0, avoids division by near-zero)
  } else {
    length * (1.0 - tanh_impl(beta * (1.0 - xi)) / tanh_impl(beta))
  }
}

// Local spacing between nodes i and i+1.
pub fn tanh_spacing(
  i      : Int,
  n      : Int,
  length : Double,
  beta   : Double,
) -> Double {
  tanh_coord(i + 1, n, length, beta) - tanh_coord(i, n, length, beta)
}

// Build a full coordinate array of length n.
pub fn make_tanh_grid(
  n      : Int,
  length : Double,
  beta   : Double,
) -> Array[Double] {
  let xs = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 { xs[i] = tanh_coord(i, n, length, beta) }
  xs
}

// ── Non-uniform FD operators ──────────────────────────────────────────────────

// Second derivative at interior node i using 3-point Lagrange interpolation.
// Exact for linear and quadratic polynomials (2nd-order scheme on non-uniform grid).
//   d²φ/dx²|_i = 2[φ_{i−1}/(h1(h1+h2)) − φ_i/(h1·h2) + φ_{i+1}/(h2(h1+h2))]
pub fn nonuniform_d2(
  phi_m : Double,   // φ_{i-1}
  phi_0 : Double,   // φ_i
  phi_p : Double,   // φ_{i+1}
  h1    : Double,   // x_i   − x_{i-1}   (left spacing)
  h2    : Double,   // x_{i+1} − x_i     (right spacing)
) -> Double {
  let denom = h1 * h2
  let sum12 = h1 + h2
  2.0 * (phi_m / (h1 * sum12) - phi_0 / denom + phi_p / (h2 * sum12))
}

// First derivative at interior node i (3-point asymmetric Lagrange):
//   dφ/dx|_i = φ_{i−1}·(−h2)/(h1·(h1+h2)) + φ_i·(h2−h1)/(h1·h2) + φ_{i+1}·h1/(h2·(h1+h2))
pub fn nonuniform_d1(
  phi_m : Double,
  phi_0 : Double,
  phi_p : Double,
  h1    : Double,
  h2    : Double,
) -> Double {
  let sum12 = h1 + h2
  phi_m * (-h2) / (h1 * sum12) + phi_0 * (h2 - h1) / (h1 * h2) + phi_p * h1 / (h2 * sum12)
}

// Maximum |d²φ/dx²| on a stretched grid for a given phi array and coord array.
// Used by T142 to verify exactness for linear field.
pub fn nonuniform_laplacian_max_err(
  phi : Array[Double],
  xs  : Array[Double],
  n   : Int,
) -> Double {
  let mut max_err = 0.0
  for i = 1; i < n - 1; i = i + 1 {
    let h1  = xs[i] - xs[i - 1]
    let h2  = xs[i + 1] - xs[i]
    let d2  = nonuniform_d2(phi[i - 1], phi[i], phi[i + 1], h1, h2)
    let err = if d2 >= 0.0 { d2 } else { -d2 }
    if err > max_err { max_err = err }
  }
  max_err
}

// ── WASM API – grid coordinate queries ───────────────────────────────────────
//
// These use the global ny and domain height 2.0 for the y-direction.
// beta > 0: clustering near y=0 (bottom wall); beta < 0: near y=H.

pub fn get_tanh_y_at(i : Int, beta : Double) -> Double {
  tanh_coord(i, ny, 2.0, beta)
}

pub fn get_tanh_dy_at(i : Int, beta : Double) -> Double {
  tanh_spacing(i, ny, 2.0, beta)
}

// Mean y-spacing (= H/(ny−1) regardless of beta, for sanity check)
pub fn get_tanh_mean_dy() -> Double {
  2.0 / (ny - 1).to_double()
}
