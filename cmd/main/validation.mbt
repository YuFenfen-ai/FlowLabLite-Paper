// FlowLabLite — Ghia (1982) Benchmark Validation
//
// Reference: Ghia, U., Ghia, K.N., Shin, C.T. (1982).
//   "High-Re solutions for incompressible flow using the Navier-Stokes equations
//    and a multigrid method." J. Comput. Phys. 48, 387–411.
//
// Domain mapping: Ghia uses unit square [0,1]×[0,1].
//   Our solver uses [0,2]×[0,2] → multiply Ghia coordinates by 2.0.
//
// Data embedded: Re=100, Table 1 (u on vertical centerline x=0.5),
//                        Table 2 (v on horizontal centerline y=0.5).
// Both tables have 17 reference points.

// ─────────────────────────────────────────────────────────────────────────────
// Ghia Re=100 reference data  (17 points each)
// ─────────────────────────────────────────────────────────────────────────────

// Vertical centerline (x/H = 0.5): y/H coordinates in Ghia's [0,1] domain.
// Multiply by 2.0 to map to our [0,2] domain.
let ghia_re100_yn : Array[Double] = [
  0.0000, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813, 0.4531, 0.5000,
  0.6172, 0.7344, 0.8516, 0.9531, 0.9609, 0.9688, 0.9766, 1.0000,
]

// u/U values at those y/H coordinates (Re=100, Table 1).
// u=0 at y=0 (bottom wall), u=1 at y=1 (top lid), non-dimensional.
let ghia_re100_u : Array[Double] = [
   0.00000, -0.03717, -0.04192, -0.04775, -0.06434, -0.10150, -0.15662,
  -0.21090, -0.20581, -0.13641,  0.00332,  0.23151,  0.68717,  0.73722,
   0.78871,  0.84123,  1.00000,
]

// Horizontal centerline (y/H = 0.5): x/H coordinates in Ghia's [0,1] domain.
let ghia_re100_xn : Array[Double] = [
  0.0000, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266, 0.2344, 0.5000,
  0.8047, 0.8594, 0.9063, 0.9453, 0.9531, 0.9609, 0.9688, 1.0000,
]

// v/U values at those x/H coordinates (Re=100, Table 2).
// v=0 at x=0 (left wall) and x=1 (right wall).
let ghia_re100_v : Array[Double] = [
   0.00000,  0.09233,  0.10091,  0.10890,  0.12317,  0.16077,  0.17507,
   0.17527,  0.05454, -0.24533, -0.22445, -0.16914, -0.10313, -0.08864,
  -0.07391, -0.05906,  0.00000,
]

// ─────────────────────────────────────────────────────────────────────────────
// Linear interpolation helper
// ─────────────────────────────────────────────────────────────────────────────

/// Linear interpolation of a 1-D field `vals` sampled at uniform grid
/// y[i] = y0 + i * h, i = 0..n-1.  Clamps to boundary if y_ref is outside.
fn lerp_uniform(vals : Array[Double], y0 : Double, h : Double, y_ref : Double) -> Double {
  let n = vals.length()
  if n == 0 { return 0.0 }
  // Fractional index
  let frac_idx = (y_ref - y0) / h
  let i0 = frac_idx.to_int()
  // Clamp to valid range
  if i0 < 0         { return vals[0] }
  if i0 >= n - 1    { return vals[n - 1] }
  let t = frac_idx - i0.to_double()
  vals[i0] * (1.0 - t) + vals[i0 + 1] * t
}

// ─────────────────────────────────────────────────────────────────────────────
// Ghia comparison result
// ─────────────────────────────────────────────────────────────────────────────

struct GhiaResult {
  l2_u  : Double   // L2 error on u vertical centerline
  max_u : Double   // max |error| on u vertical centerline
  l2_v  : Double   // L2 error on v horizontal centerline
  max_v : Double   // max |error| on v horizontal centerline
  n_pts : Int      // number of reference points used
}

// ─────────────────────────────────────────────────────────────────────────────
// Main comparison function
// ─────────────────────────────────────────────────────────────────────────────

/// Compare a solver's centerline velocities against Ghia (1982) Re=100 data.
///
/// Interpolates the solver solution linearly to each of the 17 Ghia reference
/// points.  Ghia coordinates are mapped from [0,1] to our [0,2] domain by ×2.
///
/// Returns L2 and max errors for both centerlines.
/// Only meaningful for Re=100 runs; for other Re the comparison is returned
/// but the errors will be larger (no reference data for other Re is embedded).
fn ghia_compare(s : SolverExport) -> GhiaResult {
  let ng  = ghia_re100_yn.length()   // 17 reference points
  let jmid = s.nc / 2                // vertical centerline column
  let imid = s.nc / 2                // horizontal centerline row

  // Extract our u along vertical centerline into a 1-D array
  let u_vcl : Array[Double] = Array::make(s.nc, 0.0)
  for i = 0; i < s.nc; i = i + 1 {
    u_vcl[i] = s.u[i][jmid]
  }

  // Extract our v along horizontal centerline into a 1-D array
  let v_hcl : Array[Double] = Array::make(s.nc, 0.0)
  for j = 0; j < s.nc; j = j + 1 {
    v_hcl[j] = s.v[imid][j]
  }

  // Compute L2 and max errors for u (vertical centerline)
  let mut sum_sq_u = 0.0
  let mut max_u    = 0.0
  for k = 0; k < ng; k = k + 1 {
    let y_ref = ghia_re100_yn[k] * 2.0   // map [0,1] → [0,2]
    let u_our = lerp_uniform(u_vcl, s.y0, s.h, y_ref)
    let err   = (u_our - ghia_re100_u[k]).abs()
    sum_sq_u = sum_sq_u + err * err
    if err > max_u { max_u = err }
  }
  let l2_u = (sum_sq_u / ng.to_double()).sqrt()

  // Compute L2 and max errors for v (horizontal centerline)
  let mut sum_sq_v = 0.0
  let mut max_v    = 0.0
  for k = 0; k < ng; k = k + 1 {
    let x_ref = ghia_re100_xn[k] * 2.0   // map [0,1] → [0,2]
    let v_our = lerp_uniform(v_hcl, s.x0, s.h, x_ref)
    let err   = (v_our - ghia_re100_v[k]).abs()
    sum_sq_v = sum_sq_v + err * err
    if err > max_v { max_v = err }
  }
  let l2_v = (sum_sq_v / ng.to_double()).sqrt()

  { l2_u, max_u, l2_v, max_v, n_pts: ng }
}

/// Format a GhiaResult as a human-readable summary string for CLI output.
fn format_ghia_result(
  result      : GhiaResult,
  solver_name : String,
  re          : Double,
  steps       : Int,
) -> String {
  let buf = StringBuilder::new()
  buf.write_string("Ghia (1982) validation | solver=\{solver_name} | Re=\{re} | step=\{steps}\n")
  buf.write_string("  Reference: Ghia, Ghia, Shin (1982) J. Comput. Phys. 48, 387–411\n")
  buf.write_string("  Points compared: \{result.n_pts} (Re=100 data)\n")
  buf.write_string("\n")
  buf.write_string("  u-velocity (vertical centerline, x=L/2):\n")
  buf.write_string("    L2  error : \{result.l2_u}\n")
  buf.write_string("    max error : \{result.max_u}\n")
  buf.write_string("\n")
  buf.write_string("  v-velocity (horizontal centerline, y=L/2):\n")
  buf.write_string("    L2  error : \{result.l2_v}\n")
  buf.write_string("    max error : \{result.max_v}\n")
  buf.to_string()
}
