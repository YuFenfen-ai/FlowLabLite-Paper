# PCG Pressure Solver — Design Document

**Project**: FlowLabLite  
**Author**: Fenfen Yu (余芬芬)  
**Date**: 2026-04-16  
**Status**: Phase 1 implemented (PCG on collocated grid); Phase 2 (staggered MAC) planned.

---

## 1  Motivation

The existing pressure Poisson solver uses **Gauss-Seidel (GS)** iteration with a
**fixed** iteration count `nit = 50`.  This has two weaknesses:

| Issue | Impact |
|---|---|
| Fixed 50 sweeps — no convergence criterion | May under- or over-solve; wastes CPU when already converged |
| O(N²) convergence for elliptic problems | Slow for finer grids |
| Non-optimal residual reduction | 50 GS sweeps ≈ first-order accuracy only |

**Conjugate Gradient (CG)** and its preconditioned variant **PCG** are the
industry-standard choices for symmetric positive definite (SPD) linear systems
such as the pressure Poisson equation:

- Optimal convergence in exact arithmetic (converges in at most N iterations for an N-DOF system)
- Residual-based stopping criterion → no wasted iterations
- O(N log N) iterations with a good preconditioner

---

## 2  Mathematical Formulation

### 2.1  Pressure Poisson equation

The fractional-step / projection method (Chorin 1968) requires solving, at each
time step, the pressure Poisson equation:

```
∇²p = b
```

where the right-hand side `b` encodes the velocity divergence:

```
b_{i,j} = ρ/dt · (∂u/∂x + ∂v/∂y) − ρ[(∂u/∂x)² + 2(∂u/∂y)(∂v/∂x) + (∂v/∂y)²]
```

On a uniform grid with spacing `dx`, `dy`, the 5-point finite-difference
discretisation of ∇²p gives:

```
[p_{i,j+1} − 2p_{i,j} + p_{i,j-1}] / dx²
  + [p_{i+1,j} − 2p_{i,j} + p_{i-1,j}] / dy²  =  b_{i,j}
```

This is written in matrix form as **Ap = b**, where:

```
A is (N_int × N_int),  N_int = (nx−2)(ny−2) interior unknowns
Diagonal:            a_{ii} = −(2/dx² + 2/dy²)
Horizontal neighbour: a_{i,j±1} = 1/dx²
Vertical neighbour:   a_{i±1,j} = 1/dy²
```

Boundary conditions for the lid-driven cavity:

| Boundary | Condition | Type |
|---|---|---|
| Left/right walls (x=0, x=Lx) | ∂p/∂x = 0 | Neumann |
| Bottom wall (y=0) | ∂p/∂y = 0 | Neumann |
| Top lid (y=Ly) | p = 0 | Dirichlet |

With at least one Dirichlet BC, the matrix **A is symmetric negative definite**
(equivalently, −A is SPD) → CG applies directly.

### 2.2  PCG algorithm (Jacobi preconditioner)

```
Given:  A x = b,  A is SPD
        M = diag(A)  (Jacobi preconditioner)

x₀ = initial guess (current pressure field)
r₀ = b − A x₀
z₀ = M⁻¹ r₀
d₀ = z₀
ρ₀ = r₀ᵀ z₀

for k = 0, 1, 2, ..., K_max:
    q_k   = A d_k                       ← matrix-free stencil application
    α_k   = ρ_k / (d_k · q_k)          ← step length
    x_{k+1} = x_k + α_k d_k            ← update solution
    r_{k+1} = r_k − α_k q_k            ← update residual
    if ‖r_{k+1}‖ / ‖b‖ < tol  →  stop
    z_{k+1} = M⁻¹ r_{k+1}
    ρ_{k+1} = r_{k+1}ᵀ z_{k+1}
    β_k   = ρ_{k+1} / ρ_k
    d_{k+1} = z_{k+1} + β_k d_k        ← update search direction
```

**Jacobi preconditioner** for the 5-point Laplacian on a uniform interior node:

```
M_{ii} = a_{ii} = 2/dx² + 2/dy²
M⁻¹ r_i = r_i / (2/dx² + 2/dy²)
```

**Theoretical convergence** for a 41×41 grid (1521 interior unknowns):

| Solver | Iterations | Convergence |
|---|---|---|
| Gauss-Seidel (fixed 50) | 50 | O(1/n²) spectral radius |
| CG (no preconditioner) | ≤ N_int = 1521 | exact in ≤ N steps |
| PCG (Jacobi) | O(√κ) ≈ 60–80 | residual tol = 1e-5 |
| PCG (IC(0) or SSOR) | O(√κ/ω) ≈ 20–40 | (future work) |
| PCG (Multigrid) | O(1) ≈ 10–20 | near-optimal |

For the 41×41 cavity at Re=20, GS-50 and PCG-Jacobi give comparable accuracy,
but PCG adapts to the difficulty of each time step while GS uses a fixed budget.

### 2.3  Matrix-free operator

The matrix A is **never assembled**.  At each PCG iteration, the matrix-vector
product Av is computed directly via the 5-point stencil:

```moonbit
// For interior node (i, j):
(Av)[i][j] = (v[i][j+1] + v[i][j-1]) / dx²
           + (v[i+1][j] + v[i-1][j]) / dy²
           - v[i][j] * (2/dx² + 2/dy²)
```

Boundary values of `v` must be set according to the pressure BCs before each
stencil application (Neumann ghost-cell reflection + top Dirichlet zero).

---

## 3  Grid Analysis: Collocated vs Staggered

### 3.1  Collocated grid (current implementation)

All variables `u`, `v`, `p` stored at the **same node positions** `(i·dy, j·dx)`:

```
  p  u  v         (each quantity lives on the same node)
  p  u  v
  p  u  v
```

Grid size: `ny × nx` for each of `u`, `v`, `p`.

**Advantages:**
- Simple, uniform indexing — same `(i,j)` for all variables
- Easy to extend to multi-physics (temperature, concentration)
- Industry standard: OpenFOAM, Fluent (with Rhie-Chow fix)
- Small memory footprint

**Disadvantages:**
- **Checkerboard instability**: pressure-velocity decoupling  
  The central-difference pressure gradient `(p[j+1] − p[j−1]) / 2dx` couples
  cell `j` to cells `j+2` and `j-2`, but not `j+1`.  This allows a
  `±` oscillating pressure mode that produces zero gradient but is not zero pressure.
- **Requires Rhie-Chow interpolation** for stability at higher Re — the current
  code does NOT implement this, limiting reliable operation to low Re.
- Second-order central differences for convective terms — can alias at coarse grids.

### 3.2  Staggered MAC grid (Harlow & Welch 1965)

Variables are stored on **different sub-grids**:

```
For a domain with ny×nx cells:
  p[i][j]  at cell centre   ((i+0.5)dy, (j+0.5)dx)     ny × nx
  u[i][j]  at x-face centre  ((i+0.5)dy,  j·dx)         ny × (nx+1)
  v[i][j]  at y-face centre   (i·dy,     (j+0.5)dx)    (ny+1) × nx
```

Schematic (one cell, pressure at ×, u-faces at →, v-faces at ↑):

```
 v[i+1][j]
    ↑
u[i][j] → × p[i][j]  → u[i][j+1]
    ↑
 v[i][j]
```

**Advantages:**
- **No checkerboard** — pressure and velocity are naturally coupled through
  face-centred differences:  ∂p/∂x at u-face = (p[j] − p[j-1]) / dx  (first difference, no skip)
- **Exact discrete divergence-free condition** after each projection step — by
  construction, the corrected velocity satisfies ∇·u = 0 to machine precision
- More robust for incompressible flows, especially at higher Re
- Recommended for educational/research implementations

**Disadvantages:**
- More complex indexing — u, v, p have different array dimensions
- Requires separate BC arrays for each variable
- Interpolation needed when u and v values are needed at the same point
- More memory: `ny×(nx+1) + (ny+1)×nx + ny×nx` vs `3×ny×nx`
  For 41×41: 1722 + 1722 + 1681 = 5125 vs 5043 — negligible difference

### 3.3  Recommendation and roadmap

| Phase | Grid | Pressure solver | Status |
|---|---|---|---|
| 0 (existing) | Collocated, node-centred | Gauss-Seidel (50 fixed) | ✅ done |
| 1 (this doc) | Collocated (unchanged) | **PCG-Jacobi** (adaptive) | ✅ implemented |
| 2 (planned) | **Staggered MAC** | **PCG-Jacobi** | 🔲 design below |
| 3 (future) | Staggered MAC | PCG + ILU(0) or SSOR | 🔲 planned |

---

## 4  Staggered MAC Grid — Discretisation (Phase 2 Design)

### 4.1  Grid layout

```
nx = ny = 40 cells  (note: cells not nodes; compare with current 41 nodes)
dx = Lx / nx = 2 / 40 = 0.05
dy = Ly / ny = 2 / 40 = 0.05

p [i][j]   i=0..ny-1,  j=0..nx-1
u [i][j]   i=0..ny-1,  j=0..nx     (j=0: left wall face, j=nx: right wall face)
v [i][j]   i=0..ny,    j=0..nx-1   (i=0: bottom wall face, i=ny: top lid face)
```

### 4.2  Boundary conditions (staggered)

```
Left wall  (j=0 u-face):  u[i][0] = 0        for all i
Right wall (j=nx u-face): u[i][nx] = 0       for all i
Bottom     (i=0 v-face):  v[0][j] = 0        for all j
Top lid    (i=ny v-face): v[ny][j] = 0       for all j  (v=0, no normal velocity)
Top lid    u-component:   u[ny-1][j] = U_lid = 1.0  (half-cell below lid face)
  → actually, for MAC the lid condition is set via ghost cells:
    u[ny][j] = 2*U_lid - u[ny-1][j]   (ghost above top face)
  The interior momentum equation at i=ny-1 sees this ghost.
```

### 4.3  Momentum equations (staggered Chorin projection)

**Prediction step** (u*, v* without pressure):

For `u*[i][j]`, i=0..ny-1, j=1..nx-1 (interior u-faces):

```
u*[i][j] = u[i][j]
          − dt * (u²[i][j+1/2] − u²[i][j-1/2]) / dx      ← u-advection in x
          − dt * (u·v[i+1/2][j] − u·v[i-1/2][j]) / dy     ← u-advection in y
          + ν·dt * (u[i][j+1] − 2u[i][j] + u[i][j-1]) / dx²  ← viscosity x
          + ν·dt * (u[i+1][j] − 2u[i][j] + u[i-1][j]) / dy²  ← viscosity y
```

where face-centred velocities require linear interpolation from neighbouring faces.

**Divergence source**:

```
b[i][j] = ρ/dt · ((u*[i][j+1] − u*[i][j]) / dx + (v*[i+1][j] − v*[i][j]) / dy)
```

**PCG pressure solve**: same as Phase 1.

**Correction step**:

```
u[i][j] = u*[i][j] − dt/ρ · (p[i][j] − p[i][j-1]) / dx   ← face-centred gradient
v[i][j] = v*[i][j] − dt/ρ · (p[i][j] − p[i-1][j]) / dy
```

### 4.4  Why the staggered correction guarantees zero divergence

The divergence of the corrected velocity at cell (i,j):

```
D = (u[i][j+1] − u[i][j]) / dx + (v[i+1][j] − v[i][j]) / dy
  = (u*[i][j+1] − u*[i][j]) / dx
    − dt/ρ · (p[i][j+1] − p[i][j]) / dx² + dt/ρ · (p[i][j] − p[i][j-1]) / dx²
  + (v*[i+1][j] − v*[i][j]) / dy
    − dt/ρ · (p[i+1][j] − p[i][j]) / dy² + dt/ρ · (p[i][j] − p[i-1][j]) / dy²
  = ∇·u* − dt/ρ · ∇²p
  = ∇·u* − dt/ρ · b       ← = 0 by construction of b
```

This is exact — not approximate as in the collocated case.

---

## 5  Implementation Steps

### Phase 1 (this commit): PCG on collocated grid

1. **`apply_pressure_bcs(p)`** — extracted helper (used by both GS and PCG)
2. **`laplacian_apply(v, lv, inv_dx2, inv_dy2, a_diag)`** — matrix-free A·v
3. **`pressure_poisson_pcg(p, dx, dy, b) → Int`** — PCG solve, returns iteration count
4. **`cavity_flow_pcg(nt, u, v, dt, dx, dy, p, rho, nu, pcg_iters) → (u,v,p)`** — Chorin + PCG
5. Global state: `g_u_pcg`, `g_v_pcg`, `g_p_pcg`, `g_steps_pcg`, `g_pcg_last_iters`
6. WASM API: `init_chorin_pcg`, `run_chorin_pcg_n_steps`, `get_pcg_*` getters
7. Tests 27–32

### Phase 2 (next session): Staggered MAC grid + PCG

1. New grid arrays: `g_u_mac[ny][nx+1]`, `g_v_mac[ny+1][nx]`, `g_p_mac[ny][nx]`
2. Staggered BCs and momentum predictor
3. Staggered divergence, pressure Poisson (PCG), velocity correction
4. WASM API for MAC solver
5. Validation vs Ghia et al. (1982) benchmark data

### Phase 3 (future): Better preconditioners

- SSOR (symmetric SOR) preconditioner for PCG
- Incomplete Cholesky IC(0)
- Target: ≤ 20 PCG iterations per pressure solve

---

## 6  Test Plan

| Test | Description | Pass Criterion |
|---|---|---|
| T27 | PCG on zero RHS | 0 iterations, p unchanged |
| T28 | PCG BCs satisfied | Neumann/Dirichlet met after solve |
| T29 | Chorin-PCG wall BCs | u=0 on walls, u=1 on lid after step |
| T30 | Chorin-PCG lid velocity | u[ny-1][j] = 1.0 maintained |
| T31 | Chorin-PCG independence | Does not modify Chorin/SIMPLE state |
| T32 | PCG iteration bound | iters ≤ pcg_max_iter = 200 |
| T33 | PCG vs GS pressure | ‖p_pcg − p_gs‖/‖p_gs‖ < 0.05 after 1 step |
| T34 | Chorin-PCG vortex | Negative u at centre after 200 steps |

---

## 7  References

1. Chorin, A.J. (1968). *Numerical solution of the Navier-Stokes equations.*
   Math. Comp., 22(104), 745–762.
2. Harlow, F.H. & Welch, J.E. (1965). *Numerical calculation of time-dependent
   viscous incompressible flow.* Phys. Fluids, 8(12), 2182.
3. Hestenes, M.R. & Stiefel, E. (1952). *Methods of conjugate gradients for
   solving linear systems.* J. Res. NIST, 49(6), 409–436.
4. Patankar, S.V. & Spalding, D.B. (1972). *A calculation procedure for heat,
   mass and momentum transfer in three-dimensional parabolic flows.*
   Int. J. Heat Mass Transfer, 15(10), 1787–1806.
5. Ferziger, J.H., Perić, M. & Street, R.L. (2020). *Computational Methods for
   Fluid Dynamics.* 4th ed., Springer.
6. Ghia, U., Ghia, K.N. & Shin, C.T. (1982). *High-Re solutions for
   incompressible flow using the Navier-Stokes equations and a multigrid method.*
   J. Comput. Phys., 48(3), 387–411.
