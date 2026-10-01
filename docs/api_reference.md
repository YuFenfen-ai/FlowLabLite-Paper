# FlowLabLite — WASM API Reference

All 121 exported functions are callable from JavaScript after loading the WASM module.
Functions are grouped by solver. Out-of-bounds indices return `0` (never throw).

---

## Common conventions

```javascript
// Load the module
const importObject = /* see INSTALL.md */;
const { instance } = await WebAssembly.instantiate(wasmBytes, importObject);
const w = instance.exports;

// Coordinate convention
//   i = row index,  0 = bottom row,  ny-1 = top (lid) row
//   j = col index,  0 = left col,    nx-1 = right col
// Out-of-bounds: i < 0 || i >= ny || j < 0 || j >= nx → returns 0.0
```

---

## Chorin Solver (27 functions)

### Lifecycle

#### `init_simulation() → void`
Zero-fills all Chorin global arrays (`g_u`, `g_v`, `g_p`) and resets the step counter to 0.
Must be called before the first `run_n_steps`.

```javascript
w.init_simulation();
```

#### `run_n_steps(n: i32) → void`
Advance the Chorin projection solver by `n` time steps (each step = dt = 0.001 s).
Internally calls `cavity_flow_array()` once per step.

```javascript
w.run_n_steps(50);   // advance 50 steps
```

#### `run_all_steps() → void`
Run the full `nt = 500` steps in one call. Equivalent to `run_n_steps(nt)`.

---

### Grid constants

#### `get_nx() → i32`  Returns 41 (grid nodes in x-direction).
#### `get_ny() → i32`  Returns 41 (grid nodes in y-direction).
#### `get_nt() → i32`  Returns 500 (total time steps).
#### `get_nit() → i32` Returns 50 (Gauss-Seidel iterations per pressure solve).

---

### Physical constants

#### `get_re() → f64`   Returns 20.0 (Reynolds number = U·L/ν).
#### `get_dx() → f64`   Returns 0.05 (grid spacing in x, metres).
#### `get_dy() → f64`   Returns 0.05 (grid spacing in y, metres).
#### `get_dt() → f64`   Returns 0.001 (time step, seconds).
#### `get_rho() → f64`  Returns 1.0 (fluid density, kg/m³).
#### `get_nu() → f64`   Returns 0.1 (kinematic viscosity, m²/s).

---

### State accessors

#### `get_step_count() → i32`
Number of time steps executed since the last `init_simulation()`.

#### `get_u_at(i: i32, j: i32) → f64`
x-velocity at node (i, j). Returns 0 if out of range.

#### `get_v_at(i: i32, j: i32) → f64`
y-velocity at node (i, j). Returns 0 if out of range.

#### `get_p_at(i: i32, j: i32) → f64`
Pressure at node (i, j). Returns 0 if out of range.

#### `get_velocity_magnitude_at(i: i32, j: i32) → f64`
`sqrt(u² + v²)` at node (i, j). Returns 0 if out of range.

#### `get_max_velocity_magnitude() → f64`
Maximum velocity magnitude over all grid nodes. Always ≥ 0.

#### `get_u_center() → f64`  u at the centre node (ny/2, nx/2). Equivalent to `get_u_at(20, 20)`.
#### `get_v_center() → f64`  v at the centre node.
#### `get_p_center() → f64`  Pressure at the centre node.

---

### Diagnostic functions

#### `get_divergence_norm() → f64`
Mean |∇·u| over all interior cells. Should approach 0 as the solver converges.
Defined as: `Σ |∂u/∂x + ∂v/∂y|_ij / N_interior`.

#### `get_max_u() → f64`  Maximum u over all grid nodes.
#### `get_max_v() → f64`  Maximum v over all grid nodes.
#### `get_max_p() → f64`  Maximum pressure over all grid nodes.
#### `get_min_p() → f64`  Minimum pressure (typically negative for the Chorin/SND convention).

---

## SIMPLE Solver (11 functions)

The SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) solver uses
independent global state `g_u_s / g_v_s / g_p_s`. Running SIMPLE never
modifies Chorin state.

### Lifecycle

#### `init_simple() → void`
Zero-fills SIMPLE fields and resets iteration counter.

#### `run_simple_n_iter(n: i32) → void`
Run `n` SIMPLE iterations. Each iteration: momentum predictor → pressure correction → velocity update.

---

### State accessors

#### `get_simple_step_count() → i32`  Iterations since last `init_simple()`.
#### `get_simple_residual() → f64`    Last pressure-correction residual (≥ 0).
#### `get_u_simple_at(i: i32, j: i32) → f64`  x-velocity at (i, j). Returns 0 if OOB.
#### `get_v_simple_at(i: i32, j: i32) → f64`  y-velocity at (i, j). Returns 0 if OOB.
#### `get_p_simple_at(i: i32, j: i32) → f64`  Pressure at (i, j). Returns 0 if OOB.
#### `get_max_u_simple() → f64`               Maximum u over all SIMPLE nodes.
#### `get_simple_divergence_norm() → f64`      Mean |∇·u| over SIMPLE interior cells.

---

## Chorin-PCG Solver (13 functions)

Same fractional-step algorithm as Chorin but replaces Gauss-Seidel with
**Preconditioned Conjugate Gradient** (PCG, Jacobi preconditioner).
Tolerance: 1×10⁻⁵; max iterations: 200. Independent state `g_u_pcg / g_v_pcg / g_p_pcg`.

### Lifecycle

#### `init_chorin_pcg() → void`        Zero-fills PCG fields and resets counter.
#### `run_chorin_pcg_n_steps(n: i32) → void`  Advance PCG solver by n steps.

---

### State accessors

#### `get_pcg_step_count() → i32`                  Steps since last init.
#### `get_pcg_last_iters() → i32`                  PCG iterations used in the most recent pressure solve. Typically 50–130.
#### `get_u_pcg_at(i: i32, j: i32) → f64`         x-velocity. Returns 0 if OOB.
#### `get_v_pcg_at(i: i32, j: i32) → f64`         y-velocity. Returns 0 if OOB.
#### `get_p_pcg_at(i: i32, j: i32) → f64`         Pressure. Returns 0 if OOB.
#### `get_velocity_magnitude_pcg_at(i: i32, j: i32) → f64`  sqrt(u²+v²). Returns 0 if OOB.
#### `get_max_u_pcg() → f64`   Maximum u over all PCG nodes.
#### `get_max_v_pcg() → f64`   Maximum v over all PCG nodes.
#### `get_max_p_pcg() → f64`   Maximum pressure.
#### `get_min_p_pcg() → f64`   Minimum pressure.
#### `get_pcg_divergence_norm() → f64`  Mean |∇·u| over interior cells. < 1×10⁻⁴ after 50 steps.

---

## MAC Staggered Solver (14 functions)

Harlow-Welch MAC staggered grid. Pressure at cell centres, velocities at face centres.
Grid: `mac_nc = 40` cells per direction. Independent state `g_u_mac / g_v_mac / g_p_mac`.

The staggered arrangement eliminates the pressure-velocity decoupling problem and
guarantees divergence-free interior cells by construction.

### Grid layout

```
p[i][j]   : cell centre  (mac_nc × mac_nc = 40×40)
u_mac[i][j]: x-face centre (mac_nc × (mac_nc+1) = 40×41)  — i-th row, j-th x-face
v_mac[i][j]: y-face centre ((mac_nc+1) × mac_nc = 41×40)  — i-th y-face, j-th col
```

`get_u_mac_at(i,j)` / `get_v_mac_at(i,j)` / `get_p_mac_at(i,j)` all use
**cell-centre indices** (i ∈ [0, mac_nc), j ∈ [0, mac_nc)). The getter
performs face-to-centre interpolation internally where needed.

### Lifecycle

#### `init_mac() → void`                      Zero-fills MAC fields and resets counter.
#### `run_mac_n_steps(n: i32) → void`          Advance MAC solver by n steps.

---

### State accessors

#### `get_mac_step_count() → i32`              Steps since last init.
#### `get_mac_last_iters() → i32`              PCG iterations in most recent pressure solve.
#### `get_mac_nc() → i32`                      Returns 40 (cells per direction).
#### `get_u_mac_at(i: i32, j: i32) → f64`     x-velocity at cell (i,j). Returns 0 if OOB.
#### `get_v_mac_at(i: i32, j: i32) → f64`     y-velocity at cell (i,j). Returns 0 if OOB.
#### `get_p_mac_at(i: i32, j: i32) → f64`     Pressure at cell (i,j). Returns 0 if OOB.
#### `get_velocity_magnitude_mac_at(i: i32, j: i32) → f64`  sqrt(u²+v²) at cell (i,j).
#### `get_max_u_mac() → f64`    Maximum u (cell-centre) over all MAC cells.
#### `get_max_v_mac() → f64`    Maximum v.
#### `get_max_p_mac() → f64`    Maximum pressure.
#### `get_min_p_mac() → f64`    Minimum pressure.
#### `get_mac_divergence_norm() → f64`  Mean |∇·u| over PCG-interior cells (i,j = 1..38). Guaranteed < 1×10⁻⁴ after 50 steps.

---

## JavaScript Usage Examples

### Example 1 — Run Chorin and read the velocity field

```javascript
const { instance } = await WebAssembly.instantiate(wasmBytes, importObject);
const w = instance.exports;

// Initialise and run 100 steps
w.init_simulation();
w.run_n_steps(100);

const nx = w.get_nx();   // 41
const ny = w.get_ny();   // 41

// Read all u-velocities into a typed array
const u = new Float64Array(nx * ny);
for (let i = 0; i < ny; i++) {
  for (let j = 0; j < nx; j++) {
    u[i * nx + j] = w.get_u_at(i, j);
  }
}

console.log('Centre u:', w.get_u_at(20, 20));   // ≈ −0.10 after 100 steps
console.log('Divergence norm:', w.get_divergence_norm());
```

### Example 2 — Switch solvers and compare

```javascript
// Run both Chorin and SIMPLE independently
w.init_simulation();   w.run_n_steps(500);
const u_chorin = w.get_u_at(20, 20);

w.init_simple();       w.run_simple_n_iter(100);
const u_simple = w.get_u_simple_at(20, 20);

console.log('Chorin u_centre:', u_chorin);   // ≈ −0.10
console.log('SIMPLE u_centre:', u_simple);   // similar sign, different magnitude
```

### Example 3 — Monitor PCG convergence

```javascript
w.init_chorin_pcg();
for (let chunk = 0; chunk < 10; chunk++) {
  w.run_chorin_pcg_n_steps(10);
  console.log(
    `Step ${w.get_pcg_step_count()}: ` +
    `last_iters=${w.get_pcg_last_iters()}, ` +
    `div_norm=${w.get_pcg_divergence_norm().toExponential(2)}`
  );
}
```

---

## Error Handling

| Situation | Behaviour |
|---|---|
| `get_*_at(i, j)` with OOB index | Returns `0.0` — no exception |
| `run_n_steps(0)` | No-op |
| Running without calling `init_*` first | Fields contain zeros (zero-initialised by MoonBit) |
| PCG fails to converge in 200 iterations | Solver returns best current solution; `get_*_last_iters()` = 200 |
