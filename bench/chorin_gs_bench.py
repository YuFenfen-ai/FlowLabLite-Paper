#!/usr/bin/env python3
"""
FlowLabLite — Chorin + Gauss-Seidel benchmark  (Python)

Benchmarks two implementations:
  1. numpy  — vectorised pressure Poisson and velocity update
  2. pure   — plain Python loops (no numpy)

Usage:
  python chorin_gs_bench.py              # default: small grid, 3 repeats
  python chorin_gs_bench.py 81 500 3    # grid_n steps repeats
  python chorin_gs_bench.py --all        # run all grid sizes

Reference implementation matches FlowLabLite v0.0.1 Chorin solver
(tag v0.0.1, cmd/main/main.mbt cavity_flow_array).

Physical parameters:
  domain  2×2,  rho=1.0,  nu=0.1 (Re=20),  dt=0.001
  BCs: top lid u=1, all walls no-slip
"""

import sys
import time
import math

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers
# ─────────────────────────────────────────────────────────────────────────────

def wall_clock():
    return time.perf_counter()

def make_grid(n):
    """Return n×n zero arrays for u, v, p."""
    try:
        import numpy as np
        return np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    except ImportError:
        z = [[0.0] * n for _ in range(n)]
        return ([row[:] for row in z],
                [row[:] for row in z],
                [row[:] for row in z])

# ─────────────────────────────────────────────────────────────────────────────
# NumPy implementation
# ─────────────────────────────────────────────────────────────────────────────

def run_numpy(grid_n, nt, rho=1.0, nu=0.1, dt=0.001, nit=50):
    import numpy as np
    L   = 2.0
    dx  = L / (grid_n - 1)
    dy  = L / (grid_n - 1)
    dx2 = dx * dx
    dy2 = dy * dy
    coeff     = dx2 * dy2 / (2.0 * (dx2 + dy2))
    inv_denom = 1.0 / (2.0 * (dx2 + dy2))

    u = np.zeros((grid_n, grid_n))
    v = np.zeros((grid_n, grid_n))
    p = np.zeros((grid_n, grid_n))

    for _ in range(nt):
        un = u.copy()
        vn = v.copy()

        # ── Build RHS b ───────────────────────────────────────────────────
        b = np.zeros((grid_n, grid_n))
        b[1:-1, 1:-1] = (rho * (
            (un[1:-1, 2:] - un[1:-1, :-2]) / (2*dx) +
            (vn[2:, 1:-1] - vn[:-2, 1:-1]) / (2*dy)
        ) / dt)

        # ── Gauss-Seidel pressure Poisson ─────────────────────────────────
        for _ in range(nit):
            pn = p.copy()
            p[1:-1, 1:-1] = (
                (pn[1:-1, 2:] + pn[1:-1, :-2]) * dy2 +
                (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx2
            ) * inv_denom - coeff * b[1:-1, 1:-1]
            # BCs
            p[:, -1] = p[:, -2]
            p[0, :]  = p[1, :]
            p[:, 0]  = p[:, 1]
            p[-1, :] = 0.0

        # ── Velocity update ───────────────────────────────────────────────
        u[1:-1, 1:-1] = (
            un[1:-1, 1:-1]
            - un[1:-1, 1:-1] * dt/dx * (un[1:-1, 1:-1] - un[1:-1, :-2])
            - vn[1:-1, 1:-1] * dt/dy * (un[1:-1, 1:-1] - un[:-2,  1:-1])
            - dt / (2*rho*dx) * (p[1:-1, 2:] - p[1:-1, :-2])
            + nu * dt/dx2 * (un[1:-1, 2:] - 2*un[1:-1, 1:-1] + un[1:-1, :-2])
            + nu * dt/dy2 * (un[2:,  1:-1] - 2*un[1:-1, 1:-1] + un[:-2,  1:-1])
        )
        v[1:-1, 1:-1] = (
            vn[1:-1, 1:-1]
            - un[1:-1, 1:-1] * dt/dx * (vn[1:-1, 1:-1] - vn[1:-1, :-2])
            - vn[1:-1, 1:-1] * dt/dy * (vn[1:-1, 1:-1] - vn[:-2,  1:-1])
            - dt / (2*rho*dy) * (p[2:,  1:-1] - p[:-2, 1:-1])
            + nu * dt/dx2 * (vn[1:-1, 2:] - 2*vn[1:-1, 1:-1] + vn[1:-1, :-2])
            + nu * dt/dy2 * (vn[2:,  1:-1] - 2*vn[1:-1, 1:-1] + vn[:-2,  1:-1])
        )
        # BCs
        u[0, :] = 0.0;  u[:, 0] = 0.0;  u[:, -1] = 0.0
        u[-1, :] = 1.0
        v[0, :] = 0.0;  v[-1, :] = 0.0
        v[:, 0] = 0.0;  v[:, -1] = 0.0

    return u, v, p


# ─────────────────────────────────────────────────────────────────────────────
# Pure-Python implementation (no numpy)
# ─────────────────────────────────────────────────────────────────────────────

def run_pure(grid_n, nt, rho=1.0, nu=0.1, dt=0.001, nit=50):
    L  = 2.0
    dx = L / (grid_n - 1)
    dy = L / (grid_n - 1)
    dx2 = dx * dx
    dy2 = dy * dy
    inv_denom = 1.0 / (2.0 * (dx2 + dy2))
    coeff     = dx2 * dy2 * inv_denom

    u = [[0.0] * grid_n for _ in range(grid_n)]
    v = [[0.0] * grid_n for _ in range(grid_n)]
    p = [[0.0] * grid_n for _ in range(grid_n)]

    for _ in range(nt):
        # Save old velocities
        un = [row[:] for row in u]
        vn = [row[:] for row in v]

        # Build b
        b = [[0.0] * grid_n for _ in range(grid_n)]
        for i in range(1, grid_n - 1):
            for j in range(1, grid_n - 1):
                b[i][j] = rho / dt * (
                    (un[i][j+1] - un[i][j-1]) / (2*dx) +
                    (vn[i+1][j] - vn[i-1][j]) / (2*dy)
                )

        # Gauss-Seidel
        for _ in range(nit):
            pn = [row[:] for row in p]
            for i in range(1, grid_n - 1):
                for j in range(1, grid_n - 1):
                    lap = (pn[i][j+1] + pn[i][j-1]) * dy2 + (pn[i+1][j] + pn[i-1][j]) * dx2
                    p[i][j] = lap * inv_denom - coeff * b[i][j]
            for i in range(grid_n): p[i][-1] = p[i][-2]
            for j in range(grid_n): p[0][j]  = p[1][j]
            for i in range(grid_n): p[i][0]  = p[i][1]
            for j in range(grid_n): p[-1][j] = 0.0

        # Velocity update
        for i in range(1, grid_n - 1):
            for j in range(1, grid_n - 1):
                cux = un[i][j] * dt/dx * (un[i][j] - un[i][j-1])
                cuy = vn[i][j] * dt/dy * (un[i][j] - un[i-1][j])
                pgx = dt / (2*rho*dx) * (p[i][j+1] - p[i][j-1])
                vux = nu * dt/dx2 * (un[i][j+1] - 2*un[i][j] + un[i][j-1])
                vuy = nu * dt/dy2 * (un[i+1][j] - 2*un[i][j] + un[i-1][j])
                u[i][j] = un[i][j] - cux - cuy - pgx + vux + vuy

                cvx = un[i][j] * dt/dx * (vn[i][j] - vn[i][j-1])
                cvy = vn[i][j] * dt/dy * (vn[i][j] - vn[i-1][j])
                pgy = dt / (2*rho*dy) * (p[i+1][j] - p[i-1][j])
                vvx = nu * dt/dx2 * (vn[i][j+1] - 2*vn[i][j] + vn[i][j-1])
                vvy = nu * dt/dy2 * (vn[i+1][j] - 2*vn[i][j] + vn[i-1][j])
                v[i][j] = vn[i][j] - cvx - cvy - pgy + vvx + vvy

        # BCs
        for j in range(grid_n): u[0][j] = 0.0
        for i in range(grid_n): u[i][0] = 0.0; u[i][-1] = 0.0
        for j in range(grid_n): u[-1][j] = 1.0
        for j in range(grid_n): v[0][j] = 0.0; v[-1][j] = 0.0
        for i in range(grid_n): v[i][0] = 0.0; v[i][-1] = 0.0

    return u, v, p


# ─────────────────────────────────────────────────────────────────────────────
# Benchmark runner
# ─────────────────────────────────────────────────────────────────────────────

def bench_one(label, fn, grid_n, nt, repeats):
    times = []
    for r in range(repeats):
        t0 = wall_clock()
        fn(grid_n, nt)
        t1 = wall_clock()
        times.append((t1 - t0) * 1000)
    best = min(times)
    avg  = sum(times) / len(times)
    print(f"  {label:20s}  n={grid_n:4d}  nt={nt:5d}  "
          f"best={best:9.1f} ms  avg={avg:9.1f} ms  ({repeats} runs)")
    return best

SIZES = [
    (41,  500, 5, "small"),
    (81,  500, 3, "medium"),
    (161, 500, 3, "large"),
]

def main():
    args = sys.argv[1:]
    run_all = "--all" in args

    print("=" * 70)
    print("FlowLabLite — Chorin+GS Python benchmark")
    print("=" * 70)

    # Detect numpy
    try:
        import numpy as np
        has_numpy = True
        print(f"NumPy {np.__version__} available")
    except ImportError:
        has_numpy = False
        print("NumPy not available — running pure-Python only")
    print()

    if run_all:
        sizes = SIZES
    elif len(args) >= 1 and args[0] != "--all":
        n  = int(args[0]) if len(args) > 0 else 41
        nt = int(args[1]) if len(args) > 1 else 500
        rp = int(args[2]) if len(args) > 2 else 3
        sizes = [(n, nt, rp, f"{n}×{n}")]
    else:
        sizes = SIZES

    results = []
    for grid_n, nt, repeats, label in sizes:
        print(f"Grid {grid_n}×{grid_n}, {nt} steps:")
        row = {"grid": f"{grid_n}×{grid_n}", "nt": nt}
        if has_numpy:
            row["numpy_ms"] = bench_one("numpy", run_numpy, grid_n, nt, repeats)
        row["pure_ms"]  = bench_one("pure-Python", run_pure, grid_n, nt, repeats)
        results.append(row)
        print()

    # TSV output for run_bench.sh to capture
    print("=== TSV_RESULTS ===")
    header = "grid\tnt\tnumpy_ms\tpure_ms"
    print(header)
    for r in results:
        np_ms = f"{r.get('numpy_ms', -1):.1f}"
        pu_ms = f"{r.get('pure_ms', -1):.1f}"
        print(f"{r['grid']}\t{r['nt']}\t{np_ms}\t{pu_ms}")

if __name__ == "__main__":
    main()
