/*
 * FlowLabLite — Chorin + Gauss-Seidel benchmark  (C)
 *
 * Matches FlowLabLite v0.0.1 Chorin solver (tag v0.0.1,
 * cmd/main/main.mbt cavity_flow_array).
 *
 * Compile & run:
 *   gcc -O2 -o chorin_gs_bench chorin_gs_bench.c -lm
 *   ./chorin_gs_bench              # all grid sizes
 *   ./chorin_gs_bench 41 500 5    # n nt repeats
 *
 * Physical parameters: domain 2×2, rho=1.0, nu=0.1 (Re=20), dt=0.001
 * BCs: top lid u=1, all walls no-slip
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define RHO  1.0
#define NU   0.1
#define DT   0.001
#define NIT  50
#define L    2.0

/* Wall-clock timer in milliseconds */
static double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec * 1000.0 + ts.tv_nsec / 1.0e6;
}

/* Allocate n×n double array (flat row-major) */
static double* alloc2d(int n) {
    double* a = (double*)calloc(n * n, sizeof(double));
    if (!a) { fprintf(stderr, "OOM\n"); exit(1); }
    return a;
}

#define IDX(n, i, j)  ((i)*(n) + (j))
#define U(i,j)  u[IDX(n,i,j)]
#define V(i,j)  v[IDX(n,i,j)]
#define P(i,j)  p[IDX(n,i,j)]
#define UN(i,j) un[IDX(n,i,j)]
#define VN(i,j) vn[IDX(n,i,j)]
#define PN(i,j) pn[IDX(n,i,j)]
#define B(i,j)  b[IDX(n,i,j)]

static void solve(int n, int nt) {
    double dx  = L / (n - 1);
    double dy  = L / (n - 1);
    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double invD = 1.0 / (2.0 * (dx2 + dy2));
    double coef = dx2 * dy2 * invD;

    double *u  = alloc2d(n), *v  = alloc2d(n), *p  = alloc2d(n);
    double *un = alloc2d(n), *vn = alloc2d(n), *pn = alloc2d(n);
    double *b  = alloc2d(n);

    for (int step = 0; step < nt; step++) {

        /* Save old velocities */
        memcpy(un, u, n*n*sizeof(double));
        memcpy(vn, v, n*n*sizeof(double));

        /* Build RHS b */
        for (int i = 1; i < n-1; i++)
            for (int j = 1; j < n-1; j++)
                B(i,j) = RHO/DT * (
                    (UN(i,j+1) - UN(i,j-1)) / (2*dx) +
                    (VN(i+1,j) - VN(i-1,j)) / (2*dy)
                );

        /* Gauss-Seidel pressure Poisson */
        for (int q = 0; q < NIT; q++) {
            memcpy(pn, p, n*n*sizeof(double));
            for (int i = 1; i < n-1; i++)
                for (int j = 1; j < n-1; j++) {
                    double lap = (PN(i,j+1)+PN(i,j-1))*dy2
                               + (PN(i+1,j)+PN(i-1,j))*dx2;
                    P(i,j) = lap * invD - coef * B(i,j);
                }
            /* BCs */
            for (int i = 0; i < n; i++) P(i,n-1) = P(i,n-2);
            for (int j = 0; j < n; j++) P(0,j)   = P(1,j);
            for (int i = 0; i < n; i++) P(i,0)   = P(i,1);
            for (int j = 0; j < n; j++) P(n-1,j) = 0.0;
        }

        /* Velocity update */
        for (int i = 1; i < n-1; i++)
            for (int j = 1; j < n-1; j++) {
                double cux = UN(i,j)*DT/dx*(UN(i,j)-UN(i,j-1));
                double cuy = VN(i,j)*DT/dy*(UN(i,j)-UN(i-1,j));
                double pgx = DT/(2*RHO*dx)*(P(i,j+1)-P(i,j-1));
                double vux = NU*DT/dx2*(UN(i,j+1)-2*UN(i,j)+UN(i,j-1));
                double vuy = NU*DT/dy2*(UN(i+1,j)-2*UN(i,j)+UN(i-1,j));
                U(i,j) = UN(i,j) - cux - cuy - pgx + vux + vuy;

                double cvx = UN(i,j)*DT/dx*(VN(i,j)-VN(i,j-1));
                double cvy = VN(i,j)*DT/dy*(VN(i,j)-VN(i-1,j));
                double pgy = DT/(2*RHO*dy)*(P(i+1,j)-P(i-1,j));
                double vvx = NU*DT/dx2*(VN(i,j+1)-2*VN(i,j)+VN(i,j-1));
                double vvy = NU*DT/dy2*(VN(i+1,j)-2*VN(i,j)+VN(i-1,j));
                V(i,j) = VN(i,j) - cvx - cvy - pgy + vvx + vvy;
            }

        /* BCs */
        for (int j = 0; j < n; j++) U(0,j) = 0.0;
        for (int i = 0; i < n; i++) { U(i,0) = 0.0; U(i,n-1) = 0.0; }
        for (int j = 0; j < n; j++) U(n-1,j) = 1.0;
        for (int j = 0; j < n; j++) { V(0,j) = 0.0; V(n-1,j) = 0.0; }
        for (int i = 0; i < n; i++) { V(i,0) = 0.0; V(i,n-1) = 0.0; }
    }

    /* Prevent optimiser from eliminating computation */
    if (isnan(U(n/2, n/2))) fprintf(stderr, "NaN\n");

    free(u); free(v); free(p);
    free(un); free(vn); free(pn); free(b);
}

int main(int argc, char *argv[]) {
    typedef struct { int n, nt, reps; const char *label; } Case;
    Case cases[] = {
        { 41,  500, 5, "small  (41×41)"  },
        { 81,  500, 3, "medium (81×81)"  },
        { 161, 500, 3, "large (161×161)" },
    };
    int ncases = 3;

    if (argc >= 3) {
        cases[0].n    = atoi(argv[1]);
        cases[0].nt   = atoi(argv[2]);
        cases[0].reps = argc >= 4 ? atoi(argv[3]) : 3;
        cases[0].label = "custom";
        ncases = 1;
    }

    printf("%-70s\n", "FlowLabLite — Chorin+GS C benchmark (gcc -O2)");
    printf("%-70s\n\n", "=");

    printf("=== TSV_RESULTS ===\n");
    printf("grid\tnt\tc_ms\n");

    for (int ci = 0; ci < ncases; ci++) {
        Case *c = &cases[ci];
        double best = 1e18, total = 0;
        for (int r = 0; r < c->reps; r++) {
            double t0 = now_ms();
            solve(c->n, c->nt);
            double t1 = now_ms();
            double elapsed = t1 - t0;
            if (elapsed < best) best = elapsed;
            total += elapsed;
        }
        double avg = total / c->reps;
        printf("  %-22s  n=%4d  nt=%5d  best=%9.1f ms  avg=%9.1f ms\n",
               c->label, c->n, c->nt, best, avg);
        printf("%d×%d\t%d\t%.1f\n", c->n, c->n, c->nt, best);
    }
    return 0;
}
