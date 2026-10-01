/**
 * FlowLabLite — Chorin + Gauss-Seidel benchmark  (Java)
 *
 * Matches FlowLabLite v0.0.1 Chorin solver (tag v0.0.1,
 * cmd/main/main.mbt cavity_flow_array).
 *
 * Physical parameters: domain 2×2, rho=1.0, nu=0.1 (Re=20), dt=0.001
 * BCs: top lid u=1, all walls no-slip
 *
 * Build & run:
 *   javac ChorinGSBench.java
 *   java -server ChorinGSBench            # default: all grid sizes
 *   java -server ChorinGSBench 41 500 5   # n nt repeats
 */
public class ChorinGSBench {

    static final double RHO = 1.0;
    static final double NU  = 0.1;
    static final double DT  = 0.001;
    static final int    NIT = 50;
    static final double L   = 2.0;

    // ── Core solver ─────────────────────────────────────────────────────────

    static void solve(int n, int nt) {
        double dx  = L / (n - 1);
        double dy  = L / (n - 1);
        double dx2 = dx * dx;
        double dy2 = dy * dy;
        double invD = 1.0 / (2.0 * (dx2 + dy2));
        double coef = dx2 * dy2 * invD;

        double[][] u  = new double[n][n];
        double[][] v  = new double[n][n];
        double[][] p  = new double[n][n];
        double[][] un = new double[n][n];
        double[][] vn = new double[n][n];
        double[][] b  = new double[n][n];
        double[][] pn = new double[n][n];

        for (int step = 0; step < nt; step++) {

            // Save old velocities
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++) { un[i][j] = u[i][j]; vn[i][j] = v[i][j]; }

            // Build RHS b
            for (int i = 1; i < n - 1; i++)
                for (int j = 1; j < n - 1; j++)
                    b[i][j] = RHO / DT * (
                        (un[i][j+1] - un[i][j-1]) / (2*dx) +
                        (vn[i+1][j] - vn[i-1][j]) / (2*dy)
                    );

            // Gauss-Seidel pressure Poisson
            for (int q = 0; q < NIT; q++) {
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) pn[i][j] = p[i][j];

                for (int i = 1; i < n - 1; i++)
                    for (int j = 1; j < n - 1; j++) {
                        double lap = (pn[i][j+1] + pn[i][j-1]) * dy2
                                   + (pn[i+1][j] + pn[i-1][j]) * dx2;
                        p[i][j] = lap * invD - coef * b[i][j];
                    }
                // BCs
                for (int i = 0; i < n; i++) p[i][n-1] = p[i][n-2];
                for (int j = 0; j < n; j++) p[0][j]   = p[1][j];
                for (int i = 0; i < n; i++) p[i][0]   = p[i][1];
                for (int j = 0; j < n; j++) p[n-1][j] = 0.0;
            }

            // Velocity update
            for (int i = 1; i < n - 1; i++)
                for (int j = 1; j < n - 1; j++) {
                    double cux = un[i][j] * DT/dx * (un[i][j] - un[i][j-1]);
                    double cuy = vn[i][j] * DT/dy * (un[i][j] - un[i-1][j]);
                    double pgx = DT / (2*RHO*dx) * (p[i][j+1] - p[i][j-1]);
                    double vux = NU * DT/dx2 * (un[i][j+1] - 2*un[i][j] + un[i][j-1]);
                    double vuy = NU * DT/dy2 * (un[i+1][j] - 2*un[i][j] + un[i-1][j]);
                    u[i][j] = un[i][j] - cux - cuy - pgx + vux + vuy;

                    double cvx = un[i][j] * DT/dx * (vn[i][j] - vn[i][j-1]);
                    double cvy = vn[i][j] * DT/dy * (vn[i][j] - vn[i-1][j]);
                    double pgy = DT / (2*RHO*dy) * (p[i+1][j] - p[i-1][j]);
                    double vvx = NU * DT/dx2 * (vn[i][j+1] - 2*vn[i][j] + vn[i][j-1]);
                    double vvy = NU * DT/dy2 * (vn[i+1][j] - 2*vn[i][j] + vn[i-1][j]);
                    v[i][j] = vn[i][j] - cvx - cvy - pgy + vvx + vvy;
                }

            // BCs
            for (int j = 0; j < n; j++) u[0][j] = 0.0;
            for (int i = 0; i < n; i++) { u[i][0] = 0.0; u[i][n-1] = 0.0; }
            for (int j = 0; j < n; j++) u[n-1][j] = 1.0;
            for (int j = 0; j < n; j++) { v[0][j] = 0.0; v[n-1][j] = 0.0; }
            for (int i = 0; i < n; i++) { v[i][0] = 0.0; v[i][n-1] = 0.0; }
        }
        // Prevent JIT from eliminating the computation
        if (Double.isNaN(u[n/2][n/2])) System.err.println("NaN detected");
    }

    // ── Benchmark runner ─────────────────────────────────────────────────────

    record GridCase(int n, int nt, int repeats, String label) {}

    static final GridCase[] CASES = {
        new GridCase(41,  500, 5, "small  (41×41)"),
        new GridCase(81,  500, 3, "medium (81×81)"),
        new GridCase(161, 500, 3, "large (161×161)"),
    };

    public static void main(String[] args) {
        System.out.println("=".repeat(70));
        System.out.println("FlowLabLite — Chorin+GS Java benchmark (JDK " +
            System.getProperty("java.version") + ")");
        System.out.println("=".repeat(70));

        GridCase[] cases;
        if (args.length >= 2) {
            int n  = Integer.parseInt(args[0]);
            int nt = Integer.parseInt(args[1]);
            int rp = args.length >= 3 ? Integer.parseInt(args[2]) : 3;
            cases  = new GridCase[]{ new GridCase(n, nt, rp, n + "×" + n) };
        } else {
            cases = CASES;
        }

        // JVM warmup (run smallest case twice before timing)
        System.out.println("JVM warmup...");
        solve(cases[0].n(), cases[0].nt());
        solve(cases[0].n(), cases[0].nt());
        System.out.println();

        System.out.println("=== TSV_RESULTS ===");
        System.out.println("grid\tnt\tjava_ms");

        for (GridCase c : cases) {
            double[] times = new double[c.repeats()];
            for (int r = 0; r < c.repeats(); r++) {
                long t0 = System.nanoTime();
                solve(c.n(), c.nt());
                long t1 = System.nanoTime();
                times[r] = (t1 - t0) / 1_000_000.0;
            }
            double best = Double.MAX_VALUE, sum = 0;
            for (double t : times) { if (t < best) best = t; sum += t; }
            double avg = sum / times.length;
            System.out.printf("  %-22s  n=%4d  nt=%5d  best=%9.1f ms  avg=%9.1f ms%n",
                c.label(), c.n(), c.nt(), best, avg);
            System.out.printf("%s\t%d\t%.1f%n", c.n() + "×" + c.n(), c.nt(), best);
        }
    }
}
