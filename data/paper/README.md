# paper_data — data package for "Theory-Code-Test Co-Design for AI-Assisted Scientific Software"

Every number in paper sections 4.3–4.6 traces to a file in this folder.
Run `python verify_paper_numbers.py` (exit code 0 = all checks pass); the recorded run is
`verification_report.txt` (74/74 checks passed, 2026-09-27).

## Map: paper item → data file

| Paper item | Data file(s) here | Origin |
|---|---|---|
| Table 3 (baseline centerline errors, Re=100 row) | `sec4_3/ghia_re100_g129_centerline_s40000.csv`, `sec4_4_grid_convergence/g129/time_convergence.json` | FlowLabLite validation harness (`validation/run_ghia_study.py`, Chorin-GS, 129×129, dt=0.001, 40,000 steps) |
| Table 4 (grid convergence) + §4.4 text (div norms, u_min, checkpoint series, profile changes, diffusion numbers) | `sec4_4_grid_convergence/g{65,129,257}/time_convergence.json` (+ `time_convergence.csv`, `pointwise_errors.csv`, `centerline_s40000.csv`, `run_manifest.json` per grid) | same harness; 65²/129² at dt=0.001, 257² at dt=0.0005 (explicit-diffusion stability limit), 40,000 steps each |
| Figure 8 (log–log error vs grid spacing) | `sec4_4_grid_convergence/fig8_grid_convergence.png`; data = the three `time_convergence.json` files | matplotlib; regenerated from the JSONs by the audit script's data |
| Table 5 (flow-structure validation; FlowLabLite column) | `sec4_5_flow_structure/field_re100_g129.csv` (full u,v,p field, 129×129, Re=100, 40,000 steps); ψ by trapezoidal integration of u upward from the bottom wall; primary vortex = argmin ψ; BL1/BR1 = local ψ maxima in the bottom corners | `moon run cmd/main --release --target wasm -- --format csv --solver chorin --grid 129 --re 100 --steps 40000` (header label bug fixed 2026-09-30, see caveat 1) |
| Table 5 (Ghia column) | `reference_values/ghia1982_tableV_re100.txt` | transcribed from Ghia, Ghia & Shin (1982) Table V (printed pp. 408–409) and Tables I/II (underscored minima) |
| Table 5 (Marchi column) | `reference_values/marchi2009_re100.txt` | transcribed from Marchi, Suero & Araki (2009), Tables 8/11/12 "Present" rows |
| Table 6 (Boussinesq natural convection; FlowLabLite column) | `sec4_6_natural_convection/bous_ra1e{3,4,5}.csv` (final plateau rows; u*=20u, v*=20v) | Boussinesq module CLI (see `sec4_6_natural_convection/commands.txt`): 41×41 built-in grid, dt=0.001, ν=0.071, α=0.1, β=0.8875·Ra/10³ (Ra = β·L³/(να); Pr=0.71) |
| Table 6 (benchmark column) | `reference_values/devahldavis1983_benchmark.txt` | De Vahl Davis (1983) benchmark solutions |
| §4.4 text: 257² divergence at dt=0.001 | not a data file — computed: νΔt/Δx² = 0.02·0.001/(2/256)² = 0.32768 ("0.33" in text); 0.16 at dt=0.0005 | analytic stability argument |

## Known caveats recorded in the paper

1. RESOLVED: the export-header label bug (header printed "Re=20" for Re=100 runs because
   `io_formats.mbt` read the compile-time `nu` instead of the runtime `g_nu_runtime[0]`) was fixed in
   all four export headers (CSV/VTK/Tecplot/VTI) and `field_re100_g129.csv` was regenerated. The new
   header reads Re=100 and all 16,643 data lines are byte-identical to the originally archived file,
   so every published number is unchanged. Regression-safe: without `--re` the header still shows the
   default configuration.
2. All FlowLabLite §4.5 values sit at the 40,000-step horizon, not at steady state; the ~10% vortex-strength
   deficit against both steady references is reported and explained in §4.4/§4.5.
3. Ra=10⁵ natural-convection case is reported as a documented resolution limit (thermal boundary layer
   spans 1–2 grid cells), not a benchmark passage.
4. De Vahl Davis Ra=10⁵ u_max 34.81 is the original paper's value; later benchmarks give ≈34.74.

## Reproduction environment

- Solver: FlowLabLite worktree `high-re-projection-tvd` (D:\yufenfen\MoonBit\project\FlowLabLite\.worktrees\...),
  MoonBit `moon` CLI, target wasm; 176 tests passing at data-generation time.
- Grid-convergence runs: `python -m validation.run_ghia_study ...` (see each `run_manifest.json` for the
  exact command, platform, Python version, and SHA-256 of the centerline files).
- Audit: `python verify_paper_numbers.py` (pure Python 3 standard library).
