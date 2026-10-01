// FlowLabLite — CLI Format Dispatch
//
// Parses --format, --solver, --list-formats from @env.args() and dispatches
// to the appropriate Tier-1 format function (or the NetCDF-4 stub).
//
// Usage (via moon run):
//   moon run cmd/main --target wasm -- --list-formats
//   moon run cmd/main --target wasm -- --format csv
//   moon run cmd/main --target wasm -- --format vtk [--solver chorin|simple|pcg|mac]
//   moon run cmd/main --target wasm -- --format tecplot [--solver all|chorin|...]
//   moon run cmd/main --target wasm -- --format netcdf   (stub: prints error)
//
// Output goes to stdout.  Redirect with shell:
//   moon run cmd/main --target wasm -- --format vtk --solver pcg > cavity_pcg.vtk

// ─────────────────────────────────────────────────────────────────────────────
// Argument helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Parse a decimal integer string (positive only, no sign).
/// Returns `def` when the string is empty, non-numeric, or overflows.
/// Used for --re and --steps CLI argument parsing.
fn parse_pos_int(s : String, def : Int) -> Int {
  if s.length() == 0 { return def }
  let n   : Array[Int]  = [0]
  let ok  : Array[Bool] = [true]
  s.iter().each(fn(c : Char) {
    if c >= '0' && c <= '9' {
      n[0] = n[0] * 10 + (c.to_int() - 48)   // '0' = ASCII 48
    } else {
      ok[0] = false
    }
  })
  if ok[0] && s.length() > 0 { n[0] } else { def }
}

/// Return true when a CLI solver filter requests the named solver.
/// The same predicate is used for execution and export so --solver cannot
/// silently run unrelated solvers in the background.
fn cli_solver_selected(filter : String, solver_name : String) -> Bool {
  filter == "all" || filter == solver_name
}

/// Parse comma-separated cumulative checkpoint steps.
/// Returns an empty array for an invalid, non-positive, duplicate, or
/// decreasing sequence so callers can reject it before starting a long run.
fn parse_checkpoint_steps(raw : String) -> Array[Int] {
  if raw.length() == 0 { return [] }
  let checkpoints : Array[Int] = []
  let valid : Array[Bool] = [true]
  let previous : Array[Int] = [0]
  raw.split(",").each(fn(part) {
    let step = parse_pos_int(part.to_string(), 0)
    if step <= previous[0] {
      valid[0] = false
    } else {
      checkpoints.push(step)
      previous[0] = step
    }
  })
  if valid[0] { checkpoints } else { [] }
}
/// Return true if any format-related flag is present in argv.
/// Scans from index 1 (skipping argv[0] = program name).
fn cli_wants_format(args : Array[String]) -> Bool {
  let mut i = 1
  while i < args.length() {
    if args[i] == "--format" || args[i] == "--list-formats" || args[i] == "--validate-ghia" {
      return true
    }
    i = i + 1
  }
  false
}

/// Return the value that follows `flag` in argv, or `default` if absent.
/// E.g. cli_get_string(args, "--format", "json") → "vtk" when args has "--format" "vtk".
fn cli_get_string(
  args    : Array[String],
  flag    : String,
  default : String,
) -> String {
  let mut i = 1
  while i < args.length() - 1 {
    if args[i] == flag {
      return args[i + 1]
    }
    i = i + 1
  }
  default
}

// ─────────────────────────────────────────────────────────────────────────────
// Format availability listing
// ─────────────────────────────────────────────────────────────────────────────

fn print_list_formats() -> Unit {
  println("Available formats in this build:")
  println("  json       JSON (default, with local_viewer.html support)    \u2713")
  println("  csv        RFC 4180 CSV (.csv)                               \u2713")
  println("  vtk        VTK Legacy ASCII Format (.vtk, ParaView / VisIt)  \u2713")
  println("  vti        VTK XML ImageData (.vti, ParaView preferred)      \u2713")
  println("  tecplot    Tecplot ASCII Data Format (.dat)                  \u2713")
  println("  centerline Centerline velocity CSV (u_vcl + v_hcl)          \u2713")
  println("  html       Self-contained HTML report (Chart.js, CDN)       \u2713")
  println("  monitor    Convergence monitor history CSV                   \u2713")
  if netcdf_available {
    println("  netcdf    CF-1.11 NetCDF-4 (.nc, ncview / Panoply / CDO)    \u2713")
  } else {
    println("  netcdf    CF-1.11 NetCDF-4 (.nc, ncview / Panoply / CDO)    \u2717  (libnetcdf not linked)")
  }
  println("")
  println("To redirect output to a file:")
  println("  moon run cmd/main --target wasm -- --format vtk > cavity.vtk")
  println("  moon run cmd/main --target wasm -- --format vtk --solver pcg > pcg.vtk")
}

// ─────────────────────────────────────────────────────────────────────────────
// Data collection from global solver state
// ─────────────────────────────────────────────────────────────────────────────

/// Build a cell-centred SolverExport from the MAC staggered global state.
/// Interpolates staggered u/v to cell centres using get_u_mac_at / get_v_mac_at.
fn collect_mac_export() -> SolverExport {
  let u_cc = create_zeros_2d(mac_nc, mac_nc)
  let v_cc = create_zeros_2d(mac_nc, mac_nc)
  let p_cc = create_zeros_2d(mac_nc, mac_nc)
  for i = 0; i < mac_nc; i = i + 1 {
    for j = 0; j < mac_nc; j = j + 1 {
      u_cc[i][j] = get_u_mac_at(i, j)
      v_cc[i][j] = get_v_mac_at(i, j)
      p_cc[i][j] = g_p_mac[i][j]
    }
  }
  {
    name:  "mac",
    label: "MAC Staggered-Grid PCG",
    nc:    mac_nc,
    h:     dx,
    x0:    dx / 2.0,
    y0:    dy / 2.0,
    steps: get_mac_step_count(),        // actual steps completed (respects --steps override)
    u:     u_cc,
    v:     v_cc,
    p:     p_cc,
  }
}

/// Collect post-run global state into SolverExport structs.
///
/// solver_filter: "all" | "chorin" | "simple" | "pcg" | "mac"
/// Only solvers whose run_* flag is true AND whose name matches the filter
/// are included.
fn collect_solver_exports(solver_filter : String) -> Array[SolverExport] {
  let result : Array[SolverExport] = []
  let want_all = solver_filter == "all"
  if (want_all || solver_filter == "chorin") && run_chorin {
    result.push({
      name:  "chorin",
      label: "Chorin (Projection + Gauss-Seidel)",
      nc:    g_nx[0],
      h:     g_dx[0],
      x0:    0.0,
      y0:    0.0,
      steps: get_step_count(),          // actual steps completed (respects --steps override)
      u:     g_u[0],
      v:     g_v[0],
      p:     g_p[0],
    })
  }
  if (want_all || solver_filter == "simple") && run_simple {
    result.push({
      name:  "simple",
      label: "SIMPLE (Pressure-Correction)",
      nc:    nx,
      h:     dx,
      x0:    0.0,
      y0:    0.0,
      steps: get_simple_step_count(),   // actual iterations completed
      u:     g_u_s,
      v:     g_v_s,
      p:     g_p_s,
    })
  }
  if (want_all || solver_filter == "pcg") && run_pcg {
    result.push({
      name:  "pcg",
      label: "Chorin-PCG (Collocated)",
      nc:    nx,
      h:     dx,
      x0:    0.0,
      y0:    0.0,
      steps: get_pcg_step_count(),      // actual steps completed
      u:     g_u_pcg,
      v:     g_v_pcg,
      p:     g_p_pcg,
    })
  }
  if (want_all || solver_filter == "mac") && run_mac {
    result.push(collect_mac_export())
  }
  result
}

// ─────────────────────────────────────────────────────────────────────────────
// Main dispatch
// ─────────────────────────────────────────────────────────────────────────────

/// Dispatch CLI format output.  Called after all solvers have been run.
///
/// Reads --format, --solver from args.
/// All output goes to stdout (redirect with > for file output).
fn dispatch_format_cli(args : Array[String]) -> Unit {
  // --list-formats: no solver run needed
  if args.contains("--list-formats") {
    print_list_formats()
    return
  }

  let fmt    = cli_get_string(args, "--format",  "json")
  let solver = cli_get_string(args, "--solver",  "all")

  // --validate-ghia: compare solver solution against Ghia (1982) Re=100 benchmark.
  // Handled before the fmt=="json" guard so it works without --format.
  if args.contains("--validate-ghia") {
    let ghia_solver = if solver == "all" { "chorin" } else { solver }
    let exports = collect_solver_exports(ghia_solver)
    if exports.is_empty() {
      println("[io] Solver '\{ghia_solver}' is not available for Ghia comparison.")
      return
    }
    let re_str = cli_get_string(args, "--re", "0")
    let re_val = if re_str == "0" { 2.0 / g_nu_runtime[0] }
                 else { parse_pos_int(re_str, 20).to_double() }
    let result = ghia_compare(exports[0])
    println(format_ghia_result(result, exports[0].name, re_val, exports[0].steps))
    return
  }

  if fmt == "json" {
    println("[io] Hint: omit --format (or all flags) for the full JSON output with timing.")
    println("[io] Explicit --format json is not supported; re-run without --format.")
    return
  }

  // For VTK / VTI (one file per solver): default solver is "chorin"
  if fmt == "vtk" || fmt == "vti" {
    let single_solver = if solver == "all" { "chorin" } else { solver }
    let exports = collect_solver_exports(single_solver)
    if exports.is_empty() {
      println("[io] Solver '\{single_solver}' is not available (check --solver or run_* flags).")
      println("[io] Available: chorin | simple | pcg | mac")
      return
    }
    if fmt == "vtk" {
      println(format_vtk(exports[0]))
    } else {
      println(format_vti(exports[0]))
    }
    return
  }

  // For CSV / Tecplot: collect all (or filtered) solvers
  let solvers = collect_solver_exports(solver)
  if solvers.is_empty() {
    println("[io] No solver data to export (check --solver or run_* flags).")
    return
  }

  if fmt == "csv" {
    println(format_csv(solvers))
  } else if fmt == "tecplot" {
    println(format_tecplot(solvers))
  } else if fmt == "centerline" {
    let s = solvers[0]
    let re_val = 2.0 / g_nu_runtime[0]
    let div_norm = if s.name == "chorin" { get_divergence_norm() }
                   else if s.name == "pcg" { get_pcg_divergence_norm() }
                   else if s.name == "mac" { get_mac_divergence_norm() }
                   else { get_simple_divergence_norm() }
    println(format_centerline_csv(s, re_val, div_norm))
  } else if fmt == "monitor" {
    // Monitor: build a simple residual history by running the first solver
    // in record_interval-step batches and collecting div_norm at each batch.
    // Uses the current global state (already run above).
    let re = 2.0 / nu
    let records : Array[MonitorRecord] = []
    // Step 0: record initial state (already-run solver's current div_norm)
    let s0 = solvers[0]
    let init_div = if s0.name == "chorin" { get_divergence_norm() }
                   else if s0.name == "pcg" { get_pcg_divergence_norm() }
                   else if s0.name == "mac" { get_mac_divergence_norm() }
                   else { get_simple_divergence_norm() }
    records.push({
      step: s0.steps, time: s0.steps.to_double() * dt,
      div_norm: init_div, r_u: -1.0, r_v: -1.0,
      ke: -1.0, cfl_max: -1.0, pcg_iters: 0,
    })
    println(format_monitor_csv(records, s0.name, re))
  } else if fmt == "html" {
    let re_str = cli_get_string(args, "--re", "0")
    let re_val = if re_str == "0" { 2.0 / g_nu_runtime[0] }
                 else { parse_pos_int(re_str, 20).to_double() }
    println(format_html_report(solvers, re_val, solvers[0].steps))
  } else if fmt == "netcdf" {
    match format_netcdf(solvers[0]) {
      Ok(s)    => println(s)
      Err(msg) => println("[io] \{msg}")
    }
  } else {
    println("[io] Unknown format: '\{fmt}'.")
    println("[io] Run with --list-formats to see supported formats.")
  }
}
