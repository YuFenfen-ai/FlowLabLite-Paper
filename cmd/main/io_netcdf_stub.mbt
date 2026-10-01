// FlowLabLite — NetCDF-4 Stub  (Tier-2, always compiles)
//
// This file provides the NetCDF-4 interface that the rest of the package uses.
// It is a *stub*: all calls return Err with a human-readable explanation.
//
// To enable real NetCDF-4 output (Tier-2, optional):
//   1. Install libnetcdf >= 4.4:
//        apt install libnetcdf-dev        (Debian/Ubuntu)
//        brew install netcdf              (macOS)
//        conda install netcdf4            (all platforms)
//   2. Implement C FFI bindings via moonbit-c-binding (see io.md §4.2).
//   3. Replace this file with the FFI implementation.
//
// Compiling *without* libnetcdf: always succeeds — only this stub is compiled,
// and attempting --format netcdf at runtime prints a clear error message.

/// Whether NetCDF-4 output is available in this build.
/// Always false for the stub; set to true in the FFI implementation.
let netcdf_available : Bool = false

/// Export solver data in CF-1.11 NetCDF-4 format.
///
/// This stub always returns Err.  See module comment for upgrade instructions.
fn format_netcdf(_ : SolverExport) -> Result[String, String] {
  Err(
    "NetCDF-4 output is not available in this build.\n" +
    "Requires libnetcdf >= 4.4.  See io.md §4.2 for installation instructions."
  )
}
