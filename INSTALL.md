# FlowLabLite — Installation & Deployment Guide

---

## Prerequisites

| Tool | Minimum | How to install |
|---|---|---|
| `moon` (MoonBit build tool) | 0.1.20260309 | https://www.moonbitlang.com/download/ |
| `moonc` (compiler, bundled with moon) | 0.8.3 | bundled |
| `wasmtime` (test/run runtime, bundled) | any | bundled |
| `node` (export verification only) | 18+ | https://nodejs.org/ |
| Chrome / Edge / Brave | 115+ | required for wasm-gc stringref |
| HTTP server | any | see options below |

---

## Platform-Specific Notes

### Windows
```powershell
# Run all shell scripts via Git Bash (not PowerShell):
bash build_wasm.sh release
bash run_local.sh
```
Git Bash is included with [Git for Windows](https://git-scm.com/download/win).

### macOS
```bash
# No special steps. Homebrew Python provides the HTTP server:
brew install python   # if not already installed
```

### Linux
```bash
sudo apt install python3   # or equivalent
```

---

## Option A — Browser Visualisation (Recommended)

### Step 1 — Clone

```bash
git clone https://github.com/YuFenfen-ai/FlowLabLite.git
cd FlowLabLite
```

### Step 2 — Build the WASM binary

```bash
bash build_wasm.sh release
```

This produces `_build/wasm-gc/release/build/cmd/main/main.wasm`  
(62 exported API functions + `_start`).

> **Why not `moon build` alone?**  
> moon 0.1.20260309 does not forward the export list to the WASM linker.
> `build_wasm.sh` calls `moonc link-core` directly with `-exported_functions`.

### Step 3 — Serve over HTTP

The WASM file must be served over HTTP (not opened as a local `file://` URL).

```bash
# Python (any platform)
python3 -m http.server 8080

# Node.js alternative
npx serve . -p 8080

# VS Code Live Server extension (right-click cmd/main/main.html → Open with Live Server)
```

### Step 4 — Open in browser

```
http://localhost:8080/cmd/main/main.html
```

Click **"Run 50 Steps"** to start the simulation. The velocity heatmap,
pressure heatmap, and **streamline overlay** render automatically.
Toggle streamlines with the **"Streamlines"** button in the toolbar.

---

## Option B — Local JSON Workflow (All Four Solvers)

This option runs all four solvers natively (via Wasmtime) and saves
full-field results to `results.json`, then views them in a browser.

### Step 1 — Run

```bash
bash run_local.sh               # saves to results.json
bash run_local.sh my_run.json   # custom output filename
```

Output:
```
[run_local] Simulation completed. Extracting JSON...
[run_local] Valid JSON saved to: results.json
[run_local]   Chorin grid  : 1681 points
[run_local]   SIMPLE grid  : 1681 points
```

### Step 2 — View

Open `cmd/main/local_viewer.html` in any browser (local `file://` works here).  
Drag and drop `results.json` onto the page.

Four solver tabs appear: **Chorin / SIMPLE / PCG / MAC**.  
Each tab shows velocity magnitude, u/v components, and pressure heatmaps
with statistics.

---

## Option C — MoonBit Native Mode (Command-Line Only)

Run the solver without a browser, printing statistics to stdout:

```bash
moon run cmd/main --target wasm    # runs via Wasmtime, prints JSON
```

To run on native (non-WASM) target for faster local profiling:

```bash
moon run cmd/main                  # native binary (no browser export)
```

> Native mode does not produce a WASM file; it is for local debugging only.

---

## Option D — Running the Test Suite

```bash
moon test --target wasm            # 176 tests, all should pass
```

Expected output:
```
Total tests: 176, passed: 176, failed: 0.
```

Run a subset:
```bash
moon test --target wasm --filter gamg    # run all GAMG-related tests
moon test --target wasm --filter T71     # run test named T71
```

---

## Troubleshooting

### WASM file not found in browser

**Symptom**: `Status: Error loading WASM` in the page header.

**Fix**: Build with `bash build_wasm.sh release` first, then serve over HTTP.
`main.html` tries two paths:
```
../../_build/wasm-gc/release/build/cmd/main/main.wasm   ← release
../../_build/wasm-gc/debug/build/cmd/main/main.wasm     ← debug fallback
```

### `LinkError: WebAssembly.instantiate` in browser console

**Symptom**: Module loads but functions throw `LinkError`.

**Fix**: The import object must include `spectest / env / moonbit:ffi` namespaces.
`main.html` provides a `Proxy`-based catch-all import object automatically.
This error only appears if you are loading the WASM outside of `main.html`.

### Chrome shows blank canvas

**Symptom**: Page loads, buttons work, but canvas stays grey.

**Cause**: Chrome < 115 does not support wasm-gc stringref (type code `0x77`).

**Fix**: Update Chrome to 115 or later, or use Edge/Brave 115+.

### `moon test` fails on native target

**Symptom**: Some tests fail when running without `--target wasm`.

**Explanation**: All tests are written for wasm-gc semantics. The `--target wasm`
flag is **required**. Native target is not officially supported for this test suite.

### Wasmtime not found

**Symptom**: `moon run` or `moon test` fails with "wasmtime: command not found".

**Fix**: Wasmtime is bundled with moon. Reinstall moon from
https://www.moonbitlang.com/download/ to restore the bundled runtime.

---

## Build Script Reference

### `build_wasm.sh`

```bash
bash build_wasm.sh           # debug build (faster, includes source map)
bash build_wasm.sh release   # release build (optimised, smaller)
```

The script:
1. Runs `moon build --target wasm-gc [--release]`
2. Re-links the `.core` file with `moonc link-core -exported_functions <list>`
3. Prints the export count (should be 64 = 62 API + _start + memory)

### `run_local.sh`

```bash
bash run_local.sh [output_filename.json]
```

Builds (debug), runs via `wasmtime`, extracts the JSON block from stdout,
and saves to `results.json` (or custom filename).
