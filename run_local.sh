#!/usr/bin/env bash
# run_local.sh – Build and run FlowLabLite locally.
#
# Usage:
#   bash run_local.sh                             # JSON output → results.json (default)
#   bash run_local.sh my_output.json              # JSON, custom filename
#   bash run_local.sh --format vtk                # VTK  → cavity_chorin.vtk
#   bash run_local.sh --format vtk --solver pcg   # VTK  → cavity_pcg.vtk
#   bash run_local.sh --format tecplot            # Tecplot → cavity_all.dat
#   bash run_local.sh --format csv                # CSV  → cavity_all.csv
#   bash run_local.sh --list-formats              # show available formats
#
# JSON mode: extracts between ===JSON_DATA_START=== / ===JSON_DATA_END=== markers
#            (compatible with cmd/main/local_viewer.html).
# Other formats: raw stdout written directly to file; no markers needed.
#
# Solver configuration: edit constants at top of cmd/main/main.mbt.

set -euo pipefail

# ── Parse arguments ────────────────────────────────────────────────────────────
FORMAT="json"
SOLVER="all"
OUTPUT_FILE=""
PASS_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --format)     FORMAT="$2";      PASS_ARGS+=("--format" "$2"); shift 2 ;;
    --solver)     SOLVER="$2";      PASS_ARGS+=("--solver" "$2"); shift 2 ;;
    --list-formats)                 PASS_ARGS+=("--list-formats"); FORMAT="list"; shift ;;
    --*)          PASS_ARGS+=("$1"); shift ;;
    *)
      # Positional arg = output filename (legacy: bash run_local.sh output.json)
      OUTPUT_FILE="$1"; shift ;;
  esac
done

# Default output filename per format
if [[ -z "$OUTPUT_FILE" ]]; then
  case "$FORMAT" in
    vtk)     OUTPUT_FILE="cavity_${SOLVER}.vtk" ;;
    vti)        OUTPUT_FILE="cavity_${SOLVER}.vti" ;;
    tecplot)    OUTPUT_FILE="cavity_all.dat" ;;
    csv)        OUTPUT_FILE="cavity_all.csv" ;;
    centerline) OUTPUT_FILE="centerline_${SOLVER}.csv" ;;
    monitor)    OUTPUT_FILE="monitor_${SOLVER}.csv" ;;
    html)       OUTPUT_FILE="report_${SOLVER}.html" ;;
    list)    OUTPUT_FILE="" ;;
    *)       OUTPUT_FILE="results.json" ;;
  esac
fi

echo "[run_local] Building and running FlowLabLite locally..."
echo "[run_local] Target: wasm (Wasmtime runtime)"
if [[ "$FORMAT" != "json" && "$FORMAT" != "list" ]]; then
  echo "[run_local] Format : $FORMAT  |  Solver: $SOLVER  |  Output: $OUTPUT_FILE"
fi

# ── Run the solver ─────────────────────────────────────────────────────────────
if [[ "${#PASS_ARGS[@]}" -gt 0 ]]; then
  # Format mode: stdout is clean format output; stderr has progress (suppressed by --quiet)
  FULL_OUTPUT=$(moon run cmd/main --target wasm -- "${PASS_ARGS[@]}" 2>/dev/null)
else
  # JSON mode: capture all output including timing (JSON extraction below)
  FULL_OUTPUT=$(moon run cmd/main --target wasm 2>&1)
fi

# ── Process output ─────────────────────────────────────────────────────────────
if [[ "$FORMAT" == "list" ]]; then
  echo "$FULL_OUTPUT"
  exit 0
fi

if [[ "$FORMAT" == "json" ]]; then
  echo "[run_local] Simulation completed. Extracting JSON..."
  echo "$FULL_OUTPUT" \
    | sed -n '/===JSON_DATA_START===/,/===JSON_DATA_END===/p' \
    | grep -v '===JSON_DATA' \
    > "$OUTPUT_FILE"

  if node -e "JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'))" 2>/dev/null; then
    CHORIN_PTS=$(node -e "
      const d=JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'));
      console.log((d.grid||[]).length);
    " 2>/dev/null || echo "?")
    SIMPLE_PTS=$(node -e "
      const d=JSON.parse(require('fs').readFileSync('$OUTPUT_FILE','utf8'));
      console.log((d.simple_grid||[]).length);
    " 2>/dev/null || echo "0")
    echo "[run_local] Valid JSON saved to: $OUTPUT_FILE"
    echo "[run_local]   Chorin grid  : $CHORIN_PTS points"
    echo "[run_local]   SIMPLE grid  : $SIMPLE_PTS points"
  else
    echo "[run_local] WARNING: JSON validation failed. File written: $OUTPUT_FILE"
  fi
  echo "[run_local] Open cmd/main/local_viewer.html in a browser and load $OUTPUT_FILE"
else
  # VTK / Tecplot / CSV: write stdout directly
  echo "$FULL_OUTPUT" > "$OUTPUT_FILE"
  local_lines=$(wc -l < "$OUTPUT_FILE" | tr -d ' ')
  echo "[run_local] Saved to: $OUTPUT_FILE  ($local_lines lines)"
fi
