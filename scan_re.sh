#!/usr/bin/env bash
# scan_re.sh – Parameter sweep over Reynolds numbers.
#
# Runs FlowLabLite at Re = 20, 100, 400 (default) or any list you pass.
# For each Re it produces:
#   - centerline_re<N>.csv   u/v centerline data
#   - report_re<N>.html      self-contained HTML report (Chart.js)
# If Re=100 is in the sweep it also runs --validate-ghia and writes:
#   - validation_summary.txt  Ghia L2/max errors per solver
#
# Usage:
#   bash scan_re.sh                        # default: Re=20,100,400
#   bash scan_re.sh 100 400 1000           # custom Re list
#   bash scan_re.sh --steps 2000           # override step count
#   bash scan_re.sh --solver pcg           # choose solver
#
# Output goes to the working directory.

set -euo pipefail

# ── Parse arguments ────────────────────────────────────────────────────────────
RE_LIST=()
STEPS_ARGS=()
SOLVER_ARGS=()
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --steps)  STEPS_ARGS=("--steps" "$2"); shift 2 ;;
    --solver) SOLVER_ARGS=("--solver" "$2"); shift 2 ;;
    --*)      EXTRA_ARGS+=("$1"); shift ;;
    *)        RE_LIST+=("$1"); shift ;;
  esac
done

if [[ ${#RE_LIST[@]} -eq 0 ]]; then
  RE_LIST=(20 100 400)
fi

SOLVER="${SOLVER_ARGS[1]:-all}"

echo "======================================================================"
echo " FlowLabLite — Re parameter sweep"
echo " Reynolds numbers : ${RE_LIST[*]}"
echo " Solver           : $SOLVER"
echo " Steps override   : ${STEPS_ARGS[*]:-<default>}"
echo "======================================================================"

VALIDATION_OUT="validation_summary.txt"
echo "FlowLabLite Ghia (1982) Validation Summary" > "$VALIDATION_OUT"
echo "Generated: $(date)"                          >> "$VALIDATION_OUT"
echo ""                                            >> "$VALIDATION_OUT"

for RE in "${RE_LIST[@]}"; do
  echo ""
  echo "── Re = $RE ──────────────────────────────────────────────────"

  # ── Centerline CSV ──────────────────────────────────────────────────────
  CL_FILE="centerline_re${RE}.csv"
  echo "[scan] Writing centerline → $CL_FILE"
  moon run cmd/main --target wasm -- \
    --re "$RE" "${STEPS_ARGS[@]}" "${SOLVER_ARGS[@]}" \
    --format centerline 2>/dev/null > "$CL_FILE"
  CL_LINES=$(wc -l < "$CL_FILE" | tr -d ' ')
  echo "[scan] $CL_FILE written ($CL_LINES lines)"

  # ── HTML report ─────────────────────────────────────────────────────────
  HTML_FILE="report_re${RE}.html"
  echo "[scan] Writing HTML report → $HTML_FILE"
  moon run cmd/main --target wasm -- \
    --re "$RE" "${STEPS_ARGS[@]}" "${SOLVER_ARGS[@]}" \
    --format html 2>/dev/null > "$HTML_FILE"
  echo "[scan] $HTML_FILE written"

  # ── Ghia validation (Re=100 only) ───────────────────────────────────────
  if [[ "$RE" == "100" ]]; then
    echo "[scan] Running Ghia (1982) validation at Re=100…"
    GHIA_OUT=$(moon run cmd/main --target wasm -- \
      --re 100 "${STEPS_ARGS[@]}" "${SOLVER_ARGS[@]}" \
      --validate-ghia 2>/dev/null)
    echo "$GHIA_OUT"
    echo "=== Re=100 Ghia Validation ===" >> "$VALIDATION_OUT"
    echo "$GHIA_OUT"                       >> "$VALIDATION_OUT"
    echo ""                                >> "$VALIDATION_OUT"
  fi
done

echo ""
echo "======================================================================"
echo " Sweep complete."
echo " Files generated:"
for RE in "${RE_LIST[@]}"; do
  echo "   centerline_re${RE}.csv   report_re${RE}.html"
done
if printf '%s\n' "${RE_LIST[@]}" | grep -q '^100$'; then
  echo "   $VALIDATION_OUT"
fi
echo "======================================================================"
