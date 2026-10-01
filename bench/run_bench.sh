#!/usr/bin/env bash
# run_bench.sh — FlowLabLite multi-language benchmark runner
#
# Languages: Python (numpy + pure), Java, C (gcc -O2), MoonBit (wasm)
# Grids: 41×41 (small), 81×81 (medium), 161×161 (large), each 500 steps
# Output: console + docs/multilang_ben_results.tsv
#
# Usage:
#   bash bench/run_bench.sh               # run all available languages
#   bash bench/run_bench.sh --quick       # 41×41 only (fast preview)
#   bash bench/run_bench.sh --lang java   # run only Java

set -euo pipefail
cd "$(dirname "$0")/.."   # project root

QUICK=0
LANG_FILTER=""
for arg in "$@"; do
  case "$arg" in
    --quick) QUICK=1 ;;
    --lang)  shift; LANG_FILTER="$1" ;;
  esac
done

BENCH_DIR="bench"
OUT_TSV="docs/multilang_ben_results.tsv"

echo "======================================================================"
echo " FlowLabLite — multi-language Chorin+GS benchmark"
echo " $(date)"
echo "======================================================================"
echo ""

# ── Helper: detect tool availability ─────────────────────────────────────────
have() { command -v "$1" &>/dev/null; }

# ── Results accumulator ───────────────────────────────────────────────────────
declare -A RESULTS   # key = "grid_lang"  value = "best_ms"
GRIDS=("41x41" "81x81" "161x161")
NTS=(500 500 500)
REPS=(5 3 3)

if [[ "$QUICK" == "1" ]]; then
  GRIDS=("41x41"); NTS=(500); REPS=(5)
fi

record() {
  local grid="$1" lang="$2" ms="$3"
  RESULTS["${grid}_${lang}"]="$ms"
}

# ── Java ─────────────────────────────────────────────────────────────────────
if [[ -z "$LANG_FILTER" || "$LANG_FILTER" == "java" ]]; then
  if have javac && have java; then
    echo "── Java (JDK $(java -version 2>&1 | head -1)) ──────────────────"
    javac "$BENCH_DIR/ChorinGSBench.java" -d "$BENCH_DIR/"
    for idx in "${!GRIDS[@]}"; do
      G="${GRIDS[$idx]}"; N="${G%x*}"; NT="${NTS[$idx]}"; RP="${REPS[$idx]}"
      echo -n "  $G / ${NT}st / ${RP}x ... "
      OUT=$(java -server -cp "$BENCH_DIR" ChorinGSBench "$N" "$NT" "$RP" 2>/dev/null)
      # Extract best time from TSV block
      MS=$(echo "$OUT" | grep -A9999 "TSV_RESULTS" | grep "^${N}x" | awk -F'\t' '{print $3}' | head -1)
      echo "${MS} ms (best)"
      record "$G" "java" "$MS"
    done
    echo ""
  else
    echo "⚠  Java not found — skipping"
    echo ""
  fi
fi

# ── Python ───────────────────────────────────────────────────────────────────
if [[ -z "$LANG_FILTER" || "$LANG_FILTER" == "python" ]]; then
  PY=""
  have python3 && PY=python3
  have python  && PY=python
  if [[ -n "$PY" ]]; then
    echo "── Python ($($PY --version 2>&1)) ───────────────────────────────"
    for idx in "${!GRIDS[@]}"; do
      G="${GRIDS[$idx]}"; N="${G%x*}"; NT="${NTS[$idx]}"; RP="${REPS[$idx]}"
      echo -n "  $G numpy ... "
      OUT=$($PY "$BENCH_DIR/chorin_gs_bench.py" "$N" "$NT" "$RP" 2>/dev/null)
      NP_MS=$(echo "$OUT" | grep -A9999 "TSV_RESULTS" | grep "^${N}x" | awk -F'\t' '{print $3}' | head -1)
      PU_MS=$(echo "$OUT" | grep -A9999 "TSV_RESULTS" | grep "^${N}x" | awk -F'\t' '{print $4}' | head -1)
      echo "numpy=${NP_MS}ms  pure=${PU_MS}ms"
      record "$G" "python_numpy" "$NP_MS"
      record "$G" "python_pure"  "$PU_MS"
    done
    echo ""
  else
    echo "⚠  Python not found — skipping"
    echo ""
  fi
fi

# ── C (gcc -O2) ──────────────────────────────────────────────────────────────
if [[ -z "$LANG_FILTER" || "$LANG_FILTER" == "c" ]]; then
  if have gcc; then
    echo "── C (gcc $(gcc --version | head -1)) ─────────────────────────"
    gcc -O2 -o "$BENCH_DIR/chorin_gs_bench_c" "$BENCH_DIR/chorin_gs_bench.c" -lm
    for idx in "${!GRIDS[@]}"; do
      G="${GRIDS[$idx]}"; N="${G%x*}"; NT="${NTS[$idx]}"; RP="${REPS[$idx]}"
      echo -n "  $G ... "
      OUT=$("$BENCH_DIR/chorin_gs_bench_c" "$N" "$NT" "$RP" 2>/dev/null)
      MS=$(echo "$OUT" | grep -A9999 "TSV_RESULTS" | grep "^${N}x" | awk -F'\t' '{print $3}' | head -1)
      echo "${MS} ms (best)"
      record "$G" "c" "$MS"
    done
    echo ""
  else
    echo "⚠  gcc not found — skipping"
    echo ""
  fi
fi

# ── MoonBit (wasm via wasmtime) ───────────────────────────────────────────────
if [[ -z "$LANG_FILTER" || "$LANG_FILTER" == "moonbit" ]]; then
  if have moon; then
    echo "── MoonBit ($(moon version 2>&1 | head -1), wasm target) ──────"
    # MoonBit runs the full solver (41×41 default). We time the moon run call.
    for idx in "${!GRIDS[@]}"; do
      G="${GRIDS[$idx]}"; NT="${NTS[$idx]}"; RP="${REPS[$idx]}"
      # Only 41×41 is natively supported (nx/ny compile-time constants)
      if [[ "$G" != "41x41" ]]; then
        echo "  $G — MoonBit grid size fixed at 41×41, skip"
        record "$G" "moonbit" "N/A"
        continue
      fi
      echo -n "  $G ($RP runs) ... "
      BEST_MS=999999
      for ((r=0; r<RP; r++)); do
        T0=$(date +%s%3N)
        moon run cmd/main --target wasm -- --format csv --solver chorin > /dev/null 2>&1
        T1=$(date +%s%3N)
        ELAPSED=$((T1-T0))
        [[ $ELAPSED -lt $BEST_MS ]] && BEST_MS=$ELAPSED
      done
      echo "${BEST_MS} ms (best, includes moon startup)"
      record "$G" "moonbit" "$BEST_MS"
    done
    echo ""
  else
    echo "⚠  moon not found — skipping"
    echo ""
  fi
fi

# ── Write TSV results ─────────────────────────────────────────────────────────
mkdir -p docs
{
  echo "# FlowLabLite benchmark results — $(date)"
  echo "grid	nt	java_ms	python_numpy_ms	python_pure_ms	c_gcc_O2_ms	moonbit_wasm_ms"
  for idx in "${!GRIDS[@]}"; do
    G="${GRIDS[$idx]}"; NT="${NTS[$idx]}"
    J="${RESULTS[${G}_java]:-N/A}"
    PN="${RESULTS[${G}_python_numpy]:-N/A}"
    PP="${RESULTS[${G}_python_pure]:-N/A}"
    CC="${RESULTS[${G}_c]:-N/A}"
    MB="${RESULTS[${G}_moonbit]:-N/A}"
    echo "$G	$NT	$J	$PN	$PP	$CC	$MB"
  done
} > "$OUT_TSV"

echo "======================================================================"
echo " Results written to: $OUT_TSV"
echo "======================================================================"
cat "$OUT_TSV"
