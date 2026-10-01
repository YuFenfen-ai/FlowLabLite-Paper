#!/usr/bin/env bash
# build_wasm.sh – Build FlowLabLite WASM with all exported functions.
#
# Usage:
#   bash build_wasm.sh          # debug build (faster, with source-map)
#   bash build_wasm.sh release  # release/optimised build
#
# Why this script exists:
#   moon 0.1.20260309 passes "-pkg-config-path ./cmd/main/moon.pkg" to
#   moonc link-core.  The old-format moon.pkg does not carry a link/exports
#   section, so the compiled WASM only exports "_start".  This script re-runs
#   the link step with an explicit "-exported_functions" flag so all 27 getter/
#   runner functions are visible to the browser.
export PATH="$HOME/.moon/bin:$PATH"
set -euo pipefail

MOON_HOME="$HOME/.moon"
CORE_BUNDLE="$MOON_HOME/lib/core/_build/wasm-gc/release/bundle"
EXPORTS="init_simulation,run_all_steps,run_n_steps,get_nx,get_ny,get_nt,get_nit,get_re,get_dx,get_dy,get_dt,get_rho,get_nu,get_step_count,get_u_at,get_velocity_magnitude_at,get_v_at,get_p_at,get_u_center,get_v_center,get_p_center,get_divergence_norm,get_max_velocity_magnitude,get_max_u,get_max_v,get_max_p,get_min_p,init_simple,run_simple_n_iter,get_simple_step_count,get_simple_residual,get_u_simple_at,get_v_simple_at,get_p_simple_at,get_max_u_simple,get_simple_divergence_norm,init_chorin_pcg,run_chorin_pcg_n_steps,get_pcg_step_count,get_pcg_last_iters,get_u_pcg_at,get_v_pcg_at,get_p_pcg_at,get_velocity_magnitude_pcg_at,get_max_u_pcg,get_max_v_pcg,get_max_p_pcg,get_min_p_pcg,get_pcg_divergence_norm,init_mac,run_mac_n_steps,get_mac_step_count,get_mac_last_iters,get_mac_nc,get_u_mac_at,get_v_mac_at,get_p_mac_at,get_velocity_magnitude_mac_at,get_max_u_mac,get_max_v_mac,get_max_p_mac,get_min_p_mac,get_mac_divergence_norm,init_rk3,run_rk3_n_steps,get_rk3_step_count,get_u_rk3_at,get_v_rk3_at,get_p_rk3_at,get_rk3_divergence_norm,init_fvm,run_fvm_n_iter,get_fvm_step_count,get_u_fvm_at,get_v_fvm_at,get_p_fvm_at,get_fvm_divergence_norm,init_scalar,set_scalar_alpha,get_scalar_alpha,run_scalar_n_steps,get_scalar_step_count,get_scalar_at,get_scalar_mean,init_boussinesq,set_bous_nu,set_bous_alpha,set_bous_beta,get_bous_nu,get_bous_alpha,get_bous_beta,run_boussinesq_n_steps,get_bous_step_count,get_u_bous_at,get_v_bous_at,get_p_bous_at,get_phi_bous_at,get_bous_nusselt,init_channel,set_channel_fx,set_channel_nu,get_channel_fx,run_channel_n_steps,get_channel_step_count,get_u_channel_at,get_v_channel_at,get_p_channel_at,get_channel_analytic_u,set_grid_size,get_tanh_y_at,get_tanh_dy_at,get_tanh_mean_dy,set_curv_c,get_curv_c,get_curv_x_at,get_curv_y_at,get_curv_jacobian_at,get_curv_nonortho_at,get_curv_g11_at,get_curv_g12_at,get_curv_g22_at"

MODE="${1:-debug}"
if [[ "$MODE" == "release" ]]; then
  BUILD_DIR="./_build/wasm-gc/release/build/cmd/main"
  OPT_FLAGS=""
  echo "[build_wasm] Building RELEASE wasm-gc..."
  moon build --target wasm-gc --release 2>&1 | grep -v "^$"
else
  BUILD_DIR="./_build/wasm-gc/debug/build/cmd/main"
  OPT_FLAGS="-O0 -g -source-map"
  echo "[build_wasm] Building DEBUG wasm-gc..."
  moon build --target wasm-gc 2>&1 | grep -v "^$"
fi

CORE_FILE="$BUILD_DIR/main.core"
OUT_FILE="$BUILD_DIR/main.wasm"

echo "[build_wasm] Re-linking $CORE_FILE -> $OUT_FILE with exported_functions..."
moonc link-core \
  "$CORE_BUNDLE/abort/abort.core" \
  "$CORE_BUNDLE/core.core" \
  "$CORE_FILE" \
  -main "YuFenfen-ai/FlowLabLite/cmd/main" \
  -o "$OUT_FILE" \
  -exported_functions "$EXPORTS" \
  -export-memory-name memory \
  -pkg-config-path "./cmd/main/moon.pkg" \
  -pkg-sources "YuFenfen-ai/FlowLabLite/cmd/main:./cmd/main" \
  -pkg-sources "moonbitlang/core:$MOON_HOME/lib/core" \
  -target wasm-gc \
  $OPT_FLAGS

echo "[build_wasm] Done: $OUT_FILE"
echo "[build_wasm] Checking exports in WASM binary..."
node - <<'JSEOF'
const fs = require('fs');
const paths = [
  './_build/wasm-gc/release/build/cmd/main/main.wasm',
  './_build/wasm-gc/debug/build/cmd/main/main.wasm'
];
for (const p of paths) {
  if (!fs.existsSync(p)) continue;
  const buf = Buffer.from(fs.readFileSync(p));
  let off = 8;
  const exports = [];
  while (off < buf.length) {
    const sid = buf[off++];
    let sz = 0, sh = 0, b;
    do { b = buf[off++]; sz |= (b & 0x7f) << sh; sh += 7; } while (b & 0x80);
    if (sid === 7) {
      let p2 = off;
      let cnt = 0, sh2 = 0, b2;
      do { b2 = buf[p2++]; cnt |= (b2 & 0x7f) << sh2; sh2 += 7; } while (b2 & 0x80);
      for (let i = 0; i < cnt; i++) {
        let nl = 0, sh3 = 0, b3;
        do { b3 = buf[p2++]; nl |= (b3 & 0x7f) << sh3; sh3 += 7; } while (b3 & 0x80);
        const name = buf.slice(p2, p2 + nl).toString('utf8');
        p2 += nl; const kind = buf[p2++];
        let idx = 0, sh4 = 0, b4;
        do { b4 = buf[p2++]; idx |= (b4 & 0x7f) << sh4; sh4 += 7; } while (b4 & 0x80);
        exports.push(name);
      }
      console.log('[build_wasm] WASM exports (' + exports.length + '):', exports.join(', '));
      break;
    }
    off += sz;
  }
  break;
}
JSEOF
