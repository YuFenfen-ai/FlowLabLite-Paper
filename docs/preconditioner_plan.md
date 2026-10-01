# PCG 预条件器实现规划文档

**项目**: FlowLabLite  
**作者**: Fenfen Yu（余芬芬）  
**日期**: 2026-04-16  
**状态**: 规划完成，待实现  

---

## 1  总体策略

### 1.1  设计原则

1. **不改变现有功能**：`pressure_poisson_pcg`（Jacobi）函数保持完全不变
2. **向后兼容**：所有新函数均为内部函数（无 `pub`），不影响 WASM 导出 API
3. **共享基础设施**：DILU 和 DIC 共享同一修正对角计算函数 `build_modified_diag`
4. **统一入口**：通过 `pressure_poisson_pcg_prec(p, dx, dy, b, precond_type)` 选择预条件器
5. **逐步验证**：每个预条件器独立实现、独立测试，再集成

### 1.2  预条件器类型编号

| 编号 | 名称 | 函数 |
|---|---|---|
| 0 | Jacobi（已有） | `pressure_poisson_pcg`（不变） |
| 1 | DILU | `pressure_poisson_pcg_dilu` |
| 2 | DIC | `pressure_poisson_pcg_dic` |
| 3 | GAMG | `pressure_poisson_pcg_gamg` |

---

## 2  代码结构规划

### 2.1  新增函数一览（插入位置：`pressure_poisson_pcg` 之后）

```
cmd/main/main.mbt
  ...（现有代码，不变）...
  pressure_poisson_pcg()              ← 现有，不改动
  ─── 新增区域 [Layer 2 ext] ────────────────────────────────────
  build_modified_diag()               ← DILU/DIC 共用修正对角
  apply_dilu_precond()                ← DILU apply（对角缩放）
  apply_dic_precond()                 ← DIC apply（3 步扫描）
  gamg_nc (常量)                      ← GAMG 粗网格大小 = 21
  apply_pressure_bcs_n()             ← 通用 n×n 压力 BC
  gamg_restrict()                     ← 注入限制算子
  gamg_prolongate()                   ← 双线性插值算子
  gamg_smooth_fine()                  ← 细网格阻尼 Jacobi 光滑
  gamg_coarse_solve()                 ← 粗网格 Jacobi 迭代
  apply_gamg_precond()                ← GAMG V-cycle apply
  pressure_poisson_pcg_prec()         ← 通用 PCG，按 precond_type 分支
  pressure_poisson_pcg_dilu()         ← 薄封装（precond_type=1）
  pressure_poisson_pcg_dic()          ← 薄封装（precond_type=2）
  pressure_poisson_pcg_gamg()         ← 薄封装（precond_type=3）
  ─── 以下现有代码，不改动 ─────────────────────────────────────
  cavity_flow_pcg()
  ...
```

### 2.2  新增测试一览（`cmd/main/main_wbtest.mbt` 末尾追加）

测试编号 43–58，共 16 个新测试：

| 编号 | 名称 | 验证内容 |
|---|---|---|
| 43 | `dilu_zero_rhs` | 零 RHS → 0 次迭代 |
| 44 | `dilu_pressure_bcs` | 求解后 BC 满足 Neumann/Dirichlet |
| 45 | `dilu_matches_jacobi` | 压力场与 Jacobi PCG 相差 < 1% |
| 46 | `dilu_modified_diag_positive` | 修正对角 d > 0（无零除风险） |
| 47 | `dilu_convergence_bounded` | 迭代次数 ≤ pcg_max_iter |
| 48 | `dilu_iters_vs_jacobi` | DILU 迭代数 ≤ Jacobi 迭代数 × 1.5 |
| 49 | `dic_zero_rhs` | 零 RHS → 0 次迭代 |
| 50 | `dic_pressure_bcs` | 求解后 BC 满足 |
| 51 | `dic_matches_jacobi` | 压力场与 Jacobi PCG 相差 < 1% |
| 52 | `dic_modified_diag_positive` | 修正对角 d > 0 |
| 53 | `dic_convergence_bounded` | 迭代次数 ≤ pcg_max_iter |
| 54 | `dic_iters_not_worse_than_jacobi` | DIC 迭代数 ≤ Jacobi × 1.5 |
| 55 | `gamg_zero_rhs` | 零 RHS → 0 次迭代 |
| 56 | `gamg_pressure_bcs` | 求解后 BC 满足 |
| 57 | `gamg_matches_jacobi` | 压力场与 Jacobi PCG 相差 < 1% |
| 58 | `gamg_convergence_bounded` | 迭代次数 ≤ pcg_max_iter |

**回归测试**：运行现有 42 个测试（T1–T42），全部保持通过。

---

## 3  分步实现计划

### 阶段 0：文档（已完成）

- [x] `docs/preconditioner_theory.md` — 理论推导、公式、OpenFOAM 参考
- [x] `docs/preconditioner_plan.md` — 本文件

---

### 阶段 1：DILU 预条件器

**目标**：实现最简单的改进预条件器，验证通用 PCG 框架可用。

**Step 1.1** — 实现共享基础函数

```moonbit
fn build_modified_diag(inv_dx2, inv_dy2, a_diag) -> Array[Array[Double]]
```

- 行主序扫描内部节点 (1,1)→(ny−2,nx−2)
- 公式：`d[i][j] = a_diag - inv_dx2²/d[i][j-1] - inv_dy2²/d[i-1][j]`
- 安全检查：若 d[i][j] ≤ 1e-15，回退到 a_diag

**Step 1.2** — 实现 DILU apply

```moonbit
fn apply_dilu_precond(r, z, d_inv)
```

- 逐元素：`z[i][j] = r[i][j] * d_inv[i][j]`

**Step 1.3** — 实现通用 PCG 框架

```moonbit
fn pressure_poisson_pcg_prec(p, dx_, dy_, b, precond_type) -> Int
```

- precond_type=0: Jacobi（与现有等效）
- precond_type=1: DILU
- precond_type=2/3: 占位，后续填充

**Step 1.4** — 添加 DILU 封装

```moonbit
fn pressure_poisson_pcg_dilu(p, dx_, dy_, b) -> Int
```

**Step 1.5** — 添加并运行 DILU 测试（T43–T48）

验证标准：
- T43: iters == 0 for zero b
- T44: Neumann/Dirichlet BCs satisfied after solve
- T45: `‖p_dilu - p_jacobi‖/‖p_jacobi‖ < 0.01`
- T46: all d[i][j] > 0
- T47: iters ≤ 200
- T48: iters_dilu ≤ iters_jacobi × 1.5

---

### 阶段 2：DIC 预条件器

**目标**：在 DILU 修正对角基础上增加前向/后向三角扫描。

**Step 2.1** — 实现 DIC apply（3 步）

```moonbit
fn apply_dic_precond(r, z, d, inv_dx2, inv_dy2)
```

- 步骤 1：前向扫描（行主序 1→ny-2, 1→nx-2）
  - `w[i][j] = r[i][j] + (inv_dx2/d[i][j-1])*w[i][j-1] + (inv_dy2/d[i-1][j])*w[i-1][j]`
- 步骤 2：对角缩放
  - `z[i][j] = w[i][j] / d[i][j]`
- 步骤 3：后向扫描（逆行主序 ny-2→1, nx-2→1）
  - `z[i][j] += (inv_dx2/d[i][j])*z[i][j+1] + (inv_dy2/d[i][j])*z[i+1][j]`
  - 使用 `ri = ny-2-…` 辅助索引实现递减（避免可能的 MoonBit 递减循环问题）

**Step 2.2** — 在 `pressure_poisson_pcg_prec` 中接入 DIC（precond_type=2）

**Step 2.3** — 添加 DIC 封装

```moonbit
fn pressure_poisson_pcg_dic(p, dx_, dy_, b) -> Int
```

**Step 2.4** — 添加并运行 DIC 测试（T49–T54）

---

### 阶段 3：GAMG 预条件器

**目标**：实现 2-level V-cycle，最复杂但收敛最快。

**Step 3.1** — 添加 GAMG 常量和 BC 辅助函数

```moonbit
let gamg_nc : Int = (nx + 1) / 2   // = 21
fn apply_pressure_bcs_n(p, n)      // n×n 网格通用 BC
```

**Step 3.2** — 实现限制算子

```moonbit
fn gamg_restrict(r_fine, r_coarse)  // 注入：r_c[i][j] = r_f[2i][2j]
```

**Step 3.3** — 实现插值算子

```moonbit
fn gamg_prolongate(e_coarse, e_fine)
```

4 个 for 循环，分别处理：偶偶、偶奇、奇偶、奇奇 节点类型。

**Step 3.4** — 实现阻尼 Jacobi 光滑（细网格）

```moonbit
fn gamg_smooth_fine(z, r, inv_dx2, inv_dy2, a_diag, nu)
```

- ω = 2/3，nu 次迭代
- 每次：apply_pressure_bcs → laplacian_apply → z += (2/3)/a_diag*(r-Az)

**Step 3.5** — 实现粗网格 Jacobi 求解

```moonbit
fn gamg_coarse_solve(e_c, r_c, inv_dx2_c, inv_dy2_c, a_diag_c, n_iters)
```

- 同样用阻尼 Jacobi，n_iters=20

**Step 3.6** — 实现 GAMG V-cycle

```moonbit
fn apply_gamg_precond(r, z, inv_dx2, inv_dy2, a_diag)
```

V-cycle 6 个步骤（见理论文档）。

**Step 3.7** — 在 `pressure_poisson_pcg_prec` 中接入 GAMG（precond_type=3）

**Step 3.8** — 添加 GAMG 封装

```moonbit
fn pressure_poisson_pcg_gamg(p, dx_, dy_, b) -> Int
```

**Step 3.9** — 添加并运行 GAMG 测试（T55–T58）

---

### 阶段 4：集成测试与提交

**Step 4.1** — 运行全量测试

```bash
moon test --target wasm
```

预期结果：58 个测试全部通过（原 42 个 + 新增 16 个）。

**Step 4.2** — 本地求解器验证

```bash
bash run_local.sh
```

检查：PCG 各预条件器均产生正确的流场统计（max_u ≈ 1, min_p < 0 等）。

**Step 4.3** — 更新 CLAUDE.md

测试总数从 42 更新为 58，添加预条件器测试范围描述。

**Step 4.4** — git commit

按模块分两次提交：
1. 文档：`docs/preconditioner_theory.md`, `docs/preconditioner_plan.md`
2. 代码：`cmd/main/main.mbt`, `cmd/main/main_wbtest.mbt`, `CLAUDE.md`

---

## 4  风险与注意事项

### 4.1  MoonBit 递减循环

DIC 后向扫描需要从高下标到低下标遍历。
若 `for i = n; i >= 1; i = i - 1` 在 MoonBit 中不支持，使用辅助计数器：

```moonbit
// 辅助计数器方法（规避递减）
for ri = 0; ri < ny - 2; ri = ri + 1 {
  let i = ny - 2 - ri
  for ci = 0; ci < nx - 2; ci = ci + 1 {
    let j = nx - 2 - ci
    // process (i, j)
  }
}
```

### 4.2  修正对角数值稳定性

对 Re=20、41×41 均匀网格，修正对角 d[i][j] 预期在 [1400, 1600] 之间，
不会出现零除问题。若扩展到高 Re 或非均匀网格，需要更健壮的稳定化。

### 4.3  GAMG 内存分配

`apply_gamg_precond` 每次调用分配约 7600 个 Double：
- 细网格工作数组：3 × 1681 ≈ 5043
- 粗网格数组：4 × 441 ≈ 1764

在 PCG 内循环（每步最多 200 次）中，这会增加 GC 压力。
若性能成为问题，可将工作数组提升为全局变量（后续优化）。

### 4.4  GAMG 不改变 WASM API

所有新函数均为内部函数（无 `pub` 修饰），
`cmd/main/moon.pkg.json` 中的 WASM 导出列表无需修改，
`build_wasm.sh` 也无需更改。

---

## 5  验收标准

| 验收项 | 标准 |
|---|---|
| 回归测试 | 原 42 个测试全部通过 |
| 新增测试 | 新增 16 个测试（T43–T58）全部通过 |
| Jacobi 不变性 | `pressure_poisson_pcg` 行为与修改前完全一致 |
| DILU 正确性 | `‖p_dilu - p_jacobi‖/‖p_jacobi‖ < 0.01` |
| DIC 正确性 | 同上 |
| GAMG 正确性 | 同上 |
| WASM API | 导出函数数量不变（仍为 62 个） |
| 编译 | `moon build` 和 `bash build_wasm.sh` 无错误 |
