# FlowLabLite — solver_monitor.md
# 求解器监控参数规范

版本：v0.1 · 2026-04-18
作者：Fenfen Yu（余芬芬）

---

## 1. 稳态 vs 瞬态求解器分类

### 1.1 稳态求解器（SIMPLE）

| 特征 | 说明 |
|---|---|
| 时间 | 无时间步长，只有迭代步 k |
| 收敛判据 | 前后两迭代步差异（残差）< 阈值，或方程残差 < 阈值 |
| 直接法 | 当前 FlowLabLite **不存在**稳态直接法（SIMPLE 是唯一稳态求解器，且为迭代法） |
| 典型阈值 | 方程残差 < 1e-4（流动），< 1e-6（高精度） |

**参考**：OpenFOAM simpleFoam `residualControl`；Fluent SIMPLE 残差监控

### 1.2 瞬态求解器（Chorin / Chorin-PCG / MAC）

| 特征 | 说明 |
|---|---|
| 外层循环 | 时间步 n，步长 dt = 0.001，总时间 T = nt × dt |
| 内层循环 | 每时间步压力 Poisson 求解（GS/PCG） |
| 每步近似 | 时间步内可视为准稳态子问题 |
| 全局监控 | 质量守恒（div_norm）、动能（KE）、CFL 数 |
| 内层监控 | 压力求解残差 r_p、PCG 迭代次数 |

**参考**：OpenFOAM pimpleFoam `residualControl`；Fluent PISO 残差；Chorin (1968)

---

## 2. 监控参数定义

### 2.1 散度范数（质量守恒）— 所有求解器 ★已实现★

**含义**：速度场不可压缩性的量化指标。对理想不可压流，值恒为 0。
是最重要的全局守恒指标，对应 Fluent "continuity" 残差、OpenFOAM continuityErrors。

#### Chorin / SIMPLE / PCG（同位网格 41×41，中心差分）

```
div_norm = (1 / N_int) * Σ_{i=1}^{ny-2} Σ_{j=1}^{nx-2}
             |(u[i][j+1] - u[i][j-1]) / (2·dx) + (v[i+1][j] - v[i-1][j]) / (2·dy)|
```

- N_int = (ny-2)×(nx-2) 内部节点数
- 单位：s⁻¹（速度/长度）
- 代码：`get_divergence_norm()` 行 241；`get_simple_divergence_norm()` 行 803；`get_pcg_divergence_norm()` 行 1909

#### MAC（交错网格 40×40，单侧差分）

```
div_norm_MAC = sqrt( (1/N_int) * Σ_{i=1}^{nc-2} Σ_{j=1}^{nc-2}
                     ((u_mac[i][j+1] - u_mac[i][j]) / dx
                    + (v_mac[i+1][j] - v_mac[i][j]) / dy)^2 )
```

- MAC 使用 RMS（均方根）而非平均绝对值，因交错网格散度在单元中心精确成立
- 代码：`get_mac_divergence_norm()` 行 2384

**收敛准则**：

| 状态 | 阈值 | 说明 |
|---|---|---|
| 良好 | < 1e-4 | 质量守恒满足 |
| 警告 | 1e-4 ~ 1e-2 | 建议增加压力迭代次数 |
| 发散 | > 1e-2 | 求解器不稳定 |

---

### 2.2 压力 Poisson 求解残差（内层迭代）— PCG 求解器 ★已实现（iters）★

**含义**：每时间步内，压力方程 A·p = b 的求解精度。
对应 OpenFOAM p 残差、Fluent pressure 残差。

```
r_p^(k) = ||r^(k)||_2 / ||b||_2,   r^(k) = b - A·p^(k)
```

- 停止准则（代码行 942）：r_p < 1e-5（pcg_tol）或 k ≥ 200（pcg_max_iter）
- 当前已导出：`get_pcg_last_iters()`、`get_mac_last_iters()`（迭代次数）

**RHS b 公式**：

Chorin / PCG（含交叉导数项，代码行 350-376）：
```
b[i][j] = rho * [ (1/dt) * (∂u/∂x + ∂v/∂y) - (∂u/∂x)² - 2·(∂u/∂y)·(∂v/∂x) - (∂v/∂y)² ]
```

MAC（仅散度项，无交叉导数，代码行 2167）：
```
b[i][j] = (rho/dt) * ((u_s[i][j+1] - u_s[i][j]) / dx + (v_s[i+1][j] - v_s[i][j]) / dy)
```

---

### 2.3 SIMPLE 方程残差（稳态迭代）— 两种模式 ★Step M 新增★

#### 模式 A：简化残差（默认，`--simple-residual div`）

当前实现（代码行 773）：
```
R_simple = div_norm(u^(k), v^(k))
```
即 SIMPLE 残差 = 速度场散度范数，迭代结束后缓存到 `g_simple_residual[0]`。

#### 模式 B：严格方程残差（OpenFOAM 风格，`--simple-residual eq`）

归一化动量方程不平衡：
```
R_u = Σ_P |a_P·u_P* - Σ_{nb} a_nb·u_nb* - b_P| / (Σ_P a_P·|u_ref|)

其中：
  a_P = 1/dt·rho·dx·dy（时间项系数）
  a_nb = 对流 + 扩散系数（相邻节点）
  b_P  = 压力梯度 + 源项
  u_ref = max(|u*|) over domain（归一化参考值）
```

实现要求：在 `simple_one_iter` 完成后，纯读取 g_u_s/g_v_s/g_p_s 数组，计算不平衡，不回写任何值。

---

### 2.4 速度变化量（瞬态收敛趋势）— ★Step M 新增★

**含义**：相邻时间步速度的最大相对变化。瞬态趋近稳态时趋向 0，也可检测发散。
对应 Fluent "x-velocity"/"y-velocity" 残差（瞬态模式）。

```
R_u^(n) = max_{i,j} |u^(n+1)[i][j] - u^(n)[i][j]| / (max_{i,j} |u^(n+1)[i][j]| + 1e-10)
R_v^(n) = max_{i,j} |v^(n+1)[i][j] - v^(n)[i][j]| / (max_{i,j} |v^(n+1)[i][j]| + 1e-10)
```

ε = 1e-10 防止零除。

**实现**：io_cli.mbt 在每 record_interval 步调用 compute_ru/rv，保存 u_prev 副本（约 2×41×41×8 B ≈ 27 KB 额外内存，可接受）。

---

### 2.5 CFL 数（显式稳定性）— ★Step M 新增★

**含义**：Courant-Friedrichs-Lewy 数，显式格式稳定性指标。
CFL < 1 为 Chorin 显式格式的必要稳定条件（对流主导）。

```
CFL_max = max_{i,j} (|u[i][j]|·dt/dx + |v[i][j]|·dt/dy)
```

当前参数（Re=20，dt=0.001，dx=0.05）：最大速度 ≈ 1.0 → CFL_max ≈ 0.02，远小于 1。

| CFL | 状态 | 说明 |
|---|---|---|
| < 0.5 | 安全 | 显式格式稳定 |
| 0.5 ~ 1.0 | 警告 | 接近稳定边界 |
| ≥ 1.0 | 错误 | 数值发散风险，需减小 dt 或 Re |

**WASM 导出**：`get_cfl_max()` → Step M 新增（读取 g_u/g_v 全局）

---

### 2.6 动能（流场整体能量）— ★Step M 新增★

**含义**：全场总动能，稳态时趋于常数，发散时无界增长。
对应 OpenFOAM `fieldAverage` 函数对象。

```
KE = 0.5 * Σ_{i=0}^{ny-1} Σ_{j=0}^{nx-1} (u[i][j]² + v[i][j]²) · dx · dy
```

**WASM 导出**：`get_kinetic_energy()` → Step M 新增（读取 g_u/g_v 全局）

---

### 2.7 速度变化量 WASM 导出 — ★Step M 新增★

WASM 导出 `get_ru()` / `get_rv()` 供 main.html 实时监控曲线使用。
这两个函数读取内部缓存值（在 `run_n_steps` 最后一批完成后更新）。

---

### 2.8 中心线速度（基准验证）— ★已有 API★

u 沿竖向中心线（x=1.0，j=nx/2=20）：
```
u_vcl[i] = u[i][20],   y_vcl[i] = i · dy,   i = 0..ny-1
```

v 沿横向中心线（y=1.0，i=ny/2=20）：
```
v_hcl[j] = v[20][j],   x_hcl[j] = j · dx,   j = 0..nx-1
```

**Ghia (1982) L2 误差**（Step D 实现）：
```
L2 = sqrt( (1/N_G) * Σ_{k=1}^{N_G} (u_our(y_k_G) - u_Ghia(y_k_G))² )
```

坐标映射：y_our = y_Ghia × 2.0（Ghia 域 [0,1] → 本项目域 [0,2]）。

---

## 3. 监控输出格式

### 3.1 统一监控 CSV（`--format monitor`）

```
# FlowLabLite monitor | solver=chorin | Re=20 | dt=0.001 | record_interval=50
step,time,div_norm,r_u,r_v,ke,cfl_max,pcg_iters
0,0.000,0.0000e+00,-,-,0.0000e+00,0.000,0
50,0.050,2.3e-04,1.2e-02,8.3e-03,1.23e-02,0.024,12
100,0.100,1.8e-04,6.1e-03,4.2e-03,1.45e-02,0.026,10
```

说明：
- `-` 表示该步无有效值（如第 0 步无上一步用于计算 R_u）
- `pcg_iters` 对 Chorin（GS 固定 50 次）填 50；对 PCG/MAC 填实际 PCG 次数

### 3.2 稳态监控 CSV（SIMPLE `--format monitor`）

```
# FlowLabLite monitor | solver=simple | Re=20
iter,time,div_norm,r_u,r_v,ke,cfl_max,pcg_iters
0,0.000,0.0000e+00,-,-,0.0000e+00,-,0
10,0.010,2.3e-03,5.6e-02,4.3e-02,1.12e-02,-,50
```

说明：SIMPLE 无时间概念，time = iter × dt（伪时间）；CFL 对稳态填 `-`。

---

## 4. MonitorRecord 数据结构

```moonbit
// cmd/main/monitor.mbt
struct MonitorRecord {
  step      : Int     // 时间步编号（瞬态）或迭代步（稳态）
  time      : Double  // 物理时间 = step × dt（稳态为伪时间）
  div_norm  : Double  // 散度范数（质量守恒指标）
  r_u       : Double  // u 方向速度变化残差（-1.0 表示无效/不适用）
  r_v       : Double  // v 方向速度变化残差（-1.0 表示无效/不适用）
  ke        : Double  // 动能
  cfl_max   : Double  // CFL 最大值（稳态填 -1.0）
  pcg_iters : Int     // 压力迭代次数（GS 固定 50；PCG 填实际次数）
}
```

---

## 5. 实现隔离原则

**所有 monitor.mbt 中的函数必须通过参数接受数组，不直接读取全局变量**（`g_u`、`g_v` 等）。

**理由**：
1. **可测试性**：用 `make_test_export()` 合成数据直接单元测试，无需运行完整模拟
2. **解耦**：全局变量名变更不影响监控函数
3. **复用性**：`compute_cfl_max(u, v, dt, dx, dy)` 对四套求解器完全通用
4. **无副作用**：函数不持有状态，不可能意外修改数值场
5. **OpenFOAM 类比**：`functionObjects` 接收 `fvMesh&` 引用而非访问求解器内部状态

调用点 `io_cli.mbt` 负责从全局数组取值并传入监控函数。

---

## 6. 当前实现状态（截至 2026-04-18）

| 监控量 | Chorin | SIMPLE | PCG | MAC | 实现步骤 |
|---|---|---|---|---|---|
| div_norm | ✓ | ✓ | ✓ | ✓ | 已有 |
| p 求解残差/iters | GS(固定 50) | GS(固定 50) | ✓ iters | ✓ iters | 已有 |
| SIMPLE 简化残差 | N/A | ✓(=div) | N/A | N/A | 已有 |
| SIMPLE 方程残差 | N/A | ✗ | N/A | N/A | **Step M6** |
| R_u, R_v（速度变化） | ✗ | ✗ | ✗ | ✗ | **Step M1** |
| CFL_max | ✗ | N/A | ✗ | ✗ | **Step M2** |
| KE 动能 | ✗ | ✗ | ✗ | ✗ | **Step M3** |
| 监控历史 CSV | ✗ | ✗ | ✗ | ✗ | **Step M5** |
| 中心线速度 | ✓(WASM) | ✓ | ✓ | ✓ | 已有 |
| Ghia L2 误差 | ✗ | ✗ | ✗ | ✗ | **Step D** |

---

## 7. 参考权威代码

| 代码 | 求解器类型 | 相关监控 |
|---|---|---|
| OpenFOAM simpleFoam | 稳态 SIMPLE | `residualControl {U 1e-4; p 1e-4;}`，continuityErrors |
| OpenFOAM pimpleFoam | 瞬态 PIMPLE | 同上 + 时间步 Co 数，CFL 自适应 dt |
| Fluent | 稳态/瞬态 | Residuals Monitor（continuity, x/y-velocity）|
| CFX | 稳态/瞬态 | MAX/RMS residuals，全局质量不平衡 < 1e-4 |

---

## 8. 稳定性参考（Re 扩展）

| Re | ν | CFL_max (dt=0.001, U≈1) | 建议步数 |
|---|---|---|---|
| 20 | 0.100 | ≈ 0.02 ✓ | 500 |
| 100 | 0.020 | ≈ 0.02 ✓ | 2000 |
| 400 | 0.005 | ≈ 0.02 ✓ | 5000 |
| 1000 | 0.002 | ≈ 0.02 ✓ | 10000+ |

dt=0.001，dx=0.05 时，黏性稳定条件 dt ≤ dx²/(4ν) 在所有以上 Re 下均满足。
