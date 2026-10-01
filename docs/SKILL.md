# FlowLabLite — SKILL.md
# 重要决策、Bug 解决过程与架构思想

版本：v0.1 · 2026-04-18
作者：Fenfen Yu（余芬芬）

本文档记录项目开发过程中的非显而易见决策、重要 Bug 的根因分析和修复过程、架构思想。
供未来维护者和 AI 助手参考，避免重蹈覆辙。

---

## 1. 数值算法关键决策

### 1.1 GAMG Richardson 迭代符号：减而非加

**问题背景**：实现 GAMG V-cycle 外循环时，更新公式写为 `p += z`（加法），导致求解器发散。

**根因分析**：
GAMG V-cycle 求解 B·z = r，其中 B = −A（A 为 SND Laplacian）。
更新后误差递推为 `e_{k+1} = (I + T·A)·e_k`，谱半径 = 2 → 发散。

**正确公式**：
```moonbit
p[i][j] = p[i][j] - z[i][j]   // 正确：谱半径 → 0，收敛
// p[i][j] = p[i][j] + z[i][j]   // 错误：谱半径 = 2，发散
```

数学证明：T = (−A)⁻¹，误差递推 e_{k+1} = (I + T·A)·e_k = (I − I)·e_k = 0。

**教训**：更新符号取决于预条件器求解的是 A·z = r 还是 (−A)·z = r，必须与 Laplacian 符号约定一致。

---

### 1.2 GAMG 限制算子：必须除以 4

**问题背景**：粗网格修正量 4× 放大，导致 α ≈ −0.25，求解器停滞。

**根因分析**：
`gamg_restrict` 使用权重注入（1, 0.5, 0.25）。内部节点原始权重之和 = 4，
若不除以 4，粗网格 RHS 被 4× 放大，修正量过大。

**修复**：
```moonbit
// 所有注入权重相加后，输出值除以 4
coarse_b[ci][cj] = sum_of_weighted_values / 4.0
```

**教训**：任何加权平均算子都必须验证权重归一化（总和 = 1）。

---

### 1.3 GAMG 粗网格求解器：CG 而非阻尼 Jacobi

**问题背景**：粗网格用阻尼 Jacobi（谱半径 ≈ 1）无法消除低频误差，GAMG 效率极低。

**决策**：粗网格（21×21）使用共轭梯度法（CG）直接求解。
- Jacobi 对光滑（低频）模式谱半径接近 1，限制算子传递的正是这些模式
- CG 在小系统上 O(n) 次迭代收敛，21×21 系统约 40 次足够

**教训**：多重网格的粗网格求解器必须是"精确"求解器（或足够精确的迭代法），
不能用同一个光滑器（Jacobi）——这是多重网格设计的基本原则。

---

### 1.4 GAMG 外循环：稳态 Richardson 迭代，而非 PCG

**问题背景**：尝试将 GAMG V-cycle 作为 PCG 的预条件器，导致 PCG 正交性被破坏。

**根因分析**：
PCG 要求预条件器 M 对称正定且**固定**（线性算子）。
GAMG V-cycle 内的粗网格 CG 迭代次数可变（非线性），使 M 成为可变算子，
破坏 PCG 的共轭性 → 收敛失败。

**决策**：`pressure_poisson_pcg_gamg` 实现为**稳态 Richardson 迭代**（固定点迭代），
而非 Flexible-PCG。稳态迭代无需 M 对称，收敛可靠。

**教训**：可变预条件器（GAMG、AMG）只能用于 Flexible-CG 或稳态迭代，不能直接用于标准 PCG。

---

### 1.5 SND 符号约定：压力对正 RHS 为负值

**问题背景**：测试断言 `p[i][j] > 0`，但实际解为负，导致误报为 bug。

**根因分析**：
离散 Laplacian A 是对称负定（SND）矩阵（对角线 = −a_diag = −1600 < 0）。
A·p = b，b > 0 → p < 0（数学必然）。

**决策**：不写"p > 0"的测试；测试应验证物理现象（如顶盖附近存在高压/低压区的相对关系）。

**教训**：数值求解器的符号约定必须从离散矩阵出发分析，不能依赖物理直觉（物理上压力可正可负，但离散解的符号取决于矩阵定义）。

---

## 2. MoonBit 语言特性与陷阱

### 2.1 不支持行首 `+` 续行

**问题**：多行表达式在新行以 `+` 开头，编译报错。

**正确写法**：
```moonbit
// 错误
let val = a * b
        + c * d   // ← 编译错误

// 正确
let part1 = a * b
let part2 = c * d
let val = part1 + part2
```

---

### 2.2 Double 没有超越函数（sin/cos/exp/tanh 等）

**问题**：`Double.sin()` / `Double.exp()` / `Double.tanh()` 等均不存在（编译错误 [4015]）。
MoonBit 标准库（moonbitlang/core）未提供浮点超越函数。

**工程影响**：
- 测试 RHS 无法用 sin 构造 → 改用多项式：`(fi*(ny1-fi)) * (fj*(nx1-fj))`
- tanh 拉伸格式（Step 9）需自行实现 exp/tanh

**解决方案 A（测试场）**：用多项式构造光滑、边界为零的 RHS：
```moonbit
b[i][j] = (fi * (ny1 - fi)) * (fj * (nx1 - fj))   // 光滑，边界为零
```

**解决方案 B（生产代码，Step 9）**：在同文件内实现私有 Taylor 级数近似：
```moonbit
fn exp_impl(x : Double) -> Double {   // 25 项，|x|≤20 误差 < 1e-14
  if x < 0.0 { return 1.0 / exp_impl(-x) }
  let mut term = 1.0; let mut s = 1.0
  for k = 1; k <= 25; k = k + 1 { term = term * x / k.to_double(); s = s + term }
  s
}
fn tanh_impl(x : Double) -> Double {
  if x >= 20.0 { return 1.0 }; if x <= -20.0 { return -1.0 }
  let e2x = exp_impl(2.0 * x); (e2x - 1.0) / (e2x + 1.0)
}
```
未来若 MoonBit 引入 `@math` 包，只需删除 `_impl` 函数并替换调用处。

---

### 2.3 MoonBit `Double` 渲染整数为无小数点字符串

**问题**：`0.0.to_string()` 返回 `"0"`（而非 `"0.0"`）；`0.05.to_string()` 返回 `"0.05"`。
VTK/Tecplot 解析器均接受无小数点的整数形式，不影响 ParaView 读取。

**测试影响**：T90 原测试检查坐标字符串 `"0.0 0.0 0.0"` 失败（实际为 `"0 0 0.0"`）。
修正为测试结构关键字（POINT_DATA、SCALARS、VECTORS）而非坐标值。

---

### 2.4 全局可变状态模拟

MoonBit 没有可变全局变量。用单元素数组模拟：
```moonbit
let g_quiet_mode    : Array[Bool]   = [false]   // g_quiet_mode[0] = true 开启安静模式
let g_nu_runtime    : Array[Double] = [nu]       // Step C：运行时 Re 覆盖
let g_steps_override: Array[Int]    = [0]        // Step C：运行时步数覆盖（0=使用默认）
```

---

### 2.5 bare `return` 在 Unit 函数中合法

经验证，MoonBit Unit 函数中可用 `return` 提前退出，无需任何特殊语法。

---

### 2.6 大写字母标识符被解析为类型名（Step J 新增，2026-04-19）

**问题**：MoonBit 的标识符解析规则与大多数语言不同——**大写字母开头的标识符**在值上下文中被解析为**类型构造器**，不能作为 `let` 值绑定。

```moonbit
// 错误：L 被当作类型名，编译报 Error [4021]: L is unbound
fn rhs_u_dt(...) -> Array[Array[Double]] {
  let L = create_zeros_2d(nr, nc)   // ← 编译器以为这是类型声明上下文
  L[i][j] = ...                     // ← L 找不到值绑定
  L                                 // ← 同上
}
```

**正确做法**：所有值绑定变量必须以**小写字母**或 `_` 前缀开头。
对应数学符号的替换约定：

| 数学符号（大写）| MoonBit 变量名 |
|---|---|
| `L(u)·Δt`（RHS 算子）| `rhs`, `rh1u`, `rh2u` |
| `A`（系数矩阵）| `mat_a`, `coef` |
| `B`（右端向量）| `rhs_b`, `bvec` |
| `U`, `V`, `P`（场量）| `vel_u`, `vel_v`, `press` |

**影响范围**：所有数值 CFD 函数中以数学约定命名的单字母或短大写变量，均需改名。

---

### 2.7 `for _` 通配符不能作为循环计数器（Step J 新增，2026-04-19）

**问题**：`_` 是 MoonBit 的通配符模式，只适用于模式匹配。在 C 风格 for 循环中用作计数器会触发语法混淆（编译器尝试将其解析为 `for _ in ...` 函数式迭代）。

```moonbit
// 错误：Error [3002] Parse error, unexpected token `=`
for _ = 0; _ < n; _ = _ + 1 { ... }

// 正确：以 _ 为前缀的合法标识符
for _k = 0; _k < n; _k = _k + 1 { ... }
```

**规则**：需要忽略循环计数器时，使用 `_k`、`_i`、`_step` 等带 `_` 前缀的标识符，而非裸 `_`。

---

## 3. 架构设计决策

### 3.1 监控函数隔离原则（Pure Function Monitor）

**决策**：monitor.mbt 中所有监控函数通过参数接收数组，不直接读取全局变量（`g_u`, `g_v` 等）。

**理由**：
1. **可测试性**：合成数组直接单元测试，无需运行完整模拟
2. **解耦**：全局变量重命名不影响监控函数
3. **复用性**：`compute_cfl_max(u, v, dt, dx, dy)` 对四套求解器通用
4. **无副作用**：纯函数不可能意外修改数值场
5. **OpenFOAM 类比**：`functionObjects` 接收 `fvMesh&` 而非访问求解器内部

---

### 3.2 格式输出策略：计算在 MoonBit，展示在 JS/HTML

**决策**：所有数值计算、格式序列化（VTK/CSV/Tecplot/HTML 生成）在 MoonBit 中实现；
HTML/JS 只负责触发（按钮点击）和展示（绘图）。

**理由**：
- 保证单一事实来源（MoonBit 是唯一计算层）
- JS 格式生成代码难以测试（无 `moon test` 覆盖）
- MoonBit → WASM 导出字符串格式函数，JS 通过 WASM 调用获取完整字符串，触发下载

**注意**：当前 local_viewer.html 部分下载按钮由 JS 生成（历史遗留），Step M 后逐步迁移到 WASM 调用。

---

### 3.3 安静模式（Quiet Mode）设计

**问题**：`--format vtk` 输出时，求解器进度 println 污染 stdout，破坏重定向。

**方案**：`g_quiet_mode[0] = true` 在 `--format` 触发时设置，抑制三处进度输出（Chorin/PCG/MAC）。
格式输出走 stdout，进度可通过 stderr 输出（未来改进方向）。

---

### 3.4 NetCDF 降级策略（Tier-1/Tier-2 格式）

**问题**：libnetcdf 不是所有环境都有，编译失败会阻止其他格式使用。

**方案**：`io_netcdf_stub.mbt` 永远编译通过，提供相同函数签名但返回 `Err`。
`--list-formats` 用 ✓/✗ 标注可用性。升级到真实 NetCDF 时替换该文件，其余代码不变。

---

### 3.5 边界条件顺序：侧壁优先，顶盖最后

**问题**：若顶盖 BC 先设，侧壁循环会把角点 `u[ny-1][0]` 和 `u[ny-1][nx-1]` 覆盖为 0。

**正确顺序**（已固化为项目约定，在 CLAUDE.md 中也有记录）：
```moonbit
// 1. 侧壁
for i = 0; i < ny; i = i + 1 { u[i][0] = 0.0; u[i][nx-1] = 0.0 }
// 2. 顶盖（最后）
for j = 0; j < nx; j = j + 1 { u[ny-1][j] = 1.0 }
```

---

### 3.6 SIMPLE 残差：简化版 vs 严格版

**历史**：SIMPLE 的 `get_simple_residual()` 直接返回 div_norm（散度范数），是简化实现。
与 OpenFOAM 的归一化方程残差含义不同。

**决策（Step M6）**：提供两种选项：
- `--simple-residual div`：简化版（默认，=div_norm，已有）
- `--simple-residual eq`：严格版（归一化方程残差，Step M6 新增）

用户可通过 CLI 选择，两种都是只读 post-processing，不影响数值路径。

---

## 4. 已知问题与处理方式

### 4.1 moon 0.1.20260309 不传 `-exported_functions`

**现象**：`moon build` 生成的 WASM 只有 `_start`，浏览器无法调用求解器函数。

**处理**：`build_wasm.sh` 在 moon build 后重新调用链接器（`wasm-ld` / `wasm-gc-ld`），
显式传入所有 66 个导出函数名。

**状态**：持续监控 moon 版本更新，一旦修复则移除 workaround。

---

### 4.2 Tecplot 多 ZONE 顺序

**现象**：输出 4 个 ZONE（Chorin/SIMPLE/PCG/MAC），其中 MAC 使用 nc=40，其余 nc=41。
`I=40 J=40` 与其他 ZONE 的 `I=41 J=41` 不同。

**处理**：`format_tecplot` 中每个 ZONE 使用 `s.nc` 而非硬编码 41，动态适配。

---

*本文档在重要 Bug 修复、架构变更后更新。*
