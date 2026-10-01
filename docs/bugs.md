# FlowLabLite — Bug 记录文档

**项目**: FlowLabLite  
**作者**: Fenfen Yu（余芬芬）  
**日期**: 2026-04-17  
**覆盖范围**: 从初始提交（2026-03-09）至当前 HEAD（ba1505e）

---

## Bug 类型分类

| 类型代码 | 类型名称 | 说明 |
|---|---|---|
| **SYN** | 语言语法错误 | 违反 MoonBit 语法规则，编译报错 |
| **SEM** | 语言语义错误 | 语法合法但语义错误，编译期类型/约束报错 |
| **API** | API 误用 | 错误使用语言标准库或工具链 API |
| **CFG** | 配置/构建错误 | 构建工具、包配置文件问题 |
| **RUN** | 运行时错误 | 程序运行时崩溃或越界 |
| **NUM** | 数值逻辑错误 | 公式、算法、符号错误，产生错误的计算结果 |
| **TST** | 测试设计错误 | 测试断言逻辑错误，不能正确验证目标行为 |
| **INT** | 集成/接口错误 | 跨模块/跨语言（MoonBit↔JS↔WASM）接口不匹配 |

严重程度：**P1**（崩溃/完全不可用）、**P2**（功能错误/数值错误）、**P3**（警告/次要缺陷）

---

## Bug #1 — `func` 关键字在 MoonBit 中无效

| 属性 | 内容 |
|---|---|
| **类型** | SYN |
| **严重程度** | P1 |
| **发现时机** | 初始开发阶段，编译报错 |
| **问题提交** | 早期 class-based 实现（`9a92585` 之前） |
| **修复提交** | `9c7a9fa` fix(main): fix the func identifier to fn. |

**Bug 现象**

初始代码从 Python/OpenFOAM 参考实现翻译时，使用了 `func` 作为函数定义关键字：

```moonbit
// 错误：func 不是 MoonBit 关键字
func init(nX : int, nY : int, length : float) { ... }
```

编译器报错，程序完全无法构建。

**定位与解决**

MoonBit 使用 `fn` 定义函数，而非 `func`（Go 风格）或 `def`（Python 风格）。

```moonbit
// 修复：
fn init(nX : int, nY : int, length : float) { ... }
```

---

## Bug #2 — 类型标识符使用错误大小写（`int`/`float` vs `Int`/`Double`）

| 属性 | 内容 |
|---|---|
| **类型** | SYN / SEM |
| **严重程度** | P1 |
| **发现时机** | 初始开发阶段，编译报错 |
| **问题提交** | 早期实现 |
| **修复提交** | `f98b514` fix: fix the identifier case；`1db1c0d` fix: fix the identifier case Float |

**Bug 现象**

从参考代码（Python/C++/OpenFOAM 风格）直译时，使用了小写 `int`/`float`，MoonBit 编译器找不到这些类型：

```moonbit
var Nx : int   // 错误：MoonBit 内置类型首字母大写
var L  : float // 错误：MoonBit 中浮点数是 Double，无 Float
```

**定位与解决**

MoonBit 基本类型：整数用 `Int`，浮点用 `Double`（无 `Float`）。

```moonbit
var Nx : Int
var L  : Double
```

---

## Bug #3 — struct 字段错误使用 `mut` 修饰符

| 属性 | 内容 |
|---|---|
| **类型** | SEM |
| **严重程度** | P1 |
| **发现时机** | 初始开发阶段，编译报错 |
| **问题提交** | `a5261ca`–`a487ada` 阶段 |
| **修复提交** | `577e748`、`f819f18`、`efe153d`、`7018f81`、`6bf556e`（5 次连续修复） |

**Bug 现象**

在 struct 定义中为数组字段加了 `mut`，但 MoonBit 中 `Array[T]` 是引用类型，字段本身指向可变数组，不需要也不允许对字段本身标记 `mut`：

```moonbit
struct Field[T] {
  mut data : Array[Array[T]]  // 编译错误
}
struct CavitySolver {
  mut u : Field[Double]  // 编译错误
  mut v : Field[Double]  // 编译错误
  mut p : Field[Double]  // 编译错误
}
```

**定位与解决**

MoonBit 的 `Array[T]` 天然可变（内容可原地修改），字段声明不需要 `mut`。连续修复 5 个 struct 的字段定义，去除多余的 `mut`：

```moonbit
struct Field[T] {
  data : Array[Array[T]]  // 正确：数组内容可变，字段指针不可变
}
```

---

## Bug #4 — 2D 数组行共享问题（`Array::make` 浅拷贝陷阱）

| 属性 | 内容 |
|---|---|
| **类型** | SEM / RUN |
| **严重程度** | P2 |
| **发现时机** | 运行时，所有行数据相同 |
| **问题提交** | `8d12a24` 附近（添加网格函数时） |
| **修复提交** | `8453218` fix: The value identifiers X, Y are unbound（含修复行共享问题） |

**Bug 现象**

使用 `Array::make(ny, Array::make(nx, 0.0))` 创建 2D 数组时，所有行指向同一个内部数组（浅拷贝），修改任意一行等效于修改所有行：

```moonbit
// 错误：ny 行全部共享同一个 Array::make(nx, 0.0) 实例
let X = Array::make(ny, Array::make(nx, 0.0))
```

**定位与解决**

需要在循环中为每行分别分配：

```moonbit
// 正确：每行独立分配
let x = Array::make(ny, Array::make(nx, 0.0))
for i = 0; i < ny; i = i + 1 {
  x[i] = Array::make(nx, 0.0)   // 重新分配，覆盖浅拷贝行
  for j = 0; j < nx; j = j + 1 {
    x[i][j] = x_coords[j]
  }
}
```

---

## Bug #5 — 全局 `let mut` 变量不被允许

| 属性 | 内容 |
|---|---|
| **类型** | SEM |
| **严重程度** | P1 |
| **发现时机** | 编译报错 |
| **问题提交** | `4aa1f77`（添加步骤计数器时） |
| **修复提交** | `0e26bb7` fix(mbt): replace let mut global with Array[Int] counter |

**Bug 现象**

在模块顶层声明可变全局变量失败：

```moonbit
// 错误：MoonBit 全局作用域不允许 let mut
let mut g_steps_done : Int = 0
```

**定位与解决**

MoonBit 全局变量不支持 `let mut`（区别于局部作用域）。解决方案：使用单元素数组作为可变容器：

```moonbit
// 修复：用 Array[Int] 绑定全局可变状态
let g_steps : Array[Int] = [0]
// 读取：g_steps[0]
// 写入：g_steps[0] = g_steps[0] + n
```

---

## Bug #6 — 类型推导失败（全局常量缺少类型标注）

| 属性 | 内容 |
|---|---|
| **类型** | SEM |
| **严重程度** | P1 |
| **发现时机** | 编译报错 `type cannot be inferred` |
| **问题提交** | `b1454af` 修复前的版本 |
| **修复提交** | `b1454af` fix: type cannot be inferred error |

**Bug 现象**

```moonbit
// 错误：编译器无法推导 dx/dy 的类型（Int 除法还是 Double 除法？）
let dx = 2.0 / (nx - 1).to_double()
let dy = 2.0 / (ny - 1).to_double()
```

**定位与解决**

在全局常量定义处加显式类型标注：

```moonbit
let dx : Double = 2.0 / (nx - 1).to_double()
let dy : Double = 2.0 / (ny - 1).to_double()
```

---

## Bug #7 — moon.pkg 旧格式导致 `UnexpectedToken` 构建失败

| 属性 | 内容 |
|---|---|
| **类型** | CFG |
| **严重程度** | P1 |
| **发现时机** | `moon build` 报 parse error |
| **问题提交** | 早期版本使用 DSL 格式的 `moon.pkg` |
| **修复提交** | `b40fb42` fix(lib): replace old-format moon.pkg with moon.pkg.json |

**Bug 现象**

`lib/moon.pkg` 使用了旧 DSL 格式（非 JSON），moon 工具链版本升级后不再识别，报 `UnexpectedToken` 错误，整个工程无法构建。

**定位与解决**

删除 `lib/moon.pkg`，替换为标准 JSON 格式 `lib/moon.pkg.json`：

```json
{
  "import": []
}
```

同期 `aab0ce5` 修复了 `bench-import` 字段引起的类似问题。

---

## Bug #8 — `moon build` 生成的 WASM 只导出 `_start`，无 API 函数

| 属性 | 内容 |
|---|---|
| **类型** | CFG / INT |
| **严重程度** | P1 |
| **发现时机** | 浏览器中 WASM 加载后所有 API 调用均为 `undefined` |
| **问题提交** | 所有早期版本 |
| **修复提交** | `5ce2af4` fix(wasm): add build_wasm.sh to export all 27 API functions |

**Bug 现象**

`moon build --target wasm-gc` 生成的 WASM 只包含 `_start` 入口，`init_simulation`、`run_n_steps` 等所有 API 函数均不在导出段。浏览器 JS 调用时得到 `undefined`。

**根本原因**

moon 0.1.20260309 向 `moonc link-core` 传递 `-pkg-config-path ./cmd/main/moon.pkg`，而旧格式 `moon.pkg`（及当时的 `moon.pkg.json`）在 link 阶段未正确传递 `-exported_functions`，导致链接器不知道要导出哪些函数。

**定位与解决**

新建 `build_wasm.sh`，在 `moon build` 完成核心编译后，手动调用 `moonc link-core` 并显式传入 `-exported_functions` 标志：

```bash
moonc link-core \
  "$CORE_FILE" \
  -main FlowLabLite/cmd/main \
  -exported_functions "$EXPORTS" \
  -target wasm-gc \
  -o "$OUT_FILE" \
  ...
```

---

## Bug #9 — wasm-gc 字符串不能直接跨 JS 读取

| 属性 | 内容 |
|---|---|
| **类型** | INT |
| **严重程度** | P2 |
| **发现时机** | JS 端调用字符串返回函数时读到乱码/崩溃 |
| **问题提交** | 早期设计了返回 String 的 WASM 导出（`get_simulation_config` 等） |
| **修复提交** | `1eef28b` fix(pkg): drop string-returning exports — wasm-gc strings are GC refs |

**Bug 现象**

wasm-gc 中 `String` 是托管 GC 引用（ref type），无法通过线性内存指针在 JS 端直接读取，调用后崩溃或返回无意义数据。

**定位与解决**

删除所有返回 `String` 的 WASM 导出函数（`get_simulation_config`、`run_cavity_simulation`、`run_simple_simulation`）。改为仅导出数值型 getter（`Int`、`Double`），所有字符串操作在 JS 端完成：

```moonbit
// 删除：pub fn get_simulation_config() -> String { ... }
// 改为：数值型 getter
pub fn get_nx() -> Int { nx }
pub fn get_re() -> Double { rho * 1.0 / nu }
```

---

## Bug #10 — 浏览器加载 WASM 时缺少 import object（`LinkError`）

| 属性 | 内容 |
|---|---|
| **类型** | INT |
| **严重程度** | P1 |
| **发现时机** | 浏览器控制台 `LinkError: WebAssembly.instantiate` 失败 |
| **问题提交** | 早期 HTML 直接调用 `WebAssembly.instantiate(wasmBytes)` 无参数 |
| **修复提交** | `20009d9` fix(html): provide Proxy-based MoonBit wasm-gc import object |

**Bug 现象**

MoonBit wasm-gc 编译产物需要 `spectest`、`env`、`moonbit:ffi` 等命名空间的导入，直接实例化时 `WebAssembly.instantiate` 因找不到这些导入而抛出 `LinkError`。

**定位与解决**

构造带有所需命名空间的 import object，使用 `Proxy` 将未知导入自动映射为 no-op：

```javascript
const noopFn = () => {};
const makeNoopModule = () => new Proxy({}, { get: () => noopFn });
const importObject = new Proxy({
    spectest: { print_char: noopFn, print_i32: noopFn, print_f64: noopFn },
    env: { println: (s) => console.log('[mbt]', s) },
    'moonbit:ffi': { println: (s) => console.log('[mbt]', s) },
}, { get: (t, k) => k in t ? t[k] : makeNoopModule() });
const wasmModule = await WebAssembly.instantiate(wasmBytes, importObject);
```

---

## Bug #11 — 多行表达式首行以 `+` 开头导致编译错误

| 属性 | 内容 |
|---|---|
| **类型** | SYN |
| **严重程度** | P1 |
| **发现时机** | 编写 Navier-Stokes 离散式时，编译报错 |
| **问题提交** | PCG 和 MAC 求解器开发过程中多次出现 |
| **修复提交** | 每次出现均在同一提交中修复（见 `2be1135`、`3b67401`） |

**Bug 现象**

翻译 CFD 多项物理量叠加公式时，自然地将长表达式拆行并以 `+` 开头：

```moonbit
// 错误：MoonBit 不支持行首 + 续行
let adv_u = -u * (u - u_w) / dx
           + (-v * (u - u_s) / dy)   // ← 编译报错
```

**定位与解决**

MoonBit 将行首 `+` 解析为正号（一元运算符），而非续行二元加法。改用中间变量：

```moonbit
// 正确：拆分为独立中间变量
let adv_u_x = -u * (u - u_w) / dx
let adv_u_y = -v * (u - u_s) / dy
let adv_u   = adv_u_x + adv_u_y
```

该规则已写入 `CLAUDE.md` 作为永久约定。

---

## Bug #12 — 边界条件顺序错误导致顶盖角点被覆盖

| 属性 | 内容 |
|---|---|
| **类型** | NUM |
| **严重程度** | P2 |
| **发现时机** | 调试边界条件时发现角点 u=0 而非 u=1 |
| **问题提交** | 早期边界条件实现 |
| **修复提交** | 在 `4424836`（大重写）中确立正确顺序，`CLAUDE.md` 中记录 |

**Bug 现象**

顶盖驱动腔流要求顶盖 `u[ny-1][j] = 1.0`。若先设顶盖再设侧壁，侧壁循环（`i` 从 0 到 ny-1）会将角点 `u[ny-1][0]` 和 `u[ny-1][nx-1]` 覆盖回 0：

```moonbit
// 错误顺序：侧壁循环覆盖了顶盖角点
for j = 0; j < nx; j = j + 1 { u[ny-1][j] = 1.0 }  // 先设顶盖
for i = 0; i < ny; i = i + 1 { u[i][0] = 0.0; u[i][nx-1] = 0.0 }  // 覆盖角点
```

**定位与解决**

侧壁 BC 必须**先于**顶盖 BC 设置，顶盖作为最后一道赋值：

```moonbit
// 正确顺序
for i = 0; i < ny; i = i + 1 { u[i][0] = 0.0; u[i][nx-1] = 0.0 }  // 侧壁先
for j = 0; j < nx; j = j + 1 { u[ny-1][j] = 1.0 }                  // 顶盖最后
```

---

## Bug #13 — GAMG V-cycle 外层迭代符号错误（`p += z` 发散）

| 属性 | 内容 |
|---|---|
| **类型** | NUM |
| **严重程度** | P1（求解器完全不收敛） |
| **发现时机** | T57 测试失败，rel_err ≥ 1%；调试发现迭代发散 |
| **问题提交** | GAMG 第一版实现（`06e8896` 之前的开发过程） |
| **修复提交** | `06e8896` feat(solver): DILU, DIC, and GAMG pressure preconditioners |

**Bug 现象**

GAMG V-cycle 输出 `z ≈ (−A)⁻¹·r`，外层 Richardson 迭代最初写成加法：

```moonbit
// 错误：p += z 导致谱半径 = 2，迭代发散
p[i][j] = p[i][j] + z[i][j]
```

每步误差翻倍，200 次迭代后压力场与参考解误差达几百倍，T57 失败。

**定位与解决**

数学推导：V-cycle 求解的系统是 `B·z = r`，其中 `B = −A`（SPD）。因此 `z ≈ (−A)⁻¹·r = −A⁻¹·r`。  
Richardson 迭代格式：`p_{k+1} = p_k − T·r_k`，其中 `T = (−A)⁻¹`。

误差递推：  
`e_{k+1} = (I + T·A)·e_k = (I + (−A)⁻¹·A)·e_k = (I − I)·e_k = 0`

因此应**减去** z：

```moonbit
// 正确：p -= z，谱半径 → 0
p[i][j] = p[i][j] - z[i][j]
```

---

## Bug #14 — GAMG 限制算子缺少 /4 归一化

| 属性 | 内容 |
|---|---|
| **类型** | NUM |
| **严重程度** | P2（求解器停滞，迭代数达 200 但不收敛） |
| **发现时机** | GAMG 收敛测试中发现粗网格修正过度放大，PCG 步长 α ≈ −0.25 |
| **问题提交** | GAMG 早期实现 |
| **修复提交** | `06e8896` feat(solver): DILU, DIC, and GAMG pressure preconditioners |

**Bug 现象**

`gamg_restrict` 使用注入加权平均（权重 1/2/4 邻域分别 1、0.5、0.25），对于内部节点原始权重之和 = 4，未除以 4 导致粗网格残差被放大 4 倍。PCG 中搜索方向步长计算为 α ≈ −1/4 而非 α ≈ −1，求解器停滞。

**定位与解决**

在限制算子中增加 `/4.0` 归一化：

```moonbit
fn gamg_restrict(r_fine, r_coarse) {
  // ...（构建加权和后）
  r_coarse[ci][cj] = weighted_sum / 4.0   // 必须除以 4
}
```

---

## Bug #15 — GAMG 粗网格使用阻尼 Jacobi（谱半径 ≈ 1，无法消除低频误差）

| 属性 | 内容 |
|---|---|
| **类型** | NUM |
| **严重程度** | P2（收敛极慢） |
| **发现时机** | GAMG 粗网格迭代后低频模式误差未减少 |
| **问题提交** | GAMG 中间开发版本 |
| **修复提交** | `06e8896` feat(solver): DILU, DIC, and GAMG pressure preconditioners |

**Bug 现象**

粗网格阻尼 Jacobi（ω = 2/3）对光滑（低频）模式谱半径约为 1，无法将误差从粗网格向下传递，造成 V-cycle 整体无效。

**定位与解决**

改用 CG（共轭梯度）在 21×21 SPD 系统（`−A` 在粗网格上）上精确求解：

```moonbit
fn gamg_coarse_cg(e_c, r_c, inv_dx2_c, inv_dy2_c, a_diag_c) {
  // CG 迭代直到收敛，而非固定次数的 Jacobi
}
```

---

## Bug #16 — GAMG 以 PCG 作为外层求解器（非线性算子破坏正交性）

| 属性 | 内容 |
|---|---|
| **类型** | NUM |
| **严重程度** | P2（偶发不收敛） |
| **发现时机** | 将 V-cycle 作为 PCG 预条件器时，某些 RHS 下正交条件被破坏 |
| **问题提交** | 早期方案：`pressure_poisson_pcg_prec(p, ..., precond_type=3)` 调用标准 PCG |
| **修复提交** | `06e8896` feat(solver): DILU, DIC, and GAMG pressure preconditioners |

**Bug 现象**

标准 PCG 要求预条件器为对称正定的**线性**算子。V-cycle 内部的 CG 粗网格求解迭代次数随残差变化（非线性），导致 V-cycle 本身不是线性算子，PCG 的 Lanczos 正交性条件被破坏，出现偶发不收敛。

**定位与解决**

改为**平稳 Richardson 迭代**（stationary Richardson iteration），不依赖预条件器的线性性和对称性：

```moonbit
fn pressure_poisson_pcg_gamg(p, dx_, dy_, b) -> Int {
  // 外层不用 PCG，改用平稳 Richardson
  for iter = 0; iter < pcg_max_iter; iter = iter + 1 {
    // r = b - A*p
    apply_gamg_precond(r, z, ...)   // z ≈ (−A)⁻¹·r
    // p = p - z  （减号，见 Bug #13）
  }
}
```

---

## Bug #17 — T57 断言"存在正压力节点"违反 SND 符号约定

| 属性 | 内容 |
|---|---|
| **类型** | TST |
| **严重程度** | P2（测试本身逻辑错误，即使求解器正确也必然失败） |
| **发现时机** | T57 测试失败分析 |
| **问题提交** | 初始 T57 测试设计 |
| **修复提交** | `06e8896` feat(solver): DILU, DIC, and GAMG pressure preconditioners |

**Bug 现象**

T57 原始断言：

```moonbit
let any_positive = // 检查是否有 p[i][j] > 0
assert_true(any_positive)   // 期望某些压力为正
assert_true(iters_gamg < iters_jac)  // 期望 GAMG 比 Jacobi 快
```

两个断言均错误：
1. 对于 SND 离散 Laplacian `A`，`b > 0` 时 `A·p = b` 的解 `p* < 0` 处处成立（见 decides.md #1）。
2. 平稳 Richardson 迭代与 PCG 迭代数无可比性。

**定位与解决**

改为相对误差检验，与 T45（DILU）、T51（DIC）保持一致：

```moonbit
let rel_err = diff_sq.sqrt() / ref_sq.sqrt()
assert_true(rel_err < 0.01)
```

---

## Bug #18 — `approx_eq(a, b, 0.0)` 永远返回 false

| 属性 | 内容 |
|---|---|
| **类型** | TST |
| **严重程度** | P2（T77/T78 确定性测试总是失败） |
| **发现时机** | T77 测试失败 |
| **问题提交** | T77/T78 初始实现 |
| **修复提交** | `f6084ac` test: add 25 extended white-box tests (T59–T83) |

**Bug 现象**

```moonbit
// 错误：approx_eq 实现为 (a - b).abs() < tol
// tol = 0.0 时，即使 a == b，0.0 < 0.0 也为 false
assert_true(approx_eq(u1, u2, 0.0))
```

**定位与解决**

对于确定性验证（要求**精确相等**），直接使用 MoonBit 的 `==` 运算符：

```moonbit
assert_true(u1 == u2)   // Double 精确相等判断
```

---

## Bug #19 — T79 断言 SIMPLE 残差单调递减（实际非单调）

| 属性 | 内容 |
|---|---|
| **类型** | TST |
| **严重程度** | P2（测试预期不符合物理现实） |
| **发现时机** | T79 失败，`r3 > r2 * 1.05` |
| **问题提交** | T79 初始实现 |
| **修复提交** | `f6084ac` test: add 25 extended white-box tests (T59–T83) |

**Bug 现象**

SIMPLE 求解器的散度范数（质量守恒残差）并非单调递减——随着流速增大，单步绝对残差会上升，整体呈"先升后降"模式：

```moonbit
// 错误假设：散度范数单调递减
assert_true(r3 <= r2 * 1.05)   // 实测 r3 > r2 * 1.05
```

**定位与解决**

改为验证最终残差小于绝对阈值（100 步后质量守恒到一定精度）：

```moonbit
// 正确：验证收敛到阈值，不要求单调
let div_100 = get_simple_divergence_norm()
assert_true(div_100 < 0.02)   // 实测约 0.007
```

---

## Bug #20 — MoonBit `Double` 没有 `.sin()` 方法

| 属性 | 内容 |
|---|---|
| **类型** | API |
| **严重程度** | P1（编译报错 [4015]） |
| **发现时机** | T75 初始实现，编译报 method not found |
| **问题提交** | T75 初始实现 |
| **修复提交** | `f6084ac` test: add 25 extended white-box tests (T59–T83) |

**Bug 现象**

```moonbit
// 错误：MoonBit Double 类型无 .sin() 方法
let bval = (pi * i.to_double() / (ny - 1).to_double()).sin()
// 编译错误 [4015]: Value sin not found in type Double
```

**定位与解决**

使用多项式代替三角函数生成非均匀、边界为零的测试 RHS：

```moonbit
// 正确：多项式 RHS，光滑、在边界为零
let fi = i.to_double()
let fj = j.to_double()
b[i][j] = (fi * (ny1 - fi)) * (fj * (nx1 - fj))
```

---

## Bug #21 — debug `println` 语句残留在测试中

| 属性 | 内容 |
|---|---|
| **类型** | TST |
| **严重程度** | P3（测试通过但输出噪声） |
| **发现时机** | 测试输出中发现 `DEBUG T57:` 行 |
| **问题提交** | T57 调试过程中临时加入 |
| **修复提交** | `06e8896` feat(solver): DILU, DIC, and GAMG pressure preconditioners |

**Bug 现象**

```moonbit
println("DEBUG T57: iters_gamg=\{iters_gamg}, iters_jac=\{iters_jac}")
```

测试通过后 println 未及时删除，污染测试输出。

**定位与解决**

删除所有调试输出语句，同时删除 `pressure_poisson_pcg_prec` 中 `precond_type==3` 的两段 println 调试块。

---

## 汇总表

| # | 类型 | 严重程度 | 简述 | 修复提交 |
|---|---|---|---|---|
| 1 | SYN | P1 | `func` 应为 `fn` | `9c7a9fa` |
| 2 | SYN/SEM | P1 | `int`/`float` 应为 `Int`/`Double` | `f98b514`、`1db1c0d` |
| 3 | SEM | P1 | struct 字段错误使用 `mut` | `577e748`–`6bf556e` |
| 4 | SEM/RUN | P2 | 2D 数组行共享浅拷贝 | `8453218` |
| 5 | SEM | P1 | 全局不支持 `let mut` | `0e26bb7` |
| 6 | SEM | P1 | 全局常量缺类型标注 | `b1454af` |
| 7 | CFG | P1 | `moon.pkg` 旧格式 | `b40fb42` |
| 8 | CFG/INT | P1 | WASM 仅导出 `_start` | `5ce2af4` |
| 9 | INT | P2 | wasm-gc String 不可跨 JS 读取 | `1eef28b` |
| 10 | INT | P1 | WASM 实例化缺 import object | `20009d9` |
| 11 | SYN | P1 | 多行表达式行首 `+` 无效 | `2be1135`、`3b67401` |
| 12 | NUM | P2 | BC 顺序错误，角点被覆盖 | `4424836` |
| 13 | NUM | P1 | GAMG `p += z` 发散 | `06e8896` |
| 14 | NUM | P2 | GAMG 限制算子缺 /4 | `06e8896` |
| 15 | NUM | P2 | GAMG 粗网格 Jacobi 谱半径≈1 | `06e8896` |
| 16 | NUM | P2 | GAMG 以 PCG 外层破坏正交性 | `06e8896` |
| 17 | TST | P2 | T57 断言正压力违反 SND 约定 | `06e8896` |
| 18 | TST | P2 | `approx_eq(a,b,0.0)` 永远 false | `f6084ac` |
| 19 | TST | P2 | SIMPLE 残差假设单调递减 | `f6084ac` |
| 20 | API | P1 | `Double.sin()` 不存在 | `f6084ac` |
| 21 | TST | P3 | debug println 残留 | `06e8896` |
| 22 | INT | P1 | `--validate-ghia` 被 `fmt=="json"` 提前拦截 | todo2 Step D |
| 23 | SEM | P1 | `MonitorRecord` 在 monitor.mbt 重复定义 | todo2 Step M |
| 24 | CFG | P1 | Java benchmark GBK 编码错误（需 -encoding UTF-8）| todo2 bench |

---

## Bug #22 — `--validate-ghia` 被 `fmt=="json"` 提前拦截

| 属性 | 内容 |
|---|---|
| **类型** | INT |
| **严重程度** | P1 |
| **发现时机** | todo2 Step D 实施，功能无输出 |
| **修复提交** | todo2 Step D |

**Bug 现象**

`dispatch_format_cli` 在 `--validate-ghia` 检测之前先检查 `fmt == "json"`（无 `--format` 时默认值），提前返回导致 Ghia 验证功能永远不执行。

**修复**

将 `--validate-ghia` 检测逻辑移至 `fmt == "json"` 提前返回之前。

---

## Bug #23 — `MonitorRecord` 在 monitor.mbt 重复定义

| 属性 | 内容 |
|---|---|
| **类型** | SEM |
| **严重程度** | P1 |
| **发现时机** | todo2 Step M 实施，编译失败 |
| **修复提交** | todo2 Step M |

**Bug 现象**

Step M 新建 `monitor.mbt` 时再次定义了 `MonitorRecord` struct 和 `format_monitor_csv`，而 Step B 已在 `io_formats.mbt` 中定义，导致编译错误（重复定义）。

**修复**

从 `monitor.mbt` 删除重复定义，添加注释指向 `io_formats.mbt`。

---

## Bug #24 — Java benchmark GBK 编码错误

| 属性 | 内容 |
|---|---|
| **类型** | CFG |
| **严重程度** | P1 |
| **发现时机** | todo2 bench 阶段，javac 编译失败 |
| **修复提交** | todo2 bench |

**Bug 现象**

Windows 下 `javac ChorinGSBench.java` 报 `unmappable character (0x97) for encoding GBK`，因 Java 文件注释含 `—`、`×` 等 Unicode 字符。

**修复**

```bash
javac -encoding UTF-8 bench/ChorinGSBench.java -d bench/
```

`bench/run_bench.sh` 已更新，所有 javac 调用添加 `-encoding UTF-8`。

---

## Bug #25 — MoonBit 大写单字母变量被解析为类型名

| 属性 | 内容 |
|---|---|
| **类型** | MoonBit 语法陷阱 |
| **严重程度** | P1（编译失败） |
| **发现时机** | Step J SSP-RK3 实现，`rhs_u_dt` 函数 |
| **修复提交** | feat(solver): add SSP-RK3 time integration (Step J) |

**Bug 现象**

```moonbit
fn rhs_u_dt(...) -> Array[Array[Double]] {
  let L = create_zeros_2d(nr, nc)   // ← 编译报 Error [4021]: L is unbound
  ...
  L[i][j] = ...
  L
}
```

编译器报 `Error [4021]: The value identifier L is unbound`，同理 `L1u`, `L2u`, `L3u` 也全部报错。

**根因**

MoonBit 将单个大写字母或以大写字母开头的短标识符识别为**类型构造器**，而非值绑定。
`let L = ...` 中，`L` 被解析为类型上下文，随后 `L[i][j] = ...` 找不到对应的值绑定。

**修复**

将 `L` → `rhs`，`L1u/L1v` → `rh1u/rh1v`，`L2u/L2v` → `rh2u/rh2v`，`L3u/L3v` → `rh3u/rh3v`，全部改为小写开头标识符。

**教训**

MoonBit 中变量命名规范：
- **值绑定**：必须以**小写字母**或 `_` 开头
- **类型/构造器**：以**大写字母**开头
- 即使在局部 `let` 绑定中，大写开头也会被解析为类型名，导致后续同名引用找不到值

---

## Bug #26 — `for _ = 0` 循环变量语法无效

| 属性 | 内容 |
|---|---|
| **类型** | MoonBit 语法陷阱 |
| **严重程度** | P1（编译失败） |
| **发现时机** | Step J `run_rk3_n_steps` 函数 |
| **修复提交** | feat(solver): add SSP-RK3 time integration (Step J) |

**Bug 现象**

```moonbit
pub fn run_rk3_n_steps(n : Int) -> Unit {
  for _ = 0; _ < n; _ = _ + 1 {   // ← Parse error: unexpected token `=`
    rk3_step(...)
  }
}
```

编译器报 `Error [3002]: Parse error, unexpected token '=', you may expect 'in'`。

**根因**

MoonBit 的 C 风格 `for` 循环要求初始化变量必须是合法标识符。`_` 是**通配符模式**，不能作为 for 循环计数器使用——编译器将 `for _ = 0` 解析为 `for _ in 0`（函数式迭代语法），触发语法混淆。

**修复**

使用 `_k` 替代 `_`：

```moonbit
for _k = 0; _k < n; _k = _k + 1 {
  rk3_step(...)
}
```

`_k` 以下划线开头，是合法的"有意忽略"命名约定，同时满足 for 循环的标识符要求。

**教训**

MoonBit 中纯 `_` 只能用于**模式匹配**（如 `match x { _ => ... }`）和**参数占位**（如 `fn f(_ : Int) -> Unit`），不能作为 for 循环的计数变量。
需要忽略循环变量时，用 `_k`、`_i`、`_step` 等带下划线前缀的名称。

---

## Bug #27 — TVD 限流器变量名与 MoonBit 命名规则冲突（预防性记录）

| 属性 | 内容 |
|---|---|
| **类型** | MoonBit 语法约束（预防） |
| **严重程度** | P2（潜在问题） |
| **发现时机** | Step I TVD 实现代码审查 |

**背景**

TVD 限流器在数学文献中常以 `ψ(r)` 表示，辅助变量习惯命名为 `r`（平滑比）、`L`（限流值）。
在 MoonBit 中，单字母大写如 `L`、`R` 会被解析为类型名（同 Bug #25）；`r` 等小写单字母值绑定可用。

**处理**

Step I 已使用：
- `r` → `r`（可用，小写）
- 限流器返回值直接使用 `psi`
- 辅助数组 `L` → 未使用，直接内联计算

Step J 修复记录见 Bug #25。

---

## Bug #28 — tanh 坐标公式方向相反（Step 9，2026-04-19）

**现象**：T141 测试期望 `first_spacing < 0.5 * mean_spacing`（在 x=0 处聚簇），实际 `first_spacing ≈ 0.41 >> mean = 0.2`。

**根因**：初始公式 `x_i = L·tanh(β·ξ_i)/tanh(β)` 中，tanh 为次线性函数（tanh(x)<x），节点在 ξ 较大（x=L 端）处密集，而非 x=0。

**修正**：改用 `x_i = L·(1 − tanh(β·(1−ξ_i))/tanh(β))`，镜像后在 x=0 处聚簇（β>0）。

**教训**：tanh 拉伸格式有多种约定（Anderson 1984 vs Vinokur 1983），需明确指定聚簇位置再选公式，不能仅凭"直觉"写 tanh(β·ξ)。
