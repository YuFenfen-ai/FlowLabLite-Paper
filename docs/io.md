# FlowLabLite — I/O 模块规划

版本：v0.1-draft · 日期：2026-04-18 · 状态：待确认

---

## 1. 需求总结

在不影响现有功能的前提下，增加多格式数据导出能力：

| 场景 | 目标 |
|---|---|
| 命令行 | 支持 JSON / CSV / VTK / Tecplot / NetCDF4 输出，默认 JSON |
| 浏览器 (main.html) | 增加 VTK / Tecplot / CSV 下载按钮 |
| 编译可靠性 | 无论有没有安装可选依赖，`moon build` 均能成功；缺少依赖时给出明确提示 |

---

## 2. 格式规范名称（官方/业界标准）

| 格式 | 规范全称 | 推荐扩展名 | 主要读取软件 |
|---|---|---|---|
| JSON | — | `.json` | main.html, local_viewer.html, Node.js |
| CSV | RFC 4180 Comma-Separated Values | `.csv` | Excel, Python/pandas, MATLAB |
| VTK Legacy ASCII | VTK Legacy Format v4.2 (ASCII) | `.vtk` | ParaView, VisIt, VTK Python |
| Tecplot ASCII | Tecplot Data Format Guide (ASCII .dat) | `.dat` | Tecplot 360, ParaView (Tecplot reader) |
| NetCDF-4 | CF Conventions v1.11 / NetCDF-4 Classic | `.nc` | ncview, Panoply, CDO, xarray, MATLAB |

> **注意**：VTK XML 格式（`.vts` `.vtu`）比 Legacy 更现代，但 Legacy ASCII 无依赖、可直接文本生成，优先选用。
> Tecplot 二进制（`.plt`）需官方库，本规划仅实现 ASCII `.dat`。

---

## 3. 依赖分析与实现分层

### Tier 1 — 纯 MoonBit，零外部依赖（永远可编译）

| 格式 | 实现原理 |
|---|---|
| JSON | 已实现，`output_json()` in main.mbt |
| CSV | 纯字符串拼接，RFC 4180 |
| VTK Legacy ASCII | 纯字符串拼接，固定文本协议 |
| Tecplot ASCII | 纯字符串拼接，ZONE 结构化格式 |

### Tier 2 — 可选 C FFI，缺失时自动降级

| 格式 | 外部依赖 | 最低版本 | 安装示例 |
|---|---|---|---|
| NetCDF-4 | `libnetcdf` + `libhdf5` | netcdf ≥ 4.4 | `apt install libnetcdf-dev` / Conda |

> **结论**：CSV / VTK / Tecplot 均为纯文本格式，**不需要任何外部库**，
> 只有 NetCDF-4 是真正的可选依赖。

---

## 4. 降级编译机制

### 4.1 设计方案：永远编译的 stub + 可选 FFI 文件

```
cmd/main/
  io_formats.mbt          # Tier-1：JSON/CSV/VTK/Tecplot 纯 MoonBit 实现
  io_netcdf_stub.mbt      # 默认包含：所有 netcdf 函数返回 Err 并打印提示
  io_netcdf_ffi.mbt       # 可选：C FFI 实现（用户手动替换 stub）
  io_dispatch.mbt         # 格式分发 + --list-formats 查询
  io_cli.mbt              # 命令行参数解析（--format / --output / --solver 等）
```

### 4.2 编译期提示策略

- `io_netcdf_stub.mbt` 顶部注释中列出如何启用真正的 NetCDF4 支持
- 运行时：当用户指定 `--format netcdf` 但当前使用 stub 时，输出：
  ```
  [io] NetCDF-4 输出不可用：当前构建未包含 libnetcdf FFI。
       如需启用，请参阅 io.md §4.2 安装说明。
  ```
- 将来若 MoonBit 支持条件编译标志（`#[cfg(feature = "netcdf")]`），迁移为标准方式

### 4.3 暂不实现 NetCDF-4（阶段 1 范围外）

阶段 1 仅交付 Tier-1 格式（JSON / CSV / VTK / Tecplot）+ CLI 框架 + stub 桩。
NetCDF-4 FFI 实现列入阶段 2，待 `moonbit-c-binding` 插件接口稳定后推进。

---

## 5. 模块设计

### 5.1 数据结构（`io_formats.mbt`）

```moonbit
// 单求解器的完整场数据（41×41 或 40×40）
struct SolverExport {
  name    : String           // "chorin" | "simple" | "pcg" | "mac"
  label   : String           // 显示名称
  nc      : Int              // 网格节点数（nx = ny = nc）
  dx      : Double
  steps   : Int
  u       : Array[Array[Double]]
  v       : Array[Array[Double]]
  p       : Array[Array[Double]]
  // 统计摘要
  max_u   : Double
  max_v   : Double
  div_norm: Double
}

// 导出入口
fn export_csv(s : SolverExport) -> String
fn export_vtk(s : SolverExport) -> String
fn export_tecplot(s : SolverExport) -> String
fn export_all_tecplot(solvers : Array[SolverExport]) -> String  // 多 ZONE

// stub（阶段 1）
fn export_netcdf(s : SolverExport) -> Result[String, String]
  // = Err("NetCDF-4 not available in this build.")
```

### 5.2 VTK Legacy ASCII 格式草稿

```
# vtk DataFile Version 2.0
FlowLabLite chorin Re=20 step=500
ASCII
DATASET STRUCTURED_GRID
DIMENSIONS 41 41 1
POINTS 1681 double
0.0  0.0  0.0
0.05 0.0  0.0
...
POINT_DATA 1681
SCALARS pressure double 1
LOOKUP_TABLE default
<p[0][0]> <p[0][1]> ...
VECTORS velocity double
<u[0][0]> <v[0][0]> 0.0
...
SCALARS vel_magnitude double 1
LOOKUP_TABLE default
<mag[0][0]> ...
```

### 5.3 Tecplot ASCII 格式草稿

```
TITLE = "FlowLabLite chorin Re=20 step=500"
VARIABLES = "X", "Y", "U", "V", "P", "VelMag"
ZONE T="chorin", I=41, J=41, DATAPACKING=POINT, SOLUTIONTIME=0.5
0.0  0.0  0.0 0.0  <p> <mag>
0.05 0.0  <u> <v>  <p> <mag>
...
```

多求解器时追加额外 `ZONE` 块（SOLUTIONTIME 依次递增），一个 `.dat` 文件包含全部求解器结果。

### 5.4 CSV 格式草稿

```
# FlowLabLite export v1.0 | solver=chorin | Re=20 | steps=500
solver,i,j,x,y,u,v,p,vel_mag
chorin,0,0,0.000000,0.000000,0.000000,0.000000,...
...
simple,0,0,...
```

所有求解器数据追加在同一文件（solver 字段区分），或 `--output-dir` 模式下每求解器一个文件。

---

## 6. CLI 接口设计

### 6.1 用法

```bash
# 默认：等同现有行为，输出 JSON
moon run cmd/main --target wasm

# 指定格式
moon run cmd/main --target wasm -- --format json
moon run cmd/main --target wasm -- --format csv
moon run cmd/main --target wasm -- --format vtk
moon run cmd/main --target wasm -- --format tecplot
moon run cmd/main --target wasm -- --format netcdf   # stub: 打印提示并退出

# 自定义输出文件（不指定则打印到 stdout）
moon run cmd/main --target wasm -- --format vtk --output cavity.vtk

# 多格式同时输出（写到 ./output/ 目录）
moon run cmd/main --target wasm -- --format all --output-dir ./output/

# 查询当前构建支持哪些格式
moon run cmd/main --target wasm -- --list-formats

# 仅运行指定求解器（不影响 WASM 导出）
moon run cmd/main --target wasm -- --solver chorin --format vtk
```

### 6.2 --list-formats 输出示例

```
Available formats in this build:
  json      RFC-compatible JSON (default)          ✓
  csv       RFC 4180 CSV                           ✓
  vtk       VTK Legacy ASCII (.vtk)                ✓
  tecplot   Tecplot ASCII Data Format (.dat)       ✓
  netcdf    CF-1.11 NetCDF-4 (.nc)                 ✗  (libnetcdf not linked)
```

### 6.3 参数解析实现

MoonBit 目前访问 WASI argv 的方式（`moon run --target wasm` 走 wasmtime）：

```moonbit
// 使用 moonbitlang/core 的 @sys.argv() 或 @env.args()
// 实现时确认可用接口；备选方案：从环境变量 FLOWLAB_ARGS 读取
```

> 实现阶段需先验证 `moonbitlang/core` 提供哪个 argv API（wasm-wasi vs native 有差异）。

---

## 7. 浏览器端设计（main.html）

浏览器端**不**引入额外 WASM 导出——格式生成在 JavaScript 侧完成
（数据已通过现有 `get_*_at()` WASM 函数访问，JS 直接格式化字符串）。

原因：WASM-GC 跨边界传字符串复杂；格式生成纯文本逻辑可在 JS 中等价实现；MoonBit 格式函数只用于 CLI 路径。

### 新增下载按钮（工具栏）

```
[Run N Steps] [Pause] [Step] [Reset]  |  [⬇ VTK] [⬇ Tecplot] [⬇ CSV]  |  [Streamlines] [Show u/v]
```

- 点击后生成文本 Blob，触发浏览器下载
- 文件名格式：`flowlab_<solver>_step<N>.<ext>`
- NetCDF 按钮**不**添加（二进制格式，浏览器 JS 不易生成；用户需要 NetCDF 请走命令行）

---

## 8. 文件变更清单

| 文件 | 操作 | 说明 |
|---|---|---|
| `cmd/main/io_formats.mbt` | 新建 | CSV / VTK / Tecplot 生成函数 |
| `cmd/main/io_netcdf_stub.mbt` | 新建 | NetCDF stub，返回 Err + 提示 |
| `cmd/main/io_dispatch.mbt` | 新建 | 格式枚举、可用性检查、分发逻辑 |
| `cmd/main/io_cli.mbt` | 新建 | argv 解析，`--format` / `--output` 等 |
| `cmd/main/main.mbt` | 修改 | `fn main` 分支：默认 JSON，有参数走 io_dispatch |
| `cmd/main/main.html` | 修改 | 工具栏加 VTK / Tecplot / CSV 下载按钮 |
| `run_local.sh` | 修改（可选） | 增加 `--format` 透传支持 |
| `cmd/main/moon.pkg.json` | 不变 | 无新 WASM 导出（格式函数仅 CLI 使用） |

---

## 9. 测试计划

| 测试编号 | 测试内容 | 通过标准 |
|---|---|---|
| T84 | CSV 输出包含所有 4 个求解器的行 | 行数 = 4 × 41² |
| T85 | VTK 输出 DIMENSIONS 字段为 `41 41 1` | 字符串包含检查 |
| T86 | VTK POINTS 数量 = 41×41 = 1681 | 计数验证 |
| T87 | Tecplot 输出含 4 个 ZONE 块 | 字符串包含 4 × `ZONE T=` |
| T88 | `--format netcdf` 不崩溃，返回错误提示 | 返回 `Err(...)` |
| T89 | `--list-formats` 输出包含 5 行 | stdout 检查 |
| T90 | VTK 压力值与 `get_p_at(0,0)` 一致 | 数值比对 |
| T91 | Tecplot 第一行 X/Y 坐标与期望网格坐标一致 | 数值比对 |

月测试套件目标：83 → **91 tests, all passing**

---

## 10. 实施步骤

经确认后，按以下顺序推进：

### Step A — 核心格式函数（`io_formats.mbt`）
- 实现 `export_csv()` / `export_vtk()` / `export_tecplot()`
- 单元测试 T84–T88, T90–T91
- `moon test --target wasm` 全绿后提交

### Step B — CLI 参数解析（`io_cli.mbt` + `io_dispatch.mbt`）
- 实现 argv 解析，`--format` / `--output` / `--list-formats`
- 修改 `fn main` 入口
- 集成测试：`moon run cmd/main -- --format vtk --output /tmp/test.vtk` 验证文件可被 ParaView 读取
- 通过后提交

### Step C — 浏览器下载按钮（`main.html`）
- 工具栏添加 VTK / Tecplot / CSV 下载按钮（JS 格式化，不改 WASM）
- 手动验证：下载 `.vtk` 文件，ParaView 正确显示速度场
- 通过后提交

### Step D — `run_local.sh` 扩展（可选，较小改动）
- 支持 `bash run_local.sh --format vtk` 透传参数
- 提交

---

## 11. 未决问题（需实现阶段确认）

1. **MoonBit argv API**：`moonbitlang/core` 在 wasm-wasi target 下的命令行参数访问接口（`@sys.args()`？`@env.args()`？）——需查阅最新文档或测试
2. **NetCDF-4 FFI 时间线**：`moonbit-c-binding` 插件当前能力是否足以封装 `libnetcdf`？——列入阶段 2 评估
3. **run_local.sh 格式透传**：JSON 以外的格式无需 `===JSON_DATA_START===` 标记提取，需调整 shell 脚本逻辑
4. **MAC 网格坐标系**：MAC 使用 cell-center 坐标（40×40），VTK/Tecplot ZONE 应标注为 `cell-centered`；与其他 41×41 node-centered 求解器分开输出或统一插值？——建议：各求解器独立 ZONE/文件，保留原始坐标

---

*确认后进入 Step A 实施。*
