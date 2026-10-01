# FlowLabLite — 功能清单（Feature Inventory）

**版本：** v0.2 · 2026-04-18  
**作者：** Fenfen Yu（余芬芬）  
**用途：** 作为回归基线，防止新功能引入时破坏已验证功能

---

## 1. 求解器

### 1.1 Chorin 投影法（瞬态）

| 项目 | 说明 |
|---|---|
| 算法 | 分步投影法（Fractional-Step），显式时间推进 |
| 压力求解 | Gauss-Seidel（固定 50 次内迭代） |
| 全局状态 | `g_u`, `g_v`, `g_p`（41×41 节点）；`g_steps`（步数计数） |
| 时间步 | nt=500，dt=0.001 |
| 参数 | Re=20（ν=0.1，ρ=1.0，U_lid=1，L=2） |
| 代码位置 | `cmd/main/main.mbt` 行 440–520 |

**WASM 导出（27 个）**：
`init_simulation`, `run_all_steps`, `run_n_steps`,
`get_nx`, `get_ny`, `get_nt`, `get_nit`, `get_re`, `get_dx`, `get_dy`, `get_dt`, `get_rho`, `get_nu`,
`get_step_count`, `get_u_at`, `get_v_at`, `get_p_at`,
`get_velocity_magnitude_at`, `get_max_velocity_magnitude`,
`get_u_center`, `get_v_center`, `get_p_center`,
`get_divergence_norm`, `get_max_u`, `get_max_v`, `get_max_p`, `get_min_p`

### 1.2 SIMPLE 压力修正法（稳态迭代）

| 项目 | 说明 |
|---|---|
| 算法 | Semi-Implicit Method for Pressure-Linked Equations |
| 类型 | **稳态迭代**（无时间步，只有迭代步；无直接法） |
| 全局状态 | `g_u_s`, `g_v_s`, `g_p_s`（41×41）；`g_simple_iters`, `g_simple_residual` |
| 欠松弛 | α_p=0.3（压力），α_u=0.7（速度） |
| 残差 | 简化版：div_norm（已实现）；严格版：方程残差（Step M6 新增） |
| 代码位置 | `cmd/main/main.mbt` 行 822–913 |

**WASM 导出（11 个）**：
`init_simple`, `run_simple_n_iter`,
`get_simple_step_count`, `get_simple_residual`,
`get_u_simple_at`, `get_v_simple_at`, `get_p_simple_at`,
`get_max_u_simple`, `get_simple_divergence_norm`

### 1.3 Chorin-PCG 同位网格求解器（瞬态）

| 项目 | 说明 |
|---|---|
| 算法 | Chorin 投影法 + 共轭梯度压力求解 |
| PCG 预条件器 | Jacobi / DILU / DIC / GAMG（4 种，运行时选择） |
| PCG 配置 | 容差 1e-5，最大迭代 200 |
| 全局状态 | `g_u_pcg`, `g_v_pcg`, `g_p_pcg`（41×41）；`g_steps_pcg`, `g_pcg_last_iters` |
| 代码位置 | `cmd/main/main.mbt` 行 1473–1781 |

**WASM 导出（13 个）**：
`init_chorin_pcg`, `run_chorin_pcg_n_steps`,
`get_pcg_step_count`, `get_pcg_last_iters`,
`get_u_pcg_at`, `get_v_pcg_at`, `get_p_pcg_at`,
`get_velocity_magnitude_pcg_at`, `get_max_u_pcg`, `get_max_v_pcg`,
`get_max_p_pcg`, `get_min_p_pcg`, `get_pcg_divergence_norm`

### 1.4 MAC 交错网格求解器（瞬态）

| 项目 | 说明 |
|---|---|
| 算法 | Marker-and-Cell（Harlow-Welch），交错网格 + PCG 压力 |
| 网格 | 压力 40×40 单元中心；u_mac 40×41 x 面；v_mac 41×40 y 面 |
| 特性 | 交错格式保证散度机器精度守恒 |
| 全局状态 | `g_u_mac`(40×41), `g_v_mac`(41×40), `g_p_mac`(40×40)；mac_nc=40 |
| 代码位置 | `cmd/main/main.mbt` 行 2080–2396 |

**WASM 导出（14 个）**：
`init_mac`, `run_mac_n_steps`,
`get_mac_step_count`, `get_mac_last_iters`, `get_mac_nc`,
`get_u_mac_at`, `get_v_mac_at`, `get_p_mac_at`,
`get_velocity_magnitude_mac_at`,
`get_max_u_mac`, `get_max_v_mac`, `get_max_p_mac`, `get_min_p_mac`,
`get_mac_divergence_norm`

### 1.5 PCG 压力求解预条件器（供 PCG/MAC 使用）

| 预条件器 | 类型 | 代码行 |
|---|---|---|
| Jacobi | 对角缩放 | 行 982 |
| DILU | 修正不完全 LU（非对称） | 行 1126 |
| DIC | IC(0) 不完全 Cholesky（SPD） | 行 1155 |
| GAMG | 两级几何多重网格 V-cycle（Richardson 迭代） | 行 1408 |

---

## 2. I/O 格式与 CLI

### 2.1 输出格式（已实现）

| 格式 | CLI 选项 | 输出文件 | 函数 | 代码位置 |
|---|---|---|---|---|
| JSON（默认） | （无参数） | results.json | `output_json()` | main.mbt 行 640 |
| CSV（RFC 4180） | `--format csv` | cavity_all.csv | `format_csv()` | io_formats.mbt |
| VTK Legacy ASCII | `--format vtk` | cavity_{solver}.vtk | `format_vtk()` | io_formats.mbt |
| Tecplot ASCII | `--format tecplot` | cavity_all.dat | `format_tecplot()` | io_formats.mbt |
| NetCDF-4（Stub） | `--format netcdf` | — | `format_netcdf()` | io_netcdf_stub.mbt |

**VTK 格式规范**：STRUCTURED_GRID v2.0，DIMENSIONS nx ny 1，POINT_DATA 包含 SCALARS pressure、VECTORS velocity、SCALARS vel_magnitude。

**Tecplot 格式规范**：多 ZONE（每求解器一个），DATAPACKING=POINT，VARIABLES X Y U V P VelMag，SOLUTIONTIME=step×dt。

**CSV 列**：`solver,i,j,x,y,u,v,p,vel_mag`（RFC 4180 头行，每求解器 nc² 行）。

**NetCDF Stub**：永远编译通过，返回 `Err`；`--list-formats` 中显示 ✗。

### 2.2 CLI 参数（已实现）

| 参数 | 值 | 说明 |
|---|---|---|
| `--format` | `json\|csv\|vtk\|tecplot\|netcdf` | 输出格式 |
| `--solver` | `all\|chorin\|simple\|pcg\|mac` | 筛选求解器 |
| `--list-formats` | — | 打印格式列表（含 ✓/✗ 状态） |

**安静模式（Quiet Mode）**：`--format` 触发时，`g_quiet_mode[0]=true` 抑制三处进度 println（Chorin/PCG/MAC 的 "Time step N/500 completed"），保证 stdout 为纯格式数据可直接重定向。

### 2.3 数据收集（io_cli.mbt）

| 函数 | 说明 |
|---|---|
| `collect_mac_export()` | 从 MAC 交错网格状态插值到 40×40 单元中心，构造 SolverExport |
| `collect_solver_exports(filter)` | 按过滤条件收集一个或多个求解器数据 |
| `dispatch_format_cli(args)` | 根据 args 分发到对应格式函数并输出 |

### 2.4 SolverExport 结构体

```moonbit
struct SolverExport {
  name   : String               // "chorin" | "simple" | "pcg" | "mac"
  label  : String               // 人类可读描述
  nc     : Int                  // 网格点数（41 或 40）
  h      : Double               // 网格间距 dx
  x0, y0 : Double               // 第一个网格点坐标
  steps  : Int                  // 完成的时间步数
  u, v, p : Array[Array[Double]] // 场数据（y×x 索引）
}
```

---

## 3. 浏览器可视化

### 3.1 main.html（实时求解与可视化）

**功能**：

| 功能 | 说明 |
|---|---|
| 实时 WASM 求解 | 四套求解器在浏览器中运行，无服务器依赖 |
| 速度幅值热力图 | jet 色映射（蓝→绿→红→黄），20×20 Canvas |
| 压力热力图 | 蓝红发散色映射 |
| 流线叠加 | RK4 积分，自动调整密度和颜色 |
| 求解器选择 | 四个标签页：Chorin / SIMPLE / PCG / MAC |
| 控制按钮 | Run N Steps / Pause / Reset / Step |
| 统计面板 | Re, dx, dt, 步数, div_norm, max_u/v, PCG iters 等 |
| u/v 分量场 | 独立热力图（Step 2 已实现） |
| 颜色映射控制 | min/max 范围调节（Step 2 已实现） |
| 下载按钮 | ⬇ VTK / ⬇ Tecplot / ⬇ CSV（通过 WASM 数据，JS 生成文件） |

**技术约束**：Chrome 115+（wasm-gc stringref 0x77）；必须 HTTP 服务（不能 file://）。

### 3.2 local_viewer.html（离线 JSON 结果查看）

**功能**：

| 功能 | 说明 |
|---|---|
| 拖放加载 | 接受 .json 文件（results.json） |
| 四求解器标签 | Chorin / SIMPLE / PCG / MAC |
| 热力图 | 速度幅值、压力，与 main.html 相同的渲染 |
| 颜色范围控制 | 手动 min/max 或自动缩放 |
| 背景色选择 | dark / black / white / gray / navy |
| 统计面板 | Re、步数、网格参数等 |
| 下载按钮 | VTK / CSV / Tecplot（JS 端格式生成） |

---

## 4. 构建与运行脚本

### 4.1 build_wasm.sh

| 功能 | 说明 |
|---|---|
| debug 构建 | `bash build_wasm.sh`（含源码映射） |
| release 构建 | `bash build_wasm.sh release`（优化） |
| 绕过 moon 链接限制 | 显式传 `-exported_functions`（moon 0.1.20260309 不传） |
| 导出验证 | Node.js 脚本解析 WASM export 段并打印函数名 |
| 输出路径 | `./_build/wasm-gc/{debug,release}/build/cmd/main/main.wasm` |

### 4.2 run_local.sh

| 用法 | 说明 |
|---|---|
| `bash run_local.sh` | JSON → results.json（默认） |
| `bash run_local.sh output.json` | 自定义 JSON 文件名 |
| `bash run_local.sh --format vtk` | VTK → cavity_chorin.vtk |
| `bash run_local.sh --format vtk --solver pcg` | VTK → cavity_pcg.vtk |
| `bash run_local.sh --format tecplot` | Tecplot → cavity_all.dat |
| `bash run_local.sh --format csv` | CSV → cavity_all.csv |
| `bash run_local.sh --list-formats` | 显示支持格式 |
| JSON 验证 | Node.js 验证 JSON 语法，打印 Chorin/SIMPLE 点数 |

---

## 5. 测试套件

**运行命令**：`moon test --target wasm`

**当前状态**：**165 tests, passed: 165, failed: 0**

| 文件 | 范围 | 数量 | 说明 |
|---|---|---|---|
| `cmd/main/main_wbtest.mbt` | T1–T58 | 58 | Chorin, SIMPLE, PCG, MAC, preconditioners |
| `cmd/main/main_ext_wbtest.mbt` | T59–T152 | 94 | Extended: I/O, monitors, QUICK, SOR/SSOR, RK3, FVM, scalar, Boussinesq, channel, grids |
| `cmd/main3d/main3d_wbtest.mbt` | 13 basic | 13 | 3D: Chorin GS, PCG, SIMPLE |
| **合计** | T1–T152 + 3D | **165** | — |

### 测试分组详情

| 分组 | 范围 | 内容 |
|---|---|---|
| Chorin 求解器 | T1–T16 | 初始化、边界条件、散度、涡旋、状态隔离 |
| SIMPLE 求解器 | T17–T26 | 迭代收敛、压力修正、独立状态 |
| Chorin-PCG 求解器 | T27–T34 | PCG 收敛、压力精度、状态隔离 |
| MAC 交错网格 | T35–T42 | 散度消去、涡旋形成、精度验证 |
| DILU 预条件器 | T43–T48 | 修正对角线、BC 保持、Jacobi 匹配 |
| DIC 预条件器 | T49–T54 | IC(0) Cholesky、BC 保持、Jacobi 匹配 |
| GAMG 预条件器 | T55–T58 | 两级 V-cycle、BC 保持、Jacobi 匹配 |
| 数组工具 | T59–T60 | 深拷贝隔离、行独立性 |
| Laplacian 算子 | T61–T63 | x²/y² 精确性、调和函数零值 |
| 边界条件 | T64–T65 | 幂等性、粗网格 BC |
| RHS 源项 | T66–T67 | 无散零值、公式验证 |
| GAMG 子组件 | T68–T70 | 插值保常数、限制归一化、光滑器下降 |
| PCG 残差 | T71–T74 | 四种预条件器均满足停止准则 |
| 跨预条件器一致性 | T75 | 解的相对差 < 0.1% |
| 系统物理 | T76–T79 | SND 符号约定、确定性、步数可加性、SIMPLE 质量守恒 |
| 回归 | T80–T83 | 修正对角范围、MAC 散度衰减、Laplacian 精度、全预条件器涡旋符号 |
| I/O 格式 | T84–T91 | CSV 头行、行数、VTK 头部/点数、Tecplot ZONE/SOLUTIONTIME、NetCDF stub |
| 高级 I/O & 监控 | T92–T108 | VTI 格式、MonitorRecord、Ghia 基准、QUICK 对流、SOR/SSOR、HTML 报告、收敛监控 |
| TVD 对流格式 | T109–T113 | Van Leer、Superbee、均匀/线性场精确性 |
| RK3 3D 求解器 | T114–T117 | 初始化、步数计数、涡旋检测 |
| FVM 求解器 | T118–T122 | 初始化、BC、散度、涡旋、状态隔离 |
| 标量运输 | T123–T127 | 初始化、BC、扩散、步数、隔离 |
| Boussinesq 热浮力 | T128–T132 | 初始化、速度/标量 BC、浮力环流、隔离 |
| 通道流周期 BC | T133–T137 | Neumann/周期 LR、初始化、Poiseuille 解、隔离 |
| 非均匀网格 | T138–T142 | tanh 坐标、单调性、均匀极限、拉伸、Laplacian 精度 |
| 曲线网格 | T143–T147 | 恒等映射、Jacobian、边界保持、非正交性、位移 |
| 运行时网格调整 | T148–T152 | 网格尺寸更新、小网格运行、散度有限、重置、钳制 |

---

## 6. 文档结构

| 文件 | 内容 |
|---|---|
| `README.md` | 项目概述、快速开始 |
| `CLAUDE.md` | Claude Code 项目简报，关键约束与决策 |
| `io.md` | I/O 模块规划（Tier-1/2 格式、CLI 接口、依赖分析） |
| `solver_monitor.md` | 求解器监控参数规范（本文档配套） |
| `features.md` | 已实现功能清单（本文档，回归基线） |
| `todo.md` | 第一阶段待办（已完成） |
| `todo2.md` | 第二阶段实施计划（当前进行中） |
| `docs/api_reference.md` | WASM API 参考（62 函数签名、约定、示例） |
| `docs/arch.md` | 三层架构、模块职责、设计决策 |
| `docs/preconditioner_theory.md` | DILU/DIC/GAMG 数学基础 |
| `docs/test_report_20260417.md` | 165 测试套件验证报告（T1–T83 原始 + T84–T152 扩展） |
| `docs/test_report_io_20260418.md` | I/O 模块集成验证（含 T84–T91 + 完整 165 test 覆盖） |
| `docs/SKILL.md` | 重要 Bug 解决过程、架构思想记录 |

---

## 7. 已知约束与限制

| 约束 | 影响 | 解决方向 |
|---|---|---|
| `nu = 0.1` 编译期常量 | 无法 CLI 指定 Re | Step C：引入 `g_nu_runtime` |
| `nx = ny = 41` 编译期常量 | 不能改变网格尺寸 | Step K（后续，极高难度） |
| 无文件 I/O | 读入外部网格数据困难 | Step P（FFI + wasi-fs） |
| moon 0.1.20260309 不传导出函数 | 必须用 build_wasm.sh | 等待 moon 修复 |
| wasm-gc stringref | Chrome 115+，Node.js 不可直接加载 | 无替代方案 |
| NetCDF-4 Stub | 无法实际输出 .nc 文件 | 需链接 libnetcdf ≥ 4.4 |

---

## 8. WASM 导出函数完整列表（121 个）

**Chorin（27）**：
`init_simulation`, `run_all_steps`, `run_n_steps`,
`get_nx`, `get_ny`, `get_nt`, `get_nit`, `get_re`, `get_dx`, `get_dy`, `get_dt`, `get_rho`, `get_nu`,
`get_step_count`, `get_u_at`, `get_v_at`, `get_p_at`,
`get_velocity_magnitude_at`, `get_max_velocity_magnitude`,
`get_u_center`, `get_v_center`, `get_p_center`,
`get_divergence_norm`, `get_max_u`, `get_max_v`, `get_max_p`, `get_min_p`

**SIMPLE（11）**：
`init_simple`, `run_simple_n_iter`,
`get_simple_step_count`, `get_simple_residual`,
`get_u_simple_at`, `get_v_simple_at`, `get_p_simple_at`,
`get_max_u_simple`, `get_simple_divergence_norm`

**Chorin-PCG 共位网格（13）**：
`init_chorin_pcg`, `run_chorin_pcg_n_steps`,
`get_pcg_step_count`, `get_pcg_last_iters`,
`get_u_pcg_at`, `get_v_pcg_at`, `get_p_pcg_at`,
`get_velocity_magnitude_pcg_at`, `get_max_u_pcg`, `get_max_v_pcg`,
`get_max_p_pcg`, `get_min_p_pcg`, `get_pcg_divergence_norm`

**MAC 交错网格（14）**：
`init_mac`, `run_mac_n_steps`,
`get_mac_step_count`, `get_mac_last_iters`, `get_mac_nc`,
`get_u_mac_at`, `get_v_mac_at`, `get_p_mac_at`,
`get_velocity_magnitude_mac_at`,
`get_max_u_mac`, `get_max_v_mac`, `get_max_p_mac`, `get_min_p_mac`,
`get_mac_divergence_norm`

**RK3 3D（7）**：
`init_rk3`, `run_rk3_n_steps`, `get_rk3_step_count`,
`get_u_rk3_at`, `get_v_rk3_at`, `get_p_rk3_at`, `get_rk3_divergence_norm`

**FVM（8）**：
`init_fvm`, `run_fvm_n_iter`, `get_fvm_step_count`,
`get_u_fvm_at`, `get_v_fvm_at`, `get_p_fvm_at`, `get_fvm_divergence_norm`

**标量运输（7）**：
`init_scalar`, `set_scalar_alpha`, `get_scalar_alpha`,
`run_scalar_n_steps`, `get_scalar_step_count`, `get_scalar_at`, `get_scalar_mean`

**Boussinesq 热浮力（14）**：
`init_boussinesq`, `set_bous_nu`, `set_bous_alpha`, `set_bous_beta`,
`get_bous_nu`, `get_bous_alpha`, `get_bous_beta`,
`run_boussinesq_n_steps`, `get_bous_step_count`,
`get_u_bous_at`, `get_v_bous_at`, `get_p_bous_at`, `get_phi_bous_at`, `get_bous_nusselt`

**通道流周期 BC（8）**：
`init_channel`, `set_channel_fx`, `set_channel_nu`, `get_channel_fx`,
`run_channel_n_steps`, `get_channel_step_count`,
`get_u_channel_at`, `get_v_channel_at`, `get_p_channel_at`, `get_channel_analytic_u`

**非均匀 tanh 网格（3）**：
`set_grid_size`, `get_tanh_y_at`, `get_tanh_dy_at`, `get_tanh_mean_dy`

**曲线网格（7）**：
`set_curv_c`, `get_curv_c`, `get_curv_x_at`, `get_curv_y_at`,
`get_curv_jacobian_at`, `get_curv_nonortho_at`, `get_curv_g11_at`, `get_curv_g12_at`, `get_curv_g22_at`

---

*本文档随功能迭代更新，每次 todo 完成后同步修订版本号和测试数量。*
