# FlowLabLite — 数值方法理论公式

> 本文档记录 FlowLabLite 项目已实现及计划实现的所有数值方法的理论基础和离散公式。
> 章节编号对应 `todo3.md` 的实施步骤。
> 变量约定：$i$（行，y 方向），$j$（列，x 方向），$\Delta x = \Delta y = h$（正方形网格）。

---

## §1 数值方法综述（已实现）

### 1.1 已实现求解器列表

| 求解器 | 压力求解 | 时间推进 | 网格类型 |
|---|---|---|---|
| Chorin | Gauss-Seidel（50 迭代/步）| 显式 Euler | 同位 |
| SIMPLE | Gauss-Seidel + 欠松弛 | 稳态迭代 | 同位 |
| Chorin-PCG | PCG（4 种预条件子）| 显式 Euler | 同位 |
| MAC | PCG | 显式 Euler | 交错（Harlow-Welch）|

### 1.2 Chorin 投影法（基准算法）

**步骤 1：** 求中间速度场 $\tilde{\mathbf{u}}$（忽略压力梯度）

$$\tilde{u}_{i,j} = u^n_{i,j} - \Delta t \left(
  u^n \frac{\partial u^n}{\partial x} + v^n \frac{\partial u^n}{\partial y}
\right) + \nu \Delta t \left(
  \frac{\partial^2 u^n}{\partial x^2} + \frac{\partial^2 u^n}{\partial y^2}
\right)$$

一阶迎风离散对流项：
$$u \frac{\partial u}{\partial x} \approx u_{i,j} \frac{u_{i,j} - u_{i,j-1}}{\Delta x}, \quad
  v \frac{\partial u}{\partial y} \approx v_{i,j} \frac{u_{i,j} - u_{i-1,j}}{\Delta y}$$

**步骤 2：** 求解压力 Poisson 方程
$$\frac{\partial^2 p^{n+1}}{\partial x^2} + \frac{\partial^2 p^{n+1}}{\partial y^2}
= \frac{\rho}{\Delta t} \left(
  \frac{\partial \tilde{u}}{\partial x} + \frac{\partial \tilde{v}}{\partial y}
\right)$$

离散 RHS（Gauss-Seidel 内迭代）：
$$b_{i,j} = \frac{\rho}{\Delta t} \left(
  \frac{\tilde{u}_{i,j+1} - \tilde{u}_{i,j-1}}{2\Delta x} +
  \frac{\tilde{v}_{i+1,j} - \tilde{v}_{i-1,j}}{2\Delta y}
\right)$$

**步骤 3：** 速度修正
$$u^{n+1}_{i,j} = \tilde{u}_{i,j} - \frac{\Delta t}{2\rho \Delta x}(p_{i,j+1} - p_{i,j-1})$$

### 1.3 已实现对流格式

**一阶迎风（Upwind）：**
$$(\phi u)_e \approx \phi_{i,j} \cdot u_e \quad (u_e > 0)$$

**三阶 QUICK（Leonard, 1979，Step E）：**
$$\phi_e = \frac{3}{8}\phi_{i,j+1} + \frac{6}{8}\phi_{i,j} - \frac{1}{8}\phi_{i,j-1} \quad (u_e > 0)$$

---

## §2 TVD 格式（Van Leer, Superbee）— todo3 Step I

### 2.1 总变差不增（TVD）条件

**定义：** 数值格式在一维中满足 TVD 若：
$$\text{TV}(u^{n+1}) \leq \text{TV}(u^n), \quad \text{TV}(u) = \sum_i |u_{i+1} - u_i|$$

**Sweby 图：** 限制函数 $\psi(r)$ 需落在 TVD 区域内（$0 \leq \psi(r) \leq \min(2r, 2)$，$r>0$）。

**比例因子：**
$$r_{i+1/2} = \frac{\phi_i - \phi_{i-1}}{\phi_{i+1} - \phi_i} \quad (u > 0)$$

### 2.2 Van Leer 限制函数

$$\psi_{\text{VL}}(r) = \frac{r + |r|}{1 + |r|}$$

特性：
- $\psi(0) = 0$（一阶迎风），$\psi(1) = 1$（二阶中心），$\psi(\infty) \to 2$（最大 TVD 界）
- $C^1$ 光滑，无振荡，适合光滑流动
- 精度：$O(\Delta x^2)$（光滑区），$O(\Delta x)$（极值处）

### 2.3 Superbee 限制函数

$$\psi_{\text{SB}}(r) = \max\left(0, \min(2r, 1), \min(r, 2)\right)$$

特性：
- 更激进：在 Sweby 图中贴近 TVD 区域上界
- 对激波和接触间断有更小的数值扩散
- 光滑区可能产生轻微"过激"（边界处接触 $\psi=2$）

### 2.4 TVD 对流项离散（以 $x$ 方向为例）

$$\frac{\partial(\phi u)}{\partial x} \approx \frac{F_{i+1/2} - F_{i-1/2}}{\Delta x}$$

其中东面通量（$u_e > 0$）：
$$F_{i+1/2} = u_e \left[
  \phi_i + \frac{1}{2} \psi(r_{i+1/2})(\phi_i - \phi_{i-1})
\right]$$

### 2.5 TVD 与 QUICK 精度对比

| 格式 | 光滑区精度 | 激波处 | 色散误差 |
|---|---|---|---|
| 一阶迎风 | $O(h)$ | 无振荡 | 高（数值扩散大）|
| QUICK | $O(h^3)$ | **有振荡**（非单调）| 低 |
| Van Leer TVD | $O(h^2)$ | 无振荡 | 中 |
| Superbee TVD | $O(h^2)$ | 无振荡 | 低（接近 QUICK）|

---

## §3 时间推进格式 — todo3 Step J

### 3.1 三阶 Runge-Kutta（RK3）

**Williamson (1980) 紧凑格式（低存储 SSP-RK3）：**

设右端项 $L(u)$（对流 + 粘性 + 压力梯度之和）：

$$\begin{aligned}
u^{(1)} &= u^n + \alpha_1 \Delta t \, L(u^n) \\
u^{(2)} &= u^{(1)} + \alpha_2 \Delta t \, L(u^{(1)}) \\
u^{n+1} &= u^{(2)} + \alpha_3 \Delta t \, L(u^{(2)})
\end{aligned}$$

SSP-RK3（Shu-Osher，1988）系数：
$$\alpha_1 = 1, \quad \alpha_2 = \frac{1}{4}, \quad \alpha_3 = \frac{2}{3}$$

或等价的 Butcher 表：
$$A = \begin{pmatrix} 0 & 0 & 0 \\ 1 & 0 & 0 \\ 1/4 & 1/4 & 0 \end{pmatrix}, \quad
b = (1/6, 1/6, 2/3)^T, \quad c = (0, 1, 1/2)^T$$

**与 Chorin 投影法的结合：** 每个 RK 子步需要求解一次压力 Poisson 方程，计算量为显式 Euler 的 3×。

### 3.2 隐式 Euler（线性化）

$$\frac{u^{n+1} - u^n}{\Delta t} = L(u^{n+1})$$

线性化（Picard 迭代，以已知 $u^n$ 冻结系数）：
$$\left(\frac{1}{\Delta t} - \nu \nabla^2\right) u^{n+1} = \frac{u^n}{\Delta t} - (u^n \cdot \nabla) u^n$$

优点：无条件稳定（粘性步约束解除），适合高粘性流（低 Re）。
缺点：每时间步需求解线性系统，本项目用 Gauss-Seidel 迭代近似求解。

### 3.3 时间精度与稳定性

| 格式 | 精度 | 稳定性（CFL 限制）| 成本（压力 Poisson）|
|---|---|---|---|
| 显式 Euler | $O(\Delta t)$ | CFL < 1 | 1×/步 |
| RK3（SSP）| $O(\Delta t^3)$ | CFL < $\sqrt{3}$ | 3×/步 |
| 隐式 Euler | $O(\Delta t)$ | 无条件稳定 | 1×/步（但需内迭代）|

**CFL 数定义（当前项目）：**
$$\text{CFL} = \max_{i,j} \left( \frac{|u_{i,j}| \Delta t}{\Delta x} + \frac{|v_{i,j}| \Delta t}{\Delta y} \right)$$

当前参数（Re=20，dt=0.001，dx=0.05）：CFL $\approx 0.02$（远低于 1，安全）。

---

## §4 有限体积法（FVM）基础 — todo3 Step 1

### 4.1 守恒形式的 N-S 方程

动量方程积分守恒形式：
$$\frac{d}{dt} \int_V \rho \mathbf{u} \, dV + \oint_{\partial V} \rho \mathbf{u} (\mathbf{u} \cdot \mathbf{n}) \, dS
= -\oint_{\partial V} p \mathbf{n} \, dS + \oint_{\partial V} \mu \nabla \mathbf{u} \cdot \mathbf{n} \, dS$$

### 4.2 控制体积与面通量离散

以正方形控制体积 $[x_{j-1/2}, x_{j+1/2}] \times [y_{i-1/2}, y_{i+1/2}]$ 为例：

**面积分（中点法则）：**
$$\oint_{\partial V} f \, dS \approx f_e \Delta y - f_w \Delta y + f_n \Delta x - f_s \Delta x$$

**对流通量（东面，$u_e > 0$）：**
$$F^{\text{conv}}_e = \rho u_e \phi_e \Delta y$$

其中 $\phi_e$ 用插值格式（迎风/中心/TVD）给出，$u_e = (u_{i,j} + u_{i,j+1})/2$。

**扩散通量（东面）：**
$$F^{\text{diff}}_e = \mu \frac{\phi_{i,j+1} - \phi_{i,j}}{\Delta x} \Delta y$$

### 4.3 SIMPLE 算法在 FVM 上的实现

**速度预测步（动量方程，以 $u$ 分量为例）：**

$$a_P^u u_P^* = \sum_{\text{nb}} a_{\text{nb}}^u u_{\text{nb}}^n + b_P^u - \frac{(p_e - p_w)\Delta y}{\Delta x}$$

系数（中心差分扩散，迎风对流）：
$$a_E = D_e + \max(-F_e, 0), \quad a_W = D_w + \max(F_w, 0)$$
$$D_e = \mu \Delta y / \Delta x, \quad F_e = \rho u_e \Delta y$$

**压力修正方程（连续性约束）：**
$$a_P^p p'_P = a_E^p p'_E + a_W^p p'_W + a_N^p p'_N + a_S^p p'_S - b^p_P$$

其中 $b^p_P$ 是速度预测场的质量通量不平衡：
$$b^p_P = -(F_e^* - F_w^* + F_n^* - F_s^*)$$

**速度修正步：**
$$u_P^{n+1} = u_P^* - \frac{\Delta y}{a_P^u}(p'_e - p'_w)$$

**欠松弛：**
$$p^{n+1} = p^n + \alpha_p p', \quad u^{n+1} \leftarrow \alpha_u u^{n+1} + (1-\alpha_u) u^n$$

### 4.4 FVM 与 FDM 精度对比

| 方面 | FDM（当前）| FVM（计划）|
|---|---|---|
| 守恒性 | 不严格（微分形式）| 严格（积分形式）|
| 非规则网格 | 困难 | 自然（多面体）|
| 2阶精度 | $O(h^2)$（均匀网格）| $O(h^2)$（均匀网格）|
| 激波处理 | 需特殊格式 | TVD 通量自然适配 |

---

## §5 3D 顶盖驱动流 — todo3 Step 5

### 5.1 三维 N-S 方程（Chorin 投影法）

**步骤 1：** 中间速度场（含 $w$ 分量）

$$\tilde{u}_{i,j,k} = u^n_{i,j,k} - \Delta t \left[
  u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} + w \frac{\partial u}{\partial z}
\right] + \nu \Delta t \nabla^2 u^n_{i,j,k}$$

（$v$, $w$ 分量类似）

**步骤 2：** 三维压力 Poisson
$$\frac{\partial^2 p}{\partial x^2} + \frac{\partial^2 p}{\partial y^2} + \frac{\partial^2 p}{\partial z^2}
= \frac{\rho}{\Delta t} \nabla \cdot \tilde{\mathbf{u}}$$

离散（7 点模板）：
$$b_{i,j,k} = \frac{\rho}{\Delta t}\left(
  \frac{\tilde{u}_{i,j+1,k} - \tilde{u}_{i,j-1,k}}{2\Delta x} +
  \frac{\tilde{v}_{i+1,j,k} - \tilde{v}_{i-1,j,k}}{2\Delta y} +
  \frac{\tilde{w}_{i,j,k+1} - \tilde{w}_{i,j,k-1}}{2\Delta z}
\right)$$

Gauss-Seidel 更新（7 点模板）：
$$p^{(q+1)}_{i,j,k} = \frac{
  (p_{i,j+1,k} + p_{i,j-1,k})\Delta y^2 \Delta z^2 +
  (p_{i+1,j,k} + p_{i-1,j,k})\Delta x^2 \Delta z^2 +
  (p_{i,j,k+1} + p_{i,j,k-1})\Delta x^2 \Delta y^2 -
  b_{i,j,k} \Delta x^2 \Delta y^2 \Delta z^2
}{2(\Delta y^2 \Delta z^2 + \Delta x^2 \Delta z^2 + \Delta x^2 \Delta y^2)}$$

### 5.2 3D 边界条件（6 个面）

| 面 | $u$ | $v$ | $w$ | 压力 |
|---|---|---|---|---|
| 顶盖（$z=L$）| 1 | 0 | 0 | $\partial p/\partial z = 0$ |
| 底面（$z=0$）| 0 | 0 | 0 | $\partial p/\partial z = 0$ |
| 左壁（$x=0$）| 0 | 0 | 0 | $\partial p/\partial x = 0$ |
| 右壁（$x=L$）| 0 | 0 | 0 | $\partial p/\partial x = 0$ |
| 前壁（$y=0$）| 0 | 0 | 0 | $\partial p/\partial y = 0$ |
| 后壁（$y=L$）| 0 | 0 | 0 | $\partial p/\partial y = 0$ |

### 5.3 参考数据

**Albensoeder & Kuhlmann (2002)** 给出 Re=100-1000 三维方腔流的基准数据：
- $u_{\max}$ 沿 $z$-中心线（$y=z=L/2$）
- 一次涡和角涡（Taylor-Görtler 涡）的位置

---

## §6 被动标量输运 — todo3 Step 6

### 6.1 标量对流扩散方程

$$\frac{\partial \phi}{\partial t} + u \frac{\partial \phi}{\partial x} + v \frac{\partial \phi}{\partial y}
= \alpha \left(\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2}\right) + S_\phi$$

其中 $\phi$ 为被动标量（如浓度、温度在 Pe 数较小时），$\alpha = \lambda / (\rho c_p)$ 为扩散系数，$S_\phi$ 为源项。

### 6.2 离散格式

**显式 Euler + 一阶迎风对流：**

$$\phi^{n+1}_{i,j} = \phi^n_{i,j} - \Delta t \left[
  u^n_{i,j} \frac{\phi^n_{i,j} - \phi^n_{i,j-1}}{\Delta x} +
  v^n_{i,j} \frac{\phi^n_{i,j} - \phi^n_{i-1,j}}{\Delta y}
\right] + \alpha \Delta t \left[
  \frac{\phi^n_{i,j+1} - 2\phi^n_{i,j} + \phi^n_{i,j-1}}{\Delta x^2} +
  \frac{\phi^n_{i+1,j} - 2\phi^n_{i,j} + \phi^n_{i-1,j}}{\Delta y^2}
\right]$$

### 6.3 无量纲参数

**Péclet 数：** $Pe = u L / \alpha$

- $Pe \to 0$：纯扩散（椭圆型），稳态解为调和函数
- $Pe \to \infty$：纯对流（双曲型），需要迎风格式

**扩散稳定性：** $\alpha \Delta t / \Delta x^2 \leq 1/2$（1D）；2D 需 $\alpha \Delta t (1/\Delta x^2 + 1/\Delta y^2) \leq 1/2$

### 6.4 验证方案

| Pe 极限 | 解析解 | 验证指标 |
|---|---|---|
| $Pe \to 0$（$u=v=0$，$\alpha$ 大）| $\phi_{ss} = $ 调和函数（$\nabla^2 \phi = 0$）| $\|\phi - \phi_{\text{ref}}\|_2 / \|\phi_{\text{ref}}\|_2 < 10^{-3}$ |
| $Pe \to \infty$（$\alpha=0$）| 纯对流（$\phi$ 沿流线不变）| 质量守恒 $\sum \phi \Delta x \Delta y = \text{const}$ |

---

## §7 能量方程与 Boussinesq 近似 — todo3 Step 7

### 7.1 Boussinesq 近似

在浮力驱动流中，密度仅在重力项中随温度变化（其他处视为常数）：
$$\rho = \rho_0 [1 - \beta (T - T_0)]$$

其中 $\beta$ 为体积热膨胀系数，$T_0$ 为参考温度。

### 7.2 无量纲温度方程

以参考温差 $\Delta T = T_H - T_C$ 无量纲化，令 $\theta = (T - T_C)/\Delta T \in [0, 1]$：

$$\frac{\partial \theta}{\partial t} + u \frac{\partial \theta}{\partial x} + v \frac{\partial \theta}{\partial y}
= \frac{1}{Pr \cdot Re} \left(\frac{\partial^2 \theta}{\partial x^2} + \frac{\partial^2 \theta}{\partial y^2}\right)$$

**或用 Ra/Pr 形式（自然对流）：**

$$\frac{\partial \theta}{\partial t} + \mathbf{u} \cdot \nabla \theta = \frac{1}{\sqrt{Ra \cdot Pr^{-1}}} \nabla^2 \theta$$

### 7.3 Boussinesq 体积力（浮力项）

在 $v$ 动量方程（垂直方向）中添加浮力源项：
$$\frac{\partial v}{\partial t} + \ldots = \ldots + \frac{\rho_0 \beta g \Delta T}{\rho_0} \theta
= \frac{Ra \cdot Pr^{-1}}{Re^2} \theta \quad (\text{无量纲形式})$$

**或直接用 Ra-Pr 参数化（自然对流，无强迫流动）：**
$$f_y = Ra \cdot Pr^{-1} \cdot \theta \quad (\text{仅重力方向})$$

### 7.4 无量纲参数

| 参数 | 定义 | 物理含义 |
|---|---|---|
| $Ra = g \beta \Delta T L^3 / (\nu \alpha)$ | 瑞利数 | 浮力 vs 粘性×扩散 |
| $Pr = \nu / \alpha$ | 普朗特数 | 动量 vs 热量扩散 |
| $Nu = \bar{h} L / \lambda$ | 努塞尔数 | 对流热传递强度 |

### 7.5 验证：De Vahl Davis (1983) 方腔自然对流

**基准工况：** 左壁热（$\theta=1$），右壁冷（$\theta=0$），上下绝热（$\partial\theta/\partial y=0$），$Pr=0.71$

| Ra | $Nu_{\text{avg}}$（参考）| 参考文献 |
|---|---|---|
| $10^3$ | 1.118 | De Vahl Davis (1983) |
| $10^4$ | 2.243 | De Vahl Davis (1983) |
| $10^5$ | 4.519 | De Vahl Davis (1983) |
| $10^6$ | 8.800 | De Vahl Davis (1983) |

验证指标：$|Nu_{\text{sim}} - Nu_{\text{ref}}| / Nu_{\text{ref}} < 1\%$（41×41 网格期望 3-5%）

---

## §8 边界条件库 — todo3 Step 8

### 8.1 Dirichlet 边界（速度/温度给定）

已实现（顶盖驱动流）：
$$\phi_{i,j} = \phi_{\text{BC}} \quad \text{at boundary}$$

### 8.2 Neumann 边界（零梯度，出口）

$$\frac{\partial \phi}{\partial n} = 0 \implies \phi_{\text{ghost}} = \phi_{\text{interior}}$$

离散实现（右壁出口，$j=nc-1$）：
$$\phi_{i, nc-1} = \phi_{i, nc-2}$$

**已实现的 Neumann 压力 BC（所有壁面）：** 已在压力 Poisson 求解中实现。

### 8.3 对称边界（法向分量为零）

在 $x=0$ 对称面：
$$u = 0 \quad \text{（法向）}, \quad \frac{\partial v}{\partial x} = 0 \quad \text{（切向）}$$

$$\implies v_{i, 0} = v_{i, 1}$$

### 8.4 周期性边界（通道流）

$$\phi_{i, 0} = \phi_{i, nc-1}, \quad \phi_{i, nc} = \phi_{i, 1}$$

**应用场景：** 无穷长通道的展向周期性，圆柱绕流（上下周期）。

---

## §9 拉伸网格 — todo3 Step 9

### 9.1 tanh 压缩函数

在 $x$ 方向对节点位置进行非均匀分布（壁面加密）：

$$x_j = \frac{L}{2} \left[1 + \frac{\tanh\left[\beta_s\left(\frac{j}{N-1} - \frac{1}{2}\right)\right]}{\tanh(\beta_s/2)}\right]$$

其中 $\beta_s$ 为拉伸因子（$\beta_s \to 0$：均匀；$\beta_s = 2$：适度加密；$\beta_s = 4$：强加密）。

**局部网格间距：**
$$\Delta x_j = x_{j+1} - x_j, \quad j = 0, \ldots, N-2$$

### 9.2 非等距网格离散系数

**一阶导数（二阶精度，非均匀间距）：**
$$\left.\frac{\partial \phi}{\partial x}\right|_j \approx
\frac{\phi_{j+1} \Delta x_{j-1}^2 - \phi_{j-1} \Delta x_j^2 + \phi_j(\Delta x_j^2 - \Delta x_{j-1}^2)}
{\Delta x_j \Delta x_{j-1}(\Delta x_j + \Delta x_{j-1})}$$

**二阶导数（非均匀间距）：**
$$\left.\frac{\partial^2 \phi}{\partial x^2}\right|_j \approx
\frac{2}{\Delta x_j + \Delta x_{j-1}} \left(
  \frac{\phi_{j+1} - \phi_j}{\Delta x_j} - \frac{\phi_j - \phi_{j-1}}{\Delta x_{j-1}}
\right)$$

---

## §10 曲线网格基础 — todo3 Step 10

### 10.1 坐标变换（$\xi$-$\eta$ 参数化）

物理坐标 $(x, y)$ → 计算坐标 $(\xi, \eta)$（均匀网格）的映射：
$$x = x(\xi, \eta), \quad y = y(\xi, \eta)$$

### 10.2 雅可比矩阵

$$J = \frac{\partial(x, y)}{\partial(\xi, \eta)} = \begin{pmatrix} x_\xi & x_\eta \\ y_\xi & y_\eta \end{pmatrix},
\quad |J| = x_\xi y_\eta - x_\eta y_\xi$$

逆变基矢量（计算空间度量系数）：
$$\xi_x = \frac{y_\eta}{|J|}, \quad \xi_y = -\frac{x_\eta}{|J|}, \quad
  \eta_x = -\frac{y_\xi}{|J|}, \quad \eta_y = \frac{x_\xi}{|J|}$$

### 10.3 曲线坐标系中的扩散算子

$$\nabla^2 \phi = \frac{1}{|J|}\left[
  \frac{\partial}{\partial \xi}\left(\frac{g^{22}}{|J|} \frac{\partial\phi}{\partial\xi}\right)
  - \frac{\partial}{\partial \xi}\left(\frac{g^{12}}{|J|} \frac{\partial\phi}{\partial\eta}\right)
  - \frac{\partial}{\partial \eta}\left(\frac{g^{12}}{|J|} \frac{\partial\phi}{\partial\xi}\right)
  + \frac{\partial}{\partial \eta}\left(\frac{g^{11}}{|J|} \frac{\partial\phi}{\partial\eta}\right)
\right]$$

其中度量张量 $g^{11} = x_\eta^2 + y_\eta^2$，$g^{22} = x_\xi^2 + y_\xi^2$，$g^{12} = -(x_\xi x_\eta + y_\xi y_\eta)$。

### 10.4 FVM 在曲线网格上的实现

控制体积面积：$A_P = |J| \Delta\xi \Delta\eta$

面法向量（东面，$\xi = \text{const}$）：
$$\mathbf{S}_e = (y_\eta, -x_\eta) \Delta\eta$$

面通量：
$$F^{\text{conv}}_e = \rho (\mathbf{u} \cdot \mathbf{S}_e) \phi_e$$

---

## 参考文献

1. **Chorin (1968)** — "Numerical solution of the Navier-Stokes equations," *Math. Comp.* 22, 745-762.
2. **Leonard (1979)** — "A stable and accurate convective modelling procedure based on quadratic upstream interpolation," *Comp. Meth. Appl. Mech. Eng.* 19, 59-98.
3. **Ghia, Ghia & Shin (1982)** — "High-Re solutions for incompressible flow using the Navier-Stokes equations and a multigrid method," *J. Comput. Phys.* 48, 387-411.
4. **Shu & Osher (1988)** — "Efficient implementation of essentially non-oscillatory shock-capturing schemes," *J. Comput. Phys.* 77, 439-471.
5. **De Vahl Davis (1983)** — "Natural convection of air in a square cavity: a bench mark numerical solution," *Int. J. Numer. Methods Fluids* 3, 249-264.
6. **Patankar (1980)** — *Numerical Heat Transfer and Fluid Flow*. Hemisphere Publishing.
7. **Albensoeder & Kuhlmann (2002)** — "Three-dimensional centrifugal-flow instabilities in the lid-driven cavity problem," *Phys. Fluids* 14, 2018-2029.
8. **Sweby (1984)** — "High resolution schemes using flux limiters for hyperbolic conservation laws," *SIAM J. Numer. Anal.* 21, 995-1011.

---

*作者：Fenfen Yu（余芬芬），AI 协作：Claude Sonnet 4.6*
*日期：2026-04-18*
