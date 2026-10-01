# PCG 预条件器理论文档

**项目**: FlowLabLite  
**作者**: Fenfen Yu（余芬芬）  
**日期**: 2026-04-16  
**状态**: 设计完成，待实现  

---

## 1  背景与动机

FlowLabLite 的压力 Poisson 方程求解目前使用 **Jacobi（对角）预条件共轭梯度法（PCG）**。
Jacobi 预条件器简单、稳健，但对高雷诺数、细网格问题的收敛加速有限。

本文档设计并推导以下四种预条件器：

| 预条件器 | 类型 | 实现复杂度 | 收敛加速 | 适用场景 |
|---|---|---|---|---|
| **Diagonal (Jacobi)** | 对角缩放 | ★☆☆☆ | 基准 | 均匀网格，低 Re |
| **DILU** | 修正对角 ILU | ★★☆☆ | ~2× Jacobi | 均匀/非均匀网格 |
| **DIC** | 不完全 Cholesky | ★★★☆ | ~3–4× Jacobi | SPD 系统 |
| **GAMG** | 几何多重网格 | ★★★★ | ~5–10× Jacobi | 大规模、多层次 |

---

## 2  压力 Poisson 方程的矩阵形式

### 2.1  离散化

在 41×41 节点均匀网格（dx = dy = 0.05）上，5 点有限差分格式：

$$\frac{p_{i,j+1} - 2p_{i,j} + p_{i,j-1}}{dx^2} + \frac{p_{i+1,j} - 2p_{i,j} + p_{i-1,j}}{dy^2} = b_{i,j}$$

写成矩阵形式 **Ap = b**（对内部节点排序，共 N = 39×39 = 1521 个未知量）：

$$A_{ij,ij} = -\left(\frac{2}{dx^2} + \frac{2}{dy^2}\right) = -a_{\text{diag}}$$

$$A_{ij,ij\pm1} = \frac{1}{dx^2} = \text{inv\_dx2} > 0 \quad \text{（东西相邻）}$$

$$A_{ij\pm1,ij} = \frac{1}{dy^2} = \text{inv\_dy2} > 0 \quad \text{（南北相邻）}$$

矩阵 A 是**对称负定（SND）**，−A 是**对称正定（SPD）**。

### 2.2  边界条件影响

| 边界 | 类型 | 矩阵影响 |
|---|---|---|
| 左/右壁（x=0, x=2） | Neumann dp/dx=0 | 幽灵节点反射，有效对角不变 |
| 底壁（y=0） | Neumann dp/dy=0 | 同上 |
| 顶盖（y=2） | Dirichlet p=0 | 移入 RHS，对角保持 a_diag |

### 2.3  代码中的符号约定

现有代码的 `laplacian_apply` 计算的是正 Laplacian（= −(−∇²)p = ∇²p），
其矩阵的对角元为 **−a_diag < 0**，非对角元为 **+inv_dx2 > 0**。

PCG 代码以 **a_diag（正数）** 作为 Jacobi 预条件器参数，等价于用 **−diag(A)** 近似 **−A**（SPD）。
所有新预条件器遵循同一符号约定：以 a_diag 为"正对角"，以 inv_dx2 为"正偏对角"。

---

## 3  Diagonal（Jacobi）预条件器

### 3.1  原理

最简单的预条件器，用矩阵对角 M = diag(−A) 近似 −A：

$$M_{\text{Jac}} = a_{\text{diag}} \cdot I, \quad M^{-1}_{\text{Jac}} r = \frac{r}{a_{\text{diag}}}$$

对均匀网格，a_diag 是常数，因此每次 apply 仅需一次标量除法。

### 3.2  公式

$$z_{i,j} = \frac{r_{i,j}}{a_{\text{diag}}} = \frac{r_{i,j}}{\frac{2}{dx^2} + \frac{2}{dy^2}}$$

### 3.3  收敛分析

条件数：$\kappa(M^{-1}A) \approx O(1/h^2)$，PCG 迭代次数 $O(1/h)$。  
对 41×41，理论迭代数 ~40 步；实测 20–80 步（取决于时间步 RHS 变化幅度）。

### 3.4  代码状态

**已实现**（`pressure_poisson_pcg`，见 `cmd/main/main.mbt`）。不需修改。

---

## 4  DILU（对角不完全 LU）预条件器

### 4.1  原理

DILU（Diagonal-based Incomplete LU）是对 Jacobi 的改进：  
在 LU 分解中只保留**对角元素的修正**，不做三角求解。等价于用**修正对角**替代原对角。

对称矩阵情况下，DILU 与 DIC 的**分解步骤完全相同**，但 apply 步骤只用对角（无前向/后向扫描）。

### 4.2  建立（Setup）阶段

按**行主序**（row-major）从 (1,1) 到 (ny−2, nx−2) 更新修正对角 d：

$$d[i][j] = a_{\text{diag}} - \frac{(\text{inv\_dx2})^2}{d[i][j-1]} - \frac{(\text{inv\_dy2})^2}{d[i-1][j]}$$

初始条件（无西邻/南邻时对应项为 0）：

$$d[1][1] = a_{\text{diag}}, \quad d[1][j] = a_{\text{diag}} - \frac{(\text{inv\_dx2})^2}{d[1][j-1]}, \quad \text{等}$$

**稳定性保证**：对 SPD 矩阵，d[i][j] > 0（无需担忧零除）。若极端情况 d 接近零，回退到 a_diag。

### 4.3  应用（Apply）阶段

$$z[i][j] = r[i][j] \cdot d^{-1}[i][j]$$

只做逐元素缩放（O(N) 操作），与 Jacobi 复杂度相同但收敛效果更好。

### 4.4  与 Jacobi 的差异

对 dx=dy=0.05，第一个非角点 (1,2)：

$$d[1][2] = 1600 - \frac{400^2}{1600} = 1600 - 100 = 1500$$

修正幅度约 6%。对边角节点（累积效应）修正可达 15–20%，有效降低条件数。

### 4.5  OpenFOAM 参考

**文件**: `src/OpenFOAM/matrices/lduMatrix/preconditioners/DILUPreconditioner/DILUPreconditioner.C`  
**仓库**: https://github.com/OpenFOAM/OpenFOAM-dev  
**许可**: GPL-3.0  

关键代码段（修正对角，对非对称矩阵 upper ≠ lower）：

```cpp
// src/OpenFOAM/matrices/lduMatrix/preconditioners/
// DILUPreconditioner/DILUPreconditioner.C  (OpenFOAM-dev, GPL-3.0)
void Foam::DILUPreconditioner::calcReciprocalD
(
    scalarField& rD,
    const lduMatrix& matrix
)
{
    scalar* __restrict__ rDPtr = rD.begin();
    const label* const __restrict__ uPtr = matrix.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix.upper().begin();
    const scalar* const __restrict__ lowerPtr = matrix.lower().begin();

    const label nFaces = matrix.upper().size();
    for (label face=0; face<nFaces; face++)
    {
        rDPtr[uPtr[face]] -= upperPtr[face]*lowerPtr[face]/rDPtr[lPtr[face]];
    }
    // Invert diagonal
    forAll(rD, cell)
    {
        rDPtr[cell] = 1.0/rDPtr[cell];
    }
}
```

对称矩阵时 upper = lower，因此 `upperPtr[face]*lowerPtr[face] = upper²`，
对应我们的公式 `d[i][j] -= inv_dx2² / d[west]`。

### 4.6  参考文献

- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, 2nd ed. SIAM, §10.3.
- OpenFOAM Foundation. *DILUPreconditioner* source code, v2406.

---

## 5  DIC（对角不完全 Cholesky）预条件器

### 5.1  原理

DIC 是 SPD 矩阵的不完全 Cholesky 分解 IC(0)，以 **L·D·L^T** 形式（D 为对角矩阵，L 为单位下三角）近似 −A：

$$-A \approx L \cdot D \cdot L^T$$

与 DILU 相比，DIC 在 apply 阶段额外执行**前向扫描**和**后向扫描**，
将偏对角的耦合效应传播到整个预条件器向量，显著降低条件数。

### 5.2  建立（Setup）阶段

修正对角公式与 DILU 相同：

$$d[i][j] = a_{\text{diag}} - \frac{(\text{inv\_dx2})^2}{d[i][j-1]} - \frac{(\text{inv\_dy2})^2}{d[i-1][j]}$$

L 的单位下三角元素（行主序中 west 和 south 邻居先于 (i,j)）：

$$L_{ij,\,ij-1} = \frac{-\text{inv\_dx2}}{d[i][j-1]}, \quad L_{ij,\,i-1j} = \frac{-\text{inv\_dy2}}{d[i-1][j]}$$

（符号为负，因为 −A 的偏对角为 −inv_dx2）

### 5.3  应用（Apply）阶段

给定残差 r，计算 z = M^{-1} r（求解 L·D·L^T·z = r）：

**步骤 1：前向扫描（L·w = r，按行主序 (1,1)→(ny−2,nx−2)）**

$$w[i][j] = r[i][j] + \frac{\text{inv\_dx2}}{d[i][j-1]} \cdot w[i][j-1] + \frac{\text{inv\_dy2}}{d[i-1][j]} \cdot w[i-1][j]$$

（边界节点邻居贡献为 0）

**步骤 2：对角缩放（D·y = w）**

$$y[i][j] = \frac{w[i][j]}{d[i][j]}$$

**步骤 3：后向扫描（L^T·z = y，按行主序逆序 (ny−2,nx−2)→(1,1)）**

$$z[i][j] = y[i][j] + \frac{\text{inv\_dx2}}{d[i][j]} \cdot z[i][j+1] + \frac{\text{inv\_dy2}}{d[i][j]} \cdot z[i+1][j]$$

（边界节点邻居贡献为 0）

### 5.4  符号推导验证

L^T 的元素：

$$(L^T)_{ij,\,ij+1} = L_{ij+1,\,ij} = \frac{(-\text{inv\_dx2})}{d[i][j]} = -\frac{\text{inv\_dx2}}{d[i][j]}$$

后向扫描展开：

$$z[i][j] = y[i][j] - (L^T)_{ij,ij+1} \cdot z[i][j+1] - (L^T)_{ij,i+1j} \cdot z[i+1][j]$$
$$= y[i][j] + \frac{\text{inv\_dx2}}{d[i][j]} \cdot z[i][j+1] + \frac{\text{inv\_dy2}}{d[i][j]} \cdot z[i+1][j] \quad \checkmark$$

### 5.5  OpenFOAM 参考

**文件**: `src/OpenFOAM/matrices/lduMatrix/preconditioners/DICPreconditioner/DICPreconditioner.C`  
**仓库**: https://github.com/OpenFOAM/OpenFOAM-dev  

关键代码段（apply，即前向+后向扫描）：

```cpp
// DICPreconditioner.C  (OpenFOAM-dev, GPL-3.0)
void Foam::DICPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction
) const
{
    scalar* __restrict__ wAPtr = wA.begin();
    const scalar* __restrict__ rAPtr = rA.begin();
    const scalar* __restrict__ rDPtr = rD_.begin();   // reciprocal of d
    const label* const __restrict__ uPtr = matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr = matrix_.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr = matrix_.upper().begin();

    wA = rA;
    const label nFaces = matrix_.upper().size();

    // Forward sweep (lower triangle)
    for (label face=0; face<nFaces; face++)
    {
        wAPtr[uPtr[face]] -= upperPtr[face]*rDPtr[lPtr[face]]*wAPtr[lPtr[face]];
    }
    // Multiply by reciprocal diagonal
    for (label cell=0; cell<wA.size(); cell++)
    {
        wAPtr[cell] *= rDPtr[cell];
    }
    // Backward sweep (upper triangle)
    for (label face=nFaces-1; face>=0; face--)
    {
        wAPtr[lPtr[face]] -= upperPtr[face]*rDPtr[uPtr[face]]*wAPtr[uPtr[face]];
    }
}
```

对应我们的公式：  
- 前向：`wA[upper] -= upper_coef * rD[lower] * wA[lower]`  
  = `w[i][j] -= (-inv_dx2) / d[west] * w[west]`  
  = `w[i][j] += (inv_dx2/d[west]) * w[west]` ✓

### 5.6  参考文献

- Kershaw, D.S. (1978). The incomplete Cholesky–conjugate gradient method for the iterative solution of systems of linear equations. *J. Comput. Phys.*, 26(1), 43–65.
- Meijerink, J.A. & van der Vorst, H.A. (1977). An iterative solution method for linear systems of which the coefficient matrix is a symmetric M-matrix. *Math. Comp.*, 31(137), 148–162.
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, §10.3.2.

---

## 6  GAMG（几何代数多重网格）预条件器

### 6.1  原理

多重网格（Multigrid, MG）是解椭圆方程的**最优**方法，理论上每次 V-cycle 收敛因子与网格大小无关。
用作 PCG 预条件器时，一次 V-cycle 替代 M⁻¹ 应用，称为 GAMG-PCG。

### 6.2  网格层次

针对 41×41 节点的均匀网格实现 **2 层 V-cycle**：

| 层次 | 节点 | 内部未知量 | 间距 |
|---|---|---|---|
| 细网格（Level 0） | 41×41 = 1681 | 39×39 = 1521 | dx = dy = 0.05 |
| 粗网格（Level 1） | 21×21 = 441 | 19×19 = 361 | dx_c = dy_c = 0.10 |

粗化方式：隔点粗化（2:1），细网格节点 (2i, 2j) 对应粗网格节点 (i, j)。

### 6.3  算子关系

粗网格 Poisson 算子（5 点差分，间距加倍）：

$$a_{\text{diag},c} = \frac{2}{dx_c^2} + \frac{2}{dy_c^2} = \frac{a_{\text{diag}}}{4}$$
$$\text{inv\_dx2\_c} = \frac{1}{dx_c^2} = \frac{\text{inv\_dx2}}{4}$$

粗网格边界条件与细网格相同类型（Neumann/Dirichlet）。

### 6.4  限制算子（Restriction）

采用**注入（Injection）**：将细网格残差直接采样到粗网格：

$$r_c[i][j] = r_f[2i][2j]$$

计算量：O(N_c)，N_c = 21×21 = 441。

### 6.5  插值算子（Prolongation）

采用**双线性插值**，将粗网格误差估计插回细网格（4 种点类型）：

| 细网格点类型 | 位置 | 插值公式 |
|---|---|---|
| 粗网格重合点（偶,偶） | (2ci, 2cj) | $e_f = e_c[ci][cj]$ |
| 水平中点（偶,奇） | (2ci, 2cj+1) | $e_f = \frac{e_c[ci][cj] + e_c[ci][cj+1]}{2}$ |
| 垂直中点（奇,偶） | (2ci+1, 2cj) | $e_f = \frac{e_c[ci][cj] + e_c[ci+1][cj]}{2}$ |
| 单元中心（奇,奇） | (2ci+1, 2cj+1) | $e_f = \frac{e_c[ci][cj] + e_c[ci][cj+1] + e_c[ci+1][cj] + e_c[ci+1][cj+1]}{4}$ |

### 6.6  V-cycle 算法

给定细网格残差 r，计算修正 z ≈ A⁻¹ r：

```
V-cycle(A, r, z):
  1. z ← 0
  2. Pre-smooth:   z ← z + ω·M_J⁻¹·(r − A·z)   [ν₁=2 次，ω=2/3]
  3. 计算残差:     res_f ← r − A·z
  4. 限制:         r_c ← R·res_f                  [注入]
  5. 粗网格求解:   A_c·e_c ← r_c                  [20 次 Jacobi]
  6. 插值+修正:    z ← z + P·e_c                  [双线性]
  7. Post-smooth:  z ← z + ω·M_J⁻¹·(r − A·z)   [ν₂=2 次，ω=2/3]
```

**阻尼 Jacobi 光滑**（ω = 2/3 为最优）：

$$z \leftarrow z + \frac{2}{3 \cdot a_{\text{diag}}} \cdot (r - Az)$$

参数 ω = 2/3 使谱半径从 1 降至约 1/3，有效消除高频误差。

### 6.7  收敛分析

理论 V-cycle 收敛因子（误差缩减比）：

$$\rho_{\text{V}} \approx \frac{1}{1 + \gamma}$$

其中 γ 取决于光滑次数 ν₁+ν₂ 和粗网格求解精度。对 ν₁=ν₂=2，典型 γ ≈ 2–5，即每次 V-cycle 误差缩减约 50–80%。

用作 PCG 预条件器时，PCG 外层迭代数通常 5–15 次（vs Jacobi PCG 的 30–80 次）。

### 6.8  OpenFOAM 参考

**文件**: `src/OpenFOAM/matrices/lduMatrix/solvers/GAMG/GAMGPreconditioner/GAMGPreconditioner.C`  
**仓库**: https://github.com/OpenFOAM/OpenFOAM-dev  

OpenFOAM GAMG 的限制/插值接口：

```cpp
// GAMGPreconditioner.C  (OpenFOAM-dev, GPL-3.0)
void Foam::GAMGPreconditioner::precondition
(
    scalarField& wA,
    const scalarField& rA,
    const direction cmpt
) const
{
    wA = 0;
    GAMG.Vcycle(interfaces, wA, rA, ...);
}
```

本实现为简化版 2-level V-cycle，与 OpenFOAM 的多层次代数 GAMG 等价于其最简退化情形。

### 6.9  参考文献

- Brandt, A. (1977). Multi-level adaptive solutions to boundary-value problems. *Math. Comp.*, 31(138), 333–390.
- Briggs, W.L., Henson, V.E. & McCormick, S.F. (2000). *A Multigrid Tutorial*, 2nd ed. SIAM.
- Trottenberg, U., Oosterlee, C.W. & Schüller, A. (2001). *Multigrid*. Academic Press.
- Wesseling, P. (1991). *An Introduction to Multigrid Methods*. Wiley.

---

## 7  收敛性对比（理论与实测预期）

对 41×41 均匀网格，Re=20，压力 Poisson 方程（dx=dy=0.05），容差 pcg_tol=1e-5：

| 预条件器 | 每步迭代数（预期） | 条件数改善 | 单次 apply 代价 |
|---|---|---|---|
| Jacobi（现有） | 30–80 | κ/κ_J ≈ 1 | O(N) 乘法 |
| DILU | 15–40 | ~2× | O(N) 乘法 |
| DIC | 10–25 | ~3–4× | O(N) + 2次三角扫描 |
| GAMG | 5–15 | ~5–10× | O(N log N) V-cycle |

注：GAMG 每次 PCG 迭代的 apply 代价比 DILU/DIC 高约 10× ，但迭代次数少约 5–10×，
总体 FLOPs 相近或略少。V-cycle 在细网格时优势更明显（N > 10⁴）。

---

## 8  实现对应关系

| 数学符号 | 代码变量 | 位置 |
|---|---|---|
| a_diag = 2/dx² + 2/dy² | `a_diag` | `pressure_poisson_pcg_prec` |
| inv_dx2 = 1/dx² | `inv_dx2` | 同上 |
| d[i][j]（修正对角） | `d_mod[i][j]` | `build_modified_diag` 返回值 |
| d⁻¹[i][j] | `d_inv[i][j]` | DILU apply 预处理 |
| w（前向扫描中间量） | `w` | `apply_dic_precond` 局部变量 |
| R（限制算子） | `gamg_restrict` | 函数 |
| P（插值算子） | `gamg_prolongate` | 函数 |
| e_c（粗网格误差） | `e_c` | `apply_gamg_precond` 局部变量 |
| ω = 2/3 | `omega = 2.0/3.0` | `gamg_smooth_fine` |
| N_c = 21 | `gamg_nc` | 模块级常量 |

---

## 9  参考文献汇总

1. Hestenes, M.R. & Stiefel, E. (1952). Methods of conjugate gradients for solving linear systems. *J. Res. NIST*, 49(6), 409–436.
2. Kershaw, D.S. (1978). The incomplete Cholesky–conjugate gradient method. *J. Comput. Phys.*, 26(1), 43–65.
3. Meijerink, J.A. & van der Vorst, H.A. (1977). An iterative solution method for symmetric M-matrices. *Math. Comp.*, 31(137), 148–162.
4. Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, 2nd ed. SIAM. §10.3.
5. Brandt, A. (1977). Multi-level adaptive solutions to boundary-value problems. *Math. Comp.*, 31(138), 333–390.
6. Briggs, W.L., Henson, V.E. & McCormick, S.F. (2000). *A Multigrid Tutorial*, 2nd ed. SIAM.
7. Trottenberg, U., Oosterlee, C.W. & Schüller, A. (2001). *Multigrid*. Academic Press.
8. OpenFOAM Foundation. *OpenFOAM-dev source code* (GPL-3.0). https://github.com/OpenFOAM/OpenFOAM-dev
   - `DILUPreconditioner.C`
   - `DICPreconditioner.C`
   - `GAMGPreconditioner.C`
