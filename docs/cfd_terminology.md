# CFD 术语层次：压力-速度耦合算法 vs 线性方程组求解器

> 记录日期：2026-04-15  
> 来源：文献调研 + OpenFOAM 源码/文档核查

---

## 1. 两层架构概述

求解不可压缩 Navier-Stokes 方程时，存在一个公认的**两层分离结构**。
OpenFOAM 的 `fvSolution` 配置文件直接体现了这一设计：

```
# OpenFOAM  fvSolution 文件
solvers {                             ← 第二层：线性方程组求解器
    p  { solver PCG;  preconditioner DIC; }
    U  { solver PBiCGStab; preconditioner DILU; }
}
SIMPLE {                              ← 第一层：压力-速度耦合算法
    nNonOrthogonalCorrectors 0;
    convergence 1e-5;
}
```

| 层次 | OpenFOAM 术语 | 学术标准英文术语 | 中文术语 | 代表方法 |
|---|---|---|---|---|
| **第一层** | Solution algorithm | Pressure-velocity coupling algorithm | 压力-速度耦合算法 | SIMPLE, PISO, PIMPLE, Chorin projection |
| **第二层** | Linear solver | Linear equation solver | 线性方程组求解器 | Gauss-Seidel, PCG, BiCGStab, GMRES, LU |

---

## 2. 第一层——压力-速度耦合算法

### 2.1 Chorin 投影法（FlowLabLite v0.0.1+）

**算法思路（算子分裂 / operator splitting）**

1. **预测步**（忽略压力梯度）：  
   计算中间速度 ũ，仅处理对流和粘性项。  
   此步为**显式有限差分更新**——用已知量直接代入公式，不产生线性系统。

2. **投影步**（Helmholtz 分解）：  
   求解压力 Poisson 方程 ∇²p = b，  
   再用压力梯度修正 ũ → u，使速度满足连续方程 ∇·u = 0。

**特点**：瞬态求解；每个时间步均需调用线性方程组求解器处理压力方程。

**原始文献：**
> Chorin, A. J. (1968). Numerical solution of the Navier-Stokes equations.  
> *Mathematics of Computation*, **22**(104), 745–762.  
> DOI: [10.2307/2004755](https://doi.org/10.2307/2004755)

同期独立工作（Temam 分步法）：
> Temam, R. (1969). Sur l'approximation de la solution des équations de Navier-Stokes par la méthode des pas fractionnaires.  
> *Archiv for Rational Mechanics and Analysis*, **33**, 377–385.  
> DOI: [10.1007/BF00247696](https://doi.org/10.1007/BF00247696)

---

### 2.2 SIMPLE 算法（FlowLabLite v0.0.2+）

**算法思路（预测-修正迭代 / predictor-corrector iteration）**

1. **动量预测**（Momentum predictor）：  
   用当前压力场 p\* 显式预测速度 u\*, v\*，施加欠松弛（α_u = 0.7）。  
   此步同样为**显式更新**，不产生线性系统。

2. **压力修正源项**：  
   由 ∇·u\* 构造压力修正方程的右端项 b。

3. **压力修正 Poisson 求解**：  
   ∇²p' = b，调用线性方程组求解器（Gauss-Seidel，50 次内迭代）。

4. **场量更新**：  
   p ← p\* + α_p · p'（α_p = 0.3），u/v 用 p' 梯度修正。

5. **边界条件 + 收敛检查**：残差 = mean|∂u/∂x + ∂v/∂y|。

**特点**：稳态迭代；每次 SIMPLE 扫描只需一次 Poisson 求解。

**原始文献：**
> Patankar, S. V., & Spalding, D. B. (1972). A calculation procedure for heat, mass and momentum transfer in three-dimensional parabolic flows.  
> *International Journal of Heat and Mass Transfer*, **15**(10), 1787–1806.  
> DOI: [10.1016/0017-9310(72)90054-3](https://doi.org/10.1016/0017-9310(72)90054-3)

教科书标准参考（必读）：
> Patankar, S. V. (1980). *Numerical Heat Transfer and Fluid Flow.*  
> Hemisphere Publishing. ISBN: 978-0891165224

---

### 2.3 其他常见耦合算法（供对比）

| 算法 | 提出者 | 文献 | 适用场景 |
|---|---|---|---|
| **SIMPLEC** | Van Doormaal & Raithby (1984) | *Num. Heat Transfer*, 7(2), 147–163 | 稳态，收敛更快 |
| **PISO** | Issa (1986) | *J. Comput. Phys.*, 62(1), 40–65 | 瞬态隐式 |
| **PIMPLE** | OpenFOAM 混合算法 | — | 大时间步瞬态 |

---

## 3. 第二层——线性方程组求解器

离散化之后，无论 Chorin 还是 SIMPLE，压力（修正）Poisson 方程都归结为代数方程组 **Ax = b**，
此时才需要选择线性方程组求解器：

| 类别 | 代表方法 | OpenFOAM 对应 | FlowLabLite |
|---|---|---|---|
| **直接法** | LU 分解、高斯消元 | 无（大规模太贵） | 无 |
| **定常迭代法** | Gauss-Seidel, Jacobi, SOR | `smoothSolver` | ✅ 压力 Poisson（`nit=50` 次） |
| **Krylov 子空间法** | CG, PCG, GMRES, BiCGStab | `PCG`, `PBiCGStab` | 无 |
| **多重网格法** | AMG, GMG | `GAMG` | 无 |

**线性求解器参考文献：**
> Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.).  
> SIAM. ISBN: 978-0898715347  
> — CG / GMRES 的权威参考

> Ferziger, J. H., Perić, M., & Street, R. L. (2020). *Computational Methods for Fluid Dynamics* (4th ed.).  
> Springer. ISBN: 978-3319996929  
> — 涵盖两层结构的标准 CFD 教材

---

## 4. FlowLabLite 两层结构全景

```
┌─────────────────────────────────────────────────────────────────┐
│              第一层：压力-速度耦合算法                            │
│  ┌─────────────────────────────┬──────────────────────────────┐ │
│  │  Chorin 投影法（v0.0.1+）   │  SIMPLE（v0.0.2+）           │ │
│  │  外层：显式时间推进          │  外层：预测-修正迭代          │ │
│  │  500 时间步                 │  100 次迭代                  │ │
│  └─────────────────────────────┴──────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
                          ↓ 压力(修正) Poisson 方程 Ax = b
┌─────────────────────────────────────────────────────────────────┐
│              第二层：线性方程组求解器                             │
│  ┌─────────────────────────────┬──────────────────────────────┐ │
│  │  速度动量方程：              │  速度预测步：                │ │
│  │  显式代入（无线性系统）       │  显式代入（无线性系统）       │ │
│  │                             │                              │ │
│  │  压力 Poisson ∇²p = b：      │  压力修正 ∇²p' = b：         │ │
│  │  Gauss-Seidel × 50 次       │  Gauss-Seidel × 50 次        │ │
│  └─────────────────────────────┴──────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────┘
```

两个耦合算法共享同一个线性求解器实现（`pressure_poisson_array`，Gauss-Seidel）。

---

## 5. 关于"显式更新"与"直接法"的区别

CFD 中容易混淆的两个概念：

| 概念 | 正确英文术语 | 含义 | FlowLabLite 使用情况 |
|---|---|---|---|
| **直接法** | Direct method | 线性代数意义：LU/高斯消元，有限步精确解 Ax=b | 未使用 |
| **显式更新** | Explicit scheme / explicit update | 数值格式意义：用已知时间层的值直接计算新值，不产生线性系统 | ✅ 速度动量方程 |

**FlowLabLite 速度更新是显式格式（explicit scheme），不是直接法（direct method）。**

显式格式的限制：时间步长受 CFL 条件约束（稳定性要求 dt 不能太大）。
若改用隐式格式（Crank-Nicolson 等），则速度方程也会产生 Av = b，需要线性求解器。

---

## 6. OpenFOAM 对应类名一览

| 层次 | OpenFOAM 类 / 关键字 |
|---|---|
| 耦合算法控制 | `SIMPLEControl`, `PIMPLEControl`, `PISOControl` |
| 线性求解器（压力） | `PCG`（预处理共轭梯度）, `GAMG`（代数多重网格） |
| 线性求解器（速度） | `PBiCGStab`（预处理双共轭梯度稳定化）, `smoothSolver` |
| 光顺子程序 | `GaussSeidel`, `symGaussSeidel`, `DIC`, `DILU` |

`GAMG` 使用 `GaussSeidel` 作为光顺子程序（smoother），
与 FlowLabLite 直接用 Gauss-Seidel 迭代 50 次的方式等价（但 AMG 效率更高）。

---

## 7. 完整术语体系（Ferziger & Perić 分类框架）

```
不可压缩 N-S 方程求解
│
├── 第一层：压力-速度耦合算法（Pressure-velocity coupling algorithm）
│   ├── 投影法 / 分步法 ─── Chorin (1968) ────── 瞬态显式时间推进
│   ├── SIMPLE ──────────── Patankar & Spalding (1972) ── 稳态迭代
│   ├── SIMPLEC ─────────── Van Doormaal & Raithby (1984)
│   └── PISO ────────────── Issa (1986) ─────── 瞬态隐式
│
└── 第二层：线性方程组求解器（Linear equation solver）
    ├── 直接法：高斯消元、LU 分解
    ├── 定常迭代法：Gauss-Seidel、Jacobi、SOR（逐次超松弛）
    ├── Krylov 子空间法：CG、GMRES、BiCGStab
    └── 多重网格法：AMG（代数多重网格）、GMG（几何多重网格）
```
