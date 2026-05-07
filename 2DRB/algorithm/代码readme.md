# 热对流 LBM 算法文献与代码分支对应说明

本文档整理三篇算法文献和代码分支的对应关系，供后续逐行检查
`2DRBOpenmp.F90`、`2DRBOpenacc.F90`、`3DRBOpenmp.F90`、
`3DRBOpenacc.F90` 和 `3DRBOpenaccMpi.F90` 时使用。这里暂时不判定代码一定正确，只建立检查基准。

## 文献与分支

| 文献 | 主要用途 | 对应代码分支 |
| --- | --- | --- |
| `Interpolation-supplemented lattice Boltzmann simulation of thermal convection on non-uniform meshes.pdf` | 给出热对流 DDF-MRT 框架，包含 2D/3D、D2Q9/D3Q19 流场、D2Q5/D3Q7 温度场，以及非均匀网格 ISLBM 插值步骤 | `EnableLegacyThermalScheme`；当前代码采用其中均匀网格退化形式，不采用非均匀网格插值 |
| `Lattice Boltzmann simulations of three-dimensional thermal convective flows at high Rayleigh number.pdf` | 3D 高 Rayleigh 数热对流算法，流场 D3Q19 MRT，温度场 D3Q7 MRT，高阶误差各向同性约束 | `EnableLegacyThermalScheme` 的 3D 版本 |
| `Multiple-relaxation-time lattice Boltzmann method for the Navier-Stokes and nonlinear convection-diffusion equations_ Modeling, analysis, and elements.pdf` | 主要参考温度场 CDE 的 MRT 统一框架，尤其是辅助分布函数 `G_i` 和 `B = uT` 的历史通量修正 | `EnableUseG`；流场仍沿用前两篇的 Guo forcing / D2Q9 或 D3Q19 MRT |

需要注意的是，`EnableUseG` 分支不是严格复现 Xu/ISLBM 系列文献中的温度 MRT 模型，而是只沿用其流场 D2Q9/D3Q19 MRT、Guo forcing、速度边界和热对流物理模型；温度场则替换为 Chai-Shi NCDE MRT 框架中的辅助分布函数 `G_i` 修正模型。

因此后续检查代码时不能把 `EnableLegacyThermalScheme` 和 `EnableUseG` 的温度参数混用：

- legacy 分支使用 `aT`、`Qk = 3 - sqrt(3)`、`Qe = Qnu = 4*sqrt(3) - 6`。
- UseG 分支使用固定温度权重、`csT^2`、`taug = 1/2 + kappa/csT^2` 和历史通量修正 `dB`。

## 共同物理模型

三篇文献的热对流基本模型都可以按 Boussinesq 近似理解：

```text
div(u) = 0
du/dt + u · grad(u) = -grad(p) + nu laplacian(u) + g beta (T - Tref) e_g
dT/dt + u · grad(T) = kappa laplacian(T)
```

无量纲参数：

```text
Ra = g beta DeltaT L^3 / (nu kappa)
Pr = nu / kappa
```

代码检查时要把这些量一起核对：`Rayleigh`、`Prandtl`、`Mach`、`tauf`、
`viscosity`、`diffusivity`、`gBeta`、`Tref`、`lengthUnit`、`timeUnit`、
`velocityUnit`。

`Tref` 或 `T0` 主要用于浮力项中的温度偏差：

```text
F = rho * gBeta * (T - Tref) * e_g
```

以及无量纲温度的定义：

```text
T* = (T - Tref) / DeltaT
```

但 Nusselt 数和壁面热流的核心是温度梯度和壁温差：

```text
DeltaT = Thot - Tcold
```

不能把 `Tref` 当成 `DeltaT` 使用。计算 Nu 时应检查：

- 温度梯度方向是否与冷热壁方向一致。
- 分母是否使用 `DeltaT`。
- `Tref` 是否只进入浮力项或无量纲温度定义，而不是被误用为壁温差。

## 流场算法

流场在两个温度分支中应保持同一个算法族：

- 2D 使用 D2Q9 MRT。
- 3D 使用 D3Q19 MRT。
- 流场采用低 Mach 数弱可压缩 MRT-LBM。在 `Ma` 足够小的条件下，该模型近似恢复不可压缩 Boussinesq 热对流方程。
- 当前代码按全密度 `rho` 编写，不是按密度波动 `delta rho = rho - rho0` 编写。
- 浮力项按 Guo forcing 思路进入碰撞项，力为 `rho * gBeta * (T - Tref)` 乘以重力方向。
- 宏观速度恢复需要包含半步力修正。

当前检查基准应按全密度 D2Q9/D3Q19 MRT-LBM 理解：`rho` 是参与平衡矩、力项和速度恢复的完整密度变量。Xu 系列文献会报告 density fluctuation 和 velocity divergence，但这里的 density fluctuation 是诊断量，不是代码中的主变量。如果以后另建 density-fluctuation equilibrium 或 incompressible-equilibrium 版本，需要单独注明。

宏观量恢复的基准形式：

```text
rho = sum_i f_i
u = (sum_i e_i f_i + F / 2) / rho
```

流场压力通常由 LBM 状态方程近似给出，例如：

```text
p = rho * cs^2
```

因此密度涨落虽然很小，但并不严格为零。检查代码时需要确认：

- `rho = sum_i f_i` 是否参与平衡矩。
- 速度恢复是否包含半步力修正。
- 不要把当前代码里的 `rho` 误读成 `delta rho`。

流场松弛参数的核心关系：

```text
nu = cs^2 * (tau_f - 1/2) * dt
Snu = 1 / tau_f
Sq, Sm = 8 * (2*tau_f - 1) / (8*tau_f - 1)
```

注意守恒模的处理。Xu 系列 MRT 文献中，密度和动量守恒模通常取：

```text
s_rho = 0
s_j   = 0
```

这表示碰撞不改变守恒矩。若代码中守恒模取 `1` 或其他数值，只要碰撞实现中守恒矩的 post-collision 值被强制保持为 equilibrium/conserved value，宏观守恒性通常仍可成立；但这已经不是文献逐项复现的写法。

因此查代码时要区分两件事：

1. 守恒矩是否真的被碰撞改变。
2. 松弛矩阵中守恒模的数值是否和文献一致。

后续查代码时，流场重点不是 `EnableLegacyThermalScheme` 或 `EnableUseG` 的差异，而是确认两个温度分支没有意外改变流场碰撞、力项符号、速度恢复和边界处理。

## `EnableLegacyThermalScheme` 温度场

这个分支对应 Xu/ISLBM 系列和 2019 三维高 Rayleigh 数热对流文献中的温度 MRT 模型。

### 2D D2Q5

温度分布函数为 `g_i`，宏观温度为：

```text
T = sum_i g_i
```

D2Q5 温度矩可以按下面的结构检查：

```text
n = [T, uT, vT, aT*T, 0]^T
aT = 20 * sqrt(3) * kappa - 4
```

权重应为：

```text
omega_0 = (1 - aT) / 5
omega_1..4 = (4 + aT) / 20
```

D2Q5 温度松弛矩阵建议按下面形式检查：

```text
Q = diag(0, Qk, Qk, Qe, Qnu)
```

其中：

```text
Qk  = 3 - sqrt(3)
Qe  = 4 * sqrt(3) - 6
Qnu = 4 * sqrt(3) - 6
```

这里：

- 第 0 个矩 `T` 是守恒温度矩，不参与松弛。
- 第 1、2 个矩 `uT`、`vT` 是温度通量矩，其松弛率 `Qk` 决定热扩散率。
- 第 3 个矩 `aT*T` 的松弛率应记为 `Qe`。
- 第 4 个高阶矩的松弛率记为 `Qnu`。
- 如果代码只定义一个 `Qnu` 并同时用于第 3、4 个矩，那么需要确认它确实等价于 `Qe = Qnu`，而不是漏掉第 3 个矩。

恒温边界采用 half-way anti-bounce-back：

```text
g_opp(x_f, t + dt) = -g_i^+(x_f, t) + (4 + aT) / 10 * T_wall
```

绝热边界采用 half-way bounce-back：

```text
g_opp(x_f, t + dt) = g_i^+(x_f, t)
```

### 3D D3Q7

D3Q7 是 2D D2Q5 的维度拓展，温度矩应为：

```text
n = [T, uT, vT, wT, aT*T, 0, 0]^T
aT = 42 * sqrt(3) * kappa - 6
```

权重应为：

```text
omega_0 = (1 - aT) / 7
omega_1..6 = (6 + aT) / 42
```

D3Q7 温度松弛矩阵应按下面形式检查：

```text
Q = diag(0, Qk, Qk, Qk, Qe, Qnu, Qnu)
```

其中：

```text
Qk  = 3 - sqrt(3)
Qe  = 4 * sqrt(3) - 6
Qnu = 4 * sqrt(3) - 6
```

这里：

- `T` 是守恒温度矩。
- `uT`、`vT`、`wT` 是三个方向的温度通量矩。
- `aT*T` 对应 `Qe`。
- 最后两个二阶各向异性矩对应 `Qnu`。
- 如果代码中 `Qe` 没有单独命名，而是直接复用 `Qnu`，需要确认第 4、5、6 个非通量矩的松弛率都符合文献中的 `Qe = Qnu`。

恒温边界采用：

```text
g_opp(x_f, t + dt) = -g_i^+(x_f, t) + (6 + aT) / 21 * T_wall
```

绝热边界仍采用 half-way bounce-back。

### legacy 温度边界系数

legacy 恒温边界系数来自 half-way anti-bounce-back 的统一形式：

```text
g_opp = -g_i^+ + 2 * omega_axis * T_wall
```

对于 legacy 温度权重，D2Q5 为：

```text
omega_axis = (4 + aT) / 20
2 * omega_axis = (4 + aT) / 10
```

D3Q7 为：

```text
omega_axis = (6 + aT) / 42
2 * omega_axis = (6 + aT) / 21
```

因此 `(4+aT)/10` 和 `(6+aT)/21` 只适用于 `EnableLegacyThermalScheme`，不能直接用于 `EnableUseG`。

### 均匀网格退化

`Interpolation-supplemented...` 这篇文献的非均匀网格部分主要在 streaming 后加入虚拟格点和二次插值。当前代码如果仍采用均匀网格，则应退化为标准 LBM：

```text
collision: g_i^+(x, t)
streaming: g_i(x + e_i dt, t + dt) = g_i^+(x, t)
boundary: half-way bounce-back / anti-bounce-back
```

当前代码采用均匀网格时，不需要寻找 ISLBM 的二次插值、非均匀坐标、非均匀加权平均或非均匀壁面导数。

### 网格设置与 half-way 边界位置

文献中的 half-way bounce-back / anti-bounce-back 设置意味着：数组中的内部格点是流体节点，真实物理固壁不在第一层流体节点上，而是在第一层流体节点外侧半个格距的位置。

对当前均匀网格代码，建议按下面几何关系理解：

```text
left/bottom wall   : x = 0
first fluid node   : x = 0.5 * dx_lu
fluid node i       : x_i = i - 0.5    in lattice units
last fluid node    : x = N - 0.5
right/top wall     : x = N
```

如果无量纲化后令壁面距离为 1，则坐标可写成：

```text
x_wall_left  = 0
x_wall_right = 1
x_i          = (i - 0.5) / N,   i = 1, ..., N
dx           = 1 / N
```

这对应代码中的典型设置：

```text
xp(0)    = 0
xp(nx+1) = nx
xp(i)    = i - 0.5
xp       = xp / lengthUnit
```

`yp`、`zp` 同理。因此：

- `nx`、`ny`、`nz` 是流体节点数，不是包含边界节点后的总点数。
- 物理域长度按 `N` 个 lattice unit 理解，而不是 `N-1`。
- `lengthUnit = nx` 或 `ny` 是和 half-way 边界位置配套的。
- `xp(0)`、`xp(nx+1)`、`yp(0)`、`yp(ny+1)` 是物理壁面坐标，不是流体节点。
- `xp(1)`、`yp(1)` 是第一层流体节点，距离对应物理壁面半个格距。
- 初始化温度线性场、中心线插值、壁面热流和 Nu 计算都必须使用这个几何关系。

但这不代表壁面热流可以随意用普通中心差分。由于速度和温度边界采用 half-way bounce-back / anti-bounce-back，第一层流体节点距离物理壁面是半个格距。因此计算壁面温度梯度和 Nusselt 数时，仍应使用与 half-way 边界一致的均匀网格壁面导数公式。

例如 bottom/hot wall 的温度梯度不能简单写成：

```text
(T(2) - T(1)) / dx
```

而应根据壁面位置和第一、第二层流体节点位置构造 one-sided 二阶公式。若代码中采用文献给出的均匀网格 half-way 公式，需要确认：

```text
dT/dn_wall = (-8*Twall + 9*T1 - T2) / (3*dx)
```

或对应方向上的符号版本。这里 `T1` 是离壁面最近的第一层流体节点，`T2` 是第二层流体节点。

## `EnableUseG` 温度场

这个分支对应 Chai-Shi 2020 文献中用于非线性对流扩散方程的 MRT 框架。流场仍沿用前两篇文献的 Navier-Stokes MRT 和 Guo forcing；主要变化在温度场。

目标温度方程可写成：

```text
partial_t T + div(B) = div(kappa grad(T))
B = u*T
```

在不可压缩条件 `div(u)=0` 下，这等价于常见的：

```text
partial_t T + u · grad(T) = kappa laplacian(T)
```

Chai-Shi 框架中的关键点是引入辅助分布函数 `G_i`，用于消除离散化带来的附加项。对当前温度方程，常用理解是：

```text
g_i^eq = omega_i * T * (1 + e_i · u / csT^2)
G_i    = omega_i * e_i · [(I - S1/2) partial_t B] / csT^2
```

其中：

```text
B = u*T
```

代码中通常不显式存储完整 `G_i`，而是在温度一阶矩碰撞后加入历史通量修正。

如果代码中定义：

```text
dB = B(t) - B(t-dt)
```

则因为演化方程中实际加入的是 `dt * G_i`，在 lattice unit `dt = 1` 下，moment-space 的一阶通量修正可以写成：

```text
n_flux^post = n_flux - Qk * (n_flux - B) + (1 - 0.5*Qk) * dB
```

如果代码中定义：

```text
dBdt = [B(t) - B(t-dt)] / dt
```

则修正项必须写成：

```text
n_flux^post = n_flux - Qk * (n_flux - B) + dt * (1 - 0.5*Qk) * dBdt
```

因此检查代码时要先确认变量名 `dB`、`dBx`、`dBy`、`dBz` 存的是“差值”还是“时间导数”。不能同时除以 `dt` 又省略前面的 `dt`，否则会多乘或少乘一个时间步因子。

若 `EnableUseG` 关闭，则 `dB = 0`，该辅助项消失。

### 2D D2Q5 `EnableUseG`

建议按下面结构理解：

```text
omega_0 = 1 / 3
omega_1..4 = 1 / 6
csT^2 = 1 / 3
taug = 1/2 + kappa / csT^2
Qk = 1 / taug
Qnu = 1
```

这里的 `Qnu = 1` 不是由热扩散率唯一决定的，而是当前实现中对非守恒高阶矩的快速松弛选择。Chai-Shi NCDE 框架主要约束一阶通量矩的松弛率 `S1`，从而确定热扩散率：

```text
kappa = dt * beta * csT^2 * (1/Qk - 1/2)
```

在当前简化温度方程中取：

```text
beta = 1
dt = 1
```

于是：

```text
kappa = csT^2 * (1/Qk - 1/2)
taug  = 1/Qk = 1/2 + kappa/csT^2
```

D2Q5 UseG 分支若采用：

```text
omega_0 = 1/3
omega_axis = 1/6
csT^2 = 1/3
```

则平衡分布可以写为：

```text
g_i^eq = omega_i * T * (1 + e_i · u / csT^2)
```

在常见 D2Q5 矩阵下，对应平衡矩为：

```text
n_eq = [T, uT, vT, -2T/3, 0]^T
```

一阶温度通量矩 `uT`、`vT` 加入：

```text
Bx = u*T
By = v*T

dBx = Bx(t) - Bx(t-dt)
dBy = By(t) - By(t-dt)

n_xflux^post = n_xflux - Qk*(n_xflux - Bx) + (1 - 0.5*Qk)*dBx
n_yflux^post = n_yflux - Qk*(n_yflux - By) + (1 - 0.5*Qk)*dBy
```

### 3D D3Q7 `EnableUseG` 的实现性拓展

Chai-Shi 文献给出的是一般维度的 NCDE MRT 框架。当前代码若在三维温度场中采用 D3Q7，则可以按 `DdQ(2d+1)` 的思想从 D2Q5 拓展：

```text
omega_0 = 1/4
omega_1..6 = 1/8
csT^2 = 1/4
```

热扩散率由一阶通量矩松弛率决定：

```text
kappa = csT^2 * (1/Qk - 1/2)
taug  = 1/Qk = 1/2 + kappa/csT^2
Qk    = 1/taug
```

高阶非守恒矩可采用快速松弛，例如：

```text
Qnu = 1
```

在与 legacy D3Q7 类似的矩阵结构下，平衡矩应为：

```text
n_eq = [T, uT, vT, wT, -3T/4, 0, 0]^T
```

三个一阶温度通量矩分别加入历史通量修正：

```text
Bx = u*T
By = v*T
Bz = w*T

dBx = Bx(t) - Bx(t-dt)
dBy = By(t) - By(t-dt)
dBz = Bz(t) - Bz(t-dt)

n_xflux^post = n_xflux - Qk*(n_xflux - Bx) + (1 - 0.5*Qk)*dBx
n_yflux^post = n_yflux - Qk*(n_yflux - By) + (1 - 0.5*Qk)*dBy
n_zflux^post = n_zflux - Qk*(n_zflux - Bz) + (1 - 0.5*Qk)*dBz
```

检查代码时必须确认：

- `wT` 是否位于 D3Q7 温度矩的第 3 个通量矩。
- `dBz` 是否加到 `wT` 对应的矩上。
- `ex/ey/ez` 的编号和 `N`、`N^{-1}`、`opp` 方向编号完全一致。
- 重启文件中是否保存和恢复了 `Bx_prev`、`By_prev`、`Bz_prev`。

### `EnableUseG` 温度边界条件

`EnableUseG` 分支采用固定温度权重，不再使用 `aT`。因此恒温边界的 half-way anti-bounce-back 仍为统一形式：

```text
g_opp = -g_i^+ + 2 * omega_axis * T_wall
```

但具体系数变为：

#### 2D D2Q5 `EnableUseG`

```text
omega_axis = 1/6
2 * omega_axis = 1/3

g_opp = -g_i^+ + (1/3) * T_wall
```

#### 3D D3Q7 `EnableUseG`

```text
omega_axis = 1/8
2 * omega_axis = 1/4

g_opp = -g_i^+ + (1/4) * T_wall
```

检查代码时必须确认：

- `EnableLegacyThermalScheme` 使用 `(4+aT)/10` 或 `(6+aT)/21`。
- `EnableUseG` 使用 `1/3` 或 `1/4`。
- 两个分支不能共用同一个 ABB 系数，除非代码内部已经按 `2*omega_i` 自动计算。

### D3Q7 矩阵和速度编号

三维温度场最容易出错的地方不是公式，而是编号和符号。D3Q7 的速度编号、反方向编号、矩阵 `N` 和逆矩阵 `N^{-1}` 必须逐项一致。特别是第三个方向 `ez` 的正负号、`wT` 所在矩、以及 `dBz` 的符号要重点核对。只要 `ez` 行符号和速度编号不一致，三维温度通量就会出现方向性错误。

如果代码采用如下速度编号：

```text
e0 = ( 0,  0,  0)
e1 = ( 1,  0,  0)
e2 = (-1,  0,  0)
e3 = ( 0,  1,  0)
e4 = ( 0, -1,  0)
e5 = ( 0,  0,  1)
e6 = ( 0,  0, -1)
```

则自然的反方向关系是：

```text
opp(0)=0
opp(1)=2
opp(2)=1
opp(3)=4
opp(4)=3
opp(5)=6
opp(6)=5
```

需要确认 `N` 矩阵的 `wT` 行和这个速度编号一致。如果文献矩阵或代码矩阵中 `ez` 行写成相反符号，只要 `N^{-1}` 和平衡矩、streaming 方向保持自洽，算法仍可运行；但若只有一处符号改变，就会出错。

## 边界条件

速度边界：

- 固壁采用 half-way bounce-back，对应无滑移边界。
- 对 D2Q9/D3Q19，需确认 `opp` 方向编号和 streaming 顺序一致。

温度边界：

- 恒温边界采用 half-way anti-bounce-back。
- 绝热边界采用 half-way bounce-back。

恒温边界的统一写法是：

```text
g_opp = -g_i^+ + 2 * omega_axis * T_wall
```

因此不同温度分支的系数不同。

### legacy 温度模型

D2Q5 legacy:

```text
omega_axis = (4 + aT)/20
g_opp = -g_i^+ + (4 + aT)/10 * T_wall
```

D3Q7 legacy:

```text
omega_axis = (6 + aT)/42
g_opp = -g_i^+ + (6 + aT)/21 * T_wall
```

### UseG 温度模型

D2Q5 UseG:

```text
omega_axis = 1/6
g_opp = -g_i^+ + (1/3) * T_wall
```

D3Q7 UseG:

```text
omega_axis = 1/8
g_opp = -g_i^+ + (1/4) * T_wall
```

绝热边界统一为：

```text
g_opp = g_i^+
```

后续查代码时要按算例宏区分：

- `SideHeatedCell`: 左右壁恒温，其他壁绝热，浮力方向按文献/代码设置核对。
- `RayleighBenardCell`: 下热上冷，水平或侧向周期/绝热设置要与代码宏一致。

边界检查不能只看分布函数公式，还要同时检查：

- 边界后的 `T` 是否由 `sum_i g_i` 恢复。
- `opp` 方向编号是否正确。
- 热流和 Nusselt 数计算是否使用同一个壁面方向。
- `Thot`、`Tcold`、`Tref` 和 `DeltaT` 是否被一致使用。
- `EnableLegacyThermalScheme` 和 `EnableUseG` 是否使用了各自正确的 ABB 系数。

## 后处理基准：2D 热驱动方腔稳态文献

后处理主要参考 `Lattice-Boltzmann simulations of the thermally driven 2D square cavity at high Rayleigh numbers.pdf`。这篇文献对应的是 2D side-heated square cavity：左壁热、右壁冷、上下绝热、四壁无滑移，典型参数是 `Pr = 0.71` 和高 `Ra`。

需要特别注意：这篇文献的后处理只针对稳态解。它没有给出非稳态或湍流 Rayleigh-Benard 算例的时间平均、窗口平均、统计稳态判据或采样策略。因此它可以作为 `steadyFlow` 和稳态 side-heated benchmark 的 Nu/Re/streamfunction 检查基准；当前代码中的稳态 RB 后处理也沿用这一套稳态热流诊断思路，只是把主热流方向从 `x` 方向旋转到 `y` 方向。它不能直接作为 `unsteadyFlow` 后处理依据。

### 稳态判据

文献中使用两个稳态检查：

```text
||u(t+1000)-u(t)||_2 / ||u(t+1000)||_2 < 1e-12
||theta(t+1000)-theta(t)||_infty < 1e-6
```

这里的 `theta` 是无量纲温度，文献采用 `theta_h = 0.5`、`theta_c = -0.5`，所以 `Delta theta = 1`。代码检查时不要求阈值必须逐字相同，但要确认它们表达的是稳态收敛，而不是非稳态统计收敛。

### Nusselt 数定义

对 side-heated cavity，热通量方向为 `x`，局部水平热流为：

```text
q_x = u * theta - dtheta/dx
```

文献比较三个 Nusselt 数：

```text
Nu      = 1/(L*H*Delta theta) * int int q_x dx dy
Nu_0    = 1/(H*Delta theta)   * int q_x(x=0,y) dy
Nu_1/2  = 1/(H*Delta theta)   * int q_x(x=1/2,y) dy
```

稳态且能量守恒一致时，`Nu`、`Nu_0` 和 `Nu_1/2` 应趋于同一个值。三者差异可以作为热量守恒和后处理一致性的诊断量。

在内部流体节点，文献使用中心差分：

```text
q_x(i,j) = u(i,j)*theta(i,j) + [theta(i-1,j)-theta(i+1,j)]/(2*dx)
```

第一层流体节点距离左壁半个格距，因此邻近热壁的热流使用 half-way 几何下的一侧二阶公式：

```text
q_x(1,j) = u(1,j)*theta(1,j)
         + [4*theta_h - 3*theta(1,j) - theta(2,j)]/(3*dx)
```

这一条公式来自左壁、第一层流体节点和第二层流体节点的 half-way 二次插值，并与 `q_x = u*theta - dtheta/dx` 一致。文献 OCR/排版中如果出现 `[4*theta(1,j)-theta(2,j)-3*theta_h]/(3*dx)`，按线性温度场检验会给出错误符号和量级，应作为疑似笔误一并警惕。

### 稳态 RB 的方向替换

对于稳态 Rayleigh-Benard 算例，当前代码仍可按同一套稳态热流后处理框架理解，只是主热流方向从 side-heated cavity 的 `x` 方向换成竖直 `y` 方向。也就是说：

```text
side-heated cavity: q_x = u * theta - dtheta/dx
Rayleigh-Benard  : q_y = v * theta - dtheta/dy
```

内部节点的对应形式为：

```text
q_y(i,j) = v(i,j)*theta(i,j) + [theta(i,j-1)-theta(i,j+1)]/(2*dy)
```

如果下壁为热壁、上壁为冷壁，且正的 Nusselt 数表示热量向上输运，则 bottom/hot wall 的 half-way 壁面热流应按：

```text
q_y(y=0,i) = 2 * [theta_h - theta(i,1)] / dy
```

top/cold wall 若也输出为正的向上热流，则对应：

```text
q_y(y=1,i) = 2 * [theta(i,ny) - theta_c] / dy
```

中心线 Nu 则从 side-heated 的 `x = 1/2` 竖直中心线平均，换成 RB 的 `y = 1/2` 水平中心线平均。后续查代码时，`SideHeatedCell` 和 `RayleighBenardCell` 的 Nu 公式应只在热流方向、速度分量和冷热壁位置上不同，而不应混入不同的温差尺度或不同的 `Tref` 用法。

### 热壁 Nu 笔误检查

文献计算左热壁 `Nu_0` 时给出的壁面热流公式印成了：

```text
q_x(x=0,j) = 2 * [theta(1,j) - theta_h] / dx
```

按文献前面定义的 `q_x = u*theta - dtheta/dx`，左热壁处 `u=0`，温度从热壁向腔内降低，所以 `dtheta/dx < 0`，热流 `q_x` 应为正。与这个符号约定一致的公式应写成：

```text
q_x(x=0,j) = 2 * [theta_h - theta(1,j)] / dx
```

因此后续检查代码时，应把文献 Eq. (44) 的符号当作疑似排版笔误处理。若代码输出的是正的热壁 Nusselt 数，左热壁应使用 `theta_h - theta(1,j)` 这个方向；如果代码对冷壁也输出正值，则还要确认冷壁的符号约定是“从热壁到冷壁的总热流”为正，还是按局部外法向热流为正。

### 角点、极值和流函数

文献对左壁上下角点的处理也和 half-way 边界有关。由于上下壁是绝热边界，角点附近的壁面温度不是直接给定，而是用靠近角点的四个温度点沿竖直方向拟合：

```text
theta(y) = c0 + c2 * (y - y0)^2
```

然后再用壁面热流公式计算角点附近的 `q_x`。如果当前代码在 wall-average Nu 中排除了角点，或者只使用第一层流体节点平均，需要在后续代码检查时单独标明。

文献还报告热壁上的局部极值：

```text
Nu_max, y_max
Nu_min, y_min
```

`Nu_max` 的位置用热壁上五个相邻点做抛物线最小二乘插值；`Nu_min` 靠近左上角，文献说明只使用极小点之前一侧的点做插值。代码如果输出 `Nu_hot_max/min`，应检查它是否采用相同的插值思想，还是仅取离散网格点极值。

流动强度方面，文献报告：

```text
psi_mid
|psi|_max and its location
u_max on the vertical centerline
v_max on the horizontal centerline
```

流函数定义为：

```text
psi(x,y) = int_0^y u(x,mu) dmu = - int_0^x v(eta,y) deta
```

若要严格复现文献，`|psi|_max` 需要先对 LBM 数据做插值，再在更细网格上寻找极值。若当前代码只输出网格点上的流函数极值，后续 benchmark 对比时应降低为“近似诊断”，不要视为逐项复现。

## 稳态与非稳态检查基准

稳态文献算例常用速度和温度的相对变化作为收敛准则，典型检查间隔是若干千步。代码中 `steadyFlow` 应重点核对：

- `errorU` 和 `errorT` 是否分别对应速度场和温度场相对变化。
- 检查间隔是否和输出/采样间隔混在一起。
- 达到稳态后最终输出、Nu/Re 计算、reload 文件写出是否完整。

非稳态或湍流 Rayleigh-Benard 算例不能用单次稳态误差作为终止依据。更合理的检查基准是：

- 先经过足够长的初始瞬态。
- 再进入统计采样。
- 采样窗口内检查 Nu、Re 或其他全局量的时间平均是否收敛。
- 输出的 Nu/Re 应明确是瞬时值、窗口平均值还是累计平均值。

## OpenMP / OpenACC 并行实现检查

并行版本的检查原则是：先确认它和串行/CPU 基准保持同一个算法顺序，再检查并行语义是否改变了数据依赖。`OpenMP` 和 `OpenACC` 都不应该改变下面这个主时间步结构：

```text
collision
streaming
boundary treatment
macro reconstruction
diagnostics / output / reload
```

如果某个并行版本为了性能调整了循环组织，也必须保证这些阶段之间的依赖关系不被打乱。尤其不能让下一阶段读取还没有完成更新的 `f/g/T/u/v/w/rho`，也不能让后处理读取 device 上尚未同步回 host 的数据。

### OpenMP 检查

OpenMP 版本重点检查线程作用域和 reduction：

- `parallel do` 或 `do` 循环内的局部临时变量必须是私有的，例如矩空间数组、平衡矩、碰撞后矩、source term、局部速度/温度中间量。
- 写数组时要确认每个线程只写自己的格点或自己的方向分量，不产生多个线程写同一个元素的竞争。
- 如果 streaming 采用 pull 方式，应确认读写数组分离或读写方向不会相互覆盖；如果采用 push 方式，更要检查是否有同一目标格点被多个源格点写入。
- `errorU`、`errorT`、`Nu`、`Re`、体平均热流、最大/最小值等全局诊断量必须使用正确的 `reduction` 或先线程局部汇总再合并。
- 极值及其位置不能只对数值做 `max/min reduction` 后随意更新位置；位置需要二次扫描、临界区或明确的 tie-breaking 规则。
- 如果同一个 `parallel` 区域里连续使用多个 `omp do`，要确认依赖阶段之间没有错误使用 `nowait`。
- `Bx_prev`、`By_prev`、`Bz_prev`、旧温度/旧速度等历史量的更新顺序必须和串行逻辑一致，不能在某些线程仍需要旧值时提前覆盖。

OpenMP 结果和 CPU 基准允许有很小的浮点 reduction 顺序差异，但不应出现系统性的 Nu 符号变化、壁面热流方向变化、稳态残差停不下来或不同线程数给出明显不同的物理量。

### OpenACC 检查

OpenACC 版本重点检查 device 数据生命周期、`present` 关系和 host/device 同步：

- 进入主时间步前，应明确哪些数组常驻 device，例如 `f`、`g`、`rho`、`u/v/w`、`T`、`f_post/g_post`、`Bx_prev/By_prev/Bz_prev` 等。
- GPU kernel 中使用的数组应在数据区内用 `present(...)` 或等价方式确认已在 device 上，避免每步隐式拷贝或意外 host fallback。
- 复杂碰撞 kernel 中的矩空间临时数组和 source term 必须是线程私有，不能被所有 gang/vector 共享。
- `Nu`、`Re`、`errorU`、`errorT`、体平均热流等 scalar 诊断量需要显式 `reduction`，输出到 host 前需要 `copyout` 或 `update self`。
- 如果后处理在 CPU 上做，必须先把需要的 `T/u/v/w/rho` 或相关诊断 buffer 从 device 更新回 host。
- 如果 reload、边界参数、初始条件或宏观场在 host 上被修改，继续进入 GPU kernel 前必须把对应数据更新到 device。
- 使用 `async` 时，要按真实依赖加 `wait`：collision 完成后才能 streaming，streaming 完成后才能 boundary，boundary 完成后才能 macro，macro 完成后才能诊断。
- 文件输出、restart 写出、Tecplot 输出和 CPU 端最终诊断前，应有明确的 host 同步点；不要把普通 Fortran 文件 I/O 放进 GPU kernel。
- `EnableUseG` 分支中的历史通量 `B(t)-B(t-dt)` 更需要检查 device 上的旧值更新顺序，避免 host 和 device 各保存一份不同步的历史场。

### 数组布局与访问模式检查

OpenMP 和 OpenACC 对数组的要求不完全一样。OpenMP 主要受 CPU cache、线程私有变量、内存带宽和 false sharing 影响；OpenACC 主要受 GPU 全局内存带宽、coalesced access、kernel launch 成本、device 常驻数据和线程私有资源影响。因此同一个数组在 OpenMP 中能跑，不代表它适合直接搬到 OpenACC 的大规模寻址模式。

Fortran 是列主序，第一维连续。若数组类似：

```text
f(0:nx+1, 0:ny+1, 0:q-1)
T(0:nx+1, 0:ny+1)
u(0:nx+1, 0:ny+1)
```

则 `i` 或 `x` 方向是连续内存方向。检查循环顺序时应优先确认最内层循环是否沿第一维走，或者 GPU 的 vector 维是否主要覆盖连续的第一维。否则 OpenMP 会更容易 cache miss，OpenACC 会更容易出现非合并访存。

数组可按用途分成几类检查：

- 主场大数组：`f/g/rho/u/v/w/T/Bx_prev/By_prev/Bz_prev` 等规则网格数组适合大规模并行，但应常驻 device，避免每步 host-device 拷贝。
- 碰撞后或迁移缓冲数组：`f_post/g_post` 等如果存在，应确认读写分离，适合 GPU 常驻；若 streaming 原地写，要特别查数据竞争和覆盖顺序。
- 小型常量数组：`ex/ey/ez/opp/omega/N/Ninv/S` 等适合只初始化或 copyin 一次，不应在每个 kernel 里反复传输。
- 局部 scratch 数组：`m/meq/m_post/source/n/neq/n_post` 等在 OpenMP 中通常作为线程私有小数组可以接受；在 OpenACC 中每个 GPU 线程私有一份可能消耗大量 register/local memory，必要时要检查是否应拆成标量、缩小作用域或避免在深层循环里创建过大的私有数组。
- 边界和角点数组：壁面平均、角点拟合、极值插值等通常数据量小、分支多、算术强度低，不一定适合单独开 GPU kernel；如果只在输出时做，CPU 后处理加一次 `update self` 往往更清楚。
- 诊断汇总数组：`Nu`、`Re`、局部 wall Nu、centerline Nu、error 数组或临时平均数组要区分“参与全局 reduction 的量”和“最终只给 host 输出的量”，后者常需要 `copyout` 或 `update self`，不能只写 `present`。
- 间接索引或非规则表：如果某些后处理依赖索引列表、插值表、稀疏访问或复杂分支，它们在 OpenMP 中可能只是慢一点，在 OpenACC 中可能导致访存发散和 kernel 低效，应考虑留在 CPU 或改成规则化的独立 kernel。

LBM 分布函数的方向维也要单独检查。如果数组采用 `f(i,j,k)` 或 `f(i,j,k,q)`，每个方向平面在空间上连续，适合按方向或空间块扫；但单个格点一次读取所有方向时，方向之间可能跨较大 stride。因为 D2Q9/D3Q19/D2Q5/D3Q7 的方向数很小，这通常可以接受，但后续若出现 OpenACC 性能或访存问题，应把“方向维放在最后还是最前”作为性能检查项，而不是在没有 benchmark 的情况下贸然重排数组。

判断一个数组是否适合 GPU 大规模寻址时，可以用下面几个问题：

```text
每个线程访问的位置是否规则、相邻线程是否访问相邻地址？
这个数组是否在每个时间步都被大量访问？
它是否能在整个主循环中常驻 device？
它是否只在边界/输出/最终统计时使用？
它是只读常量、读写主场、还是线程私有 scratch？
它是否需要 reduction、atomic、critical 或二次扫描？
```

满足“规则网格、大量访问、可常驻、读写责任清楚”的数组，通常适合 OpenACC 大规模并行；满足“小规模、分支多、只在输出时使用、需要复杂位置回填”的数组，不一定值得 GPU 化，先保持 CPU 后处理反而更可靠。

OpenACC 语法检查只能说明预处理和 Fortran 语法大体可通过，不能证明 GPU runtime 正确。真正的 OpenACC 检查需要和 OpenMP/CPU 基准在同一宏组合、同一网格、同一时间步数下比较 `errorU/errorT`、`Nu/Re`、壁面 Nu、中心线 Nu、关键极值和最终场。

### OpenMP 与 OpenACC 对照

后续检查 `2DRBOpenmp.F90` 和 `2DRBOpenacc.F90` 时，建议先固定同一组宏：

```text
case macro: SideHeatedCell or RayleighBenardCell
thermal branch: EnableLegacyThermalScheme or EnableUseG
flow mode: steadyFlow or unsteadyFlow
grid / Ra / Pr / Mach / output cadence
```

然后分三层对照：

1. 算法层：碰撞、迁移、边界、宏观量恢复、后处理公式是否同源。
2. 并行层：私有变量、reduction、同步、device 数据传输是否正确。
3. 输出层：`errorU/errorT`、`NuVolAvg/ReVolAvg`、壁面 Nu、中线 Nu、reload 时间和文件输出 cadence 是否表达同一个物理量。

若 OpenMP 和 OpenACC 结果不一致，优先检查并行层：OpenACC 的 `private/present/reduction/update/wait` 和 OpenMP 的 `private/reduction/nowait`。只有排除并行语义问题后，再回头怀疑算法公式或物理参数。

## 下一步代码检查清单

后续从 `2DRBOpenmp.F90` 和 `2DRBOpenacc.F90` 开始，建议按这个顺序查：

1. 顶部宏：`EnableLegacyThermalScheme`、`EnableUseG`、`steadyFlow`、`unsteadyFlow`、算例宏是否组合合理。
2. 参数：`Ra`、`Pr`、`Mach`、`tauf`、`diffusivity`、`gBeta`、`Tref`、`lengthUnit` 是否符合目标文献。
3. 流场：D2Q9/D3Q19 速度集、MRT 矩、平衡矩、松弛参数、Guo forcing、半步力速度恢复。
4. 温度场：按当前分支检查 D2Q5/D3Q7 的矩变换、平衡矩、松弛参数、`g_post` 逆变换。
5. `EnableUseG`: 检查 `Bx_prev`、`By_prev`、`Bz_prev` 的初始化、碰撞更新、重启后恢复，以及 `dB` 是否只在该分支生效。
6. streaming：均匀网格应为普通邻点迁移，不应混入非均匀网格插值。
7. 边界：速度 bounce-back、温度 anti-bounce-back / bounce-back 与算例几何一致。
8. 宏观量：`rho`、`u/v/w`、`T` 的恢复顺序是否支持下一步碰撞和力项计算。
9. 诊断：Nu、Re、壁面热流、体平均热流、中心线热流是否和当前算例方向及尺度一致。
10. 稳态 side-heated / RB 后处理：`Nu`、壁面 Nu、中心线 Nu、热壁局部极值、流函数和中心线速度是否按同一套稳态热流诊断理解；RB 只把热流方向从 `x/u` 换为 `y/v`。
11. OpenMP：循环私有变量、reduction、阶段间同步和历史场更新是否保持串行算法含义。
12. OpenACC：device 数据生命周期、`present`、`private`、`reduction`、`update self/device` 和 `async/wait` 是否完整。
13. OpenMP/OpenACC 对照：同一宏组合下比较 `errorU/errorT`、`Nu/Re`、壁面 Nu、中线 Nu 和最终场。
14. 数组布局：主场大数组、常量表、局部 scratch、边界后处理数组是否分别适合当前 OpenMP/OpenACC 访问模式。
15. 稳态/非稳态：停止准则、采样间隔、输出间隔、重启间隔和最终输出是否各自独立且含义清楚。

## 需要特别确认的点

1. 如果目标是严格复现 Xu/ISLBM 系列 legacy 温度算法，应关闭 `EnableUseG`。否则历史通量修正会引入 Chai-Shi NCDE 算法成分。

2. 如果目标是 `EnableUseG` 分支，则 legacy 中的 `aT`、`Qk = 3 - sqrt(3)`、`Qe = Qnu = 4*sqrt(3) - 6` 不应混入当前温度扩散系数定义。

3. `EnableUseG` 的温度边界系数必须跟随固定权重改变：

```text
D2Q5 UseG: ABB coefficient = 1/3
D3Q7 UseG: ABB coefficient = 1/4
```

不能继续使用 legacy 的 `(4+aT)/10` 或 `(6+aT)/21`。

4. Chai-Shi 中的 `G_i` 对应 `partial_t B`。如果代码变量 `dB` 存的是差值 `B(t)-B(t-dt)`，则修正项可写成 `(1-S1/2)*dB`；如果存的是时间导数，则必须乘回 `dt`。

5. 3D `EnableUseG` 是从 Chai-Shi 一般 CDE 框架拓展到 D3Q7 的实现，不是文献逐项列出的专门 3D 热对流算例。必须重点检查三维新增的 `wT`、`dBz`、D3Q7 矩阵、`N^{-1}` 和 `opp` 方向编号。

6. 当前使用均匀网格时，ISLBM 的二次插值和非均匀平均公式只作为背景文献，不作为代码应具备的功能。但 Nu 壁面导数仍应采用与 half-way 边界一致的均匀网格壁面导数公式。

7. 网格坐标应按 half-way 边界理解：物理壁面在 `xp(0)`、`xp(nx+1)` 等边界坐标上，第一层流体节点 `xp(1)`、`xp(nx)` 分别距离壁面半个格距。`lengthUnit` 应对应物理壁面间距，而不是首末流体节点间距。

8. 流场是低 Mach 弱可压缩 LBM，而不是严格 projection incompressible solver。当前实现使用全密度 `rho`，不是密度波动变量。需要检查 `Ma`、密度涨落诊断、速度散度、半步力修正和浮力项符号。

9. `Tref` 主要进入浮力项 `T-Tref` 和无量纲温度定义；Nu 的尺度应使用 `DeltaT = Thot - Tcold`，不要把 `Tref` 误当成温差。

10. 2D 热驱动方腔文献只给出稳态后处理；它不能为非稳态时间平均或统计稳态窗口提供依据。稳态 RB 当前只是沿用同一套稳态热流诊断并旋转热流方向。

11. 该文献 Eq. (44) 的左热壁热流符号应作为疑似笔误处理。与 `q_x = u*theta - dtheta/dx` 一致的正热壁 Nu 应使用 `2*(theta_h-theta(1,j))/dx`。

12. 并行实现的第一怀疑对象通常不是物理公式，而是作用域、reduction、同步和数据传输。特别是 OpenACC 结果异常时，应先查 `private`、`present`、`reduction`、`update self/device` 和 `wait`。

13. 任何 OpenMP/OpenACC 对照都必须使用同一套宏和参数。不要拿 `SideHeatedCell` 的后处理去比较 `RayleighBenardCell`，也不要把 `EnableLegacyThermalScheme` 和 `EnableUseG` 的温度场结果直接混比。

14. 数组是否适合 OpenACC 不能只看“能不能放到 device 上”。规则大数组适合 GPU 常驻和大规模寻址；小边界数组、复杂插值数组、极值位置回填和一次性输出数组不一定适合 GPU 化。OpenMP 和 OpenACC 可以保留不同的数组处理策略，只要最终算法含义一致。

## 最简检查表

```text
[ ] Legacy D2Q5: Q = diag(0,Qk,Qk,Qe,Qnu)
[ ] Legacy D3Q7: Q = diag(0,Qk,Qk,Qk,Qe,Qnu,Qnu)
[ ] Legacy: Qk=3-sqrt(3), Qe=Qnu=4sqrt(3)-6
[ ] Legacy D2Q5 ABB coefficient = (4+aT)/10
[ ] Legacy D3Q7 ABB coefficient = (6+aT)/21

[ ] UseG D2Q5: omega0=1/3, omega_axis=1/6, csT^2=1/3
[ ] UseG D3Q7: omega0=1/4, omega_axis=1/8, csT^2=1/4
[ ] UseG: kappa=csT^2*(1/Qk-1/2)
[ ] UseG: Qk=1/(1/2+kappa/csT^2)
[ ] UseG D2Q5 ABB coefficient = 1/3
[ ] UseG D3Q7 ABB coefficient = 1/4
[ ] UseG: dB is clearly defined as difference or derivative
[ ] UseG 3D: dBz is initialized, updated, saved, and reloaded consistently

[ ] Flow field remains same between two temperature branches
[ ] Flow field uses full rho, not density-fluctuation as the primary variable
[ ] Guo forcing uses half-step velocity correction
[ ] buoyancy force uses T-Tref, not DeltaT directly
[ ] Nu uses DeltaT and correct wall-normal direction
[ ] Uniform grid does not include ISLBM interpolation
[ ] Grid coordinates use half-way boundaries: first fluid node is half a lattice spacing from the wall
[ ] Wall derivative for Nu respects half-way boundary position
[ ] Side-heated postprocessing follows the steady-only 2D square-cavity definitions
[ ] Steady RB postprocessing uses the same steady heat-flux diagnostics with x/u replaced by y/v
[ ] Hot-wall Nu treats the printed Eq. (44) sign as a likely typo
[ ] Nu, wall Nu, and middle-line Nu are compared only as steady consistency diagnostics
[ ] Unsteady postprocessing is not inferred from the 2D square-cavity postprocessing paper

[ ] OpenMP loop temporaries are private
[ ] OpenMP global diagnostics use correct reductions
[ ] OpenMP stages do not use nowait across real dependencies
[ ] OpenACC resident arrays have a clear enter/update/exit data lifetime
[ ] OpenACC kernels use present/private/reduction clauses intentionally
[ ] OpenACC host-side output and CPU diagnostics are preceded by update self
[ ] OpenACC async kernels have waits at real collision/streaming/boundary/macro dependencies
[ ] OpenMP and OpenACC are compared under the same macro, grid, Ra, Pr, and output cadence
[ ] Main field arrays use access patterns suitable for contiguous/cache/coalesced memory access
[ ] Small constants are initialized or copied once, not repeatedly transferred per kernel
[ ] Large scratch arrays are not accidentally shared in OpenMP or over-privatized in OpenACC
[ ] Boundary-only and output-only arrays are evaluated before forcing them onto GPU
```
