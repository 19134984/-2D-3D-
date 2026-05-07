# 热对流 LBM 算法说明文档修正建议

本文档用于修正原始 md 文档：

```text
热对流 LBM 算法文献与代码分支对应说明
```

目标不是推翻原文，而是在原文基础上补充和纠正若干容易影响后续代码检查的细节。原文整体分支逻辑是合理的：

- `EnableLegacyThermalScheme` 对应 Xu 系列热对流 DDF-MRT 算法，包括均匀网格退化形式和 ISLBM 非均匀网格拓展；
- `EnableUseG` 对应 Chai-Shi MRT-LB for NCDE 框架在温度对流扩散方程上的使用；
- `EnableUseG` 分支只替换温度场算法，流场仍沿用 Xu 系列 D2Q9/D3Q19 MRT + Guo forcing 思路；
- 当前代码没有启用非均匀网格，因此不需要引入 ISLBM 的二次插值 streaming，但仍要保留 half-way 边界下的壁面导数处理。

下面给出需要修正和补充的内容。建议 Codex 按照本文件逐项修改原 md。

---

## 1. 文献与代码分支对应关系：保留原逻辑，但补充“不是完全同一算法”的说明

原文的分支对应关系可以保留。需要补充的是：`EnableUseG` 并不是对前两篇文献的简单复现，而是“流场沿用前两篇，温度场替换为第三篇 Chai-Shi 框架”。

建议在“文献与分支”表格后增加：

```markdown
需要注意的是，`EnableUseG` 分支不是严格复现 Xu/ISLBM 系列文献中的温度 MRT 模型，而是只沿用其流场 D2Q9/D3Q19 MRT、Guo forcing、速度边界和热对流物理模型；温度场则替换为 Chai-Shi NCDE MRT 框架中的辅助分布函数 `G_i` 修正模型。

因此后续检查代码时不能把 `EnableLegacyThermalScheme` 和 `EnableUseG` 的温度参数混用：

- legacy 分支使用 `aT`、`Qk=3-sqrt(3)`、`Qe=Qnu=4sqrt(3)-6`；
- UseG 分支使用固定温度权重、`csT^2`、`taug=1/2+kappa/csT^2` 和历史通量修正 `dB`。
```

---

## 2. 流场算法表述：把“不可压缩形式”改成“低 Mach 弱可压缩形式”

原文写法：

```markdown
平衡分布或平衡矩采用低 Mach 数不可压缩形式。
```

这个说法偏强。Xu 系列文献中的 D3Q19 流场模型本质上仍是低 Mach 弱可压缩 LBM；它在小 Mach 数条件下近似恢复不可压缩 Boussinesq 方程。文献中还专门检查了 density fluctuation 和 velocity divergence，这说明算法本身并不是严格投影法意义上的 incompressible solver。

建议替换为：

```markdown
流场采用低 Mach 数弱可压缩 MRT-LBM。在 `Ma` 足够小的条件下，该模型近似恢复不可压缩 Boussinesq 热对流方程。

如果代码采用 density-fluctuation equilibrium 或 incompressible equilibrium，则需要单独注明；否则默认应按弱可压缩 D2Q9/D3Q19 MRT-LBM 理解。
```

同时建议在流场部分增加：

```markdown
流场压力通常由 LBM 状态方程近似给出，例如 `p = rho * cs^2`。因此密度涨落虽然很小，但并不严格为零。检查代码时需要确认：

- `rho = sum_i f_i` 是否参与平衡矩；
- 速度恢复是否包含半步力修正；
- 是否存在专门的 density-fluctuation 写法。
```

---

## 3. 流场守恒模松弛率：补充文献写法与代码写法的区别

原文写了：

```markdown
nu = cs^2 * (tau_f - 1/2) * dt
Snu = 1 / tau_f
Sq, Sm = 8 * (2 tau_f - 1) / (8 tau_f - 1)
```

这几条是正确的，但容易漏掉守恒模。Xu 系列 D3Q19 MRT 文献中，守恒矩对应的松弛率通常写为：

```text
s_rho = 0
s_j   = 0
```

而非守恒二阶动量通量相关模写为：

```text
s_e = s_epsilon = s_nu = s_pi = 1 / tau_f
s_q = s_m = 8 * (2*tau_f - 1) / (8*tau_f - 1)
```

建议在流场松弛参数段落后加入：

```markdown
注意守恒模的处理：

在 Xu 系列 MRT 文献中，密度和动量守恒模通常取

```text
s_rho = 0
s_j   = 0
```

这表示碰撞不改变守恒矩。若代码中守恒模取 `1` 或其他数值，只要碰撞实现中守恒矩的 post-collision 值被强制保持为 equilibrium/conserved value，宏观守恒性通常仍可成立；但这已经不是文献逐项复现的写法。

因此查代码时要区分两件事：

1. 守恒矩是否真的被碰撞改变；
2. 松弛矩阵中守恒模的数值是否和文献一致。
```

---

## 4. `EnableLegacyThermalScheme` 的 D2Q5 温度松弛率：补上 `Qe`

原文 D2Q5 legacy 部分写了：

```markdown
Qk  = 3 - sqrt(3)
Qnu = 4 * sqrt(3) - 6
```

应补充 `Qe`。D2Q5 温度矩通常可以写成：

```text
n = [T, uT, vT, aT*T, 0]^T
```

其中第 4 个矩 `aT*T` 对应的松弛率也需要明确。建议替换为：

```markdown
D2Q5 温度松弛矩阵建议按下面形式检查：

```text
Q = diag(0, Qk, Qk, Qe, Qnu)
```

其中

```text
Qk  = 3 - sqrt(3)
Qe  = 4 * sqrt(3) - 6
Qnu = 4 * sqrt(3) - 6
```

这里：

- 第 0 个矩 `T` 是守恒温度矩，不参与松弛；
- 第 1、2 个矩 `uT`、`vT` 是温度通量矩，其松弛率 `Qk` 决定热扩散率；
- 第 3 个矩 `aT*T` 的松弛率应记为 `Qe`；
- 第 4 个高阶矩的松弛率记为 `Qnu`。
```

---

## 5. `EnableLegacyThermalScheme` 的 D3Q7 温度松弛率：同样补上 `Qe`

原文 D3Q7 legacy 部分写了：

```markdown
Qk  = 3 - sqrt(3)
Qnu = 4 * sqrt(3) - 6
```

建议替换为：

```markdown
D3Q7 温度矩为：

```text
n = [T, uT, vT, wT, aT*T, 0, 0]^T
```

对应松弛矩阵应按下面形式检查：

```text
Q = diag(0, Qk, Qk, Qk, Qe, Qnu, Qnu)
```

其中

```text
Qk  = 3 - sqrt(3)
Qe  = 4 * sqrt(3) - 6
Qnu = 4 * sqrt(3) - 6
```

这里：

- `T` 是守恒温度矩；
- `uT`、`vT`、`wT` 是三个方向的温度通量矩；
- `aT*T` 对应 `Qe`；
- 最后两个二阶各向异性矩对应 `Qnu`。
```

这一点很重要，因为后续检查代码时，如果只查 `Qk` 和 `Qnu`，很容易漏掉 `Qe` 是否正确赋值。

---

## 6. legacy 温度边界：保留原式，但注明它只适用于 legacy 权重

原文 legacy 恒温边界写法是正确的：

```text
D2Q5: g_opp = -g_i^+ + (4 + aT)/10 * T_wall
D3Q7: g_opp = -g_i^+ + (6 + aT)/21 * T_wall
```

但需要强调：这两个系数来自 legacy 模型的温度权重，不能用于 `EnableUseG` 分支。

建议在 legacy 边界公式后加入：

```markdown
上述系数来自 half-way anti-bounce-back 的统一形式：

```text
g_opp = -g_i^+ + 2 * omega_axis * T_wall
```

对于 legacy 温度权重：

D2Q5:

```text
omega_axis = (4 + aT) / 20
2 * omega_axis = (4 + aT) / 10
```

D3Q7:

```text
omega_axis = (6 + aT) / 42
2 * omega_axis = (6 + aT) / 21
```

因此 `(4+aT)/10` 和 `(6+aT)/21` 只适用于 `EnableLegacyThermalScheme`，不能直接用于 `EnableUseG`。
```

---

## 7. `EnableUseG` 的恒温边界系数必须单独写出

这是原 md 中最需要补充的地方。`EnableUseG` 分支不再使用 `aT` 权重，而是采用固定权重：

D2Q5:

```text
omega_0 = 1/3
omega_1..4 = 1/6
csT^2 = 1/3
```

D3Q7:

```text
omega_0 = 1/4
omega_1..6 = 1/8
csT^2 = 1/4
```

因此恒温 anti-bounce-back 的系数也要改变。

建议在 `EnableUseG` 温度场部分增加一个小节：

```markdown
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

- `EnableLegacyThermalScheme` 使用 `(4+aT)/10` 或 `(6+aT)/21`；
- `EnableUseG` 使用 `1/3` 或 `1/4`；
- 两个分支不能共用同一个 `ABB` 系数，除非代码内部已经按 `2*omega_i` 自动计算。
```

---

## 8. `EnableUseG` 中 `G_i` 或历史通量修正项：明确 `dB` 是差值还是时间导数

原文写法：

```markdown
B      = u * T
dB     = B(t) - B(t - dt)
SG     = 1 - 0.5 * Qk
n_flux^post = n_flux - Qk * (n_flux - B) + SG * dB
```

这个写法在 `dt=1` 且 `dB` 表示差值时是正确的。但 Chai-Shi 文献中的 `G_i` 是通过 `partial_t B` 定义的，而不是单纯 `B(t)-B(t-dt)`。因此必须在文档中说清楚：

```markdown
Chai-Shi 框架中的辅助分布函数为：

```text
G_i = omega_i * e_i · [(I - S1/2) partial_t B] / csT^2
```

其中

```text
B = u*T
```

若代码中定义：

```text
dB = B(t) - B(t-dt)
```

则因为演化方程中实际加入的是 `dt * G_i`，在 lattice unit `dt=1` 下，moment-space 的一阶通量修正可以写成：

```text
n_flux^post = n_flux - Qk * (n_flux - B) + (1 - 0.5*Qk) * dB
```

若代码中定义：

```text
dBdt = [B(t) - B(t-dt)] / dt
```

则修正项必须写成：

```text
n_flux^post = n_flux - Qk * (n_flux - B) + dt * (1 - 0.5*Qk) * dBdt
```

因此检查代码时要先确认变量名 `dB`、`dBx`、`dBy`、`dBz` 存的是“差值”还是“时间导数”。不能同时除以 `dt` 又省略前面的 `dt`，否则会多乘或少乘一个时间步因子。
```

---

## 9. `EnableUseG` 的 D2Q5 温度模型：保留原式，但标成当前实现选择

原文写法：

```text
omega_0 = 1 / 3
omega_1..4 = 1 / 6
csT^2 = 1 / 3
taug = 1/2 + kappa / csT^2
Qk = 1 / taug
Qnu = 1

n_eq = [T, uT, vT, -2T/3, 0]^T
```

这个结构是自洽的。建议补充解释：

```markdown
这里的 `Qnu=1` 不是由热扩散率唯一决定的，而是当前实现中对非守恒高阶矩的快速松弛选择。Chai-Shi NCDE 框架主要约束一阶通量矩的松弛率 `S1`，从而确定热扩散率：

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

D2Q5 UseG 分支若采用

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
```

---

## 10. `EnableUseG` 的 D3Q7 温度模型：明确这是三维实现性拓展

原文写法：

```text
omega_0 = 1 / 4
omega_1..6 = 1 / 8
csT^2 = 1 / 4
taug = 1/2 + kappa / csT^2
Qk = 1 / taug
Qnu = 1

n_eq = [T, uT, vT, wT, -3T/4, 0, 0]^T
```

这个结构是合理的，但必须明确：这是根据 Chai-Shi 一般 `DdQ(2d+1)` CDE-MRT 框架做的三维 D3Q7 拓展，不是文献中逐项给出的 3D 热对流算例。

建议将该段改为：

```markdown
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

- `wT` 是否位于 D3Q7 温度矩的第 3 个通量矩；
- `dBz` 是否加到 `wT` 对应的矩上；
- `ex/ey/ez` 的编号和 `N`、`N^{-1}`、`opp` 方向编号完全一致；
- 重启文件中是否保存和恢复了 `Bx_prev`、`By_prev`、`Bz_prev`。
```

---

## 11. D3Q7 矩阵和速度编号：增加专门检查项

三维温度场最容易出错的地方不是公式，而是编号和符号。建议在原 md 的“需要特别确认的点”中加入：

```markdown
- D3Q7 的速度编号、反方向编号、矩阵 `N` 和逆矩阵 `N^{-1}` 必须逐项一致。特别是第三个方向 `ez` 的正负号、`wT` 所在矩、以及 `dBz` 的符号要重点核对。只要 `ez` 行符号和速度编号不一致，三维温度通量就会出现方向性错误。
```

如果原代码采用如下速度编号：

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

---

## 12. 均匀网格退化：保留“不需要 ISLBM 插值”，但补充 Nu 壁面导数仍需 half-way 修正

原文写：

```markdown
后续检查代码时不需要寻找 ISLBM 的二次插值、非均匀坐标、非均匀加权平均或非均匀壁面导数。
```

建议改成：

```markdown
当前代码采用均匀网格时，不需要寻找 ISLBM 的二次插值、非均匀坐标、非均匀加权平均或非均匀壁面导数。

但这不代表壁面热流可以随意用普通中心差分。由于速度和温度边界采用 half-way bounce-back / anti-bounce-back，第一层流体节点距离物理壁面是半个格距。因此计算壁面温度梯度和 Nusselt 数时，仍应使用与 half-way 边界一致的均匀网格壁面导数公式。
```

例如 bottom/hot wall 的温度梯度不能简单写成：

```text
(T(2)-T(1))/dx
```

而应根据壁面位置和第一、第二层流体节点位置构造 one-sided 二阶公式。若代码中采用文献给出的均匀网格 half-way 公式，需要确认：

```text
dT/dn_wall = (-8*Twall + 9*T1 - T2) / (3*dx)
```

或对应方向上的符号版本。这里 `T1` 是离壁面最近的第一层流体节点，`T2` 是第二层流体节点。

---

## 13. `Tref` 的表述：避免把参考温度误写进 Nu

原文写：

```markdown
Tref 只影响浮力项和热流后处理中的参考温度，不应和壁温差 DeltaT 混淆。
```

建议改成：

```markdown
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

- 温度梯度方向是否与冷热壁方向一致；
- 分母是否使用 `DeltaT`；
- `Tref` 是否只进入浮力项或无量纲温度定义，而不是被误用为壁温差。
```

---

## 14. 建议重写原 md 中的“边界条件”小节

原边界条件部分比较简略，建议替换为下面版本：

```markdown
## 边界条件

速度边界：

- 固壁采用 half-way bounce-back，对应无滑移边界；
- 对 D2Q9/D3Q19，需确认 `opp` 方向编号和 streaming 顺序一致。

温度边界：

- 恒温边界采用 half-way anti-bounce-back；
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

边界检查不能只看分布函数公式，还要同时检查：

- 边界后的 `T` 是否由 `sum_i g_i` 恢复；
- `opp` 方向编号是否正确；
- 热流和 Nusselt 数计算是否使用同一个壁面方向；
- `Thot`、`Tcold`、`Tref` 和 `DeltaT` 是否被一致使用；
- `EnableLegacyThermalScheme` 和 `EnableUseG` 是否使用了各自正确的 ABB 系数。
```

---

## 15. 建议重写原 md 中的“需要特别确认的点”

建议把原来的“需要特别确认的点”替换或扩展为：

```markdown
## 需要特别确认的点

1. 如果目标是严格复现 Xu/ISLBM 系列 legacy 温度算法，应关闭 `EnableUseG`。否则历史通量修正会引入 Chai-Shi NCDE 算法成分。

2. 如果目标是 `EnableUseG` 分支，则 legacy 中的 `aT`、`Qk=3-sqrt(3)`、`Qe=Qnu=4sqrt(3)-6` 不应混入当前温度扩散系数定义。

3. `EnableUseG` 的温度边界系数必须跟随固定权重改变：
   - D2Q5 UseG: ABB 系数为 `1/3`；
   - D3Q7 UseG: ABB 系数为 `1/4`。
   不能继续使用 legacy 的 `(4+aT)/10` 或 `(6+aT)/21`。

4. Chai-Shi 中的 `G_i` 对应 `partial_t B`。如果代码变量 `dB` 存的是差值 `B(t)-B(t-dt)`，则修正项可写成 `(1-S1/2)*dB`；如果存的是时间导数，则必须乘回 `dt`。

5. 3D `EnableUseG` 是从 Chai-Shi 一般 CDE 框架拓展到 D3Q7 的实现，不是文献逐项列出的专门 3D 热对流算例。必须重点检查三维新增的 `wT`、`dBz`、D3Q7 矩阵、`N^{-1}` 和 `opp` 方向编号。

6. 当前使用均匀网格时，ISLBM 的二次插值和非均匀平均公式只作为背景文献，不作为代码应具备的功能。但 Nu 壁面导数仍应采用与 half-way 边界一致的均匀网格壁面导数公式。

7. 流场是低 Mach 弱可压缩 LBM，而不是严格 projection incompressible solver。需要检查 `Ma`、密度涨落、速度散度、半步力修正和浮力项符号。

8. `Tref` 主要进入浮力项 `T-Tref` 和无量纲温度定义；Nu 的尺度应使用 `DeltaT=Thot-Tcold`，不要把 `Tref` 误当成温差。
```

---

## 16. 给 Codex 的执行要求

请 Codex 按下面要求修正原 md：

1. 保留原文总体结构，不要重写成完全不同的文档。
2. 优先修改以下段落：
   - 文献与分支；
   - 流场算法；
   - `EnableLegacyThermalScheme` 温度场；
   - `EnableUseG` 温度场；
   - 边界条件；
   - 均匀网格退化；
   - 需要特别确认的点。
3. 所有公式都保持 plain Markdown + code block 风格，便于后续继续复制到 Fortran 注释或代码检查说明中。
4. 不要把 `EnableLegacyThermalScheme` 和 `EnableUseG` 的温度参数混用。
5. 对所有新增内容保持“检查基准”的语气，不要声称代码已经正确。
6. 最终输出一份修正后的完整 md 文档，而不是只输出 diff。

---

## 17. 最简检查表

后续检查代码时，可按下面 checklist 快速判断：

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
[ ] Guo forcing uses half-step velocity correction
[ ] buoyancy force uses T-Tref, not DeltaT directly
[ ] Nu uses DeltaT and correct wall-normal direction
[ ] Uniform grid does not include ISLBM interpolation
[ ] Wall derivative for Nu respects half-way boundary position
```
