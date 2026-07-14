# Ginzburg–d’Humières 2003 证据账本

## 论文身份与证据边界

- 作者：Irina Ginzburg, Dominique d’Humières。
- 标题：*Multireflection boundary conditions for lattice Boltzmann models*。
- 载体：*Physical Review E* 68 (2003), 066614。
- 本地 PDF：`pdf/Multireflection boundary conditions for lattice Boltzmann models.pdf`。
- 页码口径：PDF p. `n` 对应文章页 `066614-n`。
- 证据边界：以本地 PDF 黑色印刷正文为证据；`ginzburg2003.txt` 只用于定位，p. 3 的黄色高亮等批注层不作为原文证据。

## 原符号与统一符号

| 原符号 | 原文含义 | 统一符号 | 精确映射或备注 |
| --- | --- | --- | --- |
| `λ_ν` | 中心反演下不变的偶/黏性通道碰撞率 | `r_even` | `r_even=λ_ν`，`H_even=1/r_even-1/2` |
| `λ_2` | 其他模态的奇/第二通道碰撞率 | `r_odd` | `r_odd=λ_2`，`H_odd=1/r_odd-1/2` |
| `Λ²` | Eq. (41) 定义的归一化偶奇 Hénon 乘积 | `Λ²` | `Λ²=(4/3)H_even H_odd` |
| `H_{qeff}`, `H_q` | 沿链接 `c_q` 的有效通道宽度与规定宽度 | `W_{q,eff}`, `W_q` | 仅为避免与本教程 Hénon 参数 `H` 混淆而改名；二者都是格子长度 |
| `D_eff`, `D_{1/2}` | 有效通道直径；以截断链接中点计算的直径 | 保持不变 | 格子长度，不是扩散系数 |
| `I_f` | Eq. (9) 动量定义 `J=j-I_f F` 中的 forcing gauge | `I_f` | 原文标准定义取 `0`；本文使用的修正定义取 `-1/2` |
| `δ_q` | 边界沿链接的无量线性位置 | 保持不变 | halfway bounce-back 中 `δ_q=1/2` |
| `Θ_q` | `c_q` 在墙面法向上的投影（格子单位） | 保持不变 | 主坐标轴强迫下的相关链接取 `Θ_q²=1` |
| `κ_{±1}`, `\barκ_{−1}`, `\barκ_{−2}` | 边界闭合系数 | 保持原记号 | 不是温度扩散率 `κ` |

所有 `H=1/r-1/2` 换算均是定义和精确代数，要求相应 `r_even≠0` 与 `r_odd≠0`。本账本不从这一定义自动假定 `0<r<2`。

## 核心证据

### 1. 两率简化模型、forcing gauge 与线性 Stokes 语境

- **原文位置：** PDF p. 3，文章页 066614-3，Eqs. (6), (9)–(12) 及 §II B。
- **原符号：** 简化模型只使用 `λ_ν` 与 `λ_2` 两个碰撞特征值；Eq. (9) 定义

  \[
  \mathbf J=\mathbf j-I_f\mathbf F,
  \]

  并明说标准定义为 `I_f=0`，本文使用的修正定义为 `I_f=-1/2`。线性平衡 Eq. (11) 导向稳态不可压 Stokes Eq. (12)。
- **统一符号：** `λ_ν→r_even`，`λ_2→r_odd`；`I_f` 保持不变。
- **允许主张：** `I_f=-1/2` 是该文实际采用的动量/forcing gauge；后续 Eq. (42) 的简化明确依赖此值。
- **不支持 / 不可外推：** `I_f=-1/2` 不是所有 forcing 实现共有的普适定义；本页的 Stokes 结论使用线性平衡且限于稳态不可压语境，不支持直接换成非线性或非稳态 Navier–Stokes 结论。

### 2. 稳态平行平板 Poiseuille 剖面（Eqs. (38)–(39)）

- **原文位置：** PDF p. 6，文章页 066614-6，§IV D，Eqs. (38)–(39) 及紧邻定义。
- **原文前提：** 两个静止、平行平面墙关于原点对称，均匀力 `\mathbf F` 平行于墙面，所考察的链接方向与坐标是 `c_q,x_q`，流动为稳态抛物线剖面。
- **原公式：**

  \[
  j_q=j_{0q}\left(1-\frac{4x_q^2}{H_{q\mathrm{eff}}^2}\right),
  \tag{38}
  \]

  \[
  F_q=j_{0q}\frac{8\nu}{\Theta_q^2H_{q\mathrm{eff}}^2}.
  \tag{39}
  \]

- **统一符号：** 为避免与 Hénon 参数冲突，将原文几何宽度 `H_{qeff}` 重记为 `W_{q,eff}`；因此 Eq. (38) 的分母是 `W_{q,eff}²`，不是 `H_even` 或 `H_odd`。
- **允许主张：** Eqs. (38)–(39) 定义了本节中有效通道宽度与抛物线幅值、均匀体力和黏性的关系。
- **不支持 / 不可外推：** 该剖面不是任意几何、任意时变 forcing 或任意非线性流的通用解；`H_{qeff}` 也不是碰撞模态的 Hénon 参数。

### 3. 通用边界系数下的有效宽度（Eq. (40)）

- **原文位置：** PDF p. 6，文章页 066614-6，Eq. (40) 及其紧前一句。
- **原文前提：** 取 `F_{\bar q}^{pc}=0`，使用 Eq. (10) 和 Eqs. (35)–(38)，采用线性平衡 Eq. (11) 并舍去非线性项。
- **原公式以统一几何宽度 `W` 重记：**

  \[
  \begin{aligned}
  W_{q,\mathrm{eff}}^2={}&W_q^2+6\Lambda^2-4\delta_q^2
  -\frac{8\nu}{\Theta_q^2}\left(\frac1{\lambda_2}+I_f\right)\\
  &+\frac{4(\kappa_{-1}-\bar\kappa_{-2})}
  {1-\bar\kappa_{-1}-\bar\kappa_{-2}}\\
  &+\frac{12\nu}{1-\bar\kappa_{-1}-\bar\kappa_{-2}}
  \left[1-2\delta_q+(1+2\delta_q)\bar\kappa_{-1}
  +(3+2\delta_q)\bar\kappa_{-2}\right].
  \end{aligned}
  \tag{40}
  \]

  原文左端与第一项分别写作 `H_{qeff}²` 与 `H_q²`。式中所有分母必须非零，特别是 `Θ_q≠0`, `λ_2≠0` 与 `1-\barκ_{-1}-\barκ_{-2}≠0`。
- **允许主张：** Eq. (40) 明示有效宽度同时受几何位置、碰撞参数、forcing gauge 和边界系数影响；不是只由某一个松弛时间决定。
- **不支持 / 不可外推：** 脱离上述线性平衡、稳态 Poiseuille、均匀平行力和本节边界闭合系数时，Eq. (40) 不能作为通用壁面位置公式。

### 4. 归一化偶奇乘积（Eq. (41)）

- **原文位置：** PDF p. 6，文章页 066614-6，Eq. (41)。
- **原公式：**

  \[
  \Lambda^2=\frac43
  \left(\frac1{\lambda_\nu}-\frac12\right)
  \left(\frac1{\lambda_2}-\frac12\right).
  \tag{41}
  \]

- **统一符号：**

  \[
  \Lambda^2=\frac43H_{\mathrm{even}}H_{\mathrm{odd}},
  \quad
  H_{\mathrm{even}}=\frac1{r_{\mathrm{even}}}-\frac12,
  \quad
  H_{\mathrm{odd}}=\frac1{r_{\mathrm{odd}}}-\frac12.
  \]

- **允许主张：** 论文的 `Λ²` 包含因子 `4/3`；因此 `Λ²` 不等于未归一化的 `H_even H_odd`。
- **不支持 / 不可外推：** Eq. (41) 是本文两率简化模型的定义，不能把 `Λ²=1/4` 直接改写为 `H_even H_odd=1/4`；也不能与 D2Q5 温度体相的 `H_κH_ν=1/6` 混同。

### 5. halfway bounce-back 的有效直径（Eq. (42)）

- **原文位置：** PDF p. 7，文章页 066614-7，Eq. (42) 及紧前后段落。
- **原文前提：** bounce-back 系数取 Eq. (23)，`δ_q=1/2`，且 `I_f=-1/2`；它是 Eq. (40) 在这些特殊条件下的简化。
- **原公式：**

  \[
  D_{\mathrm{eff}}^2=D_{1/2}^2+4\Lambda^2-1,
  \tag{42}
  \]

  其中 `D_eff` 是有效直径，`D_{1/2}` 是以截断链接中点定位两面墙时的直径。
- **统一符号：** `D` 不改名；只将 Eq. (41) 中的两个率换成 `r_even,r_odd`。
- **允许主张：** 在上述受限场景中，`Λ²=1/4` 精确使 `D_eff²=D_{1/2}²`；取正直径时得 `D_eff=D_{1/2}`。
- **不支持 / 不可外推：** Eq. (42) 不支持任意方向、曲面、非稳态、非线性流动、其他 forcing gauge 或温度 ABB 也必然有精确 halfway 壁面；它是 D2Q9/D3Q15 流场边界分析，不是 D2Q5 体相四阶系数。

### 6. `Λ²=1/4`、率关系与 `3/16`（Eq. (43)）

- **原文位置：** PDF p. 7，文章页 066614-7，Eq. (43) 及其前后句。
- **原公式与统一符号：** 论文说 `Λ²=1/4` 或等价地

  \[
  \lambda_2(\lambda_\nu)
  =8\frac{2-\lambda_\nu}{8-\lambda_\nu},
  \tag{43}
  \]

  即

  \[
  r_{\mathrm{odd}}
  =8\frac{2-r_{\mathrm{even}}}{8-r_{\mathrm{even}}}.
  \]

  由 Eq. (41) 与 Hénon 参数定义进行精确代数换写，得

  \[
  \Lambda^2=\frac14
  \quad\Longleftrightarrow\quad
  \frac43H_{\mathrm{even}}H_{\mathrm{odd}}=\frac14
  \quad\Longleftrightarrow\quad
  H_{\mathrm{even}}H_{\mathrm{odd}}=\frac{3}{16}.
  \]

  反过来求解碰撞率就得 Eq. (43)。等价链的代数定义域要求 `r_even≠0`, `r_odd≠0`；Eq. (43) 还要求 `r_even≠8`。在这个定义域内，`r_even=2` 会导致 `r_odd=0`，因而不属于可用的 Hénon 换算点。
- **允许主张：** 对 Eq. (42) 的受限 halfway bounce-back Poiseuille 场景，Eq. (43) 与 `Λ²=1/4` 是同一条率约束；原文据此说明墙面位于 `±(x_{qb}+1/2)` 时的抛物线剖面与精确 Poiseuille 解一致。
- **不支持 / 不可外推：** `3/16` 是本边界场景的未归一化偶奇 Hénon 乘积，不是 `Λ²` 本身，不是 D2Q5 温度各向同性的 `1/6`，也不是任意边界、任意 forcing 或任意格子的普适 magic 常数。

### 7. 原文自带的限制与不可外推边界

- **原文位置：** PDF p. 6 的 Eq. (40) 前后；PDF p. 7 的 Eqs. (42)–(43) 后两段；PDF p. 26，文章页 066614-26，§VII Conclusion。
- **平壁与方向限制：** 本账本的 Eq. (42) halfway 链条针对平行平壁、网格对齐链接与主坐标轴 forcing。p. 7 明说，对任意倾斜 Poiseuille 流，在恒定 `λ_ν,λ_2,I_f` 下不能对所有 `δ_q` 消掉误差。
- **线性、稳态与 forcing 限制：** 精确剖面的取消对实际算法细节敏感；原文后文特别要求 forcing 项与 Eq. (2) 一致，动量按 Eq. (9) 及 `I_f=-1/2` 定义，而非线性平衡会引入非零高阶导数并破坏这一精确解。p. 26 还明说本文理论分析是对稳态流做的，非稳态扩展仍待完成。
- **证明等级限制：** p. 7 关于复杂流中保持常数 `Λ²` 可使渗透率与黏性无关的说法，作者明确标注为尚未理论证明、仅被数值模拟强烈支持。不得将它改写为任意复杂几何的已证严格定理。
- **多反射方法限制：** p. 26 说明 multireflection 沿一条链接至少需要三个流体节点；少于三点时的替代做法未被全面研究。这条是 multireflection 方法限制，不应倒推为 Eq. (42) 的三点前提。
- **总体不可外推：** 本论文不支持把 `Λ²=1/4 ⇔ H_evenH_odd=3/16` 宣称为曲面、移动边界、粗糙几何、非稳态/非线性流、其他 forcing gauge、温度 ABB 或整个热对流求解器的普适精度条件。

## 目视核对记录

| 页图 | 文章页 | 核对内容 | 观察结果与批注层处理 |
| --- | --- | --- | --- |
| `rendered/ginzburg2003-pdf-p3.png` | 066614-3 | 两率简化模型；Eqs. (9)–(12) 的 `I_f`、黏性与线性 Stokes 语境 | 印刷正文可读；黄色高亮为批注/标记层，全部排除于证据之外 |
| `rendered/ginzburg2003-pdf-p6.png` | 066614-6 | Eqs. (38)–(41)：抛物线剖面、有效宽度通式和 `Λ²` 定义 | 式号、上下文与分母清晰；未见遮盖公式的批注 |
| `rendered/ginzburg2003-pdf-p7.png` | 066614-7 | Eqs. (42)–(43)；`I_f=-1/2`、`δ_q=1/2`、有效直径、倾斜流限制和“未理论证明”声明 | 公式与限制段落清晰；未见遮盖正文的批注 |
| `rendered/ginzburg2003-pdf-p26.png` | 066614-26 | 结论中的稳态适用范围、非稳态待扩展与 multireflection 三节点限制 | 正文清晰，未见遮盖性批注 |

四张页图均以 180 dpi 从本地 PDF 渲染，并在本次接管中重新目视查看。布局文本只承担定位功能。
