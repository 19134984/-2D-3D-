# Dubois–Lallemand 2009 证据账本

## 论文身份

- 作者：François Dubois, Pierre Lallemand。
- 标题：*Towards higher order lattice Boltzmann schemes*。
- 载体：*Journal of Statistical Mechanics: Theory and Experiment* (2009), P06006；本地 PDF 首页标明 `arXiv:0811.0599v3`。
- 本地文件：`pdf/Towards higher order lattice Boltzmann schemes.pdf`。
- 页码口径：本地 PDF 页与印刷页一致。

## 原符号与统一符号

| 原符号 | 原文含义 | 统一符号 | 精确映射 |
| --- | --- | --- | --- |
| `s_k` | 第 `k` 个非守恒矩的碰撞率 | `r_k` | `r_k=s_k` |
| `σ_k` | Eq. (16) 定义的 Hénon shift | `H_k` | `H_k=σ_k=1/s_k-1/2` |
| `σ_1,σ_3,σ_4` | D2Q5 扩散、能量、高阶模态的 Hénon shift | `H_κ,H_e,H_ν` | `H_1=H_κ`, `H_3=H_e`, `H_4=H_ν` |
| `α` | D2Q5 平衡能量矩参数 | `a` | `a=α` |
| `λ=Δx/Δt` | 格子速度尺度 | `λ` | 保持不变；应用论文取格子单位 `λ=1` |

> 本文的裸 `σ_k` 是 Hénon 参数，不是 Wang/Contrino 中作为碰撞率的 `σ_κ,σ_e,σ_ν`。

## 核心证据

### 1. Hénon 参数定义（Eq. (16)）

- **原文位置：** PDF/印刷 p. 5，Eq. (16)。
- **原公式：** `σ_k ≡ 1/s_k - 1/2`。
- **统一符号：**

  \[
  H_k\equiv\frac1{r_k}-\frac12=\tau_k-\frac12,
  \qquad r_k=s_k,\quad \tau_k=1/r_k.
  \]

- **允许主张：** 本教程的 `H` 与 Dubois 的 `σ` 完全同义，是定义和精确代数，不是渐近近似。
- **不支持 / 不可外推：** 不得把 Dubois `σ_k` 当成碰撞率；也不能仅由定义推出 `0<r_k<2`，该区间若使用必须另列模型假设。

### 2. D2Q5 矩阵（Appendix Eq. (79)）

- **原文位置：** PDF/印刷 p. 38，Appendix B，Eq. (79)。
- **原公式：**

  \[
  M=\begin{pmatrix}
  1&1&1&1&1\\
  0&\lambda&0&-\lambda&0\\
  0&0&\lambda&0&-\lambda\\
  -4&1&1&1&1\\
  0&1&-1&1&-1
  \end{pmatrix}.
  \]

- **统一符号：** 格子单位 `λ=1` 时与 Wang/Contrino 的 `N` 完全一致；列顺序为 rest、`+x,+y,-x,-y`。
- **允许主张：** Appendix Eq. (79) 支持应用论文 D2Q5 矩阵的来源和模态编号。
- **不支持 / 不可外推：** `λ≠1` 时第二、三行含速度尺度；不得不经单位和列序检查就声称两个矩阵逐元素相同。

### 3. D2Q5 四阶等效方程与系数（Eqs. (39)–(42)）

- **原文位置：** PDF/印刷 p. 11，Eqs. (39)–(42)。
- **原文结构：** Eq. (39) 给纯热 D2Q5 碰撞矩阵 `Ψ`；Eq. (40) 给守恒律至四阶展开，其中四阶空间导数项为

  \[
  \kappa_{40}(\partial_x^4\rho+\partial_y^4\rho)
  +\kappa_{22}\partial_x^2\partial_y^2\rho.
  \]

  Eqs. (41)–(42) 给

  \[
  \kappa_{40}=8-3\alpha+12(\alpha+4)\sigma_1^2
  -12(1-\alpha)\sigma_1\sigma_3-60\sigma_1\sigma_4,
  \]

  \[
  \kappa_{22}=-6(\alpha+4)+24(\alpha+4)\sigma_1^2
  -24(1-\alpha)\sigma_1\sigma_3+120\sigma_1\sigma_4.
  \]

- **统一符号：** `α→a`, `σ_1→H_κ`, `σ_3→H_e`, `σ_4→H_ν`；`κ_40,κ_22` 保持不变。
- **允许主张：** Eqs. (40)–(42) 明确分开轴向四阶导数和混合四阶导数，并给出二者系数。
- **不支持 / 不可外推：** 这些是纯热 D2Q5 体相等效方程系数；不能自动替代有对流速度、边界、D2Q9 或完整热流耦合的误差分析。Eq. (40) 的 `O(Δt^4)` 也不能脱离其渐近展开语境解释为任意量的全局误差。

### 4. 由 Eqs. (41)–(42) 独立复算各向同性

- **原文位置：** PDF/印刷 p. 11，Eqs. (41)–(42)；以下差分为本账本的精确代数复算。
- **复算：**

  \[
  \begin{aligned}
  \kappa_{22}-2\kappa_{40}
  &=[-6(\alpha+4)-2(8-3\alpha)]\\
  &\quad+[24(\alpha+4)-24(\alpha+4)]\sigma_1^2\\
  &\quad+[-24(1-\alpha)+24(1-\alpha)]\sigma_1\sigma_3\\
  &\quad+[120+120]\sigma_1\sigma_4\\
  &=40(6\sigma_1\sigma_4-1).
  \end{aligned}
  \]

  要使四阶组合与双调和算子

  \[
  (\partial_x^2+\partial_y^2)^2
  =\partial_x^4+2\partial_x^2\partial_y^2+\partial_y^4
  \]

  成比例，需要 `κ_22=2κ_40`，故 `σ_1σ_4=1/6`，即 `H_κH_ν=1/6`。
- **允许主张：** 该复算支持 Wang/Contrino Eq. (17) 是 D2Q5 四阶各向同性关系。
- **不支持 / 不可外推：** `κ_22=2κ_40` 只规定方向比例，并未要求 `κ_40=κ_22=0`；各向同性不等于零误差或全局四阶。

### 5. 完整四阶消除（Eq. (55)）

- **原文位置：** PDF/印刷 p. 16，Eq. (55) 及前一句。实际页码是 16，不采用 Contrino 批注层所写的 “Page 15”。
- **原公式：** 原文明确说令 Eqs. (41)–(42) 中 `κ_40=0`、`κ_22=0`，得到

  \[
  \sigma_3=\sigma_1\frac{\alpha+4}{1-\alpha}
  -\frac{1}{12\sigma_1}\frac{2+3\alpha}{1-\alpha},
  \qquad
  \sigma_4=\frac{1}{6\sigma_1}.
  \]

- **统一符号：**

  \[
  H_e=H_\kappa\frac{a+4}{1-a}
  -\frac{2+3a}{12H_\kappa(1-a)},
  \qquad H_\nu=\frac1{6H_\kappa}.
  \]

  分母条件：`H_κ≠0`、`a≠1`。转换到碰撞率时还要求各 `r_j≠0`。
- **允许主张：** Eq. (55) 支持在该纯热 D2Q5 等效方程中同时令两个四阶系数为零；第二式同时包含各向同性乘积。
- **不支持 / 不可外推：** Eq. (55) 不证明三阶、边界或耦合流场误差为零，不支持把整个热对流程序称作全局四阶，也不适用于 `a=1` 或 `H_κ=0` 的直接代入。

### 6. TRT 特殊点

- **原文位置：** PDF/印刷 p. 16，Eq. (55) 后紧接段落。
- **原文声明：** 施加 intermediate TRT `σ_3=σ_4`，原文给

  \[
  \sigma_1=\frac1{\sqrt{12}},\qquad \sigma_3=\sigma_4=\frac1{\sqrt3}.
  \]

- **统一符号与率：**

  \[
  H_\kappa=\frac1{\sqrt{12}},\quad H_e=H_\nu=\frac1{\sqrt3},
  \]

  \[
  r_\kappa=3-\sqrt3,\qquad r_e=r_\nu=4\sqrt3-6.
  \]

  对 `a≠1` 回代 Eq. (55)：`1/(12H_κ)=H_κ`，从而第一式右端为 `2H_κ=1/√3`；第二式给 `1/(6H_κ)=1/√3`。
- **允许主张：** 该正 Hénon 特殊点分别把 `κ_40`、`κ_22` 置零，并与 Contrino Eqs. (19)–(20) 一致。
- **不支持 / 不可外推：** 原文此段给出的正分支不自动排除其他代数分支或增加稳定性结论；参数可用性还需物理扩散率和稳定性条件另行检查。

### 7. 平面波验证及其证据等级

- **原文位置：** PDF/印刷 p. 17，Eq. (55) 后的第一种验证方法及 Eqs. (56)–(57) 开端。
- **原文要点：** 对周期 homogeneous plane wave 构造放大矩阵，论文说数值实验确认色散关系的理论展开；随后讨论有限格点特征模态。
- **允许主张：** p. 17 提供对体相色散/特征值预测的数值核对思路。
- **不支持 / 不可外推：** 平面波核对不是对 Contrino 方腔、BB/ABB 边界或任意非线性热对流的数值证明；它不能替代 Eq. (55) 的适用条件。

## 目视核对记录

| 页图 | 印刷页 | 核对内容 | 观察结果 |
| --- | ---: | --- | --- |
| `rendered/dubois2009-pdf-p5.png` | 5 | Eq. (16) Hénon 参数定义 | 正文清晰，无可见批注层 |
| `rendered/dubois2009-pdf-p11.png` | 11 | Eqs. (39)–(42) 与四阶系数 | 公式完整，无截断或覆盖 |
| `rendered/dubois2009-pdf-p16.png` | 16 | Eq. (55)、TRT 条件和特殊点 | 实际印刷页 16；正文清晰 |
| `rendered/dubois2009-pdf-p17.png` | 17 | 平面波/特征模态核对上下文 | 验证对象为体相色散，不是应用论文边界 |
| `rendered/dubois2009-pdf-p38.png` | 38 | Appendix Eq. (79) D2Q5 矩阵 | 矩阵与 `λ` 因子清晰，无可见批注层 |

所有页图均以 180 dpi 从本地 PDF 渲染；布局文本只用于定位。
