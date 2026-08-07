---
title: 高 Rayleigh 数热对流 LBM 稳定性增强算法——完整推导、参数区域与验证边界
status: 理论与符号推导完成，待数值验证
version: 1.0
date: 2026-07-30
completed_date: 2026-07-30
language: zh-CN
---

# 0. 文件用途

本文件整理当前已经讨论并确认的算法思路，目标是作为后续使用 Codex 继续进行符号推导、等效方程分析、边界条件推导和程序实现的统一入口。

本文件严格区分三类内容：

- **已确认**：可直接作为后续推导的前提；
- **工作假设**：目前合理，但必须通过推导验证；
- **待推导**：不得在论文或代码说明中写成已证明结论。

当前优先路线为：

1. 保留现有双分布函数 MRT/TRT 热 LBM 框架；
2. 流场加入 LBM-CDE 型应力修正源项；
3. 温度场优先保留 D2Q5，采用“主分支 A”；
4. 不预设原流场魔法参数 \(3/16\) 在修正后仍成立；
5. 不预设原温度四阶关系 \(1/6\) 在修正后仍成立；
6. 先重新推导 D2Q5 是否存在可接受参数解；
7. 若 D2Q5 确实无法兼顾精度与稳定性，再升级温度格子至 D2Q9。

---

# 1. 研究目标

研究对象为高 Rayleigh 数下的自然对流和 Rayleigh–Bénard 对流。

在高 \(Ra\) 下，物理黏度 \(\nu\) 与热扩散率 \(\kappa\) 很小，传统热 LBM 容易出现：

- 物理输运模态的净阻尼减弱；
- 某些高阶或 ghost 模态因魔法参数约束趋于弱阻尼；
- 温度平衡分布权重退化；
- 棋盘格、边界层振荡、短波误差累积；
- 热边界层和羽流区局部失稳；
- 半程反弹边界位置可能随参数变化；
- 高阶精度条件与最优稳定性条件发生冲突。

本研究的目标不是完全替换现有算法，而是在尽可能保留原模型优势的前提下：

\[
\boxed{
\text{将物理小输运系数与基础数值碰撞参数部分解耦}
}
\]

并检查能否同时保持：

- 正确的 Navier–Stokes 方程和温度对流扩散方程；
- Boussinesq 浮力耦合；
- 二阶时间精度；
- 原有无滑移边界精度；
- D2Q5 温度模型的高阶误差性质；
- 高 \(Ra\) 下更好的稳定性；
- GPU/OpenACC 友好的局部计算结构。

---

# 2. 基线算法

## 2.1 二维模型

流场：

\[
\mathrm{D2Q9}
\]

温度场：

\[
\mathrm{D2Q5}
\]

流场采用 MRT/TRT 碰撞，温度场采用 MRT/TRT 碰撞。

## 2.2 三维模型

流场：

\[
\mathrm{D3Q19}
\]

温度场：

\[
\mathrm{D3Q7}
\]

二维算法与三维算法的基本思想一致，三维版本是二维版本的推广。

## 2.3 宏观控制方程

采用 Boussinesq 近似：

\[
\nabla\cdot \boldsymbol u=0,
\]

\[
\frac{\partial \boldsymbol u}{\partial t}
+\boldsymbol u\cdot\nabla\boldsymbol u
=
-\frac{1}{\rho_0}\nabla p
+\nu\nabla^2\boldsymbol u
+\boldsymbol F_b,
\]

\[
\frac{\partial T}{\partial t}
+\boldsymbol u\cdot\nabla T
=
\kappa\nabla^2T,
\]

其中浮力为

\[
\boldsymbol F_b
=
\rho_0 g\beta_T(T-T_0)\hat{\boldsymbol g}.
\]

---

# 3. 统一记号

为了避免松弛率、松弛时间、Hénon 参数和碰撞特征值混淆，统一采用：

\[
s=\frac{1}{\tau},
\]

\[
\sigma=\frac{1}{s}-\frac12=\tau-\frac12,
\]

\[
\lambda=-s.
\]

其中：

- \(s\in(0,2)\)：松弛率；
- \(\tau>1/2\)：通常 LBM 文献中的松弛时间；
- \(\sigma>0\)：Hénon 参数；
- \(\lambda\in(-2,0)\)：部分 TRT 文献中的碰撞特征值。

碰撞后某非守恒矩的非平衡部分满足：

\[
m^{neq,*}=(1-s)m^{neq}.
\]

因此单步放大因子为

\[
G=1-s.
\]

三种极限：

\[
s\rightarrow 1
\quad\Longrightarrow\quad
G\rightarrow 0,
\]

表示该模态在一次碰撞中被完全消除；

\[
s\rightarrow 2
\quad\Longrightarrow\quad
G\rightarrow -1,
\]

表示该模态几乎不衰减，只发生符号交替；

\[
s\rightarrow 0
\quad\Longrightarrow\quad
G\rightarrow 1,
\]

表示该模态几乎原样保留。

所以 \(s\to 0\) 和 \(s\to 2\) 都属于弱阻尼极限。

---

# 4. 高 \(Ra\) 下的原始问题

## 4.1 流场物理黏度模态

标准 LBM 中：

\[
\nu
=
c_s^2
\left(
\frac{1}{s_\nu}-\frac12
\right)\delta t
=
c_s^2\sigma_\nu\delta t.
\]

当

\[
\nu\rightarrow 0
\]

时：

\[
\sigma_\nu\rightarrow 0,
\qquad
s_\nu\rightarrow 2,
\qquad
1-s_\nu\rightarrow -1.
\]

这意味着剪切应力非平衡模态只发生近似符号翻转，而不能被有效衰减。

这部分弱阻尼与物理小黏度本身有关，不能通过简单“强行增大耗散”彻底消除，否则会改变目标物理黏度。

## 4.2 原流场魔法参数引起的伴随弱阻尼

原算法采用类似：

\[
\Lambda_f
=
\sigma_\nu\sigma_q
=
\frac{3}{16}.
\]

若仍以物理黏度对应的

\[
\sigma_\nu\rightarrow 0
\]

代入，则：

\[
\sigma_q
=
\frac{3}{16\sigma_\nu}
\rightarrow \infty,
\]

从而：

\[
s_q
=
\frac{1}{\sigma_q+1/2}
\rightarrow 0.
\]

于是同时出现：

\[
s_\nu\to 2,
\qquad
s_q\to 0.
\]

即：

- 物理剪切应力模态接近 \(-1\) 放大；
- 与边界精度耦合的高阶模态接近 \(+1\) 放大；
- 两个模态均弱阻尼；
- 误差可能分别以交替形式和同号形式累积。

---

# 5. D2Q5 温度模型的平衡分布

## 5.1 D2Q5 离散速度

采用：

\[
\boldsymbol e_0=(0,0),
\]

\[
\boldsymbol e_1=(1,0),\quad
\boldsymbol e_2=(0,1),\quad
\boldsymbol e_3=(-1,0),\quad
\boldsymbol e_4=(0,-1).
\]

## 5.2 温度矩

采用矩向量：

\[
\boldsymbol n
=
\begin{bmatrix}
n_0&n_x&n_y&n_e&n_\nu
\end{bmatrix}^{\mathrm T}.
\]

平衡矩为：

\[
\boldsymbol n^{eq}
=
\begin{bmatrix}
T&u_xT&u_yT&aT&0
\end{bmatrix}^{\mathrm T}.
\]

矩变换矩阵取：

\[
N=
\begin{bmatrix}
1&1&1&1&1\\
0&1&0&-1&0\\
0&0&1&0&-1\\
-4&1&1&1&1\\
0&1&-1&1&-1
\end{bmatrix}.
\]

## 5.3 平衡分布的显式形式

由

\[
\boldsymbol g^{eq}=N^{-1}\boldsymbol n^{eq}
\]

得到：

\[
g_0^{eq}
=
\frac{1-a}{5}T,
\]

\[
g_{+x}^{eq}
=
\frac{4+a}{20}T+\frac12Tu_x,
\]

\[
g_{-x}^{eq}
=
\frac{4+a}{20}T-\frac12Tu_x,
\]

\[
g_{+y}^{eq}
=
\frac{4+a}{20}T+\frac12Tu_y,
\]

\[
g_{-y}^{eq}
=
\frac{4+a}{20}T-\frac12Tu_y.
\]

定义：

\[
c_T^2
=
\frac{4+a}{10},
\]

则移动方向的平衡分布可写成：

\[
\boxed{
g_{\pm\alpha}^{eq}
=
\frac{T}{2}
\left(
c_T^2\pm u_\alpha
\right)
},
\qquad
\alpha=x,y.
\]

---

# 6. D2Q5 平衡分布非负条件

这里必须区分“平衡分布非负”和“线性稳定”。

假定：

\[
T>0.
\]

要求：

\[
g_{+\alpha}^{eq}\ge 0,
\qquad
g_{-\alpha}^{eq}\ge 0,
\]

得到：

\[
c_T^2+u_\alpha\ge0,
\]

\[
c_T^2-u_\alpha\ge0.
\]

所以：

\[
\boxed{
|u_\alpha|\le c_T^2
}.
\]

这不是近似，而是在当前 D2Q5 一阶速度平衡分布下的精确非负条件。

静止分布的非负要求：

\[
g_0^{eq}\ge0
\quad\Longrightarrow\quad
a\le1.
\]

移动方向对称权重非负要求：

\[
c_T^2\ge0
\quad\Longrightarrow\quad
a\ge-4.
\]

因此完整参数范围为：

\[
\boxed{
-4\le a\le1
}
\]

以及：

\[
\boxed{
\max(|u_x|,|u_y|)\le c_T^2
}.
\]

若只用速度模长表示，一个对任意方向都成立的充分条件为：

\[
|\boldsymbol u|^2\le (c_T^2)^2.
\]

注意这里右端是 \((c_T^2)^2\)，不是 \(c_T^2\)。

---

# 7. 非负条件与长波必要稳定条件的区别

## 7.1 未补偿数值扩散时的有效扩散张量

对于最小速度集 D2Q5，若平衡分布中没有完整的 \(u_\alpha u_\beta\) 二阶补偿，则长波等效扩散张量含有：

\[
\boldsymbol D_{\mathrm{eff}}
=
c_T^2\boldsymbol I
-
\boldsymbol u\boldsymbol u^{\mathrm T}.
\]

在二维中：

\[
\boldsymbol D_{\mathrm{eff}}
=
\begin{bmatrix}
c_T^2-u_x^2&-u_xu_y\\
-u_xu_y&c_T^2-u_y^2
\end{bmatrix}.
\]

要求该张量半正定，得到：

\[
\boxed{
|\boldsymbol u|^2\le c_T^2
}.
\]

这是长波极限下的必要稳定条件，而不是平衡分布非负条件。

## 7.2 两个条件的比较

平衡分布非负条件：

\[
\boxed{
\max_\alpha |u_\alpha|\le c_T^2
}
\]

或保守地写成：

\[
|\boldsymbol u|^2\le(c_T^2)^2.
\]

长波必要稳定条件：

\[
\boxed{
|\boldsymbol u|^2\le c_T^2
}.
\]

当：

\[
0<c_T^2\ll1
\]

时：

\[
(c_T^2)^2\ll c_T^2.
\]

所以平衡分布非负条件通常比长波必要稳定条件更严格。

平衡分布出现负值，并不自动意味着 von Neumann 失稳；反之，满足平衡分布非负也不能代替完整稳定性分析。

## 7.3 温度变量的平移问题

如果程序直接演化：

\[
\theta\in[-0.5,0.5],
\]

则分布函数会随 \(\theta\) 改变符号，此时不能简单使用人口非负性作为物理判据。

更适合的做法是演化：

\[
T\in[0,1]
\]

或其他整体为正的温度变量，而浮力使用：

\[
T-T_0.
\]

---

# 8. 原 D2Q5 中小 \(\kappa\) 的实现方式

原模型中：

\[
\kappa
=
c_T^2
\left(
\frac{1}{q_\kappa}-\frac12
\right)\delta t.
\]

原 Luo 类算法固定：

\[
\sigma_\kappa
=
\frac{1}{q_\kappa}-\frac12
=
\frac{\sqrt3}{6},
\]

因此：

\[
q_\kappa
=
3-\sqrt3.
\]

于是小热扩散率主要通过：

\[
c_T^2\rightarrow0
\]

实现。

由于：

\[
c_T^2=\frac{4+a}{10},
\]

所以：

\[
\kappa\rightarrow0
\quad\Longrightarrow\quad
a\rightarrow-4.
\]

这会导致移动方向的对称平衡权重：

\[
\frac12c_T^2T
\]

趋于零，而对流反对称部分：

\[
\pm\frac12Tu_\alpha
\]

不随 \(c_T^2\) 一起消失。

因此高 \(Ra\)、小 \(\kappa\) 时容易出现：

- 平衡分布负值；
- 有效扩散张量半正定条件变严；
- 对流项相对扩散项过强；
- 棋盘格和短波振荡；
- 温度边界层附近误差放大。

---

# 9. LBM-CDE 型输运系数解耦思路

## 9.1 流场目标关系

引入应力修正参数 \(\chi_s\)，使物理黏度为：

\[
\boxed{
\nu
=
(1-\chi_s)
c_s^2
\sigma_\nu^{(0)}
\delta t
}
\]

其中：

\[
\sigma_\nu^{(0)}
=
\frac{1}{s_\nu^{(0)}}-\frac12
\]

为基础碰撞算子中的数值参数。

定义物理等效 Hénon 参数：

\[
\boxed{
\sigma_{\nu,\mathrm{eff}}
=
(1-\chi_s)\sigma_\nu^{(0)}
=
\frac{\nu}{c_s^2\delta t}
}.
\]

## 9.2 温度场目标关系

引入扩散修正参数 \(\chi_\kappa\)，使物理热扩散率为：

\[
\boxed{
\kappa
=
(1-\chi_\kappa)
c_T^2
\sigma_\kappa^{(0)}
\delta t
}
\]

其中：

\[
\sigma_\kappa^{(0)}
=
\frac{1}{q_\kappa^{(0)}}-\frac12.
\]

定义：

\[
\boxed{
\sigma_{\kappa,\mathrm{eff}}
=
(1-\chi_\kappa)\sigma_\kappa^{(0)}
=
\frac{\kappa}{c_T^2\delta t}
}.
\]

这样可以固定一个适中的 \(c_T^2\)，不再依靠 \(a\to-4\) 实现小 \(\kappa\)。

---

# 10. 需要正确理解“稳定性提高”的来源

LBM-CDE 型方法不能把物理剪切应力或物理热通量模态任意强阻尼，因为这些模态的净耗散必须与目标 \(\nu\) 和 \(\kappa\) 一致。

其主要作用是：

1. 基础碰撞时间不必直接进入极端值；
2. 高阶非水动力模态不必全部跟随物理输运参数；
3. 平衡权重不必因小 \(\kappa\) 退化至零；
4. 可独立调整基础碰撞谱和物理输运系数；
5. 可为 ghost 模态选择更合理的阻尼；
6. 为后续正则化和选择性过滤提供空间。

但物理应力或通量的净有效松弛仍会随：

\[
\nu\to0
\quad\text{或}\quad
\kappa\to0
\]

进入弱耗散极限。

---

# 11. 流场 \(3/16\) 关系为什么必须重新推导

原关系：

\[
\boxed{
\sigma_\nu\sigma_q=\frac{3}{16}
}
\]

是在无附加应力修正源项的碰撞—传播—边界系统中推导的。

加入 LBM-CDE 型源项后，存在两个不同参数：

基础碰撞参数：

\[
\sigma_\nu^{(0)}
\]

物理等效参数：

\[
\sigma_{\nu,\mathrm{eff}}
=
(1-\chi_s)\sigma_\nu^{(0)}.
\]

如果简单采用：

\[
\sigma_{\nu,\mathrm{eff}}\sigma_q
=
\frac{3}{16},
\]

则当：

\[
\nu\to0
\]

时：

\[
\sigma_{\nu,\mathrm{eff}}\to0,
\]

所以：

\[
\sigma_q\to\infty,
\qquad
s_q\to0.
\]

即 \(q\) 模态仍然弱阻尼。

如果改用：

\[
\sigma_\nu^{(0)}\sigma_q
=
\frac{3}{16},
\]

则虽然可以让 \(s_q\) 保持适中，但不能直接保证原无滑移壁面位置仍然正确。

因此不能预先选择任一关系。

必须重新推导：

\[
\boxed{
\mathcal B
\left(
\sigma_\nu^{(0)},
\chi_s,
\sigma_q
\right)=0
}
\]

其中 \(\mathcal B=0\) 表示源项修正后的半程反弹无滑移闭合条件。

要求：

\[
\mathcal B
\left(
\sigma_\nu^{(0)},0,\sigma_q
\right)
=
\sigma_\nu^{(0)}\sigma_q-\frac{3}{16}.
\]

## 11.1 需要检查的三种可能结果

### 情形 1

新边界条件主要依赖：

\[
\sigma_{\nu,\mathrm{eff}}.
\]

则：

\[
s_q\to0
\]

可能仍不可避免。

### 情形 2

新边界条件主要依赖：

\[
\sigma_\nu^{(0)}.
\]

则可能同时保持：

- 小物理黏度；
- 适中的基础剪切松弛率；
- 适中的 \(q\) 模态松弛率。

### 情形 3

标准 halfway bounce-back 无法同时满足目标，需要增加源项一致的边界修正项。

当前不得把“修正后仍直接采用 \(3/16\)”写成既定算法。

---

# 12. 温度主分支 A

当前确认优先采用：

\[
\boxed{
\text{D2Q5 + LBM-CDE 型扩散修正源项}
}
\]

而不是立即升级 D2Q9。

主分支 A 的目标是：

- 固定适中的 \(a\)；
- 固定适中的 \(c_T^2\)；
- 使用 \(\chi_\kappa\) 实现小 \(\kappa\)；
- 重新推导四阶误差条件；
- 检查 D2Q5 是否存在可接受参数窗口；
- 只有 D2Q5 确实失败后才升级格子。

---

# 13. 原 D2Q5 的 \(1/6\) 关系

设：

\[
\sigma_1
=
\frac{1}{q_\kappa}-\frac12
\]

对应两个温度通量矩；

设：

\[
\sigma_3
\]

对应能量型偶矩；

设：

\[
\sigma_4
\]

对应另一个高阶偶矩。

无源纯扩散 D2Q5 的四阶完全消除条件为：

\[
\boxed{
\sigma_4
=
\frac{1}{6\sigma_1}
}
\]

即：

\[
\boxed{
\sigma_1\sigma_4=\frac16
}
\]

以及：

\[
\boxed{
\sigma_3
=
\sigma_1\frac{a+4}{1-a}
-
\frac{2+3a}
{12\sigma_1(1-a)}
}.
\]

这里的 \(1/6\) 不是一个孤立条件。它与 \(\sigma_3\) 的关系共同用于消除两个独立四阶误差系数。

加入扩散修正源项后，这两个关系都需要重新推导。

---

# 14. 纯扩散下有效温度通量松弛率

这是后续推导必须复核的关键中间结果。

## 14.1 基础通量松弛率

定义：

\[
q_\kappa^{(0)}
\]

以及：

\[
\sigma_\kappa^{(0)}
=
\frac{1}{q_\kappa^{(0)}}-\frac12.
\]

## 14.2 扩散修正源项

温度扩散修正源项的一阶矩写成：

\[
R_\alpha
=
\chi_\kappa c_T^2\partial_\alpha T.
\]

离散源项需要包含半步修正因子。

在纯扩散、线性问题中，将局部温度梯度表达式代回通量矩碰撞，可以得到一个净有效松弛率：

\[
\boxed{
q_{\kappa,\mathrm{eff}}
=
\frac{
2q_\kappa^{(0)}
}{
2(1-\chi_\kappa)
+
\chi_\kappa q_\kappa^{(0)}
}
}
\]

对应：

\[
\boxed{
\sigma_{\kappa,\mathrm{eff}}
=
\frac{1}{q_{\kappa,\mathrm{eff}}}-\frac12
=
(1-\chi_\kappa)
\left(
\frac{1}{q_\kappa^{(0)}}-\frac12
\right)
}.
\]

所以：

\[
\boxed{
\kappa
=
c_T^2
\sigma_{\kappa,\mathrm{eff}}
\delta t
}.
\]

这说明，在纯扩散线性极限中，扩散修正源项的净作用等价于改变温度通量模态的有效 Hénon 参数。

---

# 15. \(1/6\) 关系加入修正源项后的核心矛盾

如果重新推导后，四阶条件仍然只依赖有效通量参数，并保持类似：

\[
\boxed{
\sigma_{\kappa,\mathrm{eff}}\sigma_4
=
\frac16
},
\]

则当：

\[
\kappa\to0
\]

且 \(c_T^2\) 固定时：

\[
\sigma_{\kappa,\mathrm{eff}}
=
\frac{\kappa}{c_T^2\delta t}
\to0.
\]

于是：

\[
\sigma_4
=
\frac{1}{6\sigma_{\kappa,\mathrm{eff}}}
\to\infty,
\]

从而：

\[
q_4
=
\frac{1}{\sigma_4+1/2}
\to0.
\]

即高阶温度模态仍然会进入弱阻尼极限。

因此 LBM-CDE 虽然能避免：

\[
a\to-4,
\]

却不自动保证原四阶精度关系与高阶模态强阻尼可以同时满足。

这正是必须重新推导 \(1/6\) 的原因。

---

# 16. D2Q5 是否可行的判定方式

不能仅根据 \(\kappa\to0\) 极限直接判定 D2Q5 失败，因为实际计算中的 \(\kappa\) 是有限小量。

应给定高阶松弛率下限：

\[
q_{4,\min}>0.
\]

要求：

\[
q_4\ge q_{4,\min}.
\]

由于：

\[
q_4
=
\frac{1}{\sigma_4+1/2},
\]

得到：

\[
\sigma_4
\le
\frac{1}{q_{4,\min}}-\frac12.
\]

若四阶条件为：

\[
\sigma_{\kappa,\mathrm{eff}}\sigma_4=\frac16,
\]

则：

\[
\sigma_{\kappa,\mathrm{eff}}
\ge
\frac{
1
}{
6\left(
1/q_{4,\min}-1/2
\right)
}.
\]

于是最小允许热扩散率为：

\[
\boxed{
\kappa_{\min}
=
\frac{
c_T^2\delta t
}{
6\left(
1/q_{4,\min}-1/2
\right)
}
}.
\]

将目标 \(Ra\)、\(Pr\)、网格尺度、Ma 和无量纲化方式代入，即可判断 D2Q5 是否进入不可接受区域。

---

# 17. D2Q5 重新推导的正确顺序

## A1. 建立完整 D2Q5 源项模型

明确：

\[
\boldsymbol n^{eq}
=
[T,u_xT,u_yT,aT,0]^{\mathrm T}.
\]

明确松弛矩阵：

\[
Q
=
\operatorname{diag}
\left(
0,
q_\kappa^{(0)},
q_\kappa^{(0)},
q_3,
q_4
\right).
\]

明确温度源项：

\[
\boldsymbol\Psi
=
\left(
I-\frac{Q}{2}
\right)
\boldsymbol R.
\]

必须给出 \(\boldsymbol R\) 在矩空间中的显式分量。

## A2. 先做无对流纯扩散展开

令：

\[
\boldsymbol u=0.
\]

推导到四阶的等效方程：

\[
\partial_tT
=
\kappa\nabla^2T
+
\Delta t^2\mathcal E_3
+
\Delta t^3\mathcal E_4
+
O(\Delta t^4).
\]

重点得到两个四阶系数：

\[
\widetilde\kappa_{40},
\qquad
\widetilde\kappa_{22}.
\]

它们应写成：

\[
\widetilde\kappa_{40}
=
\widetilde\kappa_{40}
\left(
a,
q_\kappa^{(0)},
q_3,
q_4,
\chi_\kappa
\right),
\]

\[
\widetilde\kappa_{22}
=
\widetilde\kappa_{22}
\left(
a,
q_\kappa^{(0)},
q_3,
q_4,
\chi_\kappa
\right).
\]

## A3. 求两类条件

四阶各向同性：

\[
\boxed{
\widetilde\kappa_{22}
=
2\widetilde\kappa_{40}
}
\]

四阶完全消除：

\[
\boxed{
\widetilde\kappa_{40}=0,
\qquad
\widetilde\kappa_{22}=0
}
\]

## A4. 检查无源极限

当：

\[
\chi_\kappa=0
\]

时，结果必须退化到原 D2Q5 条件：

\[
\sigma_1\sigma_4=\frac16
\]

以及原 \(\sigma_3\) 关系。

## A5. 加入恒定对流速度

再令：

\[
\boldsymbol u=\text{constant}.
\]

推导完整对流扩散等效方程，检查：

- \(u_x^2\partial_{xx}T\)；
- \(u_y^2\partial_{yy}T\)；
- \(u_xu_y\partial_{xy}T\)；
- 三阶色散项；
- 四阶各向异性项；
- 源项与对流项的交叉误差。

## A6. 判断 D2Q5 是否足够

只有满足以下条件时才保留 D2Q5：

\[
-4<a<1,
\]

\[
0<q_\kappa^{(0)},q_3,q_4<2,
\]

并最好满足：

\[
q_i\in[\varepsilon,2-\varepsilon].
\]

还需保证：

\[
|\boldsymbol u|^2\le c_T^2
\]

的长波必要稳定范围在目标参数区间内可接受。

---

# 18. 何时升级到 D2Q9

D2Q5 无法局部表示交叉二阶矩：

\[
u_xu_y.
\]

因此无法完整补偿：

\[
-u_xu_y
\]

形式的数值交叉扩散。

只有在以下任一情况明确发生时才升级到 D2Q9：

1. 四阶消除方程无可接受实数解；
2. 解要求某个松弛率趋于 \(0\) 或 \(2\)；
3. 目标 \(Ra\) 下 D2Q5 的稳定窗口过窄；
4. 交叉扩散项无法控制；
5. 对角方向出现系统性棋盘格或色散误差；
6. 温度边界层和羽流方向与格线不一致时误差过大。

D2Q9 的理论优势是：

- 含对角速度；
- 可表示 \(u_xu_y\)；
- 可表示 \(\partial_{xy}T\)；
- 可构造完整二阶数值扩散补偿；
- 有更多高阶矩自由度。

但代价是：

- 分布函数数量增加；
- 内存带宽开销增加；
- GPU 性能下降；
- 边界处理更复杂；
- 参数空间更大。

所以当前不优先升级。

---

# 19. 流场源项的待推导结构

流场源项分为：

\[
\boldsymbol\Phi
=
\boldsymbol\Phi^{F}
+
\boldsymbol\Phi^{\nu}.
\]

其中：

\[
\boldsymbol\Phi^F
\]

为 Boussinesq 浮力项，

\[
\boldsymbol\Phi^\nu
\]

为应力修正项。

连续层面的应力修正形式应满足：

\[
S_f^{(\nu)}
\propto
\rho_0
\left(
\frac{c_{i\alpha}c_{i\beta}}{c_s^2}
-\delta_{\alpha\beta}
\right)
\left[
\chi_s S_{\alpha\beta}
+
\frac{\chi_b-\chi_s}{D}
S_{\gamma\gamma}\delta_{\alpha\beta}
\right].
\]

离散源项统一采用：

\[
\boxed{
\boldsymbol\Phi
=
\left(
I-\frac{S}{2}
\right)
\boldsymbol S_f
}.
\]

不得把 Guo 浮力源项和 LBM-CDE 应力源项分别做不一致的时间离散。

---

# 20. D3Q19 偏应力矩说明

三维无迹对称应力张量有五个独立分量：

- 两个独立无迹对角分量；
- 三个剪切分量 \(xy,xz,yz\)。

所以 D3Q19 中应力修正应投影到五个独立偏应力矩。

这五个矩共同对应一个物理剪切黏度，并不表示存在五个不同黏度。

体积黏度或迹部分应投影到独立的能量/体积模态。

后续必须基于代码中实际使用的 \(M\) 矩阵，逐行确定：

- 哪两行为无迹对角应力；
- 哪三行为剪切应力；
- 哪一行为体积应力；
- 每个分量的归一化系数。

---

# 21. 边界条件

第一阶段仍采用：

流场：

\[
\text{halfway bounce-back}
\]

定温边界：

\[
\text{anti-bounce-back}
\]

绝热边界：

\[
\text{bounce-back}
\]

但加入源项后，必须检查：

- 未知入射分布是否需要源项修正；
- 边界温度是否仍位于半格点；
- 流场滑移误差是否仍与原魔法参数一致；
- 温度边界位置是否依赖 \(\chi_\kappa\)；
- Nusselt 数壁面梯度公式是否需要修改；
- 浮力源项和边界反弹次序是否一致。

多重反射边界目前作为理论参考和后续曲面/非均匀网格扩展，不作为第一版方腔算法的必选项。

---

# 22. 正则化和过滤的定位

LBM-CDE 型解耦不是唯一稳定化手段。

若源项解耦后仍存在高频噪声，可进一步采用：

## 22.1 矩空间正则化

只保留所需的低阶非平衡矩：

流场保留：

- 剪切应力；
- 必要的体积应力。

温度场保留：

- 温度通量矩。

删除不需要的高阶 ghost 非平衡分量。

## 22.2 选择性 ghost 模态阻尼

对不影响物理输运系数的高阶矩选择：

\[
s_g\approx 1
\]

或其他具有明显阻尼的参数。

## 22.3 显式空间滤波

只作为最后手段，因为会：

- 引入人工扩散；
- 改变 \(Nu-Ra\) 标度；
- 模糊边界层；
- 抑制真实羽流；
- 难以区分物理耗散与数值耗散。

推荐顺序：

\[
\boxed{
\text{输运系数解耦}
\rightarrow
\text{矩空间正则化}
\rightarrow
\text{选择性 ghost 阻尼}
\rightarrow
\text{必要时显式滤波}
}
\]

---

# 23. Codex 原始推导任务清单（完成结果见第 29--41 节）

## 任务 1：D2Q5 源项的矩空间表达

给出：

\[
\boldsymbol R
\]

及：

\[
\boldsymbol\Psi
=
\left(
I-\frac{Q}{2}
\right)\boldsymbol R.
\]

要求明确每个分量。

## 任务 2：D2Q5 纯扩散四阶等效方程

推导：

\[
\widetilde\kappa_{40},
\qquad
\widetilde\kappa_{22}.
\]

必须通过 \(\chi_\kappa=0\) 的退化检验。

## 任务 3：D2Q5 恒定对流等效方程

检查：

- 数值扩散张量；
- 交叉扩散；
- 三阶色散；
- 四阶各向异性；
- 源项交叉项。

## 任务 4：D2Q5 稳定参数区域

在以下参数空间中搜索：

\[
(a,\chi_\kappa,q_\kappa^{(0)},q_3,q_4).
\]

约束：

\[
-4<a<1,
\]

\[
0<q_i<2,
\]

\[
q_i\in[\varepsilon,2-\varepsilon]
\]

以及目标 \(\kappa\) 和目标速度范围。

## 任务 5：流场源项一致的边界推导

推导：

\[
\mathcal B
\left(
\sigma_\nu^{(0)},
\chi_s,
\sigma_q
\right)=0.
\]

检查是否能避免：

\[
s_q\to0.
\]

## 任务 6：D2Q9/D3Q19 流场源项投影

基于实际 \(M\) 矩阵，给出完整矩空间源项。

## 任务 7：应变率局部表达式

推导：

\[
S_{\alpha\beta}
=
\mathcal F
\left(
\sum_i c_{i\alpha}c_{i\beta}
(f_i-f_i^{eq}),
\boldsymbol F,
\chi_s
\right).
\]

## 任务 8：温度梯度局部表达式

推导：

\[
\nabla T
=
\mathcal G
\left(
\sum_i\boldsymbol c_i
(g_i-g_i^{eq}),
\boldsymbol F,
Q,
\chi_\kappa
\right).
\]

要求避免显式隐式循环。

## 任务 9：边界条件源项修正

分别推导：

- 流场 halfway bounce-back；
- 温度 anti-bounce-back；
- 温度绝热 bounce-back。

## 任务 10：必要时升级 D2Q9

只有在 D2Q5 无可接受解时启动。

---

# 24. Codex 推导要求

Codex 必须遵循以下规则：

1. 所有符号必须先定义；
2. 不得混用 \(s,\tau,\sigma,\lambda\)；
3. 每个矩阵的行列顺序必须固定；
4. 所有源项必须区分连续形式和离散形式；
5. 必须保留半步源项修正；
6. 必须检查守恒矩；
7. 必须检查 \(\chi=0\) 极限；
8. 必须检查 \(u=0\) 极限；
9. 必须检查纯扩散极限；
10. 必须检查 \(\nu,\kappa\to0\) 极限；
11. 必须输出无法闭合的步骤，不得用“显然”跳过；
12. 所有高阶展开尽量使用符号代数复核；
13. 推导结果必须给出可直接转写为 Fortran 的形式；
14. 最终需同时给出碰撞空间和矩空间实现；
15. 不得把工作假设写成已证明结论。

---

# 25. 数值验证计划

## 25.1 纯扩散测试

验证：

- 扩散系数；
- 二阶精度；
- 四阶误差系数；
- 各向同性；
- 高阶松弛率影响；
- \(\chi_\kappa\) 影响。

## 25.2 恒定速度对流扩散

验证：

- 不同速度方向；
- 对角方向对流；
- 数值扩散；
- 相位误差；
- 棋盘格模式；
- 最大稳定速度。

## 25.3 热 Poiseuille / Couette

验证：

- 源项时间离散；
- 温度梯度局部公式；
- 应变率公式；
- 壁面精度；
- 参数变化下的误差。

## 25.4 二维自然对流方腔

建议参数：

\[
Ra=10^3,10^4,10^5,10^6,10^7,10^8.
\]

比较：

- 平均 Nusselt 数；
- 最大局部 Nusselt 数；
- 中心线速度；
- 温度等值线；
- 边界层厚度；
- 稳态/非稳态；
- 最大稳定 \(Ra\)。

## 25.5 Rayleigh–Bénard 对流

比较：

- 临界 \(Ra\)；
- \(Nu-Ra\) 标度；
- \(Re-Ra\) 标度；
- 羽流结构；
- 热耗散；
- 动能耗散；
- 长时间统计量。

## 25.6 稳定窗口比较

至少比较：

1. 原始 D2Q5；
2. D2Q5 + \(\chi_\kappa\)；
3. D2Q5 + \(\chi_\kappa\) + 正则化；
4. 必要时 D2Q9；
5. 流场原 \(3/16\)；
6. 流场新边界关系；
7. 仅改流场；
8. 仅改温度场；
9. 流场和温度场同时修改。

---

# 26. 参数选择原则

## 26.1 温度平衡参数

优先固定：

\[
a
\]

在远离：

\[
-4
\]

和：

\[
1
\]

的位置。

相应：

\[
c_T^2=\frac{4+a}{10}
\]

应具有足够正权重裕量。

## 26.2 基础松弛率

优先避免：

\[
q_\kappa^{(0)}\approx0,
\qquad
q_\kappa^{(0)}\approx2,
\]

\[
q_3\approx0,
\qquad
q_3\approx2,
\]

\[
q_4\approx0,
\qquad
q_4\approx2.
\]

同理流场高阶模态也应避免接近 \(0\) 或 \(2\)。

## 26.3 \(\chi\) 参数

虽然理论上可取：

\[
\chi<1,
\]

但：

\[
\frac{1}{1-\chi}
\]

过大会放大源项和截断误差。

必须通过精度和稳定性共同确定推荐范围。

---

# 27. 推导前的阶段性结论

## 27.1 已确认

1. D2Q5 平衡分布非负条件为：

\[
|u_\alpha|\le c_T^2.
\]

2. 长波必要稳定条件为：

\[
|\boldsymbol u|^2\le c_T^2.
\]

3. 这两个条件不同，不能混用。

4. 流场加入应力修正源项后，原 \(3/16\) 必须重新推导。

5. 温度加入扩散修正源项后，原 \(1/6\) 必须重新推导。

6. 温度优先保留 D2Q5。

7. 只有 D2Q5 无法给出可接受参数解时才升级 D2Q9。

## 27.2 工作假设

在纯扩散线性极限中，扩散修正源项可等效为：

\[
\sigma_{\kappa,\mathrm{eff}}
=
(1-\chi_\kappa)\sigma_\kappa^{(0)}.
\]

该结果应由 Codex 用完整离散方程重新复核。

## 27.3 当时待证明（完成状态见第 41 节）

1. 修正后的流场边界关系；
2. 修正后的温度四阶关系；
3. D2Q5 是否存在兼顾稳定性和高阶精度的参数窗口；
4. 是否需要边界源项修正；
5. 是否需要 D2Q9；
6. 是否需要额外正则化。

---

# 28. 参考文献与源文件

1. `Multireflection boundary conditions for lattice Boltzmann models.pdf`  
   Ginzburg, I.; d’Humières, D.  
   DOI: `10.1103/PhysRevE.68.066614`

2. `LBM-CDE(1).pdf`  
   用于小黏度、小扩散率下通过源项调节物理输运系数的逆向设计模型。

3. `Lattice Boltzmann simulations of three-dimensional thermal convective flows at high Rayleigh number.pdf`  
   DOI: `10.1016/j.ijheatmasstransfer.2019.06.002`

4. `Accelerated lattice Boltzmann simulation using GPU and OpenACC with data management.pdf`  
   DOI: `10.1016/j.ijheatmasstransfer.2017.02.032`

5. `Lattice Boltzmann simulations of thermal convective flows in two dimensions(1).pdf`  
   DOI: `10.1016/j.camwa.2012.07.001`

6. `[Luo2014]_JCP_LB simulations of the thermally driven 2D square cavity at high Rayleigh numbers(1).pdf`  
   DOI: `10.1016/j.jcp.2014.06.047`

7. `Towards higher order lattice Boltzmann schemes.pdf`  
   DOI: `10.1088/1742-5468/2009/06/P06006`

8. `Optimal Stability of Advection-Diffusion Lattice Boltzmann Models with Two Relaxation Times for Positive_Negative Equilibrium.pdf`  
   DOI: `10.1007/s10955-010-9969-9`

---

# 29. Codex 完成推导：适用范围与最终结论

> 本节是对前述“任务 1--10”的完成稿。推导以
> `N=[1,1,1,1,1; ...]` 所定义的 D2Q5 矩顺序为准，并以仓库
> `均匀网格/2DRBOpenmp.F90`、`均匀网格/3DRBOpenmp.F90`
> 中的实际 D2Q9/D3Q19 矩顺序完成流场源项投影。

先给出结论：

1. 对“主分支 A”的无热源、恒定速度线性问题，局部温度梯度代回碰撞后，温度通量矩的更新**严格等价**于使用
   \[
   q_{\kappa,\mathrm{eff}}
   =
   \frac{2q_\kappa^{(0)}}
   {2(1-\chi_\kappa)+\chi_\kappa q_\kappa^{(0)}}.
   \]
   因而所有纯扩散四阶公式只需把原 \(\sigma_1\) 换成
   \(\sigma_{\kappa,\mathrm{eff}}\)，不是另一个独立的四阶模型。
2. D2Q5 四阶各向同性仍要求
   \[
   \sigma_{\kappa,\mathrm{eff}}\sigma_4=\frac16.
   \]
   四阶完全消除还要求原来的 \(\sigma_3\) 关系，只是其中
   \(\sigma_1\) 换成 \(\sigma_{\kappa,\mathrm{eff}}\)。
3. 当 \(\kappa\to0\) 且 \(c_T^2\) 固定时，
   \(q_4\sim6\sigma_{\kappa,\mathrm{eff}}\to0\)。因此
   “小扩散率 + 严格四阶各向同性 + 所有 ghost 松弛率远离 0”
   三者不能同时满足。这个结论与 \(q_\kappa^{(0)}\) 的选择无关。
4. 若放弃严格四阶条件，只保留二阶一致性并独立优化 \(q_3,q_4\)，
   D2Q5 仍存在 von Neumann 稳定参数区，因此目前**不能判定必须升级
   D2Q9**。
5. 对流场应力修正，在标准 halfway bounce-back、源项先并入
   post-collision 分布的约定下，精确壁面关系为
   \[
   \boxed{
   \sigma_{\nu,\mathrm{eff}}\sigma_q=\frac{3}{16}
   },
   \qquad
   \sigma_{\nu,\mathrm{eff}}
   =(1-\chi_s)\sigma_\nu^{(0)}.
   \]
   所以标准 halfway bounce-back 不能避免
   \(s_q\to0\)。若要避免，只能改变边界算子，而不是把
   \(3/16\) 强行施加在基础参数上。

以下推导使用两套源项记号：

- \(\boldsymbol R_c\)：连续时间源项，量纲为“每单位时间”；
- \(\widehat{\boldsymbol R}=\delta t\,\boldsymbol R_c\)：一个时间步的原始源增量；
- 真正加入碰撞的离散源为
  \[
  \boldsymbol\Psi
  =
  \left(I-\frac Q2\right)\widehat{\boldsymbol R}.
  \]

这样可避免把连续源和离散源混写。

---

# 30. 任务 1：D2Q5 源项的完整矩空间表达

## 30.1 晶格、权重与矩顺序

速度顺序固定为

\[
\boldsymbol e_0=(0,0),\quad
\boldsymbol e_1=(1,0),\quad
\boldsymbol e_2=(0,1),\quad
\boldsymbol e_3=(-1,0),\quad
\boldsymbol e_4=(0,-1).
\]

由前述平衡分布可得

\[
w_0=\frac{1-a}{5},
\qquad
w_1=w_2=w_3=w_4=\frac{4+a}{20}
=\frac{c_T^2}{2},
\]

\[
\sum_iw_i=1,
\qquad
\sum_iw_i e_{i\alpha}=0,
\qquad
\sum_iw_i e_{i\alpha}e_{i\beta}
=c_T^2\delta_{\alpha\beta}.
\]

矩顺序固定为

\[
\boldsymbol n
=
[T,j_x,j_y,e,\nu_T]^{\mathrm T}
=N\boldsymbol g.
\]

## 30.2 完整弱可压缩源与主分支 A

为与 `LBM-CDE.pdf` 的式 (24)、(28)、(35) 对齐，定义

\[
\boldsymbol A
=
\frac{p\nabla T+T\boldsymbol F}{\rho_0}.
\]

允许体热源 \(Q_T\) 时，一阶 Hermite 截断的 D2Q5 原始分布源可写成

\[
\boxed{
\widehat R_i
=
\delta t\,w_i
\left[
Q_T\left(
1+\frac{\boldsymbol e_i\cdot\boldsymbol u}{c_T^2}
\right)
+\frac{\boldsymbol e_i\cdot\boldsymbol A}{c_T^2}
+\chi_\kappa\boldsymbol e_i\cdot\nabla T
\right]
}.
\]

其矩空间投影为

\[
\boxed{
\widehat{\boldsymbol R}^{\,m}
=
\delta t
\begin{bmatrix}
Q_T\\
A_x+u_xQ_T+\chi_\kappa c_T^2\partial_xT\\
A_y+u_yQ_T+\chi_\kappa c_T^2\partial_yT\\
aQ_T\\
0
\end{bmatrix}
}.
\]

因此离散源的五个分量为

\[
\boxed{
\boldsymbol\Psi^m
=
\begin{bmatrix}
\delta t\,Q_T\\[2mm]
\left(1-\dfrac{q_\kappa^{(0)}}2\right)\widehat R_x^m\\[2mm]
\left(1-\dfrac{q_\kappa^{(0)}}2\right)\widehat R_y^m\\[2mm]
\left(1-\dfrac{q_3}2\right)a\delta t\,Q_T\\[2mm]
0
\end{bmatrix}
}.
\]

主分支 A 用于本次四阶推导的闭合选择是

\[
Q_T=0,\qquad
\boldsymbol A=\boldsymbol0,
\]

即只保留扩散修正：

\[
\boxed{
\widehat{\boldsymbol R}^{\,m}_{A}
=
\delta t
\begin{bmatrix}
0\\
\chi_\kappa c_T^2\partial_xT\\
\chi_\kappa c_T^2\partial_yT\\
0\\
0
\end{bmatrix}
}.
\]

相应分布空间源为

\[
\boxed{
\widehat{\boldsymbol R}_{A}
=
\frac{\chi_\kappa c_T^2\delta t}{2}
\begin{bmatrix}
0\\
\partial_xT\\
\partial_yT\\
-\partial_xT\\
-\partial_yT
\end{bmatrix}
}.
\]

守恒检查直接给出

\[
\sum_i\widehat R_{A,i}=0,
\]

所以扩散修正不改变 \(T=\sum_i g_i\)。

## 30.3 局部梯度代回后的精确有效松弛率

定义通量非平衡矩

\[
j_\alpha^{neq}
=
\sum_i e_{i\alpha}(g_i-g_i^{eq})
=j_\alpha-u_\alpha T.
\]

在主分支 A 中，局部梯度公式为

\[
\boxed{
\partial_\alpha T
=
-\frac{q_{\kappa,\mathrm{eff}}}
{c_T^2\delta t}
j_\alpha^{neq}
},
\]

其中

\[
\boxed{
q_{\kappa,\mathrm{eff}}
=
\frac{2q_\kappa^{(0)}}
{2(1-\chi_\kappa)+
\chi_\kappa q_\kappa^{(0)}}
}.
\]

证明如下。通量矩的碰撞为

\[
j_\alpha^*
=
j_\alpha
-q_\kappa^{(0)}j_\alpha^{neq}
+\left(1-\frac{q_\kappa^{(0)}}2\right)
\chi_\kappa c_T^2\delta t\,\partial_\alpha T.
\]

把局部梯度代入：

\[
j_\alpha^*
=
j_\alpha
-\left[
q_\kappa^{(0)}
+\left(1-\frac{q_\kappa^{(0)}}2\right)
\chi_\kappa q_{\kappa,\mathrm{eff}}
\right]j_\alpha^{neq}.
\]

而

\[
q_\kappa^{(0)}
+\left(1-\frac{q_\kappa^{(0)}}2\right)
\chi_\kappa q_{\kappa,\mathrm{eff}}
=q_{\kappa,\mathrm{eff}},
\]

故

\[
\boxed{
j_\alpha^*
=
j_\alpha
-q_{\kappa,\mathrm{eff}}
(j_\alpha-u_\alpha T)
}.
\]

因此这不是仅在 Chapman--Enskog 主阶成立的近似，而是在当前线性局部闭合下的逐格点代数恒等式。

对应 Hénon 参数为

\[
\boxed{
h
\equiv
\sigma_{\kappa,\mathrm{eff}}
=
\frac1{q_{\kappa,\mathrm{eff}}}-\frac12
=(1-\chi_\kappa)
\left(
\frac1{q_\kappa^{(0)}}-\frac12
\right)
}.
\]

于是

\[
\kappa=c_T^2h\delta t.
\]

这也说明：对 MRT D2Q5 而言，主分支 A 的 \(\chi_\kappa\) 在常系数线性体节点中不会产生一个额外独立的稳定参数；它的净作用就是把通量矩改成 \(q_{\kappa,\mathrm{eff}}\)。它的主要价值是允许固定 \(a\)，避免平衡权重随 \(\kappa\) 退化，而不是把物理温度通量模态任意强阻尼。

---

# 31. 任务 2：D2Q5 纯扩散四阶等效方程

令

\[
\boldsymbol u=\boldsymbol0,
\qquad
h=\sigma_{\kappa,\mathrm{eff}},
\qquad
r=\sigma_3,
\qquad
t=\sigma_4,
\]

\[
\lambda=\frac{\delta x}{\delta t}.
\]

对碰撞--传播矩阵作四阶 Taylor 展开，得到

\[
\boxed{
\begin{aligned}
\partial_tT
&-
\frac{\lambda^2\delta t}{10}
h(a+4)\nabla^2T\\
&+
\frac{\lambda^4\delta t^3}{1200}
h(a+4)
\left[
K_{40}
\left(
\partial_x^4T+\partial_y^4T
\right)
+K_{22}\partial_x^2\partial_y^2T
\right]
=O(\delta t^4),
\end{aligned}
}
\]

其中

\[
\boxed{
K_{40}
=
8-3a
+12(a+4)h^2
-12(1-a)hr
-60ht
},
\]

\[
\boxed{
K_{22}
=
-6(a+4)
+24(a+4)h^2
-24(1-a)hr
+120ht
}.
\]

因为

\[
c_T^2=\frac{a+4}{10},
\]

二阶项正好给出

\[
\kappa=\lambda^2\delta t\,c_T^2h.
\]

若把方程写成

\[
\partial_tT
=
\kappa\nabla^2T
+\widetilde\kappa_{40}
(\partial_x^4+\partial_y^4)T
+\widetilde\kappa_{22}
\partial_x^2\partial_y^2T
+O(\delta t^4),
\]

则有

\[
\boxed{
\widetilde\kappa_{40}
=
-\frac{\lambda^4\delta t^3}{1200}
h(a+4)K_{40}
},
\]

\[
\boxed{
\widetilde\kappa_{22}
=
-\frac{\lambda^4\delta t^3}{1200}
h(a+4)K_{22}
}.
\]

## 31.1 四阶各向同性

四阶各向同性要求

\[
K_{22}=2K_{40}.
\]

两式相减：

\[
K_{22}-2K_{40}
=
-40+240ht.
\]

所以

\[
\boxed{
ht=\frac16
},
\]

即

\[
\boxed{
\sigma_{\kappa,\mathrm{eff}}\sigma_4=\frac16
}.
\]

## 31.2 四阶完全消除

再要求

\[
K_{40}=K_{22}=0,
\]

得到

\[
\boxed{
t=\frac1{6h}
},
\]

\[
\boxed{
r
=
h\frac{a+4}{1-a}
-
\frac{2+3a}{12h(1-a)}
}.
\]

当

\[
\chi_\kappa=0
\]

时

\[
h=\sigma_\kappa^{(0)},
\]

上述公式严格退化为 `Towards higher order lattice Boltzmann schemes.pdf`
式 (40)--(42)、(55)。符号脚本
`tools/derive_d2q5_cde.py` 已逐项验证该退化。

## 31.3 小扩散率下的无解结论

四阶各向同性要求

\[
q_4
=
\frac1{t+1/2}
=
\boxed{
\frac{6h}{1+3h}
}.
\]

所以

\[
h\to0
\quad\Longrightarrow\quad
q_4\sim6h\to0.
\]

进一步考察 \(r\)：

\[
r
=
\frac{
h(a+4)
-(2+3a)/(12h)
}{1-a}.
\]

对固定 \(a<1\)：

- 若 \(a>-2/3\)，则 \(r\to-\infty\)，松弛率 \(q_3\) 不可接受；
- 若 \(a<-2/3\)，则 \(r\to+\infty\)，所以 \(q_3\to0\)；
- 若 \(a=-2/3\)，则 \(r=2h\to0\)，所以 \(q_3\to2\)。

因此在 \(h\to0\) 极限下，即使忽略 \(q_4\)，\(q_3\) 也会趋向不可接受区域或失去正 Hénon 参数。

结论是

\[
\boxed{
\text{D2Q5 在小扩散率极限不能同时保持严格四阶条件和有界 ghost 阻尼。}
}
\]

---

# 32. 任务 3：D2Q5 恒定对流等效方程

## 32.1 精确 Fourier 放大矩阵

令

\[
\boldsymbol u=(u,v)=\text{constant},
\]

并采用上节已经证明的有效通量松弛率

\[
q_1=q_{\kappa,\mathrm{eff}}.
\]

源项闭合后的矩空间碰撞矩阵为

\[
C=
\begin{bmatrix}
1&0&0&0&0\\
q_1u&1-q_1&0&0&0\\
q_1v&0&1-q_1&0&0\\
q_3a&0&0&1-q_3&0\\
0&0&0&0&1-q_4
\end{bmatrix}.
\]

对 Fourier 模式

\[
\boldsymbol g(\boldsymbol x,t)
=
\widehat{\boldsymbol g}(t)
\exp(i\boldsymbol k\cdot\boldsymbol x),
\]

一步放大矩阵为

\[
\boxed{
A(\boldsymbol k)
=
D(\boldsymbol k)N^{-1}CN
},
\]

\[
D(\boldsymbol k)
=
\operatorname{diag}
\left(
1,
e^{-ik_x},
e^{-ik_y},
e^{ik_x},
e^{ik_y}
\right).
\]

该矩阵既用于等效方程，也用于完整 von Neumann 稳定性搜索。

## 32.2 二阶数值扩散张量

设 \(A\) 的水动力特征值为 \(\Lambda_h(\boldsymbol k)\)，并定义

\[
\log\Lambda_h
=z_1+z_2+z_3+z_4+O(|\boldsymbol k|^5),
\]

其中 \(z_n\) 是 \(\boldsymbol k\) 的 \(n\) 次齐次多项式。

一阶项为

\[
z_1=-i(uk_x+vk_y).
\]

二阶项为

\[
\boxed{
z_2
=
-h
\left[
c_T^2(k_x^2+k_y^2)
-(uk_x+vk_y)^2
\right]
}.
\]

因此长波等效扩散张量为

\[
\boxed{
\boldsymbol D_{\mathrm{num}}
=
h\delta t
\left(
c_T^2I-\boldsymbol u\boldsymbol u^{\mathrm T}
\right)
}.
\]

展开后

\[
D_{xx}=h\delta t(c_T^2-u^2),
\]

\[
D_{yy}=h\delta t(c_T^2-v^2),
\]

\[
D_{xy}=D_{yx}=-h\delta t\,uv.
\]

所以 D2Q5 的交叉扩散项并未被 \(\chi_\kappa\) 独立消除；它只通过
\(h=\sigma_{\kappa,\mathrm{eff}}\) 缩放。长波必要条件仍为

\[
\boxed{
u^2+v^2\le c_T^2
}.
\]

## 32.3 三阶色散和四阶误差

令

\[
h=\sigma_{\kappa,\mathrm{eff}},
\qquad
r=\sigma_3,
\qquad
t=\sigma_4.
\]

三阶项写成

\[
z_3
=
i\left(
d_{30}k_x^3+d_{21}k_x^2k_y
+d_{12}k_xk_y^2+d_{03}k_y^3
\right),
\]

其中

\[
\begin{aligned}
d_{30}
=\frac{u}{120}\big[
&-24ah^2-12ahr+3a
+240h^2u^2-96h^2\\
&+12hr+60ht-20u^2+2
\big],
\end{aligned}
\]

\[
\begin{aligned}
d_{21}
=\frac{v}{40}\big[
&-8ah^2-4ahr+a
+240h^2u^2-32h^2\\
&+4hr-20ht-20u^2+4
\big],
\end{aligned}
\]

\[
\begin{aligned}
d_{12}
=\frac{u}{40}\big[
&-8ah^2-4ahr+a
+240h^2v^2-32h^2\\
&+4hr-20ht-20v^2+4
\big],
\end{aligned}
\]

\[
\begin{aligned}
d_{03}
=\frac{v}{120}\big[
&-24ah^2-12ahr+3a
+240h^2v^2-96h^2\\
&+12hr+60ht-20v^2+2
\big].
\end{aligned}
\]

四阶项写成

\[
z_4
=
e_{40}k_x^4+e_{31}k_x^3k_y
+e_{22}k_x^2k_y^2
+e_{13}k_xk_y^3+e_{04}k_y^4.
\]

五个 \(e_{ij}\) 的完整展开见本文件附录 A。它们满足以下复核：

\[
u=v=0
\quad\Longrightarrow\quad
e_{40}=e_{04}
=-\frac{h(a+4)}{1200}K_{40},
\]

\[
e_{22}
=-\frac{h(a+4)}{1200}K_{22},
\qquad
e_{31}=e_{13}=0.
\]

因此常对流展开与纯扩散四阶式完全相容。

在当前 \(Q_T=\boldsymbol F=\boldsymbol0\)、常 \(\boldsymbol u\) 的闭合中，
所有显式 \(\chi_\kappa\) 交叉项都已经代数合并进 \(h\)；
不存在另一个独立的“源项--对流”系数。若恢复
\(\boldsymbol A\)、\(Q_T\) 或空间变化速度，则本节矩阵 \(C\) 不再完整，
必须按式 (35) 的全源公式重新展开，不能把本节结果外推到该情形。

---

# 33. 任务 4：D2Q5 稳定参数区域

## 33.1 参数空间降维

给定目标扩散率后，

\[
h
=
\frac{\kappa}{c_T^2\delta t},
\qquad
q_{\kappa,\mathrm{eff}}
=
\frac1{h+1/2}.
\]

基础参数满足

\[
\sigma_\kappa^{(0)}
=
\frac1{q_\kappa^{(0)}}-\frac12,
\]

\[
\boxed{
\chi_\kappa
=
1-\frac{h}{\sigma_\kappa^{(0)}}
}.
\]

所以 \((q_\kappa^{(0)},\chi_\kappa)\) 在主分支 A 的线性体节点中并非两个独立稳定参数。固定 \(h\) 后，它们只决定如何把同一个有效碰撞拆成“基础碰撞 + 源项”。

完整可行集合定义为

\[
\mathcal P
=
\left\{
(a,q_\kappa^{(0)},\chi_\kappa,q_3,q_4):
\begin{array}{l}
-4<a<1,\\
q_i\in[\varepsilon,2-\varepsilon],\\
0\le\chi_\kappa<1,\\
\kappa=c_T^2(1-\chi_\kappa)
\sigma_\kappa^{(0)}\delta t,\\
\max_{\boldsymbol k\in[-\pi,\pi]^2}
\rho[A(\boldsymbol k)]\le1
\end{array}
\right\}.
\]

还必须同时检查

\[
\max_\alpha|u_\alpha|\le c_T^2
\]

的人口非负条件，以及

\[
|\boldsymbol u|^2\le c_T^2
\]

的长波必要条件；二者不能互相替代。

## 33.2 严格四阶子集的解析下界

若要求 \(q_4\ge\varepsilon\)，由

\[
q_4=\frac{6h}{1+3h}
\]

得到

\[
\boxed{
h
\ge
\frac{\varepsilon}{3(2-\varepsilon)}
}.
\]

相应最小扩散率为

\[
\boxed{
\kappa_{\min}^{(4)}
=
c_T^2\delta t
\frac{\varepsilon}{3(2-\varepsilon)}
}.
\]

例如

\[
\varepsilon=0.05
\]

时

\[
h_{\min}^{(4)}
\approx8.547\times10^{-3}.
\]

只要目标 \(h\) 低于该值，严格四阶各向同性子集就必为空，与
\(q_\kappa^{(0)}\) 和 \(\chi_\kappa\) 如何拆分无关。

还应检查

\[
\sigma_{\min}
=
\frac1{2-\varepsilon}-\frac12,
\qquad
\sigma_{\max}
=
\frac1\varepsilon-\frac12,
\]

\[
\sigma_{\min}\le r(a,h)\le\sigma_{\max}.
\]

这给出 \(a\) 的第二个可行性约束。

## 33.3 代表性 von Neumann 扫描

为判断“放弃严格四阶后 D2Q5 是否仍有稳定区域”，进行了如下可复现实验：

- \(a=0\)，故 \(c_T^2=0.4\)；
- \(|\boldsymbol u|=0.1\)，同时检查轴向和 \(45^\circ\) 对流；
- \(q_3,q_4=0.05,0.10,\ldots,1.95\)，共 \(39^2=1521\) 组；
- \(k_x,k_y\in[0,\pi]\)，每方向 25 点；
- 判据为所有特征值满足
  \(\rho[A(\boldsymbol k)]\le1+2\times10^{-10}\)。

结果如下：

| \(h\) | \(q_{\kappa,\mathrm{eff}}\) | 稳定组数 / 1521 | 稳定点的 \(q_3\) 外包络 | 稳定点的 \(q_4\) 外包络 |
|---:|---:|---:|---:|---:|
| \(10^{-1}\) | 1.666667 | 1521 | 0.05--1.95 | 0.05--1.95 |
| \(10^{-2}\) | 1.960784 | 1166 | 0.05--1.95 | 0.05--1.95 |
| \(10^{-3}\) | 1.996008 | 102 | 0.95--1.95 | 1.65--1.95 |

“外包络”只表示稳定点的最小/最大坐标，不表示包络中的每一对参数都稳定。

同一 \(a=0\) 下，严格四阶公式给出：

| \(h\) | \(q_3^{(4)}\) | \(q_4^{(4)}\) | 结论 |
|---:|---:|---:|---|
| \(10^{-1}\) | -1.30435 | 0.46154 | \(q_3\) 非法 |
| \(10^{-2}\) | -0.06201 | 0.05825 | \(q_3\) 非法 |
| \(10^{-3}\) | -0.00602 | 0.00598 | \(q_3\) 非法且 \(q_4\) 近 0 |

这组扫描支持两个不同结论：

1. 严格四阶 D2Q5 在小 \(h\) 下不可行；
2. 二阶 D2Q5 仍可能通过独立选择 ghost 松弛率获得稳定区。

该扫描不是高 \(Ra\) 实际算例的运行验证。正式参数必须把目标
\(Ra,Pr,Ma,N,\delta x,\delta t\) 映射到 \(h\)，再按实际速度上界重扫。

---

# 34. 任务 5：流场源项一致的 halfway bounce-back 推导

## 34.1 应力源使剪切矩成为有效碰撞

以 D2Q9 的 \(xy\) 剪切矩为例，忽略体力交叉项时

\[
\Pi_{xy}^{neq}
=
\sum_i c_{ix}c_{iy}
(f_i-f_i^{eq}).
\]

`LBM-CDE.pdf` 式 (31) 给出

\[
S_{xy}
=
-\frac{\Pi_{xy}^{neq}}
{\rho_0(2\nu+c_s^2\delta t)}.
\]

应力修正源投影到该矩的原始分量为

\[
\widehat R_{xy}^{\nu}
=
2\rho_0c_s^2\chi_sS_{xy}\delta t.
\]

设基础剪切松弛率为 \(s_\nu^{(0)}\)，则

\[
\Pi_{xy}^{neq,*}
=
\left(1-s_\nu^{(0)}\right)\Pi_{xy}^{neq}
+\left(1-\frac{s_\nu^{(0)}}2\right)
\widehat R_{xy}^{\nu}.
\]

利用

\[
\nu
=(1-\chi_s)c_s^2
\sigma_\nu^{(0)}\delta t,
\]

可化简为

\[
\boxed{
\Pi_{xy}^{neq,*}
=
\left(1-s_{\nu,\mathrm{eff}}\right)
\Pi_{xy}^{neq}
},
\]

\[
\boxed{
s_{\nu,\mathrm{eff}}
=
\frac1{
(1-\chi_s)\sigma_\nu^{(0)}+1/2
}
}.
\]

即应力修正源没有让物理剪切矩保持基础放大因子；它把该矩的实际放大因子变为有效值。

## 34.2 标准 halfway bounce-back 的闭合

标准 halfway bounce-back 只看到完整 post-collision 分布，因此其偶模态参数必须使用实际的

\[
\sigma_{\nu,\mathrm{eff}}
=
\frac1{s_{\nu,\mathrm{eff}}}-\frac12.
\]

对当前 D2Q9/TRT 矩顺序，原精确半格点关系遂变为

\[
\boxed{
\mathcal B
\left(
\sigma_\nu^{(0)},\chi_s,\sigma_q
\right)
=
(1-\chi_s)\sigma_\nu^{(0)}\sigma_q
-\frac{3}{16}
=0
}.
\]

无源极限检查：

\[
\chi_s=0
\quad\Longrightarrow\quad
\sigma_\nu^{(0)}\sigma_q=\frac3{16}.
\]

由闭合关系

\[
\sigma_q
=
\frac{3}{16\sigma_{\nu,\mathrm{eff}}},
\]

\[
\boxed{
s_q
=
\frac{16\sigma_{\nu,\mathrm{eff}}}
{3+8\sigma_{\nu,\mathrm{eff}}}
}.
\]

因此

\[
\nu\to0
\quad\Longrightarrow\quad
s_q\sim\frac{16}{3}\sigma_{\nu,\mathrm{eff}}\to0.
\]

最终判定是：

\[
\boxed{
\text{标准 halfway bounce-back 不能借助 }\chi_s
\text{ 避免 }s_q\to0.
}
\]

若改用

\[
\sigma_\nu^{(0)}\sigma_q=\frac3{16},
\]

则得到的是“基础碰撞的旧壁面关系”，而不是含源实际传播算子的闭合；它一般会产生参数相关壁面偏移或滑移。

要同时保持适中 \(s_q\) 与精确壁面，必须引入边界专用自由度，例如 multireflection、插值反弹或源项一致的 link correction。仅由体节点控制方程不能唯一确定这个额外边界项，因此本文件不把某个未经 Poiseuille 解析解验证的修正写成已证明公式。

---

# 35. 任务 6：实际 D2Q9/D3Q19 流场源项投影

## 35.1 统一定义

定义

\[
A_{\alpha\beta}
=
\chi_sS_{\alpha\beta}
+\frac{\chi_b-\chi_s}{D}
S_{\gamma\gamma}\delta_{\alpha\beta}.
\]

分布空间的应力原始源为

\[
\widehat R_i^\nu
=
\delta t\,w_i\rho_0
\left(
\frac{c_{i\alpha}c_{i\beta}}{c_s^2}
-\delta_{\alpha\beta}
\right)A_{\alpha\beta}.
\]

总原始矩源为

\[
\widehat{\boldsymbol R}^{\,m}
=M
\left(
\widehat{\boldsymbol R}^{\,F}
+\widehat{\boldsymbol R}^{\,\nu}
\right),
\]

真正加入碰撞的是

\[
\boxed{
\boldsymbol\Phi^m
=
\left(I-\frac S2\right)
\widehat{\boldsymbol R}^{\,m}
}.
\]

以下给出的都是乘半步因子之前的原始矩源，且省略公共的 \(\delta t\)。

## 35.2 D2Q9

仓库 D2Q9 矩顺序为

\[
[\rho,e,\epsilon,j_x,q_x,j_y,q_y,p_{xx},p_{xy}]^{\mathrm T}.
\]

令

\[
u_F=u_xF_x+u_yF_y.
\]

Guo/Boussinesq 力源投影为

\[
\widehat{\boldsymbol R}^{\,F}_{D2Q9}
=
\begin{bmatrix}
0\\
6u_F\\
-6u_F\\
F_x\\
-F_x\\
F_y\\
-F_y\\
2(u_xF_x-u_yF_y)\\
u_xF_y+u_yF_x
\end{bmatrix}.
\]

应力修正投影为

\[
\boxed{
\widehat{\boldsymbol R}^{\,\nu}_{D2Q9}
=
\rho_0
\begin{bmatrix}
0\\
2(A_{xx}+A_{yy})\\
-2(A_{xx}+A_{yy})\\
0\\
0\\
0\\
0\\
\dfrac23(A_{xx}-A_{yy})\\
\dfrac23A_{xy}
\end{bmatrix}
}.
\]

所以 Fortran 中每一行应实现

\[
\Phi_j
=
\left(1-\frac{s_j}{2}\right)
\left(
\widehat R_j^F+\widehat R_j^\nu
\right).
\]

不能把 BGK 的单一公共 prefactor 直接复制到 MRT 的全部矩。

## 35.3 D3Q19

仓库 D3Q19 矩顺序为

\[
\begin{aligned}
[
&\rho,e,\epsilon,
j_x,q_x,j_y,q_y,j_z,q_z,\\
&p_{xx},\pi_{xx},p_{ww},\pi_{ww},
p_{xy},p_{yz},p_{xz},
m_x,m_y,m_z
]^{\mathrm T}.
\end{aligned}
\]

令

\[
u_F=u_xF_x+u_yF_y+u_zF_z.
\]

力源投影为

\[
\widehat{\boldsymbol R}^{\,F}_{D3Q19}
=
\begin{bmatrix}
0\\
38u_F\\
-11u_F\\
F_x\\
-\frac23F_x\\
F_y\\
-\frac23F_y\\
F_z\\
-\frac23F_z\\
4u_xF_x-2u_yF_y-2u_zF_z\\
-2u_xF_x+u_yF_y+u_zF_z\\
2u_yF_y-2u_zF_z\\
-u_yF_y+u_zF_z\\
u_xF_y+u_yF_x\\
u_yF_z+u_zF_y\\
u_xF_z+u_zF_x\\
0\\
0\\
0
\end{bmatrix}.
\]

应力修正投影为

\[
\boxed{
\widehat{\boldsymbol R}^{\,\nu}_{D3Q19}
=
\rho_0
\begin{bmatrix}
0\\
\dfrac{38}{3}(A_{xx}+A_{yy}+A_{zz})\\
-\dfrac{11}{3}(A_{xx}+A_{yy}+A_{zz})\\
0\\
0\\
0\\
0\\
0\\
0\\
\dfrac23(2A_{xx}-A_{yy}-A_{zz})\\
-\dfrac13(2A_{xx}-A_{yy}-A_{zz})\\
\dfrac23(A_{yy}-A_{zz})\\
-\dfrac13(A_{yy}-A_{zz})\\
\dfrac23A_{xy}\\
\dfrac23A_{yz}\\
\dfrac23A_{xz}\\
0\\
0\\
0
\end{bmatrix}
}.
\]

五个独立偏应力分量对应矩 9、11、13、14、15；矩 10、12 是当前非正交基下与两个对角应力组合伴随的高阶矩，因此实际投影不为零。忽略 10、12 会破坏“先在分布空间构造源，再乘实际 \(M\)”的一致性。

---

# 36. 任务 7：应变率的局部表达式

定义

\[
\Pi_{\alpha\beta}^{neq}
=
\sum_i c_{i\alpha}c_{i\beta}
(f_i-f_i^{eq}),
\]

\[
\mu=\rho_0\nu,
\qquad
\mu_B=\rho_0\nu_B.
\]

对非对角分量 \(\alpha\ne\beta\)，`LBM-CDE.pdf` 式 (31) 给出

\[
\boxed{
S_{\alpha\beta}
=
-
\frac{
2\Pi_{\alpha\beta}^{neq}
+\delta t(u_\alpha F_\beta+u_\beta F_\alpha)
}{
4\mu+2\delta t\rho_0c_s^2
}
}.
\]

速度散度为

\[
\boxed{
S_{\gamma\gamma}
=
-
\frac{
\Pi_{\gamma\gamma}^{neq}
+\delta t\,u_\gamma F_\gamma
}{
D\mu_B+\delta t\rho_0c_s^2
}
}.
\]

对角分量为

\[
\boxed{
S_{\alpha\alpha}
=
-
\frac{
\Pi_{\alpha\alpha}^{neq}
+\delta t\,u_\alpha F_\alpha
}{
2\mu+\delta t\rho_0c_s^2
}
-
\frac{
2\mu-D\mu_B
}{
D(2\mu+\delta t\rho_0c_s^2)
}
S_{\gamma\gamma}
}.
\]

这里使用的是物理黏度

\[
\mu
=
\rho_0(1-\chi_s)c_s^2
\sigma_\nu^{(0)}\delta t,
\]

不是未修正的基础黏度。

这些公式已经把应力源对非平衡矩的贡献包含在分母中，所以计算顺序可以直接写成

1. 由当前 \(f_i\)、\(\rho\)、\(\boldsymbol u\) 计算
   \(\Pi_{\alpha\beta}^{neq}\)；
2. 用上式一次性得到 \(S_{\alpha\beta}\)；
3. 构造 \(A_{\alpha\beta}\)；
4. 投影并加入 \((I-S/2)\widehat{\boldsymbol R}^m\)。

不需要 fixed-point 迭代，也没有显式--隐式循环。

---

# 37. 任务 8：温度梯度的局部表达式

定义

\[
\boldsymbol j_T^{neq}
=
\sum_i\boldsymbol e_i(g_i-g_i^{eq}).
\]

若采用完整弱可压缩源

\[
\boldsymbol A
=
\frac{p\nabla T+T\boldsymbol F}{\rho_0}
\]

并允许 \(Q_T\ne0\)，由 `LBM-CDE.pdf` 式 (35) 得

\[
\boxed{
\partial_\alpha T
=
-
\frac{
2j_{T,\alpha}^{neq}
+\delta t
\left(
TF_\alpha/\rho_0+u_\alpha Q_T
\right)
}{
\delta t
\left\{
c_T^2
\left[
\dfrac{2(1-\chi_\kappa)}
{q_\kappa^{(0)}}
+\chi_\kappa
\right]
+p/\rho_0
\right\}
}
}.
\]

主分支 A 令

\[
p/\rho_0\ \text{耦合项}=0,
\qquad
\boldsymbol F=\boldsymbol0,
\qquad
Q_T=0,
\]

于是

\[
\boxed{
\nabla T
=
-\frac{q_{\kappa,\mathrm{eff}}}
{c_T^2\delta t}
\boldsymbol j_T^{neq}
}.
\]

Fortran 可直接实现为

```fortran
jneq_x = (g(1)-geq(1)) - (g(3)-geq(3))
jneq_y = (g(2)-geq(2)) - (g(4)-geq(4))
gradT_x = -qk_eff*jneq_x/(cT2*dt)
gradT_y = -qk_eff*jneq_y/(cT2*dt)
```

若使用完整源，则必须按完整分母和力/热源分子实现，不能把上面的主分支简式用于自然对流力耦合项已经显式加入温度源的版本。

---

# 38. 任务 9：边界条件的源项一致性

## 38.1 统一时间推进次序

采用

\[
\boxed{
\text{宏观量}
\to
\text{局部梯度/应变率}
\to
\text{碰撞并加入半步源}
\to
\text{传播}
\to
\text{link 边界}
}.
\]

记

\[
f_i^*
\]

和

\[
g_i^*
\]

为已经包含 \((I-S/2)\widehat R\) 的完整 post-collision 分布。边界式必须反射这个完整值。

## 38.2 流场 halfway bounce-back

对从壁面进入流体的未知方向 \(i\)：

\[
\boxed{
f_i(\boldsymbol x_b,t+\delta t)
=
f_{\bar i}^*(\boldsymbol x_b,t)
+2w_i\rho_0
\frac{\boldsymbol c_i\cdot\boldsymbol u_w}{c_s^2}
}.
\]

静止壁面 \(\boldsymbol u_w=\boldsymbol0\) 时就是完整 post-collision 反弹。

不需要再手工增加第二个“源项反弹项”；若碰撞例程先生成无源
\(f^{coll}\)，而源项在别处另加，则必须保证被反弹的是

\[
f_{\bar i}^{coll}+\Phi_{\bar i},
\]

否则时间离散不一致。

壁面精确位置条件已经在第 34 节得到：

\[
\sigma_{\nu,\mathrm{eff}}\sigma_q=\frac3{16}.
\]

## 38.3 定温 anti-bounce-back

对静止定温壁面 \(T_w\)：

\[
\boxed{
g_i(\boldsymbol x_b,t+\delta t)
=
-g_{\bar i}^*(\boldsymbol x_b,t)
+2w_iT_w
}.
\]

D2Q5 轴向权重为

\[
w_i=\frac{c_T^2}{2},
\]

所以

\[
\boxed{
g_i
=
-g_{\bar i}^*
+c_T^2T_w
}.
\]

这正是

\[
\frac{4+a}{10}T_w
\]

形式。

对一维稳态线性温度

\[
T_j=T_w+G(j-1/2),
\]

体节点解满足

\[
j_x^{neq}
=
-\frac{c_T^2}{q_{\kappa,\mathrm{eff}}}G.
\]

代入 anti-bounce-back 可严格得到

\[
T_1-T_w=\frac G2,
\]

所以线性纯导热的壁面仍在半格点，且不单独依赖
\(\chi_\kappa\) 或 \(q_\kappa^{(0)}\)，只通过已合并的
\(q_{\kappa,\mathrm{eff}}\) 进入体节点分布。

该结论不等价于“任意曲率温度场都无边界误差”。当
\(\partial_n^2T\ne0\)、\(Q_T\ne0\) 或使用完整
\(p\nabla T+T\boldsymbol F\) 源时，仍需用 manufactured solution
测量温度滑移。

## 38.4 绝热 bounce-back

绝热壁面采用

\[
\boxed{
g_i(\boldsymbol x_b,t+\delta t)
=
g_{\bar i}^*(\boldsymbol x_b,t)
}.
\]

在一维稳态情形代入通量矩更新可得

\[
(2-q_{\kappa,\mathrm{eff}})j_n^{neq}=0.
\]

只要

\[
q_{\kappa,\mathrm{eff}}\ne2,
\]

即得到

\[
j_n^{neq}=0
\quad\Longleftrightarrow\quad
\partial_nT=0.
\]

当 \(q_{\kappa,\mathrm{eff}}\to2\) 时，该约束退化，这与小扩散率下通量矩弱阻尼是同一个极限问题。

## 38.5 Nusselt 数

物理定义不变：

\[
Nu_w
=
-\frac{L}{\Delta T}
\left.\frac{\partial T}{\partial n}\right|_w.
\]

但后处理必须与碰撞使用同一个局部梯度公式。不能在碰撞中使用
\(q_{\kappa,\mathrm{eff}}\) 的非平衡梯度，却在壁面 \(Nu\) 中继续使用旧
\(q_\kappa^{(0)}\) 系数。

---

# 39. 任务 10：是否升级 D2Q9 的最终判定

本轮推导得到的是条件判定，而不是无条件升级：

\[
\boxed{
\begin{array}{l}
\text{若目标是稳定的二阶高 }Ra\text{ 算法，先保留 D2Q5；}\\
\text{若必须同时满足严格四阶各向同性且 }
h<h_{\min}^{(4)},\text{ 则升级 D2Q9。}
\end{array}
}
\]

D2Q9 的启动条件为以下任一项：

1. 目标 \(h\) 低于
   \[
   \varepsilon/[3(2-\varepsilon)]
   \]
   且仍要求严格四阶各向同性；
2. D2Q5 的交叉扩散
   \[
   -h\,uv\,\partial_{xy}T
   \]
   在对角对流测试中不可接受；
3. 完整谱扫描在目标速度范围内不存在满足
   \(q_3,q_4\in[\varepsilon,2-\varepsilon]\) 的稳定点；
4. ABB/BB 边界误差在目标参数下无法通过网格收敛控制。

D2Q9 回退分支应至少采用能表示完整二阶速度张量的平衡分布：

\[
g_i^{eq}
=
w_iT
\left[
1
+\frac{\boldsymbol c_i\cdot\boldsymbol u}{c_T^2}
+\frac{
(\boldsymbol c_i\cdot\boldsymbol u)^2
-c_T^2|\boldsymbol u|^2
}{
2c_T^4
}
\right],
\]

从而局部表示 \(u_xu_y\) 并消除 D2Q5 缺失对角速度导致的结构性限制。

是否采用该分支必须由第 25 节的恒定速度对流扩散和自然对流方腔验证决定；当前推导没有证明 D2Q5 的二阶稳定分支失败，因此第一版实现仍应保留 D2Q5。

---

# 40. 可直接转写为 Fortran 的主分支 A

## 40.1 参数预处理

```fortran
cT2    = (4.0d0 + paraA)/10.0d0
sigma0 = 1.0d0/qk0 - 0.5d0
hEff   = diffusivity/(cT2*dt)
chiK   = 1.0d0 - hEff/sigma0
qkEff  = 1.0d0/(hEff + 0.5d0)
```

必须检查

```fortran
if (cT2 <= 0.0d0) error stop "cT2 must be positive"
if (qk0 <= 0.0d0 .or. qk0 >= 2.0d0) error stop "invalid qk0"
if (chiK < 0.0d0 .or. chiK >= 1.0d0) error stop "invalid chiK"
```

## 40.2 逐格点碰撞

```fortran
! n = N*g
n0 = g0 + g1 + g2 + g3 + g4
n1 = g1 - g3
n2 = g2 - g4
n3 = -4.0d0*g0 + g1 + g2 + g3 + g4
n4 = g1 - g2 + g3 - g4

neq0 = T
neq1 = T*u
neq2 = T*v
neq3 = paraA*T
neq4 = 0.0d0

jneqX = n1 - neq1
jneqY = n2 - neq2
gradTx = -qkEff*jneqX/(cT2*dt)
gradTy = -qkEff*jneqY/(cT2*dt)

r1 = chiK*cT2*dt*gradTx
r2 = chiK*cT2*dt*gradTy

npost0 = n0
npost1 = n1 - qk0*(n1-neq1) + (1.0d0-0.5d0*qk0)*r1
npost2 = n2 - qk0*(n2-neq2) + (1.0d0-0.5d0*qk0)*r2
npost3 = n3 - q3*(n3-neq3)
npost4 = n4 - q4*(n4-neq4)
```

上述 `npost1/npost2` 可用

```fortran
npost1 = n1 - qkEff*(n1-neq1)
npost2 = n2 - qkEff*(n2-neq2)
```

逐点验证；两者在主分支 A 中应只相差舍入误差。这个对照是最便宜的实现自检。

## 40.3 物理与数值不变量

实现后必须逐项检查：

\[
\sum_i\Psi_i=0,
\]

\[
\chi_\kappa=0
\Rightarrow
q_{\kappa,\mathrm{eff}}=q_\kappa^{(0)},
\]

\[
\boldsymbol u=0
\Rightarrow
\partial_tT=\kappa\nabla^2T+O(\delta x^4),
\]

\[
\kappa\to0
\Rightarrow
q_{\kappa,\mathrm{eff}}\to2,
\]

以及

\[
\max_\alpha|u_\alpha|\le c_T^2,
\qquad
|\boldsymbol u|^2\le c_T^2.
\]

---

# 41. 完成状态、验证边界与文献依据

## 41.1 已由代数和符号计算闭合

- D2Q5 源项的分布/矩空间分量；
- 半步源离散；
- 局部梯度与 \(q_{\kappa,\mathrm{eff}}\) 的代数等价；
- 纯扩散四阶 \(K_{40},K_{22}\)；
- \(\chi_\kappa=0\) 退化；
- 四阶各向同性和完全消除条件；
- 常对流二至四阶 Fourier 展开；
- D2Q9/D3Q19 实际矩阵的完整源项投影；
- 应变率和温度梯度局部表达；
- 标准 HBB 下的有效 \(3/16\) 关系；
- 线性纯导热 ABB 和绝热 BB 的半格点/零通量闭合。

## 41.2 尚需数值验证，不能写成已验证性能

- 非线性自然对流中完整 \(p\nabla T+T\boldsymbol F\) 温度源的四阶误差；
- 曲率温度场的 ABB 隐藏误差；
- 高 \(Ra\) 实际速度范围内的全谱稳定窗口；
- 新流场源投影在 D3Q19 上的 Poiseuille/Couette 壁面误差；
- OpenACC/MPI 实现的设备正确性和运行稳定性；
- \(Nu\)、\(Re\) 与基准数据的一致性。

本轮没有修改求解器，也没有完成 GPU/高 \(Ra\) 运行。因此这里的“完成”指推导任务已经闭合，并不等于新算法已经通过算例验证。

## 41.3 主要原始文献

1. `pdf/LBM-CDE.pdf`：式 (24)--(29) 的源项与输运系数映射；式 (31)--(35) 的局部应变率/梯度；式 (36)--(39) 的边界形式。
2. `pdf/Towards higher order lattice Boltzmann schemes.pdf`：D2Q5 纯扩散式 (40)--(42) 与四阶参数式 (55)。
3. `pdf/Multireflection boundary conditions for lattice Boltzmann models.pdf`：link 边界、偶/奇碰撞参数与壁面位置的关系。
4. `pdf/Optimal Stability of Advection-Diffusion Lattice Boltzmann Models with Two Relaxation Times for Positive_Negative Equilibrium.pdf`：人口正负、长波条件和全谱稳定性之间的区别。
5. `均匀网格/2DRBOpenmp.F90`：本推导采用的实际 D2Q9/D2Q5 矩顺序。
6. `均匀网格/3DRBOpenmp.F90`：本推导采用的实际 D3Q19 矩顺序。
7. `文档/high_ra_derivation_source_matrix.md`：逐任务区分文献直接支持、本文新推导和仍需数值验证的结论。

---

# 附录 A：恒定对流四阶系数完整展开

本附录采用

\[
h=\sigma_{\kappa,\mathrm{eff}},
\qquad
r=\sigma_3,
\qquad
t=\sigma_4.
\]

## A.1 \(e_{40}\)

\[
\begin{aligned}
e_{40}=\frac1{1200}\big[
&-12a^2h^3-12a^2h^2r+3a^2h
+720ah^3u^2-96ah^3\\
&+360ah^2ru^2-36ah^2r+60ah^2t
+120ahr^2u^2\\
&-150ahu^2+4ah-30aru^2
-6000h^3u^4+2880h^3u^2\\
&-192h^3-360h^2ru^2+48h^2r
-1800h^2tu^2+240h^2t\\
&-120hr^2u^2-600ht^2u^2
+900hu^4-100hu^2-32h\\
&+30ru^2+150tu^2
\big].
\end{aligned}
\]

## A.2 \(e_{31}\)

\[
\begin{aligned}
e_{31}
=\frac{uv}{60}\big[
&72ah^3+36ah^2r+12ahr^2-15ah-3ar\\
&-1200h^3u^2+288h^3-36h^2r-12hr^2\\
&+180hu^2-35h+3r
\big].
\end{aligned}
\]

## A.3 \(e_{22}\)

\[
\begin{aligned}
e_{22}=\frac1{200}\big[
&-4a^2h^3-4a^2h^2r+a^2h\\
&+120ah^3(u^2+v^2)-32ah^3
+60ah^2r(u^2+v^2)-12ah^2r\\
&-20ah^2t
+20ahr^2(u^2+v^2)
-25ah(u^2+v^2)+8ah\\
&-5ar(u^2+v^2)
-6000h^3u^2v^2
+480h^3(u^2+v^2)-64h^3\\
&-60h^2r(u^2+v^2)+16h^2r
+300h^2t(u^2+v^2)-80h^2t\\
&-20hr^2(u^2+v^2)
+100ht^2(u^2+v^2)
+900hu^2v^2\\
&-100h(u^2+v^2)+16h
+5r(u^2+v^2)-25t(u^2+v^2)
\big].
\end{aligned}
\]

## A.4 \(e_{13}\)

\[
\begin{aligned}
e_{13}
=\frac{uv}{60}\big[
&72ah^3+36ah^2r+12ahr^2-15ah-3ar\\
&-1200h^3v^2+288h^3-36h^2r-12hr^2\\
&+180hv^2-35h+3r
\big].
\end{aligned}
\]

## A.5 \(e_{04}\)

\[
e_{04}=e_{40}\big|_{u\rightarrow v}.
\]

这些系数由

\[
\log\Lambda_h
=
-i(uk_x+vk_y)
-h[c_T^2|\boldsymbol k|^2-(\boldsymbol u\cdot\boldsymbol k)^2]
+z_3+z_4+O(|\boldsymbol k|^5)
\]

定义。符号生成和纯扩散退化断言位于
`tools/derive_d2q5_cde.py`。
