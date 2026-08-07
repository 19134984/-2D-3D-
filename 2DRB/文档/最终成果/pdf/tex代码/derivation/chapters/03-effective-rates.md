# 局部反馈消元与物理模态有效率

本章把 Task 2 的源项代回同一个 transformed TRT 碰撞矩，不从论文的输运系数公式反猜松弛率。结论分成四级：原论文直接给出的闭合式、项目内精确代数、冻结系数局部模态解释、以及明确不作的高阶或边界主张。

## 四种不可混用的量

对任一名义松弛率 $s$，定义 Hénon 平移

$$
\sigma=\frac1s-\frac12,
\qquad
1-\frac{s}{2}=s\sigma.
$$

算法输入仍是 Task 2 的四个名义率

$$
(s_f^+,s_f^-),\qquad (s_g^+,s_g^-).
$$

局部应变或温度梯度反馈消元后，某个物理矩块可表现为

$$
\sigma_{\rm eff},\qquad
s_{\rm eff}=\frac1{\sigma_{\rm eff}+1/2}.
$$

$s_{\rm eff}$ 是局部闭合后的物理模态系数，不是第五个或第六个 TRT 输入。输运系数还要乘格子尺度与矩归一化；例如剪切黏度对应 $c_s^2\Delta t\,\sigma_{s,\rm eff}$，而 $D$ 维体黏度多一个 $2/D$。

## 流场二阶矩的直接代回

记 $s=s_f^+$、$\sigma=1/s-1/2$，并把三个物理块统一写成

$$
\mathcal P^*=\mathcal P
-s\mathcal P
+\Delta t\left(1-\frac{s}{2}\right)
\left(a+2\rho_0c_s^2\chi\,\mathcal S\right).
$$

这里 $\mathcal P$ 是 transformed 二阶非平衡矩，$a$ 是 $uF$ 仿射项，$\mathcal S$ 是相应应变分量。原始 PDF Eq. (30) 及 Eqs. (31)-(33) 可统一改写为

$$
\mathcal S
=-\frac{2\mathcal P+\Delta t\,a}
{2\Delta t\rho_0c_s^2[1+2(1-\chi)\sigma]}.
$$

三个块的具体选择是：

| 物理块 | $\mathcal P$ | $a$ | $\mathcal S$ | $\chi$ |
| --- | --- | --- | --- | --- |
| 非对角剪切 | $\Pi^{neq}_{\alpha\beta}$，$\alpha\ne\beta$ | $u_\alpha F_\beta+u_\beta F_\alpha$ | $S_{\alpha\beta}$ | $\chi_s$ |
| 对角偏差 | $\Pi^{neq}_{xx}-\Pi^{neq}_{yy}$ | $2(u_xF_x-u_yF_y)$ | $S_{xx}-S_{yy}$ | $\chi_s$ |
| 迹/体积 | $\Pi^{neq}_{\gamma\gamma}$ | $2u_\gamma F_\gamma$ | $S_{\gamma\gamma}$ | $\chi_b$ |

迹块中使用了

$$
D\mu^B=2\rho_0c_s^2\Delta t(1-\chi_b)\sigma,
$$

所以它与剪切块具有相同的代数分母，但物理体黏度仍保留论文的 $2/D$ 归一化。

把闭合式代回碰撞增量，并用 $1-s/2=s\sigma$ 化简，可得

$$
\Delta\mathcal P
=-s_{\rm eff}\mathcal P
+\Delta t\left(1-\frac{s_{\rm eff}}2\right)a,
$$

其中

$$
\boxed{\sigma_{s,\rm eff}=(1-\chi_s)\sigma_f^+},
\qquad
\boxed{\sigma_{b,\rm eff}=(1-\chi_b)\sigma_f^+}.
$$

因此

$$
\boxed{\nu=c_s^2\Delta t\,\sigma_{s,\rm eff}},
\qquad
\boxed{\nu^B=\frac{2}{D}c_s^2\Delta t\,\sigma_{b,\rm eff}}.
$$

这一步同时说明 $u_\alpha F_\beta+u_\beta F_\alpha$ 没有被吸收入率；它仍以有效物理块自己的梯形因子 $1-s_{\rm eff}/2$ 出现。

## 标量奇通量的直接代回

令

$$
\pi=\frac{p}{\rho_0},\qquad
\boldsymbol a_T=\frac{T\boldsymbol F}{\rho_0}+Q\boldsymbol u,
$$

并记 transformed 标量奇通量非平衡矩为 $\boldsymbol j^{neq}$。Task 2 的 Eq. (24) 源矩是

$$
\boldsymbol J
=(\pi+\chi_\kappa c_s^2)\nabla T+\boldsymbol a_T.
$$

原始 PDF Eq. (35) 给出局部闭合

$$
\nabla T
=-\frac{2\boldsymbol j^{neq}+\Delta t\boldsymbol a_T}
{2\kappa+\Delta t(\pi+c_s^2)},
$$

其中论文 Eq. (29) 对名义奇率的 TRT 推广写成

$$
\kappa=(1-\chi_\kappa)c_s^2\Delta t\,\sigma_g^-.
$$

标量奇矩的实际碰撞增量为

$$
\Delta\boldsymbol j
=-s_g^-\boldsymbol j^{neq}
+\Delta t\left(1-\frac{s_g^-}{2}\right)
\left[(\pi+\chi_\kappa c_s^2)\nabla T+\boldsymbol a_T\right].
$$

设

$$
A=c_s^2+\pi,\qquad B=(1-\chi_\kappa)c_s^2.
$$

直接代回后，$\boldsymbol j^{neq}$ 的系数为

$$
-s_g^-\left[1+\frac{2\sigma_g^-(A-B)}{A+2B\sigma_g^-}\right]
=-\frac{2A}{A+2B\sigma_g^-},
$$

仿射项系数为

$$
\Delta t\frac{2B\sigma_g^-}{A+2B\sigma_g^-}.
$$

两者正好组成

$$
\Delta\boldsymbol j
=-s_{\rm flux,eff}\boldsymbol j^{neq}
+\Delta t\left(1-\frac{s_{\rm flux,eff}}2\right)\boldsymbol a_T,
$$

且

$$
\boxed{
\sigma_{\rm flux,eff}
=\frac{(1-\chi_\kappa)c_s^2}{c_s^2+\pi}\sigma_g^-,
\qquad
s_{\rm flux,eff}
=\frac1{\sigma_{\rm flux,eff}+1/2}.
}
$$

冻结压力下物理扩散系数中的压力严格消去：

$$
\boxed{
\kappa=(c_s^2+\pi)\Delta t\,\sigma_{\rm flux,eff}
=(1-\chi_\kappa)c_s^2\Delta t\,\sigma_g^-.
}
$$

$T\boldsymbol F/\rho_0$ 与 $Q\boldsymbol u$ 始终保留为仿射通量；它们不是 $\sigma_{\rm flux,eff}$ 的一部分。

## TRT 骨架内的受限块 MRT

反馈只作用于特定物理矩。对同一奇偶子空间内的 ghost，源矩没有相应闭合，因此仍使用名义率。

| 场/奇偶块 | 物理或 ghost 模式 | 使用的 Hénon 平移 |
| --- | --- | --- |
| flow even | 剪切、对角偏差 | $(1-\chi_s)\sigma_f^+$ |
| flow even | 迹/体积 | $(1-\chi_b)\sigma_f^+$ |
| flow even | 其余 even ghost | 名义 $\sigma_f^+$ |
| flow odd | odd ghost | 名义 $\sigma_f^-$ |
| scalar odd | 一阶物理通量 | $\dfrac{(1-\chi_\kappa)c_s^2}{c_s^2+\pi}\sigma_g^-$ |
| scalar odd | 其余 odd ghost | 名义 $\sigma_g^-$ |
| scalar even | 非守恒 even ghost | 名义 $\sigma_g^+$ |

所以完整算法仍是 TRT 奇偶骨架，但局部消元后的物理解释已经在同一 parity 内分裂；准确说法是“源反馈诱导的受限块 MRT 有效算子”。这不意味着生产碰撞矩阵新增了可独立输入的 MRT 参数。

## 极限与结论等级

- 当 $\chi_s=\chi_b=\chi_\kappa=0$ 且 $\pi=0$ 时，各物理移位回到相应名义移位；若再令四个名义率相等，则回到 Task 2 的逐分量 transformed BGK 极限。
- 固定名义 $\sigma=1/2$（即 $s=1$）并令 $\chi\to1^-$ 时，物理输运趋于零，而名义碰撞率不必趋于 $2$；局部物理块的 $s_{\rm eff}$ 才趋于 $2$。
- $\chi\le1$ 是论文用来避免负输运系数的物理限制；数值稳定性仍需由后续求解器和边界验证，不由局部率公式自动保证。

以上剪切、偏差、迹和标量通量等式属于“原论文闭合式 + 项目内精确代回”的代数结论。标量率中的 $\pi$ 解释只属于冻结系数局部模态。若在参考点 $\pi_0$ 冻结该率而实际压力为 $\pi_0+\delta\pi(\boldsymbol x)$，首个乘积导数残差是

$$
\Delta t\,\sigma_{\rm flux,eff}(\pi_0)
\left[\nabla\delta\pi\cdot\nabla T
+\delta\pi\,\nabla^2T\right],
$$

一般不为零。本章不把局部率式宣称为变系数四阶等效方程、边界 magic parameter 或普适稳定性判据。

在第四章采用的 $\nabla=\epsilon\nabla_1$ 标度下，若 $\delta\pi$ 是慢变量上的 $O(1)$ 系数变化，则 $\nabla\delta\pi\cdot\nabla T$ 与 $\delta\pi\nabla^2T$ 都是 $O(\epsilon^2)$。它们不是常系数模型的 $O(\epsilon^3)$ 首个遗漏项。要保持所述二阶目标方程，必须把 $\pi$ 在 CE2 内视为常数/冻结系数，或把该乘积导数显式加入目标方程；仅有空间光滑性不会提高它的 CE 阶数。
