# 04 梯形 transformed LBE 的二阶宏观恢复

本章从 Task 1 的平衡矩和 Task 2 的源矩出发，对 transformed 离散 LBE 作二阶 Chapman--Enskog（CE）展开。目标方程只用于最后比对，不作为符号残差的输入夹具。

## 1. 两套互不替代的标度

### 1.1 CE 慢变量标度

采用标准两时间展开

$$
\nabla=\epsilon\nabla_1,
\qquad
\partial_t=\epsilon\partial_{t_1}
+\epsilon^2\partial_{t_2},
$$

$$
\widetilde h=h^{(0)}+\epsilon h^{(1)}+\epsilon^2h^{(2)}+O(\epsilon^3),
\qquad
q=\epsilon q^{(1)}+\epsilon^2q^{(2)}+O(\epsilon^3).
$$

$\boldsymbol F$、$Q$ 和 $\nabla T$ 都按慢尺度进入相应源矩；这只说明它们在 CE 空间/时间展开中的位置，不决定 Mach 幅值。

### 1.2 低 Mach 幅值标度

为与原始 PDF Eqs. (11)、(22) 的删项等级一致，另用

$$
\boldsymbol u=O(\mathrm{Ma}),\qquad
\delta\rho/\rho_0=O(\mathrm{Ma}^2),
$$

$$
p=p_0+O(\mathrm{Ma}^2),\qquad
\boldsymbol F^{(1)}=O(\mathrm{Ma}),\qquad
Q^{(1)}=O(\mathrm{Ma}^0).
$$

于是 $uF$ 是 $O(\mathrm{Ma}^2)$，$TF/\rho_0$ 与 $Qu$ 是 $O(\mathrm{Ma})$，而 $p_0\nabla T/\rho_0$ 可在标量通量的 $O(\mathrm{Ma}^0)$ 层保留。这里 $p_0$ 的梯度为零；空间变化部分才产生后文单列的变系数残差。

论文 Eq. (11) 中标记为 $O(\mathrm{Ma}^3)$ 及更高的流场项、Eq. (22) 中的 $-Tu_\alpha\partial_\beta u_\beta=O(\mathrm{Ma}^3)$ 不参与二阶目标闭合，但会在第 8 节明确列出。

## 2. transformed 离散方程与三层人口方程

对流场或标量场统一记

$$
\Lambda=s_+P_++s_-P_-,
\qquad
\Sigma=\Lambda^{-1}-\frac12I
=\sigma_+P_++\sigma_-P_-.
$$

Task 2 的碰撞--迁移方程是

$$
\widetilde h_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
-\widetilde h_i(\boldsymbol x,t)
=-\left[\Lambda(\widetilde h-h^{eq})\right]_i
+\Delta t\left[\left(I-\frac{\Lambda}{2}\right)q\right]_i.
$$

定义

$$
D_{1i}=\partial_{t_1}+c_{i\alpha}\partial_{\alpha1}.
$$

对左端 Taylor 展开至 $O(\epsilon^2)$，逐阶得到

$$
\epsilon^0:\qquad
h^{(0)}=h^{eq},
$$

$$
\epsilon^1:\qquad
\Delta t D_1 h^{(0)}
=-\Lambda h^{(1)}
+\Delta t\left(I-\frac{\Lambda}{2}\right)q^{(1)},
$$

$$
\epsilon^2:\qquad
\Delta t\partial_{t_2}h^{(0)}
+\Delta t D_1h^{(1)}
+\frac{\Delta t^2}{2}D_1^2h^{(0)}
=-\Lambda h^{(2)}
+\Delta t\left(I-\frac{\Lambda}{2}\right)q^{(2)}.
$$

这里每个源块自动使用本 parity 的因子：flow even/odd 分别是 $1-s_f^+/2$、$1-s_f^-/2$，scalar even/odd 分别是 $1-s_g^+/2$、$1-s_g^-/2$。

## 3. 梯形半源对时间离散项的消去

由一阶人口方程

$$
h^{(1)}
=-\Delta t\Lambda^{-1}D_1h^{(0)}
+\Delta t\left(\Lambda^{-1}-\frac12I\right)q^{(1)},
$$

所以

$$
\boxed{
h^{(1)}+\frac{\Delta t}{2}D_1h^{(0)}
=-\Delta t\Sigma\left(D_1h^{(0)}-q^{(1)}\right).
}
$$

二阶方程中的后两项恰好组合为

$$
D_1h^{(1)}+\frac{\Delta t}{2}D_1^2h^{(0)}
=-\Delta t D_1\Sigma(D_1h^{(0)}-q^{(1)}).
$$

对守恒矩取矩时，Task 2 的

$$
\rho_0\boldsymbol u
=\sum_i\boldsymbol c_i\widetilde f_i+\frac{\Delta t}{2}\boldsymbol F,
\qquad
T=\sum_i\widetilde g_i+\frac{\Delta t}{2}Q
$$

又使 $h^{(2)}$ 的半源守恒矩与右端 $q^{(2)}$ 的 $-\Lambda/2$ 部分配对。因而守恒矩的二阶通式是

$$
\boxed{
\partial_{t_2}\langle\psi h^{(0)}\rangle
-\Delta t\left\langle
\psi D_1\Sigma(D_1h^{(0)}-q^{(1)})
\right\rangle
=\langle\psi q^{(2)}\rangle.
}
$$

若把某个显式源因子写成任意 $b$，它进入构成关系的归一化系数是

$$
r=\frac{b}{s\sigma}.
$$

正确梯形因子满足 $b=1-s/2=s\sigma$，故 $r=1$。例如把 $b$ 错换成 $1$，就得到 $r=2/(2-s)$，不再能消去相应源矩；这也是后续扰动测试的代数入口。

对守恒力/热源还需计入 transformed half-source nonequilibrium 的松弛贡献。其 Euler 层净系数不是 $r$，而是

$$
c=b+\frac{s}{2}.
$$

正确 $b$ 同样给出 $c=1$。前者控制 even stress 或 odd scalar flux 的非守恒构成矩，后者控制 $\boldsymbol F$、$Q$ 进入守恒方程的净幅值。

## 4. 连续性与动量方程

Task 1/原论文 Eq. (6) 的流场平衡矩与 Task 2 的源矩是

$$
\sum_i f_i^{eq}=\delta\rho,
\qquad
\sum_i c_{i\alpha}f_i^{eq}=\rho_0u_\alpha,
$$

$$
\sum_i c_{i\alpha}c_{i\beta}f_i^{eq}
=\delta\rho c_s^2\delta_{\alpha\beta}+\rho_0u_\alpha u_\beta,
$$

$$
\sum_i S_i=0,
\qquad
\sum_i c_{i\alpha}S_i=F_\alpha,
$$

$$
\sum_i c_{i\alpha}c_{i\beta}S_i
=u_\alpha F_\beta+u_\beta F_\alpha
+2\rho_0c_s^2
\left[\chi_sS_{\alpha\beta}
+\frac{\chi_b-\chi_s}{D}S_{\gamma\gamma}\delta_{\alpha\beta}
\right].
$$

零阶矩在 $O(\epsilon)$ 给出

$$
\partial_{t_1}\delta\rho
+\partial_{\alpha1}(\rho_0u_\alpha)=0.
$$

一阶矩在 $O(\epsilon)$ 给出 Euler 层

$$
\partial_{t_1}(\rho_0u_\alpha)
+\partial_{\beta1}
(\rho_0u_\alpha u_\beta+p\delta_{\alpha\beta})
=F_\alpha,
$$

其中 $p$ 的常数基准部分不贡献梯度。把该 Euler 层和连续性层代入二阶非平衡应力核，平衡矩时间导数产生的
$u_\alpha F_\beta+u_\beta F_\alpha$ 与 even 源二阶矩中的同名项相消。保留到论文 Eq. (11) 的 $O(\mathrm{Ma}^2)$ 后，三个 even 物理块分别得到

$$
\Pi^{(1),dev}_{\alpha\beta}
=-2\rho_0c_s^2\Delta t
(1-\chi_s)\sigma_f^+S^{dev}_{\alpha\beta},
$$

$$
\Pi^{(1),trace}
=-2\rho_0c_s^2\Delta t
(1-\chi_b)\sigma_f^+S_{\gamma\gamma}.
$$

合并 $O(\epsilon)$ 与 $O(\epsilon^2)$ 后恢复

$$
\partial_t\delta\rho+\partial_\alpha(\rho_0u_\alpha)
=O(\epsilon^3),
$$

$$
\partial_t(\rho_0u_\alpha)
+\partial_\beta(\rho_0u_\alpha u_\beta)
=-\partial_\alpha p+\partial_\beta\tau_{\alpha\beta}+F_\alpha
+O(\epsilon^3)+O(\mathrm{Ma}^3),
$$

$$
\tau_{\alpha\beta}
=2\rho_0\nu
\left(S_{\alpha\beta}-\frac1D S_{\gamma\gamma}\delta_{\alpha\beta}\right)
+\rho_0\nu^B S_{\gamma\gamma}\delta_{\alpha\beta},
$$

$$
\nu=(1-\chi_s)c_s^2\Delta t\sigma_f^+,
\qquad
\nu^B=\frac2D(1-\chi_b)c_s^2\Delta t\sigma_f^+.
$$

黏性应力是偶二阶矩，故 $s_f^-$ 不进入上述二阶体相输运系数。这个结论不排除 $s_f^-$ 影响三阶误差、ghost 阻尼或壁面位置。

## 5. 标量 CDE 的分量消去

原论文 Eq. (17) 和 Task 2 源矩给出

$$
\sum_i g_i^{eq}=T,
\qquad
\sum_i c_{i\alpha}g_i^{eq}=Tu_\alpha,
$$

$$
\sum_i c_{i\alpha}c_{i\beta}g_i^{eq}
=T(c_s^2+\pi)\delta_{\alpha\beta}+Tu_\alpha u_\beta,
$$

$$
\sum_i R_i=Q,
$$

$$
\sum_i c_{i\alpha}R_i
=\pi\partial_\alpha T
+\frac{TF_\alpha}{\rho_0}
+Qu_\alpha
+\chi_\kappa c_s^2\partial_\alpha T.
$$

零阶矩的一阶方程是

$$
\partial_{t_1}T+\partial_{\alpha1}(Tu_\alpha)=Q.
$$

在奇一阶非平衡通量核
$\langle c_\alpha(D_1g^{eq}-R)\rangle$ 中，利用上述标量守恒层与流场 Euler 层，来自平衡矩导数的各项是

$$
\pi\partial_\alpha T
+\frac{TF_\alpha}{\rho_0}
+Qu_\alpha
+c_s^2\partial_\alpha T
-Tu_\alpha\partial_\beta u_\beta.
$$

逐项减去 odd 源一阶矩后：

- $p\,\partial_\alpha T/\rho_0$ 精确相消；
- $TF_\alpha/\rho_0$ 精确相消；
- $Qu_\alpha$ 精确相消；
- 剩余的目标扩散通量是 $(1-\chi_\kappa)c_s^2\partial_\alpha T$；
- $-Tu_\alpha\partial_\beta u_\beta$ 按论文 Eq. (22) 记为首个 $O(\mathrm{Ma}^3)$ 余项，不强行置零。

因此二阶恢复为

$$
\partial_tT+\partial_\alpha(Tu_\alpha)
=\partial_\alpha(\kappa\partial_\alpha T)+Q
+O(\epsilon^3)+O(\mathrm{Ma}^3),
$$

$$
\boxed{
\kappa=(1-\chi_\kappa)c_s^2\Delta t\sigma_g^-.
}
$$

物理扩散通量是一阶奇矩，故 $s_g^+$ 不进入二阶体相扩散系数；它仍控制 scalar even ghost、半源重构配对以及高阶/边界行为。

## 6. 命名残差不是硬编码零

`second_order_residual_table()` 以两组数据类保存原论文/Tasks 1--2 的平衡矩和源矩系数，再从四个实际率计算

$$
r_f^+=\frac{b_f^+}{s_f^+\sigma_f^+},\quad
r_f^-=\frac{b_f^-}{s_f^-\sigma_f^-},\quad
r_g^+=\frac{b_g^+}{s_g^+\sigma_g^+},\quad
r_g^-=\frac{b_g^-}{s_g^-\sigma_g^-}.
$$

其中真正进入非守恒构成矩的是 $r_f^+$ 与 $r_g^-$。守恒源另定义

$$
c_F=b_f^-+\frac{s_f^-}{2},
\qquad
c_Q=b_g^++\frac{s_g^+}{2}.
$$

残差按下式生成：

$$
\mathcal R_{p\nabla T}
=E_{g,2p}-r_g^-R_{g,1p},
$$

$$
\mathcal R_{TF}
=E_{g,1u}\,c_FR_{f,1F}-r_g^-R_{g,1TF},
$$

$$
\mathcal R_{Qu}
=E_{g,1u}\,c_QR_{g,0Q}-r_g^-R_{g,1Qu},
$$

$$
\mathcal R_{uF}
=E_{f,2uu}\,c_FR_{f,1F}-r_f^+R_{f,2uF},
$$

$$
\mathcal R_{\partial_tF}
=\frac{b_f^-}{s_f^-}-\sigma_f^-,
\qquad
\mathcal R_{\partial_tQ}
=\frac{b_g^+}{s_g^+}-\sigma_g^+.
$$

所有 $E$、$R$ 的默认值都来自命名矩约束且等于 1，所有正确 $b$ 都是各自的 $1-s/2$，所以默认残差经 SymPy 化简为零。这个零是输入矩关系和 CE 系数的结果，而不是零字典。

扰动测试提供反证：删去任一 $R_{g,1p}$、$R_{g,1TF}$、$R_{g,1Qu}$ 或 $R_{f,2uF}$ 后，相应残差变成 1。若把非守恒 scalar odd 或 flow even 因子从 $1-s/2$ 错换成 1，相关构成残差含

$$
1-\frac{2}{2-s}=\frac{s}{s-2}\ne0,
$$

若把守恒 flow odd 或 scalar even 因子错换成 1，则 $c-1=s/2\ne0$；相应时间中心化残差直接变成

$$
\frac1s-\left(\frac1s-\frac12\right)=\frac12\ne0.
$$

代码还对输运表达式直接求偏导：

$$
\frac{\partial\nu}{\partial s_f^-}=0,
\qquad
\frac{\partial\kappa}{\partial s_g^+}=0,
$$

同时验证 $\partial\nu/\partial s_f^+\ne0$、$\partial\kappa/\partial s_g^-\ne0$。这只证明二阶体相输运的 parity 归属。

## 7. 冻结压力与变系数余项

局部冻结 $\pi=p/\rho_0$ 时，第三章的

$$
(c_s^2+\pi)\Delta t\sigma_{\rm flux,eff}
=(1-\chi_\kappa)c_s^2\Delta t\sigma_g^-
$$

严格成立。若只在参考点 $\pi_0$ 冻结碰撞模态，而 $\pi=\pi_0+\delta\pi(\boldsymbol x)$ 变化，则对冻结通量取散度产生

$$
\mathcal R_{\pi,\rm var}
=\Delta t\sigma_{\rm flux,eff}(\pi_0)
\left[
\nabla\delta\pi\cdot\nabla T
+\delta\pi\nabla^2T
\right].
$$

该项一般非零。因此“压力从物理 $\kappa$ 中消去”是局部冻结系数模态结论；若要主张空间变率、边界层或更高阶 modified equation 仍完全消去，必须另做变系数 Taylor 分析。

## 8. 首个遗漏项族与结论边界

本章首个 CE 遗漏阶是 $O(\epsilon^3)$，包括：

- Taylor 展开的 $\Delta t^3D_1^3h^{(0)}/6$、$D_1D_2$ 混合项和 $h^{(2)}$ 通量；
- $\partial_t\boldsymbol F$、$\partial_tQ$ 在二阶 midpoint 配对之后留下的更高时间导数；
- 非均匀 $s$、$\chi$、$p/\rho_0$ 所产生的乘积导数；
- Task 2 已保留的三、四阶 D2Q9 源矩进入的色散/各向异性项；
- ghost 率和边界闭合对三阶及更高误差的影响。

首个 Mach 遗漏阶是 $O(\mathrm{Ma}^3)$。原论文 Eq. (11) 明列的代表项为

$$
-c_s^2(u_\alpha\partial_\beta\delta\rho
+u_\beta\partial_\alpha\delta\rho),
$$

$$
-\rho_0(u_\alpha u_\gamma\partial_\gamma u_\beta
+u_\beta u_\gamma\partial_\gamma u_\alpha),
$$

标量 Eq. (22) 的代表项为

$$
-Tu_\alpha\partial_\beta u_\beta.
$$

Eq. (11) 还含 $O(\mathrm{Ma}^4)$ 的
$-2\rho_0u_\alpha u_\beta\partial_\gamma u_\gamma$。原论文 Eq. (6) 后也把截断 Maxwellian 的首个平衡矩误差标为 $O(\mathrm{Ma}^3)$。

因此本章允许的最终主张是：在上述 CE 与低 Mach 标度、正确的 Task 1--2 矩约束、parity-specific 梯形因子以及局部光滑系数下，连续性、动量和 CDE 恢复到二阶，并得到第三章的物理输运系数。它不是无条件精确格式、变系数四阶结果或边界一致性证明。
