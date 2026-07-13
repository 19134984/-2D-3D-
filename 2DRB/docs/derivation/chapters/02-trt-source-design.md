# 02 奇偶分解源项与算子梯形 TRT 碰撞

本章只建立精确离散代数，不修改生产 Fortran 求解器。所有格点常数均沿用 Task 1 的 D2Q9 精确有理数；`Xs/` 只可作为既有实现对照，不能作为本章证明依据。

## 1. 记号与直接文献约定

对任意九分量向量 $h$，以 Task 1 的相反方向映射定义反转算子 $R$ 和投影算子

$$
P_+=\frac{I+R}{2},\qquad P_-=\frac{I-R}{2},\qquad
h^+=P_+h,\quad h^-=P_-h.
$$

它们满足 $P_+^2=P_+$、$P_-^2=P_-$、$P_+P_-=0$ 和 $P_++P_-=I$。流动分布和标量分布分别使用独立松弛率对
$(s_f^+,s_f^-)$ 与 $(s_g^+,s_g^-)$。

### 1.1 Eq. (13) 的符号与归一化

已直接目视核对原始 `LBM-CDE.pdf` 第 4 页 Eq. (2) 和第 6 页 Eq. (13)。论文定义

$$
S_{\alpha\beta}=\frac12
(\partial_\alpha u_\beta+\partial_\beta u_\alpha),
\qquad S_{\gamma\gamma}=\nabla\cdot\boldsymbol u.
$$

本项目选择无量纲二阶 Hermite 张量

$$
H_{i,\alpha\beta}^{(2)}
=\frac{c_{i\alpha}c_{i\beta}}{c_s^2}-\delta_{\alpha\beta},
$$

并定义二维张量

$$
A_{\alpha\beta}
=\chi_s S_{\alpha\beta}
+\frac{\chi_b-\chi_s}{2}
S_{\gamma\gamma}\delta_{\alpha\beta}.
$$

因此 Eq. (13) 的收缩项必须写成

$$
\boxed{+\rho_0 w_i H_i^{(2)}:A}.
$$

这里是正号；$1/c_s^2$ 已包含在 $H_i^{(2)}$ 中，$A$ 自身不再除以 $c_s^2$。若改用有量纲张量 $c_ic_i-c_s^2I$，才应把相应的 $1/c_s^2$ 移到 $A$ 一侧。两种记法不可混用。

## 2. 源项的奇偶分解

记

$$
\boldsymbol u=(u_x,u_y),\quad
\boldsymbol F=(F_x,F_y),\quad
\boldsymbol G_T=\nabla T,
$$

并定义标量源的一阶目标通量

$$
\boldsymbol J
=\frac{p\boldsymbol G_T+T\boldsymbol F}{\rho_0}
+Q\boldsymbol u+\chi_\kappa c_s^2\boldsymbol G_T.
$$

由离散速度反转 $\boldsymbol c_{\bar i}=-\boldsymbol c_i$ 可直接分出偶、奇项。

### 2.1 流场源

$$
S_i^- = w_i\frac{\boldsymbol c_i\cdot\boldsymbol F}{c_s^2},
$$

$$
S_i^+ = w_i\left[
\frac{(\boldsymbol c_i\cdot\boldsymbol u)
(\boldsymbol c_i\cdot\boldsymbol F)}{c_s^4}
-\frac{\boldsymbol u\cdot\boldsymbol F}{c_s^2}
+\rho_0 H_i^{(2)}:A
\right].
$$

这里第一行是奇函数，第二行的每一项都是偶函数，故代码仍通过 $P_\pm$ 重新投影并逐分量验证，而不是只依赖目测分类。

### 2.2 标量源

原始 PDF 第 8 页 Eq. (24) 给出

$$
R_i^+=w_iQ,
$$

$$
R_i^-=w_i\left[
\frac{\boldsymbol c_i\cdot(p\boldsymbol G_T+T\boldsymbol F)}
{\rho_0c_s^2}
+Q\frac{\boldsymbol c_i\cdot\boldsymbol u}{c_s^2}
+\chi_\kappa\boldsymbol c_i\cdot\boldsymbol G_T
\right]
=w_i\frac{\boldsymbol c_i\cdot\boldsymbol J}{c_s^2}.
$$

## 3. 精确 raw moment 表

定义

$$
M_{pq}[h]=\sum_i c_{ix}^p c_{iy}^q h_i.
$$

再记流场偶源的二阶矩张量

$$
E_{\alpha\beta}
=u_\alpha F_\beta+u_\beta F_\alpha
+2\rho_0c_s^2A_{\alpha\beta}.
$$

### 3.1 零至二阶必需矩

| 分量 | $M_{00}$ | $(M_{10},M_{01})$ | 二阶矩 |
| --- | --- | --- | --- |
| $S^+$ | $0$ | $(0,0)$ | $M_{20}=E_{xx},\ M_{11}=E_{xy},\ M_{02}=E_{yy}$ |
| $S^-$ | $0$ | $(F_x,F_y)$ | 全部为 $0$ |
| $R^+$ | $Q$ | $(0,0)$ | $M_{20}=M_{02}=c_s^2Q,\ M_{11}=0$ |
| $R^-$ | $0$ | $(J_x,J_y)$ | 全部为 $0$ |

特别地，Eq. (13) 的正号和无量纲 $H_i^{(2)}$ 归一化给出

$$
\sum_i c_{i\alpha}c_{i\beta}S_i^+
=u_\alpha F_\beta+u_\beta F_\alpha
+2\rho_0c_s^2A_{\alpha\beta}.
$$

若二阶 Hermite 收缩的符号或 $c_s^2$ 位置错误，这个测试会分别得到错误的负号或多/少一个 $c_s^2$。

### 3.2 三、四阶非零矩

D2Q9 上存在 $c_x^3=c_x$、$c_y^3=c_y$ 的格点别名关系。以下表格列出 Task 5 必须保留的全部三、四阶非零项；表格没有对五阶及以上作任何消失声明。

| 分量 | 非零三阶矩 | 非零四阶矩 |
| --- | --- | --- |
| $S^-$ | $M_{30}=F_x,\ M_{21}=c_s^2F_y,\ M_{12}=c_s^2F_x,\ M_{03}=F_y$ | 无；本阶五个分量均为零 |
| $S^+$ | 无 | $M_{40}=E_{xx},\ M_{31}=M_{13}=E_{xy},\ M_{04}=E_{yy}$；另有 $M_{22}=2c_s^2(\boldsymbol u\cdot\boldsymbol F)+2\rho_0c_s^4\operatorname{tr}A$ |
| $R^-$ | $M_{30}=J_x,\ M_{21}=c_s^2J_y,\ M_{12}=c_s^2J_x,\ M_{03}=J_y$ | 无；本阶五个分量均为零 |
| $R^+$ | 无 | $M_{40}=M_{04}=3c_s^4Q,\ M_{22}=c_s^4Q$ |

对 $R^+$ 还有 $M_{31}=M_{13}=0$。这些值来自 Eq. (13)、Eq. (24) 与精确 D2Q9 权重的逐项求和，不是论文直接列出的高阶闭合条件。

## 4. 从算子梯形规则推导 TRT

对流场或标量场统一写成

$$
\partial_t h+\boldsymbol c\cdot\nabla h
=\Omega(h),\qquad
\Omega(h)=-K(h-h^{eq})+q.
$$

沿特征线采用梯形规则，并定义

$$
\widetilde h=h-\frac{\Delta t}{2}\Omega(h)
=h+\frac{\Delta t}{2}K(h-h^{eq})-\frac{\Delta t}{2}q.
$$

这一符号与原始 PDF Appendix A Eq. (A.2)/(A.14) 完全一致。端点隐式项消去后有

$$
\widetilde h^*-\widetilde h=\Delta t\,\Omega(h).
$$

由变换式解出 $h-h^{eq}$，定义离散松弛算子

$$
\mathsf S
=\Delta t K\left(I+\frac{\Delta t}{2}K\right)^{-1}
=s_+P_++s_-P_-.
$$

于是

$$
\Delta t\,\Omega(h)
=-\mathsf S(\widetilde h-h^{eq})
+\Delta t\left(I-\frac{\mathsf S}{2}\right)q.
$$

把 $q=q^++q^-$ 投影，得到要求的显式 TRT 碰撞：

$$
\boxed{
\begin{aligned}
\widetilde h^*={}&\widetilde h
-s_+P_+(\widetilde h-h^{eq})
-s_-P_-(\widetilde h-h^{eq})\\
&+\Delta t\left[
\left(1-\frac{s_+}{2}\right)q^+
+\left(1-\frac{s_-}{2}\right)q^-
\right].
\end{aligned}}
$$

这说明两个源前因子是算子梯形变换的结果，不是把论文 BGK 的 $1-s/2$ 经验复制两次。

## 5. 半源宏观重构

原始 PDF Eqs. (A.6b)-(A.7b) 与 (A.18)-(A.19) 给出

$$
\boxed{\rho_0\boldsymbol u
=\sum_i\boldsymbol c_i\widetilde f_i
+\frac{\Delta t}{2}\boldsymbol F},
$$

$$
\boxed{T=\sum_i\widetilde g_i+\frac{\Delta t}{2}Q}.
$$

因此 transformed nonequilibrium 的守恒矩并不为零，而是

$$
\sum_i\boldsymbol c_i(\widetilde f_i-f_i^{eq})
=-\frac{\Delta t}{2}\boldsymbol F,
$$

$$
\sum_i(\widetilde g_i-g_i^{eq})
=-\frac{\Delta t}{2}Q.
$$

## 6. 一次碰撞的净源证明

### 6.1 动量

动量属于奇子空间。松弛项对动量的贡献为

$$
-s_f^-\left(-\frac{\Delta t}{2}\boldsymbol F\right)
=\frac{s_f^-\Delta t}{2}\boldsymbol F.
$$

显式奇源贡献为

$$
\Delta t\left(1-\frac{s_f^-}{2}\right)\boldsymbol F.
$$

两者相加，对任意符号松弛率 $s_f^-$ 都严格得到

$$
\boxed{\Delta(\rho_0\boldsymbol u)=\Delta t\boldsymbol F}.
$$

### 6.2 标量

标量零阶守恒矩属于偶子空间。松弛项贡献

$$
\frac{s_g^+\Delta t}{2}Q,
$$

显式偶源贡献

$$
\Delta t\left(1-\frac{s_g^+}{2}\right)Q.
$$

所以对任意 $s_g^+$

$$
\boxed{\Delta T=\Delta t Q}.
$$

忽略半源 nonequilibrium 守恒矩会错误地只剩 $(1-s/2)$ 倍净源。

## 7. 逐分量 BGK 极限

令 $s_+=s_-=s$。利用 $P_++P_-=I$ 和 $q^++q^-=q$，对每个离散方向 $i$ 有

$$
\widetilde h_i^*
=\widetilde h_i-s(\widetilde h_i-h_i^{eq})
+\Delta t\left(1-\frac{s}{2}\right)q_i.
$$

这与 transformed BGK 碰撞逐分量相同，而不只是若干低阶矩相同。自动测试使用任意九分量符号向量比较 TRT 与 BGK 的每一个分量。

## 8. 实现与验证映射

- `tools/derivation/sources.py`：构造 $S,R$，用 Task 1 投影器得到奇偶分量，并精确生成零至四阶 raw moment 表。
- `tools/derivation/collision.py`：实现算子梯形 TRT、逐分量 BGK 和两种半源重构。
- `tests/derivation/test_sources_collision.py`：直接锁定 Eq. (13) 的正号/归一化、全部必需源矩、高阶活动矩、一次碰撞净源和 BGK 极限。

本任务没有改变格点权重、平衡分布或生产求解器。物理不变量是：流场净动量源严格为 $\Delta t\boldsymbol F$，标量净源严格为 $\Delta tQ$，且二者不依赖各自的 TRT 松弛率。
