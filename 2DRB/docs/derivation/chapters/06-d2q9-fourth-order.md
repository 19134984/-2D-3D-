# D2Q9 LBM-CDE-TRT 冻结系数四阶等效方程

## 1. 范围与约定

本章只讨论 $u=F=Q=0$ 且 $p/\rho_0$、松弛参数和
$\chi_\kappa$ 在展开过程中冻结的线性温度模型。结果不包含边界误差，
也不外推到空间变压力、变松弛率或变 $\chi_\kappa$。

全文固定 Fourier 约定

$$
g_i(\boldsymbol x)=\widehat g_i\exp(+\mathrm i\boldsymbol k\cdot\boldsymbol x),
\qquad
E_{ii}(\boldsymbol k)=\exp(-\mathrm i\Delta t\,\boldsymbol k\cdot\boldsymbol c_i),
$$

$$
\nabla T=+\mathrm i\boldsymbol kT,
\qquad
\Gamma(\boldsymbol k)=-\frac{\log z_h(\boldsymbol k)}{\Delta t}.
$$

因此扩散支满足
$\lim_{q\to0}\operatorname{Re}\Gamma(q,0)/q^2=\kappa>0$。这与
Dubois--Lallemand D2Q5 平面波段采用的正相位矩阵 $B$ 相反；两者只有在
$\boldsymbol k\mapsto-\boldsymbol k$ 后才能比较，不能混用源项奇次项的符号。

## 2. 冻结 D2Q9 模型

记

$$
\pi=\frac p{\rho_0},\qquad
a=c_s^2+\pi,\qquad
b=(1-\chi_\kappa)c_s^2,\qquad
d=a-b,
$$

其中标准 D2Q9 的 $c_s^2=1/3$。令 $\ell^T=(1,\ldots,1)$，

$$
e_i=w_i+\frac{\pi}{c_s^2}\lambda_i,\qquad
G=e\ell^T,\qquad
H_{i\alpha}=w_ic_{i\alpha}\frac d{c_s^2},
\qquad
J_{\alpha i}=c_{i\alpha}.
$$

以 $P_+$、$P_-$ 表示速度反演的偶、奇投影，原始 TRT 碰撞为

$$
S=s_eP_+ +s_oP_-,\qquad
C_0=I-S(I-G),\qquad
\sigma_j=\frac1{s_j}-\frac12.
$$

三个实际方案分别是：

$$
C_{\rm base}=C_0\quad(\pi=\chi_\kappa=0),
$$

$$
C_{\rm ext}(\boldsymbol k)
=C_0+\mathrm i\Delta t\left(1-\frac{s_o}{2}\right)
H\boldsymbol k\ell^T,
$$

$$
C_{\rm fb}=C_0-
\frac{2(1-s_o/2)}{a+2b\sigma_o}HJ(I-G).
$$

放大矩阵一律为 $A=EC$。外梯度方案使用论文给定的精确梯度；反馈方案
使用

$$
\nabla T=-\frac{2J(I-G)\widetilde g}
{\Delta t\,(a+2b\sigma_o)}.
$$

特别地，$\chi_\kappa=0$、$\pi\ne0$ 时 $d=\pi\ne0$，论文的
$p\nabla T$ 标量源仍然存在；它不是无源基线。

## 3. 物理通量块与奇 ghost 块

直接精确求和给出

$$
JH=dI_2,\qquad (HJ)^2=dHJ,\qquad \operatorname{rank}(HJ)=2
\quad(d\ne0).
$$

实现中始终构造

$$
K_{i\alpha}=\frac{w_ic_{i\alpha}}{c_s^2},\qquad
P_{\rm flux}=KJ,\qquad
P_{\rm odd\ ghost}=P_--P_{\rm flux},
$$

而不计算 $HJ/d$。所以 $d=0$ 时也没有除零；$HJ/d=P_{\rm flux}$ 只是一条
$d\ne0$ 分支上的证明恒等式。

一般三块生成器允许独立的 $\sigma_f$（两个物理通量模）、$\sigma_o$
（两个奇 ghost 模）和 $\sigma_e$（四个偶 ghost 模）。实际方案随后才专门化为

$$
\sigma_f=\sigma_o\quad\text{(基线、外梯度)},
\qquad
\sigma_f=\frac ba\sigma_o\quad\text{(局部反馈)}.
$$

任意满足 $P_-v=v$、$Jv=0$ 的奇 ghost 向量仍满足
$C_{\rm fb}v=(1-s_o)v$；反馈只改变物理通量块，不能把整个奇投影都换成
$s_f$。

`three_block_amplification()` 是显式的一般三块代数入口。其中任意
$\sigma_f$ 的 `external` 只是一种推广模型，只有取
$\sigma_f=\sigma_o$ 后才是本章的实际 TRT 外梯度方案；`homogeneous` 在取
$\sigma_f=(b/a)\sigma_o$ 后与实际反馈碰撞等价。

## 4. 不得删去的高阶矩

对本项目的 D2Q9 $\lambda_i$，精确 raw moments 为

$$
L_{40}=L_{04}=\frac13,\qquad L_{22}=\frac19,\qquad
L_{31}=L_{13}=0.
$$

因此平衡态高阶矩不是零：

$$
M_{40}^{\rm eq}=M_{04}^{\rm eq}=aT,\qquad
M_{22}^{\rm eq}=c_s^2aT.
$$

外梯度源同样保留由 D2Q9 四阶权重矩诱导的三阶矩：

$$
\sum_i c_{i\beta}c_{i\gamma}c_{i\delta}H_{i\alpha}
=dc_s^2
(\delta_{\alpha\beta}\delta_{\gamma\delta}
+\delta_{\alpha\gamma}\delta_{\beta\delta}
+\delta_{\alpha\delta}\delta_{\beta\gamma}).
$$

两条推导路线均从九个离散速度逐项求这些矩，没有把模型简化为只登记
二阶平衡矩和一阶源矩。

## 5. 两条独立推导路线

Route A 直接构造原始 $9\times9$ 矩阵，用形式递推跟踪
$z_h(0)=1$ 的守恒根，再从 $-\log z_h$ 取到总次数四。符号答案的生成过程
不使用浮点特征值，也不调用 Route B。

Route B 从物理空间恒等式

$$
\exp(\Delta t\,\partial_t)g_i
=\exp(-\Delta t\,\boldsymbol c_i\cdot\nabla)g_i^*
$$

出发，按阶消去非守恒块并形成 PDE 系数。它不构造放大矩阵、特征多项式、
特征值、$z_h$ 或 $\log z_h$，也不调用 Route A。测试同时扫描公开入口及内部
helper 的禁用词，并双向 monkeypatch 另一条路线及其 helper；六个一般有理点
（三种实际方案各两个）上的全部精确残差均为零。

内部递推取 $\Delta t=1$ 以控制符号规模，出口恢复
$\kappa\propto\Delta t$ 和四阶项 $\propto\Delta t^3$。

## 6. 二阶与四阶系数

写等效方程为

$$
\partial_tT=\kappa\nabla^2T
+\Delta t^3\left[
C_{40}(\partial_x^4+\partial_y^4)+C_{22}\partial_x^2\partial_y^2
\right]T+O(\Delta t^4).
$$

一般齐次三块模型给出 $\kappa=a\sigma_f\Delta t$；一般外梯度三块模型给出
$\kappa=b\sigma_f\Delta t$。专门化到三个实际方案后，

$$
\kappa_{\rm base}=\frac{\sigma_o\Delta t}{3},
\qquad
\kappa_{\rm ext}=\kappa_{\rm fb}=b\sigma_o\Delta t.
$$

三个实际方案都精确满足

$$
C_{22}=2C_{40},
$$

即各向同性残差 $C_{22}-2C_{40}$ 恒为零。四阶完全消除仍要求
$C_{40}=C_{22}=0$，不能把“各向同性”误写成“完全消除”。

### 6.1 真正无源基线

$$
C_{40}^{\rm base}
=\frac{\sigma_o(8\sigma_e\sigma_o-4\sigma_o^2-1)}{36},
\qquad C_{22}^{\rm base}=2C_{40}^{\rm base}.
$$

在 $\sigma_o\ne0$ 的主分支，生成的消除条件是

$$
\sigma_e=\frac{4\sigma_o^2+1}{8\sigma_o}.
$$

它是一族条件，不是硬编码的单个 TRT 点。附加
$\sigma_e=2\sigma_o$ 后才得到
$\sigma_o=1/\sqrt{12}$、$\sigma_e=1/\sqrt3$。

### 6.2 精确外梯度

$$
\begin{aligned}
C_{40}^{\rm ext}=-\frac{\sigma_o}{12}\bigl(&
12ab\sigma_e\sigma_o-3ab+a+12b^2\sigma_o^2\\
&-12b\sigma_e\sigma_o+b\bigr),
\qquad C_{22}^{\rm ext}=2C_{40}^{\rm ext}.
\end{aligned}
$$

在 $b\sigma_o(a-1)\ne0$ 的主分支，

$$
\sigma_e=
\frac{3ab-a-12b^2\sigma_o^2-b}{12b\sigma_o(a-1)}.
$$

### 6.3 局部非平衡反馈

先保留一般齐次物理通量位移，得到

$$
C_{40}^{\rm hom}=-\frac{a\sigma_f}{12}
\left[12a\sigma_e\sigma_f+12a\sigma_f^2-3a
-12\sigma_e\sigma_f+2\right],
\qquad C_{22}^{\rm hom}=2C_{40}^{\rm hom}.
$$

实际反馈取 $\sigma_f=(b/a)\sigma_o$。在
$b\sigma_o(a-1)\ne0$ 的主分支，其消除条件为

$$
\sigma_e=
\frac{3a^2-2a-12b^2\sigma_o^2}{12b\sigma_o(a-1)}.
$$

外梯度与反馈虽有相同的二阶扩散率，却一般不具有相同的四阶系数。例如
测试点 $a=4/9$、$b=1/4$、$\sigma_o=2/7$、$\sigma_e=3/11$ 给出

$$
C_{40}^{\rm ext}-C_{40}^{\rm fb}=-\frac{1823}{465696},
\qquad
C_{22}^{\rm ext}-C_{22}^{\rm fb}=-\frac{1823}{232848}.
$$

## 7. Gamma 与 PDE 方向系数

由固定 Fourier 约定，

$$
\Gamma(\boldsymbol k)=\kappa(k_x^2+k_y^2)
-\Delta t^3\left[C_{40}(k_x^4+k_y^4)+C_{22}k_x^2k_y^2\right]
+O(|\boldsymbol k|^6).
$$

所以 $\Gamma$ 的四阶符号与 PDE 的 $C_{40},C_{22}$ 相反。对轴向
$\boldsymbol k=(q,0)$，$K_{\rm axis}=C_{40}$；对等模对角方向
$\boldsymbol k=(q/\sqrt2,q/\sqrt2)$，

$$
K_{\rm diag}=\frac{C_{40}}2+\frac{C_{22}}4,
\qquad C_{22}=4K_{\rm diag}-2K_{\rm axis}.
$$

代码中的 `gamma_qq4` 明确表示 $\boldsymbol k=(q,q)$，因此它等于
$-\Delta t^3(2C_{40}+C_{22})$；旧兼容名 `gamma_diagonal4` 也采用这一
$(q,q)$ 含义。数值拟合使用的 `equal_diagonal` 才是
$(q/\sqrt2,q/\sqrt2)$，两者不可混名。

## 8. Dubois--Lallemand 印刷式审计

原 PDF 第 12 页在 Eq. (44) 后印出的 D2Q9 系数被逐字符编码为

$$
\kappa_{40}^{\rm printed}=\sigma_1\left[
2\sigma_5(\sigma_7-\sigma_3)(a_4-4)
+6\xi\{1-\sigma_1\sigma_7-5\sigma_1\sigma_3
+2\sigma_5(\sigma_7-\sigma_3)\}\right],
$$

$$
\kappa_{22}^{\rm printed}=2A(a_4-4)+12\xi B,
$$

其中

$$
A=\sigma_1+\sigma_5
-2\sigma_1\sigma_5(\sigma_3+\sigma_7+4\sigma_8),
$$

$$
B=\sigma_5+3\sigma_1
-2\sigma_1\sigma_5(\sigma_3+\sigma_7)
-2\sigma_1\sigma_3\sigma_5
-8\sigma_1\sigma_8(\sigma_1+\sigma_5)
+\sigma_1^2\sigma_7.
$$

按同页印刷的 TRT 关系
$\sigma_1=\sigma_5$、$\sigma_3=\sigma_4=\sigma_7=\sigma_8$，在

$$
\sigma_1=\frac1{\sqrt{12}},\qquad
\sigma_3=\frac1{\sqrt3},\qquad
\xi=a_4=\frac13
$$

处逐项代入得到精确审计残差

$$
\kappa_{40}^{\rm printed}=0,
\qquad
\kappa_{22}^{\rm printed}=\frac1{\sqrt3}\ne0.
$$

另一方面，本章与之匹配的零源 D2Q9 模型由两条独立路线都得到
$C_{40}=C_{22}=0$。这里保留不一致，不猜测原文哪个符号有误，也不把印刷式
改写后再称为原文。该论文使用的一般矩模型和归一化只作为外部审计对象，
没有作为本章 LBM-CDE 四阶生成器的输入。

## 9. 80 位定向 Fourier 验证

数值验证直接在每个波数上以至少 80 位精度构造精确矩阵并求靠近 1 的
hydrodynamic 根；它不参与符号公式生成。

一般基线点 $\sigma_o=1/5$、$\sigma_e=4/13$ 的精确
$\Gamma$ 轴向四阶系数为 $217/58500$。在
$q=(1/50,1/100,1/200,1/400)$ 上，轴向和等模对角拟合的末项分别为
$0.00370940280055184$ 和 $0.00370940423498509$，相对相邻波数的观测阶
均收敛到 $4$，且 $\operatorname{Re}\Gamma/q^2>0$。

在零源 TRT 消除点，改用
$q=(1/25,1/50,1/100,1/200)$ 后，扣除二阶扩散项的残差观测阶由 $q^4$
转为 $q^6$；轴向和等模对角末级观测阶均在 $6\pm0.1$ 内。这一递减序列
验证的是阶数变化，不是用单个很小的 $q^4$ 数字代替消除证明。

## 10. 奇异与退化分支

- $a=0$ 是奇异平衡闭合，拒绝进入生成器。
- $a+2b\sigma_o=0$ 只使局部反馈梯度闭合无定义；反馈分支必须拒绝，
  但不含该分母的 baseline/external 分支仍然定义，不能被这个奇点误伤。
- $d=0$ 时 $H=0$，源与反馈均消失，且 $a=b$、$\sigma_f=\sigma_o$；实现不形成 $HJ/d$。
- $b=0$ 时扩散率为零。外梯度仍一般有 $C_{40}=-a\sigma_o/12$；反馈有 $\sigma_f=0$，因而 $C_{40}=C_{22}=0$。
- $a=1$ 时主分支消除公式不能除以 $a-1$。在非零通量位移分支，外梯度的兼容条件变为 $1-2b+12b^2\sigma_o^2=0$，反馈的兼容条件变为 $12(b\sigma_o)^2=1$；零位移分支另由消失的整体前因子判断。求解器将这些分支区分为 `identically_satisfied` 或 `incompatible`。
- $\chi_\kappa\to1$ 或 $\sigma_o\to0$ 是零扩散或边缘松弛极限，不是默认物理分支。
- $\pi=\chi_\kappa=0$ 且 $\sigma_o=\sigma_e$ 恢复标准 BGK D2Q9。

Wang/Luo 的 D2Q5 四阶关系和边界 magic 常数没有进入任何 D2Q9 体相
四阶计算；它们不能作为本章系数或消除条件的输入。
