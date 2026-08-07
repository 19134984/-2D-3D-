# 基于 LBM-CDE 的 TRT 推广、边界 magic 与温度四阶误差推导

## 摘要

本文从用户提供的 `LBM-CDE.pdf` 的梯形 transformed 离散方程出发，重新完成
奇偶分解、TRT 源项、局部反馈消元、二阶宏观恢复、源感知边界残差与 D2Q9
四阶等效方程推导。目标不是把 BGK 公式中的单一松弛率机械替换成两个数，
而是回答三个可验证的问题：TRT 推广是否保持原 LBM-CDE 的宏观方程；所谓
magic 参数在加入局部反馈、压力修正和热源后还剩多大适用域；体相四阶消除、
低扩散率与边界标定能否同时满足。

全文所有主要结论分为四级：`strictly_proved` 表示已由精确代数与至少一条独立
路线闭环；`restricted_model` 表示只在明确的平壁、稳态、冻结系数或二次场假设
下成立；`numerical_evidence` 表示高精度 Fourier 或有限算例抽查；`unresolved`
表示一般源项/边界 jet 尚未闭合。`Xs/` 仅是此前复现代码的实现比较对象，不作为
本文理论证明的依据。

## 主要结论

1. **BGK 推广到 TRT 是可行的，但源项必须按奇偶子空间分别乘因子。**
   对 transformed 分布，正确碰撞项包含
   $(1-s_+/2)\Phi_i^+$ 与 $(1-s_-/2)\Phi_i^-$。当 $s_+=s_-$ 时逐分量退化回
   LBM-CDE 的 BGK 形式；因此这是算子梯形规则的推广，而非经验替换。

2. **局部反馈改变物理块的有效 Hénon shift，不会自动改变同奇偶性的 ghost。**
   冻结系数、CE2 范围内，剪切、体积和温度通量块分别满足
   $\sigma_{\rm sh}^{\rm eff}=(1-\chi_s)\sigma_f^+$、
   $\sigma_{\rm bu}^{\rm eff}=(1-\chi_b)\sigma_f^+$ 与
   $\sigma_{\rm q}^{\rm eff}=b\sigma_g^-/a$，其中
   $a=c_s^2+\pi$、$b=(1-\chi_\kappa)c_s^2$。名义 ghost shift 仍是碰撞输入本身。

3. **不存在脱离边界实现与源规范的普适 magic 常数。**
   $\sigma_f^+\sigma_f^-=3/16$ 只对应平直、格点对齐、halfway、稳态 Stokes、
   均匀体力、half-force 重构与无 feedback 下的受限标定；压力边界驱动文献中的 $3/8$ 是另一
   gauge。一般压力、力、时间和法/切向导数 jets 仍要求显式边界修正。

4. **D2Q9 温度四阶条件必须针对实际 external 或 local-feedback 方案重推。**
   两种方案二阶扩散率都为 $\kappa=b\sigma_o\Delta t$，但其四阶系数一般不同。
   本文保留 $\lambda_i$ 的非零四阶矩，并用放大矩阵流体根与物理空间 Taylor/矩
   消元两条路线得到 $C_{22}=2C_{40}$ 及各自的完整消除条件。

5. **体相四阶消除与受限 ABB 标定会锁定扩散尺度。**
   令 $K=\kappa/\Delta t$，与受限 ABB 条件联立后，local-feedback 与 external
   分别给出
   $$
   K^2=\frac{a(3a+1)}{48},\qquad
   K^2=\frac{12ab+5a-4b-9a^2}{48}.
   $$
   这里的 ABB 条件仅指 R8 子域：平直、格点对齐的 halfway 壁，冻结压力、稳态
   一维二次温度场、half-source 重构、完整 pressure-wall equilibrium term 与匹配的
   external/local population chain、CDE-consistent 均匀 $Q$，且无流/力/切向 jets；
   它不是一般边界恒等式。
   在 $\pi=0$ 的正支，两式都退化为 $K=1/\sqrt{72}$。因此“任意低 $\kappa$、
   体相四阶完全消除、同一个受限 ABB magic”三项通常不可同时强制。

6. **低扩散率下的诚实实现选择有三类。**
   一条已识别但仍待推导、实现和验证的方向，是保持体相四阶条件并加入显式受限 ABB 修正；
   也可保留受限 ABB 标定而
   接受非零 $q^4$ 误差；若要用 split-even MRT 同时释放两项约束，必须另行证明
   模态归属与参数 Jacobian 满秩，本文只把它列为待推导候选。

## 记号与目标方程

速度格点为 $\boldsymbol c_i$，反向指标为 $\bar i$，奇偶部分定义为
$$
h_i^+=\frac{h_i+h_{\bar i}}2,\qquad
h_i^-=\frac{h_i-h_{\bar i}}2.
$$
本文 D2Q9 速度顺序固定为
$(0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)$；
相应 LBM-CDE 压力修正向量由低阶矩约束显式得到
$$
\boldsymbol\lambda_t=
\left(-\frac59,\frac19,\frac19,\frac19,\frac19,
\frac1{36},\frac1{36},\frac1{36},\frac1{36}\right)^{\mathsf T}.
$$
Hénon shift 与松弛率的关系为
$$
\sigma=\frac1s-\frac12,\qquad s=\frac1{\sigma+1/2}.
$$
本文区分名义碰撞 shift、局部反馈消元后的物理有效 shift、受限边界乘积和一般
边界 residual；四者不能互换。

目标低 Mach、弱可压宏观方程直接沿用第四章的恢复形式：
$$
\partial_t\delta\rho+\partial_\alpha(\rho_0u_\alpha)=O(\epsilon^3),
$$
$$
\partial_t(\rho_0u_\alpha)+\partial_\beta(\rho_0u_\alpha u_\beta)
=-\partial_\alpha p+\partial_\beta\tau_{\alpha\beta}+F_\alpha
+O(\epsilon^3)+O(\mathrm{Ma}^3),
$$
$$
\tau_{\alpha\beta}
=2\rho_0\nu\left(S_{\alpha\beta}-\frac1D S_{\gamma\gamma}
\delta_{\alpha\beta}\right)
+\rho_0\nu^B S_{\gamma\gamma}\delta_{\alpha\beta},
$$
$$
\nu=(1-\chi_s)c_s^2\Delta t\,\sigma_f^+,
\qquad
\nu^B=\frac2D(1-\chi_b)c_s^2\Delta t\,\sigma_f^+,
$$
以及
$$
\partial_tT+\nabla\!\cdot(T\boldsymbol u)
=\nabla\!\cdot(\kappa\nabla T)+Q.
$$
二阶恢复章节会逐项说明梯形半源重构、压力修正和冻结系数余项；这些方程不被
直接假定为结果。

## transformed TRT 离散骨架

以流场为例，一步碰撞--迁移写成
$$
\widetilde f_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
-\widetilde f_i(\boldsymbol x,t)
=\Omega_i^{\rm TRT},
$$

$$
\Omega_i^{\rm TRT}
=-s_f^+\bigl(\widetilde f_i^+-f_i^{\mathrm{eq},+}\bigr)
-s_f^-\bigl(\widetilde f_i^--f_i^{\mathrm{eq},-}\bigr)
+\Delta t\,\Psi_i,
$$

其中 parity-specific 梯形源增量为

$$
\Psi_i=\left(1-\frac{s_f^+}{2}\right)\Phi_i^+
+\left(1-\frac{s_f^-}{2}\right)\Phi_i^-.
$$
温度分布 $\widetilde g_i$ 具有同一算子结构，只需把平衡态、源矩与松弛率替换为
标量对应量。后续章节从这一离散式推导有效率、宏观方程和高阶误差。

## 正文组织

完整 Markdown 与 PDF 按下列顺序装配 reviewed chapter：

1. 奇偶源项与算子梯形 TRT 碰撞；
2. 局部反馈消元与物理模态有效率；
3. 二阶 Chapman--Enskog/等效方程恢复与 BGK 极限；
4. 源感知速度、ABB、绝热与混合角点边界；
5. D2Q9 温度四阶等效方程；
6. 参数兼容性、不可兼容性与低扩散率选择；
7. 实现顺序、伪代码与运行时不变量；
8. D2Q5 文献公式复算附录；
9. 证据矩阵、限制与参考文献。

下文直接给出全部详细公式、适用假设、程序入口与测试映射，不省略边界残差和
四阶推导细节。各分章文件同时保留为可独立审查的生成源。

# 奇偶分解源项与算子梯形 TRT 碰撞

本章只建立精确离散代数，不修改生产 Fortran 求解器。所有格点常数均沿用 Task 1 的 D2Q9 精确有理数；`Xs/` 只可作为既有实现对照，不能作为本章证明依据。

## 记号与直接文献约定

对任意九分量向量 $h$，以 Task 1 的相反方向映射定义反转算子 $R$ 和投影算子

$$
P_+=\frac{I+R}{2},\qquad P_-=\frac{I-R}{2},\qquad
h^+=P_+h,\quad h^-=P_-h.
$$

它们满足 $P_+^2=P_+$、$P_-^2=P_-$、$P_+P_-=0$ 和 $P_++P_-=I$。流动分布和标量分布分别使用独立松弛率对
$(s_f^+,s_f^-)$ 与 $(s_g^+,s_g^-)$。

### Eq. (13) 的符号与归一化

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

## 源项的奇偶分解

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

### 流场源

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

### 标量源

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

## 精确 raw moment 表

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

### 零至二阶必需矩

- $S^+$：$M_{00}=0$、一阶矩为零，二阶矩为
  $M_{20}=E_{xx}$、$M_{11}=E_{xy}$、$M_{02}=E_{yy}$。
- $S^-$：$M_{00}=0$、一阶矩为 $(F_x,F_y)$，二阶矩全部为零。
- $R^+$：$M_{00}=Q$、一阶矩为零，二阶矩为
  $M_{20}=M_{02}=c_s^2Q$、$M_{11}=0$。
- $R^-$：$M_{00}=0$、一阶矩为 $(J_x,J_y)$，二阶矩全部为零。

特别地，Eq. (13) 的正号和无量纲 $H_i^{(2)}$ 归一化给出

$$
\sum_i c_{i\alpha}c_{i\beta}S_i^+
=u_\alpha F_\beta+u_\beta F_\alpha
+2\rho_0c_s^2A_{\alpha\beta}.
$$

若二阶 Hermite 收缩的符号或 $c_s^2$ 位置错误，这个测试会分别得到错误的负号或多/少一个 $c_s^2$。

### 三、四阶非零矩

D2Q9 上存在 $c_x^3=c_x$、$c_y^3=c_y$ 的格点别名关系。以下表格列出 Task 5 必须保留的全部三、四阶非零项；表格没有对五阶及以上作任何消失声明。

- $S^-$ 的非零三阶矩为 $M_{30}=F_x$、$M_{21}=c_s^2F_y$、
  $M_{12}=c_s^2F_x$、$M_{03}=F_y$；本阶五个四阶分量均为零。
- $S^+$ 的三阶矩为零；四阶矩为 $M_{40}=E_{xx}$、
  $M_{31}=M_{13}=E_{xy}$、$M_{04}=E_{yy}$，以及
  $M_{22}=2c_s^2(\boldsymbol u\cdot\boldsymbol F)
  +2\rho_0c_s^4\operatorname{tr}A$。
- $R^-$ 的非零三阶矩为 $M_{30}=J_x$、$M_{21}=c_s^2J_y$、
  $M_{12}=c_s^2J_x$、$M_{03}=J_y$；四阶矩为零。
- $R^+$ 的三阶矩为零；四阶矩为
  $M_{40}=M_{04}=3c_s^4Q$、$M_{22}=c_s^4Q$。

对 $R^+$ 还有 $M_{31}=M_{13}=0$。这些值来自 Eq. (13)、Eq. (24) 与精确 D2Q9 权重的逐项求和，不是论文直接列出的高阶闭合条件。

## 从算子梯形规则推导 TRT

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

## 半源宏观重构

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

## 一次碰撞的净源证明

### 动量

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

### 标量

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

## 逐分量 BGK 极限

令 $s_+=s_-=s$。利用 $P_++P_-=I$ 和 $q^++q^-=q$，对每个离散方向 $i$ 有

$$
\widetilde h_i^*
=\widetilde h_i-s(\widetilde h_i-h_i^{eq})
+\Delta t\left(1-\frac{s}{2}\right)q_i.
$$

这与 transformed BGK 碰撞逐分量相同，而不只是若干低阶矩相同。自动测试使用任意九分量符号向量比较 TRT 与 BGK 的每一个分量。

## 实现与验证映射

- `tools/derivation/sources.py`：构造 $S,R$，用 Task 1 投影器得到奇偶分量，并精确生成零至四阶 raw moment 表。
- `tools/derivation/collision.py`：实现算子梯形 TRT、逐分量 BGK 和两种半源重构。
- `tests/derivation/test_sources_collision.py`：直接锁定 Eq. (13) 的正号/归一化、全部必需源矩、高阶活动矩、一次碰撞净源和 BGK 极限。

本任务没有改变格点权重、平衡分布或生产求解器。物理不变量是：流场净动量源严格为 $\Delta t\boldsymbol F$，标量净源严格为 $\Delta tQ$，且二者不依赖各自的 TRT 松弛率。

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

# 梯形 transformed LBE 的二阶宏观恢复

本章从 Task 1 的平衡矩和 Task 2 的源矩出发，对 transformed 离散 LBE 作二阶 Chapman--Enskog（CE）展开。目标方程只用于最后比对，不作为符号残差的输入夹具。

## 两套互不替代的标度

### CE 慢变量标度

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

### 低 Mach 幅值标度

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

## transformed 离散方程与三层人口方程

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

## 梯形半源对时间离散项的消去

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

以下把连续性、动量和 CDE 的首个离散遗漏写为 $O(\epsilon^3)$ 时，采用精确假设

$$
\boxed{\texttt{coefficient\_assumption = frozen\_through\_ce2}}.
$$

即 $s_f^\pm,s_g^\pm,\chi_s,\chi_b,\chi_\kappa$ 与 $\pi=p/\rho_0$ 在 CE 二阶内为常数/冻结系数。若这些系数在慢变量上变化，它们的乘积导数必须在二阶方程中保留。

## 连续性与动量方程

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

## 标量 CDE 的分量消去

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

## 命名残差不是硬编码零

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

## 冻结压力与二阶变系数余项

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

该项一般非零。在 $\nabla=\epsilon\nabla_1$ 且 $\delta\pi=O(1)$ 的系数计数下，两项都是 $O(\epsilon^2)$；光滑性只保证导数存在，不会把它们提升到 $O(\epsilon^3)$。因此“压力从物理 $\kappa$ 中消去”是局部冻结系数模态结论。若系数不冻结，必须把该余项以及空间变 $s$、$\chi$ 引起的同类乘积导数纳入二阶目标方程。

残差元数据据此分离为

```text
coefficient_variation_epsilon_order = 2
coefficient_assumption = frozen_through_ce2
constant_coefficient_first_omitted_epsilon_order = 3
```

旧键 `first_omitted_epsilon_order=3` 只保留为常系数路径的兼容值，不能用它给变系数余项定阶。

## 常/冻结系数下的首个遗漏项族与结论边界

在 `frozen_through_ce2` 假设下，首个 CE 遗漏阶是 $O(\epsilon^3)$，包括：

- Taylor 展开的 $\Delta t^3D_1^3h^{(0)}/6$、$D_1D_2$ 混合项和 $h^{(2)}$ 通量；
- $\partial_t\boldsymbol F$、$\partial_tQ$ 在二阶 midpoint 配对之后留下的更高时间导数；
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

因此本章允许的最终主张是：在上述 CE 与低 Mach 标度、正确的 Task 1--2 矩约束、parity-specific 梯形因子，并且 $s$、$\chi$、$p/\rho_0$ 在 CE2 内为常数/冻结系数时，连续性、动量和 CDE 恢复到二阶，并得到第三章的物理输运系数。仅有局部光滑性不足以满足这个假设。它不是无条件精确格式、一般变系数二阶闭合、变系数四阶结果或边界一致性证明。

# 源感知边界残差与 magic 分类

本章不把某个松弛率乘积先验称为“magic”。判定顺序固定为：先写 transformed population 的碰撞后链，保留独立墙面 Taylor jets，再检查同一条件能否同时消去所有非零系数。`classify_magic()` 同时返回系数表、假设、参数映射、未满足 jets 和已显式相消的零行。对于 general velocity，返回值还包含 `rate_compatibility_status`：它只审计已解析的 shear/bulk 两行，不能代替完整墙面闭合状态。

## 记号、率映射和六个适用域

Hénon 平移量统一记为

$$
\sigma=\frac1s-\frac12.
$$

任何数值结论必须同时附带 normalization、驱动/源实现、宏观量 gauge 和几何。Task 6 分开处理六个适用域：

- **classical velocity：** D2Q9 两率、均匀体力、half-force momentum、平直
  halfway wall、无反馈时，
  $\Lambda^2=1/4\Leftrightarrow\sigma_f^+\sigma_f^-=3/16$；分类为
  `restricted_calibration`。
- **flow feedback：** 一维剪切、均匀体力、冻结反馈时，候选条件为
  $(1-\chi_s)\sigma_f^+\sigma_f^-=3/16$；仍是 `restricted_calibration`。
- **scalar external-gradient：** D2Q9 温度 ABB、显式精确梯度源、稳态一维
  二次场下，$\sigma_{g,\mathrm{flux}}\sigma_g^+=3/16$。
- **scalar local-feedback：** homogeneous physical-flux block 在同一稳态一维
  二次场产生相同 ABB 行，但不得外推到一般 source/time/force jets。
- **general source-aware wall：** velocity BB、ABB 或绝热 BB 的压力、热源、
  时间、法/切向导数和力 jets 独立；未闭合项标为 nonzero unresolved，顶层恒为
  `boundary_correction_required`，velocity 另报 rate 子系统状态。
- **mixed corner：** 一个共享对角 population 同时受 Dirichlet ABB 与绝热 BB
  约束；一般 $\operatorname{rank}A=1<\operatorname{rank}[A|b]=2$，分类为
  `corner_closure_conflict`。

压力边界驱动 Poiseuille 的 $3/8$ 是独立子域，不属于均匀体力一行。Wang/Luo 的 D2Q5 各向同性 $1/6$ 和特定稳定参数 $1/4$ 也不属于本章 D2Q9 温度边界标定。

## Classical velocity：1/4 与 3/16 是同一受限标定

Ginzburg--d'Humières 的 normalization 是

$$
\Lambda^2=\frac43
\left(\frac1{\lambda_\nu}-\frac12\right)
\left(\frac1{\lambda_2}-\frac12\right),
$$

$$
D_{\mathrm{eff}}^2=D_{1/2}^2+4\Lambda^2-1.
$$

在 Wang/Luo 的 D2Q9 矩排序中，黏性应力是中心对称 even block，三阶伴随通量是 odd block，因此

$$
\lambda_\nu=s_f^+,\qquad
\lambda_2=s_f^-,\qquad
\Lambda^2=\frac43\sigma_f^+\sigma_f^-.
$$

于是精确 halfway channel 条件

$$
\Lambda^2=\frac14
$$

与

$$
\sigma_f^+\sigma_f^-=\frac3{16}
$$

完全等价。这里的完整假设是：D2Q9 两率碰撞、平直静止 grid-aligned halfway wall、稳态线性 Stokes 平衡、均匀体力、论文的 half-force momentum gauge、无 LBM-CDE feedback。它不是 universal constant。

Luo 的 pressure-boundary drive 使用另一套 forcing/gauge；代码把

$$
\sigma_f^+\sigma_f^-=\frac38
$$

保留为独立 `pressure_boundary` 分支，绝不与均匀体力的 $3/16$ 合并。

直接制造解不调用 `velocity_bb_residual()`。线性 Couette 只验证 halfway 位置和符号；hydrostatic rest 使用

$$
\rho_0\boldsymbol u=\sum_i\boldsymbol c_i\widetilde f_i+\frac12\boldsymbol F.
$$

验证 transformed momentum 的 $-\boldsymbol F/2$ 被 half-force reconstruction 恢复，不产生伪滑移。Uniform-force 的 $3/16$ 与 pressure-drive 的 $3/8$ 当前只编码为 `source_evidence_only`：它们保留各自论文 normalization/gauge 的系数映射，但 `probes_quadratic_slip=False`，没有被描述为由本项目独立二次场或人口 stencil 重新生成。二者仍分支记录，避免把两个 gauge 混为一谈。

## Flow feedback：物理块与 nominal ghost 不可混用

Task 3 的映射是

$$
\sigma_{f,\mathrm{shear}}^{\mathrm{phys}}
=(1-\chi_s)\sigma_f^+,\qquad
\sigma_{f,\mathrm{bulk}}^{\mathrm{phys}}
=(1-\chi_b)\sigma_f^+.
$$

未参与该物理闭合的 ghost 仍是

$$
\sigma_{f,\mathrm{even\ ghost}}=\sigma_f^+,\qquad
\sigma_{f,\mathrm{odd\ ghost}}=\sigma_f^-.
$$

仅在稳态一维剪切、均匀体力、无 bulk/time/tangential jets 时，二次 slip 行给出候选

$$
(1-\chi_s)\sigma_f^+\sigma_f^-=\frac3{16}.
$$

一般场同时保留

$$
C_{\mathrm{shear}}
=\frac{16}{3}(1-\chi_s)\sigma_f^+\sigma_f^- -1,
$$

$$
C_{\mathrm{bulk}}
=\frac{16}{3}(1-\chi_b)\sigma_f^+\sigma_f^- -1.
$$

这里必须区分两层状态。

- 完整 general velocity 表还保留法/切向速度一、二阶导数，wall-time 导数，压力与压力梯度，力与力梯度，以及 velocity/source/time jets；这些尚未由完整 population chain 闭合的系数均以 nonzero unresolved symbol 保留。因此顶层 `status` 恒为 `boundary_correction_required`。
- `rate_compatibility_status` 只检查上面的 shear/bulk 两行。当 $\chi_s=\chi_b$ 时，两行共享一个零点，子系统为 `restricted_calibration`，且 `rate_compatibility_conditions` 返回该 product 的唯一精确解；当 $\chi_s\ne\chi_b$ 时，同一 nominal product 不能同时令两行归零，子系统为 `no_single_magic`，条件映射为空。

即使 $\chi_s=\chi_b$，也只能说“已解析的 rate 子系统兼容”，不能说完整边界已经闭合。

## Temperature Dirichlet ABB：两条独立 population 路线

### 墙面平衡项与 transformed source

每条 Dirichlet ABB link 使用完整墙面项

$$
2\left(w_i+\lambda_{t,i}\frac{\pi_w}{c_s^2}\right)T_w,
\qquad \pi=\frac p{\rho_0},\qquad c_s^2=\frac13.
$$

温度宏观重构与源因子是

$$
T=\sum_i\widetilde g_i+\frac Q2,
$$

$$
\left(1-\frac{s_g^+}{2}\right)R_i^+
+\left(1-\frac{s_g^-}{2}\right)R_i^-.
$$

物理 scalar-flux Hénon shift 与 nominal ghosts 分开记录：

$$
\sigma_{g,\mathrm{flux}}^{\mathrm{phys}}
=\frac{(1-\chi_\kappa)c_s^2}{c_s^2+\pi}\sigma_g^-,
$$

$$
\sigma_{g,\mathrm{odd\ ghost}}=\sigma_g^-,\qquad
\sigma_{g,\mathrm{even\ ghost}}=\sigma_g^+.
$$

### 显式精确梯度源链

把三个法向 D2Q9 links 聚合为 D1Q3，权重为 $W_0=2/3$、$W_+=W_-=1/6$，并代入

$$
T(x)=T_0+T_1x+T_2x^2.
$$

求解 collision、streaming、$\sum\widetilde g=T-Q/2$ 和 ABB link 的有限多项式方程，得到

$$
\frac{R_{\mathrm{ABB}}}{T_2}
=\frac{16(1-\chi_\kappa)\sigma_g^+\sigma_g^-
-3(1+3\pi)}{36}.
$$

用物理通量 shift 改写为

$$
\frac{R_{\mathrm{ABB}}}{T_2}
=\frac{1+3\pi}{36}
\left(16\sigma_{g,\mathrm{flux}}^{\mathrm{phys}}\sigma_g^+-3\right).
$$

只要 $1+3\pi\ne0$，受限稳态二次行的零条件就是

$$
\sigma_{g,\mathrm{flux}}^{\mathrm{phys}}\sigma_g^+=\frac3{16}.
$$

### 局部非平衡反馈消元链

第二条路线不使用显式梯度源。它把 local nonequilibrium feedback 消去为 homogeneous physical-flux block，使用

$$
s_{g,\mathrm{flux}}
=\frac1{\sigma_{g,\mathrm{flux}}^{\mathrm{phys}}+1/2},
$$

但 even heat source 仍使用原始 D1Q3 权重 $(2/3,1/6,1/6)$，不能换成 pressure-modified equilibrium weights。直接解同一二次 population recurrence 得到

$$
\frac Q{T_2}=-2(c_s^2+\pi)\sigma_{g,\mathrm{flux}}^{\mathrm{phys}},
$$

$$
\frac{R_{\mathrm{ABB}}^{\mathrm{local}}}{T_2}
=\frac{1+3\pi}{36}
\left(16\sigma_{g,\mathrm{flux}}^{\mathrm{phys}}\sigma_g^+-3\right).
$$

因此两条路线只在“冻结压力、稳态一维二次温度、CDE-consistent uniform $Q$、无流/力/切向/墙时变”这一行相同。一般 source/time/force/variable-pressure 系数必须分别推导，不能从该相同结果外推。

### 主生成路线与一般 ABB 独立 jets

主 API 不再直接填写最终系数表。`_temperature_abb_population_chain()` 先生成有限 D1Q3 population、collision、streaming、宏观重构与 ABB 方程，再从解中提取 steady quadratic 行；affine $T\pi$ 墙面平衡乘积由独立有限差分方程生成。`temperature_abb_residual()` 只消费该生成结果。测试通过调用追踪锁定这条依赖关系，而直接制造解仍使用独立 helper，不复用主生成器。

Local-feedback 路线在本任务中只独立证明了上一节 steady 1D quadratic ABB 行相同；一般 pressure/source/time/force/tangential wall chain 仍未闭合，结论是 unresolved 且需要 boundary correction，不能读取为已经从受限行完整外推。

External-gradient 一般表显式保留

$$
T,T_n,T_\tau,T_{nn},T_{n\tau},T_{\tau\tau},
T_t,T_{tt},T_{tn},T_{t\tau},
$$

以及 $\pi,\pi_n,\pi_\tau,F_n,F_\tau,Q,Q_n,Q_\tau,u_n,u_\tau$ 和 $p\nabla T$、$T\boldsymbol F$、$Q\boldsymbol u$ 源 jets。当前已由生成方程解析的是 steady quadratic $T_{nn}$ 行和 affine $T\pi$ 的 $\pi_n$ 行；其余未闭合项不是零，而是显式 nonzero unresolved coefficient，并全部进入 `unsatisfied_jets`。

两个最小链说明一个 product 不够：

- 非定常 ABB 的独立制造入口会消费 wall-time 输入，但在本任务冻结稳态假设之外，明确返回 `unsupported_unresolved`，不伪造一个已验证系数；
- affine $T$、affine $\pi$ 且 $F_n/\rho_0=\pi_n$ 的直接 $T\pi$ 展开给出率无关残差

$$
R_{\pi_n}=-\frac14T_n\pi_n.
$$

因此 general ABB 是 `boundary_correction_required`，而不是 universal magic。

## Adiabatic BB：先求 kinetic odd flux

`_adiabatic_population_chain()` 不预填法向或对角系数：它在 `flat_grid_aligned_halfway`、`half_source`、`transformed_cde_chain`、D2Q9 $c_s^2=1/3$ 下，分别求解 affine 法向人口对和两条 quadratic tangential diagonal population 的 collision、streaming、odd-source 与 reflection 有限方程，再从求解后的 reflected defect 提取法向 kinetic odd flux 与切向曲率。方程显式含 nominal even/odd rate、$c_s^2$ 与压力比；改变 even rate 或压力会改变中间方程，而对称消元后的两个受限系数保持不变。`adiabatic_bb_residual()` 只消费该生成结果。压力、力、热源导数、wall-time 和法向曲率等一般行仍显式保留为 nonzero unresolved。

绝热 BB 首先给出 reflected odd-flux condition，而不是先代入 $T_n=0$。下述一维 affine 法向链还要求稳态、常压力比与零切向 jets，并产生

$$
R_{\mathrm{odd}}
=-\frac{1-\chi_\kappa}{3}\sigma_g^-T_n.
$$

这里没有可调的 even--odd product；$T_n=0$ 是边界数据，不是 magic relation。

施加 smooth homogeneous Neumann data 后，两条 diagonal links 的切向一阶项相消，但二阶项相加为

$$
R_{\tau\tau}
=\frac{1-\chi_\kappa}{9}\sigma_g^-T_{\tau\tau}.
$$

独立 source probes 还必须分开：

$$
R_F=\sigma_g^-T(F_n/\rho_0-\pi_n),
$$

所以 force-only 非零，而 hydrostatic pair $F_n/\rho_0=\pi_n$ 精确相消。Uniform $Q$ probe 把同一个 $(1-s_g^+/2)w_iQ$ 增量实际加入正、负人口，再由两者之差证明 $Q$ 不进入 homogeneous odd-flux 行；推导对象保留 $Q$，不是用常数零替代输入。绝热 wall-time/法向曲率制造入口在冻结稳态范围之外明确返回 `unsupported_unresolved`。一般绝热墙由此分类为 `boundary_correction_required`。

## Mixed corner：覆盖顺序不是同时满足

在采用 `grid_aligned_right_angle_corner`、half-source 规范的 Dirichlet/adiabatic 直角角点，同一个 diagonal incoming population $z$ 被两条墙方程同时使用：

$$
z=b_D=-g_{\bar i}^*+2g_{i,w}^{eq},
\qquad
z=b_A=g_{\bar i}^*.
$$

因此

$$
A=\begin{bmatrix}1\\1\end{bmatrix},
\qquad
A z=\begin{bmatrix}b_D\\b_A\end{bmatrix}.
$$

一般 $b_D\ne b_A$ 时

$$
\operatorname{rank}A=1,
\qquad
\operatorname{rank}[A|b]=2.
$$

对制造场

$$
T(x,y)=T_w+T_{xy}xy
$$

在半对角 fluid point $(h/2,h/2)$ 上，直接平衡 link 给出两种覆盖顺序之差的首个 $h^2$ 系数

$$
\frac{b_D-b_A}{T_{xy}h^2}
=-\frac{1+3\pi}{72}.
$$

这项不含 $\sigma_g^+$、$\sigma_g^-$ 或 $\chi_\kappa$，无法靠 rate tuning 消除。`_corner_population_chain()` 从两条共享对角人口赋值生成方程、秩、距离和源计数；naive 计数等于赋值序列长度，共享计数等于人口标识去重后的长度，因此得到两次覆盖与一次共享源，而不是在结果表中写死 `2/1`。独立 corner 制造解另行构造同样的有限赋值序列。分类必须是 `corner_closure_conflict`。

## 结论边界

- `universal_magic`：本章没有找到满足该定义的条件。
- `restricted_calibration`：只用于假设完整写出的 classical、1D shear 或 1D scalar quadratic 行。
- `no_single_magic`：在本章作为 `rate_compatibility_status`，表示已解析 shear/bulk 子系统要求不兼容条件；它不覆盖完整 general wall 状态。
- `boundary_correction_required`：率无关或不受同一 product 控制的 source/time/pressure/tangential 行仍在。
- `corner_closure_conflict`：共享 unknown 的方程秩不兼容或 sequential overwrite 顺序相关。

特别地，$1/4$、$3/16$、$3/8$ 均不得脱离各自 normalization、驱动/gauge 和 geometry 单独引用；D2Q5 bulk constants 也不得移植为 D2Q9 boundary magic。

# D2Q9 LBM-CDE-TRT 冻结系数四阶等效方程

## 范围与约定

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

## 冻结 D2Q9 模型

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

## 物理通量块与奇 ghost 块

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

## 不得删去的高阶矩

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

## 两条独立推导路线

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

## 二阶与四阶系数

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

### 真正无源基线

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

### 精确外梯度

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

### 局部非平衡反馈

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

## Gamma 与 PDE 方向系数

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

## Dubois--Lallemand 印刷式审计

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

## 80 位定向 Fourier 验证

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

## 奇异与退化分支

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

# TRT 参数兼容性与不可兼容性

本章把第五章的受限壁面残差与第六章的冻结系数四阶体相条件放进同一个精确求解器。目标不是选一个文献参数点，而是回答：给定输运系数、名义 TRT 率、物理有效块、体相四阶条件和壁面标定后，这组约束是否仍有解。

## 量、约束与报告接口

统一记

$$
\sigma=\frac1s-\frac12,
\qquad
a=c_s^2+\pi,
\qquad
b=(1-\chi_\kappa)c_s^2,
\qquad
K=\frac\kappa{\Delta t}.
$$

标量名义奇 ghost、名义偶 ghost 和反馈后的物理通量平移分别为

$$
\sigma_o=\frac Kb,
\qquad
\sigma_e,
\qquad
\sigma_{\rm flux}=\frac ba\sigma_o=\frac Ka.
$$

这里的 $\sigma_{\rm flux}$ 不能替换整个奇子空间的 $\sigma_o$。碰撞输入仍是名义奇偶率；物理通量率只是局部反馈消元后该物理块的等效系数。

`ParameterReport` 不返回裸率元组，而是同时给出：

- `status`：`feasible_exact`、`feasible_restricted`、`no_feasible_solution`、`degenerate_branch`、`boundary_correction_required` 或 `mrt_extension_required`；
- `exact_substitutions`：精确有理式、根式和极限；
- `collision_rates` 与逐块 `open_interval_checks`；
- `violated_constraints`、全部推导假设和 `minimal_extension`；
- `consumed_evidence`：Task 5 的规范 `D2Q9EquivalentCoefficients` / `QuarticConditionSystem` 对象和 Task 6 的实际 residual/classification 对象，不用字符串冒充推导证据；
- `is_conditional` 与 `feasibility_conditions`：自由符号输入尚不能判定正性或开区间时，显式保存待满足条件。

Task 5 的 `canonical_quartic_condition()` 从受审闭式输入构造保留因子的四阶多项式，并在多个精确参数点同时与 amplification route 和 Taylor route 核对。参数层直接消费该规范对象；包括 $a=1$ 在内都从同一个未除以 $a-1$ 的多项式分类，不再复制一套 reviewed family，也不把某个有理专门化点冒充一般证据。

## Task 5 体相条件

`recover_baseline_quartic_family()` 对 `amplification_route(order=4)` 的输出调用 `quartic_condition_system()`，恢复完整基线族

$$
\sigma_e=\frac{4\sigma_o^2+1}{8\sigma_o},
$$

而不是直接返回某个熟知 TRT 点。

对外置精确梯度源和局部非平衡反馈，Task 5 的受审主分支分别为

$$
\sigma_e^{\rm ext}
=\frac{3ab-a-b-12b^2\sigma_o^2}
{12b\sigma_o(a-1)},
$$

$$
\sigma_e^{\rm fb}
=\frac{3a^2-2a-12b^2\sigma_o^2}
{12b\sigma_o(a-1)}.
$$

它们只适用于 $a\ne1$、$b\ne0$ 和非零主平移。求解器先分类 $a=0$、$b=0$ 与 $a=1$，不会把含 $a-1$ 的表达式称为普适公式。

## 与受限标量 ABB 的精确消元

Task 6 在以下假设内给出唯一可用的标量 ABB 行：D2Q9 两率标量碰撞、平直 grid-aligned halfway wall、稳态一维二次温度场、半热源重构、完整压力平衡墙项，以及对应的 external-gradient 或 local-feedback population chain。在这个适用域内，

$$
\sigma_{\rm flux}\sigma_e=\frac3{16}.
$$

代入 $\sigma_{\rm flux}=K/a$ 得

$$
\sigma_e=\frac{3a}{16K}.
$$

### 局部反馈

把 $\sigma_o=K/b$ 与上式代入 Task 5 反馈主分支，精确约去分母后得到

$$
\boxed{K^2=\frac{a(3a+1)}{48}}.
$$

### 外置精确梯度源

相同消元给出

$$
\boxed{
K^2=\frac{12ab+5a-4b-9a^2}{48}
}.
$$

这两式是“冻结体相四阶完全消除 + 受限 ABB 二次行 + 给定输运”同时成立的条件，不是推荐扩散率。

当 $\pi=0$、$a=c_s^2=1/3$ 时，两式都退化为正物理分支

$$
\frac\kappa{\Delta t}=\frac1{\sqrt{72}}.
$$

相应物理通量与名义偶平移固定为

$$
\sigma_{\rm flux}=\frac{\sqrt2}{4},
\qquad
\sigma_e=\frac{3\sqrt2}{8},
$$

但名义奇平移仍为

$$
\sigma_o=\frac1{2\sqrt2(1-\chi_\kappa)}.
$$

因此必须给出 $\chi_\kappa$ 才能生成名义奇率。正物理分支还要求 $a\ne0$、$b>0$、$K>0$、三个相关 Hénon 平移为正，并逐一检查名义奇率、名义偶率和物理通量率是否属于开区间 $(0,2)$。该区间只表示形式可接受，不表示条件数良好或数值稳健。

## 低扩散率与诚实的约束选择

若 $\Delta t=1$ 且目标 $\kappa=10^{-3}$，一般不能同时满足上述体相四阶条件与受限 ABB 条件。求解器返回 `no_feasible_solution`，并把 `bulk_quartic_and_restricted_abb` 列入违反项；它不会静默保留其中一个条件。

若明确选择“保留体相四阶消除，并把显式墙面修正列为待完成扩展”，低扩散率仍可有可接受的名义率。这里仅证明名义率在形式上可接受，并未给出修正公式或修正后 residual。以反馈分支

$$
a=\frac13,
\qquad
b=\frac14,
\qquad
\Delta t=1,
\qquad
\kappa=\frac1{1000}
$$

为例，Task 5 体相条件给出

$$
\sigma_o=\frac1{250},
\qquad
\sigma_e=\frac{250009}{6000},
\qquad
\sigma_{\rm flux}=\frac3{1000},
$$

$$
s_o=\frac{125}{63},
\qquad
s_e=\frac{6000}{253009},
\qquad
s_{\rm flux}=\frac{1000}{503}.
$$

三个率都在 $(0,2)$ 内，但状态仍是 `boundary_correction_required`。这里的显式修正只释放已经审计过的受限 ABB product，不能自动消除一般 pressure/source/time/force/tangential jets。

工程上有三种彼此独立、不能偷换的选择：

1. 低 $\kappa$ + 体相四阶消除 + 显式墙面修正；
2. 低 $\kappa$ + 受限 ABB 标定 + 保留体相四阶误差；
3. 当两者都必须满足时，把显式受限 ABB 边界修正列为已识别的结构性最小扩展；它与 split-even MRT 一样，仍需另行完成推导、实现和 residual 验证。

显式边界修正与 split-even MRT 不是同一种机制。前者是当前根据违反的受限 ABB product 识别出的结构性最小扩展，但仓库尚无 correction formula 或 corrected residual，因而不能称为充分方案；后者只有在矩空间推导明确“哪个偶模态进入 ABB、哪个独立偶模态进入 $C_{40}$”，并证明双约束 Jacobian 满秩后，才能升级为可行结论。当前二者都不能标成 `feasible_exact`；split-even 只报告 `candidate_requiring_mode_jacobian_derivation`。

## 特殊与退化分支

求解顺序固定为先分类、后除法：

- 数值 $a\le0$：平衡/物理块非正，返回 `no_feasible_solution` 与
  `a_not_positive`。
- 数值 $b\le0$：输运块非正，返回 `no_feasible_solution` 与
  `b_not_positive`。
- 数值 $K\le0$：目标扩散输运非正，返回 `no_feasible_solution`。
- $a=1$、反馈：直接使用未除式条件 $12K^2=1$，不除以 $a-1$。
- $a=1$、外置梯度：直接使用 $1-2b+12K^2=0$；正支要求 $b>1/2$。
- 兼容根号为零：零 Hénon 平移不作为可行 ABB 极限。
- 兼容根号为负：不存在实正物理分支。
- 反馈 $a+2K=0$ 且其余输入为正：局部梯度闭合奇异，返回
  `degenerate_branch`。

数值可行性优先于代数边界分类：形式交点 $a=1/9$、$K=-1/18$ 虽满足 $a+2K=0$，但因 $K<0$ 返回 `no_feasible_solution`。任一实际碰撞率不在 $(0,2)$ 也采用相同优先级。只有自由符号无法判定时才保留条件报告，而不是宣称已经可行。

## 流动侧兼容性

剪切输运给出

$$
\nu=c_s^2\Delta t(1-\chi_s)\sigma_f^+.
$$

同一名义偶 shift 还产生独立体积块

$$
\sigma_{b,\mathrm{eff}}=(1-\chi_b)\sigma_f^+,
\qquad
\nu^B=\frac2D c_s^2\Delta t\,\sigma_{b,\mathrm{eff}}.
$$

二维参数报告因此同时返回 sigma_bulk_physical、physical_bulk 实际率与
nu_bulk_2d，并在任何壁面分类之前检查 $1-\chi_b>0$、体积 shift、$\nu^B$
及体积实际率 $0<s_{b,\mathrm{eff}}<2$。剪切率可接受不能替代这组门禁。

只保留 Task 6 的稳态一维均匀体力剪切行时，

$$
(1-\chi_s)\sigma_f^+\sigma_f^-=\frac3{16}.
$$

对给定正 $\nu$ 精确解得

$$
\sigma_f^+=\frac{\nu}{c_s^2\Delta t(1-\chi_s)},
\qquad
\sigma_f^-=\frac{3c_s^2\Delta t}{16\nu},
$$

$$
s_f^-=\frac{16\nu}{3c_s^2\Delta t+8\nu}
\longrightarrow0^+
\quad(\nu\to0^+).
$$

所以任意正 $\nu$ 下的开区间检查不等于小黏性时具有良好阻尼或条件数。

若保留 trace jets，求解器改为消费 Task 6 的 general residual/classification。$\chi_b\ne\chi_s$ 时，已解析 shear/bulk 两行的 `rate_compatibility_status` 为 `no_single_magic`，返回 `mrt_extension_required`；若选择独立剪切/体积率，还必须同时加入一般速度边界修正，不能只解剪切/体积两行后宣称通用壁面标定。另一条独立充分路线是直接采用一般速度边界修正。即使 $\chi_b=\chi_s$，完整一般壁面表仍保留未闭合 jets，因此顶层状态是 `boundary_correction_required`，而不是 universal magic。

## 结论边界

- $1/\sqrt{72}$ 只是特定约束集合在 $\pi=0$ 的兼容点，不是通用推荐值。
- $3/16$ 只在 Task 6 明示的受限墙面、source/gauge、几何和场假设内使用。
- 冻结系数四阶条件不外推到变系数四阶精度。
- 名义 ghost 率、物理有效率、输运系数和边界校准量始终分栏报告。
- 任何不可兼容的约束集合都返回非空 `minimal_extension`，不静默删约束。

# LBM-CDE-TRT 单步算法规范

本章把前七章的代数结果整理为一个可直接实现的时间步。流程顺序是规范的一部分：宏观重构、反馈量、源项、碰撞和边界不能为方便而交换。生产实现仍应从仓库相同维度、相同后端的最近基线增量修改。

## 输入与启动检查

每个算例必须显式给出

$$
(s_f^+,s_f^-),\qquad(s_g^+,s_g^-),
\qquad
(\chi_s,\chi_b,\chi_\kappa),
\qquad
\Delta t, c_s^2, \rho_0,
$$

以及压力修正、力、热源、速度边界、温度边界和 mixed-corner 策略。启动时执行：

1. 将每个名义率转换为 $\sigma=1/s-1/2$，分别检查 $0<s<2$；剪切、bulk/trace 和标量通量的任一实际数值率越界，或相应物理输运非正，即为物理不可行；
2. 检查 $a=c_s^2+p/\rho_0$、标量输运块 $b$、目标 $K$ 或 $\nu$ 的正性与反馈闭合分母；自由符号无法判定时必须保留条件报告；
3. 用第七章求解器验证请求的体相/壁面约束集合，若状态为 `no_feasible_solution` 或 `degenerate_branch` 则停止；
4. 记录所有 restricted 假设、显式边界修正和 corner 策略，不能只记录数值率。

开区间检查不是稳定性证明；实际网格、Mach 数、浮力、边界和并行实现仍需独立验证。

## 一个时间步的固定顺序

### 用半力、半热源重构 rho、u、T

从 transformed populations $\widetilde f_i,\widetilde g_i$ 重构

$$
\delta\rho=\sum_i\widetilde f_i,
$$

$$
\rho_0\boldsymbol u
=\sum_i\boldsymbol c_i\widetilde f_i
+\frac{\Delta t}{2}\boldsymbol F,
$$

$$
T=\sum_i\widetilde g_i+\frac{\Delta t}{2}Q.
$$

随后更新状态方程、$p/\rho_0$、浮力 $\boldsymbol F$ 与热源 $Q$。若 $\boldsymbol F$ 或 $Q$ 依赖刚重构的宏观量，必须在本步内按已声明的显式/迭代策略统一更新，不能一部分用旧时刻、一部分用新时刻。

### 构造含压力修正的平衡态

构造 $f_i^{eq}$ 与 $g_i^{eq}$，并逐点保证至少满足

$$
\sum_i f_i^{eq}=\delta\rho,
\qquad
\sum_i c_{i\alpha}f_i^{eq}=\rho_0u_\alpha,
$$

$$
\sum_i c_{i\alpha}c_{i\beta}f_i^{eq}
=\delta\rho c_s^2\delta_{\alpha\beta}
+\rho_0u_\alpha u_\beta,
$$

$$
\sum_i g_i^{eq}=T,
\qquad
\sum_i c_{i\alpha}g_i^{eq}=Tu_\alpha,
$$

$$
\sum_i c_{i\alpha}c_{i\beta}g_i^{eq}
=T(c_s^2+p/\rho_0)\delta_{\alpha\beta}
+Tu_\alpha u_\beta.
$$

温度 Dirichlet ABB 的墙面平衡项必须使用完整压力修正权重；不能只在体相平衡态中保留压力、在墙面项中删除它。

### 从局部非平衡重构应变与温度梯度

先形成 transformed nonequilibrium

$$
f_i^{neq}=\widetilde f_i-f_i^{eq},
\qquad
g_i^{neq}=\widetilde g_i-g_i^{eq}.
$$

流动侧把非平衡二阶矩分成非对角剪切、对角偏差和迹三个块。每一块使用第三章的局部闭合

$$
\mathcal S
=-\frac{2\mathcal P+\Delta t\,a_{uF}}
{2\Delta t\rho_0c_s^2[1+2(1-\chi)\sigma_f^+]},
$$

其中剪切/对角偏差用 $\chi_s$，迹块用 $\chi_b$，相应 $uF$ 仿射项不能省略。

标量侧使用

$$
\nabla T
=-\frac{2\boldsymbol j^{neq}+\Delta t\boldsymbol a_T}
{2\kappa+\Delta t(c_s^2+p/\rho_0)},
\qquad
\boldsymbol a_T=\frac{T\boldsymbol F}{\rho_0}+Q\boldsymbol u.
$$

在冻结、静止的第四阶审计子域中，这与

$$
\nabla T
=-\frac{2J(I-G)\widetilde g}
{\Delta t(a+2b\sigma_g^-)}
$$

一致。若 $a+2b\sigma_g^-=0$，反馈分支必须立即停止，不能以极小数替代分母。

### 分别构造 raw even/odd 源

先按物理公式构造 raw source，再用 $P_\pm$ 投影复核奇偶性。流场源为

$$
S_i^-=w_i\frac{\boldsymbol c_i\cdot\boldsymbol F}{c_s^2},
$$

$$
S_i^+=w_i\left[
\frac{(\boldsymbol c_i\cdot\boldsymbol u)
(\boldsymbol c_i\cdot\boldsymbol F)}{c_s^4}
-\frac{\boldsymbol u\cdot\boldsymbol F}{c_s^2}
+\rho_0H_i^{(2)}:A
\right],
$$

$$
A_{\alpha\beta}
=\chi_sS_{\alpha\beta}
+\frac{\chi_b-\chi_s}{2}
S_{\gamma\gamma}\delta_{\alpha\beta}.
$$

标量源为

$$
R_i^+=w_iQ,
$$

$$
R_i^-=w_i\left[
\frac{\boldsymbol c_i\cdot(p\nabla T+T\boldsymbol F)}
{\rho_0c_s^2}
+Q\frac{\boldsymbol c_i\cdot\boldsymbol u}{c_s^2}
+\chi_\kappa\boldsymbol c_i\cdot\nabla T
\right].
$$

必须维护独立数组或独立表达式 `source_even`、`source_odd`；不能先用一个率乘总源，再试图恢复奇偶因子。

### 用名义奇偶率和各自源因子碰撞

对 $h=\widetilde f$ 或 $\widetilde g$，统一执行

$$
\begin{aligned}
\widetilde h_i^*={}&\widetilde h_i
-s_+P_+(\widetilde h-h^{eq})_i
-s_-P_-(\widetilde h-h^{eq})_i\\
&+\Delta t\left[
\left(1-\frac{s_+}{2}\right)q_i^+
+\left(1-\frac{s_-}{2}\right)q_i^-
\right].
\end{aligned}
$$

这里始终使用四个名义输入率 $(s_f^+,s_f^-,s_g^+,s_g^-)$。反馈后不同于名义 parity ghost 的物理块仅为

$$
\sigma_{f,\mathrm{shear}}^{\rm phys}
=(1-\chi_s)\sigma_f^+,
\qquad
\sigma_{f,\mathrm{bulk}}^{\rm phys}
=(1-\chi_b)\sigma_f^+,
$$

$$
\sigma_{g,\mathrm{flux}}^{\rm phys}
=\frac{(1-\chi_\kappa)c_s^2}{c_s^2+p/\rho_0}\sigma_g^-.
$$

流动 even/odd ghosts 仍使用 $\sigma_f^+,\sigma_f^-$；标量 even/odd ghosts 仍使用 $\sigma_g^+,\sigma_g^-$。物理有效率不是额外碰撞输入，也不能覆盖整个 parity 子空间。

BGK 调试开关只做

$$
s_f^+=s_f^-,
\qquad
s_g^+=s_g^-,
$$

于是每套分布逐分量退化为

$$
\widetilde h_i^*=\widetilde h_i-s(\widetilde h_i-h_i^{eq})
+\Delta t(1-s/2)q_i.
$$

BGK 开关不自动令 $\chi_s=\chi_b=\chi_\kappa=0$；若要验证无反馈 BGK 基线，必须另行关闭反馈参数。

### Stream

执行 pull 或 push streaming，但全程序只能选择一种索引约定：

$$
\widetilde h_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
=\widetilde h_i^*(\boldsymbol x,t).
$$

MPI/OpenACC 版本必须在同一逻辑位置完成 halo/device 同步；同步只改变数据可见性，不改变碰撞—stream—边界的数学顺序。

### 依次施加速度 BB、温度 ABB、绝热 BB 和显式 corner rule

1. 对固壁速度 link 应用所选 halfway BB。只有平直、格点对齐的 halfway 壁、稳态一维 Stokes 剪切、均匀体力、half-force gauge、冻结 feedback 且无 bulk/time/tangential jets 的受限子域，才可使用 $(1-\chi_s)\sigma_f^+\sigma_f^-=3/16$ 标定；pressure-driven gauge 的 $3/8$ 不与它混用。
2. 对 Dirichlet 温度 link 应用 ABB，墙项保留

   $$
   2\left(w_i+\lambda_i\frac{p_w/\rho_0}{c_s^2}\right)T_w.
   $$

3. 对绝热 link 应用 homogeneous Neumann/BB；不要把 uniform heat source 重复写入 odd-flux 边界行。已解析行限定于平直、格点对齐的 halfway 壁、half-source、transformed CDE chain 和 D2Q9 $c_s^2=1/3$；其中 affine-normal 还要求稳态、常压力比和零切向 jets。它们不代表一般闭合：一般 wall-time、normal-curvature 和 source rows 仍为 `unresolved`，需要另行边界修正。
4. 在格点对齐直角、half-source 规范的 mixed Dirichlet/adiabatic corner，对共享 diagonal population 只赋值一次。默认显式策略为：轴向 links 各自执行所属边界，共享 diagonal 采用 Dirichlet-priority ABB，跳过随后对同一 population 的 adiabatic overwrite，并把未满足的绝热 corner residual 记入诊断。若项目选择别的单一闭合，必须替换整条 corner 策略，而不是交换两次 overwrite 的顺序。

共享 corner source increment 只计一次，diagonal wall distance 使用 $h/\sqrt2$。一般 mixed corner 的两条方程秩不兼容，因此这里明确选择工程闭合，不宣称同时精确满足两类边界。

### 计算诊断与不变量

每步或按固定间隔检查：

- 质量、动量和标量总量，以及施加源后的预算残差；
- $\rho,T,p,\boldsymbol u$ 和所有 populations 是否有限；
- $a$、$b$ 与反馈闭合分母是否跨过奇点；
- 名义率、物理有效率、输运系数是否仍与启动参数映射一致；
- 边界 residual、corner overwrite 计数和共享 source 计数；
- `errorU`、`errorT`、Nu、Re、壁面极值及算例规定的收敛判据。

冻结系数四阶条件只在其审计子域内作为诊断；若 $p/\rho_0$、松弛率或 $\chi$ 空间变化，必须另行记录乘积导数残差，不能继续标记为四阶完全消除。

## 三种标量运行策略

配置层必须显式选择以下一种，不允许求解器暗中替换：

### 低 kappa + 体相四阶消除 + 显式边界修正

名义率满足 Task 5 四阶条件；求解器把 Dirichlet 墙面的显式修正报告为结构性最小扩展。本文尚未给出 correction formula 或 corrected residual，因此该分支状态仍是 `boundary_correction_required`，不是已实现方案，更不是一般 source/time/pressure jets 的通用修复。

### 低 kappa + 受限 ABB 标定 + 保留体相 q4 误差

仅在平直、格点对齐的 halfway 壁、稳态一维二次温度场、half-source 重构、完整 pressure-wall equilibrium term，以及相匹配的 external-gradient 或 local-feedback population chain 下，用 $\sigma_{\rm flux}\sigma_e=3/16$ 生成名义偶率。随后计算并报告非零 $C_{40},C_{22}$；不得把“各向同性残差为零”写成“四阶误差为零”。

### 两项都强制时采用独立新推导

显式受限 ABB 边界修正是已识别的结构性最小扩展，但尚待 correction formula、实现与 corrected-residual 验证。split-even MRT 是增加独立偶模态的另一候选；在明确 ABB 与 $C_{40}$ 分别依赖哪个偶模态、并证明两约束 Jacobian 满秩前，只能标记 `split_even_mrt_candidate_requiring_mode_jacobian_derivation`。两条路线当前都不是已证明充分方案，也不能预先标记 `feasible_exact`。

## 最小伪代码

```text
validate_parameter_report()
for each time step:
    reconstruct_rho_u_T_with_half_sources()
    update_pressure_force_heat()
    build_pressure_corrected_equilibria()
    reconstruct_strain_and_temperature_gradient()
    build_raw_even_odd_sources_separately()
    collide_with_nominal_parity_rates_and_source_factors()
    stream_and_exchange_halos()
    apply_velocity_BB()
    apply_temperature_ABB()
    apply_adiabatic_BB()
    apply_one_explicit_corner_closure()
    compute_diagnostics_and_invariant_checks()
```

该顺序同时适用于 CPU、OpenMP、OpenACC 和 MPI+OpenACC 后端；后端可以改变循环组织和数据移动，但不能改变物理步骤的依赖关系。

# 参数化 D2Q5 四阶参考验证器

本附录只验证 Dubois--Lallemand 的 D2Q5 纯扩散模型。它为后续 D2Q9 推导校验符号工具，但任何 D2Q5 四阶系数都不能直接迁移到 D2Q9。

## 原始模型与参数化平衡态

速度按静止、东、北、西、南排列。Dubois Appendix Eq. (79) 的矩阵为

$$
M=\begin{pmatrix}
1&1&1&1&1\\
0&\lambda&0&-\lambda&0\\
0&0&\lambda&0&-\lambda\\
-4&1&1&1&1\\
0&1&-1&1&-1
\end{pmatrix}.
$$

纯热模型只有 $m_0=\rho$ 守恒，平衡矩是

$$
m^{eq}=\rho(1,0,0,\alpha,0)^{\mathsf T}.
$$

直接计算 $f^{eq}=M^{-1}m^{eq}$，得到一般参数权重

$$
w_0=\frac{1-\alpha}{5},
\qquad
w_1=w_2=w_3=w_4=\frac{4+\alpha}{20}.
$$

因此

$$
\sum_i w_i=1,
\qquad
\sum_iw_i c_{i\alpha}=0,
\qquad
\sum_iw_i c_{i\alpha}c_{i\beta}
=c_e\lambda^2\delta_{\alpha\beta},
$$

其中

$$
\boxed{c_e=\frac{4+\alpha}{10}}.
$$

Task 1 的固定 D2Q5 权重只对应 $\alpha=-2/3$：此时
$w_0=1/3$、$w_{1\ldots4}=1/6$、$c_e=1/3$。验证器始终先保留一般 $\alpha$，不会用这个固定点代替参数化模型。

Dubois Eq. (39) 的碰撞矩阵是

$$
\Psi=\begin{pmatrix}
1&0&0&0&0\\
0&1-s_1&0&0&0\\
0&0&1-s_1&0&0\\
\alpha s_3&0&0&1-s_3&0\\
0&0&0&0&1-s_4
\end{pmatrix},
\qquad
\sigma_j=\frac1{s_j}-\frac12.
$$

## 论文 Eq. (40) 的符号与尺度

原文 Eq. (40) 为

$$
\begin{aligned}
\partial_t\rho
&-\frac{\lambda^2\Delta t}{10}\sigma_1(4+\alpha)
(\partial_x^2+\partial_y^2)\rho\\
&+\frac{\Delta t^3\lambda^4}{1200}\sigma_1(4+\alpha)
\left[
\kappa_{40}(\partial_x^4+\partial_y^4)
+\kappa_{22}\partial_x^2\partial_y^2
\right]\rho
=O(\Delta t^4).
\end{aligned}
$$

所以扩散率和两个有量纲四阶系数分别是

$$
D=c_e\lambda^2\Delta t\,\sigma_1,
$$

$$
C_{40}=\frac{\Delta t^3\lambda^4}{1200}
\sigma_1(4+\alpha)\kappa_{40},
\qquad
C_{22}=\frac{\Delta t^3\lambda^4}{1200}
\sigma_1(4+\alpha)\kappa_{22}.
$$

代码内部先取 $\lambda=\Delta t=1$ 完成精确递推，再只在出口恢复
$D\mapsto\lambda^2\Delta tD$、$C_4\mapsto\lambda^4\Delta t^3C_4$。这是代数尺度分离，不是修改模型或把论文闭式代入递推。

## Route A：放大矩阵的流体根级数

原文第 17 页采用平面波 $f\propto\exp(i k_xx+i k_yy)$，并明确给出正相位

$$
B=\operatorname{diag}
(1,e^{iq_x},e^{iq_y},e^{-iq_x},e^{-iq_y}),
\qquad
q_\alpha=\lambda\Delta t\,k_\alpha,
$$

$$
G=B M^{-1}\Psi M.
$$

这里需要区分原文内部的两套符号写法：Dubois Eq. (15) 的空间迁移符号与第 17 页印出的正相位 $B$ 的空间位移符号相反。本验证器固定采用第 17 页的正相位约定；若改用 Eq. (15) 的约定，等价于作 $\boldsymbol k\mapsto-\boldsymbol k$。当前问题是纯扩散，本文核对的偶数阶系数因此不变。Task 5 一旦加入源项或平流项，相关奇数阶与源项符号不能由这里直接沿用，必须在其自行固定的 Fourier 约定下重新推导。

验证器不调用浮点特征值求解器生成符号答案。令

$$
G(\epsilon)=\sum_{n=0}^4\epsilon^nG_n,
\quad
z_h(\epsilon)=1+\sum_{n=1}^4\epsilon^nz_n,
\quad
v(\epsilon)=v_0+\sum_{n=1}^4\epsilon^nv_n,
$$

其中 $v_0=f^{eq}/\rho$，并选规范
$\boldsymbol 1^{\mathsf T}v_n=0$。第 $n$ 阶由同一个增广系统得到：

$$
\begin{pmatrix}
G_0-I&-v_0\\
\boldsymbol1^{\mathsf T}&0
\end{pmatrix}
\begin{pmatrix}v_n\\z_n\end{pmatrix}
{}={}
\begin{pmatrix}
-\displaystyle\sum_{r=1}^nG_rv_{n-r}
+\displaystyle\sum_{r=1}^{n-1}z_rv_{n-r}\\
0
\end{pmatrix}.
$$

最后展开

$$
\Gamma=-\frac1{\Delta t}\log z_h.
$$

奇数总次数严格为零；二次和四次齐次项分别给出 $D,C_{40},C_{22}$。

## Route B：物理空间 Taylor/矩递推

第二条路线不构造 $G$，不计算特征多项式、特征向量、$z_h$ 或对数。它直接从与正相位约定一致的离散迁移恒等式

$$
f_i(\boldsymbol x,t+\Delta t)
=f_i^*(\boldsymbol x+\boldsymbol c_i\Delta t,t)
$$

作物理空间 Taylor 展开：

$$
\sum_{p=0}^4\frac{\Delta t^p}{p!}\partial_t^pf_i
{}={}
\sum_{p=0}^4\frac{\Delta t^p}{p!}
(\boldsymbol c_i\cdot\nabla)^pf_i^*+O(\Delta t^5).
$$

写成 differential-jet 级数

$$
m=m^{eq}+\sum_{n=1}^4\epsilon^nm^{[n]},
\qquad
\partial_t\rho
=\sum_{r=0}^3\epsilon^rL_r(\partial_x,\partial_y)\rho.
$$

已知低阶 jet 后，第 $n$ 阶 Taylor 残差记为 $R_n$。未知的非守恒矩修正和守恒等效方程系数通过

$$
\boxed{
(I-\Psi)m^{[n]}+m^{eq}L_{n-1}=-R_n,
\qquad m_0^{[n]}=0
}
$$

逐阶消去。结果是

$$
L_0=0,
\qquad
L_1=c_e\sigma_1(\partial_x^2+\partial_y^2),
\qquad
L_2=0,
$$

而 $L_3$ 给出 Eq. (40) 的两个四阶项。Route B 只与 Route A 共享 $M$、$\Psi$、速度和通用截断多项式工具；测试还把 Route A 两个入口替换成抛异常函数，Route B 仍独立生成二阶系数。

## 两路线对 Eqs. (41)--(42) 的精确复现

两条递推的差经 SymPy 化简严格为零，并分别得到

$$
\boxed{
\kappa_{40}=8-3\alpha
+12(\alpha+4)\sigma_1^2
-12(1-\alpha)\sigma_1\sigma_3
-60\sigma_1\sigma_4
}
$$

以及

$$
\boxed{
\kappa_{22}=-6(\alpha+4)
+24(\alpha+4)\sigma_1^2
-24(1-\alpha)\sigma_1\sigma_3
+120\sigma_1\sigma_4.
}
$$

这两个式子是生成器输出与论文闭式的比较结果；Route A 和 Route B 的递推本身没有调用它们。

## 各向同性不等于完整消除

四阶算子与 $\nabla^4$ 各向同性一致只要求
$\kappa_{22}=2\kappa_{40}$。精确残差是

$$
\boxed{
\kappa_{22}-2\kappa_{40}
=40(6\sigma_1\sigma_4-1).
}
$$

所以各向同性关系为

$$
\sigma_1\sigma_4=\frac16.
$$

它只使两个四阶方向系数满足 $1:2$，一般并不使
$\kappa_{40}$ 或 $\kappa_{22}$ 为零。完整四阶消除必须同时解
$\kappa_{40}=\kappa_{22}=0$，得到 Dubois Eq. (55)：

$$
\boxed{
\sigma_3=\sigma_1\frac{\alpha+4}{1-\alpha}
-\frac{2+3\alpha}{12\sigma_1(1-\alpha)},
\qquad
\sigma_4=\frac1{6\sigma_1}.
}
$$

若再施加中间 TRT 约束 $\sigma_3=\sigma_4$，则

$$
\boxed{
\sigma_1=\frac1{\sqrt{12}},
\qquad
\sigma_3=\sigma_4=\frac1{\sqrt3}.
}
$$

代回后两个系数对任意允许的 $\alpha$ 都为零。

## Wang/Luo 符号映射

Wang Eq. (12) 与 Luo Eq. (14) 在 $\lambda=1$ 时使用相同的 D2Q5 矩阵；两文的自由参数 $a$ 对应本附录的 $\alpha$。两文的

$$
Q=\operatorname{diag}(\cdots,\sigma_\kappa,\sigma_\kappa,
\sigma_e,\sigma_\nu)
$$

中的 $\sigma_\kappa,\sigma_e,\sigma_\nu$ 是碰撞率，不是 Dubois 的 Hénon 平移。映射必须写成

$$
s_1=\sigma_\kappa,
\qquad s_3=\sigma_e,
\qquad s_4=\sigma_\nu,
$$

$$
\sigma_1=\frac1{\sigma_\kappa}-\frac12,
\quad
\sigma_3=\frac1{\sigma_e}-\frac12,
\quad
\sigma_4=\frac1{\sigma_\nu}-\frac12.
$$

Wang/Luo 的各向同性式因此正好对应
$\sigma_1\sigma_4=1/6$。该关系单独使用时不是完整四阶消除条件，更不是 D2Q9 条件。

## 精确与高精度抽查

一般符号比较之外，测试还使用三组允许范围内的有理参数：

| $\alpha$ | $\sigma_1$ | $\sigma_3$ | $\sigma_4$ | 轴向 $q^4$ 误差 | 对角 $q^4$ 误差 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| $-1$ | $1/5$ | $2/7$ | $3/11$ | $2.10\times10^{-16}$ | $1.41\times10^{-16}$ |
| $0$ | $1/3$ | $1/4$ | $1/6$ | $1.10\times10^{-15}$ | $3.79\times10^{-16}$ |
| $1/2$ | $2/5$ | $3/8$ | $4/9$ | $4.84\times10^{-16}$ | $9.91\times10^{-16}$ |

抽查使用 80 位精度、$q=10^{-6}$。轴向比较

$$
\frac{\Gamma(q,0)-Dq^2}{q^4}\to C_{40},
$$

对角比较

$$
\frac{\Gamma(q,q)-2Dq^2}{q^4}
\to2C_{40}+C_{22}.
$$

数值特征值只用于独立抽查，没有参与符号系数的生成。

## 结论边界

- 本附录验证的是参数化 D2Q5 纯扩散模型。
- D2Q5 不是通过删除 D2Q9 对角速度得到的连续参数退化。
- Wang/Luo 各向同性式不等于完整四阶消除。
- Eq. (55) 和 TRT 特殊点均为 D2Q5-only 结果，不能移植到 Task 5 的 D2Q9 目标模型。

# 需求—证据—验证矩阵

本矩阵把每个主要主张连接到报告位置、原始论文依据、项目内生成入口和可执行
测试。论文只提供低阶定义或受限边界结果时，项目内代数不会反向冒充论文原式。

## R1：D2Q9 与 LBM-CDE 压力修正向量

- 结论等级：`strictly_proved`。
- 报告位置：总报告“记号”直接列出九分量向量；完整离散求和见 evidence ledger“本任务采用的精确离散常数”；第 06 章第 4 节使用该向量。
- 论文依据：LBM-CDE Eqs. (16)--(17)、Eq. (36)；Wang Eq. (1)。
- 项目生成：`d2q9()`、`raw_moment()`。
- 可执行验证：`test_d2q9_exact_definition`、`test_d2q9_lambda_t_constraints`、
  `test_d2q9_lambda_t_fourth_raw_moments`。
- 限制：$\lambda_i$ 的四阶矩是离散求和结果，不是论文声明为零的量。

## R2：奇偶源与 transformed TRT 碰撞

- 结论等级：`strictly_proved`。
- 报告位置：第 02 章第 2--6 节。
- 论文依据：LBM-CDE Eqs. (13)、(24)、(26)、(28)，Appendix A.1--A.5。
- 项目生成：`flow_source()`、`scalar_source()`、`trt_collision()`。
- 可执行验证：`test_sources_are_exact_even_odd_projections`、
  `test_trt_collision_matches_operator_trapezoidal_formula`、
  `test_one_flow_collision_has_exact_net_momentum_source`、
  `test_one_scalar_collision_has_exact_net_heat_source`。
- 限制：TRT 双因子是项目内投影算子推导，不是原文 BGK 分量式的逐字抄录。

## R3：逐分量 BGK 极限

- 结论等级：`strictly_proved`。
- 报告位置：第 02 章第 7 节；第 08 章 BGK switch。
- 论文依据：LBM-CDE Eqs. (26)、(28)。
- 项目生成：TRT 碰撞对象在 $s_+=s_-$ 下的符号差。
- 可执行验证：`test_equal_rates_recover_bgk_componentwise`、
  `test_equal_nominal_rates_and_zero_feedback_recover_bgk_blocks`。

## R4：局部反馈后的物理有效 shift

- 结论等级：`strictly_proved`（冻结系数、CE2）。
- 报告位置：第 03 章；第 04 章第 4--7 节。
- 论文依据：LBM-CDE Eqs. (30)--(35)。
- 项目生成：`shear_effective_rate()`、`bulk_effective_rate()`、`scalar_flux_effective_rate()`、`effective_operator_blocks()` 与直接碰撞代回。
- 可执行验证：`test_scalar_flux_feedback_is_reduced_from_direct_collision`、
  `test_off_diagonal_and_deviatoric_diagonal_use_shear_block`、
  `test_trace_mode_uses_bulk_block_and_dimension_normalization`。
- 限制：物理 block shift 不替换同奇偶 ghost shift；变系数乘积导数另列 residual。

## R5：二阶宏观恢复与首个遗漏项族

- 结论等级：`strictly_proved`（声明的 CE/Mach 阶次内）。
- 报告位置：第 04 章。
- 论文依据：LBM-CDE Eqs. (1)--(14)、(18)--(29) 与 Appendix A。
- 项目生成：`second_order_residual_table()` 及逐矩源扰动。
- 可执行验证：`test_default_named_moments_generate_all_second_order_cancellations`、
  `test_removing_each_inverse_design_source_moment_exposes_residual`、
  `test_wrong_parity_specific_trapezoidal_factors_expose_residuals`。
- 限制：$O(Ma^3)$、$O(Ma^4)$ 与空间变系数余项不是零。

## R6：classical halfway velocity magic

- 结论等级：`restricted_model`。
- 报告位置：第 05 章第 2 节。
- 论文依据：Ginzburg--d'Humières Eq. (41)--(43)。
- 项目生成：`velocity_bb_residual()` 的 multireflection 映射。
- 可执行验证：`test_multireflection_quarter_maps_exactly_to_three_sixteenths`、
  `test_quarter_and_three_sixteenths_are_not_naked_constants`。
- 限制：平壁、格点对齐、halfway、稳态 Stokes、均匀体力、half-force gauge、无反馈。

## R7：一般速度边界无单一 magic

- 结论等级：`unresolved`（完整 wall table）；`restricted_model`（shear/bulk 子系统）。
- 报告位置：第 05 章第 3、7 节。
- 论文依据：LBM-CDE 流源与局部应变闭合；边界 kinetic rule。
- 项目生成：`velocity_bb_residual(force_mode="general_source_aware")`、
  `classify_magic()`。
- 可执行验证：`test_general_velocity_table_retains_independent_unresolved_wall_jets`、
  `test_equal_shear_and_bulk_feedback_has_one_restricted_rate_calibration`、
  `test_distinct_shear_and_bulk_feedback_has_no_single_magic`。
- 限制：Poiseuille 的 $3/16$、$3/8$ 仅为分别标注 gauge 的 `source_evidence_only`。

## R8：温度 ABB 受限二次行

- 结论等级：`restricted_model`。
- 报告位置：第 05 章第 4 节。
- 论文依据：LBM-CDE Eqs. (22)--(24)、(34)--(38)。
- 项目生成：external D1Q3 population solve 与 local-feedback homogeneous solve。
- 可执行验证：`test_quadratic_abb_chain_maps_physical_flux_times_nominal_even_ghost`、
  `test_direct_quadratic_temperature_stencil_is_independent_and_agrees`。
- 限制：平直、格点对齐的 halfway 壁，冻结压力、稳态一维二次场、half-source 重构、完整 pressure-wall equilibrium term、匹配的 external/local population chain、CDE-consistent 均匀 $Q$，且无流/力/切向 jets。

## R9：绝热与混合角点

- 结论等级：`restricted_model`（已解析行）与 `unresolved`（一般 time/normal-curvature/source rows）。
- 报告位置：第 05 章第 5--6 节。
- 论文依据：LBM-CDE Eq. (39) 及 transformed source 定义。
- 项目生成：affine normal 与 quadratic diagonal population solve；finite corner assignments。
- 可执行验证：`test_adiabatic_primary_api_consumes_executable_population_chain`、
  `test_direct_adiabatic_diagonal_pair_exposes_tangential_curvature`、
  `test_one_diagonal_unknown_two_wall_equations_are_rank_incompatible`。
- 限制：绝热已解析行要求 `flat_grid_aligned_halfway`、`half_source`、`transformed_cde_chain`、D2Q9 $c_s^2=1/3$；affine-normal 还要求稳态、常压力比、零切向 jets。一般 wall-time、normal-curvature 与 source rows 保持 `unresolved`。mixed corner 要求 `grid_aligned_right_angle_corner` 与 half-source 规范，采用单一显式闭合并记录未满足 residual。

## R10：D2Q9 二阶扩散与三种实际方案

- 结论等级：`strictly_proved`（冻结系数）。
- 报告位置：第 06 章第 2、5--7 节。
- 论文依据：LBM-CDE Eqs. (16)--(24)、(34)--(35)。
- 项目生成：`amplification_route()` 与 `taylor_moment_route()`。
- 可执行验证：`test_both_routes_recover_k2_transport_in_all_actual_cases`、
  `test_external_and_feedback_match_at_k2_after_rate_specialization`。

## R11：D2Q9 四阶系数与消除条件

- 结论等级：`strictly_proved`（冻结系数 modified equation）。
- 报告位置：第 06 章第 5--9 节。
- 论文依据：论文只供低阶模型定义；Dubois--Lallemand D2Q9 印刷式仅外部审计。
- 项目生成：独立放大矩阵根级数与物理空间 Taylor/矩递推；
  `canonical_quartic_condition()` 保存未除式条件和 provenance。
- 可执行验证：`test_a_generic_rational_points_have_zero_exact_route_residual`、
  `test_z_bidirectional_monkeypatch_keeps_each_route_independent`、
  `test_canonical_undivided_conditions_match_both_routes_at_exact_points`。
- 数值抽查：80 位 `high_precision_directional_fit()`；消除点 residual 从 $q^4$ 转为 $q^6$。

## R12：D2Q5 文献公式复算

- 结论等级：`strictly_proved`（D2Q5-only）。
- 报告位置：附录 A。
- 论文依据：Dubois--Lallemand Eqs. (39)--(42)、(55)、Appendix Eq. (79)；
  Wang/Luo D2Q5 参数式。
- 项目生成：D2Q5 放大矩阵与独立 Taylor 递推。
- 可执行验证：`test_routes_match_eq41_eq42_symbolically`、
  `test_isotropy_is_distinct_from_complete_cancellation`、
  `test_eq55_and_intermediate_trt_point_cancel_both_coefficients`。
- 限制：任何 D2Q5 条件都不作为 D2Q9 系数输入。

## R13：体相四阶消除与受限 ABB 的兼容恒等式

- 结论等级：`strictly_proved`（两组受限条件的精确消元）。
- 报告位置：第 07 章第 2--3 节。
- 论文依据：Task5/Task6 已审对象；无论文直接给出 LBM-CDE--TRT 联立式。
- 项目生成：`derive_scalar_compatibility()` 消费 Task5 canonical 未除式对象和 Task6
  boundary residual/classification。
- 可执行验证：`test_external_and_feedback_identities_are_generated_by_exact_elimination`、
  `test_derive_retains_the_exact_task5_canonical_objects_it_consumes`。
- 限制：消元中的 ABB product 只来自 R8 所列的平直 halfway、稳态一维二次、half-source、完整墙面平衡项与匹配 source-chain 子域。

## R14：pi=0 的 K=1/sqrt(72) 与低扩散率冲突

- 结论等级：`strictly_proved`（正支、受限 ABB）；不是普适推荐值。
- 报告位置：总报告主要结论 5；第 07 章第 3--4 节。
- 项目生成：canonical 未除式 quartic condition 与 ABB product 的消元。
- 可执行验证：`test_pi_zero_positive_compatible_point_fixes_physical_and_even_shifts`、
  `test_low_positive_kappa_is_rejected_when_both_constraints_are_mandatory`。
- 限制：$1/\sqrt{72}$ 是体相四阶条件与 R8 受限 ABB product 联立的正支，不外推到一般边界或任意低扩散率。

## R15：参数可行性、退化分支与实现顺序

- 结论等级：`strictly_proved`（代数分类）；稳定性只作必要率区间检查。
- 报告位置：第 07 章第 4--7 节；第 08 章。
- 项目生成：`solve_scalar_parameters()`、`solve_flow_parameters()` 与结构化
  `ParameterReport`。
- 可执行验证：`test_negative_b_is_rejected_by_derive_and_target_solver`、
  `test_negative_nu_precedes_general_trace_classification`、
  `test_a_one_special_branches_are_derived_without_dividing_by_a_minus_one`、
  `test_low_kappa_bulk_quartic_with_explicit_wall_correction_is_feasible`、
  `test_a_one_split_even_extension_is_explicitly_only_a_candidate`、
  `test_unequal_shear_and_bulk_feedback_rejects_one_magic_with_trace_jets`、
  `test_negative_bulk_transport_precedes_restricted_or_general_wall_status`。
- 限制：名义率、剪切、bulk/trace 与标量通量实际率必须分别通过门禁；$0<s<2$ 不等于稳定性或良好条件数。显式 ABB 修正只是已识别但尚未推导、实现和验证的结构性最小扩展；split-even MRT 也是待推导候选。一般速度 wall jets 仍要求独立边界修正。

## 交叉审查状态

- TRT/source/CE 代数：Reviewer A 最终签核，`critical=0`、`important=0`、`minor=0`；bulk/trace 负输运反例已纳入 25 项参数回归。
- 边界假设、常数与分类：Reviewer B 最终签核，`critical=0`、`important=0`；全部过度主张与适用域遗漏已修复并复审。
- D2Q5/D2Q9 高阶、兼容性与程序输出：Reviewer C 最终签核，`critical=0`、`important=0`、`minor=0`；独立回归为 D2Q5 19/19、D2Q9 38/38、参数 25/25。
- 三条审查线的全部 `critical`、`important` 与 `minor` 发现均已修复并复审；最终 PDF 另经 60/60 页渲染检查。

# 参考文献与来源说明

1. LBM-CDE 原始论文，用户提供文件 `LBM-CDE.pdf`。本文使用 Eqs. (1)--(35)、
   Eq. (36)--(39) 及 Appendix A 的 transformed trapezoidal LBE；逐页符号映射见
   `derivation/evidence-ledger.md`。

2. Wang 等，*Lattice Boltzmann simulations of thermal convective flows in two
   dimensions*。本文只把其 D2Q5 参数化温度矩与边界 kinetic rules 用作限定来源，
   不把 D2Q5 四阶条件外推到 D2Q9。

3. Contrino 等 / Luo，*Lattice Boltzmann simulations of the thermally driven
   2D square cavity at high Rayleigh numbers*，Journal of Computational Physics
   (2014)。本文审计其 D2Q5 四阶参数与 pressure-boundary Poiseuille 的 $3/8$
   标定，并与 uniform-body-force 的 $3/16$ 分开记录。

4. Dubois 与 Lallemand，*Towards higher order lattice Boltzmann schemes*。
   本文复算其 D2Q5 Eqs. (39)--(42)、(55) 与 Appendix Eq. (79)，并逐字审计
   第 12 页 D2Q9 印刷式；这些公式不作为 LBM-CDE--D2Q9 待求系数的输入。

5. Ginzburg 与 d'Humières，*Multireflection boundary conditions for lattice
   Boltzmann models*，Physical Review E 68, 066614 (2003)。本文使用其
   Eqs. (41)--(43) 建立 $\Lambda^2=1/4$ 与受限 Hénon product $3/16$ 的映射。

所有页码、方程号、允许主张与限制均在证据账本逐条登记。项目内 exact algebra
的入口与测试名见需求—证据—验证矩阵；`Xs/` 不在理论来源列表中。
