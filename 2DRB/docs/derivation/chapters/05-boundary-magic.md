# 05 源感知边界残差与 magic 分类

本章不把某个松弛率乘积先验称为“magic”。判定顺序固定为：先写 transformed population 的碰撞后链，保留独立墙面 Taylor jets，再检查同一条件能否同时消去所有非零系数。`classify_magic()` 同时返回系数表、假设、参数映射、未满足 jets 和已显式相消的零行。对于 general velocity，返回值还包含 `rate_compatibility_status`：它只审计已解析的 shear/bulk 两行，不能代替完整墙面闭合状态。

## 1. 记号、率映射和六个适用域

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

## 2. Classical velocity：1/4 与 3/16 是同一受限标定

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

## 3. Flow feedback：物理块与 nominal ghost 不可混用

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

## 4. Temperature Dirichlet ABB：两条独立 population 路线

### 4.1 墙面平衡项与 transformed source

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

### 4.2 显式精确梯度源链

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

### 4.3 局部非平衡反馈消元链

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

### 4.4 主生成路线与一般 ABB 独立 jets

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

## 5. Adiabatic BB：先求 kinetic odd flux

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

## 6. Mixed corner：覆盖顺序不是同时满足

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

## 7. 结论边界

- `universal_magic`：本章没有找到满足该定义的条件。
- `restricted_calibration`：只用于假设完整写出的 classical、1D shear 或 1D scalar quadratic 行。
- `no_single_magic`：在本章作为 `rate_compatibility_status`，表示已解析 shear/bulk 子系统要求不兼容条件；它不覆盖完整 general wall 状态。
- `boundary_correction_required`：率无关或不受同一 product 控制的 source/time/pressure/tangential 行仍在。
- `corner_closure_conflict`：共享 unknown 的方程秩不兼容或 sequential overwrite 顺序相关。

特别地，$1/4$、$3/16$、$3/8$ 均不得脱离各自 normalization、驱动/gauge 和 geometry 单独引用；D2Q5 bulk constants 也不得移植为 D2Q9 boundary magic。
