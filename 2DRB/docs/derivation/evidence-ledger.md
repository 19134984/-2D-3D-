# LBM-CDE/TRT 推导证据账本

本账本只登记已在原始 PDF 中逐页目视核对的内容。`tmp/pdfs/` 中的文本提取只可用于定位，不作为证据；`Xs/` 是既有复现代码，只可用于实现对照，不能替代论文中的理论依据。

核对日期：2026-07-13。

## 证据使用规则

- “允许主张”限定了后续推导可以从该条证据引用的结论；超出该栏的结论必须另找证据或给出显式代数推导。
- “D2Q5-only” 条目不得直接外推到 D2Q9。
- 论文中的自由参数、连续权函数约束和本项目选定的离散系数必须分开记录。
- 页码优先给出论文印刷页；同时给出 PDF 页序，避免期刊页码与文件页序混淆。

## 1. LBM-CDE.pdf

原始 PDF：`D:\桌面\代码\代码\热对流\2DRB\pdf\LBM-CDE.pdf`

| 印刷页 / PDF 页 | 方程或位置 | 符号映射 | 允许主张 | 限制 |
| --- | --- | --- | --- | --- |
| 第 4 页 / PDF 第 4 页 | Eq. (2) 及其后定义 | 论文 `S_alpha_beta` 对应 `flow_source(..., strain=S)`；`S_gamma_gamma=tr(S)`；`D=2` | 应变率定义为 `S_alpha_beta=(partial_alpha u_beta+partial_beta u_alpha)/2`，迹为速度散度；这是 Eq. (13) 中张量 `A` 的输入约定。 | 该条只固定应变率的正号和 `1/2` 归一化，不提供离散源项。 |
| 第 4 页 / PDF 第 4 页 | Eq. (5) | 论文的连续速度 `v` 对应离散速度 `c_i`；`omega(v)` 对应离散权重 `w_i`；`c_s^2` 对应 `Lattice.cs2` | 目标权函数满足零阶归一、一阶为零、二阶各向同性以及四阶各向同性矩关系。 | Eq. (5) 是连续高斯积分约束，不能单独证明某组离散权重。 |
| 第 6 页 / PDF 第 6 页 | Eq. (13) | `v` 对应 `c_i`；`omega` 对应 `w_i`；`H_i^(2)=c_i c_i/c_s^2-I`；`A=chi_s S+(chi_b-chi_s)tr(S)I/D` | 流场源项的二阶 Hermite 收缩是正号 `+rho_0 w_i H_i^(2):A`。本项目把 `/c_s^2` 放在无量纲 `H_i^(2)` 内，因此 `A` 不再含额外 `/c_s^2`。 | Eq. (13) 是连续 Hermite 形式；离散 D2Q9 三、四阶 raw moments 必须另行精确求和，不能假定消失。 |
| 第 7 页 / PDF 第 7 页 | Eqs. (16)-(17) 及 Eq. (16) 后正文 | 论文 `lambda(v)` 对应 `lambda_t[i]`；`g^eq` 对应温度平衡分布；`T p/(rho_0 c_s^2)` 是该修正的二阶矩系数 | 在论文列出的低阶矩构造中，`lambda` 的零阶、一阶矩为零，二阶矩为 `c_s^2 delta_alpha_beta`。正文说明离散化时相对 `omega` 修改零速度权因子。 | “二阶矩修正”只限定论文规定的低阶作用，不表示所选离散向量的高阶 raw moments 为零。论文也没有逐项印出本项目的 9 元 `lambda_t` 向量；该向量必须由这些约束与标准 D2Q9 权重显式推导。 |
| 第 7 页 / PDF 第 7 页 | Eq. (18) | 论文 `R` 对应标量源 `scalar_source().raw`；`Q` 对应 `heat_source` | 标量源的目标零阶矩是 `integral R dv=Q`。 | 这是连续矩约束；离散奇偶分量的矩需结合 Eq. (24) 和 D2Q9 权重验证。 |
| 第 8 页 / PDF 第 8 页 | Eqs. (22)-(24) | `p partial_alpha T+T F_alpha` 对应 `p gradT+T F`；`chi_kappa` 对应 `chi_kappa`；`R` 对应 `R^++R^-` | Eq. (24) 固定标量源：偶部为 `w_i Q`，奇部的一阶通量为 `(p gradT+T F)/rho_0+Q u+chi_kappa c_s^2 gradT`。 | 只使用 Eq. (22) 中保留到论文声明阶数的源设计；未指定的高阶离散矩不宣称为零。 |
| 第 9 页 / PDF 第 9 页 | Eqs. (26) 与 (28) | `Phi,Psi` 对应显式 LBE 源；`1/tau_fL,1/tau_gL` 对应 BGK 离散松弛率 | 在单松弛率 BGK 情形，梯形变换给出源前因子 `1-s/2`。 | 论文未给出双松弛率 TRT 分量式；本项目的 `(1-s_+/2)R^++(1-s_-/2)R^-` 必须从 Appendix A 的算子变换推导，不能说成 Eq. (26)/(28) 的逐字公式。 |
| 第 11 页 / PDF 第 11 页 | Eq. (36) 及本页离散化说明 | `bar(i)` 对应 `opposite[i]`；论文 2D 离散速度集对应 `d2q9()` | `bar(i)` 表示与 `i` 方向相反的离散速度；论文的二维流动场和标量场数值算例都采用 D2Q9。该页还说明标量平衡分布最高保留到二阶矩，而流动平衡分布涉及三阶矩。 | 该段关于求积阶数的说明是模型充分性说明，不是 D2Q5 公式。 |
| 第 33 页 / PDF 第 33 页 | Eqs. (A.1)-(A.4) | `tilde(f)=f-dt Omega/2` 对应本项目的一般算子梯形变量；`Phi` 对应显式源增量 | 沿特征线使用梯形规则并引入变换后，可消除隐式端点；BGK 特例显式碰撞和 `1-s/2` 源因子由此产生。 | TRT 需要把标量 BGK 算子推广为 `s_+P_++s_-P_-`，推广步骤属于本项目代数推导。 |
| 第 34 页 / PDF 第 34 页 | Eqs. (A.6b)-(A.7b) | `tilde(f)` 对应 `f_tilde`；`F_alpha` 对应 `force` | 变换分布的宏观动量重构为 `rho_0 u=sum_i c_i f_tilde_i+dt F/2`。 | 半源项不能在一次碰撞净源证明中遗漏；它使守恒动量的 transformed nonequilibrium moment 为 `-dt F/2`。 |
| 第 35 页 / PDF 第 35 页 | Eqs. (A.14)-(A.19) | `tilde(g)` 对应 `g_tilde`；`Psi` 对应显式标量源；`Q` 对应 `heat_source` | 标量场使用相同梯形变换；宏观量重构为 `T=sum_i g_tilde_i+dt Q/2`。 | Eq. (A.16) 是 BGK 特例；TRT 偶、奇源因子仍须按投影算子分别推导。 |

### Task 3：有效率与二阶恢复的新增核对

| 印刷页 / PDF 页 | 方程或位置 | 符号映射 | 允许主张 | 限制 |
| --- | --- | --- | --- | --- |
| 第 4 页 / PDF 第 4 页 | Eqs. (1a)-(1c)、Eq. (2) | `delta rho,u,p,T,F,Q,kappa` 对应本项目同名宏观符号；`sigma_alpha_beta` 对应恢复的黏性应力 | 目标连续性、弱可压 N-S 与 CDE 的分量形式；黏性应力分解为剪切偏差部分与 `mu^B tr(S) I`。 | 目标 PDE 只用于最终比对，不能作为 `second_order_residual_table()` 的零残差输入。 |
| 第 4-5 页 / PDF 第 4-5 页 | Eqs. (3)-(6) | `f, f^eq` 对应流场分布与平衡分布；连续矩对应 Task 1 离散矩约束 | `f^eq` 的零至三阶矩分别给出 `delta rho`、`rho0 u`、`delta rho c_s^2 I+rho0 uu` 与对称三阶速度矩；正文把截断 Maxwellian 的首误差标为 `O(Ma^3)`。 | 连续 Hermite 矩必须与 D2Q9 精确求积能力分开；不能据此消去所有离散高阶矩。 |
| 第 5-6 页 / PDF 第 5-6 页 | Eqs. (7)-(10) | `integral S dv`、`integral v S dv`、`integral vv S dv` 对应 flow source 的零、一、二阶矩 | 由零、一阶矩获得连续性和 Euler 层；Eq. (10) 展示平衡矩导数中的 `uF+Fu`，供二阶应力与源二阶矩逐项相消。 | Eq. (10) 使用连续 BGK 的领先 CE 近似；TRT 离散 Hénon 平移与半源项需另行推导。 |
| 第 6 页 / PDF 第 6 页 | Eqs. (11)-(14) | `chi_s,chi_b` 对应 shear/bulk feedback；`nu,nu^B` 对应物理输运系数 | Eq. (11) 给出源二阶矩及 Mach 分级；保留至 `O(Ma^2)` 得到 Eq. (13)，连续输运为 `nu=(1-chi_s)tau_f c_s^2`、`nu^B=2(1-chi_b)tau_f c_s^2/D`。 | Eq. (11) 的密度梯度和二次速度导数项是 `O(Ma^3)`，最后一项是 `O(Ma^4)`；不得称二阶恢复为无条件精确。 |
| 第 7 页 / PDF 第 7 页 | Eqs. (18)-(21) | 标量零阶源、平衡一/二阶矩导数对应 CDE 守恒层和非平衡通量核 | `integral R dv=Q`；Eq. (21) 在标量通量核中产生 `p grad(T)/rho0`、`T F/rho0`、`Q u` 与 `c_s^2 grad(T)`。 | Eq. (21) 仍是连续领先 CE；离散二阶系数应使用 `dt(1/s-1/2)`。 |
| 第 8 页 / PDF 第 8 页 | Eqs. (22)-(25) | `chi_kappa` 对应 scalar-flux feedback；`kappa` 对应物理扩散率 | Eq. (22) 要求源一阶矩包含 `p grad(T)/rho0+T F/rho0+Q u+(c_s^2-kappa/tau_g)grad(T)`；Eq. (24) 的参数化形式和 Eq. (25) 的连续 `kappa` 由此得到。 | Eq. (22) 的 `-T u div(u)` 明标为 `O(Ma^3)` 并被论文舍弃；本项目把它列为首个 Mach 余项而不是零。 |
| 第 9 页 / PDF 第 9 页 | Eqs. (26)-(29) | `tau_fL,tau_gL` 对应 BGK 格子松弛时间；TRT 中分别映射到物理 even/odd Hénon 平移 | 梯形离散后 `nu=(tau_fL-1/2)(1-chi_s)c_s^2 dt`、`nu^B=2(tau_fL-1/2)(1-chi_b)c_s^2dt/D`、`kappa=(tau_gL-1/2)(1-chi_kappa)c_s^2dt`。 | 原式是 BGK；TRT 中 shear/bulk 用 `sigma_f^+`、scalar flux 用 `sigma_g^-` 是项目内 parity 矩推导，不是论文逐字结论。 |
| 第 10 页 / PDF 第 10 页 | Eq. (30) | `integral vv(tilde f-f^eq)dv` 对应 transformed 二阶非平衡矩 | 给出 transformed 非平衡应力、应变、`uF+Fu` 和半时间步应变项的完整关系，可分别投影到 off-diagonal、deviatoric diagonal 与 trace。 | Eq. (30) 本身使用论文 BGK `tau_f`；有效 TRT 块率由代回 Task 2 的 parity collision 得到。 |
| 第 10 页 / PDF 第 10 页 | Eqs. (31)-(33) | off-diagonal closure、trace closure、diagonal reconstruction对应三个 flow physical blocks | Eq. (31) 闭合非对角应变；Eq. (32) 闭合 `tr(S)`；Eq. (33) 用迹重构对角应变。它们固定所有因子 2、`D`、`rho0`、`c_s^2` 与 `dt`。 | 局部应变闭合不改变未参与源反馈的 ghost 率，也不提供边界 magic parameter。 |
| 第 10 页 / PDF 第 10 页 | Eqs. (34)-(35) | `integral v(tilde g-g^eq)dv` 对应 scalar odd nonequilibrium flux；`p/rho0` 对应 `pressure_ratio` | Eq. (35) 给出 `grad(T)` 对局部非平衡通量以及 affine `T F/rho0+Q u` 的闭合；直接代回可导出冻结压力下的 scalar-flux 有效率。 | 该式支持局部闭合。空间变 `p/rho0` 的乘积导数、高阶 modified equation 和边界行为需另行分析。 |
| 第 33-35 页 / PDF 第 33-35 页 | Eqs. (A.1)-(A.5)、(A.6b)-(A.7b)、(A.14)-(A.19) | 一般梯形变换、flow/scalar 显式 LBE 与半源宏观重构 | 二阶 CE 可使用 `I-Lambda/2=Lambda Sigma`，并由半源重构消去 Taylor 时间离散项；正确 parity 源传递系数为 `(1-s/2)/(s sigma)=1`。 | Appendix 只印出 BGK 分量式；`Lambda=s_+P_++s_-P_-` 的算子推广和残差扰动属于本项目显式代数。 |

## 2. Wang et al. (2013)

原始 PDF：`D:\桌面\代码\代码\热对流\2DRB\pdf\Lattice Boltzmann simulations of thermal convective flows in two dimensions.pdf`

| 印刷页 / PDF 页 | 方程或位置 | 符号映射 | 允许主张 | 限制 |
| --- | --- | --- | --- | --- |
| 期刊页 263 / PDF 第 2 页 | Eq. (1) | 论文 `e_i` 对应 `d2q9().velocities[i]` | D2Q9 可按 `(0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)` 排序。 | 这里只固定速度及顺序，不固定温度模型的自由参数。 |
| 期刊页 264 / PDF 第 3 页 | Eq. (7) | 论文 `c_s` 对应 `sqrt(Lattice.cs2)` | 在格子速度取 1 的标准 D2Q9 标度下，`c_s^2 = 1/3`。 | 不可由此声称所有 D2Q5 参数化模型都唯一采用同一权重。 |
| 期刊页 264 / PDF 第 3 页 | Eqs. (11)-(12) | 论文 D2Q5 的 `e_i, i=0,...,4` 对应 `d2q5().velocities`；`N` 是温度分布的矩变换 | 该 D2Q5 使用 D2Q9 的前五个速度；矩基可用于从平衡矩反解静止态权重。 | **D2Q5-only。** |
| 期刊页 265 / PDF 第 4 页 | Eqs. (14)-(17) | 论文自由参数 `a` 对应 D2Q5 平衡矩第四分量的系数；`T` 对应零阶标量矩；`sigma_kappa,sigma_nu` 是两组温度非守恒矩的松弛率 | D2Q5 平衡矩为 `(T,uT,vT,aT,0)`，`a` 是自由参数；Eq. (17) 只给出四阶误差项的**各向同性条件**。 | **D2Q5-only。** Eq. (17) 本身不是完整四阶误差消除条件；Eq. (14) 也不会唯一固定 `w_0,w_1,...,w_4`，除非另行选定 `a`。 |
| 期刊页 265 / PDF 第 4 页 | Eqs. (18)-(19) | `sigma_e,sigma_nu,sigma_kappa` 是 D2Q5 温度模型的非守恒矩松弛率 | 若要进一步主张完整四阶误差消除，还必须结合额外松弛率关系和这些特殊参数值，不能只引用 Eq. (17)。 | **D2Q5-only。** 这些关系不得导入 D2Q9。 |

## 3. Contrino et al. / Luo (2014)

原始 PDF：`D:\桌面\代码\代码\热对流\2DRB\pdf\[Luo2014]_JCP_LB simulations of the thermally driven 2D square cavity at high Rayleigh numbers.pdf`

| 印刷页 / PDF 页 | 方程或位置 | 符号映射 | 允许主张 | 限制 |
| --- | --- | --- | --- | --- |
| 期刊页 258 / PDF 第 2 页 | Eq. (1) | 论文 `e_i` 对应 D2Q9 离散速度 | 独立复核标准 D2Q9 速度集合。 | 紧凑的正负号写法不单独指定本项目的索引顺序。 |
| 期刊页 259 / PDF 第 3 页 | Eq. (9) | 论文 `c_s` 对应 `sqrt(Lattice.cs2)` | 独立复核标准格子标度下 `c_s^2 = 1/3`。 | 该式位于流动 D2Q9 模型部分。 |
| 期刊页 260 / PDF 第 4 页 | Eqs. (13)-(17) | 论文 `n^eq=(T,uT,vT,aT,0)^T` 对应 D2Q5 平衡矩；`a` 是自由参数 | 复核 D2Q5 的参数化平衡矩、松弛矩阵和四阶误差条件。 | **D2Q5-only。** `a` 的存在再次说明 D2Q5 权重不是由速度集唯一决定。 |
| 期刊页 260 / PDF 第 4 页 | Eqs. (18)-(21) | 论文 `sigma_k=1/s_k-1/2` 对应后续 TRT 的 Hénon 平移参数 | 这些特定参数关系可消除该 D2Q5 热模型的完整四阶误差。 | **D2Q5-only。** 不得作为 D2Q9 TRT 系数的证据。 |

## 4. Dubois & Lallemand, Towards higher order lattice Boltzmann schemes

原始 PDF：`D:\桌面\代码\代码\热对流\2DRB\pdf\Towards higher order lattice Boltzmann schemes.pdf`

| 印刷页 / PDF 页 | 方程或位置 | 符号映射 | 允许主张 | 限制 |
| --- | --- | --- | --- | --- |
| 第 5 页 / PDF 第 5 页 | Eq. (16) | 论文 `sigma_k` 对应后续代码中的 Hénon 参数；`s_k` 对应矩空间松弛率 | 可使用 `sigma_k = 1/s_k - 1/2` 在松弛率与平移参数之间转换。 | 这只是参数定义，不自动给出 TRT magic parameter。 |
| 第 11 页 / PDF 第 11 页 | Eqs. (40)-(42) | 论文 D2Q5 的 `kappa_40,kappa_22` 对应四阶等效方程系数 | 可引用 D2Q5 四阶等效方程及其两个误差系数。 | **D2Q5-only。** |
| 第 16 页 / PDF 第 16 页 | Eq. (55) | 论文 `sigma_1,sigma_3` 对应 D2Q5 两组非守恒矩的 Hénon 参数 | 可引用该 D2Q5 模型完整四阶消除的一组参数值。 | **D2Q5-only。** 不得迁移到 D2Q9。 |

## 5. Ginzburg & d'Humières (2003)

原始 PDF：`D:\桌面\代码\代码\热对流\2DRB\pdf\Multireflection boundary conditions for lattice Boltzmann models.pdf`

| 印刷页 / PDF 页 | 方程或位置 | 符号映射 | 允许主张 | 限制 |
| --- | --- | --- | --- | --- |
| 文章页 066614-2 / PDF 第 2 页 | Table I 及 Eq. (1) 前的中心对称定义 | `c_q` 对应 `velocities[q]`；`c_bar(q)=-c_q` 对应 `opposite[q]`；表中 `t_0,t_1,t_2` 对应 D2Q9 静止、轴向、对角权重 | 在 `c_s^2=1/3` 时，标准 D2Q9 权重是 `4/9,1/9,1/36`；非零速度按中心对称成对，因而定义相反方向映射。 | Table I 使用 `t_1^*=3t_1`、`t_2^*=3t_2`，换算时必须保留这个因子。 |
| 文章页 066614-6 / PDF 第 6 页 | Eq. (41) | 论文 `lambda_nu,lambda_2` 对应偶、奇（或黏性/伴随）松弛率；括号项对应 Hénon 平移参数 | 可把论文的边界参数写成 `Lambda^2=(4/3)(1/lambda_nu-1/2)(1/lambda_2-1/2)`。 | 这是该论文的归一化；与别处的 `Lambda` 命名比较时必须先做符号映射。 |
| 文章页 066614-7 / PDF 第 7 页 | Eqs. (42)-(43) 及其后正文 | `D_eff` 是有效通道宽度；`D_1/2` 是半格点壁面位置对应宽度 | 当 `Lambda^2=1/4` 时，Eq. (42) 的宽度偏移项消失；Eq. (43) 给出相应的两松弛率关系。 | 该结论属于反弹边界位置分析，不是体相 D2Q9 权重或 LBM-CDE `lambda_t` 的来源。 |

## 本任务采用的精确离散常数

### D2Q9

- 速度顺序：`(0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)`。
- 权重：`(4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36)`。
- `c_s^2=1/3`。
- 相反方向：`(0,3,4,1,2,7,8,5,6)`。
- LBM-CDE 修正：`lambda_t=(-5/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36)`。

最后一项不是论文逐项抄录值，而是可复查的离散推导：按 `LBM-CDE.pdf` 第 7 页的构造保留全部移动方向的标准 D2Q9 权重，只修改零速度项；再施加 Eq. (16) 后的零阶约束，移动方向之和为 `5/9`，故零速度项必须为 `-5/9`。由 D2Q9 权重矩关系可立即验证其一阶矩为零、二阶矩为 `c_s^2 delta_alpha_beta`、三阶矩为零。

该低阶构造不会消去离散高阶矩。对上述 `lambda_t` 精确求和得到 `M40=M04=1/3`、`M22=1/9`、`M31=M13=0`。Task 5 的等效方程推导必须保留这些非零四阶贡献，不能把论文的“二阶矩修正”简化成“只有二阶矩非零”。

Task 2 的源项三、四阶 raw moment 表是把 Eqs. (13)、(24) 代入本节 D2Q9 权重后逐项精确求和得到的项目内代数结果，不是论文直接列出的高阶闭合条件。账本只允许在已显式计算到的四阶范围内使用这些值；更高阶矩不得据此宣称为零。

### D2Q5 验证格点

本任务的通用代数验证器固定采用
`w=(1/3,1/6,1/6,1/6,1/6)`、`c_s^2=1/3`，速度是 D2Q9 的前五项，相反方向为 `(0,3,4,1,2)`。这是一个明确选定的验证点，不是 Wang 或 Contrino/Luo 论文唯一规定的权重。

由 Wang 的 Eqs. (12) 与 (14) 在静止态反解可得
`w_0=(1-a)/5`、四个移动方向 `w_i=(4+a)/20`。选择 `a=-2/3` 才得到本任务的 `1/3,1/6` 权重以及 `c_s^2=1/3`。因此，后续若研究论文中的可变 `a` 模型，必须显式构造参数化 D2Q5，不能把本验证格点当作一般结论。
