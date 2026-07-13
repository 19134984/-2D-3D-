# 08 LBM-CDE-TRT 单步算法规范

本章把前七章的代数结果整理为一个可直接实现的时间步。流程顺序是规范的一部分：宏观重构、反馈量、源项、碰撞和边界不能为方便而交换。生产实现仍应从仓库相同维度、相同后端的最近基线增量修改。

## 1. 输入与启动检查

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

## 2. 一个时间步的固定顺序

### 第 1 步：用半力、半热源重构 rho、u、T

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

### 第 2 步：构造含压力修正的平衡态

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

### 第 3 步：从局部非平衡重构应变与温度梯度

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

### 第 4 步：分别构造 raw even/odd 源

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

### 第 5 步：用名义奇偶率和各自源因子碰撞

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

### 第 6 步：stream

执行 pull 或 push streaming，但全程序只能选择一种索引约定：

$$
\widetilde h_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
=\widetilde h_i^*(\boldsymbol x,t).
$$

MPI/OpenACC 版本必须在同一逻辑位置完成 halo/device 同步；同步只改变数据可见性，不改变碰撞—stream—边界的数学顺序。

### 第 7 步：依次施加速度 BB、温度 ABB、绝热 BB 和显式 corner rule

1. 对固壁速度 link 应用所选 halfway BB。只有平直、格点对齐的 halfway 壁、稳态一维 Stokes 剪切、均匀体力、half-force gauge、冻结 feedback 且无 bulk/time/tangential jets 的受限子域，才可使用 $(1-\chi_s)\sigma_f^+\sigma_f^-=3/16$ 标定；pressure-driven gauge 的 $3/8$ 不与它混用。
2. 对 Dirichlet 温度 link 应用 ABB，墙项保留

   $$
   2\left(w_i+\lambda_i\frac{p_w/\rho_0}{c_s^2}\right)T_w.
   $$

3. 对绝热 link 应用 homogeneous Neumann/BB；不要把 uniform heat source 重复写入 odd-flux 边界行。已解析行限定于平直、格点对齐的 halfway 壁、half-source、transformed CDE chain 和 D2Q9 $c_s^2=1/3$；其中 affine-normal 还要求稳态、常压力比和零切向 jets。它们不代表一般闭合：一般 wall-time、normal-curvature 和 source rows 仍为 `unresolved`，需要另行边界修正。
4. 在格点对齐直角、half-source 规范的 mixed Dirichlet/adiabatic corner，对共享 diagonal population 只赋值一次。默认显式策略为：轴向 links 各自执行所属边界，共享 diagonal 采用 Dirichlet-priority ABB，跳过随后对同一 population 的 adiabatic overwrite，并把未满足的绝热 corner residual 记入诊断。若项目选择别的单一闭合，必须替换整条 corner 策略，而不是交换两次 overwrite 的顺序。

共享 corner source increment 只计一次，diagonal wall distance 使用 $h/\sqrt2$。一般 mixed corner 的两条方程秩不兼容，因此这里明确选择工程闭合，不宣称同时精确满足两类边界。

### 第 8 步：计算诊断与不变量

每步或按固定间隔检查：

- 质量、动量和标量总量，以及施加源后的预算残差；
- $\rho,T,p,\boldsymbol u$ 和所有 populations 是否有限；
- $a$、$b$ 与反馈闭合分母是否跨过奇点；
- 名义率、物理有效率、输运系数是否仍与启动参数映射一致；
- 边界 residual、corner overwrite 计数和共享 source 计数；
- `errorU`、`errorT`、Nu、Re、壁面极值及算例规定的收敛判据。

冻结系数四阶条件只在其审计子域内作为诊断；若 $p/\rho_0$、松弛率或 $\chi$ 空间变化，必须另行记录乘积导数残差，不能继续标记为四阶完全消除。

## 3. 三种标量运行策略

配置层必须显式选择以下一种，不允许求解器暗中替换：

### A. 低 kappa + 体相四阶消除 + 显式边界修正

名义率满足 Task 5 四阶条件；求解器把 Dirichlet 墙面的显式修正报告为结构性最小扩展。本文尚未给出 correction formula 或 corrected residual，因此该分支状态仍是 `boundary_correction_required`，不是已实现方案，更不是一般 source/time/pressure jets 的通用修复。

### B. 低 kappa + 受限 ABB 标定 + 保留体相 q4 误差

仅在平直、格点对齐的 halfway 壁、稳态一维二次温度场、half-source 重构、完整 pressure-wall equilibrium term，以及相匹配的 external-gradient 或 local-feedback population chain 下，用 $\sigma_{\rm flux}\sigma_e=3/16$ 生成名义偶率。随后计算并报告非零 $C_{40},C_{22}$；不得把“各向同性残差为零”写成“四阶误差为零”。

### C. 两项都强制时采用独立新推导

显式受限 ABB 边界修正是已识别的结构性最小扩展，但尚待 correction formula、实现与 corrected-residual 验证。split-even MRT 是增加独立偶模态的另一候选；在明确 ABB 与 $C_{40}$ 分别依赖哪个偶模态、并证明两约束 Jacobian 满秩前，只能标记 `split_even_mrt_candidate_requiring_mode_jacobian_derivation`。两条路线当前都不是已证明充分方案，也不能预先标记 `feasible_exact`。

## 4. 最小伪代码

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
