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
