# Task 02 实施报告：奇偶源项与算子梯形 TRT 碰撞

## 结果

- 分支：`codex/lbm-cde-trt-derivation`
- 起始提交：`d0d2b29de6e4d847a77a586413aa6b82de367e60`
- 提交信息：`feat: derive parity-resolved TRT sources`
- 最终提交哈希：随父任务交接提供；Git 提交无法在自身内容中嵌入不会再次改变的自身对象 ID。
- 状态：Task 02 数学合同已由精确 SymPy 测试闭合；没有未决的符号或归一化约定。

## 交付文件

- `tools/derivation/sources.py`
- `tools/derivation/collision.py`
- `tests/derivation/test_sources_collision.py`
- `docs/derivation/chapters/02-trt-source-design.md`
- `docs/derivation/evidence-ledger.md`
- `.superpowers/sdd/task-02-report.md`

未修改生产 Fortran、原始 PDF、Task 1 格点常数或 `Xs/`。

## Eq. (13) 直接核对结论

已目视核对原始 `LBM-CDE.pdf` 第 4 页 Eq. (2) 和第 6 页 Eq. (13)：

- `S_alpha_beta=(partial_alpha u_beta+partial_beta u_alpha)/2`；
- 采用无量纲 `H_i^(2)=c_i c_i/c_s^2-I`；
- `A=chi_s S+(chi_b-chi_s)tr(S)I/2`，`A` 内不再额外除以 `c_s^2`；
- 收缩项是正号 `+rho_0 w_i H_i^(2):A`；
- 对称 `A` 的离散二阶矩贡献为 `+2 rho_0 c_s^2 A`。

测试 `test_eq13_positive_dimensionless_hermite_contraction` 会同时捕获错误负号和多/少一个 `c_s^2` 的归一化错误。

## 实现摘要

- `flow_source()` 与 `scalar_source()` 按 D2Q9 相反方向算子精确生成 raw、even、odd 三个向量及零至四阶 raw moment 表。
- 源矩测试覆盖 brief 中全部零至二阶约束，并保留后续 Taylor 分析所需的全部非零三、四阶矩；没有对未计算的更高阶矩作消失声明。
- `trt_collision()` 从一般算子梯形变换得到两个独立源因子 `1-s_plus/2` 与 `1-s_minus/2`。
- `reconstruct_momentum()` 和 `reconstruct_scalar()` 实现半源宏观重构。
- 一次碰撞测试分别拆出 half-source nonequilibrium 的松弛贡献和显式源贡献，证明对任意符号松弛率总增量严格为 `dt*F` 与 `dt*Q`。
- `s_plus=s_minus=s` 后，用任意九分量符号向量逐分量验证 TRT 与 transformed BGK 完全相同。

## RED

先只创建 delayed-import 测试，未创建 `sources.py` 或 `collision.py`。命令：

```text
python -m unittest tests.derivation.test_sources_collision -v
```

全部失败测试：

```text
test_equal_rates_recover_bgk_componentwise ... FAIL
test_half_source_macroscopic_reconstructions ... FAIL
test_one_flow_collision_has_exact_net_momentum_source ... FAIL
test_one_scalar_collision_has_exact_net_heat_source ... FAIL
test_trt_collision_matches_operator_trapezoidal_formula ... FAIL
test_eq13_positive_dimensionless_hermite_contraction ... FAIL
test_flow_source_required_low_order_moments ... FAIL
test_nonzero_third_and_fourth_raw_source_moments_are_retained ... FAIL
test_scalar_source_required_low_order_moments ... FAIL
test_sources_are_exact_even_odd_projections ... FAIL

AssertionError: "ModuleNotFoundError: No module named
'tools.derivation.collision'" is not None : source/collision API is
unavailable: ModuleNotFoundError: No module named 'tools.derivation.collision'

Ran 10 tests in 0.002s
FAILED (failures=10)
```

这是预期 assertion `FAIL`；没有 import `ERROR`、语法错误或测试收集错误。

## GREEN

最小 API 完成后的 targeted 命令：

```text
python -m unittest tests.derivation.test_sources_collision -v
```

最终 targeted 输出：

```text
test_equal_rates_recover_bgk_componentwise ... ok
test_half_source_macroscopic_reconstructions ... ok
test_one_flow_collision_has_exact_net_momentum_source ... ok
test_one_scalar_collision_has_exact_net_heat_source ... ok
test_trt_collision_matches_operator_trapezoidal_formula ... ok
test_eq13_positive_dimensionless_hermite_contraction ... ok
test_flow_source_required_low_order_moments ... ok
test_nonzero_third_and_fourth_raw_source_moments_are_retained ... ok
test_scalar_source_required_low_order_moments ... ok
test_sources_are_exact_even_odd_projections ... ok

Ran 10 tests in 9.619s
OK
```

首次实现运行时有 5 项 `assertEqual` 因 SymPy 展开式与因式分解式结构不同而失败。逐项计算 `simplify(actual-expected)` 均为零；只把符号表达式断言改为差值化简，没有修改任何物理系数。随后得到上述 targeted GREEN。

完整回归命令：

```text
python -m unittest discover -s tests/derivation -v
```

结果：Task 1 与 Task 2 共 21 项全部 `ok`，`Ran 21 tests in 9.740s`，`OK`。

## 原始 PDF 目视核对记录

原始文件：`D:\桌面\代码\代码\热对流\2DRB\pdf\LBM-CDE.pdf`。

- 第 4 页：Eq. (2)，应变率、迹和黏性应力定义；Eq. (5)，连续权函数矩。
- 第 5-6 页：Eqs. (7a)-(13)，流场源的零、一、二阶目标矩和最终 Hermite 形式。
- 第 7-8 页：Eqs. (15)-(25)，标量平衡矩、标量源零/一阶目标矩及 Eq. (24) 最终源式。
- 第 9 页：Eqs. (26)-(29)，BGK 离散源因子 `1-s/2`。
- 第 33 页：Eqs. (A.1)-(A.5)，沿特征线梯形积分、transformed distribution 与显式 BGK 碰撞。
- 第 34 页：Eqs. (A.6b)-(A.7b)，动量半力重构。
- 第 35 页：Eqs. (A.14)-(A.19)，标量梯形变换、显式源与 `dt*Q/2` 重构。

论文只直接给出 BGK 离散式；TRT 的奇偶源因子来自 Appendix A 变换的一般算子推广，已在中文推导章中逐步推导，不把该推广冒充成论文原式。

## 物理与数值边界

- 物理不变量保持为一次碰撞净动量源 `dt*F`、净标量源 `dt*Q`，且与各自的 TRT 松弛率无关。
- 流场与标量场松弛率对彼此独立；没有把 D2Q5-only 的高阶消除关系导入 D2Q9。
- 非零三、四阶源矩被明确保留给 Task 5；未指定高阶矩不宣称为零。
- 本任务只实现推导工具，不改变生产碰撞、边界条件、Nu/Re 后处理或基准输出，因此不需要运行 Fortran benchmark。
