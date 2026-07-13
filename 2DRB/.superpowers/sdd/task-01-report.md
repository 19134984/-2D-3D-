# Task 01 实施报告：精确 D2Q 格点代数与证据账本

## 结果

- 分支：`codex/lbm-cde-trt-derivation`
- 起始提交：`bc1299a`
- 提交信息：`test: establish exact D2Q lattice algebra`
- 最终提交哈希：随父任务交接报告提供。Git 对象 ID 由包括本报告在内的提交内容计算，无法把自身最终哈希原样嵌入同一个提交而不再次改变该哈希。
- 状态：Task 01 范围内实现与验证通过；没有会改变格点系数的未决歧义。

## 交付文件

- `tools/derivation/__init__.py`
- `tools/derivation/lattice.py`
- `tests/derivation/__init__.py`
- `tests/derivation/test_lattice.py`
- `docs/derivation/evidence-ledger.md`
- `.superpowers/sdd/task-01-report.md`

未修改任何生产 Fortran 文件、原始 PDF 或 `tmp/` 文本提取。

## 实现摘要

- 使用 SymPy `Integer`/`Rational` 定义 D2Q9 和 D2Q5 的速度、精确权重、`c_s^2`、相反方向，以及 LBM-CDE D2Q9 `lambda_t`。
- 实现任意二维幂次的精确 raw moment、0 至 3 阶 Cartesian Hermite moment 和由相反方向置换构造的偶/奇 parity projector。
- 测试覆盖相反方向对合、无浮点常数、D2Q9 二/四阶权重矩、D2Q5 四阶局限、`lambda_t` 的零至三阶约束、Hermite 矩和投影算子代数。
- D2Q5 验证格点明确固定为 Wang/Contrino-Luo 参数化模型的 `a=-2/3` 点；没有把这一点误写成文献唯一规定。
- D2Q9 `lambda_t` 明确标记为依据 `LBM-CDE.pdf` Eq. (16) 后约束和“只修改零速度项”的文字构造推导，而非论文逐项印出的向量。

## RED

先只创建测试包和测试文件，未创建 `tools/derivation` 实现。命令：

```text
python -m unittest tests.derivation.test_lattice -v
```

输出：

```text
test_d2q5_exact_definition (tests.derivation.test_lattice.LatticeDefinitionTests.test_d2q5_exact_definition) ... FAIL
test_d2q9_exact_definition (tests.derivation.test_lattice.LatticeDefinitionTests.test_d2q9_exact_definition) ... FAIL
test_lattice_constants_contain_no_floating_values (tests.derivation.test_lattice.LatticeDefinitionTests.test_lattice_constants_contain_no_floating_values) ... FAIL
test_opposite_maps_are_involutions (tests.derivation.test_lattice.LatticeDefinitionTests.test_opposite_maps_are_involutions) ... FAIL
test_d2q5_weight_moments_and_fourth_order_limitation (tests.derivation.test_lattice.LatticeMomentTests.test_d2q5_weight_moments_and_fourth_order_limitation) ... FAIL
test_d2q9_lambda_t_constraints (tests.derivation.test_lattice.LatticeMomentTests.test_d2q9_lambda_t_constraints) ... FAIL
test_d2q9_weight_moments (tests.derivation.test_lattice.LatticeMomentTests.test_d2q9_weight_moments) ... FAIL
test_raw_and_hermite_moments_are_exact (tests.derivation.test_lattice.LatticeMomentTests.test_raw_and_hermite_moments_are_exact) ... FAIL
test_projector_algebra (tests.derivation.test_lattice.ParityProjectorTests.test_projector_algebra) ... FAIL
test_projectors_match_pairwise_even_odd_definition (tests.derivation.test_lattice.ParityProjectorTests.test_projectors_match_pairwise_even_odd_definition) ... FAIL

======================================================================
FAIL: test_d2q5_exact_definition (tests.derivation.test_lattice.LatticeDefinitionTests.test_d2q5_exact_definition)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "...\\tests\\derivation\\test_lattice.py", line 29, in setUp
    self.assertIsNone(
AssertionError: "ModuleNotFoundError: No module named 'tools.derivation'" is not None : exact lattice API is unavailable: ModuleNotFoundError: No module named 'tools.derivation'

[其余 9 个 FAIL 具有相同的 setUp assertion traceback；上方已保留全部失败测试名。]

----------------------------------------------------------------------
Ran 10 tests in 0.002s

FAILED (failures=10)
```

这是预期的 assertion failure；测试收集和 Python 语法均正常，也没有以 `ERROR` 结束。

## GREEN

实现最小精确代数 API 后运行同一命令：

```text
python -m unittest tests.derivation.test_lattice -v
```

输出：

```text
test_d2q5_exact_definition (tests.derivation.test_lattice.LatticeDefinitionTests.test_d2q5_exact_definition) ... ok
test_d2q9_exact_definition (tests.derivation.test_lattice.LatticeDefinitionTests.test_d2q9_exact_definition) ... ok
test_lattice_constants_contain_no_floating_values (tests.derivation.test_lattice.LatticeDefinitionTests.test_lattice_constants_contain_no_floating_values) ... ok
test_opposite_maps_are_involutions (tests.derivation.test_lattice.LatticeDefinitionTests.test_opposite_maps_are_involutions) ... ok
test_d2q5_weight_moments_and_fourth_order_limitation (tests.derivation.test_lattice.LatticeMomentTests.test_d2q5_weight_moments_and_fourth_order_limitation) ... ok
test_d2q9_lambda_t_constraints (tests.derivation.test_lattice.LatticeMomentTests.test_d2q9_lambda_t_constraints) ... ok
test_d2q9_weight_moments (tests.derivation.test_lattice.LatticeMomentTests.test_d2q9_weight_moments) ... ok
test_raw_and_hermite_moments_are_exact (tests.derivation.test_lattice.LatticeMomentTests.test_raw_and_hermite_moments_are_exact) ... ok
test_projector_algebra (tests.derivation.test_lattice.ParityProjectorTests.test_projector_algebra) ... ok
test_projectors_match_pairwise_even_odd_definition (tests.derivation.test_lattice.ParityProjectorTests.test_projectors_match_pairwise_even_odd_definition) ... ok

----------------------------------------------------------------------
Ran 10 tests in 0.004s

OK
```

## 原始 PDF 目视核对记录

以下页面均从原始 PDF 直接渲染后逐页目视检查；文本提取只用于搜索。

- `LBM-CDE.pdf`：印刷/PDF 第 4 页 Eq. (5)；第 7 页 Eqs. (15)-(21)；第 11 页 Eqs. (36)-(40) 及离散化说明。
- `Lattice Boltzmann simulations of thermal convective flows in two dimensions.pdf`：期刊页 263-265（PDF 第 2-4 页），Eqs. (1), (7), (11)-(17)。
- `[Luo2014]_JCP_LB simulations of the thermally driven 2D square cavity at high Rayleigh numbers.pdf`：期刊页 258-260（PDF 第 2-4 页），Eqs. (1), (9), (13)-(21)。
- `Towards higher order lattice Boltzmann schemes.pdf`：第 5 页 Eq. (16)，第 11 页 Eqs. (39)-(43)，第 16 页 Eq. (55)。
- `Multireflection boundary conditions for lattice Boltzmann models.pdf`：文章页 066614-2 的 Table I 与中心对称定义、066614-6 Eq. (41)、066614-7 Eqs. (42)-(43)。

任务计划总约束中的一处文字写“four PDFs”，但更具体的 Task 01 brief 明确列出五篇本地论文并要求账本包含五篇；本任务按 brief 核对了全部五篇。

## 歧义与边界

- **无系数变更歧义：** D2Q9 权重、`c_s^2`、相反方向和 `lambda_t` 均可由已核对证据与精确矩测试闭合验证。
- **D2Q5 参数边界：** 文献的 D2Q5 平衡矩含自由参数 `a`；本任务只把 `a=-2/3` 的标准二阶各向同性格点作为代数验证器。后续参数化 D2Q5 工作必须另行建模。
- **符号边界：** Ginzburg & d'Humières 使用 `Lambda^2` 和自身松弛率归一化；账本只记录原式与映射，尚未把它当作后续 D2Q9 TRT magic parameter 的无条件同义词。
- PDF 渲染器报告个别字体边界框/缺失字体警告，但上述公式、表格和页码在渲染图中清晰可辨，不影响证据读取。
