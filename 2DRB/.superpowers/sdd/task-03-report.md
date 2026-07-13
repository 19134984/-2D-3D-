# Task 03 实施报告：有效率与二阶宏观恢复

## 结果

- 分支：`codex/lbm-cde-trt-derivation`
- 起始提交：`c56f022ff670f74c704e2f670c300ca72110ab6d`
- 提交信息：`feat: derive effective rates and second-order recovery`
- 最终提交哈希：随父任务交接提供；提交不能在自身内容中稳定嵌入自己的对象 ID。
- 范围：只新增符号推导工具、测试、中文推导章、证据账本和本报告；未修改生产 Fortran、原始 PDF 或 `Xs/`。

## 交付文件

- `tools/derivation/effective_rates.py`
- `tests/derivation/test_effective_rates.py`
- `docs/derivation/chapters/03-effective-rates.md`
- `docs/derivation/chapters/04-second-order-recovery.md`
- `docs/derivation/evidence-ledger.md`
- `.superpowers/sdd/task-03-report.md`

## 原始 PDF 式号核对

原始证据文件为 `D:\桌面\代码\代码\热对流\2DRB\pdf\LBM-CDE.pdf`，没有使用 `Xs/` 作为理论证据。

- PDF 第 4-6 页：Eqs. (1)-(14)，目标宏观方程、流场平衡矩、逆向源矩、Mach 分级和连续 shear/bulk 输运系数。
- PDF 第 7-8 页：Eqs. (18)-(25)，标量守恒矩、平衡通量核、`p grad(T)`、`T F`、`Q u` 及 `chi_kappa` 源设计。
- PDF 第 9 页：Eqs. (26)-(29)，BGK 梯形离散源因子与离散输运系数。
- PDF 第 10 页：Eqs. (30)-(35)，transformed flow 二阶非平衡矩、off-diagonal/trace/diagonal 应变闭合和 scalar gradient 闭合。
- PDF 第 33-35 页：Eqs. (A.1)-(A.7b)、(A.14)-(A.19)，梯形变量、显式 LBE 和两种半源宏观重构。

账本把论文直接结论与项目内 TRT parity 推广分开记录；论文只印出 BGK 分量式。

## 局部反馈消元结论

定义名义 Hénon 平移 `sigma=1/s-1/2`。把 Eqs. (31)-(35) 直接代回 Task 2 的实际 collision moment，得到

```text
sigma_shear_eff = (1-chi_s) sigma_f_plus
sigma_bulk_eff  = (1-chi_b) sigma_f_plus
sigma_flux_eff  = (1-chi_kappa) cs2/(cs2+p/rho0) sigma_g_minus
s_eff           = 1/(sigma_eff+1/2)
```

物理输运系数为

```text
nu    = (1-chi_s) cs2 dt sigma_f_plus
nu_B  = (2/D)(1-chi_b) cs2 dt sigma_f_plus
kappa = (1-chi_kappa) cs2 dt sigma_g_minus
```

`T F/rho0`、`Q u` 和 flow 的 `uF+Fu` 均保留为仿射项，没有被吸收入率。有效物理块与同 parity ghost 分裂：flow odd ghost、scalar odd ghost 以及未反馈 even ghost 仍使用名义移位；所以算法是 TRT 骨架，局部解释是受限块 MRT，不是新增名义输入。

标量压力消去只按 frozen `pi=p/rho0` 局部模态主张。若在 `pi0` 冻结率但实际有 `delta pi(x)`，首个乘积导数残差是

```text
dt sigma_flux_eff(pi0)
  [grad(delta pi).grad(T) + delta pi laplacian(T)]
```

该式一般非零。

## 二阶 CE 与结论等级

离散推导采用

```text
grad = epsilon grad_1
d_t  = epsilon d_t1 + epsilon^2 d_t2
h_tilde = h_eq + epsilon h_1 + epsilon^2 h_2 + O(epsilon^3)
```

并独立采用与论文 Eqs. (11)、(22) 一致的低 Mach 幅值计数：`u=O(Ma)`、`delta_rho/rho0=O(Ma^2)`、`p=p0+O(Ma^2)`、`F^(1)=O(Ma)`、`Q^(1)=O(Ma^0)`。

章节逐式写出 transformed LBE 的 `epsilon^0/1/2` 人口方程。恒等式

```text
I-Lambda/2 = Lambda Sigma
```

与 `dt*F/2`、`dt*Q/2` 宏观重构一起消去二阶 Taylor 时间离散项。流场 `uF`、标量 `p grad(T)`、`T F`、`Q u` 均由命名平衡矩和源矩逐项相消，不把目标 PDE 输入残差表。

首个遗漏 CE 家族为 `O(epsilon^3)`；首个 Mach 家族为 `O(Ma^3)`。报告章明确保留论文 Eq. (11) 的密度梯度/二次速度导数项、Eq. (22) 的 `-T u div(u)`、变系数乘积导数、`d_t F/d_t Q` 的更高时间导数、ghost/边界和 Task 2 三四阶矩贡献。

## 残差生成与扰动证据

`second_order_residual_table()` 从 `EquilibriumMomentConstraints`、`SourceMomentConstraints` 和四个 parity 梯形因子生成：

```text
p_grad_T, T_F, Q_u, u_F, d_t_F, d_t_Q,
s_f_minus_transport, s_g_plus_transport,
first_omitted_epsilon_order, first_omitted_mach_order
```

非守恒构成矩的源传递系数为 `r=b/(s sigma)`；守恒 force/heat 因 half-source nonequilibrium 的松弛贡献，Euler 净系数为 `c=b+s/2`。正确 `b=1-s/2` 时两者都为 1。

默认六个命名残差均由 SymPy 化简为零。扰动的精确结果是：

- 删除 `p grad(T)`、`T F`、`Q u` 或 `uF` 源矩时，对应残差为 1。
- scalar odd 因子错设为 1 时，`p_grad_T=T_F=Q_u=s_g_minus/(s_g_minus-2)`。
- flow even 因子错设为 1 时，`u_F=s_f_plus/(s_f_plus-2)`。
- flow odd 因子错设为 1 时，`T_F=u_F=s_f_minus/2`，`d_t_F=1/2`。
- scalar even 因子错设为 1 时，`Q_u=s_g_plus/2`，`d_t_Q=1/2`。
- 直接求导验证 `partial(nu)/partial(s_f_minus)=0`、`partial(kappa)/partial(s_g_plus)=0`，而对正确物理 parity 率的导数非零。

## TDD 记录

### 初始 RED

先只创建 delayed-import 测试，没有创建 `effective_rates.py`：

```text
python -m unittest tests.derivation.test_effective_rates -v

Ran 11 tests in 0.002s
FAILED (failures=11)
AssertionError: "ModuleNotFoundError: No module named
'tools.derivation.effective_rates'" is not None
```

全部是捕获导入异常后的 assertion `FAIL`，没有测试收集、解析或 import `ERROR`。

### 初次 GREEN 与真实符号失败

最小实现首轮为 10/11；唯一失败来自变压力残差辅助路径把 `c_s^2=1/3` 写死，而测试使用任意符号 `c_s^2`。将 `c_s^2,dt` 随有效模态结果传递后，11/11 通过。随后增加 ghost block 与 equal-rate BGK 测试，targeted 达到 13/13。

### CE 守恒源系数返工 RED/GREEN

数学自审发现默认正确因子下 `r=c=1` 会掩盖守恒/非守恒源传递的差别。新增精确错因子断言后先得到

```text
test_wrong_conserved_source_factor_has_exact_half_source_residual ... FAIL
AssertionError: -s_f_minus**2/(2*s_f_minus - 4) != 0

Ran 1 test in 0.189s
FAILED (failures=1)
```

将 Euler force/heat 改为 `c=b+s/2` 后，该目标测试 1/1 通过，最终 Task 3 targeted 14/14 通过。

## 验证结果

```text
python -m unittest tests.derivation.test_effective_rates -v
Ran 14 tests
OK

python -m unittest discover -s tests/derivation -v
Ran 35 tests
OK

git diff --check
exit 0
```

Markdown 扫描覆盖 03/04 两章和证据账本：除 tab/LF/CR 外无控制字符；display/inline dollar delimiter 数为偶数；LaTeX 花括号总数配平。最终提交前清除 `__pycache__`、`.pyc` 和 `.superpowers/sdd/pdf-checks`。

## 不作的主张

- 不把 `sigma_eff` 当成新的 nominal TRT 输入。
- 不把 frozen-pressure modal identity 外推成变系数四阶等效方程。
- 不声称 `s_f_minus`、`s_g_plus` 对边界或高阶误差没有影响。
- 不声称二阶恢复消除了列出的 `O(epsilon^3)`、`O(Ma^3)` 或更高项。
- 未运行 Fortran benchmark，因为本任务没有修改生产数值求解器。
