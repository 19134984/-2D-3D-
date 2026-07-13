# Task 04 实施报告：D2Q5 四阶参考生成器

## 结果

- 起始提交：`5b8d2456adcdbd0a699404ed2563785f78ad3985`
- 提交信息：`test: reproduce D2Q5 fourth-order reference`
- 范围：仅新增 D2Q5 推导工具、测试、附录、证据账本和本报告；未修改生产 Fortran、D2Q9 目标公式或原始 PDF。

## 交付文件

- `tools/derivation/series.py`
- `tools/derivation/d2q5_reference.py`
- `tests/derivation/test_d2q5_reference.py`
- `docs/derivation/chapters/A-d2q5-verifier.md`
- `docs/derivation/evidence-ledger.md`
- `.superpowers/sdd/task-04-report.md`

## Minor 审查修复

在 `1dab741` 基础上补强了两个审查项，未修改任何推导系数或生产代码：

- 附录明确区分 Dubois Eq. (15) 与第 17 页正相位 $B$ 的相反空间位移符号。本验证器固定采用第 17 页约定；两者以 $\boldsymbol k\mapsto-\boldsymbol k$ 对应，所以当前纯扩散偶数阶系数不变。Task 5 的源项或平流符号必须在其自行固定的 Fourier 约定下重新推导。
- Route B 禁词门禁同时扫描 `taylor_moment_route` 与 `_streaming_taylor_residual`，防止独立性只在公开入口表面成立。
- monkeypatch 门禁改为一般有理点上的 `order=4` 调用，并精确断言 `diffusion=7/100`、`kappa40=26483/3850`、`kappa22=-25317/1925`。

## 原始证据核对

直接提取并渲染检查了本地 Dubois--Lallemand PDF：

- 第 11 页：Eq. (39) 的 `Psi` 和 Eqs. (40)-(42) 的扩散/四阶归一化；
- 第 16 页：Eq. (55) 和中间 TRT 特殊点；
- 第 17 页：平面波 `exp(i k.x)`、正相位 `B`、`G=B M^-1 Psi M`；
- 第 38 页：Appendix Eq. (79) 的 D2Q5 矩阵和速度顺序。

另核对 Wang PDF 第 4 页 Eqs. (14)-(19) 与 Luo PDF 第 4 页 Eqs. (13)-(21)。两文的 `sigma_kappa,sigma_e,sigma_nu` 是碰撞率，映射到 Dubois `s1,s3,s4`；对应的 Hénon 参数必须再作 `1/s-1/2` 转换。

## 参数化模型

验证器没有复用 Task 1 固定 `alpha=-2/3` 权重作为一般模型。它从

`m_eq = rho*(1,0,0,alpha,0)`

和 Eq. (79) 的矩阵精确反解：

- `w0=(1-alpha)/5`；
- `w1=w2=w3=w4=(4+alpha)/20`；
- `c_e=(4+alpha)/10`。

测试逐项验证 `M*M.inv()==I`、Eq. (39) 的碰撞矩阵、平衡人口重构以及 `alpha=-2/3` 才回到 Task 1 的 `1/3,1/6` 固定点。

## 两条独立路线

### Route A：放大矩阵级数

Route A 构造原文正相位 `G=B M^-1 Psi M`，但不使用浮点特征值生成符号答案。它展开 `G_n,z_n,v_n`，固定守恒规范后，以一个 6x6 增广系统逐阶跟随 `z_h(0)=1` 的根，最后从 `-log(z_h)/dt` 取得二阶与四阶齐次项。

### Route B：物理空间 Taylor/矩递推

Route B 不构造 `G`，也不使用特征值、特征向量、特征多项式、`z_h` 或对数。它从离散 streaming 的物理空间 Taylor 展开建立 differential jets，逐阶解

`(I-Psi) m^[n] + m_eq L_(n-1) = -R_n,  m^[n]_0=0`

以消去非守恒矩，得到 `L0=0`、二阶扩散、`L2=0` 和四阶空间算子。

两路只共享矩阵、速度和通用截断多项式工具。测试还 monkeypatch `amplification_matrix` 和 `amplification_route` 为抛异常函数，Route B 仍能独立返回一般有理点上的精确二阶与四阶系数；公开入口和内部 streaming Taylor helper 也同时接受禁词扫描。

## 精确结果

两路线的符号差严格为零，并分别复现：

`kappa40 = 8 - 3*alpha + 12*(alpha+4)*sigma1**2 - 12*(1-alpha)*sigma1*sigma3 - 60*sigma1*sigma4`

`kappa22 = -6*(alpha+4) + 24*(alpha+4)*sigma1**2 - 24*(1-alpha)*sigma1*sigma3 + 120*sigma1*sigma4`

有量纲系数恢复为：

- `D=(alpha+4)*lambda**2*dt*sigma1/10`；
- `C40=lambda**4*dt**3*sigma1*(alpha+4)*kappa40/1200`；
- `C22=lambda**4*dt**3*sigma1*(alpha+4)*kappa22/1200`。

内部无量纲推导只用于减少符号规模；出口恢复上述 `lambda,dt` 尺度，未改变物理归一化。

## 各向同性与完整消除

精确计算得到

`kappa22-2*kappa40 = 40*(6*sigma1*sigma4-1)`。

所以 `sigma1*sigma4=1/6` 只保证四阶误差各向同性。测试选取满足该关系但 `kappa40` 非零的点，明确防止把各向同性误写成完整消除。

同时解 `kappa40=kappa22=0` 才得到 Eq. (55)：

- `sigma3=sigma1*(alpha+4)/(1-alpha) - (2+3*alpha)/(12*sigma1*(1-alpha))`；
- `sigma4=1/(6*sigma1)`。

加上 `sigma3=sigma4` 后得到 `sigma1=1/sqrt(12)`、`sigma3=sigma4=1/sqrt(3)`；测试验证它对符号 `alpha` 同时消掉两个系数。

## 高精度数值抽查

在三组一般有理参数点上，以 80 位精度和 `q=1e-6` 直接计算放大矩阵的 hydrodynamic 特征值。轴向检查 `(Gamma(q,0)-D*q**2)/q**4`，对角检查 `(Gamma(q,q)-2*D*q**2)/q**4`。

最大轴向差约 `1.10e-15`，最大对角差约 `9.91e-16`。数值特征值只作抽查，没有参与符号答案生成。

## TDD 记录

### 初始 delayed-import RED

`python -m unittest tests.derivation.test_d2q5_reference -v`

- `Ran 1 test in 0.001s`
- `FAILED (failures=1)`
- 捕获值：`ModuleNotFoundError: No module named 'tools.derivation.d2q5_reference'`

异常被 delayed import 捕获后表现为 assertion `FAIL`，没有测试收集、解析或 import `ERROR`。添加最小 callable 骨架后该测试 1/1 通过。

### 参数化模型与级数工具 RED/GREEN

参数化模型测试先因骨架不接受 `lattice_speed` 产生 4 个捕获后的 assertion `FAIL`；实现 Eq. (79)、Eq. (39) 与平衡重构后 4/4 通过。

级数工具先因 `tools.derivation.series` 不存在产生 2 个 assertion `FAIL`；实现总次数、齐次分量和尺度截断后 2/2 通过。

### 二阶路线 RED/GREEN

Route A/B 与相位测试先因新入口不存在产生 4 个 assertion `FAIL`。实现后 4/4 通过；增加 monkeypatch 独立性门禁后为 5/5。

### 四阶路线 RED/GREEN

四阶组最初为 7 个 assertion `FAIL`，原因是两路线明确只支持 `order=2`。扩展递推后，公式、Eq. (55)、TRT 点、三有理点和高精度检查均已通过。

首次四阶运行唯一剩余失败为测试夹具的 SymPy 替换顺序：同一字典同时替换 `sigma4=1/(6*sigma1)` 和 `sigma1=1/5`，留下 `-40+8/sigma1`。这不是两路线的非零残差；改为顺序替换后 7/7 通过。

一般符号递推最初约 112 秒。把内部计算无量纲化、仅在出口恢复物理尺度后，四阶组降到约 54 秒，targeted 全组为 49.793 秒。

## 验证结果

`python -m unittest tests.derivation.test_d2q5_reference -v`

- `Ran 19 tests in 52.000s`
- `OK`

`python -m unittest discover -s tests/derivation -v`

- `Ran 55 tests in 62.081s`
- `OK`

最终还执行 diff check、控制字符扫描、Markdown dollar delimiter 和 LaTeX 花括号平衡检查；结果记录在提交前验证输出中。

## 不作的主张

- 不把 D2Q5 描述成删除 D2Q9 对角速度的连续退化。
- 不把 Wang/Luo 各向同性关系描述成完整四阶消除。
- 不把 Eq. (55) 或 TRT 特殊点迁移到 D2Q9。
- 不把高精度数值抽查作为符号公式的来源。
- 未运行 Fortran benchmark，因为本任务没有修改生产求解器。
