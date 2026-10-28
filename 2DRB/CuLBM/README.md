# 二维均匀网格 CuLBM，单 GPU OpenACC

本目录实现 Molla（2025）[Development and numerical validations of a two-dimensional cumulant-driven thermal lattice Boltzmann method](https://doi.org/10.1016/j.rineng.2025.106225) 的牛顿流体热对流路线：流场 D2Q9、温度 D2Q5。侧热和 RB 使用项目现有的封闭腔体、均匀格心网格及单 GPU OpenACC 风格。

原文存在归一化、力项和逆变换的不一致。本实现按累积量生成函数及矩变换定义进行数学一致的补全，详见 `FORMULATION_CN.md`。当前没有作者原代码或勘误来确认这些补全与作者实际程序完全一致，因此交付不等于逐式照抄印刷公式或已经复现全部数值表格。

## 文件与框架

- `2DRBOpenaccCuLBM.F90`：独立完整求解器。
- `FORMULATION_CN.md`：累积量、力、正温度约束及原文差异。
- `verify_local.py` / `validation_local.json`：局部数学、编译、边界、小网格运行及严格续算。
- `verify_transport.py` / `transport_validation_local.json`：实际流场/温度例程的周期小幅波输运系数检查。
- `PROVENANCE.json`：来源、源码校验值及验证范围。

二维结构和诊断依据 `均匀网格/2DRBOpenmp.F90` 核对；实际 GPU 数据区、异步循环和输出风格参照同维 `均匀网格/2DRBOpenacc.F90`。主文件复用 CLBM 已整理的双分布框架，并移植 HybridRLB4 中稳态参考场续算、Nu 标度分离等修复；两套碰撞和初始化均重新实现为真正累积量方法。原有求解器均保持不变。

## 算例和参数

默认 `steadyFlow + SideHeatedCell`：左热右冷、上下绝热、四壁无滑移。RB 改选 `unsteadyFlow + RayleighBenardCell`，边界自动变为下热上冷、左右绝热、四壁无滑移；初始温度为导热态加确定性扰动。两组宏各启用一个。

参数直接在 `commondata` 修改，保留项目命名及中文注释：

| 参数 | 默认值或含义 |
| --- | --- |
| `nx,ny` | `256,256`，均匀格心网格，格距和格子步长均为1 |
| `Rayleigh,Prandtl` | `1e6,0.71` |
| `Mach` | `sqrt(3)*0.01`，自由落体格子速度 `Uff=0.01` |
| `Thot,Tcold` | `2,1`，两壁和内部载体温度均要求为正 |
| `cs2,thermalCs2` | 固定格子常数 `1/3,1/5` |
| `cuViscousRate` | `1/tauf`，三个流场二阶累积量均采用此率 |
| `cuThirdRate,cuFourthRate` | 论文经验值 `0.6,1.4`，不是任意算例的最优保证 |
| `cuThermalFluxRate` | `1/taug`，`taug=0.5+kappa/thermalCs2` |
| `cuThermalSecondRate` | `1`，按论文设置，不套用 CLBM 的其他松弛关系 |
| `steadyCheckIntervalItc` | `2000`，按绝对累计步数检查 |
| `epsU,epsT` | `1e-7`；速度相对 L2 变化、平均温度 L1 变化/温差 |
| `itc_max` | 稳态总目标上限 `20000000`，是否收敛单独记录 |
| `unsteadyRunDuration` | 非稳态总目标 `1000 t_ff`，续算只补剩余时长 |

侧热特征长度沿用 `L=nx`，RB 为 `L=ny`。黏度、扩散率、浮力和自由落体时间由 Ra、Pr、Ma 及温差耦合推导。

这里的温度归一化需要 **严格正的载体 T**。不能直接把其他程序的 `Thot=1,Tcold=0` 或负温标搬入本版本。启动时检查 `Thot>Tcold>0`，设备碰撞检查内部 `T>0`。失败时设置持久错误标记，在检查及输出边界明确停止；原样复制分布只是避免设备除零的错误占位，不是温度裁剪或恢复稳定。程序不要求每个 `g_i>0`，也不保证温度自动保持在两壁值之间。

## 推进、输出和续算

两套碰撞读取同一时刻的宏观场；两套迁移和边界后，先恢复新温度，再恢复含新浮力半步修正的物理速度。流场高阶累积量先围绕实际 raw 均速提取，碰撞后按完整受力均速逆变换，保留完整动量增量。没有将 CLBM 的固定四阶中心矩冒充累积量，也没有叠加旧 MRT、useG 或滤波。

保留稳态误差曲线、Nu/Re、全场/冷热壁/中线 Nu、壁面及中心线极值、流函数/涡量、Tecplot、二进制宏观快照、非稳态历史及窗口平均。输出前缀含 `2DOpenaccCuLBM`。Nu 按真实温差及 `L/kappa` 归一化，与速度显示参数 `velocityScaleCompare` 分离。

`loadInitField=1` 读取本版本配套的 latest `.meta` 和 `.bin`，并要求 Nu/Re 历史对齐。检查点标识 `CULBM2D_D2Q9_D2Q5_V1`，元数据版本 `901`，保存17个物理/时间轴参数。二进制包含 `f/g`；稳态另存 `up/vp/Tp` 和误差，写出前同步设备参考场。恢复时在执行设备上重新由 `f/g` 求宏观量。

非采样点续算保持原采样相位，旧历史保留并追加；已达到目标或已经满足当前稳态阈值的状态再次启动，不重复写同一时刻历史。目标时长及阈值可以调整，算法、网格、物理参数和采样定义必须匹配。Fortran sequential unformatted 格式不保证跨编译器兼容。

## 验证与限制

单 GPU 编译示例，本机尚未执行 NVHPC：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccCuLBM.F90 -o culbm
```

本地脚本依赖 Python/NumPy 和 GNU Fortran，实际启用 OpenACC host、数组边界和浮点异常检查。临时构建目录须是可写 ASCII 路径，例如：

```text
python verify_local.py --compiler C:/msys64/ucrt64/bin/gfortran.exe --temp-root C:/Temp
python verify_transport.py --compiler C:/msys64/ucrt64/bin/gfortran.exe --temp-root C:/Temp
```

以 JSON 的 `passed`、源码 SHA-256 和具体测试参数为准。真正热累积量含非线性的热流平方项，默认热松弛率下，线性导热一般不是机器精度的离散不动点。边界精确导热核验会单独采用 `cuThermalFluxRate=1` 的配置并记录，不能把其结果当成默认参数已被验证。程序不承诺有限网格下任意温标平移完全等价。

本次 `64×6`、2000步、小幅周期波独立测试测得黏度误差约 `0.0249%`、热扩散率误差约 `0.0281%`。周期链只由测试驱动补齐，生产算例仍为封闭腔体。

本地 `12×10`、48步的背景温标观察中，将两壁同时加 `3.25` 而保持温差不变，侧热/RB 的最大 Nu 变化分别约 `0.00528/0.00556`；并非任意温标平移精确不变。相同小网格、600步、常规热率的导热观察也出现冷热壁 Nu 偏离1，数值已保留在报告的 `default_rate_observation` 中；这些未收敛短程观察不作为稳态基准或网格无关结论。整体温度及分布同倍缩放、相应调整浮力系数的等价性检查则通过。

这些是实现和输运层面的开发检查。尚未验证 NVHPC/NVIDIA GPU、论文对流基准、高 Ra 长时间稳定性、长期 RB 统计或网格无关性；本次没有访问服务器或提交任务。多块网格工作保持搁置。
