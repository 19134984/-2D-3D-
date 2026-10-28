# 二维均匀网格混合 RLB-4，单 GPU OpenACC

本目录按 Nee（2025）[Evaluation of hybrid lattice Boltzmann models for laminar and turbulent natural convection](https://doi.org/10.1016/j.icheatmasstransfer.2025.108635) 实现高阶正则化 D2Q9 流场与宏观温度方程的混合方法。当前交付二维均匀网格代码和本地开发核验，服务器和实际 NVIDIA GPU 测试留待后续。

## 文件与算例

- `2DRBOpenaccHybridRLB4.F90`：完整单文件求解器。
- `FORMULATION_CN.md`：公式、原文缺项和本实现补全的说明。
- `verify_local.py`、`validation_local.json`：局部公式、温度离散、编译、小网格运行及续算检查。
- `verify_transport.py`、`transport_validation_local.json`：使用实际流场例程核验剪切波衰减黏度。
- `PROVENANCE.json`：论文与参考代码的 SHA-256 来源记录。

主框架、参数和后处理组织按仓库二维核心文件核对；实际 GPU 数据管理与循环参照 `均匀网格/2DRBOpenacc.F90`。同时复用 FPC 衍生文件已经整理的输出和严格续算框架，重新实现本方法的碰撞、速度定义及整个温度推进部分。原始均匀网格、FPC、CLBM 和 USTLBM 文件保持不变。

默认启用 `steadyFlow + SideHeatedCell`：左热右冷、上下绝热、四壁无滑移。切换 RB 时启用 `unsteadyFlow + RayleighBenardCell`，其边界自动切为下热上冷、左右绝热、四壁无滑移。两组宏分别只能选一个。RB 从线性导热态加确定性小扰动启动。

## 参数与时间推进

主要参数直接在 `commondata` 修改，并保留中文注释：

| 参数 | 默认值和含义 |
| --- | --- |
| `nx,ny` | `256,256`，格心均匀网格，格距与流场时间步均为 1 |
| `Rayleigh,Prandtl` | `1e6,0.71` |
| `Mach` | `sqrt(3)*0.01`，对应自由落体格子速度 `Uff=0.01` |
| `Thot,Tcold` | `1,0`；浮力按 `T-Tref`，Nu 按真实温差归一化 |
| `thermalSubsteps` | `1`；每个流场步内进行多少个完整温度 RK4 子步 |
| `steadyCheckIntervalItc` | `2000`，以绝对累计步数检查稳态误差 |
| `epsU,epsT` | `1e-7,1e-7`，分别为相对速度 L2 变化和温差归一化的平均温度 L1 变化 |
| `itc_max` | 稳态总目标上限 `20000000`，达到上限不等于收敛 |
| `unsteadyRunDuration` | 非稳态总目标 `1000 t_ff`，续算只补足到总目标 |
| `outputSnapshotInterval` | 稳态 `10 t_ff`；非稳态 `0.5 t_ff` |
| `unsteadyAverageStartTf/EndTf` | 非稳态后处理窗口，默认后半段 |

侧热特征长度沿用项目定义 `L=nx`，RB 为 `L=ny`。若以后改成长方形并与论文高腔比较，需要先核对论文的 Ra 特征长度、纵横比和时间尺度；默认方腔不是论文所有高 Ra 算例的同一几何。

主循环为：流场碰撞、迁移、无滑移反弹、恢复 raw 动量速度，然后冻结这一步速度完成温度 RK4。温度整步状态是 `T`，没有第二套 `g` 分布。检查点保存 `f/T`，临时 RK 数组不保存。

本程序的“4”指四阶 Hermite 重构及温度经典 RK4；温度内区空间离散补齐为中心二阶。原文没有给齐 RLB-4 的空间模板、阶段耦合和壁面细节，不能称整个流热耦合已达到四阶时间或空间精度。补全选择见 `FORMULATION_CN.md`。

`thermalSubsteps` 仅细分温度步，不改变流场时间步、输出时间或物理参数。启动时按二次 ghost 闭合的保守扩散谱界检查 `10*kappa/thermalSubsteps<=2.7`；这只是纯扩散限制，不保证包含对流项的全方法稳定。当前没有额外滤波、温度裁剪、上风耗散或亚格子模型。

## Nu/Re、后处理与续算

保留原项目的稳态误差曲线、全场/冷热壁/中线 Nu、壁面与中心线极值、流函数和涡量、Tecplot、二进制宏观快照、Nu/Re 历史及窗口统计。输出前缀含 `2DOpenaccHybridRLB4`。

`NuVolAvg=1+L*<u_heat*(T-Tref)>/(kappa*DeltaT)`，侧热 `u_heat=u`，RB 为 `v`；`ReVolAvg=L*sqrt(<u²+v²>)/nu`。这里的速度按论文式(8)取 raw 动量，未额外加半步力。`velocityScaleCompare` 只影响速度显示，Nu 的物理热扩散标度固定为 `L/kappa`。

结束时另写 `ThermalBudget_2DOpenaccHybridRLB4.txt`，报告当前瞬时温度方程的总 RHS、向外壁热流、两者之和及离散速度散度。原文的对流形式在离散速度不严格无散时不保证全离散总热量守恒；该残差用于量化这一问题，不自动判定失稳或统计收敛。

续算时设置 `loadInitField=1`，保留本模型配套的 `.bin`、latest `.meta` 和 Nu/Re 历史。标识为 `HYBRID_RLB4_D2Q9_T_V1`，元数据版本 `801`，含14个物理/时间轴参数。稳态检查点还保存 `up/vp/Tp` 和已有误差，读取后恢复原检查相位；写出前同步设备上的参考场。目标总时长和收敛阈值可调整，模型、网格、物理参数和采样间隔必须匹配。非采样点续算不改变采样相位；旧历史保留并追加，不接受缺失或超前的历史。

已达到目标或已满足稳态阈值的检查点再次启动时不重复写同一时刻历史。检查点采用 Fortran sequential unformatted，不保证跨编译器格式兼容。

## 编译与验证边界

单 GPU NVHPC 编译示例，本机尚未执行：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccHybridRLB4.F90 -o hybrid_rlb4
```

本地开发验证用 GNU Fortran 的 OpenACC host 后端，开启数组边界及浮点异常检查。验证脚本需要 Python 和 NumPy，仅改 ASCII 临时目录里的小网格副本，不改正式源参数。可运行：

```text
python verify_local.py --compiler C:/msys64/ucrt64/bin/gfortran.exe --temp-root C:/Temp
python verify_transport.py --compiler C:/msys64/ucrt64/bin/gfortran.exe --temp-root C:/Temp
```

以报告中的 `passed` 和源码 SHA-256 为实际核验记录。局部核、纯导热、有限步输出及严格续算通过，不代表论文高 Ra 基准通过；温度单独 RK4 时间检验不代表全耦合时间阶。NVHPC 编译、NVIDIA GPU 执行、论文对流基准、高 Ra 长时稳定性、RB 统计收敛和网格无关性仍待验证。本次没有访问或提交服务器任务。

本次独立周期剪切波测试专门设 `64×6`、`nu=0.05`、速度幅值 `1e-5`，推进2000步，拟合得到 `nu=0.05001968248`，相对误差约 `0.0394%`。周期填充只存在于测试驱动，生产算例仍为封闭腔体。
