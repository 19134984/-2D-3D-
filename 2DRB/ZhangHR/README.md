# Zhang 2022 HR 温度 LBM：二维均匀网格、单 GPU OpenACC

主程序：`2DRBOpenaccZhangHR.F90`。

论文为 Zhang、Li、Wu（2022），*Advection–Diffusion Lattice Boltzmann Method With and Without Dynamical Filter*，Frontiers in Physics 10:875628，[出版社全文](https://doi.org/10.3389/fphy.2022.875628)。正式原文位于 References 的 `05_论文原文PDF/Zhihong Zhang-Yunke Wu-2022-Advection–Diffusion Lattice Boltzmann Method With and Without Dynamical Filter.pdf`；仓库副本为 `../pdf/算法/fphy-10-875628.pdf`。这篇与同作者另一篇 *A hybrid regularized lattice Boltzmann model for convection–diffusion equation* 是两篇文章。

本程序复现论文 **HR 标量碰撞及可选动态滤波**，并按本项目要求组合 **Luo/Wang D2Q9 MRT 流场**，扩展到侧热和 RB。原文流场是 RR、算例是周期被动标量；这里不声称完整重现原文全部流场和图表。公式、原文缺项及项目适配见 `FORMULATION_CN.md`。

数值定义参考 `均匀网格/2DRBOpenmp.F90`，流场和 GPU 风格参考 `均匀网格/2DRBOpenacc.F90`；已经核验的 `EntropicADE/` 只作为同框架的设备常驻、输出与严格续算外壳来源。新的温度核没有保留 Feng 熵碰撞或正支回退。原有求解器未修改。

## 算例及主要参数

文件顶部两组宏各开启一个，修改后重新编译：

- `steadyFlow` + `SideHeatedCell`：默认侧热稳态，左热右冷、上下绝热。
- `unsteadyFlow` + `RayleighBenardCell`：RB 非稳态，下热上冷、左右绝热。

四壁无滑移，均匀半格壁面；当前不提供生产周期侧壁、多块网格、MPI 或多 GPU。

参数仍直接写在 `commondata`，沿用项目命名与中文注释：

| 参数 | 默认值和含义 |
| --- | --- |
| `nx,ny` | 256×256 格心网格 |
| `Rayleigh,Prandtl` | 1e6、0.71 |
| `Mach` | `sqrt(3)*0.01`，自由落体格子速度0.01 |
| `Thot,Tcold` | 0.5、−0.5，沿用核心温标；本模型允许带符号的温度及分布 |
| `taug` | `0.5+3*diffusivity`，对应物理热扩散率 |
| `zhangThermalModel` | 2：HR；1：PR1st；0：带相同时间修正的 BGK，便于对照 |
| `zhangBeta` | 0.95；PR 权重，`1-beta` 为差分估计权重；只影响 HR |
| `zhangEnableFilter` | 0，默认不加动态滤波；1 启用论文宏观滤波 |
| `zhangFilterSigma0,zhangFilterEpsilon` | 0.2、0.1，滤波经验参数 |
| `steadyCheckIntervalItc` | 2000，稳态误差与场检查间隔 |
| `epsU,epsT` | 1e−7、1e−7，速度相对L2变化、按实际温差归一化的平均温度L1变化 |
| `rbInitialPerturbation` | 1e−3，确定性初温扰动相对壁面温差 |
| `unsteadyRunDuration` | 1000个自由落体时间；默认后半段统计 |

Ra、Pr、Ma、特征长度、黏度、扩散率和浮力系数联动，侧热长度取 `nx`、RB 取 `ny`。默认参数不是已完成长期验证的高 Ra 工况。与 Feng 对比时还应统一温标，不能只改算法名。

## 方法及动态滤波开关

温度使用一阶平衡态 `geq=w*T*(1+c.u/cs2)`，从现有分布提取 PR 通量，与 `-cs2*taug*grad(T)-0.5*dt(Tu)` 混合；`grad(T)` 用二阶中心差分，`dt(Tu)` 用当前与上一时层的格子一步差分。初始化历史设为初始 `T*u/T*v`，故首步时间差分为零。温度无物理热源；该时间修正不是额外加热。

两碰撞使用同一时刻宏观场，再迁移和施加壁面，先重建温度，再重建含新温度浮力半步的物理速度。初始化流场含 `−F/2` 补偿，与已有 Feng 版保持一致；这是相对核心旧推进的适配，非原文 RR 流场。

默认关闭滤波，适合作为 HR 的 DNS 精度对照。开启后，在宏观温度恢复之后按原文三点模板同时滤两个方向，传感器采用前向一阶梯度，尺度按正文取 `epsilon*(Tmax-Tmin)`；恒定场关闭滤波。恒温墙外按半格温度延拓，绝热墙外作镜面延拓，这是原文周期算例之外的项目适配。

滤波只改变宏观 `T`，不修改各方向 `g`；此时 `T` 与 `sum(g)` 可不同，下一步以滤后温度计算平衡态、PR/FD 通量及浮力。原文没有提供 population 同步规则，本实现明确保留宏观滤波的接口，不补入未声明的分布修正；独立温度和通量历史均进入检查点。

动态系数滤波一般不严格保持总标量，且增加数值耗散。不能凭更平滑或不发散就称为 DNS 精度改善；需要比较滤波开关、参数敏感性、热收支和网格收敛。本程序没有 LES 模型，也未声称任意工况保正或稳定；不裁剪温度和分布来掩盖失稳。PR2nd 因原文表达存在歧义，本版不加入该模式。

## 输出与严格续算

保留项目的稳态误差、体积／冷热壁／中线 Nu、壁面极值、中心线速度极值、流函数／涡量、Tecplot、u/v/T/rho 快照、Nu/Re 时间序列及非稳态窗口统计。文件前缀全部为 `ZhangHR`。

- `Nu=1+L*<u_heat*(T-Tref)>/(kappa*DeltaT)`，侧热取 u、RB 取 v。
- `Re=L*sqrt(<u²+v²>)/nu`；非稳态时间 RMS Re 对瞬时 Re 先平方，按实际时间积分平均后开方。
- `ZhangHistory2DOpenaccZhangHR.txt` 的列为步数、实际 t_ff、当前 sigma 最大值、全域累计滤波净温度改变量、累计绝对改变量、当前 `sum(g)`、当前 `sum(T)`、T最小值、T最大值、流场质量。

滤波开启时，常规 Nu 仍按输入的物理 kappa 定义，尚未把滤波算子换算成壁面热通量；需结合单独热量诊断，不能以该 Nu 单独证明完整热预算闭合。

续算时设置 `loadInitField=1`，配套保留 latest `.meta`、对应 `.bin`、Nu/Re 和 Zhang 诊断历史。检查点保存 f/g、独立宏观 T、上一步 Cu、逐格累计滤波改变量及 sigma；稳态还保存上次检查参考场和误差。标识为 `ZHANG_HR_LUO2D_V1`，metadata版本1002，不能读取其他模型检查点。网格、物理参数、模型、beta、滤波开关／参数及采样轴必须匹配；可以延长总时长、调整收敛阈值。缺失、超前或错位历史会拒绝续算。快照不是续算文件。

## 编译及本地核验

NVIDIA 编译示例，需具备 NVHPC 环境：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccZhangHR.F90 -o zhang_hr
```

本地开发核验脚本使用支持 OpenACC 的 GNU Fortran host 后端：

```text
python3 verify_local.py --compiler /path/to/gfortran
python3 verify_transport.py --compiler /path/to/gfortran
python3 verify_benchmark.py --compiler /path/to/gfortran
```

核验仅修改隔离副本，JSON 摘要保存被测生产源 SHA-256。实际通过项目及指标以报告为准。尚未在 NVIDIA GPU/NVHPC 实跑，未完成高 Ra RB 稳定窗口、长期统计或网格无关性验证，不能宣称优于现有代码。本轮未登录服务器，未暂存、提交或推送 Git。

2026-10-29 本地结果：

- 公式、beta退化、模型切换、生产梯度、时间历史、导热 Nu=1、宏观滤波接口及四种算例宏组合通过检查。侧热／RB无滤波／RB有滤波的48步连续与32+16步续算检查点和共同快照逐字一致；18个破坏历史用例被拒绝。
- 周期测量的黏度／热扩散率相对误差均小于0.055%；平流温度波最大相位误差约0.00259弧度。这是无滤波、小振幅检查。
- 原文400²、200步 Gaussian 的 Pe=1／1000，BGK、PR1、HR beta=0.5／0共8组均完成。各组 L2 与表1所报 L2 数值的相对差小于本地1%验收阈值；HR beta=0.5的 L2 为0.003428583／0.003107507，对照原表0.0034286／0.0031075。Pe=1 BGK的 L2 为0.02274578，对照0.022573，相差约0.765%；原因尚未确定，不能声称逐位完全重现。其他7组相对差不超过0.041%。
- 32²、Ra=1000、Pr=0.71、Ma=0.05侧热在第16000步收敛，体积 Nu=1.11813451，冷热壁均约1.11825228。相对 Wang/Luo（2013）表2／9的1.1178，体积 Nu 偏差约0.03%。这里只验证一个低 Ra 网格。

没有发现上述差异由后处理或并行装配导致的证据；Pe=1 BGK 与原表的剩余差别可能涉及原文未完整交代的初始化或离散细节，当前只能保留为未确定原因。
