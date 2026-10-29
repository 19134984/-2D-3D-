# 二维均匀网格 Feng 熵稳定温度 LBM，单 GPU OpenACC

主文件为 `2DRBOpenaccEntropicADE.F90`。温度算法对应 Feng、Leng、Jiang、Chu（2026）*Entropic lattice Boltzmann method for general anisotropic advection–diffusion*，arXiv:2605.01774v1，重点为式(7)–(41)。原文在 `../pdf/算法/2605.01774v1.pdf`。

按用户确认，流场采用已经复现的 **Luo/Wang D2Q9 MRT**；温度采用 **Feng D2Q9 标量熵碰撞**。当前只做各向同性热扩散 `D=kappa I`，二维均匀格心网格、单 GPU OpenACC。两套分布分别为流场 `f` 和温度 `g`，`T=sum(g)`；不能与 `EntropicLBM/` 的 Frapolli 总能量分布混同。

数值和二维诊断核对 `均匀网格/2DRBOpenmp.F90`，流场碰撞从同维 `均匀网格/2DRBOpenacc.F90` 保留。设备常驻、输出和严格续算外壳复用其同框架后续实现 `FPC/2DRBOpenaccFPC.F90`，这里重新编写全部温度碰撞；没有修改原有求解器。差异和原文缺项见 `FORMULATION_CN.md`。

## 算例与参数

在文件顶部两组宏中各开启一项：

- 侧壁差温稳态：`steadyFlow` + `SideHeatedCell`，这是默认配置；左热右冷、上下绝热、四壁无滑移。
- RB 非稳态：`unsteadyFlow` + `RayleighBenardCell`；下热上冷、左右绝热、四壁无滑移，初始线性温度加确定性扰动。

壁面宏按算例自动匹配。生产版本暂不提供周期侧壁或多块网格；原文二维 RB 是水平周期、长高比2，与这里的封闭方腔不能直接混作同一基准。

主要参数直接修改 `commondata`，沿用原有中文注释和命名：

| 参数 | 默认值与含义 |
| --- | --- |
| `nx,ny` | `256,256`，全部 `1:nx,1:ny` 为流体格心，壁面在半格距离处 |
| `Rayleigh,Prandtl` | `1e6,0.71` |
| `Mach` | `sqrt(3)*0.01`，自由落体格子速度 `Uff=0.01` |
| `Thot,Tcold` | `2,1`；本模型需要严格正的温度载体及温度分布 |
| `Snu,Sq` | 沿用 Luo 流场率，二阶能量／四阶矩取 `Snu`，奇矩取 `Sq` |
| `fengFluxRate` | `1/(0.5+3*kappa)`，控制目标热扩散率 |
| `fengPositivityFloor` | `1e-14`，相对于本格最小正分布的正支阈值 |
| `fengGhostTolerance` | `64*epsilon(real8)`，只保护归一化 ghost 的机器零范数 |
| `steadyCheckIntervalItc` | `2000`，稳态误差和场检查间隔 |
| `rbInitialPerturbation` | `1e-3`，相对于壁面温差 |
| `epsU,epsT` | `1e-7,1e-7`；相对速度 L2 变化、温差归一化平均温度 L1 变化 |
| `unsteadyRunDuration` | 非稳态总目标 `1000 t_ff`；默认后半段为平均窗口 |

侧热特征长度为 `nx`，RB 为 `ny`；黏度、热扩散率、浮力和自由落体时间按 Ra/Pr/Ma/真实温差联动。默认网格和时长只是可编辑设置，不代表已经完成该工况的长时间验证。

## 碰撞与耦合

主循环保留项目的 `initial`、`collision`、`collisionT`、`streaming`、`bounceback`、`streamingT`、`bouncebackT`、`macroT`、`macro` 组织。两个碰撞都读取同一时刻的宏观场；迁移和边界之后先恢复新温度，再恢复含新浮力半步修正的物理速度。

流场是原 Luo MRT 碰撞及矩空间 Guo 型力；没有额外流场熵稳定器、FPC、正则化或 LES。初始化给 `f` 半步力补偿，使初始物理速度为零；迁移后的速度使用新温度计算半步浮力。这两项时间一致性处理与核心旧推进次序存在差别，不能把本实现与旧程序的差异全部归于温度碰撞。浮力保留项目的 `rho*gBeta*(T-Tref)`，与原文使用导热参考剖面的写法差异在公式说明中公开记录。

温度核 `feng_thermal_collision` 保留式(38)的 ADE 修正项，不把它替换为 Frapolli 的 Newton 熵求解，也不把 lambda 裁剪到固定区间。正支回退对完整碰撞增量施加同一 alpha。它可改变局部有效扩散，因此增加回退次数、最小 alpha 和 lambda 范围的诊断，不能把“回退后不发散”等同于目标热扩散率始终不变。

半程恒温使用 ABB，绝热使用保留切向来源的镜面反射。碰撞保正不能保证 ABB 之后仍为正；检查到非正温度分布时会明确报错，不进行温度／分布裁剪来掩盖失稳。正温标也不能保证所有参数稳定。

## 输出、Nu/Re 与续算

输出前缀均带 `EntropicADE`，与原有算法分开。保留稳态误差曲线、体积／冷热壁／中线 Nu、壁面 Nu 极值、中心线速度极值、流函数／涡量、Tecplot、`u/v/T/rho` 二进制快照、Nu/Re 历史和非稳态窗口统计。

- `NuVolAvg=1+L*<u_heat*(T-Tref)>/(kappa*DeltaT)`；侧热取 `u_heat=u`，RB 取 `v`。
- `ReVolAvg=L*sqrt(<u^2+v^2>)/nu`；时间窗口的 RMS Re 为瞬时空间 RMS Re 先平方、时间积分后开方。
- `velocityScaleCompare` 只改变速度类显示单位，Nu 始终按实际 `kappa` 和壁温差归一化。
- `EntropyHistory2DOpenaccEntropicADE.txt` 随每个 Nu/Re 样本输出累计回退次数／比例、全程最小 alpha、当前 lambda 极值、最小温度／温度分布及质量。它用于检查稳定化的活跃程度，不是额外的物理耗散模型。

续算设置 `loadInitField=1`，保留同目录的 latest `.meta`、对应 `.bin`、Nu/Re 及熵诊断历史。历史缺失、错位或超过检查点会拒绝续算，需先恢复对应检查点的文件组合。二进制保存 `f/g`、熵诊断状态；稳态还保存上次检查的 `up/vp/Tp` 和 `errorU/errorT`。算法标识为 `FENG_ADE_LUO2D_V1`，metadata版本1001；不同模型、网格、参数或采样时间轴不能混用。允许增加总目标时长、调整收敛阈值；非采样点续算保持原采样相位。快照仅供后处理，不能用于严格续算。

## 本地核验

`verify_local.py` 使用实际生产源中的碰撞例程及完整推进流程，检查公式、守恒、两方向导热、宏组合、输出及严格续算。`verify_transport.py` 在独立周期驱动中使用生产碰撞核测输运系数；周期条件只用于核验，不改变封闭生产算例。`verify_benchmark.py` 使用实际完整求解器检查低 Ra 侧热稳态。

```text
python3 verify_local.py --compiler /path/to/gfortran
python3 verify_transport.py --compiler /path/to/gfortran
python3 verify_benchmark.py --compiler /path/to/gfortran
```

临时构建和小网格参数只作用于隔离副本；摘要保存在验证 JSON 中并带生产源 SHA-256。NVIDIA 编译示例：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccEntropicADE.F90 -o entropic_ade
```

2026-10-29 本机核验通过：四种宏组合、局部公式和矩、正支回退、两方向导热、小网格完整推进和严格续算；周期输运的黏度／热扩散率相对误差均小于 0.04%。32×32、Ra=1000、Pr=0.71、Ma=0.05 侧热在第16000步满足两项稳态阈值，体积 Nu=1.11813，对照 Wang/Luo（2013）表2／9的1.1178，偏差约0.03%；冷热壁 Nu 分别为1.11829和1.11829。各项设置和原始指标见三个验证 JSON。

上述检查使用 GNU Fortran OpenACC host 后端。它们不等于 NVHPC/NVIDIA GPU 执行验证，也不等于高 Ra DNS、长期 RB 统计收敛或网格无关性验证；目前没有实测证明其高 Ra 稳定性优于已有代码。本轮不登录服务器，不提交 Git 或推送远端。
