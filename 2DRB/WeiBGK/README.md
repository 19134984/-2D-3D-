# Wei 2015：二维均匀网格双 BGK，单 GPU OpenACC

主程序 `2DRBOpenaccWeiBGK.F90` 对应 Wei、Wang、Yang、Dou、Qian（2015），*A simple lattice Boltzmann model for turbulence Rayleigh–Bénard thermal convection*，[DOI](https://doi.org/10.1016/j.compfluid.2015.06.003)。正式原文位于 References 的 `05_论文原文PDF/Yikun Wei-Hua-Shu Dou and Yuehong Qian-2015-A simple lattice Boltzmann model for turbulence Rayleigh-Benard thermal convection.pdf`。

按用户要求，论文已定义的流场直接采用论文模型：**D2Q9 双 BGK、线性显式浮力、原始动量速度**。没有替换为 Luo MRT、Guo forcing、中心矩、正则化或熵碰撞；也没有加滤波、LES 或温度裁剪。

重要复现边界：魏文只写两套平衡态为多项式，没有给出完整系数、物性换算及 population 级壁面规则。本实现借助原文引用的 Shan（1997）二阶平衡态与标准零相互作用扩散关系补全。具体证据、排版修复和墙面适配见 `FORMULATION_CN.md`，不能把本实现视为已恢复作者全部原始代码。

数值和诊断采用项目 `均匀网格/2DRBOpenmp.F90`、同维 `均匀网格/2DRBOpenacc.F90` 的组织；设备常驻、输出与续算使用已核验 `ZhangHR/` 同框架外壳，并移除其全部 HR／通量历史／滤波状态。已有源文件均未修改。

## 参数和算例

顶部两组宏各开一个，修改后重新编译：

- 默认 `steadyFlow` + `SideHeatedCell`：左热右冷，上下绝热，四壁无滑移。
- `unsteadyFlow` + `RayleighBenardCell`：下热上冷，左右绝热，四壁无滑移。

当前是单块均匀格心方腔，不增加 MPI、多 GPU 或多块网格。魏文第一算例为水平周期通道，第二个为封闭方腔但未交代侧温条件；此处封闭 RB 的绝热侧壁、侧热算例、半程反弹／ABB／镜面温度规则是项目补全，不是逐字复制原文全部几何。

参数直接修改 `commondata`，保留中文注释和原有命名：

| 参数 | 默认值及含义 |
| --- | --- |
| `nx,ny` | 256×256，壁面在流体格心外半格处 |
| `Rayleigh,Prandtl` | 1e6、0.71 |
| `Mach` | `sqrt(3)*0.01`，自由落体格子速度0.01 |
| `Thot,Tcold` | 1、0，沿用魏文温标 |
| `Tref` | 固定0；浮力精确使用 T，不减平均温度 |
| `tauf,taug` | `0.5+3*viscosity`、`0.5+3*diffusivity` |
| `steadyCheckIntervalItc` | 2000，稳态误差及有限性检查 |
| `epsU,epsT` | 1e−7、1e−7，速度相对L2变化、温差归一化平均温度L1变化 |
| `rbInitialPerturbation` | 1e−3，相对壁面温差的确定性初温扰动 |
| `unsteadyRunDuration` | 1000个自由落体时间，默认统计后半段 |

Ra、Pr、Ma、黏度、热扩散率、浮力和时间尺度联动；侧热长度为 nx、RB 为 ny。默认参数不代表已完成该工况高 Ra／长期验证。

## 碰撞和速度定义

原文式(1)–(4)按 BGK 语义实现为：

```text
fpost_i = f_i + (feq_i-f_i)/tauf + 3*w_i*rho*gBeta*T*c_iy
gpost_i = g_i + (geq_i-g_i)/taug
rho=sum(f),  u=sum(c_x*f)/rho,  v=sum(c_y*f)/rho,  T=sum(g)
```

没有 Guo 的 `1-1/(2*tauf)` 前因子，速度没有 `+F/2`，初态 f 也没有 `−F/2` 补偿。两套碰撞读取同一时间层，迁移和边界后恢复温度及原始动量速度。代码 y 轴向上，正浮力对应热流体向上；保留式(3)的 T，因此常量背景浮力也由密度／压力场平衡，不静默改成 T−平均温度。

补全的两套平衡态均为标准二阶 D2Q9 多项式，以 rho 或 T 为相应零阶矩。它没有额外的时间导数修正或热源。极端参数下可以失稳；遇非有限状态或非正密度明确停止，不用裁剪伪装稳定。

## 输出与续算

保留稳态误差、体积／冷热壁／中线 Nu、墙面 Nu 极值、中心线速度极值、流函数／涡量、Tecplot、u/v/T/rho 快照、Nu/Re 时间序列和非稳态窗口统计，前缀均为 `WeiBGK`。

- `Nu=1+L*<u_heat*T>/(kappa*DeltaT)`，侧热取 u，RB 取 v。
- `Re=L*sqrt(<u²+v²>)/nu`；时间 RMS Re 先平方、按实际时间平均，再开方。
- 快照速度为后处理无量纲值，不能用作严格续算。

严格续算设置 `loadInitField=1`，保留配套 latest `.meta`、对应 `.bin` 和 Nu/Re 历史。方案标识 `WEI_BGK_2D_V1`，metadata版本1004、14个物性／采样参数；保存 f/g，稳态还保存上次检查参考场和误差。密度、温度和无半步修正的速度都从 f/g 重建。不读取其他模型检查点，缺失或超前的 Nu/Re 历史会拒绝续算。可延长总目标或修改稳态阈值，保持绝对检查／采样相位。

## 编译和核验

NVIDIA 环境示例：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccWeiBGK.F90 -o wei_bgk
```

本地脚本使用 GNU Fortran OpenACC host 后端，只改隔离副本：

```text
python3 verify_local.py --compiler /path/to/gfortran
python3 verify_transport.py --compiler /path/to/gfortran
python3 verify_benchmark.py --compiler /path/to/gfortran
```

2026-10-29 已通过：四宏组合、逐方向 BGK／显式浮力公式、质量／动量矩、无半步速度、两方向导热 Nu=1、完整侧热和 RB 短程、非对齐严格续算、十个破坏历史／参数／二进制用例。48步连续与32+16步续算检查点和共同快照逐字一致。

周期测得黏度／热扩散率相对误差均小于0.074%。32²、Ra=1000、Pr=0.71、Ma=0.05 侧热第16000步达到稳态阈值，体积 Nu=1.11562330，冷热壁 Nu≈1.11577；相对 Wang/Luo（2013）表2／9的1.1178偏约0.195%。这不是魏文未报告的原始 Nu 数值。原始指标及生产源 SHA256 在三个 JSON 中。

尚未验证 NVHPC／NVIDIA GPU实跑、高 Ra RB 稳定窗口、长时间统计、网格无关性或 DNS 分辨率；也未证明优于原有算法。当前低 Ra 余差更可能与有限网格／显式强迫离散有关，但没有据此确定误差来源；检查未显示后处理或并行装配错误。本轮不登录服务器、不暂存、不提交、不推送 Git。
