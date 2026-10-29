# Fei 2018 双中心矩热 LBM：二维均匀网格、单 GPU OpenACC

主文件 `2DRBOpenaccFeiCLBM.F90`，对应 Fei、Luo 等（2018），*Modeling incompressible thermal flows using a central-moments-based lattice Boltzmann method*，International Journal of Heat and Mass Transfer 120,624–634，[DOI](https://doi.org/10.1016/j.ijheatmasstransfer.2017.12.052)。正式PDF位于 References 的 `05_论文原文PDF/Linlin Fei-Kai Hong Luo-2018-Modeling incompressible thermal flows using a central-moments-based lattice Boltzmann method.pdf`。

按用户确认，**流场和温度场都采用原文模型**：D2Q9 中心矩流场及其一致力源、D2Q5 连续温度中心矩。没有替换成 Luo MRT，也没有叠加 HR、熵回退、滤波、LES 或温度裁剪。二阶流体中心矩20／02通过耦合块松弛，温度保留原始平衡矩中的 `T*ux²`、`T*uy²`。

二维网格／无量纲定义／Nu/Re核对 `均匀网格/2DRBOpenmp.F90`，单GPU组织参考同维 `均匀网格/2DRBOpenacc.F90`。设备常驻、I/O和严格续算外壳复用已核验 `ZhangHR/` 的同一项目框架，已删去其通量历史和滤波状态。只创建本文件夹，没有修改已有求解器。

## 算例与主要参数

文件顶部两组宏各启用一项，修改后重新编译：

- 默认 `steadyFlow` + `SideHeatedCell`：左热右冷、上下绝热、四壁无滑移。
- `unsteadyFlow` + `RayleighBenardCell`：下热上冷、左右绝热、四壁无滑移。

均为二维均匀格心网格，壁面在半格处。侧热是本项目扩展；封闭RB、温标和初扰也不与原文 Table1 混作同一设置。当前不做多块网格、MPI、多GPU或三维。

| 参数 | 默认值及含义 |
| --- | --- |
| `nx,ny` | 256×256 |
| `Rayleigh,Prandtl` | 1e6、0.71 |
| `Mach` | `sqrt(3)*0.01`，自由落体格子速度0.01 |
| `Thot,Tcold,Tref` | 0.5、−0.5、0，沿用项目温标，允许有符号T与g |
| `cs2` | 流场1/3 |
| `thermalCs2` | **温度1/2**；原文 `cT=sqrt(1/2)`，不能误写为cT²=1/4 |
| `Snu` | `(0.5+3*nu)^(-1)`，流体剪切率 |
| `feiBulkRate` | 1，20／02块的体积率；与Snu共同形成s± |
| `feiThirdRate,feiFourthRate` | 1、1，高阶流体中心矩率 |
| `feiThermalRate` | `(0.5+kappa/thermalCs2)^(-1)`，温度通量率 |
| `feiThermalHighRate` | 1，两个二阶温度中心矩率 |
| `steadyCheckIntervalItc` | 2000，稳态误差和场检查间隔 |
| `epsU,epsT` | 1e−7、1e−7，相对速度L2变化及按壁温差归一化的平均温度L1变化 |
| `rbInitialPerturbation` | 1e−3，相对温差的确定性温度扰动 |
| `unsteadyRunDuration` | 1000个自由落体时间，默认后半程统计 |

Ra、Pr、Ma、长度、黏度、扩散率和浮力联动。侧热特征长度为nx，RB为ny。源项按原文参考密度 `rho0=1` 归一化力密度 `Fy=gBeta*(T-Tref)`，不乘瞬时rho；速度仍为 `(sum(c*f)+F/2)/rho`。初始化、碰撞、宏观恢复与重启使用同一定义。

## 边界及算法适配

论文正文采用非平衡反弹速度边界／非平衡外推温度边界，本版按项目统一框架采用半程无滑移反弹、恒温ABB和绝热镜面。这是壁面适配，不能声称逐项重现原文的边界算例。D2Q5只有轴向链，默认静止温度权重为零；负温度或负温度分布本身不是错误，不能套用Feng正支条件。

两个碰撞都读取当前宏观场，迁移和壁面之后先恢复温度，再恢复含新温度力项的物理速度。初始化f含−F/2补偿，使物理初速为零。宏观层的参数／诊断保持项目定义，论文原始算法与项目几何等差异详见 `FORMULATION_CN.md`。

## 输出与续算

输出前缀带 `FeiCLBM`。保留稳态误差、体积／冷热壁／中线Nu、壁面极值、中心线速度极值、流函数／涡量、Tecplot、u/v/T/rho快照、Nu/Re历史和非稳态窗口统计。

- `Nu=1+L*<u_heat*(T-Tref)>/(kappa*DeltaT)`，侧热用u，RB用v。
- `Re=L*sqrt(<u²+v²>)/nu`；窗口时间RMS Re对瞬时Re平方按时间积分，再开方。
- 显示速度标度由 `velocityScaleCompare` 控制，不改变Nu和物理参数。

严格续算设置 `loadInitField=1`，保留 latest `.meta`、对应 `.bin` 与Nu/Re历史。检查点写入完整D2Q9 f、D2Q5 g，稳态另存上次检查参考场及误差；不需要额外的通量时间历史。标识 `FEI_CLBM_2D_V1`，metadata版本1003，检查网格、物理参数、松弛率、采样轴和算例，拒绝混用其他模型。可以延长时长或调整收敛阈值；快照不能代替检查点。

## 编译与核验

NVIDIA/NVHPC编译示例：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccFeiCLBM.F90 -o fei_clbm
```

开发核验使用支持OpenACC的GNU Fortran host后端：

```text
python3 verify_local.py --compiler /path/to/gfortran
python3 verify_transport.py --compiler /path/to/gfortran
python3 verify_benchmark.py --compiler /path/to/gfortran
```

每项只修改隔离副本，JSON报告绑定生产源SHA-256。local核验包含直接投影全部中心矩的独立oracle、非零双向力、20／02耦合及非默认高阶率、参考密度力、两方向纯导热Nu=1、宏组合、完整短算和不对齐检查相位的严格续算。transport使用实际生产核，重现论文L=100、kappa=0.05、Ma=0.3移动温度波至t*约2，并测量流场黏度。

2026-10-29上述三个JSON均通过：

- 公式／力矩、两向导热Nu=1、4种宏编译、侧热和RB完整短算通过；48步连续与32+16步续算的检查点及共同快照逐字一致，8种错误历史／模型组合被拒绝。
- 目标nu=kappa=0.05，测得nu=0.05000806，误差0.0161%；静止温度波kappa=0.05000428，误差0.00855%；原文Ma=0.3移动温度波kappa=0.05000691，误差0.0138%。移动波t*=1.99998的完整波形相对误差约0.772%，最大相位误差0.00769弧度；不能将扩散率接近视为波形完全相同。
- 32²、Ra=1000、Pr=0.71、Ma=0.05侧热在16000步收敛，体积Nu=1.11686967，冷热壁Nu均约1.11702312。对照Wang/Luo2013的1.1178，体积偏差约0.0832%；这里仅验证一个低Ra网格。差异更符合有限网格和壁面离散误差，未发现后处理或并行装配缺陷的证据。

尚未运行NVHPC／NVIDIA GPU，没有完成高Ra稳定窗口、长时RB统计、网格无关性或DNS分辨率验证；不根据原文高Ra快照宣称当前程序已达到该精度或优于已有方法。本轮没有登录服务器，没有暂存、提交或推送Git。
