# 二维均匀网格 FPC 热对流

本目录以单 GPU OpenACC 复现 Schupbach、Premnath（2025）《Fokker-Planck Model Based Central Moment Lattice Boltzmann Method for Effective Simulations of Thermal Convective Flows》的流场 D2Q9 和温度 D2Q9 option A。主求解循环在单一 OpenACC 设备执行；文件输出及部分最终诊断在同步后由主机完成。当前完成代码与本地检查，论文基准、长时间高瑞利数稳定性及网格无关性留待服务器测试。

## 框架与算例

- 按仓库要求核对 `../均匀网格/2DRBOpenmp.F90` 的二维框架、变量及诊断组织；实际执行框架、OpenACC 数据常驻、并行循环和文件输出沿用同维度 `../均匀网格/2DRBOpenacc.F90`。本实现只有单 GPU OpenACC 执行路径。
- 默认侧壁差温稳态：左热右冷，上下绝热，四壁无滑移。
- RB：下热上冷，左右绝热，四壁无滑移；初始线性温度叠加可直接设置幅度的确定性扰动。
- 均为单块均匀网格。代码保留原项目的 `commondata`、子程序顺序、中文说明、Nu/Re、壁面 Nu、中心线极值、流函数/涡量、快照和 Tecplot 输出风格。

在主文件顶部选且只选一个流动模式和一个算例。侧热使用 `steadyFlow` + `SideHeatedCell`；RB 使用 `unsteadyFlow` + `RayleighBenardCell`。边界随算例自动配套。所有参数直接在 `commondata` 修改，默认 `256×256`、`Ra=10^6`、`Pr=0.71`、`Uff=0.01`、冷热壁温度 `2/1`。这些是可启动的默认参数，不代表该网格已经通过精度验证。

## 碰撞与边界处理

流场采用中心矩 FPC 碰撞：二阶迹/偏差分别松弛；四阶吸引子为碰撞后二阶矩的组合 `(k20*k02+2*k11^2)/rho`。温度采用论文 option A，三阶、四阶吸引子依赖碰撞后的低阶温度矩。没有另加滤波或正则化。

两套碰撞读取同一时刻的宏观场。迁移及边界后先恢复温度，再以新温度恢复浮力和带半步力的物理速度。初始流场分布包含半步力修正，保证初始物理速度为零。

恒温采用半格距反反弹；D2Q9 绝热边界采用仅反射法向、保留切向上游位置的镜面反射，角点恒温条件优先。这一处理已通过两方向纯导热检查。

文中式(37)的 `sigma11` 下标及 `sigma12` 系数存在排印不一致，本实现按通式(29)构造高阶中心力源。连续一阶中心矩取零，不将存储分布的 `-F/2` 代回二阶源。温度矩按中心矩二项式定义实现，以避免附录 B 中个别符号混用。

## 输出和续算

输出前缀统一含 `2DOpenaccFPC`。非稳态 Nu/Re 记录绝对自由落体时间；体平均 Nu 使用温差归一化，Re 使用全域均方根速度。结束时根据完整历史重建运行均值和指定统计窗口的前、后半段平均。窗口比较仅提供统计诊断，不自动证明统计收敛。

`loadInitField=1` 读取 `reloadFile2DOpenaccFPC-latest.meta` 指向的检查点。请配套保存 `.bin`、`.meta`、Nu/Re 历史；程序检查算法标识、网格、两套 Q 数、工况、物理参数、采样间隔和历史计数，不猜测缺失元数据，也不接收旧 D2Q5 检查点。采样始终锚定原算例时间零点，非采样时刻的重启不改变采样相位。稳态续算保留旧 Nu/Re 行并追加新结果。

`itc_max` / `unsteadyRunDuration` 表示整个算例的总目标，可在续算时提高；目标已达到时不重复推进。二进制使用 Fortran sequential unformatted，跨编译器或不同记录标记格式的兼容性未验证。

## 本地检查

`verify_fpc_local.py` 编译实际主源的临时小网格副本，使用 GNU Fortran 的 OpenACC host 后端，并开启边界检查及浮点异常捕获。它检查独立矩恒等式、碰撞守恒与松弛、力增量、初值、600 步纯导热、48 步侧热稳态/RB 非稳态，以及 32+16 步续算；采样间隔为 7 步，特意覆盖非采样点续算。报告与所测主源 SHA-256 写入 `validation_local.json`。

可在 ASCII 临时目录运行，例如：

```text
python verify_fpc_local.py --compiler C:/msys64/ucrt64/bin/gfortran.exe --temp-root C:/Temp
```

完整主程序的本地编译形式为：

```text
gfortran -cpp -fopenacc -ffree-line-length-none -O2 2DRBOpenaccFPC.F90 -o fpc.exe
```

后续单 GPU 正式运行使用 NVHPC，编译形式如下；本机没有 `nvfortran`，尚未实际核验这一步：

```text
nvfortran -O3 -acc -Mpreprocess -Mextend 2DRBOpenaccFPC.F90 -o fpc
```

本地检查不能替代 NVHPC 编译及 NVIDIA GPU 执行，也未验证论文 Nu/Re 基准、长时间高 Ra 稳定性、统计收敛或网格无关性。当前未向服务器提交任务。

## 原文件保护记录

原文件 SHA-256 与开始实施前一致：

| 原文件 | SHA-256 |
| --- | --- |
| `均匀网格/2DRBOpenmp.F90` | `83A17F30CD4715BD14E9FA2D06AA43F026CA4FF57A11154CBCDBBBCF6B7A713F` |
| `均匀网格/2DRBOpenacc.F90` | `68DEBF94F75F3A70AF6BF4F7205E18AFC14618D48DD227CCF1D1A8E73D8AFD48` |
| `后处理/postprocess_statAvg_2DOpenacc.F90` | `89FC4708CBC489DDF528175B6C35EA35D5442F08655D9F7C88CD9A5453A1F310` |

论文原件位于 `../pdf/算法/William Schupbach-Kannan Premnath-2025-Fokker-Planck Model Based Central Moment Lattice Boltzmann Method for Effective Simulations of Thermal Convective Flows.pdf`。
