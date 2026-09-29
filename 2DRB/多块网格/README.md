# 2DRBOpenaccMultiblock：中心粗网格与四套紧凑细数组

当前版本（2026-09-29）直接使用 `_coarse、_left、_right、_bottom、_top` 五套普通数组，不再通过块编号选择数据。本次简化参数与积分表达，保留数组布局和接口交换；基准为 `均匀网格/2DRBOpenacc.F90` 与 `ISLBM/2DRBOpenaccISLBM.F90`，两份基准源码未修改。变量和计算次序见 [代码逻辑说明](代码逻辑说明.md)。

## 参数

物理参数、网格数、运行及输出开关继续沿用基准命名；多块网格只需额外设置粗细比 `refineRatio` 和四个分区位置 `fineLayerCellsLeft/Right/Bottom/Top`。交换层宽度只保留内部常量 `interfaceSkin=2`：粗、细网格各延伸两个本地格距；换成细格距单位，粗网格外延为 `interfaceSkin*refineRatio`，细网格内延为 `interfaceSkin`。

细网格直接复用 `Snu、Sq、Qk、Qnu`，粗网格保留按格距换算后的松弛率。前一轮已删除 34 个重复的全局量：4 个细网格松弛率、16 个细区积分起止下标、4 个细区积分外框及 10 个一维积分权重数组。

全文件检查以保留基准参数约定为前提，主要清理多块新增的重复量和临时别名。`gBeta1/gBeta`、`velocityUnit`、`velocityScaleCompare`、`thermalGeqCoeff` 和 `unsteadyAverageStartTf/EndTf/MidTf` 保留原名及公式；多块子程序也使用 `Snu、Sq、Qk、Qnu、gBeta` 传参。初始速度为零，当前温度初始化直接写成 `omegaT*T`；快照仍保存格子速度，因此不为使用参数而改变输出单位。

主程序沿用基准的 `string`、`numAccDevices`、`unsteadyItcRemaining`，以及稳态在循环入口判断收敛、非稳态按剩余步数推进的组织。续算参数集中在 `loadInitField/reloadFileNum` 后。多块子程序统一使用 `integer(kind=4)`、`real(kind=8)` 和 `character(len=...)` 的声明风格；基准已有变量不因使用次数少而删除。

多块新增的 `gBetaCoarse`、`coarseOverlapCells`、`fineOverlapCells` 不再另存，使用处直接写 `dxCoarse*gBeta` 或 `interfaceSkin`。其余清理包括换算宽度、循环边界、中线位置和输出中转量；重叠宽度重复检查已删除，实际网格尺寸和插值模板检查保留。

格距、粗松弛率、几何尺寸、运行状态、累加器和多次使用的复杂表达式仍保留。两个接口插值子程序的 `timeWeight0/1/2` 是 NVHPC 标量捕获所需，保留其 `firstprivate` 写法。

## 数组与坐标

每套场数组的空间下标从 1 开始。例如左细网格：

```fortran
x = xOffsetLeft+(i-0.5d0)*dxFine
y = yOffsetLeft+(j-0.5d0)*dxFine
```

左右细区贯穿整个高度，上下细区只填左右之间的部分，角点不重复存储。默认参数下：

| 数组后缀 | 本地尺寸 | xOffset | yOffset | dx |
|---|---:|---:|---:|---:|
| _coarse | 389×389 | 122.5 | 122.5 | 2 |
| _left | 130×1024 | 0 | 0 | 1 |
| _right | 131×1024 | 893 | 0 | 1 |
| _bottom | 763×130 | 130 | 0 | 1 |
| _top | 763×131 | 130 | 893 | 1 |

细数组共 466,407 个节点，原完整细矩形为 1,048,576 个节点；节点存储减少约 55.5%。另有少量 f_post/g_post 迁移外圈。单个数组的一维仍可能为 1024，但没有完整的 1024×1024 细场数组，也不再分配中心空区。

内部交换层使粗块向分区线外延伸两个粗格距、细网格向中心延伸两个细格距；两侧各重建两层接收节点。默认左侧分区线为 127.5，粗接收节点为 123.5、125.5，细接收节点为 128.5、129.5，四边均采用这一规则。节点数与本次简化前一致；相比更早的对称延伸布局减少 6,088 个细节点（约 1.29%），未测量 GPU 加速比。

当前紧凑布局要求细存储的中心空区宽高为正；参数检查会拒绝重叠层填满中心空区的配置。`refineRatio=1` 仍只使用 _coarse 一套全域均匀数组。

## 计算与交换

- 保留原 D2Q9 流场、D2Q5 温度场、墙面处理及粗细接口尺度变换。
- 细区之间直接复制碰撞后的 f_post/g_post 外圈，含 D2Q9 对角方向；不做空间或时间插值。
- 每个细步先完成四区流场碰撞，交换 f_post，再迁移、处理墙面和恢复速度；随后完成四区温度碰撞，交换 g_post，再推进温度。
- 粗→细仍使用四点空间 Lagrange 和三时间层插值，初始粗步线性启动；细→粗使用同步时刻的共址数据。
- rhoHistory、uHistory、vHistory、THistory、FxHistory、FyHistory、flowNeqHistory、thermalNeqHistory 各有五套。粗历史末下标为 0:2，四套细历史为 0:0。
- 原 h 已统一命名为 dx。本程序最细格子单位下 dt=dx；力历史仍为 Fx/dx、Fy/dx，接收端再乘接收网格的 dx。
- 积分直接累加中心粗矩形与左、右、下、上四个细物理矩形，不再采用细存储区积分后扣除中心面积的表达。重叠存储只用于推进和交换，不另计物理面积；分区线上的节点仍按控制面积取半权，角点取相应乘积权重。中线模板跨细数组接缝时按全局坐标取值。

源码仍为 module commondata → program main → 外部子程序。没有 target、pointer、派生类型、contains、nBlocks 或块编号选择器。通用计算子程序显式接收数组、尺寸和所需几何参数。

## 续算和输出

计步沿用基准约定：`itc` 是本次运行已完成的细步数，每次启动为 0；`restartItcOffset` 是检查点中已有的累计细步数，新算例为 0。总步数为 `restartItcOffset+itc`，用于采样、输出、物理时间和时间插值启动判断。MLUPS 只计本次 `itc` 对应的工作量，无需额外的 `itcStart`。

`reloadDimensionlessTime=restartItcOffset/timeUnit` 记录续算前的无量纲时间；`reloadMetadataLoaded` 初始为假，检查点元数据和场量完整读取并核对成功后置真。两者均写入设置日志。当前 `.meta` 仅保存检查点文件名，真实时间元数据来自 `.bin`；手动指定检查点读取成功时，该标记同样为真。输出时间统一用累计步数除以 `timeUnit`，保持连续运行和续算的时间计算一致。

稳态 `itc_max` 限制本次运行步数；非稳态 `itc_max` 对应整个算例的总目标，续算只补剩余步数。两者都向上对齐粗细同步步。检查点和快照中的步数字段仍写累计值；读取检查点时将其放入 `restartItcOffset`，并把 `itc` 清零，文件布局不变。

本次简化不改变检查点与快照格式。检查点直接从网格参数开始，不写入 16 字节的文本格式标识；按 coarse、left、right、bottom、top 固定顺序保存数组及历史，保留物理配置检查、输出时钟和稳态检查场。

带有旧文本格式标识的文件不能直接读取。更早的对称布局会被几何检查拒绝，需使用匹配旧布局的源码续算，不能只修改当前的 `interfaceSkin`。下面的 v10 转换输出当前不带文本标识的检查点，但保留旧布局。

旧连通细环 **v10** 文件可转换；输入不会修改，输出路径必须不存在：

```text
python convert_restart_v10.py old-v10.bin converted-current.bin
```

转换后保留与检查点匹配的 NuRe/收敛历史，将 `reloadFile2DOpenaccMultiblock-latest.meta` 内容设为转换后的文件名，再设置 `loadInitField=1`。转换不改变计算时刻、物理参数、输出计数或历史时间层。v9 及更早格式不支持该工具。

快照直接从区域数量开始，不写入 16 字节的文本格式标识。区域数量为 1 或 5；每区含几何、一维权重、显式二维积分面积及 u/v/T/rho。一维权重仅在输出时计算，无需常驻数组。读取器应按文件记录的区域数量循环。Tecplot 区域采用具名标题，零积分面积节点需在物理区域绘图时屏蔽。

稳态误差检查、非稳态 Nu/Re 窗口统计、独立输出间隔和开关均保留。原均匀代码中此前尚未移植的流函数/涡量等后处理不属于此次新增功能。

## 验证

最新风格对齐记录见 `baseline_style_verification.json`：75 个子程序除类型声明写法外保持一致；检查实际主程序的稳态多次续算、本次计数、连续与续算对照、非稳态剩余时长、输出开关和最终输出。此前数值与计步检查按下面的报告哈希追溯。

计步约定验证见 `restart_counter_verification.json`：语法检查、7 组短程前后对照、旧检查点读取、实际稳态/非稳态主程序连续与续算对比均通过；多次续算的本次步数与偏移量、非稳态已达目标时零步退出、累计输出时刻和本次 MLUPS 计数均检查。场量、历史和诊断对照未发现差异。以下其他报告为此前版本，按源码哈希区分。

```text
python verify_simplified.py path/to/pre-simplification.F90
python verify_simplified.py path/to/pre-inline.F90 --report inline_verification.json
python verify_simplified.py path/to/pre-conventions.F90 --report baseline_conventions_verification.json
python verify_asymmetric.py
python verify_asymmetric.py path/to/pre-change-compact.F90
```

`simplified_verification.json` 保留前一轮参数与积分简化的结果：7 个算例各运行 40 个细步后通过检查，场量及交换历史最大绝对差约 `9.02e-16`，Nu/Re 最大差约 `2.23e-16`，稳态误差记录一致；快照几何和逐节点积分面积精确相同；新源码读取旧检查点后的场量及交换历史逐位相同。该报告以记录的源码哈希为准。

上一轮保留基准参数约定的结果见 `baseline_conventions_verification.json`：7 组各 40 个细步的场量、交换历史、Nu/Re、稳态误差及快照字段前后差值均为 0，几何与逐节点面积精确一致，旧检查点读取逐位一致。解析时间统计也通过：两个历史时刻的 Re 为 1、3，整窗 RMS 为 `sqrt(5)`，两半窗为 `sqrt(3)`、`sqrt(7)`；其他 Nu 均值、窗口覆盖及相对半窗差符合解析值。旧对称布局测试副本的场量逐位一致、旧几何检查点拒绝和绘图参数读取也已检查。

`inline_verification.json` 保存此前表达式内联及实际主程序连续/续算、输出开关、最终输出和数组释放的检查结果，按各报告的源码哈希区分版本。

`verify_asymmetric.py` 检查粗细比 1/2/4/8 的小网格运行、全部细接收点的三次多项式插值（含角点）、快照面积及有限值、稳态/非稳态精确续算和独立输出开关。提供旧紧凑数组源码时，还记录相同算例的 Nu/Re 变化、检查恢复旧宽度后场量逐位一致，并验证旧几何检查点被拒绝。构建和运行输出位于系统临时目录。

`verify_asymmetric.py` 的报告为 `asymmetric_verification.json`，包含源码哈希和实际检查内容；应按报告哈希区分验证版本。`compact_verification.json` 及旧连通细环逐位对照属于此前对称布局的验证记录；非对称布局改变了接口位置，不要求它与旧布局场量逐位相同。历史报告保留。

这些是 **gfortran OpenACC host** 短程检查，不是当前源码的 nvfortran/P100 GPU 验证，也不证明长期统计收敛。此前的 P100 报告、ring/plain/named_history 报告和 transient_compare.py 对应历史版本，不能作为本次参数与积分简化的验证结果。

## 当前图示

[图示目录](图示/README.md) 中的总体布局、接口放大、同级交换和时间推进图已按当前五套数组与非对称重叠更新。运行 `python 图示/draw_layout.py` 可从源码参数重新生成 PNG、SVG 和带源码哈希的几何记录。图示是代码说明，不代表数值验证。
