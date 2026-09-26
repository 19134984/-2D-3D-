# 2DRBOpenaccMultiblock：中心粗网格与四套紧凑细数组

当前版本（2026-09-20）直接使用 `_coarse、_left、_right、_bottom、_top` 五套普通数组，不再通过块编号选择数据。原均匀网格和 ISLBM 源码未修改。变量和计算次序见 [代码逻辑说明](代码逻辑说明.md)。

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

`coarseOverlapCells=2` 控制粗块向分区线外延伸两个粗格距；`fineOverlapCells=2` 控制细网格向中心延伸两个细格距。`interfaceSkin=2` 仍表示各自重建的两层节点。默认左侧分区线为 127.5，粗接收节点为 123.5、125.5，细接收节点为 128.5、129.5，四边均采用这一规则。相比此前对称延伸的 472,495 个细节点，减少 6,088 个（约 1.29%）；未测量 GPU 加速比。

当前紧凑布局要求细存储的中心空区宽高为正；参数检查会拒绝重叠层填满中心空区的配置。`refineRatio=1` 仍只使用 _coarse 一套全域均匀数组。

## 计算与交换

- 保留原 D2Q9 流场、D2Q5 温度场、墙面处理及粗细接口尺度变换。
- 细区之间直接复制碰撞后的 f_post/g_post 外圈，含 D2Q9 对角方向；不做空间或时间插值。
- 每个细步先完成四区流场碰撞，交换 f_post，再迁移、处理墙面和恢复速度；随后完成四区温度碰撞，交换 g_post，再推进温度。
- 粗→细仍使用四点空间 Lagrange 和三时间层插值，初始粗步线性启动；细→粗使用同步时刻的共址数据。
- rhoHistory、uHistory、vHistory、THistory、FxHistory、FyHistory、flowNeqHistory、thermalNeqHistory 各有五套。粗历史末下标为 0:2，四套细历史为 0:0。
- 原 h 已统一命名为 dx。本程序最细格子单位下 dt=dx；力历史仍为 Fx/dx、Fy/dx，接收端再乘接收网格的 dx。
- 积分使用原物理分区，细区扣除中心面积。中线模板跨细数组接缝时按全局坐标取值。

源码仍为 module commondata → program main → 外部子程序。没有 target、pointer、派生类型、contains、nBlocks 或块编号选择器。通用计算子程序显式接收数组、尺寸和所需几何参数。

## 续算和输出

检查点直接从网格参数开始，不再写入 16 字节的文本格式标识；按 coarse、left、right、bottom、top 固定顺序保存数组及历史，保留物理配置检查、输出时钟和稳态检查场。

带有旧文本格式标识的文件不能直接读取。缩小细数组后，旧对称布局也会被几何检查拒绝；沿用旧布局时需要设置 `fineOverlapCells=coarseOverlapCells*refineRatio`，保持数组尺寸和坐标一致。下面的 v10 转换输出当前不带文本标识的检查点，但保留旧布局，续算同样需要恢复旧宽度。

旧连通细环 **v10** 文件可转换；输入不会修改，输出路径必须不存在：

```text
python convert_restart_v10.py old-v10.bin converted-current.bin
```

转换后保留与检查点匹配的 NuRe/收敛历史，将 `reloadFile2DOpenaccMultiblock-latest.meta` 内容设为转换后的文件名，再设置 `loadInitField=1`。转换不改变计算时刻、物理参数、输出计数或历史时间层。v9 及更早格式不支持该工具。

快照直接从区域数量开始，不再写入 16 字节的文本格式标识。区域数量为 1 或 5；每区含几何、一维权重、显式二维积分面积及 u/v/T/rho。读取器应按文件记录的区域数量循环。Tecplot 区域采用具名标题，零积分面积节点需在物理区域绘图时屏蔽。

稳态误差检查、非稳态 Nu/Re 窗口统计、独立输出间隔和开关均保留。原均匀代码中此前尚未移植的流函数/涡量等后处理不属于此次新增功能。

## 验证

```text
python verify_asymmetric.py
python verify_asymmetric.py path/to/pre-change-compact.F90
```

检查粗细比 1/2/4/8 的小网格运行、全部细接收点的三次多项式插值（含角点）、快照面积及有限值、稳态/非稳态精确续算和独立输出开关。提供旧紧凑数组源码时，还记录相同算例的 Nu/Re 变化、检查恢复旧宽度后场量逐位一致，并验证旧几何检查点被拒绝。构建和运行输出位于系统临时目录。

当前报告为 `asymmetric_verification.json`，包含源码哈希和实际检查内容。`compact_verification.json` 及旧连通细环逐位对照属于此前对称布局的验证记录；非对称布局改变了接口位置，不要求它与旧布局场量逐位相同。

这些是 **gfortran OpenACC host** 检查，不是当前源码的 nvfortran/P100 GPU 验证。此前的 P100 报告、ring/plain/named_history 报告和 transient_compare.py 对应历史版本，不能作为本次五套数组布局的验证结果。

## 当前图示

[图示目录](图示/README.md) 中的总体布局、接口放大、同级交换和时间推进图已按当前五套数组与非对称重叠更新。运行 `python 图示/draw_layout.py` 可从源码参数重新生成 PNG、SVG 和带源码哈希的几何记录。图示是代码说明，不代表数值验证。
