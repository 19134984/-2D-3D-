# 2DRBOpenaccMultiblock：中心粗网格与四套紧凑细数组

当前版本（2026-09-18）直接使用 `_coarse、_left、_right、_bottom、_top` 五套普通数组，不再通过块编号选择数据。原均匀网格和 ISLBM 源码未修改。变量和计算次序见 [代码逻辑说明](代码逻辑说明.md)。

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
| _left | 132×1024 | 0 | 0 | 1 |
| _right | 133×1024 | 891 | 0 | 1 |
| _bottom | 759×132 | 132 | 0 | 1 |
| _top | 759×133 | 132 | 891 | 1 |

细数组共 472,495 个节点，原完整细矩形为 1,048,576 个节点；节点存储减少约 54.9%。另有少量 f_post/g_post 迁移外圈。单个数组的一维仍可能为 1024，但没有完整的 1024×1024 细场数组，也不再分配中心空区。

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

新检查点为 **v11**：按 coarse、left、right、bottom、top 固定顺序保存数组及历史，保留物理配置检查、输出时钟和稳态检查场。

旧连通细环 **v10** 文件可转换；输入不会修改，输出路径必须不存在：

```text
python convert_restart_v10.py old-v10.bin converted-v11.bin
```

转换后保留与检查点匹配的 NuRe/收敛历史，将 `reloadFile2DOpenaccMultiblock-latest.meta` 内容设为转换后的文件名，再设置 `loadInitField=1`。转换不改变计算时刻、物理参数、输出计数或历史时间层。v9 及更早格式不支持该工具。

快照仍为 **v3**，区域数量为 1 或 5；每区含几何、一维权重、显式二维积分面积及 u/v/T/rho。读取器应按文件记录的区域数量循环。Tecplot 区域采用具名标题，零积分面积节点需在物理区域绘图时屏蔽。

稳态误差检查、非稳态 Nu/Re 窗口统计、独立输出间隔和开关均保留。原均匀代码中此前尚未移植的流函数/涡量等后处理不属于此次新增功能。

## 验证

```text
python verify_multiblock.py
python verify_compact.py path/to/pre-change-ring.F90
```

第一条自动运行紧凑数组版本的本地检查；第二条另外与指定的改动前连通细环源码逐节点对照。小网格参数只修改临时源码，构建和运行输出位于系统临时目录。

`compact_verification.json` 记录源码哈希和实际检查内容：粗细比 1/2/4/8、流场/温度/接口历史逐位对照、稳态和非稳态实际主程序续算、v10 转换、中线模板跨接缝、快照面积以及独立输出开关。

这些是 **gfortran OpenACC host** 检查，不是当前源码的 nvfortran/P100 GPU 验证。此前的 P100 报告、ring/plain/named_history 报告、transient_compare.py 和图示目录中的两区域图对应历史版本，不能作为本次五套数组布局的验证结果。
