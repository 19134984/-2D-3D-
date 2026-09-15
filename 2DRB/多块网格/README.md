# 2DRBOpenaccMultiblock：连通细网格环

当前版本（2026-09-15）为 **中心粗矩形＋外围一个连通细环**，不是原来的五个独立矩形。修改基于用户当前版本，原 `均匀网格/2DRBOpenacc.F90` 未修改。

先读 [代码逻辑说明](代码逻辑说明.md)，查看 [布局](图示/multiblock-layout.png)、[粗细节点](图示/multiblock-interface.png)、[连通迁移](图示/multiblock-samelevel.png) 和 [时间顺序](图示/multiblock-timestep.png)。

## 当前实现

- `nBlocks=2`：1 是粗网格，2 是一套连通细网格数组。只建立粗→细和细→粗两条连接，原上/下/左/右细区之间不再打包或交换。
- 细数组为 `nx×ny`，中心空区由 `fine_active` 排除在碰撞、迁移、宏观更新之外，保持直接二维邻点寻址。
- 保留块内 D2Q9 MRT、D2Q5 MRT、`EnableUseG`/旧温度分支、真实壁面以及流热推进顺序；只在细空区增加跳过逻辑。
- 保留 Huang 含力矩重标定、四点空间 Lagrange、三层时间插值，以及 Chen 碰撞前缓冲交换。
- 参数仍在 `module commondata`；默认四边节点编号 Left/Right/Bottom/Top 为 `128/129/128/129`，基准分界为 `127.5/895.5`，粗细比 2，单侧延伸 4 个细格距。
- 面积和中线统计显式扣除细数组中的中心物理区域，零面积处不参与极值或积分。

源码仍为 `module commondata → program main → 外部子程序`；普通数组保存区域属性，连续存储配合 `select_block` 指针选择区域。没有引入派生类型或 MPI。

## 存储取舍与输出

细数组保留中心空区的占位存储，因此不是内存最省的实现。默认两级合计分配 1,199,897 个节点位置，其中活动节点 623,816 个。以后若采用紧凑环形存储，需要另做索引和内核改造。

检查点为 **v9**，旧五块检查点不兼容，需重新初始化。快照为 **v3**，在每区一维权重后新增二维实际面积，后面才是 `u,v,T,rho`。读取器必须适配；完整顺序见 [逻辑说明](代码逻辑说明.md)。

Tecplot 输出中的 `integration_area/L^2=0` 表示不属于该区域物理积分范围，绘图必须屏蔽这些节点，不能显示中心占位数据。输出文件继续使用 `Multiblock` 前缀。

原均匀网格的流函数/涡量、壁面极值拟合等结束后处理尚未全部移植；这次保留现有功能范围。五块历史记录保留，但不作为本版验证依据。

## 本地检查

在本目录运行：

```text
python verify_multiblock.py
python transient_compare.py
python 图示/draw_layout.py
```

第一条在当前两区域版本自动调用 `verify_ring.py`。构建及运行输出写入系统临时目录，验证参数只在临时源码中修改。

- `ring_verification_results.json`：五种宏配置语法、2/4/8 粗细比几何、原接缝直接迁移、空区排除、线性导热、单块退化和精确重启。
- `ring_transient_comparison.json`：96²、Ra=1e4 侧壁差温 2000 细步瞬态，与原均匀代码比较，并比较不同粗细重叠宽度。
- `图示/layout_parameters.json`：当前源码哈希、几何与活动节点计数。

以上为本机 gfortran OpenACC host 验证，不代表 nvfortran/P100 GPU、高 Ra 长期稳定性、严格守恒或网格/时间收敛阶验证。
