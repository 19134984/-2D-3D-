# 当前多块网格图示

2026-09-20：中心粗数组与左、右、下、上四套细数组，采用非对称重叠。

- [总体布局](multiblock-layout.png)：五套数组的尺寸、积分分区和细数组接缝。
- [接口放大](multiblock-interface.png)：粗边缘 123.5、分区线 127.5、细边缘 129.5，以及两侧接收层。
- [同级细网格交换](multiblock-samelevel.png)：四套细数组之间先交换碰撞后外圈，再迁移。
- [时间推进](multiblock-timestep.png)：粗步预测、两个细子步、双向同步及时间插值。

每张图都有同名 SVG 矢量版本。粉色虚线表示积分分区线；人工边缘是各自最外计算节点，不是物理墙壁。接口放大图的粉框为粗接收节点，绿点为细接收节点。

生成脚本读取上级目录的 `2DRBOpenaccMultiblock.F90`，不修改求解器：

```powershell
python .\draw_layout.py
```

需要 Pillow 和 Windows 微软雅黑字体。当前时间推进图针对 `refineRatio=2`，其他比值会明确报错以避免错误时序。`layout_parameters.json` 记录源码 SHA-256、五区节点尺寸和坐标、两侧延伸量、积分面积；生成时检查总积分面积等于 `nx*ny`。
