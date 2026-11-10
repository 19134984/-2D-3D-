# Palabos 四叉树与八叉树加密及通信分析

调查日期：2026-11-09。证据版本：本地 `Palabos` 仓库提交 `8f8ecd277204f9d24019834068906128565c9346`。本报告根据当前源码作静态追踪，未编译运行算例，未修改 Palabos 源码。文中的源码路径均相对于本报告所在的 `文档/`。

## 一、先区分二维矩形局部加密与三维八叉树

在当前仓库全文搜索 `quadtree` / `quad.tree` 没有命中；`src` 中名称带 tree 的网格实现集中在 `gridRefinement/octree*`。这只是搜索结果，不能单独作为“不存在”的证明；更有力的实现证据是：二维公开管理类使用 `std::vector<std::vector<Box2D>> bulks` 保存各层矩形块，另存粗细接口、方向和 MPI 进程号，没有 parent / child[4] 结构。二维示例直接调用 `management.refine(level, Box2D)` 添加矩形区域。因此，对当前版本应称其为“二维二倍尺度的层级矩形块加密”，不宜直接称为“二维四叉树实现”。

证据：[multiGridManagement2D.h](../Palabos/src/multiGrid/multiGridManagement2D.h) 第 48–76、111–134 行；[dipole.cpp](../Palabos/examples/showCases/gridRefinement2d/dipole.cpp) 第 301–336 行。

二维主要入口在 `src/multiGrid/`；三维显式八叉树另在 `src/gridRefinement/`。这两条实现路径必须分开追踪，尤其不能把三维的过滤和插值规则直接用于解释二维。

## 二、二维加密区域如何生成

`MultiGridManagement2D::initialize` 从最粗层包围盒逐级放大坐标，只有 `referenceLevel` 初始具有整域 bulk。随后 `refine(coarseLevel, coarseDomain)` 做以下处理：

1. `fineLevel = coarseLevel + 1`；输入区域使用 coarseLevel 的坐标。
2. 将区域裁剪到该层外边界，并记录是否接触四条外边界。
3. `reducedCoarseDomain = coarseDomain.enlarge(-1)`：从原粗网格删除内部，留下接口所需的粗节点。
4. `extendedCoarseDomain = coarseDomain.enlarge(overlapWidth)`：向外扩展重叠区，再 `multiply(2)` 生成细层矩形块。
5. 在没有接触物理外边界的四边分别建立粗侧、细侧耦合接口，并保存方向与朝向；接触物理外边界时调整扩展区且不生成相应的内部粗细接口。

证据：[multiGridManagement2D.cpp](../Palabos/src/multiGrid/multiGridManagement2D.cpp) 第 96–109、154–258 行。逆向的几何操作由同文件第 265 行开始的 `coarsen` 提供；这不等于示例在运行时自动依据误差估计动态粗化。

例如不接触外边界时，粗层区域取 `[10,20]×[10,20]`，`overlapWidth=1`：粗网格删除 `[11,19]×[11,19]`，新增细矩形为 `[18,42]×[18,42]`。粗侧左接口位于粗坐标 `x=10`，细侧左边界在粗单位中为 `x=9`，实际细坐标为 `x=18`。粗细网格之间保留重叠带，因此它不是把所有父格点简单删除后只留下四个子格点。

当前二维示例 `gridRefinement2d/dipole.cpp` 第 309–316 行手工指定两块区域，再分别加密 level 0 与 level 1；第 335–336 行选择三次空间插值和细回粗过滤。没有在该示例主循环中发现基于流场误差判据重新构造加密区域的逻辑。

## 三、二维交换的是什么量

`MultiGridLattice2D::createInterfaces` 将每层的粗侧接口与下一细层、细侧接口与上一粗层配对，并安装 data processor。二维默认 generator 使用 `scalingOrder=0`，源码注释称其为分解成 rho、u、fneq；进一步查看实际等温 / BGK dynamics，向量实际保存的是 `rhoBar`、`j`、各方向 `fneq`，随后还有 descriptor 的 external scalars。因此描述当前 BGK 实现时应精确写成 `rhoBar / j / fneq`。

证据：[multiGridLattice2D.hh](../Palabos/src/multiGrid/multiGridLattice2D.hh) 第 236–269 行；[multiGridGenerator2D.hh](../Palabos/src/multiGrid/multiGridGenerator2D.hh) 第 238–240、334–340 行；[isoThermalDynamics.hh](../Palabos/src/basicDynamics/isoThermalDynamics.hh) 第 124–135、366–410 行。

它不是直接对粗细层原始 `f_i` 作完全相同的复制，而是先分解、做尺度变换和插值，再重构 `f_i = f_i^eq + f_i^neq`。

### 3.1 粗 → 细：尺度变换、空间插值、两个时间槽

`ConvectiveRescaleEngine::scaleCoarseFine` 首先调用粗格点 dynamics 的 `decompose`。对非 advection-diffusion dynamics，计算

\[
\tau_f = 2\tau_c-\frac12,
\qquad \tau_c=1/\omega_c,
\]

先将非平衡分量乘 `tauFine/tauCoarse`，然后调用 dynamics 的 `rescale(xDxInv=2, xDt=1/2)`。对当前等温 `order=0` 实现，密度不变，`j` 的倍率为 `xDt*xDxInv=1`，`fneq` 还要再乘 `xDt=1/2`，所以总倍率为

\[
f^{neq}_f = \frac{\tau_f}{2\tau_c}\, f^{neq}_c.
\]

这里的“总倍率”针对上述等温 dynamics；不能无条件推广给所有 dynamics。源码的 `isAdvectionDiffusion()` 分支会跳过这段 tau/rescale 处理；external field 的尺度由具体 dynamics 负责。

证据：[gridRefinement.hh](../Palabos/src/multiGrid/gridRefinement.hh) 第 89–126 行；[isoThermalDynamics.hh](../Palabos/src/basicDynamics/isoThermalDynamics.hh) 第 260–279 行。

粗细坐标重合的细格点接收尺度变换后的粗向量；位于两粗点之间的细格点进行插值。二维有两种可选算法：

- 线性：`(g0+g1)/2`。
- 中心三次：`9/16*(g0+g1)-1/16*(g[-1]+g2)`。

这些权重逐分量用于分解后的向量 `g`。接口端部改用非中心二次插值 `3/8*g0+3/4*g1-1/8*g2`，并有角点附近 helper；不要把接口所有位置都理解成同一个四点模板。

证据：[gridRefinement.hh](../Palabos/src/multiGrid/gridRefinement.hh) 第 217–255 行；[fineGridProcessors2D.hh](../Palabos/src/multiGrid/fineGridProcessors2D.hh) 第 89–172、327–424、736–775、829–863 行；[multiGridGenerator2D.hh](../Palabos/src/multiGrid/multiGridGenerator2D.hh) 第 276–320 行。

粗到细处理器把结果写到细边界 dynamics 的 `getDecomposedValues(1)`，即 `T1` 槽；`FineGridBoundaryDynamics` 另保存 `T0`。细时间处于周期起点时直接用 T0 重构，否则按 `relativeTime/numTimeSteps` 在线性时间插值后重构。`numTimeSteps=2` 时，中间时间使用 `(T0+T1)/2`。

证据：[fineGridProcessors2D.hh](../Palabos/src/multiGrid/fineGridProcessors2D.hh) 第 106、145–154、347、397–403 行；[gridRefinementDynamics.hh](../Palabos/src/multiGrid/gridRefinementDynamics.hh) 第 132–163 行。

### 3.2 细 → 粗：重合点限制与可选非平衡过滤

`CopyFineToCoarse2D::process` 仅在 `fineTime % 2 == 1` 执行：从粗点坐标找到重合细点，调用 `scaleFineCoarse`，再用粗 dynamics 重构。对上述等温 dynamics，逆变换为

\[
\tau_c = \frac{\tau_f+1/2}{2},
\qquad f^{neq}_c = \frac{2\tau_c}{\tau_f}\,f^{neq}_f.
\]

证据：[coarseGridProcessors2D.hh](../Palabos/src/multiGrid/coarseGridProcessors2D.hh) 第 95–130 行；[gridRefinement.hh](../Palabos/src/multiGrid/gridRefinement.hh) 第 130–161 行。

有过滤版本 `CopyFineToCoarseWithFiltering2D`：先取重合细点，再遍历 `iPop=1..q-1` 的邻点，各点先尺度变换，再把 `start=1+d` 到向量末尾的分量作等权平均，除以 `q`。因此：

- 中心点已在循环前计入，过滤包含中心点及 q−1 个邻点；D2Q9 为 3×3 共九点。
- `rhoBar` 与 `j` 不参加这里的空间平均，仍来自重合细点。
- `fneq` 以及该向量尾部的 external scalars 参加这段平均。
- 这不是将九个完整 `f_i` 平均后原样覆盖粗点，也不是将四个子节点统一求平均。

证据：[coarseGridProcessors2D.hh](../Palabos/src/multiGrid/coarseGridProcessors2D.hh) 第 190–238 行。二维 generator 将细回粗处理器安装在 processor level 2，粗到细在 coarse processor level 1，细时间槽轮换在 fine processor level 0，见 [multiGridGenerator2D.hh](../Palabos/src/multiGrid/multiGridGenerator2D.hh) 第 272–320、325–360 行。

## 四、二维一粗步 / 两细步的实际调度

`iterateMultiGrid(level)` 的顺序是：

```text
当前层 collideAndStream 一次（包含该层已安装的接口处理器）
若有下一细层：
    iterateMultiGrid(level+1)
    iterateMultiGrid(level+1)
    当前粗层 duplicateOverlaps(staticVariables)
```

所以 level 0 走一步，level 1 走两步，level 2 走四步。粗层先推进得到下一粗时刻，再给细边界提供 T1，细层两步利用 T0 / T1 补齐边界。

证据：[multiGridLattice2D.hh](../Palabos/src/multiGrid/multiGridLattice2D.hh) 第 388–407 行。

时间计数的细节很重要：单层 `MultiBlockLattice2D::collideAndStream` 先碰撞迁移、再执行内部处理器、最后才 `incrementTime`。因此处理器看到 `fineTime%2==1` 时，正在完成的是本周期第二个细步。此时 level 0 的 `Copy_t1_to_t0_2D` 先把 T1 复制到 T0 并重新补齐边界分布，level 2 的细回粗随后执行；接着才增长计数。

证据：[multiBlockLattice2D.hh](../Palabos/src/multiBlock/multiBlockLattice2D.hh) 第 331–337 行；[multiBlock2D.cpp](../Palabos/src/multiBlock/multiBlock2D.cpp) 第 493–518 行；[fineGridProcessors2D.hh](../Palabos/src/multiGrid/fineGridProcessors2D.hh) 第 879–897 行。初始化时还会专门调用一次 Copy_t1_to_t0，避免最初的 T0 为空，见 [multiGridLattice2D.hh](../Palabos/src/multiGrid/multiGridLattice2D.hh) 第 204–226 行。

## 五、二维粗细交换和 MPI 通信的衔接

这里要区分两层操作：上文的分解、尺度变换、插值是粗细层之间的数值耦合；`duplicateOverlaps` 是块与块之间同步 envelope 的通信。

1. `generateLattices` 为每一加密层构造独立的 `MultiBlockLattice2D`，按矩形块登记 blockId 和 MPI rank，交给该层 communicator；见 [multiGridGenerator2D.hh](../Palabos/src/multiGrid/multiGridGenerator2D.hh) 第 47–77 行。
2. 示例采用 `ParallellizeBySquares2D`。它在最细层划分 MPI 区域，再逐层将同一批区域坐标除以二，沿用同一 MPI 分配，使相同物理位置的粗细块落在相应 rank；见 [parallelizer2D.cpp](../Palabos/src/multiGrid/parallelizer2D.cpp) 第 128–175 行和 [dipole.cpp](../Palabos/examples/showCases/gridRefinement2d/dipole.cpp) 第 321–328 行。
3. `MultiProcessing2D` 只枚举各 multiblock 的本 rank blocks，将其他层的几何缩放到第一层坐标后求交，再给本地 AtomicBlock 安装耦合处理器。因此这条二维耦合入口本身不包含每个粗细节点直接发送 MPI 消息的代码，也不能据此宣称任意粗细块跨 rank 布局都能自动耦合；见 [multiBlockOperations2D.hh](../Palabos/src/multiBlock/multiBlockOperations2D.hh) 第 105–183 行和 [multiBlockOperations2D.cpp](../Palabos/src/multiBlock/multiBlockOperations2D.cpp) 第 139–167 行。
4. data processor 注册时声明哪些块、哪些变量被修改，从而安排处理器执行后的 envelope 同步。粗到细声明细块为 `allVariables`，因为 T0/T1 位于 boundary dynamics；细到粗声明粗块为 `staticVariables`。见 [fineGridProcessors2D.h](../Palabos/src/multiGrid/fineGridProcessors2D.h) 第 77–81 行、[coarseGridProcessors2D.h](../Palabos/src/multiGrid/coarseGridProcessors2D.h) 第 69–77、106–114 行。
5. `executeInternalProcessors(level)` 本地运行处理器后调用 `duplicateOverlapsInModifiedMultiBlocks`，再交给 block communicator。细回粗只写 bulk，整个两细步完成后 multigrid 递归例程还明确对粗层做一次 staticVariables overlap 复制。见 [multiBlock2D.cpp](../Palabos/src/multiBlock/multiBlock2D.cpp) 第 509–518、593–634 行、[multiGridLattice2D.hh](../Palabos/src/multiGrid/multiGridLattice2D.hh) 第 395–398 行。

由此，二维实现可以概括为：**同一 MPI 区域内部完成相邻粗细层的尺度变换与时空插值，层内多块 communicator 负责跨 rank 的 envelope 同步，再进入下一阶段计算。** 对通信缓冲区和实际 MPI send/recv 的实现，还需继续阅读底层 parallel block communicator，不能只停在 `refine` 函数。

## 六、 三维八叉树：网格生成与邻居关系

### 6.1 真正的树结构和分裂

`OctreeNode<DataT>` 保存 `level/isLeaf/parent/child[8]`。一片叶子调用 `splitOctreeNode` 后创建 8 个孩子，`AssignOctreeData` 用 `getOctant` 分别赋予子长方体，再把父节点的 `isLeaf` 设为 false。它是块级八叉树：叶子对应一块有许多格点的规则网格，不是每个叶子只存一个 LBM 格点。

证据：[octree.h](../Palabos/src/gridRefinement/octree.h) 第 41–50 行；[octree.hh](../Palabos/src/gridRefinement/octree.hh) 第 97–109 行；[octreeGridGenerator.hh](../Palabos/src/gridRefinement/octreeGridGenerator.hh) 第 133–142、658–675 行。

`RefineOctree::operator()` 根据给定的网格密度函数决定目标层级：
\[
L_\mathrm{target}=\operatorname{round}\{L_\min+(L_\max-L_\min)\,g_\max\}.
\]
其中 `g_max` 是当前块内采样点或体积查询得到的最大密度值，乘 `gridDensityScaleFactor` 后以上限 1 截断；低于最小层级的块直接分裂，当前层级低于目标层级也继续分裂。例如 `Lmin=1,Lmax=4,gmax=0.8`，目标为 `round(3.4)=3`。这里的密度是网格分辨率需求函数，不是流体密度。

随后 `CalibrateOctree` 检查面、棱、角邻居；若细叶子的层级比粗邻居大 2 级或更多，就继续分裂粗邻居。生成器循环执行这个检查，因此相邻有效叶子满足 2:1 平衡，即层级差不超过 1。

证据：[octreeGridGenerator.hh](../Palabos/src/gridRefinement/octreeGridGenerator.hh) 第 358–428、448–462、1479–1504 行。

### 6.2 邻居如何找到、如何保存

`getOctreeGreaterEqualFaceNeighbor` 的步骤是：若目标面越过父节点边界，就递归到父节点找邻居；找到相邻祖先后，若它不是叶子，再通过 `OctreeTables::reflect` 找对应孩子。棱和角有类似函数。`gatherOctreeGreaterEqualNeighbors` 汇集 6 个面、12 条棱、8 个角共 26 个方向。

生成后的 `OctreeGridStructure` 保存 `blockId -> bulk/level/processId/isOverlap` 及 `blockId -> 26方向邻接表`：
- 同级或更粗的叶子：记录其 blockId；
- 更细的邻域：记录 `OctreeTables::smaller()` 标记；
- 域外：记录 `border()`；
- 未分配的叶子：记录 `undef()`。

所以邻接表并不为一个粗面的所有细块逐个列出 ID；后续 `CouplingInterfaces3D` 根据这些邻居类型和规则盒几何构造粗到细、细到粗的面、棱、角接口，再合并可合并的接口。

证据：[octree.hh](../Palabos/src/gridRefinement/octree.hh) 第 269–288、350–389 行；[octreeGridGenerator.hh](../Palabos/src/gridRefinement/octreeGridGenerator.hh) 第 1117–1154 行；[octreeGridStructure.h](../Palabos/src/gridRefinement/octreeGridStructure.h) 第 114–124 行；[couplingInterfaceGenerator3D.cpp](../Palabos/src/gridRefinement/couplingInterfaceGenerator3D.cpp) 第 56–100 行。

### 6.3 树如何变成求解网格

在粗细交界处，生成器还把部分非叶父块标记为 `isOverlap`，保留交界所需的粗层计算块。它与每个 MultiBlock 的 envelope/halo 概念不同：前者是跨层重叠块，后者是同层规则块周边的通信缓冲格点。

生成器最后构造 `OctreeGridStructure`，释放临时树 `freeOctree(root)`。求解器按层调用 `getMultiBlockManagement`，把同层同进程的可合并规则块合并，形成 `SparseBlockStructure3D` 与进程归属。运行时主要操作规则块数组和接口处理器。

示例 `examples/showCases/gridRefinement3d/offLatticeExternalFlow.cpp` 在主循环前从 `gridDensityFunctionFile` 生成树和固定层级，然后调用 `MultiLevelCoupling3D::collideAndStream` 推进。这条示例路径没有在每一步自动重新细化/粗化；`refinementCriteria3D` 中存在指标计算函数，不能据此认定该示例已经实现完整动态 AMR 重建。

证据：[octreeGridGenerator.hh](../Palabos/src/gridRefinement/octreeGridGenerator.hh) 第 584–602、2082–2107 行；[octreeGridStructure.cpp](../Palabos/src/gridRefinement/octreeGridStructure.cpp) 第 371–422 行；[offLatticeExternalFlow.cpp](../Palabos/examples/showCases/gridRefinement3d/offLatticeExternalFlow.cpp) 第 321–352、771、793–798、1044–1045 行。

## 七、 三维八叉树：粗细网格具体交换什么

### 7.1 交换的是分解并换算后的状态

`DecomposeAndRescaleFunctional3D` 调用 Engine 的 `decomposeAndRescale`，把状态写入辅助 `MultiNTensorField3D`。接收侧 `RecomposeFunctional3D` 调用目标格点动力学的 `recompose`，恢复目标层的分布函数。

对于等温基本动力学，源码实际布局是：
- `order=0`：`[rhoBar, j_x, j_y, j_z, fneq_0 ... fneq_(q-1), external...]`；
- `order=1`：`[rhoBar, j_x, j_y, j_z, PiNeq的6个独立分量, external...]`。

这里 `rhoBar` 是 descriptor 使用的密度变量，`j` 是动量；源码若干注释将其笼统称作 density/velocity，分析时应以实际赋值为准。`order` 控制状态分解方式，不是空间插值阶数。示例先写 `order=0`，随后在真正构造 lattice 前又赋值为 1，因此该示例实际采用 `order=1`。

证据：[gridRefinementFunctional3D.hh](../Palabos/src/gridRefinement/gridRefinementFunctional3D.hh) 第 1367–1404、1442–1467 行；[isoThermalDynamics.hh](../Palabos/src/basicDynamics/isoThermalDynamics.hh) 第 124–148、178–215 行；[offLatticeExternalFlow.cpp](../Palabos/examples/showCases/gridRefinement3d/offLatticeExternalFlow.cpp) 第 784–798 行。

以示例的 `ConvectiveNoForceRescaler` 为例，令 `r=xDt=目标层dt/源层dt`：
\[
\omega_t=\frac{r\omega_s}{1-\tfrac12\omega_s(1-r)},\qquad
q^{neq}_t=\frac{\omega_s}{\omega_t}r\,q^{neq}_s.
\]
密度量、动量在这种对流缩放中不缩放，非平衡部分按上述比例缩放。粗到细 `r=0.5`，细到粗 `r=2`。该关系针对这里的无力缩放器，其他动力学、外力或温度方程需要分别检查，不能机械套用。

证据：[rescaleEngine.hh](../Palabos/src/gridRefinement/rescaleEngine.hh) 第 55–77、85–103 行；[couplingInterfaceGenerator3D.hh](../Palabos/src/gridRefinement/couplingInterfaceGenerator3D.hh) 第 588–595、633–636 行。

### 7.2 粗到细：空间和时间都要补点

每对相邻层建立 4 个辅助状态场：`decomposed_t0`、`decomposed_t1`、`decomposed_t12`、`decomposed_fine`。粗层辅助场采用 2 格 envelope，细层辅助场采用 1 格 envelope。它们只覆盖接口所需的区域，按规则块分布对齐。

证据：[couplingInterfaceGenerator3D.hh](../Palabos/src/gridRefinement/couplingInterfaceGenerator3D.hh) 第 97–247 行，尤其 `216-247`。

时间插值逐分量计算：
\[
q_c^{n+1/2}=\tfrac12(q_c^n+q_c^{n+1}).
\]
空间插值则把粗层接口上的换算后状态送到细层：
- 重合位置直接复制；
- 一个切向方向的半格位置用中心 4 点；
- 接口面内两个切向方向都位于半格的位置用 16 点；
- 靠近棱、角、物理边界采用专门的缩减/偏置模板，不能一律使用内部 16 点模板。

中心 4 点公式为
\[
q_{i+1/2}=\frac{-q_{i-1}+9q_i+9q_{i+1}-q_{i+2}}{16}.
\]
这是对分解状态的插值，而不是单独插值速度后忽略非平衡信息。

证据：[gridRefinementFunctional3D.hh](../Palabos/src/gridRefinement/gridRefinementFunctional3D.hh) 第 87–126、134–300、763–799 行；[fdStencils1D.h](../Palabos/src/finiteDifference/fdStencils1D.h) 第 121–127 行；[fdStencils2D.h](../Palabos/src/finiteDifference/fdStencils2D.h) 第 47–67 行。

### 7.3 一个粗步中的实际调度

`MultiLevelCoupling3D::collideAndStream` 的实际执行顺序如下（再细一层仍按同一方式递归）：

1. 粗层碰撞迁移，从 `t_n` 推进到 `t_(n+1)`。
2. 分解、缩放粗层新状态，与保存的旧状态平均，得到半时刻状态。
3. 细层递归推进第 1 个子步。
4. 用粗层半时刻状态做空间插值，重构细层粗细接口。
5. 细层递归推进第 2 个子步。
6. 用粗层新时刻状态做空间插值，再重构细层接口。
7. 细层状态分解并按 `xDt=2` 换回粗层尺度，过滤/复制后重构粗层接口。
8. 交换粗层旧、新状态缓冲，为下一粗步保留数据。

注意：源码在每个细子步推进之后更新接口，阅读或移植时应保持与碰撞迁移、处理器顺序一致，不能仅凭常见 AMR 示意图把这些步骤随意前移。两层满足 `dx_f=dx_c/2,dt_f=dt_c/2`；三层时一个最粗步对应中层 2 步、最细层 4 步。

证据：[couplingInterfaceGenerator3D.hh](../Palabos/src/gridRefinement/couplingInterfaceGenerator3D.hh) 第 674–818 行，尤其第一个子步 `720-735`、第二个子步 `752-766`、回传 `780-818`。

### 7.4 细到粗：回传与过滤

细层推进两个子步后，`CopyAndFilterNonEquilibriumFunctional3D` 从重合细节点及其 26 个邻点取值，对选中的分解分量做 27 点等权平均，再重构粗层接口：
\[
q_c(I,J,K)=\frac1{27}\sum_{a,b,c=-1}^{1}q_{f\to c}(2I+a,2J+b,2K+c).
\]
实际局部坐标还要减去 block 的 location。

虽然类名包含 NonEquilibrium，构造参数 `filterAll` 才决定过滤范围：
- 默认 `filterAll=true`：从分量 0 开始，密度量、动量、非平衡分量都参与平均；
- `filterAll=false`：密度量和动量从重合细节点复制，只过滤非平衡部分；
- external scalars 在此实现中直接复制，不参与这段 27 点平均；
- 某些边界棱使用 `CopyFunctional3D` 直接复制，因为完整过滤模板不可用。

27 点是过滤模板，不表示流体碰撞模型必须使用 D3Q27。这与二维 `multiGrid` 的默认过滤范围不同。

证据：[couplingInterfaceGenerator3D.h](../Palabos/src/gridRefinement/couplingInterfaceGenerator3D.h) 第 244–248 行；[couplingInterfaceGenerator3D.hh](../Palabos/src/gridRefinement/couplingInterfaceGenerator3D.hh) 第 629–666 行；[gridRefinementFunctional3D.hh](../Palabos/src/gridRefinement/gridRefinementFunctional3D.hh) 第 1167–1215、1311–1338 行。

## 八、 树加密和 MPI 通信的关系

### 8.1 跨层算术主要在本地处理，MPI 同步规则块的 envelope

`MultiProcessing3D::subdivideGenerator` 收集各参与 MultiBlock 的本 rank 块，把不同层坐标缩放到同一参考层后求交；`addInternalProcessor` 提取本地 AtomicBlock 指针，安装粗细耦合处理器。因此要保证粗细耦合区域在进程划分上兼容。八叉树生成器的 `EnforceOctreeProcessAssignmentConstraints` 会把交界处部分细孩子分配到 overlap 父块所在 rank；其余块再按负载分配。

这解释了实现层次：树确定空间拓扑和进程归属；处理器在本地执行分解、插值、过滤、重构；MultiBlock 通信器同步 lattice 和辅助状态场的 envelope。不能把它描述成每个四叉树/八叉树父子节点逐个调用 MPI 发送整个子树。

证据：[multiBlockOperations3D.hh](../Palabos/src/multiBlock/multiBlockOperations3D.hh) 第 105–183 行；[multiBlockOperations3D.cpp](../Palabos/src/multiBlock/multiBlockOperations3D.cpp) 第 133–159 行；[octreeGridGenerator.hh](../Palabos/src/gridRefinement/octreeGridGenerator.hh) 第 759–803、880–895 行。

### 8.2 通信如何自动接到处理器后面

`MultiBlock3D::executeInternalProcessors(level,communicate)` 先遍历本地块执行处理器；若 `communicate=true`，则按处理器声明的修改范围调用 `duplicateOverlapsInModifiedMultiBlocks`。例如分解处理器声明张量的 `staticVariables` 被修改；重构处理器声明 lattice 的 `staticVariables` 被修改。

`couplingActionsGenerator3D.hh` 还提供 Actions 形式的调度，它显式写出 `addCommunication`：粗层分解后同步、插值后的细张量同步、细层回传分解后同步，以及粗层过滤后的张量同步。它是辅助理解通信时机的另一实现入口；上文示例实际使用的是 `MultiLevelCoupling3D` 路径。

证据：[multiBlock3D.cpp](../Palabos/src/multiBlock/multiBlock3D.cpp) 第 537–553、631–649 行；[gridRefinementFunctional3D.hh](../Palabos/src/gridRefinement/gridRefinementFunctional3D.hh) 第 1420–1424、1478–1488 行；[couplingActionsGenerator3D.hh](../Palabos/src/gridRefinement/couplingActionsGenerator3D.hh) 第 260–318 行。

### 8.3 一次跨块通信的具体过程

`CommunicationStructure3D` 把相交区域转换为源、目标块的局部坐标，记录 `fromBlockId/toBlockId/fromDomain/toDomain/fromProcessId/toProcessId`，按归属分为本地复制、发送、接收三类包。

`ParallelBlockCommunicator3D::communicate` 依次执行：

1. 对已知长度的数据预挂非阻塞接收。
2. 源块 `getDataTransfer().send` 把指定接口/envelope区域打包，交给发送池进行非阻塞发送。
3. 同 rank 的块之间调用 `attribute` 直接复制。
4. 等待远程接收完成，用 `getDataTransfer().receive` 解包写入目标区域。
5. 等待本轮发送结束。

`SendPoolCommunicator` 将发往同一 rank 的多条小消息合并到同一个字节向量，再调用 `iSend`；静态消息用预先知道的长度配合 `iRecv`。动态消息还先传长度，接收路径不同。这里的“staticVariables”指固定结构的数据种类，不表示流场数值在时间上不变。

因此通信内容由被同步对象决定：lattice 同步它的格点静态数据；辅助 NTensor 同步分解、换算后的状态分量。MPI 传的是指定局部区域的数据缓冲，插值和尺度换算由耦合处理器完成。

证据：[parallelBlockCommunicator3D.cpp](../Palabos/src/parallelism/parallelBlockCommunicator3D.cpp) 第 58–116、207–234、253–309 行；[sendRecvPool.cpp](../Palabos/src/parallelism/sendRecvPool.cpp) 第 101–165 行。非 MPI 编译时选择串行通信器，MPI 编译时可选择并行或 blocking 实现，见 [defaultMultiBlockPolicy3D.h](../Palabos/src/multiBlock/defaultMultiBlockPolicy3D.h) 第 66–76 行。

## 九、 建议的阅读顺序与适用边界

- 需要二维粗细网格耦合：先读 `examples/showCases/gridRefinement2d/dipole.cpp`，再读 `multiGridManagement2D.cpp`、`multiGridLattice2D.hh`、`fineGridProcessors2D.hh`、`coarseGridProcessors2D.hh` 和 `gridRefinementDynamics.hh`。
- 需要树状网格生成、平衡和邻居查找：先读 `gridRefinement/octree.hh`、`octreeGridGenerator.hh`，再读 `octreeGridStructure.cpp`。
- 需要完整三维耦合调度：读 `couplingInterfaceGenerator3D.hh` 的 `integrateProcessingFunctionals` 与 `collideAndStream`，配合 `gridRefinementFunctional3D.hh` 和 `rescaleEngine.hh`。
- 需要 MPI 数据如何真正搬运：读 `multiBlockOperations3D.hh` 的本地求交、`multiBlock3D.cpp` 的处理器后通信，以及 `parallelBlockCommunicator3D.cpp`、`sendRecvPool.cpp` 的打包收发。

这是基于当前本地提交的静态源码分析，未编译运行示例，也没有据此验证插值守恒性、稳定性、计算精度或热对流模型适用性。上游 Palabos 源码保持不变。

## 十、与 Lagrava 等人 JCP 2012 论文的对应关系

2026-11-09 补充核对：二维的三次插值加非平衡过滤路径，对应 D. Lagrava、O. Malaspinas、J. Latt、B. Chopard 的 *Advances in multi-domain lattice Boltzmann grid refinement*，Journal of Computational Physics **231**(14), 4808–4822 (2012)，[DOI: 10.1016/j.jcp.2012.03.015](https://doi.org/10.1016/j.jcp.2012.03.015)。[出版方摘要](https://www.sciencedirect.com/science/article/pii/S002199911200157X)明确说明该方法已在 Palabos 中实现；[日内瓦大学收录的作者稿](https://archive-ouverte.unige.ch/unige:34397)第 3.3 节式 (31)–(33) 给出细到粗非平衡缩放和过滤，第 3.4–3.6 节说明粗到细重构、两细步耦合和空间插值。

当前本地代码的直接证据是 [dipole.cpp](../Palabos/examples/showCases/gridRefinement2d/dipole.cpp) 第 334–336 行调用 `createRefinedLatticeCubicInterpolationFiltering`，其二维 generator 安装 `CopyCoarseToFineCubicInterp2D` 与 `CopyFineToCoarseWithFiltering2D`，见 [multiGridGenerator2D.hh](../Palabos/src/multiGrid/multiGridGenerator2D.hh) 第 238–278、325–340 行。细到粗只对非平衡分量进行九点平均、保留重合点的密度量和动量，与论文的选择相符；external scalars 是代码通用接口的额外分量，不能据此当作原论文已验证内容。

这里要区分“粗细接口数值耦合方法”与“树结构、块管理、MPI 数据搬运”。Lagrava 论文支持前者的来源；八叉树平衡、邻接表和通信池属于本地实现的组织机制。三维路径沿用相关耦合思想，但当前默认 `filterAll=true` 还过滤密度量和动量，而原论文明确选择只过滤非平衡部分，因此不能将当前三维全部默认设置视作该二维论文的逐项复现。三维这一差异见本报告第 7.4 节。


