# 2DRBOpenaccMultiblock

`2DRBOpenaccMultiblock.F90` 从 `../均匀网格/2DRBOpenacc.F90` 复制派生，原文件未修改。

参考文献：Huang, R. and Wu, H., *Multiblock approach for the passive scalar thermal lattice Boltzmann method*, Physical Review E **89**, 043303 (2014), DOI: [10.1103/PhysRevE.89.043303](https://doi.org/10.1103/PhysRevE.89.043303)。本地全文为 `../pdf/Multiblock approach for the passive scalar thermal lattice Boltzmann method.pdf`，主要依据第 III 节式 (19)-(22)。

## 保留的原算法

- D2Q9 流场 MRT：原矩基、平衡矩、Guo 浮力源、宏观速度半步力修正。
- D2Q5 温度 MRT：默认 `EnableUseG`，保留 `T*u/T*v` 的历史差分修正；也保留原 `EnableLegacyThermalScheme` 分支。
- 每块执行原顺序：`collision → streaming → bounceback → macro → collisionT → streamingT → bouncebackT → macroT`。
- 保留原流体节点距物理壁面半个格距的布置以及物理壁面的 BB/ABB。
- 默认仍为原程序的封闭 Rayleigh-Bénard 算例：`Ra=1e7, Pr=0.7, Ma=0.1`，上下恒温，侧壁绝热，四壁无滑移，初始温度扰动与原程序相同。

论文使用 D2Q9 温度模型。本文件把其一般 MRT 接口关系应用于原 D2Q5 矩基，未把块内温度模型换成论文的 D2Q9。文献的局部热流后处理公式也没有直接套入本文件。

## 网格参数

参数集中在文件开头的 `module commondata`。

| 参数 | 默认值 | 含义 |
|---|---:|---|
| `nx,ny` | 1024,1024 | 最细格距下的等效全域分辨率；也支持 `NX_OVERRIDE/NY_OVERRIDE` |
| `refineRatio` | 2 | 粗细格距及时间步之比；支持 2，设 1 时退化为原单块均匀网格 |
| `wallCellsX,wallCellsY` | `nx/8,ny/8` | 左右、上下细网格层的物理厚度，以最细格距计 |
| `overlapCells` | 4 | 从接口参考节点向外延伸的粗格距数；每侧延伸 `4×2=8` 个细格距 |
| `interfaceSkin` | 2 | 每次推进后由相邻块重建的人工边界节点层数，以各块自身格距计 |
| `loadInitField` | 0 | 设 1 读取本程序的 `latest.meta` 并精确续算 |
| `unsteadyRunDuration` | 1000 | 非稳态绝对目标时间，单位 `t_ff` |

默认五块的**物理积分分区**如下；分区面积互不重叠。各块计算节点可以重合，统计按裁剪后的控制面积加权，避免面积重复计入。

```text
  +-----------------------------------+
  |             3: 上细块             |
  +--------+-----------------+--------+
  |        |                 |        |
  |4: 左细 |    1: 中心粗    |5: 右细 |
  |        |                 |        |
  +--------+-----------------+--------+
  |             2: 下细块             |
  +-----------------------------------+
```

每个块内部是等距正方形格子，全域通过不同格距的块实现局部加密。默认中心 `h=2`、四周 `h=1`；最细格距为长度单位 1。

`refineRatio=2` 时，`nx,ny,wallCellsX,wallCellsY` 必须为偶数；壁面层厚度至少为 `2*overlapCells+2` 个细格距，中心所有权区至少容纳 `8×8` 个粗格子。默认 `nx/8` 的设置适合从 `96×96` 开始做小规模检查。

### 2026-09-12：改为粗细节点对齐

细格点保持原程序的半格距位置，粗格点成为细格点的子集：

\[
x_f=0.5+k,\qquad x_c=0.5+2k,\quad k\in\mathbb Z,
\]

`y` 同理。代码内部仍用 `x=x0+(i-0.5)*h`，但 `x0=xb+0.5-h/2`；因此粗块的 `x0` 是虚拟控制区的左侧面，不能再当作物理积分分区边界。真实壁面附近仅使用细块，第一、最后节点仍为 `0.5`、`N-0.5`，原 BB/ABB 不变。

以左细块与中心粗块为例：物理积分分界是 `x=128`，对应的共同节点参考列是 `x=128.5`。粗块从 `120.5` 开始，细块延伸至 `136.5`，重叠节点范围的跨度为 `16=2×overlapCells×refineRatio`。其中 `120.5,122.5,…,136.5` 与共同的粗节点行相交处，都是粗细重合点。两级仍各自存储状态，交换时作矩重标定。`interfaceSkin=2` 对应粗块首两列 `120.5,122.5`、细块末两列 `135.5,136.5`。

| 块 | `h` | 计算节点 x 范围 | 计算节点 y 范围 | 实际推进节点数 | 正面积积分节点数 |
|---|---:|---|---|---|---|
| 中心粗块 | 2 | 120.5–904.5 | 120.5–904.5 | 393×393 | 385×385 |
| 下细块 | 1 | 0.5–1023.5 | 0.5–136.5 | 1024×137 | 1024×128 |
| 上细块 | 1 | 0.5–1023.5 | 888.5–1023.5 | 1024×136 | 1024×128 |
| 左细块 | 1 | 0.5–136.5 | 120.5–904.5 | 137×785 | 128×768 |
| 右细块 | 1 | 888.5–1023.5 | 120.5–904.5 | 136×785 | 128×768 |

计算节点共 648306 个，正面积积分样本共 606977 个；后者不是去重后的坐标数。上/下、左/右计算节点数相差一列/行，是共同半格距相位与人工端点包含规则所致，物理分区仍对称。接口历史数组和插值增加了开销，不能仅凭节点数推断 GPU 加速比。

图示见 `图示/multiblock-layout.png` 和 `图示/multiblock-interface.png`；对应 SVG 可缩放。运行 `python 图示/draw_layout.py` 可重绘，并生成包含源码 SHA-256 的 `layout_parameters.json`。放大图中的蓝色方框包围橙色圆点，表示同一坐标上的粗细节点。

当前是静态两级布局。多块模式支持原封闭 RB 和侧壁差温两类物理边界；周期侧边需要专门的周期块连接，目前多块模式会明确拒绝此配置。`refineRatio=1` 保留原周期宏路径。

## 粗细块的物理单位与松弛率

设块 (b) 的格距、时间步为

\[
\Delta x_b=\Delta t_b=h_b,\qquad c=1.
\]

原文件计算的 `viscosity,diffusivity,gBeta,Snu,Sq,Qk,Qnu,paraA` 均作为最细格距下的参数；全域共用原来的 `lengthUnit,timeUnit`。每个非守恒矩的松弛率满足论文式 (20)：

\[
h_b\left(\frac1{s_b}-\frac12\right)
=\left(\frac1{s_f}-\frac12\right),\qquad
s_b=\left[\frac12+\frac{1/s_f-1/2}{h_b}\right]^{-1}.
\]

流场的 `Snu` 与 `Sq` 均作这一换算；**不能只换剪切矩，也不能在粗块重新套细块的 `Sq(tauf)` 魔术参数公式**。温度场的 `Qk,Qnu` 同样换算。`paraA` 及热平衡矩在两级间保持相同。

因此下式在各块相同：

\[
\nu=h_b\frac13\left(\frac1{S_{\nu,b}}-\frac12\right),\qquad
\kappa=h_b\frac{4+A_T}{10}\left(\frac1{Q_{k,b}}-\frac12\right).
\]

当前 D2Q5 分支 (A_T=-2/3)，旧分支 (A_T=\mathrm{paraA})。块内每步的浮力增量使用 `h_b*gBeta`，`SideHeatedHa` 分支的每步磁场源系数也乘 `h_b`。各块的速度单位相同。

## 接口交换

代码使用**迁移后的分布函数**构造宏观量和矩。设 `Fm` 为尚未乘 `(I-S/2)` 的力矩增量，论文式 (19) 可写为块间公共量

\[
K_a=\frac{s_{a,b}}{h_b}\left(m_{a,b}-m^{eq}_{a,b}
                    +\frac12 F^{LB}_{m,a,b}\right).
\]

插值传递 `rho,u,v,T,Fx/h,Fy/h,Kf,Kg`。目标块按

\[
m_{a,d}=m^{eq}_{a,d}+\frac{h_d}{s_{a,d}}K_a
                         -\frac12F^{LB}_{m,a,d}
\]

重建全部非守恒矩。守恒矩直接恢复为

\[
m_0=\rho,\qquad m_3=\rho u-F_x^{LB}/2,\qquad
m_5=\rho v-F_y^{LB}/2,\qquad n_0=T.
\]

`EnableUseG` 的两个非守恒热流矩也有源修正。接口保留原离散历史，通过当前 (B=(uT,vT)) 与上次 `collisionT` 保存的 `B_prev` 定义

\[
D_B=(B-B_{prev})/h_b,\quad
R^{LB}_T=(0,h_bD_{B,x},h_bD_{B,y},0,0),\quad
B_{prev,d}=B_d-h_dD_B.
\]

`Kg` 的定义相应包含 `R_T^LB/2`。这是对原 D2Q5 离散修正及其分裂推进时序的适配；论文没有直接给出这个 `B_prev` 接口实现。纯导热和有限时间对照已覆盖基础正确性，更严格的时间/空间收敛阶仍需单独验证。

空间上区分两类目标：**重合节点直接读取同坐标来源节点的交换量**；非重合细节点使用四点 Lagrange 插值，二维采用张量积的 `4×4` 模板。随后均按目标块松弛率重建矩，不能直接复制跨级的原始 `f,g`。粗块从细块获取的人工边界节点必须全部重合，建连接时对此作显式检查。同级重合点也走直接读取路径。

来源节点或插值模板严格避开来源块两层人工边界。人工边界两层在完整的流动/温度推进后统一重建，覆盖分裂推进可能污染的节点；真实物理壁面仍执行原 BB/ABB。

时间上先得到粗块 (t+\Delta t_c) 的预测，细块执行两个子步。对 (θ=(t_f-t)/\Delta t_c)，三时间层 (t-\Delta t_c,t,t+\Delta t_c) 的权重为

\[
w_-={\theta(\theta-1)\over2},\quad w_0=1-\theta^2,
\quad w_+={\theta(\theta+1)\over2}.
\]

首个粗步缺少负时间层，采用线性启动；其后使用上述三时间层公式。两个细步完成后，再从细块修复粗块人工边界，最后滚动粗块历史。插值先读固定的来源快照，再更新全部目标节点，避免块处理次序改变结果。

## 统计、输出与重启

所有体积分使用 `ownedBox` 定义的物理分区。每个节点以自身为中心的 `h×h` 控制区与该分区取交集，得到 `A_ij=dxWeight(i)*dyWeight(j)`。内部权重仍为 `h²`，分区边缘使用裁剪后的面积：

\[
Nu_V=1+\frac{L}{\kappa\Delta T}
\frac{\sum_b\sum_{A_{ij}>0}u_\parallel T A_{ij}}{L_xL_y},
\quad Re_V=\frac{L}{\nu}
\sqrt{\frac{\sum_b\sum_{A_{ij}>0}(u^2+v^2)A_{ij}}{L_xL_y}}.
\]

例如中心块在 `[128,896]` 上使用 `128.5,130.5,…,896.5` 共 385 个积分节点；首、末一维权重分别为 `1.5` 和 `0.5`，其余为 `2`，总长度为 `768`，线性坐标积分也精确。节点 `896.5` 同时存在于右细块中，但粗块代表 `[895.5,896]`，右细块代表 `[896,897]`，两者积分面积不重叠。不能仅按坐标去重，也不能把每个样本都乘 `h²`。

热壁、冷壁沿用原半步长壁面温度梯度公式；壁面和中线沿切向按裁剪长度加权。中线值及梯度在真实中线坐标处由四点多项式求得。稳态收敛误差也按面积加权。

- `NuRe_2DOpenaccMultiblock.dat`：实际 `t_ff,NuVolAvg,ReVolRMS,Nu_hot,Nu_cold,Nu_middle,mass,meanT,Tmin,Tmax,rhoMin,rhoMax`。
- `NuReStatistics_2DOpenaccMultiblock.dat`：指定时间窗和前/后半窗的平均值、相对差异；对 `Re²` 时间积分后开根号。时间窗覆盖不完整时明确输出 `INCOMPLETE`。
- `...Tecplot-*.dat`：ASCII 多 zone 格式，只写正面积积分节点；含物理坐标、局部格距和第八列 `integration_area/L^2`。原单块 `.plt` 二进制读取器不能直接读取这个格式。
- `...Snapshot-*.bin`：多块 stream 格式，16 字节 `MB2DSNAPSHOT0002` 标识（空格补齐），随后为 4 个 int32 `nBlocks,nx,ny,itc` 和 2 个 float64 `t_ff,lengthUnit`。每块依次写 2 个 int32 `ni,nj`，7 个 float64 `first_node_x,first_node_y,h,xmin,xmax,ymin,ymax`，再写一维 `dxWeight(ni),dyWeight(nj)`，最后按 Fortran 顺序写 `u,v,T,rho` 四个二维数组。坐标原点现在是首个积分节点，物理分区单独记录。
- `reloadFile2DOpenaccMultiblock-*.bin` 与 `-latest.meta`：完整保存各块 `f,g`、宏观量、力、`B_prev`、粗时间层、采样编号和输出编号。仅在粗细同步时写入，网格、算法、参数和历史样本数不匹配会拒绝重启。当前标识为 `MB2DRESTART0002`，增加物理积分分区元数据。

**旧版错位网格的检查点和快照不兼容本版，不能继续旧算例；需新建运行目录从初场开始。** 本目录两个验证 JSON 和图示已按对齐后的源码重新生成。

输出时刻和最终目标时刻向上对齐到粗细同步步，时间列写实际时刻。时间间隔至少为一个粗步。检查点完整写完后才更新 `latest.meta`，旧编号检查点保留。续算可延长总目标时间，需同时保留对应的 Nu/Re 历史。

原程序只适用于单块均匀网格的流函数/涡量、壁面极值五点拟合后处理未直接搬入本文件；当前提供上述适用于原生多块网格的输出。所有输出使用新的 `Multiblock` 前缀。

## 编译与检查

在运行目录编译，避免输出落入源文件目录。NVIDIA 编译命令示例（本地尚未执行）：

```sh
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccMultiblock.F90 -o 2DRBOpenaccMultiblock
./2DRBOpenaccMultiblock
```

本地 OpenACC **host** 检查：

```sh
gfortran -cpp -fopenacc -ffree-line-length-none -O2 2DRBOpenaccMultiblock.F90 -o 2DRBOpenaccMultiblock
```

可复现的轻量验证（所有临时改参源码、模块、可执行文件和运行结果均写入系统临时目录，原文件不会修改）：

```sh
python verify_multiblock.py
python transient_compare.py
```

结果记录在 `verification_results.json`、`transient_comparison.json`，并包含源文件 SHA-256。

验证范围包括五种配置的 OpenACC 语法、方形/矩形域节点对齐与积分、接口直接传值和时空插值、单块退化回归、两种 D2Q5 分支的精确线性导热、完整状态及实际主程序输出历史的断点续算、旧检查点拒绝和新版快照读取。另以 `96², Ra=1e4` 的侧壁差温有限时间瞬态，与原程序作同时间对比。

本次本机结果：

| 检查 | 结果 |
|---|---|
| 五种宏配置的 `gfortran -fopenacc` 语法 | 通过 |
| `96×96` 与 `128×96` 的节点和积分几何 | 粗节点均属细节点子集；物理壁面位置、总面积与线性坐标积分通过 |
| 重合点直接传值、非重合点时空插值 | 指定多项式的最大误差为 0；包含数组末端直接读取的边界检查 |
| 单块退化，40 个细步，RB `EnableUseG` | 完整状态最大绝对差 `8.88e-16` |
| 单块退化，40 个细步，侧壁差温旧 D2Q5 | 完整状态逐值一致 |
| 五块线性导热，300 个粗步 | 最大温度误差不超过 `1.73e-15`；体平均/热壁/冷壁/中线 Nu 均为 1 至舍入误差 |
| 40 个细步在第 20 步重启 | 与连续运行完整状态逐值一致，包括粗块历史量 |
| 实际主程序运行及续算 | 10 个 Nu/Re 样本逐值一致，统计窗完整覆盖 |
| 文件格式保护 | 旧版检查点被拒绝；v2 快照坐标和积分面积通过 |
| `Ra=1e4` 侧壁差温，2000 个细步，`t_ff=1.202813` | 相对均匀网格：Nu 差 `0.0613%`，Re 差 `0.0731%`，热壁/冷壁 Nu 差 `0.00756% / 0.00887%` |

最后一行是**瞬态一致性检查**，不是稳态基准。比较时先把均匀网格场映射到多块原生节点，再按不重复物理面积比较；`u,v` 相对 L2 差分别约 `0.0966%,0.0983%`，温度相对 L2 差约 `0.0302%`。

这些是本机 CPU 上的检查；尚未进行 `nvfortran/P100` GPU 编译运行、长期稳态或统计收敛验证，也未确认高 `Ra` 稳定性和网格/时间收敛阶。重叠插值不是有限体积的严格通量守恒修正，因此仍需检查长时间质量漂移和热通量一致性。
