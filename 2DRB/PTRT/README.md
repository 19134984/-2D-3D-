# P-TRT 流场移植记录

## 当前状态

2026-09-18：已从 `../均匀网格/2DRBOpenacc.F90` 复制 `2DRBOpenacc.F90`，
并按用户确认的方案完成 **P-TRT D2Q9 流场 + Luo/Wang D2Q5 温度场**。
采用守恒一致的 Guo/MRT 体力处理。均匀网格原文件未修改。

默认仍为四壁无滑移 Rayleigh–Bénard 非稳态算例，`1024×1024`、`Ra=1e7`、
`Pr=0.7`、`Ma=0.1`、目标 `1000 t_ff`；仅测试临时副本使用小网格和短时长。
没有提交服务器作业。

文献：Yuan Yu et al., *Purified Two-Relaxation-Time Lattice Boltzmann Method:
Removing Ghost Modes from TRT for Enhanced Stability*, arXiv:2602.06686v1,
2026-02-06，https://arxiv.org/abs/2602.06686v1 。依据用户上传的本地 PDF。

## 可以直接采用的 D2Q9 部分

取格子单位 `dt=dx=1`、`cs2=1/3`，速度编号与现有程序完全一致。

- 黏性松弛率 `Snu=1/tauf`，`viscosity=(tauf-0.5)/3`。
- 保留现有奇松弛率 `Sq=8*(2*tauf-1)/(8*tauf-1)`，满足
  `(1/Snu-0.5)*(1/Sq-0.5)=3/16`。
- 式(3)采用三阶 Hermite 平衡分布。D2Q9 中可写成
  `feq_i=wi*rho*[1+3*eu+4.5*eu^2-1.5*u2+4.5*eu^3-4.5*eu*u2]`。
- `h_i=f_i-feq_i`；`hplus_i=(h_i+h_opp)/2`；`hminus_i=(h_i-h_opp)/2`。
- 式(17)、(24)：`phi_i=(ex_i^2-1/3)*(ey_i^2-1/3)`，
  `SG=sum(phi_i*h_i)`，`PG_i=(81/4)*wi*phi_i*SG`。
  `PG` 是唯一的四阶偶幽灵模态；`R2=hplus-PG` 是二阶应力投影。
- 无体力且宏观量与分布一致时，式(30)为
  `fpost=feq+(1-Snu)*R2+(1-Sq)*hminus`。

## 论文印刷体力公式的离散矩核对

这是对上传 v1 公式的直接代数检查，不是对作者实际计算程序的判断。
已对照 PDF 第4、10页图像，排除文本抽取造成的括号或下标错误。

令 `a=1-Snu/2`、`b=1-Sq/2`。式(6)原文等价于

`S_i=wi*[a*(e_i.u)*(e_i.F)/cs2^2+b*(e_i.F-u.F)/cs2]`。

利用 D2Q9 求积恒等式直接得到：

`sum(S_i)=(a-b)*(u.F)/cs2=(Sq-Snu)*(u.F)/(2*cs2)`。

该量一般不为零，而浮力不应增加质量。

同时，式(9)要求 `sum(e_i*h_i)=-F/2`。式(29)把奇非平衡部分的一阶矩去掉，
因此式(30)中两个正则化项的一阶矩均为零。结合式(6)有

`p_post=rho*u+(1-Sq/2)*F`，

而 `p_pre=rho*u-F/2`，故

`p_post-p_pre=(3/2-Sq/2)*F`，一般不等于 `F`。

例：`tauf=0.6`，`Snu=5/3`，`Sq=8/19`，实际动量增量为 `49/38*F`。
这不仅是平衡态三阶截断误差，而是局部碰撞的零阶/一阶矩不一致。
换用式(7)本身也不能自动消除一阶矩问题，故不以另一印刷体力式替换了事。

复现：`python PTRT/verification/audit_paper_force.py`，使用精确有理数计算。

## 已采用的守恒一致方案

保留现有 MRT 对守恒动量矩 `s_j=0` 的处理；对二阶矩使用 `Snu`，
对三阶矩使用 `Sq`，对四阶幽灵矩使用完全消除（松弛率为1）。
不引入论文 Scheme II 的速度三次项梯度修正。

令

`Jneq=sum(e_i*h_i)`，`P1_i=wi*(e_i.Jneq)/cs2`，`R3=hminus-P1`，

`A_i=wi*(e_i.F)/cs2`，

`B_i=wi*[(e_i.u)*(e_i.F)/cs2^2-(u.F)/cs2]`。

实际碰撞采用

`fpost=feq+P1+(1-Snu)*R2+(1-Sq)*R3+A+(1-Snu/2)*B`。

这样不碰撞守恒动量，并且源项严格满足 `sum(source)=0`、
`sum(e_i*source)=F`；无体力时退化为论文的幽灵模态滤除形式。
二阶源矩为 `(1-Snu/2)*(uF+Fu)`，与现有 Guo/MRT 应力源矩一致。
显式计算 `Jneq` 还避免假设推进中刚更新的温度浮力与上一时刻的宏观半步力完全相同。

本方案是 **P-TRT 幽灵模态滤除 + 守恒一致 Guo/MRT 体力**，
不能描述成对论文有力公式(6)+(30)的逐字复现。

## 已选用的温度分支

复制源默认是 `EnableUseG`。之前的 Luo/Wang D2Q5 温度算法对应
`EnableLegacyThermalScheme`：

`Qk=3-sqrt(3)`，`Qnu=4*sqrt(3)-6`，
`paraA=20*sqrt(3)*diffusivity-4`。

经用户确认选用后者。温度碰撞与边界函数沿用基线，
不向 Luo 分支混入相邻时间步 `T*u/T*v` 修正，也不引入 Xs 的 chi 参数。

## 修改范围

- 流场 `collision()`：三阶平衡分布、奇偶分解、四阶幽灵投影扣除、
  一阶守恒矩保留、二阶/三阶松弛和守恒体力。
- `initial()`：流场初始化改为相同三阶平衡分布；记录流场模型、体力方案与松弛率。
- 增加 D2Q9 反方向表 `opposite`，纳入 OpenACC 数据区；临时标量和数组均为格点私有。
- 默认关闭 `EnableUseG`，开启 `EnableLegacyThermalScheme`。
- 黏度、扩散率、浮力定义、推进顺序、边界、温度碰撞和后处理函数保持基线实现。
  流场 `streaming/bounceback/macro` 与温度 `collisionT/streamingT/bouncebackT/macroT`
  已由脚本逐段比较，内容相同。

保留基线输出命名约定，运行时使用单独结果目录，避免覆盖其他算法的数据。
续算应使用本算法、相同参数产生的检查点；尚未验证跨算法检查点转换，
也尚未对本算法执行完整的中断/续算对照验证。

## 已完成的验证

环境：Windows，GNU Fortran 15.2.0，`-fopenacc`，`ACC_DEVICE_TYPE=host`。
未使用 NVIDIA GPU 或 NVHPC；主机 OpenACC 测试不代表 GPU 验证。

1. 完整原尺寸主程序通过 OpenACC 语法检查及 `-O2` 编译。
2. `verification/verify_ptrt.py` 提取实际源码中的模块和碰撞函数编译执行，
   独立参考使用完整二阶、三阶 Hermite 张量收缩，不复用生产代码的幽灵扣除公式。
   三组黏性参数（测试网格16，对应 `Ra=1e3/1e6/1e10`）、浮力及倾斜磁场双分量体力、
   串行与 OpenACC 主机两种执行，共12组，每组256个局部状态：
   - 与独立投影的最大绝对差 `1.1102230246251565e-16`；
   - 质量守恒、动量增量 `F`、二阶源矩、三阶松弛率检查通过；
   - 四阶幽灵消除及注入任意幽灵模态后输出不变检查通过；
   - 串行/OpenACC 主机结果对照通过。
   高 Ra 在这里仅用于生成松弛率，不是高 Ra 热对流运行验证。
3. 调用实际初始化、设备数据管理和完整流热推进子程序，运行 `32×32、Ra=1e4`，
   2000格子步（约 `3.60844 t_ff`），同时对照原基线的 Luo 温度分支：

| 短程末态量 | 原基线 + Luo 温度场 | P-TRT + Luo 温度场 |
|---|---:|---:|
| 质量相对漂移 | -1.1169e-13 | -3.1408e-13 |
| 最小密度 | 0.99924648 | 0.99924589 |
| 最大密度 | 1.00040808 | 1.00040854 |
| 最大速度 | 7.20740e-5 | 7.20090e-5 |
| 瞬时体积 Nu | 0.99663254 | 0.99665860 |
| 瞬时体积 RMS Re | 0.09026034 | 0.09021962 |

分布和场量有限，密度为正；每套算法的串行/OpenACC 主机结果一致。
两种算法的场量最大绝对差约 `5.91e-7`，这是短程回归记录，不要求不同碰撞模型逐位相等。
没有发现守恒或主机执行差异；上述细小差异与碰撞模型改变相容，不能据此判定长期精度。
表中不是稳态或时间平均结果。

复现（仓库根目录）：

```powershell
python PTRT/verification/audit_paper_force.py
python PTRT/verification/verify_ptrt.py
```

编译/运行中间文件放在 ASCII 临时目录。验证结果、源码哈希和编译器信息保存在
`verification/results.json`，短程 Nu/Re 与场范围见 `verification/smoke_*.log`。

源码 SHA-256：

- 原基线：`68debf94f75f3a70af6bf4f7205e18afc14618d48dd227ccf1d1a8e73d8afd48`
- P-TRT：`390f8297ee60eec7e8b515e44d33570e77f0b583e512537d4a2562b97765d1cc`

实际检查：局部碰撞代数、OpenACC 主机执行、完整编译、小网格流热耦合短程回归。
仍未验证：P100/NVHPC 执行、高 Ra 长期稳定性、空间收敛、论文完整基准复现及严格续算。
问题分类：论文印刷体力项的不一致属于数值模型的离散矩问题；本次未发现主机并行执行差异，
后处理函数未改动，尚不能用短程结果断言整体稳定性提高。
