# 二维均匀网格耦合 DUGKS

`2DRBOpenaccDUGKS.F90` 实现 Wang、Tao、Guo 的 **A coupled discrete unified gas-kinetic scheme for Boussinesq flows**，*Computers & Fluids* 120 (2015) 70–81，DOI：`10.1016/j.compfluid.2015.07.012`。对应原文位于 `pdf/算法/`。

主程序结构、`commondata` 直接参数、算例宏、稳态误差和输出组织参照 `均匀网格/2DRBOpenmp.F90`；单 GPU 数据区和异步计算参照 `均匀网格/2DRBOpenacc.F90`。既有基线没有改动。DUGKS 是有限体积方法，因此原来的 LBM 碰撞—整格迁移被替换为论文的碰撞—特征线面重构—通量差推进，保持同一时刻的流场和温度耦合。

## 算例与直接参数

默认开启 `steadyFlow` 与 `SideHeatedCell`：左热右冷、上下绝热、四壁无滑移。RB 开启 `unsteadyFlow` 与 `RayleighBenardCell`：下热上冷、左右绝热、四壁无滑移。两组宏分别只能开启一个。只需在文件开头注释/启用相应宏，温度边界会跟随算例自动匹配。

RB 可选 `VerticalWallsPeriodical`，同时使速度和温度横向周期。默认保持项目的封闭 RB；论文第 3.2 节则采用横向周期、长高比 2。若要对应论文 RB，需同时选择周期宏和 `nx*dx=2*ny*dy`。不能把封闭方腔与该周期长方形的结果直接混作同一基准。

主要参数直接位于 `commondata`，带中文注释：

| 参数 | 默认值/用途 |
| --- | --- |
| `nx, ny` | `256,256`，单元中心均匀网格 |
| `Rayleigh, Prandtl, Mach` | `1e6, 0.71, 0.01` |
| `Thot, Tcold` | `0.5,-0.5`；浮力采用 `T-Tref` |
| `dx, dy` | `1,1`；壁面位于最外层单元的外侧半步位置 |
| `Courant` | `0.5`，时间步由离散速度最大模长确定 |
| `loadInitField` | `0` 新算，`1` 从最新检查点续算 |
| `steadyCheckIntervalItc` | 相隔 `2000` 个 DUGKS 步比较稳态误差 |
| `epsU, epsT` | `1e-7`；连续通过 3 次才报告稳态收敛 |
| `unsteadyRunDuration` | 非稳态总目标 `1000 t_ff`；续算只补足到总目标 |
| `unsteadyAverageStartTf/EndTf` | 非稳态时间平均窗口，默认后半段 |
| `outputSnapshotInterval` | 稳态 `10 t_ff`，非稳态 `0.5 t_ff` |

`cs2=1/3` 是与整数 D2Q9 速度绑定的格子常数，不作为独立可调声速。若修改物理尺度，保持速度集与声速一致。

## 论文方程与代码对应

两套分布均使用 D2Q9，平衡分布是论文式 (26)、(40) 的完整二阶速度展开；温度平衡没有替换成 D2Q5 或一阶展开。`equilibrium` 同时供两套分布使用。流场力项是式 (12) 的离散式 `F_i=((c_i-u)·a)/cs2 * f_eq_i`，浮力向上。

- `f/g` 保存式 (15)、(31) 的 **tilde 变换分布**，不是物理原始分布。
- `collision` 按式 (24)、(38) 计算 `bar+`，两套碰撞读取同一个 `t^n` 宏观场。式 (34) 中 `bar_g+` 的负号与式 (33)、(38) 不一致，代码采用由梯形积分推导得到的正号定义及式 (38)。
- `reconstruction_slopes` 使用内部中心斜率，壁面邻层使用单边斜率。`face_flux_x/y` 沿特征线追溯到 `x_face-h*c_i`，实施式 (19)–(21)、(35) 的二阶线性重构。切向粒子在实体边界也外推到物理壁面，不把壁面位置当作单元中心。
- 在界面先由 `bar_f/bar_g` 恢复密度、温度和含 `h*a/2` 的速度，再按式 (23)、(37) 恢复物理面分布。式 (23) 中 `tau_v*h*F/(2*tau_v+h)` 没有遗漏。
- `streaming` 以式 (13)、(30) 的共享面通量差更新单元平均值，使用式 (25)、(39) 的 `tilde+=4*bar+/3-tilde/3`。没有整格 streaming。
- `macro` 在两套分布同时推进后先求新温度，再按式 (46) 加 `dt*a(T_new)/2` 恢复物理速度。

采用 `nu=tau_v*cs2`、`kappa=tau_c*cs2`，**没有 LBM 的 `tau-1/2` 黏性修正**。`Uff=Mach*sqrt(cs2)`，`nu=Uff*L*sqrt(Pr/Ra)`，`gBeta=Uff^2/(L*DeltaT)`，`t_ff=L/Uff`。时间步为 `dt=Courant*min(dx,dy)/sqrt(2)`，所以物理时间为 `itc*dt`，所有输出、最大步数和采样均相应换算；不能将 DUGKS 步数直接当作格子时间。

无滑移边界在 **bar 分布** 上使用式 (42)，静止壁也保留 `u_wall'= -h*a/2`。代码从反弹关系自洽解出 `rho_wall`，再代入力修正，避免论文中 `rho_wall≈rhoRef` 的近不可压缩简化造成极小质量泄漏。恒温采用式 (43) 的反反弹，绝热采用式 (44) 的全方向反弹，然后恢复物理面分布。没有通过全局重归一化强行修正质量。

## Nu、Re、后处理

保留项目常用的 `NuVolAvg=1+L*<u_n*(T-Tref)>/(kappa*DeltaT)` 和 `ReVolAvg=L*sqrt(<u^2+v^2>)/nu`。侧热 `u_n=u`，RB `u_n=v`。用于文献速度比较的中心线输出默认采用 `UL/kappa` 标度。

`Nu_global`、冷热壁 `Nu`、中截面 `Nu_middle` 直接由实际 DUGKS 的 `sum(c_n*g_face)` 热通量计算，冷热壁采用同一个从热壁指向冷壁的正方向。全场数值热通量在对应方向对面值作梯形积分；没有用 LBM 分布或节点差分替代 DUGKS 的通量。

宏观 `NuVolAvg/ReVolAvg` 标在采样步 `t^n`；面通量诊断由这个状态重新构造，是 **`t^n+h` 的半步预测通量**。因此非稳态下两类诊断有半步时间差；壁通量文件和最终诊断中也注明这一点。重新构造不改变推进分布或采样场。

输出包括：

- `NuRe_2DOpenaccDUGKS.dat`：实际累计步数、真实 `t_ff`、体平均 Nu/Re、数值全场/冷热壁/中面 Nu、质量、最低温度。
- `Convergence_2DOpenaccDUGKS.dat`：相隔指定步数的场相对 L2 误差。速度为两分量联合范数，温度以 `T-Tref` 范数归一化。
- `FinalDiagnostics_2DOpenaccDUGKS.dat`、`WallNu*.dat`、`Centerline*.dat`：最终标量、壁热流及极值、中心线速度。极值位置是单元中心采样位置，未冒充基线五点抛物线/细样条插值的极值精度。
- `*Snapshot-*.bin`：有版本头、网格、累计步数、时间、坐标及 `u/v/T/rho` 的二进制快照。
- `*Tecplot-*.plt`：可直接读取的 ASCII Tecplot，包括速度、温度、密度、沿竖向积分的流函数和差分涡量；涡量内部采用中心差分，边缘单元采用单边差分。流函数采用项目封闭腔体的下壁零基准。
- `UnsteadyStatistics_2DOpenaccDUGKS.dat`：对实际时间戳作分段线性积分得到 Nu/Re 时间均值、时间标准差，以及窗口前后半段均值；窗口覆盖不足会标记 `INCOMPLETE_WINDOW`，不会据此声称统计收敛。

采样按每个绝对目标 `k*outputSnapshotInterval` 单独换算为整数步，避免先舍入采样步长再累计造成时间漂移。段末不在采样点的最终诊断单独输出，不追加进固定时间序列，因此改变续算分段不会改变历史或统计。

## 检查点续算

保留对应的编号 `.bin` 与 `reloadFile2DOpenaccDUGKS-latest.meta`，把 `loadInitField` 改为 `1`。检查点保存两套 **tilde 分布**、上次收敛比较参考场、误差、连续通过次数、输出编号及已有 Nu/Re 历史。临时面分布和斜率会从保存状态重新计算。

读取会检查版本头/尾、网格、算例、稳态/非稳态模式、周期选项、`Ra/Pr/Mach`、温标、`dt/CFL/cs2`、输出间隔等，并拒绝非有限分布及不完整文件。允许增加总目标步数/时长，不能把不同模型或不同时间步的检查点混用。Nu/Re 历史从检查点缓存重建；收敛曲线保留检查点以前的有效行，裁掉超前尾部。若旧收敛曲线缺失，会明确提示后从检查点继续记录，推进状态不受影响。二进制先完整关闭，再更新最新文件指针；若指针损坏可从保留的旧编号文件恢复指针，代码不会猜测一个文件继续运行。

## 编译与本地检查

单 GPU 编译示例（需要已安装 NVIDIA HPC SDK）：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccDUGKS.F90 -o 2DRBOpenaccDUGKS
```

代码也可以用 GNU Fortran 忽略 OpenACC 指令作本地小规模检查；这不等于验证 NVIDIA 设备执行。可重复本地验证脚本为 `verify_local.py`，例如：

```text
python verify_local.py --build-dir C:/msys64/tmp/2drb-dugks-check --compiler C:/msys64/ucrt64/bin/gfortran.exe
```

脚本只修改 ASCII 临时构建目录中的 16×16 测试副本。实际检查记录和被测主源 SHA-256 保存在 `validation_local.json`。已检查侧热/RB 与稳态/非稳态组合、周期 RB、双方向无浮力导热、带浮力质量守恒、37+43 步与 80 步连续推进的分布/历史一致性、错误参数续算拒绝，以及 GNU OpenACC 语法和 host 设备执行。

本地导热测试保留论文全反弹绝热边界；`NuVolAvg=1`、`Re=0`，数值面平均 Nu 误差约 `3.6e-6`，不将其称为逐分布精确导热不动点。带浮力 600 步测试相对质量偏差约 `1.3e-13`。这些是实现层面的检查，**没有进行服务器提交、nvfortran 编译、NVIDIA GPU 运行、论文基准复核、网格无关性或统计收敛验证**；默认 256×256 生产参数也没有长跑。
