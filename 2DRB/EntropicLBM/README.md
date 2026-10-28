# 二维均匀网格熵稳定双分布热 LBM

主文件为 `2DRBOpenaccEntropic.F90`。文献是 Frapolli、Chikatamarla、Karlin，*Entropic lattice Boltzmann simulation of thermal convective turbulence*，Computers & Fluids 175 (2018) 2–19，DOI：10.1016/j.compfluid.2018.08.021。对应仓库 `pdf/算法/` 中的同名 PDF。

这是该文热双分布 KBC 模型的**二维扩展实现**，以 `均匀网格/2DRBOpenmp.F90` 的算例组织、推进和诊断方式为主框架，参照同维 OpenACC 版本的设备常驻方式；检查点、快照及后处理辅助代码沿用现有 `FPC/2DRBOpenaccFPC.F90` 的项目风格。原文件均未修改。当前为单 GPU OpenACC、二维均匀网格、封闭腔体，未加入多块网格、旋转或三维求解。

## 直接修改的参数和开关

文件开头两组选项各启用一项：`steadyFlow/unsteadyFlow` 和 `SideHeatedCell/RayleighBenardCell`。默认 `SideHeatedCell + steadyFlow`；RB 通常改为 `RayleighBenardCell + unsteadyFlow`。两种算例均可编译为稳态或非稳态模式。

`commondata` 中可直接修改 `nx, ny, Rayleigh, Prandtl, Mach, Thot, Tcold`、稳态阈值/检查周期、非稳态运行时长/平均窗口和输出间隔；均有中文说明。默认 `256×256, Ra=1e6, Pr=0.71, Uff=0.01`，其中 `Mach=sqrt(3)*Uff`。温度使用 `Tcold=T0=1/3, Thot=T0+0.01`，符合本文固定参考压力和小温差前提；不能直接把其他被动温度模型的任意温标搬过来而认为仍是同一热模型。`cs2=thermalCs2=1/3` 是固定 D2Q9 格子常数。

侧热的特征长度为 `nx`，RB 为 `ny`；`gBeta=Ra*nu*kappa/(L^3*DeltaT)`，`t_ff=L/Uff`。黏度和热扩散率按 `nu=Uff*L*sqrt(Pr/Ra)`、`kappa=nu/Pr` 耦合确定。

## 方程与代码的对应

将论文三维速度张量积降为二维 D2Q9，采用同一 D1Q3 的 `{-1,0,1}` 和权重 `{1/6,2/3,1/6}`，设 `D=2, cv=D/2=1`。方向排列与本仓库其他二维程序相同。

| 论文部分 | 实现位置和含义 |
|---|---|
| 式 (2.1)–(2.3) | `macro`：`rho=sum(f)`，`u_raw=sum(c*f)/rho`，`T=(sum(g)/rho-|u_raw|^2)/D` |
| 式 (2.9)–(2.12) | `entropic_moment_population`、`entropic_equilibria`：用零、一、二阶矩生成二阶 Hermite 平衡态，保留能量通量和 `R_eq` |
| Table 1，`Pr<=1` | `fstar=feq`；`gstar` 一阶矩为 `q-2*u·(P-Peq)`，二阶矩为 `R_eq` |
| Table 1，`Pr>1` | `fstar` 保留完整压力张量 `P`；`gstar` 一阶矩为 `qeq+2*u·(P-Peq)`，二阶矩为 `R_eq` |
| 式 (2.13)–(2.14) | 自然原始矩顺序为 `00,10,01,20,02,11,21,12,22`；`k` 承载前三个守恒矩，`s` 仅承载 `20,02,11`，`h` 承载 `21,12,22` |
| 式 (2.4)–(2.5) | `collision/collisionT`：KBC 流场碰撞及准平衡双松弛总能量碰撞，读取同一时刻状态 |
| 式 (2.6)、第 2.3 节的极值定义 | `entropic_relax`：沿实际碰撞路径极小化 `H=sum(f*log(f/W))` |
| 式 (2.17)–(2.18) | `collision`：只对 `f` 施加 EDM 浮力，`a_y=gBeta*(T-Tref)` |
| 式 (2.22)–(2.24) | `beta=1/(2*min(tauf,taug))`，`beta1=1/(2*max(tauf,taug))`，分别恢复两侧 Pr 分支的 `nu/kappa` |

这里的 `g` 是**总能量分布**：`sum(g)=2*rho*E=D*rho*T+rho*|u_raw|²`，不是 `sum(g)=T` 的被动温度分布。没有叠加旧 D2Q5、`useG` 或 FPC 温度碰撞。仅 `f` 要求正性；`g` 不要求逐方向为正，例如 `T=T0,u=0` 时静止方向的 `g_eq` 可以为零。

## 原文内部不一致与本实现的选择

这是本次与用户确认后采用的处理：**按熵极小值定义计算，并记录与原文印刷公式的差异**，不能把本程序表述为逐字照抄式 (2.15)/(2.16)。尚无作者勘误作为依据。

原文式 (2.15) 打印的高阶部分是 `gamma_p*h+(1-gamma_p)*heq`，而通常 KBC 松弛记号写成 `(1-gamma)*h+gamma*heq`，两者须取 `gamma=1-gamma_p` 才表示同一路径。本代码使用后一个记号，实际碰撞为：

```text
fpost = feq + (1-2*beta)*ds + (1-beta*gamma)*dh
        + 2*(beta-beta1)*(fstar-feq)
```

从这条实际路径推导二次熵驻点得到：

```text
gamma = 1/beta - (2-1/beta)*<ds|dh>/<dh|dh>
        + 2*(beta-beta1)/beta * <fstar-feq|dh>/<dh|dh>
<a|b> = sum(a*b/feq)
```

这与论文式 (2.16) 的 `beta1` 分母及第三项符号不能直接对上；`Pr<=1` 时最后一项消失，仍存在 `beta/beta1` 的不一致。本实现默认 `entropyGammaStrategy=2`，绕过这一近似式，直接对上面的实际 `fpost` 求 `H` 极小值。`strategy=1` 仅供比较上述重新推导的二次熵近似，也不等于论文印刷式 (2.16)。

算法先求 `fpost>0` 的可行区间，再用有界 Newton/二分法；终止由 `entropyTolerance` 控制。若无正分布可行区、求解未达到容差或 EDM 后出现非正 `f`，持久记录 `collisionStatus`，在有限性检查/采样/终止时明确失败，不靠裁剪分布掩盖失稳。熵极小值是在固定其他碰撞分量的路径上求解，不代表在任意粗网格、任意 Mach 下都能保证稳定。

平衡分布沿用论文式 (2.12) 的二阶多项式近似，并非完整的熵极小平衡分布。因此，上述路径求极小值不构成离散 H 定理或全局熵单调性的证明。

`dh` 从独立的高阶自然矩差反演，禁止用 `pop-feq-ds` 混入低阶舍入残量；熵求解使用归一化方向及实际人口修正振幅，避免近平衡时巨大 `gamma` 放大舍入误差。此项已用周期衰减和近平衡守恒回归检查。

## 力、温度和输出速度

数值碰撞始终从 `f` 恢复论文原始动量速度 `u_raw`。温度也按该速度从总能量扣除动能。输出及 Nu/Re 诊断的 `u/v` 使用 EDM 的物理半步速度 `u_phys=u_raw+a/2`；两种速度没有在碰撞中混用。初始 `f` 为 `u_raw=0` 的平衡态，因此若初始浮力非零，最初的诊断速度为 `a/2`。

按论文明确给出的式 (2.17)，EDM 只作用于 `f`，没有自行增加 `g` 的体力功源项。因此均匀受力、无能量输运时，`sum(g)` 守恒而 `u_raw` 改变，会有 `T_new-T_old=-(2*u_raw·a+|a|²)/D`。这是需要保留并解释的低 Mach 总能量约定；不能声称强迫过程中温度严格不变。论文第 2.5 节也以压缩功及黏性耗散在所考察 Boussinesq 情形下可忽略为前提。

## 边界、诊断与输出

默认四壁无滑移。侧热：左热右冷、上下绝热；RB：下热上冷、左右绝热。流场用项目半程 bounce-back；恒温边界对**总能量平衡分布**做半程 anti-bounce-back，壁面密度用相邻流体密度，壁速为零；绝热边界反射法向离散速度、保留切向来源。角点恒温条件优先。原文采用三维、水平周期和 Grad 壁面处理，本程序的二维封闭几何和半程壁面是对用户现有算例风格的适配，不能直接当作原文三维湍流表格的同一验证算例。生产源当前不提供周期几何切换；周期条件仅用于独立输运系数检查。

保留项目的稳态 `errorU/errorT`、Nu 全场/冷热壁/中线及极值、中心线速度极值、流函数/涡量、二进制宏观快照、Tecplot、Nu/Re 历史和非稳态窗口统计。`NuVolAvg=1+L*<u_heat*(T-Tref)>/(kappa*DeltaT)`；`ReVolAvg=L*sqrt(<u_phys²+v_phys²>)/nu`。所有 Nu 的对流与导热项均按真实 `DeltaT` 归一化，Nu 的热扩散标度与可编辑的速度显示标度 `velocityScaleCompare` 分离。稳态温度误差为 `sum(abs(T-Tprev))/(nx*ny*DeltaT)`，不受温标平移影响。

`itc_max` 到达与稳态收敛分别记录；非稳态运行结束、窗口样本覆盖、统计收敛和网格无关性也分别看待。平均窗口必须满足 `0<=start<mid<end<=duration`；若实际历史没有覆盖窗口，输出 `INCOMPLETE` 并不生成名义完整窗口平均值。

`loadInitField=1` 从同工作目录的 `reloadFile2DOpenaccEntropic-latest.meta` 续算，需配套保存最新 `.bin`、Nu/Re 历史和稳态收敛历史。V2 检查点含 `f/g`，稳态还含上一检查时刻的 `up/vp/Tp` 及 `errorU/errorT`，恢复绝对检查相位；网格、算例、流动模式、热/流体参数、采样间隔和稳态检查间隔必须匹配，目标总时长及误差阈值可以调整。历史超前或缺失会拒绝续算。GNU/NVHPC 的顺序非格式化文件不承诺跨编译器兼容，旧模型/旧版本检查点也不兼容。

## 本地检查及当前验证边界

单 GPU 编译示例（需要 NVIDIA HPC SDK；本次未执行该编译器）：

```text
nvfortran -O3 -acc -Minfo=accel -Mpreprocess -Mextend 2DRBOpenaccEntropic.F90 -o 2DRBOpenaccEntropic
```

在本目录用 Python 运行 `verify_local.py`；可传入 GNU Fortran 编译器绝对路径。它不改生产源，变体和输出保存在 `_local_checks/`，摘要为 `validation_local.json`。检查涵盖四种算例/流动宏组合的 GNU OpenACC 编译、两个主要组合的 16×16 短运行、Pr 两侧矩/熵/能量契约、近平衡守恒、连续与分段推进的逐字节检查点/历史一致性、非检查相位的稳态续算，以及不完整平均窗口的明确标记。GNU OpenACC 运行设备为 host。

`verify_transport.py` 提取实际生产碰撞核，在独立周期驱动上测试 `Pr=0.7,1,7` 的剪切波/温度波衰减；测试专门覆盖为 `64×1, nu=0.05`、关闭浮力，各运行 5000 步，不改变生产网格/边界。结果保存于 `transport_validation_local.json`：本次测得黏度相对偏差约 `0.0394%`，热扩散率相对偏差约 `0.0731%–0.0802%`，属于有限波数离散误差量级。摘要中保存被检验源文件的 SHA-256。

`verify_conduction.py` 用实际初始化、OpenACC 推进、边界和 Nu 例程检查双向纯导热：测试专门覆盖为 `16×12`、关闭浮力及 RB 初始扰动，各推进 600 步，结果为 `conduction_validation_local.json`，临时输出位于 `_local_checks/conduction/`。侧热与 RB 的温度最大绝对误差分别约 `1.55e-15` 和 `2.33e-15`，速度保持零，总质量最大绝对误差约 `2.84e-14`；全域、冷热壁及中线 Nu 与 1 的最大偏差小于 `1.8e-12`。这检验半程恒温/绝热壁和 Nu 归一化的一致性，不构成有浮力对流基准验证。

已检查：局部矩与守恒契约、熵路径、线性输运系数、双向纯导热边界、本地主机短运行、输出及续算。尚未检查：NVHPC/P100 实际 GPU 编译与运行、封闭侧热稳态基准精度、RB 非稳态长期统计、二维网格无关性，以及原文三维湍流结果。当前没有证据将剩余误差归因为并行装配；后续基准若不符，需要区分二维/壁面适配、物理模型假设与后处理误差。本次没有登录服务器或提交任务。
