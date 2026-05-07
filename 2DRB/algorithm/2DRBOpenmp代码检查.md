# 2DRBOpenmp.F90 代码检查记录

检查对象：`2DRBOpenmp.F90`

检查依据：`algorithm/代码readme.md`

检查日期：2026-04-29

## 当前宏组合

当前文件顶部实际启用的组合是：

```text
flow mode      : unsteadyFlow
case           : RayleighBenardCell
thermal branch : EnableUseG
legacy branch  : off
velocity BC    : horizontal/vertical no-slip
thermal BC     : horizontal constant T, vertical adiabatic
```

对应代码位置：

- `steadyFlow` 被注释，`unsteadyFlow` 启用：`2DRBOpenmp.F90:11-12`
- RB 宏启用，side-heated 宏注释：`2DRBOpenmp.F90:28-39`
- `EnableUseG` 启用，`EnableLegacyThermalScheme` 注释：`2DRBOpenmp.F90:44-46`

结论：正式源文件当前应按“非稳态 RB + UseG 温度分支”理解；但本报告同时检查稳态和非稳态两条编译路径。稳态路径通过临时副本 `2DRBOpenmp_steady_check.F90` 将 `steadyFlow/unsteadyFlow` 对调后做 syntax-only 检查，检查完成后已删除临时副本。本轮按检查结果修正了正式 `2DRBOpenmp.F90` 中的几处问题：`calNuRe()` 的参考温度、reload 独立编号、OpenMP 文件中的无效 OpenACC 注释，以及快照/重启/plt 输出控制命名。

## 双模式检查范围

本轮按两层做检查：

```text
unsteadyFlow: 读取正式 2DRBOpenmp.F90 当前宏组合，重点检查非稳态时间推进、采样、NuVolAvg/ReVolAvg、reload/output cadence。
steadyFlow  : 使用临时副本切换宏，重点检查 steady-only 的 check()、最终 Nu/Re、壁面 Nu、中线 Nu、极值和流函数后处理是否能通过预处理与语法检查。
```

需要强调：稳态路径目前做的是代码阅读和 syntax-only 验证，还没有跑实际稳态算例到收敛。

## 1. 宏与参数

已检查：

- 宏互斥保护存在：
  - `steadyFlow` 和 `unsteadyFlow` 同开或同关会触发 `#error`。
  - `RayleighBenardCell` 和 `SideHeatedCell` 同开或同关会触发 `#error`。
  - `VerticalWallsNoslip` 和 `VerticalWallsPeriodicalU` 同开或同关会触发 `#error`。
  - `EnableUseG` 和 `EnableLegacyThermalScheme` 同开或同关会触发 `#error`。
- RB 下 `lengthUnit = ny`，side-heated 下 `lengthUnit = nx`，与主热流方向一致。
- 当前 `Thot=0.5`、`Tcold=-0.5`、`Tref=0`，因此当前配置下 `T` 与 `T-Tref` 数值一致。
- `heatFluxScale = lengthUnit/diffusivity`，与 `Nu = 1 + L/kappa * <convective flux>` 的体平均写法一致。
- UseG 下 `taug = 0.5 + (tauf-0.5)/Prandtl`，等价于 `0.5 + diffusivity/csT^2`，其中 D2Q5 `csT^2=1/3`。
- UseG 下 `Qk=1/taug`、`Qnu=1`、`thermalGeqCoeff=3`，符合 README 里的 D2Q5 UseG 检查基准。

需要注意：

- 当前正式文件默认并不是稳态，而是 `unsteadyFlow`。这不是算法错误，但后续做稳态 benchmark 时需要切换宏或使用独立稳态配置。
- `calNuRe()` 已改为使用 `T-Tref`，因此未来即使壁温中心不为 0，体平均对流热通量也保持和浮力/稳态后处理相同的参考温度约定。

## 2. 网格和初始化

已检查：

- half-way 网格坐标正确：

```text
xp(0)=0, xp(nx+1)=nx, xp(i)=i-0.5, xp=xp/lengthUnit
yp(0)=0, yp(ny+1)=ny, yp(j)=j-0.5, yp=yp/lengthUnit
```

- `nx/ny` 表示流体节点数，物理壁面位于 `0` 和 `N/lengthUnit`。
- RB 初始温度使用上下壁线性场，低 Ra 时才加小扰动；当前 `Rayleigh=1e7` 时扰动会跳过。
- UseG 历史通量 `Bx_prev/By_prev` 在非重启初始化和重启重构后都会设置为当前 `u*T`、`v*T`。

结论：网格位置与 README 的 half-way 边界解释一致。

## 3. 流场算法

已检查：

- D2Q9 速度集和权重是标准方向编号。
- `rho=sum_i f_i`，当前代码使用全密度，不是密度波动变量。
- D2Q9 MRT 矩、平衡矩和松弛参数结构符合常见 MRT-LBM：
  - 守恒矩 `rho, jx, jy` 的松弛率为 0。
  - 黏性相关矩使用 `Snu=1/tauf`。
  - 能量通量相关矩使用 `Sq=8*(2*tauf-1)/(8*tauf-1)`。
- 浮力项为 `rho*gBeta*(T-Tref)`，当前 RB 方向进入 `Fy`。
- Guo forcing 的半步力修正体现在 `fSource` 和 `macro()` 的 `+0.5*F` 速度恢复中。
- 主时间步顺序是：

```text
collision -> streaming -> bounceback -> macro
collisionT -> streamingT -> bouncebackT -> macroT
```

结论：流场主算法和 README 基准一致，未见 UseG/legacy 分支意外改变流场。

## 4. 温度场算法

已检查：

- 当前启用 UseG，不启用 legacy。
- UseG D2Q5 权重为 `omegaT(0)=1/3`、`omegaT(1:4)=1/6`。
- 温度矩：

```text
n = [T, x-flux, y-flux, high-order, anisotropic]^T
neq = [T, uT, vT, -2T/3, 0]^T
```

- `dBx = Bx(t)-Bx_prev`、`dBy = By(t)-By_prev`，修正项是 `(1-0.5*Qk)*dB`。
- `Bx_prev/By_prev` 在同一个格点碰撞中先读旧值、再写当前值；OpenMP 下每个线程只处理自己的格点，没有交叉写。
- 恒温 ABB 系数在 UseG 下通过 `2*omegaT(axis)` 得到 `1/3`，没有误用 legacy 的 `(4+aT)/10`。
- 绝热边界仍采用 half-way bounce-back。

结论：当前 D2Q5 UseG 温度场实现和 README 中的算法基准一致。

## 5. 后处理与 Nu

### 非稳态体平均 Nu/Re

当前 `unsteadyFlow` 中每 `0.5 t_ff` 调用一次 `calNuRe()`：

- side-heated 使用 `u*(T-Tref)`。
- RB 使用 `v*(T-Tref)`。
- 体平均 Nu 写为 `1 + lengthUnit/diffusivity * <convective flux>`。
- 当前 `Tref=0` 时数值上和旧写法等价，但新写法对非对称壁温更稳妥。

### 稳态 RB 后处理

稳态 RB 的全场和壁面 Nu 已按 side-heated 的方向旋转实现：

```text
q_y = v*theta - dtheta/dy
```

检查结果：

- `RBcalc_Nu_global()` 内部节点使用 `-dT/dy`，边界相邻节点使用 half-way 二阶公式。
- bottom/hot wall 使用正的上向热流：

```text
(8*Thot - 9*T(i,1) + T(i,2))/(3*dy)
```

- top/cold wall 输出为同一正向热流：

```text
(-8*Tcold + 9*T(i,ny) - T(i,ny-1))/(3*dy)
```

- 中线 Nu 对偶数 `ny` 使用 `jB/jT` 两层夹住 `y=1/2`，方向从 `x/u` 正确替换为 `y/v`。

结论：稳态 RB 后处理与 README 的“方向替换”原则一致。

### side-heated 后处理

side-heated 的 wall Nu 使用 half-way 壁面二阶公式：

```text
q_x(x=0) = (8*Thot - 9*T(1,j) + T(2,j))/(3*dx)
```

这和 README 中 Eq.(44) 符号笔误修正后的方向一致。

额外发现：README 之前写的第一层流体节点近壁公式 `4*theta(1)-theta(2)-3*theta_h` 与 half-way 几何不一致。代码使用的是 `4*theta_h-3*theta(1)-theta(2)` 对应的形式，按线性温度场检验是正确的。我已同步修正 README。

## 6. OpenMP 并行检查

已检查：

- 主要计算 kernel 使用 `default(none)`。
- `collision()` 的 `m/meq/m_post/s/fSource` 是私有数组。
- `collisionT()` 的 `n/neq/n_post/q` 和 `Bx/By/dBx/dBy` 是私有变量。
- `streaming()` 和 `streamingT()` 是 pull streaming，每个线程写自己的 `(i,j,alpha)`，未见目标写竞争。
- 边界 bounce-back 循环按单壁方向写入，角点处重复写入的分布函数对应同一个反弹值，未见不一致赋值。
- `calNuRe()`、稳态 Nu、Re 和残差计算使用 reduction。
- 极值位置没有在 OpenMP reduction 中直接更新，而是先得到局部 wall Nu 数组，再串行寻找极值和做五点拟合，避免了“max 值和位置不同步”的典型问题。
- 未发现 `nowait` 跨越真实依赖阶段。

数组访问模式：

- `u/v/T/rho` 以 `i` 为内层循环，符合 Fortran 第一维连续访问。
- `f/g` 方向维放在第一维，碰撞和宏观量恢复时读取一个格点的所有方向是连续的；OpenMP 下这是可接受布局。
- 若以后移植/对照 OpenACC，需要重新评估方向维在第一维对 GPU coalescing 的影响，但这不是当前 OpenMP 版本的正确性问题。

## 7. 输出、重启和非稳态控制

已检查：

- 非稳态采样时间点用目标 `dimensionlessTime+1` 和 `outputSnapshotInterval` 每次重新换算到 `itc`，避免累积漂移。
- `output_ReloadFile()` 只写 `f/g`，重启后通过 `reconstruct_macro_from_fg()` 重建 `rho/u/v/T/F/B_prev`，逻辑自洽。
- `output_SnapshotFile()` 是后处理快照，写出无量纲速度、`T` 和 `rho`，不用于严格重启。
- 输出命名已整理成三条线：
  - 快照：`outputSnapshotInterval`、`outputSnapshotIntervalItc`、`outputSnapshotFile`、`snapshotFileNum`、`snapshotFilePrefix`、`output_SnapshotFile()`。
  - 重启：`loadInitField`、`reloadDimensionlessTime`、`reloadFileNum`、`reloadFileInterval`、`reloadFileIntervalItc`、`outputReloadFile`、`reloadFilePrefix`、`output_ReloadFile()`。
  - Tecplot：`outputPltFileInterval`、`outputPltFileIntervalItc`、`outputPltFile`、`pltFileNum`、`pltFolderPrefix`、`output_Tecplot()`。
- `output_SnapshotFile()`、`output_ReloadFile()` 和 `output_Tecplot()` 的文件编号规则保持一致：
  - `steadyFlow`：文件名使用当前 `itc`。
  - `unsteadyFlow`：快照使用 `snapshotFileNum`，重启使用 `reloadFileNum`，Tecplot 使用 `pltFileNum`。
  - 文件编号统一按 `i12.12` 写成 12 位补零格式，重启读取文件名也按同样格式生成。
- 后处理数据文件保持双精度输出：
  - `convergence.plt` 和 `convergence_all.plt` 中的 `errorU/errorT` 使用 `ES24.16E3`。
  - `Nu_VolAvg.dat` 和 `Re_VolAvg.dat` 中的时间列与数值列使用 `ES24.16E3`。
  - 二进制 snapshot/reload/Tecplot 主场输出写入 `real(kind=8)` 数据。
  - `Nu_global` 的屏幕输出和 `settingsFile` 输出也使用 `ES24.16E3`。

稳态路径已检查：

- `check()` 只在 `steadyFlow` 下编译，使用 `up/vp/Tp` 计算速度和温度残差，并在同一循环末尾更新旧场。
- 稳态主循环按 `errorU/errorT` 和 `itc_max` 停止，不走非稳态固定时长逻辑。
- 稳态主循环中已按 `outputSnapshotFile/outputSnapshotIntervalItc` 增加可选周期快照输出。
- 稳态结束后会输出最终 Tecplot/二进制快照，并调用 RB 或 side-heated 对应的最终后处理。
- 当前 RB 稳态最终后处理包括 `RBcalc_Nu_global()`、`RBcalc_Nu_wall_avg()`、`RBcalc_umid_max()`、`RBcalc_vmid_max()`、`calc_psi_vort_and_output()` 和最终一次 `calNuRe()`。

已修正：

- `output_ReloadFile()` 的非稳态文件编号已改为独立的 `reloadFileNum`，不再依赖 `snapshotFileNum` 或快照输出是否开启；稳态文件名使用 `itc`。
- 快照输出已由 `outputFrequency/outputBinFile/binFileNum/binFolderPrefix/output_binary()` 语义整理为 `outputSnapshotInterval/outputSnapshotFile/snapshotFileNum/snapshotFilePrefix/output_SnapshotFile()`。
- `output_Tecplot()` 中的 `!$acc update self(u,v,T)` 已删除，改为普通注释说明 OpenMP 版本数据常驻 host。

需要注意：

- `calc_psi_vort_and_output()` 当前在 steady/unsteady 结束后都会执行。对于非稳态，这只是最终场诊断，不应解释为文献中的稳态后处理基准。

## 8. 本轮发现汇总

未发现的严重问题：

- 未发现 UseG 和 legacy 温度参数混用。
- 未发现 RB 后处理方向仍停留在 `x/u` 的问题。
- 未发现 OpenMP 主要 kernel 的明显数据竞争。
- 未发现 half-way 网格坐标与边界热流公式冲突。

需要确认或后续处理的点：

1. 当前正式文件启用 `unsteadyFlow`；稳态路径已通过临时副本做 syntax-only，但若要检查“稳态 RB 后处理运行结果”，仍需要切到 `steadyFlow` 后实际跑到收敛。
2. `calNuRe()` 已改为使用 `T-Tref`。
3. `output_ReloadFile()` 已改为使用独立的 `reloadFileNum`。
4. OpenMP 文件中的 `!$acc update self(u,v,T)` 已删除。
5. README 中第一层流体节点近壁热流公式已修正，代码本身使用的公式更符合 half-way 几何。

## 9. 验证

已执行 syntax-only 检查：

```text
powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile 2DRBOpenmp.F90 -SyntaxOnly
```

结果：通过。

稳态临时副本 syntax-only 检查：

```text
powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile 2DRBOpenmp_steady_check.F90 -SyntaxOnly
```

结果：通过。临时副本已删除。

说明：这验证了当前非稳态宏组合和临时稳态宏组合下的预处理和 Fortran/OpenMP 语法；没有做长时间数值运行，也没有验证 legacy 温度分支或 side-heated 宏组合。
