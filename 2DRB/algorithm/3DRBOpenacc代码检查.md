# 3DRBOpenacc.F90 代码检查记录

检查对象：`3DRBOpenacc.F90`

检查依据：`algorithm/代码readme.md`

检查日期：2026-04-29

## 当前宏组合

当前文件顶部实际启用的组合是：

```text
flow mode      : steadyFlow
case           : SideHeatedCell
thermal branch : EnableUseG
legacy branch  : off
velocity BC    : horizontal/vertical/spanwise no-slip
thermal BC     : vertical constant T, horizontal/spanwise adiabatic
```

对应代码位置：

- `steadyFlow` 启用，`unsteadyFlow` 注释：`3DRBOpenacc.F90:10-11`
- `SideHeatedCell` 启用，`RayleighBenardCell` 注释：`3DRBOpenacc.F90:31-46`
- `EnableUseG` 启用，`EnableLegacyThermalScheme` 注释：`3DRBOpenacc.F90:51-55`

结论：正式源文件当前应按“稳态 side-heated closed cavity + UseG 温度分支”理解；但本报告同时检查 `steadyFlow` 和 `unsteadyFlow` 两条编译/控制路径。`unsteadyFlow` 路径通过临时源码流对调 `steadyFlow/unsteadyFlow` 后做 syntax-only 检查，正式 `3DRBOpenacc.F90` 未改动。

## 双模式检查范围

本轮按两层做检查：

```text
steadyFlow  : 读取正式 3DRBOpenacc.F90 当前宏组合，重点检查残差停止、check3d()、最终 Nu/Re、壁面 Nu、中线 Nu、极值和最终输出。
unsteadyFlow: 临时切换宏后检查固定时长推进、采样点换算、NuVolAvg/ReVolAvg、bin/reload/output cadence，以及 OpenACC host/device 同步。
```

需要强调：两条路径目前做的是代码阅读和 syntax-only 验证，没有跑实际 GPU 算例，也没有跑稳态到收敛或非稳态到 `1000 t_ff`。

## 1. 宏与参数

已检查：

- `steadyFlow` 和 `unsteadyFlow` 有互斥保护；同开或同关会触发 `#error`。
- side-heated 下 `lengthUnit = nx`，RB 下 `lengthUnit = ny`，与主热流方向一致。
- 当前 `Thot=0.5`、`Tcold=-0.5`、`Tref=0`，因此当前配置下 `T` 与 `T-Tref` 数值一致。
- `tauf`、`viscosity`、`diffusivity`、`gBeta` 由 `Rayleigh/Prandtl/Mach/lengthUnit` 一起决定，符合 README 中要求成组核对的关系。
- UseG D3Q7 下 `csT^2=1/4`、`taug=0.5+diffusivity/csT^2`、`Qk=1/taug`、`Qnu=1`、`thermalGeqCoeff=4`。
- legacy 分支的 `Qk=3-sqrt(3)`、`Qnu=4*sqrt(3)-6` 和 `paraA` 只在 `EnableLegacyThermalScheme` 下生效，未混入当前 UseG 参数。

需要注意：

- 当前正式文件默认是 `steadyFlow`。检查非稳态时必须切换宏，不能直接拿默认编译路径代表非稳态。
- 如果未来把 `Thot/Tcold` 改成不关于 0 对称，`calNuRe3d()` 里体平均 Nu 是否应使用 `T-Tref` 而不是 `T` 需要重新确认。

## 2. 网格和初始化

已检查：

- half-way 网格坐标正确：

```text
xp(0)=0, xp(nx+1)=nx, xp(i)=i-0.5, xp=xp/lengthUnit
yp(0)=0, yp(ny+1)=ny, yp(j)=j-0.5, yp=yp/lengthUnit
zp(0)=0, zp(nz+1)=nz, zp(k)=k-0.5, zp=zp/lengthUnit
```

- `nx/ny/nz` 表示流体节点数，物理壁面位于 `0` 和 `N/lengthUnit`。
- side-heated 初始温度使用左右壁线性场；RB 路径使用上下壁线性场，低 Ra 时才加小扰动。
- 非重启初始化时，UseG 的 `Bx_prev/By_prev/Bz_prev` 设置为当前 `u*T/v*T/w*T`。
- 重启时由 `f/g` 重构 `T/rho/u/v/w`，并把 `Bx_prev/By_prev/Bz_prev` 重新设为当前 `u*T/v*T/w*T`。

风险点：

- 对 `EnableUseG` 来说，严格重启最好保存和恢复真实的 `Bx_prev/By_prev/Bz_prev`。当前 reload 文件只写 `f/g`，重启后把历史通量重置为当前值，因此第一步 `dB = B(t)-B(t-dt)` 不是严格连续的。

## 3. 流场算法

已检查：

- D3Q19 速度集、反向索引和权重结构一致。
- `rho=sum_i f_i`，当前代码使用全密度，不是密度波动变量。
- D3Q19 MRT 矩、平衡矩和松弛参数结构符合当前基准：
  - 守恒矩 `rho, jx, jy, jz` 的松弛率为 0。
  - 黏性相关矩使用 `Snu=1/tauf`。
  - 能量通量相关矩使用 `Sq=8*(2*tauf-1)/(8*tauf-1)`。
- 浮力项为 `rho*gBeta*(T-Tref)`，方向进入 `Fy`。
- Guo forcing 的半步力修正体现在矩空间 source 和 `macro3d()` 的 `+0.5*Fy` 速度恢复中。
- 主时间步顺序是：

```text
collision3d -> streaming3d -> bounceback3d -> macro3d
collisionT3d -> streamingT3d -> bouncebackT3d -> macroT3d
```

结论：流场主算法和 README 基准一致，未见 UseG/legacy 分支意外改变流场。

## 4. 温度场算法

已检查：

- 当前启用 UseG，不启用 legacy。
- UseG D3Q7 权重为 `omegaT(0)=1/4`、`omegaT(1:6)=1/8`。
- 温度矩结构为：

```text
n = [T, uT, vT, wT, high-order, anisotropic-1, anisotropic-2]^T
neq = [T, uT, vT, wT, -3T/4, 0, 0]^T
```

- `dBx/dBy/dBz` 定义为差值 `B(t)-B(t-dt)`，修正项是 `(1-0.5*Qk)*dB`，没有额外除以 `dt`。
- `wT` 位于第 3 个一阶通量矩，`dBz` 加到 `n_post(3)`，和 D3Q7 的 `ezT` 编号一致。
- `Bx_prev/By_prev/Bz_prev` 在同一个格点碰撞中先读旧值、再写当前值；OpenACC 下每个并行迭代只写自己的格点。
- 恒温 ABB 系数在 UseG 下通过 `2*omegaT(axis)` 得到 `1/4`，没有误用 legacy 的 `(6+aT)/21`。
- 绝热边界仍采用 half-way bounce-back。

结论：当前 D3Q7 UseG 温度场实现和 README 中的算法基准一致；三维新增的 `wT/dBz` 已对上。

## 5. 边界和 streaming

已检查：

- 流场和温度场都采用 pull streaming；每个格点从上游 `post` 缓冲读取，写回自己的目标格点。
- 当前速度边界为三个方向全无滑移，使用 half-way bounce-back。
- 当前温度边界为左右壁恒温、上下和前后绝热。
- 温度恒温边界按当前分支选择系数：
  - legacy: `(6+paraA)/21`
  - UseG: `2*omegaT(alpha)`，当前轴向方向为 `1/4`
- 温度绝热边界为 `g = g_post(oppT)`，符合 half-way bounce-back。

结论：当前 side-heated 配置下，边界方向、系数和 README 的 UseG D3Q7 检查表一致。

## 6. 后处理与 Nu

### 非稳态体平均 Nu/Re

`unsteadyFlow` 中采样逻辑为：

```text
nextSampleItc = round((dimensionlessTime+1) * outputFrequency * timeUnit)
call calNuRe3d()
```

已检查：

- 每个采样点都从目标 `dimensionlessTime+1` 重新换算到 `itc`，避免累积漂移。
- side-heated 使用 `u*T`，RB 使用 `v*T`。
- 体平均 Nu 写为 `1 + lengthUnit/diffusivity * <convective flux>`。
- Re 使用三维 RMS 速度并乘以 `lengthUnit/viscosity`。
- 采样后如果 `outputBinFile=1`，先 `update_host_all_3d_openacc()`，再输出 `u/v/w/T/rho` 快照。

需要注意：

- 当前 `Tref=0`，所以 `u*T` 与 `u*(T-Tref)` 等价；未来若壁温中心改变，需要重新确认。
- `unsteadyFlow` 里 `outputPltFile` 当前不参与循环内周期 Tecplot 输出，文件只在结束后强制输出一次。若以后希望非稳态周期 Tecplot，应补上对应分支。

### 稳态最终后处理

`steadyFlow` 中：

- 主循环用 `errorU/errorT` 和 `itc_max` 停止。
- `check3d()` 每 2000 步计算残差并更新 `up/vp/wp/Tp`。
- 结束后先 `update_host_all_3d_openacc()`，再输出最终 Tecplot/bin 快照。
- side-heated 路径调用 `SideHeatedcalc_Nu_global3d()`、`SideHeatedcalc_Nu_wall_avg3d()` 和中面速度极值。
- RB 路径调用 `RBcalc_Nu_global3d()`、`RBcalc_Nu_wall_avg3d()` 和中面速度极值。

已检查：

- side-heated 壁面 Nu 使用 half-way 二阶公式：

```text
q_x(x=0) = (8*Thot - 9*T(1,j,k) + T(2,j,k))/(3*dx)
```

- RB 壁面 Nu 是 side-heated 公式的方向替换，bottom/hot wall 使用正的上向热流：

```text
q_y(y=0) = (8*Thot - 9*T(i,1,k) + T(i,2,k))/(3*dy)
```

- Nu 分母使用 `DeltaT = Thot - Tcold`，没有把 `Tref` 当成温差。
- 稳态最终后处理在 CPU 上读数组前已经同步了 device 数据。

风险点：

- `steadyFlow` 结束后调用一次 `calNuRe3d()`，但写文件时间轴使用 `dimensionlessTime*outputFrequency`，不是实际 `itc/timeUnit`。因此稳态最终写入 `Nu_VolAvg_3D.dat` / `Re_VolAvg_3D.dat` 的时间可能不代表真实结束时刻。

## 7. OpenACC 并行检查

已检查：

- `enter_data_3d_openacc()` 把主场、常量表、`Bx_prev/By_prev/Bz_prev` 常驻 device，并创建 `f_post/g_post`。
- `update_host_all_3d_openacc()` 在文件输出、reload 写出和最终 CPU 后处理前同步 `u/v/w/T/rho/f/g/B*_prev`。
- `exit_data_3d_openacc()` 释放 steady-only 旧场和主场 device 映射。
- 主要 kernel 使用 `default(none)` 和显式 `present(...)`。
- `collision3d()` 的 `m/meq/m_post/s/fSource` 是 `private`。
- `collisionT3d()` 的 `n/neq/q/n_post/B/dB` 是 `private`。
- `check3d()`、`calNuRe3d()`、RB 全局 Nu、RB 壁面统计中的全局求和使用 `reduction` 或 `copyout`。
- 所有时间推进 kernel 放在同一个 `async(1)` 队列，主循环在检查/输出前有 `!$acc wait(1)`。
- 没有把普通 Fortran 文件 I/O 放进 GPU kernel。

需要注意：

- `SideHeatedcalc_Nu_global3d()` 和 `SideHeatedcalc_Nu_wall_avg3d()` 目前是 CPU 后处理，依赖主程序结束后的 `update_host_all_3d_openacc()`。这在当前调用顺序下是安全的；如果未来把这些子程序移到循环中调用，必须先同步或改写成 device 后处理。
- `gfortran -fopenacc -fsyntax-only` 只能说明预处理、Fortran 语法和部分 directive 写法可通过，不能证明 GPU runtime 正确。

## 8. 输出、重启和控制流

已检查：

- `steadyFlow` 和 `unsteadyFlow` 的输出频率、reload 频率、Tecplot 频率和 `itc_max` 分别放在宏分支中。
- `unsteadyFlow` 的 `itc_max` 由 `unsteadyRunDuration*timeUnit` 换算。
- `output_binary3d()` 是后处理快照，输出无量纲速度、`T` 和 `rho`，不用于严格重启。
- `writeReloadFile3d()` 当前只写 `f/g`，重启后由 `reconstruct_macro_from_fg3d()` 恢复宏观量。

需要注意：

- `writeReloadFile3d()` 在 `unsteadyFlow` 下用 `binFileNum` 命名。如果未来设置 `outputBinFile=0` 但 `outputReloadFile=1`，可能出现 reload 文件名不递增或覆盖风险。当前默认 `outputBinFile=1`，不会触发。
- `EnableUseG` 下严格重启缺少 `B*_prev`，这一点比普通宏观量重构更重要，建议优先修正。

## 9. 本轮发现汇总

未发现的严重问题：

- 未发现 UseG 和 legacy 温度参数混用。
- 未发现 D3Q7 的 `wT/dBz/ezT/oppT` 编号错位。
- 未发现流场浮力项误用 `DeltaT` 或漏掉半步力修正。
- 未发现当前 OpenACC 主 kernel 的明显 shared scratch 或 reduction 漏项。
- 未发现 half-way 网格坐标与壁面 Nu 公式冲突。
- 稳态和非稳态两条宏路径均通过 syntax-only 检查。

需要修正或确认的点：

1. `EnableUseG` 严格重启应保存并恢复 `Bx_prev/By_prev/Bz_prev`，否则重启后第一步历史通量不连续。
2. `steadyFlow` 最终一次 `calNuRe3d()` 写出的时间轴不是真实 `itc/timeUnit`。
3. `unsteadyFlow` 的 reload 文件名依赖 `binFileNum`；若关闭 bin 输出但开启 reload，存在命名风险。
4. `unsteadyFlow` 当前只在结束时强制输出 Tecplot；如果希望周期 Tecplot，`outputPltFile` 分支需要补齐。
5. GPU runtime 正确性仍需要用 OpenACC 编译器和同宏组 CPU/OpenMP 基准对比 `errorU/errorT`、`Nu/Re`、壁面 Nu、中线 Nu 和最终场。

## 10. 验证

已执行默认稳态路径 syntax-only 检查：

```text
powershell -ExecutionPolicy Bypass -File .\build_in_ascii_path.ps1 -SourceFile 3DRBOpenacc.F90 -SyntaxOnly -ParallelFlag -fopenacc
```

结果：通过。

已执行临时非稳态路径 syntax-only 检查：

```text
Get-Content 3DRBOpenacc.F90
临时把 #define steadyFlow 改为 !#define steadyFlow
临时把 !#define unsteadyFlow 改为 #define unsteadyFlow
通过 stdin 输入 gfortran -x f95-cpp-input -cpp -fopenacc -fsyntax-only -
```

结果：通过。

说明：以上验证只覆盖预处理、Fortran 语法和部分 OpenACC directive 语法；本机未检测到 `nvfortran`，因此不能视为 GPU runtime 验证。
