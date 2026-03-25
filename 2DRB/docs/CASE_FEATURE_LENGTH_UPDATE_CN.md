# 2DRB 特征长度自动切换修改说明

## 修改目的

本次修改的目标是把 `2DRB` 相关代码中的特征长度 `lengthUnit` 改成按当前算例自动切换，避免在不同 case 之间手工修改 `Nx` 或 `Ny`。

当前约定如下：

- 对于左右壁面差温的 `SideHeatedCell`，特征长度取 `Nx`
- 对于上下壁面差温的 `RayleighBenardCell`，特征长度取 `Ny`

也就是说，整套无量纲长度仍然只使用一个统一的特征长度，但这个特征长度会随着 case 自动变化。

## 主要修改内容

### 1. 主程序 `2DRB.F90`

在 `commondata` 中将 `lengthUnit` 改成按宏自动选择：

- `SideHeatedCell` 时使用 `dble(nx)`
- 其他情况使用 `dble(ny)`

同时把坐标构造和步长定义恢复为统一按 `lengthUnit` 缩放：

- `xp(i) = (i - 0.5) / lengthUnit`
- `yp(j) = (j - 0.5) / lengthUnit`
- `dx = 1.0d0 / lengthUnit`
- `dy = 1.0d0 / lengthUnit`

### 2. 边界位置不再写死为 `1.0`

由于 `lengthUnit` 在侧壁差温时改成了 `Nx`，这时：

- 右壁位置仍然是 `xp(nx+1) = nx / lengthUnit = 1`
- 顶壁位置则变成 `yp(ny+1) = ny / lengthUnit = ny / nx`

因此原来代码中把顶壁或右壁写死成 `1.0d0` 的地方，在这次修改中统一改成使用真实边界坐标：

- 右壁使用 `xp(nx+1)`
- 顶壁使用 `yp(ny+1)`

这样在 `Nx /= Ny` 时，边界位置与当前无量纲定义仍然保持一致。

### 3. `2DRB_ISLBM.F90`

对 ISLBM 版本做了同样的自动切换处理，但由于该版本支持非均匀网格，修改比均匀网格版本多一步：

- `SideHeatedCell` 时，通过首层 `x` 向中心坐标反推 `lengthUnit`
- `RayleighBenardCell` 时，通过首层 `y` 向中心坐标反推 `lengthUnit`

随后增加了统一的坐标重缩放步骤：

- `xp = xp * (nx / lengthUnit)`
- `yp = yp * (ny / lengthUnit)`

这样非均匀网格生成逻辑保持不变，但最终输出到求解和后处理的坐标，始终与当前 case 采用的特征长度一致。

此外还同步修正了以下依赖边界坐标的逻辑：

- 有效长度 `Lx_eff` / `Ly_eff`
- 对称性检查
- 右壁、顶壁一侧的单边导数距离
- 样条扩展末端坐标 `xExt(nXExt)` / `yExt(nYExt)`

### 4. 参考代码与诊断代码同步

以下侧壁差温参考/诊断代码同步改成了 `lengthUnit = dble(nx)`：

- `2DRB/ChaiPRE2020.F90`
- `2DRB/diagnostics/ChaiPRE2020_41.F90`
- `2DRB/diagnostics/ChaiPRE2020_81.F90`
- `2DRB/diagnostics/ChaiPRE2020_sweep.F90`

同时也把这些文件里写死为 `1.0d0` 的顶壁坐标改成了 `yp(ny+1)`。

## 本次修改的直接效果

- 切换 `SideHeatedCell` 与 `RayleighBenardCell` 时，不再需要手工改特征长度
- `Nx /= Ny` 时，侧壁差温和上下差温都会自动使用各自正确的特征长度
- 与边界位置相关的插值、拟合、导数计算不再依赖固定的 `1.0d0`

## 编译检查说明

本次修改后：

- `2DRB_ISLBM.F90` 已通过最小编译检查
- `2DRB.F90` 仍存在一个原有的编译问题：`SG = 1.0d0 - 0.5d0*s_j` 处把非常量写成了 `parameter`，该问题不是本次修改引入的
- `ChaiPRE2020.F90` 本身也存在原有的独立编译问题，本次修改未额外扩大该问题范围

## 涉及文件

- `2DRB/2DRB.F90`
- `2DRB/2DRB_ISLBM.F90`
- `2DRB/ChaiPRE2020.F90`
- `2DRB/diagnostics/ChaiPRE2020_41.F90`
- `2DRB/diagnostics/ChaiPRE2020_81.F90`
- `2DRB/diagnostics/ChaiPRE2020_sweep.F90`
- `2DRB/docs/CASE_FEATURE_LENGTH_UPDATE_CN.md`
