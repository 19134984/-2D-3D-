# D2Q9-TRT + D2Q9-TRT-a-chi_kappa LBM-CDE

主程序：`2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90`

- 程序框架来自当前二维 OpenACC 基线，流场与温度场均采用 D2Q9。
- 流场使用 TRT LBM-CDE：手动设置 `chi_nu`，由目标运动黏度换算基础松弛时间 `tauf`。
- 流场奇模态可在源码顶部切换 original/effective 两种 magic 尺度，默认使用 original magic。
- 温度场使用 D2Q9-TRT-`a`-`chi_kappa` LBM-CDE。与流场的 `chi_nu` 一样，
  `chi_kappa` 是可手动调节的稳定性参数；代码在保持目标热扩散率不变的条件下反算 `ce`，
  再由 `paraA=6*ce-4` 得到 `paraA`。`ce` 与流场声速平方
  `cs2=1/3` 相互独立，压力仍严格使用 `p=cs2*(rho-rho0)`。
- 温度 TRT 固定松弛率为：

  ```text
  Qe  = 4*sqrt(3)-6     ! nT(1:2)：能量型矩 e、epsilon
  Qk  = 3-sqrt(3)       ! nT(3:6)：温度通量矩 jx、qx、jy、qy
  Qnu = 4*sqrt(3)-6     ! nT(7:8)：二阶各向异性矩 pxx、pxy
  ```

  当前固定 `Qe=Qnu`；程序仍在矩空间分别给能量矩和二阶各向异性矩赋值，便于以后独立调节。

- D2Q9 参数满足 `paraA=6*ce-4`，并决定基础平衡态权重：

  ```text
  W0          = -(5*paraA+2)/18
  W1:4        = (paraA+4)/18
  W5:8        = (paraA+4)/72
  ```

  这些表达式分别等价于 `1-5*ce/3`、`ce/3` 和 `ce/12`。

  参数链为

  ```text
  ce = diffusivity/[(1-chi_kappa)*(1/Qk-1/2)*dt]
  paraA = 6*ce-4
  diffusivity = (1-chi_kappa)*ce*(1/Qk-1/2)*dt
  ```

  同时要求 `chi_kappa>=0` 和基础权重严格为正时，必须满足

  ```text
  0 <= chi_kappa < 1-5*diffusivity/[3*(1/Qk-1/2)*dt]
  diffusivity/[(1/Qk-1/2)*dt] <= ce < 3/5
  ```

  当前固定 `Qk=3-sqrt(3)`，即 `1/Qk-1/2=1/sqrt(12)`，且格子时间步 `dt=1`，
  因此 `chi_kappa` 上限为

  ```text
  chi_kappa < 1-10*sqrt(3)*diffusivity/3
  sqrt(12)*diffusivity <= ce < 3/5
  ```

  上限是严格不等号；若反算得到 `ce>=3/5`，静止方向基础权重将不再严格为正，
  因此代码会在启动阶段停止。`ce=1/3`、`paraA=-2` 时恢复标准 D2Q9 基础权重；同时静态平衡矩
  `nT1_eq=-2*T`、`nT2_eq=T`。对流 Hermite 项仍使用标准权重 `omega` 和 `cs2=1/3`。
- 目标热扩散率始终采用 `diffusivity=viscosity/Prandtl`；调节 `chi_kappa` 时，代码通过改变
  `ce`、`paraA` 和基础权重保持该目标热扩散率不变。

- 温度平衡态由 `W_i(paraA)` 基础项、标准二阶 Hermite 速度项和 `lambda_i` 压力项组成。
  在代码采用的矩顺序中，两个静态自由矩严格满足

  ```text
  nT1_eq = [paraA+3*(u^2+v^2)+6*p/rho0]*T
  nT2_eq = [-(3*paraA+4)/2-3*(u^2+v^2)-9*p/rho0]*T
  ```

  因此 `nT2_eq` 的常数项不是另一个可独立调节的参数。
- 温度梯度由该平衡态的非平衡一阶矩重构；温度源项同步包含
  `(pressure/rho0+chi_kappa*ce)*grad(T)+T*F/rho0`，并作为纯奇源项使用 `1-Qk/2`。
- 碰撞采用与流场一致的正交 D2Q9 矩阵，矩顺序为
  `(T,e,epsilon,jx,qx,jy,qy,pxx,pxy)`；对应松弛率为
  `nT(0)->0`、`nT(1:2)->Qe`、`nT(3:6)->Qk`、`nT(7:8)->Qnu`。
- 恒温壁采用包含 `2*W_i(paraA)*Twall` 和动态压力项的 anti-bounce-back，绝热壁采用 bounce-back；D2Q9 对角方向和角点均显式处理。
- 初始化时一次性检查权重严格为正、权重和与二阶矩、九个平衡矩、源项矩、固定 TRT 参数、
  `paraA=-2` 标准极限、`chi_kappa->ce->paraA` 参数链以及热扩散率重构误差。
- 非稳态输出包含 Nu/Re、原始动能与热耗散、RB 精确耗散关系、前后半统计窗口比较、温度 RMS 边界层和网格/时间分辨率指标。

## 快速回归测试

在仓库根目录运行：

```text
python Xs/D2Q9TRT_D2Q9TRT/test_d2q9_trt_ce_chikappa.py
```

测试会检查 `paraA=-2` 标准 D2Q9 极限、九个平衡矩、纯奇源项矩、矩空间 TRT 及其与等价奇偶形式的一致性、
`chi_kappa` 合法范围与热扩散率重构、三个合法 `chi_kappa` 下保持同一目标热扩散率的周期扩散，
以及 halfway-ABB 一维稳态导热。

## 验证边界

固定 TRT 松弛率只对应特定体相误差条件，不能自动推出包含压力、浮力、非稳态反馈和 halfway BB/ABB 后的全局高阶精度。正式使用仍需进行 OpenACC 编译、同参数数值对照、耗散关系检查和网格收敛验证。
