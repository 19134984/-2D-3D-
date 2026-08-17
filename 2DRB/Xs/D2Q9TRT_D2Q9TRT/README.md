# D2Q9-TRT + D2Q9-TRT LBM-CDE

主程序：`2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90`

- 程序框架来自当前二维 OpenACC 基线，流场与温度场均采用 D2Q9。
- 流场使用 TRT LBM-CDE：手动设置 `chi_nu`，由目标运动黏度换算基础松弛时间 `tauf`。
- 流场奇模态可在源码顶部切换 original/effective 两种 magic 尺度，默认使用 original magic。
- 温度场使用 D2Q9-TRT-a-`chi_kappa` LBM-CDE，固定松弛率为：

  ```text
  Qk  = 3-sqrt(3)       ! 奇模态/温度通量
  Qnu = 4*sqrt(3)-6     ! 偶模态
  ```

- `paraA` 是手动设置的温度平衡态/稳定性设计参数：

  ```text
  ce          = (4+paraA)/10
  W0          = (2-paraA)/6
  W1:4        = (4+paraA)/30
  W5:8        = (4+paraA)/120
  ```

  标准 D2Q9 Hermite 权重 `omega` 和 `cs2=1/3` 不随 `paraA` 改变。
- 目标热扩散率保持为现有物理参数给出的 `diffusivity=viscosity/Prandtl`，代码自动反算

  ```text
  chi_kappa = 1-20*sqrt(3)*diffusivity/(4+paraA)
  diffusivityReconstructed = (1-chi_kappa)*ce*(1/Qk-1/2)
  ```

- 温度梯度由新平衡态的非平衡一阶矩重构；温度源项同步包含
  `(pressure/rho0+chi_kappa*ce)*grad(T)+T*F/rho0`，并作为纯奇源项使用 `1-Qk/2`。
- 恒温壁采用包含 `2*W_i(paraA)*Twall` 和动态压力项的 anti-bounce-back，绝热壁采用 bounce-back；D2Q9 对角方向和角点均显式处理。
- 初始化时一次性检查 `W_i` 归一性、平衡态零/一/二阶矩及热扩散率重构误差。
- 非稳态输出包含 Nu/Re、原始动能与热耗散、RB 精确耗散关系、前后半统计窗口比较、温度 RMS 边界层和网格/时间分辨率指标。

## 验证边界

固定 TRT 松弛率只对应特定体相误差条件，不能自动推出包含压力、浮力、非稳态反馈和 halfway BB/ABB 后的全局高阶精度。正式使用仍需进行 OpenACC 编译、同参数数值对照、耗散关系检查和网格收敛验证。
