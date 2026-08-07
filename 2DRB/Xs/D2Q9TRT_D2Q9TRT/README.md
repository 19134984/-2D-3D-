# D2Q9 TRT + D2Q9 TRT LBM-CDE

主程序：`2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9TRT.F90`

- 流场：D2Q9 TRT，保留 LBM-CDE 的 `chi_s/chi_b` 偶应力修正。
- 流场奇模态：固定按 `chi_s` 修正后的有效偶尺度满足
  `H_even,eff*H_odd=3/16`；无需额外编译开关。
- 温度场：保留原 LBM-CDE D2Q9 平衡态、`chi_kappa`、压力/力/梯度源项、
  局部温度梯度闭合和 D2Q9 BB/ABB 边界。
- 温度碰撞：由 BGK 改为速度反演奇偶 TRT；每一部分分别使用
  `(1-s_even/2)` 和 `(1-s_odd/2)` 源项前因子。
- 固定名义温度率：

  ```text
  sigma_odd  = 1/sqrt(12),  s_odd  = 3-sqrt(3)
  sigma_even = 1/sqrt(3),   s_even = 4*sqrt(3)-6
  ```

  它们满足冻结、无源、`p=0` 的 D2Q9 体相四阶完全消除条件。目标热扩散率
  通过

  ```text
  kappa = (1-chi_kappa)*cT2*sigma_odd
  ```

  自动反算 `chi_kappa`，因此该变体禁止 `CHI_KAPPA_OVERRIDE` 和
  `BASE_TAUG_OVERRIDE`。

## 适用边界

固定名义率的四阶结论不能外推为完整高 Rayleigh 数自然对流的全局四阶精度。
加入 `chi_kappa` 局部反馈、空间变化压力、浮力、非稳态项以及 halfway BB/ABB
后，体相和壁面仍需分别验证。该文件是一个明确的数值对照分支，不是已经完成
网格收敛或 GPU 基准验证的生产结论。
