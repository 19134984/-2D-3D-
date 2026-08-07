# D2Q9 TRT + D2Q9 BGK

主程序：`2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9BGK.F90`

- 流场：D2Q9 TRT，偶模态保留 LBM-CDE 的 `chi_s/chi_b` 修正。
- 温度场：原 LBM-CDE D2Q9 BGK，保留 `chi_kappa`、压力/力/梯度源项。
- 流场魔法参数固定使用 `chi_s` 修正后的有效偶尺度：
  `H_even,eff*H_odd=3/16`；无需额外编译开关。
