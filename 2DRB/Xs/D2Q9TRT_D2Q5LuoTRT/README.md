# D2Q9 TRT + Luo/Wang D2Q5 TRT

主程序：`2DRBOpenaccLBMCDE_D2Q9TRT_D2Q5LuoTRT.F90`

- 流场：D2Q9 TRT，偶模态保留 LBM-CDE 的 `chi_s/chi_b` 修正。
- 温度场：Luo/Wang D2Q5 矩空间 TRT，不包含相邻时刻 `T*u/T*v` 修正。
- 温度参数：`s_k=3-sqrt(3)`，`s_e/s_nu` 按 Dubois Eq. (55) 设置；该结论消除冻结纯扩散体算子的四阶误差，不等同于对流和壁面全局四阶。
- 流场魔法参数固定使用 `chi_s` 修正后的有效偶尺度：
  `H_even,eff*H_odd=3/16`；无需额外编译开关。

温度梯度仅用于 Nu、耗散和诊断输出，采用二阶有限差分，不参与 D2Q5 碰撞。
