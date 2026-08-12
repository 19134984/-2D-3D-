# D2Q9 TRT + D2Q9 BGK

主程序：`2DRBOpenaccLBMCDE_D2Q9TRT_D2Q9BGK.F90`

- 流场：D2Q9 TRT，在矩空间完成碰撞。输入参数为 `chi_nu` 和 `chi_b`，代码由目标运动黏度反算基础松弛时间 `tauf`；守恒矩不松弛，其余偶矩使用 `Snu`，非守恒奇矩使用 `Sq`。
- 温度场：保留原 LBM-CDE D2Q9-BGK，不替换为 Luo D2Q5。温度碰撞也写在矩空间中：零阶温度矩守恒，其余 8 个非守恒矩统一使用 `Sg=1/taug`，因此与速度空间 BGK 数学等价。
- `taug` 在 `commondata` 中直接设置，`chi_kappa` 由目标热扩散率自动换算。压力、浮力、温度梯度和 `chi_kappa*grad(T)` 耦合项均保留。
- 流场 TRT 魔法尺度由文件开头的 `FLOW_ODD_ORIGINAL_MAGIC`（默认，基于 `tauf-0.5`）和 `FLOW_ODD_EFFECTIVE_MAGIC`（基于有效黏性尺度）手动二选一；魔法参数值在 `commondata` 的 `flowMagicParameter` 中直接修改。
- 非稳态统计频率 `statisticSampleInterval` 与全场快照频率 `outputSnapshotInterval` 相互独立。Nu、Re、耗散、密度波动、速度散度和温度剖面只在 `unsteadyAverageStartTf--unsteadyAverageEndTf` 内累计，快照、Tecplot 和重启文件可从计算开始输出。
- 续算时通过 `latest.meta` 恢复推进位置，并重读统计历史，避免 Nu/Re 与耗散、可压缩性和温度剖面使用不同的累计区间。
- 后处理同时输出瞬时/运行 Re、体对流 Nu、壁面 Nu、由动能耗散和热耗散精确关系换算的两个 Nu，以及相对于体对流 Nu 的误差；还包括 `theta_rms(z)`、整体热边界层估计、`N_BL` 和网格分辨率指标。
