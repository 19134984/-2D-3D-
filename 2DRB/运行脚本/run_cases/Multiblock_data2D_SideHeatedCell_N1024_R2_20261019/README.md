# 多块网格二维侧壁稳态：node07 串行作业

PBS 作业号：`7265.master`。只提交一个作业，顺序为 `Ra1e6 → Ra1e7 → Ra1e8`，每个求解器以前台方式结束后才启动下一个，没有后台并行或额外 qsub。独立算例若失败，记录退出码后继续下一项；任一失败使整批最终退出码为非零。整批 walltime 为 96 小时。

远端根目录：`/data2/XLLi/Multiblock/data2D/SideHeatedCell/`。每个 Ra 下的 `source/` 保存源码快照与来源记录，`results/` 保存数据、编译及运行状态；`serial_node07/` 保存唯一 PBS 脚本、总日志与逐算例执行状态。`RB/` 仅在上级 data2D 下预留，没有提交非稳态。

三个算例均为 nx=ny=1024、refineRatio=2、四周细层参数128/129/128/129、Pr=0.71、Mach=0.1、epsU=epsT=1e-7，从零开始。最大2000万细步，稳态场检查间隔2000细步；Nu/Re和快照每10 t_ff，检查点与Tecplot每100 t_ff。保留 ReVolRMS 和 ReVolMean，最终文献对比使用实际收敛时刻的 ReVolMean。

提交时的本地母版校验值见 campaign.json；参数化差异见各 source/master_to_case.diff；case source SHA-256 在 source/provenance.json；PBS实际编译前在 results/runtime_source.sha256 记录并核对。未改变本地主程序的默认算例设置。此目录是已提交快照，不随母版后续修改更新。
