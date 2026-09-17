# Ra=10⁷ 提交记录

2026-09-16，用户要求启动下一算例。作业 **6980.master**，PBS batch，node05:ppn=1，walltime=96:00:00。

稳态侧壁差温，Ra=1e7、Pr=0.71、Ma=0.1，最细等效网格 512×512，粗细比 2，外围细环左/下接口编号 64、右/上 65。从零启动，最大 20000000 细步；每 2000 步检查 errorU/errorT，均不超过 1e-7 即停止。每 10 tff 采样，每 100 tff 保存检查点与 Tecplot。

使用通过 Ra=10⁶ 验证的修复版本：接口插值时间权重按三个标量传入 GPU，无 EnableUseG 分支，旧 D2Q5 MRT 温度算法。

提交时再次读取本地母版，其 SHA-256 与待提交副本一致；重新上传原母版快照、参数化源码、运行脚本并检查 LF、bash -n 与远端哈希。

- 本地母版：`多块网格/2DRBOpenaccMultiblock.F90`
- 母版 SHA-256：`5eff362a8eceff71d30d59457a78d22dc2d6850d97e56cf274f4d096d6187bb3`
- 参数化源码 SHA-256：`90e5702cea5255daade80cd9bed0afce1fb2f0c688e59d71a869b5a252d08c8e`
- PBS 实际编译哈希：见 `results/runtime_source.sha256`，运行脚本强制与上述参数化哈希一致。
- 远端：`/data2/XLLi/Multiblock/Ra1e7/ratio2/sideheated_N512_20260916_v3/`

提交前确认无同算例排队或运行记录，results 为空。Ra=10⁶ 作业已完成，本次没有与其并行。

node05 提交前已有其他用户的 Python 与 particle_simulation GPU 进程，约 2236 MiB 显存、99% GPU 利用率。沿用用户先前允许共享 node05 的授权；本次墙钟耗时不作为独占 GPU 性能基准。没有取消或修改其他用户作业。

提交后核查文件保存在 results；运行中不等于达到稳态，需完成后另行验收历史、最终检查点与 Nu/Re。
