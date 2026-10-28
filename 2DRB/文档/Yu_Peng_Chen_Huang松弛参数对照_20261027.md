# Yu、Peng、Chen、Huang：粗细网格松弛参数对照

核查日期：2026-10-27。以本地 `pdf/多块网格/` 下四篇完整PDF为准，核对原公式和关键页图像，并用出版社或作者机构链接确认文献身份。本次未修改求解器、运行或提交算例。

## 统一记号

s是松弛率，τ=1/s。r=Δxc/Δxf=Δtc/Δtf。本文将剪切率统一记sν，体黏性率se，四阶能量型率sε，三阶能流型率sq；不同文章的原始编号不同。

按原始声学粗细转换，对需要匹配的率：

\[
\frac1{s_i^f}-\frac12=r\left(\frac1{s_i^c}-\frac12\right).
\]

不是把s直接乘r。正向或反向计算是同一关系。

## 松弛参数策略总表

| 本地文献 | 模型 | 剪切sν | 体黏se | 四阶sε | 三阶sq | 温度 |
|---|---|---|---|---|---|---|
| Yu、Mei、Shyy 2002 | BGK | 唯一τ按r转换 | 不独立设置 | 不独立设置 | 不独立设置 | 无独立温度场 |
| Peng等2006 | D2Q9 MRT | 按r转换 | 按r转换 | 原文也要求按r转换 | 原文也要求按r转换 | 无独立温度场 |
| Chen等2016 | D2Q9 MRT | 按r转换 | 按r转换 | 作者选两级相同 | 作者选两级相同 | 无独立温度场 |
| Huang与Wu2014 | 流场与温度均D2Q9 MRT | 按r转换，原文sp^u | 按r转换 | 按r转换 | 按r转换 | 全部非守恒温度率按r转换 |

“按r转换”描述的是相应论文的方案；不表示所有MRT方法必须采用同一选法。守恒矩的松弛率不承担物理输运匹配，不能将非守恒缩放公式机械套到s=0的守恒项。

## 原文定位与身份

1. **Dazhi Yu, Renwei Mei, Wei Shyy (2002)**. *A multi-block lattice Boltzmann method for viscous fluid flows*. IJNMF 39(2),99–120。[Wiley](https://onlinelibrary.wiley.com/doi/10.1002/fld.280)。p101式(8)/(10)为BGK；p103式(13)给唯一τ的转换。它不是Huang引文中的Yu博士论文，不能把那篇博士论文可能采用的MRT处理移到当前这份PDF。
2. **Y. Peng, C. Shu, Y.T. Chew, X.D. Niu, X.Y. Lu (2006)**. *Application of multi-block approach in the immersed boundary–lattice Boltzmann method for viscous fluid flows*. JCP218,460–478。[出版社](https://www.sciencedirect.com/science/article/abs/pii/S0021999106000933)，[DOI](https://doi.org/10.1016/j.jcp.2006.02.017)。p468式(47)先列s2、s8、s9，紧接正文要求其他松弛率也按该式转换。s2=se、s3=sε、s5/s7=sq、s8/s9=sν。
3. **Songying Chen, Cheng Peng, Yihua Teng, Lian-Ping Wang, Kun Zhang (2016)**. *Improving lattice Boltzmann simulation of moving particles in a viscous flow using local grid refinement*. Computers & Fluids136,228–246。[作者全文](https://research.me.udel.edu/lwang/reprints/Chen_etal_CaF_2016.pdf)，[DOI](https://doi.org/10.1016/j.compfluid.2016.06.009)。p231式(16)/(17)规定剪切、体黏率的转换；p232明确sq、sε是自由率，并为方便取两级相同。
4. **Rongzong Huang, Huiying Wu (2014)**. *Multiblock approach for the passive scalar thermal lattice Boltzmann method*. PRE89,043303。[APS](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.89.043303)。p4式(20)要求全部非守恒率满足逐矩转换；式(21)排除守恒项。其热场是D2Q9，不能直接称作本项目D2Q5 Luo热场的逐项复现。

## 参数转换和非平衡矩转换必须分开

对无外力的碰前非平衡矩，Chen的逐矩关系可写为：

\[
m_i^{neq,c}=r\frac{s_i^f}{s_i^c}m_i^{neq,f}.
\]

这给出两种容易混淆的情况：

- **Chen固定sq、sε**：这些率两级相同，但其非平衡矩仍乘r，不是直接复制。
- **Peng全部非守恒率转换**：其式(43)的碰前矩传递矩阵却仅对两个偏应力项设置显式r sf/sc，其余项系数为1。其式(45)/(46)传递碰后量还含两级(I−S)因子，不能将“碰前系数1”说成“碰后也不转换”。
- **Huang**：式(19)将所有非守恒矩和半步外力一起转换；式(20)另规定松弛矩阵的转换。这两条是不同条件。
- **Yu**：式(22)是碰前分布非平衡部分的比例；式(24)/(25)是碰后比例。不能把这些比例当成τ的粗细变换。

Chen在§5.2.4、p244式(47)比较的是ghost非平衡矩映射因子r和1，而不是比较sq、sε的固定与缩放。再把e的非平衡映射设成1复现Peng处理后，界面附近速度散度振荡更大。这是其所测粒子流结果，不是高Ra热对流的普遍稳定性判据。

## Chen的具体参数示例

p240 Couette圆柱算例，r=2：

| 参数 | 粗格 | 细格 |
|---|---:|---:|
| sν | 1.2 | 6/7 |
| se | 1.5 | 1.2 |
| sq | 1.5 | 1.5 |
| sε | 1.4 | 1.4 |

该例清楚展示了“物理输运率转换、自由高阶率固定”的实际选法。不应自行替换成未经原文核对的1.64/1.54/1.9等参数。

## Huang第7页的固定参数表述

热Couette网格收敛研究在不同整体网格试验之间保持松弛率、通过降低速度保持Re。它不表示同一个多块算例的粗细两级松弛率相同，不能用它否定第4页式(20)。

## 与本项目联合反算的区别

现有联合反算：sν按r转换；se、sε仍与sν绑定而随之变化；sq由魔法参数反算。D2Q5温度Qk按r转换，Qe/Qnu由四阶误差关系反算。它既不是Chen的“独立se转换、sq/sε固定”，也不是Huang全部非守恒率按同一关系缩放的原样方案。

这些文献支持存在不同参数策略，并不证明任意一种移入当前热对流程序后就能稳定。必须同时保持碰撞、平衡态、体力和接口映射一致。

逐篇核查记录：

- [Yu](D:/桌面/代码/代码/热对流/2DRB/文档/多块网格松弛参数核查_20261027_Yu.md)
- [Peng](D:/桌面/代码/代码/热对流/2DRB/文档/多块网格松弛参数核查_20261027_Peng.md)
- [Chen](D:/桌面/代码/代码/热对流/2DRB/文档/多块网格松弛参数核查_20261027_Chen.md)
- [Huang](D:/桌面/代码/代码/热对流/2DRB/文档/多块网格松弛参数核查_20261027_Huang.md)
