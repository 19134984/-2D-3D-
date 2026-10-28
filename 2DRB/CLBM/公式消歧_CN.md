# Sharma 2018 CTLBM 的公式核对

## 1. 本地原文与补充来源

主要依据为用户提供的 2018 年论文 PDF，尤其式(5)–(7)、(16)–(19)、(22)、(27)–(30)、(47)–(48)。这些公式已检查高分辨率页面，未仅依赖文本提取。PDF 路径为原工作目录 `pdf/算法/Natural convection heat transfer modeling by the cascaded thermal lattice Boltzmann method.pdf`。

另检索了原文所引的两篇文献：

1. Sharma 等，2017，*New Cascaded Thermal Lattice Boltzmann Method for simulations of advection-diffusion and convective heat transfer*，[出版社记录](https://doi.org/10.1016/j.ijthermalsci.2017.04.020)。本次可核对出版社摘要、作者和模型结构，正文访问返回 403；没有把摘要当作已读全文。作者的[公开博士论文记录](https://pantheon.ufrj.br/handle/11422/12788)也已找到，但其 PDF 文件链接返回 404，未用未读正文支持公式结论。
2. Cui 等，2016，*Discrete effect on the halfway bounce-back boundary condition of multiple-relaxation-time lattice Boltzmann model for convection-diffusion equations*，[APS 原始记录](https://doi.org/10.1103/PhysRevE.93.043311)。该文摘要确认可以通过二阶松弛参数消除离散壁面滑移；本次正文接口返回 401，不能声称逐式读到了该文。

对于具体松弛关系，核对了另一份可访问的原始研究全文，即 Sharma 2018 的参考文献 [46]：Fei 与 Luo，2018，*Cascaded lattice Boltzmann method for incompressible thermal flows with heat sources and general thermal boundary conditions*，[作者所在大学公开全文](https://discovery.ucl.ac.uk/id/eprint/10043940/1/1-s2.0-S0045793018300288-main.pdf)，DOI 10.1016/j.compfluid.2018.01.020。其第 4 页明确以松弛率写出关系，并给出两组数值测试。本实现的通式还经过下面的独立推导与实际离散程序核验，并非依靠无法读取的引用补齐。

## 2. 参数 a 与温度格子声速

Sharma 的式(22)在零速度时给出静止权重 `1-a`，每个轴向权重 `a/4`，因此二阶格子矩为

`chi = cT^2 = sum(w_i*e_ix^2) = a/2`。

式(29)也给出 `alpha = (a/2)*(tauFlux-1/2)`，两处相互一致。第 4 页文字把 `a` 设为 `sqrt(2/5)`；已核对 PDF 矢量根号横线覆盖整个 `2/5`。程序因此保留：

`thermalA = sqrt(2/5)`，`thermalCs2 = thermalA/2 ≈ 0.3162277660`。

Fei 与 Luo 的另一篇文章使用 `cT=sqrt(2/5)`，即 `cT^2=0.4`，但它与这里的 `a` 不是同一个参数定义。不能为了套用一条松弛关系，把 Sharma 的参数悄悄改为 `a=0.8`。

## 3. 式(30)为何不能直接当松弛时间

Sharma 式(26)–(29)把 `tau` 定义为松弛时间，碰撞乘子为 `1/tau`。然而式(30)的发表形式为

`tauEven = 12*(1-2*tauFlux)/(1-12*tauFlux)`。

例如 `tauFlux=0.6` 时右端约为 `0.3871`；若把它作为时间再取倒数，二阶松弛率就超过 2。

Fei 与 Luo 的公开全文用清楚的松弛率记号给出 `sEven=12*(sFlux-2)/(sFlux-12)`。代入 `sFlux=1/tauFlux` 后，右端恰好成为 Sharma 式(30)的右端，说明其左端至少混用了时间与率。Fei 的两组检查为 `(sFlux,sEven)=(4/3,3/4)` 和 `(3/4,4/3)`，可逐项代回，而不是凭数值直觉取倒数。

进一步换为时间可得该专用关系的乘积 `(tauFlux-1/2)*(tauEven-1/2)=5/24`。它适用于那篇文章的 `chi=0.4`，不能在改变权重后继续称为同一个零滑移条件。

## 4. 按实际权重独立推导零滑移关系

考虑速度为零的一维稳态扩散，方向法向为 x，均匀热源为 Q，半程 ABB 壁面在 x=0，首个流体节点位于 x=1/2。记 `theta=T-Q/2` 为存储分布的零阶矩，`j=g_+-g_-`、`epsilon=g_++g_-`，`chi=a/2`。热源使用二阶时间积分，其二阶矩为 `chi*Q`。

由离散稳态碰撞与迁移，若 `T=A*x^2+B*x+C`，可以得到：

- `Q=-2*alpha*A`，`alpha=chi*(tauFlux-1/2)`；
- `j=-chi*tauFlux*T'(x)`；
- `epsilon=chi*T+[chi*(tauEven-1/2)-tauEven]*Q`。

半程 ABB 为 `g_+(1/2)=-g_-^*(1/2)+chi*Tw`。代入上式，壁面温度滑移为

`delta = A * [1/4 - 2*(1-chi)*(tauFlux-1/2)*(tauEven-1/2)]`。

故零滑移条件为

`magic = (tauFlux-1/2)*(tauEven-1/2) = 1/[8*(1-chi)]`，

`sEven = 1/[1/2 + magic/(tauFlux-1/2)]`。

当 `chi=0.4` 时，此式还原 `magic=5/24` 和 Fei 的率公式。对于程序保留的 Sharma 字面参数，则应使用 `chi=sqrt(2/5)/2` 对应的乘积。代码采用该通式，并把原文关系与采用关系的差别明确保留在中文注释中。

本次核验不是只比较这条代数式：测试驱动使用独立推导的抛物温度和稳态分布，调用实际 `collisionT/streamingT/bouncebackT/macroT`，在测试驱动中加入均匀热源的二阶积分，检查侧热和 RB 两种方向的解析抛物解保持不变。正式求解器仍是无内部热源模型。

## 5. 温度 ABB 系数与碰撞目标

式(48)的高分辨率原图是 `g_in=-g_out+(a/2)*Tw`。文本提取容易把分式次序弄乱。它与静止轴向权重的 `2*w_axis*Tw` 完全一致，因此程序采用 `2*omegaT(alpha)*Tw`，没有使用 `2*a*Tw`。

温度碰撞直接按式(27)实施。将其转换为中心矩后，二阶目标分别为 `(chi-u^2)*T`、`(chi-v^2)*T`；这与式(22)的线性平衡分布相符。Fei 的模型使用 `chi*T` 作为二阶中心目标，属于另一种模型。本实现没有用后一目标替换前者，也不因文献的稳定性宣传就略过这个区别。

## 6. 浮力时序与去混叠

原文式(7)的第三阶源投影满足 `S21=u^2*Fy+2*u*v*Fx`、`S12=v^2*Fx+2*u*v*Fy`。按二项式换到中心矩后，它们均为零；对应第四阶中心源也为零。中心源因此只有一阶的 Fx、Fy。程序按该定义更新守恒动量，二阶及更高阶独立松弛到式(16)目标。四阶目标固定为 `rho*cs^4`，没有 FPC 的后碰撞二阶矩乘积。

宏观恢复统一为 `rho*u=sum(f*ex)+Fx/2`、`rho*v=sum(f*ey)+Fy/2`。初值的存储一阶中心矩为 `-F/2`，第三/四阶保持原文去混叠定义，避免只套用另一个模型的 Hermite 力分布。两套碰撞读取同一时刻的速度、温度；两套迁移完成后先更新 T，再更新 rho、F、u、v。

## 7. 当前结论的边界

已经完成源文件实现、主机后端数学检查、纯导热和有限步重启检查；未声称已经重现原文表格。尚未完成原文高 Ra 精度、GPU 执行、RB 长时间统计或网格收敛验证。原始引用的全文可获取性限制和上述参数歧义均保留，以便后续与作者代码或新增原文交叉核对。
