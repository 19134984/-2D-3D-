# 高 Ra 热对流 LBM 推导：来源--任务矩阵

本文件区分三类证据：

- **文献直接给出**：可按论文公式号复核；
- **本文重新推导**：以文献模型为起点，但结论不是论文原句；
- **仍需数值验证**：理论闭合不等于当前求解器已经通过算例。

| 任务 | 直接来源 | 文献直接支持 | 本文新增推导 | 仍需验证 |
|---|---|---|---|---|
| 1. D2Q5 源项 | `LBM-CDE.pdf` Eqs. (24), (28) | 连续/离散标量源和半步因子 | 用当前 D2Q5 权重、\(N\) 矩阵投影出五个矩分量 | 完整 \(p\nabla T+T\boldsymbol F\) 分支 |
| 2. 纯扩散四阶式 | `Towards higher order lattice Boltzmann schemes.pdf` Eqs. (40)--(42), (55) | 无源 D2Q5 的 \(K_{40},K_{22}\) 和 quartic 条件 | 证明局部源碰撞等价于 \(q_{\kappa,\mathrm{eff}}\)，再独立作四阶矩阵展开 | manufactured diffusion 的实测系数 |
| 3. 恒定对流 | 同上作为矩方法参考；`Optimal Stability...pdf` 作为谱分析参考 | 对流扩散谱分析框架 | 当前 D2Q5 的 \(z_1\)--\(z_4\) 全展开 | 不同速度方向、相位误差 |
| 4. 稳定区域 | `Optimal Stability...pdf` | 非负、必要条件和全谱条件必须区分 | 当前实际 \(N,C,D(\boldsymbol k)\) 的扫描与严格四阶可行性下界 | 目标 \(Ra,Pr,Ma,N\) 的实际参数扫描 |
| 5. 流场边界 | `Multireflection boundary conditions...pdf` | 无源 link 边界和偶/奇参数控制壁面位置 | 源修正剪切矩的实际放大率为 \(s_{\nu,\mathrm{eff}}\)，故标准 HBB 使用有效参数 | Poiseuille 壁面滑移；新 link correction |
| 6. D2Q9/D3Q19 投影 | `LBM-CDE.pdf` Eq. (26)；仓库两个核心 OpenMP 文件 | 张量源形式、Guo 力半步离散 | 逐行乘当前实际 \(M\) 得到全部矩分量 | D3Q19 编译/运行与矩级单元测试 |
| 7. 应变率 | `LBM-CDE.pdf` Eqs. (31)--(33) | 非平衡二阶矩的局部公式 | 映射到当前黏度、力和维数记号 | 自然对流中的机器精度/误差 |
| 8. 温度梯度 | `LBM-CDE.pdf` Eq. (35) | 含 \(TF_\alpha/\rho_0+u_\alpha Q\) 的局部公式 | 主分支 A 的 \(q_{\kappa,\mathrm{eff}}\) 简式 | 完整热源/力耦合分支 |
| 9. 边界源一致性 | `LBM-CDE.pdf` Eqs. (36)--(39) | HBB、ABB、绝热 BB 的基础 link 形式 | 线性导热下 ABB 半格点与 BB 零通量证明 | 曲率温度场的隐藏误差 |
| 10. D2Q9 判定 | `Towards higher order...pdf` 的 D2Q9 章节；二维热对流基线论文 | D2Q9 有对角速度和更多矩自由度 | 给出基于 \(h_{\min}^{(4)}\)、交叉扩散和全谱空集的触发条件 | D2Q5 与 D2Q9 成本/精度比较 |

## 已核对的文献身份

1. `LBM-CDE.pdf` 是输运系数源项解耦的核心论文。其原始实现以 BGK、D2Q9/D3Q27 为主，不能直接当作当前 MRT D2Q9/D3Q19 的逐矩公式。
2. `[Luo2014]_JCP_LB simulations of the thermally driven 2D square cavity at high Rayleigh numbers.pdf` 实际对应 Contrino 等人的高 \(Ra\) 二维热方腔论文，DOI `10.1016/j.jcp.2014.06.047`。
3. `Lattice Boltzmann simulations of thermal convective flows in two dimensions.pdf`，DOI `10.1016/j.camwa.2012.07.001`，用于二维 D2Q9/D2Q5 热模型和基准。
4. `Lattice Boltzmann simulations of three-dimensional thermal convective flows at high Rayleigh number.pdf`，DOI `10.1016/j.ijheatmasstransfer.2019.06.002`，用于 D3Q19/D3Q7 高 \(Ra\) 基线与验证。
5. `Multireflection boundary conditions for lattice Boltzmann models.pdf`，DOI `10.1103/PhysRevE.68.066614`，用于 link 边界理论；它没有直接证明加入当前应力源后仍使用原 \(3/16\)。
6. `Towards higher order lattice Boltzmann schemes.pdf`，DOI `10.1088/1742-5468/2009/06/P06006`，直接给出无源 D2Q5 四阶公式；含源结论由本文重新展开，而不是把论文公式未经证明地改名。
7. `Optimal Stability of Advection-Diffusion Lattice Boltzmann Models with Two Relaxation Times for Positive_Negative Equilibrium.pdf`，DOI `10.1007/s10955-010-9969-9`，提供稳定性分析框架；具体“最优参数”不能脱离其平衡权重和 TRT 假设直接移植。
8. `Accelerated lattice Boltzmann simulation using GPU and OpenACC with data management.pdf`，DOI `10.1016/j.ijheatmasstransfer.2017.02.032`，主要支持 GPU/OpenACC 数据组织，不是新稳定化公式的理论来源。

## 不得越过的结论边界

- \(3/16\) 的含源版本是本文在“标准 HBB 反射完整 post-collision 分布”假设下的推导，仍需 Poiseuille 壁面滑移验证。
- \(1/6\) 的含源版本之所以成立，是因为主分支 A 在常系数线性体节点中与有效通量碰撞代数等价，并已用独立四阶符号展开验证；该结论不能外推到 \(Q\ne0\)、变速度或完整压力/力耦合源。
- D3Q19 投影来源是实际代码矩阵乘法，不是原论文直接列出的 MRT 行向量；必须做矩级数值单元测试。
- 参数谱扫描只是代表性切片，不是目标高 \(Ra\) 算例的运行证明。
