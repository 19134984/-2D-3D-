# A comparative study of the simplified thermal lattice Boltzmann method with different lattice velocity models

## 元数据

- 题名: A comparative study of the simplified thermal lattice Boltzmann method with different lattice velocity models
- 作者: Zheng-Wei He, Chang Shu, Zhen Chen
- 期刊: Physics of Fluids
- 年份: 2025
- 卷期页码: 37, 113601
- DOI: 10.1063/5.0300977
- 投稿/接收/发表: 2025-09-05 / 2025-10-16 / 2025-11-03
- 建议标签: LBM, STLBM, thermal convection, stability analysis, lattice velocity models
- 期刊分区/IF: 未在本文 PDF 中给出，未设置

## 主要工作

这篇文章比较的是简化热格子 Boltzmann 方法（STLBM）中温度场采用不同热格子速度模型时的表现，核心对象是二维热场的 D2Q4、D2Q5 和 D2Q9。作者关心的不是传统 LBM 中流场模型怎么选，而是在 STLBM 这个“直接演化宏观变量”的框架下，温度分布函数的速度集合减少以后，会不会影响收敛、精度、稳定性和计算效率。

文章的主结论很实用：只要数值稳定能够保证，D2Q4、D2Q5、D2Q9 对 STLBM 的精度影响很小；真正明显的差别在稳定性和总耗时上。D2Q4 每步最省，且在一个高 Ra 稳态自然对流测试中总时间最低，但它在高松弛参数的细网格下最先发散。D2Q9 稳定性最好，但每步成本最高。D2Q5 是折中方案，在低 Re/Ra 或可能非定常的热流问题中更稳健。

## 这篇文章真正解决的问题

STLBM 已被用于高 Reynolds/Rayleigh 数热流、纳米流体、相变、浸没边界等问题，但既有 STLBM 工作大多沿用 D2Q9 作为二维温度场格子模型。这个选择有历史惯性：传统 LBM 里 D2Q4/D2Q5 常被认为稳定性不如 D2Q9，而已有 STLBM 稳定性理论也主要基于 D2Q9。本文补上的空白是：在 STLBM 的预测-校正宏观变量演化方式下，热格子速度模型是否真的需要 D2Q9，还是可以用更少方向换取效率。

这对热对流代码很直接。你的当前代码里温度场的离散速度数越多，`g` 或温度相关分布函数数组越大，每个格点的流式/重构/累加循环越重。如果模型精度差不多，D2Q5 或 D2Q4 的吸引力就在于减少温度场方向数。但如果在高松弛参数或细网格下稳定性变差，省下的每步成本可能被发散、重启、调参或更慢收敛抵消。

## 从引言值得学的背景

作者把 STLBM 放在传统 LBM 的三个痛点之后：一是 BGK LBM 在高 Re/Ra 下稳定性不足，二是分布函数数量大导致内存开销高，三是物理边界条件的实现有时不直接。STLBM 的思路是保留 LBM 的显式和格子结构，但把演化对象转为密度、速度、温度等宏观变量，避免在所有时间层都保存完整分布函数。

这条发展线很重要：STLBM 并不是简单把传统 DDF 热 LBM 的 `f_i, g_i` 改名，而是在每一步用平衡分布函数从邻居位置重构宏观变量，再用校正项补上粘性和热扩散效应。因此，传统 LBM 中“速度集合阶数越高越稳/越准”的经验不能直接套用。本文正是用 von Neumann 稳定性分析和两个 benchmark 去判断这个经验在 STLBM 中是否仍成立。

引言还强调，STLBM 已经用于 `Ra = 10^6` 的相变问题，以及 `Ra = 10^6` 到 `10^9` 的高强度对流问题。这说明作者不是在讨论一个只适合教材算例的小模型，而是在为高 Rayleigh 数热流应用选择更经济的热格子模型。

## 核心模型

本文采用双分布函数背景下的 STLBM。流场仍用 D2Q9 以保证质量和动量守恒，温度场比较 D2Q4、D2Q5 和 D2Q9。STLBM 每步由 predictor 和 corrector 两部分组成。

预测步公式为：

```text
rho* = sum_alpha f_alpha^eq(r - e_alpha dt, t - dt)                         (1)
(rho u)* = sum_alpha e_alpha f_alpha^eq(r - e_alpha dt, t - dt)              (2)
T* = sum_alpha g_alpha^eq(r - e_alpha dt, t - dt)                            (3)
```

含义是：不用保存并推进完整的非平衡分布函数，而是从上一步邻居位置的平衡分布函数重构当前的密度、动量和温度预测值。这一点和常规 pull-streaming 的数据访问形态相似，但物理变量是宏观量。

校正步公式为：

```text
rho = rho*                                                                 (4)
rho u = (rho u)* + (tau_v - 1) [ sum_alpha e_alpha f_alpha^eq(r + e_alpha dt, t)
        - rho(r, t - dt) u(r, t - dt) ] + F_E dt                            (5)
T = T* + (tau_c - 1) [ sum_alpha e_alpha g_alpha^eq(r + e_alpha dt, t)
        - T(r, t - dt) ]                                                     (6)
```

这三式说明密度直接采用预测值；动量和温度则通过松弛参数相关的校正项恢复宏观方程中的粘性和扩散行为。`tau_v` 控制动量扩散，`tau_c` 控制热扩散。

浮力项采用 Boussinesq 近似：

```text
F_E = [0, -rho g beta (T - T0)]^T                                            (7)
```

这里 `g` 是重力加速度，`beta` 是热膨胀系数，`T0` 是平均温度。本文默认低 Mach 数不可压缩热流，并采用被动标量传热模型，所以忽略压缩功和粘性热。

松弛参数和物性关系为：

```text
nu  = c_s^2 (tau_v - 0.5) dx^2 / dt                                          (8)
chi = c_s^2 (tau_c - 0.5) dx^2 / dt                                          (9)
```

`nu` 是运动粘度，`chi` 是热扩散率。这个关系对你的代码很关键：当物性固定、网格加密或 `dx/dt` 变化时，`tau_v` 和 `tau_c` 会移动；本文后面发现，高松弛参数场景下不同热格子模型稳定性差异会变得明显。

流场 D2Q9 平衡分布函数为：

```text
f_alpha^eq(r,t) = w_alpha rho [1 + (e_alpha . u)/c_s^2
    + (e_alpha . u)^2/(2 c_s^4) - u^2/(2 c_s^2)]                             (10)
```

D2Q9 速度集合为：

```text
e_alpha = (0,0), alpha = 0
e_alpha = (+/-1,0), (0,+/-1), alpha = 1..4
e_alpha = (+/-1,+/-1), alpha = 5..8                                          (11)
```

D2Q9 权重和声速为：

```text
w_0 = 1/9, w_1..w_4 = 4/9, w_5..w_8 = 1/36                                  (12a)
c_s = 1/sqrt(3)                                                              (12b)
```

注意：这里的权重按 PDF 页面显示记录。若用于实际复现，应回到原文公式或作者代码再次核对，因为常见 D2Q9 权重通常为 `4/9, 1/9, 1/36`，本文 PDF 页面中 `w_0` 和轴向权重的排版可能需要谨慎确认。

温度场比较的速度集合为：

```text
D2Q4: e_alpha = (+/-1,0), (0,+/-1), alpha = 1..4                             (13a)
D2Q5: e_alpha = (0,0), alpha = 0;
      e_alpha = (+/-1,0), (0,+/-1), alpha = 1..4                             (13b)
```

热平衡分布函数可以写成二阶形式：

```text
g_alpha^eq(r,t) = w_alpha T [1 + (e_alpha . u)/c_s^2
    + (e_alpha . u)^2/(2 c_s^4) - u^2/(2 c_s^2)]                             (14a)
```

也可以写成线性形式：

```text
g_alpha^eq(r,t) = w_alpha T [1 + (e_alpha . u)/c_s^2]                         (14b)
```

本文为了简单和效率采用线性形式 (14b)。这意味着温度场平衡函数只保留速度的一阶耦合项，适合本文低 Mach 数、被动标量热传输的假设。

D2Q4/D2Q5 的热权重和声速为：

```text
D2Q4: w_alpha = 1/4, alpha = 1..4                                             (15a)
D2Q5: w_0 = 1/3, w_1..w_4 = 1/6                                              (15b)
c_s = 1/sqrt(3)                                                              (15c)
```

## 稳定性分析主线

作者先在等温 SLBM 中复习 von Neumann 分析。把变量定义为：

```text
y1 = rho, y2 = rho u, y3 = rho v                                             (16a)
y1* = rho*, y2* = (rho u)*, y3* = (rho v)*                                   (16b)
Y = (y1, y2, y3)^T                                                           (17a)
Y* = (y1*, y2*, y3*)^T                                                       (17b)
dY^{n+1} = G dY^n                                                           (18)
lambda_max = f(rho, u, v, theta_x, theta_y, tau_v)                           (19)
```

`G` 是特征矩阵，`|lambda_max| <= 1` 是线性稳定判据。物理上，`|lambda_max|` 表示一个数值扰动经过一个时间步后的最坏放大因子。小于 1 说明扰动会衰减；等于 1 是临界；大于 1 会增长。

已有分析常把平均速度方向和波数方向都简化为水平方向。本文指出，这个一维波数假设信息不足。因为二维离散速度集合可能有方向各向异性，扰动也可能沿斜向、边缘或角落频率传播。作者因此放开波数向量，让 `k_x` 和 `k_y` 都在 `[-pi, pi]` 中变化，并用二维灰度图观察不同频率方向上的最大特征值。

热流条件下变量扩展为：

```text
y1 = rho, y2 = rho u, y3 = rho v, y4 = T                                     (20a)
y1* = rho*, y2* = (rho u)*, y3* = (rho v)*, y4* = T*                         (20b)
Y = (y1, y2, y3, y4)^T                                                       (21a)
Y* = (y1*, y2*, y3*, y4*)^T                                                  (21b)
```

然后把从邻居位置来的平衡分布函数记为：

```text
f_alpha^eq = f_alpha^eq(r - e_alpha dt, t - dt)
g_beta^eq  = g_beta^eq(r - e_beta dt, t - dt)
y_{j,alpha/beta} = y_j(r - e_{alpha/beta} dt, t - dt)                       (22)
```

将流场和热场平衡函数代入后，得到用 `y1..y4` 表示的 `f_alpha^eq` 和 `g_beta^eq`：

```text
f_alpha^eq = w_alpha^v y1
  + w_alpha^v/c_s^2 (e_x alpha y2 + e_y alpha y3)
  + w_alpha^v/(2 c_s^4) (e_x alpha y2 + e_y alpha y3)^2 / y1
  - w_alpha^v/(2 c_s^2) (y2^2 + y3^2) / y1                                  (23a)

g_beta^eq = w_beta^c y4
  + w_beta^c/c_s^2 y4 (e_x beta y2 + e_y beta y3) / y1                       (23b)
```

这两式是稳定性推导的核心，因为不同 D2Qq 热模型改变的是 `g_beta^eq` 中的 `e_beta` 和 `w_beta^c`，也就是温度扰动如何沿不同方向被传播和放大。

预测步可写成非线性映射：

```text
Y* = F(y1^n, y2^n, y3^n, y4^n) = F(Y^n)                                     (24)
```

线性化后得到：

```text
dY* = E dY^n                                                                 (25a)
dy_j* = E_jl dy_l^n, j,l = 1..4                                              (25b)
E_jl = partial y_j* / partial y_l^n                                          (26)
```

校正步同理得到：

```text
dY^{n+1} = G dY^n                                                            (27a)
dy_j^{n+1} = G_jl dy_l^n, j,l = 1..4                                         (27b)
```

热流条件下最大特征值写为：

```text
lambda_max = f(rho, u, v, T, theta_x, theta_y, tau_v, tau_c)                 (28)
```

`theta_x = k_x dt`，`theta_y = k_y dt`。这说明稳定性不仅由 `tau_v/tau_c` 控制，也和扰动方向有关。本文最关键的理论发现正来自这一点：如果只看单方向波数，D2Q4、D2Q5、D2Q9 看起来差不多；一旦看二维波数平面，高频角落区域的差异就暴露出来。

## 复杂点展开：为什么 D2Q4 会更容易在高松弛参数下失稳

为什么重要：D2Q4 方向最少，单步成本最低；如果它总是稳定，那它会是最有吸引力的热模型。但本文数值实验显示，在 `Ra = 10^4, Pr = 0.71` 的自然对流细网格中，D2Q4 在 `301 x 301` 就发散，而 D2Q5 和 D2Q9 仍收敛。

机制解释：二维 von Neumann 图中，角落区域对应 `|k_x|, |k_y| -> pi` 的高频扰动。D2Q4 只有四个轴向速度，没有静止粒子和对角方向。它对某些高频组合的抑制能力弱，所以最大特征值区域会向 1 甚至超过 1 逼近。D2Q5 加入静止方向，D2Q9 进一步加入对角方向，能更好地压制这些高频扰动。

论文证据：Sec. III 的 Fig. 5 显示，在 `tau_v = tau_c = 0.5` 的极低松弛参数分析中，D2Q5 和 D2Q9 对高频扰动呈现类似等温 SLBM 的抑制特征，而 D2Q4 的高频角落区域接近临界。Sec. IV B 2 的 Table VIII 进一步给出实际发散顺序：`301 x 301` 时 D2Q4 发散，D2Q5/D2Q9 收敛；`361 x 361` 时 D2Q4 和 D2Q5 发散，D2Q9 收敛；`401 x 401` 时三者都发散。Fig. 16 解释了这些发散：D2Q4 在 `301 x 301` 的最大特征值可达 1.456，D2Q5/D2Q9 在更细网格时才形成不稳定区。

对代码的含义：如果你要在当前热对流 LBM 代码中尝试减少温度场方向数，不应只看每步方向数。需要同时检查 `tau_c` 随网格、`Ra/Pr`、`dt/dx` 的变化。如果某组参数把 `tau_c` 推到较大区间，D2Q4 的风险会比 D2Q5/D2Q9 高。

## 数值实验和可复用 benchmark

### 收敛准则

所有模拟使用相对变化作为收敛判据：

```text
max(V_err, T_err) < 1e-8                                                     (29)
```

其中 `V_err` 是速度场相邻时间步相对 L2 变化，`T_err` 是温度场相邻时间步相对 L2 变化。这个准则适合稳态问题，与你当前代码中检查 `errorU/errorT` 的思路相近。

### Porous plate problem

这个算例用于验证精度，因为有解析解。物理设置是完全发展通道流：上冷板以水平速度 `u0` 移动，下热板以竖直速度 `v0` 注入，上板以相同速度抽出；左右边界周期；几何为 `L x H`，且 `L = 2H`。

解析解：

```text
u = u0 [exp(Re y/H) - 1] / [exp(Re) - 1]                                    (30a)
T = Th - DeltaT [exp(Pr Re y/H) - 1] / [exp(Pr Re) - 1]                     (30b)
```

其中：

```text
Re = v0 H / nu
Pr = nu / chi
DeltaT = Th - Tc
Ra = g beta DeltaT H^3 / (nu chi)
```

本文参数：

```text
Th = 1.0, Tc = 0.0, nu = 0.1, u0 = 0.1, L = 2.0, H = 1.0, Ra = 100
```

网格：

```text
41 x 21, 81 x 41, 121 x 61, 161 x 81, 201 x 101
```

误差定义为相对 L2 误差：

```text
Err_V = sqrt(sum_i [u_num(x_i) - u_ana(x_i)]^2)
        / sqrt(sum_i [u_ana(x_i)]^2)                                        (31a)

Err_T = sqrt(sum_i [T_num(x_i) - T_ana(x_i)]^2)
        / sqrt(sum_i [T_ana(x_i)]^2)                                        (31b)
```

关键结果：Fig. 7 显示 D2Q4、D2Q5、D2Q9 的速度误差收敛阶均约 1.93，温度误差收敛阶约 1.99、2.00、2.00。也就是说三种热速度模型在该解析算例中都接近二阶精度。Figs. 8-11 进一步用 `Pr = 0.71, Re = 5/10/20` 和 `Re = 10, Pr = 0.2/0.8/1.5` 的中心线速度/温度剖面验证，数值结果与解析解吻合。

对你代码的可复用价值：这个算例适合用来检查温度场离散速度模型切换是否破坏二阶收敛，尤其适合在实现 D2Q5 或 D2Q4 热场后做单元级验证。

### 2D square cavity natural convection

这个算例更接近热对流主问题。设置为二维方腔 `L x L`，左壁热 `Th`，右壁冷 `Tc`，上下壁绝热。Prandtl 数固定为 `Pr = 0.71`，特征速度 `Vc = sqrt(g beta DeltaT L)` 取相关参数设置，考虑多个 Rayleigh 数。

控制参数定义：

```text
Ra = g beta DeltaT L^3 / (nu chi) = Vc^2 L^2 / (nu chi)                     (32)
Pr = nu / chi                                                               (33)
```

Nusselt 数定义：

```text
Nu_volume = L/(chi DeltaT L^2) integral [u T - chi partial_x T] dX          (34)
Nu0_hotwall = 1/DeltaT integral (partial_x T)|_{x=x0} dy                    (35)
Nu_local = L/DeltaT partial_x T                                             (36)
```

无量纲速度定义：

```text
u' = u L / chi                                                              (37)
```

精度测试使用 `Ra = 10^4, 10^5, 10^6`，并比较最大水平速度 `u_max`、最大竖直速度 `v_max`、极值位置，以及体平均/壁面平均/局部最大最小 Nusselt 数。参考解包括 DQ、TLBM、FDM、CSLBM。

关键结果：

- `Ra = 10^4`：在 `301 x 301` 网格下 D2Q4 发散，D2Q5 和 D2Q9 收敛；D2Q5/D2Q9 的 `u_max` 约 16.186/16.187，`v_max` 约 19.639/19.640，与 DQ 参考值 16.190/19.638 很接近。
- `Ra = 10^5`：三种模型在 `301 x 301` 均收敛，`u_max` 约 34.754/34.746/34.757，`v_max` 约 68.614/68.611/68.620，和 DQ 参考值 34.736/68.640 接近。
- `Ra = 10^6`：三种模型在 `301 x 301` 均收敛，`u_max` 约 64.618/64.584/64.704，`v_max` 约 219.83/219.77/219.85，和 DQ 参考值 64.775/220.64 接近。
- 热量传输指标方面，Tables IV-VI 显示三种模型的 Nusselt 数基本处在 DQ、TLBM、FDM、CSLBM 参考范围内。`Ra = 10^6, 301 x 301` 时，三种模型的体平均 `Nu` 约 8.807/8.801/8.801，壁面平均 `Nu0` 约 8.802/8.802/8.799。
- Figs. 13-14 显示流线和温度等值线大部分重合。差异主要出现在 `Ra = 10^6` 顶部附近，D2Q4/D2Q5 相对 D2Q9 有轻微偏离。

对你代码的可复用价值：这组方腔自然对流可以作为热 LBM 版本切换模型后的主 benchmark。应记录 `Ra, Pr, grid, errorU/errorT, Nu, Nu0, Nu_local, u_max/v_max`，同时观察是否在细网格高 `tau_c` 情况下出现模型相关发散。

## 稳定性测试结果

低松弛参数测试用非常粗的 `11 x 11` 网格和高 Rayleigh 数 `Ra = 10^7, 10^8, 10^9`，使 `tau_v/tau_c` 接近 0.5：

```text
Ra = 1e7: tau_v = 0.500799, tau_c = 0.501126
Ra = 1e8: tau_v = 0.500253, tau_c = 0.500356
Ra = 1e9: tau_v = 0.500080, tau_c = 0.500113
```

Table VII 显示 D2Q4、D2Q5、D2Q9 全部收敛。这说明 STLBM 在接近传统 LBM 容易失稳的低松弛参数区间仍然很稳。

高松弛参数测试通过细网格获得：

```text
301 x 301: tau_v = 1.258353, tau_c = 1.568103
361 x 361: tau_v = 1.410024, tau_c = 1.781724
401 x 401: tau_v = 1.511138, tau_c = 1.924138
```

Table VIII 的收敛/发散顺序为：

```text
301 x 301: D2Q4 diverge, D2Q5 converge, D2Q9 converge
361 x 361: D2Q4 diverge, D2Q5 diverge, D2Q9 converge
401 x 401: D2Q4 diverge, D2Q5 diverge, D2Q9 diverge
```

这个结果是全文最值得抓住的稳定性证据：热速度模型方向数越多，高松弛参数下的稳定边界越宽。D2Q9 最稳，D2Q5 次之，D2Q4 最脆弱。

## 计算效率结果

效率测试选择 `Ra = 10^6, 301 x 301` 的二维自然对流。以 D2Q9 总耗时为 1.00，Table IX 给出：

```text
D2Q4: total time = 0.65, iterative steps = 285000, time per 1000 dt = 0.84
D2Q5: total time = 0.85, iterative steps = 366000, time per 1000 dt = 0.85
D2Q9: total time = 1.00, iterative steps = 366000, time per 1000 dt = 1.00
```

两个层次要分开看：

- 每 1000 步耗时：D2Q4/D2Q5 比 D2Q9 少约 16%/15%，这是热格子方向减少带来的直接收益。
- 总耗时：D2Q4 比 D2Q9 少约 35%，不仅因为每步便宜，还因为该算例中 D2Q4 需要的迭代步数更少。

作者也提醒，效率和程序实现策略有关，尤其并行代码中还会受内存访问、数组布局、通信和 GPU 并行粒度影响。本文所有效率测试都在串行代码中完成。因此不能直接把 35% 当成 OpenACC/MPI 代码中的收益上限，但它能说明模型方向数确实会影响温度场循环成本。

## 作者给出的模型选择指南

作者的选择建议是：

- 高 Reynolds/Rayleigh 数且期望稳态的热流问题：D2Q4 是首选，因为在本文典型算例中收敛效率最好。
- 低 Reynolds/Rayleigh 数或预期非定常的热流问题：D2Q5 更合适，因为它每步成本接近 D2Q4，但比 D2Q4 有更宽的稳定参数范围。
- D2Q9 更稳，但成本最高；当稳定性优先、参数极端、或者需要保守实现时可以作为稳健基准。

我对这条建议的理解是：D2Q4 适合“已知会稳定收敛的稳态问题”的性能优化；D2Q5 适合“还不确定是否会遇到稳定性问题”的默认轻量模型；D2Q9 适合做基准、排错和高风险参数区间。

## 对当前热对流代码的直接启发

如果你想把当前 2D/3D 热对流代码中的温度场速度集合从 D2Q9/D3Q? 简化，应该先把“精度”和“稳定性”分开验证。本文显示精度差异很小，但稳定性差异很大，尤其在高松弛参数细网格下。

代码层面最值得检查的是：

- `tau_c` 或热扩散率相关参数随 `Ra, Pr, dx, dt` 的变化范围；
- 温度场分布函数/平衡函数循环是否可以独立切换速度集合；
- Nusselt 数计算是否只依赖宏观温度梯度，避免被速度集合切换影响后处理；
- 方腔自然对流 benchmark 中是否同时比较 `Nu`、`u_max/v_max`、极值位置和流线/等温线；
- OpenACC/MPI 版本里，减少热方向数是否真的减少 device memory traffic，还是被 halo exchange、host-device update 或其他数组访问掩盖。

对你的 `3DRBOpenaccMpi.F90` 更具体地说，本文不是 MPI/OpenACC 性能论文，不能直接给出 GPU 并行收益。但它提供了一个模型选择逻辑：如果温度场方向数减少，GPU 上每格点读写和循环次数会下降；不过在 MPI+OpenACC 中还要重新评估 halo 打包、边界层访问和内存布局。D2Q4/D2Q5 在串行代码中省时间，不保证在多 GPU 上按同样比例省时间。

## 局限性

作者承认本文的稳定性理论是线性 von Neumann 分析，默认均匀和周期性背景。它不能完整覆盖以下情况：

- buoyancy-coupled operator 的非正规效应导致的瞬态增长；
- 具体壁面边界条件与 D2Q4 等简化热模型之间的相互作用；
- 强各向异性或强剪切热流中的完全非线性失稳发展；
- 复杂几何中的高保真热流模拟。

因此，本文的建议适合作为模型选择的初始准则，但不能替代你在目标边界条件、目标 Ra/Pr、目标网格和目标并行实现下的实际验证。

## 原始证据地图

- 题名、作者、期刊、DOI：PDF 首页，Physics of Fluids 37, 113601 (2025), DOI 10.1063/5.0300977。
- STLBM predictor-corrector：Sec. II, Eqs. (1)-(6)。
- Boussinesq 浮力和松弛参数-物性关系：Sec. II, Eqs. (7)-(9)。
- D2Q9 流场平衡函数、速度集合、权重：Sec. II, Eqs. (10)-(12)。
- D2Q4/D2Q5 热速度集合和热平衡函数：Sec. II, Eqs. (13)-(15)。
- von Neumann 稳定性变量和特征矩阵：Sec. III, Eqs. (16)-(28), Figs. 2-5。
- porous plate 精度测试：Sec. IV A, Eqs. (29)-(31), Figs. 6-11。
- 方腔自然对流定义和 Nu 指标：Sec. IV B, Eqs. (32)-(37), Tables I-VI, Figs. 12-14。
- 低/高松弛参数稳定性：Sec. IV B 2, Tables VII-VIII, Figs. 15-16。
- 计算效率：Sec. IV B 3, Table IX。
- 模型选择建议和局限：Sec. V, Fig. 17。
