# 公式对应与补全选择

来源为 Alexander Nee（2025），*Evaluation of hybrid lattice Boltzmann models for laminar and turbulent natural convection*，DOI：10.1016/j.icheatmasstransfer.2025.108635。已核对本地 PDF 第3–5页公式及第2.4节算法流程。本实现只采用其中 RLB-4，不同时加入 BGK-1、TRT-2。

## 流场：式(6)–(8)、(13)–(16)

格距与流场时间步为1，`cs2=1/3`：

```text
rho = sum(f)
rho*u = sum(ex*f), rho*v = sum(ey*f)
F = (0,rho*gBeta*(T-Tref))
source_i = w_i*(e_i dot F)/cs2
f_post = feq4 + (1-omega)*fneq_regularized + source
omega = 1/(0.5+nu/cs2)
```

初始化为原始平衡态，宏观恢复、诊断和重启都使用同一 raw 动量定义。没有 FPC 的变换分布、Guo 半步力或额外力源修正；因此不能把此处强迫与原框架的半步物理速度混用，也不应先验声称整体二阶强迫精度。

式(14)以两个一维二阶 Hermite 因子相乘实现，并保留全部乘积：

```text
Hxx = ex²-cs2, Hyy = ey²-cs2
feq4 = rho*w*(1+ux*ex/cs2+ux²*Hxx/(2*cs2²))
            *(1+uy*ey/cs2+uy²*Hyy/(2*cs2²))
```

式(16)从当前分布提取 `xi_ab=sum(e_a*e_b*(f-feq4))`。式(15)首项写 `tau_ab`，其余写 `xi_ab`，但没有另一套独立 `tau_ab` 定义；本实现将其按同一非平衡二阶张量处理，不自行乘黏性或 `(1-omega/2)` 比例。

高阶非平衡量仅保留由低阶张量重构的部分：

```text
a21 = uy*xiXX + 2*ux*xiXY
a12 = ux*xiYY + 2*uy*xiXY
a22 = uy²*xiXX + 4*ux*uy*xiXY + ux²*xiYY

fneqR = w*((xiXX*Hxx+2*xiXY*ex*ey+xiYY*Hyy)/(2*cs2²)
         +(a21*Hxx*ey+a12*ex*Hyy)/(2*cs2³)
         +a22*Hxx*Hyy/(4*cs2⁴))
```

程序不提取原始分布的三阶和四阶非平衡矩作为独立演化量。该正则化会去除不能由上述低阶张量重构的高阶扰动。

## 温度：式(2)及第2.4节

用格子单位写原文的对流–扩散方程：

```text
L(T) = -u*Dx(T)-v*Dy(T)+kappa*(Dxx(T)+Dyy(T))
```

原文给了 BGK-1 的详细离散式，但没有给齐 RLB-4 的空间差分、RK阶段耦合或壁面公式。这里明确补齐如下：

1. 内区一阶导数采用中心二阶 `(T_plus-T_minus)/2`，二阶导数采用标准三点差分。
2. 先完成流场，冻结新恢复的速度；温度每个子步采用经典 `k1,k2,k3,k4`，权重 `1,2,2,1`，每个阶段重新补边界。
3. 不用保守通量形式替代原文的 `u·grad(T)`。这使离散速度不严格无散时仍可保持常温场；全离散总热量守恒则需要额外检查。
4. `thermalSubsteps` 将一个流场步等分为若干温度子步。默认1，不修改流场推进或时间尺度。

经典 RK4 对固定速度的半离散温度 ODE 是四阶时间方法。整个耦合还包含显式流场力与顺序分步，不能把“RK4”外推为整体四阶精度。生产实现的温度内区为空间二阶。

## 半格壁面与热流一致性

延续本项目格心几何：最外层流体格心距物理壁面半个格距，不把首末流体层覆盖成壁面值。恒温 ghost 采用二次延拓：

```text
Tghost = (8*Twall-6*T1+T2)/3
```

其数值壁面梯度 `(T1-Tghost)/h` 等于项目三点半格距公式 `(-8*Twall+9*T1-T2)/(3*h)`，因此温度扩散算子与壁面 Nu 使用同一壁面数值热流。绝热 ghost 取相邻流体值，法向热流为零。角点 halo 不参与五点算子中的跨对角计算。

二次恒温闭合改变了扩散矩阵边界谱，不能沿用纯内区的 `8*kappa` 上界。启动检查取保守的 `10*kappa/thermalSubsteps<=2.7`，为经典 RK4 负实轴稳定区留裕量；它不代替含对流及耦合的稳定性验证。

## 热预算与已有后处理

结束时在当前场上重算 `L(T)`，写出：

```text
dMdt_FD = sum(L(T))        # 格子单元面积为1
Qout = sum(diffusive heat directed out of all walls)
Rheat = dMdt_FD+Qout
```

扩散部分可望远镜闭合；残差反映 `-sum(u*Dx(T)+v*Dy(T))`，可能来自非零离散速度散度及对流离散。输出散度 RMS、最大值及以 `Uff/L` 归一化的 RMS，供后续与低 Mach、网格收敛联合判断。该残差不是稳态判据，也不证明统计收敛。

Nu 对流和导热项按同一真实温差归一化；修改速度显示标度不改变 Nu。稳态温度误差按 `mean(abs(T-Tprev))/DeltaT`，避免温标平移改变停止标准。

## 与论文结果的比较边界

原文的 `Ra=10^12` 稳定性专项采用纵横比4的二维侧热高腔、`3201×801` 网格、`omega=1.9795`，稳定性对照时段为 `10^5` 次迭代。作者将其定位为 pDNS，不能解释为任意粗网格上充分分辨 DNS 或长期统计收敛的保证。本程序默认方腔，RB 为按用户要求增加的封闭腔体扩展。

当前交付属于公式与项目边界组织的可核验实现。只有完成一致几何、特征尺度、低 Mach/时间步及网格检查后，才能评价其是否重现论文的 Nu、速度和统计量。
