# Zhang, Zhou & Sun (2017) 二维 RB 耗散统计：算例与代码输出要求核对

核对日期：2026-08-03  
核对范围：仅使用仓库内原始论文 `pdf/Statistics of kinetic and thermal energy dissipation rates in two-dimensional turbulent Rayleigh–Bénard convection.pdf`。  
文献：Y. Zhang, Q. Zhou & C. Sun, *Journal of Fluid Mechanics* **814** (2017), 165-184, DOI `10.1017/jfm.2017.19`。

> 页码说明：下文 `p.` 指论文印刷页码；印刷 p. 165 对应 PDF 第 1 页。

## 1. 适用边界：论文只能直接核验下壁加热 RB 分支

论文算例是二维、均匀网格、单位宽高比方腔：

$$
\Gamma=\frac{L}{H}=1.
$$

底板为恒温热壁 $\theta=+0.5$，顶板为恒温冷壁 $\theta=-0.5$；左右侧壁绝热；全部固壁均不可穿透、无滑移。论文采用高度 $H$、温差 $\Delta$ 和自由落体速度 $U=\sqrt{\beta g\Delta H}$ 无量纲化，因此

$$
\nu=\sqrt{\frac{Pr}{Ra}},
\qquad
\kappa=\frac{1}{\sqrt{PrRa}}.
$$

证据：p. 168，§2，Eqs. (2.1)-(2.4) 前后的几何、边界和无量纲化说明。

**适用性结论：**该论文可直接作为“下热上冷、侧壁绝热”分支的物理统计基准；它没有侧壁差温加热算例。Eqs. (1.5)、(1.6) 的精确全局关系依赖这里的 RB 热边界和定义，不能不加推导地宣称也适用于侧壁加热方腔。

## 2. 全部模拟参数和最低 `Ra`

论文对每个 $Pr=0.7$ 和 $Pr=5.3$ 都计算了以下 9 个 Rayleigh 数：

$$
Ra\in\left\{
10^6,\ 3\times10^6,\ 10^7,\ 3\times10^7,\ 10^8,\ 3\times10^8,
\ 10^9,\ 3\times10^9,\ 10^{10}
\right\}.
$$

两组均使用上述单位宽高比二维方腔和均匀网格。网格分辨率如下（p. 169，Table 1）：

| $Ra$ | $N_x\times N_z$，$Pr=0.7$ | $N_x\times N_z$，$Pr=5.3$ |
|---:|---:|---:|
| $1\times10^6$ | $129\times129$ | $129\times129$ |
| $3\times10^6$ | $193\times193$ | $193\times193$ |
| $1\times10^7$ | $257\times257$ | $257\times257$ |
| $3\times10^7$ | $385\times385$ | $385\times385$ |
| $1\times10^8$ | $513\times513$ | $513\times513$ |
| $3\times10^8$ | $769\times769$ | $769\times769$ |
| $1\times10^9$ | $1025\times1025$ | $1025\times1025$ |
| $3\times10^9$ | $1537\times1537$ | $1537\times1537$ |
| $1\times10^{10}$ | $2049\times2049$ | $3073\times3073$ |

最低 Rayleigh 数是 $Ra=10^6$，但它对应两个同等最低的算例，而不是唯一算例：

| $Pr$ | $Ra$ | $N_x\times N_z$ | $Nu$ | $Re$ | $R_u$ | $R_\theta$ | $N_{BL}$ | $\Delta_g/\eta$ | $\Delta_g/\eta_B$ | $\Delta t/\tau_\eta$ |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.7 | $10^6$ | $129\times129$ | 6.30 | 279 | 1.0041 | 1.0023 | 11 | 0.45 | 0.37 | 0.0069 |
| 5.3 | $10^6$ | $129\times129$ | 6.87 | 38 | 1.0026 | 1.0019 | 10 | 0.17 | 0.38 | 0.0042 |

这里 $R_u$、$R_\theta$ 是下面第 3 节定义的精确关系比值。若 P100 只提交一个最低 `Ra` 算例，选择 `Pr=0.7`（Table 1 第一行）还是 `Pr=5.3` 属于实现决定；论文没有给出二选一的优先级。证据：p. 169，Table 1。

## 3. Table 1 每一列的含义

Table 1 从左到右共 11 列（p. 169）：

1. $Pr$。
2. $Ra$。
3. 水平和竖直方向分辨率 $N_x\times N_z$。
4. 全局 Nusselt 数 $Nu$。
5. Reynolds 数
   $$
   Re=\frac{U_{\mathrm{rms}}H}{\nu},
   \qquad
   U_{\mathrm{rms}}=\sqrt{\langle u^2+w^2\rangle_{V,t}}.
   $$
6. 直接积分的动能耗散与精确关系预测值之比
   $$
   R_u=
   \frac{\langle\epsilon_u\rangle_{V,t}}
        {(Nu-1)/\sqrt{RaPr}}.
   $$
7. 直接积分的热能耗散与精确关系预测值之比
   $$
   R_\theta=
   \frac{\langle\epsilon_\theta\rangle_{V,t}}
        {Nu/\sqrt{RaPr}}.
   $$
8. 热边界层内的网格点数 $N_{BL}$。
9. 网格间距与全局 Kolmogorov 长度之比 $\Delta_g/\eta$，其中
   $$
   \eta=\frac{HPr^{1/2}}{[Ra(Nu-1)]^{1/4}}.
   $$
10. 网格间距与 Batchelor 长度之比 $\Delta_g/\eta_B$，其中
    $$
    \eta_B=\eta Pr^{-1/2}.
    $$
11. 时间步长与 Kolmogorov 时间尺度之比 $\Delta t/\tau_\eta$，其中
    $$
    \tau_\eta=\sqrt{\frac{\nu}{\langle\epsilon_u\rangle}}
    =\sqrt{\frac{Pr}{Nu-1}}.
    $$

论文还报告所有算例满足：热边界层至少 10 个网格点、$\Delta_g\lesssim0.57\eta$、$\Delta_g\lesssim0.48\eta_B$、$\Delta t/\tau_\eta<0.01$，且 CFL 数不大于 0.3。证据：p. 168，§2；p. 169，Table 1 及表下注释。

## 4. 局部耗散定义和精确全局关系

### 4.1 局部定义

动能耗散率为（p. 166，Eq. (1.3)）

$$
\epsilon_u(\boldsymbol{x},t)
=
\frac{\nu}{2}\sum_{i,j}
\left[
\frac{\partial u_j}{\partial x_i}
+
\frac{\partial u_i}{\partial x_j}
\right]^2.
$$

在二维 $(x,z)$、速度 $(u,w)$ 下可展开为

$$
\epsilon_u
=
2\nu\left[\left(\frac{\partial u}{\partial x}\right)^2
+\left(\frac{\partial w}{\partial z}\right)^2\right]
+\nu\left(\frac{\partial u}{\partial z}+\frac{\partial w}{\partial x}\right)^2.
$$

热能耗散率为（p. 166，Eq. (1.4)）

$$
\epsilon_\theta(\boldsymbol{x},t)
=
\kappa\sum_i
\left[\frac{\partial\theta}{\partial x_i}\right]^2
=
\kappa\left[
\left(\frac{\partial\theta}{\partial x}\right)^2
+\left(\frac{\partial\theta}{\partial z}\right)^2
\right].
$$

二维展开式是由 Eq. (1.3)/(1.4) 直接展开，论文正文给出的是求和形式。

### 4.2 有量纲精确关系

空间-时间（或系综）平均满足（p. 166，Eqs. (1.5), (1.6)）

$$
\langle\epsilon_u\rangle_{V,t}
=
\frac{\nu^3}{H^4}(Nu-1)RaPr^{-2},
$$

$$
\langle\epsilon_\theta\rangle_{V,t}
=
\kappa\frac{\Delta^2}{H^2}Nu.
$$

在论文自由落体无量纲化下，它们化为（p. 179，Eqs. (3.10a,b)）

$$
\langle\epsilon_u\rangle_{V,t}
=\frac{Nu-1}{\sqrt{RaPr}},
\qquad
\langle\epsilon_\theta\rangle_{V,t}
=\frac{Nu}{\sqrt{RaPr}}.
$$

论文的全局输运定义为（p. 171，Eqs. (3.1), (3.2)）

$$
Nu=1+\sqrt{PrRa}\,\langle w\theta\rangle_{V,t},
$$

$$
Re=\frac{U_{\mathrm{rms}}H}{\nu},
\qquad
U_{\mathrm{rms}}=\sqrt{\langle u^2+w^2\rangle_{V,t}}.
$$

## 5. 全体积、边界层、主体、羽流和背景的口径

### 5.1 边界层/主体分解

论文把两种全局耗散拆成（p. 175，Eqs. (3.4), (3.5)）

$$
\langle\epsilon_u\rangle_{V,t}
=\langle\epsilon_u\rangle_{V_{BL},t}
+\langle\epsilon_u\rangle_{V_{bulk},t},
$$

$$
\langle\epsilon_\theta\rangle_{V,t}
=\langle\epsilon_\theta\rangle_{V_{BL},t}
+\langle\epsilon_\theta\rangle_{V_{bulk},t}.
$$

等号右侧不是未加权的区域条件均值。论文明确说明：区域平均已经乘以对应区域的体积百分比，所以每一项是对全局总耗散的加权贡献。

黏性/热边界层厚度 $\delta_{u,\mathrm{rms}}$、$\delta_{\theta,\mathrm{rms}}$ 定义为“从墙面到 r.m.s. 速度/温度达到最大值位置的距离”。$\delta_{\theta,\mathrm{rms}}$ 接近 $H/(2Nu)$。由于四壁无滑移，速度边界层包括上下板和两侧壁；由于两侧壁绝热，热边界层仅包括上下板。证据：p. 176，§3.3。

### 5.2 羽流/背景分解

热耗散还被拆成（p. 177，Eq. (3.6)）

$$
\langle\epsilon_\theta\rangle_{V,t}
=\langle\epsilon_\theta\rangle_{V_{pl},t}
+\langle\epsilon_\theta\rangle_{V_{bg},t}.
$$

$V_{pl}$ 包括边界层和主体区中检测到的羽流，$V_{bg}$ 是背景区。主体内某点必须同时满足（p. 177，Eqs. (3.7a,b)）

$$
|\theta(x,z,t)-\langle\theta\rangle_{x,t}|
>c\,\theta_{\mathrm{rms}},
$$

$$
\sqrt{PrRa}\,|w(x,z,t)\theta(x,z,t)|
>c\,Nu,
$$

论文取 $c=1.2$。第二个判据取局部对流热流的绝对值，因为二维流中羽流可能产生负的局部热输运。

## 6. 论文实际统计的量

除 Table 1 外，论文图 1-11 还统计或展示：

- $Nu$、$Re$ 以及瞬时 $\theta,u,w,\epsilon_u,\epsilon_\theta$ 场；
- 全场 $\epsilon_u$、$\epsilon_\theta$ 的 PDF，并分别以各自 r.m.s. 值归一化（Figs. 3(a,b), 4(a,b)）；
- $\log_{10}\epsilon_u$、$\log_{10}\epsilon_\theta$ 的 PDF、均值 $\mu$ 和标准差 $\sigma$，并与对数正态分布比较（Figs. 3(c,d), 4(c,d)）；
- PDF 尾部的拉伸指数拟合（p. 173，Eq. (3.3)）
  $$
  p(Y)=\frac{C}{\sqrt{Y}}\exp(-mY^\alpha),
  \qquad Y=X-X_{mp},
  $$
  其中正文以 $X=\epsilon_\theta/(\epsilon_\theta)_{\mathrm{rms}}$ 解释符号，并将该形式用于两种耗散率；
- 水平和时间平均的垂向剖面 $\langle\epsilon_u\rangle_{x,t}$、$\langle\epsilon_\theta\rangle_{x,t}$（Fig. 5）；
- 两种耗散的 BL/bulk 加权贡献和百分比（Figs. 6, 7）；
- 热耗散的 plume/background 加权贡献和百分比（Figs. 8, 9）；
- 全场平均 $\langle\epsilon_u\rangle_{V,t}$、$\langle\epsilon_\theta\rangle_{V,t}$ 随 $Ra$ 的变化（Fig. 10）；
- 各分区耗散归一化后随 $Re$ 的变化（Fig. 11，Eqs. (3.11)-(3.16)）。

这里所谓“非稳态统计”不是从静止初场到稳态的瞬态响应表；论文研究的是时间依赖湍流在统计状态下的空间-时间统计。

## 7. 时间窗口和独立样本：论文没有给出具体数值

论文没有报告以下信息：

- 进入统计状态的起始时刻；
- 总统计时长或自由落体时间数；
- 采样间隔；
- 快照总数；
- 独立样本数或独立性判据。

论文只提供两项时间相关核查：

1. 对全部算例，$\Delta t/\tau_\eta<0.01$（p. 168；p. 169，Table 1）。
2. 分别用每次模拟的前半段和后半段计算 $Nu$、$Re$，所得相对差均小于 1%（p. 171，§3.1）。

因此，代码新增统计起止步、采样间隔和样本数时，必须把它们记录为新的实现参数，不能称为论文复现值。论文也没有给出可直接照抄的“独立样本数要求”。

## 8. 最低 `Ra` 算例必须输出的最小完备指标

### 8.1 可复核 Table 1 与两条精确关系的最小集合

- 元数据：`HeatingMode`、$Pr$、$Ra$、$\Gamma$、$N_x$、$N_z$、$\Delta t$、边界条件和无量纲口径；
- 全局量：$Nu$、$U_{\mathrm{rms}}$、$Re$；
- 直接积分量：$\langle\epsilon_u\rangle_{V,t}$、$\langle\epsilon_\theta\rangle_{V,t}$；
- 精确关系预测值：$(Nu-1)/\sqrt{RaPr}$、$Nu/\sqrt{RaPr}$；
- 校验量：$R_u$、$R_\theta$ 及各自相对误差；
- 分辨率量：$N_{BL}$、$\Delta_g/\eta$、$\Delta_g/\eta_B$、$\Delta t/\tau_\eta$、CFL；
- 收敛量：前半段和后半段各自的 $Nu$、$Re$ 及相对差；
- 统计元数据：实际统计起止步、采样间隔、累计样本数。

其中 `HeatingMode` 不是论文 Table 1 列，但同一代码同时支持侧壁/下壁加热时必须显式输出，防止错误套用 RB 精确关系。

### 8.2 若要声称完成论文的耗散统计，还必须保留

- 两种局部耗散的 r.m.s. 与全场直方图/PDF 数据；
- 两种 $\log_{10}\epsilon$ 的 $\mu$、$\sigma$ 和直方图；
- $\langle\epsilon_u\rangle_{x,t}$、$\langle\epsilon_\theta\rangle_{x,t}$ 的垂向剖面；
- $\delta_{u,\mathrm{rms}}$、$\delta_{\theta,\mathrm{rms}}$；
- $\epsilon_u$、$\epsilon_\theta$ 的 BL/bulk 加权贡献和百分比；
- $\epsilon_\theta$ 的 plume/background 加权贡献和百分比；
- 羽流判据所需的 $\langle\theta\rangle_{x,t}$、$\theta_{\mathrm{rms}}$ 和瞬时 $w\theta$。

最低 `Ra` 的数值锚点应从第 2 节两行中按实际选择的 $Pr$ 取值。仅输出 Table 1 风格的全局量，只能声称完成全局与分辨率核验；没有 PDF、垂向剖面和分区贡献，就不能声称完整复现了论文的非稳态耗散统计。
