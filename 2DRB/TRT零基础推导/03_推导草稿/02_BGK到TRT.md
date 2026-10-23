# 第 2 章：从“向平衡靠拢”一步一步走到 TRT

> 本章目标不是背公式，而是回答五个最基础的问题：什么叫松弛？BGK 的一行碰撞式从哪里来？为什么真正进入输运系数的是 `tau-1/2`？偶部和奇部为什么各有一个 `1/2`？TRT 又怎样在两个碰撞率相等时退化回 BGK？

## 0. 先固定本章记号

上一章已经介绍 population、碰撞、迁移和反向方向。本章统一使用下面的记号：

| 记号 | 本章含义 | 不是什么 |
| --- | --- | --- |
| $h_i$ | 碰撞前、第 $i$ 个离散方向的 population | 不是单颗真实粒子 |
| $h_i^{\mathrm{eq}}$ | 由当前宏观量算出的第 $i$ 个方向平衡值 | 不是整个计算结束后的最终答案 |
| $h_i^*$ | 在同一格点完成碰撞、尚未迁移的值 | 星号不是乘方 |
| $\bar i$ | 与 $i$ 方向完全相反的方向编号 | 不是复共轭 |
| $r$ | 离散碰撞率，英文 relaxation rate | 不是连续时间本身 |
| $\tau=1/r$ | 以一个格子时间步为单位的离散松弛时间 | 不是最终有效输运时间尺度 |
| $H=\tau-1/2=1/r-1/2$ | Hénon 参数，本文也称 Hénon shift | 不是几何高度 |
| $r_+$、$r_-$ | TRT 的偶通道、奇通道碰撞率 | 下标正负不是率的正负号 |

这里的 $i$ 只是方向编号。例如 D1Q3 中可取 $i=0,1,2$；D2Q5 中可取 $i=0,1,2,3,4$。

本章先使用通用 population 字母 $h$。后面回到四篇论文时：流场 population 通常写 $f_i$，温度 population 通常写 $g_i$。换字母不改变本章代数。

---

## 1. 什么叫“松弛到平衡”

**推导身份：** 本节是主稿从一阶松弛 ODE 所作的独立教学推演，不是四篇论文的逐式原文或历史推导；它只为建立“偏差—平衡—离散率”的直觉。

### 1.1 先不用 LBM，只看一个数

设有一个随时间变化的量 $q(t)$：

- $t$ 是连续时间；
- $q^{\mathrm{eq}}$ 是它希望靠近的平衡值；
- $\theta>0$ 是连续松弛时间。

最简单的一阶松弛方程是

\[
\frac{\mathrm dq}{\mathrm dt}
=-\frac{1}{\theta}\left(q-q^{\mathrm{eq}}\right).
\tag{2.1}
\]

先逐项读式 (2.1)。

1. $q-q^{\mathrm{eq}}$ 是“当前值减平衡值”，也就是偏离平衡的量。
2. 前面的负号表示：若 $q$ 高于平衡值，导数为负，$q$ 应下降；若 $q$ 低于平衡值，导数为正，$q$ 应上升。
3. $\theta$ 越大，右端绝对值越小，靠近平衡越慢。

取具体数字 $q=0.20$、$q^{\mathrm{eq}}=0.50$、$\theta=2$，则

\[
q-q^{\mathrm{eq}}=0.20-0.50=-0.30,
\]

\[
\frac{\mathrm dq}{\mathrm dt}
=-\frac12(-0.30)=0.15.
\]

导数为正，符合“当前值低于 0.50，所以应向上走”的直觉。

### 1.2 把连续方程解出来

为避免每次都写长串，定义偏差

\[
d(t)\equiv q(t)-q^{\mathrm{eq}}.
\tag{2.2}
\]

本节暂时假设 $q^{\mathrm{eq}}$ 不随时间变化，所以

\[
\frac{\mathrm dd}{\mathrm dt}=\frac{\mathrm dq}{\mathrm dt}.
\]

将式 (2.2) 代入式 (2.1)：

\[
\frac{\mathrm dd}{\mathrm dt}=-\frac d\theta.
\tag{2.3}
\]

先分两种初值：

- 若 $d(0)=0$，式 (2.3) 给出 $\mathrm dd/\mathrm dt=0$，因此 $d(t)=0$ 是一直停在平衡态的合法解；
- 以下处理 $d(0)\ne0$ 的分支，此时才做除以 $d$ 的操作。

对非零分支，两边除以 $d$ 再乘 $\mathrm dt$：

\[
\frac{\mathrm dd}{d}=-\frac{\mathrm dt}{\theta}.
\]

从初始时刻 $0$ 积分到时刻 $t$：

\[
\ln|d(t)|-\ln|d(0)|=-\frac t\theta.
\]

其中 $|d|$ 表示偏差的绝对值，$\ln$ 表示自然对数。对上式两边取指数，先得到的是

\[
|d(t)|=|d(0)|\exp\left(-\frac t\theta\right).
\]

令 $s$ 表示初始偏差的符号：$d(0)>0$ 时 $s=+1$，$d(0)<0$ 时 $s=-1$。指数因子始终为正，所以该连续衰减分支保持初始符号。于是

\[
d(t)=s|d(t)|
=s|d(0)|\exp\left(-\frac t\theta\right)
=d(0)\exp\left(-\frac t\theta\right).
\]

因此可以写成

\[
d(t)=d(0)\exp\left(-\frac t\theta\right).
\tag{2.4}
\]

其中 $\exp(x)=e^x$，$e\approx2.71828$。把 $d=q-q^{\mathrm{eq}}$ 放回去：

\[
q(t)=q^{\mathrm{eq}}
+\left[q(0)-q^{\mathrm{eq}}\right]
\exp\left(-\frac t\theta\right).
\tag{2.5}
\]

因为 $\exp(-t/\theta)$ 对有限 $t$ 始终为正，若初始偏差 $d(0)\ne0$，连续解的偏差只会缩小而不会改变正负号。后面离散更新在 $r>1$ 时可以越过平衡，正是离散一步更新与连续指数衰减的一个区别。

仍取 $q(0)=0.20$、$q^{\mathrm{eq}}=0.50$、$\theta=2$。在 $t=2$ 时，

\[
q(2)=0.50-0.30e^{-1}
\approx0.50-0.11036=0.38964.
\]

它确实从 0.20 向 0.50 靠近，但还没有一下跳到 0.50。

### 1.3 从连续时间变成一步离散更新

计算机每次只能前进一步。先把一维时间 Taylor 展开写出：

\[
q(t+\Delta t)
=q(t)+\Delta t\,q_t(t)
+\frac{\Delta t^2}{2}q_{tt}(t)
+O(\Delta t^3).
\]

其中 $q_t=\mathrm dq/\mathrm dt$，$q_{tt}=\mathrm d^2q/\mathrm dt^2$。两边减 $q(t)$，再除以非零的 $\Delta t$：

\[
\frac{q(t+\Delta t)-q(t)}{\Delta t}
=q_t(t)+\frac{\Delta t}{2}q_{tt}(t)
+O(\Delta t^2).
\]

若把 $q^*$ 记作这一步 ODE 松弛更新后的 $q(t+\Delta t)$，丢掉首个修正项 $\Delta t\,q_{tt}/2=O(\Delta t)$ 及更高项，就得到最简单的一阶前向差分：

\[
\frac{\mathrm dq}{\mathrm dt}
\approx\frac{q^*-q}{\Delta t}.
\tag{2.6}
\]

这里：

- $q$ 是这一步更新前的值；
- $q^*$ 是完成这次“松弛更新”后的值；
- $\Delta t$ 是一个时间步的长度。

把式 (2.6) 代入式 (2.1)：

\[
\frac{q^*-q}{\Delta t}
=\frac{q^{\mathrm{eq}}-q}{\theta}.
\]

两边乘 $\Delta t$：

\[
q^*-q=\frac{\Delta t}{\theta}
\left(q^{\mathrm{eq}}-q\right).
\]

定义离散碰撞率

\[
r\equiv\frac{\Delta t}{\theta},
\tag{2.7}
\]

就得到

\[
q^*=q+r\left(q^{\mathrm{eq}}-q\right).
\tag{2.8}
\]

在格子单位中通常令 $\Delta t=1$。再定义

\[
\tau\equiv\frac1r,
\tag{2.9}
\]

则 $\tau$ 表示以一个格子时间步为单位的离散松弛时间。

### 1.4 为什么 $r>1$ 会“越过”平衡值

从式 (2.8) 两边减去 $q^{\mathrm{eq}}$：

\[
q^*-q^{\mathrm{eq}}
=q-q^{\mathrm{eq}}
+r(q^{\mathrm{eq}}-q).
\]

由于 $q^{\mathrm{eq}}-q=-(q-q^{\mathrm{eq}})$，所以

\[
q^*-q^{\mathrm{eq}}
=(1-r)(q-q^{\mathrm{eq}}).
\tag{2.10}
\]

这一个式子把松弛分成三种直观情况：

| $r$ | 偏差倍率 $1-r$ | 一步发生什么 |
| ---: | ---: | --- |
| $0<r<1$ | $0<1-r<1$ | 从原来一侧靠近平衡，不越过 |
| $r=1$ | $1-r=0$ | 一步正好到平衡 |
| $1<r<2$ | $-1<1-r<0$ | 越过平衡，但偏差绝对值变小 |

例如 $q=0.20$、$q^{\mathrm{eq}}=0.50$、$r=1.2$：

\[
q^*=0.20+1.2(0.50-0.20)=0.56.
\]

新值 0.56 越过了 0.50。偏差由 $-0.30$ 变为 $+0.06$，而

\[
+0.06=(1-1.2)(-0.30)=(-0.2)(-0.30).
\]

这就是 over-relaxation，中文常译“过松弛”。“过”指越过平衡，不表示计算必然出错。

对单个标量松弛模态，若希望每一步偏差绝对值缩小，需要

\[
|1-r|<1.
\]

逐步解不等式：

\[
-1<1-r<1,
\]

\[
-2<-r<0,
\]

乘以 $-1$ 时不等号反向：

\[
0<r<2.
\tag{2.11}
\]

Dubois 与 Lallemand 在其一般矩方法中也把 $0<s_k<2$ 列为松弛参数的稳定性约束之一。[Dubois–Lallemand 2009，PDF/印刷 p. 4，Eq. (12) 前]

但是必须分清：式 (2.11) 只证明一个孤立线性模态的局部放大因子小于 1。完整 LBM 还包括迁移、边界、多个模态、非线性、离散速度和浮点误差，所以 $0<r<2$ 不是“任意算例都稳定”的充分证明。

---

## 2. 从一个数推广到每个 population：BGK

**推导身份：** 本节是主稿把标量更新推广到 populations 后所作的独立教学演算；Dubois–Lallemand 的对应论文证据会在具体公式后标出。

### 2.1 BGK 的名字与作用

BGK 来自 Bhatnagar、Gross、Krook 三位作者的姓氏。对初学者，本章只需要抓住它的离散作用：**所有需要松弛的 population 都使用同一个碰撞率 $r$。**

把式 (2.8) 中的标量 $q$ 换成每个方向的 $h_i$：

\[
h_i^*=h_i+r\left(h_i^{\mathrm{eq}}-h_i\right).
\tag{2.12}
\]

等价地，也可以写成

\[
h_i^*=h_i-r\left(h_i-h_i^{\mathrm{eq}}\right).
\tag{2.13}
\]

两式完全相同，因为

\[
h_i^{\mathrm{eq}}-h_i
=-(h_i-h_i^{\mathrm{eq}}).
\]

### 2.2 用三个 population 完整算一次

取某格点的三个 population：

\[
(h_0,h_1,h_2)=(0.60,0.30,0.10).
\]

假设由当前宏观量构造出的平衡值为

\[
(h_0^{\mathrm{eq}},h_1^{\mathrm{eq}},h_2^{\mathrm{eq}})
=(0.50,0.25,0.25).
\]

先检查更新前和平衡态的总和：

\[
0.60+0.30+0.10=1.00,
\]

\[
0.50+0.25+0.25=1.00.
\]

二者总和相同，这一条件将保证碰撞不改变该总和。

取 $r=1.2$，逐分量计算：

\[
h_0^*=0.60+1.2(0.50-0.60)=0.48,
\]

\[
h_1^*=0.30+1.2(0.25-0.30)=0.24,
\]

\[
h_2^*=0.10+1.2(0.25-0.10)=0.28.
\]

碰撞后总和为

\[
0.48+0.24+0.28=1.00.
\]

再检查每个方向相对平衡态的偏差：

\[
(h_0-h_0^{\mathrm{eq}},h_1-h_1^{\mathrm{eq}},h_2-h_2^{\mathrm{eq}})
=(0.10,0.05,-0.15),
\]

\[
(h_0^*-h_0^{\mathrm{eq}},h_1^*-h_1^{\mathrm{eq}},h_2^*-h_2^{\mathrm{eq}})
=(-0.02,-0.01,0.03).
\]

后一行正好是前一行乘 $1-r=-0.2$。三个分量都越过平衡，但偏差绝对值缩小到原来的五分之一。

### 2.3 一般地证明碰撞守恒总和

对式 (2.12) 的所有方向求和：

\[
\sum_i h_i^*
=\sum_i h_i
+r\sum_i\left(h_i^{\mathrm{eq}}-h_i\right).
\tag{2.14}
\]

若平衡态的构造满足

\[
\sum_i h_i^{\mathrm{eq}}=\sum_i h_i,
\tag{2.15}
\]

则

\[
\sum_i\left(h_i^{\mathrm{eq}}-h_i\right)=0.
\]

代回式 (2.14)：

\[
\boxed{\sum_i h_i^*=\sum_i h_i}.
\tag{2.16}
\]

所以“BGK 碰撞守恒”不是因为碰撞项凭空消失，而是因为平衡态和当前态被设计为具有相同的守恒矩。

在 Dubois 与 Lallemand 的矩表述中，非守恒矩使用

\[
m_k^*=m_k+s_k\left(m_k^{\mathrm{eq}}-m_k\right),
\]

而若所有非守恒矩的 $s_k$ 都相同，矩碰撞就退化为单松弛率 BGK。[Dubois–Lallemand 2009，PDF/印刷 p. 4，Eqs. (12)–(14) 及其后说明]

---

## 3. 为什么输运系数里总出现 $\tau-1/2$

这是初学者最容易被一句“众所周知”跳过去的地方。本节不把 $-1/2$ 当作经验修正，而从“碰撞后还要整步迁移”的二阶 Taylor 展开中把它算出来。

**推导身份：** 下面是主稿基于碰撞—迁移定义所作的独立 Taylor 教学推演，用来复核 Dubois–Lallemand Eq. (16) 的结构；它不是该论文的逐字或逐式历史推导。

先把本节推导边界写清楚：

1. 暂时没有外力项和热源项；若加入源项，它自己的时间离散也可能带来半步修正，必须重新分析。
2. 使用一个守恒标量和空间、时间上均取常数的统一碰撞率 $r$；本节只解释 Hénon shift 的时间离散来源。
3. 假设 population 沿每条特征方向至少有连续且有界的四阶导数，使显式三阶项和 $O(\Delta t^4)$ 余项有意义。
4. 这是一段体相、光滑区域的等效方程推导，不包含壁面闭合，因此不能从本节直接推出 boundary magic。

还有一个与前面 ODE 类比不同的计阶约定：从式 (2.17) 开始，我们研究一族形式上的 LBM 离散格式，令 $\Delta t$ 变小，同时把无量纲碰撞率 $r=O(1)$ 和格子速度比例保持固定。若硬要换回第 1 节的 ODE 类比，则 $\theta=\Delta t/r$ 会随网格一起变化；这里不是“固定连续 $\theta$ 的前向 Euler 收敛证明”。

### 3.1 先写完整碰撞—迁移方程

把碰撞和迁移合在一行：

\[
h_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
-h_i(\boldsymbol x,t)
=r\left[h_i^{\mathrm{eq}}(\boldsymbol x,t)-h_i(\boldsymbol x,t)\right].
\tag{2.17}
\]

左边不是只在时间上前进一步，还沿 $\boldsymbol c_i$ 在空间上移动一步。

为缩短记号，定义沿第 $i$ 条特征方向的微分算子

\[
D_i\equiv\partial_t+\boldsymbol c_i\cdot\boldsymbol\nabla.
\tag{2.18}
\]

逐项解释：

- $\partial_t$ 看时间变化；
- $\boldsymbol c_i\cdot\boldsymbol\nabla$ 看沿第 $i$ 个速度方向的空间变化；
- $D_i$ 把二者合在一起。

$D_i^2h_i$ 表示把同一个算子连续作用两次，即 $D_i(D_i h_i)$；$D_i^3h_i=D_i(D_i^2h_i)$ 表示连续作用三次。它们不是普通数字 $D_i$ 的平方或立方。

在一维、速度 $c_i=c$ 时，

\[
D_i=\partial_t+c\partial_x.
\]

### 3.2 左端 Taylor 展开到二阶

假设 $h_i$ 对时间和空间足够光滑。沿 $(\boldsymbol c_i\Delta t,\Delta t)$ 展开。完整写到下一项是

\[
h_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
=h_i+\Delta tD_i h_i
+\frac{\Delta t^2}{2}D_i^2h_i
+\frac{\Delta t^3}{6}D_i^3h_i
+O(\Delta t^4).
\]

因为本节只保留到二阶，第一个被舍弃的显式项是 $\Delta t^3D_i^3h_i/6$。把它和更高项合写为 $O(\Delta t^3)$，得到

\[
h_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
=h_i
+\Delta t D_i h_i
+\frac{\Delta t^2}{2}D_i^2h_i
+O(\Delta t^3).
\tag{2.19}
\]

式中的 $1/2$ 来自 Taylor 展开的 $2!=2$，不是碰撞模型另行指定的参数。

将式 (2.19) 代入式 (2.17)，左边的 $h_i$ 抵消：

\[
\Delta t D_i h_i
+\frac{\Delta t^2}{2}D_i^2h_i
+O(\Delta t^3)
=r(h_i^{\mathrm{eq}}-h_i).
\tag{2.20}
\]

### 3.3 把 population 写成“平衡态 + 小偏差”

当 $\Delta t$ 足够小，引入与 $\Delta t$ 无关并在极限中保持有界的系数场 $h_i^{(1)},h_i^{(2)},h_i^{(3)}$。上标括号只是阶数编号，不是幂。若 $\Delta t$ 保留时间单位，则 $\Delta t h_i^{(1)}$、$\Delta t^2h_i^{(2)}$ 必须与 $h_i$ 同量纲，所以 $h_i^{(1)}$、$h_i^{(2)}$ 的量纲分别是 population/时间、population/时间平方。

完整多写一项是

\[
h_i=h_i^{\mathrm{eq}}
+\Delta t h_i^{(1)}
+\Delta t^2h_i^{(2)}
+\Delta t^3h_i^{(3)}
+O(\Delta t^4).
\]

本节保留到二阶，把首个舍弃项 $\Delta t^3h_i^{(3)}$ 与更高项合并，得到

\[
h_i=h_i^{\mathrm{eq}}
+\Delta t h_i^{(1)}
+\Delta t^2 h_i^{(2)}
+O(\Delta t^3).
\tag{2.21}
\]

式 (2.21) 中的 $O(\Delta t^3)$ 就包含刚才显式写出的 $\Delta t^3h_i^{(3)}$。

由式 (2.21)，

\[
h_i^{\mathrm{eq}}-h_i
=-\Delta t h_i^{(1)}
-\Delta t^2 h_i^{(2)}
+O(\Delta t^3).
\tag{2.22}
\]

### 3.4 比较一阶项，求第一阶非平衡量

把式 (2.21)–(2.22) 放入式 (2.20)。只收集含一个 $\Delta t$ 的项：

\[
\Delta t D_i h_i^{\mathrm{eq}}
=-r\Delta t h_i^{(1)}.
\]

两边除以非零的 $\Delta t$：

\[
D_i h_i^{\mathrm{eq}}=-r h_i^{(1)}.
\]

若 $r\ne0$，则

\[
h_i^{(1)}=-\frac1rD_i h_i^{\mathrm{eq}}.
\tag{2.23}
\]

这一步产生了 $1/r=\tau$，但还没有产生 $-1/2$。后者来自下一步保留的二阶迁移项。

### 3.5 对守恒方程求和，$-1/2$ 出现

设本节考察的守恒量为

\[
\rho=\sum_i h_i=\sum_i h_i^{\mathrm{eq}}.
\tag{2.24}
\]

平衡态承担全部守恒量，因此展开中的非平衡修正满足

\[
\sum_i h_i^{(1)}=0,
\qquad
\sum_i h_i^{(2)}=0.
\tag{2.25}
\]

对式 (2.20) 的所有 $i$ 求和。右边为

\[
r\sum_i(h_i^{\mathrm{eq}}-h_i)=0.
\]

左边除以 $\Delta t$，并保留到一阶修正：

\[
0=\sum_iD_i h_i^{\mathrm{eq}}
+\Delta t\sum_iD_i h_i^{(1)}
+\frac{\Delta t}{2}\sum_iD_i^2h_i^{\mathrm{eq}}
+O(\Delta t^2).
\tag{2.26}
\]

将式 (2.23) 代入中间一项：

\[
\sum_iD_i h_i^{(1)}
=-\frac1r\sum_iD_i^2h_i^{\mathrm{eq}}.
\]

于是式 (2.26) 变成

\[
0=\sum_iD_i h_i^{\mathrm{eq}}
+\Delta t\left(\frac12-\frac1r\right)
\sum_iD_i^2h_i^{\mathrm{eq}}
+O(\Delta t^2).
\tag{2.27}
\]

把括号中的负号提出：

\[
\frac12-\frac1r
=-\left(\frac1r-\frac12\right).
\]

因此自然定义

\[
\boxed{
H\equiv\frac1r-\frac12
=\tau-\frac12
}.
\tag{2.28}
\]

式 (2.27) 就是

\[
0=\sum_iD_i h_i^{\mathrm{eq}}
-\Delta t H\sum_iD_i^2h_i^{\mathrm{eq}}
+O(\Delta t^2).
\tag{2.29}
\]

现在可以看清两个来源：

- $1/r$ 来自碰撞产生的第一阶非平衡量，见式 (2.23)；
- $1/2$ 来自整步迁移 Taylor 展开的二阶项，见式 (2.19)。

两者相减，才是进入二阶等效输运项的时间尺度。

### 3.6 用最小 D1Q2 扩散例把 $H$ 变成扩散系数

为了只观察 $-1/2$，取一维两速度模型：

\[
c_+=+c,\qquad c_-=-c.
\]

取对称平衡态

\[
h_+^{\mathrm{eq}}=\frac\rho2,
\qquad
h_-^{\mathrm{eq}}=\frac\rho2.
\tag{2.30}
\]

两个方向算子为

\[
D_+=\partial_t+c\partial_x,
\qquad
D_-=\partial_t-c\partial_x.
\]

先算式 (2.29) 的第一项：

\[
\begin{aligned}
\sum_{i=+,-}D_i h_i^{\mathrm{eq}}
&=(\partial_t+c\partial_x)\frac\rho2
+(\partial_t-c\partial_x)\frac\rho2\\
&=\frac12\partial_t\rho+\frac c2\partial_x\rho
+\frac12\partial_t\rho-\frac c2\partial_x\rho\\
&=\partial_t\rho.
\end{aligned}
\tag{2.31}
\]

再把平方算子完全展开：

\[
D_+^2
=(\partial_t+c\partial_x)^2
=\partial_t^2+2c\partial_t\partial_x+c^2\partial_x^2,
\]

\[
D_-^2
=(\partial_t-c\partial_x)^2
=\partial_t^2-2c\partial_t\partial_x+c^2\partial_x^2.
\]

作用于两个 $\rho/2$ 并相加：

\[
\begin{aligned}
\sum_{i=+,-}D_i^2 h_i^{\mathrm{eq}}
&=\frac12(D_+^2+D_-^2)\rho\\
&=\partial_t^2\rho+c^2\partial_x^2\rho.
\end{aligned}
\tag{2.32}
\]

其中两个混合项 $+2c\partial_t\partial_x$ 与 $-2c\partial_t\partial_x$ 正好抵消。

把式 (2.31)–(2.32) 代入式 (2.29)：

\[
\partial_t\rho
=\Delta t H
\left(\partial_t^2\rho+c^2\partial_x^2\rho\right)
+O(\Delta t^2).
\tag{2.33}
\]

为什么可以暂时略去 $\partial_t^2\rho$？引入一个只用于标记阶数的无量纲小参数 $\varepsilon$，把 $\Delta t$ 记作 $O(\varepsilon)$，再定义慢时间 $T=\varepsilon t$。链式法则给

\[
\partial_t=\varepsilon\partial_T,
\qquad
\partial_t^2=\varepsilon^2\partial_T^2.
\]

所以 $\partial_t\rho=O(\varepsilon)$，而 $\partial_t^2\rho=O(\varepsilon^2)$；后者在式 (2.33) 右端又乘一个 $\Delta t=O(\varepsilon)$，贡献为 $O(\varepsilon^3)$。这正是本小节采用的扩散慢时间计阶，并不是说任何快速瞬态都能删除 $\partial_t^2\rho$。保留最低非零扩散阶：

\[
\partial_t\rho
=\kappa\partial_x^2\rho+\text{更高阶项},
\]

其中

\[
\boxed{\kappa=c^2\Delta t H
=c^2\Delta t\left(\frac1r-\frac12\right)}.
\tag{2.34}
\]

这里 $\kappa$ 是这个 D1Q2 教学模型在慢时间极限下的扩散系数。

这个 D1Q2 只是为展示 $-1/2$ 的来源而搭建的最小课堂模型。D2Q5、D2Q9 的具体扩散率或黏性系数还会乘各自平衡矩和速度集给出的因子，不能把式 (2.34) 的 $c^2$ 原样复制过去；但 Hénon 组合 $1/r-1/2$ 会再次出现。

Dubois 与 Lallemand 的 Eq. (16) 正是定义

\[
\sigma_k\equiv\frac1{s_k}-\frac12,
\]

并在后文写出二阶输运系数与 $\Delta t\sigma_k$ 成正比。[Dubois–Lallemand 2009，PDF/印刷 p. 5，Eqs. (16)–(17) 后]

### 3.7 $r$、$\tau$、$H$ 三联换算

三条定义是

\[
\tau=\frac1r,
\qquad
H=\tau-\frac12,
\qquad
r=\frac1{H+1/2}=\frac2{2H+1}.
\tag{2.35}
\]

最后一条逐步来自

\[
H=\frac1r-\frac12,
\]

\[
H+\frac12=\frac1r,
\]

\[
r=\frac1{H+1/2}.
\]

倒数换算要求 $r\ne0$；从 $H$ 反求 $r$ 时还要求 $H\ne-1/2$，否则分母为零。

数值检查：

| $r$ | $\tau=1/r$ | $H=\tau-1/2$ | 偏差倍率 $1-r$ |
| ---: | ---: | ---: | ---: |
| $1/2$ | $2$ | $3/2$ | $1/2$ |
| $1$ | $1$ | $1/2$ | $0$ |
| $6/5$ | $5/6$ | $1/3$ | $-1/5$ |
| $3/2$ | $2/3$ | $1/6$ | $-1/2$ |

例如 $H=1/3$ 时，

\[
r=\frac1{1/3+1/2}
=\frac1{5/6}=\frac65=1.2,
\]

与表中一致。

若只考虑式 (2.11) 的 $0<r<2$，则

\[
\tau=\frac1r>\frac12,
\qquad H>0.
\]

但后续读论文时不能看见任意一个 $\sigma$ 就自动把它当 Hénon 参数；不同论文复用了同一个希腊字母，含义可能相反。

---

## 4. 偶部和奇部：TRT 的两个抽屉

**推导身份：** 本节从和与差的定义独立推演偶奇投影，目的是让每个 $1/2$ 和反向变号都可手算，不把它冒充四篇论文的原文段落。

### 4.1 先认清反向方向

对一个非静止方向 $i$，用 $\bar i$ 表示它的反向方向：

\[
\boldsymbol c_{\bar i}=-\boldsymbol c_i.
\tag{2.36}
\]

反向操作做两次必须回到原方向：

\[
\bar{\bar i}=i.
\]

例如 D2Q5 中，若 $i=1$ 是东向，则 $\bar i=3$ 是西向；若 $i=2$ 是北向，则 $\bar i=4$ 是南向。

静止方向满足 $\boldsymbol c_0=\boldsymbol0$，所以它反向后仍是自己，通常写 $\bar0=0$。

### 4.2 用两个具体数字算偶部和奇部

取一对反向 population：

\[
h_i=0.30,
\qquad h_{\bar i}=0.10.
\]

先算和与差：

\[
h_i+h_{\bar i}=0.40,
\qquad
h_i-h_{\bar i}=0.20.
\]

定义偶部为和的一半：

\[
h_i^+\equiv\frac{h_i+h_{\bar i}}2
=\frac{0.40}{2}=0.20.
\tag{2.37}
\]

定义奇部为差的一半：

\[
h_i^-\equiv\frac{h_i-h_{\bar i}}2
=\frac{0.20}{2}=0.10.
\tag{2.38}
\]

这里的 $1/2$ 只因为要从“和、差”恢复原来的两个数。它不是 Hénon 参数中的 $-1/2$，两者来源完全不同。

重构第 $i$ 个方向：

\[
h_i^++h_i^-=0.20+0.10=0.30=h_i.
\tag{2.39}
\]

重构反向方向：

\[
h_i^+-h_i^-=0.20-0.10=0.10=h_{\bar i}.
\tag{2.40}
\]

因此一般地

\[
\boxed{h_i=h_i^++h_i^-},
\qquad
\boxed{h_{\bar i}=h_i^+-h_i^-}.
\tag{2.41}
\]

### 4.3 为什么叫“偶”和“奇”

交换 $i$ 与 $\bar i$。反向方向的偶部是

\[
h_{\bar i}^+
=\frac{h_{\bar i}+h_i}{2}
=h_i^+.
\tag{2.42}
\]

它在反向操作下不变，所以叫偶部。

反向方向的奇部是

\[
h_{\bar i}^-
=\frac{h_{\bar i}-h_i}{2}
=-\frac{h_i-h_{\bar i}}2
=-h_i^-.
\tag{2.43}
\]

它在反向操作下变号，所以叫奇部。

静止方向是自己的反向方向。令 $\bar0=0$ 代入定义：

\[
h_0^+=\frac{h_0+h_0}{2}=h_0,
\qquad
h_0^-=\frac{h_0-h_0}{2}=0.
\]

所以静止 population 只有偶部，不能再把它与自己重复组成一个“双元素方向对”。

平衡态也用完全相同的方法分解：

\[
h_i^{\mathrm{eq},+}
=\frac{h_i^{\mathrm{eq}}+h_{\bar i}^{\mathrm{eq}}}{2},
\qquad
h_i^{\mathrm{eq},-}
=\frac{h_i^{\mathrm{eq}}-h_{\bar i}^{\mathrm{eq}}}{2}.
\tag{2.44}
\]

下图把“两个原数—和差各除以 2—再相加减恢复”的关系画在一起。

![图 2-1：一对反向 population 的偶奇分解与重构](../04_图表/偶奇分解.svg)

**图 2-1。** 上半部从 $h_i,h_{\bar i}$ 生成 $h_i^+,h_i^-$；下半部用 $h_i^+\pm h_i^-$ 恢复原来的两个方向。读图时不要把投影定义中的 $1/2$ 与 Hénon shift 中的 $-1/2$ 混在一起。

---

## 5. 给两个抽屉不同碰撞率：TRT

**推导身份：** 本节的分量式、两个数字例和 BGK 极限证明是主稿基于偶奇定义所作的独立教学复算；论文原符号将在第 6 节单独映射。

### 5.1 从 BGK 的一个率变成两个率

BGK 对全部偏差使用同一个 $r$。TRT 的全称是 two-relaxation-time，即“两松弛时间”。等价地，它使用两个碰撞率：

- 偶部偏差使用 $r_+$，也记作 $r_{\mathrm{even}}$；
- 奇部偏差使用 $r_-$，也记作 $r_{\mathrm{odd}}$。

先分别写碰撞后的偶部和奇部：

\[
(h_i^+)^*
=h_i^+
+r_+\left(h_i^{\mathrm{eq},+}-h_i^+\right),
\tag{2.45}
\]

\[
(h_i^-)^*
=h_i^-
+r_-\left(h_i^{\mathrm{eq},-}-h_i^-\right).
\tag{2.46}
\]

碰撞后原方向等于碰撞后偶部加奇部：

\[
h_i^*=(h_i^+)^*+(h_i^-)^*.
\]

把式 (2.45)–(2.46) 代入：

\[
\boxed{
h_i^*=h_i
+r_+\left(h_i^{\mathrm{eq},+}-h_i^+\right)
+r_-\left(h_i^{\mathrm{eq},-}-h_i^-\right)
}.
\tag{2.47}
\]

这就是本教程统一使用的 TRT 分量式。

也常写成“减去当前偏差”的形式：

\[
h_i^*=h_i
-r_+\left(h_i^+-h_i^{\mathrm{eq},+}\right)
-r_-\left(h_i^--h_i^{\mathrm{eq},-}\right).
\tag{2.48}
\]

式 (2.47) 与式 (2.48) 完全相同。

### 5.2 反向分量为什么奇通道前会变号

定义

\[
\Delta^+\equiv h_i^{\mathrm{eq},+}-h_i^+,
\qquad
\Delta^-\equiv h_i^{\mathrm{eq},-}-h_i^-.
\]

则第 $i$ 个方向是

\[
h_i^*=h_i+r_+\Delta^++r_-\Delta^-.
\tag{2.49}
\]

对反向方向，偶部不变，奇部变号：

\[
h_{\bar i}^{\mathrm{eq},+}-h_{\bar i}^+=\Delta^+,
\]

\[
h_{\bar i}^{\mathrm{eq},-}-h_{\bar i}^-
=(-h_i^{\mathrm{eq},-})-(-h_i^-)
=-\Delta^-.
\]

所以

\[
h_{\bar i}^*
=h_{\bar i}+r_+\Delta^+-r_-\Delta^-.
\tag{2.50}
\]

现在把式 (2.49) 与式 (2.50) 相加：

\[
h_i^*+h_{\bar i}^*
=h_i+h_{\bar i}+2r_+\Delta^+.
\tag{2.51}
\]

奇通道修正互相抵消，所以偶率控制一对 population 的“和”。

再相减：

\[
h_i^*-h_{\bar i}^*
=h_i-h_{\bar i}+2r_-\Delta^-.
\tag{2.52}
\]

偶通道修正互相抵消，所以奇率控制一对 population 的“差”。这正是 TRT 可以分别调节两类性质的代数原因。

这里还要避免一个误解：式 (2.51) 表明**单独一对** population 的和一般可以被偶通道改变，只有当该对的 $\Delta^+=0$ 时才单独守恒。真正要求的是所有方向的守恒总和。对完整、反向成对的速度集，奇部在每一对中相消，而偶部在全体方向上的和等于原 populations 的和，因此

\[
\sum_i(h_i^{\mathrm{eq},-}-h_i^-)=0,
\]

\[
\sum_i(h_i^{\mathrm{eq},+}-h_i^+)
=\sum_i(h_i^{\mathrm{eq}}-h_i)=0.
\]

后一等号使用了平衡态与当前态共享同一个守恒总和。于是对式 (2.47) 求和后仍有

\[
\sum_i h_i^*=\sum_i h_i.
\]

所以 TRT 守恒的是规定的全体方向总和，而不是无条件守恒每一对的和。

### 5.3 完整算一对 TRT 碰撞

仍用

\[
h_i=0.30,
\qquad h_{\bar i}=0.10,
\]

因此

\[
h_i^+=0.20,
\qquad h_i^-=0.10.
\]

取平衡 population

\[
h_i^{\mathrm{eq}}=0.24,
\qquad h_{\bar i}^{\mathrm{eq}}=0.16.
\]

平衡偶部和奇部为

\[
h_i^{\mathrm{eq},+}
=\frac{0.24+0.16}{2}=0.20,
\]

\[
h_i^{\mathrm{eq},-}
=\frac{0.24-0.16}{2}=0.04.
\]

取

\[
r_+=1,
\qquad r_-=1.5.
\]

偶修正为

\[
r_+(h_i^{\mathrm{eq},+}-h_i^+)
=1(0.20-0.20)=0.
\]

奇修正为

\[
r_-(h_i^{\mathrm{eq},-}-h_i^-)
=1.5(0.04-0.10)=-0.09.
\]

所以

\[
h_i^*=0.30+0-0.09=0.21.
\]

反向方向的奇修正变号：

\[
h_{\bar i}^*=0.10+0+0.09=0.19.
\]

检查一对 population 的和：

\[
0.21+0.19=0.40=0.30+0.10.
\]

新偶部为

\[
\frac{0.21+0.19}{2}=0.20,
\]

与平衡偶部 0.20 相同。

新奇部为

\[
\frac{0.21-0.19}{2}=0.01.
\]

原奇部相对平衡奇部的偏差为

\[
0.10-0.04=0.06.
\]

碰撞后偏差为

\[
0.01-0.04=-0.03.
\]

而

\[
-0.03=(1-r_-)(0.06)=(1-1.5)(0.06).
\]

这再次验证：奇部偏差的倍率是 $1-r_-$。

刚才特意取了 $h_i^+=h_i^{\mathrm{eq},+}=0.20$，所以偶修正为零。它只完整测试了奇通道。下面再用一个 D1Q3 例子让偶、奇两个通道同时工作，并检查全体方向总和。

取静止方向 $0$ 和一对反向方向 $i,\bar i$：

\[
(h_0,h_i,h_{\bar i})
=\left(\frac15,\frac7{10},\frac3{10}\right),
\]

平衡态取

\[
(h_0^{\mathrm{eq}},h_i^{\mathrm{eq}},h_{\bar i}^{\mathrm{eq}})
=\left(\frac25,\frac9{20},\frac7{20}\right).
\]

先检查总和：

\[
\frac15+\frac7{10}+\frac3{10}=\frac65,
\]

\[
\frac25+\frac9{20}+\frac7{20}=\frac65.
\]

当前方向对的偶、奇部为

\[
h_i^+=\frac{7/10+3/10}{2}=\frac12,
\qquad
h_i^-=\frac{7/10-3/10}{2}=\frac15.
\]

平衡方向对的偶、奇部为

\[
h_i^{\mathrm{eq},+}
=\frac{9/20+7/20}{2}=\frac25,
\qquad
h_i^{\mathrm{eq},-}
=\frac{9/20-7/20}{2}=\frac1{20}.
\]

所以

\[
\Delta^+=\frac25-\frac12=-\frac1{10},
\qquad
\Delta^-=\frac1{20}-\frac15=-\frac3{20}.
\]

取

\[
r_+=\frac65,
\qquad
r_-=\frac45.
\]

第 $i$ 个方向为

\[
\begin{aligned}
h_i^*
&=\frac7{10}+\frac65\left(-\frac1{10}\right)
+\frac45\left(-\frac3{20}\right)\\
&=\frac7{10}-\frac3{25}-\frac3{25}
=\frac{23}{50}.
\end{aligned}
\]

反向方向的奇修正变号：

\[
\begin{aligned}
h_{\bar i}^*
&=\frac3{10}+\frac65\left(-\frac1{10}\right)
-\frac45\left(-\frac3{20}\right)\\
&=\frac3{10}-\frac3{25}+\frac3{25}
=\frac3{10}.
\end{aligned}
\]

静止方向只有偶部，因此

\[
h_0^*
=\frac15+\frac65\left(\frac25-\frac15\right)
=\frac{11}{25}.
\]

全体方向总和为

\[
\frac{11}{25}+\frac{23}{50}+\frac3{10}
=\frac{22+23+15}{50}
=\frac65.
\]

再从碰撞后数值反算两个投影：

\[
\frac{23/50+3/10}{2}=\frac{19}{50}
=\frac12+\frac65\left(-\frac1{10}\right),
\]

\[
\frac{23/50-3/10}{2}=\frac2{25}
=\frac15+\frac45\left(-\frac3{20}\right).
\]

偶部从 $1/2$ 变成 $19/50$，奇部从 $1/5$ 变成 $2/25$；两条通道和全体守恒都得到数字验证。

### 5.4 逐项证明 TRT 在两率相等时回到 BGK

令

\[
r_+=r_-=r.
\tag{2.53}
\]

从 TRT 式 (2.47) 出发：

\[
h_i^*=h_i
+r(h_i^{\mathrm{eq},+}-h_i^+)
+r(h_i^{\mathrm{eq},-}-h_i^-).
\]

提出共同因子 $r$：

\[
h_i^*=h_i
+r\left[
h_i^{\mathrm{eq},+}+h_i^{\mathrm{eq},-}
-h_i^+-h_i^-
\right].
\]

根据重构式，

\[
h_i^{\mathrm{eq},+}+h_i^{\mathrm{eq},-}=h_i^{\mathrm{eq}},
\]

\[
h_i^++h_i^-=h_i.
\]

因此

\[
h_i^*=h_i+r(h_i^{\mathrm{eq}}-h_i),
\tag{2.54}
\]

这正是 BGK 式 (2.12)。证明完成。

静止方向也没有漏掉：因为 $h_0^-=h_0^{\mathrm{eq},-}=0$ 且 $h_0^+=h_0$，两率相等时式 (2.47) 给

\[
h_0^*=h_0+r(h_0^{\mathrm{eq}}-h_0),
\]

同样是 BGK。上面的 $i$ 又是任意移动方向，所以它的反向方向也已经包含在证明中。

注意：这不是“两个方法碰巧给出相近结果”，而是每个方向的代数表达式完全相同。

### 5.5 TRT 和 MRT 的关系

MRT 是 multiple-relaxation-time，即“多松弛时间”方法。它先把 populations 变到矩空间，再允许不同非守恒矩使用不同碰撞率。

TRT 可以看成一种有结构的两率限制：

1. 先按中心反演把成对方向分成偶部和奇部；
2. 一组偶模态共享 $r_+$；
3. 一组奇模态共享 $r_-$；
4. 守恒模态仍不被碰撞改变。

所以可以写成关系链：

\[
\text{BGK：一个率}
\subset
\text{TRT：偶、奇两个率}
\subset
\text{MRT：多个模态率}.
\]

这里的“包含”表示碰撞参数自由度的包含，不表示任意一套矩基或任意非线性碰撞都能不加条件地互换。

---

## 6. 四篇论文的符号翻译表

这一节非常重要，因为四篇论文对 $\sigma$、$\lambda$ 和 $\Lambda$ 的用法并不统一。

**表格身份：** 下表是主稿对四个本地论文证据账本的综合翻译，不是任意一篇论文原有的统一符号表。

| 论文 | 论文原符号 | 原符号类型 | 本教程统一记号 | 精确换算 | 原文位置 |
| --- | --- | --- | --- | --- | --- |
| Wang et al. 2013 | $s_\nu,s_q$ | D2Q9 流场碰撞率 | $r_\nu^{(f)},r_q^{(f)}$ | $H_j^{(f)}=1/r_j^{(f)}-1/2$ | PDF p. 3, Eqs. (8)–(10) |
| Wang et al. 2013 | $\sigma_\kappa,\sigma_e,\sigma_\nu$ | D2Q5 温度碰撞率 | $r_\kappa,r_e,r_\nu$ | $H_j=1/r_j-1/2$ | PDF pp. 3–4, Eqs. (13)–(17) |
| Contrino et al. 2014 | $s_\nu,s_q$ | D2Q9 流场碰撞率 | $r_\nu^{(f)},r_q^{(f)}$ | $H_j^{(f)}=1/r_j^{(f)}-1/2$ | PDF p. 3, Eqs. (7a)–(7b) |
| Contrino et al. 2014 | $\sigma_\kappa,\sigma_e,\sigma_\nu$ | D2Q5 温度碰撞率 | $r_\kappa,r_e,r_\nu$ | $H_j=1/r_j-1/2$ | PDF pp. 3–4, Eqs. (13)–(17) |
| Dubois–Lallemand 2009 | $s_k$ | 第 $k$ 个非守恒矩的碰撞率 | $r_k$ | $r_k=s_k$ | PDF/印刷 p. 4, Eq. (12) |
| Dubois–Lallemand 2009 | $\sigma_k$ | Hénon 参数 | $H_k$ | $H_k=\sigma_k=1/s_k-1/2$ | PDF/印刷 p. 5, Eq. (16) |
| Ginzburg–d'Humières 2003 | $\lambda_\nu,\lambda_2$ | 两个碰撞率 | $r_{\mathrm{even}},r_{\mathrm{odd}}$ | $H_j=1/\lambda_j-1/2$ | PDF p. 3, Eqs. (6), (10)–(12) |
| Ginzburg–d'Humières 2003 | $\Lambda^2$ | 归一化偶奇 Hénon 乘积 | 保留 $\Lambda^2$ | $\Lambda^2=(4/3)H_{\mathrm{even}}H_{\mathrm{odd}}$ | PDF p. 6, Eq. (41) |

逐条读表会发现：

1. Wang/Contrino 的温度 $\sigma_\kappa$ 是碰撞率；
2. Dubois 的裸 $\sigma_k$ 却已经是 $1/s_k-1/2$；
3. Ginzburg 的小写 $\lambda$ 是碰撞率；
4. Ginzburg 的大写 $\Lambda^2$ 不是碰撞率，也不等于未经缩放的 $H_{\mathrm{even}}H_{\mathrm{odd}}$。

另外，第二篇应用论文的本地文件名前缀虽然是 `[Luo2014]`，正文第一作者实际是 **Dario Contrino**，所以本教程写作 Contrino et al. 2014。

因此，看到两个公式都写 $\sigma$，不能直接令它们相等。正确流程要先判断“原符号的类型”，再分三条路：

1. 若原符号本身是碰撞率，例如 Wang/Contrino 的 $\sigma_\kappa$ 或 Ginzburg 的 $\lambda_\nu$，才使用

\[
\text{论文原符号}
\longrightarrow
\text{碰撞率 }r
\longrightarrow
\tau=1/r
\longrightarrow
H=\tau-1/2
\]

这一整条率—时间—Hénon 翻译链。

2. Dubois 的 $\sigma_k$ 已经是 Hénon 参数，所以直接写

\[
\sigma_k\longrightarrow H_k,
\qquad
r_k=\frac1{H_k+1/2}.
\]

3. Ginzburg 的 $\Lambda^2$ 是两个 Hénon 参数的归一化乘积，保留为

\[
\Lambda^2=\frac43H_{\mathrm{even}}H_{\mathrm{odd}},
\]

不能把它塞进单一的 $r\to\tau\to H$ 链。

---

## 7. 本章最容易混淆的五件事

### 7.1 投影的 $1/2$ 与 Hénon 的 $-1/2$ 不是同一个来源

- $h_i^\pm$ 中的 $1/2$ 来自两个数的和差反解；
- $H=1/r-1/2$ 中的 $-1/2$ 来自二阶 Taylor 迁移项。

它们数值上都出现二分之一，但物理和代数来源不同。

### 7.2 $\tau$ 不是最终进入黏性或扩散率的全部尺度

$\tau=1/r$ 只做率—时间的倒数换算。二阶等效输运系数通常含 $H=\tau-1/2$，还要乘速度集和平衡矩给出的模型系数。

### 7.3 $1<r<2$ 不等于“不稳定”

此时偏差越过平衡，但局部绝对值仍缩小。完整稳定性必须结合具体模型、边界和参数判断。

### 7.4 TRT 的“两个率”不是随意挑两个矩

它们与中心反演下的偶、奇投影对应。若没有先定义反向方向和投影，直接写两个率就失去 TRT 的结构含义。

### 7.5 后面的 magic 常数不是一个普适数字

下一章会看到：在 Ginzburg 的受限流场边界分析中，$\Lambda^2=1/4$ 与 $H_{\mathrm{even}}H_{\mathrm{odd}}=3/16$ 是**同一个边界约束**的归一化与未归一化写法；温度 D2Q5 体相四阶各向同性的 $H_\kappa H_\nu=1/6$ 才是另一个对象。前两者也只在稳态、格线对齐平壁、halfway bounce-back、线性/Stokes、均匀力及 $I_f=-1/2$ 等条件下由 Ginzburg p. 7, Eqs. (42)–(43) 得到，不能当成任意边界的普适常数。

---

## 8. 读完本章应能独立完成的检查

1. 给定 $r=1.25$，算出 $\tau=0.8$、$H=0.3$，并算出偏差倍率 $1-r=-0.25$。
2. 给定 $h_i=0.7$、$h_{\bar i}=0.3$，算出 $h_i^+=0.5$、$h_i^-=0.2$，再重构两个原数。
3. 交换 $i$ 与 $\bar i$，验证偶部不变、奇部变号。
4. 从式 (2.47) 令 $r_+=r_-=r$，不跳步地恢复 BGK。
5. 解释 Hénon 参数中的 $1/r$ 和 $-1/2$ 分别来自哪一步。
6. 看到 Wang 的 $\sigma_\kappa$ 时，先把它识别为 $r_\kappa$；看到 Dubois 的 $\sigma_1$ 时，把它识别为 $H_1$。

如果这六项中有任何一项还不能独立完成，先回到对应小节，不要直接进入 magic 参数。

---

## 9. 本章来源边界

- 松弛矩更新、$0<s_k<2$、BGK 作为统一率特例：Dubois & Lallemand, *Towards higher order lattice Boltzmann schemes*, PDF/印刷 p. 4, Eqs. (12)–(14)。
- Hénon 参数 $\sigma_k=1/s_k-1/2$ 及二阶输运系数与之成正比：同文 PDF/印刷 p. 5, Eq. (16) 及后文说明。
- Wang/Contrino 的 D2Q9 与 D2Q5 碰撞率符号：两篇应用论文 PDF pp. 3–4 的相应碰撞矩阵与参数式。
- Ginzburg 的 $\lambda_\nu,\lambda_2,\Lambda^2$ 定义：*Multireflection boundary conditions for lattice Boltzmann models*, PDF pp. 3, 6，尤其 Eq. (41)；$\Lambda^2=1/4$、对应 $H_{\mathrm{even}}H_{\mathrm{odd}}=3/16$ 及 halfway 有效宽度结论来自 PDF p. 7, Eqs. (42)–(43)，并依赖稳态平行平板 Poiseuille、格线对齐 halfway BB、线性/Stokes、均匀力和 $I_f=-1/2$ 等限定条件。
- 本章 D1Q2 Taylor 展开是对式 (2.17) 的教学性逐阶演算，用于解释 $-1/2$ 的来源；它不替代四篇论文中 D2Q5/D2Q9 的具体高阶等效方程。

下一章将在这里建立的 $r\leftrightarrow\tau\leftrightarrow H$ 与偶奇投影之上，严格限定场景，推导边界滑移和 boundary magic。
