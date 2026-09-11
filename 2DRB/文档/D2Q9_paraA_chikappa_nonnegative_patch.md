# D2Q9–TRT 温度模型：`paraA` 参数化修改说明（\(\chi_\kappa\ge0\) 版本）

> 目的：交给 Codex 修改现有 D2Q9–TRT–LBM-CDE 温度代码。  
> 本文件只处理 `paraA` / \(a\) 的定义、\(c_e\) 映射、平衡态权重、\(\chi_\kappa\) 反算和参数合法区间。其余算法保持原设计。

## 1. `paraA` 的定义

继续使用

```text
paraA
```

作为唯一主动调节的温度平衡态参数，并定义

\[
a\equiv \texttt{paraA}.
\]

在新的 D2Q9 模型中，直接规定 D2Q9 的第三个非守恒标量平衡矩的静态部分为

\[
\boxed{m_3^{eq,0}=aT}.
\]

D2Q9 该矩满足

\[
m_3=3(\Pi_{xx}+\Pi_{yy})-4T.
\]

基础各向同性二阶矩定义为

\[
\Pi_{\alpha\beta}^{eq,0}=c_eT\delta_{\alpha\beta},
\]

因此

\[
m_3^{eq,0}=(6c_e-4)T.
\]

与 \(m_3^{eq,0}=aT\) 对比：

\[
\boxed{a=6c_e-4},
\qquad
\boxed{c_e=\frac{a+4}{6}}.
\]

代码必须使用

```text
ce = (paraA + 4.0) / 6.0
```

不要再使用 Luo D2Q5 的

\[
c_e=(4+a)/10.
\]

后者只属于 D2Q5 的矩基。

---

## 2. 标准 D2Q9 基准点

标准 D2Q9 有

\[
c_e=c_s^2=\frac13.
\]

因此

\[
\frac{a+4}{6}=\frac13
\]

给出

\[
\boxed{a=-2}.
\]

所以

```text
paraA = -2.0
```

必须严格恢复标准 D2Q9。

---

## 3. 第一版基础平衡态权重

第一版继续固定轴向与对角方向的基础权重比例

\[
\boxed{W_c:W_d=4:1}.
\]

由

\[
W_0+4W_c+4W_d=1,
\]

\[
2W_c+4W_d=c_e
\]

得到

\[
W_0=1-\frac53c_e,\qquad
W_c=\frac{c_e}{3},\qquad
W_d=\frac{c_e}{12}.
\]

代入 \(c_e=(a+4)/6\)：

\[
\boxed{W_0=-\frac{5a+2}{18}},
\]

\[
\boxed{W_c=\frac{a+4}{18}},
\]

\[
\boxed{W_d=\frac{a+4}{72}}.
\]

代码：

```text
W0 = -(5.0 * paraA + 2.0) / 18.0
Wc =  (paraA + 4.0) / 18.0
Wd =  (paraA + 4.0) / 72.0
```

必须检查

\[
W_0+4W_c+4W_d=1,
\]

\[
2W_c+4W_d=c_e.
\]

---

## 4. \(m_4\) 不再作为独立参数

第一版固定 \(W_c:W_d=4:1\) 后，

\[
a_4^{(0)}=4-9c_e.
\]

代入 \(c_e=(a+4)/6\)：

\[
\boxed{a_4^{(0)}=-\frac{3a+4}{2}}.
\]

所以不要再把 \(a_4\) 作为独立输入。

完整平衡矩中：

\[
\boxed{
m_3^{eq}
=
\left[
a+3(u^2+v^2)+\frac{6p}{\rho_0}
\right]T
}
\]

以及

\[
\boxed{
m_4^{eq}
=
\left[
-\frac{3a+4}{2}
-3(u^2+v^2)
-\frac{9p}{\rho_0}
\right]T.
}
\]

其余保持

\[
m_0^{eq}=T,
\]

\[
m_1^{eq}=uT,\qquad m_2^{eq}=vT,
\]

\[
m_5^{eq}=-uT,\qquad m_6^{eq}=-vT,
\]

\[
m_7^{eq}=(u^2-v^2)T,
\qquad
m_8^{eq}=uvT.
\]

---

## 5. 速度空间平衡态中的 \(a\)

完整平衡态仍写成

\[
\boxed{
\begin{aligned}
g_i^{eq}
={}&W_i(a)T\\
&+\omega_iT
\left[
\frac{\mathbf c_i\cdot\mathbf u}{c_s^2}
+\frac{(\mathbf c_i\cdot\mathbf u)^2}{2c_s^4}
-\frac{u^2+v^2}{2c_s^2}
\right]\\
&+\lambda_i\frac{Tp}{\rho_0c_s^2}.
\end{aligned}
}
\]

其中固定

\[
c_s^2=\frac13,
\]

标准 D2Q9 权重

\[
\omega_0=\frac49,\qquad
\omega_{1-4}=\frac19,\qquad
\omega_{5-8}=\frac1{36},
\]

压力权重

\[
\lambda_0=-\frac59,\qquad
\lambda_{1-4}=\frac19,\qquad
\lambda_{5-8}=\frac1{36}.
\]

基础温度权重使用

\[
W_0=-\frac{5a+2}{18},
\qquad
W_{1-4}=\frac{a+4}{18},
\qquad
W_{5-8}=\frac{a+4}{72}.
\]

注意：

- \(c_e\) 不是新的格子声速；
- 速度项仍使用 \(c_s^2=1/3\)；
- 压力仍使用流场定义，例如 \(p=c_s^2\delta\rho\)；
- 不允许使用 \(p=c_e\delta\rho\)；
- \(\chi_\kappa\) 不进入 \(g_i^{eq}\)。

---

## 6. `paraA` 与基础热扩散率

无 \(\chi_\kappa\) 修正时

\[
\kappa_0
=
c_e
\left(
\frac1{s_-}-\frac12
\right)\Delta t.
\]

定义

\[
\sigma_-=\frac1{s_-}-\frac12.
\]

则

\[
\boxed{
\kappa_0=\frac{a+4}{6}\sigma_-\Delta t.
}
\]

第一版固定

\[
\boxed{\sigma_-=\frac1{\sqrt{12}}},
\qquad
\boxed{\sigma_+=\frac1{\sqrt3}},
\]

即

\[
s_-=3-\sqrt3,
\qquad
s_+=4\sqrt3-6.
\]

---

## 7. \(\chi_\kappa\) 的自动反算

最终目标扩散率定义为

\[
\boxed{
\kappa_{\rm target}
=
(1-\chi_\kappa)c_e\sigma_-\Delta t.
}
\]

代入 \(c_e=(a+4)/6\)：

\[
\boxed{
\kappa_{\rm target}
=
(1-\chi_\kappa)
\frac{a+4}{6}
\sigma_-\Delta t.
}
\]

因此

\[
\boxed{
\chi_\kappa
=
1-
\frac{6\kappa_{\rm target}}
{(a+4)\sigma_-\Delta t}.
}
\]

固定 \(\sigma_-=1/\sqrt{12}\) 后：

\[
\boxed{
\chi_\kappa
=
1-
\frac{12\sqrt3\,\kappa_{\rm target}}
{(a+4)\Delta t}.
}
\]

格子单位 \(\Delta t=1\)：

\[
\boxed{
\chi_\kappa
=
1-
\frac{12\sqrt3\,\kappa_{\rm target}}
{a+4}.
}
\]

代码中：

```text
ce = (paraA + 4.0) / 6.0
sigmaMinus = 1.0 / sMinus - 0.5

chiKappa = 1.0
         - kappaTarget / (ce * sigmaMinus * dt)
```

`chiKappa` 是从属量，不是独立输入参数。

---

## 8. 本版本强制要求 \(\chi_\kappa\ge0\)

要求

\[
\boxed{\chi_\kappa\ge0}.
\]

由

\[
\chi_\kappa
=
1-
\frac{6\kappa_{\rm target}}
{(a+4)\sigma_-\Delta t}
\]

得到

\[
\boxed{
a
\ge
\frac{6\kappa_{\rm target}}
{\sigma_-\Delta t}
-4.
}
\]

固定

\[
\sigma_-=\frac1{\sqrt{12}}
\]

后：

\[
\boxed{
a
\ge
\frac{12\sqrt3\,\kappa_{\rm target}}
{\Delta t}
-4.
}
\]

格子单位：

\[
\boxed{
a\ge12\sqrt3\,\kappa_{\rm target}-4.
}
\]

该下界对应

\[
\chi_\kappa=0.
\]

当 \(a\) 增大时，\(c_e\) 和默认扩散率增大，\(\chi_\kappa>0\) 用来把默认扩散率修正回 \(\kappa_{\rm target}\)。

---

## 9. 第一版 `paraA` 的完整允许范围

第一版还要求基础平衡态权重严格为正。

由

\[
W_c>0,\qquad W_d>0
\]

得到

\[
a>-4.
\]

由

\[
W_0>0
\]

得到

\[
a<-\frac25.
\]

因此仅由基础权重得到

\[
-4<a<-\frac25.
\]

再叠加 \(\chi_\kappa\ge0\)：

\[
\boxed{
\max\left(
-4,\,
\frac{6\kappa_{\rm target}}
{\sigma_-\Delta t}-4
\right)
\le a<-\frac25.
}
\]

因为 \(\kappa_{\rm target}>0\)，实际 \(\chi_\kappa\) 下界自然高于 \(-4\)，所以可以直接使用

\[
\boxed{
\frac{6\kappa_{\rm target}}
{\sigma_-\Delta t}-4
\le a<-\frac25.
}
\]

固定 \(\sigma_-=1/\sqrt{12}\)：

\[
\boxed{
\frac{12\sqrt3\,\kappa_{\rm target}}
{\Delta t}-4
\le a<-\frac25.
}
\]

格子单位：

\[
\boxed{
12\sqrt3\,\kappa_{\rm target}-4
\le a<-\frac25.
}
\]

如果下界大于等于 \(-2/5\)，说明在当前固定松弛率、基础权重正值和 \(\chi_\kappa\ge0\) 三个要求下没有可用参数区间，代码应直接报错。

---

## 10. 参数检查代码

建议在初始化阶段加入：

```text
sigmaMinus = 1.0 / sMinus - 0.5

aMin = 6.0 * kappaTarget / (sigmaMinus * dt) - 4.0
aMax = -0.4

if paraA < aMin:
    error("Invalid paraA: chiKappa would be negative")

if paraA >= aMax:
    error("Invalid paraA: W0 would be non-positive")

ce = (paraA + 4.0) / 6.0

chiKappa = 1.0
         - kappaTarget / (ce * sigmaMinus * dt)

W0 = -(5.0 * paraA + 2.0) / 18.0
Wc =  (paraA + 4.0) / 18.0
Wd =  (paraA + 4.0) / 72.0
```

检查：

```text
assert(ce > 0)
assert(chiKappa >= 0)
assert(chiKappa < 1)

assert(W0 > 0)
assert(Wc > 0)
assert(Wd > 0)

assert(abs(W0 + 4*Wc + 4*Wd - 1) < tol)
assert(abs(2*Wc + 4*Wd - ce) < tol)
```

下界处由于浮点舍入可能得到极小负数，例如 `-1e-15`；只允许对机器舍入量使用容差，不要把明显为负的 `chiKappa` 强行截断成 0。

---

## 11. 源项中需要同步更新

当前 \(Q=0\) 时：

\[
\boxed{
\mathbf A
=
\left(
\frac p{\rho_0}
+
\chi_\kappa c_e
\right)\nabla T
+
\frac{T\mathbf F}{\rho_0}.
}
\]

现在必须使用

\[
c_e=\frac{a+4}{6},
\]

所以

\[
\boxed{
\mathbf A
=
\left[
\frac p{\rho_0}
+
\chi_\kappa\frac{a+4}{6}
\right]\nabla T
+
\frac{T\mathbf F}{\rho_0}.
}
\]

速度空间源项保持

\[
R_i
=
\frac{\omega_i}{c_s^2}
\mathbf c_i\cdot\mathbf A.
\]

不要把 \(\chi_\kappa c_e\) 改成 \(\chi_\kappa c_s^2\)。

---

## 12. 局部温度梯度同步更新

若代码使用局部非平衡分布计算温度梯度，则

\[
\boxed{
\nabla T
=
-
\frac{
2\sum_i\mathbf c_i(\widetilde g_i-g_i^{eq})
+
\Delta t\,T\mathbf F/\rho_0
}{
2\kappa_{\rm target}
+
\Delta t
\left[
p/\rho_0+(a+4)/6
\right]
}.
}
\]

注意整个分子前有总负号。

---

## 13. 等温边界同步更新

静止等温壁的反弹回格式：

\[
g_{\bar i}
=
-g_i^*
+
2W_i(a)T_w
+
2\lambda_i\frac{T_wp_w}{\rho_0c_s^2}.
\]

基础温度部分：

轴向方向

\[
\boxed{
2W_cT_w=\frac{a+4}{9}T_w,
}
\]

对角方向

\[
\boxed{
2W_dT_w=\frac{a+4}{36}T_w.
}
\]

不要再使用 D2Q5 的 \((4+a)/10\) 作为 D2Q9 边界修正系数。

绝热边界保持

\[
g_{\bar i}=g_i^*.
\]

---

## 14. 必须增加的回归测试

### 测试 1：标准 D2Q9 点

令

```text
paraA = -2.0
```

必须得到

\[
c_e=\frac13,
\]

\[
W_0=\frac49,\qquad
W_c=\frac19,\qquad
W_d=\frac1{36}.
\]

并且静态时

\[
m_3^{eq,0}=-2T,
\qquad
m_4^{eq,0}=T.
\]

### 测试 2：扩散率一致性

对任意合法 `paraA`，重新计算

\[
\kappa_{\rm check}
=
(1-\chi_\kappa)c_e\sigma_-\Delta t.
\]

必须满足

\[
\boxed{
\kappa_{\rm check}=\kappa_{\rm target}
}
\]

至机器精度。

### 测试 3：\(\chi_\kappa\) 下界

令

\[
a=a_{\min}
=
\frac{6\kappa_{\rm target}}
{\sigma_-\Delta t}-4.
\]

应得到

\[
\chi_\kappa=0.
\]

若 \(a<a_{\min}\)，程序必须拒绝运行。

---

## 15. Codex 最终修改清单

1. 保留 `paraA` 作为唯一主动温度平衡态参数。
2. 改为
   ```text
   ce = (paraA + 4) / 6
   ```
3. 删除 D2Q9 中所有
   ```text
   ce = (paraA + 4) / 10
   ```
   的旧写法。
4. 更新 \(W_0,W_c,W_d\)。
5. 更新 \(m_3^{eq}\) 和 \(m_4^{eq}\)。
6. 不允许独立输入 \(a_4\)。
7. 根据 `paraA` 自动计算 `chiKappa`。
8. 强制 `chiKappa >= 0`。
9. 使用
   \[
   6\kappa_{\rm target}/(\sigma_-\Delta t)-4
   \le a< -2/5
   \]
   作为第一版合法区间。
10. 固定特殊松弛率时使用
    \[
    12\sqrt3\,\kappa_{\rm target}/\Delta t-4
    \le a< -2/5.
    \]
11. 同步修改源项、局部温度梯度和等温边界中的 \(c_e\)。
12. 保持 \(c_s^2=1/3\) 与流场压力定义不变。
13. `chiKappa` 不进入 \(g_i^{eq}\)。
14. 增加 `paraA=-2` 标准 D2Q9 回归测试。
15. 增加 \(\kappa_{\rm check}=\kappa_{\rm target}\) 和 \(\chi_\kappa\ge0\) 测试。

---

## 16. 最终参数链

本版本代码中的参数关系必须严格保持为

\[
\boxed{
\texttt{paraA}=a
\longrightarrow
c_e=\frac{a+4}{6}
\longrightarrow
\{W_i,m_k^{eq}\}
}
\]

并由

\[
\boxed{
\kappa_{\rm target}
=
(1-\chi_\kappa)c_e\sigma_-\Delta t
}
\]

自动反算

\[
\boxed{\chi_\kappa\ge0}.
\]

因此：

- `paraA`：唯一主动调节参数；
- `ce`：由 `paraA` 得到；
- \(a_4\)：由 `paraA` 和第一版 \(4:1\) 权重封闭自动得到；
- `chiKappa`：由 `paraA`、目标 \(\kappa\) 和固定松弛率自动得到。
