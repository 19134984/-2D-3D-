# D2Q9–TRT–\(a\)–\(\chi_\kappa\) 温度场代码修改说明

> 用途：让 Codex 在现有自然对流代码中实现/修改一套 **D2Q9–TRT–\(a\)–\(\chi_\kappa\)** 温度场算法，用于与 Luo D2Q5–TRT 和原始 LBM-CDE D2Q9–BGK 温度模型对比。
>
> 本文件只描述温度场。**不要修改现有流场算法、Ra/Pr 定义、浮力模型和流场边界条件。**

---

## 1. 目标方程

温度场求解

\[
\frac{\partial T}{\partial t}
+\nabla\cdot(T\mathbf u)
=
\nabla\cdot(\kappa\nabla T)+Q.
\]

当前自然对流问题取

\[
\boxed{Q=0}.
\]

目标热扩散率 \(\kappa_{\rm target}\) 由现有物理参数设置得到，例如

\[
\kappa_{\rm target}=\frac{\nu}{Pr}.
\]

不同温度算法比较时必须保持完全相同的物理

\[
\boxed{\kappa_{\rm target}}.
\]

---

## 2. D2Q9 与固定格子参数

温度场使用标准 D2Q9：

\[
\mathbf c_0=(0,0),
\]

\[
\mathbf c_{1-4}
=
(1,0),(0,1),(-1,0),(0,-1),
\]

\[
\mathbf c_{5-8}
=
(1,1),(-1,1),(-1,-1),(1,-1).
\]

反向映射：

\[
\bar 0=0,\;
\bar1=3,\;
\bar2=4,\;
\bar3=1,\;
\bar4=2,\;
\bar5=7,\;
\bar6=8,\;
\bar7=5,\;
\bar8=6.
\]

标准格子声速保持

\[
\boxed{c_s^2=\frac13}.
\]

标准 D2Q9 Hermite 权重保持

\[
\boxed{
\omega_0=\frac49,\qquad
\omega_{1-4}=\frac19,\qquad
\omega_{5-8}=\frac1{36}.
}
\]

注意：

\[
\boxed{
a\text{ 不修改 }c_s^2\text{ 和 }\omega_i.
}
\]

\(a\) 只修改温度 equilibrium 的基础扩散权重。

---

## 3. 参数 \(a\) 与基础温度权重

沿用 Luo 的参数化

\[
c_e=\frac{4+a}{10}.
\]

第一版固定 standard D2Q9 weight family：

\[
\boxed{
W_0(a)=\frac{2-a}{6},
}
\]

\[
\boxed{
W_{1-4}(a)=\frac{4+a}{30},
}
\]

\[
\boxed{
W_{5-8}(a)=\frac{4+a}{120}.
}
\]

这些权重满足

\[
\sum_iW_i(a)=1,
\]

以及

\[
\boxed{
\sum_iW_i(a)c_{i\alpha}c_{i\beta}
=
\frac{4+a}{10}\delta_{\alpha\beta}.
}
\]

标准 D2Q9 对应

\[
\frac{4+a}{10}=\frac13,
\]

所以

\[
\boxed{
a_{\rm std}=-\frac23.
}
\]

此时

\[
W_0=\frac49,\qquad
W_{1-4}=\frac19,\qquad
W_{5-8}=\frac1{36}.
\]

为了至少保证基础静态权重非负：

\[
\boxed{-4\le a\le2.}
\]

这只是权重非负条件，不代表新模型完整稳定区已经证明。

---

## 4. 压力修正基向量 \(\lambda_i\)

保留原始 LBM-CDE 的压力修正：

\[
\boxed{
\lambda_0=-\frac59,
}
\]

\[
\boxed{
\lambda_{1-4}=\frac19,
}
\]

\[
\boxed{
\lambda_{5-8}=\frac1{36}.
}
\]

满足

\[
\sum_i\lambda_i=0,
\]

\[
\sum_i\mathbf c_i\lambda_i=0,
\]

\[
\sum_i c_{i\alpha}c_{i\beta}\lambda_i
=
c_s^2\delta_{\alpha\beta}.
\]

---

## 5. 新温度平衡态

将 equilibrium 分为四部分：

\[
g_i^{eq}
=
g_i^{(0)}
+
g_i^{(u)}
+
g_i^{(uu)}
+
g_i^{(p)}.
\]

### 5.1 基础扩散部分

\[
\boxed{
g_i^{(0)}=W_i(a)T.
}
\]

### 5.2 一阶对流部分

仍使用标准 D2Q9 Hermite 项：

\[
\boxed{
g_i^{(u)}
=
\omega_iT
\frac{\mathbf c_i\cdot\mathbf u}{c_s^2}.
}
\]

### 5.3 二阶速度部分

\[
\boxed{
g_i^{(uu)}
=
\omega_iT
\left[
\frac{(\mathbf c_i\cdot\mathbf u)^2}{2c_s^4}
-
\frac{u^2}{2c_s^2}
\right].
}
\]

该项必须保留，用来恢复完整

\[
Tu_\alpha u_\beta,
\]

特别是

\[
Tu_xu_y.
\]

### 5.4 压力修正部分

\[
\boxed{
g_i^{(p)}
=
\lambda_i
\frac{Tp}{\rho_0c_s^2}.
}
\]

### 5.5 完整 equilibrium

\[
\boxed{
\begin{aligned}
g_i^{eq}
={}&
W_i(a)T
+
\omega_iT
\frac{\mathbf c_i\cdot\mathbf u}{c_s^2}
\\
&+
\omega_iT
\left[
\frac{(\mathbf c_i\cdot\mathbf u)^2}{2c_s^4}
-
\frac{u^2}{2c_s^2}
\right]
+
\lambda_i
\frac{Tp}{\rho_0c_s^2}.
\end{aligned}
}
\]

---

## 6. equilibrium 的矩约束

代码必须满足

\[
\boxed{
\sum_i g_i^{eq}=T,
}
\]

\[
\boxed{
\sum_i c_{i\alpha}g_i^{eq}=Tu_\alpha,
}
\]

以及

\[
\boxed{
\sum_i
c_{i\alpha}c_{i\beta}g_i^{eq}
=
\frac{4+a}{10}T\delta_{\alpha\beta}
+
Tu_\alpha u_\beta
+
\frac{Tp}{\rho_0}\delta_{\alpha\beta}.
}
\]

建议在 debug/unit-test 模式中直接验证上述三个矩。

---

## 7. TRT 奇偶分解

\[
g_i^+
=
\frac{g_i+g_{\bar i}}2,
\qquad
g_i^-=
\frac{g_i-g_{\bar i}}2.
\]

\[
g_i^{eq,+}
=
\frac{g_i^{eq}+g_{\bar i}^{eq}}2,
\qquad
g_i^{eq,-}
=
\frac{g_i^{eq}-g_{\bar i}^{eq}}2.
\]

定义

\[
n_i^+=g_i^+-g_i^{eq,+},
\qquad
n_i^-=g_i^--g_i^{eq,-}.
\]

---

## 8. TRT 松弛参数

第一版固定

\[
\boxed{
\sigma_-=
\frac1{s_-}-\frac12
=
\frac1{\sqrt{12}},
}
\]

\[
\boxed{
\sigma_+=
\frac1{s_+}-\frac12
=
\frac1{\sqrt3}.
}
\]

因此

\[
\boxed{
s_-=3-\sqrt3
\approx1.267949192431123,
}
\]

\[
\boxed{
s_+=4\sqrt3-6
\approx0.928203230275509.
}
\]

第一版中 \(s_-\)、\(s_+\) 都固定，不由物理 \(\kappa\) 决定。

---

## 9. 基础热扩散率

当

\[
\chi_\kappa=0
\]

时：

\[
\boxed{
\kappa_0(a)
=
\frac{4+a}{10}
\left(
\frac1{s_-}-\frac12
\right)\Delta t.
}
\]

采用固定 quartic 参数后：

\[
\boxed{
\kappa_0(a)
=
\frac{\sqrt3(4+a)}{60}\Delta t.
}
\]

格子单位 \(\Delta t=1\)：

\[
\boxed{
\kappa_0(a)
=
\frac{\sqrt3(4+a)}{60}.
}
\]

---

## 10. \(\chi_\kappa\) 与目标物理扩散率

定义

\[
\boxed{
\kappa_{\rm target}
=
(1-\chi_\kappa)\kappa_0(a).
}
\]

所以

\[
\boxed{
\kappa_{\rm target}
=
(1-\chi_\kappa)
\frac{4+a}{10}
\left(
\frac1{s_-}-\frac12
\right)\Delta t.
}
\]

固定 \(s_-\) 后：

\[
\boxed{
\kappa_{\rm target}
=
(1-\chi_\kappa)
\frac{\sqrt3(4+a)}{60}\Delta t.
}
\]

第一版的参数策略：

\[
\boxed{
a=\text{主动数值设计参数},
}
\]

\[
\boxed{
\chi_\kappa=\text{由 }a\text{ 和 }\kappa_{\rm target}\text{ 自动反算}.
}
\]

即

\[
\boxed{
\chi_\kappa(a)
=
1-
\frac{60\kappa_{\rm target}}
{\sqrt3(4+a)\Delta t}
}
\]

或

\[
\boxed{
\chi_\kappa(a)
=
1-
\frac{20\sqrt3\,\kappa_{\rm target}}
{(4+a)\Delta t}.
}
\]

格子单位：

\[
\boxed{
\chi_\kappa(a)
=
1-
\frac{20\sqrt3\,\kappa_{\rm target}}
{4+a}.
}
\]

注意：

\[
\boxed{
a,\chi_\kappa,s_-\text{ 不是三个独立自由参数。}
}
\]

固定 \(\kappa_{\rm target}\) 和 \(s_-\) 后，只主动选择 \(a\)，\(\chi_\kappa\) 必须同步更新。

---

## 11. 温度源项的一阶矩

按照修改后的 equilibrium 二阶矩重新进行 LBM-CDE inverse design：

\[
\boxed{
R_\alpha
=
\left[
\frac p{\rho_0}
+
\chi_\kappa\frac{4+a}{10}
\right]
\partial_\alpha T
+
\frac{TF_\alpha}{\rho_0}
+
u_\alpha Q.
}
\]

当前 \(Q=0\)：

\[
\boxed{
R_\alpha
=
\left[
\frac p{\rho_0}
+
\chi_\kappa\frac{4+a}{10}
\right]
\partial_\alpha T
+
\frac{TF_\alpha}{\rho_0}.
}
\]

---

## 12. 离散速度空间源项

当前 \(Q=0\)：

\[
\boxed{
R_i
=
\frac{\omega_i}{c_s^2}
\mathbf c_i\cdot
\left[
\left(
\frac p{\rho_0}
+
\chi_\kappa\frac{4+a}{10}
\right)\nabla T
+
\frac{T\mathbf F}{\rho_0}
\right].
}
\]

如果代码存储的是加速度

\[
\mathbf b=\frac{\mathbf F}{\rho_0},
\]

使用

\[
\boxed{
R_i
=
\frac{\omega_i}{c_s^2}
\mathbf c_i\cdot
\left[
\left(
\frac p{\rho_0}
+
\chi_\kappa\frac{4+a}{10}
\right)\nabla T
+
T\mathbf b
\right].
}
\]

必须确认现有流场变量到底是 \(\mathbf F\) 还是 \(\mathbf F/\rho_0\)，避免重复除以 \(\rho_0\)。

---

## 13. TRT 离散源项

当前 \(Q=0\) 时

\[
R_{\bar i}=-R_i,
\]

所以 source 为纯奇源项。

使用

\[
\boxed{
\Psi_i
=
\left(
1-\frac{s_-}{2}
\right)R_i.
}
\]

---

## 14. TRT 碰撞

\[
\boxed{
\begin{aligned}
g_i^*
={}&
g_i
-
s_+
\left(
g_i^+-g_i^{eq,+}
\right)
\\
&-
s_-
\left(
g_i^--g_i^{eq,-}
\right)
+
\Delta t
\left(
1-\frac{s_-}{2}
\right)R_i.
\end{aligned}
}
\]

随后 streaming：

\[
\boxed{
g_i(\mathbf x+\mathbf c_i\Delta t,t+\Delta t)
=
g_i^*(\mathbf x,t).
}
\]

如果现有 LBM-CDE BGK 代码中的 `g` 已经是梯形源项变换后的演化变量，则继续沿用其数据结构和 convention，不要混用原始 distribution 与 transformed distribution。

---

## 15. 宏观温度

当前 \(Q=0\)：

\[
\boxed{
T=\sum_i g_i.
}
\]

---

## 16. 局部温度梯度

梯度仍然来自一阶非平衡矩。

因为新的 equilibrium 二阶矩包含

\[
\frac{4+a}{10}T\delta_{\alpha\beta},
\]

所以 generalized local-gradient formula 中自然出现 \(a\)。

当前 \(Q=0\)：

\[
\boxed{
\nabla T
=
-
\frac{
2\displaystyle\sum_i
\mathbf c_i
(g_i-g_i^{eq})
+
\Delta t\,T\mathbf F/\rho_0
}{
2\kappa_{\rm target}
+
\Delta t
\left[
p/\rho_0+(4+a)/10
\right]
}.
}
\]

如果代码存储的是加速度 \(\mathbf b=\mathbf F/\rho_0\)：

\[
\boxed{
\nabla T
=
-
\frac{
2\displaystyle\sum_i
\mathbf c_i
(g_i-g_i^{eq})
+
\Delta t\,T\mathbf b
}{
2\kappa_{\rm target}
+
\Delta t
\left[
p/\rho_0+(4+a)/10
\right]
}.
}
\]

### 符号检查

分子是整体负号：

\[
-\left[
2\sum_i\mathbf c_i(g_i-g_i^{eq})
+\Delta t\,T\mathbf F/\rho_0
\right].
\]

不要写成

\[
-2\sum_i\mathbf c_i(g_i-g_i^{eq})
+\Delta t\,T\mathbf F/\rho_0.
\]

---

## 17. 等温 Dirichlet 边界

使用 halfway anti-bounce-back。

通用形式：

\[
\boxed{
g_{\bar i}(\mathbf x_f,t+\Delta t)
=
-g_i^*(\mathbf x_f,t)
+
2g_i^{eq}(T_w,\mathbf u_w=0,p_w;a).
}
\]

静止壁面：

\[
g_i^{eq}(T_w,0,p_w;a)
=
W_i(a)T_w
+
\lambda_i\frac{T_wp_w}{\rho_0c_s^2}.
\]

所以

\[
\boxed{
g_{\bar i}
=
-g_i^*
+
2W_i(a)T_w
+
2\lambda_i\frac{T_wp_w}{\rho_0c_s^2}.
}
\]

### 17.1 轴向方向

\[
W_{1-4}(a)=\frac{4+a}{30},
\qquad
\lambda_{1-4}=\frac19.
\]

因此

\[
\boxed{
g_{\bar i}
=
-g_i^*
+
\left[
\frac{4+a}{15}
+
\frac{2p_w}{3\rho_0}
\right]T_w.
}
\]

若当前代码关闭/忽略 scalar pressure correction：

\[
\boxed{
g_{\bar i}
=
-g_i^*
+
\frac{4+a}{15}T_w.
}
\]

### 17.2 对角方向

\[
W_{5-8}(a)=\frac{4+a}{120},
\qquad
\lambda_{5-8}=\frac1{36}.
\]

因此

\[
\boxed{
g_{\bar i}
=
-g_i^*
+
\left[
\frac{4+a}{60}
+
\frac{p_w}{6\rho_0}
\right]T_w.
}
\]

若关闭 pressure correction：

\[
\boxed{
g_{\bar i}
=
-g_i^*
+
\frac{4+a}{60}T_w.
}
\]

### 17.3 对现有 LBM-CDE 边界代码的要求

如果原代码已经存在

\[
2\omega_iT_w
+
2\lambda_i\frac{T_wp_w}{\rho_0c_s^2}
\]

这种 wall-equilibrium 处理，则保留现有 `p_w` 处理方式不变，只把

\[
2\omega_iT_w
\]

替换为

\[
\boxed{
2W_i(a)T_w.
}
\]

不要为这次修改额外发明新的 wall-pressure 插值规则。

---

## 18. 绝热 Neumann 边界

继续使用 halfway bounce-back：

\[
\boxed{
g_{\bar i}(\mathbf x_f,t+\Delta t)
=
g_i^*(\mathbf x_f,t).
}
\]

---

## 19. 初始化

初始化温度后，从现有流场取得

\[
\mathbf u(\mathbf x,0),\qquad
p(\mathbf x,0),
\]

使用新的 equilibrium：

\[
\boxed{
g_i(\mathbf x,0)
=
g_i^{eq}
(T,\mathbf u,p;a).
}
\]

保持现有 LBM-CDE transformed-distribution convention。

---

## 20. 推荐配置结构

```text
temperature_model = D2Q9_TRT_A_CHIKAPPA

a_temperature = -2.0/3.0

s_minus_temperature = 3 - sqrt(3)
s_plus_temperature  = 4*sqrt(3) - 6

kappa_target = existing_physical_kappa
chi_kappa_mode = AUTO_FROM_A
```

每次参数初始化：

```text
ce = (4 + a_temperature) / 10

W0 = (2 - a_temperature) / 6
Wc = (4 + a_temperature) / 30
Wd = (4 + a_temperature) / 120

sigma_minus = 1/s_minus_temperature - 0.5

kappa0 = ce * sigma_minus * dt

chi_kappa = 1 - kappa_target / kappa0
```

固定 quartic \(s_-\) 时等价于：

```text
chi_kappa =
    1 - 20*sqrt(3)*kappa_target /
        ((4 + a_temperature)*dt)
```

---

## 21. 每步计算建议顺序

1. 从当前温度 populations 计算 \(T\)；
2. 从流场读取 \(\mathbf u,p,\mathbf F\)；
3. 用 \(T,\mathbf u,p,a\) 计算新的 \(g_i^{eq}\)；
4. 用 local nonequilibrium relation 计算 \(\nabla T\)；
5. 根据 \(\nabla T,p,T,\mathbf F,\chi_\kappa,a\) 计算 \(R_i\)；
6. 做 TRT 奇偶分解；
7. 做 TRT collision，并加入
   \[
   (1-s_-/2)R_i;
   \]
8. streaming；
9. 温度边界：
   - 等温：anti-bounce-back；
   - 绝热：bounce-back；
10. 进入下一时间步。

如果现有代码 collision/streaming/boundary 的执行顺序不同，尽量保持现有框架，只确保数学依赖关系正确。

---

## 22. 必须增加的 consistency checks

### 22.1 基础权重

检查

\[
\left|\sum_iW_i(a)-1\right|.
\]

### 22.2 equilibrium 零阶矩

\[
\left|
\sum_i g_i^{eq}-T
\right|.
\]

### 22.3 equilibrium 一阶矩

\[
\left|
\sum_i c_{ix}g_i^{eq}-Tu_x
\right|,
\]

\[
\left|
\sum_i c_{iy}g_i^{eq}-Tu_y
\right|.
\]

### 22.4 equilibrium 二阶矩

检查

\[
\Pi_{xx}^{eq}
=
\frac{4+a}{10}T+Tu_x^2+\frac{Tp}{\rho_0},
\]

\[
\Pi_{yy}^{eq}
=
\frac{4+a}{10}T+Tu_y^2+\frac{Tp}{\rho_0},
\]

\[
\boxed{
\Pi_{xy}^{eq}=Tu_xu_y.
}
\]

### 22.5 扩散率重构

运行前输出：

```text
a
ce
s_minus
s_plus
sigma_minus
kappa0
chi_kappa
kappa_reconstructed
```

其中

\[
\boxed{
\kappa_{\rm reconstructed}
=
(1-\chi_\kappa)
\frac{4+a}{10}
\sigma_-\Delta t
}
\]

必须与

\[
\kappa_{\rm target}
\]

一致。

---

## 23. 与已有两套温度算法的关系

### Luo D2Q5–TRT

固定 quartic relaxation 后：

\[
\kappa
=
\frac{\sqrt3(4+a)}{60}.
\]

所以目标物理 \(\kappa\) 直接锁定 \(a\)。

### 原始 LBM-CDE D2Q9–BGK

固定 \(\tau_g\) 后：

\[
\kappa
=
(1-\chi_\kappa)c_s^2
\left(
\tau_g-\frac12
\right)\Delta t.
\]

所以目标物理 \(\kappa\) 锁定 \(\chi_\kappa\)。

### 新 D2Q9–TRT–\(a\)–\(\chi_\kappa\)

\[
\boxed{
\kappa
=
(1-\chi_\kappa)
\frac{\sqrt3(4+a)}{60}\Delta t.
}
\]

因此在固定物理 \(\kappa\) 和固定 TRT relaxation 后：

\[
\boxed{
a\text{ 可以作为 equilibrium/stability 的设计参数，}
}
\]

而

\[
\boxed{
\chi_\kappa(a)\text{ 自动补偿，保证物理 }\kappa\text{ 不变。}
}
\]

---

## 24. 当前理论状态与本次代码修改范围

本次只实现上述第一版模型。

当前可依赖的设计逻辑：

- equilibrium 的零阶、一阶、二阶矩满足目标；
- 二阶 Chapman–Enskog / inverse-design 层面恢复目标 CDE；
- \(\chi_\kappa\) 补偿 \(a\) 对基础扩散率的影响；
- 局部梯度和等温 ABB 与新的 \(a\)-dependent equilibrium 同步修改。

当前**不要假定**

\[
\chi_\kappa\neq0
\]

以后仍严格保持 source-free D2Q9 TRT 原有的全部 quartic cancellation。

因此本次修改：

1. 不新增额外四阶修正项；
2. 不重新优化 \(s_-\)、\(s_+\)；
3. 不修改 flow solver；
4. 先完成数值测试；
5. 后续再单独进行含 \(\chi_\kappa\) source 的高阶等效方程与稳定性分析。

---

## 25. Codex 修改时最核心的五项

### A. 基础 equilibrium weights

\[
\boxed{
W_0=\frac{2-a}{6},\qquad
W_{1-4}=\frac{4+a}{30},\qquad
W_{5-8}=\frac{4+a}{120}.
}
\]

### B. 一阶、二阶速度项仍使用标准 D2Q9

\[
\boxed{
c_s^2=\frac13,\qquad
\omega_i=(4/9,1/9,1/36).
}
\]

### C. TRT relaxation 固定

\[
\boxed{
s_-=3-\sqrt3,\qquad
s_+=4\sqrt3-6.
}
\]

### D. 给定 \(a\) 后自动反算

\[
\boxed{
\chi_\kappa
=
1-
\frac{20\sqrt3\,\kappa_{\rm target}}
{(4+a)\Delta t}.
}
\]

### E. source、gradient、isothermal boundary 必须同步使用新的 \(a\)-dependent equilibrium

\[
\boxed{
R_i
=
\frac{\omega_i}{c_s^2}
\mathbf c_i\cdot
\left[
\left(
\frac p{\rho_0}
+
\chi_\kappa\frac{4+a}{10}
\right)\nabla T
+
\frac{T\mathbf F}{\rho_0}
\right],
}
\]

\[
\boxed{
\nabla T
=
-
\frac{
2\sum_i\mathbf c_i(g_i-g_i^{eq})
+
\Delta t\,T\mathbf F/\rho_0
}{
2\kappa_{\rm target}
+
\Delta t[p/\rho_0+(4+a)/10]
},
}
\]

\[
\boxed{
g_{\bar i}
=
-g_i^*
+
2W_i(a)T_w
+
2\lambda_i\frac{T_wp_w}{\rho_0c_s^2}.
}
\]
