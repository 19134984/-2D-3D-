# 三维热对流 MRT-LBM 算法实现说明（依据文献复现）

## 1. 目标

本文档用于指导按文献实现三维热对流的双分布函数 MRT-LBM 程序。  
文献采用：

- **流场：D3Q19-MRT**
- **温度场：D3Q7-MRT**
- **耦合方式：Boussinesq 浮力项**
- **外力处理：Guo forcing 的 MRT 形式**
- **流场边界：half-way bounce-back**
- **温度定温边界：half-way anti-bounce-back**
- **温度绝热边界：half-way bounce-back**

实现目标是恢复如下宏观方程：

\[
\nabla \cdot \mathbf{u} = 0
\]

\[
\frac{\partial \mathbf{u}}{\partial t}
+ \mathbf{u} \cdot \nabla \mathbf{u}
=
-\frac{1}{\rho_0}\nabla p
+ \nu \nabla^2 \mathbf{u}
+ g \beta_T (T-T_0)\hat{\mathbf{z}}
\]

\[
\frac{\partial T}{\partial t}
+ \mathbf{u}\cdot \nabla T
=
\kappa \nabla^2 T
\]

---

## 2. 程序结构建议

建议将代码分为如下模块：

### 2.1 参数模块
包含：

- 网格数：`nx, ny, nz`
- 时间步长：`dt = 1`
- 格距：`dx = 1`
- Rayleigh 数 `Ra`
- Prandtl 数 `Pr`
- Mach 数 `Ma`
- 黏性系数 `nu`
- 热扩散率 `kappa`
- 流场松弛参数
- 温度场松弛参数

### 2.2 速度集模块
分别定义：

- 流场 D3Q19：`exf, eyf, ezf, wf`
- 温度 D3Q7：`ext, eyt, ezt, wt`

### 2.3 矩阵模块
分别定义：

- 流场变换矩阵：`M(19,19)`，`Minv(19,19)`
- 温度变换矩阵：`N(7,7)`，`Ninv(7,7)`

### 2.4 宏观量恢复模块
恢复：

- `rho`
- `ux, uy, uz`
- `T`

### 2.5 平衡态模块
计算：

- `feq(i)`
- `meq(k)`
- `geq(i)`
- `neq(k)`

### 2.6 外力项模块
计算：

- `Fx, Fy, Fz`
- `Fbar`
- `Fp = M^{-1}(I-S/2)M Fbar`

### 2.7 碰撞模块
- 流场 MRT 碰撞
- 温度 MRT 碰撞

### 2.8 迁移模块
- `stream_f`
- `stream_g`

### 2.9 边界条件模块
- 流场无滑移
- 温度定温
- 温度绝热

### 2.10 主循环模块
每一步依次执行：

1. 恢复宏观量
2. 计算浮力
3. 流场碰撞
4. 流场迁移
5. 流场边界
6. 温度碰撞
7. 温度迁移
8. 温度边界
9. 输出统计量

---

## 3. 流场模型：D3Q19-MRT

## 3.1 D3Q19 速度集

取 \(c=\delta x/\delta t\)，文献中 \(\delta x=\delta t=1\)，因此 \(c=1\)。

D3Q19 速度集为：

\[
\mathbf{e}_0=(0,0,0)
\]

\[
\mathbf{e}_1=(1,0,0),\quad
\mathbf{e}_2=(-1,0,0),\quad
\mathbf{e}_3=(0,1,0),\quad
\mathbf{e}_4=(0,-1,0),\quad
\mathbf{e}_5=(0,0,1),\quad
\mathbf{e}_6=(0,0,-1)
\]

\[
\mathbf{e}_7=(1,1,0),\quad
\mathbf{e}_8=(-1,1,0),\quad
\mathbf{e}_9=(1,-1,0),\quad
\mathbf{e}_{10}=(-1,-1,0)
\]

\[
\mathbf{e}_{11}=(1,0,1),\quad
\mathbf{e}_{12}=(-1,0,1),\quad
\mathbf{e}_{13}=(1,0,-1),\quad
\mathbf{e}_{14}=(-1,0,-1)
\]

\[
\mathbf{e}_{15}=(0,1,1),\quad
\mathbf{e}_{16}=(0,-1,1),\quad
\mathbf{e}_{17}=(0,1,-1),\quad
\mathbf{e}_{18}=(0,-1,-1)
\]

程序中可写为：

```fortran
exf = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 1,-1, 1,-1, 0, 0, 0, 0]
eyf = [0, 0, 0, 1,-1, 0, 0, 1, 1,-1,-1, 0, 0, 0, 0, 1,-1, 1,-1]
ezf = [0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1,-1,-1, 1, 1,-1,-1]
```

---

## 3.2 D3Q19 权重

\[
\omega_0=\frac{1}{3},\qquad
\omega_{1\sim 6}=\frac{1}{18},\qquad
\omega_{7\sim 18}=\frac{1}{36}
\]

程序中可写为：

```fortran
wf(0) = 1.0d0 / 3.0d0
wf(1:6) = 1.0d0 / 18.0d0
wf(7:18) = 1.0d0 / 36.0d0
```

---

## 3.3 流场声速

\[
c_s = \frac{c}{\sqrt{3}},\qquad c_s^2=\frac{1}{3}
\]

程序中通常设置：

```fortran
cs2_f = 1.0d0 / 3.0d0
```

---

## 3.4 流场演化方程

流场分布函数 \(f_i\) 满足 MRT 格式：

\[
f_i(\mathbf{x}+\mathbf{e}_i\delta_t,t+\delta_t)-f_i(\mathbf{x},t)
=
-(M^{-1}S)_{ij}\left(m_j-m_j^{(eq)}\right)
+\delta_t F_i'
\]

其中：

\[
m = M f
\]

\[
m^{(eq)} = M f^{(eq)}
\]

\[
F' = M^{-1}\left(I-\frac{S}{2}\right)M\tilde{F}
\]

更适合编程的矩空间碰撞写法是：

\[
m^+ = m - S(m-m^{(eq)}) + \delta_t \left(I-\frac{S}{2}\right)M\tilde{F}
\]

\[
f^+ = M^{-1} m^+
\]

然后再进行迁移。

---

## 3.5 流场平衡态分布函数 \(f_i^{eq}\)

采用标准二阶平衡态：

\[
f_i^{(eq)}
=
\omega_i \rho
\left[
1+\frac{\mathbf{e}_i\cdot\mathbf{u}}{c_s^2}
+\frac{(\mathbf{e}_i\cdot\mathbf{u})^2}{2c_s^4}
-\frac{|\mathbf{u}|^2}{2c_s^2}
\right]
\]

其中：

\[
|\mathbf{u}|^2=u^2+v^2+w^2
\]

实际实现时，更推荐直接使用矩空间平衡态 \(m^{eq}\)。

---

## 3.6 流场变换矩阵 \(M\)

文献使用的是 D3Q19 MRT 的正交矩阵基。  
基底顺序为：

\[
\langle 1 \rangle,
\langle 19e^2-30 \rangle,
\left\langle \frac{21}{2}e^4-\frac{53}{2}e^2+12 \right\rangle,
\langle e_x \rangle,
\langle (5e^2-9)e_x \rangle,
\langle e_y \rangle,
\langle (5e^2-9)e_y \rangle,
\langle e_z \rangle,
\langle (5e^2-9)e_z \rangle,
\]

\[
\langle 3e_x^2-e^2 \rangle,
\langle (3e^2-5)(3e_x^2-e^2) \rangle,
\langle e_y^2-e_z^2 \rangle,
\langle (3e^2-5)(e_y^2-e_z^2) \rangle,
\langle e_x e_y \rangle,
\langle e_y e_z \rangle,
\langle e_x e_z \rangle,
\]

\[
\langle (e_y^2-e_z^2)e_x \rangle,
\langle (e_z^2-e_x^2)e_y \rangle,
\langle (e_x^2-e_y^2)e_z \rangle
\]

### 实现建议
- **不要自己重新推导矩阵元素**
- 直接将论文给出的 19×19 数值矩阵写入 `M`
- 然后数值求逆得到 `Minv`
- 或直接把 `Minv` 也写死

若后续要交给 Codex，实现时应明确要求：  
**流场矩阵 M 采用论文原矩阵，不改动基底顺序。**

---

## 3.7 流场平衡矩 \(m^{eq}\)

文献给出的平衡矩为：

\[
m^{(eq)}=
\rho
\begin{bmatrix}
1 \\
-11+19|\mathbf{u}|^2 \\
3-\frac{11}{2}|\mathbf{u}|^2 \\
u \\
-\frac{2}{3}u \\
v \\
-\frac{2}{3}v \\
w \\
-\frac{2}{3}w \\
2u^2-v^2-w^2 \\
-\frac{1}{2}(2u^2-v^2-w^2) \\
v^2-w^2 \\
-\frac{1}{2}(v^2-w^2) \\
uv \\
vw \\
uw \\
0 \\
0 \\
0
\end{bmatrix}
\]

推荐代码直接按此公式构造 `meq(0:18)`。

---

## 3.8 流场松弛矩阵 \(S\)

\[
S=\mathrm{diag}(s_\rho,s_e,s_\varepsilon,s_j,s_q,s_j,s_q,s_j,s_q,
s_\nu,s_\pi,s_\nu,s_\pi,s_\nu,s_\nu,s_\nu,s_m,s_m,s_m)
\]

文献中取：

\[
s_\rho=s_j=0
\]

\[
s_e=s_\varepsilon=s_\nu=s_\pi=\frac{1}{\tau_f}
\]

\[
s_q=s_m=\frac{8(2\tau_f-1)}{8\tau_f-1}
\]

### 说明
若严格复现论文，就按上述方式设定。  
注意这与某些实现中“守恒模态取 1”不同，本文献原始做法是守恒模态取 0。

---

## 3.9 黏性与 \(\tau_f\) 的关系

\[
\nu = c_s^2 \left(\tau_f-\frac{1}{2}\right)\delta_t
\]

由于 \(\delta_t=1\)、\(c_s^2=1/3\)，因此：

\[
\nu = \frac{1}{3}\left(\tau_f-\frac{1}{2}\right)
\]

反过来：

\[
\tau_f = \frac{1}{2}+3\nu
\]

程序中通常流程为：

1. 先由无量纲参数换算得到格子单位下的 `nu`
2. 再计算 `tau_f`
3. 再计算各流场松弛率

---

## 3.10 外力项处理

文献采用 Guo forcing 的 MRT 形式：

\[
F' = M^{-1}\left(I-\frac{S}{2}\right)M\tilde{F}
\]

其中浮力为：

\[
\mathbf{F}=\rho g \beta_T (T-T_0)\hat{\mathbf{z}}
\]

即：

\[
F_x=0,\qquad F_y=0,\qquad F_z=\rho g\beta_T(T-T_0)
\]

文献直接给出 \(M\tilde{F}\) 的形式：

\[
M\tilde{F}=
\begin{bmatrix}
0 \\
38(\mathbf{u}\cdot\mathbf{F}) \\
-11(\mathbf{u}\cdot\mathbf{F}) \\
F_x \\
-\frac{2}{3}F_x \\
F_y \\
-\frac{2}{3}F_y \\
F_z \\
-\frac{2}{3}F_z \\
4uF_x-2vF_y-2wF_z \\
-2uF_x+vF_y+wF_z \\
2vF_y-2wF_z \\
-vF_y+wF_z \\
uF_y+vF_x \\
vF_z+wF_y \\
uF_z+wF_x \\
0 \\
0 \\
0
\end{bmatrix}
\]

### 编码建议
每个格点依次做：

1. 先由当前 `rho, ux, uy, uz, T` 计算 `Fx, Fy, Fz`
2. 构造 `MF(0:18)`
3. 乘以 `(I - S/2)`
4. 得到碰撞修正项
5. 加入矩空间碰撞方程

---

## 3.11 流场宏观量恢复

密度定义为：

\[
\rho = \sum_{i=0}^{18} f_i
\]

速度定义为：

\[
\mathbf{u}
=
\frac{1}{\rho}
\left(
\sum_{i=0}^{18}\mathbf{e}_i f_i+\frac{1}{2}\mathbf{F}
\right)
\]

即：

\[
u_x = \frac{\sum_i e_{ix} f_i + \frac{1}{2}F_x}{\rho}
\]

\[
u_y = \frac{\sum_i e_{iy} f_i + \frac{1}{2}F_y}{\rho}
\]

\[
u_z = \frac{\sum_i e_{iz} f_i + \frac{1}{2}F_z}{\rho}
\]

### 注意
这里的 \(+\frac{1}{2}\mathbf{F}\) 不能漏掉。

---

## 3.12 流场边界条件

文献对无滑移壁面采用 **half-way bounce-back**。

即对于指向流体内的反向分布：

\[
\bar{f}_i(\mathbf{x}_f,t+\delta t)=f_i^+(\mathbf{x}_f,t)
\]

其中：

- \(f_i^+\) 为碰撞后分布
- \(\bar{f}_i\) 为对向速度方向上的分布

---

## 4. 温度模型：D3Q7-MRT

## 4.1 D3Q7 速度集

D3Q7 速度集为：

\[
\mathbf{e}_0=(0,0,0)
\]

\[
\mathbf{e}_1=(1,0,0),\quad
\mathbf{e}_2=(-1,0,0),\quad
\mathbf{e}_3=(0,1,0),\quad
\mathbf{e}_4=(0,-1,0),\quad
\mathbf{e}_5=(0,0,1),\quad
\mathbf{e}_6=(0,0,-1)
\]

程序中可写为：

```fortran
ext = [0, 1,-1, 0, 0, 0, 0]
eyt = [0, 0, 0, 1,-1, 0, 0]
ezt = [0, 0, 0, 0, 0, 1,-1]
```

---

## 4.2 温度变换矩阵 \(N\)

文献的基底为：

\[
\langle 1 \rangle,\quad
\langle e_x \rangle,\quad
\langle e_y \rangle,\quad
\langle e_z \rangle,\quad
\langle -6+7e^2 \rangle,\quad
\langle 3e_x^2-e^2 \rangle,\quad
\langle e_y^2-e_z^2 \rangle
\]

对应数值矩阵可写为：

\[
N=
\begin{pmatrix}
1&1&1&1&1&1&1\\
0&1&-1&0&0&0&0\\
0&0&0&1&-1&0&0\\
0&0&0&0&0&1&-1\\
-6&1&1&1&1&1&1\\
0&2&2&-1&-1&-1&-1\\
0&0&0&1&1&-1&-1
\end{pmatrix}
\]

### 注意
论文文本若有 OCR 误差，应以基底定义为准。  
特别是第 4 行必须对应 \(\langle e_z\rangle\)，即 `[0,0,0,0,0,1,-1]`。

---

## 4.3 D3Q7 权重

文献温度平衡态权重不是固定常数，而是与参数 \(a_T\) 有关：

\[
\omega_0=\frac{1-a_T}{7},\qquad
\omega_{1\sim 6}=\frac{6+a_T}{42}
\]

因此温度场的权重初始化顺序必须是：

1. 先确定 `kappa`
2. 再确定 `a_T`
3. 最后计算 `wt(0:6)`

---

## 4.4 温度演化方程

温度分布函数 \(g_i\) 满足：

\[
g_i(\mathbf{x}+\mathbf{e}_i\delta_t,t+\delta_t)-g_i(\mathbf{x},t)
=
-(N^{-1}Q)_{ij}(n_j-n_j^{(eq)})
\]

矩空间写法：

\[
n = N g
\]

\[
n^+ = n - Q(n-n^{(eq)})
\]

\[
g^+ = N^{-1}n^+
\]

然后做迁移。

---

## 4.5 温度平衡态分布函数 \(g_i^{eq}\)

文献采用：

\[
g_i^{(eq)}
=
\omega_i T
\left[
1+\frac{7}{6+a_T}\frac{\mathbf{e}_i\cdot \mathbf{u}}{c_s^2}
\right]
\]

这里实际只保留到速度一阶项。

---

## 4.6 温度平衡矩 \(n^{eq}\)

直接写为：

\[
n^{(eq)}=
\begin{bmatrix}
T\\
uT\\
vT\\
wT\\
a_T T\\
0\\
0
\end{bmatrix}
\]

编程时推荐直接构造 `neq(0:6)`。

---

## 4.7 温度松弛矩阵 \(Q\)

\[
Q=\mathrm{diag}(0,q_\kappa,q_\kappa,q_\kappa,q_e,q_\nu,q_\nu)
\]

即：

- 零阶守恒模态：0
- 三个热流模态：`q_kappa`
- 一个能量模态：`q_e`
- 两个高阶模态：`q_nu`

---

## 4.8 热扩散率公式

文献中热扩散率为：

\[
\kappa
=
\frac{6+a_T}{21}
\left(
\frac{1}{q_\kappa}-\frac{1}{2}
\right)
\]

这是温度模型最核心的参数关系。

---

## 4.9 四阶误差各向同性条件

文献采用 Dubois 等人的条件，使 D3Q7 温度模型的四阶截断误差满足各向同性。

条件为：

\[
\left(
\frac{1}{q_\kappa}-\frac{1}{2}
\right)
\left(
\frac{1}{q_e}-\frac{1}{2}
\right)=\frac{1}{6}
\]

\[
\frac{1}{q_\nu}-\frac{1}{2}
=
\frac{a_T+6}{1-a_T}
\left(
\frac{1}{q_\kappa}-\frac{1}{2}
\right)
-
\frac{4+3a_T}{12(1-a_T)}
\left(
\frac{1}{q_\kappa}-\frac{1}{2}
\right)^{-1}
\]

并可写成：

\[
q_\nu
=
\frac{
6(1-a_T)(2-q_\kappa)q_\kappa
}{
(11+3a_T)(q_\kappa-6)q_\kappa + 12(a_T+6)
}
\]

---

## 4.10 文献最终采用的温度参数

文献进一步选取：

\[
\frac{1}{q_\kappa}-\frac{1}{2}=\frac{\sqrt{3}}{6}
\]

于是：

\[
q_\kappa = 3-\sqrt{3}
\]

\[
q_e=q_\nu=4\sqrt{3}-6
\]

\[
a_T = 42\sqrt{3}\,\kappa - 6
\]

### 推荐初始化方式
严格复现文献时，按以下顺序设置：

1. 给定目标热扩散率 `kappa`
2. 设置  
   \[
   q_\kappa=3-\sqrt{3}
   \]
3. 设置  
   \[
   q_e=q_\nu=4\sqrt{3}-6
   \]
4. 设置  
   \[
   a_T = 42\sqrt{3}\,\kappa - 6
   \]
5. 再由 `a_T` 求温度权重

---

## 4.11 温度宏观量恢复

温度定义为：

\[
T=\sum_{i=0}^{6} g_i
\]

温度场没有外力半步修正项。

---

## 4.12 温度边界条件

### 4.12.1 定温边界

文献使用 **half-way anti-bounce-back**：

\[
\bar{g}_i(\mathbf{x}_f,t+\delta t)
=
- g_i^+(\mathbf{x}_f,t)
+
\frac{6+a_T}{21}T_w
\]

其中 \(T_w\) 为壁面给定温度。

### 4.12.2 绝热边界

文献使用 **half-way bounce-back**：

\[
\bar{g}_i(\mathbf{x}_f,t+\delta t)=g_i^+(\mathbf{x}_f,t)
\]

---

## 5. 无量纲参数到格子参数的使用顺序

文献计算时给定：

- \(Ra\)
- \(Pr\)
- \(Ma\)

并通常固定：

\[
Ma=0.1
\]

实现中建议采用以下初始化顺序：

1. 设定 `Ra, Pr, Ma`
2. 确定特征长度 `L0`
3. 确定格点数 `nx, ny, nz`
4. 通过无量纲换算得到格子单位下的 `nu, kappa`
5. 用 `nu` 计算 `tau_f`
6. 用 `tau_f` 计算流场松弛率
7. 用 `kappa` 计算 `a_T`
8. 设置温度松弛率与权重

> 注：本文档只整理文献所需算法框架。若后续需要，还可以单独再写一份  
> **Ra-Pr-Ma 到格子单位参数的完整换算说明**。

---

## 6. 单步推进流程

下面给出每个时间步最推荐的执行顺序。

### Step 1. 恢复宏观量

#### 流场
\[
\rho = \sum_i f_i
\]

\[
\mathbf{u}
=
\frac{1}{\rho}
\left(
\sum_i \mathbf{e}_i f_i+\frac{1}{2}\mathbf{F}
\right)
\]

#### 温度
\[
T = \sum_i g_i
\]

---

### Step 2. 计算浮力

\[
F_x=0,\qquad
F_y=0,\qquad
F_z=\rho g\beta_T(T-T_0)
\]

---

### Step 3. 计算流场平衡矩

按公式构造 `meq(0:18)`。

---

### Step 4. 计算流场外力矩项

构造 `MF(0:18)`，即 \(M\tilde{F}\)。

---

### Step 5. 流场碰撞

\[
m = M f
\]

\[
m^+ = m - S(m-m^{eq}) + \left(I-\frac{S}{2}\right)M\tilde{F}
\]

\[
f^+ = M^{-1}m^+
\]

---

### Step 6. 流场迁移

\[
f_i(\mathbf{x}+\mathbf{e}_i,t+\delta t)=f_i^+(\mathbf{x},t)
\]

---

### Step 7. 流场边界

无滑移壁面施加 half-way bounce-back。

---

### Step 8. 计算温度平衡矩

按公式构造：

\[
n^{eq}=[T,uT,vT,wT,a_TT,0,0]^T
\]

---

### Step 9. 温度碰撞

\[
n = N g
\]

\[
n^+ = n - Q(n-n^{eq})
\]

\[
g^+ = N^{-1}n^+
\]

---

### Step 10. 温度迁移

\[
g_i(\mathbf{x}+\mathbf{e}_i,t+\delta t)=g_i^+(\mathbf{x},t)
\]

---

### Step 11. 温度边界

- 定温：anti-bounce-back
- 绝热：bounce-back

---

## 7. 代码实现时必须注意的事项

### 7.1 关于流场守恒模态
如果目标是**严格按文献复现**，则必须取：

\[
s_\rho=s_j=0
\]

不要擅自改成 1。

### 7.2 关于速度恢复
速度恢复式中的半步力修正项必须保留：

\[
\mathbf{u}
=
\frac{1}{\rho}
\left(
\sum_i \mathbf{e}_i f_i+\frac{1}{2}\mathbf{F}
\right)
\]

### 7.3 关于 D3Q7 权重
温度场权重依赖 \(a_T\)，不是固定常数。

### 7.4 关于温度参数
如果目的是复现论文，最稳妥做法是直接采用：

\[
q_\kappa=3-\sqrt{3},\qquad
q_e=q_\nu=4\sqrt{3}-6
\]

然后由目标 \(\kappa\) 求 \(a_T\)。

### 7.5 关于矩阵 \(N\)
若论文 OCR 有局部错误，应以基底定义为准，尤其 \(\langle e_z\rangle\) 那一行。

### 7.6 关于实现方式
建议碰撞统一在**矩空间**完成，而不是在分布空间直接展开矩阵乘法。

---

## 8. 最短复现清单

如果只需要最短的“照文献敲代码”版本，可以按下列固定方案：

### 流场
- D3Q19-MRT
- \(c_s^2=1/3\)
- \(s_\rho=s_j=0\)
- \(s_e=s_\varepsilon=s_\nu=s_\pi=1/\tau_f\)
- \(s_q=s_m=8(2\tau_f-1)/(8\tau_f-1)\)
- \(\nu=(1/3)(\tau_f-1/2)\)

### 温度场
- D3Q7-MRT
- \(q_\kappa=3-\sqrt{3}\)
- \(q_e=q_\nu=4\sqrt{3}-6\)
- \(a_T=42\sqrt{3}\kappa-6\)
- \(\omega_0=(1-a_T)/7\)
- \(\omega_{1\sim 6}=(6+a_T)/42\)

### 耦合
- \(\mathbf{F}=\rho g\beta_T(T-T_0)\hat z\)

### 边界
- 流场：bounce-back
- 温度定温：anti-bounce-back
- 温度绝热：bounce-back

---

## 9. 推荐的 Codex 任务拆分方式

后续如果把这份文档交给 Codex，建议按下面顺序拆任务：

1. 建立公共参数模块 `commondata`
2. 写 D3Q19 和 D3Q7 速度集与权重初始化
3. 写 `M, Minv, N, Ninv` 初始化
4. 写流场宏观恢复子程序
5. 写温度宏观恢复子程序
6. 写流场 `meq` 与外力项子程序
7. 写温度 `neq` 子程序
8. 写流场 MRT 碰撞子程序
9. 写温度 MRT 碰撞子程序
10. 写流场迁移与边界
11. 写温度迁移与边界
12. 写主时间推进循环
13. 写输出与后处理模块

---

## 10. 结论

这篇文献的可执行算法本质上是：

- 用 **D3Q19-MRT** 求解流场
- 用 **D3Q7-MRT** 求解温度场
- 用 **Boussinesq 浮力 + Guo forcing** 实现热-流耦合
- 用特别选取的 D3Q7 参数保证温度模型四阶误差各向同性
- 通过标准的 bounce-back / anti-bounce-back 完成边界处理

因此，只要按本文档中的：

- 速度集
- 权重
- 平衡态
- 平衡矩
- 松弛率
- 外力项
- 宏观恢复
- 边界处理
- 时间推进顺序

逐项实现，就可以完成文献算法的代码复现。
