# 参数化 D2Q5 四阶参考验证器

本附录只验证 Dubois--Lallemand 的 D2Q5 纯扩散模型。它为后续 D2Q9 推导校验符号工具，但任何 D2Q5 四阶系数都不能直接迁移到 D2Q9。

## 原始模型与参数化平衡态

速度按静止、东、北、西、南排列。Dubois Appendix Eq. (79) 的矩阵为

$$
M=\begin{pmatrix}
1&1&1&1&1\\
0&\lambda&0&-\lambda&0\\
0&0&\lambda&0&-\lambda\\
-4&1&1&1&1\\
0&1&-1&1&-1
\end{pmatrix}.
$$

纯热模型只有 $m_0=\rho$ 守恒，平衡矩是

$$
m^{eq}=\rho(1,0,0,\alpha,0)^{\mathsf T}.
$$

直接计算 $f^{eq}=M^{-1}m^{eq}$，得到一般参数权重

$$
w_0=\frac{1-\alpha}{5},
\qquad
w_1=w_2=w_3=w_4=\frac{4+\alpha}{20}.
$$

因此

$$
\sum_i w_i=1,
\qquad
\sum_iw_i c_{i\alpha}=0,
\qquad
\sum_iw_i c_{i\alpha}c_{i\beta}
=c_e\lambda^2\delta_{\alpha\beta},
$$

其中

$$
\boxed{c_e=\frac{4+\alpha}{10}}.
$$

Task 1 的固定 D2Q5 权重只对应 $\alpha=-2/3$：此时
$w_0=1/3$、$w_{1\ldots4}=1/6$、$c_e=1/3$。验证器始终先保留一般 $\alpha$，不会用这个固定点代替参数化模型。

Dubois Eq. (39) 的碰撞矩阵是

$$
\Psi=\begin{pmatrix}
1&0&0&0&0\\
0&1-s_1&0&0&0\\
0&0&1-s_1&0&0\\
\alpha s_3&0&0&1-s_3&0\\
0&0&0&0&1-s_4
\end{pmatrix},
\qquad
\sigma_j=\frac1{s_j}-\frac12.
$$

## 论文 Eq. (40) 的符号与尺度

原文 Eq. (40) 为

$$
\begin{aligned}
\partial_t\rho
&-\frac{\lambda^2\Delta t}{10}\sigma_1(4+\alpha)
(\partial_x^2+\partial_y^2)\rho\\
&+\frac{\Delta t^3\lambda^4}{1200}\sigma_1(4+\alpha)
\left[
\kappa_{40}(\partial_x^4+\partial_y^4)
+\kappa_{22}\partial_x^2\partial_y^2
\right]\rho
=O(\Delta t^4).
\end{aligned}
$$

所以扩散率和两个有量纲四阶系数分别是

$$
D=c_e\lambda^2\Delta t\,\sigma_1,
$$

$$
C_{40}=\frac{\Delta t^3\lambda^4}{1200}
\sigma_1(4+\alpha)\kappa_{40},
\qquad
C_{22}=\frac{\Delta t^3\lambda^4}{1200}
\sigma_1(4+\alpha)\kappa_{22}.
$$

代码内部先取 $\lambda=\Delta t=1$ 完成精确递推，再只在出口恢复
$D\mapsto\lambda^2\Delta tD$、$C_4\mapsto\lambda^4\Delta t^3C_4$。这是代数尺度分离，不是修改模型或把论文闭式代入递推。

## Route A：放大矩阵的流体根级数

原文第 17 页采用平面波 $f\propto\exp(i k_xx+i k_yy)$，并明确给出正相位

$$
B=\operatorname{diag}
(1,e^{iq_x},e^{iq_y},e^{-iq_x},e^{-iq_y}),
\qquad
q_\alpha=\lambda\Delta t\,k_\alpha,
$$

$$
G=B M^{-1}\Psi M.
$$

这里需要区分原文内部的两套符号写法：Dubois Eq. (15) 的空间迁移符号与第 17 页印出的正相位 $B$ 的空间位移符号相反。本验证器固定采用第 17 页的正相位约定；若改用 Eq. (15) 的约定，等价于作 $\boldsymbol k\mapsto-\boldsymbol k$。当前问题是纯扩散，本文核对的偶数阶系数因此不变。Task 5 一旦加入源项或平流项，相关奇数阶与源项符号不能由这里直接沿用，必须在其自行固定的 Fourier 约定下重新推导。

验证器不调用浮点特征值求解器生成符号答案。令

$$
G(\epsilon)=\sum_{n=0}^4\epsilon^nG_n,
\quad
z_h(\epsilon)=1+\sum_{n=1}^4\epsilon^nz_n,
\quad
v(\epsilon)=v_0+\sum_{n=1}^4\epsilon^nv_n,
$$

其中 $v_0=f^{eq}/\rho$，并选规范
$\boldsymbol 1^{\mathsf T}v_n=0$。第 $n$ 阶由同一个增广系统得到：

$$
\begin{pmatrix}
G_0-I&-v_0\\
\boldsymbol1^{\mathsf T}&0
\end{pmatrix}
\begin{pmatrix}v_n\\z_n\end{pmatrix}
{}={}
\begin{pmatrix}
-\displaystyle\sum_{r=1}^nG_rv_{n-r}
+\displaystyle\sum_{r=1}^{n-1}z_rv_{n-r}\\
0
\end{pmatrix}.
$$

最后展开

$$
\Gamma=-\frac1{\Delta t}\log z_h.
$$

奇数总次数严格为零；二次和四次齐次项分别给出 $D,C_{40},C_{22}$。

## Route B：物理空间 Taylor/矩递推

第二条路线不构造 $G$，不计算特征多项式、特征向量、$z_h$ 或对数。它直接从与正相位约定一致的离散迁移恒等式

$$
f_i(\boldsymbol x,t+\Delta t)
=f_i^*(\boldsymbol x+\boldsymbol c_i\Delta t,t)
$$

作物理空间 Taylor 展开：

$$
\sum_{p=0}^4\frac{\Delta t^p}{p!}\partial_t^pf_i
{}={}
\sum_{p=0}^4\frac{\Delta t^p}{p!}
(\boldsymbol c_i\cdot\nabla)^pf_i^*+O(\Delta t^5).
$$

写成 differential-jet 级数

$$
m=m^{eq}+\sum_{n=1}^4\epsilon^nm^{[n]},
\qquad
\partial_t\rho
=\sum_{r=0}^3\epsilon^rL_r(\partial_x,\partial_y)\rho.
$$

已知低阶 jet 后，第 $n$ 阶 Taylor 残差记为 $R_n$。未知的非守恒矩修正和守恒等效方程系数通过

$$
\boxed{
(I-\Psi)m^{[n]}+m^{eq}L_{n-1}=-R_n,
\qquad m_0^{[n]}=0
}
$$

逐阶消去。结果是

$$
L_0=0,
\qquad
L_1=c_e\sigma_1(\partial_x^2+\partial_y^2),
\qquad
L_2=0,
$$

而 $L_3$ 给出 Eq. (40) 的两个四阶项。Route B 只与 Route A 共享 $M$、$\Psi$、速度和通用截断多项式工具；测试还把 Route A 两个入口替换成抛异常函数，Route B 仍独立生成二阶系数。

## 两路线对 Eqs. (41)--(42) 的精确复现

两条递推的差经 SymPy 化简严格为零，并分别得到

$$
\boxed{
\kappa_{40}=8-3\alpha
+12(\alpha+4)\sigma_1^2
-12(1-\alpha)\sigma_1\sigma_3
-60\sigma_1\sigma_4
}
$$

以及

$$
\boxed{
\kappa_{22}=-6(\alpha+4)
+24(\alpha+4)\sigma_1^2
-24(1-\alpha)\sigma_1\sigma_3
+120\sigma_1\sigma_4.
}
$$

这两个式子是生成器输出与论文闭式的比较结果；Route A 和 Route B 的递推本身没有调用它们。

## 各向同性不等于完整消除

四阶算子与 $\nabla^4$ 各向同性一致只要求
$\kappa_{22}=2\kappa_{40}$。精确残差是

$$
\boxed{
\kappa_{22}-2\kappa_{40}
=40(6\sigma_1\sigma_4-1).
}
$$

所以各向同性关系为

$$
\sigma_1\sigma_4=\frac16.
$$

它只使两个四阶方向系数满足 $1:2$，一般并不使
$\kappa_{40}$ 或 $\kappa_{22}$ 为零。完整四阶消除必须同时解
$\kappa_{40}=\kappa_{22}=0$，得到 Dubois Eq. (55)：

$$
\boxed{
\sigma_3=\sigma_1\frac{\alpha+4}{1-\alpha}
-\frac{2+3\alpha}{12\sigma_1(1-\alpha)},
\qquad
\sigma_4=\frac1{6\sigma_1}.
}
$$

若再施加中间 TRT 约束 $\sigma_3=\sigma_4$，则

$$
\boxed{
\sigma_1=\frac1{\sqrt{12}},
\qquad
\sigma_3=\sigma_4=\frac1{\sqrt3}.
}
$$

代回后两个系数对任意允许的 $\alpha$ 都为零。

## Wang/Luo 符号映射

Wang Eq. (12) 与 Luo Eq. (14) 在 $\lambda=1$ 时使用相同的 D2Q5 矩阵；两文的自由参数 $a$ 对应本附录的 $\alpha$。两文的

$$
Q=\operatorname{diag}(\cdots,\sigma_\kappa,\sigma_\kappa,
\sigma_e,\sigma_\nu)
$$

中的 $\sigma_\kappa,\sigma_e,\sigma_\nu$ 是碰撞率，不是 Dubois 的 Hénon 平移。映射必须写成

$$
s_1=\sigma_\kappa,
\qquad s_3=\sigma_e,
\qquad s_4=\sigma_\nu,
$$

$$
\sigma_1=\frac1{\sigma_\kappa}-\frac12,
\quad
\sigma_3=\frac1{\sigma_e}-\frac12,
\quad
\sigma_4=\frac1{\sigma_\nu}-\frac12.
$$

Wang/Luo 的各向同性式因此正好对应
$\sigma_1\sigma_4=1/6$。该关系单独使用时不是完整四阶消除条件，更不是 D2Q9 条件。

## 精确与高精度抽查

一般符号比较之外，测试还使用三组允许范围内的有理参数：

| $\alpha$ | $\sigma_1$ | $\sigma_3$ | $\sigma_4$ | 轴向 $q^4$ 误差 | 对角 $q^4$ 误差 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| $-1$ | $1/5$ | $2/7$ | $3/11$ | $2.10\times10^{-16}$ | $1.41\times10^{-16}$ |
| $0$ | $1/3$ | $1/4$ | $1/6$ | $1.10\times10^{-15}$ | $3.79\times10^{-16}$ |
| $1/2$ | $2/5$ | $3/8$ | $4/9$ | $4.84\times10^{-16}$ | $9.91\times10^{-16}$ |

抽查使用 80 位精度、$q=10^{-6}$。轴向比较

$$
\frac{\Gamma(q,0)-Dq^2}{q^4}\to C_{40},
$$

对角比较

$$
\frac{\Gamma(q,q)-2Dq^2}{q^4}
\to2C_{40}+C_{22}.
$$

数值特征值只用于独立抽查，没有参与符号系数的生成。

## 结论边界

- 本附录验证的是参数化 D2Q5 纯扩散模型。
- D2Q5 不是通过删除 D2Q9 对角速度得到的连续参数退化。
- Wang/Luo 各向同性式不等于完整四阶消除。
- Eq. (55) 和 TRT 特殊点均为 D2Q5-only 结果，不能移植到 Task 5 的 D2Q9 目标模型。
