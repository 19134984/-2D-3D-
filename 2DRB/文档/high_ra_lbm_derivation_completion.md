# 29. Codex 完成推导：适用范围与最终结论

> 本节是对前述“任务 1--10”的完成稿。推导以
> `N=[1,1,1,1,1; ...]` 所定义的 D2Q5 矩顺序为准，并以仓库
> `均匀网格/2DRBOpenmp.F90`、`均匀网格/3DRBOpenmp.F90`
> 中的实际 D2Q9/D3Q19 矩顺序完成流场源项投影。

先给出结论：

1. 对“主分支 A”的无热源、恒定速度线性问题，局部温度梯度代回碰撞后，温度通量矩的更新**严格等价**于使用
   \[
   q_{\kappa,\mathrm{eff}}
   =
   \frac{2q_\kappa^{(0)}}
   {2(1-\chi_\kappa)+\chi_\kappa q_\kappa^{(0)}}.
   \]
   因而所有纯扩散四阶公式只需把原 \(\sigma_1\) 换成
   \(\sigma_{\kappa,\mathrm{eff}}\)，不是另一个独立的四阶模型。
2. D2Q5 四阶各向同性仍要求
   \[
   \sigma_{\kappa,\mathrm{eff}}\sigma_4=\frac16.
   \]
   四阶完全消除还要求原来的 \(\sigma_3\) 关系，只是其中
   \(\sigma_1\) 换成 \(\sigma_{\kappa,\mathrm{eff}}\)。
3. 当 \(\kappa\to0\) 且 \(c_T^2\) 固定时，
   \(q_4\sim6\sigma_{\kappa,\mathrm{eff}}\to0\)。因此
   “小扩散率 + 严格四阶各向同性 + 所有 ghost 松弛率远离 0”
   三者不能同时满足。这个结论与 \(q_\kappa^{(0)}\) 的选择无关。
4. 若放弃严格四阶条件，只保留二阶一致性并独立优化 \(q_3,q_4\)，
   D2Q5 仍存在 von Neumann 稳定参数区，因此目前**不能判定必须升级
   D2Q9**。
5. 对流场应力修正，在标准 halfway bounce-back、源项先并入
   post-collision 分布的约定下，精确壁面关系为
   \[
   \boxed{
   \sigma_{\nu,\mathrm{eff}}\sigma_q=\frac{3}{16}
   },
   \qquad
   \sigma_{\nu,\mathrm{eff}}
   =(1-\chi_s)\sigma_\nu^{(0)}.
   \]
   所以标准 halfway bounce-back 不能避免
   \(s_q\to0\)。若要避免，只能改变边界算子，而不是把
   \(3/16\) 强行施加在基础参数上。

以下推导使用两套源项记号：

- \(\boldsymbol R_c\)：连续时间源项，量纲为“每单位时间”；
- \(\widehat{\boldsymbol R}=\delta t\,\boldsymbol R_c\)：一个时间步的原始源增量；
- 真正加入碰撞的离散源为
  \[
  \boldsymbol\Psi
  =
  \left(I-\frac Q2\right)\widehat{\boldsymbol R}.
  \]

这样可避免把连续源和离散源混写。

---

# 30. 任务 1：D2Q5 源项的完整矩空间表达

## 30.1 晶格、权重与矩顺序

速度顺序固定为

\[
\boldsymbol e_0=(0,0),\quad
\boldsymbol e_1=(1,0),\quad
\boldsymbol e_2=(0,1),\quad
\boldsymbol e_3=(-1,0),\quad
\boldsymbol e_4=(0,-1).
\]

由前述平衡分布可得

\[
w_0=\frac{1-a}{5},
\qquad
w_1=w_2=w_3=w_4=\frac{4+a}{20}
=\frac{c_T^2}{2},
\]

\[
\sum_iw_i=1,
\qquad
\sum_iw_i e_{i\alpha}=0,
\qquad
\sum_iw_i e_{i\alpha}e_{i\beta}
=c_T^2\delta_{\alpha\beta}.
\]

矩顺序固定为

\[
\boldsymbol n
=
[T,j_x,j_y,e,\nu_T]^{\mathrm T}
=N\boldsymbol g.
\]

## 30.2 完整弱可压缩源与主分支 A

为与 `LBM-CDE.pdf` 的式 (24)、(28)、(35) 对齐，定义

\[
\boldsymbol A
=
\frac{p\nabla T+T\boldsymbol F}{\rho_0}.
\]

允许体热源 \(Q_T\) 时，一阶 Hermite 截断的 D2Q5 原始分布源可写成

\[
\boxed{
\widehat R_i
=
\delta t\,w_i
\left[
Q_T\left(
1+\frac{\boldsymbol e_i\cdot\boldsymbol u}{c_T^2}
\right)
+\frac{\boldsymbol e_i\cdot\boldsymbol A}{c_T^2}
+\chi_\kappa\boldsymbol e_i\cdot\nabla T
\right]
}.
\]

其矩空间投影为

\[
\boxed{
\widehat{\boldsymbol R}^{\,m}
=
\delta t
\begin{bmatrix}
Q_T\\
A_x+u_xQ_T+\chi_\kappa c_T^2\partial_xT\\
A_y+u_yQ_T+\chi_\kappa c_T^2\partial_yT\\
aQ_T\\
0
\end{bmatrix}
}.
\]

因此离散源的五个分量为

\[
\boxed{
\boldsymbol\Psi^m
=
\begin{bmatrix}
\delta t\,Q_T\\[2mm]
\left(1-\dfrac{q_\kappa^{(0)}}2\right)\widehat R_x^m\\[2mm]
\left(1-\dfrac{q_\kappa^{(0)}}2\right)\widehat R_y^m\\[2mm]
\left(1-\dfrac{q_3}2\right)a\delta t\,Q_T\\[2mm]
0
\end{bmatrix}
}.
\]

主分支 A 用于本次四阶推导的闭合选择是

\[
Q_T=0,\qquad
\boldsymbol A=\boldsymbol0,
\]

即只保留扩散修正：

\[
\boxed{
\widehat{\boldsymbol R}^{\,m}_{A}
=
\delta t
\begin{bmatrix}
0\\
\chi_\kappa c_T^2\partial_xT\\
\chi_\kappa c_T^2\partial_yT\\
0\\
0
\end{bmatrix}
}.
\]

相应分布空间源为

\[
\boxed{
\widehat{\boldsymbol R}_{A}
=
\frac{\chi_\kappa c_T^2\delta t}{2}
\begin{bmatrix}
0\\
\partial_xT\\
\partial_yT\\
-\partial_xT\\
-\partial_yT
\end{bmatrix}
}.
\]

守恒检查直接给出

\[
\sum_i\widehat R_{A,i}=0,
\]

所以扩散修正不改变 \(T=\sum_i g_i\)。

## 30.3 局部梯度代回后的精确有效松弛率

定义通量非平衡矩

\[
j_\alpha^{neq}
=
\sum_i e_{i\alpha}(g_i-g_i^{eq})
=j_\alpha-u_\alpha T.
\]

在主分支 A 中，局部梯度公式为

\[
\boxed{
\partial_\alpha T
=
-\frac{q_{\kappa,\mathrm{eff}}}
{c_T^2\delta t}
j_\alpha^{neq}
},
\]

其中

\[
\boxed{
q_{\kappa,\mathrm{eff}}
=
\frac{2q_\kappa^{(0)}}
{2(1-\chi_\kappa)+
\chi_\kappa q_\kappa^{(0)}}
}.
\]

证明如下。通量矩的碰撞为

\[
j_\alpha^*
=
j_\alpha
-q_\kappa^{(0)}j_\alpha^{neq}
+\left(1-\frac{q_\kappa^{(0)}}2\right)
\chi_\kappa c_T^2\delta t\,\partial_\alpha T.
\]

把局部梯度代入：

\[
j_\alpha^*
=
j_\alpha
-\left[
q_\kappa^{(0)}
+\left(1-\frac{q_\kappa^{(0)}}2\right)
\chi_\kappa q_{\kappa,\mathrm{eff}}
\right]j_\alpha^{neq}.
\]

而

\[
q_\kappa^{(0)}
+\left(1-\frac{q_\kappa^{(0)}}2\right)
\chi_\kappa q_{\kappa,\mathrm{eff}}
=q_{\kappa,\mathrm{eff}},
\]

故

\[
\boxed{
j_\alpha^*
=
j_\alpha
-q_{\kappa,\mathrm{eff}}
(j_\alpha-u_\alpha T)
}.
\]

因此这不是仅在 Chapman--Enskog 主阶成立的近似，而是在当前线性局部闭合下的逐格点代数恒等式。

对应 Hénon 参数为

\[
\boxed{
h
\equiv
\sigma_{\kappa,\mathrm{eff}}
=
\frac1{q_{\kappa,\mathrm{eff}}}-\frac12
=(1-\chi_\kappa)
\left(
\frac1{q_\kappa^{(0)}}-\frac12
\right)
}.
\]

于是

\[
\kappa=c_T^2h\delta t.
\]

这也说明：对 MRT D2Q5 而言，主分支 A 的 \(\chi_\kappa\) 在常系数线性体节点中不会产生一个额外独立的稳定参数；它的净作用就是把通量矩改成 \(q_{\kappa,\mathrm{eff}}\)。它的主要价值是允许固定 \(a\)，避免平衡权重随 \(\kappa\) 退化，而不是把物理温度通量模态任意强阻尼。

---

# 31. 任务 2：D2Q5 纯扩散四阶等效方程

令

\[
\boldsymbol u=\boldsymbol0,
\qquad
h=\sigma_{\kappa,\mathrm{eff}},
\qquad
r=\sigma_3,
\qquad
t=\sigma_4,
\]

\[
\lambda=\frac{\delta x}{\delta t}.
\]

对碰撞--传播矩阵作四阶 Taylor 展开，得到

\[
\boxed{
\begin{aligned}
\partial_tT
&-
\frac{\lambda^2\delta t}{10}
h(a+4)\nabla^2T\\
&+
\frac{\lambda^4\delta t^3}{1200}
h(a+4)
\left[
K_{40}
\left(
\partial_x^4T+\partial_y^4T
\right)
+K_{22}\partial_x^2\partial_y^2T
\right]
=O(\delta t^4),
\end{aligned}
}
\]

其中

\[
\boxed{
K_{40}
=
8-3a
+12(a+4)h^2
-12(1-a)hr
-60ht
},
\]

\[
\boxed{
K_{22}
=
-6(a+4)
+24(a+4)h^2
-24(1-a)hr
+120ht
}.
\]

因为

\[
c_T^2=\frac{a+4}{10},
\]

二阶项正好给出

\[
\kappa=\lambda^2\delta t\,c_T^2h.
\]

若把方程写成

\[
\partial_tT
=
\kappa\nabla^2T
+\widetilde\kappa_{40}
(\partial_x^4+\partial_y^4)T
+\widetilde\kappa_{22}
\partial_x^2\partial_y^2T
+O(\delta t^4),
\]

则有

\[
\boxed{
\widetilde\kappa_{40}
=
-\frac{\lambda^4\delta t^3}{1200}
h(a+4)K_{40}
},
\]

\[
\boxed{
\widetilde\kappa_{22}
=
-\frac{\lambda^4\delta t^3}{1200}
h(a+4)K_{22}
}.
\]

## 31.1 四阶各向同性

四阶各向同性要求

\[
K_{22}=2K_{40}.
\]

两式相减：

\[
K_{22}-2K_{40}
=
-40+240ht.
\]

所以

\[
\boxed{
ht=\frac16
},
\]

即

\[
\boxed{
\sigma_{\kappa,\mathrm{eff}}\sigma_4=\frac16
}.
\]

## 31.2 四阶完全消除

再要求

\[
K_{40}=K_{22}=0,
\]

得到

\[
\boxed{
t=\frac1{6h}
},
\]

\[
\boxed{
r
=
h\frac{a+4}{1-a}
-
\frac{2+3a}{12h(1-a)}
}.
\]

当

\[
\chi_\kappa=0
\]

时

\[
h=\sigma_\kappa^{(0)},
\]

上述公式严格退化为 `Towards higher order lattice Boltzmann schemes.pdf`
式 (40)--(42)、(55)。符号脚本
`tools/derive_d2q5_cde.py` 已逐项验证该退化。

## 31.3 小扩散率下的无解结论

四阶各向同性要求

\[
q_4
=
\frac1{t+1/2}
=
\boxed{
\frac{6h}{1+3h}
}.
\]

所以

\[
h\to0
\quad\Longrightarrow\quad
q_4\sim6h\to0.
\]

进一步考察 \(r\)：

\[
r
=
\frac{
h(a+4)
-(2+3a)/(12h)
}{1-a}.
\]

对固定 \(a<1\)：

- 若 \(a>-2/3\)，则 \(r\to-\infty\)，松弛率 \(q_3\) 不可接受；
- 若 \(a<-2/3\)，则 \(r\to+\infty\)，所以 \(q_3\to0\)；
- 若 \(a=-2/3\)，则 \(r=2h\to0\)，所以 \(q_3\to2\)。

因此在 \(h\to0\) 极限下，即使忽略 \(q_4\)，\(q_3\) 也会趋向不可接受区域或失去正 Hénon 参数。

结论是

\[
\boxed{
\text{D2Q5 在小扩散率极限不能同时保持严格四阶条件和有界 ghost 阻尼。}
}
\]

---

# 32. 任务 3：D2Q5 恒定对流等效方程

## 32.1 精确 Fourier 放大矩阵

令

\[
\boldsymbol u=(u,v)=\text{constant},
\]

并采用上节已经证明的有效通量松弛率

\[
q_1=q_{\kappa,\mathrm{eff}}.
\]

源项闭合后的矩空间碰撞矩阵为

\[
C=
\begin{bmatrix}
1&0&0&0&0\\
q_1u&1-q_1&0&0&0\\
q_1v&0&1-q_1&0&0\\
q_3a&0&0&1-q_3&0\\
0&0&0&0&1-q_4
\end{bmatrix}.
\]

对 Fourier 模式

\[
\boldsymbol g(\boldsymbol x,t)
=
\widehat{\boldsymbol g}(t)
\exp(i\boldsymbol k\cdot\boldsymbol x),
\]

一步放大矩阵为

\[
\boxed{
A(\boldsymbol k)
=
D(\boldsymbol k)N^{-1}CN
},
\]

\[
D(\boldsymbol k)
=
\operatorname{diag}
\left(
1,
e^{-ik_x},
e^{-ik_y},
e^{ik_x},
e^{ik_y}
\right).
\]

该矩阵既用于等效方程，也用于完整 von Neumann 稳定性搜索。

## 32.2 二阶数值扩散张量

设 \(A\) 的水动力特征值为 \(\Lambda_h(\boldsymbol k)\)，并定义

\[
\log\Lambda_h
=z_1+z_2+z_3+z_4+O(|\boldsymbol k|^5),
\]

其中 \(z_n\) 是 \(\boldsymbol k\) 的 \(n\) 次齐次多项式。

一阶项为

\[
z_1=-i(uk_x+vk_y).
\]

二阶项为

\[
\boxed{
z_2
=
-h
\left[
c_T^2(k_x^2+k_y^2)
-(uk_x+vk_y)^2
\right]
}.
\]

因此长波等效扩散张量为

\[
\boxed{
\boldsymbol D_{\mathrm{num}}
=
h\delta t
\left(
c_T^2I-\boldsymbol u\boldsymbol u^{\mathrm T}
\right)
}.
\]

展开后

\[
D_{xx}=h\delta t(c_T^2-u^2),
\]

\[
D_{yy}=h\delta t(c_T^2-v^2),
\]

\[
D_{xy}=D_{yx}=-h\delta t\,uv.
\]

所以 D2Q5 的交叉扩散项并未被 \(\chi_\kappa\) 独立消除；它只通过
\(h=\sigma_{\kappa,\mathrm{eff}}\) 缩放。长波必要条件仍为

\[
\boxed{
u^2+v^2\le c_T^2
}.
\]

## 32.3 三阶色散和四阶误差

令

\[
h=\sigma_{\kappa,\mathrm{eff}},
\qquad
r=\sigma_3,
\qquad
t=\sigma_4.
\]

三阶项写成

\[
z_3
=
i\left(
d_{30}k_x^3+d_{21}k_x^2k_y
+d_{12}k_xk_y^2+d_{03}k_y^3
\right),
\]

其中

\[
\begin{aligned}
d_{30}
=\frac{u}{120}\big[
&-24ah^2-12ahr+3a
+240h^2u^2-96h^2\\
&+12hr+60ht-20u^2+2
\big],
\end{aligned}
\]

\[
\begin{aligned}
d_{21}
=\frac{v}{40}\big[
&-8ah^2-4ahr+a
+240h^2u^2-32h^2\\
&+4hr-20ht-20u^2+4
\big],
\end{aligned}
\]

\[
\begin{aligned}
d_{12}
=\frac{u}{40}\big[
&-8ah^2-4ahr+a
+240h^2v^2-32h^2\\
&+4hr-20ht-20v^2+4
\big],
\end{aligned}
\]

\[
\begin{aligned}
d_{03}
=\frac{v}{120}\big[
&-24ah^2-12ahr+3a
+240h^2v^2-96h^2\\
&+12hr+60ht-20v^2+2
\big].
\end{aligned}
\]

四阶项写成

\[
z_4
=
e_{40}k_x^4+e_{31}k_x^3k_y
+e_{22}k_x^2k_y^2
+e_{13}k_xk_y^3+e_{04}k_y^4.
\]

五个 \(e_{ij}\) 的完整展开见本文件附录 A。它们满足以下复核：

\[
u=v=0
\quad\Longrightarrow\quad
e_{40}=e_{04}
=-\frac{h(a+4)}{1200}K_{40},
\]

\[
e_{22}
=-\frac{h(a+4)}{1200}K_{22},
\qquad
e_{31}=e_{13}=0.
\]

因此常对流展开与纯扩散四阶式完全相容。

在当前 \(Q_T=\boldsymbol F=\boldsymbol0\)、常 \(\boldsymbol u\) 的闭合中，
所有显式 \(\chi_\kappa\) 交叉项都已经代数合并进 \(h\)；
不存在另一个独立的“源项--对流”系数。若恢复
\(\boldsymbol A\)、\(Q_T\) 或空间变化速度，则本节矩阵 \(C\) 不再完整，
必须按式 (35) 的全源公式重新展开，不能把本节结果外推到该情形。

---

# 33. 任务 4：D2Q5 稳定参数区域

## 33.1 参数空间降维

给定目标扩散率后，

\[
h
=
\frac{\kappa}{c_T^2\delta t},
\qquad
q_{\kappa,\mathrm{eff}}
=
\frac1{h+1/2}.
\]

基础参数满足

\[
\sigma_\kappa^{(0)}
=
\frac1{q_\kappa^{(0)}}-\frac12,
\]

\[
\boxed{
\chi_\kappa
=
1-\frac{h}{\sigma_\kappa^{(0)}}
}.
\]

所以 \((q_\kappa^{(0)},\chi_\kappa)\) 在主分支 A 的线性体节点中并非两个独立稳定参数。固定 \(h\) 后，它们只决定如何把同一个有效碰撞拆成“基础碰撞 + 源项”。

完整可行集合定义为

\[
\mathcal P
=
\left\{
(a,q_\kappa^{(0)},\chi_\kappa,q_3,q_4):
\begin{array}{l}
-4<a<1,\\
q_i\in[\varepsilon,2-\varepsilon],\\
0\le\chi_\kappa<1,\\
\kappa=c_T^2(1-\chi_\kappa)
\sigma_\kappa^{(0)}\delta t,\\
\max_{\boldsymbol k\in[-\pi,\pi]^2}
\rho[A(\boldsymbol k)]\le1
\end{array}
\right\}.
\]

还必须同时检查

\[
\max_\alpha|u_\alpha|\le c_T^2
\]

的人口非负条件，以及

\[
|\boldsymbol u|^2\le c_T^2
\]

的长波必要条件；二者不能互相替代。

## 33.2 严格四阶子集的解析下界

若要求 \(q_4\ge\varepsilon\)，由

\[
q_4=\frac{6h}{1+3h}
\]

得到

\[
\boxed{
h
\ge
\frac{\varepsilon}{3(2-\varepsilon)}
}.
\]

相应最小扩散率为

\[
\boxed{
\kappa_{\min}^{(4)}
=
c_T^2\delta t
\frac{\varepsilon}{3(2-\varepsilon)}
}.
\]

例如

\[
\varepsilon=0.05
\]

时

\[
h_{\min}^{(4)}
\approx8.547\times10^{-3}.
\]

只要目标 \(h\) 低于该值，严格四阶各向同性子集就必为空，与
\(q_\kappa^{(0)}\) 和 \(\chi_\kappa\) 如何拆分无关。

还应检查

\[
\sigma_{\min}
=
\frac1{2-\varepsilon}-\frac12,
\qquad
\sigma_{\max}
=
\frac1\varepsilon-\frac12,
\]

\[
\sigma_{\min}\le r(a,h)\le\sigma_{\max}.
\]

这给出 \(a\) 的第二个可行性约束。

## 33.3 代表性 von Neumann 扫描

为判断“放弃严格四阶后 D2Q5 是否仍有稳定区域”，进行了如下可复现实验：

- \(a=0\)，故 \(c_T^2=0.4\)；
- \(|\boldsymbol u|=0.1\)，同时检查轴向和 \(45^\circ\) 对流；
- \(q_3,q_4=0.05,0.10,\ldots,1.95\)，共 \(39^2=1521\) 组；
- \(k_x,k_y\in[0,\pi]\)，每方向 25 点；
- 判据为所有特征值满足
  \(\rho[A(\boldsymbol k)]\le1+2\times10^{-10}\)。

结果如下：

| \(h\) | \(q_{\kappa,\mathrm{eff}}\) | 稳定组数 / 1521 | 稳定点的 \(q_3\) 外包络 | 稳定点的 \(q_4\) 外包络 |
|---:|---:|---:|---:|---:|
| \(10^{-1}\) | 1.666667 | 1521 | 0.05--1.95 | 0.05--1.95 |
| \(10^{-2}\) | 1.960784 | 1166 | 0.05--1.95 | 0.05--1.95 |
| \(10^{-3}\) | 1.996008 | 102 | 0.95--1.95 | 1.65--1.95 |

“外包络”只表示稳定点的最小/最大坐标，不表示包络中的每一对参数都稳定。

同一 \(a=0\) 下，严格四阶公式给出：

| \(h\) | \(q_3^{(4)}\) | \(q_4^{(4)}\) | 结论 |
|---:|---:|---:|---|
| \(10^{-1}\) | -1.30435 | 0.46154 | \(q_3\) 非法 |
| \(10^{-2}\) | -0.06201 | 0.05825 | \(q_3\) 非法 |
| \(10^{-3}\) | -0.00602 | 0.00598 | \(q_3\) 非法且 \(q_4\) 近 0 |

这组扫描支持两个不同结论：

1. 严格四阶 D2Q5 在小 \(h\) 下不可行；
2. 二阶 D2Q5 仍可能通过独立选择 ghost 松弛率获得稳定区。

该扫描不是高 \(Ra\) 实际算例的运行验证。正式参数必须把目标
\(Ra,Pr,Ma,N,\delta x,\delta t\) 映射到 \(h\)，再按实际速度上界重扫。

---

# 34. 任务 5：流场源项一致的 halfway bounce-back 推导

## 34.1 应力源使剪切矩成为有效碰撞

以 D2Q9 的 \(xy\) 剪切矩为例，忽略体力交叉项时

\[
\Pi_{xy}^{neq}
=
\sum_i c_{ix}c_{iy}
(f_i-f_i^{eq}).
\]

`LBM-CDE.pdf` 式 (31) 给出

\[
S_{xy}
=
-\frac{\Pi_{xy}^{neq}}
{\rho_0(2\nu+c_s^2\delta t)}.
\]

应力修正源投影到该矩的原始分量为

\[
\widehat R_{xy}^{\nu}
=
2\rho_0c_s^2\chi_sS_{xy}\delta t.
\]

设基础剪切松弛率为 \(s_\nu^{(0)}\)，则

\[
\Pi_{xy}^{neq,*}
=
\left(1-s_\nu^{(0)}\right)\Pi_{xy}^{neq}
+\left(1-\frac{s_\nu^{(0)}}2\right)
\widehat R_{xy}^{\nu}.
\]

利用

\[
\nu
=(1-\chi_s)c_s^2
\sigma_\nu^{(0)}\delta t,
\]

可化简为

\[
\boxed{
\Pi_{xy}^{neq,*}
=
\left(1-s_{\nu,\mathrm{eff}}\right)
\Pi_{xy}^{neq}
},
\]

\[
\boxed{
s_{\nu,\mathrm{eff}}
=
\frac1{
(1-\chi_s)\sigma_\nu^{(0)}+1/2
}
}.
\]

即应力修正源没有让物理剪切矩保持基础放大因子；它把该矩的实际放大因子变为有效值。

## 34.2 标准 halfway bounce-back 的闭合

标准 halfway bounce-back 只看到完整 post-collision 分布，因此其偶模态参数必须使用实际的

\[
\sigma_{\nu,\mathrm{eff}}
=
\frac1{s_{\nu,\mathrm{eff}}}-\frac12.
\]

对当前 D2Q9/TRT 矩顺序，原精确半格点关系遂变为

\[
\boxed{
\mathcal B
\left(
\sigma_\nu^{(0)},\chi_s,\sigma_q
\right)
=
(1-\chi_s)\sigma_\nu^{(0)}\sigma_q
-\frac{3}{16}
=0
}.
\]

无源极限检查：

\[
\chi_s=0
\quad\Longrightarrow\quad
\sigma_\nu^{(0)}\sigma_q=\frac3{16}.
\]

由闭合关系

\[
\sigma_q
=
\frac{3}{16\sigma_{\nu,\mathrm{eff}}},
\]

\[
\boxed{
s_q
=
\frac{16\sigma_{\nu,\mathrm{eff}}}
{3+8\sigma_{\nu,\mathrm{eff}}}
}.
\]

因此

\[
\nu\to0
\quad\Longrightarrow\quad
s_q\sim\frac{16}{3}\sigma_{\nu,\mathrm{eff}}\to0.
\]

最终判定是：

\[
\boxed{
\text{标准 halfway bounce-back 不能借助 }\chi_s
\text{ 避免 }s_q\to0.
}
\]

若改用

\[
\sigma_\nu^{(0)}\sigma_q=\frac3{16},
\]

则得到的是“基础碰撞的旧壁面关系”，而不是含源实际传播算子的闭合；它一般会产生参数相关壁面偏移或滑移。

要同时保持适中 \(s_q\) 与精确壁面，必须引入边界专用自由度，例如 multireflection、插值反弹或源项一致的 link correction。仅由体节点控制方程不能唯一确定这个额外边界项，因此本文件不把某个未经 Poiseuille 解析解验证的修正写成已证明公式。

---

# 35. 任务 6：实际 D2Q9/D3Q19 流场源项投影

## 35.1 统一定义

定义

\[
A_{\alpha\beta}
=
\chi_sS_{\alpha\beta}
+\frac{\chi_b-\chi_s}{D}
S_{\gamma\gamma}\delta_{\alpha\beta}.
\]

分布空间的应力原始源为

\[
\widehat R_i^\nu
=
\delta t\,w_i\rho_0
\left(
\frac{c_{i\alpha}c_{i\beta}}{c_s^2}
-\delta_{\alpha\beta}
\right)A_{\alpha\beta}.
\]

总原始矩源为

\[
\widehat{\boldsymbol R}^{\,m}
=M
\left(
\widehat{\boldsymbol R}^{\,F}
+\widehat{\boldsymbol R}^{\,\nu}
\right),
\]

真正加入碰撞的是

\[
\boxed{
\boldsymbol\Phi^m
=
\left(I-\frac S2\right)
\widehat{\boldsymbol R}^{\,m}
}.
\]

以下给出的都是乘半步因子之前的原始矩源，且省略公共的 \(\delta t\)。

## 35.2 D2Q9

仓库 D2Q9 矩顺序为

\[
[\rho,e,\epsilon,j_x,q_x,j_y,q_y,p_{xx},p_{xy}]^{\mathrm T}.
\]

令

\[
u_F=u_xF_x+u_yF_y.
\]

Guo/Boussinesq 力源投影为

\[
\widehat{\boldsymbol R}^{\,F}_{D2Q9}
=
\begin{bmatrix}
0\\
6u_F\\
-6u_F\\
F_x\\
-F_x\\
F_y\\
-F_y\\
2(u_xF_x-u_yF_y)\\
u_xF_y+u_yF_x
\end{bmatrix}.
\]

应力修正投影为

\[
\boxed{
\widehat{\boldsymbol R}^{\,\nu}_{D2Q9}
=
\rho_0
\begin{bmatrix}
0\\
2(A_{xx}+A_{yy})\\
-2(A_{xx}+A_{yy})\\
0\\
0\\
0\\
0\\
\dfrac23(A_{xx}-A_{yy})\\
\dfrac23A_{xy}
\end{bmatrix}
}.
\]

所以 Fortran 中每一行应实现

\[
\Phi_j
=
\left(1-\frac{s_j}{2}\right)
\left(
\widehat R_j^F+\widehat R_j^\nu
\right).
\]

不能把 BGK 的单一公共 prefactor 直接复制到 MRT 的全部矩。

## 35.3 D3Q19

仓库 D3Q19 矩顺序为

\[
\begin{aligned}
[
&\rho,e,\epsilon,
j_x,q_x,j_y,q_y,j_z,q_z,\\
&p_{xx},\pi_{xx},p_{ww},\pi_{ww},
p_{xy},p_{yz},p_{xz},
m_x,m_y,m_z
]^{\mathrm T}.
\end{aligned}
\]

令

\[
u_F=u_xF_x+u_yF_y+u_zF_z.
\]

力源投影为

\[
\widehat{\boldsymbol R}^{\,F}_{D3Q19}
=
\begin{bmatrix}
0\\
38u_F\\
-11u_F\\
F_x\\
-\frac23F_x\\
F_y\\
-\frac23F_y\\
F_z\\
-\frac23F_z\\
4u_xF_x-2u_yF_y-2u_zF_z\\
-2u_xF_x+u_yF_y+u_zF_z\\
2u_yF_y-2u_zF_z\\
-u_yF_y+u_zF_z\\
u_xF_y+u_yF_x\\
u_yF_z+u_zF_y\\
u_xF_z+u_zF_x\\
0\\
0\\
0
\end{bmatrix}.
\]

应力修正投影为

\[
\boxed{
\widehat{\boldsymbol R}^{\,\nu}_{D3Q19}
=
\rho_0
\begin{bmatrix}
0\\
\dfrac{38}{3}(A_{xx}+A_{yy}+A_{zz})\\
-\dfrac{11}{3}(A_{xx}+A_{yy}+A_{zz})\\
0\\
0\\
0\\
0\\
0\\
0\\
\dfrac23(2A_{xx}-A_{yy}-A_{zz})\\
-\dfrac13(2A_{xx}-A_{yy}-A_{zz})\\
\dfrac23(A_{yy}-A_{zz})\\
-\dfrac13(A_{yy}-A_{zz})\\
\dfrac23A_{xy}\\
\dfrac23A_{yz}\\
\dfrac23A_{xz}\\
0\\
0\\
0
\end{bmatrix}
}.
\]

五个独立偏应力分量对应矩 9、11、13、14、15；矩 10、12 是当前非正交基下与两个对角应力组合伴随的高阶矩，因此实际投影不为零。忽略 10、12 会破坏“先在分布空间构造源，再乘实际 \(M\)”的一致性。

---

# 36. 任务 7：应变率的局部表达式

定义

\[
\Pi_{\alpha\beta}^{neq}
=
\sum_i c_{i\alpha}c_{i\beta}
(f_i-f_i^{eq}),
\]

\[
\mu=\rho_0\nu,
\qquad
\mu_B=\rho_0\nu_B.
\]

对非对角分量 \(\alpha\ne\beta\)，`LBM-CDE.pdf` 式 (31) 给出

\[
\boxed{
S_{\alpha\beta}
=
-
\frac{
2\Pi_{\alpha\beta}^{neq}
+\delta t(u_\alpha F_\beta+u_\beta F_\alpha)
}{
4\mu+2\delta t\rho_0c_s^2
}
}.
\]

速度散度为

\[
\boxed{
S_{\gamma\gamma}
=
-
\frac{
\Pi_{\gamma\gamma}^{neq}
+\delta t\,u_\gamma F_\gamma
}{
D\mu_B+\delta t\rho_0c_s^2
}
}.
\]

对角分量为

\[
\boxed{
S_{\alpha\alpha}
=
-
\frac{
\Pi_{\alpha\alpha}^{neq}
+\delta t\,u_\alpha F_\alpha
}{
2\mu+\delta t\rho_0c_s^2
}
-
\frac{
2\mu-D\mu_B
}{
D(2\mu+\delta t\rho_0c_s^2)
}
S_{\gamma\gamma}
}.
\]

这里使用的是物理黏度

\[
\mu
=
\rho_0(1-\chi_s)c_s^2
\sigma_\nu^{(0)}\delta t,
\]

不是未修正的基础黏度。

这些公式已经把应力源对非平衡矩的贡献包含在分母中，所以计算顺序可以直接写成

1. 由当前 \(f_i\)、\(\rho\)、\(\boldsymbol u\) 计算
   \(\Pi_{\alpha\beta}^{neq}\)；
2. 用上式一次性得到 \(S_{\alpha\beta}\)；
3. 构造 \(A_{\alpha\beta}\)；
4. 投影并加入 \((I-S/2)\widehat{\boldsymbol R}^m\)。

不需要 fixed-point 迭代，也没有显式--隐式循环。

---

# 37. 任务 8：温度梯度的局部表达式

定义

\[
\boldsymbol j_T^{neq}
=
\sum_i\boldsymbol e_i(g_i-g_i^{eq}).
\]

若采用完整弱可压缩源

\[
\boldsymbol A
=
\frac{p\nabla T+T\boldsymbol F}{\rho_0}
\]

并允许 \(Q_T\ne0\)，由 `LBM-CDE.pdf` 式 (35) 得

\[
\boxed{
\partial_\alpha T
=
-
\frac{
2j_{T,\alpha}^{neq}
+\delta t
\left(
TF_\alpha/\rho_0+u_\alpha Q_T
\right)
}{
\delta t
\left\{
c_T^2
\left[
\dfrac{2(1-\chi_\kappa)}
{q_\kappa^{(0)}}
+\chi_\kappa
\right]
+p/\rho_0
\right\}
}
}.
\]

主分支 A 令

\[
p/\rho_0\ \text{耦合项}=0,
\qquad
\boldsymbol F=\boldsymbol0,
\qquad
Q_T=0,
\]

于是

\[
\boxed{
\nabla T
=
-\frac{q_{\kappa,\mathrm{eff}}}
{c_T^2\delta t}
\boldsymbol j_T^{neq}
}.
\]

Fortran 可直接实现为

```fortran
jneq_x = (g(1)-geq(1)) - (g(3)-geq(3))
jneq_y = (g(2)-geq(2)) - (g(4)-geq(4))
gradT_x = -qk_eff*jneq_x/(cT2*dt)
gradT_y = -qk_eff*jneq_y/(cT2*dt)
```

若使用完整源，则必须按完整分母和力/热源分子实现，不能把上面的主分支简式用于自然对流力耦合项已经显式加入温度源的版本。

---

# 38. 任务 9：边界条件的源项一致性

## 38.1 统一时间推进次序

采用

\[
\boxed{
\text{宏观量}
\to
\text{局部梯度/应变率}
\to
\text{碰撞并加入半步源}
\to
\text{传播}
\to
\text{link 边界}
}.
\]

记

\[
f_i^*
\]

和

\[
g_i^*
\]

为已经包含 \((I-S/2)\widehat R\) 的完整 post-collision 分布。边界式必须反射这个完整值。

## 38.2 流场 halfway bounce-back

对从壁面进入流体的未知方向 \(i\)：

\[
\boxed{
f_i(\boldsymbol x_b,t+\delta t)
=
f_{\bar i}^*(\boldsymbol x_b,t)
+2w_i\rho_0
\frac{\boldsymbol c_i\cdot\boldsymbol u_w}{c_s^2}
}.
\]

静止壁面 \(\boldsymbol u_w=\boldsymbol0\) 时就是完整 post-collision 反弹。

不需要再手工增加第二个“源项反弹项”；若碰撞例程先生成无源
\(f^{coll}\)，而源项在别处另加，则必须保证被反弹的是

\[
f_{\bar i}^{coll}+\Phi_{\bar i},
\]

否则时间离散不一致。

壁面精确位置条件已经在第 34 节得到：

\[
\sigma_{\nu,\mathrm{eff}}\sigma_q=\frac3{16}.
\]

## 38.3 定温 anti-bounce-back

对静止定温壁面 \(T_w\)：

\[
\boxed{
g_i(\boldsymbol x_b,t+\delta t)
=
-g_{\bar i}^*(\boldsymbol x_b,t)
+2w_iT_w
}.
\]

D2Q5 轴向权重为

\[
w_i=\frac{c_T^2}{2},
\]

所以

\[
\boxed{
g_i
=
-g_{\bar i}^*
+c_T^2T_w
}.
\]

这正是

\[
\frac{4+a}{10}T_w
\]

形式。

对一维稳态线性温度

\[
T_j=T_w+G(j-1/2),
\]

体节点解满足

\[
j_x^{neq}
=
-\frac{c_T^2}{q_{\kappa,\mathrm{eff}}}G.
\]

代入 anti-bounce-back 可严格得到

\[
T_1-T_w=\frac G2,
\]

所以线性纯导热的壁面仍在半格点，且不单独依赖
\(\chi_\kappa\) 或 \(q_\kappa^{(0)}\)，只通过已合并的
\(q_{\kappa,\mathrm{eff}}\) 进入体节点分布。

该结论不等价于“任意曲率温度场都无边界误差”。当
\(\partial_n^2T\ne0\)、\(Q_T\ne0\) 或使用完整
\(p\nabla T+T\boldsymbol F\) 源时，仍需用 manufactured solution
测量温度滑移。

## 38.4 绝热 bounce-back

绝热壁面采用

\[
\boxed{
g_i(\boldsymbol x_b,t+\delta t)
=
g_{\bar i}^*(\boldsymbol x_b,t)
}.
\]

在一维稳态情形代入通量矩更新可得

\[
(2-q_{\kappa,\mathrm{eff}})j_n^{neq}=0.
\]

只要

\[
q_{\kappa,\mathrm{eff}}\ne2,
\]

即得到

\[
j_n^{neq}=0
\quad\Longleftrightarrow\quad
\partial_nT=0.
\]

当 \(q_{\kappa,\mathrm{eff}}\to2\) 时，该约束退化，这与小扩散率下通量矩弱阻尼是同一个极限问题。

## 38.5 Nusselt 数

物理定义不变：

\[
Nu_w
=
-\frac{L}{\Delta T}
\left.\frac{\partial T}{\partial n}\right|_w.
\]

但后处理必须与碰撞使用同一个局部梯度公式。不能在碰撞中使用
\(q_{\kappa,\mathrm{eff}}\) 的非平衡梯度，却在壁面 \(Nu\) 中继续使用旧
\(q_\kappa^{(0)}\) 系数。

---

# 39. 任务 10：是否升级 D2Q9 的最终判定

本轮推导得到的是条件判定，而不是无条件升级：

\[
\boxed{
\begin{array}{l}
\text{若目标是稳定的二阶高 }Ra\text{ 算法，先保留 D2Q5；}\\
\text{若必须同时满足严格四阶各向同性且 }
h<h_{\min}^{(4)},\text{ 则升级 D2Q9。}
\end{array}
}
\]

D2Q9 的启动条件为以下任一项：

1. 目标 \(h\) 低于
   \[
   \varepsilon/[3(2-\varepsilon)]
   \]
   且仍要求严格四阶各向同性；
2. D2Q5 的交叉扩散
   \[
   -h\,uv\,\partial_{xy}T
   \]
   在对角对流测试中不可接受；
3. 完整谱扫描在目标速度范围内不存在满足
   \(q_3,q_4\in[\varepsilon,2-\varepsilon]\) 的稳定点；
4. ABB/BB 边界误差在目标参数下无法通过网格收敛控制。

D2Q9 回退分支应至少采用能表示完整二阶速度张量的平衡分布：

\[
g_i^{eq}
=
w_iT
\left[
1
+\frac{\boldsymbol c_i\cdot\boldsymbol u}{c_T^2}
+\frac{
(\boldsymbol c_i\cdot\boldsymbol u)^2
-c_T^2|\boldsymbol u|^2
}{
2c_T^4
}
\right],
\]

从而局部表示 \(u_xu_y\) 并消除 D2Q5 缺失对角速度导致的结构性限制。

是否采用该分支必须由第 25 节的恒定速度对流扩散和自然对流方腔验证决定；当前推导没有证明 D2Q5 的二阶稳定分支失败，因此第一版实现仍应保留 D2Q5。

---

# 40. 可直接转写为 Fortran 的主分支 A

## 40.1 参数预处理

```fortran
cT2    = (4.0d0 + paraA)/10.0d0
sigma0 = 1.0d0/qk0 - 0.5d0
hEff   = diffusivity/(cT2*dt)
chiK   = 1.0d0 - hEff/sigma0
qkEff  = 1.0d0/(hEff + 0.5d0)
```

必须检查

```fortran
if (cT2 <= 0.0d0) error stop "cT2 must be positive"
if (qk0 <= 0.0d0 .or. qk0 >= 2.0d0) error stop "invalid qk0"
if (chiK < 0.0d0 .or. chiK >= 1.0d0) error stop "invalid chiK"
```

## 40.2 逐格点碰撞

```fortran
! n = N*g
n0 = g0 + g1 + g2 + g3 + g4
n1 = g1 - g3
n2 = g2 - g4
n3 = -4.0d0*g0 + g1 + g2 + g3 + g4
n4 = g1 - g2 + g3 - g4

neq0 = T
neq1 = T*u
neq2 = T*v
neq3 = paraA*T
neq4 = 0.0d0

jneqX = n1 - neq1
jneqY = n2 - neq2
gradTx = -qkEff*jneqX/(cT2*dt)
gradTy = -qkEff*jneqY/(cT2*dt)

r1 = chiK*cT2*dt*gradTx
r2 = chiK*cT2*dt*gradTy

npost0 = n0
npost1 = n1 - qk0*(n1-neq1) + (1.0d0-0.5d0*qk0)*r1
npost2 = n2 - qk0*(n2-neq2) + (1.0d0-0.5d0*qk0)*r2
npost3 = n3 - q3*(n3-neq3)
npost4 = n4 - q4*(n4-neq4)
```

上述 `npost1/npost2` 可用

```fortran
npost1 = n1 - qkEff*(n1-neq1)
npost2 = n2 - qkEff*(n2-neq2)
```

逐点验证；两者在主分支 A 中应只相差舍入误差。这个对照是最便宜的实现自检。

## 40.3 物理与数值不变量

实现后必须逐项检查：

\[
\sum_i\Psi_i=0,
\]

\[
\chi_\kappa=0
\Rightarrow
q_{\kappa,\mathrm{eff}}=q_\kappa^{(0)},
\]

\[
\boldsymbol u=0
\Rightarrow
\partial_tT=\kappa\nabla^2T+O(\delta x^4),
\]

\[
\kappa\to0
\Rightarrow
q_{\kappa,\mathrm{eff}}\to2,
\]

以及

\[
\max_\alpha|u_\alpha|\le c_T^2,
\qquad
|\boldsymbol u|^2\le c_T^2.
\]

---

# 41. 完成状态、验证边界与文献依据

## 41.1 已由代数和符号计算闭合

- D2Q5 源项的分布/矩空间分量；
- 半步源离散；
- 局部梯度与 \(q_{\kappa,\mathrm{eff}}\) 的代数等价；
- 纯扩散四阶 \(K_{40},K_{22}\)；
- \(\chi_\kappa=0\) 退化；
- 四阶各向同性和完全消除条件；
- 常对流二至四阶 Fourier 展开；
- D2Q9/D3Q19 实际矩阵的完整源项投影；
- 应变率和温度梯度局部表达；
- 标准 HBB 下的有效 \(3/16\) 关系；
- 线性纯导热 ABB 和绝热 BB 的半格点/零通量闭合。

## 41.2 尚需数值验证，不能写成已验证性能

- 非线性自然对流中完整 \(p\nabla T+T\boldsymbol F\) 温度源的四阶误差；
- 曲率温度场的 ABB 隐藏误差；
- 高 \(Ra\) 实际速度范围内的全谱稳定窗口；
- 新流场源投影在 D3Q19 上的 Poiseuille/Couette 壁面误差；
- OpenACC/MPI 实现的设备正确性和运行稳定性；
- \(Nu\)、\(Re\) 与基准数据的一致性。

本轮没有修改求解器，也没有完成 GPU/高 \(Ra\) 运行。因此这里的“完成”指推导任务已经闭合，并不等于新算法已经通过算例验证。

## 41.3 主要原始文献

1. `pdf/LBM-CDE.pdf`：式 (24)--(29) 的源项与输运系数映射；式 (31)--(35) 的局部应变率/梯度；式 (36)--(39) 的边界形式。
2. `pdf/Towards higher order lattice Boltzmann schemes.pdf`：D2Q5 纯扩散式 (40)--(42) 与四阶参数式 (55)。
3. `pdf/Multireflection boundary conditions for lattice Boltzmann models.pdf`：link 边界、偶/奇碰撞参数与壁面位置的关系。
4. `pdf/Optimal Stability of Advection-Diffusion Lattice Boltzmann Models with Two Relaxation Times for Positive_Negative Equilibrium.pdf`：人口正负、长波条件和全谱稳定性之间的区别。
5. `均匀网格/2DRBOpenmp.F90`：本推导采用的实际 D2Q9/D2Q5 矩顺序。
6. `均匀网格/3DRBOpenmp.F90`：本推导采用的实际 D3Q19 矩顺序。
7. `文档/high_ra_derivation_source_matrix.md`：逐任务区分文献直接支持、本文新推导和仍需数值验证的结论。

---

# 附录 A：恒定对流四阶系数完整展开

本附录采用

\[
h=\sigma_{\kappa,\mathrm{eff}},
\qquad
r=\sigma_3,
\qquad
t=\sigma_4.
\]

## A.1 \(e_{40}\)

\[
\begin{aligned}
e_{40}=\frac1{1200}\big[
&-12a^2h^3-12a^2h^2r+3a^2h
+720ah^3u^2-96ah^3\\
&+360ah^2ru^2-36ah^2r+60ah^2t
+120ahr^2u^2\\
&-150ahu^2+4ah-30aru^2
-6000h^3u^4+2880h^3u^2\\
&-192h^3-360h^2ru^2+48h^2r
-1800h^2tu^2+240h^2t\\
&-120hr^2u^2-600ht^2u^2
+900hu^4-100hu^2-32h\\
&+30ru^2+150tu^2
\big].
\end{aligned}
\]

## A.2 \(e_{31}\)

\[
\begin{aligned}
e_{31}
=\frac{uv}{60}\big[
&72ah^3+36ah^2r+12ahr^2-15ah-3ar\\
&-1200h^3u^2+288h^3-36h^2r-12hr^2\\
&+180hu^2-35h+3r
\big].
\end{aligned}
\]

## A.3 \(e_{22}\)

\[
\begin{aligned}
e_{22}=\frac1{200}\big[
&-4a^2h^3-4a^2h^2r+a^2h\\
&+120ah^3(u^2+v^2)-32ah^3
+60ah^2r(u^2+v^2)-12ah^2r\\
&-20ah^2t
+20ahr^2(u^2+v^2)
-25ah(u^2+v^2)+8ah\\
&-5ar(u^2+v^2)
-6000h^3u^2v^2
+480h^3(u^2+v^2)-64h^3\\
&-60h^2r(u^2+v^2)+16h^2r
+300h^2t(u^2+v^2)-80h^2t\\
&-20hr^2(u^2+v^2)
+100ht^2(u^2+v^2)
+900hu^2v^2\\
&-100h(u^2+v^2)+16h
+5r(u^2+v^2)-25t(u^2+v^2)
\big].
\end{aligned}
\]

## A.4 \(e_{13}\)

\[
\begin{aligned}
e_{13}
=\frac{uv}{60}\big[
&72ah^3+36ah^2r+12ahr^2-15ah-3ar\\
&-1200h^3v^2+288h^3-36h^2r-12hr^2\\
&+180hv^2-35h+3r
\big].
\end{aligned}
\]

## A.5 \(e_{04}\)

\[
e_{04}=e_{40}\big|_{u\rightarrow v}.
\]

这些系数由

\[
\log\Lambda_h
=
-i(uk_x+vk_y)
-h[c_T^2|\boldsymbol k|^2-(\boldsymbol u\cdot\boldsymbol k)^2]
+z_3+z_4+O(|\boldsymbol k|^5)
\]

定义。符号生成和纯扩散退化断言位于
`tools/derive_d2q5_cde.py`。
