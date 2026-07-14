# Boundary magic 独立复算

## 1. 输入身份

论文直接输入：

- Ginzburg Eq. (23) 的 BB 系数；
- Eq. (9) 的 $I_f=-1/2$；
- Eq. (10) 的黏性关系；
- Eq. (40) 的通用有效宽度式；
- Eq. (41) 的 $\Lambda^2$ 定义。

本文件只做独立代回和代数复算，不声称重新证明 Eq. (40)。

## 2. Eq. (40) 五项代回

使用

\[
\delta_q=\frac12,
\quad I_f=-\frac12,
\quad\Theta_q^2=1,
\]

\[
\kappa_{-1}=\bar\kappa_{-1}=\bar\kappa_{-2}=0,
\]

\[
\nu=\frac13H_{\mathrm{even}},
\quad
H_{\mathrm{odd}}=\frac1{\lambda_2}-\frac12,
\quad
\Lambda^2=\frac43H_{\mathrm{even}}H_{\mathrm{odd}}.
\]

五项分别为

\[
6\Lambda^2,
\]

\[
-4\delta_q^2=-1,
\]

\[
\begin{aligned}
-\frac{8\nu}{\Theta_q^2}
\left(\frac1{\lambda_2}+I_f\right)
&=-8\nu H_{\mathrm{odd}}\\
&=-\frac83H_{\mathrm{even}}H_{\mathrm{odd}}\\
&=-2\Lambda^2,
\end{aligned}
\]

\[
\frac{4(0-0)}{1-0-0}=0,
\]

\[
\frac{12\nu}{1-0-0}
\left[1-2\left(\frac12\right)+2(0)+4(0)\right]=0.
\]

总和：

\[
6\Lambda^2-1-2\Lambda^2+0+0
=4\Lambda^2-1.
\]

因此

\[
D_{\mathrm{eff}}^2
=D_{1/2}^2+4\Lambda^2-1.
\]

## 3. $1/4$ 与 $3/16$

正直径条件下：

\[
D_{\mathrm{eff}}=D_{1/2}
\Longleftrightarrow
4\Lambda^2-1=0
\Longleftrightarrow
\Lambda^2=\frac14.
\]

再由 Eq. (41)：

\[
\Lambda^2=\frac43H_{\mathrm{even}}H_{\mathrm{odd}}
=\frac14
\Longleftrightarrow
H_{\mathrm{even}}H_{\mathrm{odd}}=\frac3{16}.
\]

反向检查：

\[
\frac43\cdot\frac3{16}=\frac14.
\]

## 4. 率关系双向复算

记

\[
a=r_{\mathrm{even}},
\qquad b=r_{\mathrm{odd}}.
\]

正向：

\[
\frac{2-a}{2a}\frac{2-b}{2b}=\frac3{16},
\]

\[
4(2-a)(2-b)=3ab,
\]

\[
16-8a-8b+ab=0,
\]

\[
b(a-8)=8(a-2),
\]

\[
\boxed{b=8\frac{2-a}{8-a}}.
\]

反向：设

\[
b=8\frac{2-a}{8-a}.
\]

则

\[
\begin{aligned}
H_{\mathrm{odd}}
&=\frac1b-\frac12\\
&=\frac{8-a}{8(2-a)}-\frac12\\
&=\frac{8-a-4(2-a)}{8(2-a)}\\
&=\frac{3a}{8(2-a)}.
\end{aligned}
\]

乘 $H_{\mathrm{even}}=(2-a)/(2a)$：

\[
H_{\mathrm{even}}H_{\mathrm{odd}}
=\frac{2-a}{2a}\frac{3a}{8(2-a)}
=\frac3{16}.
\]

定义域：

\[
a\ne0,
\quad b\ne0,
\quad a\ne8,
\quad a\ne2.
\]

$a=2$ 被单独排除，因为率公式给 $b=0$，Hénon 参数无定义。

## 5. 表观滑移复算

\[
u(x)=u_0\left(1-\frac{4x^2}{D_{\mathrm{eff}}^2}\right),
\quad
u_0=\frac{FD_{\mathrm{eff}}^2}{8\nu}.
\]

在名义壁面 $x=D_{1/2}/2$：

\[
\begin{aligned}
u_{\mathrm{slip}}
&=u_0\left(1-\frac{D_{1/2}^2}{D_{\mathrm{eff}}^2}\right)\\
&=\frac{F}{8\nu}
\left(D_{\mathrm{eff}}^2-D_{1/2}^2\right)\\
&=\frac{F}{8\nu}(4\Lambda^2-1).
\end{aligned}
\]

在 $F\ne0$、$\nu\ne0$ 的非平凡受迫剖面中：

\[
u_{\mathrm{slip}}=0
\Longleftrightarrow
\Lambda^2=\frac14.
\]

## 6. 范围检查

此复算不覆盖曲/斜壁、非稳态、非线性、其他 forcing gauge、温度 ABB、一般 Knudsen layer、全局收敛阶或稳定性。
