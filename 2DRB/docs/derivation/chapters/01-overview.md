# 基于 LBM-CDE 的 TRT 推广、边界 magic 与温度四阶误差推导

## 摘要

本文从用户提供的 `LBM-CDE.pdf` 的梯形 transformed 离散方程出发，重新完成
奇偶分解、TRT 源项、局部反馈消元、二阶宏观恢复、源感知边界残差与 D2Q9
四阶等效方程推导。目标不是把 BGK 公式中的单一松弛率机械替换成两个数，
而是回答三个可验证的问题：TRT 推广是否保持原 LBM-CDE 的宏观方程；所谓
magic 参数在加入局部反馈、压力修正和热源后还剩多大适用域；体相四阶消除、
低扩散率与边界标定能否同时满足。

全文所有主要结论分为四级：`strictly_proved` 表示已由精确代数与至少一条独立
路线闭环；`restricted_model` 表示只在明确的平壁、稳态、冻结系数或二次场假设
下成立；`numerical_evidence` 表示高精度 Fourier 或有限算例抽查；`unresolved`
表示一般源项/边界 jet 尚未闭合。`Xs/` 仅是此前复现代码的实现比较对象，不作为
本文理论证明的依据。

## 主要结论

1. **BGK 推广到 TRT 是可行的，但源项必须按奇偶子空间分别乘因子。**
   对 transformed 分布，正确碰撞项包含
   $(1-s_+/2)\Phi_i^+$ 与 $(1-s_-/2)\Phi_i^-$。当 $s_+=s_-$ 时逐分量退化回
   LBM-CDE 的 BGK 形式；因此这是算子梯形规则的推广，而非经验替换。

2. **局部反馈改变物理块的有效 Hénon shift，不会自动改变同奇偶性的 ghost。**
   冻结系数、CE2 范围内，剪切、体积和温度通量块分别满足
   $\sigma_{\rm sh}^{\rm eff}=(1-\chi_s)\sigma_f^+$、
   $\sigma_{\rm bu}^{\rm eff}=(1-\chi_b)\sigma_f^+$ 与
   $\sigma_{\rm q}^{\rm eff}=b\sigma_g^-/a$，其中
   $a=c_s^2+\pi$、$b=(1-\chi_\kappa)c_s^2$。名义 ghost shift 仍是碰撞输入本身。

3. **不存在脱离边界实现与源规范的普适 magic 常数。**
   $\sigma_f^+\sigma_f^-=3/16$ 只对应平直、格点对齐、halfway、稳态 Stokes、
   均匀体力、half-force 重构与无 feedback 下的受限标定；压力边界驱动文献中的 $3/8$ 是另一
   gauge。一般压力、力、时间和法/切向导数 jets 仍要求显式边界修正。

4. **D2Q9 温度四阶条件必须针对实际 external 或 local-feedback 方案重推。**
   两种方案二阶扩散率都为 $\kappa=b\sigma_o\Delta t$，但其四阶系数一般不同。
   本文保留 $\lambda_i$ 的非零四阶矩，并用放大矩阵流体根与物理空间 Taylor/矩
   消元两条路线得到 $C_{22}=2C_{40}$ 及各自的完整消除条件。

5. **体相四阶消除与受限 ABB 标定会锁定扩散尺度。**
   令 $K=\kappa/\Delta t$，与受限 ABB 条件联立后，local-feedback 与 external
   分别给出
   $$
   K^2=\frac{a(3a+1)}{48},\qquad
   K^2=\frac{12ab+5a-4b-9a^2}{48}.
   $$
   这里的 ABB 条件仅指 R8 子域：平直、格点对齐的 halfway 壁，冻结压力、稳态
   一维二次温度场、half-source 重构、完整 pressure-wall equilibrium term 与匹配的
   external/local population chain、CDE-consistent 均匀 $Q$，且无流/力/切向 jets；
   它不是一般边界恒等式。
   在 $\pi=0$ 的正支，两式都退化为 $K=1/\sqrt{72}$。因此“任意低 $\kappa$、
   体相四阶完全消除、同一个受限 ABB magic”三项通常不可同时强制。

6. **低扩散率下的诚实实现选择有三类。**
   一条已识别但仍待推导、实现和验证的方向，是保持体相四阶条件并加入显式受限 ABB 修正；
   也可保留受限 ABB 标定而
   接受非零 $q^4$ 误差；若要用 split-even MRT 同时释放两项约束，必须另行证明
   模态归属与参数 Jacobian 满秩，本文只把它列为待推导候选。

## 记号与目标方程

速度格点为 $\boldsymbol c_i$，反向指标为 $\bar i$，奇偶部分定义为
$$
h_i^+=\frac{h_i+h_{\bar i}}2,\qquad
h_i^-=\frac{h_i-h_{\bar i}}2.
$$
本文 D2Q9 速度顺序固定为
$(0,0),(1,0),(0,1),(-1,0),(0,-1),(1,1),(-1,1),(-1,-1),(1,-1)$；
相应 LBM-CDE 压力修正向量由低阶矩约束显式得到
$$
\boldsymbol\lambda_t=
\left(-\frac59,\frac19,\frac19,\frac19,\frac19,
\frac1{36},\frac1{36},\frac1{36},\frac1{36}\right)^{\mathsf T}.
$$
Hénon shift 与松弛率的关系为
$$
\sigma=\frac1s-\frac12,\qquad s=\frac1{\sigma+1/2}.
$$
本文区分名义碰撞 shift、局部反馈消元后的物理有效 shift、受限边界乘积和一般
边界 residual；四者不能互换。

目标低 Mach、弱可压宏观方程直接沿用第四章的恢复形式：
$$
\partial_t\delta\rho+\partial_\alpha(\rho_0u_\alpha)=O(\epsilon^3),
$$
$$
\partial_t(\rho_0u_\alpha)+\partial_\beta(\rho_0u_\alpha u_\beta)
=-\partial_\alpha p+\partial_\beta\tau_{\alpha\beta}+F_\alpha
+O(\epsilon^3)+O(\mathrm{Ma}^3),
$$
$$
\tau_{\alpha\beta}
=2\rho_0\nu\left(S_{\alpha\beta}-\frac1D S_{\gamma\gamma}
\delta_{\alpha\beta}\right)
+\rho_0\nu^B S_{\gamma\gamma}\delta_{\alpha\beta},
$$
$$
\nu=(1-\chi_s)c_s^2\Delta t\,\sigma_f^+,
\qquad
\nu^B=\frac2D(1-\chi_b)c_s^2\Delta t\,\sigma_f^+,
$$
以及
$$
\partial_tT+\nabla\!\cdot(T\boldsymbol u)
=\nabla\!\cdot(\kappa\nabla T)+Q.
$$
二阶恢复章节会逐项说明梯形半源重构、压力修正和冻结系数余项；这些方程不被
直接假定为结果。

## transformed TRT 离散骨架

以流场为例，一步碰撞--迁移写成
$$
\widetilde f_i(\boldsymbol x+\boldsymbol c_i\Delta t,t+\Delta t)
-\widetilde f_i(\boldsymbol x,t)
=\Omega_i^{\rm TRT},
$$

$$
\Omega_i^{\rm TRT}
=-s_f^+\bigl(\widetilde f_i^+-f_i^{\mathrm{eq},+}\bigr)
-s_f^-\bigl(\widetilde f_i^--f_i^{\mathrm{eq},-}\bigr)
+\Delta t\,\Psi_i,
$$

其中 parity-specific 梯形源增量为

$$
\Psi_i=\left(1-\frac{s_f^+}{2}\right)\Phi_i^+
+\left(1-\frac{s_f^-}{2}\right)\Phi_i^-.
$$
温度分布 $\widetilde g_i$ 具有同一算子结构，只需把平衡态、源矩与松弛率替换为
标量对应量。后续章节从这一离散式推导有效率、宏观方程和高阶误差。

## 正文组织

完整 Markdown 与 PDF 按下列顺序装配 reviewed chapter：

1. 奇偶源项与算子梯形 TRT 碰撞；
2. 局部反馈消元与物理模态有效率；
3. 二阶 Chapman--Enskog/等效方程恢复与 BGK 极限；
4. 源感知速度、ABB、绝热与混合角点边界；
5. D2Q9 温度四阶等效方程；
6. 参数兼容性、不可兼容性与低扩散率选择；
7. 实现顺序、伪代码与运行时不变量；
8. D2Q5 文献公式复算附录；
9. 证据矩阵、限制与参考文献。

下文直接给出全部详细公式、适用假设、程序入口与测试映射，不省略边界残差和
四阶推导细节。各分章文件同时保留为可独立审查的生成源。
