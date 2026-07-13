# LBM-CDE 源反馈 TRT 推导与验证设计

**日期：** 2026-07-13
**状态：** 已由用户批准
**主模型：** D2Q9 流场分布 \(f_i\) + D2Q9 温度分布 \(g_i\)
**辅助模型：** D2Q5 温度模型仅用于符号程序退化校验，不作为 D2Q9 参数来源

## 1. 研究目标

本研究从 pdf/LBM-CDE.pdf 的逆向设计路线出发，不直接移植其他 TRT 模型的现成参数，完成以下目标：

1. 将 LBM-CDE 的双 BGK 碰撞严格推广为流场、温度场各自独立的 TRT 碰撞骨架。
2. 从连续源项矩约束、奇偶投影和梯形离散重新推出 TRT 版离散源项。
3. 证明新模型在二阶、弱可压缩尺度下恢复原目标 N-S/C-D 方程，并给出物理黏度、热扩散率与参数的映射。
4. 消去局部非平衡应变率和温度梯度反馈，辨识真正作用于物理矩的有效松弛率。
5. 针对当前 D2Q9 平衡态、压力耦合、源项和边界规则重新推导速度与温度的边界残差及 magic 条件。
6. 针对当前 D2Q9 温度模型重新生成四阶等效方程，区分四阶各向同性与四阶项全消除。
7. 联立物理输运率、边界条件、四阶条件和稳定区间，判断源反馈 TRT 是否存在可接受参数族。
8. 形成完整中文推导报告、可复现符号验证程序、最小数值验证及实现伪代码。

研究不得预设“纯 TRT 一定能同时满足全部条件”，也不得预设“必须升级为显式 MRT”。兼容性必须由最终方程和验证结果决定。

## 2. 权威资料和使用边界

理论依据按以下优先级使用：

1. pdf/LBM-CDE.pdf：目标宏观方程、逆向源项设计、梯形离散、局部应变率/梯度恢复和原始边界规则。
2. pdf/Multireflection boundary conditions for lattice Boltzmann models.pdf：链式边界 Taylor 分析和壁面位置误差方法。
3. pdf/Towards higher order lattice Boltzmann schemes.pdf：递归等效方程、D2Q5/D2Q9 高阶分析方法。
4. Wang 2013 与 Luo 2014 JCP：参数关系和热对流实现的交叉检查。

Xs/2DRBOpenmpLBMCDE.F90 与 Xs/2DRBOpenaccLBMCDE.F90 是此前按论文复现的实现，只能用于核对公式落地和设计验证接口，不能作为理论证明依据。

## 3. 统一记号与模型结构

对每一反向速度对 \(\mathbf c_{\bar i}=-\mathbf c_i\)，定义

\[
h_i^+=\frac{h_i+h_{\bar i}}2,\qquad
h_i^-=\frac{h_i-h_{\bar i}}2,
\qquad h\in\{f,g,S,R\}.
\]

定义投影算子 \(P^+h=h^+\)、\(P^-h=h^-\)，以及 Hénon 参数

\[
\sigma_h^\pm=\frac1{s_h^\pm}-\frac12.
\]

流场和标量场使用互不绑定的四个名义松弛率：

\[
(s_f^+,s_f^-),\qquad (s_g^+,s_g^-).
\]

连续候选方程为

\[
D_tf=-\left(\frac{P^+}{\tau_f^+}+\frac{P^-}{\tau_f^-}\right)
(f-f^{eq})+S,
\]

\[
D_tg=-\left(\frac{P^+}{\tau_g^+}+\frac{P^-}{\tau_g^-}\right)
(g-g^{eq})+R.
\]

流场二阶应力为偶矩，与 \(s_f^+\) 关联；温度一阶扩散通量为奇矩，与 \(s_g^-\) 关联。守恒矩必须通过平衡态矩和半源重构严格配平。

## 4. 源项奇偶分解

### 4.1 流场

\[
S_i^-=w_i\frac{\mathbf c_i\cdot\mathbf F}{c_s^2},
\]

\[
S_i^+=w_i\left[
\frac{(\mathbf c_i\cdot\mathbf u)(\mathbf c_i\cdot\mathbf F)}{c_s^4}
-\frac{\mathbf u\cdot\mathbf F}{c_s^2}
+\rho_0H_i^{(2)}:A
\right],
\]

\[
A_{\alpha\beta}=\chi_sS_{\alpha\beta}
+\frac{\chi_b-\chi_s}{D}S_{\gamma\gamma}\delta_{\alpha\beta}.
\]

外力的一阶矩在奇源项；Guo 型 \(uF\) 修正、剪切和体黏度调节均在偶源项。

### 4.2 温度

\[
R_i^+=w_iQ,
\]

\[
R_i^-=w_i\left[
\frac{\mathbf c_i\cdot(p\nabla T+T\mathbf F)}{\rho_0c_s^2}
+Q\frac{\mathbf c_i\cdot\mathbf u}{c_s^2}
+\chi_\kappa\mathbf c_i\cdot\nabla T
\right].
\]

\(R_i^+\) 除零阶矩 \(Q\) 外还具有二阶矩 \(c_s^2Q\delta_{\alpha\beta}\)，高阶展开必须保留其产生的 \(\nabla^2Q\) 等项。

## 5. 算子化梯形离散

不得把 BGK 的单一 sourcePref 直接沿用到 TRT。离散碰撞必须为

\[
\widetilde f_i^*=\widetilde f_i
-s_f^+(\widetilde f_i^+-f_i^{eq,+})
-s_f^-(\widetilde f_i^--f_i^{eq,-})
+\Delta t\,\Phi_i,
\]

\[
\Phi_i=\left(1-\frac{s_f^+}{2}\right)S_i^+
+\left(1-\frac{s_f^-}{2}\right)S_i^-,
\]

\[
\widetilde g_i^*=\widetilde g_i
-s_g^+(\widetilde g_i^+-g_i^{eq,+})
-s_g^-(\widetilde g_i^--g_i^{eq,-})
+\Delta t\,\Psi_i,
\]

\[
\Psi_i=\left(1-\frac{s_g^+}{2}\right)R_i^+
+\left(1-\frac{s_g^-}{2}\right)R_i^-.
\]

宏观量继续采用半源重构：

\[
\delta\rho=\sum_i\widetilde f_i,\qquad
\rho_0\mathbf u=\sum_i\mathbf c_i\widetilde f_i
+\frac{\Delta t}{2}\mathbf F,
\]

\[
T=\sum_i\widetilde g_i+\frac{\Delta t}{2}Q.
\]

一次碰撞更新的净动量源和净标量源必须严格等于 \(\Delta t\mathbf F\) 和 \(\Delta t Q\)，且不得依赖松弛率。

## 6. 源反馈与有效松弛率

名义 TRT 碰撞率不能直接当作边界与高阶分析中的物理矩松弛率。局部非平衡应变率和温度梯度代回源项后，源反馈只改变选定物理矩，导致同一奇偶扇区内的有效率分裂。

纯热、常压力退化问题中，

\[
\sigma_{1,\mathrm{eff}}=(1-\chi_\kappa)\sigma_g^-,
\]

而三阶奇 ghost 矩仍由名义 \(\sigma_g^-\) 控制，偶 ghost 矩由 \(\sigma_g^+\) 控制。非零压力时，有效率还可能含

\[
\sigma_{1,\mathrm{eff}}
=\frac{(1-\chi_\kappa)c_s^2}{c_s^2+p/\rho_0}\sigma_g^-.
\]

流场剪切、体积应力也会因 \(\chi_s,\chi_b\) 反馈而分别形成有效率。因此算法称为“LBM-CDE-TRT 源反馈模型”；从等效碰撞算子看，它是 TRT 骨架上的受限块-MRT。

二阶物理输运映射为

\[
\nu=(1-\chi_s)c_s^2\Delta t\,\sigma_f^+,
\qquad
\kappa=(1-\chi_\kappa)c_s^2\Delta t\,\sigma_g^-.
\]

这些关系必须由完整二阶离散展开和局部闭合重新证明。

## 7. 证明任务

### 7.1 二阶内点恢复

离散 Chapman-Enskog 或等价 Taylor 展开必须证明：

- 恢复连续方程、N-S 方程和 C-D 方程至论文声明的 \(O(Ma^2)\)；
- \(\nu\) 中不出现 \(s_f^-\)，\(\kappa\) 中不出现 \(s_g^+\)；
- \(p\nabla T\)、\(T\mathbf F\)、\(Q\mathbf u\)、\(uF\) 项正确抵消；
- 时变 \(F,Q\) 的半步离散项处理正确；
- 明确列出未消除的 \(O(\epsilon^3)\)、\(O(Ma^3)\) 项。

### 7.2 BGK 退化

当

\[
s_f^+=s_f^-=\tau_{fL}^{-1},\qquad
s_g^+=s_g^-=\tau_{gL}^{-1},
\]

模型必须逐项恢复 LBM-CDE 的 BGK 离散方程、源项前因子、宏观重构和局部闭合公式。

### 7.3 边界

分别推导：

1. 速度 halfway bounce-back；
2. 温度 Dirichlet anti-bounce-back；
3. 温度绝热 bounce-back；
4. 混合角点。

一般残差保留

\[
\mathcal E_w=C_2\partial_n^2\phi
+C_p\partial_np+C_FF_n+C_Q\partial_nQ
+C_t\partial_t\phi+C_{nt}\partial_n\partial_t\phi+\cdots.
\]

只有所有独立系数能被同一条件消除时，才称为完整 magic 条件。一维候选式

\[
(1-\chi_s)\sigma_f^+\sigma_f^-=\frac{3}{16},
\qquad
(1-\chi_\kappa)\sigma_g^-\sigma_g^+=\frac{3}{16}
\]

只能作为退化制造解校准目标。

### 7.4 D2Q9 温度四阶误差

\[
\partial_tT+\nabla\cdot(T\mathbf u)
=\kappa\nabla^2T+Q
+\Delta t^3\left[
C_{40}(\partial_x^4+\partial_y^4)T
+C_{22}\partial_x^2\partial_y^2T
+\mathcal E_{p,u,F,Q}
\right]+\cdots.
\]

\[
C_{22}=2C_{40}
\]

只表示四阶各向同性，而

\[
C_{40}=C_{22}=0
\]

才表示四阶项全消除。

Dubois-Lallemand 的 D2Q9 印刷版 \(\kappa_{22}\) 与后续 TRT 参数存在代数不一致，必须从放大矩阵和递归 Taylor 展开重新生成，不能静默修正或直接引用。

## 8. 双证据验证

### 8.1 常系数离散放大矩阵

\[
A(\mathbf k)=E(\mathbf k)C,\qquad
E_{ii}=\exp(-i\mathbf k\cdot\mathbf c_i\Delta t).
\]

对温度水动力特征值 \(z_h\) 展开

\[
\Gamma(\mathbf k)=-\frac1{\Delta t}\log z_h
=\kappa(k_x^2+k_y^2)
+C_{40}(k_x^4+k_y^4)+C_{22}k_x^2k_y^2+O(k^6).
\]

轴向 \((q,0)\) 与等模长对角向 \((q/\sqrt2,q/\sqrt2)\) 必须同时验证。全消除后残差应由 \(q^4\) 转为 \(q^6\)。

### 8.2 递归 Taylor 生成器

以 \(D_i=\partial_t+\mathbf c_i\cdot\nabla\) 展开 streaming，递归消去非守恒矩。每个系数标明来自 equilibrium、\(R^+\)、\(R^-\)、梯度反馈或时间导数替换。

常系数交集上 Taylor 系数必须与放大矩阵逐项一致；变压力、变速度、变外力和变热源由 Taylor 推导为主，以侧带 Fourier 块矩阵抽查。

### 8.3 D2Q5 校验

符号生成器先精确复现 D2Q5 已知四阶等效方程、各向同性和全消除关系，再用于 D2Q9。D2Q5 的 \(1/6\) 等参数不得移植到 D2Q9。

## 9. 参数兼容性与决策

最终联立：

- \(0<s_f^\pm,s_g^\pm<2\)；
- \(\nu>0,\kappa>0\)；
- 目标 \(\nu,\kappa\)；
- 速度、温度边界残差；
- D2Q9 四阶条件；
- 大 \(\chi\)、小松弛率的误差和稳定性限制。

优先级为：守恒与目标方程、正确物性、稳定区间、边界一致性、四阶优化。

若可行域非空，给出参数族和推荐选择；若为空，给出无解证明，并只释放满足缺失约束所需的最少额外矩松弛率。

## 10. 多智能体交叉审稿

三名子智能体分别负责：

- A：连续 TRT、源项投影、梯形离散、二阶宏观恢复；
- B：边界递推、magic 参数、角点闭合；
- C：D2Q5/D2Q9 四阶等效方程、放大矩阵和参数可行域。

每个主要推导经过主推导、敌对审查、独立复算三轮。主智能体统一符号、复核程序输出并裁决冲突。任何单一子智能体的结论都不能直接进入最终报告。

## 11. 最终交付物

1. docs/lbm-cde-trt-derivation.md：完整中文推导报告。
2. output/pdf/lbm-cde-trt-derivation.pdf：排版并经页面渲染检查的 PDF。
3. tools/derivation/：符号矩、放大矩阵、Taylor 递归和参数求解程序。
4. tests/derivation/：D2Q5、D2Q9、BGK 退化、Fourier 与边界制造解检查。
5. 报告中的实现伪代码。

本阶段不修改生产 Fortran 求解器。只有理论条件、参数可行域和最小验证通过后，才建议创建新求解器变体。

## 12. 完成标准

- 完成 BGK 到源反馈 TRT 的逐步推导；
- 二阶宏观恢复和 BGK 退化有符号或逐分量证据；
- 边界结论区分通用、限模型和无通用条件；
- D2Q9 四阶系数由两条独立证据链一致给出；
- Dubois 印刷式问题被明确复核；
- D2Q5 仅作为校验；
- 参数兼容性有可行解、无解证明或最小扩展结论；
- 报告、程序、测试和 PDF 页面均经过最新验证；
- 所有未证明命题和适用范围显式列出。
