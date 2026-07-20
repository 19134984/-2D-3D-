# 参考文献与来源说明

1. LBM-CDE 原始论文，用户提供文件 `LBM-CDE.pdf`。本文使用 Eqs. (1)--(35)、
   Eq. (36)--(39) 及 Appendix A 的 transformed trapezoidal LBE；逐页符号映射见
   `docs/derivation/evidence-ledger.md`。

2. Wang 等，*Lattice Boltzmann simulations of thermal convective flows in two
   dimensions*。本文只把其 D2Q5 参数化温度矩与边界 kinetic rules 用作限定来源，
   不把 D2Q5 四阶条件外推到 D2Q9。

3. Contrino 等 / Luo，*Lattice Boltzmann simulations of the thermally driven
   2D square cavity at high Rayleigh numbers*，Journal of Computational Physics
   (2014)。本文审计其 D2Q5 四阶参数与 pressure-boundary Poiseuille 的 $3/8$
   标定，并与 uniform-body-force 的 $3/16$ 分开记录。

4. Dubois 与 Lallemand，*Towards higher order lattice Boltzmann schemes*。
   本文复算其 D2Q5 Eqs. (39)--(42)、(55) 与 Appendix Eq. (79)，并逐字审计
   第 12 页 D2Q9 印刷式；这些公式不作为 LBM-CDE--D2Q9 待求系数的输入。

5. Ginzburg 与 d'Humières，*Multireflection boundary conditions for lattice
   Boltzmann models*，Physical Review E 68, 066614 (2003)。本文使用其
   Eqs. (41)--(43) 建立 $\Lambda^2=1/4$ 与受限 Hénon product $3/16$ 的映射。

所有页码、方程号、允许主张与限制均在证据账本逐条登记。项目内 exact algebra
的入口与测试名见需求—证据—验证矩阵；`Xs/` 不在理论来源列表中。
