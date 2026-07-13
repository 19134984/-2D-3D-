# 07 TRT 参数兼容性与不可兼容性

本章把第五章的受限壁面残差与第六章的冻结系数四阶体相条件放进同一个精确求解器。目标不是选一个文献参数点，而是回答：给定输运系数、名义 TRT 率、物理有效块、体相四阶条件和壁面标定后，这组约束是否仍有解。

## 1. 量、约束与报告接口

统一记

$$
\sigma=\frac1s-\frac12,
\qquad
a=c_s^2+\pi,
\qquad
b=(1-\chi_\kappa)c_s^2,
\qquad
K=\frac\kappa{\Delta t}.
$$

标量名义奇 ghost、名义偶 ghost 和反馈后的物理通量平移分别为

$$
\sigma_o=\frac Kb,
\qquad
\sigma_e,
\qquad
\sigma_{\rm flux}=\frac ba\sigma_o=\frac Ka.
$$

这里的 $\sigma_{\rm flux}$ 不能替换整个奇子空间的 $\sigma_o$。碰撞输入仍是名义奇偶率；物理通量率只是局部反馈消元后该物理块的等效系数。

`ParameterReport` 不返回裸率元组，而是同时给出：

- `status`：`feasible_exact`、`feasible_restricted`、`no_feasible_solution`、`degenerate_branch`、`boundary_correction_required` 或 `mrt_extension_required`；
- `exact_substitutions`：精确有理式、根式和极限；
- `collision_rates` 与逐块 `open_interval_checks`；
- `violated_constraints`、全部推导假设和 `minimal_extension`；
- `consumed_evidence`：Task 5 的规范 `D2Q9EquivalentCoefficients` / `QuarticConditionSystem` 对象和 Task 6 的实际 residual/classification 对象，不用字符串冒充推导证据；
- `is_conditional` 与 `feasibility_conditions`：自由符号输入尚不能判定正性或开区间时，显式保存待满足条件。

Task 5 的 `canonical_quartic_condition()` 从受审闭式输入构造保留因子的四阶多项式，并在多个精确参数点同时与 amplification route 和 Taylor route 核对。参数层直接消费该规范对象；包括 $a=1$ 在内都从同一个未除以 $a-1$ 的多项式分类，不再复制一套 reviewed family，也不把某个有理专门化点冒充一般证据。

## 2. Task 5 体相条件

`recover_baseline_quartic_family()` 对 `amplification_route(order=4)` 的输出调用 `quartic_condition_system()`，恢复完整基线族

$$
\sigma_e=\frac{4\sigma_o^2+1}{8\sigma_o},
$$

而不是直接返回某个熟知 TRT 点。

对外置精确梯度源和局部非平衡反馈，Task 5 的受审主分支分别为

$$
\sigma_e^{\rm ext}
=\frac{3ab-a-b-12b^2\sigma_o^2}
{12b\sigma_o(a-1)},
$$

$$
\sigma_e^{\rm fb}
=\frac{3a^2-2a-12b^2\sigma_o^2}
{12b\sigma_o(a-1)}.
$$

它们只适用于 $a\ne1$、$b\ne0$ 和非零主平移。求解器先分类 $a=0$、$b=0$ 与 $a=1$，不会把含 $a-1$ 的表达式称为普适公式。

## 3. 与受限标量 ABB 的精确消元

Task 6 在以下假设内给出唯一可用的标量 ABB 行：D2Q9 两率标量碰撞、平直 grid-aligned halfway wall、稳态一维二次温度场、半热源重构、完整压力平衡墙项，以及对应的 external-gradient 或 local-feedback population chain。在这个适用域内，

$$
\sigma_{\rm flux}\sigma_e=\frac3{16}.
$$

代入 $\sigma_{\rm flux}=K/a$ 得

$$
\sigma_e=\frac{3a}{16K}.
$$

### 3.1 局部反馈

把 $\sigma_o=K/b$ 与上式代入 Task 5 反馈主分支，精确约去分母后得到

$$
\boxed{K^2=\frac{a(3a+1)}{48}}.
$$

### 3.2 外置精确梯度源

相同消元给出

$$
\boxed{
K^2=\frac{12ab+5a-4b-9a^2}{48}
}.
$$

这两式是“冻结体相四阶完全消除 + 受限 ABB 二次行 + 给定输运”同时成立的条件，不是推荐扩散率。

当 $\pi=0$、$a=c_s^2=1/3$ 时，两式都退化为正物理分支

$$
\frac\kappa{\Delta t}=\frac1{\sqrt{72}}.
$$

相应物理通量与名义偶平移固定为

$$
\sigma_{\rm flux}=\frac{\sqrt2}{4},
\qquad
\sigma_e=\frac{3\sqrt2}{8},
$$

但名义奇平移仍为

$$
\sigma_o=\frac1{2\sqrt2(1-\chi_\kappa)}.
$$

因此必须给出 $\chi_\kappa$ 才能生成名义奇率。正物理分支还要求 $a\ne0$、$b>0$、$K>0$、三个相关 Hénon 平移为正，并逐一检查名义奇率、名义偶率和物理通量率是否属于开区间 $(0,2)$。该区间只表示形式可接受，不表示条件数良好或数值稳健。

## 4. 低扩散率与诚实的约束选择

若 $\Delta t=1$ 且目标 $\kappa=10^{-3}$，一般不能同时满足上述体相四阶条件与受限 ABB 条件。求解器返回 `no_feasible_solution`，并把 `bulk_quartic_and_restricted_abb` 列入违反项；它不会静默保留其中一个条件。

若明确选择“保留体相四阶消除，墙面使用显式修正”，低扩散率仍可有可接受的名义率。以反馈分支

$$
a=\frac13,
\qquad
b=\frac14,
\qquad
\Delta t=1,
\qquad
\kappa=\frac1{1000}
$$

为例，Task 5 体相条件给出

$$
\sigma_o=\frac1{250},
\qquad
\sigma_e=\frac{250009}{6000},
\qquad
\sigma_{\rm flux}=\frac3{1000},
$$

$$
s_o=\frac{125}{63},
\qquad
s_e=\frac{6000}{253009},
\qquad
s_{\rm flux}=\frac{1000}{503}.
$$

三个率都在 $(0,2)$ 内，但状态仍是 `boundary_correction_required`。这里的显式修正只释放已经审计过的受限 ABB product，不能自动消除一般 pressure/source/time/force/tangential jets。

工程上有三种彼此独立、不能偷换的选择：

1. 低 $\kappa$ + 体相四阶消除 + 显式墙面修正；
2. 低 $\kappa$ + 受限 ABB 标定 + 保留体相四阶误差；
3. 当两者都必须满足时，采用已验证的显式受限 ABB 边界修正；split-even MRT 只能列作有待另行推导的候选。

显式边界修正与 split-even MRT 不是同一种机制。前者是当前用于释放受限 ABB product 的充分扩展；后者只有在矩空间推导明确“哪个偶模态进入 ABB、哪个独立偶模态进入 $C_{40}$”，并证明双约束 Jacobian 满秩后，才能升级为可行结论。当前 split-even 只报告 `candidate_requiring_mode_jacobian_derivation`，不能与显式修正并列为已证明充分方案，也不能标成 `feasible_exact`。

## 5. 特殊与退化分支

求解顺序固定为先分类、后除法：

| 分支 | 精确结论 | 报告 |
| --- | --- | --- |
| 数值 $a\le0$ | 平衡/物理块非正 | `no_feasible_solution`，列出 `a_not_positive` |
| 数值 $b\le0$ | 输运块非正 | `no_feasible_solution`，列出 `b_not_positive` |
| 数值 $K\le0$ | 目标扩散输运非正 | `no_feasible_solution`，列出相应非正约束 |
| $a=1$，反馈 | 直接条件 $12K^2=1$ | 不使用主分支除以 $a-1$ |
| $a=1$，外置梯度 | $1-2b+12K^2=0$，正支要求 $b>1/2$ | 不使用主分支除以 $a-1$ |
| 兼容根号为零 | 零 Hénon 平移不作为可行 ABB 极限 | `no_feasible_solution` |
| 兼容根号为负 | 无实正物理分支 | `no_feasible_solution` |
| 反馈 $a+2K=0$ 且其余输入为正 | 局部梯度闭合奇异 | `degenerate_branch` |

数值可行性优先于代数边界分类：形式交点 $a=1/9$、$K=-1/18$ 虽满足 $a+2K=0$，但因 $K<0$ 返回 `no_feasible_solution`。任一实际碰撞率不在 $(0,2)$ 也采用相同优先级。只有自由符号无法判定时才保留条件报告，而不是宣称已经可行。

## 6. 流动侧兼容性

剪切输运给出

$$
\nu=c_s^2\Delta t(1-\chi_s)\sigma_f^+.
$$

只保留 Task 6 的稳态一维均匀体力剪切行时，

$$
(1-\chi_s)\sigma_f^+\sigma_f^-=\frac3{16}.
$$

对给定正 $\nu$ 精确解得

$$
\sigma_f^+=\frac{\nu}{c_s^2\Delta t(1-\chi_s)},
\qquad
\sigma_f^-=\frac{3c_s^2\Delta t}{16\nu},
$$

$$
s_f^-=\frac{16\nu}{3c_s^2\Delta t+8\nu}
\longrightarrow0^+
\quad(\nu\to0^+).
$$

所以任意正 $\nu$ 下的开区间检查不等于小黏性时具有良好阻尼或条件数。

若保留 trace jets，求解器改为消费 Task 6 的 general residual/classification。$\chi_b\ne\chi_s$ 时，已解析 shear/bulk 两行的 `rate_compatibility_status` 为 `no_single_magic`，返回 `mrt_extension_required`；若选择独立剪切/体积率，还必须同时加入一般速度边界修正，不能只解剪切/体积两行后宣称通用壁面标定。另一条独立充分路线是直接采用一般速度边界修正。即使 $\chi_b=\chi_s$，完整一般壁面表仍保留未闭合 jets，因此顶层状态是 `boundary_correction_required`，而不是 universal magic。

## 7. 结论边界

- $1/\sqrt{72}$ 只是特定约束集合在 $\pi=0$ 的兼容点，不是通用推荐值。
- $3/16$ 只在 Task 6 明示的受限墙面、source/gauge、几何和场假设内使用。
- 冻结系数四阶条件不外推到变系数四阶精度。
- 名义 ghost 率、物理有效率、输运系数和边界校准量始终分栏报告。
- 任何不可兼容的约束集合都返回非空 `minimal_extension`，不静默删约束。
