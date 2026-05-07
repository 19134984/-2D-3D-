# Section Patterns

Use these patterns only as selection guidance. Do not copy all sections blindly.

## Baseline pattern

- `# 主要工作`
- `# 这篇文章真正解决了什么问题`
- `# 核心模型/方法`
- `# 数值实验与关键结论`
- `# 对我当前工作的直接启发`
- `# 局限`
- `# 一句话总结`

## Model comparison paper

Use when the paper compares schemes such as `D2Q5 vs D2Q9`, `SRT vs MRT`, or different thermal models.

- `# 主要工作`
- `# 这篇文章真正解决了什么问题`
- `# 文中比较的模型到底是什么`
- `# 理论分析主线`
- `# 数值实验与作者最后的判断`
- `# 对我当前工作的直接启发`
- `# 一句话总结`

## Boundary-condition paper

Use when the paper focuses on `Dirichlet`, `Neumann`, curved boundaries, interface treatment, or flux recovery.

- `# 主要工作`
- `# 这篇文章真正解决了什么问题`
- `# 边界构造的核心思路`
- `# 边界条件部分最关键的发现`
- `# 数值验证到底验证了什么`
- `# 对我当前工作的直接启发`
- `# 局限`
- `# 一句话总结`

## GPU / MPI / performance paper

Use when the paper focuses on acceleration, scalability, OpenACC, CUDA, MPI, overlap, or memory layout.

- `# 主要工作`
- `# 这篇文章真正解决了什么问题`
- `# 算法与并行对象`
- `# 并行实现与优化细节`
- `# 性能结果到底说明了什么`
- `# 对我当前工作的直接启发`
- `# 一句话总结`

## Writing rules for this workspace

- Do not stop at abstract translation.
- Always explain why the paper matters for the user's current thermal or LBM work.
- Prefer explicit judgments such as “真正困难的部分是...” or “这篇文章最值得抓住的一点是...”.
- When a result only holds under specific conditions, state those conditions clearly.
- If the paper compares methods, say when each one wins instead of only saying one is “better”.
