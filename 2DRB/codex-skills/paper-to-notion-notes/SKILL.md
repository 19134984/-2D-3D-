---
name: paper-to-notion-notes
description: Turn a paper in the local detail/ folder into a detailed Chinese reading note in the user's Notion literature database. Use when the user asks for the same operation again, wants a detailed reading note, or wants a local PDF turned into a Notion literature page.
---

# Paper To Notion Notes

Use this skill when the user wants a literature note in Notion from a paper stored locally, especially in `detail/`.

## Target Notion database

- Default literature database page: `2eb69c3fe8438043a961cce26a1f2620`
- Default data source: `collection://2eb69c3f-e843-8031-b0b6-000b68182560`
- Always fetch the database first if the schema may have changed or a prior create or update call failed.

## What "same operation" means in this workspace

1. Read the paper from `detail/` using the PDF or an existing extracted `.txt`.
2. Extract verified metadata: title, authors, year, venue or journal, journal quartile / impact factor when reliably sourced, and DOI if visible.
3. Use the fixed note framework below. The framework is calibrated from the user's completed Notion notes; do not improvise a different structure just because a new paper has a different topic.
4. Build a paper evidence ledger before drafting: claims, formulas, figures/tables, benchmark cases, difficult points, and where each appears in the paper.
5. Write a detailed Chinese note, not just an abstract translation; difficult points should be unpacked until the user can connect the paper's claim to its formulas, assumptions, mechanisms, evidence, and relevance to their own LBM / thermal-convection work.
6. Run the quality gate below before creating or updating Notion.
7. Create or update the page in the Notion literature database, then set the domain tags, then fetch to verify.
8. Reply with the Notion link and a short summary of what the note emphasizes.

## Trigger phrases

- "same operation"
- "write the paper in detail/ into Notion"
- "make a detailed reading note"
- "follow my existing Notion literature notes"
- "Use $paper-to-notion-notes ..."

## Workflow

### 1) Resolve the paper

- If the user gives an exact title, find the matching PDF in `detail/`.
- If the user says they just uploaded a paper, prefer the newest PDF in `detail/`.
- If multiple files are equally plausible and the choice is risky, ask one concise question.
- Reuse an existing extracted text file in `detail/` when the filename clearly matches; otherwise run `pdftotext -layout`.
- If formula, figure, or table text is missing or mangled in extraction, inspect the relevant PDF pages directly before writing. Do not reconstruct formulas, benchmark tables, or result values from memory.

### 2) Avoid duplicates

- Search Notion for the exact paper title first.
- If a matching page exists, update that page instead of creating a duplicate unless the user explicitly asked for a new one.

### 3) Use the fixed house framework

- The main framework is fixed from the user's completed notes. Use it by default instead of only sampling 1-2 pages ad hoc.
- If calibration is needed, fetch pages in the literature database whose reading status is `已完成`, especially manually enriched notes, and use them only to confirm depth and tone, not to replace the fixed framework.
- Mirror the established style:
  - Chinese prose in the final note
  - structured sections
  - detailed interpretation rather than abstract translation
  - mechanism explanation, concrete example, code-level implication, and practical takeaway
  - no separate "one-line summary" section
- See `references/section-patterns.md` when choosing a section layout.

### 4) Build a paper evidence ledger

Before drafting the note, create an internal evidence ledger from the paper. Do not paste the ledger into Notion unless it would help readability, but use it to control quality.

Minimum ledger columns:

- `Claim / difficult point`
- `Original evidence`: section, equation number, figure, table, page, benchmark value, parameter range, or explicit author discussion
- `What the evidence proves`
- `Formula / data / algorithm detail to include`
- `Connection to current LBM / thermal-convection work`
- `Confidence`: high when directly supported; uncertain when extraction is incomplete

Rules:

- Every nontrivial claim in the note must trace back to a ledger row.
- If a difficult point has no original-paper evidence, skip it or mark it uncertain. Do not fill gaps with generic topic knowledge.
- For survey or comparison papers, organize evidence by method/theme, consensus, disagreement, and gap instead of only listing papers one by one.
- For papers with user annotations, completed Notion notes, or manually enriched notes, treat those user insights as high-priority reading signals, but still verify technical claims against the paper.

### 5) Extract the paper's real contribution

Prioritize:

- what problem the paper is actually solving
- model or method choice
- algorithm, boundary treatment, or optimization details
- what the numerical experiments really prove
- where the method is stronger or weaker
- direct relevance to the user's thermal or LBM code
- every important quantitative or qualitative result explicitly emphasized by the paper
- smaller but still meaningful results, side findings, or secondary discussion points that help explain the full picture
- ideas in the introduction that are worth learning, not just the technical core
- the historical development context when the introduction or related work provides a clear lineage
- algorithm evolution, model evolution, or why the authors moved from older methods to the current one
- new concepts, new terms, or unfamiliar method names that the paper introduces or relies on
- the paper-specific difficult points that are easy to hand-wave; identify them from the paper itself before writing, then expand those points in detail
- the causal chain behind performance or accuracy numbers, not just the numbers themselves
- how the paper's formulas, model assumptions, algorithm, implementation choices, or diagnostics would map onto the user's LBM / thermal-convection code when relevant
- benchmark cases that can be reused for comparison, including the physical configuration, boundary conditions, parameter values such as `Ra`, `Pr`, `Re`, `Le`, `Da`, grid resolution, and the reported comparison metrics

Do not omit an important result just because it appears in the introduction, discussion, or conclusion instead of the main methods section.
If the paper frames its contribution through historical comparison or method evolution, make that storyline explicit in the note instead of reducing it to a flat list of techniques.
In the results-and-discussion part, do not only summarize the biggest headline findings; also surface smaller but informative observations, secondary comparisons, boundary-case remarks, and author discussions that materially help the reader understand the paper.
Write the note so that a beginner can learn from it: when the paper mentions a new concept, new term, or specialized method name, briefly explain what it means and why it matters instead of assuming prior knowledge.
Every formula that appears in the paper must appear in the note and receive interpretation. Include displayed equations, important inline definitions, dimensionless numbers, closure relations, relaxation-time formulas, boundary formulas, diagnostic definitions, and performance metrics. Preserve equation numbers when visible. For each formula, explain the symbols, physical/numerical meaning, and why it matters for the paper or the user's code. Do not paste formulas without explanation, and do not omit formulas because they look standard.
The introduction is always worth learning from. Always include a detailed introduction-learning section that explains the motivation, historical context, method lineage, problem framing, and why the authors position the paper the way they do, even if the introduction is not a formal survey.
When the paper contains benchmark or validation cases, write them as reusable comparison targets. For each benchmark, state:
- what the case is, such as lid-driven cavity, side-heated natural convection, Rayleigh-Benard convection, channel flow, Taylor-Green vortex, thermal wave, porous cavity, or particle-laden convection
- boundary conditions and geometry
- control parameters, especially `Ra`, `Pr`, `Re`, `Ma`, `Le`, `Da`, aspect ratio, particle/porosity parameters, or other paper-specific parameters
- grid resolution, time-step or statistical window when reported
- what quantities are compared, such as `Nu`, `Re`, `Sh`, velocity extrema, streamfunction/vorticity extrema, centerline profiles, spectra, convergence order, stability limits, MLUPS/GLUPS, or parallel efficiency
- which references or tables/figures the paper compares against
- how this benchmark could be used to validate the user's current or future LBM / thermal-convection code
Before building the note, decide what the paper's own hard points are. They may be GPU/OpenACC details, but they may also be physical modeling assumptions, MRT/SRT collision structure, source terms, boundary conditions, stability analysis, benchmark design, error/convergence analysis, MPI/domain decomposition, turbulence statistics, phase-change physics, porous-media treatment, nonuniform grids, or another paper-specific issue.
Hard-point expansions must be evidence-backed from the original paper, not inferred from topic knowledge alone. For every selected complex point, identify the paper evidence first: equation number, figure, table, section, benchmark value, parameter range, algorithm step, or an explicitly stated author discussion. If the extract does not provide enough evidence, mark the point as uncertain or skip it rather than filling with generic knowledge.
For the selected complex points, use a "why -> how -> paper evidence -> implication" expansion when useful:
- why the issue matters physically, numerically, or architecturally
- how the mechanism works in plain Chinese, with formulas or a small pseudo-code / Fortran-style sketch when it clarifies the point
- what evidence the paper gives, including equation numbers, figures, tables, benchmark values, parameter scans, convergence behavior, or qualitative observations
- what this means for the user's current work when there is a real connection
Do not compress a hard section into a single confident sentence. If a reader would still wonder "why is that true?" or "how would I implement/check that?", add the missing bridge.

### 6) Build the note

Default section pattern:

- main contribution
- what problem the paper really solves
- introduction: what is worth learning from the motivation, history, and problem framing
- historical background or method-evolution context
- core model or method
- complete formula walkthrough: list and interpret every formula from the paper
- explanations of new concepts or terms that a beginner would otherwise not understand
- evidence-backed deep-dive subsections for the paper's own difficult points, chosen from the actual content rather than from a fixed topic template
- original evidence map when the paper is complex: briefly name the key equations, figures, tables, or sections that support the note's main claims
- important results and what they mean
- smaller results or discussion points worth keeping
- benchmark / validation cases and reusable comparison conditions
- experiments and key conclusions
- direct relevance to the user's current work
- limitations

Adapt sections to the paper:

- For model-heavy papers, split out a theory-analysis section.
- For boundary papers, split out the most important boundary-condition findings.
- For GPU or MPI papers, split out implementation and optimization details.
- For comparison papers, explicitly explain why one model wins and under what conditions.
- For GPU/OpenACC papers specifically, discuss CPU-vs-GPU memory behavior, data layout such as AoS/SoA, device data regions, loop fusion, gang/worker/vector choices, and whether performance is limited by bandwidth, latency hiding, PCIe transfers, or arithmetic.
- For non-GPU papers, do not force GPU/OpenACC discussion. Instead, give the same level of detail to that paper's real difficult points, such as collision modeling, forcing schemes, boundary reconstruction, stability limits, validation metrics, statistical averaging, physical assumptions, or algorithmic tradeoffs.
- For papers with formulas or models, include every formula and explain the role of each term or parameter rather than only pasting the equation.
- For papers with pseudo-code or workflow diagrams, translate the workflow into the user's likely code-level operations: which arrays are read/written, where data lives, which loop is parallelized, and which diagnostics or restart variables might be affected.
- Always add a detailed section on what is worth learning from the introduction, regardless of whether the introduction is a strong survey.
- For papers that present a historical or algorithmic lineage, add a dedicated section such as "development context", "method evolution", or "why this method emerged".
- When multiple headline results appear in the paper, mention all of them briefly even if only one or two are explored in depth.
- In the results discussion, include smaller but meaningful observations instead of compressing everything into only the top conclusions.
- When the paper uses beginner-unfriendly terminology, add short explanations in plain Chinese instead of only repeating the original jargon.

### 7) Quality gate before Notion write

Before creating or updating a Notion page, check:

- Formula gate: every displayed equation and important inline definition found in the paper appears in the note with interpretation. If extraction missed a formula, inspect the PDF page or mark the formula coverage incomplete.
- Evidence gate: every difficult-point explanation cites paper evidence internally, such as equation, figure, table, section, parameter range, benchmark value, or author discussion.
- Benchmark gate: every validation case reports reusable comparison conditions: geometry, boundary conditions, `Ra`/`Pr`/`Re`/other parameters, grid/resolution, compared quantities, and reference figure/table/source.
- Claim gate: no claim is based only on general knowledge when the paper itself did not support it.
- Relevance gate: code/work relevance is specific to the user's LBM / thermal-convection work and not generic encouragement.
- Structure gate: use the fixed framework and do not add a separate one-line summary section.
- Metadata gate: title, authors, year, venue, DOI, and journal quartile / impact factor are verified from the paper or a trusted source; unknown fields are omitted or marked uncertain rather than guessed.

If any gate fails, revise the note before writing to Notion. If the source PDF/extract is insufficient, say what could not be verified instead of producing a polished but unsupported note.

### 8) Create the page in Notion

- Confirm the current schema by fetching the database when needed.
- Use the exact current property names returned by the fetched database schema.
- The common fields in this workspace correspond to:
  - title
  - authors
  - publication year
  - venue or journal
  - journal quartile / impact factor
  - importance
  - reading status
  - domain tags
- If the database has a `期刊分区/IF` rich-text property, set it only when the journal quartile or impact factor has been reliably verified from the paper, journal site, JCR, CAS journal ranking, or another trusted source. Include the source year in the value, for example `JCR Q1 / IF 6.2 (2024 JCR)` or `中科院 1区 / IF 6.2 (2024)`.
- If the journal quartile or impact factor cannot be verified, leave `期刊分区/IF` blank. Do not guess or infer the value from the journal name, topic, or memory, because quartile and IF values are year-dependent.
- Common values:
  - For newly created literature pages, set importance to the lowest level by default (`低`).
  - For newly created literature pages, set reading status to the not-started state by default (`未开始`).
  - For newly created literature pages, set the page icon / leading emoji to a paper document emoji (`📄`).
  - Only use a different importance, reading status, or icon when the user explicitly asks for it in the current request.
- Create the page first with the basic properties and content.
- Then set the domain tags.
- Fetch the page to verify the properties and content.

### 9) Known Notion quirk in this workspace

- The domain tag field is a multi-select field.
- If direct setting is awkward with the current tool behavior, use `notion_update_page` and pass a JSON-style string such as:
  - `"[\"LBM\",\"GPU\"]"`
- Verify by fetching the page and checking that the domain tag field is an actual array in the fetched properties.

## Category guidance

Choose the smallest accurate set from the current database options.

- Always prefer the exact current tag values returned by the fetched database schema.
- Common themes in this workspace include LBM, MRT, SRT, convection-diffusion, thermal convection, side-heated cases, RB convection, 3D, GPU, CPU, Navier-Stokes, and higher-moment tuning.

## Metadata rules

- Never hallucinate authors, venue, year, or DOI.
- If a metadata field cannot be verified from the paper text or a trusted local extract, leave it out or mark it as uncertain in the body rather than guessing.
- Prefer the paper's first page and header/footer text for journal and DOI.
- When the PDF text extraction is messy, cross-check multiple snippets from the extract before creating the page.

## Local command hints

Useful local commands:

- Find PDFs: `Get-ChildItem -LiteralPath .\\detail -Filter *.pdf`
- Reuse or inspect extracts: `Get-ChildItem -LiteralPath .\\detail -Filter *.txt`
- Extract text: `pdftotext -layout <input.pdf> <output.txt>`
- Inspect first page text: `Get-Content <txt> -TotalCount 120`
- Search key phrases: `rg -n "Abstract|Keywords|Introduction|Conclusions|doi|D2Q5|D2Q9|Neumann|Dirichlet|OpenACC|MPI" <txt>`

## Final response

Keep the final reply short:

- give the Notion link
- say it has been added or updated
- mention 3-5 emphasized points
- mention the chosen domain tags if useful
