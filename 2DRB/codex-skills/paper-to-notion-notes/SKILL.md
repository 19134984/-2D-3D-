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
2. Extract verified metadata: title, authors, year, venue or journal, and DOI if visible.
3. Read 1-2 existing Notion literature notes to match the user's style and depth.
4. Write a detailed Chinese note, not just an abstract translation.
5. Create or update the page in the Notion literature database, then set the domain tags, then fetch to verify.
6. Reply with the Notion link and a short summary of what the note emphasizes.

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

### 2) Avoid duplicates

- Search Notion for the exact paper title first.
- If a matching page exists, update that page instead of creating a duplicate unless the user explicitly asked for a new one.

### 3) Match the house style

- Fetch 1-2 existing pages from the literature database in the same topic area when possible.
- Mirror the established style:
  - Chinese prose in the final note
  - structured sections
  - detailed interpretation rather than abstract translation
  - explicit "direct relevance to my current work" and "one-line summary" style endings when appropriate
- See `references/section-patterns.md` when choosing a section layout.

### 4) Extract the paper's real contribution

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

Do not omit an important result just because it appears in the introduction, discussion, or conclusion instead of the main methods section.
If the paper frames its contribution through historical comparison or method evolution, make that storyline explicit in the note instead of reducing it to a flat list of techniques.
In the results-and-discussion part, do not only summarize the biggest headline findings; also surface smaller but informative observations, secondary comparisons, boundary-case remarks, and author discussions that materially help the reader understand the paper.
Write the note so that a beginner can learn from it: when the paper mentions a new concept, new term, or specialized method name, briefly explain what it means and why it matters instead of assuming prior knowledge.

### 5) Build the note

Default section pattern:

- main contribution
- what problem the paper really solves
- historical background or method-evolution context when the paper provides it
- core model or method
- explanations of new concepts or terms that a beginner would otherwise not understand
- important results and what they mean
- smaller results or discussion points worth keeping
- experiments and key conclusions
- direct relevance to the user's current work
- limitations
- one-line summary

Adapt sections to the paper:

- For model-heavy papers, split out a theory-analysis section.
- For boundary papers, split out the most important boundary-condition findings.
- For GPU or MPI papers, split out implementation and optimization details.
- For comparison papers, explicitly explain why one model wins and under what conditions.
- For papers with a strong introduction survey, add a dedicated section on what is worth learning from the introduction.
- For papers that present a historical or algorithmic lineage, add a dedicated section such as "development context", "method evolution", or "why this method emerged".
- When multiple headline results appear in the paper, mention all of them briefly even if only one or two are explored in depth.
- In the results discussion, include smaller but meaningful observations instead of compressing everything into only the top conclusions.
- When the paper uses beginner-unfriendly terminology, add short explanations in plain Chinese instead of only repeating the original jargon.

### 6) Create the page in Notion

- Confirm the current schema by fetching the database when needed.
- Use the exact current property names returned by the fetched database schema.
- The common fields in this workspace correspond to:
  - title
  - authors
  - publication year
  - venue or journal
  - importance
  - reading status
  - domain tags
- Common values:
  - importance: usually the highest level for papers the user explicitly asked to read in detail
  - reading status: usually the completed state
- Create the page first with the basic properties and content.
- Then set the domain tags.
- Fetch the page to verify the properties and content.

### 7) Known Notion quirk in this workspace

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
