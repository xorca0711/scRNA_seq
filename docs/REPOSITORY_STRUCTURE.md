# Repository structure and label scope

Updated 26 September 2026 (A0 registration). This contract follows the existing shared-analysis
and paper-study layout. It defines where current material belongs; dated
protocols, original trial names and immutable run records remain historical evidence.

| Material | Canonical location | Rule |
|---|---|---|
| Repository-wide scientific questions | [RESEARCH_QUESTIONS.md](../RESEARCH_QUESTIONS.md) | One current register, A0–A14; link to source studies rather than maintaining a second register |
| Question-specific plans and analyses | `RQ_Specified/A<id>_<topic>/` | Prospective plan, source metadata, configuration, scripts, tables, reports and gallery; reference the root question register |
| Shared question figures | [analysis/figures/rq/README.md](../analysis/figures/rq/README.md) | Curated gallery/captions; assets, tables and provenance alongside; fresh script-16 runs under `renders/` |
| Cross-question measurement checks | [RQ_MEASUREMENT_CONTRACTS.md](RQ_MEASUREMENT_CONTRACTS.md) | Reusable checks and legacy-ID crosswalk; biological decisions stay in the root register |
| Shared plotting and pipeline entrypoints | `analysis/scripts/` | Numbered entrypoints; reusable utilities in `analysis/lib/` |
| Palette, reusable settings and contracts | `analysis/config/` | Shared configuration, independent of a paper's local trial ID |
| Paper notes, contracts, inference and galleries | `Research Article/<paper>/` | Preserve the paper-specific source evidence, scripts, tables and figure gallery |
| Cross-study specificity analysis | `Research Article/epithelial_state_specificity/` | Existing substantive analysis module, referenced by the global questions |
| Corrective scientific analyses | `analysis/corrections/` | Keep original results and corrections distinguishable |
| Current scientific findings and evidence ledger | [FINDINGS.md](../FINDINGS.md), [CLAIMS.md](../CLAIMS.md), [NEGATIVE_RESULTS.md](../NEGATIVE_RESULTS.md) | Claims change only with supporting evidence; a layout migration does not promote them |
| Methods, structure, portfolio and communication drafts | `docs/` | Each page identifies whether it describes execution, reference methods or proposed work |
| Reading order and paper status | [Research Article/ROADMAP.json](../Research%20Article/ROADMAP.json), [Research Article/README.md](../Research%20Article/README.md) | Stable paper identifiers differ from the reading sequence; update the two views together |
| Current handoff and operating context | [PROGRESS.md](../PROGRESS.md), [AI_CONTEXT.md](../AI_CONTEXT.md) | Current section first; older dated checkpoints remain historical |
| Raw inputs, caches and private reading annotations | Ignored local directories | Do not copy them into public figure or documentation directories |

## Identifier namespaces

- `A0`–`A14` identify repository-wide questions; figures use the associated
  question ID, with panel/group suffixes where needed.
- Paper-local IDs require paper context: `Niethamer/W1`, `Niethamer/S1`,
  `Sikkema/S1`, `Choi/D1`, `Yu/N1`, `Yu/U5`, `Yu/F01`. Equal short labels do
  not mean equal analyses. Existing historical scripts are not renamed solely
  to make all short labels globally unique.
- `W1` is the historical Wagner-branch myeloid pseudobulk trial label. It is
  neither a statistical evidence grade nor evidence that Compass flux modelling ran.
- The former Yu follow-ups `RQ1`–`RQ4` are now `A11`–`A14`.
  Display labels D1, D2a–c, D0 and D4 map to A11, A12a–c, A12-S1 and A14.
  Historical D0/D1/D2 CSV basenames are retained and documented by the
  [relocation manifest](../analysis/figures/rq/il1b_context/relocation_manifest.json).
- `C1` onward identifies claims. Reading-order paper numbers, trial IDs,
  global questions, figure groups and claim IDs are different namespaces.

## Paths and historical provenance

On 25 September 2026, `Thesis/` was renamed to `Research Article/` and
`RQ_Specified/` was introduced at the owner's request. Paper folder identifiers,
scientific results and historical run identities are unchanged. The
[migration record](migrations/2026-09-25-research-layout/README.md) lists preserved
hashes and archives original bytes of updated code and documentation. Historical
paths are resolved by `analysis/lib/repository_paths.py`; this is a record check,
not evidence that the updated scripts have rerun. Quote paths containing
`Research Article` in commands; Markdown links encode the space as `%20`.

`ROADMAP.json` paper `folder` values are relative to `Research Article/`, without a
trailing slash; other artifact fields are repository-relative. The owner-selected
`gate2_C3_yu_lee_choi_min_2026/` is branch 2C item 3, stable paper 13.
This documented exception does not renumber the roadmap.

Shared IL-1 figures now use scripts 19–21. Their recorded UMAP/PCA preparation
was executed before relocation; the original script bytes and identities
remain archived alongside the [current figure report](../analysis/figures/rq/il1b_context/REPORT.md).
Original preparation input paths are study-relative; current rendering input
and output paths are repository-relative. No execution date or original hash
is rewritten to make an old run appear new. Future render records archive
their predecessor; visual review hashes must match the current PNGs.

Current summaries take status from stage reports and validation artifacts.
Older plans and first-batch reports do not override a completed continuation.
“Feasible analyses complete” does not mean that unavailable biological
identities, independent validation, spatial regions or causal endpoints are resolved.

The main README links to galleries without embedding a partial selection.
PI fit, contact preferences and personal outreach planning stay in the owner's
Notion workspace; public documents contain scientific questions and portfolio text.
