# A10: epithelial and niche programmes against measured organoid growth

Question-specific work for [A10](../../RESEARCH_QUESTIONS.md#a10). The canonical
hypothesis stays in the register.

**Status: identity audit, original fit and revised specification complete; PRs #76/#78 merged.** Read
[RATIONALE.md](RATIONALE.md) for the biology and the argument,
[PLAN.md](PLAN.md) for the analysis structure, and
[the audit](reports/STAGE1_IDENTITY_AUDIT.md) for what the design can and cannot
support.

**Result after the revised specification: the epithelial growth programmes clear the
declared margin for within-group association; the fibroblast block does not.** Read
[the stage 4 report](reports/STAGE4_REVISED_REPORT.md) first, then
[stage 3](reports/STAGE3_FIT_REPORT.md) for how the first specification failed and why.

- **The supported cell** pairs the within-unit baseline with growth-nominated Hallmark
  programmes: an increment of 0.0426 against a 0.02 margin, clearing it in all four
  mandatory sensitivity checks.
- **Both changes were needed.** Growth programmes with the old baseline give 0.0136, and
  the old repair modules with the new baseline give 0.0036. Neither alone suffices.
- **The earlier fibroblast increment is specification-dependent.** It turns negative
  in the revised growth model. This does not establish absence of a niche role.
- **It is uneven.** The increment is positive in 8 of 15 units, and plates 2 and 4 still
  have baselines that do not predict. The pooled result is carried by two of four plates.

Target-transcript reduction supplies a limited consistency check. Biological units
remain unresolved, the revised specification follows inspection of the original fit,
and the within-unit metric uses held-out outcome means. Everything here remains
a within-screen descriptive association. See the
[logical review](../../docs/LOGICAL_RATIONALE_REVIEW.md#a10-corrections) for the
scope of the checks and differences between the written plan and implementation.

## Why this question is worth the effort

A10 relates compartment-level RNA to a measured imaging endpoint: day-14 mean
organoid area conditional on day-7 mean area. It complements A1's existing lineage
and microscopy outcomes. Mean area measures organoid size and morphology; it does
not by itself measure total tissue production, mature AT1 function or lung repair.

Fibroblast signals can affect type 2 stemness, motivating measurement of both
compartments. This deposit cultures mouse type 2 cells with human fibroblasts and
assigns reads by species. The perturbation is epithelial, but both endpoint RNA
profiles may reflect feedback, cell abundance and shared culture conditions.
Species assignment distinguishes RNA origin, not cell-autonomous versus niche-mediated
causal effects. The [rationale](RATIONALE.md) explains the biological premise.

## Where the work stands

- **Retrieved and processed.** Metadata, the count workbook and frozen gene sets
  were used for the completed fits. Local input caches are ignored; hashes and
  output records are tracked in `tables/`. A clean clone does not contain raw inputs.
- **Done.** Identity and join audit, count extraction, original nested fits,
  revised four-cell grid and required sensitivities. The recorded public-data
  search found no eligible external validation cohort.
- **Resolved since stage 1.** The count workbook names its two sheets by species, so
  the epithelial and fibroblast compartments no longer have to be inferred.
- **Still unresolved.** The biological unit. With no external cohort available, that is
  now the binding constraint on how strongly any result can be read.

## The gate in one sentence

The deposit has 886 RNA libraries across 203 perturbation targets, but a well is not
an animal, and the deposit does not say how many independent preparations there are.

## Layout

| Path | Contents |
|---|---|
| [RATIONALE.md](RATIONALE.md) | Biological context, hypothesis, logical flow, limits, cross-links, governance |
| [PLAN.md](PLAN.md) | Four-stage structure, outcome, covariates, tests, holdout, decision rules |
| [config/a10_outcome_contract.json](config/a10_outcome_contract.json) | The same, machine-readable, with cached file hashes |
| [reports/STAGE1_IDENTITY_AUDIT.md](reports/STAGE1_IDENTITY_AUDIT.md) | What joined, why the unit is unresolved, and what stage 2 must assume |
| [reports/STAGE4_REVISED_REPORT.md](reports/STAGE4_REVISED_REPORT.md) | The four-cell grid, the supported cell, and the limitation that qualifies it |
| [reports/STAGE3_FIT_REPORT.md](reports/STAGE3_FIT_REPORT.md) | The first specification, its inconclusive result and the per-unit failure that pointed at the fix |
| [reports/PUBLIC_DATA_SEARCH.md](reports/PUBLIC_DATA_SEARCH.md) | Whether any public cohort could validate A10, and why none can |
| `scripts/01_audit_identities.py` | Stage 1: joins and units, metadata only |
| `scripts/02_extract_counts.py` | Stage 3a: stream the two count sheets, cache the declared gene set |
| `scripts/03_fit_outcome_models.py` | Stage 3b: nested models, leave-one-unit-out, plus assay validation |
| `tables/` | Stage outputs and run records |
| `cache/` | Local metadata/count intermediates; ignored, hashes tracked |

## Two things a later session must not do

- **Do not write a study note on the source paper.** Its reading by the owner is
  recorded as not started, and an earlier session's note for an unread paper was
  rejected (DEVELOPMENT decision 21). This work uses the deposit, not the paper.
- **Do not split wells at random for evaluation.** Wells may share preparation and
  handling. Completed holdouts are entire deposited plate-replicate groups;
  whether those groups are independent preparations is unresolved.
