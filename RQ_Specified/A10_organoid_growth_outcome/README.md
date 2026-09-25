# A10: epithelial and niche programmes against measured organoid growth

Question-specific work for [A10](../../RESEARCH_QUESTIONS.md#a10). The canonical
hypothesis stays in the register.

**Status: stages 1 and 3 complete; the analysis has run.** Read
[RATIONALE.md](RATIONALE.md) for the biology and the argument,
[PLAN.md](PLAN.md) for the analysis structure, and
[the audit](reports/STAGE1_IDENTITY_AUDIT.md) for what the design can and cannot
support.

**Result: both increments are inconclusive.** The epithelial increment fails the
declared margin, and the fibroblast increment clears it only in the primary fit, then
falls below it once the one over-replicated target is removed. The baseline itself does
not transfer to a held-out unit, with a median within-unit predictive R squared near
zero against a pooled 0.586. Read
[the stage 3 report](reports/STAGE3_FIT_REPORT.md) first.

The assay validation passes convincingly, so this is a negative result about the
hypothesis rather than a broken pipeline. The biological unit also remains unresolved,
which would have capped any positive result at within-screen description.

## Why this question is worth the effort

Every other question in the register compares RNA states with other RNA states. A10
is the first whose answer can be set against something the experiment measured about
tissue building: how much organoid a well actually produced.

The assay also captures the biology that makes alveolar repair hard to read from
epithelium alone. Type 2 stemness is maintained by fibroblast signals, so a gene
knocked out in the epithelium can change repair either by changing that cell or by
changing what its niche does in response. This deposit cultures mouse type 2 cells
with human fibroblasts and assigns reads by species, so one well gives both sides of
that exchange, with the perturbation made only in the epithelium.

## Where the work stands

- **Cached.** Three small metadata tables: the perturbation design, the imaging
  outcomes and the species-assignment quality control. Hashes are in the contract.
- **Deliberately not fetched.** The 303 MB count table. It is stage 3 input, and
  fetching it before the design gate would prejudge the gate.
- **Done.** Stage 1, the identity and join audit. A public data search, which found no
  cohort that could validate this screen externally.
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
| [reports/STAGE3_FIT_REPORT.md](reports/STAGE3_FIT_REPORT.md) | The fitted result, its sensitivity checks and the per-unit failure |
| [reports/PUBLIC_DATA_SEARCH.md](reports/PUBLIC_DATA_SEARCH.md) | Whether any public cohort could validate A10, and why none can |
| `scripts/01_audit_identities.py` | Stage 1: joins and units, metadata only |
| `scripts/02_extract_counts.py` | Stage 3a: stream the two count sheets, cache the declared gene set |
| `scripts/03_fit_outcome_models.py` | Stage 3b: nested models, leave-one-unit-out, plus assay validation |
| `tables/` | Stage outputs and run records |
| `cache/` | The three fetched metadata files; ignored, hashes tracked |

## Two things a later session must not do

- **Do not write a study note on the source paper.** Its reading by the owner is
  recorded as not started, and an earlier session's note for an unread paper was
  rejected (DEVELOPMENT decision 21). This work uses the deposit, not the paper.
- **Do not split wells at random for evaluation.** Sibling wells share isolation,
  fibroblast lot and plate handling. Holdouts are whole preparations.
