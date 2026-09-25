# A10: epithelial and niche programmes against measured organoid growth

Question-specific work for [A10](../../RESEARCH_QUESTIONS.md#a10). The canonical
hypothesis stays in the register.

**Status: structured; stage 1 authorized; nothing fitted.** Read
[RATIONALE.md](RATIONALE.md) for the biology and the argument, then
[PLAN.md](PLAN.md) for the analysis structure.

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
- **Next.** Stage 1, the identity and join audit. It decides whether a model is
  possible at all, or whether results are restricted to description.

## The gate in one sentence

The deposit has 886 RNA libraries across 203 perturbation targets, but a well is not
an animal, so the biological unit has to be established from the files before any
model is fitted.

## Layout

| Path | Contents |
|---|---|
| [RATIONALE.md](RATIONALE.md) | Biological context, hypothesis, logical flow, limits, cross-links, governance |
| [PLAN.md](PLAN.md) | Four-stage structure, outcome, covariates, tests, holdout, decision rules |
| [config/a10_outcome_contract.json](config/a10_outcome_contract.json) | The same, machine-readable, with cached file hashes |
| `scripts/` | Stage 1 audit; later stages are written only after the owner reviews stage 1 |
| `tables/` | Stage 1 outputs and its run record |
| `cache/` | The three fetched metadata files; ignored, hashes tracked |

## Two things a later session must not do

- **Do not write a study note on the source paper.** Its reading by the owner is
  recorded as not started, and an earlier session's note for an unread paper was
  rejected (DEVELOPMENT decision 21). This work uses the deposit, not the paper.
- **Do not split wells at random for evaluation.** Sibling wells share isolation,
  fibroblast lot and plate handling. Holdouts are whole preparations.
