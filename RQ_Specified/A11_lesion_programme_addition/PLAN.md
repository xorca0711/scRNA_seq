# A11 pre-registered test in the Kim 2020 cohort

25 September 2026. **Pre-registered; not scored.** The owner chose this design at
stage 4 of the [shared component contract](../A5_A11_shared_component_contract/README.md)
(DEVELOPMENT decision 37). The machine-readable version is
[`config/kim2020_test_contract.json`](config/kim2020_test_contract.json). Scoring
waits for the owner to retain this plan.

## The question

[A11](../../RESEARCH_QUESTIONS.md#a11) asks whether lesion-associated epithelial
programmes add to a shared plasticity component. The test here is narrower. Does
the frozen lesion-specific module rise in lung adenocarcinoma epithelium, relative
to the same patient's normal type 2 cells, in a cohort that has never been scored
with it?

## Why this cohort, and why not the obvious one

The lesion-specific module from the contract is identical, gene for gene, to a
module the completed human run already scored in 23 patients. So that cohort can
only be the discovery. Scoring it again would reproduce a reported number.

| Discovery, already reported | Value |
|---|--:|
| Patients | 23 |
| Mean paired difference, log2 CPM | 0.323 |
| Standard deviation of differences | 0.385 |
| Patients positive | 19 |

The evaluation cohort is Kim et al. 2020
([10.1038/s41467-020-16164-1](https://doi.org/10.1038/s41467-020-16164-1)), held
locally. Ten patients have both a tumour sample and a normal-lung sample. The
pairing comes from the deposit's own patient identifier column, not from matching
sample names. The cohort was opened once before, for EGFR ligand RNA, and never
scored with any module in this plan.

## Populations

- **Normal arm.** Cells the authors labelled type 2 in the patient's normal-lung
  sample.
- **Lesion arm.** All tumour epithelial cells in the patient's tumour sample, author
  subtypes tS1, tS2 and tS3 pooled. No subtype is picked, to avoid a post hoc choice.

This differs from the discovery, and the difference matters. The discovery used
atlas-mapped type 2-like cells in both arms, which in lesions included lesional
cells that resemble type 2. Here the lesion arm is all tumour epithelium. So this
is a stricter transport test: a positive result survives a change of cell
definition, and a negative one may reflect that change. The author annotation has
no copy-number malignancy label, so the lesion arm may also include some
non-malignant epithelium.

If the owner prefers the discovery's own cell definition, the alternative is to
rerun its atlas mapping on Kim before any scoring. That costs a reference-mapping
run and changes nothing else in this plan.

## Eligibility

The unit is the patient. A patient is eligible with at least 50 cells in each arm.
That floor is inherited from the discovery's primary configuration, not chosen for
this cohort. At least three eligible patients are required. The metadata show
eight eligible patients, before any score exists; two fall short on tumour
epithelium.

Each module must also reach an assayed source fraction of 0.7 in Kim's gene index,
computed as in the contract. A module that fails is not scored. If the primary
fails, the test does not run.

## Scoring uses the discovery's own instrument

The discovery scores come from a tracked R script,
`u5_human_paired_pathways.R`. This test uses the identical procedure:

1. Sum raw counts over each arm's cells for each patient.
2. Normalize all eligible pseudobulks together with edgeR TMM.
3. Take log2 CPM with a prior count of 1.
4. Score a module as the mean log2 CPM of its genes present.
5. Take the lesion arm minus the normal arm for each patient.

Before any Kim score, the same code must reproduce the discovery's tracked
per-patient differences for this module from the cached discovery pseudobulks,
within 1e-6. If it cannot, Kim is not scored and the discrepancy is reported.

## Primary test

| Item | Choice |
|---|---|
| Module | lesion-specific, 91 mouse genes |
| Estimand | mean within-patient difference, lesion minus normal type 2 |
| Test | exact two-sided Wilcoxon signed-rank, alpha 0.05 |
| Estimate | Hodges-Lehmann with exact 95 percent interval |
| Sensitivity | paired t 95 percent interval |
| Smallest effect of interest | 0.10 log2 CPM, about a third of the discovery mean |

Decision rules, fixed now:

- **Transports** if p is below 0.05 and the estimate is above zero.
- **Contradicts** if p is below 0.05 and the estimate is below zero.
- **Precisely absent** if the upper bound of the exact interval is below 0.10. This
  retires the specified addition only, not every lesion programme.
- **Inconclusive** otherwise.

## Secondary tests

One family, Benjamini-Hochberg at q below 0.05, same test as the primary:

1. **Beyond shared.** Each patient's lesion-specific difference minus their
   shared-remodelling difference. This is the literal form of "adds to a shared
   component".
2. **Stress excluded.** The lesion-specific module with every gene of the three
   Hallmark control sets removed. Nineteen of its 91 genes sit in those sets,
   leaving 72, so a rise driven only by stress would fail here.
3. **Injury-lesion pair.** The seven-gene component shared by injury and lesions,
   five of them human-mappable.

Reported without a test: the shared module; type 1 and type 2 identity axes; the
three Hallmark axes; the leave-one-patient-out range of the primary; and the primary
without the one double-primary patient.

## Power, stated before any score

From a logged calculation using only the discovery's tracked differences:

| Scenario | Exact Wilcoxon | Paired t |
|---|--:|--:|
| Discovery effect, 8 patients | 0.46 | 0.53 |
| Half the discovery effect | 0.15 | 0.18 |

So even if the discovery effect transports in full, this test is close to a coin
flip. First estimates are often inflated, and at half the effect the chance is
about one in six. A null result will therefore be inconclusive unless the
precise-absence rule is met, which is unlikely with eight patients. A positive
result would be meaningful precisely because the test is hard to pass.

## What this cannot show

- A null is not a refutation unless the precise-absence rule is met.
- Tumour epithelium is not confirmed malignant.
- Same-patient pairing controls the donor, not stage, smoking or driver mutation,
  which all vary across these patients.
- A positive result supports transport of this one addition. It does not show that
  no other lesion programme exists.

## Order of work

1. This plan and its contract are committed before the gates run.
2. The eligibility gates run and their tables are committed.
3. The owner retains or revises this plan.
4. Only then does scoring run, and its first step is the instrument check.
