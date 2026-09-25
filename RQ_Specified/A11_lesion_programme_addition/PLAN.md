# A11 pre-registered test in the Kim 2020 cohort

25 September 2026. **Pre-registered; not scored.** The owner chose this design at
stage 4 of the [shared component contract](../A5_A11_shared_component_contract/README.md)
(DEVELOPMENT decision 37). The machine-readable version is
[`config/kim2020_test_contract.json`](config/kim2020_test_contract.json). The owner now authorized revisions and execution after review. This amendment
is committed before new scores. See the shared [biological rationale](../A5_A11_shared_component_contract/BIOLOGICAL_LOGIC.md).

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
| Mean paired difference, log2 CPM | 0.327 |
| Standard deviation of differences | 0.384 |
| Patients positive | 19 |

These figures use the discovery run's own primary, its broad type 2 compartment
label, which the register also quotes. Its narrower type 2 label gives a mean of
0.323 with the same 19 positive patients and serves as a sensitivity.

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
is a changed population contrast, not automatically a stricter test. Either
direction may partly reflect different cellular composition. The author annotation has
no copy-number malignancy label, so the lesion arm may also include some
non-malignant epithelium.

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
within 1e-6, under both the broad and the narrower type 2 labels. If it cannot, Kim
is not scored and the discrepancy is reported.

For this check reproduce the original normalization scope: all retained
patient/histology units within each label/floor view, including histologies outside
the final LUAD-normal contrast. The original R script normalizes these together
before forming the 23 paired differences. Use its original human gene membership
for discovery; Kim has its own gene-index coverage. Hash inputs and record versions.

## Primary test

| Item | Choice |
|---|---|
| Module | lesion-specific, 91 mouse genes |
| Estimand | Location shift / pseudomedian of paired differences, lesion minus normal type 2 |
| Test | exact two-sided Wilcoxon signed-rank, alpha 0.05 |
| Estimate | Hodges-Lehmann with exact 95 percent interval |
| Mean sensitivity | arithmetic mean and paired t 95 percent interval, a different estimand |
| Smallest effect of interest | 0.10 log2 CPM, about a third of the discovery mean |

Signed-rank location inference assumes independent patients and a symmetric
distribution of paired differences. The Hodges–Lehmann estimator targets the
pseudomedian, not generally the arithmetic mean. See the
[R documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html).
If zeros/ties or interval failure prevent the installed implementation from giving
exact inference, report it unavailable; do not promote an approximate fallback.

The old overlapping categories are replaced by **two separate axes**:

- Direction: positive if exact p < 0.05 and HL > 0, negative if p < 0.05 and HL < 0,
  otherwise unresolved. Failed exact inference is unavailable.
- Magnitude: meaningful positive shift supported if exact CI lower > 0.10;
  positive shift of 0.10 ruled out if exact CI upper < 0.10; otherwise unresolved.
  An unavailable interval gives unavailable magnitude inference.

A small positive effect can have a supported direction while falling below the
planning margin without contradictory labels. The 0.10 margin is pragmatic and
based on discovery, not clinically validated. Ruling out this positive margin
is not two-sided equivalence to zero and does not retire all lesion programmes.

## Secondary tests

One family, Benjamini-Hochberg at q below 0.05, same test as the primary:

1. **Beyond shared.** Each patient's lesion-specific difference minus their
   shared-remodelling difference. This is relative expression activation, not
   conditional adjustment, independent mechanism or incremental prediction.
2. **Stress excluded.** The lesion-specific module with every gene of the three
   Hallmark control sets removed. Nineteen of its 91 genes sit in those sets,
   leaving 72. This addresses those memberships, not every stress/cycling
   response; a smaller module also has less sensitivity.
3. **Injury-lesion pair.** The seven-gene component shared by injury and lesions,
   five of them human-mappable.

An unavailable secondary enters the fixed three-test BH family as p=1 and is
marked unavailable. A positive primary means replicated lesion association.
Positive primary plus positive BH-significant beyond-shared and stress-excluded
results supports relative activation exceeding the nominated shared score and
surviving those exclusions. If the shared score decreases, do not describe this
as two co-activated programmes. None of these tests establishes specificity to
neoplasia without a non-neoplastic injury comparator.

Reported without a test: the shared module; type 1 and type 2 identity axes; the
three Hallmark axes; the leave-one-patient-out range of the primary; and the primary
without the one double-primary patient.

## Power, stated before any score

From a logged calculation using only the discovery's tracked differences:

| Scenario | Exact Wilcoxon | Paired t |
|---|--:|--:|
| Discovery effect, 8 patients | 0.47 | 0.54 |
| Half the discovery effect | 0.15 | 0.18 |

So even if the discovery effect transports in full, this test is close to a coin
flip. First estimates are often inflated, and at half the effect the chance is
about one in six. A null result will therefore be inconclusive unless the
positive-margin interval resolves magnitude. These are normal-location-model
planning numbers, not guaranteed power. A significant result is not inherently
more credible because power is low.

## What this cannot show

- Unresolved direction is not refutation; the magnitude margin is not equivalence to zero.
- Tumour epithelium is not confirmed malignant.
- Same-patient pairing controls the donor, not stage, smoking or driver mutation,
  which all vary across these patients.
- A positive primary supports lesion association in this defined population.
  Relative activation requires secondary evidence; neoplasia specificity is untested.

## Order of work

1. This plan and its contract are committed before the gates run.
2. The eligibility gates run and their tables are committed.
3. The owner authorized revision and execution; commit the amended plan before new scores.
4. Only then does scoring run, and its first step is the instrument check.

## Gate outcome

The gates ran on the pre-registered configuration and read gene names and cell
labels only. Kim's gene index holds 29,634 genes.

| Module | Role | Assayed fraction | Gate 1 |
|---|---|--:|---|
| Lesion-specific | primary | 0.802 | pass |
| Lesion-specific without stress genes | secondary | 0.778 | pass |
| Injury-lesion pair | secondary | 0.714 | pass |
| Shared remodelling | secondary baseline | 0.833 | pass |

Every identity and control axis also passes. For gate 2, eight patients meet the
50-cell floor in both arms, as the metadata indicated; P0008 and P0009 fall short on
tumour epithelium. The test is eligible; owner-authorized scoring first requires the instrument check.

## Correction made before scoring

The version committed in `e6b9dc9` cited the discovery run's narrower type 2 label
as its reference, at a mean of 0.323. The discovery run's own primary is its broad
type 2 compartment label, at 0.327, which is also the closer analogue of Kim's
pooled tumour epithelium. The reference, the power calculation and the instrument
check now use the broad label, with the narrower one kept as a sensitivity. The
first power output is preserved under `tables/superseded_power_label_AT2/`.

This was corrected before any Kim score existed. The Kim populations, test,
margin and decision rules did not change.


## Current amendment before scoring

The owner requested proceeding with review revisions and biological rationale.
The amendments align the estimand with Wilcoxon, separate direction/magnitude,
correct the population interpretation and narrow the biological claim. Module
membership, populations, test, margin, secondary family and original gate outputs
are retained. The 23-patient result remains discovery. Descriptive patient omissions
use already-normalized paired scores and do not constitute new tests or a separately
renormalized analysis. No favorable subtype or alternative test is selected after
viewing Kim.
