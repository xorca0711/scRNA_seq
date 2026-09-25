# Biological logic linking A5 and A11

25 September 2026. This is the rationale for the revised tests, not a claim that
their hypotheses have been established. The original frozen partition is retained.

## Start from the biological problem

AT2 cells produce surfactant; thin AT1 cells provide the epithelial surface for gas
exchange. Building an alveolus at birth and restoring one after injury both demand
changes in epithelial identity, shape and interaction with the surrounding tissue.
Some transcriptional responses may therefore recur. Recurrence does not mean that
an adult injured cell has become a fetal progenitor.

The developmental source is specifically the Guo **postnatal-day-1 AT1/AT2
population**, not all lung development. Birth also brings oxygen exposure and
unfolded-protein stress. These are plausible explanations for an injury-like RNA
signature that do not require reuse of a developmental differentiation programme.
Guo provided spatial co-detection of alveolar markers, stronger context than the
list-overlap check alone; neither establishes the fate of every sequenced cell.
[Guo 2019](https://doi.org/10.1038/s41467-018-07770-1).

Adult injury induces Krt8-positive transitional cells between epithelial states.
Strunz reported convergence from alveolar and airway sources and persistence of
related cells in fibrosis. The same paper already reported poor correspondence
between developmental signatures (including Guo's) and its injury signature. Thus
A5 tests **partial recruitment of a developmental source programme**, with an
explicit prior counter-observation; it does not begin by assuming global identity.
[Strunz 2020, discussion](https://doi.org/10.1038/s41467-020-17358-3).

Neonatal and adult repair can also follow different lineage directions. Penkala's
neonatal injury work reported AT1-to-AT2 reprogramming, whereas adult alveolar repair
can involve AT2-to-AT1 differentiation. A shared RNA response is compatible with
different starting cells and destinations.
[Penkala 2021](https://doi.org/10.1016/j.stem.2021.04.026).

In tumours, epithelial cells experience altered identity, stress and tissue
organization. A lesion-derived gene list can therefore rise during ordinary repair
too. The existing 91-gene score did so in repair/development and IPF, as well as
19/23 LUAD pairs. A11 must separate **replication of lesion association**, **relative
activation beyond a nominated shared response**, and **specificity to neoplasia**.
These are three different claims. The present Kim experiment addresses the first
two, with limitations, and has no non-neoplastic injury arm for the third.

## Logical sequence and why each step is needed

| Step | Biological question | Analysis | What the result can establish |
|---|---|---|---|
| 1. Define sources | Which genes characterize the nominated developmental, injury and lesion populations? | Preserve published lists and operational-gene exclusions; report overlap | A reproducible candidate definition, not a universal programme |
| 2. Check independence | Did the proposed test data help choose the genes? | Trace lists to the original experiment, including genes removed | Whether evaluation is external to module selection |
| 3. Make a within-unit contrast | Do transitional or lesion-associated cells differ from a biologically relevant reference? | Same-mouse A5 and same-patient A11 contrasts | A paired expression association; cells are not replicates |
| 4. Challenge the explanation | Is the difference just AT1 maturation, generic injury activation or stress? | External identity/stress exclusions and reference sensitivity in A5; shared-response contrast and stress exclusions in A11 | How much of the nominated association survives specific rivals |
| 5. Quantify direction and size | Is there a positive effect, and is its magnitude resolved? | Estimates, intervals, fixed multiplicity and omission diagnostics | Supported association, bounded magnitude, or uncertainty |
| 6. Connect to outcomes | Does programme activity change actual repair or transformation? | Interpret alongside existing lineage/perturbation evidence; require linked outcomes for new causal claims | Context for a mechanism hypothesis, not causation from an RNA score |

## What the shared partition means

The 12-gene union contains 5 development–injury and 7 injury–lesion genes, with no
three-way intersection. Four of the first five overlap AT1 identity; the second
pair contains several stress-associated genes. Absence from a capped marker list
does not establish absence of expression. Keep `development_specific` and
`lesion_specific` as historical machine identifiers; in prose they mean
**source-list-exclusive candidates**, not proven biological specificity.

The common contract coordinates definitions and provenance. It does not require
one shared biological baseline, pool mouse and human units, or imply a sequence
from development through repair to cancer.

## A5: remove test-cohort selection from the primary gene set

The cached Strunz supplement description explicitly assigns Supplementary Data 3
to the high-resolution epithelial experiment. The contract reads its `cell_types_2`
sheet for ADI, AT1 and AT2 markers. Those are the proposed test data. Both positive
selection of the 5 shared genes and negative filtering of the 94/51-gene variants
therefore use that cohort. Negative filtering may reduce an expected effect, but
does not create held-out validation.

The revised primary uses all 99 Guo AT1/AT2 genes remaining after the already fixed
operational exclusions. Its identity-excluded sensitivity uses **Guo's own AT1 and
AT2 lists**, leaving 57 genes, followed by externally defined stress/cycling
exclusions. None of these definitions uses Strunz marker ranks. The old 94/51 and
pairwise sets remain descriptive references. This is an independently sourced
signature test in an existing, author-annotated cohort; it is not a new untouched
experiment or a validation of the authors' clustering.

Activated AT2 is the primary reference because both arms experience injury. A
resting-AT2 comparison asks a broader question and is a declared sensitivity.
Raw-count depth matching and one predeclared day window prevent technical depth
or a selected peak day from defining the answer. See [A5 plan](../A5_developmental_programme_reuse/PLAN.md).

## A11: preserve the distinction between association and addition

The 23-patient cohort is discovery. Kim supplies eight eligible patient pairs,
using all author tumour epithelial subtypes against normal AT2. This is a changed
population contrast, not automatically a stricter version of the discovery test.
It neither confirms malignant identity nor removes within-sample composition.

The primary tests the location shift of the lesion-derived candidate. A separate
comparison asks whether its paired change exceeds that of the shared union.
This is a **relative score contrast**, not proof of an independent mechanism or
incremental predictive information. Stress-gene exclusion only addresses the
listed controls; it cannot eliminate every stress response. A positive primary
alone is reported as replicated lesion association. See [A11 plan](../A11_lesion_programme_addition/PLAN.md).

## How the completed A1 evidence constrains both questions

A1's AP-1 analysis found opposite regional directions for HOPX acquisition after
perturbation; it did not measure completed AT1 function. Its TP53 comparison found
109 of 868 shared significant genes changing in opposite directions across
origins. These results oppose a universal interpretation that a shared marker
means the same state, or that reducing a transitional programme always improves
repair. New RNA results should nominate a context-dependent programme for a later
mechanistic test, not stand in for lineage, regulatory mediation or function.
See [A1 regulatory/fate report](../A1_transitional_epithelial_state_distinction/reports/REGULATORY_FATE_REPORT.md).

## Evidence and correction record

- Existing numerical precedents: ES1; the human U5/U6 specificity run; shared
  freeze and coverage tables; A1 regulatory/fate report. Reuse their biological
  units and distinguish previously observed results from fresh tests.
- Strunz archive: `raw_data/literature/PMC7366678_supplementary.zip`, member
  `41467_2020_17358_MOESM3_ESM.pdf`, description of Supplementary Data 3. The
  source marker workbook is `41467_2020_17358_MOESM6_ESM.xlsx`.
- The current A5 provenance wording supersedes the original audit's claim of no
  circularity; the current A11 plan supersedes overlapping outcome categories.
  Original outputs remain unchanged and the revisions precede new scores.
