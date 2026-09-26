# Erratum to the rival-2 stage 1 audit and its addendum

27 September 2026. Two statements in the committed stage 1 material are **wrong about the
deposit**, and one is wrong in a way that changed the design. Both were found by the
adversarial review of the draft freeze, before anything was scored, and both are corrected
here rather than by editing the committed records, whose bytes stay verifiable. The
[withdrawal report](RIVAL2_FREEZE_V1_WITHDRAWN.md) records the consequences for the freeze;
this page records the factual corrections.

The source of truth is the series summary of GSE190821, retrieved from the Entrez
`gds` index and quoted here in the parts that matter.

## Correction 1: the deposit does state genotype, dose and schedule

[`tables/rival2/stage1b_addendum_run.json`](../tables/rival2/stage1b_addendum_run.json)
lists, under `what_the_deposit_does_not_state`, both "per-mouse genotype" and "antibody dose,
schedule or measured target engagement". The second is wrong and the first is imprecise.

The deposit states:

- **Genotype, at cohort level.** RiboTag mice crossed to ShhCre, both on C57BL/6, with
  HA-tagged Rpl22 expressed conditionally in the epithelium. What the deposit does not give is
  a **per-mouse** genotype confirmation, so constancy across arms remains an assumption; the
  cohort construction itself is stated.
- **Bleomycin dose and route.** Intranasal, 3 U/kg, on day 1.
- **The 3G9 schedule and dose.** 10 mg/kg on day 1 and day 4.
- **The control antibody.** Axum8, inert, in saline.
- **The harvest.** Day 7 after bleomycin, flash-frozen, homogenised later under native
  conditions with cycloheximide.

Only **measured target engagement** is genuinely absent, and the corrected list separates that
from dose and schedule, which are recorded.

## Correction 2: immunoprecipitation purity is computable, not unverifiable

This is the correction that changed the design. The withdrawn freeze stated that the epithelial
compartment's "purity is not verified in the deposit" and treated composition as an
unresolvable caveat inherited from A10. That is factually wrong about this deposit.

The series summary states that **an aliquot of each mouse's homogenate was removed as the
input**, which is the library the deposit labels Whole Lung, and that **the remainder was
anti-HA immunoprecipitated**, which is the library labelled Epithelium. Every one of the 24
mice, and therefore all 11 in the arms used here, has a paired input from the same homogenate.

That is the standard RiboTag enrichment control, and it makes three things computable that the
withdrawn freeze had disclaimed:

1. **Per-mouse immunoprecipitation purity and enrichment**, as immunoprecipitation-minus-input
   log2 CPM for lineage-restricted genes: epithelial enrichment for `Epcam`, `Cdh1` and
   `Nkx2-1`, and de-enrichment for `Ptprc`, `Pecam1`, `Col1a1` and `Col3a1`.
2. **Per-mouse injury severity**, from the same input.
3. **A perturbation-engagement control**, from a whole-lung collagen and myofibroblast score,
   whose direction under this antibody is published.

All three are declared in
[`config/a15_rival2_freeze_v2.json`](../config/a15_rival2_freeze_v2.json) before any value is
read. The third is the most consequential: A15's own gate condition 1 requires a recorded
validation that the perturbation took effect, the withdrawn freeze had dropped it, and the
paired input is what lets v2 put it back.

## What remains genuinely unstated, corrected list

- harvest time of day, warm ischaemia time and homogenisation order, for every mouse. This
  matters because a chaperone and immediate-early signal is the standard readout of handling
  time, and the demoted A0 covariate is largely made of such genes.
- immunoprecipitation yield, input RNA mass and RNA integrity number, so **efficiency** stays
  unknown even though **purity** is computable. The two are different and the withdrawn freeze
  conflated them.
- per-mouse genotype confirmation, cage, litter and weight.
- measured target engagement, and how mice were allocated to antibody, so randomisation is not
  established.
- epithelial subtype composition. The paired input constrains this and does not remove it.

## What is not corrected, because it was right

The design audit's substantive findings stand and are not re-run: the arm membership, the
single-batch and sex matching, the bijective join and its independent agreement with A1's ten
recorded assignments, the two naming hazards, the depth audit, the exact power floors, the
endpoint coverage figures and the precedent finding that A1 used this deposit for a different
contrast. Stage 1 passed all ten of its stop rules and this erratum does not change that.
