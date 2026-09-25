# Measurement contracts supporting the research questions

25 September 2026. This is an index of decision-relevant checks, not a new
analysis protocol or a second question register. The biological hypotheses live
in [RESEARCH_QUESTIONS.md](../RESEARCH_QUESTIONS.md). Existing trial contracts,
run records and [claim grades](../CLAIMS.md) remain authoritative for past work.
No eligibility threshold is relaxed by this reorganisation.

Complete a check for its specified input and estimand, record the decision, then
stop. Repeat only when those inputs/definitions change or a concrete discrepancy
could change the conclusion. Additional sensitivity searches do not substitute
for an independent endpoint or missing biological replication.

<a id="mc1"></a>

## MC1. Biological units, confounding and endpoints

Cells describe heterogeneity; animals, donors, patients or independent culture
preparations support the corresponding inference. Preserve pairs and repeated
measures. Technical seeds, libraries from the same pool, organoids from one
preparation and multiple histologies from one patient do not increase independent n.
Require an estimable design before fitting; regression cannot separate perfectly
confounded age, injury time and processing.

- [Reproduction guide](../REPRODUCIBILITY.md) and
  [sample-level inference](../analysis/corrections/statistics/README.md): G1/W1 late
  contrasts retain confounding; reference CAMERA W1 has no significant sets.
  A seven-gene ornithine set failing the minimum size is not a negative result.
- [A1 execution report](../RQ_Specified/A1_transitional_epithelial_state_distinction/reports/FIRST_BATCH_REPORT.md):
  IRE1α is day-7 epithelial RiboTag, five KIRA8 versus five vehicle mice,
  with batch and sex terms; exclude Axum8 antibody controls. PATS source endpoints
  use measured descendants, and control 0/0 fractions remain undefined.
- [Outcome-data gate](NEXT_DATASET_GATE.md): verify well, guide, plate and biological
  preparation joins and harvest timing before A10 prediction. Hold out whole
  preparations; outcome-time RNA is concurrent association.
- [Niche report](../Research%20Article/gate2_C3_yu_lee_choi_min_2026/trials/u5_human_niche/REPORT.md):
  count exact complete paired triads before A13. IPF counts six and three remain
  below the ten-unit joint-model floor. Three-per-arm and ten-unit floors in
  their respective protocols are eligibility rules, not power guarantees.

**Decision:** freeze a meaningful effect, equivalence or useful-prediction margin
and assess precision before calling a hypothesis unsupported. Insufficient
coverage gives a coverage report. It does not authorize cell-level substitution,
an underidentified model or a claim of biological absence.

<a id="mc2"></a>

## MC2. Depth, reference definition and assay comparability

Use the existing [ligand correction](../analysis/corrections/ligand/README.md),
[ES1 definitions](../Research%20Article/epithelial_state_specificity/README.md) and
[A1 assay plan](../RQ_Specified/A1_transitional_epithelial_state_distinction/PLAN.md).
Match the control to the biological contrast; a changed reference may itself be
biology (A7), not merely a measurement nuisance.

Recorded checks that must retain their scope:

- C37 source comparison: unadjusted epithelial/myeloid medians 0.4591/0.2715,
  p=0.001953; primary 1,000-UMI 0.1224/0.0999, p=0.130859. All ten donors and
  16,064 cells remain at 500/1,000 UMIs. Only the 2,000-UMI sensitivity excludes
  1,188 cells. The 500/2,000 results do not replace the primary budget.
- C131 M3 background-centred contrasts (-37.8/-18.3/-15.5%) differ from raw
  contrasts (-31.8/-16.0/-8.5%). The uninjured comparator is a seven-week
  **Cebpa mutant**, not healthy wild type. No well passes every gate (C133);
  reference-cell split intervals are not animal confidence intervals.
- C119 common-depth neonatal control/mutant fractions (3.69%/8.07%) differ from
  the raw fractions in the A5 display. Mutant development is not normal development.
- A1 H3K4me3 peak callers differ between homeostasis and injury; peak geometry
  cannot be interpreted biologically without compatible inputs. ATAC, direct
  histone marks, methylation and 3D contacts measure different properties.

**Decision:** reuse the frozen comparison when applicable. If coverage or assay
comparability fails, narrow the estimand or hold the test. Never reinterpret an
alternative depth budget as an independent validation or infer temporal closure
from cross-sectional RNA/accessibility measurements.

<a id="mc3"></a>

## MC3. State and source identity

Retain source labels alongside reviewed alternatives; record the evidence and
version for each. Harmonize biological states before comparing them, and avoid
defining a state using the same programme being tested. Histology alone does not
prove malignant cell identity; RNA trajectories do not establish ancestry.
The [lineage audit](../RQ_Specified/A1_transitional_epithelial_state_distinction/LINEAGE_AUDIT.md)
separates measured labels/descendants from computational trajectories.

<a id="a12-s1"></a>

### A12-S1. What are the unassigned IL1B-expressing cells?

The [annotation review](../Research%20Article/gate2_C3_yu_lee_choi_min_2026/trials/u5_human_full/ANNOTATION_REVIEW.md)
and [source report](../Research%20Article/gate2_C3_yu_lee_choi_min_2026/trials/u5_human_sources/REPORT.md)
show unassigned categories carrying median 52–72% of recovered IL1B counts across
human histologies. This is source-attribution evidence, not a secretion fraction
or a newly discovered macrophage state.

Review broad identity, confidence, independent markers, mixed profiles and
assay-appropriate background. Preserve ambiguous cells and original labels.
The [source figure](../analysis/figures/rq/README.md#a12-s1) uses 34,178 display
cells capped at 75 per patient/histology/source; numerical fractions use all QC
cells. A UMAP cluster or a relaxed confidence threshold does not validate identity.

**Decision:** a coherent independently supported phenotype may motivate a new
biological state hypothesis. Until then, retain A12-S1 as an enabling question
and qualify strong macrophage-source interpretations in A12/A13.

<a id="mc4"></a>

## MC4. Ligand–receptor representation and functional specificity

The [ligand-resource report](../analysis/corrections/ligand/RESULTS.md) separates
absent database entries, expression eligibility and ranking. Exact EGFR and
EGFR_ERBB2 are different definitions. AREG–EGFR_ERBB2 is present in consensus but
scored in only 7/22 donors, below the 11-donor rule; resource absence is not the
explanation. AREG's CellChatDB donor-median rank is 10.5 among 312 retained pairs.
Shared expression data and curation do not become independent experiments when
run through several databases.

For A12 use the [recipient/target methods](../analysis/figures/rq/il1b_context/REPORT.md):
keep ligand, receptor subunits and inhibitors visible, use a frozen candidate
universe per receiver and retain alternative ligand families. The focused human
up-target figure contains 32 eligible candidates across three receivers; IL1B
fails expression eligibility in the macrophage panel. An unsigned prior fit is
neither activation nor inhibition, and ranks with different denominators are
not comparable effect sizes.

**Decision:** RNA compatibility nominates experiments. Ligand secretion,
processing, receptor complexes, activity and necessity need their appropriate
protein or perturbation endpoints. Spatial proximity alone does not supply them.
Retire a representation-only question once its source/eligibility discrepancy is
resolved; the functional A2/A9/A12 hypotheses retain separate decisions.

<a id="mc5"></a>

## MC5. Programme dependence, enrichment and validation

The [ES1 analysis](../Research%20Article/epithelial_state_specificity/README.md)
preserves complete source lists versus short panels, overlap membership and
label-excluded versions. ADI/AT1's 119 shared genes quantify definition overlap,
not developmental continuity. A frozen component must be tested outside the data
used to choose it; do not define its validation endpoint from its own genes.

The [statistical correction](../analysis/corrections/statistics/README.md) and
[human pathway report](../Research%20Article/gate2_C3_yu_lee_choi_min_2026/trials/u5_human_niche/REPORT.md)
keep the primary estimated-correlation analysis separate from fixed-0.01
sensitivity. G2's 30 frozen candidates under fixed correlation do not constitute
primary confirmation (C161); human results are 0/279 versus 148/279 at q<0.05.
Show eligible set sizes, direction, multiplicity family and sample coverage,
including ineligible sets. Do not select the model because it gives discoveries.

**Decision:** a full programme beyond selection genes can motivate biological
sharing (A5/A11), while list arithmetic alone provides weaker motivation (A8).
Freeze biological alternatives and useful-effect/prediction margins before the
next fit. New hypotheses keep their own evidence status.

<a id="mc6"></a>

## MC6. Figures, captions and provenance

The root register is authored. The [shared gallery](../analysis/figures/rq/README.md)
owns curated captions and links; scripts own image files, plotted tables and
generated panel facts in run records. Script 16 never edits the root register or
inserts missing blocks. Script 18 owns the current A2 figure; the old script-16
A2 remains historical. Explicit selection and fresh render directories prevent
an A3 title correction from rerunning every figure or overwriting old records.

UMAP locates heterogeneity; PCA shows sample structure. Neither tests mechanism.
Use donor/animal points and uncertainty where the design supports inference;
cell-level violins remain descriptive. Enrichment figures retain primary and
sensitivity results. Experimental schematics contain no invented result curves.
Use the shared palette and preserve source tables, original code and run hashes.
The [migration record](migrations/2026-09-25-rq-reframing/README.md) documents the
presentation-only A3 correction and preservation checks.

<a id="crosswalk"></a>

## ID crosswalk

All A identifiers remain in their original order. This table maps the old emphasis
to the new biological decision; it is not a count of accepted or rejected hypotheses.
A12-S1 stays an explicit enabling question. A14 has two separable hypotheses.

| ID | Previous emphasis | Biological decision retained in the register | Supporting contract |
|---|---|---|---|
| A1 | Regulatory/phenotypic distinction beyond markers | Regulatory features and independently measured response/fate | MC1, MC2 |
| A2 | AREG expression and resource ranks | Context-dependent functional source contribution | MC2–MC4 |
| A3 | Phase programmes and confounding | Late injury-associated macrophage biology beyond aging | MC1, MC5 |
| A4 | Axin2/Il1r1 overlap | Sequential Wnt maintenance and IL-1 response in a lineage | MC1–MC3 |
| A5 | Injury-specificity of signatures | Developmental programme reuse | MC1, MC2, MC5 |
| A6 | Proliferation and subtype composition | Within-state IPF change beyond abundance | MC1, MC3, MC5 |
| A7 | Shifted genotype reference | Broad versus state-selective Cebpa identity effects | MC1–MC3 |
| A8 | ADI/AT1 score overlap | Added maturation-endpoint information | MC1, MC5 |
| A9 | Receptor representation/coverage | Receptor-dependent functional AREG response | MC3, MC4 |
| A10 | Organoid outcome association | Molecular information about measured growth across preparations | MC1, MC5 |
| A11 | Shared versus neoplasia-associated scores | Shared remodelling with independently evaluated lesion additions | MC1, MC3, MC5 |
| A12 | Recipient compatibility and IL-1 specificity | Recipient context beyond ligand RNA | MC1, MC3–MC5 |
| A12-S1 | Unknown-label IL1B sources | Source identity enabling attribution; no presumed novel state | MC3 |
| A13 | Reciprocal niche associations | Added fibroblast information beyond macrophage IL1B | MC1, MC3–MC5 |
| A14 | Withdrawal and fibroblast dependence | Separate exposure-duration and recipient-specific recovery decisions | MC1, MC3, MC4 |
