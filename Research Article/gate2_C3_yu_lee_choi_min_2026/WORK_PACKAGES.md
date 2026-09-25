# Remaining analysis execution

The owner authorized all planned work on 24 September 2026 and requested a
report after each work package finishes. No renewed approval is required.
Packages finish as evaluated, partially evaluable, or scientifically
unidentifiable with the precise missing information. An access check alone
does not count as completion of an otherwise feasible analysis.

| Package | Scope | Status |
|---|---|---|
| IPF contrasts (U5) | Two-cohort donor-level LR RNA-compatibility contrasts, components and cell-floor/prior-count robustness | Completed descriptively; 277 primary contrasts checked |
| IPF source/robustness (U5) | Native LR full-cell/cap/threshold checks, second-cohort pass and source/recipient/inhibitor profiles | Completed; 586 native calls, two cohorts, donor/pair and threshold checks passed |
| Mouse perturbation (U3/U4) | Mouse annotation, per-animal coverage and eligible blockade analyses | Completed eligible early niche analyses; 613 primary RNA contrasts checked, 0/46 primary pathways at q < 0.05; primary fraction lacks three eligible treated animals and KAC classifier |
| Human lesions (U5) | Human paired niches, source profiles, pathways and targets | Completed; 75 libraries, 1,214 native LR calls, 1080 primary paired RNA contrasts and 279 primary pathway tests |
| Spatial context (U5) | Human lesion and pathological post-viral measurements | Completed 56 human sections and nine post-viral matrices; region/neighborhood questions retain missing-annotation limits |
| Ligand targets (U5) | Eligible ligand-target analysis and stability | Completed for eligible IPF and human receivers; 274 IPF and 943 human candidate rows with omission checks |
| Specificity (U6) | State specificity and repair/fibrosis/neoplasia comparison | Completed HPCS/repair, source ISR, IPF state and paired human extensions; original Han expression arm remains conditional on missing matrices |
| Figures/synthesis (U6) | Figure gallery, evidence synthesis and completion register | Completed gallery and evidence register; ready for interpretation review |

The initial batch is preserved in [INITIAL_RUN_REPORT.md](INITIAL_RUN_REPORT.md).
Scientific constraints remain those in [FINAL_REVIEW.md](FINAL_REVIEW.md) and
the two analysis contracts. Sources and outputs remain paper-local; historical
raw data and results are not overwritten. Completion reports will be added
below and linked from the gallery.

## Completed: IPF RNA-compatibility contrasts (U5)

Two cohorts, 277 eligible primary resource/subtype contrasts, both resource
views and 30/50/100-cell and 0.5/1/2-prior sensitivities. All primary summaries
were independently checked from donor values. IL1B macrophage-to-fibroblast
compatibility differs by cohort; no uniform increase or causal effect is
established. [Report and figure](trials/u5_ipf_compatibility/REPORT.md).

## Completed: IPF full-cell robustness and sources (U5)

Both complete raw matrices were streamed once; retained source counts match
the audited pseudobulks. All-cell native LIANA, two sampling seeds, three
cell floors, three detection thresholds and both resources were evaluated in
586 calls. The source/context scan covers 243,472 and 89,326 annotated
IPF/control cells. Missing edges and absent sample pairs are distinguished;
13 empty second-cohort calls were independently verified as having no eligible
edge. [Report, figure and checks](trials/u5_liana_robustness/REPORT.md).

## Completed: IPF ligand-target prioritization (U5)

Checksum-verified NicheNet v2 prior, official Pearson statistic, source/receiver
expression gates and fully refitted donor-omission target sets. The 274 primary
candidate rows are observational rankings, not causal activation evidence.
GSE136831 AT2 and GSE135893 broad fibroblasts fail the target-count gate.
IL-1 is not the leading fit to upregulated fibroblast targets; better ranks
for downregulated genes cannot establish inhibition with an unsigned prior.
[Report and figure](trials/u5_ligand_targets/REPORT.md).

## Completed: treatment-blind mouse lineage inputs (U3)

All 28,243 QC cells received a reviewed broad-lineage or explicit ambiguous
label. Raw counts were aggregated into 197 sample/subtype units; symbol
collapsing and pseudobulk totals agree. The broad CSF1R+ compartment is not
macrophage-only, and mixed alveolar cells are not assigned as KAC/DATP.
At the primary 100-alveolar-cell floor only two treated animals remain (three controls remain).
The original KAC classifier is also unavailable in audited public sources.
[Annotations and rationale](trials/u3_lineage_annotation/frozen_cluster_annotations.csv)
and [coverage](trials/u3_lineage_annotation/sample_subtype_coverage.csv).

## Completed: eligible early mouse niche analyses (U4)

All 96 native LR calls and 613 primary RNA-compatibility contrasts passed
eligibility and numerical checks. No primary estimated-correlation pathway
passes q < 0.05 (46 tests); four fixed-correlation sensitivity tests do.
No receiver has enough significant DE targets for the conditional
ligand-target step. The phenotype and later-treatment-history gates remain
unresolved; these niche measurements cannot substitute for a KAC response.
[Report, figure and validation](trials/u4_mouse_niche/REPORT.md).

## Completed: HPCS and repair specificity extension (U6)

The reduced HPCS signature rises in both author-HPCS source identifiers and
all seven evaluable pooled repair/developmental libraries. This is evidence
against treating an increased HPCS signature as specific to neoplasia; it
does not assign HPCS identity to the repair cells. Technical-seed direction
agrees in 127/130 evaluable module contrasts. The pooled wells and sorting
gates are not independent animal replication.
[Report, figure and source-defined modules](trials/u6_specificity/REPORT.md).

## Completed: source ISR signature extension (U6)

Recovered the 129-symbol Han supplementary set and evaluated the complete
set and 105-gene overlap-removed sensitivity in the reusable contexts.
All 24 eligible module/unit contrasts retain direction in the second
technical seed; several directions differ between full and overlap-removed
sets. The original Han developmental expression comparison remains
conditional because the audited processed-data sources contain code and
local matrix references, without the matrices themselves.
[Report, figure and source discrepancy audit](trials/u6_isr_extension/REPORT.md).

## Completed: spatial and post-viral context measurements (U5)

All 56 human lesion sections and nine post-viral count matrices were
processed. All 23 RNA patient labels overlap the spatial series, with
consistent deposited demographic fields. Whole-section paired patient
values, independent samples, assay coverage and measured maps are available.
Pathology-region inference remains unavailable without independent labels;
post-viral neighborhood inference additionally lacks deposited coordinates.
[Report and figures](trials/u5_spatial_context/REPORT.md).

## Completed: fibrosis epithelial-state specificity (U6)

Full raw matrices yielded exact donor/state pseudobulks in both IPF cohorts.
The reduced HPCS program is higher in KRT5-/KRT17+ than AT2 cells in all three
eligible GSE135893 paired donors. The primary GSE136831 comparison has no
complete donor pairs at 50 cells per state, and the transitional-AT2
secondary comparison has one. All individual values and lower/higher
coverage sensitivities remain available; none of these states is relabelled
as HPCS or malignant.
[Report and figure](trials/u6_ipf_specificity/REPORT.md).

## Completed: human inputs, native LR and all-QC sources (U5)

All 75 GSE308103 libraries were processed using the frozen reference mapping;
1,214 native LR calls passed independent complex-detection eligibility checks.
The all-QC source scan includes unassigned cells and passes all 150
library/confidence-cutoff count-parity checks. Unassigned cells carry median
IL1B count fractions of 52–72% across histologies. Source attribution and
subtype comparisons therefore remain conditional on confident label retention.
[Source report and figure](trials/u5_human_sources/REPORT.md).

## Completed: pathway assay-coverage sensitivity (U6)

All 52 cohort/pathway combinations have at least 92% gene coverage. The
declared 50%, 70% and 80% coverage gates select identical sets, so unchanged
CAMERA fits and multiplicity families supply all three coverage sensitivities.
No repeated model run or threshold revision is necessary.
[Coverage table](trials/u6_completion/pathway_assay_fraction_sensitivity.csv)
and [verification](trials/u6_completion/pathway_coverage_validation.json).

## Completed: paired human niches, targets and specificity (U5/U6)

All primary and declared sensitivity stages completed. The primary pathway
family contains 279 tests with 0 global q<0.05 discoveries.
The 1080 primary RNA contrasts retain individual paired patients.
Conditional target fits yield 943 expression-supported candidate rows
across overlapping source scopes and receiver views. These are not independent
replications or causal activation tests. Source-program contrasts remain
separate from the unavailable KAC/NF-κB endpoint.
[Niche report](trials/u5_human_niche/REPORT.md) ·
[target report](trials/u5_human_ligand_targets/REPORT.md) ·
[program report](trials/u6_human_specificity/REPORT.md).

## Completed: gallery and evidence register

The [paper gallery](README.md#figure-gallery) includes context, source,
compatibility, pathway, perturbation, human/spatial and specificity panels.
The root README links to paper galleries without selective figure embeds.
The [evidence review](EVIDENCE_REVIEW.md) separates supported descriptive
findings, tested results without primary discoveries, and questions that
the public data cannot identify. The original criteria and historical claims
are unchanged. [Reproduction and provenance](REPRODUCIBILITY.md).
