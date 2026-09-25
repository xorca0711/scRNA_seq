# A1 first analysis batch — 25 September 2026

**Completed:** measured-lineage source reconstruction, the first eligible
treatment count model, descriptive ATAC/RNA sample profiles and a direct-histone
input audit. The larger A1 epigenetic/state-distinction analysis remains open.
Branch: `codex/a1-state-distinction-plan`; folder PR #66 is separate.

The owner's request to challenge stages 3–4 and initiate analysis authorized
this batch. The revised [plan](../PLAN.md) and [lineage audit](../LINEAGE_AUDIT.md)
put measured ancestry/descendant endpoints ahead of computational trajectories,
and same-study functional evidence ahead of broad spatial/proteomic screening.
The targeted review covers Choi, Kobayashi, Strunz, Kathiriya, Auyeung, Chan,
Krt8 perturbation and CD44 phenotype studies. Their assay capabilities and
experimental contexts are not interchangeable.

## What ran

| Work | Real input and unit | Result | Limit |
|---|---|---|---|
| PATS lineage | Kobayashi Extended Data 4 source workbook; three explicitly named mice per marker, nested fields | 38 fields reconstructed; 18 zero-denominator control fields excluded from fraction estimates | Published source reproduction; no new tracing experiment or exit-rate measurement |
| IRE1α treatment RNA | GSE190821 epithelial RiboTag; five vehicle and five KIRA8 mice | 14,811 genes tested; four pass whole-family BH FDR < 0.05, all lower with KIRA8 | Day-7 bulk epithelial ribosome-associated RNA; raw reads and RNA-quality metrics not reprocessed |
| HPCS bulk ATAC | GSE154966: 112,729 features, eight libraries, four deposited source aliases | Integer-count audit and sample PCA completed | Pool/animal independence remains unresolved; no differential-accessibility tests |
| CD44-sorted RNA | GSE273123: 27,179 features, 16 libraries, eight deposited source aliases | Integer-count audit and sample PCA completed | R26/OG column aliases not yet mapped conclusively to GEO genotype/sample identities |
| Direct H3K4me3 | GSE141635 deposited called intervals | Different calling parameters verified in the original headers | No biological interpretation of raw interval counts, widths or overlap |

All downloaded inputs have source URLs and hashes in the
[processed-file inventory](processed_input_inventory.json) and
[workbook inventory](lineage_source_inventory.json). Four figure groups are in
the [gallery](../figures/README.md), with numerical tables and PNG/SVG versions.

## IRE1α: specified comparison and observed results

The [frozen contract](../config/ire1_kira8.json) preceded expression fitting.
GEO mouse IDs, tissue compartments and treatments map each selected count column
to one GSM in the [manifest](../tables/ire1/sample_manifest.tsv). Both treatment
arms contain three S061 and two S135 mice, with the same sex composition.
Axum8 antibody controls and whole-lung libraries are excluded. Whole-lung
libraries from the same mice would not be independent replication.

The primary model is `~ batch + sex + group`, with TMM, robust edgeR
quasi-likelihood and the frozen expression filter. Its design has rank four
and six residual degrees of freedom. Positive log2FC means KIRA8 minus vehicle.
The sex-omitted sensitivity retains batch adjustment; its all-gene log2FC
correlation with the primary model is 0.9992. This checks one model choice,
not independence from batch or sensitivity to individual mice.

None of the predefined markers passes the all-gene BH threshold:

| Marker | log2FC | BH FDR |
|---|---:|---:|
| Itgb6 | −0.811 | 0.143 |
| Krt8 | −0.566 | 0.171 |
| Krt19 | −0.604 | 0.206 |
| Cldn4 | −0.981 | 0.116 |
| Ager | +0.255 | 0.707 |
| Sftpc | +0.659 | 0.387 |
| Cdkn1a | −0.903 | 0.104 |
| Krt7 | −0.526 | 0.128 |

The published TGF-β signature is directionally lower but does not pass CAMERA
with estimated correlation (47 genes; BH FDR 0.125). The general unfolded-protein
response set also does not pass (94 genes; BH FDR 0.529). The terminal-UPR set is
held: six of eight source symbols enter the tested universe (75%), below the
frozen 80% rule. `Spa5` has no exact Ensembl match; a second mapped gene fails
expression filtering. No spelling repair or reduced threshold was applied.
Full [coverage](../tables/ire1/gene_set_coverage.tsv),
[mapping](../tables/ire1/gene_set_mapping_tested.tsv),
[pathway results](../tables/ire1/gene_set_camera.tsv) and
[gene effects](../tables/ire1/gene_effects.tsv) are retained.

PC1 explains 89.5% of unadjusted sample variance and separates the two batches.
Treatment is represented in both, so the additive contrast is estimable, but
strong batch structure limits transportability. The count model does not
separate changing cell proportions from within-cell expression or translation
changes. No effect is attributed to chromatin or irreversible state arrest.
Non-significance also does not establish equivalence or refute the paper's
different assays/statistical questions.

The published Krt19 tracing uses a day-3/4 pulse and day-14 AGER endpoint;
the RNA experiment is day-7 epithelial RiboTag. These are different measurements
and experimental units. See [Auyeung et al. primary figure legends](https://pubmed.ncbi.nlm.nih.gov/35170357/).
This reanalysis is source-informed; it is not unseen-data preregistration or
independent confirmation of the published mechanism.

## PATS: measured endpoint, correct denominator

The workbook labels the injury fields `BleoD12`. Mean within-mouse field
fractions, then averaged across three mice, reproduce 64.69% KRT8-positive and
32.84% AGER-positive among alveolar tdTomato-labelled cells. Count-pooled
within-mouse summaries are shown alongside these to expose weighting choices.
The markers are not assumed exclusive or paired across panels; their values
are not added into a state-composition model. Pulse timing still needs exact
final-paper protocol reconciliation before calculating a chase-dependent rate.

The 18 control fields have zero labelled-cell denominators despite stored
percentages of zero. Their fractions are undefined. They cannot serve as
zero-percent controls in a fate contrast. See the
[field table](../tables/lineage/pats_source_fields.tsv),
[mouse table](../tables/lineage/pats_mouse_endpoints.tsv) and
[source study](https://www.nature.com/articles/s41556-020-0542-8).

## What remains before stronger conclusions

1. **Direct epigenetic analysis remains central.** Quantify compatible common
   regions with H3/input controls, or use consistently normalized direct-mark
   tracks with explicit low-replication limits. GSE141635 homeostasis uses
   `size 1000/minDist 2000`; CTGF-positive injury uses `size 4000/minDist 4000`,
   with other parameter differences. This first batch cannot compare their
   called-region abundance biologically. Tsutsui histone profiles and the
   normal-differentiation methylation reference remain separate planned tasks.
2. **Check treatment stability before promoting a molecular conclusion.** Add
   within-batch and leave-one-mouse-out sensitivities under a new dated contract;
   retain all mice and the current primary results. Do not select a favourable
   gene-set test or combine unrelated controls to obtain significance.
3. **Resolve source identities for ATAC/CD44 inference.** Retain descriptive
   source-alias plots until pooled-animal independence and genotype crosswalks
   are verified. The public Katzen repository inspected here provides scRNA
   code but did not resolve the bulk R26/OG crosswalk.
4. **Recover exact tracing endpoints.** Audit HPCS reporter/chase/sample mapping
   and numerical descendant composition. The downloaded Chan workbook supplies
   growth/ablation panels, not its Figure 2 lineage fractions. Add trajectories
   only where observed time or lineage data can test them; spatial/protein work
   stays conditional on a specific contrast and donor/ROI map.

## Execution evidence

[Lineage run](lineage_run.json), [IRE1 run](ire1_run.json),
[descriptive run](descriptive_run.json), [figure render](figure_render.json),
[R session](../tables/ire1/R_session.txt) and
[executed-source archive](execution_sources/2026-09-25/manifest.json) record the
batch. The archive retains exact executed code before adding overwrite guards
and removing overlapping figure labels. Figure-only rendering verified that
all scientific table hashes remained unchanged. Completed numerical runs refuse
silent overwrite; archive a prior batch before rerunning it.

Seven targeted tests cover zero denominators, incorrect treatment/library
selection, duplicate animals, unresolved identities, overlapping pools and
count alignment. All passed in the scientific runtime. The four new integration
tests skip explicitly when the lightweight CI environment lacks scientific
packages. Numerical fits completed with edgeR 4.10.5 and limma 3.68.5; R emitted
locale fallback notices, retained in its log, without a model failure.

Repository checks: 2,105 validation checks and 18 claim numeric bindings passed.
The lightweight unit suite passed 13 tests and explicitly skipped the scientific
integration module; the seven targeted scientific-runtime tests passed separately.
All four figure groups were visually inspected; overlapping PCA point labels
were removed without changing coordinates or analysis tables.
