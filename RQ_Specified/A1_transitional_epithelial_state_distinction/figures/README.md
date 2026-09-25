# A1 figure gallery

Four figure groups were generated from real inputs on 25 September 2026.
Read the [batch report](../reports/FIRST_BATCH_REPORT.md) for models, results,
source hashes and limits. PNGs below have adjacent SVG versions.

## IRE1α perturbation: sample structure and treatment effects

![Epithelial RiboTag sample PCA and treatment effects](a1_ire1_rna_overview.png)

Five vehicle and five KIRA8 mice, GSE190821. PCA uses the 2,000 most variable
retained genes without batch removal; PC1 explains 89.5% and separates batches.
The volcano uses a batch + sex + treatment edgeR model and all-gene BH FDR.
Four of 14,811 tested genes pass FDR < 0.05. These are day-7 epithelial
ribosome-associated RNA measurements, not traced-cell fate or chromatin assays.
[Coordinates](../tables/ire1/sample_PCA.tsv),
[gene effects](../tables/ire1/gene_effects.tsv),
[SVG](a1_ire1_rna_overview.svg).

## Predefined molecular endpoints

![Mouse-level predefined marker expression](a1_ire1_marker_panel.png)

All eight frozen markers are shown, with one point per mouse and group means.
Shapes identify batch. Display values are unadjusted TMM log2CPM; effects and
q-values come from the adjusted model and the full tested-gene family. None
passes FDR < 0.05. Apparent directions do not establish restored cell fate.
[Plotted values](../tables/ire1/marker_logCPM.tsv),
[effects](../tables/ire1/predefined_marker_effects.tsv),
[SVG](a1_ire1_marker_panel.svg).

## Measured PATS lineage endpoint reconstruction

![Krt19 lineage endpoint fractions per mouse](a1_pats_lineage_endpoints.png)

Kobayashi Extended Data 4 source data, labelled BleoD12. Three named mice per
marker; fields remain nested within each mouse. Filled circles average field
fractions; open squares pool counts within that mouse. Control denominators
are zero, so their fractions are undefined. No between-marker pairing,
mutually exclusive composition or transition-rate estimate is assumed.
[Mouse table](../tables/lineage/pats_mouse_endpoints.tsv),
[field table](../tables/lineage/pats_source_fields.tsv),
[SVG](a1_pats_lineage_endpoints.svg).

## Descriptive ATAC and protein-sorted RNA profiles

![Deposited-source PCA for TIGIT ATAC and CD44 RNA](a1_deposited_source_PCA.png)

PCA of deposited integer counts after CPM/log2 transformation, using 2,000
variable features. Lines connect literal source aliases. ATAC pool independence
and CD44 alias-to-genotype mapping are unresolved; neither plot supplies
biological replication or an inferential contrast. Labels and coordinates:
[ATAC](../tables/descriptive/GSE154966_source_PCA_QC.tsv),
[CD44 RNA](../tables/descriptive/GSE273123_source_PCA_QC.tsv),
[SVG](a1_deposited_source_PCA.svg).

No biological peak-overlap plot was made for GSE141635: the deposited histone
calls use incompatible region settings. The
[technical audit](../tables/descriptive/H3K4me3_technical_geometry.tsv) records why.

## Remaining figure plan

The existing
[RNA/accessibility figure](../../../analysis/figures/rq/rq_a1_chromatin.png)
remains a historical descriptive result; it does not answer the expanded
histone-modification or fate questions.

| Proposed figure | Main evidence | Presentation |
|---|---|---|
| A1-1 | Sampling and modality-specific state maps | Assay/design diagram; separate RNA and ATAC UMAPs; cross-assignment matrix |
| A1-2 | Sample-level accessible regulatory programmes | PCA with paired sources joined; accessibility effects; motif heatmap; genome tracks |
| A1-3 | Direct histone marks and methylation reference | Mark-specific locus tracks and heatmaps; observed methylation fractions with coverage |
| A1-4 | Measured lineage, then testable time/topology | Mouse-level endpoints first; measured descendant matrix; trajectory only with a validation endpoint |
| A1-5 | Functional perturbation and linked phenotype | RNA response, eligible pathways and separately measured differentiation/fibroblast outcomes |
| A1-6 | Optional tissue proteins and cross-assay synthesis | Donor/ROI protein, morphology and regional effects; measured-versus-inferred evidence matrix |

Each figure requires plotted values, sample counts, an independent-unit
definition, source/code hashes and visual review. Cell-level violins are
descriptive; donor/animal points carry inference. No embedding should imply a
trajectory without independent temporal/fate evidence. The [analysis plan](../PLAN.md)
specifies figure-specific eligibility and sensitivity checks.
