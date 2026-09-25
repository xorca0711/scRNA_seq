# Fibrosis epithelial state specificity

Both full IPF matrices were streamed into donor-by-deposited-state raw-count pseudobulks. Every epithelial unit total agrees with the earlier independent full-library stream. Source-defined HPCS, ADI, alveolar, stress and ISR programs use the frozen strict human/mouse ortholog map. Exact source-list denominators remain visible; missing genes are not interpreted as unexpressed.

The primary state comparisons are Aberrant_Basaloid versus ATII in GSE136831 and KRT5-/KRT17+ versus AT2 in GSE135893. Transitional AT2 versus AT2 is a separate secondary comparison. These source labels describe different populations; similar scores are not evidence that they represent one shared state.

Primary normalization uses all assayed genes in retained epithelial donor/state pseudobulks, with edgeR TMM and log2 CPM prior 1. The 30/100-cell and 0.5/2-prior analyses remain sensitivities. Comparisons retain paired IPF donors; no pooled-cell P values or cross-cohort equivalence test is used.

## Reduced HPCS results

- GSE135893, KRT5-/KRT17+ versus AT2: 3 paired donors, mean +0.955 log2 CPM units; 3/3 individual differences positive. Declared pair-count gate: eligible.
- GSE135893, Transitional AT2 versus AT2: 1 paired donors, mean +0.104 log2 CPM units; 1/1 individual differences positive. Declared pair-count gate: not met; individual values remain descriptive.
- GSE136831, Aberrant_Basaloid versus ATII: unevaluable at the primary setting (no_complete_same_donor_state_pairs); no zero effect is imputed.

These within-state RNA-program contrasts do not establish IL-1 dependence, fibrosis-to-cancer progression or malignant transformation. The annotated AT2 IPF/control contrast is provided separately and cannot stand in for an aberrant-state contrast.

![Fibrosis state specificity](../../figures/ipf_epithelial_state_specificity.png)

[Donor/state scores](donor_state_scores.csv), [paired values](within_donor_state_values.csv), [state summaries and sensitivities](state_contrast_summary.csv), [AT2 disease contrasts](AT2_disease_descriptive_contrasts.csv), [frozen labels and questions](specification.json).
