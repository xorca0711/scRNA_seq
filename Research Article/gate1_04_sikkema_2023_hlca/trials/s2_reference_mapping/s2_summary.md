# Trial S2 output: scArches mapping of GSE178360 to the HLCA core and label transfer

Reference: 584944 HLCA core cells, embedding `X (30-dim scANVI embedding stored as X in HLCA_full_v1.1_emb.h5ad)`; transferred levels ['ann_level_1', 'ann_level_2', 'ann_level_3', 'ann_level_4', 'ann_level_5', 'ann_finest_level']; k = 50; Unknown if u > 0.2. Deviations are listed in the script docstring.

## Uncertainty per donor

| donor | n_cells | mean_uncertainty | median_uncertainty | unknown_fraction_u0.2 | unknown_fraction_u0.3 |
|---|---|---|---|---|---|
| DD046Q | 8957 | 0.1621 | 0.08 | 0.3346 | 0.2356 |
| DD047Q | 10493 | 0.0926 | 0.0 | 0.183 | 0.1243 |
| DD073R | 8279 | 0.156 | 0.0601 | 0.3255 | 0.2341 |

## Cluster-level transfer

| Cluster | Cells | Blind proposal | HLCA level 1 | Transferred finest (mode) | Mode fraction | Unknown fraction | Mean uncertainty | Verdict vs blind | S3 flat | S3 hier |
|--:|--:|---|---|---|--:|--:|--:|---|---|---|
| 0 | 1000 | Mono_Mac | Immune | Interstitial Mph perivascular | 0.354 | 0.49 | 0.2234 | agree | Interstitial Mph perivascular | Interstitial Mph perivascular |
| 1 | 377 | Mono_Mac | Immune | DC2 | 0.782 | 0.371 | 0.1815 | partial (same compartment) | DC2 | DC2 |
| 2 | 1173 | Endothelial | Endothelial | EC general capillary | 0.975 | 0.093 | 0.0609 | agree | EC general capillary | EC general capillary |
| 3 | 314 | Endothelial | Endothelial | EC arterial | 0.975 | 0.092 | 0.048 | agree | EC arterial | EC arterial |
| 4 | 1118 | CD8_T | Immune | CD8 T cells | 0.839 | 0.358 | 0.1708 | agree | CD8 T cells | CD8 T cells |
| 5 | 1326 | Fibroblast | Stroma | Alveolar fibroblasts | 0.78 | 0.336 | 0.169 | agree | Alveolar fibroblasts | Alveolar fibroblasts |
| 6 | 1231 | Fibroblast | Stroma | Peribronchial fibroblasts | 0.649 | 0.552 | 0.2526 | agree | Peribronchial fibroblasts | Peribronchial fibroblasts |
| 7 | 2096 | NK | Immune | NK cells | 0.775 | 0.148 | 0.0841 | agree | NK cells | NK cells |
| 8 | 1738 | T_cell | Immune | CD4 T cells | 0.851 | 0.253 | 0.1221 | agree | CD4 T cells | CD4 T cells |
| 9 | 1086 | Endothelial | Endothelial | EC general capillary | 0.766 | 0.302 | 0.157 | agree | EC general capillary | EC general capillary |
| 10 | 961 | Endothelial | Endothelial | EC aerocyte capillary | 0.815 | 0.267 | 0.1241 | agree | EC aerocyte capillary | EC aerocyte capillary |
| 11 | 807 | AT2 | Epithelial | AT2 | 0.648 | 0.278 | 0.1258 | agree | AT0 | AT2 |
| 12 | 1292 | Transitional | Epithelial | Basal resting | 0.682 | 0.332 | 0.1617 | agree | Basal resting | Basal resting |
| 13 | 730 | Club | Epithelial | pre-TB secretory | 0.432 | 0.677 | 0.3191 | agree | AT0 | EC aerocyte capillary |
| 14 | 893 | Neutrophil | Immune | Classical monocytes | 0.944 | 0.093 | 0.0559 | reference lacks identity | Classical monocytes | Classical monocytes |
| 15 | 446 | Platelet | Immune | DC2 | 0.269 | 0.608 | 0.3064 | reference lacks identity | EC aerocyte capillary | EC aerocyte capillary |
| 16 | 1198 | Club | Epithelial | pre-TB secretory | 0.831 | 0.407 | 0.185 | agree | AT0 | AT2 |
| 17 | 23 | Plasma | Stroma | Adventitial fibroblasts | 0.391 | 0.957 | 0.5635 | disagree | Subpleural fibroblasts | Subpleural fibroblasts |
| 18 | 486 | Pan_epithelial | Epithelial | Multiciliated (non-nasal) | 0.545 | 0.409 | 0.1925 | partial (same compartment) | AT0 | AT2 |
| 19 | 2359 | Pan_immune | Immune | Classical monocytes | 0.577 | 0.386 | 0.1848 | partial (same compartment) | Classical monocytes | Classical monocytes |
| 20 | 437 | B_cell | Immune | B cells | 0.961 | 0.037 | 0.0207 | agree | B cells | B cells |
| 21 | 281 | Proliferating | Immune | NK cells | 0.384 | 0.391 | 0.203 | partial (blind label is compartment-free) | NK cells | NK cells |
| 22 | 513 | B_cell | Immune | Mast cells | 0.998 | 0.006 | 0.0025 | partial (same compartment) | Mast cells | Mast cells |
| 23 | 205 | Lymphatic_EC | Endothelial | Lymphatic EC mature | 0.776 | 0.278 | 0.1304 | agree | Lymphatic EC mature | Lymphatic EC mature |
| 24 | 2266 | Ciliated | Epithelial | Multiciliated (non-nasal) | 0.995 | 0.005 | 0.0032 | agree | Multiciliated (non-nasal) | Multiciliated (non-nasal) |
| 25 | 481 | Ciliated | Epithelial | Multiciliated (non-nasal) | 0.96 | 0.083 | 0.0466 | agree | Multiciliated (non-nasal) | Multiciliated (non-nasal) |
| 26 | 663 | AT1 | Epithelial | AT1 | 0.949 | 0.09 | 0.0465 | agree | pre-TB secretory | AT1 |
| 27 | 1234 | Endothelial | Endothelial | EC venous systemic | 0.637 | 0.254 | 0.1253 | agree | EC venous pulmonary | EC venous pulmonary |
| 28 | 742 | SMC_pericyte | Stroma | Smooth muscle | 0.589 | 0.329 | 0.163 | agree | Smooth muscle | Smooth muscle |
| 29 | 31 | Neuroendocrine | Epithelial | Neuroendocrine | 0.387 | 0.258 | 0.1303 | agree | AT0 | Ionocyte |
| 30 | 222 | Ciliated | Epithelial | Multiciliated (non-nasal) | 0.982 | 0.018 | 0.0206 | agree | AT0 | AT2 |

Verdicts versus the blind proposals: agree: 23, disagree: 1, partial (blind label is compartment-free): 1, partial (same compartment): 4, reference lacks identity: 2.

## AT0 check per donor

| donor | strict_gate | transferred_AT0_all | transferred_AT0_u_le_0.2 | transferred_preTB_all | transferred_preTB_u_le_0.2 | gate_and_AT0 | AT0_over_gate | within_factor_2 | s3_flat_AT0 | s3_hier_AT0 |
|---|---|---|---|---|---|---|---|---|---|---|
| DD046Q | 596 | 206 | 26 | 783 | 354 | 10 | 0.044 | False | 1018 | 106 |
| DD047Q | 880 | 203 | 87 | 632 | 472 | 48 | 0.099 | False | 1631 | 167 |
| DD073R | 134 | 16 | 6 | 44 | 23 | 3 | 0.045 | False | 152 | 8 |

Epithelial subcluster 4 (existing AT0 candidate analogue, n=328): transferred AT0 75, pre-TB secretory 1, AT2 126, Unknown 117; top {'AT2': 126, 'Unknown': 117, 'AT0': 75, 'AT1': 6, 'Multiciliated (non-nasal)': 2}.

## Post hoc: our transfer versus the HLCA authors' own transfer for the same dataset (Tata_unpubl = GSE178360)

HLCA Tata_unpubl cells 21700; matched by unique 16-mer barcode: 19796 cells.

- level_3 (19648 cells labelled on both sides): agreement 0.9924 over those cells; 0.9999 when both sides are confident (fraction both confident 0.8989); mean uncertainty ours 0.046 vs HLCA 0.0488; Spearman of uncertainties 0.6882.
- level_4 (17789 cells labelled on both sides): agreement 0.9744 over those cells; 0.9999 when both sides are confident (fraction both confident 0.7487); mean uncertainty ours 0.1071 vs HLCA 0.1115; Spearman of uncertainties 0.8953.
- level_5 (5172 cells labelled on both sides): agreement 0.9778 over those cells; 0.9997 when both sides are confident (fraction both confident 0.7069); mean uncertainty ours 0.1215 vs HLCA 0.1285; Spearman of uncertainties 0.8749.
