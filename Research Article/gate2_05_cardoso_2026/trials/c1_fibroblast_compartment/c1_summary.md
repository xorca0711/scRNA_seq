# Trial C1 output: the fibroblast compartment of GSE316241

**Verdict by the frozen rule: the reprogrammed fibroblast state is NOT recovered.** Shape of the Red2Kras fibroblast score: **two states** (BIC 1 component 8.9, 2 components -495.2).

Design constraint carried from trial C0: one library per genotype, three mice pooled per
library, one time point. Every number below describes two libraries, not two groups of
animals; no P value is computed for the genotype contrast.

## QC per library

| library | genotype | barcodes | median_counts | median_genes | median_pct_mt | kept_repository_rule | kept_paper_rule | kept_by_both | expected_doublet_rate | doublet_call_method | doublets_removed | cells_analysed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Expt1_Confetti_mesenchyme | Confetti | 6435 | 4529 | 2111 | 1.18 | 6057 | 5390 | 5300 | 0.0485 | top 4.85% of scores (automatic rejected: it called 0.07%) | 306 | 5751 |
| Expt1_Red2Kras_mesenchyme | Red2Kras | 6791 | 3943 | 1824 | 1.4 | 6260 | 4797 | 4681 | 0.0501 | top 5.01% of scores (automatic rejected: it called 0.00%) | 321 | 5939 |

## Mesenchymal clusters (Leiden 0.5, labels from the paper's marker sets)

| cluster | n_cells | call | mode_fraction | confident | fraction_Red2Kras | fraction_Confetti |
|---|---|---|---|---|---|---|
| 0 | 1709 | alveolar fibroblast | 0.506 | True | 0.414 | 0.586 |
| 1 | 2017 | alveolar fibroblast | 0.973 | True | 0.223 | 0.777 |
| 2 | 1547 | adventitial fibroblast | 0.97 | True | 0.556 | 0.444 |
| 3 | 733 | mesothelium | 0.966 | True | 0.501 | 0.499 |
| 4 | 331 | mesothelial-like | 0.489 | False | 0.988 | 0.012 |
| 5 | 702 | peri-bronchial fibroblast | 0.806 | True | 0.527 | 0.473 |
| 6 | 1519 | smooth muscle | 0.982 | True | 0.504 | 0.496 |
| 7 | 201 | pericyte | 0.851 | True | 0.517 | 0.483 |
| 8 | 79 | reprogrammed fibroblast | 0.418 | False | 0.747 | 0.253 |
| 9 | 1008 | alveolar fibroblast | 0.575 | True | 0.192 | 0.808 |
| 10 | 288 | mesothelial-like | 0.41 | False | 0.951 | 0.049 |
| 11 | 184 | mesothelial-like | 0.342 | False | 0.989 | 0.011 |
| 12 | 39 | reprogrammed fibroblast | 0.385 | False | 0.436 | 0.564 |
| 13 | 86 | mesothelium | 0.419 | False | 0.523 | 0.477 |
| 14 | 937 | alveolar fibroblast | 0.394 | False | 0.996 | 0.004 |
| 15 | 167 | reprogrammed fibroblast | 0.401 | False | 0.874 | 0.126 |
| 16 | 143 | proliferating | 0.301 | False | 0.965 | 0.035 |

## Fibroblast subclusters

| cluster | n_cells | call | mode_fraction | confident | fraction_Red2Kras | fraction_Confetti | mean_reprogrammed_score | mean_alveolar_score | median_doublet_score | det_Tnc | det_Fst | det_Runx1 | det_Runx2 | det_Acta2 | det_Pdgfrb | det_Pdgfra | det_Col13a1 | det_Lcn2 | det_Saa3 | det_Cxcl12 | det_pericyte_max | det_smooth_muscle_max |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 858 | alveolar fibroblast | 0.691 | True | 0.337 | 0.663 | 0.0623 | 0.5513 | 0.0354 | 0.156 | 0.01 | 0.389 | 0.031 | 0.387 | 0.097 | 0.844 | 0.731 | 0.019 | 0.013 | 0.417 | 0.029 | 0.481 |
| 1 | 1443 | adventitial fibroblast | 0.933 | True | 0.563 | 0.437 | 0.0861 | -0.5523 | 0.0351 | 0.011 | 0.227 | 0.175 | 0.026 | 0.306 | 0.23 | 0.581 | 0.03 | 0.014 | 0.018 | 0.693 | 0.126 | 0.184 |
| 2 | 777 | adventitial fibroblast | 0.499 | False | 0.559 | 0.441 | 0.1749 | -0.251 | 0.0354 | 0.077 | 0.026 | 0.589 | 0.069 | 0.452 | 0.269 | 0.746 | 0.485 | 0.005 | 0.014 | 0.475 | 0.027 | 0.306 |
| 3 | 1802 | alveolar fibroblast | 0.95 | True | 0.235 | 0.765 | -0.0242 | 1.1626 | 0.0351 | 0.25 | 0.002 | 0.07 | 0.028 | 0.329 | 0.045 | 0.908 | 0.754 | 0.007 | 0.007 | 0.102 | 0.019 | 0.649 |
| 4 | 448 | adventitial fibroblast | 0.652 | True | 0.33 | 0.67 | 0.0681 | -0.1569 | 0.0375 | 0.02 | 0.022 | 0.384 | 0.087 | 0.411 | 0.279 | 0.783 | 0.473 | 0.009 | 0.007 | 0.554 | 0.031 | 0.27 |
| 5 | 875 | alveolar fibroblast | 0.561 | True | 0.109 | 0.891 | 0.0088 | 0.4004 | 0.0462 | 0.026 | 0.011 | 0.273 | 0.029 | 0.397 | 0.154 | 0.822 | 0.683 | 0.015 | 0.014 | 0.653 | 0.029 | 0.305 |
| 6 | 78 | alveolar fibroblast | 0.872 | True | 0.103 | 0.897 | 0.0164 | 1.1881 | 0.0152 | 0.064 | 0.0 | 0.013 | 0.0 | 0.295 | 0.051 | 0.872 | 0.628 | 0.013 | 0.0 | 0.154 | 0.013 | 0.551 |

