# Trial D5b: organoids under the corrected annotation rule, beside D5

Primary resolution 0.3; states per resolution: 0.3: AT1, DATP, cAT2, hAT2; 0.5: AT1, DATP, hAT2; 1.0: AT1, DATP, cAT2, hAT2

Readings: {'G1': False, 'G2': True, 'G3': False, 'G3_DATP_direction_only': True, 'G4': True, 'G5': {'early': True, 'late': False, 'datp': True, 'set_mean_ratios': {'AT1_early': 0.964, 'AT1_late': 0.723, 'DATP_figure7': 3.655}}, 'G5_n_at1': {'IL-1beta': 279, 'control': 202}, 'G6': {'G3_in_window': False, 'G3_DATP_direction_in_window': True, 'G5_in_window': {'early': True, 'late': True, 'datp': True, 'set_mean_ratios': {'AT1_early': 0.94, 'AT1_late': 0.693, 'DATP_figure7': 3.544}}}, 'paper_epithelial_counts': {'control': 1286, 'IL-1beta': 2584}, 'this_run_epithelial_counts': {'IL-1beta': 2423, 'control': 1193}}

## Contaminant clusters removed (cells per treatment)

| treatment | state | n |
|---|---|---|
| IL-1beta | contaminant_mesenchyme | 276 |
| control | contaminant_mesenchyme | 675 |

## Composition (fractions of epithelial cells)

| treatment | AT1_n | DATP_n | cAT2_n | hAT2_n | AT1_frac | DATP_frac | cAT2_frac | hAT2_frac | epithelial_cells |
|---|---|---|---|---|---|---|---|---|---|
| IL-1beta | 279.0 | 2105.0 | 37.0 | 2.0 | 0.1151 | 0.8688 | 0.0153 | 0.0008 | 2423.0 |
| control | 202.0 | 16.0 | 42.0 | 933.0 | 0.1693 | 0.0134 | 0.0352 | 0.7821 | 1193.0 |

## Composition inside the 2,500 to 4,500 gene window

| treatment | AT1_n | DATP_n | cAT2_n | hAT2_n | AT1_frac | DATP_frac | cAT2_frac | hAT2_frac | epithelial_cells |
|---|---|---|---|---|---|---|---|---|---|
| IL-1beta | 198.0 | 1692.0 | 25.0 | 2.0 | 0.1033 | 0.8826 | 0.013 | 0.001 | 1917.0 |
| control | 107.0 | 10.0 | 14.0 | 722.0 | 0.1254 | 0.0117 | 0.0164 | 0.8464 | 853.0 |

## AT1 cells: marker detection, IL-1beta over control

| set | gene | control | n_control | IL-1beta | n_IL-1beta | ratio_il1b_over_control |
|---|---|---|---|---|---|---|
| AT1_early | Lmo7 | 0.99 | 202 | 0.961 | 279 | 0.97 |
| AT1_early | Pdpn | 0.871 | 202 | 0.867 | 279 | 0.996 |
| AT1_early | Hopx | 0.95 | 202 | 0.882 | 279 | 0.928 |
| AT1_late | Aqp5 | 1.0 | 202 | 0.996 | 279 | 0.996 |
| AT1_late | Vegfa | 0.985 | 202 | 0.9 | 279 | 0.913 |
| AT1_late | Cav1 | 0.797 | 202 | 0.462 | 279 | 0.58 |
| AT1_late | Spock2 | 0.683 | 202 | 0.276 | 279 | 0.404 |
| DATP_figure7 | Cldn4 | 0.213 | 202 | 0.803 | 279 | 3.772 |
| DATP_figure7 | AW112010 | 0.173 | 202 | 0.925 | 279 | 5.337 |
| DATP_figure7 | Lhfp | 0.262 | 202 | 0.487 | 279 | 1.858 |

## Cluster annotation at the primary resolution

| cluster | n_cells | Sftpc | hAT2_canonical | AT2_identity | AT2_lipid | pAT2_inflammatory | cAT2 | DATP | AT1_canonical | AT1_early | AT1_late | DATP_figure7 | p53 | arrest | hypoxia | hypoxia_without_Ndrg1 | ifng_response | glycolysis | responder | ciliated | mesenchyme | immune | endothelium | club_or_airway | state | rule | AT1_reference_detection | sftpc_informative |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 586 | 1.0 | 0.982 | 0.549 | 0.788 | 0.574 | 0.09 | 0.219 | 0.336 | 0.517 | 0.37 | 0.511 | 0.529 | 0.904 | 0.476 | 0.739 | 0.522 | 0.463 | 0.444 | 0.033 | 0.701 | 0.008 | 0.01 | 0.188 | contaminant_mesenchyme | 1b | 0.794 | False |
| 1 | 935 | 1.0 | 1.0 | 0.944 | 0.856 | 0.434 | 0.039 | 0.215 | 0.119 | 0.258 | 0.446 | 0.062 | 0.42 | 0.421 | 0.325 | 0.624 | 0.278 | 0.331 | 0.027 | 0.114 | 0.175 | 0.008 | 0.007 | 0.287 | hAT2 | 5 (reference) | 0.794 | False |
| 2 | 481 | 1.0 | 0.985 | 0.469 | 0.731 | 0.474 | 0.031 | 0.481 | 0.794 | 0.918 | 0.746 | 0.519 | 0.427 | 0.571 | 0.339 | 0.601 | 0.524 | 0.466 | 0.195 | 0.148 | 0.16 | 0.004 | 0.008 | 0.141 | AT1 | 3b | 0.794 | False |
| 3 | 365 | 1.0 | 0.997 | 0.4 | 0.771 | 0.624 | 0.103 | 0.28 | 0.426 | 0.432 | 0.405 | 0.385 | 0.527 | 0.925 | 0.648 | 0.871 | 0.425 | 0.658 | 0.512 | 0.04 | 0.636 | 0.016 | 0.024 | 0.133 | contaminant_mesenchyme | 1b | 0.794 | False |
| 4 | 79 | 1.0 | 1.0 | 0.781 | 0.772 | 0.661 | 0.523 | 0.423 | 0.266 | 0.392 | 0.497 | 0.283 | 0.604 | 0.519 | 0.519 | 0.772 | 0.487 | 0.46 | 0.101 | 0.253 | 0.283 | 0.022 | 0.008 | 0.333 | cAT2 | 2 | 0.794 | False |
| 5 | 2121 | 1.0 | 1.0 | 0.757 | 0.652 | 0.621 | 0.023 | 0.422 | 0.186 | 0.299 | 0.354 | 0.291 | 0.374 | 0.403 | 0.464 | 0.726 | 0.339 | 0.442 | 0.033 | 0.055 | 0.142 | 0.005 | 0.007 | 0.237 | DATP | 4 | 0.794 | False |
