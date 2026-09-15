# Trial D4: composition and its controls

Readings: {'C1_calling_sensitive': {'PBS_AT2_Tomato': False, 'Day14_AT2_Tomato': True, 'Day28_AT2_Tomato': False}, 'C2_R4_in_window': False, 'C2_R5_in_window': False, 'C3_datp_in_day14_nontomato': np.False_}

## C1 calling sensitivity

| library | floor_only_barcodes | median_counts_extra | datp_like_among_extra | datp_fraction_D2b | datp_fraction_if_added | abs_change |
|---|---|---|---|---|---|---|
| PBS_AT2_Tomato | 2156 | 1218.0 | 0 | 0.0031 | 0.002 | 0.001 |
| Day14_AT2_Tomato | 844 | 1108.0 | 7 | 0.1825 | 0.1498 | 0.0326 |
| Day28_AT2_Tomato | 1257 | 1134.0 | 1 | 0.063 | 0.0453 | 0.0177 |

## C2 composition inside the 1,000 to 2,000 gene window

| library | AT1_n | DATP_n | cAT2_n | hAT2_n | AT1_frac | DATP_frac | cAT2_frac | hAT2_frac | cells_in_window |
|---|---|---|---|---|---|---|---|---|---|
| PBS_AT2_Tomato | 35.0 | 7.0 | 0.0 | 2681.0 | 0.0129 | 0.0026 | 0.0 | 0.9846 | 2723.0 |
| Day14_AT2_Tomato | 39.0 | 162.0 | 27.0 | 1151.0 | 0.0283 | 0.1175 | 0.0196 | 0.8347 | 1379.0 |
| Day28_AT2_Tomato | 112.0 | 44.0 | 4.0 | 1229.0 | 0.0806 | 0.0317 | 0.0029 | 0.8848 | 1389.0 |

## C3 Tomato-negative libraries

| library | cells | alveolar_cells | alveolar_fraction | hAT2_n | cAT2_n | pAT2_n | DATP_n | AT1_n | DATP_cluster_present | DATP_dispersed_fraction | DATP_dispersed_n |
|---|---|---|---|---|---|---|---|---|---|---|---|
| PBS_AT2_nonTomato | 6552 | 21 | 0.0032 | 21 | 0 | 0 | 0 | 0 | False | 0.0 | 0 |
| Day14_AT2_nonTomato | 4138 | 0 | 0.0 | 0 | 0 | 0 | 0 | 0 | False |  | 0 |
| Day28_AT2_nonTomato | 5252 | 708 | 0.1348 | 699 | 0 | 0 | 9 | 0 | True | 0.0071 | 5 |

A transcriptome cannot say whether an unlabelled DATP-like cell came from an AT2 cell the lineage did not label or from another origin; the count is reported, the origin is not.
