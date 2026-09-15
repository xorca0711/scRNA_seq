# Trial D7: attack claims

## A1

{'not_computable': 'pAT2 or hAT2 below floor', 'n_pAT2': 0, 'n_hAT2': 9546}

## A2

{'DATP_mean_scrublet': 0.05227589096976077, 'DATP_double_marker_fraction': 0.2375, 'doublet_like': False}

## A3

not attempted: needs the published dissociation gene list (van den Brink 2017, 10.1038/nmeth.4437), not on disk

## A4

{'n_cells': 3661, 'krt8_detected': 0.9008467631794591, 'cldn4_detected': 0.16716744059000274, 'double_positive': 578, 'expected_under_independence': 551.318219065829, 'ratio': 1.0483963344788088, 'permutation_p_within_depth_deciles': 0.004975124378109453, 'null_mean': 558.7, 'reading': 'no evidence of a discrete co-expressing population'}

### d7_a1_depth_by_state.csv

| state | total_counts | n_genes | pct_mt | n |
|---|---|---|---|---|
| AT1 | 3906.0 | 1618.0 | 1.8092 | 337.0 |
| DATP | 7679.0 | 2578.0 | 1.6899 | 880.0 |
| cAT2 | 13075.0 | 3093.0 | 2.7556 | 313.0 |
| hAT2 | 6447.5 | 1855.0 | 3.6807 | 9546.0 |

### d7_a2_doublet_by_state.csv

| state | n | mean_scrublet_score | score_gt_0.7_fraction | sftpc_and_2plus_AT1_fraction |
|---|---|---|---|---|
| hAT2 | 9546.0 | 0.0617 | 0.0 | 0.0119 |
| cAT2 | 313.0 | 0.0505 | 0.0 | 0.0927 |
| DATP | 880.0 | 0.0523 | 0.0 | 0.2375 |
| AT1 | 337.0 | 0.0401 | 0.0 | 0.7359 |

