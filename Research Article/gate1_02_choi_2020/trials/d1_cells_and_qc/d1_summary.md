# Trial D1: cells, quality and doublets

Tomato-positive libraries: 12865 barcodes pass the paper's filter and 11683 remain after the MAD rule and Scrublet, against the paper's 12514 captured and 12086 in Figure 1B (their counts precede this repository's doublet removal and follow Cell Ranger's caller, so the numbers are comparable in size, not identical by construction).

## Per library

| library | n_barcodes | paper_rule_cells | floor_rule_cells | floor_only_barcodes | mad_removed | binding_bounds | scrublet_method | scrublet_threshold | doublets_flagged | cells_final | median_counts | median_genes | median_pct_mt |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| PBS_AT2_Tomato | 737280 | 4418 | 6574 | 2156 | 98 | pct_mt:removed_high | scrublet automatic threshold | 0.42870180266350993 | 1 | 4319 | 4775.5 | 1568.0 | 5.417104970280347 |
| PBS_AT2_nonTomato | 737280 | 6963 | 12558 | 5595 | 410 | total_counts:removed_high;pct_mt:removed_high | scrublet automatic threshold | 0.5720757961327718 | 1 | 6552 | 4101.0 | 1977.0 | 3.462798315395414 |
| Day14_AT2_Tomato | 737280 | 4226 | 5070 | 844 | 436 | pct_mt:removed_high | scrublet automatic threshold | 0.38035022602673685 | 0 | 3790 | 7902.0 | 2267.5 | 3.1404077766686136 |
| Day14_AT2_nonTomato | 737280 | 4584 | 8977 | 4393 | 446 | pct_mt:removed_high | scrublet automatic threshold | 0.414629575796172 | 0 | 4138 | 5089.5 | 2237.0 | 3.4124105058679213 |
| Day28_AT2_Tomato | 737280 | 4221 | 5478 | 1257 | 646 | pct_mt:removed_high | scrublet automatic threshold | 0.35982188272189775 | 1 | 3574 | 7543.0 | 2089.0 | 2.6607309245263484 |
| Day28_AT2_nonTomato | 737280 | 5694 | 11120 | 5426 | 439 | pct_mt:removed_high | scrublet automatic threshold | 0.5008792714738662 | 3 | 5252 | 4764.0 | 2126.0 | 3.035377852208499 |
| Control_Organoids | 2101 | 2083 | 2101 | 18 | 207 | n_genes:removed_low;pct_mt:removed_high | scrublet automatic threshold | 0.11744013741510784 | 8 | 1868 | 15197.0 | 3672.0 | 4.102860502674471 |
| Il1b_Organoids | 3066 | 3059 | 3066 | 7 | 358 | n_genes:removed_low;pct_mt:removed_high | scrublet automatic threshold | 0.19935662335447457 | 2 | 2699 | 11970.0 | 3270.0 | 3.939034045922407 |

The organoid libraries were deposited already filtered; the paper's filter, the MAD rule and Scrublet were run on them unchanged. The paper reports 1,286 control and 2,584 IL-1beta epithelial cells after removing an EpCAM-negative stromal cluster, which is trial D5's job, not this one's.
