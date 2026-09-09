# Trial S4 output: mouse cluster 23

Verdict by the frozen rules: **low-count, ambient-like**. Low-count call True; ambient-like call True.

## QC medians

| group | total_counts | n_genes_by_counts | pct_counts_mt | pct_counts_ribo | pct_counts_in_top_20_genes | doublet_score | n_cells |
|---|---|---|---|---|---|---|---|
| cluster_23 | 1867.5 | 1095.5 | 2.964 | 16.723 | 15.364 | 0.041 | 2912.0 |
| all_other_cells | 6790.0 | 2481.0 | 2.835 | 9.146 | 21.748 | 0.028 | 159263.0 |
| random_comparator | 6774.0 | 2466.5 | 2.826 | 8.972 | 21.652 | 0.028 | 2912.0 |
| dominant_sample_EEM-scRNA-289_outside_23 | 10563.0 | 3454.5 | 2.415 | 10.36 | 17.601 | 0.034 | 1870.0 |

## Lineage-marker detection (raw counts > 0)

| group | Sftpc | Scgb1a1 | Ptprc | Pecam1 | Col1a1 | ge3_of_5 | ge2_of_5 |
|---|---|---|---|---|---|---|---|
| cluster_23 | 0.989 | 0.9948 | 0.1387 | 0.6075 | 0.7263 | 0.8427 | 0.9924 |
| random_comparator | 0.1442 | 0.1724 | 0.137 | 0.4928 | 0.2644 | 0.0446 | 0.2174 |

## Cohort split

- cluster23_cells: 2912
- cluster23_from_tracing_samples: 2280
- cluster23_from_annotated_samples: 632
- cluster23_annotated_cohort_labelled: 14
- cluster23_annotated_cohort_labelled_fraction: 0.0222
- atlas_annotated_cohort_labelled_fraction: 0.8595
- labels_of_the_labelled_cluster23_cells: {'CAP1': 7, 'Plasma_cell': 3, 'B_lymphocyte': 2, 'AT2': 1, 'Adventitial_fibroblast': 1}

## Sample composition of cluster 23

| sample_id | cells_in_cluster23 | fraction_of_cluster23 | fraction_of_that_sample | sample_median_total_counts | cohort |
|---|---|---|---|---|---|
| EEM-scRNA-289 | 2278 | 0.7823 | 0.5492 | 4141.0 | tracing |
| EEM-scRNA-181 | 389 | 0.1336 | 0.0604 | 5419.0 | annotated |
| EEM-scRNA-118 | 205 | 0.0704 | 0.043 | 7071.0 | annotated |
| EEM-scRNA-172 | 25 | 0.0086 | 0.0053 | 8464.0 | annotated |
| EEM-scRNA-113 | 4 | 0.0014 | 0.0012 | 5166.0 | annotated |
| EEM-scRNA-238 | 3 | 0.001 | 0.0006 | 6051.0 | annotated |
| EEM-scRNA-237 | 2 | 0.0007 | 0.0003 | 7005.0 | annotated |
| EEM-scRNA-116 | 2 | 0.0007 | 0.0004 | 9762.0 | annotated |
| EEM-scRNA-117 | 1 | 0.0003 | 0.0002 | 10209.0 | annotated |
| EEM-scRNA-239 | 1 | 0.0003 | 0.0001 | 6922.0 | annotated |
| EEM-scRNA-251 | 1 | 0.0003 | 0.0002 | 5085.0 | tracing |
| EEM-scRNA-292 | 1 | 0.0003 | 0.0003 | 11143.0 | tracing |

