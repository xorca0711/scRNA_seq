# Trial D0: what the Choi-2020 deposit contains

Raw 10x barcode whitelists rather than called cells: 6 of 8 libraries (PBS_AT2_Tomato, PBS_AT2_nonTomato, Day14_AT2_Tomato, Day14_AT2_nonTomato, Day28_AT2_Tomato, Day28_AT2_nonTomato). Already filtered to called cells: Control_Organoids, Il1b_Organoids.
Distinct gene spaces across the eight libraries: 1.
Non-gene features found: 0.
Any contrast with within-group replication: False.

## Accessions not parsed, and why

- **GSE144598**: ATAC-seq deposited as bigwig coverage tracks only; no peaks and no reads, so the epigenetic claim cannot be re-derived without going to SRA
- **GSE144553**: SuperSeries over the other three; holds no files of its own

## Library inventory

| accession | library | n_features | n_barcodes | barcodes_with_any_count | barcodes_over_500_counts | median_genes_per_barcode_over_500 | is_raw_whitelist | non_gene_features |
|---|---|---|---|---|---|---|---|---|
| GSE145031 | PBS_AT2_Tomato | 27998 | 737280 | 465139 | 6877 | 1061.0 | True | 0 |
| GSE145031 | PBS_AT2_nonTomato | 27998 | 737280 | 493160 | 12592 | 1196.5 | True | 0 |
| GSE145031 | Day14_AT2_Tomato | 27998 | 737280 | 525637 | 5336 | 1825.0 | True | 0 |
| GSE145031 | Day14_AT2_nonTomato | 27998 | 737280 | 480656 | 9213 | 891.0 | True | 0 |
| GSE145031 | Day28_AT2_Tomato | 27998 | 737280 | 512172 | 5693 | 1689.0 | True | 0 |
| GSE145031 | Day28_AT2_nonTomato | 27998 | 737280 | 501502 | 11177 | 968.0 | True | 0 |
| GSE144468 | Control_Organoids | 27998 | 2101 | 2101 | 2101 | 3598.0 | False | 0 |
| GSE144468 | Il1b_Organoids | 27998 | 3066 | 3066 | 3066 | 3193.0 | False | 0 |

## Gene spaces

| gene_space | n_features | n_libraries | libraries |
|---|---|---|---|
| 1 | 27998 | 8 | PBS_AT2_Tomato, PBS_AT2_nonTomato, Day14_AT2_Tomato, Day14_AT2_nonTomato, Day28_AT2_Tomato, Day28_AT2_nonTomato, Control_Organoids, Il1b_Organoids |

## The replicate rule, applied

| accession | contrast | held_constant | cells_with_more_than_one_library | has_within_group_replication |
|---|---|---|---|---|
| GSE145031 | timepoint | sort | 0 | False |
| GSE145031 | sort | timepoint | 0 | False |
| GSE144468 | treatment | none | 0 | False |

## Marker sets with any gene absent from a library

None. Every marker set in the extract is fully present in every library.
