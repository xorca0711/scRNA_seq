# Trial C2 output: the Areg deletion arm (Gate 2b)

Design constraint from trial C0: one library per genotype per sort, three mice pooled per
library, one time point. Directions are read; nothing is tested; no P value is computed.

Pre-registered directions met: **3 of 5**.

| pre-registered direction | met |
|---|---|
| reprogrammed_fibroblast_share_falls_in_hom | False |
| DATP_like_share_falls_in_hom | True |
| Cd177_positive_share_falls_in_hom | False |
| AT2_share_rises_in_hom | True |
| AM_inflammatory_falls_and_MHCII_rises_in_hom | True |

## QC per library

| library | genotype | sort | barcodes | kept_after_qc | median_genes | doublet_call_method | doublets_removed | cells_analysed | BSD_detected_fraction |
|---|---|---|---|---|---|---|---|---|---|
| Expt3_Het_niche | Areg-flox/+ | niche | 8541 | 7367 | 2615 | top 5.89% of scores (automatic rejected: it called 0.00%) | 445 | 6922 | 0.0506 |
| Expt3_Hom_niche | Areg-flox/flox | niche | 10891 | 9844 | 2750 | scrublet automatic threshold | 356 | 9488 | 0.0284 |
| Expt3_Het_RFP | Areg-flox/+ | RFP+ epithelium | 7395 | 6531 | 3874 | top 5.22% of scores (automatic rejected: it called 0.11%) | 396 | 6135 | 0.4137 |
| Expt3_Hom_RFP | Areg-flox/flox | RFP+ epithelium | 6929 | 5489 | 3545 | top 4.39% of scores (automatic rejected: it called 0.07%) | 274 | 5215 | 0.4594 |
| Expt1_Confetti_mesenchyme | Confetti | mesenchyme | 6435 | 6057 | 2111 | top 4.85% of scores (automatic rejected: it called 0.07%) | 306 | 5751 |  |
| Expt1_Red2Kras_mesenchyme | Red2Kras | mesenchyme | 6791 | 6260 | 1824 | top 5.01% of scores (automatic rejected: it called 0.00%) | 321 | 5939 |  |

## Fibroblast composition by genotype

| compartment | genotype | call | n_cells | pct_of_library |
|---|---|---|---|---|
| fibroblasts | Areg-flox/+ | adventitial fibroblast | 1150 | 50.07 |
| fibroblasts | Areg-flox/+ | alveolar fibroblast | 637 | 27.73 |
| fibroblasts | Areg-flox/flox | adventitial fibroblast | 1745 | 47.82 |
| fibroblasts | Areg-flox/flox | alveolar fibroblast | 1870 | 51.25 |

## RFP+ epithelial composition by genotype

| compartment | genotype | call | n_cells | pct_of_library |
|---|---|---|---|---|
| RFP+ epithelium | Areg-flox/+ | AT1_like | 5 | 0.08 |
| RFP+ epithelium | Areg-flox/+ | AT2 | 1794 | 29.24 |
| RFP+ epithelium | Areg-flox/+ | DATP_like | 2860 | 46.62 |
| RFP+ epithelium | Areg-flox/+ | cycling | 54 | 0.88 |
| RFP+ epithelium | Areg-flox/flox | AT1_like | 6 | 0.12 |
| RFP+ epithelium | Areg-flox/flox | AT2 | 3081 | 59.08 |
| RFP+ epithelium | Areg-flox/flox | DATP_like | 1147 | 21.99 |
| RFP+ epithelium | Areg-flox/flox | cycling | 11 | 0.21 |

## What survives ligand deletion (Part B, secondary)

Confounded by CellRanger version, sort composition and tamoxifen dose; see the run record.

| gene | absent_from_reference | Confetti | Red2Kras | Areg-flox/+ | Areg-flox/flox | verdict |
|---|---|---|---|---|---|---|
| Acta2 | False | 0.3506 | 0.3074 | 0.2786 | 0.2149 | collapses |
| Col13a1 | False | 0.5813 | 0.2715 | 0.4802 | 0.6544 | Areg-independent |
| Fst | False | 0.0518 | 0.1271 | 0.3361 | 0.1425 | collapses |
| Npnt | False | 0.7384 | 0.3619 | 0.5364 | 0.6901 | unchanged or higher |
| Pdgfra | False | 0.8325 | 0.6567 | 0.9168 | 0.9408 | Areg-independent |
| Pdgfrb | False | 0.1271 | 0.165 | 0.4519 | 0.3886 | Areg-independent |
| Runx1 | False | 0.146 | 0.177 | 0.5599 | 0.5434 | Areg-independent |
| Runx2 | False | 0.0303 | 0.0326 | 0.1754 | 0.0691 | collapses |
| Scube2 | False | 0.5285 | 0.2096 | 0.3361 | 0.5056 | unchanged or higher |
| Tcf21 | False | 0.7882 | 0.6221 | 0.5076 | 0.7257 | unchanged or higher |
| Tnc | False | 0.1118 | 0.1138 | 0.1828 | 0.1381 | collapses |

