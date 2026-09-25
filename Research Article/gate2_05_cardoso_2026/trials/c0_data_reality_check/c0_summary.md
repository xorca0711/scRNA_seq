# Trial C0 output: what the Cardoso 2026 deposits contain (Gate 0)

All 30 deposited libraries of the five accessions are public and readable: 123,807 barcodes before any quality control. Every species call agrees with the accession's stated organism. 3 distinct gene spaces are present, so integration across accessions needs an explicit intersection on Ensembl gene ID.

## What each accession is, and what it can carry

| accession | libraries | barcodes | median_genes_per_cell_range | mice_pooled_per_library | independent_animals_per_genotype | distinct_genotypes_in_characteristics | role |
|---|---|---|---|---|---|---|---|
| GSE316241 | 2 | 13226 | 1824 to 2111 | 3 | 1 library, 3 mice pooled: no within-genotype replication | 2 | mesenchymal cells, Confetti and Red2Kras lungs (paper Fig. 1b-e) |
| GSE316243 | 2 | 7836 | 677 to 1037 | 3 | 1 library, 3 mice pooled: no within-genotype replication | 0 | immune and stromal cells, Confetti and Red2Kras lungs (paper Fig. 1f-h) |
| GSE316244 | 4 | 33756 | 2615 to 3874 | 3 | 1 library, 3 mice pooled: no within-genotype replication | 2 | Areg-flox arm: niche and RFP+ mutant epithelium (paper Fig. 4d-m) |
| GSE310335 | 2 | 9408 | 3757 to 5576 | not stated | see the deposited design | 0 | human KRAS-G12D alveolar organoids (paper Fig. 5f-h) |
| GSE247505 | 20 | 59581 | 1136 to 3212 | not stated | see the deposited design | 4 | England 2025: RFP+ lineage-labelled mutant epithelium, the CellChat partner |

The entry marked *no within-genotype replication* is the binding constraint of this deposit: each mouse library is one library per genotype built from three pooled mice, so a genotype difference in these data is a difference between two libraries. It can be described; it cannot be tested, and no P value is admissible on it.

## Gene spaces

| gene_space | n_libraries | libraries |
|---|---|---|
| 1 | 6 | GSE316241:GSM9447763_Expt1_Confetti_mesenchyme; GSE316241:GSM9447764_Expt1_Red2Kras_mesenchyme; GSE316244:GSM9447779_Expt3_Het_niche; GSE316244:GSM9447780_Expt3_Hom_niche; GSE316244:GSM9447781_Expt3_Het_RFP; GSE316244:GSM9447782_Expt3_Hom_RFP |
| 2 | 22 | GSE316243:GSM9447777_Expt2_Confetti_niche; GSE316243:GSM9447778_Expt2_Red2Kras_niche; GSE247505:GSM7890829_Expt1_ConfettiRFP; GSE247505:GSM7890830_Expt1_ConfettiYFP; GSE247505:GSM7890831_Expt1_4dRFPr1; GSE247505:GSM7890832_Expt1_4dRFPr2; GSE247505:GSM7890833_Expt1_4dYFPr1; GSE247505:GSM7890834_Expt1_4dYFPr2; GSE247505:GSM7890835_Expt1_2wRFPr1; GSE247505:GSM7890836_Expt1_2wRFPr2; GSE247505:GSM7890837_Expt1_2wYFPr1; GSE247505:GSM7890838_Expt1_2wYFPr2; GSE247505:GSM7890839_Expt2_ConfettiRFPr1; GSE247505:GSM7890840_Expt2_ConfettiRFPr2; GSE247505:GSM7890841_Expt2_2wRFP_hetr1; GSE247505:GSM7890842_Expt2_2wRFP_hetr2; GSE247505:GSM7890843_Expt2_2wRFP_homr1; GSE247505:GSM7890844_Expt2_2wRFP_homr2; GSE247505:GSM7890845_Expt2_12wRFP_hetr1; GSE247505:GSM7890846_Expt2_12wRFP_hetr2; GSE247505:GSM7890847_Expt2_12wRFP_homr1; GSE247505:GSM7890848_Expt2_12wRFP_homr2 |
| 3 | 2 | GSE310335:GSE310335_SIGAE; GSE310335:GSE310335_SIGAF |

## Per-library quality, measured and not applied

| accession | library | barcodes | median_total_counts | median_genes | median_pct_mt | cells_over_10pct_mt | non_gene_features |
|---|---|---|---|---|---|---|---|
| GSE316241 | GSM9447763_Expt1_Confetti_mesenchyme | 6435 | 4529 | 2111 | 1.18 | 46 | none |
| GSE316241 | GSM9447764_Expt1_Red2Kras_mesenchyme | 6791 | 3943 | 1824 | 1.4 | 85 | none |
| GSE316243 | GSM9447777_Expt2_Confetti_niche | 4498 | 2042 | 677 | 1.88 | 12 | none |
| GSE316243 | GSM9447778_Expt2_Red2Kras_niche | 3338 | 5885 | 1037 | 1.71 | 49 | none |
| GSE316244 | GSM9447779_Expt3_Het_niche | 8541 | 7056 | 2615 | 1.83 | 306 | BSD |
| GSE316244 | GSM9447780_Expt3_Hom_niche | 10891 | 7509 | 2750 | 1.89 | 226 | BSD |
| GSE316244 | GSM9447781_Expt3_Het_RFP | 7395 | 13367 | 3874 | 1.86 | 300 | BSD |
| GSE316244 | GSM9447782_Expt3_Hom_RFP | 6929 | 12277 | 3545 | 2.76 | 920 | BSD |
| GSE310335 | GSE310335_SIGAE | 6190 | 16165 | 3757 | 7.83 | 1570 | none |
| GSE310335 | GSE310335_SIGAF | 3218 | 30586 | 5576 | 8.94 | 1231 | none |
| GSE247505 | GSM7890829_Expt1_ConfettiRFP | 5366 | 3581 | 1286 | 4.47 | 723 | none |
| GSE247505 | GSM7890830_Expt1_ConfettiYFP | 6213 | 2864 | 1143 | 3.9 | 443 | none |
| GSE247505 | GSM7890831_Expt1_4dRFPr1 | 4811 | 4306 | 1408 | 3.1 | 351 | none |
| GSE247505 | GSM7890832_Expt1_4dRFPr2 | 2144 | 4623 | 1456 | 1.39 | 7 | none |
| GSE247505 | GSM7890833_Expt1_4dYFPr1 | 2584 | 5154 | 1633 | 3.76 | 360 | none |
| GSE247505 | GSM7890834_Expt1_4dYFPr2 | 3480 | 5932 | 1705 | 2.12 | 30 | none |
| GSE247505 | GSM7890835_Expt1_2wRFPr1 | 4859 | 14451 | 3212 | 2.44 | 266 | none |
| GSE247505 | GSM7890836_Expt1_2wRFPr2 | 4429 | 11484 | 2793 | 2.04 | 304 | none |
| GSE247505 | GSM7890837_Expt1_2wYFPr1 | 5687 | 5027 | 1639 | 2.72 | 377 | none |
| GSE247505 | GSM7890838_Expt1_2wYFPr2 | 2439 | 5642 | 1642 | 1.99 | 78 | none |
| GSE247505 | GSM7890839_Expt2_ConfettiRFPr1 | 1987 | 4537 | 1597 | 4.06 | 244 | none |
| GSE247505 | GSM7890840_Expt2_ConfettiRFPr2 | 2780 | 2890 | 1136 | 3.29 | 174 | none |
| GSE247505 | GSM7890841_Expt2_2wRFP_hetr1 | 983 | 8347 | 2327 | 3.07 | 78 | none |
| GSE247505 | GSM7890842_Expt2_2wRFP_hetr2 | 973 | 8417 | 2330 | 2.75 | 108 | none |
| GSE247505 | GSM7890843_Expt2_2wRFP_homr1 | 1015 | 7200 | 2176 | 2.98 | 92 | none |
| GSE247505 | GSM7890844_Expt2_2wRFP_homr2 | 784 | 7037 | 1981 | 3.69 | 112 | none |
| GSE247505 | GSM7890845_Expt2_12wRFP_hetr1 | 2648 | 2693 | 1238 | 3.11 | 101 | none |
| GSE247505 | GSM7890846_Expt2_12wRFP_hetr2 | 2146 | 4455 | 1581 | 1.78 | 57 | none |
| GSE247505 | GSM7890847_Expt2_12wRFP_homr1 | 2255 | 4613 | 1727 | 3.11 | 88 | none |
| GSE247505 | GSM7890848_Expt2_12wRFP_homr2 | 1998 | 4523 | 1485 | 2.94 | 107 | none |

## Non-gene features

One feature in the whole deposit is not an Ensembl gene: `BSD` (blasticidin-S deaminase, a selection marker carried by the reporter construct), present in GSE316244 only. It is removed from the expression matrix by the frozen rule and carried per cell. Its detection tracks the sort:

| library | cells_detected | fraction_of_cells | total_counts |
|---|---|---|---|
| GSM9447779_Expt3_Het_niche | 426 | 0.0499 | 545 |
| GSM9447780_Expt3_Hom_niche | 311 | 0.0286 | 364 |
| GSM9447781_Expt3_Het_RFP | 2875 | 0.3888 | 4932 |
| GSM9447782_Expt3_Hom_RFP | 2757 | 0.3979 | 4727 |

## Are the populations the paper's claims rest on detectable at all

Fraction of cells with a non-zero count, mesenchymal series. A detection fraction is a presence measure: it is not abundance and it is not a cell-type call.

| panel | gene | Confetti_mesenchyme | Red2Kras_mesenchyme |
|---|---|---|---|
| EGF axis (ED Fig. 5f-h) | Areg | 0.008 | 0.062 |
| EGF axis (ED Fig. 5f-h) | Egfr | 0.36 | 0.322 |
| EGF axis (ED Fig. 5f-h) | Ereg | 0.011 | 0.072 |
| EGF axis (ED Fig. 5f-h) | Hbegf | 0.038 | 0.078 |
| EGF axis (ED Fig. 5f-h) | Tgfa | 0.022 | 0.031 |
| alveolar fibroblast (Fig. 1d) | Col13a1 | 0.424 | 0.213 |
| alveolar fibroblast (Fig. 1d) | Npnt | 0.66 | 0.417 |
| alveolar fibroblast (Fig. 1d) | Pdgfra | 0.579 | 0.415 |
| alveolar fibroblast (Fig. 1d) | Scube2 | 0.353 | 0.131 |
| alveolar fibroblast (Fig. 1d) | Tcf21 | 0.594 | 0.387 |
| fibroblast-to-epithelium (ED Fig. 5i-j) | Igf1 | 0.438 | 0.394 |
| fibroblast-to-epithelium (ED Fig. 5i-j) | Spp1 | 0.025 | 0.228 |
| fibroblast-to-epithelium (ED Fig. 5i-j) | Wnt5a | 0.169 | 0.195 |
| inflammatory fibroblast (ED Fig. 4l-m) | Cxcl12 | 0.37 | 0.437 |
| inflammatory fibroblast (ED Fig. 4l-m) | Lcn2 | 0.006 | 0.076 |
| inflammatory fibroblast (ED Fig. 4l-m) | Saa3 | 0.007 | 0.069 |
| inflammatory fibroblast (ED Fig. 4l-m) | Sfrp1 | 0.202 | 0.276 |
| reprogrammed AM (Fig. 1h) | Cdh1 | 0.007 | 0.041 |
| reprogrammed AM (Fig. 1h) | Ch25h | 0.027 | 0.066 |
| reprogrammed AM (Fig. 1h) | Cxcl2 | 0.023 | 0.103 |
| reprogrammed AM (Fig. 1h) | Ear6 | 0.0 | 0.019 |
| reprogrammed AM (Fig. 1h) | Fstl1 | 0.75 | 0.673 |
| reprogrammed AM (Fig. 1h) | Msr1 | 0.001 | 0.025 |
| reprogrammed/fibrotic fibroblast (Fig. 1d) | Acta2 | 0.489 | 0.527 |
| reprogrammed/fibrotic fibroblast (Fig. 1d) | Fst | 0.043 | 0.124 |
| reprogrammed/fibrotic fibroblast (Fig. 1d) | Pdgfrb | 0.171 | 0.204 |
| reprogrammed/fibrotic fibroblast (Fig. 1d) | Runx1 | 0.22 | 0.39 |
| reprogrammed/fibrotic fibroblast (Fig. 1d) | Runx2 | 0.03 | 0.111 |
| reprogrammed/fibrotic fibroblast (Fig. 1d) | Tnc | 0.114 | 0.2 |

Genes of the paper's own panels that are absent from a deposited reference: Siglec5.

