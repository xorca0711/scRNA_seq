# Trial C5 output: figures for the three findings, and the corrections they forced

## Ligand rankings by three notions

| library | n_DATP_like | rank_by_abundance | rank_by_enrichment_over_AT2 | rank_by_detection_ratio | Hbegf_place_abundance | Hbegf_place_enrichment |
|---|---|---|---|---|---|---|
| Expt1_4dRFPr1 | 44 | Areg > Hbegf > Ereg > Tgfa | Areg > Ereg > Hbegf > Tgfa | Ereg > Tgfa > Areg > Hbegf | 2 | 3 |
| Expt1_4dRFPr2 | 541 | Areg > Hbegf > Ereg > Tgfa | Areg > Ereg > Hbegf > Tgfa | Ereg > Tgfa > Areg > Hbegf | 2 | 3 |
| Expt1_2wRFPr1 | 818 | Areg > Hbegf > Ereg > Tgfa | Areg > Ereg > Hbegf > Tgfa | Ereg > Tgfa > Areg > Hbegf | 2 | 3 |
| Expt1_2wRFPr2 | 1058 | Areg > Hbegf > Ereg > Tgfa | Areg > Hbegf > Ereg > Tgfa | Ereg > Areg > Hbegf > Tgfa | 2 | 2 |
| Expt1_2wYFPr1 | 21 | Areg > Hbegf > Ereg > Tgfa | Areg > Ereg > Hbegf > Tgfa | Ereg > Tgfa > Hbegf > Areg | 2 | 3 |
| Expt1_2wYFPr2 | 36 | Areg > Hbegf > Ereg > Tgfa | Areg > Hbegf > Ereg > Tgfa | Ereg > Tgfa > Hbegf > Areg | 2 | 2 |

## DATP-like share per library

| library | arm | clone | timepoint | n_cells | pct_of_library |
|---|---|---|---|---|---|
| Expt1_4dRFPr1 | 4d_RFP | KrasG12D | 4 days | 44 | 1.16 |
| Expt1_4dRFPr2 | 4d_RFP | KrasG12D | 4 days | 541 | 28.16 |
| Expt1_4dYFPr1 | 4d_YFP | wild type | 4 days | 1 | 0.05 |
| Expt1_4dYFPr2 | 4d_YFP | wild type | 4 days | 3 | 0.1 |
| Expt1_2wRFPr1 | 2w_RFP | KrasG12D | 2 weeks | 818 | 20.7 |
| Expt1_2wRFPr2 | 2w_RFP | KrasG12D | 2 weeks | 1058 | 30.98 |
| Expt1_2wYFPr1 | 2w_YFP | wild type | 2 weeks | 21 | 0.49 |
| Expt1_2wYFPr2 | 2w_YFP | wild type | 2 weeks | 36 | 1.79 |

## Fibrotic gene retention after Areg deletion

| gene | flox_plus | flox_flox | c2_verdict | retention_pct | reading |
|---|---|---|---|---|---|
| Runx1 | 0.5599 | 0.5434 | Areg-independent | 97.1 | retained |
| Pdgfrb | 0.4519 | 0.3886 | Areg-independent | 86.0 | retained |
| Acta2 | 0.2786 | 0.2149 | collapses | 77.1 | falls by about a quarter |
| Tnc | 0.1828 | 0.1381 | collapses | 75.5 | falls by about a quarter |
| Fst | 0.3361 | 0.1425 | collapses | 42.4 | falls by more than half |
| Runx2 | 0.1754 | 0.0691 | collapses | 39.4 | falls by more than half |

## Cluster 11 discriminating check

| group | n_cells | pct_Epcam | pct_Cdh1 | pct_Krt8 | pct_Krt18 | pct_Sftpc | pct_Lamp3 | pct_Cldn4 | pct_Sox9 | pct_Itga2 | pct_Ager | pct_Hopx | pct_Areg | pct_Col1a1 | pct_Pdgfra | pct_Krt8_and_Col1a1 | median_doublet_score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cluster 11 | 184 | 78.3 | 71.2 | 90.8 | 92.9 | 56.5 | 43.5 | 60.3 | 47.8 | 62.5 | 54.9 | 45.7 | 88.0 | 21.2 | 10.9 | 19.0 | 0.0615 |
| cluster 12 | 39 | 84.6 | 79.5 | 89.7 | 92.3 | 23.1 | 0.0 | 25.6 | 7.7 | 0.0 | 0.0 | 82.1 | 69.2 | 23.1 | 10.3 | 23.1 | 0.069 |
| cluster 14 | 937 | 0.7 | 0.0 | 0.3 | 2.6 | 8.5 | 0.4 | 0.6 | 4.2 | 3.7 | 1.0 | 4.7 | 3.8 | 89.9 | 69.5 | 0.3 | 0.0366 |
| whole object | 11690 | 2.5 | 2.3 | 6.5 | 8.9 | 6.3 | 0.8 | 1.3 | 4.9 | 4.3 | 1.3 | 10.6 | 3.6 | 80.7 | 51.5 | 3.8 | 0.0351 |

