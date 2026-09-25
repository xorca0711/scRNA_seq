# Trial C3 output: is Areg a property of the DATP-like state (Gate 2a)

This is an expression and state-specificity result, not a communication analysis. CellChat
is R-only and unavailable here; no ligand-receptor probability is computed and nothing below
establishes that a fibroblast receives the signal.

**T1, state specificity:** Areg higher in DATP-like than in AT2 cells in all mutant RFP libraries.

| library | areg_DATP_like | areg_AT2 | higher_in_DATP |
|---|---|---|---|
| Expt1_2wRFPr1 | 2.9211 | 1.1283 | True |
| Expt1_2wRFPr2 | 2.8343 | 0.6862 | True |
| Expt1_4dRFPr1 | 2.5535 | 0.5353 | True |
| Expt1_4dRFPr2 | 2.1185 | 0.4935 | True |

**T2, ligand ranking within the DATP-like state.** What ranks second is the part the paper
does not show.

| library | ranking | top | areg_is_top | mean_Areg | mean_Ereg | mean_Hbegf | mean_Tgfa |
|---|---|---|---|---|---|---|---|
| Expt1_2wRFPr1 | Areg > Hbegf > Ereg > Tgfa | Areg | True | 2.9211 | 0.9357 | 1.3555 | 0.2448 |
| Expt1_2wRFPr2 | Areg > Hbegf > Ereg > Tgfa | Areg | True | 2.8343 | 0.4922 | 1.2782 | 0.1 |
| Expt1_4dRFPr1 | Areg > Hbegf > Ereg > Tgfa | Areg | True | 2.5535 | 0.4417 | 0.5697 | 0.1824 |
| Expt1_4dRFPr2 | Areg > Hbegf > Ereg > Tgfa | Areg | True | 2.1185 | 0.5099 | 0.6746 | 0.2398 |

**T3, the within-animal wild-type control.**

| timepoint | median_pct_DATP_RFP | median_pct_DATP_YFP | higher_in_RFP |
|---|---|---|---|
| 4 days | 14.66 | 0.08 | True |
| 2 weeks | 25.84 | 1.14 | True |

## Replicate mixing, and the batch decision it drove

| arm | same_library_kNN_enrichment |
|---|---|
| 2w_RFP | 1.7387 |
| 2w_YFP | 1.5314 |
| 4d_RFP | 1.6197 |
| 4d_YFP | 1.8029 |

Failure threshold 2.0; correction applied: False.

## Clusters and their calls

| cluster | n_cells | call | mode_fraction | confident | n_libraries | max_library_fraction |
|---|---|---|---|---|---|---|
| 0 | 11031 | AT2 | 0.999 | True | 10 | 0.352 |
| 1 | 3045 | AT2 | 0.983 | True | 10 | 0.301 |
| 2 | 7359 | AT2 | 0.998 | True | 10 | 0.369 |
| 3 | 3104 | AT2 | 0.942 | True | 10 | 0.443 |
| 4 | 421 | AT1_like | 0.487 | False | 10 | 0.304 |
| 5 | 2255 | AT2 | 0.867 | True | 10 | 0.358 |
| 6 | 2522 | DATP_like | 0.697 | True | 8 | 0.42 |
| 7 | 138 | AT2 | 0.667 | True | 8 | 0.326 |
| 8 | 421 | AT2 | 0.622 | True | 10 | 0.173 |
| 9 | 967 | Cd177_positive | 0.384 | False | 10 | 0.303 |
| 10 | 173 | AT2 | 0.618 | True | 10 | 0.243 |
| 11 | 1422 | Cd177_positive | 0.546 | True | 10 | 0.421 |
| 12 | 175 | AT1_like | 0.446 | False | 10 | 0.246 |
| 13 | 184 | AT2 | 0.538 | True | 10 | 0.451 |

## QC per library

| library | arm | replicate | clone | timepoint | barcodes | kept_after_qc | median_genes | doublet_call_method | doublets_removed | cells_analysed |
|---|---|---|---|---|---|---|---|---|---|---|
| Expt1_ConfettiRFP | Confetti_RFP | 1 | wild type | Confetti baseline | 5366 | 4309 | 1286 | top 3.45% of scores (automatic rejected: it called 0.00%) | 203 | 4106 |
| Expt1_ConfettiYFP | Confetti_YFP | 1 | wild type | Confetti baseline | 6213 | 5174 | 1143 | top 4.14% of scores (automatic rejected: it called 0.00%) | 294 | 4880 |
| Expt1_4dRFPr1 | 4d_RFP | 1 | KrasG12D | 4 days | 4811 | 3968 | 1408 | top 3.17% of scores (automatic rejected: it called 0.00%) | 171 | 3797 |
| Expt1_4dRFPr2 | 4d_RFP | 2 | KrasG12D | 4 days | 2144 | 1931 | 1456 | scrublet automatic threshold | 10 | 1921 |
| Expt1_4dYFPr1 | 4d_YFP | 1 | wild type | 4 days | 2584 | 1984 | 1633 | top 1.59% of scores (automatic rejected: it called 0.25%) | 58 | 1926 |
| Expt1_4dYFPr2 | 4d_YFP | 2 | wild type | 4 days | 3480 | 3070 | 1705 | top 2.46% of scores (automatic rejected: it called 0.03%) | 122 | 2948 |
| Expt1_2wRFPr1 | 2w_RFP | 1 | KrasG12D | 2 weeks | 4859 | 4087 | 3212 | top 3.27% of scores (automatic rejected: it called 0.05%) | 135 | 3952 |
| Expt1_2wRFPr2 | 2w_RFP | 2 | KrasG12D | 2 weeks | 4429 | 3534 | 2793 | top 2.83% of scores (automatic rejected: it called 0.25%) | 119 | 3415 |
| Expt1_2wYFPr1 | 2w_YFP | 1 | wild type | 2 weeks | 5687 | 4498 | 1639 | top 3.60% of scores (automatic rejected: it called 0.00%) | 234 | 4264 |
| Expt1_2wYFPr2 | 2w_YFP | 2 | wild type | 2 weeks | 2439 | 2019 | 1642 | scrublet automatic threshold | 11 | 2008 |

