# Trial C1c output: are fibrotic and inflammatory fibroblasts separate cells here

The paper makes fibrotic and inflammatory fibroblasts distinct populations, with the
inflammatory cells appearing from 4 weeks, sitting at the tumour periphery and lacking Tnc.
These libraries are a single 2-week time point, so the timing half of that claim cannot be
checked here at all. What can be checked is whether the two marker programmes mark the same
cells at 2 weeks.

| group | n_cells | pct_Tnc | pct_Lcn2_or_Saa3 | pct_double_positive | pct_expected_if_independent | observed_over_expected | reading |
|---|---|---|---|---|---|---|---|
| Red2Kras, C1b private clusters | 1854 | 27.78 | 29.5 | 9.65 | 8.2 | 1.178 | consistent with separate cells |
| Red2Kras, all other clusters | 4085 | 17.01 | 4.5 | 0.73 | 0.77 | 0.958 | consistent with separate cells |
| Confetti, all clusters | 5751 | 11.46 | 1.15 | 0.14 | 0.13 | 1.058 | consistent with separate cells |

## By cluster (Red2Kras cells only, clusters with at least 50)

| cluster | n_Red2Kras_cells | is_C1b_private | pct_Tnc | pct_Lcn2_or_Saa3 | pct_double_positive | observed_over_expected |
|---|---|---|---|---|---|---|
| 0 | 708 | False | 12.85 | 4.52 | 0.28 | 0.486 |
| 1 | 449 | False | 34.3 | 3.12 | 1.34 | 1.25 |
| 2 | 860 | False | 1.16 | 4.53 | 0.12 | 2.205 |
| 3 | 367 | False | 2.18 | 9.54 | 0.0 | 0.0 |
| 4 | 327 | True | 21.41 | 52.29 | 8.26 | 0.738 |
| 5 | 370 | False | 53.51 | 4.05 | 3.24 | 1.495 |
| 6 | 766 | False | 26.24 | 1.96 | 0.65 | 1.27 |
| 7 | 104 | False | 14.42 | 2.88 | 1.92 | 4.622 |
| 8 | 59 | False | 1.69 | 1.69 | 0.0 | 0.0 |
| 9 | 194 | False | 3.61 | 3.09 | 0.0 | 0.0 |
| 10 | 274 | True | 0.36 | 4.74 | 0.0 | 0.0 |
| 11 | 182 | True | 14.29 | 44.51 | 2.75 | 0.432 |
| 14 | 933 | True | 40.09 | 26.37 | 13.72 | 1.298 |
| 15 | 146 | False | 2.05 | 8.9 | 0.0 | 0.0 |
| 16 | 138 | True | 31.88 | 26.09 | 13.77 | 1.655 |

One library, three pooled mice: no P value, and the ratio is a description of one library.

