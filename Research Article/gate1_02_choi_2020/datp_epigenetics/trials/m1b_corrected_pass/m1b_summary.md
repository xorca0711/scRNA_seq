# Trial M1b: the corrected pass, and why it also refused

**Reading: NOT COMPUTABLE in all ten wells.** Not one well reached the
hundred-cell floor, and the reason is arithmetic rather than biological.

M1b repaired M1's two disclosed defects: the airway filter became a
within-cell ratio instead of a presence call, and every label was called
on RNA downsampled to one common depth instead of on raw counts. It also
inherited M1's refutation of the GSE247130 suffix map. What it could not
repair was its own rule R4b.

## The ceiling R4b could not clear

R4b required Cldn4 detection AND Krt8 at or above the 99th percentile of
the uninjured well. A 99th-percentile magnitude cut passes one per cent of
the control well by construction. Cldn4 is detected in five to nine per
cent of an infected well. Their conjunction cannot exceed a few tenths of
a per cent anywhere, and it did not:

|  | well | injured | after_filter | transitional | transitional_pct | U3_max_airway_excess |
|---|---|---|---|---|---|---|
| 0 | wildtype_PBS | False | 5492 | 3 | 0.05 | 0.226 |
| 1 | wildtype_SeV | True | 4578 | 18 | 0.39 | 0.02 |
| 2 | AP1mut_PBS | False | 8134 | 4 | 0.05 | 0.234 |
| 3 | AP1mut_SeV | True | 7444 | 11 | 0.15 | 0.053 |
| 4 | P9_control | False | 10475 | 35 | 0.33 | 0.255 |
| 5 | P9_Cebpa_mutant | False | 10252 | 68 | 0.66 | 0.007 |
| 6 | 7wk_control | False | 6483 | 8 | 0.12 | 0.088 |
| 7 | 7wk_Cebpa_mutant | False | 9016 | 15 | 0.17 | 0.086 |
| 8 | SeV_control | True | 6344 | 2 | 0.03 | 0.107 |
| 9 | SeV_Cebpa_mutant | True | 7851 | 34 | 0.43 | 0.028 |

The source paper reports a twelve per cent transitional fraction for its
infected wildtype library. The largest this rule produced anywhere was
0.66 per cent.

## The calibration run, disclosed

Two technical settings were fixed after a first run that computed labels
only and no accessibility quantity: the airway rule became relative, and
the depth budget moved from the deposit's 10th percentile to the largest
budget retaining 60 per cent of every well. The 10th percentile leaves
560 UMI against a median of about 4,700, at which Cldn4 is essentially
undetectable. The first run's numbers are kept beside this one:

|  | well | injured | airway_dropped | after_filter | transitional |
|---|---|---|---|---|---|
| 0 | wildtype_PBS | False | 610 | 6403 | 0 |
| 1 | wildtype_SeV | True | 3844 | 4068 | 2 |
| 2 | AP1mut_PBS | False | 169 | 10781 | 1 |
| 3 | AP1mut_SeV | True | 2501 | 7490 | 3 |
| 4 | P9_control | False | 630 | 11076 | 28 |
| 5 | P9_Cebpa_mutant | False | 290 | 10644 | 56 |
| 6 | 7wk_control | False | 183 | 6865 | 8 |
| 7 | 7wk_Cebpa_mutant | False | 124 | 9714 | 10 |
| 8 | SeV_control | True | 558 | 8553 | 3 |
| 9 | SeV_Cebpa_mutant | True | 389 | 8841 | 28 |

That first run removed 49 per cent of one infected well as airway, because
an absolute cut frozen from an uninjured well over-removes where Sendai
virus has raised airway transcripts globally.

## What was not computed

No peak matrix was opened for a statistic. Every well refused at the
cell-count gate before the peak pass ran.

|  | well | n_transitional | reading | why |
|---|---|---|---|---|
| 0 | wildtype_PBS | 3 | not computable | fewer than 100 transitional cells |
| 1 | wildtype_SeV | 18 | not computable | fewer than 100 transitional cells |
| 2 | AP1mut_PBS | 4 | not computable | fewer than 100 transitional cells |
| 3 | AP1mut_SeV | 11 | not computable | fewer than 100 transitional cells |
| 4 | P9_control | 35 | not computable | fewer than 100 transitional cells |
| 5 | P9_Cebpa_mutant | 68 | not computable | fewer than 100 transitional cells |
| 6 | 7wk_control | 8 | not computable | fewer than 100 transitional cells |
| 7 | 7wk_Cebpa_mutant | 15 | not computable | fewer than 100 transitional cells |
| 8 | SeV_control | 2 | not computable | fewer than 100 transitional cells |
| 9 | SeV_Cebpa_mutant | 34 | not computable | fewer than 100 transitional cells |
