# Is chromatin less sparse than the transcript, for these two genes?

A measurement, not a trial. No state is compared and no biology is reported.

|  | deposit | role | gene | cells | rna_detection | distal_peaks | any_distal_peak | promoter_peaks | any_promoter_peak | all_peaks | any_peak | atac_over_rna |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | GSE310539 | marker | Axin2 | 39849 | 0.0534 | 15 | 0.1161 | 1 | 0.1356 | 16 | 0.2205 | 4.13 |
| 1 | GSE310539 | marker | Il1r1 | 39849 | 0.249 | 37 | 0.2486 | 2 | 0.0943 | 39 | 0.2988 | 1.2 |
| 2 | GSE310539 | Wnt target | Nkd1 | 39849 | 0.3554 | 34 | 0.3481 | 1 | 0.1074 | 35 | 0.3903 | 1.1 |
| 3 | GSE310539 | Wnt target | Notum | 39849 | 0.0367 | 4 | 0.0431 | 1 | 0.0496 | 5 | 0.0871 | 2.37 |
| 4 | GSE310539 | Wnt target | Lef1 | 39849 | 0.0091 | 11 | 0.0268 | 2 | 0.0231 | 13 | 0.0458 | 5.04 |
| 5 | GSE310539 | Wnt target | Tcf7l2 | 39849 | 0.3042 | 65 | 0.2775 | 2 | 0.1204 | 67 | 0.3345 | 1.1 |
| 6 | GSE310539 | Wnt target | Wif1 | 39849 | 0.0263 | 12 | 0.0798 | 1 | 0.0884 | 13 | 0.1527 | 5.81 |
| 7 | GSE310539 | Wnt target | Lgr5 | 39849 | 0.0037 | 10 | 0.0375 | 1 | 0.0148 | 11 | 0.0504 | 13.48 |
| 8 | GSE310539 | control | Sftpc | 39849 | 0.9999 | 1 | 0.0867 | 1 | 0.0757 | 2 | 0.1453 | 0.15 |
| 9 | GSE310539 | control | Etv5 | 39849 | 0.5703 | 31 | 0.3164 | 1 | 0.2255 | 32 | 0.4363 | 0.77 |
| 10 | GSE310539 | control | Krt8 | 39849 | 0.4131 | 3 | 0.1001 | 1 | 0.0963 | 4 | 0.1762 | 0.43 |
| 11 | GSE310539 | control | Cldn4 | 39849 | 0.0297 | 8 | 0.1697 | 1 | 0.032 | 9 | 0.1883 | 6.34 |
| 12 | GSE247130_SeV | marker | Axin2 | 22409 | 0.038 | 14 | 0.2399 | 1 | 0.2742 | 15 | 0.3937 | 10.37 |
| 13 | GSE247130_SeV | marker | Il1r1 | 22409 | 0.1869 | 38 | 0.4611 | 2 | 0.2237 | 40 | 0.516 | 2.76 |
| 14 | GSE247130_SeV | Wnt target | Nkd1 | 22409 | 0.3516 | 38 | 0.5168 | 2 | 0.2631 | 40 | 0.5636 | 1.6 |
| 15 | GSE247130_SeV | Wnt target | Notum | 22409 | 0.0086 | 4 | 0.0861 | 1 | 0.0957 | 5 | 0.1615 | 18.75 |
| 16 | GSE247130_SeV | Wnt target | Lef1 | 22409 | 0.0093 | 14 | 0.1034 | 1 | 0.0266 | 15 | 0.1214 | 13.01 |
| 17 | GSE247130_SeV | Wnt target | Tcf7l2 | 22409 | 0.3053 | 74 | 0.533 | 2 | 0.268 | 76 | 0.5852 | 1.92 |
| 18 | GSE247130_SeV | Wnt target | Wif1 | 22409 | 0.007 | 12 | 0.2109 | 1 | 0.2228 | 13 | 0.3444 | 49.47 |
| 19 | GSE247130_SeV | Wnt target | Lgr5 | 22409 | 0.003 | 10 | 0.1295 | 1 | 0.044 | 11 | 0.1585 | 53.01 |
| 20 | GSE247130_SeV | control | Sftpc | 22409 | 0.9995 | 1 | 0.1931 | 1 | 0.1489 | 2 | 0.2735 | 0.27 |
| 21 | GSE247130_SeV | control | Etv5 | 22409 | 0.6156 | 43 | 0.5497 | 1 | 0.4227 | 44 | 0.6521 | 1.06 |
| 22 | GSE247130_SeV | control | Krt8 | 22409 | 0.4466 | 4 | 0.3126 | 1 | 0.2305 | 5 | 0.4159 | 0.93 |
| 23 | GSE247130_SeV | control | Cldn4 | 22409 | 0.0668 | 8 | 0.4325 | 1 | 0.1519 | 9 | 0.4678 | 7.0 |

## The reading

**Yes, for Axin2, and by a large margin.** Its transcript is detected in
3.8% to 5.3%
of cells, at roughly one molecule per positive cell, which is a Poisson coin
flip rather than a phenotype. At least one of its fifteen or sixteen linked
peaks is detected in 22.1% to
39.4% of the same cells, a gain of
4.1 and 10.4 fold.
Il1r1 gains less because its transcript is less sparse to begin with.

**The controls order correctly, so the measurement is doing what it claims.**
Sftpc has a near-universal transcript and two peaks, and the chromatin route
is worse for it. Etv5 is moderately expressed with many peaks and comes out
roughly even. Lgr5 is the sparsest transcript here and shows the largest gain.
The rule is that chromatin helps exactly where the transcript is sparse and
the gene carries many linked peaks, which is the situation Axin2 is in and
Sftpc is not.

## What this does not say

An accessible Axin2 locus is not a Wnt-responsive cell. Accessibility reports
that a locus is in a configuration permitting expression, not that the cell is
signalling now, and it is slower and more permissive than transcription. This
is a weaker proxy for the lineage reporter than the transcript would be if the
transcript worked. It is on the table because the transcript does not work.
