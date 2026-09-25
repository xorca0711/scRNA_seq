# Trial M4: the corroboration the register has been quoting

Every number on this page is formatted from the table below. None is typed.

## Why this trial exists

An eight-adversary audit of rows C116 to C150, each finding sent to an
independent verifier, returned 45 confirmed defects. The worst sits in this
branch's self-described one validated result: **the Sox9 axis the register
presents as independent corroboration was never computed.** Trial M1's
frozen panel was Cldn4, Fos and Cebpa; its write loop iterated over four
genes behind a guard that silently dropped the fourth. The figure existed
only as a hard-coded string in M1's narrative paragraph, and propagated from
there into five other files.

The same paragraph was the source of two further wrong numbers, for the same
reason: it was typed as string literals while the dataframe holding the real
values sat in scope two lines above.

## The corroboration, computed

|  | deposit | file | suffix | well | cells | total_rna_counts | cpm10k_Cldn4 | cpm10k_Fos | cpm10k_Cebpa | cpm10k_Sox9 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | GSE310539 | totalaggr | 1 | wildtype_PBS | 7340 | 38879858 | 0.0224 | 5.1199 | 0.8807 | 0.0139 |
| 1 | GSE310539 | totalaggr | 2 | wildtype_SeV | 8093 | 57932184 | 0.2363 | 5.5108 | 0.7366 | 0.0649 |
| 2 | GSE310539 | totalaggr | 3 | AP1mut_PBS | 13622 | 43681176 | 0.0227 | 1.3951 | 1.3553 | 0.0133 |
| 3 | GSE310539 | totalaggr | 4 | AP1mut_SeV | 10794 | 52034333 | 0.1199 | 1.3708 | 1.1763 | 0.0234 |
| 4 | GSE247130 | P9 | 1 | P9_mutant | 12186 | 121903574 | 0.2422 | 11.0879 | 0.9642 | 0.318 |
| 5 | GSE247130 | P9 | 2 | P9_control | 11272 | 137528687 | 0.4044 | 9.7084 | 0.0664 | 2.2429 |
| 6 | GSE247130 | 7wk | 1 | 7wk_mutant | 7589 | 56653751 | 0.0552 | 7.2521 | 0.646 | 0.0115 |
| 7 | GSE247130 | 7wk | 2 | 7wk_control | 10838 | 78181340 | 0.086 | 6.5494 | 0.0592 | 0.0325 |
| 8 | GSE247130 | SeV | 1 | SeV_mutant | 11773 | 55219121 | 0.0647 | 2.7002 | 0.2043 | 0.0246 |
| 9 | GSE247130 | SeV | 2 | SeV_control | 10636 | 74679531 | 0.2728 | 2.1053 | 0.0455 | 0.0611 |

## The ratios, computed

|  | deposit | file | gene | axis | numerator | denominator | ratio |
|---|---|---|---|---|---|---|---|
| 0 | GSE247130 | P9 | Cebpa | control over Cebpa mutant | 0.9642 | 0.0664 | 14.521 |
| 1 | GSE247130 | P9 | Sox9 | Cebpa mutant over control | 2.2429 | 0.318 | 7.053 |
| 2 | GSE247130 | P9 | Cldn4 | Cebpa mutant over control | 0.4044 | 0.2422 | 1.67 |
| 3 | GSE247130 | 7wk | Cebpa | control over Cebpa mutant | 0.646 | 0.0592 | 10.912 |
| 4 | GSE247130 | 7wk | Sox9 | Cebpa mutant over control | 0.0325 | 0.0115 | 2.826 |
| 5 | GSE247130 | 7wk | Cldn4 | Cebpa mutant over control | 0.086 | 0.0552 | 1.558 |
| 6 | GSE247130 | SeV | Cebpa | control over Cebpa mutant | 0.2043 | 0.0455 | 4.49 |
| 7 | GSE247130 | SeV | Sox9 | Cebpa mutant over control | 0.0611 | 0.0246 | 2.484 |
| 8 | GSE247130 | SeV | Cldn4 | Cebpa mutant over control | 0.2728 | 0.0647 | 4.216 |
| 9 | GSE310539 | totalaggr | Fos | wildtype over AP-1 mutant | 5.1199 | 1.3951 | 3.67 |
| 10 | GSE310539 | totalaggr | Cldn4 | Sendai over PBS, wildtype | 0.2363 | 0.0224 | 10.549 |
| 11 | GSE310539 | totalaggr | Cldn4 | Sendai over PBS, AP-1 mutant | 0.1199 | 0.0227 | 5.282 |

## What the register said and what is true

| The register said | Computed here |
|---|---|
| Cebpa 11 to 14 times higher in all three files | P9 14.5, 7wk 10.9, SeV 4.5. The direction holds in all three; the range does not, and the infected file, which every downstream trial uses, is the weakest |
| Sox9 7-fold higher in the neonatal file, confirmed independently | P9 7.1, 7wk 2.8, SeV 2.5. The neonatal figure was right, and it was never computed by the logged run that claimed it |
| Cldn4 4.5-fold higher in suffix 2 of the infected file | 4.2 |
| Cldn4 rises 12-fold with infection in the wildtype pair and 6-fold in the mutant pair | 10.5 and 5.3 |
| Fos falls from 5.12 and 5.51 to 1.40 and 1.37 | holds; the smallest wildtype over the largest mutant is 3.7 |

## What stands and what does not

**The suffix map of GSE247130 is still inverted, and that is unaffected.**
Cebpa is higher in the suffix the GEO sample order calls the knockout in
all three files (P9 14.5, 7wk 10.9, SeV 4.5), and a conditional knockout cannot carry more
of its own target than its control. Sox9 agrees in all three files
(P9 7.1, 7wk 2.8, SeV 2.5), with the neonatal file much the strongest, which is the
stage-specificity Hassan and Chen report. Cldn4 agrees in the infected file.

**What does not stand is the register's account of how well it was shown.**
One of the two axes it called independent corroboration was not measured,
and two of the magnitudes it quoted were wrong. The map was right for
reasons that were partly unlogged, which by this repository's own rule means
they were Not established at the moment they were quoted.
