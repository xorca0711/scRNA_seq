# Trial M1: closed or merely silenced, as frozen

**Reading: NOT COMPUTABLE in every well.** Two of the frozen rules were
wrong in ways that only reading the data could show. Under this
repository's convention the thresholds were not moved, the first outcome
stays here, and the corrected pass sits beside it as trial M1b. What M1
did establish, and what stands, is its rule R2.

## R2 corroborated one deposit and REFUTED the other

The barcode-suffix map is the GEO sample order and is not deposited, so
R2 required it to be checked on an independent axis before any name was
used. The check is the knockout carrying less of its own target gene.

|  | deposit | file | suffix | assumed_well | cpm10k_Cldn4 | cpm10k_Fos | cpm10k_Cebpa |
|---|---|---|---|---|---|---|---|
| 0 | GSE310539 | totalaggr | 1 | wildtype_PBS | 0.022 | 5.12 | 0.881 |
| 1 | GSE310539 | totalaggr | 2 | wildtype_SeV | 0.236 | 5.511 | 0.737 |
| 2 | GSE310539 | totalaggr | 3 | AP1mut_PBS | 0.023 | 1.395 | 1.355 |
| 3 | GSE310539 | totalaggr | 4 | AP1mut_SeV | 0.12 | 1.371 | 1.176 |
| 4 | GSE247130 | 7wk | 2 | 7wk_control | 0.086 | 6.549 | 0.059 |
| 5 | GSE247130 | 7wk | 1 | 7wk_mutant | 0.055 | 7.252 | 0.646 |
| 6 | GSE247130 | P9 | 1 | P9_mutant | 0.242 | 11.088 | 0.964 |
| 7 | GSE247130 | P9 | 2 | P9_control | 0.404 | 9.708 | 0.066 |
| 8 | GSE247130 | SeV | 1 | SeV_mutant | 0.065 | 2.7 | 0.204 |
| 9 | GSE247130 | SeV | 2 | SeV_control | 0.273 | 2.105 | 0.046 |

**GSE310539 is corroborated.** Fos falls from 5.12 and 5.51 counts per
10,000 in suffixes 1 and 2 to 1.40 and 1.37 in suffixes 3 and 4, so
suffixes 3 and 4 are the Fos/Fosb/Junb mutant as the GEO order says.
Cldn4 rises 12-fold with infection in the wildtype pair and 6-fold in the
mutant pair, in the direction and the ratio the source paper reports.
Note that only Fos discriminates: Fosb and Junb are not lower in the
mutant wells, so a genotype check built on those two would have passed
the wrong answer.

**GSE247130 is REFUTED, and the deposited order is inverted.** Cebpa is
11 to 14 times HIGHER in the suffix the GEO order calls the Cebpa
knockout, in all three files. A conditional knockout cannot carry more of
its own target than its control, so suffix 1 is the control and suffix 2
is the mutant, the reverse of the sample order. Two further axes agree
and neither was used to reach that conclusion: Cldn4 is 4.5-fold higher
in suffix 2 of the infected file, which is the expansion of transitional
cells the source paper reports for the mutant, and Sox9 is 7-fold higher
in suffix 2 of the neonatal file, which is the SOX9 reactivation the
source paper reports for the neonatal mutant. Under R2 every trial in
this folder now uses the inverted map for GSE247130.

## The two rule defects, disclosed and not repaired here

**R3 is a depth filter wearing the name of a compartment filter.** It
dropped any cell with a non-zero count of Scgb1a1, Scgb3a2, Foxj1 or
Krt5. Scgb1a1 is detected in 73 to 100 per cent of cells in every well at
a median of 1 to 20 counts, and the rule removed a median of 99.9 per
cent of cells. Worse, what it removed tracks depth: Scgb3a2 detection is
0.950 in the deepest well and 0.028 in the shallowest, a 34-fold spread
on a 2.1-fold depth difference. Presence of an abundant secreted
transcript in single-nucleus data is ambient, not identity.

**R4 froze its threshold on raw counts.** Median RNA UMI per cell varies
2.6-fold across the ten wells, so a raw Krt8 count frozen from one well
does not mean the same thing in another. The threshold should have been
frozen on a depth-normalised count, which is what M1b does.

## U7 passed, and it passed vacuously

The largest transitional group in an uninjured well was 3 cells, under the
floor of 100, so U7 reads as satisfied. That reading is worthless here,
because R3 had already removed almost every cell in those wells. A
negative control that passes because the data are gone is not a negative
control, and this is recorded so that M1b's U7 is not mistaken for a
second confirmation.

## What was NOT computed

No peak matrix was opened for a state-level quantity in this trial. The
chromatin question is untouched and moves to M1b.
