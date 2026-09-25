# Trial A2: the Il1r1 deletion cannot be verified from its own data

**Reading: NOT COMPUTABLE at both timepoints.** Rule R2 required the genotype
labels to be corroborated on Il1r1 itself before any Wnt quantity was computed.
They cannot be.

## Route A as specified was impossible, and that is the first result

The assessment proposed re-asking England et al. 2025 Figure S2L to S2N, the
only direct measurement of Axin2 and Il1r1 overlap in the primary literature,
on GSE247505 which it said was already local. **That data was never deposited.**
The superseries and both subseries contain zero mentions of ZsGreen, and the
NCBI link service returns exactly three GEO series for that paper, none holding
the sorted libraries. The only direct evidence for the overlap in the field
cannot be re-examined by anyone.

## What was asked instead, and why it was a better question

GSE247504 carries Kras-mutant AT2-derived clones on an Il1r1 heterozygous
against homozygous-deleted background, two libraries per group at two and
twelve weeks. That is a genetic perturbation of IL-1 receptor dosage with a Wnt
readout in the same cells, and it tests the direction Choi 2020 hypothesised and
nobody has tested. The reverse direction, Wnt inducing IL-1beta, is established
by Aumiller 2013 and almost never cited.

## Why it refused

|  | library | timepoint | genotype | cells | Il1r1 | wnt_module |
|---|---|---|---|---|---|---|
| 0 | GSM7890839_Expt2_ConfettiRFPr1 | Confetti | control | 1987 | 0.1308 | 1.7425 |
| 1 | GSM7890840_Expt2_ConfettiRFPr2 | Confetti | control | 2780 | 0.1313 | 1.5932 |
| 2 | GSM7890841_Expt2_2wRFP_hetr1 | 2w | het | 983 | 0.2017 | 1.3143 |
| 3 | GSM7890842_Expt2_2wRFP_hetr2 | 2w | het | 973 | 0.1619 | 0.8846 |
| 4 | GSM7890843_Expt2_2wRFP_homr1 | 2w | hom | 1015 | 0.1607 | 1.5417 |
| 5 | GSM7890844_Expt2_2wRFP_homr2 | 2w | hom | 784 | 0.1866 | 1.3463 |
| 6 | GSM7890845_Expt2_12wRFP_hetr1 | 12w | het | 2648 | 0.199 | 1.5578 |
| 7 | GSM7890846_Expt2_12wRFP_hetr2 | 12w | het | 2146 | 0.0947 | 1.4362 |
| 8 | GSM7890847_Expt2_12wRFP_homr1 | 12w | hom | 2255 | 0.2 | 1.7522 |
| 9 | GSM7890848_Expt2_12wRFP_homr2 | 12w | hom | 1998 | 0.1222 | 1.4838 |

Il1r1 per ten thousand counts is 0.202 and 0.162 in the two-week heterozygous
libraries against 0.161 and 0.187 in the homozygous ones, and 0.199 and 0.095
against 0.200 and 0.122 at twelve weeks. The groups interleave at both
timepoints. **The homozygous deletion is not visible in the deposited count**
**matrices.** Why is Not established: this page first offered a three-prime
explanation the trial had not tested, and row C149 withdrew it, because the
allele's floxed exons are 3 and 4 and its authors report that excising them
eliminates downstream expression.

A genotype contrast whose genotype cannot be verified from the data is not a
contrast, so the trial refuses. The consequence for anyone else using this
deposit is direct: **the Il1r1 genotype labels rest on the deposited metadata**
**alone and cannot be checked against the matrices.**

## What was not read, recorded so a later pass cannot present it as new

The Wnt module was computed before R2 was evaluated and is reported here unread.
At two weeks the heterozygous libraries give 1.314 and 0.885 and the homozygous
1.542 and 1.346, which is complete separation in the direction of more Wnt
without Il1r1. At twelve weeks it is 1.558 and 1.436 against 1.752 and 1.484,
which is not separation. None of this is read, because the labels generating it
are unverifiable, and because with two libraries per group the smallest
attainable permutation p is one in six and nothing here could be significant
even if the labels were sound.

## The disclosed implementation bug

A first run looked for a column named `symbol` in the features table, fell
through to the Ensembl identifier, and returned zero for every gene in every
library including the controls. That is an impossible table rather than an
outcome, so it was repaired rather than recorded as a refusal, and an assertion
now fails the trial if Sftpc comes back at zero in AT2-derived cells.
