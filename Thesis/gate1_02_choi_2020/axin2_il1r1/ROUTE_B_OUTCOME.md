# Route B: locus co-accessibility. Closed on these deposits, under a rule that was frozen first

**Outcome: NOT COMPUTABLE in every well of every pass. Route B is closed on
GSE310539 and GSE247130 under the decision rule trial A1c froze before it ran,
and there is no fourth pass.**

Route B was the one route this branch opened rather than closed. The assessment
measured that at least one of Axin2's linked peaks is detected in 22.1 and 39.4
per cent of cells against a transcript detected in 5.3 and 3.8 per cent, a gain
of 4.1 and 10.4 fold, so the chromatin readout escapes the floor that kills the
transcript readout. That was true, and it was not enough.

## Three passes, three instruments, one gate

| Trial | Instrument | Its positive control | Outcome |
|---|---|---|---|
| A1 | binary: is at least one of the locus's peaks detected | Krt8 against Krt18, 24 kb apart | **Unreadable.** The control did not clear. It was the wrong control: the vendor annotation assigns each peak to its NEAREST gene, so the two loci largely do not share peaks and the control was a biological question rather than an arithmetic identity |
| A1b | the same, with a better control | a split-locus control: two random halves of ONE gene's peak set, 20 genes, median | **Not computable.** The split-locus control did not clear either, at z +0.51 to +2.05 |
| A1c | graded: the FRACTION of the locus's peaks detected, correlated across cells | the same split-locus control | **Not computable**, at z +0.24 to +1.62, and Route B closed under A1c's own pre-registered rule |

## What the graded statistic did show, and it is the useful part

The binary form could not detect co-accessibility between anything. The graded
form can:

|  | well | pair | correlation | null_mean | null_sd | z | clears |
|---|---|---|---|---|---|---|---|
| 2 | wildtype_PBS | Krt8/Krt18 | 0.02365 | 0.00288 | 0.01619 | 1.283 | False |
| 7 | wildtype_SeV | Krt8/Krt18 | -0.00383 | 0.00121 | 0.01585 | -0.318 | False |
| 12 | AP1mut_PBS | Krt8/Krt18 | 0.00501 | 0.00226 | 0.00998 | 0.276 | False |
| 17 | AP1mut_SeV | Krt8/Krt18 | 0.00054 | 0.00385 | 0.01285 | -0.258 | False |
| 22 | 7wk_control | Krt8/Krt18 | 0.02447 | 0.00507 | 0.01378 | 1.407 | False |
| 27 | 7wk_Cebpa_mutant | Krt8/Krt18 | 0.0393 | 0.00534 | 0.01163 | 2.919 | False |
| 32 | SeV_control | Krt8/Krt18 | 0.02501 | 0.00503 | 0.01399 | 1.428 | False |
| 37 | SeV_Cebpa_mutant | Krt8/Krt18 | 0.06095 | 0.00514 | 0.01555 | 3.589 | True |
| 3 | wildtype_PBS | Etv5/Abca3 | -0.00131 | 0.00288 | 0.01619 | -0.259 | False |
| 8 | wildtype_SeV | Etv5/Abca3 | 0.03746 | 0.00121 | 0.01585 | 2.287 | False |
| 13 | AP1mut_PBS | Etv5/Abca3 | 0.02345 | 0.00226 | 0.00998 | 2.124 | False |
| 18 | AP1mut_SeV | Etv5/Abca3 | 0.00669 | 0.00385 | 0.01285 | 0.221 | False |
| 23 | 7wk_control | Etv5/Abca3 | 0.02891 | 0.00507 | 0.01378 | 1.729 | False |
| 28 | 7wk_Cebpa_mutant | Etv5/Abca3 | 0.03625 | 0.00534 | 0.01163 | 2.657 | False |
| 33 | SeV_control | Etv5/Abca3 | 0.05465 | 0.00503 | 0.01399 | 3.547 | True |
| 38 | SeV_Cebpa_mutant | Etv5/Abca3 | 0.04139 | 0.00514 | 0.01555 | 2.331 | False |

**Etv5 against Abca3 clears at z = +3.55 in one well and Krt8 against Krt18 at
z = +3.59 in another.** Those are two AT2 identity genes on different
chromosomes, and two keratins 24 kilobases apart. So the graded instrument has
real sensitivity to co-accessibility, both the kind that comes from shared
regulation and the kind that comes from physical proximity. It is not a dead
instrument; it is an instrument whose gate was set by a control that could not
pass it.

## The defect in the gate, disclosed rather than repaired

The split-locus control was introduced in A1b as one that could not fail for
biological reasons. That claim was too strong, in two ways this folder did not
see until A1c had run:

1. **It is judged against the wrong null.** The matched-pair null is built for
   a pair with Axin2's peak count on one side and Il1r1's on the other, roughly
   15 and 38. A split-locus control on a gene with 30 peaks is 15 against 15.
   Its two scores are each estimated from half as many peaks as the null's
   larger side, so its correlation is attenuated relative to the distribution
   it is being compared with.
2. **Peaks within one gene's annotation are not necessarily one regulatory
   domain.** A promoter peak and a distal peak 100 kilobases away are both
   linked to the same gene by the vendor annotation and need not be
   co-accessible at all.

So the closure of Route B is **a decision taken under a frozen rule, not a
demonstration that the approach cannot work.** The rule was frozen before the
trial ran and it is honoured here rather than revised, because revising a gate
after it refuses is how a result gets manufactured. But the register status is
Not established, never Refuted, and the next section says what a correct gate
would be.

## What the numbers were, with the caveat that none of them is readable

Recorded because a reader is entitled to see what the trial saw, and because a
later pass with a correct gate must not be able to present these as new.

|  | well | pair | correlation | null_mean | null_sd | z | percentile |
|---|---|---|---|---|---|---|---|
| 0 | wildtype_PBS | Axin2/Il1r1 | 0.01527 | 0.00288 | 0.01619 | 0.765 | 0.7833 |
| 5 | wildtype_SeV | Axin2/Il1r1 | 0.00093 | 0.00121 | 0.01585 | -0.018 | 0.4833 |
| 10 | AP1mut_PBS | Axin2/Il1r1 | -0.00644 | 0.00226 | 0.00998 | -0.872 | 0.1967 |
| 15 | AP1mut_SeV | Axin2/Il1r1 | 0.01129 | 0.00385 | 0.01285 | 0.579 | 0.72 |
| 20 | 7wk_control | Axin2/Il1r1 | -0.01287 | 0.00507 | 0.01378 | -1.302 | 0.09 |
| 25 | 7wk_Cebpa_mutant | Axin2/Il1r1 | -0.01289 | 0.00534 | 0.01163 | -1.567 | 0.05 |
| 30 | SeV_control | Axin2/Il1r1 | -0.01879 | 0.00503 | 0.01399 | -1.703 | 0.04 |
| 35 | SeV_Cebpa_mutant | Axin2/Il1r1 | -0.01284 | 0.00514 | 0.01555 | -1.156 | 0.1133 |
| 1 | wildtype_PBS | WNT_MODULE/Il1r1 | 0.00942 | 0.00288 | 0.01619 | 0.404 | 0.6633 |
| 6 | wildtype_SeV | WNT_MODULE/Il1r1 | 0.02096 | 0.00121 | 0.01585 | 1.246 | 0.89 |
| 11 | AP1mut_PBS | WNT_MODULE/Il1r1 | -0.01898 | 0.00226 | 0.00998 | -2.128 | 0.0267 |
| 16 | AP1mut_SeV | WNT_MODULE/Il1r1 | 0.01286 | 0.00385 | 0.01285 | 0.701 | 0.7633 |
| 21 | 7wk_control | WNT_MODULE/Il1r1 | 0.00195 | 0.00507 | 0.01378 | -0.226 | 0.44 |
| 26 | 7wk_Cebpa_mutant | WNT_MODULE/Il1r1 | 0.02129 | 0.00534 | 0.01163 | 1.371 | 0.9167 |
| 31 | SeV_control | WNT_MODULE/Il1r1 | -0.0022 | 0.00503 | 0.01399 | -0.517 | 0.31 |
| 36 | SeV_Cebpa_mutant | WNT_MODULE/Il1r1 | 0.00147 | 0.00514 | 0.01555 | -0.236 | 0.3767 |

In the AT2 sets the Axin2 against Il1r1 correlation is **negative in every well**
and the Wnt module against Il1r1 is near zero. In the one transitional set large
enough to compute, both flip positive: Axin2 against Il1r1 at r = +0.0697
(z = +1.61) and the Wnt module at r = +0.0917 (z = +2.18), on 773 cells.
Neither clears, both are unreadable, and the direction of a difference between
cell states is exactly the kind of thing an underpowered design invents.

## What would make Route B computable, for whoever picks it up

1. **A gate matched to the null it is judged against.** Build the split-locus
   control from genes with about 53 peaks, split 15 against 38 rather than in
   half, so the control has the same peak-count profile as the pair under test.
   Restrict the split to peaks within 50 kilobases of each other so the halves
   really are one regulatory domain.
2. **Use the fragments files.** Both deposits carry `atac_fragments.tsv.gz`,
   2.6 GB for GSE310539, which this branch has not downloaded. They allow peaks
   to be re-called on the analysis cells, per-cell quality control and TSS
   enrichment, none of which the filtered matrix supports.
3. **Do not raise the depth budget to buy significance.** The budget was held at
   the twentieth percentile through all three passes precisely so that it could
   not be tuned, and it should stay frozen for a fourth.
