# Trial M1c: the label at the depth available, and the negative control that fired

**Reading: TRIAL UNREADABLE under rule U7.** And the reason U7 fired is
worth more than the trial was.

## What U7 required and what happened

U7 required the uninjured wells to hold less of the labelled group than
the injured ones, because both source papers report the CLDN4-positive
state to be absent or Cldn4-negative without injury. They did not:

|  | well | injured | after_filter | transitional | transitional_pct |
|---|---|---|---|---|---|
| 0 | wildtype_PBS | False | 5492 | 18 | 0.33 |
| 1 | wildtype_SeV | True | 4578 | 80 | 1.75 |
| 2 | AP1mut_PBS | False | 8134 | 18 | 0.22 |
| 3 | AP1mut_SeV | True | 7444 | 55 | 0.74 |
| 4 | P9_control | False | 10475 | 387 | 3.69 |
| 5 | P9_Cebpa_mutant | False | 10252 | 827 | 8.07 |
| 6 | 7wk_control | False | 6483 | 55 | 0.85 |
| 7 | 7wk_Cebpa_mutant | False | 9016 | 140 | 1.55 |
| 8 | SeV_control | True | 6344 | 44 | 0.69 |
| 9 | SeV_Cebpa_mutant | True | 7851 | 400 | 5.09 |

The two neonatal P9 wells labelled 3.69 and
8.07 per cent, above every injured well, the
largest of which reached 5.09 per cent. Enrichment by
injury came out at 0.6-fold, below the threefold floor,
so the whole-trial negative control fired and the trial refused.

## The refusal is one of the two answers this branch was built to find

The owner's first question was whether the transition state carries a
character specific to regeneration or disease. At the transcript level,
for the markers that define it, the answer here is no. Krt8 and Cldn4 are
expressed across immature postnatal alveolar epithelium, so a
CLDN4-positive KRT8-positive call cannot separate neonatal developmental
immaturity from injury-induced transition.

Hassan and Chen's own argument predicts this and never states it as a
limitation of the marker: their whole case is that neonatal AT2 cells are
plastic and that Cebpa deletion returns mature cells toward that neonatal
state. A marker set shared with normal development is not, by itself, a
damage-associated marker set.

## What the label does do, where the stage is held fixed

Within the mature wells the label behaves as both papers describe. In
GSE310539 it rises from 0.33 and 0.22 per cent under PBS to 1.75 and 0.74
per cent after Sendai virus, and the AP-1 mutant reaches 0.42 of the
wildtype value against the 0.5 the source paper reports. In the mature
GSE247130 wells it rises from 0.69 per cent in the infected control to
5.09 per cent in the infected Cebpa mutant, which is the expansion that
paper reports. None of those comparisons is testable: each condition is
one library.
