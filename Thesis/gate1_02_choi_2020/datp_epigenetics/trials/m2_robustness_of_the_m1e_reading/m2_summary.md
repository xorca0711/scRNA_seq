# Trial M2: M1e's chromatin reading does not survive its own arithmetic

**M2 changes no rule and asks no new question. It attacks M1e's reading on
four fronts, and the chromatin half does not survive.** The RNA half does,
and comfortably.

## The stability table, four downsampling seeds

|  | well | arm | atac_clears | seeds | atac_z_min | atac_z_max | rna_clears | rna_z_min | rna_z_max | detectable_floor |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 7wk_Cebpa_mutant | AT1 | 0 | 4 | 1.4703 | 2.5661 | 0 | 0.7077 | 2.5064 | 0.0091 |
| 1 | 7wk_Cebpa_mutant | AT2_identity | 0 | 4 | -1.4987 | -0.4305 | 3 | -4.8737 | -2.9617 | 0.0043 |
| 2 | 7wk_Cebpa_mutant | transitional | 0 | 4 | -1.0703 | 1.0368 | 4 | 5.1039 | 8.0394 | 0.0076 |
| 3 | SeV_Cebpa_mutant | AT1 | 0 | 4 | -0.962 | 1.5009 | 0 | 0.5636 | 0.822 | 0.0091 |
| 4 | SeV_Cebpa_mutant | AT2_identity | 0 | 4 | -1.8786 | -1.2585 | 4 | -11.4224 | -8.7604 | 0.0046 |
| 5 | SeV_Cebpa_mutant | transitional | 2 | 4 | 1.909 | 4.2245 | 4 | 13.4021 | 18.7408 | 0.0048 |
| 6 | wildtype_SeV | AT1 | 0 | 4 | -0.1638 | 0.8501 | 0 | 0.4406 | 0.9194 | 0.012 |
| 7 | wildtype_SeV | AT2_identity | 0 | 4 | -1.1317 | 0.6808 | 4 | -7.4813 | -4.0206 | 0.0079 |
| 8 | wildtype_SeV | transitional | 3 | 4 | 2.1907 | 4.0896 | 4 | 9.335 | 13.3177 | 0.0083 |

## What this says, front by front

**1. The RNA half is solid.** The AT2 identity arm clears in every seed of
both injured wells, at z from -4.0 to -11.4, losing 12.3 to 16.7 detection
points across the seeds (12.6 and 17.1 at the single M1e budget). The transitional arm clears in every seed everywhere. The AT1 arm
clears in none, which is the negative arm behaving. Per gene the AT2 loss
is carried by at least five genes, not one: in the two injured wells Etv5
falls 0.227 and 0.239, Napsa 0.272 and 0.167, Slc34a2 0.271 and 0.189,
Abca3 0.182 and 0.212, and Lamp3 0.160 and 0.115.

**2. The chromatin positive control is not stable, which is fatal to the
null beside it.** M1e's rule R9 permits the AT2 chromatin null to be read
as retention only behind a transitional arm that fired. Across seeds that
arm clears three times in four in wildtype_SeV and twice in four in
SeV_Cebpa_mutant, five of eight in total. Under leave-one-out it fails on
two of its six genes in each injured well:

|  | well | dropped_gene | atac_z |
|---|---|---|---|
| 9 | wildtype_SeV | Ndrg1 | 2.802 |
| 12 | wildtype_SeV | Sfn | 2.716 |
| 28 | 7wk_Cebpa_mutant | Ndrg1 | 0.76 |
| 29 | 7wk_Cebpa_mutant | Sprr1a | 1.161 |
| 30 | 7wk_Cebpa_mutant | AW112010 | 1.306 |
| 31 | 7wk_Cebpa_mutant | Sfn | 0.59 |
| 32 | 7wk_Cebpa_mutant | Krt19 | 0.645 |
| 33 | 7wk_Cebpa_mutant | Lgals3 | 1.105 |
| 47 | SeV_Cebpa_mutant | Ndrg1 | 2.933 |
| 52 | SeV_Cebpa_mutant | Lgals3 | 2.912 |

Per gene, the arm is carried by Sfn and Ndrg1. In wildtype_SeV, Sfn alone
returns 0.0357 over seven peaks against an arm mean of 0.0097.

**3. The offset correction was not inert, it was ABSENT, and replacing it
changes a verdict.** M1e subtracted the genome-wide MEDIAN distal
difference, which came out exactly 0.0 in every well because more than
half of all distal peaks are detected in neither group at this budget.
That is not the same as there being nothing to remove. The mean global
distal difference is +0.00080, +0.00034 and +0.00201 in the three wells,
POSITIVE IN ALL THREE, across two deposits, two peak atlases and two
fragment budgets. A consistently signed nuisance added to two different
true values is precisely what produces arm values that disagree in sign,
so the sign disagreement between the two wells damns the analysis rather
than rescuing it. Recomputed against the mean:

|  | well | arm | offset_kind | offset | observed | z | clears |
|---|---|---|---|---|---|---|---|
| 0 | wildtype_SeV | AT2_identity | median | 0.0 | -0.00221 | -1.13173 | False |
| 1 | wildtype_SeV | AT2_identity | mean | 0.00079 | -0.00301 | -1.4936 | False |
| 2 | wildtype_SeV | AT2_identity | trimmed | 0.0001 | -0.00232 | -1.17868 | False |
| 3 | wildtype_SeV | transitional | median | 0.0 | 0.00971 | 3.54119 | True |
| 4 | wildtype_SeV | transitional | mean | 0.00079 | 0.00892 | 3.25996 | True |
| 5 | wildtype_SeV | transitional | trimmed | 0.0001 | 0.00961 | 3.5047 | True |
| 36 | 7wk_Cebpa_mutant | AT2_identity | median | 0.0 | -0.00206 | -1.49866 | False |
| 37 | 7wk_Cebpa_mutant | AT2_identity | mean | 0.00034 | -0.0024 | -1.72934 | False |
| 38 | 7wk_Cebpa_mutant | AT2_identity | trimmed | 1e-05 | -0.00208 | -1.50692 | False |
| 39 | 7wk_Cebpa_mutant | transitional | median | 0.0 | 0.00238 | 1.03684 | False |
| 40 | 7wk_Cebpa_mutant | transitional | mean | 0.00034 | 0.00205 | 0.85605 | False |
| 41 | 7wk_Cebpa_mutant | transitional | trimmed | 1e-05 | 0.00237 | 1.03037 | False |
| 72 | SeV_Cebpa_mutant | AT2_identity | median | 0.0 | -0.00247 | -1.87859 | False |
| 73 | SeV_Cebpa_mutant | AT2_identity | mean | 0.00201 | -0.00448 | -3.23099 | True |
| 74 | SeV_Cebpa_mutant | AT2_identity | trimmed | 0.00139 | -0.00386 | -2.81281 | False |
| 75 | SeV_Cebpa_mutant | transitional | median | 0.0 | 0.00515 | 4.22451 | True |
| 76 | SeV_Cebpa_mutant | transitional | mean | 0.00201 | 0.00314 | 2.63693 | False |
| 77 | SeV_Cebpa_mutant | transitional | trimmed | 0.00139 | 0.00376 | 3.12783 | True |

In SeV_Cebpa_mutant the AT2 arm moves from z = -1.88 and not clearing, to
z = -3.23 and clearing in the CLOSING direction, purely on the choice
between a median that was degenerate and a mean that was not. A reading
that turns on that choice is not a reading.

**And the direction matters, because it is not the one M1e reported.** An
independent reviewer reproduced the cell selection and pushed this further:
under a background matched on baseline accessibility the same arm reaches
z = -2.89, and under a multiplicative correction z = -3.55, which would
clear as AT2 programme CLOSED. Under that same correction the
transitional positive control falls to z = +2.53, below the frozen floor,
so rule R9 would have refused the well outright. Across the estimators,
the well reads either NOT COMPUTABLE or CLOSED. It never reads what M1e
reported. Nothing here supports the idea that AT2 chromatin stays open,
and if anything the salvageable signal points the other way.

**4. The effects and the detection floor are the same size, and the two
arms are not commensurable.** Three sham standard deviations, the smallest
value the frozen rule would have called a clearance, is 0.0043 to 0.0120
in detection fraction. The transitional arm's observed effects are 0.004
to 0.012. The design is working at its own floor.

Worse, an adversarial review found that the AT2 arm is BOUNDED and the
positive control is NOT. The reference-side baseline of the AT2 arm is
about 0.0185 per peak in GSE310539 and 0.0380 in GSE247130, so the arm
cannot fall below those values even if every distal element at all nine
AT2 loci shut completely. M1e's defence, that an effect the size of the
transitional arm's would have cleared, therefore describes a counterfactual
requiring 58 per cent of ALL distal accessibility at those loci to vanish
in fourteen days. The transitional arm has no such ceiling, because de novo
opening can drive it arbitrarily high. An unbounded statistic was being
used to certify the power of a bounded one.

The honest statement of what was measured is therefore not that nothing
moved, but that **no change was detected at a bar corresponding to a 45
per cent loss of distal accessibility at the AT2 loci in GSE310539, and a
12 per cent loss in GSE247130.** The second of those is a bound worth
having. The first excludes almost nothing.

## The reading

**The chromatin question is NOT ESTABLISHED from these matrices at this
depth.** These data cannot separate the AT2 programme's chromatin not
moving from this design being unable to see it move. M1e's reading of
"silenced but not closed" is withdrawn and stays in the record beside this
page. That phrase is the output of decision rule R10a, a label and not a
measurement, and it appears in this folder only in quotation marks and
only as the trial's own output.

**What is established** is the RNA half, and it is worth stating on its
own: in two deposits, the CLDN4-positive KRT8-positive alveolar group
loses the AT2 identity programme by 12.3 to 16.7 detection points against
a sham band, across five or more genes, with the AT1 programme flat in the
same cells. That is Descriptive only, because each well is one library
pooling two mice.

**And the two deposits are not independent, which an earlier draft of this
page got wrong.** Jichao Chen is a contributor on both GEO series and the
contact laboratory for GSE247130, so this is one laboratory with two first
authors, not two laboratories. They also share the SftpcCreER and
RosaSun1GFP lineage tools, the Sendai virus model, E-cadherin-positive
sorting, the 10x Multiome kit and cellranger-arc on mm10, and on this side
they share the code, the seeds, the vendor annotation and the budget rule.
What genuinely differs is the first author, the mouse cohort, the
institution, the year and the genotype. Agreement across them is a
consistency check, not a replication.

## What would settle the chromatin question

Not more statistics on these matrices. The binding constraints are the
size of the labelled group, 64 and 320 cells, and the fact that peaks were
called on the whole library and are therefore ascertained on the majority
population. Both deposits also carry `atac_fragments.tsv.gz` files, which
this branch did not download and which permit peaks to be re-called on the
labelled cells themselves, proper per-cell ATAC quality control, TSS
enrichment and footprinting. That is the next instrument, not another
threshold.

## Per-gene table

|  | well | arm | gene | distal_peaks | atac_delta | rna_delta |
|---|---|---|---|---|---|---|
| 0 | wildtype_SeV | AT2_identity | Etv5 | 31 | 0.0005 | -0.2272 |
| 1 | wildtype_SeV | AT2_identity | Abca3 | 11 | -0.0075 | -0.1823 |
| 2 | wildtype_SeV | AT2_identity | Cebpa | 3 | 0.0091 | -0.0837 |
| 3 | wildtype_SeV | AT2_identity | Lamp3 | 18 | -0.0004 | -0.1596 |
| 4 | wildtype_SeV | AT2_identity | Sftpb | 13 | -0.0003 | -0.1105 |
| 5 | wildtype_SeV | AT2_identity | Slc34a2 | 12 | 0.0049 | -0.2713 |
| 6 | wildtype_SeV | AT2_identity | Lyz2 | 4 | -0.0049 | -0.0692 |
| 7 | wildtype_SeV | AT2_identity | Napsa | 7 | -0.0173 | -0.2722 |
| 8 | wildtype_SeV | AT2_identity | Pon1 | 30 | -0.004 | 0.0 |
| 9 | wildtype_SeV | transitional | Ndrg1 | 10 | 0.0148 | 0.0597 |
| 10 | wildtype_SeV | transitional | Sprr1a | 1 | 0.0 | 0.0891 |
| 11 | wildtype_SeV | transitional | AW112010 | 6 | -0.0052 | 0.056 |
| 12 | wildtype_SeV | transitional | Sfn | 7 | 0.0357 | 0.2421 |
| 13 | wildtype_SeV | transitional | Krt19 | 13 | 0.0024 | 0.3421 |
| 14 | wildtype_SeV | transitional | Lgals3 | 13 | 0.0105 | 0.2318 |
| 15 | wildtype_SeV | AT1 | Ager | 1 | 0.0078 | -0.1492 |
| 16 | wildtype_SeV | AT1 | Hopx | 4 | -0.0049 | 0.1895 |
| 17 | wildtype_SeV | AT1 | Pdpn | 13 | 0.0033 | 0.0395 |
| 18 | wildtype_SeV | AT1 | Cav1 | 5 | 0.0039 | -0.0347 |
| 19 | wildtype_SeV | AT1 | Akap5 | 4 | 0.0059 | 0.0422 |
| 20 | 7wk_Cebpa_mutant | AT2_identity | Etv5 | 30 | -0.0011 | -0.1094 |
| 21 | 7wk_Cebpa_mutant | AT2_identity | Abca3 | 10 | -0.0062 | -0.1077 |
| 22 | 7wk_Cebpa_mutant | AT2_identity | Lamp3 | 18 | -0.0017 | -0.0433 |
| 23 | 7wk_Cebpa_mutant | AT2_identity | Sftpb | 10 | -0.0031 | -0.0013 |
| 24 | 7wk_Cebpa_mutant | AT2_identity | Slc34a2 | 10 | -0.0009 | -0.0077 |
| 25 | 7wk_Cebpa_mutant | AT2_identity | Lyz2 | 3 | 0.0 | -0.0531 |
| 26 | 7wk_Cebpa_mutant | AT2_identity | Napsa | 6 | -0.0052 | -0.0778 |
| 27 | 7wk_Cebpa_mutant | AT2_identity | Pon1 | 27 | 0.0018 | -0.0012 |
| 28 | 7wk_Cebpa_mutant | transitional | Ndrg1 | 7 | 0.0038 | 0.0227 |
| 29 | 7wk_Cebpa_mutant | transitional | Sprr1a | 1 | 0.0 | 0.0195 |
| 30 | 7wk_Cebpa_mutant | transitional | AW112010 | 4 | -0.0006 | 0.0203 |
| 31 | 7wk_Cebpa_mutant | transitional | Sfn | 5 | 0.0067 | 0.0576 |
| 32 | 7wk_Cebpa_mutant | transitional | Krt19 | 13 | 0.0038 | 0.0891 |
| 33 | 7wk_Cebpa_mutant | transitional | Lgals3 | 8 | 0.0006 | 0.1102 |
| 34 | 7wk_Cebpa_mutant | AT1 | Ager | 1 | 0.0223 | -0.0098 |
| 35 | 7wk_Cebpa_mutant | AT1 | Hopx | 3 | -0.0022 | 0.061 |
| 36 | 7wk_Cebpa_mutant | AT1 | Pdpn | 13 | 0.0015 | -0.0109 |
| 37 | 7wk_Cebpa_mutant | AT1 | Cav1 | 4 | -0.0017 | 0.0279 |
| 38 | 7wk_Cebpa_mutant | AT1 | Akap5 | 4 | 0.0106 | 0.0395 |
| 39 | SeV_Cebpa_mutant | AT2_identity | Etv5 | 43 | 0.0019 | -0.2391 |
| 40 | SeV_Cebpa_mutant | AT2_identity | Abca3 | 11 | -0.0065 | -0.2119 |
| 41 | SeV_Cebpa_mutant | AT2_identity | Lamp3 | 19 | -0.001 | -0.1147 |
| 42 | SeV_Cebpa_mutant | AT2_identity | Sftpb | 13 | 0.0037 | -0.0171 |
| 43 | SeV_Cebpa_mutant | AT2_identity | Slc34a2 | 18 | -0.0004 | -0.1887 |
| 44 | SeV_Cebpa_mutant | AT2_identity | Lyz2 | 5 | 0.0011 | -0.0647 |
| 45 | SeV_Cebpa_mutant | AT2_identity | Napsa | 6 | -0.0201 | -0.1667 |
| 46 | SeV_Cebpa_mutant | AT2_identity | Pon1 | 37 | 0.0015 | -0.0022 |
| 47 | SeV_Cebpa_mutant | transitional | Ndrg1 | 14 | 0.0109 | 0.0394 |
| 48 | SeV_Cebpa_mutant | transitional | Sprr1a | 2 | 0.0098 | 0.0952 |
| 49 | SeV_Cebpa_mutant | transitional | AW112010 | 6 | 0.0007 | 0.007 |
| 50 | SeV_Cebpa_mutant | transitional | Sfn | 7 | 0.0002 | 0.2271 |
| 51 | SeV_Cebpa_mutant | transitional | Krt19 | 13 | 0.0037 | 0.1654 |
| 52 | SeV_Cebpa_mutant | transitional | Lgals3 | 14 | 0.0056 | 0.1464 |
| 53 | SeV_Cebpa_mutant | AT1 | Ager | 1 | -0.0055 | -0.0175 |
| 54 | SeV_Cebpa_mutant | AT1 | Hopx | 4 | 0.0029 | 0.0607 |
| 55 | SeV_Cebpa_mutant | AT1 | Pdpn | 14 | -0.0025 | -0.0282 |
| 56 | SeV_Cebpa_mutant | AT1 | Cav1 | 10 | -0.0023 | -0.0226 |
| 57 | SeV_Cebpa_mutant | AT1 | Akap5 | 5 | -0.0025 | 0.0335 |
