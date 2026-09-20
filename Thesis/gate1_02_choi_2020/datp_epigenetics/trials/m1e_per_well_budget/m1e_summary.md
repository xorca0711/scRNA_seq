# Trial M1e: the first pass to reach a statistic, and its reading is superseded

**M1e read SILENCED BUT NOT CLOSED in two test wells. Trial M2 then showed
that the chromatin half of that reading is not stable, and it is withdrawn.**
This page records what M1e computed and what survived. The withdrawal is
in [`../m2_robustness_of_the_m1e_reading/m2_summary.md`](../m2_robustness_of_the_m1e_reading/m2_summary.md).

## What it computed

|  | well | role | n_labelled | B_atac | reading | why |
|---|---|---|---|---|---|---|
| 0 | wildtype_PBS | negative control |  |  | not computable | fewer than 50 labelled cells |
| 1 | wildtype_SeV | test | 64.0 | 2700.0 | AT2 programme SILENCED BUT NOT CLOSED | R10a: clears in RNA and not in ATAC, behind a positive control that fired |
| 2 | AP1mut_PBS | negative control |  |  | not computable | fewer than 50 labelled cells |
| 3 | AP1mut_SeV | test | 44.0 | 2368.0 | not computable | fewer than 50 labelled cells above the budget, or too small a sham pool |
| 4 | 7wk_control | negative control | 44.0 | 2768.0 | not computable | fewer than 50 labelled cells above the budget, or too small a sham pool |
| 5 | 7wk_Cebpa_mutant | negative control | 112.0 | 2177.0 | not computable | R9: the transitional-marker arm did not clear, so the ATAC side has no demonstrated sensitivity here |
| 6 | SeV_control | test |  |  | not computable | fewer than 50 labelled cells |
| 7 | SeV_Cebpa_mutant | test | 320.0 | 6985.0 | AT2 programme SILENCED BUT NOT CLOSED | R10a: clears in RNA and not in ATAC, behind a positive control that fired |

## The arms

|  | well | arm | modality | genes_used | observed | sham_sd | z | clears_band |
|---|---|---|---|---|---|---|---|---|
| 0 | wildtype_SeV | AT2_identity | atac | 9 | 0.0012958070186732 | 0.0026683690174919 | 0.6301312645003281 | False |
| 1 | wildtype_SeV | AT2_identity | rna | 9 | -0.1709822755604351 | 0.0297604486304115 | -5.92858432129269 | True |
| 2 | wildtype_SeV | transitional | atac | 6 | 0.0107130933302808 | 0.0023911189910804 | 4.402151435532966 | True |
| 3 | wildtype_SeV | transitional | rna | 6 | 0.1916806000830849 | 0.015902776159142 | 12.043414300149903 | True |
| 4 | wildtype_SeV | AT1 | atac | 5 | -0.005 | 0.0045893723884884 | -0.9838980970342636 | False |
| 5 | wildtype_SeV | AT1 | rna | 5 | 0.0096574346865932 | 0.0232630030049682 | 0.3752818964882206 | False |
| 6 | 7wk_Cebpa_mutant | AT2_identity | atac | 8 | -0.001922123015873 | 0.0020654309174224 | -0.9392466920011951 | False |
| 7 | 7wk_Cebpa_mutant | AT2_identity | rna | 8 | -0.0608136327852565 | 0.0103964807560988 | -5.713525444852438 | True |
| 8 | 7wk_Cebpa_mutant | transitional | atac | 6 | -0.005741533392203 | 0.0027917475436756 | -2.0719942175891783 | False |
| 9 | 7wk_Cebpa_mutant | transitional | rna | 6 | 0.0593011356183722 | 0.0059039660141814 | 9.82452046982027 | True |
| 10 | 7wk_Cebpa_mutant | AT1 | atac | 5 | 0.0026470924908424 | 0.0033527916975189 | 0.6527401043684219 | False |
| 11 | 7wk_Cebpa_mutant | AT1 | rna | 5 | 0.0132768221978795 | 0.0189813016696133 | 0.4912433214400209 | False |
| 12 | SeV_Cebpa_mutant | AT2_identity | atac | 8 | -0.0019135613676306 | 0.0014962201041796 | -1.5780728121142291 | False |
| 13 | SeV_Cebpa_mutant | AT2_identity | rna | 8 | -0.1264921542197687 | 0.0107219161651248 | -11.718758578105184 | True |
| 14 | SeV_Cebpa_mutant | transitional | atac | 6 | 0.0072473099816849 | 0.0012608978937113 | 5.819541563018825 | True |
| 15 | SeV_Cebpa_mutant | transitional | rna | 6 | 0.1193852208147229 | 0.007944309899362 | 15.305546611286704 | True |
| 16 | SeV_Cebpa_mutant | AT1 | atac | 5 | 0.0048046875 | 0.0026122236424466 | 1.8810730453157347 | False |
| 17 | SeV_Cebpa_mutant | AT1 | rna | 5 | 0.0066766666726563 | 0.0091180411601657 | 0.4903245182836303 | False |

## What survives M2 and what does not

**Survives.** The RNA half. In both injured wells the AT2 identity arm
falls far outside the sham band (z of -5.9 and -11.7, a loss of 12.6 to
17.1 detection points) while the AT1 arm does not move at all. M2 repeats
this at four downsampling seeds and it clears in eight of eight.

**Withdrawn.** The chromatin half. M1e read the AT2 chromatin arm's
failure to clear as evidence that the programme stays open, which is
licensed only if the positive control is solid. M2 shows it is not: the
transitional arm clears in five of eight seed-and-well combinations and
fails leave-one-out on two of its six genes in both wells. The AT2 arm's
own null also flips to a clearance in one well when the degenerate median
offset is replaced by a mean.

## The defect this leaves in the record

M1e's rule R9 was written to stop exactly this error and it did not,
because it tested the positive control ONCE, at one seed, on the whole
arm. A sensitivity control that is itself unstable does not license a
null. The rule is not moved; M2 records the correct form, which is that a
positive control has to clear under resampling and under leave-one-out
before any null beside it may be read.
