# Trial M3: one instrument, two nulls, and a direction that will not go away

M2 withdrew M1e's chromatin reading and left the question Not established.
A five-lens adversarial review of that withdrawal then found three defects
in M2 itself, and M3 implements their corrections. **This repository did
not find them**; the review did, and the run record says so.

## The three corrections

**D1, the two deposits were not running the same instrument.** Rule R12
removed Cebpa from the AT2 arm where it is genetically deleted and kept it
everywhere else. Cebpa carries only three distal peaks but the largest
positive per-gene value in wildtype_SeV, and it is the only reason that
well's AT2 arm was positive at all. M2's own leave-one-out had already
said so: dropping Cebpa gave the most negative of the nine drops. M3 drops
it everywhere, leaving one eight-gene instrument.

**D2, the sham band is the wrong null for the question.** It permutes
CELLS, so it asks whether this SPLIT is special. It cannot ask whether
these GENES are special, which is what a gene-set arm claims. M3 adds 300
random gene sets matched gene by gene on distal peak count, evaluated on
the real split, and uses their mean as the offset.

**D3, the equivalence statistic erred unsafely.** Three sham standard
deviations is the 50-per-cent-power detection floor, not a bound the data
support. M3 reports a confidence interval on the relative scale instead.

## What one instrument and the right null return

|  | well | role | arm | genes | distal_peaks | observed_raw | geneset_null_mean | geneset_null_sd | geneset_percentile | sham_sd | z_sham | z_geneset | clears_both | reference_baseline | relative_pct | relative_lo_pct | relative_hi_pct |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | wildtype_SeV | test | AT2_identity | 8 | 126 | -0.00363 | 0.00069 | 0.00158 | 0.0 | 0.002 | -2.16114 | -2.72504 | False | 0.01141 | -37.82075 | -72.15575 | -3.48575 |
| 1 | wildtype_SeV | test | transitional | 6 | 50 | 0.00971 | 0.00038 | 0.00306 | 0.99333 | 0.00269 | 3.46133 | 3.04531 | True | 0.0151 | 61.77448 | 26.75928 | 96.78967 |
| 2 | wildtype_SeV | test | AT1 | 5 | 27 | 0.0032 | 0.00022 | 0.0038 | 0.83333 | 0.00325 | 0.91452 | 0.78341 | False | 0.00779 | 38.2213 | -43.77621 | 120.21882 |
| 3 | 7wk_Cebpa_mutant | negative control | AT2_identity | 8 | 114 | -0.00206 | 0.00029 | 0.00106 | 0.00667 | 0.00151 | -1.56001 | -2.21564 | False | 0.01288 | -18.29001 | -41.29255 | 4.71252 |
| 4 | 7wk_Cebpa_mutant | negative control | transitional | 6 | 38 | 0.00238 | 0.00035 | 0.00195 | 0.88 | 0.00228 | 0.89003 | 1.04221 | False | 0.01308 | 15.53388 | -18.70855 | 49.7763 |
| 5 | 7wk_Cebpa_mutant | negative control | AT1 | 5 | 25 | 0.00611 | 0.00027 | 0.00259 | 0.97667 | 0.00349 | 1.67606 | 2.25742 | False | 0.0106 | 55.15891 | -9.40898 | 119.7268 |
| 6 | SeV_Cebpa_mutant | test | AT2_identity | 8 | 152 | -0.00247 | 0.00205 | 0.00129 | 0.0 | 0.00133 | -3.40635 | -3.50618 | True | 0.0292 | -15.48975 | -24.4114 | -6.56811 |
| 7 | SeV_Cebpa_mutant | test | transitional | 6 | 56 | 0.00515 | 0.00182 | 0.00198 | 0.96667 | 0.0017 | 1.95969 | 1.68349 | False | 0.03483 | 9.56479 | -0.01105 | 19.14063 |
| 8 | SeV_Cebpa_mutant | test | AT1 | 5 | 34 | -0.00196 | 0.00195 | 0.00297 | 0.05667 | 0.00348 | -1.12256 | -1.31612 | False | 0.0287 | -13.61636 | -37.41434 | 10.18163 |

## The reading, and it is still Not established

|  | well | role | n_labelled | reading | why | at2_relative_pct | at2_ci_lo | at2_ci_hi | at2_geneset_percentile |
|---|---|---|---|---|---|---|---|---|---|
| 0 | wildtype_SeV | test | 64 | not established | R19: the AT2 arm clears at most one null | -37.821 | -72.156 | -3.486 | 0.0 |
| 1 | 7wk_Cebpa_mutant | negative control | 112 | not computable | R9: the transitional arm does not clear both nulls | -18.29 | -41.293 | 4.713 | 0.007 |
| 2 | SeV_Cebpa_mutant | test | 320 | not computable | R9: the transitional arm does not clear both nulls | -15.49 | -24.411 | -6.568 | 0.0 |

**The AT2 identity arm is negative in all three computable wells, at
percentile 0.000, 0.007 and 0.000 of 300 matched random gene sets.** On
the relative scale that is -37.8, -18.3 and -15.5 per cent of the
reference group's own accessibility at those loci. The direction is the
same everywhere and it is the opposite of what M1e reported.

**But no well satisfies both gates, and the two wells fail different
ones.** wildtype_SeV has a positive control that clears (z 3.46 and 3.05)
and an AT2 arm that does not (z -2.16 and -2.73, short of 3).
SeV_Cebpa_mutant has an AT2 arm that clears both nulls at z -3.41 and
-3.51, in the closing direction, behind a positive control that does not
clear (z 1.96 and 1.68), so rule R9 refuses the well. Each well holds half
of what a reading needs and neither holds both. The chromatin question
stays **Not established**, now with a direction and a bound rather than
without one.

## The offset was doing more work than anyone thought

Re-centred on the gene-set null, the transitional arm in SeV_Cebpa_mutant
falls from clearing under M1e's degenerate median offset to z = 1.96.
**Most of what M1e counted as signal in that well was the uncorrected
global shift**, which is exactly what the review predicted when it said
the offset was absent rather than inert. Only wildtype_SeV has a
transitional signal that survives proper centring, at +61.8 per cent of
reference with an interval of +26.8 to +96.8.

## The one bound worth quoting

In SeV_Cebpa_mutant, AT2 distal accessibility in the labelled group is
**-15.5 per cent of reference accessibility, interval -24.4 to -6.6**, at
percentile 0 of 300 matched gene sets. It is not readable as a finding,
because its well's positive control did not fire, and the frozen rules
refuse it. It is recorded because a bound that a rule refuses is still a
number, and the next instrument should be sized against it.

## The pre-registered predictions, checked

The review predicted values before M3 ran and R22 froze them. The gene-set
null means were predicted almost exactly: +0.0009 and +0.0021 against
+0.00069 and +0.00205 observed. The AT2 direction and the percentile-zero
result were confirmed and came out stronger than predicted (z -3.51
against -2.80). The transitional arm came out **weaker** than predicted,
+3.05 and +1.68 against +3.93 and +4.46, which is why R9 now refuses a
well the review expected to pass. Agreement here confirms the review's
arithmetic; it is not an independent discovery by this repository.

## What this does not do

It does not resurrect M1e's reading, and it could not have: with one
instrument the arm is negative under every weighting. Nothing in this
folder supports the AT2 programme staying open in the transitional state.
