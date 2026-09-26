# A2 stage 1: instrument, eligibility and power

26 September 2026. Every number below comes from
[`tables/`](../tables/) and the [run record](../tables/stage1_run.json). The script
computed **no endpoint**: the human pass collected gene symbols only and read no
count, and the mouse pass read library totals and seven named gene rows.

## Verdict

**The design can carry the declared test, with three qualifications that stage 2
must absorb.** All seven stop rules pass. The qualifications are that one receptor
target is not expressed, that fibroblast sequencing depth on this plate spans four
orders of magnitude, and that the four Areg wells sit above their unit median in
that depth.

## Item 1: the plate-3 layout is as the contract states

The three cached inputs match the hashes recorded when they were downloaded, and
the 886 libraries in the deposited QC table give:

| Fact | Value |
|---|--:|
| Plate-3 units | 4 |
| Wells in every plate-3 unit | 60 |
| Wells per axis target per unit | 1 |
| Control wells per unit | 2 |
| Axis targets outside plate 3 | 0 |

Hbegf is on plate 4 at one well per unit, so the ligand-specificity contrast
crosses plates, as the plan already states. Full table:
[`stage1_layout.tsv`](../tables/stage1_layout.tsv).

## Item 2: the knockouts validate, except that one receptor is not there

The script reproduces A10's metric exactly before extending it. All 886 recomputed
mouse library totals agree with the A10 extraction record, and the five genes A10
already tested reproduce to three decimals.

| Target | Mean log2 CPM, targeted | Mean log2 CPM, other | Difference |
|---|--:|--:|--:|
| Areg | 6.226 | 7.267 | -1.042 |
| Egfr | 0.933 | 1.967 | -1.034 |
| Erbb2 | 3.147 | 5.589 | -2.441 |
| Erbb3 | 5.718 | 6.229 | -0.511 |
| Erbb4 | 0.000 | 0.367 | -0.367 |
| Itgb6 | 7.043 | 8.181 | -1.138 |
| Hbegf | 5.734 | 7.656 | -1.921 |

**Erbb4 is not expressed in this compartment.** Its non-targeted mean is 0.367
log2 CPM and its targeted mean is exactly zero, meaning no counts at all in the
four targeted libraries. There is no receptor to remove, so an Erbb4 contrast
cannot discriminate anything and stage 2 drops it.

**Egfr is near the detection floor** at 1.967 log2 CPM when not targeted, about
2.9 counts per million. Its knockout validates, but the contrast is weaker than
Erbb2 or Erbb3 and must be labelled so.

Erbb2 and Erbb4 had no recorded validation before this stage, because A10 tested
only genes that were also members of its module set. Both now have one, on the
same footing. Table: [`stage1_knockout_axis.tsv`](../tables/stage1_knockout_axis.tsv).

## Item 3: both endpoints are fully covered

| Endpoint | Members in the human index | Floor | Clears |
|---|--:|--:|---|
| HALLMARK_TGF_BETA_SIGNALING | 54 of 54 | 0.90 | yes |
| Fibroblast activation score | 5 of 5 | 1.00 | yes |

The human sheet carries 56,648 unique symbols across 58,302 rows, so the workbook
is genome-wide on the fibroblast side and neither endpoint loses a member.
Table: [`stage1_endpoint_coverage.tsv`](../tables/stage1_endpoint_coverage.tsv).

## Item 4: the covariate audit is the most consequential result

Fibroblast total counts across the 240 plate-3 wells:

| Quantile | Fibroblast counts |
|---|--:|
| Minimum | 1 |
| 1st percentile | 1,651 |
| 5th percentile | 10,674 |
| Median | 277,755 |
| Maximum | 4,409,390 |

Two facts follow, and both are covariate facts, established with no endpoint in
hand.

**The four Areg wells are all above their unit median in fibroblast depth**, at
unit percentiles 0.683, 0.700, 0.817 and 0.883. A mean log2 CPM score over a gene
set rises with depth, because an undetected gene contributes zero and deeper
libraries detect more of the set. The expected direction of that artefact is
therefore upward, which is the opposite of the declared prediction, so the
one-sided test is made conservative rather than permissive by this asymmetry. That
reasoning is an expectation, not a measurement, so stage 3 must report the
within-unit association between the endpoint and log depth next to the result.

**Three single-well contrasts sit too low to carry a 54-gene score.** Erbb2 in
unit 3-1 has 1,759 fibroblast counts, Itgb6 in 3-3 has 4,865, and Egfr in 3-2 has
51,605. At 25,000 counts a gene at 50 counts per million is expected about once,
so such a well reports detection, not programme state.

**Four of the eight plate-3 control wells are also shallow**, at 10,231, 11,417,
14,219 and 17,021 fibroblast counts. The contract gave the control wells the role
of anchoring the unperturbed location; on this plate that anchor is weak, and
stage 2 downgrades it. Table:
[`stage1_covariates.tsv`](../tables/stage1_covariates.tsv).

## Item 5: the test is still blind

No tracked table in the repository pairs either declared endpoint with per-target
or per-library rows. A10's per-well score table carries nine frozen fibroblast
modules and three Hallmark control sets, and none of them is the declared primary.
The inspected files and their hashes are in
[`stage1_precedent.json`](../tables/stage1_precedent.json).

## Item 6: what the design can and cannot reach

Exact null for the declared statistic, ranks uniform within each unit:

| Quantity | Value |
|---|--:|
| Rank space, 4 units of 60 wells | 12,960,000 |
| Smallest attainable one-sided p | 7.72e-08 |
| Critical rank sum at alpha 0.05 | 64 |
| Mean rank that requires | 16 of 60 |
| Mean percentile that requires | 0.267 |
| Exact size of the joint rule | 0.04903 |

So the design is not underpowered in the usual sense; it is demanding in a
specific way. The Areg well must average the 27th percentile of its unit or lower.
Consistency alone cannot carry it: being just below the median in all four units
does not reach the threshold. The joint rule that adds the three-of-four
consistency requirement has an exact size of 0.049, so it does not inflate the
test. Record: [`stage1_power.json`](../tables/stage1_power.json).

## Item 7: the cross-species assumption is carried, not settled

No sequence resource exists under the local `raw_data` tree, and the plan confines
this stage to local resources, so no fetch was made. The assumption that the mouse
amphiregulin EGF-like domain activates human EGFR is recorded as carried and
unverified. It can be settled by aligning the mature domains of UniProt P31955 and
P15514, which remains open work.

## What stage 2 must absorb

1. **Drop Erbb4** from the discriminating set, for want of an expressed receptor.
2. **Label Egfr as low-abundance**, so its contrast is not read as equal to Erbb2
   or Erbb3.
3. **Declare a fibroblast depth floor** for reading any single well, and declare
   depth-restricted sensitivities for the primary, because the Areg wells are
   systematically deeper than their units and three contrast wells are too shallow
   to read at all.
4. **Downgrade the control wells' anchor role**, since half of them are shallow.

None of these changes the hypothesis, the endpoints or the declared direction, and
all four follow from covariates rather than from any endpoint value.
