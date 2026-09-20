# Trials index: the DATP epigenetics branch

Every trial of this branch sits in one flat folder. This page is the index. The
full pre-registration for any trial is the docstring at the top of its script;
the narrative is [`../ANALYSIS_TRIAL_PLAN.md`](../ANALYSIS_TRIAL_PLAN.md) and
the claims table is [`../README.md`](../README.md).

**The M series** works on the two 10x multiome deposits, GSE310539 and
GSE247130. A trailing letter (M1b, M1c, M1d, M1e) marks a disclosed corrected
pass over an earlier trial, never a replacement of it: the first outcome stays
in the record with the rule that produced it.

**A warning before the table.** Six of these eight trials returned NOT
COMPUTABLE, UNREADABLE or Not established, and one had its reading withdrawn. That is the
intended behaviour of the guards, not a failure of the branch: each refusal
names a different defect, and two of them are among the branch's results.

| Trial | Question | What it returned |
|---|---|---|
| [M0](m0_deposit_reality_check.py) | What do the two deposits contain, before anything is fitted | 10 libraries in 4 files, 104,143 cells, a clean peak-to-gene join, and **not one contrast between conditions with within-group replication**. The fifth and sixth deposit in this project with that ceiling |
| [M1](m1_closed_or_merely_silenced.py) | Is the AT2 identity programme closed or merely silenced in the transitional state | **NOT COMPUTABLE**, two rule defects disclosed. But its rule R2 **refuted the GSE247130 suffix map**: Cebpa is 11 to 14 times higher in the suffix the GEO order calls the knockout, so the deposited order is inverted. Confirmed independently on Cldn4 and Sox9 |
| [M1b](m1b_corrected_pass.py) | The same, with the airway rule and the depth handling repaired | **NOT COMPUTABLE in all ten wells.** A 99th-percentile magnitude cut times a five per cent detection call cannot exceed a few tenths of a per cent by arithmetic |
| [M1c](m1c_label_at_the_depth_available.py) | The same, with the label set to what three-prime single-nucleus counting supports | **TRIAL UNREADABLE** under its own negative control, and the reason is a result: the neonatal wells label more heavily than any injured well, so **the marker set cannot separate developmental immaturity from injury-induced transition** |
| [M1d](m1d_mature_wells_only.py) | The same, with the neonatal wells out of the design | **NOT COMPUTABLE in all eight mature wells**: the budget rule retained 60 per cent of every well while the drop gate forbade losing more than 20 per cent. Two rules that were not composable |
| [M1e](m1e_per_well_budget.py) | The same, with the budget set per well so differential drop is zero | The first pass to reach a statistic. Read **silenced but not closed** in two wells. **That reading is withdrawn by M2**; the RNA half stands |
| [M2](m2_robustness_of_the_m1e_reading.py) | Does M1e's reading survive per-gene decomposition, leave-one-out, a non-degenerate offset and a change of seed | **The chromatin half does not.** The positive control clears in 5 of 8 seed-and-well combinations and fails leave-one-out on 2 of 6 genes. The RNA half clears in 8 of 8 |
| [M3](m3_one_instrument_and_the_right_null.py) | Does M2's Not established survive one common instrument and a null over gene sets rather than over cells | **The direction stops being a lean.** Cebpa was the only reason one deposit's AT2 arm was positive and it was already dropped from the other, so the two were never comparable. With one instrument the AT2 arm is negative in all three wells at percentile 0.000, 0.007 and 0.000 of 300 matched gene sets. Still Not established: no well satisfies both gates |

## Shared code

| File | Role |
|---|---|
| [`multiome_utils.py`](multiome_utils.py) | the deposits, the streaming reader for matrices too large to load, the peak-to-gene join, downsampling and the suffix corroboration |
| [`write_summaries.py`](write_summaries.py) | regenerates every trial summary from its tracked CSVs, so a summary can never drift from the numbers it describes |

## Conventions every trial in this folder follows

- Frozen rules are written to a run record **before** any data is read, and the
  record is completed afterwards with inputs, sizes, modification times,
  package versions and outputs.
- A rule is never moved after its result is seen. Where a rule was wrong, the
  first outcome stays and a corrected pass sits beside it.
- Every threshold carries a clause saying what would make its own answer
  unreadable, and a tripped clause reports "not computable", never "refuted".
- The statistical unit is the library. Every well is one library pooling two
  mice, so nothing in this folder is a test and everything is Descriptive only.
- Labels come from the RNA modality only. The labelling genes appear in no
  accessibility arm and their loci are excluded by coordinate.
- A null is readable only behind a positive control that fired, and after M2,
  only behind one that fired under resampling and leave-one-out.
- Regenerable intermediates live under `raw_data/`, which is gitignored. Every
  tracked artefact is a table, a figure or a run record.
