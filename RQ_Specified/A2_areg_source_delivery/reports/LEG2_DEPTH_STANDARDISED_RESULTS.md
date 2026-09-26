# A2 leg 2, second pass: the pair becomes readable, and reads null

27 September 2026. Computed under
[`config/a2_leg2_depth_spec.json`](../config/a2_leg2_depth_spec.json), committed before any
standardised value existed. Exploratory, and explicitly not an independent test: the
unstandardised value was seen first, so a positive here could not have been promoted. Every
number comes from [`tables/`](../tables/) and the
[run record](../tables/leg2_depth_run.json).

## Verdict

**The depth standardisation works, the C51 rule no longer refuses the pair, and the
correlation it was protecting against turns out not to be there.** Most of the nominal
result from the first pass was sequencing depth.

## What changed

| Quantity, pooled across 22 donors | Raw detection fraction | Standardised at the 1,000-molecule budget |
|---|--:|--:|
| Predictor against fibroblast depth | 0.770 | **0.232** |
| Outcome against fibroblast depth | 0.412 | **0.341** |
| Machinery composite against activation | 0.433, nominal p 0.044 | **0.293, nominal p 0.186** |

Both members now sit below the frozen 0.4 threshold, so the pair is readable by the rule
that refused it. Read, it shows no coupling: rho 0.293 at a nominal p of 0.186 in 22 donors,
18 of them IPF. The IPF stratum gives 0.393 at p 0.106. The mural stratum, the population
the cited mechanism actually used, gives 0.029 in six donors and remains underpowered.

The controls behave as the specification required. The machinery composite exceeds both the
detection-matched co-regulation control, at 0.023, and the housekeeping control, at 0.136,
so neither co-regulation nor shared detection is what remains. There is simply little left
to explain.

The abundance comparator reproduces the register. Epithelial AREG against the activation
score gives -0.184 standardised, against the -0.150 that C50 records, so the instrument is
still faithful after the change of measure.

## Why this is evidence that the first result was depth

The correlation tracks how much of each library the measure is allowed to use, in order:

| Measure | Fraction of each library used | Machinery against activation |
|---|---|--:|
| Standardised, 500 molecules | least | 0.230 |
| Standardised, 1,000 molecules (primary) | more | 0.293 |
| Standardised, 2,000 molecules | more still | 0.379 |
| Raw detection fraction | all of it | 0.433 |

The more of each library the measure consumes, the closer the answer returns to the
unstandardised one. That is what a depth-driven correlation looks like when depth is removed
by degrees, and it is the clearest internal evidence available that the frozen rule was
protecting against something real.

## The limitation this pass exposed, which is not small

The specification inherited C37's 1,000-molecule primary budget. The
[measure test](../scripts/test_depth_measure.py) shows, from the formula alone and with no
cohort value, that rate invariance holds only when the budget is small relative to the
smallest library, and degenerates as the budget approaches it, because a cell sampled in its
entirety returns its raw detection. For a gene at a constant rate the standardised value
still spreads by 0.1336 across library sizes at the primary budget, against 0.0086 at a
budget of 100.

**This deposit is floored at exactly 1,000 molecules, so the inherited budget sits on that
degenerate edge.** The standardisation applied here is therefore partial by construction,
which is consistent with the residual outcome coupling of 0.341 that remains. A budget near
100 is what invariance would require in this cohort. That was not substituted into this pass,
because the primary budget was declared before anything ran and changing it after seeing the
problem would be indistinguishable from choosing a budget that reads.

One further caveat on the 2,000-molecule sensitivity: it drops 470 of the 4,648 fibroblasts,
and for typical gene rates the measure is close to its ceiling there, so its apparent
stability is saturation rather than invariance. It is reported because it was declared, not
because it is the better estimate.

## What this establishes and what it does not

**Establishes, as a decision record.** That the C51 depth rule earned its refusal. A
correlation that was nominally significant at p 0.044 on a raw detection fraction falls to
0.293 at p 0.186 once the measure is put on a common molecule budget, and the residual
tracks how much depth the budget leaves in. The rule stopped a reading that would have been
wrong, on a pair declared before it was computed.

**Establishes, descriptively.** That in this cohort, at this measure, the fibroblast
integrin and latent-complex composite does not track the fibroblast activation programme
across donors, and neither does epithelial AREG. Both predictors are null at the donor level.

**Does not establish.** That recipient machinery is unrelated to fibroblast activation. The
standardisation is partial, 22 donors bound what any correlation here could detect, and a
donor-level correlation carries no direction and no proximity in either case.

**Does not close leg 2 for good.** The named next step is a budget near 100 molecules,
declared before it runs, which is the only version of this measure that is rate invariant in
this deposit. If that also reads null, the donor level in this cohort is exhausted for this
question and should be recorded as such.

**Does not touch C49 or C50.** Their grades stand, and this pass reproduced C50 under a
second measure rather than retesting it.

## Proposed claim wording, not graded here

Grading is the owner's decision and no row is added.

1. **Not established.** In GSE136831, a six-gene fibroblast integrin and latent-complex
   composite does not correlate with the frozen five-gene fibroblast activation score across
   22 donors when both are scored as expected detection at a common 1,000-molecule budget
   (rho 0.293, nominal p 0.186); the raw-detection version reached rho 0.433 at p 0.044 but
   was refused by the C51 depth rule.
2. **Validated as a decision record.** The C51 depth rule refused a nominally significant
   donor-level correlation which, measured on a common molecule budget, falls to
   non-significance, and the residual correlation tracks how much of each library the measure
   consumes.
