# A10 stage 3: nested models against measured organoid growth

25 September 2026. The model, its covariates and its 0.02 margin were committed in
`918238b` before any fit ran. All numbers below come from
[`tables/`](../tables/) and their run records.

## Verdict

**Both increments are inconclusive, and neither supports the hypothesis as specified.**

- The **epithelial increment fails the declared margin** in the primary fit and in
  every sensitivity check.
- The **fibroblast increment clears the margin in the primary fit but is not robust.**
  It falls below the margin in two of four checks and is unstable across units.
- A third finding undercuts both: **the baseline does not generalize to a new unit.**
  Pooled predictive performance is largely between-unit spread, not within-unit skill.
- The **assay validation passes convincingly**, so the data and the join are sound.
  The negative result is about the hypothesis, not the pipeline.

## What was fitted

885 wells across 15 plate-replicate units. Outcome is log2 day-14 mean organoid area.
Baseline holds log2 day-7 area, plate, epithelial read fraction and epithelial library
size. Evaluation leaves out one whole unit at a time.

| Model | Held-out R squared |
|---|--:|
| Baseline | 0.586 |
| Baseline plus epithelial block | 0.601 |
| Baseline plus epithelial plus fibroblast block | 0.624 |

| Increment | Value | Margin | Clears |
|---|--:|--:|---|
| Epithelial, primary | 0.0144 | 0.02 | no |
| Fibroblast, second comparison | 0.0230 | 0.02 | yes |

## Why the fibroblast increment cannot be reported as support

| Sensitivity check | Wells | Epithelial | Fibroblast | Fibroblast clears |
|---|--:|--:|--:|---|
| Primary, as fitted | 885 | 0.0144 | 0.0230 | yes |
| Add fibroblast library size to baseline | 885 | 0.0178 | 0.0231 | yes |
| Drop the extreme-replication target | 798 | 0.0163 | 0.0179 | no |
| Drop the control-like libraries | 856 | 0.0113 | 0.0222 | yes |
| Drop both | 769 | 0.0134 | 0.0159 | no |

Stage 1 flagged one target appearing in 87 libraries against four for most, and warned
it would dominate a pooled fit if left unhandled. It does. Removing it takes the
fibroblast increment from 0.023 to 0.018, below the margin, and removing it together
with the controls gives 0.016. An increment that depends on one over-replicated target
is not a programme association.

Per unit, the fibroblast increment is positive in 10 of 15 and ranges from -0.31 to
+0.51. The epithelial increment is positive in 8 of 15 and ranges from -1.08 to +0.58.

## The finding that matters most: the baseline does not transfer

Held-out R squared computed within each held-out unit, rather than pooled:

| Statistic across the 15 units | Baseline R squared |
|---|--:|
| Best unit | 0.401 |
| Median unit | -0.009 |
| Worst unit | -1.698 |
| Units with negative values | 6 of 15 |

A negative value means the model predicts that unit worse than its own mean does. So
for roughly half the units, day-7 area plus plate plus composition carries no usable
information about day-14 area at all. Plate 1 and plate 3 units behave reasonably;
plate 2 and plate 4 units do not.

**The pooled 0.586 is therefore misleading on its own.** It is inflated by differences
between unit means, which the plate terms partly capture, rather than by predicting
which well inside a unit grows. Any increment measured on top of it inherits that
problem. This is the single most important caveat in this report, and it was not
visible until the per-unit breakdown was computed.

One consequence worth stating: a block of nine fibroblast features, each of which
individually adds at most nothing, together adds 0.023. In the exploratory
single-feature scan the largest gain from any one feature is 0.0038, and every
fibroblast feature alone is negative. A block that helps only jointly, while none of
its parts helps alone, is more consistent with fitting unit-level structure than with
a coherent niche programme.

## The assay validation passes

Genes that are both module members and perturbation targets should fall in the
libraries where they were targeted.

| Check | Value |
|---|--:|
| Genes testable | 76 |
| Reduced when targeted | 58 |
| One-sided binomial p | 0.0000024 |

This supports two things at once: editing reduces the targeted transcript, and the
sheet labelled mouse is the perturbed compartment. It closes the stage 1 hold on which
assigned genome is which, now from data rather than from the sheet name alone.

So the negative result is not a broken pipeline. The joins hold, the species split is
confirmed, and the perturbations worked.

## What this does and does not establish

**Does not establish.** That epithelial programme state carries useful information
about organoid growth beyond baseline and batch. The increment is below the margin.

**Does not retire the hypothesis either.** The pre-declared retirement rule requires a
confidence bound excluding the margin, which was not computed and would not be credible
given how poorly the baseline transfers. This is inconclusive, not precise absence.

**Does establish, descriptively.** The frozen modules from the shared component
contract, scored in this assay, do not add usefully to a size-and-batch baseline for
predicting organoid growth in a held-out unit. That is a real if unwelcome fact about
those modules in this assay, and it is a fact about nine specific frozen modules, not
about all possible epithelial programmes.

## Why no stronger reading is available

The units are batches, not verified independent preparations, so even a robust
increment would have been a within-screen association. The
[public data search](PUBLIC_DATA_SEARCH.md) found no cohort that could validate this
screen externally. Expression and outcome are both from day 14, so nothing here bears
on prediction. Organoid area is growth, not mature cell fate or repair in a lung.

## What would change the answer

1. **Resolve the unit.** If replicate indices are independent preparations, the
   holdout means something; if not, this analysis was always descriptive. The deposit
   does not say, and that remains an owner decision.
2. **Model the unit heterogeneity.** The per-unit failures suggest plate and replicate
   effects are not additive. A specification allowing unit-varying slopes on day-7
   area would be a better baseline, and it should be declared before being fitted.
3. **Use a programme set chosen for this assay.** The frozen modules came from repair,
   development and lesion biology. Their failure here is informative but it is not a
   test of programmes that a screen of this kind would nominate.
4. **Handle the over-replicated target explicitly** in any future specification, rather
   than as a sensitivity check.

## Provenance

The margin and model were committed before fitting. The extraction streamed both count
sheets and kept only the declared gene set, with 1,575 of 1,613 mouse and 1,442 of
1,449 human genes found. The fit refuses to run unless the cached counts match the
extraction record. Sensitivity and per-unit tables come from a separate script that
refuses to overwrite and does not alter the primary.
