# A10 stage 4: the revised specification, and what it changes

**Interpretation addendum, 26 September 2026.** The 0.0426 epithelial increment
and -0.0154 fibroblast increment are retained, as is the declared descriptive
margin decision. The [logical review](../../../docs/LOGICAL_RATIONALE_REVIEW.md#a10-corrections)
qualifies the original narrative below: this is a same-screen adaptive
specification using held-out group outcome means, not independent confirmation.
The four-cell pattern is specification dependence, not a tested biological
interaction. A fixed R-squared margin is not automatically conservative across
different outcome scales, and no extra fibroblast score information does not
exclude niche function. Secondary outputs have no implemented BH inference and
use both primary blocks as their baseline. The original narrative and numbers
below are retained as the dated reporting record; this addendum governs current
interpretation. Proposed follow-ups now follow the review's diagnostic-first order.

## Original report, 25 September 2026

25 September 2026. The specification, its blocks and its 0.02 margin were committed in
`165a4f5` and the scripts in `8084083`, both before any fit. Numbers come from
[`tables/`](../tables/) and their run records.

## Verdict

**The epithelial growth programmes clear the declared margin, but only with the
within-unit baseline, and unevenly across plates. The fibroblast block adds nothing.**

The four-cell grid was declared precisely so the two changes could be told apart, and
it earns its keep: neither change works alone.

| Programme set | Baseline | Baseline R squared | Epithelial | Fibroblast |
|---|---|--:|--:|--:|
| Frozen repair modules | pooled | 0.586 | 0.0144 | 0.0230 |
| Frozen repair modules | within-unit | 0.377 | 0.0036 | 0.0161 |
| Growth Hallmark | pooled | 0.586 | 0.0136 | -0.0017 |
| **Growth Hallmark** | **within-unit** | **0.377** | **0.0426** | **-0.0154** |

Only the last cell clears the 0.02 margin, at more than double it. The growth
programmes with the old baseline give 0.0136, and the old modules with the new baseline
give 0.0036. **Both changes were necessary, and that is an interaction, not an
additive improvement.** It is also biologically coherent: removing between-unit variance
lets a within-unit signal show, and proliferation, cell loss and biosynthetic signalling
are the right programmes for a growth endpoint, which the repair modules never contained.

## The supported result passes every mandatory check

| Check | Wells | Epithelial | Clears margin |
|---|--:|--:|---|
| Primary | 885 | 0.0426 | yes |
| Drop the over-replicated target | 798 | 0.0414 | yes |
| Drop the control-like libraries | 856 | 0.0485 | yes |
| Drop both | 769 | 0.0460 | yes |

The declared robustness rule required the margin to be cleared in the primary fit and
after dropping both. It is, in all four. This is the rule the earlier fibroblast
increment failed, so it is not being applied leniently here.

## The earlier fibroblast finding was a specification artefact

The first analysis reported a fibroblast increment of 0.0230 that cleared the margin,
and I described the niche half as the more interesting one. The grid shows that result
does not survive either change:

| Condition | Fibroblast increment |
|---|--:|
| Old modules, pooled baseline, as first reported | 0.0230 |
| Old modules, within-unit baseline | 0.0161 |
| Growth programmes, pooled baseline | -0.0017 |
| Growth programmes, within-unit baseline | -0.0154 |

With the better specification the fibroblast block is consistently negative, meaning it
costs held-out accuracy. **So the direction of my earlier framing was wrong.** The
epithelial side carries the information about growth, and the niche block does not add
once the epithelial state and the baseline are in the model. That does not mean the
niche is irrelevant to the biology. It means these fibroblast programme scores, in bulk
species-split RNA from the same well, add nothing measurable beyond the epithelium.

## The limitation that qualifies the result

The within-unit transform did not fix the per-unit heterogeneity, it only stopped it
being hidden.

| Per-unit statistic, supported cell | Value |
|---|--:|
| Units | 15 |
| Units where the epithelial increment is positive | 8 |
| Units where the baseline R squared is negative | 7 |
| Worst unit increment | -0.45 |
| Best unit increment | +0.13 |

Plate 1 and plate 3 units behave well, with positive baselines and positive increments
almost throughout. Plate 2 and plate 4 units still have negative baseline R squared,
and two of them have strongly negative increments.

**So the pooled result is carried by two of the four plates.** It satisfies the
pre-declared rule, and I am not overturning that rule after seeing the numbers. But an
association concentrated in half the plates is weaker than the headline suggests, and
whether those plates differ biologically or technically cannot be answered from the
deposit.

## Secondary blocks: nothing added

None of the four declared secondary blocks clears the margin, in either baseline.

| Block | Pooled | Within-unit |
|---|--:|--:|
| Epithelial niche signalling | -0.0084 | -0.0081 |
| Epithelial environment | 0.0135 | 0.0030 |
| Fibroblast niche signalling | 0.0059 | -0.0052 |
| Fibroblast environment | 0.0099 | 0.0163 |

Notably the Wnt, TGF-beta, Notch and Hedgehog block does not add in either compartment,
despite Wnt being the best-established route by which fibroblasts maintain type 2
stemness. In this assay and at this resolution, that pathway block carries no extra
information about growth beyond proliferation and cell loss.

## What is now established, and what is not

**Established, descriptively.** In this screen, epithelial proliferation, cell-loss and
biosynthetic-signalling programmes carry information about which well inside a batch
grows more, beyond starting size and composition. The increment is 0.0426 and survives
every declared sensitivity check.

**Not established.** That this generalizes to an unseen batch. The metric centres each
held-out unit on its own mean, which was declared in advance, so it answers a
within-batch question. Half the units have a baseline that does not predict at all.

**Not established.** Any role for fibroblast programmes. The evidence now points the
other way within this measurement.

**Unchanged limits.** The biological unit is still unresolved, so this remains a
within-screen descriptive association. Expression and outcome are both from day 14, so
nothing here concerns prediction. Organoid area is growth, not mature cell fate.

## Why the first analysis was not wrong to report

The first report said both increments were inconclusive and named the per-unit transfer
failure as the reason. That diagnosis was correct and it is what pointed at the right
fix. Its numbers stand unchanged, and its outputs were not overwritten. The revision
found a signal the first specification could not have seen, because the first programme
set contained no proliferation and the first baseline let batch variance dominate.

## What would strengthen this next

1. **Resolve why plates 2 and 4 behave differently.** This is now the highest-value
   question about the data, above any further modelling.
2. **Resolve the biological unit.** Without it the result cannot leave descriptive
   status, whatever its size.
3. **Test the epithelial result against a held-out plate rather than a held-out unit.**
   Plate-level holdout is a harder and more honest test, and it should be declared
   before it is run.
4. **Ask whether the signal is simply proliferation.** A single-feature breakdown of the
   growth block would show whether the increment reduces to a cell-cycle score, which
   would be a plainer and more testable claim.
