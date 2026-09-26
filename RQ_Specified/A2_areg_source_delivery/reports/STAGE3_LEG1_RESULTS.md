# A2 stage 3, leg 1: the epithelial source contribution, measured

26 September 2026. Computed under
[`config/a2_stage2_freeze_v2.json`](../config/a2_stage2_freeze_v2.json), which
forbids a p-value. Every number comes from [`tables/`](../tables/) and the
[run record](../tables/stage3_run.json). The 886 recomputed human library totals
agree with the A10 extraction record, and both endpoints recovered every member,
5 of 5 and 54 of 54.

## Verdict

**The prediction fails for Areg, and the one large consistent effect in the panel
belongs to Itgb6.** Both statements are descriptive, from one well per target per
unit in four copies of a single plate layout, and neither carries a decision.

## Scale, so the effects can be read

The primary endpoint, the five-gene fibroblast activation score, has a standard
deviation of 0.482 log2 CPM across the 168 eligible plate-3 wells, with a median of
9.522 and an interquartile range of 9.227 to 9.835. The secondary Hallmark TGF-beta
score has a standard deviation of 0.271.

## The Areg contrast, against the screen's own controls

Effects are the Areg well's score minus the mean of its unit's eligible control
wells, TIGIT plus TDTOMATO, restricted to controls within a factor of four of the
Areg well's fibroblast depth.

| Unit | Depth-matched effect, activation score |
|---|--:|
| 3-1 | +0.131 |
| 3-2 | +0.414 |
| 3-3 | -0.064 |
| 3-4 | -0.059 |

| Summary, Areg | Activation score | Hallmark TGF-beta |
|---|--:|--:|
| Median effect against controls | +0.034 | +0.019 |
| Median depth-matched | +0.036 | +0.028 |
| Units lower than controls | 2 of 4 | 1 of 4 |
| Median after epithelial adjustment | +0.057 | not lower |

The predicted direction was lower. The median is about 0.07 standard deviations in
the wrong direction, two units of four sit below their controls, and one unit is
strongly positive. **No epithelial source contribution is detected on either
endpoint.**

## What that null does not mean

The freeze declared before this ran that no null here can support absence, and three
recorded facts make that binding rather than formulaic.

1. **The recipient supplies the ligand.** These fibroblasts transcribe AREG at 9.846
   mean log2 CPM in 99.2 per cent of plate-3 wells, above the mouse epithelial Areg
   at 7.267 that the knockout removes. An unremoved autocrine source at that level
   can mask an epithelial contribution entirely.
2. **The knockout is partial.** The mouse Areg transcript falls by 1.042 log2 CPM and
   remains at 6.226.
3. **The culture medium is not in the deposit.** Exogenous EGF would saturate the
   recipient's receptor and a TGF-beta receptor inhibitor would clamp the endpoint.

So the reading is: **no detectable epithelial AREG contribution on top of an
unremoved autocrine source, in this assay, with this partial perturbation.** Nothing
weaker and nothing stronger.

## The result that was not predicted: epithelial Itgb6

| Target | Median depth-matched effect | Units lower | Interpretable units | After epithelial adjustment |
|---|--:|---|---|--:|
| **Itgb6** | **-0.938** | 3 of 3 | 3 of 4 | -0.992 |
| Erbb3 | -0.187 | 4 of 4 | 4 of 4 | +0.101 |
| Areg | +0.036 | 2 of 4 | 4 of 4 | +0.057 |
| Egfr | +0.247 | 0 of 3 | 3 of 4 | +0.284 |
| Erbb4 | +0.503 | 0 of 3 | 3 of 4 | +0.385 |
| Erbb2 | -0.169 | 1 of 1 | 1 of 4 | -0.130 |

Per unit, Itgb6 gives -0.938, -1.575 and -0.680, with its fourth well ineligible at
4,865 fibroblast counts. That is about 1.9 standard deviations of the endpoint, it
holds in every unit where the well is readable, and it survives the
epithelial-fraction adjustment at -0.992, so it is not a composition artefact.

**Integrin beta-6 pairs with integrin alphaV on epithelium to activate latent
TGF-beta.** Removing it from the epithelium lowers the human fibroblast
myofibroblast and collagen programme. That is the epithelial-activation rival the
card named, and it is the arm that moved, while the ligand arm did not.

Two cautions belong with it, and they are not decorative. Erbb3 also falls in all
four units, and ERBB3 binds neuregulins rather than AREG, so part of any epithelial
perturbation effect may be generic; Itgb6's effect is about five times larger.
Egfr and Erbb4 move upward, so the panel is not uniformly downward. And one well per
target per unit in four copies of one layout cannot separate a target from its plate
position.

## The depth diagnostic the freeze required

Within the eligible wells the endpoint is essentially uncorrelated with depth, which
is what the eligibility floor was for.

| Unit | Endpoint against log depth | Endpoint against epithelial fraction | Eligible wells |
|---|--:|--:|--:|
| 3-1 | -0.049 | +0.069 | 32 |
| 3-2 | +0.043 | +0.134 | 48 |
| 3-3 | -0.029 | -0.056 | 31 |
| 3-4 | -0.109 | +0.003 | 57 |

The depth asymmetry that worried the first freeze is therefore handled by the floor
rather than by the adjustment, and no effect above rests on a depth gradient.

## What this establishes and what it proposes

**Establishes, descriptively and within this screen only.** The frozen five-gene
fibroblast activation score does not move when epithelial Areg is knocked out in this
co-culture, against the screen's own controls, in four units of one plate layout.

**Does not establish.** That epithelial AREG has no role. The recipient supplies the
ligand, the perturbation is partial, the medium is unknown, and preparation
independence is unresolved.

**Proposes, as a hypothesis and not a result.** That the epithelial contribution to
this fibroblast programme runs through epithelial integrin-mediated activation of
latent TGF-beta rather than through amphiregulin. Testing it needs a design this
deposit cannot supply: more than one well per target per unit, positions that vary,
an epithelial integrin perturbation with its own validation, and a TGF-beta
activation readout rather than a transcript proxy.

No claim row is proposed as validated, no register grade changes, and the register's
A2 readiness stays conditional and descriptive.
