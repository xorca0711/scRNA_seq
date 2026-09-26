# A2 stage 4, leg 2: refused by its own control

26 September 2026. Exploratory by declaration, under
[`config/a2_leg2_spec.json`](../config/a2_leg2_spec.json), which was committed before
this ran. Numbers come from [`tables/`](../tables/) and the
[run record](../tables/stage4_run.json).

## Verdict

**The pre-declared depth rule refuses the pair, and the refusal is the result.** The
specification said in advance that this was the likely outcome and why, so this is a
recorded prediction coming true, not a disappointment discovered afterwards.

## What was computed

One pass over the deposited sparse matrix recovered 209 of 217 requested gene vectors,
including every declared gene. Twenty-two donors carry both compartments at the 50-cell
floor, 18 of them IPF and 4 Control, with 52 to 608 fibroblasts each. All six machinery
members cleared the 0.05 detection floor, so none was dropped.

| Stratum | Pair | Donors | rho | Nominal p |
|---|---|--:|--:|--:|
| Pooled | machinery composite against activation | 22 | +0.433 | 0.044 |
| Pooled | epithelial AREG against activation | 22 | -0.150 | 0.506 |
| Pooled | co-regulation control against activation | 22 | +0.208 | 0.352 |
| Pooled | housekeeping control against activation | 22 | +0.066 | 0.770 |
| Pooled | **machinery against fibroblast depth** | 22 | **+0.770** | 0.00003 |
| Pooled | **activation against fibroblast depth** | 22 | **+0.412** | 0.057 |
| IPF only | machinery against activation | 18 | +0.453 | 0.059 |
| IPF only | machinery against fibroblast depth | 18 | +0.717 | 0.0008 |
| IPF only | activation against fibroblast depth | 18 | +0.321 | 0.194 |
| Mural | machinery against activation | 6 | -0.143 | 0.787 |

## Why the pair is refused

The inherited C51 rule refuses a pair when both members track their own compartment's
depth at an absolute rho of 0.4 or more. The predictor reaches 0.770 and the outcome
0.412, so the rule fires on the pooled set. The nominal +0.433 between them is therefore
not read, whatever it looks like.

The reason is visible in the measure. A detection fraction computed from 52 to 608 cells
per donor, over genes detected in 22 to 72 per cent of cells, across a 2.1-fold spread in
median genes per cell, is largely a statement about sequencing depth. A six-gene
composite of such fractions is more depth-coupled than any single member, which is what
the 0.770 says.

**The IPF stratum is also not read, and that is a judgement stricter than the declared
rule.** By the letter of the rule the IPF pair survives, because the outcome's depth
coupling falls to 0.321 there. But the predictor sits at 0.717 in that stratum, so the
correlation remains a depth comparison. Refusing it is a post hoc decision, made after
seeing the numbers, and it is recorded as such; it errs in the conservative direction,
which is the only direction in which a post hoc tightening is acceptable.

## Two things the leg did establish

**The instrument reproduces.** Epithelial AREG against the fibroblast activation score
gives -0.150 here, against the -0.150 the register records for C50. The reuse of the
trial E6 instrument is therefore faithful, which is what makes the refusal credible
rather than a sign of a broken pipeline.

**The co-regulation rival was addressed, and it was not the problem.** The matched
control was built by the declared rule, from oxidative-phosphorylation members outside
the TGF-beta and mesenchymal-transition sets, matched on pooled fibroblast detection
fraction. The matching is close: ITGAV 0.549 against SLC25A3 0.547, ITGB1 0.552 against
NDUFA1 0.554, ITGB8 0.224 against ATP6V1D 0.224, LTBP1 0.718 against COX4I1 0.716, TGFB1
0.289 against BAX 0.287, THBS1 0.516 against CYB5R3 0.513. The machinery composite
(+0.433) exceeds both that control (+0.208) and the housekeeping control (+0.066), so a
co-regulation or shared-detection explanation is not what defeats this leg. The depth
confound is.

## The population the mechanism actually used

The cited work acted on PDGFRB-positive pericytes. In this deposit the mural stratum,
Pericyte plus SMC, reaches the 50-cell floor in only 6 donors, with a median of 18 mural
cells per donor. Its correlation is -0.143. That stratum is underpowered and is reported
as underpowered; nothing follows from it either way.

## What would make this leg readable

The repository already owns the fix and used it for C37: put the measure on a common
molecule budget instead of a raw detection fraction. Concretely, a depth-standardized
expectation per donor, or a pseudobulk expression measure rather than a per-cell
detection fraction, with the budget declared before any correlation is computed. Until
then the donor level in this cohort cannot separate recipient state from sequencing
depth, and that is a property of the measure rather than of the biology.

## What this does and does not establish

**Establishes.** Nothing about the biology. The declared reading for a refused pair is
that the refusal is the result.

**Does establish, as a decision record.** That the C51 depth rule, applied to a new pair
in the same cohort, refuses it, and that the rule was written down before the pair was
computed. That is the second time this rule has stopped a reading in this repository, and
the first time it stopped one that would have been nominally significant at p 0.044.

**Does not establish.** That recipient machinery is unrelated to fibroblast activation.
The measure cannot see it here, and no bound on the effect is available from a refused
pair.

**Does not touch C49 or C50.** Their grades stand, and this leg reproduced C50 rather
than retesting it.
