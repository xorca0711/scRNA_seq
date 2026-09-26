# A15 rival-2 side-branch: what the GSE190821 epithelial arm shows

27 September 2026. Computed under
[`config/a15_rival2_freeze_v3.json`](../config/a15_rival2_freeze_v3.json), which was
committed before any value was read and which is the third version of this freeze; the
[first](../config/a15_rival2_freeze.json) and [second](../config/a15_rival2_freeze_v2.json)
were withdrawn on adversarial review with nothing scored. Every number comes from
[`tables/rival2/`](../tables/rival2/) and the
[run record](../tables/rival2/stage3_execute_run.json). An
[independent recomputation](../tables/rival2/stage3_verification.json) through numpy and
scipy passes 39 of 39 checks.

**This bounds A15 rival 2 only. It is not a test of A15 and may not be cited as evidence for
or against it.** The deposit fails A15 gate conditions 1, 3, 4 and 5.

## The result in one paragraph

The 3G9 anti-integrin-beta6 antibody produced a large, completely separating reduction in
the whole-lung collagen and myofibroblast programme, reproducing the published direction for
this antibody. In the **same eight mice**, with the epithelial immunoprecipitation shown to
be pure and comparable between arms, the epithelium's own identity and transitional
programme showed no detectable difference, and neither did a 15,062-gene omnibus profile.
Under the frozen rules that is a **weak bound on rival 2**: the epithelial change that rival 2
needs in order to mediate anything was not detected while the antibody was demonstrably
doing something. At four mice against four it is a weak bound and nothing stronger.

## The engagement control worked, and that is reproduction

The frozen control is the mouse orthologs of this repository's five-gene fibroblast
activation score, `Col1a1`, `Acta2`, `Postn`, `Cthrc1`, `Tnc`, measured on the paired
whole-lung **input** library and declared before any value was read to fall under 3G9,
because that is the published effect of this antibody
([Horan 2008](https://doi.org/10.1164/rccm.200706-805OC)).

| Arm | Per-mouse standardised score |
|---|---|
| 3G9 | -1.0233, -0.3892, -0.4731, -0.5992 |
| Axum8 | +0.8272, -0.1937, +0.5364, +1.3149 |

Mean difference **-1.2424**, exact two-sided p **0.028571**, which is the floor at four
against four and therefore **complete separation**: every treated mouse sits below every
control mouse. The Hodges-Lehmann shift is **-1.2583** with a 97.14 per cent
distribution-free interval of **[-2.3382, -0.1955]**, which excludes zero. The same score
separates bleomycin from saline at +1.3699, p 0.0571, the attainable two-sided floor at four
against three.

**This is reproduction, not a finding of this repository.** The direction was published in
2008 and was declared here before the value was read. What it buys is the thing A15's gate
condition 1 demands and the first freeze had dropped: evidence that the perturbation did
something in these particular mice.

Its limit is real and was recorded in advance: a whole-lung collagen programme at day 7 is a
downstream fibrotic outcome, not a direct measurement of integrin occupancy. It establishes
a pharmacodynamic effect, not target engagement in the strict sense.

## The epithelium's own programme did not separate

All contrasts are 3G9 against Axum8, both bleomycin, both batch S151, both one female and
three males, four mice each, two-sided by declaration because neither the A15 mechanism nor
rival 2 predicts a sign.

| Endpoint | Members | Mean difference | Exact two-sided p | Shift, 97.14 per cent interval |
|---|---|--:|--:|---|
| Transitional: `Krt8`, `Krt19`, `Cldn4`, `Cdkn1a`, `Krt7` | 5 of 5 measurable | -0.3599 | 0.8857 | -0.4181, [-1.9112, +1.4275] |
| Identity: `Sftpc`, `Ager` | 2 of 2 measurable | -0.2107 | 0.6857 | -0.4366, [-2.4426, +1.5693] |
| Omnibus centroid, 15,062 genes, `Itgb6` excluded | | statistic 79.667 | 0.6571 | |

The omnibus dispersion diagnostic, declared first because the statistic it guards confounds
location with dispersion, gives a within-arm distance ratio of **1.0661** against a declared
threshold of 1.5, so the omnibus is readable as a location statistic and is not being
rescued by a dispersion artefact.

Both declared sensitivities agree, so no downgrade applies: the unstandardised composites
give -0.1563 and -0.0497, and the median-of-ratios normalisation gives -0.1491 and -0.0425,
all non-separating. The demoted between-minus-within omnibus statistic gives -0.038365 at
p 0.6286.

For scale, the same transitional composite separates bleomycin from saline at **+1.5637**,
p 0.0571. The injury axis this panel indexes therefore moves in these data; what did not move
detectably is the 3G9 contrast. The two figures are standardised within their own contrasts
and are not on a common scale, so they are reported side by side and not divided.

## The immunoprecipitation was pure, and comparably so in both arms

This is the check the first freeze wrongly declared impossible. Every mouse has a paired
input from the same homogenate, so enrichment is computable per mouse.

| Marker group | Enrichment, 3G9 | Enrichment, Axum8 | Arm difference |
|---|--:|--:|--:|
| `Epcam` | +1.173 | +1.073 | +0.101 |
| `Cdh1` | +1.084 | +0.887 | +0.197 |
| `Nkx2-1` | +0.779 | +0.693 | +0.087 |
| `Ptprc` | -4.831 | -4.792 | -0.039 |
| `Pecam1` | -3.812 | -3.892 | +0.080 |
| `Col3a1` | -4.146 | -4.093 | -0.053 |
| `Lyz2` | +0.322 | +0.102 | +0.220 |
| `Cd68` | -3.692 | -3.398 | -0.294 |
| `Itgax` | -4.401 | -4.224 | -0.177 |
| `Mrc1` | -3.906 | -4.010 | +0.104 |

Epithelial markers are enriched by 0.78 to 1.17 log2 units and immune, endothelial and
mesenchymal markers are de-enriched by 3.8 to 4.8. **No arm difference exceeds 0.294 against
a declared downgrade threshold of 1.0**, and none was left uncomputed. So the epithelial
null is not a purity artefact, and a composition gradient between the arms is not what
produced it.

The A0 set, demoted from primary to a handling and stress covariate, gives -0.1632 at
p 0.4857 and does not separate, so a handling-time artefact is not visible either. The
deposit still records no harvest time, warm ischaemia time or homogenisation order, so that
covariate constrains the possibility rather than excluding it.

## What this establishes, and what it does not

**Establishes, descriptively, in one experiment of four mice against four.** That under a
dose of this antibody sufficient to separate the whole-lung collagen and myofibroblast
programme completely, no difference was detected in the epithelial compartment's own
identity or transitional programme, nor in a 15,062-gene omnibus profile, with the
immunoprecipitation pure and comparable between arms and no covariate separating.

**Does not establish.** Anything about A15. This deposit has no activated-TGF-beta readout,
no separated fibroblast compartment and no ligand arm, and the freeze forbids reading either
direction as evidence about the parent question.

**Does not establish absence.** Four against four reaches nominal significance only under
complete separation, and the shift intervals span **[-1.91, +1.43]** and **[-2.44, +1.57]**
standardised units. A null here is compatible with a true shift of roughly two to three
between-mouse standard deviations. The frozen prohibition against writing that the epithelium
did not move is deliberate and is observed here.

**Does not transfer as a fact.** A systemic antibody at day 7 in a mouse lung is not a
two-week organoid co-culture of mouse epithelium with human fibroblasts. The bound transfers
as a plausibility argument only.

## Two things worth recording against the next design

1. **The secondary is not corroboration, and the numbers confirm it.** The A5/A11 injury
   residual module gives -0.1412 at p 0.6857, and the Pearson correlation between its
   per-mouse score vector and the A0 covariate's is **0.8112**, above the 0.7 the freeze set
   as the point where corroboration language must be dropped. Two overlapping sets from one
   source study, scored on the same eight mice, agree because they are nearly the same
   instrument.
2. **The depth diagnostic is not clean.** The transitional composite correlates with log
   library total at +0.6165 and the identity composite at -0.7958 across the eight mice.
   Both were declared descriptive with no decision weight, so neither changes the reading,
   but at eight mice a correlation of that size is a caution for any future design on this
   deposit rather than a reassurance.

Two genes were declared in this freeze rather than inherited, and carry no decision: `Sfn` at
-0.342, p 0.3429, and `Hopx` at -0.0931, p 0.6857.

## Comparison against the source paper, in the declared order

The freeze required the source paper's own 3G9 results to be read **after** the computation,
not before, and that order was kept. Auyeung 2022
([doi 10.1152/ajplung.00408.2021](https://doi.org/10.1152/ajplung.00408.2021)) is the
deposit's source paper. Its abstract, which was on the record here before execution, reports
that the IRE1alpha inhibitor lowers integrin alphaVbeta6 expression with correspondingly
lower TGF-beta-induced epithelial gene expression. **The whole-lung collagen result here
reproduces the established direction for this antibody and is reported as reproduction.** No
claim of novelty attaches to any number in this report.

## Register position

No claim row is added, and the register still ends at C168. Grading is the owner's decision.
If a row is wanted, the sentence the artefacts support is:

**Descriptive only, in one experiment, bounding a rival and not a hypothesis.** In
GSE190821, 3G9 anti-integrin-beta6 lowers a frozen five-gene collagen and myofibroblast
programme in whole lung by 1.2424 standardised units with complete separation of four
treated from four control mice, while in the same mice the epithelial immunoprecipitation
shows no detectable difference in a frozen transitional panel (-0.3599, p 0.8857), a frozen
identity panel (-0.2107, p 0.6857) or a 15,062-gene omnibus centroid statistic (p 0.6571),
with epithelial enrichment differing between arms by at most 0.294 log2 units.
