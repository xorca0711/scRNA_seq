# A2 stage 5: what the question now stands on

26 September 2026. Closing report for the six declared stages. It draws only on the
logged artefacts: the [stage 1 audit](STAGE1_AUDIT.md), the
[withdrawal](STAGE2_WITHDRAWN.md) of the first freeze, the
[leg 1 results](STAGE3_LEG1_RESULTS.md) and the [leg 2 results](STAGE4_LEG2_RESULTS.md).

## The short version

The delivery hypothesis is neither supported nor refuted, and the two experiments this
repository could run have told us something else instead. The ligand arm did not move.
The epithelial integrin arm moved a lot. The donor-level leg was refused by its own
control. What A2 gained is a sharper question and one measured lead; what it lost is the
belief that this screen could settle source necessity.

## Leg 1, in the organoid screen

The five-gene fibroblast activation score, against the screen's own depth-matched
controls, across four plate-replicate units:

| Epithelial knockout | Median effect | Units in the predicted direction | Culture unchanged |
|---|--:|---|---|
| Areg | +0.036 | 2 of 4 | yes |
| **Itgb6** | **-0.938** | 3 of 3 readable | yes |
| Erbb3 | -0.187 | 4 of 4 | no |
| Egfr | +0.247 | 0 of 3 | yes |
| Erbb4 | +0.503 | 0 of 3 | yes |
| Erbb2 | -0.169 | 1 of 1 readable | no |

The endpoint's standard deviation over the 168 eligible wells is 0.482, so Areg is about
0.07 standard deviations in the wrong direction and Itgb6 about 1.9 in the right one. The
Hallmark secondary agrees on Areg. There is no p-value in this leg by declaration,
because the four units are copies of one plate layout with each target at a fixed well.

**The Areg null is bounded, not absent.** The recipient transcribes AREG at 9.846 mean
log2 CPM in 99.2 per cent of wells against the 7.267 the knockout removes; the
perturbation halves one transcript; the culture medium is not in the deposit, so
exogenous EGF or a TGF-beta receptor inhibitor cannot be excluded. The honest sentence is
that no epithelial AREG contribution is detectable on top of an unremoved autocrine
source, in this assay, with this partial perturbation.

**The Itgb6 effect survives the checks available.** It holds in every unit where the well
clears the fibroblast eligibility floor, it survives the epithelial-fraction adjustment
at -0.992, and the post hoc culture check shows organoid size within 0.022 log units of
its controls and fibroblast content at 0.93 of theirs. Erbb3, the other consistent
decrease, fails that same check with organoids 0.605 log units smaller and 2.31 times the
fibroblast material, so the generic-perturbation caution weighs on Erbb3 rather than on
Itgb6.

## Leg 2, in the fibrosis cohort

Refused by the inherited C51 depth rule. The machinery composite reaches rho 0.770
against fibroblast depth and the activation outcome 0.412, so both members exceed the
0.4 threshold and the nominal +0.433 between them is not read. The IPF stratum is also
not read, a judgement stricter than the declared rule, because the predictor still sits
at 0.717 there.

Two by-products are worth keeping. The instrument reproduced C50 exactly, at -0.150,
so the refusal is not a broken pipeline. And the co-regulation rival was addressed and
was not the problem: the machinery composite exceeded a detection-matched
oxidative-phosphorylation control (+0.208) and a housekeeping control (+0.066). Depth is
what defeats this leg, and the repository already owns the fix it needs, which is the
common-molecule-budget treatment it used for C37.

## What is established, and what is not

**Established, descriptively and within one screen.** The frozen five-gene fibroblast
activation score does not move when epithelial Areg is knocked out in this co-culture,
against the screen's own controls, in four units of one plate layout.

**Established as a decision record.** That the C51 depth rule, written down before the
pair was computed, refused a correlation that would have been nominally significant at
p 0.044. That is the rule earning its place a second time.

**Established as a measurement.** That the recipient in this co-culture carries the
machinery the cited mechanism needs: integrin alphaV at 5.835, ITGB1 at 9.487, ITGB8 at
5.254, LTBP1 at 9.757 and EGFR at 5.389 mean log2 CPM, each in nine tenths of wells or
more. That bears on C36 without settling it, because transcript is not protein.

**Not established.** That epithelial AREG has no role. That recipient machinery is
unrelated to fibroblast activation. That delivery rather than abundance governs the
response, which neither leg could test. Nothing here re-grades C37, C49 or C50.

**Proposed, as a hypothesis.** That the epithelial input to this fibroblast programme
runs through epithelial integrin-mediated activation of latent TGF-beta rather than
through amphiregulin. One well per target per unit at a fixed position cannot establish
it, and it needs its own design.

## Claim wording proposed for the register, not graded here

Grading is the owner's decision and no row is added by this analysis. If rows are wanted,
these are the sentences the artefacts support:

1. **Descriptive only.** In the GSE307112 alveolosphere screen, the frozen five-gene
   fibroblast activation score in epithelial Areg-knockout wells differs from
   depth-matched control wells by a median of +0.036 log2 CPM across four plate-replicate
   units, two of four in the predicted direction, against an endpoint standard deviation
   of 0.482. The recipient fibroblasts transcribe AREG at 9.846 mean log2 CPM, so this
   does not bound an epithelial contribution.
2. **Descriptive only, and the strongest observation in the panel.** In the same screen,
   epithelial Itgb6-knockout wells differ from depth-matched controls by a median of
   -0.938 log2 CPM on that score, in three of three readable units, with organoid size
   and fibroblast content unchanged and the effect surviving epithelial-fraction
   adjustment at -0.992.
3. **Not established, refused by a frozen control.** In GSE136831 a six-gene fibroblast
   integrin and latent-complex composite correlates with the fibroblast activation score
   at rho 0.433 across 22 donors, but the composite tracks fibroblast depth at 0.770 and
   the outcome at 0.412, so the C51 rule refuses the pair.
4. **Validated as a decision record.** The C51 depth rule refused a nominally
   significant donor-level correlation for the second time, on a pair declared before it
   was computed.

## The negative-results path, and why it is not taken yet

The plan directed this outcome pair, an inconclusive leg 1 and a leg 2 refused by its own
control, to NEGATIVE_RESULTS.md through the usual generation path. That file carries a
standing instruction not to edit it by hand: it is generated from CLAIMS.md by
`analysis/scripts/14_write_negative_results.py`, and the wording it prints is the
register's own.

So the entry cannot be written from here. Grading is the owner's decision and this
analysis adds no register row, which means there is nothing for the generator to pick up.
Once rows 1 and 3 of the wording proposed above are graded and entered, re-running that
script places them on the negative-results page without further work. Naming the gap here
rather than leaving it silent, because a result that exists only in a stage report is
easier to lose than one the generator owns.

## Other genomic layers, with feasibility verdicts

1. **A depth-standardized donor measure in the cohort already on disk.** The direct fix
   for leg 2, using the common-molecule-budget treatment this repository applied to C37,
   declared before any correlation. Verdict: **feasible now**, and the cheapest open item.
2. **Spatial transcriptomics of human fibrotic or tumour lung.** Still the only layer
   where proximity varies, so still the only one that can separate delivery from
   abundance. Verdict: **feasible, pending a cohort audit**, with a frozen proximity
   definition and the donor as the unit.
3. **An epithelial integrin perturbation with a TGF-beta activation readout.** What the
   Itgb6 lead needs: more than one well per target, positions that vary, and a measure of
   activated TGF-beta rather than a transcript proxy. Verdict: **not feasible from public
   data**; it is a bench design.
4. **A fibroblast-side ligand perturbation.** The only way to bound the autocrine source
   that leg 1 leaves in place. Verdict: **not feasible from public data**.
5. **Fibroblast surface protein by CITE-seq or equivalent.** Would convert the machinery
   and receptor transcript measurements above into the protein evidence C36 asks for.
   Verdict: **unverified feasibility**; needs a lung cohort with both modalities and
   adequate fibroblast depth.
6. **Chromatin accessibility at fibroblast TGF-beta response elements.** Tests licensing
   at the regulatory layer. Verdict: **unverified feasibility**, conditional on a public
   lung multiome with fibroblast depth.
7. **Sequence comparison of the mouse and human amphiregulin EGF-like domains.** Bounds
   the cross-species assumption leg 1 carried throughout. Verdict: **feasible and small**,
   from UniProt P31955 and P15514, and still open.

## Register readiness

A2 stays conditional and descriptive. Its decisive experiment is not in this repository's
current data, and the question now has one measured lead that points at the epithelial
integrin rather than the ligand. Whether that lead becomes its own question is the
owner's decision.
