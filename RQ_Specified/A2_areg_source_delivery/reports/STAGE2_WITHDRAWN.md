# A2 stage 2: the first freeze is withdrawn

26 September 2026. The first freeze,
[`config/a2_stage2_freeze.json`](../config/a2_stage2_freeze.json), and its report
[STAGE2_FREEZE.md](STAGE2_FREEZE.md) are **preserved unchanged and withdrawn**. The
replacement is [`config/a2_stage2_freeze_v2.json`](../config/a2_stage2_freeze_v2.json).
No endpoint was scored under either freeze, so nothing computed is being revised;
what is being withdrawn is an inference the design cannot support.

## How this was found

Three independent review lenses were run over the stage 0 to 2 documents: one on
biological and literature accuracy, one on statistical validity and
pre-registration integrity, one on internal consistency. They returned 17, 21 and
17 findings. Six were blocking and the three lenses converged on two of them. In
response I ran a fibroblast-side covariate pass,
[`scripts/01b_stage1_addendum.py`](../scripts/01b_stage1_addendum.py), which
supplied the decisive number. Every claim below was then verified directly against
the deposit or this repository's own tables before being accepted.

## The five reasons

**1. The recipient makes the ligand.** This is the one that matters.

| Compartment | AREG, mean log2 CPM | Wells with any count |
|---|--:|--:|
| Human fibroblasts, plate 3 | 9.846 | 99.2% |
| Mouse epithelium, not targeted | 7.267 | recorded by A10 |

Each figure uses its own compartment's total as the denominator, so this is a
within-compartment comparison rather than a count of molecules. Even so, the
unedited human fibroblasts transcribe AREG at a higher within-compartment level
than the mouse epithelium whose Areg the knockout removes. **Removing the
epithelial source does not remove AREG from the culture.** No necessity claim and
no sufficiency claim is available from this design, and a null is uninterpretable.
The most leg 1 can now do is bound the incremental contribution of the epithelial
source on top of an unremoved autocrine source.

The repository's own chosen citation predicted this. Zhou et al. 2012, cited in the
rationale for the fibroblast arm of the mechanism, is the paper showing that
TGF-beta1 induces amphiregulin in lung fibroblasts and that the silencing which
blocked the phenotype was directed at fibroblast AREG
([10.1074/jbc.M112.356824](https://doi.org/10.1074/jbc.M112.356824)). The rival was
in the citation and not in the rival list.

**2. There is no randomization.** Fifty of the 53 plate-3 targets occupy one fixed
well position in all four units, with Areg always at F07. Target is therefore
confounded with plate position and with guide pool, and the four units are four
copies of one layout rather than four draws. A uniform-rank null has nothing behind
it, so the exact p-value, the critical rank sum of 64 and the alpha of 0.05 are all
withdrawn.

**3. The reference set was not neutral.** Plate 3 carries AREG, EGFR, ERBB2, ERBB3,
ERBB4 and the MAP kinase cascade among its 53 targets. Ranking Areg among its
plate-mates therefore asks whether it is extreme among other perturbations of the
same pathway, which is not the estimand the freeze declared, and the claim that
same-pathway neighbours make the test conservative was asserted rather than derived.

**4. The denominator held ineligible wells.** Plate-3 fibroblast totals start at one
count. The first freeze applied its 100,000-count floor only to sensitivities, so
wells whose scores are deterministic sat in the primary rank denominator.

**5. The consistency requirement was vacuous.** Its exact size equalled the
rank-sum size, 0.04903, so it added no stringency while reading as though it did.

## What replaces it

An effect size with a direction count, and no p-value. That is what the register
card promised from the start, and the design supports it.

- **Primary endpoint promoted:** the five-gene fibroblast activation score, because
  alpha-smooth-muscle actin and collagen are what both cited mechanism papers
  measured. The Hallmark TGF-beta set is demoted to a pathway-level secondary: 16 of
  its 54 members are negative regulators or feedback genes and none of the five
  activation genes is in it, so its direction is mixed and it was never the
  mechanism-matched endpoint.
- **The screen's own control is used.** TIGIT is the in-plate control the article
  names, at six wells per unit, which A10's source design check already recorded.
  With TDTOMATO that gives eight control wells per unit instead of two.
- **The comparison** is the Areg well's endpoint minus the mean of its unit's
  eligible control wells, reported per unit with the direction count, both unmatched
  and restricted to controls within a factor of four of the Areg well's depth.
  Eligible controls per unit are 3, 6, 5 and 8; within the depth band, 3, 5, 5 and 7.
- **Epithelial fraction is demoted** to a two-sided sensitivity, because epithelial
  abundance plausibly lies on the causal path from the knockout to the fibroblast
  read, and about a fifth of reads are unassigned to either species. An effect that
  depends on that adjustment either way is reported as mediation-ambiguous.
- **The receptor arm is restricted to EGFR**, the only AREG receptor among the
  targets, with ERBB2 as its heterodimer partner. ERBB3 binds neuregulins and ERBB4
  binds neuregulins, HB-EGF, betacellulin and epiregulin, so neither is a test of
  AREG reception; both are reclassified as non-AREG-receptor perturbation controls.
- **Leg 2 is reclassified exploratory.** The logged trial E6 table already records
  the activation score against fibroblast depth at rho 0.4116 and fibroblast EGFR at
  0.4918, and both members of the new pair are fibroblast detection fractions, so
  the inherited depth gate is likely to refuse it. A co-regulation control is now
  required alongside the detection control, because TGFB1, THBS1, ITGB8 and LTBP1 are
  themselves TGF-beta inducible.

## Corrections to the stage 1 report

- It says stage 1 passed "all seven stop rules". The plan declares four stop
  conditions; seven is the number of audit items. Item 7 returned a carried
  assumption, not a pass.
- It names three sub-threshold contrast wells. At the 100,000-count floor the freeze
  then adopted, six axis wells are below it: Erbb2 in all of 3-1, 3-2 and 3-3, Erbb4
  in 3-1, Egfr in 3-2 and Itgb6 in 3-3.
- It implies Egfr is the weakest contrast. By readability Erbb2 is weaker: Erbb2
  clears the floor in one unit of four, Egfr in three, Erbb3 in all four.

## What this costs, stated plainly

A2's decisive experiment is no longer in this screen. Source necessity needs a
system where the recipient does not supply the ligand, or a fibroblast-side
perturbation, or measured extracellular ligand. The delivery-versus-abundance
contrast that makes A2 distinct was never testable in a single well with one source
compartment and no spatial variation; it belongs to the spatial layer, which the
rationale already lists as feasible pending a cohort audit. What remains here is a
bounded, descriptive, within-screen effect size, and it is worth reporting as that
and nothing more.

The hypothesis itself is untouched. The mechanism still names the recipient, the
recipient machinery is expressed in these fibroblasts (integrin alphaV at 5.835 log2
CPM in 93.8 per cent of wells, ITGB1 at 9.487, ITGB8 at 5.254, LTBP1 at 9.757), and
fibroblast EGFR is present at 5.389 in 90.8 per cent of wells. Those are the first
measurements in this repository showing the receiver is equipped, which bears on
C36 without settling it.

## Authorization

The plan placed an owner retain step between stages 1 and 2. The owner authorized
both stages in one instruction, which necessarily predates stage 1's findings, so
**that gate was not exercised against what stage 1 found.** The second freeze is
provisional pending owner review, and stage 3 is not authorized.
