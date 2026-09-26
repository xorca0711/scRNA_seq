# Why the first rival-2 freeze is withdrawn, before anything was scored

27 September 2026. The first freeze,
[`config/a15_rival2_freeze.json`](../config/a15_rival2_freeze.json), is **withdrawn and
preserved unedited** so its recorded bytes keep verifying.
[`config/a15_rival2_freeze_v2.json`](../config/a15_rival2_freeze_v2.json) is the authority
on what may now be computed. **No endpoint had been scored when it was withdrawn**, and the
execution script refused to run against it, which is the fail-closed check working rather
than a near miss.

It was withdrawn after a three-lens adversarial review of the draft, run before any value
was read. This is the second time in this repository that such a review has withdrawn a
freeze; the first was A2's rank test. Every factual claim the review made was independently
verified against the tracked files before any of it was accepted, and one recommendation was
rejected with its reason recorded below.

| Lens | Verdict | Findings | Fatal |
|---|---|--:|--:|
| Design and confounding | freeze should be withdrawn | 19 | 3 |
| Inference and statistics | freeze needs revision | 14 | 3 |
| Scope, logic and interpretation | freeze needs revision | 16 | 3 |

The three lenses converged independently on the same two defects, the endpoint's construct
validity and the missing engagement check, and each added one the others did not see. The
inference lens confirmed the audit's arithmetic, that the declared p-value floors of 0.0143
one-sided and 0.0286 two-sided at four against four, and 0.0286 one-sided at four against
three, are correct; what it rejected was the inference built on them.

## The fatal reasons

### 1. The primary endpoint could not measure what rival 2 asserts

The first freeze made the A0 frozen 50-gene transition programme the primary endpoint and
justified it as being "exactly the epithelial state change rival 2 asserts". That
justification is false, and the files say so.

Rival 2's cited support is Morris 2003, whose phenotype is Mmp12-dependent emphysema, and
Koth 2007, whose phenotype is bronchoalveolar phospholipid and collectin accumulation. The
A0 programme contains **no matrix metalloproteinase, no surfactant protein and no
collectin**. Two of its fifty members touch phospholipid handling at all, `Pla2g15` and
`Mfge8`.

This is not an accident of gene selection, it is structural. A0's own discovery
configuration excludes surfactant and identity genes from its universe by name:
`SFTPC`, `SFTPA1`, `SFTPA2`, `SFTPB`, `SFTPD`, `ABCA3`, `NAPSA`, `SLC34A2`, `AGER`, `AQP5`,
`HOPX`, `PDPN`, `CAV1`, `KRT8`, `KRT19`, `CLDN4`, `CDKN1A` and `LYZ` among others. **The A0
programme is incapable of carrying rival 2's axis by construction.** What remains in it is a
stress, chaperone, immediate-early, mitotic and matrix residual set.

The endpoint had been chosen for administrative fitness, that it was frozen in-house before
this deposit was opened, that it recovered 50 of 50 members, and that it did not contain the
perturbed gene, and the freeze then presented that administrative fitness as biological
fitness. That is the error.

### 2. The score is an abundance proxy, not a state measure

A0's fifty genes were selected to discriminate Krt8-positive alveolar differentiation
intermediate cells from activated AT2 cells, in single-cell data, in Strunz GSE141259, at
bleomycin **days 10 to 15**. Averaged over a bulk whole-epithelium immunoprecipitation at
**day 7**, a state-discriminating programme estimates how much of the compartment occupies
that state. That is injury severity and cell composition, which is A10's caution in A10's
own words, and it is the instrument's definition rather than a limitation of it.

This destroys the one direction the first freeze declared informative. A nonzero effect is
predicted a priori by the A15 mechanism itself, since less activated TGF-beta means less
transitional state and less severity, and 3G9 is published to attenuate bleomycin fibrosis.
A null means severity had not moved by day 7, which says nothing about whether the
epithelium was altered. A declared limitation cannot repair a design whose only informative
outcome that limitation invalidates. There is also a window mismatch, since the programme
was defined at days 10 to 15 and this deposit is day 7.

### 3. No perturbation-engagement check, while only the null was called informative

The first freeze declared that a small effect bounds rival 2, and established no check that
3G9 engaged integrin beta6 in these four mice. Its only pre-read gate was the
bleomycin-against-saline instrument check, which tests the endpoint's dynamic range and not
the perturbation. 3G9 blocks the protein, so the Itgb6 transcript cannot serve, and the
freeze attached no decision to it.

A15's own [PLAN.md](../PLAN.md) gate condition 1 requires "a recorded validation that the
perturbation took effect". **The side-branch dropped the parent gate's most protective
condition for the one run that actually executes.** The outcome the freeze would have called
a bound is exactly the outcome produced by an antibody that did not work, an insufficient
dose, a wrong window, or a translatome insensitive to the change, and those are not
distinguishable here.

### 4. Purity was disclaimed when the deposit makes it computable

The design lens found a factual error in the freeze, and it is the one that most improved the
replacement. The withdrawn freeze said the epithelial compartment's purity "is not verified in
the deposit". The deposit states that **an aliquot of each mouse's homogenate was taken as the
input**, which is the library labelled Whole Lung, and the remainder was immunoprecipitated.
Every mouse therefore has a paired input from the same homogenate, which is the standard
RiboTag enrichment control. Per-mouse purity, per-mouse injury severity, and a
perturbation-engagement control are all computable from it. The
[erratum](RIVAL2_STAGE1_ERRATUM.md) records that correction and the related one, that the
deposit does state genotype, bleomycin dose and the 3G9 schedule, which the stage 1 addendum
had listed as unstated.

### 5. The A0 set is sensitive to handling and to contamination, neither of which is recorded

Two further properties of the withdrawn primary endpoint, both verified. Its overlap with A1's
frozen epithelial marker panel for this same deposit is **exactly zero**. And a chaperone and
immediate-early score is the standard readout of warm ischaemia and handling time, while the
deposit records no harvest time, ischaemia time or homogenisation order for any of the 11 mice;
its members `Fn1`, `Myl9`, `Vcl` and `Mcam` additionally make it sensitive to mesenchymal
contamination of the immunoprecipitation. At four against four a single mouse handled
differently could dominate a rank statistic.

### 6. The ratio rule was quantifiably biased toward its own favourable verdict

The inference lens showed analytically, and confirmed by simulation, that the primary effect
and the injury yardstick share the Axum8 bleomycin arm, giving a correlation of **-0.463** by
construction. It further showed that gating the primary on an instrument check that only
complete separation can pass selects the denominator upward, so the compound rule declared
rival 2 "bounded" a substantial fraction of the time when the true effect sat exactly on the
freeze's own "rival stays live" threshold. A decision rule that reaches its favourable verdict
by construction is worse than no rule.

## The serious reasons that also forced changes

| Finding | Verified against the files | Change in v2 |
|---|---|---|
| The gate accounting understated the shortfall as conditions 3, 4 and 5 | Condition 1 requires recorded validation of the perturbation; nothing establishes it | v2 states conditions **1**, 3, 4 and 5 |
| The two endpoints were called independent corroboration | **11 of A0's 50 genes are members of the 386-gene injury residual module**, 22.0 per cent: `2200002D01Rik`, `Atf3`, `Dusp1`, `Edn1`, `Hspb1`, `Mfge8`, `Nop58`, `Plk2`, `Tnfrsf12a`, `Tnip3`, `Ubc`. Both trace to Strunz 2020: A0's discovery set is GSE141259 and the injury residual module is `ADI_published_400` from doi 10.1038/s41467-020-17358-3 | v2 records the overlap and the shared provenance, and relabels the secondary a **correlated sensitivity**, not corroboration |
| The ratio-to-injury margin benchmarks against the largest perturbation in the deposit | A bounded criterion of 0.25 times an injury effect of, say, 1.5 log2 units would pass 0.375 log2 units as small, which is of the same order as the -0.938 log2 effect this repository treats as A2's founding observation | v2 drops the ratio as a decision rule, keeps the injury contrast as context with no decision, and benchmarks effects in absolute log2 CPM units against A2's -0.938 explicitly |
| Decisions rested on point estimates with no interval | Numerator from 4 against 4, denominator from 4 against 3, no interval declared for either | v2 carries the decision on an exact permutation test and reports distribution-free intervals |

## One correction to the review, recorded because it changes the remedy

The reviewing lens recommended replacing the endpoint with a rival-2-specific set drawn from
Morris 2003 and Koth 2007. That recommendation is **not adopted**, for a reason the review
did not weigh: **both of those papers attribute their phenotypes to loss of TGF-beta
activation**, in their own titles. Morris 2003 is "Loss of integrin alpha(v)beta6-mediated
TGF-beta activation causes Mmp12-dependent emphysema" and Koth 2007 is "Integrin beta6
mediates phospholipid and collectin homeostasis by activation of latent TGF-beta1". A
metalloproteinase and surfactant panel would therefore measure **TGF-beta-mediated**
consequences of removing the integrin, not a TGF-beta-independent route, so it would not
separate rival 2 from the mechanism A15 proposes either.

That correction forces rival 2 to be restated, and v2 restates it: the rival is not that the
epithelium changes by some non-TGF-beta route, but that **the epithelial integrin's effect on
the fibroblast may be indirect, mediated by a change in the epithelium rather than by
TGF-beta activated at the epithelial surface acting on the fibroblast**. Under that
statement the question is whether the epithelium changes **at all**, which is an omnibus
question and needs no curated gene list. That is what v2 declares, and it is why v2's primary
is endpoint-free.

A second recommendation, an omnibus divergence statistic, **is** adopted, and becomes the
primary.

## What v2 keeps from v1

The design audit stands unchanged and is not re-run: the arms, the batch and sex matching,
the join and its independent agreement with A1's records, the naming hazards, the depth
audit, the exact power floors and the precedent finding. Stage 1 passed all ten of its stop
rules and nothing in this withdrawal touches it. What changed is the endpoint, the statistic,
the decision rule and the reading, all before any value was read.
