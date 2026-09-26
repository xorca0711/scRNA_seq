# A15 rival-2 stage 1: can GSE190821 carry the bounding comparison

27 September 2026. The side-branch authorized by the owner bounds
[rival 2](../RATIONALE.md) only: whether blocking epithelial integrin beta6 moves the
epithelium's own programme. **It is not a test of A15 and cannot become one.** Nothing is
scored in this report. Every number comes from
[`tables/rival2/`](../tables/rival2/) and the
[run record](../tables/rival2/stage1_design_run.json).

## Verdict

**The design can carry the comparison, and it is cleaner than expected.** All ten declared
stop rules pass. The antibody experiment is self-contained: the treated arm, its antibody
control and a saline control all sit in one batch, sex-matched, with one epithelial library
per mouse.

**Three findings change what stage 2 may declare**, and two of them were not anticipated
when the side-branch was proposed.

## The arms, recomputed from the deposit

| Arm | Mice | Identifiers | Batch | Sex |
|---|--:|---|---|---|
| Bleomycin, 3G9 anti-integrin-beta6 | 4 | 218, 219, 226R, 228 | S151 | 1 F, 3 M |
| Bleomycin, Axum8 inert antibody | 4 | 204, 206, 208R, 225 | S151 | 1 F, 3 M |
| Saline, Axum8 inert antibody | 3 | 207, 221, 230 | S151 | not used as a primary arm |
| Bleomycin, Vehicle | 5 | 17, 31, 45, 148, 149 | S061, S135 | A1 used this arm |
| Bleomycin, KIRA8 | 5 | 20, 32, 41, 154, 155 | S061, S135 | A1 used this arm |
| Saline, Vehicle | 3 | 146, 152, 153 | S135 | not used |

The primary arms **share one batch and are sex-matched**, so treatment is separable from
both. A saline arm on the same antibody sits in the same batch, so the injury yardstick
does not have to cross batches. Every one of the 24 mice carries exactly one epithelial and
one whole-lung library, so no mouse contributes two libraries to one arm.

Depth across the eleven eligible epithelial libraries runs from 24,751,240 to 36,057,881
counts, a ratio of 1.457, and none falls below the 100,000-count floor A1 recorded for this
deposit. This is the opposite of the four-orders-of-magnitude depth span that constrained
A2, so depth is not the binding problem here.

## Finding 1: the column names cannot identify the antibody control

The deposited counts file uses two naming conventions and three token variants. The join
had to be built by rule and verified, not assumed:

- twelve columns carry zero-padded mouse identifiers (`Bleo_017_E` for mouse 17);
- one carries an `FT` suffix (`Bleo_3G9_I_218FT` for mouse 218, whole lung);
- some columns place the compartment token before the mouse identifier and some after.

The join resolves one-to-one onto all 48 libraries with no duplicates and no unmatched
rows in either direction, and **the rule independently reproduces all ten of A1's recorded
assignments** for this deposit, which A1 built for a different contrast.

The hazard is worth stating on its own. The prefix `Bleo` covers **both**
Bleomycin/Axum8 and Bleomycin/Vehicle, so a session that grouped libraries by column-name
prefix would silently pool the inert-antibody control with the KIRA8 vehicle control and
call the result a bleomycin control arm. **Only the `treatment` field in the series
metadata separates them.** A1 recorded the reason this matters in its own words: Axum8 is
the antibody control and is not a vehicle control for KIRA8.

## Finding 2: A1 has already used this deposit, and left this contrast untouched

The precedent check fired, correctly. A1 froze an analysis of **the same deposit, the same
compartment and the same exposure** on 25 September 2026
([`config/ire1_kira8.json`](../../A1_transitional_epithelial_state_distinction/config/ire1_kira8.json)),
and the counts file this audit downloaded is **byte-identical to A1's recorded hash**,
`2c8a21fd…a0b595d1f7`. That is an input verification, not a coincidence.

What A1 scored was Vehicle against KIRA8 in the epithelium, ten libraries. It explicitly
excluded the Axum8 arm. **No prior analysis in this repository has scored 3G9 against
Axum8**, so the contrast is untouched.

A1's instrument cannot be reused. It is edgeR and limma in R, with TMM normalisation and a
quasi-likelihood fit, and **no R interpreter is available in this session**. Stage 2
therefore declares a score-level statistic of its own and must not present it as A1's
instrument. Reimplementing TMM and calling it edgeR output is refused.

## Finding 3: the source paper is Auyeung 2022, and it may already report this contrast

The deposit names PubMed 35170357: Auyeung et al., *IRE1alpha drives lung epithelial
progenitor dysfunction to establish a niche for pulmonary fibrosis*, with Dean Sheppard
among the authors ([doi 10.1152/ajplung.00408.2021](https://doi.org/10.1152/ajplung.00408.2021),
retrieved from PubMed). A1 already engaged with this paper, citing its Fig. 1K and its
Supplemental Table 1, so this is not a paper under the unread-source prohibition that
applies to roadmap paper 14.

Two consequences, and they must be carried into the freeze rather than discovered later.

1. **The 3G9 arm exists because the authors used it.** The abstract states that the IRE1alpha
   inhibitor decreases expression of integrin alphaVbeta6 and that this corresponds to
   decreased TGF-beta-induced gene expression in the epithelium and decreased collagen
   accumulation. A published 3G9 epithelial contrast is therefore likely. **Agreement with
   it would be reproduction and must be reported as reproduction, never as a finding of
   this repository.**
2. **Prior exposure is recorded, not denied.** This audit read that abstract. It states a
   direction for an epithelial TGF-beta response module. The primary endpoint declared in
   stage 2 is **not** a TGF-beta response module, so the exposure does not reveal the
   primary's direction, but it is logged here because A0's precedent is that prior exposure
   is recorded rather than assumed away.

The paper's own reported 3G9 results are deliberately **not** read before the freeze is
committed. They are compared against afterwards, in that order.

## Endpoint coverage, against the deposited gene index

The counts file carries 42,548 Ensembl mouse stable gene IDs, not symbols, so every
declared gene set was mapped through the Ensembl REST symbol endpoint, the same one A1
used, with current stable IDs only, no genomic coordinates transported to GRCm38 and **no
alias repair**. Symbols that do not resolve exactly are reported, not silently fixed.

| Gene set | Declared | Present in counts | Fraction | Unresolved symbols | Contains Itgb6 |
|---|--:|--:|--:|---|---|
| A0 frozen transition programme | 50 | 50 | 1.000 | none | no |
| A5/A11 injury residual module | 386 | 373 | 0.966 | `1500012F01Rik`, `Atp5d`, `Cyr61`, `Gnb2l1` and nine more | no |
| A5/A11 shared remodelling module | 12 | 11 | 0.917 | `Nars` | **yes** |
| A1 frozen marker panel | 8 | 8 | 1.000 | none | **yes** |

Two sets contain the perturbed gene and are therefore **barred from use as the primary
endpoint**. The shared remodelling module is excluded outright. A1's marker panel is
retained only for per-gene reporting, where Itgb6 is informative in its own right: 3G9
blocks the protein, so whether its transcript moves is a separate question from whether the
programme moves.

The A0 programme recovers all 50 members and does not contain the perturbed gene, so it is
eligible as the primary.

## Power, computed exactly before any value is read

| Contrast | Sizes | Assignments | Smallest attainable one-sided p |
|---|---|--:|--:|
| 3G9 against Axum8, bleomycin | 4 against 4 | 70 | 0.0143 |
| Axum8 bleomycin against Axum8 saline | 4 against 3 | 35 | 0.0286 |

Both floors are below 0.05, so an exact rank-sum test can reach nominal significance, **but
only under complete separation of the arms**. There is no configuration short of perfect
separation that reaches it. **Precise absence is unavailable at these sizes**, and that is
declared here rather than discovered after a null.

## What stage 2 must therefore fix

1. A primary endpoint that does not contain the perturbed gene, which means the A0 frozen
   programme, with the injury residual module as a second, larger set.
2. A score-level statistic declared in this repository's own terms, not presented as A1's
   edgeR instrument, with the absence of R recorded.
3. An injury yardstick inside the same batch and antibody, so that the size of the 3G9
   effect can be read against the size of the effect injury itself produces.
4. A reading that is **asymmetric on purpose**: a small effect bounds rival 2, while a large
   effect does not establish it, because an epithelial change could be downstream of the
   mechanism this question proposes rather than independent of it.
5. An explicit declaration that a published 3G9 contrast may exist, that agreement is
   reproduction, and that the endpoint may not change after the paper is read.
6. The cross-system limits, stated before any value is read: a systemic antibody at day 7
   in a mouse lung is not a two-week organoid co-culture with human fibroblasts, and a
   RiboTag immunoprecipitation measures ribosome-associated RNA in a compartment whose
   subtype composition can itself shift under injury.
