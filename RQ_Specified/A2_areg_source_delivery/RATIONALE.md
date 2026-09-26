# A2 rationale: why AREG source matters through delivery, not through abundance

Originally declared 26 September 2026; **revised the same day after a three-lens
review and the stage 1 addendum it prompted.** This is the current biological
argument, not a retrospective pre-registration. The first stage 2 freeze is
withdrawn and preserved; read [STAGE2_WITHDRAWN.md](reports/STAGE2_WITHDRAWN.md)
for why, and [`config/a2_stage2_freeze_v2.json`](config/a2_stage2_freeze_v2.json)
for what may now be computed. The [contract](config/a2_delivery_contract.json) stays
unedited so its recorded hash keeps verifying, and [PLAN.md](PLAN.md) keeps its
original staging.

**What the revision changed.** The recipient-derived ligand rival was missing and is
now central. Leg 1 no longer claims to test delivery against abundance. The
recipient's own EGFR is restored to the mechanism. The cited recipient is named as a
pericyte. Two citations are corrected for hedging and scope, and the C37 row now
carries the register's own wording.

## The question the register used to ask, and why it is closed

A2 previously asked which compartment is the dominant AREG source, and whether that
dominance changes with tissue context. The repository has answered the abundance
version seven times, across lung adenocarcinoma and pulmonary fibrosis cohorts, and
no depth-independent epithelial hierarchy survived. Each row states the register's
own proposition, not a conclusion drawn from it:

| Register proposition | Claim | Status and what the numbers were |
|---|---|---|
| The epithelial versus myeloid AREG detection contrast in tumour lung is sensitive to annotation and molecule-depth standardization | C37 | Descriptive only. Medians 0.4591 against 0.2715 on deposited labels, all ten donors epithelial-higher, then six of ten at the primary 1,000-UMI budget |
| Epithelium is the AREG source in human lung adenocarcinoma | C39 | Not established. CD1c-positive dendritic cells reach 0.618 against tumour states at 0.555 and 0.532 |
| Dendritic cells and monocytes are a major AREG and HBEGF source, at or above the epithelial states | C45 | Descriptive only, in two diseases and three deposits. This is a ranking in the opposite direction to an epithelial hierarchy |
| AREG or HBEGF is tumour-enriched against matched normal lung from the same donors | C40 | Refuted at compartment level |
| AREG or HBEGF is enriched in fibrotic lung against control lung | C48 | Not established |
| Epithelial AREG and fibroblast EGFR are correlated across donors | C49 | Not established. Spearman rho 0.348, p 0.112, in 22 donors; it excludes nothing weaker than rho 0.43 |
| Epithelial AREG and a fibroblast activation score are correlated across donors | C50 | Not established. Spearman rho -0.150, with the sign reversed |

C51 is the decisive methodological entry. The only significant correlation in that
donor-level trial was a control pair whose two members both tracked sequencing
depth. A frozen rule refused it. The same number on AREG would have been read as the
discovery the source paper predicted.

Two readings are possible. The first is that the axis is not real. The second is
that total ligand abundance is not the quantity the biology rate-limits, so the
abundance comparisons measured the wrong thing. The mechanism literature supports
the second, and it is testable, though not by the experiment this repository first
reached for.

## The mechanism, including the part the first draft dropped

Amphiregulin is synthesized as a membrane-anchored precursor, and its recorded action
on mesenchymal cells is to convert a store the recipient already holds.

- **The recipient's own EGFR comes first, then integrin alphaV.** In two acute injury
  models, amphiregulin released bioactive TGF-beta from latent complexes by activating
  integrin alphaV on mesenchymal stromal cells, driving their differentiation into
  myofibroblasts, with macrophages identified as a critical source
  ([Minutti et al. 2019, Immunity](https://doi.org/10.1016/j.immuni.2019.01.008)).
  The integrin step is downstream of the recipient's EGFR: the effect is abolished by
  EGFR inhibition and in pericytes lacking EGFR. The recipient in that work is a
  PDGFRB-positive pericyte, so extending it to a lung fibroblast is an assumption
  carried here, on the same footing as the cross-species assumption below.
- **The fibroblast arm needs amphiregulin, including the fibroblast's own.** Silencing
  amphiregulin or inhibiting EGFR reduced TGF-beta1-driven fibroblast proliferation,
  smooth-muscle actin and collagen, and reduced collagen accumulation in TGF-beta1
  transgenic lungs ([Zhou et al. 2012, JBC](https://doi.org/10.1074/jbc.M112.356824)).
  The amphiregulin silenced there is the fibroblast's own, and TGF-beta1 induces it in
  lung fibroblasts. That makes the recipient a source in its own right, which is the
  rival this rationale first missed and which the stage 1 addendum has now measured.
- **Leukocyte amphiregulin can be non-redundant.** Removing amphiregulin from
  regulatory T cells alone caused severe lung damage during influenza with no change
  in suppressor function, antiviral response or viral load
  ([Arpaia et al. 2015, Cell](https://doi.org/10.1016/j.cell.2015.08.021)). That
  establishes non-redundance for one leukocyte lineage. No cited source addresses
  dendritic-cell or monocyte amphiregulin function, so C45's functional standing is
  open.
- **In a Hermansky-Pudlak syndrome model**, type 2 innate lymphoid cells may
  stimulate primary lung fibroblast proliferation and differentiation partly through
  amphiregulin and EGFR, hedged in those words by the authors
  ([Sorkhdini et al. 2024, JCI Insight](https://doi.org/10.1172/jci.insight.178381)).
  It is nonetheless the closest cited support for a lung fibroblast as the recipient.
- **Both epithelial and leukocyte sources are established**, in a review
  ([Zaiss et al. 2015, Immunity](https://doi.org/10.1016/j.immuni.2015.01.020)), which
  is why C45 constrains an epithelium-only reading rather than refuting the axis.

A short-range ligand that converts a recipient-held store does not require
tissue-level ligand to be rate-limiting, so C49 and C50 do not refute the axis.
Neither row is evidence for this framing either. C49's point estimate is positive, in
the direction the abundance version predicted, and underpowered.

## The hypothesis, stated so it can fail

**AREG's contribution to a fibroblast response is set by where the ligand is released
relative to a competent recipient, rather than by how much of it a compartment
transcribes.**

Two legs follow. Neither tests the distinguishing clause, and this rationale now says
so in the open rather than in a limitations paragraph.

**Leg 1, epithelial source contribution in contact.** Removing the epithelial source
changes a frozen fibroblast activation programme, while removing the epithelium's own
receptor does not. This is a source-contribution contrast, not a delivery test: a
single well holds one source compartment and no spatial variation, so a positive
result is equally consistent with the abundance version.

**Leg 2, the post-receptor recipient layer.** Across donors, the fibroblast activation
programme tracks the recipient's integrin and latent-complex genes more closely than
it tracks epithelial ligand abundance. The receptor layer the cited mechanism makes
obligatory is not measured here, because A9 owns fibroblast EGFR, so the estimand is
the layer below the receptor rather than recipient licensing as a whole.

**The delivery-versus-abundance contrast needs the spatial layer**, where proximity
varies. That is listed under other layers below, and it is the decisive experiment for
A2 as stated.

## Rivals, with the two that were missing

1. **The recipient supplies the same ligand.** Measured, and large: see below.
2. **Another EGFR ligand carries the response.** The register ranks AREG first among
   the retained exact-EGFR ligands (C32) but records that Ereg outranks Hbegf by
   enrichment over AT2 in three of four libraries (C33), and the screen carries no
   Ereg, Tgfa or Btc well. Redundancy cannot be excluded here.
3. **The epithelial state, not its ligand**, changes the fibroblast.
4. **The epithelium activates TGF-beta itself**, through its own integrin.
5. **Well composition or read depth** moves the fibroblast profile.

## What the audit found, and what it costs

The stage 1 addendum measured the fibroblast side directly. Each figure uses that
compartment's own total as its denominator.

| Gene, human fibroblast side | Mean log2 CPM | Wells with any count |
|---|--:|--:|
| AREG | 9.846 | 99.2% |
| EREG | 7.659 | 96.3% |
| EGFR | 5.389 | 90.8% |
| ITGAV | 5.835 | 93.8% |
| ITGB1 | 9.487 | 99.2% |
| ITGB8 | 5.254 | 88.3% |
| LTBP1 | 9.757 | 99.6% |

Two things follow, and they point in opposite directions.

**Against leg 1 as first specified:** the unedited fibroblasts transcribe AREG at a
higher within-compartment level than the mouse epithelium whose Areg is removed
(7.267 log2 CPM when not targeted). Removing the epithelial source does not remove the
ligand from the culture, so no necessity or sufficiency claim is available and a null
is uninterpretable. Leg 1 can at most bound an increment on top of an unremoved
autocrine source.

**For the premise:** the recipient is equipped. Integrin alphaV, ITGB1, ITGB8 and
LTBP1 are all expressed in these fibroblasts, and so is EGFR. These are the first
measurements in this repository showing the receiver carries the machinery the
mechanism needs. That bears on C36 without settling it, since transcript is not
protein.

## What the screen cannot do, stated before it runs

- **It cannot test geometry or delivery.** Every cell in a well is in one Matrigel drop.
- **It cannot support a necessity claim**, because the recipient supplies the ligand.
- **It cannot support absence in any arm.** Every recorded reduction is partial: Areg
  -1.042, Egfr -1.034, Erbb3 -0.511, Itgb6 -1.138, Hbegf -1.921, and Erbb4 has no
  receptor to remove. The receptor arm therefore contributes only as a consistency
  check, and no null in it can exclude a role.
- **It cannot carry randomization-based inference.** Fifty of the 53 plate-3 targets
  sit at one fixed well position in all four units, Areg always at F07, so target is
  confounded with position and guide pool and the units are copies of one layout.
  Plate 3 is forced by the design, since the remaining axis targets all sit there; it
  was not selected.
- **It cannot measure secreted ligand.** Whether the fibroblast receives the signal is
  still C36, Not established.
- **It cannot be independent replication** while A10's preparation identity is
  unresolved.
- **Two assumptions are carried:** that the mouse amphiregulin EGF-like domain
  activates human EGFR, and that a lung fibroblast stands in for the pericyte of the
  cited mechanism. Neither is established here.
- **The culture medium and matrix are not in the deposit.** Exogenous EGF would
  saturate the recipient's receptor and a TGF-beta receptor inhibitor would clamp the
  endpoint, so a null is pharmacologically ambiguous until the methods are checked.

## Why leg 2 uses a cohort that already returned a null

Trial E6 tested the abundance version in GSE136831
([Adams et al. 2020](https://doi.org/10.1126/sciadv.aba1983)) and found nothing, with
a depth control that refused the one significant pair. Leg 2 reuses that instrument,
with the same donors, fibroblast label rule, 50-cell floor and depth rule, and changes
only the predictor.

The reuse is honest, but the leg is now exploratory rather than a declared test, for a
reason already in the logged table: the activation score correlates with fibroblast
depth at rho 0.4116 and fibroblast EGFR at 0.4918, and both members of the new pair
are fibroblast detection fractions, so the inherited gate is likely to refuse it. A
refusal is a result and is reported as one. A second rival must also be named: TGFB1,
THBS1, ITGB8 and LTBP1 are themselves TGF-beta inducible, so predictor and outcome are
co-regulated programme members in the same cells, and a detection control cannot bound
that. A co-regulation control is required alongside it.

## Connections and boundaries

A9 owns the fibroblast EGFR itself: its abundance, its complex composition and whether
it is engaged. Leg 2 does not measure EGFR, so it tests the layer below the receptor
and cannot settle A9, and A9 cannot settle it. The cited mechanism makes that receptor
layer obligatory, so the boundary between the two questions runs through the middle of
the mechanism, and this rationale states that rather than hiding it. A12 owns
interleukin-1 recipients and is untouched here.

A10 owns the organoid growth outcome and its own modules. This plan uses the same
deposit for a different endpoint and does not reinterpret A10's fits. It inherits
A10's prohibition in A10's own words: a fibroblast transcriptome change may not be read
as fibroblast causation, the Areg contrast included, because that needs a
fibroblast-side perturbation. Nothing here re-grades C37, C49 or C50, and no claim row
is added; grading is the owner's decision.

The screen's source paper is roadmap paper 14, DOI 10.1073/pnas.2606113123, and the
owner's reading is recorded as pending. This plan uses the deposit and the
repository's own records, writes no study note and accepts none of that paper's
claims, following DEVELOPMENT decision 21. One bounded design fact is reused from
A10's own source design check: TIGIT is the article's in-plate control.

## Other genomic layers, with feasibility verdicts

1. **Spatial transcriptomics of human fibrotic or tumour lung.** Now the decisive
   layer rather than an extension, because it is the only one where proximity varies,
   so it is the only one that can separate delivery from abundance. Ask whether
   fibroblast TGF-beta targets rise with proximity to AREG-high epithelium rather than
   with tissue-level AREG. Verdict: **feasible, pending a cohort audit**, with a frozen
   proximity definition and a donor-level unit declared before any distance is computed.
2. **A fibroblast-side perturbation.** The only way to bound the autocrine source this
   screen leaves in place. Verdict: **not feasible from public data**; it needs a
   co-culture where the recipient's own ligand can be removed.
3. **Secreted and shed ligand by targeted protein measurement.** Decisive for C36 and
   for the ligand-availability rival. Verdict: **not feasible here**; bench work.
4. **Fibroblast surface protein by CITE-seq or equivalent.** Addresses C36 and A9
   directly. Verdict: **unverified feasibility**; needs a lung cohort with both
   modalities and adequate fibroblast depth.
5. **Chromatin accessibility at fibroblast TGF-beta response elements.** Tests
   licensing at the regulatory layer. Verdict: **unverified feasibility**, conditional
   on a public lung multiome with fibroblast depth; MC2 already warns these assays are
   not interchangeable.
6. **Sequence layer, cross-species ligand activity.** Aligning the mouse and human
   amphiregulin EGF-like domains bounds one carried assumption. Verdict: **feasible and
   small**, from UniProt P31955 and P15514; stage 1 could not do it from local
   resources and it remains open.
