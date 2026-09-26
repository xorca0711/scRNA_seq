# A2 rationale: why AREG source matters through delivery, not through abundance

Declared 26 September 2026, before any endpoint value was computed in either
dataset named below. This is the biological argument for the rewritten
[A2 card](../../RESEARCH_QUESTIONS.md#a2). The machine-readable decisions are in
[`config/a2_delivery_contract.json`](config/a2_delivery_contract.json) and the
stages are in [PLAN.md](PLAN.md).

## The question the register used to ask, and why it is closed

A2 previously asked which compartment is the dominant AREG source, and whether
that dominance changes with tissue context. The repository has answered the
abundance version of that question seven times, across lung adenocarcinoma and
pulmonary fibrosis cohorts, and no hierarchy survived. Each row below states the
register's own proposition, not a conclusion drawn from it:

| Register proposition | Claim | Status and what the numbers were |
|---|---|---|
| Epithelial AREG detection exceeds myeloid detection in tumour lung, independently of annotation and molecule depth | C37 | Descriptive only. Medians 0.4591 against 0.2715, all ten donors epithelial-higher before matching, six of ten at the primary 1,000-UMI budget |
| Epithelium is the AREG source in human lung adenocarcinoma | C39 | Not established. CD1c-positive dendritic cells reach 0.618 against tumour states at 0.555 and 0.532 |
| Dendritic cells and monocytes are a major AREG and HBEGF source, at or above the epithelial states | C45 | Descriptive only, in two diseases and three deposits |
| AREG or HBEGF is tumour-enriched against matched normal lung from the same donors | C40 | Refuted at compartment level |
| AREG or HBEGF is enriched in fibrotic lung against control lung | C48 | Not established |
| Epithelial AREG and fibroblast EGFR are correlated across donors | C49 | Not established. Spearman rho 0.348, p 0.112, in 22 donors |
| Epithelial AREG and a fibroblast activation score are correlated across donors | C50 | Not established. Spearman rho -0.150, with the sign reversed |

C51 is the decisive methodological entry. The only significant correlation in
that donor-level trial was a control pair whose two members both tracked
sequencing depth. A frozen rule refused it. The same number on AREG would have
been read as the discovery the paper predicted.

Two readings of that record are possible. The first is that the axis is not
real. The second is that total ligand abundance is the wrong observable, so
every abundance comparison was measuring a quantity the biology does not use.
The mechanism literature supports the second reading, and it is testable.

## The mechanism that makes abundance the wrong observable

Amphiregulin does not act as a diffusible growth factor whose tissue
concentration sets a response. It is synthesized as a membrane-anchored
precursor, and its recorded action on mesenchymal cells is to convert a store
the recipient already holds.

- **AREG releases the recipient's own latent TGF-beta.** In two acute injury
  models, amphiregulin freed bioactive TGF-beta from latent complexes by
  activating integrin alphaV on mesenchymal stromal cells, and that activation
  drove their differentiation into myofibroblasts. Macrophages were identified
  as a critical source
  ([Minutti et al. 2019, Immunity](https://doi.org/10.1016/j.immuni.2019.01.008)).
- **The fibroblast arm of TGF-beta signalling needs AREG.** Silencing
  amphiregulin, or inhibiting EGFR, reduced TGF-beta1-driven fibroblast
  proliferation, smooth-muscle actin and collagen, and reduced collagen
  accumulation in TGF-beta1 transgenic lungs
  ([Zhou et al. 2012, JBC](https://doi.org/10.1074/jbc.M112.356824)).
- **Immune AREG is non-redundant, not merely present.** Removing amphiregulin
  from regulatory T cells alone caused severe lung damage during influenza with
  no change in suppressor function, antiviral response or viral load
  ([Arpaia et al. 2015, Cell](https://doi.org/10.1016/j.cell.2015.08.021)).
  Type 2 innate lymphoid cells drive primary lung fibroblast proliferation and
  differentiation partly through amphiregulin and EGFR
  ([Sorkhdini et al. 2024, JCI Insight](https://doi.org/10.1172/jci.insight.178381)).
- **Both epithelial and leukocyte sources are established**, which is why C45 is
  a constraint on an epithelium-only reading rather than a refutation of the
  axis ([Zaiss et al. 2015, Immunity](https://doi.org/10.1016/j.immuni.2015.01.020)).

A short-range ligand that converts a recipient-held store predicts the nulls
this repository recorded. If the rate-limiting quantity is the recipient's
latent complex and its integrin, then donor-level ligand abundance should not
predict fibroblast activation, and C49 and C50 are what the mechanism expects
rather than evidence against it. That is a prediction, not a rescue: it commits
the hypothesis to a different observable, and that observable can fail.

## The hypothesis, stated so it can fail

**AREG's contribution to a fibroblast response is set by where the ligand is
released relative to a competent recipient, not by which compartment
transcribes the most of it.**

Two legs follow, and each is tested separately because each can fail alone.

**Leg 1, delivery necessity.** Where an epithelial source and a fibroblast are
in direct contact, removing the epithelial ligand lowers a frozen fibroblast
TGF-beta response programme. Removing the epithelium's own receptors does not,
because epithelial reception is not what the fibroblast receives.

**Leg 2, recipient licensing.** Across donors, the fibroblast response tracks
the recipient's own TGF-beta activation machinery more closely than it tracks
epithelial ligand abundance.

A result that would weaken the hypothesis, stated before either test: if
removing the epithelial ligand leaves the fibroblast programme unchanged while
removing epithelial receptors moves it, the fibroblast is responding to the
epithelial state and not to delivered ligand, and the delivery framing is wrong
for this system. If leg 2 shows that ligand abundance predicts the fibroblast
response as well as recipient machinery does, the abundance rival is not
excluded and the reframing bought nothing.

## Why the organoid screen can separate the rivals

GSE307112 cultures mouse alveolar type 2 cells with human lung fibroblasts and
introduces CRISPR knockouts **in the mouse epithelium only**. The fibroblasts
are unedited, and reads are assigned by species. That asymmetry, not the number
of wells, is what makes the design informative for A2:

| Epithelial knockout | What it removes | What the delivery hypothesis predicts for the fibroblast endpoint |
|---|---|---|
| Areg | the ligand itself, for autocrine and paracrine use alike | lower |
| Egfr, Erbb2, Erbb3, Erbb4 | epithelial reception only; the ligand is still made and released | not systematically lower |
| Itgb6 | epithelial integrin-mediated TGF-beta activation | not required to be lower, because the hypothesis places the activating integrin on the recipient |
| Hbegf | a different EGFR ligand ranked below AREG in this repository (C32) | weaker or absent effect, if the axis is ligand-specific |

The autocrine rival makes the opposite prediction for the receptor row. The
epithelial-activation rival makes the opposite prediction for the Itgb6 row.
Neither rival is excluded by any result already recorded, and both are named on
the card.

Plate 3 carries Areg, Egfr, Erbb2, Erbb3, Erbb4 and Itgb6, one well each in each
of its four plate-replicate units, alongside two control wells per unit. The
comparison is therefore matched within a plate, which matters because A10 found
that plate-level behaviour differs sharply across this screen. Hbegf sits on
plate 4, so the ligand-specificity contrast crosses plates and is weaker by
construction; the plan reports it as secondary and says so.

## What the screen cannot do, stated before it runs

- **It cannot test geometry.** Every cell in a well is in the same Matrigel
  drop. The screen tests one prediction of the delivery hypothesis, that a
  contact-proximal epithelial source is necessary, and says nothing about
  distance, gradient or spatial arrangement in a tissue.
- **It cannot measure secreted ligand.** Transcript loss is not protein loss.
  Whether the fibroblast receives the signal is still C36, Not established.
- **It cannot support absence.** The recorded Areg knockout lowers the mouse
  transcript by 1.042 log2 CPM and leaves it at 6.226, so the perturbation is
  partial and pooled. A null is inconclusive by design, and the contract says so
  before the test rather than after it.
- **It cannot be independent replication.** A10 stage 1 could not establish
  whether the 15 plate-replicate labels are separate epithelial isolations and
  fibroblast lots. Until that is resolved, every result here is a within-screen
  association. That restriction is inherited, not new.
- **One assumption is carried, not established.** The ligand is mouse and the
  receptor is human. Cross-species activity of the amphiregulin EGF-like domain
  on human EGFR is assumed by this design. Stage 1 records the assumption and
  attempts a sequence-level check; if that check cannot be run offline, the
  assumption is carried explicitly into every reading.

## Why leg 2 uses a cohort that already returned a null

Trial E6 tested the abundance version in GSE136831
([Adams et al. 2020](https://doi.org/10.1126/sciadv.aba1983)) and found nothing,
with a depth control that refused the one significant pair. Leg 2 reuses that
instrument unchanged, with the same donors, the same fibroblast label rule, the
same 50-cell floor and the same depth rule, and changes only the predictor: the
recipient's own activation machinery in place of the epithelial ligand.

This is honest reuse rather than a second chance, for three reasons. The cohort
is not being re-searched for a positive result on the same variables; the
estimand is new and follows from a mechanism published before the trial ran; and
the depth rule that refused the earlier pair applies unchanged and may refuse
this one too. Both members of the new pair are measured in the same fibroblasts,
so shared detection is the dominant artefact and refusal is a plausible outcome.
A leg that reports "refused by its own control" has produced a result.

## Connections and boundaries

A9 owns the fibroblast EGFR itself: its abundance, its complex composition and
whether it is engaged. Leg 2 does not measure EGFR at all. It asks about the
machinery the recipient uses after engagement, the integrins and the latent
TGF-beta complex, so it cannot settle A9 and A9 cannot settle it. A12 owns
interleukin-1 recipients and is untouched here. A10 owns the organoid growth
outcome and its own modules; this plan uses the same deposit for a different
endpoint and does not reinterpret A10's fits. Nothing here re-grades C37, C49 or
C50, and no claim row is added; grading is the owner's decision.

The screen's source paper is roadmap paper 14, DOI 10.1073/pnas.2606113123, and
the owner's reading is recorded as pending. This plan uses the deposit and the
repository's own records, writes no study note and accepts none of that paper's
claims, following DEVELOPMENT decision 21.

## Other genomic layers, with feasibility verdicts

1. **Secreted and shed ligand, by targeted protein measurement.** This is the
   decisive measurement for the ligand-availability rival and for C36. Verdict:
   **not feasible here.** It needs bench work on conditioned medium; no public
   deposit substitutes for it.
2. **Spatial transcriptomics of human fibrotic or tumour lung.** Tests the
   geometry claim the organoid cannot: whether fibroblast TGF-beta targets rise
   with proximity to AREG-high epithelium rather than with tissue-level AREG.
   Verdict: **feasible, pending a cohort audit.** Public lung spatial cohorts
   exist; a frozen proximity definition and a donor-level unit must be declared
   before any distance is computed.
3. **Surface protein on fibroblasts, by CITE-seq or equivalent.** Addresses C36
   and A9 directly, since EGFR transcript detection is a poor proxy for receptor
   protein. Verdict: **unverified feasibility.** It requires a lung cohort with
   both modalities and enough fibroblasts per donor, which has not been searched.
4. **Chromatin accessibility at fibroblast TGF-beta response elements.** Tests
   licensing at the regulatory layer rather than at the transcript layer.
   Verdict: **unverified feasibility,** conditional on a public lung multiome
   with fibroblast depth; MC2 already warns that accessibility, histone marks
   and methylation measure different properties and are not interchangeable.
5. **Sequence layer, cross-species ligand activity.** The cheapest open item: an
   alignment of the mouse and human amphiregulin EGF-like domains bounds the
   assumption leg 1 carries. Verdict: **feasible and small,** and stage 1 owns
   it.
