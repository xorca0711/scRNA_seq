# A15 rationale: why the epithelial integrin arm is a separate question

Declared 27 September 2026, as a proposal. This is the biological argument for a
question, written after the observation that motivates it, and it says so rather than
reading as a pre-registration. Nothing here is a result of this repository. Article
details below were retrieved from PubMed; each is cited with its DOI.

## The proposition, stated so it can fail

**Where an epithelium and a fibroblast share a matrix, the epithelial contribution to the
fibroblast myofibroblast and collagen programme is carried by epithelial integrin
alphaVbeta6 converting latent TGF-beta that is already present, rather than by epithelial
supply of an EGFR ligand.**

Three predictions follow, and each can fail on its own:

1. **Direction and specificity.** Removing epithelial ITGB6 lowers the fibroblast
   activation programme; removing epithelial ligand supply does not, in a system where
   the recipient's own ligand has also been removed so that the null is bounded.
2. **The readout.** The same removal lowers **activated** TGF-beta at the recipient,
   measured as protein or as receptor-proximal signalling, not as a transcript score.
3. **Mediation.** Blocking TGF-beta receptor signalling in the recipient abolishes the
   epithelial ITGB6 effect. If it does not, the epithelial integrin acts through
   something else and the proposition is wrong even if prediction 1 holds.

Prediction 3 is the one that separates this question from a correlation. An epithelial
perturbation that changes a fibroblast is not evidence about TGF-beta activation unless
the recipient's TGF-beta receptor is shown to be the route.

## The mechanism, and how old it is

This is not a new mechanism, and the rationale states that at the top rather than in a
limitation.

- **Integrin alphaVbeta6 binds the TGF-beta1 latency-associated peptide and activates
  latent TGF-beta1 in a spatially restricted way.** Cells expressing the integrin
  activate the cytokine locally, and mice lacking the integrin are protected from
  pulmonary fibrosis
  ([Munger 1999, Cell](https://doi.org/10.1016/s0092-8674(00)80545-0)). The spatial
  restriction is the part that matters for this repository: an integrin on one cell
  activates a signal on an adjacent cell, so the mechanism is short range by construction
  ([Sheppard 2005](https://doi.org/10.1007/s10555-005-5131-6)).
- **The integrin is epithelial and injury-induced in human lung.** AlphaVbeta6 is
  overexpressed in human pulmonary fibrosis within pneumocytes lining the alveolar ducts
  and alveoli, and a blocking antibody prevents murine bleomycin fibrosis, with partial
  inhibition attenuating collagen without exacerbating inflammation
  ([Horan 2008](https://doi.org/10.1164/rccm.200706-805OC)).
- **Removing it has consequences beyond the fibroblast.** Itgb6-null mice develop
  age-related, Mmp12-dependent emphysema, rescued by transgenic active TGF-beta1
  ([Morris 2003, Nature](https://doi.org/10.1038/nature01413)), and they accumulate
  bronchoalveolar phospholipids and collectins with abnormal alveolar macrophages
  ([Koth 2007](https://doi.org/10.1165/rcmb.2006-0428OC)). An Itgb6 loss is therefore not
  a clean single-pathway perturbation of the epithelium, which is rival 2 below.
- **Two sources of active TGF-beta are proposed, not one.** The epithelial alphaVbeta6
  step is held to drive fibroblast expansion and activation, with an amplification step in
  which an alphaV integrin on the activated fibroblast itself activates further TGF-beta;
  deleting the alphaV subunit from activated fibroblasts protected mice from pulmonary,
  hepatic and renal fibrosis despite intact epithelial alphaVbeta6
  ([Sheppard 2015](https://doi.org/10.1513/AnnalsATS.201406-245MG)). This is rival 1, and
  it is the strongest rival this question has.
- **The activation step is regulated from the epithelial side.** Lysophosphatidic acid
  induces alphaVbeta6-mediated TGF-beta activation in human epithelial cells through LPA2,
  G alpha q, RhoA and Rho kinase, and both the receptor and the integrin are upregulated
  over fibrotic areas in human usual interstitial pneumonia
  ([Xu 2009](https://doi.org/10.2353/ajpath.2009.080160)). The epithelial integrin is a
  regulated node, not a constitutive one.

## What A2 measured, and exactly how far it goes

The observation is A2's, computed under a freeze that forbids a p-value. It is restated
here with its limits so that this workspace cannot be read as owning it.

| Epithelial knockout | Median depth-matched effect on the activation score | Units in the predicted direction | Culture unchanged |
|---|--:|---|---|
| **Itgb6** | **-0.938** | 3 of 3 readable | yes |
| Areg | +0.036 | 2 of 4 | yes |
| Erbb3 | -0.187 | 4 of 4 | no |
| Egfr | +0.247 | 0 of 3 | yes |
| Erbb4 | +0.503 | 0 of 3 | yes |
| Erbb2 | -0.169 | 1 of 1 readable | no |

The endpoint's standard deviation over the 168 eligible plate-3 wells is 0.482, so the
Itgb6 effect is about 1.9 standard deviations and the Areg effect about 0.07 in the wrong
direction. Per unit, Itgb6 gives -0.938, -1.575 and -0.680, with its fourth well
ineligible at 4,865 fibroblast counts. Source:
[A2 stage 3](../A2_areg_source_delivery/reports/STAGE3_LEG1_RESULTS.md).

**Five limits travel with it, and none of them is decorative.**

1. **Target is confounded with position and guide pool.** Fifty of the 53 plate-3 targets
   sit at one fixed well position in all four units. The four ITGB6 libraries are
   plate3-1 through plate3-4 at well **C02**, from the deposit's own sample titles
   (GSM9216275, GSM9216335, GSM9216395, GSM9216455). No randomization-based inference is
   available in this screen, which is why A2's first freeze was withdrawn.
2. **The ligand arm is a bounded null, not an absence.** These fibroblasts transcribe AREG
   at 9.846 mean log2 CPM in 99.2 per cent of plate-3 wells, above the 7.267 the
   epithelial knockout removes within its own compartment, and the knockout is partial.
   The comparison at the heart of this question therefore sets one measured decrease
   against one null that cannot exclude a role.
3. **Preparation independence is unresolved.** A10 stage 1 could not establish whether the
   replicate indices are independent preparations, so everything from this screen is a
   within-screen association.
4. **The knockout is partial.** The recorded mouse Itgb6 transcript reduction is 1.138
   log2 CPM. No arm of this screen can support absence.
5. **The readout is RNA.** The proposition names a protein event. A five-gene transcript
   score in the recipient is consistent with it and does not measure it.

**One check the observation does pass, and it is worth stating.** Erbb3 is the other
consistent decrease, and it travels with organoids 0.605 log units smaller and 2.31 times
the fibroblast material, which is what a viability or composition shift looks like. Itgb6
travels with organoid size within 0.022 log units of controls and fibroblast content at
0.93 of theirs. The generic-perturbation caution therefore weighs on Erbb3 rather than on
Itgb6, without being removed from either.

## Rivals

1. **Fibroblast-side amplification carries the effect.** If the recipient's own alphaV
   integrins do most of the activating, the epithelial step only seeds a loop the
   fibroblast runs, and the epithelial integrin is not the rate-limiting input
   ([Sheppard 2015](https://doi.org/10.1513/AnnalsATS.201406-245MG)). Not separable
   without a fibroblast-side perturbation.
2. **The epithelial state, not the delivered ligand.** **Restated 27 September 2026, and the
   original wording was wrong.** It read: Itgb6 loss changes the epithelium itself, so a
   changed epithelium could change a fibroblast "for reasons unrelated to TGF-beta
   activation". That does not follow from the sources cited for it, because **both of them
   attribute their phenotypes to loss of TGF-beta activation, in their own titles**:
   [Morris 2003](https://doi.org/10.1038/nature01413) is "Loss of integrin
   alpha(v)beta6-mediated TGF-beta activation causes Mmp12-dependent emphysema" and
   [Koth 2007](https://doi.org/10.1165/rcmb.2006-0428OC) is "Integrin beta6 mediates
   phospholipid and collectin homeostasis by activation of latent TGF-beta1". What they
   establish is that the epithelial integrin's TGF-beta activation has consequences **inside
   the epithelium**, not that a TGF-beta-independent route exists.

   The rival, correctly stated, is therefore: **the epithelial integrin's effect on the
   fibroblast may be indirect, mediated by a change in the epithelium rather than by TGF-beta
   activated at the epithelial surface acting on the fibroblast.** The route may still be
   TGF-beta; what is in question is what it acts on first. So the separable question is
   whether the epithelium changes at all, which is an omnibus question and needs no curated
   gene set. This restatement was forced by an adversarial review of the side-branch freeze,
   and it is recorded in
   [reports/RIVAL2_FREEZE_V1_WITHDRAWN.md](reports/RIVAL2_FREEZE_V1_WITHDRAWN.md) rather than
   applied silently. It is the one rival a public dataset can address today; see the layer
   proposals below.
3. **Generic perturbation.** Some part of any epithelial knockout effect may be generic.
   Erbb3 fell in all four units, and the panel is not uniformly downward, since Egfr and
   Erbb4 rose.
4. **Position and guide pool.** Limit 1 above is also a rival: the C02 well, or the guide
   pool assigned to it, could carry the effect.
5. **Composition and depth.** The eligibility floor and the epithelial-fraction
   sensitivity reduce this rival and do not remove it; about a fifth of reads are
   unassigned to either species.
6. **The latent pool, not the activator.** Latent TGF-beta in the well may come from
   either compartment or from the medium, which is not in the deposit. AlphaVbeta6
   activates TGF-beta1 and TGF-beta3 but not TGF-beta2
   ([Sheppard 2015](https://doi.org/10.1513/AnnalsATS.201406-245MG)), so isoform and
   source both matter and neither is measured.
7. **Cross-species activation is assumed.** The integrin is mouse and the recipient is
   human. Whether mouse alphaVbeta6 activates latent TGF-beta1 in this culture, and
   whether the activated ligand signals on human receptors, is carried and unverified.
   Human and mouse TGF-beta1 proproteins are both 390 residues
   (UniProt [P01137](https://www.uniprot.org/uniprotkb/P01137),
   [P04202](https://www.uniprot.org/uniprotkb/P04202)), and integrin beta-6 is 788 and 787
   residues (UniProt [P18564](https://www.uniprot.org/uniprotkb/P18564),
   [Q9Z0T9](https://www.uniprot.org/uniprotkb/Q9Z0T9)). No alignment was performed here,
   and equal length is not conservation of the activating motif.

## The clinical record, which bounds any answer this question could give

This belongs in the rationale rather than in a discussion section, because it changes what
a positive result would be worth.

- **An anti-alphaVbeta6 monoclonal antibody failed in idiopathic pulmonary fibrosis.**
  BG00011 was tested against placebo in a phase 2b trial that terminated early for an
  imbalance in adverse events and a lack of clinical benefit. At week 26 there was no
  significant difference in forced vital capacity change; after week 26 the treated group
  trended worse, more participants showed worsening fibrosis on imaging, and serious
  adverse events, including four deaths, were more frequent
  ([Raghu 2022](https://doi.org/10.1164/rccm.202112-2824OC)).
- **A dual alphaVbeta6 and alphaVbeta1 inhibitor was better tolerated.** Bexotegrast
  reached a favourable safety profile over 12 weeks in a phase 2a trial, with exploratory
  signals on forced vital capacity, quantitative lung fibrosis imaging and circulating
  fibrosis biomarkers ([Lancaster 2024](https://doi.org/10.1164/rccm.202403-0636OC)).

Two readings are available and this rationale commits to neither. Blocking the epithelial
activator systemically may remove a homeostatic function of TGF-beta that the lung needs,
which is what the Itgb6-null emphysema and surfactant phenotypes would predict. Or the
epithelial arm may not be the rate-limiting one in established human disease, which is
what rival 1 would predict. Either way, **a positive answer to A15 in a co-culture would
be a statement about which epithelial output moves a fibroblast programme, not about
whether blocking that output helps a patient.**

## What is genuinely open, given all of the above

The mechanism is established. What is not established, and is what this question owns, is
the **partition of the epithelial contribution**: in one system, with the recipient's own
ligand removed so that the ligand arm is bounded, and with activated TGF-beta measured,
does the integrin arm carry the epithelial input while the ligand arm does not? No cited
source answers that, because the cited work compares epithelial activation against
fibroblast activation (rival 1) or against no perturbation, not against epithelial ligand
supply. A2 attempted the comparison and could produce only one measured arm and one
unbounded null.

## Connections and boundaries

- **A2** owns amphiregulin, delivery against abundance, and the grading of the founding
  observation. A15 does not re-grade it, does not reuse its claim identifiers, and takes
  nothing from it as evidence. If A2's proposal 2 is graded, that row stays A2's.
- **A9** owns the fibroblast EGFR: its abundance, its complex composition and whether it
  is engaged. A15 measures no receptor on the fibroblast except where prediction 3
  requires the TGF-beta receptor, which is a different receptor and a different question.
- **A10** owns the organoid growth outcome and its own modules. A15 uses no A10 result as
  evidence and inherits A10's prohibition in A10's words: a fibroblast transcriptome
  change may not be read as fibroblast causation.
- **A12 and A13** own interleukin-1 recipients and fibroblast programmes beyond macrophage
  IL1B, and are untouched here.
- **The register's own rule applies to this card.** A new biological interpretation does
  not inherit a historical claim's status, so A15 begins with no evidence of its own.

## Other genomic layers, with feasibility verdicts

1. **The epithelial translatome under integrin beta6 blockade in vivo.** GSE190821 blocks
   integrin beta6 with the 3G9 antibody in bleomycin-injured mice and reads an epithelial
   RiboTag immunoprecipitation alongside a paired whole-lung input, with the mouse as the
   unit. It cannot test the hypothesis: the fibroblast compartment is not separated and there
   is no activation readout. It can address **rival 2** as restated above, by asking whether
   the epithelium changes at all when the integrin is blocked. Verdict: **authorized and
   executed, 27 September 2026**, under
   [`config/a15_rival2_freeze_v2.json`](config/a15_rival2_freeze_v2.json); read
   [the results](reports/RIVAL2_RESULTS.md). The arms are four treated mice against four
   inert-antibody mice, not four against seven: only four of the seven Axum8 mice are
   bleomycin-exposed, and the other three are the saline context arm. The caveats stand: the
   antibody is systemic rather than epithelium-restricted, bleomycin day 7 in vivo is not a
   two-week co-culture, and four against four reaches nominal significance only under
   complete separation.
2. **A whole-lung activation signature in the same deposit.** Verdict: **feasible but
   non-discriminating, and not recommended.** Whole lung confounds the fibroblast
   compartment with fibrosis extent, and reduced collagen under this antibody is already
   published ([Horan 2008](https://doi.org/10.1164/rccm.200706-805OC)), so the run would
   re-measure a known result.
3. **Spatial transcriptomics of human fibrotic lung.** Ask whether a fibroblast TGF-beta
   response rises with proximity to ITGB6-high epithelium rather than with tissue-level
   ITGB6. Verdict: **feasible, pending a cohort audit**, with a frozen proximity
   definition and the donor as the unit. This shares A2's own spatial proposal and must be
   counted once, under whichever question runs it.
4. **Activated TGF-beta in the co-culture, measured as protein.** A bioassay or
   receptor-proximal signalling measurement on the fibroblast compartment, with more than
   one well per target and positions that vary. Verdict: **not feasible from public
   data**; this is the decisive test and it is a bench design.
5. **Recipient-side TGF-beta receptor blockade in the same co-culture.** Verdict:
   **bench**. Decisive for prediction 3, and without it a positive prediction 1 stays a
   correlation.
6. **Fibroblast-side alphaV deletion in the same co-culture.** The only way to bound rival
   1. Verdict: **bench**.
7. **Chromatin accessibility at fibroblast SMAD response elements.** Tests licensing at
   the regulatory layer. Verdict: **unverified feasibility**, conditional on a public lung
   multiome with fibroblast depth; MC2 already warns that these assays are not
   interchangeable.
8. **Sequence layer, cross-species activation.** Align the mouse and human TGF-beta1
   latency-associated peptides at the integrin-binding motif, and the mouse and human
   integrin beta-6 ligand-binding regions, to bound rival 7. Verdict: **feasible and
   small**, from UniProt P01137, P04202, P18564 and Q9Z0T9; not done here and still open.
   A2 has the same item open for amphiregulin.
9. **Targeted proteomics of the latent complex in the co-culture.** Whether the latent pool
   the integrin would act on is present, and from which compartment. Verdict: **not
   feasible from public data**.
