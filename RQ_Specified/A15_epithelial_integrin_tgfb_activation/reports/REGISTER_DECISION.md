# Does the Itgb6 lead deserve its own register question

27 September 2026. A2 closed with one measured lead that its own card says needs its own
design, and left the registration to the owner. This report is the argument, not the
decision. It recommends a new identifier, states what that identifier does not buy, and
records the three alternatives that were considered and rejected, with reasons.

**Recommendation: register it as A15, narrowly scoped to the partition question, with
readiness declared blocked.** The owner retains or rejects the registration; no claim row
is added either way, and the register ends at C168.

## The criteria used

They are the repository's own, taken from the rules already in the tree rather than
invented here.

| Criterion | Source | Met? |
|---|---|---|
| It carries its own biological proposition | the register's card format | Yes. Epithelial integrin-mediated activation of latent TGF-beta, not amphiregulin |
| It carries its own decision, not a shared one | `RQ_Specified/README.md`: an entry owned by two questions is enabling work and takes no new identifier | Yes. Its decision is which epithelial output moves the fibroblast programme; A2's is delivery against abundance |
| It has its own biological units and multiplicity family | same | Yes, and they do not exist yet, which is the point of the gate in PLAN.md |
| An identifier is not an evidence grade | `RQ_Specified/README.md` | Relied on. Registering it grades nothing |
| A new interpretation does not inherit a historical claim's status | `RESEARCH_QUESTIONS.md` preamble | Relied on. A15 opens with no evidence of its own |

## The case for

1. **Keeping it inside A2 would corrupt A2's estimand.** A2 is now a question about
   amphiregulin: whether its contribution is set by delivery to a competent recipient or
   by source abundance. Integrin-mediated activation of latent TGF-beta is not
   amphiregulin biology and does not answer that question. A non-AREG mechanism riding on
   A2's card is exactly the drift the register's rewrite was meant to stop.
2. **The observation is the largest effect this repository has produced in this screen,
   and an unhomed lead gets re-cited loosely.** At about 1.9 endpoint standard deviations,
   surviving the epithelial-fraction adjustment and a culture check that the neighbouring
   Erbb3 decrease failed, it will be quoted. A card with its own prohibitions is how it
   gets quoted with its limits attached.
3. **The register already hosts questions whose test cannot run.** A3, A4, A7, A8, A12 and
   A14 all state a missing design or a missing cohort under readiness. Blocked readiness
   is a recognised state here, not a reason to refuse an identifier.
4. **The partition is genuinely open in the literature.** The mechanism is not, and the
   rationale says so. But the cited work compares epithelial activation against fibroblast
   activation, or against no perturbation. It does not compare epithelial ligand supply
   against epithelial integrin activation in one system with an autocrine-competent
   recipient. That comparison is what A15 owns.

## The case against, which is real

1. **Nothing can be run.** No public dataset supplies a TGF-beta activation readout under
   an integrin perturbation. A15 is a parked question from the day it opens, and parked
   questions consume attention at every status scan.
2. **The mechanism is 27 years old and clinically probed.** Integrin alphaVbeta6 activating
   latent TGF-beta1 has been established since 1999, a blocking antibody prevented murine
   bleomycin fibrosis in 2008, and two clinical programmes have tested the axis in
   idiopathic pulmonary fibrosis. A register question is not a place to re-ask a settled
   mechanism, and the framing must be the partition or it adds nothing.
3. **The founding observation is thin.** Three readable units, one well per target per
   unit, all at well C02, target confounded with plate position and guide pool, no
   randomization-based inference, preparation independence unresolved, the knockout
   partial, and the readout a transcript score rather than the protein event proposed.
4. **The comparison is asymmetric.** The ligand arm is a bounded null, not a measured
   absence, because the recipient transcribes AREG above the level removed. So the headline
   contrast is one measured decrease against one uninterpretable null.

## How the case is resolved

The case against argues for a **narrow scope and a blocked readiness**, not against an
identifier. Two constraints were therefore written into the card and the contract:

- The question is the **partition of the epithelial output**, not whether the integrin
  activates TGF-beta. The latter is established elsewhere and this repository does not
  re-ask it.
- Readiness is **blocked**, with the binding constraint named as the absence of a
  TGF-beta activation readout in any public deposit, and the decisive test named as a
  bench design.

Objection 3 is not resolved by scoping and is instead recorded: the founding observation
stays A2's, it is not evidence for A15, and A15 opens with none. That is why this
workspace has no `tables/` directory.

## Alternatives considered and rejected

| Alternative | Why it was rejected |
|---|---|
| **A third leg of A2** | It is not an amphiregulin question, and A2's second freeze forbids adding an arm after the endpoint has been read. A2 reported all six stages complete; reopening it to host a different mechanism would put a post hoc arm inside a closed pre-registration |
| **Fold it into A9** | A9 owns the fibroblast EGFR: receptor abundance, complex composition and engagement. A15's receptor is the recipient's TGF-beta receptor, and its perturbation is on the epithelium, not the fibroblast. The two questions share a compartment, not an estimand |
| **An enabling entry, like A12-S1** | An enabling entry supplies a shared input to questions that keep their own decisions. A15 has a biological proposition and a decision of its own, which is the register's own test for a full identifier |
| **A lead in FINDINGS.md only, with no card** | A lead with no card has no prohibitions attached. The specific risks here are naming the transcript score as a TGF-beta activation measurement and reading an unrandomized single-well contrast as a mechanism, and both are exactly what a card can forbid |

## What was deliberately not done

- **No claim row was added**, and no existing row was re-graded. A2's stage 5 already
  proposed wording for the observation as its own proposal 2; that row, if it is ever
  created, belongs to A2.
- **No file under `RQ_Specified/A2_areg_source_delivery/` was edited.** Those files are
  under the owner's review in pull request 85, and adding a cross-reference into them
  while they are being reviewed would change a document under review. A15 is discoverable
  from the register table and from the `RQ_Specified` index instead.
- **No study note was written on the screen's source paper**, roadmap paper 14
  (DOI 10.1073/pnas.2606113123), whose reading by the owner is recorded as not started.
- **The sibling deposit's series summary was not used.** GSE307128, in the same
  superseries, carries a summary that states the source paper's conclusions. It was read
  while searching GEO and is deliberately not quoted or relied on; only its structural
  facts, the accession, four samples, mouse, appear in the search report.

## What the owner is being asked to decide

1. **Register A15, or refuse it.** If refused, the honest alternative is a FINDINGS entry
   under A2 with the same limits stated, and this folder is removed.
2. **The scope wording**, if registered: the partition of the epithelial output, with the
   mechanism itself treated as established elsewhere.
3. **Whether to spend a run on layer proposal 1**, the GSE190821 epithelial translatome
   under integrin beta6 blockade, which bounds one rival and cannot test the hypothesis.
