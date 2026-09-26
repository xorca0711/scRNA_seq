# A15 analysis plan: what has to exist before anything is computed

Declared 27 September 2026. **Stage 0 only is authorized at declaration, and stage 0
computes nothing.** Stage 1 is a gate rather than an analysis, and it currently fails, so
stages 2 to 4 are not authorized and no endpoint may be computed in any dataset. The
biology is in [RATIONALE.md](RATIONALE.md), the machine-readable decisions in
[`config/a15_question_contract.json`](config/a15_question_contract.json), the registration
argument in [reports/REGISTER_DECISION.md](reports/REGISTER_DECISION.md), and the question
in the [A15 card](../../RESEARCH_QUESTIONS.md#a15).

This plan inherits the prohibitions that the [dataset gate](../../docs/NEXT_DATASET_GATE.md)
and [A10 stage 1](../A10_organoid_growth_outcome/reports/STAGE1_IDENTITY_AUDIT.md) placed
on GSE307112, and the measurement contracts
[MC1](../../docs/RQ_MEASUREMENT_CONTRACTS.md#mc1) and
[MC2 to MC4](../../docs/RQ_MEASUREMENT_CONTRACTS.md#mc2).

## Structure

Four stages. Each can stop the work, and stopping is a result.

| Stage | Work | Reads | Can it stop the work |
|---|---|---|---|
| 0 | Register the question, state the founding observation with its limits, search for data | nothing computed | Yes. The owner may refuse the registration |
| 1 | Eligibility gate on a candidate dataset | that dataset's metadata only | Yes, and it stops the work today |
| 2 | Freeze the test | nothing new | No, but it fixes what stage 3 may do |
| 3 | Execute the frozen test | the candidate dataset | Yes, on any failed input check |
| 4 | Report, and propose claim wording without grading it | own outputs | No |

## Stage 0: registration, which is what this folder contains

The card is proposed, the observation is restated with its limits, the rivals are named,
the clinical record is stated where it bounds the answer, and the public-data search is
recorded with its queries and counts. Nothing is computed and no claim row is added. The
owner retains or rejects the registration before stage 1 is attempted on any dataset.
Rejection is recorded in DEVELOPMENT.md as a rejection, not as a revision.

## Stage 1: the eligibility gate, frozen before any candidate is opened

The gate exists so that a later session cannot relax the readout requirement in order to
get a run. **A candidate dataset clears stage 1 only if all six conditions hold**, each
established from the candidate's own files rather than from its abstract or its series
summary.

1. **A perturbation of the activator on the source cell.** Epithelial ITGB6 or ITGAV is
   removed, blocked or induced, with a recorded validation that the perturbation took
   effect. An antibody dose arm counts as a perturbation of the protein and is labelled as
   such; it is not a genetic perturbation of the epithelium and may not be described as
   one.
2. **More than one independent unit per perturbation, and positions that vary.** At least
   three units per arm, where a unit is an animal, a donor or an independent preparation.
   A well is not a unit. If the design places a target at one fixed position in every
   replicate, the candidate fails this condition, which is the condition GSE307112 fails.
3. **A measure of activated TGF-beta, distinguished from the latent pool.** A reporter
   bioassay, an immunoassay specific to the active dimer, or a receptor-proximal
   signalling measurement in the recipient. **A transcript score of TGF-beta target genes
   does not satisfy this condition**, and no pre-declared alternative to it exists.
4. **A recipient readout separated from the source.** The fibroblast compartment is
   measured on its own, by species assignment, sorting, label or physical separation. A
   whole-tissue homogenate fails this condition.
5. **A ligand arm that can be bounded.** Either the recipient's own EGFR ligand is also
   removed, or the recipient is shown not to supply it. Without this, the partition that
   the question asks about cannot be formed, and the work would repeat A2's asymmetry.
6. **Independent units that are established, not assumed.** The deposit must say what its
   units are. If it cannot, the candidate may still be described, and anything computed
   from it is a within-dataset association with no claim of independent replication, which
   is the restriction A10 stage 1 already placed on this repository.

**Stop rules.** If any of conditions 1 to 5 fails, the candidate does not carry the test
and the work stops for that candidate; the failure is reported with the condition it
failed. If only condition 6 fails, stages 2 and 3 may proceed and every result is labelled
descriptive. **As of 27 September 2026 no public candidate clears conditions 3, 4 and 5
simultaneously**, which is the recorded state of this question, and
[reports/PUBLIC_DATA_SEARCH.md](reports/PUBLIC_DATA_SEARCH.md) is the evidence.

## Stage 2: what a freeze would have to fix, declared now so it cannot drift

Stage 2 cannot be written in full before a candidate exists, because the unit and the
reference set depend on the design. Six invariants are fixed now, and a freeze that
contradicts any of them is not a freeze of this question.

- **Unit.** The animal, the donor or the independent preparation. Never the well, the guide
  or the organoid.
- **Primary endpoint.** Activated TGF-beta at the recipient, as condition 3 defines it. The
  fibroblast activation programme is a **secondary** endpoint, not a substitute, and the
  five-gene score may be reused only as it was frozen elsewhere, without re-weighting.
- **The mediation arm is required, not optional.** A frozen test must include the
  recipient-side TGF-beta receptor blockade of prediction 3, or state that it does not test
  prediction 3 and therefore cannot establish the route.
- **Direction, declared in advance.** Lower activated TGF-beta and a lower fibroblast
  programme when the epithelial activator is removed. A shift in the opposite direction is
  reported and does not support the hypothesis.
- **Precise absence.** Unavailable unless the perturbation is shown to be complete. Every
  perturbation recorded in this axis so far is partial, and no bound will be reported as
  absence.
- **Multiplicity family.** Two endpoints and one mediation arm. A block of pathway sets is
  not used.

## Stage 3: execution, not authorized

No dataset may be opened for an endpoint until stage 1 passes on that dataset and stage 2
is committed for it. The A2 precedent is binding here: A2's first freeze declared an
inference its design could not support, and it had to be withdrawn before anything was
scored. The gate above is what that withdrawal bought.

## Stage 4: report and register

The report states what is established, what is not, and which rival each arm excluded.
Claim wording is proposed for the register; grading is the owner's decision and no row is
added by the assistant. If the question is still blocked, that is the report, and it
belongs in NEGATIVE_RESULTS.md through the usual generation path rather than as a silent
absence.

## A separately labelled side-branch, if the owner authorizes it

Layer proposal 1 in [RATIONALE.md](RATIONALE.md) is the only thing runnable today: asking
whether integrin beta6 blockade moves the epithelium's own programme in GSE190821. It is
**not part of the staged test above** and must never be reported as one. It bounds rival 2
only, its unit is the mouse, its arms are four treated against seven antibody-control mice,
and it would need its own pre-declared endpoint and its own stop rule before the counts are
opened. It cannot measure activated TGF-beta and it cannot separate the fibroblast
compartment.

## What this plan refuses to do

- **No transcript score is called a TGF-beta activation measurement**, under any framing,
  and no pre-declared alternative to condition 3 exists.
- **No second endpoint on GSE307112's plate-3 wells counts as a second experiment.** Same
  design, same fixed positions, same guide pools.
- **No claim row is added, and no existing row is re-graded.** The register ends at C168 and
  the founding observation's wording belongs to A2.
- **No wells, guides or organoids as biological replicates.**
- **No reading of an epithelial perturbation as fibroblast causation**, which is A10's
  prohibition inherited in A10's words.
- **No study note on the screen's source paper**, roadmap paper 14, whose reading by the
  owner is recorded as not started, and no use of the sibling deposit's series summary.
- **No therapeutic reading in either direction.** The failure of an anti-alphaVbeta6
  antibody in idiopathic pulmonary fibrosis is not a refutation of the mechanism, and a
  positive result here is not support for an intervention.
- **No re-asking of the established mechanism.** The question is the partition of the
  epithelial output; that integrin alphaVbeta6 activates latent TGF-beta is established
  elsewhere and is not this repository's to rediscover.

## Order of work

1. This plan, its rationale, the contract, the registration argument and the search report
   are committed before anything else happens.
2. The owner retains or rejects the registration, and separately decides on the
   side-branch.
3. If a candidate dataset appears, stage 1 runs on that candidate's metadata alone and its
   gate result is committed, whatever it shows.
4. Only then is stage 2 committed for that candidate and stage 3 run.
