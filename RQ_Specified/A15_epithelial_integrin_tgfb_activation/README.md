# A15: does the epithelial input to fibroblast activation run through the integrin or through the ligand

**Status, 27 September 2026: proposed, and nothing is scored here.** This workspace
holds a question, not a result. The observation that motivates it belongs to
[A2](../A2_areg_source_delivery/README.md) and is graded there, if it is graded at all.
The registration itself is the owner's decision; the argument for and against a new
register identifier is in
[reports/REGISTER_DECISION.md](reports/REGISTER_DECISION.md), and until that decision is
recorded this folder is a proposal.

Read in this order:

1. [reports/REGISTER_DECISION.md](reports/REGISTER_DECISION.md), for why this is a
   question rather than a leg of A2 and what a new identifier does not buy.
2. [RATIONALE.md](RATIONALE.md), for the biology, the founding observation with its
   limits, the rivals, and the clinical record that bounds how far any answer travels.
3. [reports/PUBLIC_DATA_SEARCH.md](reports/PUBLIC_DATA_SEARCH.md), for the search that
   found no public dataset able to test the question, and the nearest candidates.
4. [PLAN.md](PLAN.md), for the staged structure and the gate that any candidate dataset
   must clear before anything is computed.
5. [`config/a15_question_contract.json`](config/a15_question_contract.json), for the same
   decisions in machine-readable form.

## The question in one sentence

Where an epithelium and a fibroblast are in contact, is the epithelial contribution to
the fibroblast myofibroblast and collagen programme carried by epithelial integrin
alphaVbeta6 converting latent TGF-beta that the shared matrix already holds, rather than
by epithelial supply of an EGFR ligand?

## What it stands on, and what that is worth

One unplanned observation in one screen. In the GSE307112 alveolosphere screen, which
perturbs mouse alveolar type 2 cells and leaves the human lung fibroblasts unedited,
epithelial **Itgb6** knockout wells sit a median of **-0.938 log2 CPM** below their
unit's depth-matched control wells on the frozen five-gene fibroblast activation score
(COL1A1, ACTA2, POSTN, CTHRC1, TNC), in three of three units where the well clears the
fibroblast eligibility floor, against an endpoint standard deviation of 0.482. It
survives the epithelial-fraction adjustment at -0.992, and a post hoc culture check
puts organoid size within 0.022 log units of controls and fibroblast content at 0.93 of
theirs. The numbers are in
[A2 stage 3](../A2_areg_source_delivery/reports/STAGE3_LEG1_RESULTS.md) and
[A2 stage 5](../A2_areg_source_delivery/reports/STAGE5_SYNTHESIS.md).

Four things make that one observation and not a result:

- **No randomization is available.** Fifty of the 53 plate-3 targets occupy one fixed
  well position in all four units. The ITGB6 wells are all at position C02, so target is
  confounded with plate position and with guide pool.
- **The comparison is asymmetric.** The ligand arm did not move (+0.036, two of four
  units), but the recipient fibroblasts transcribe AREG at 9.846 mean log2 CPM, so that
  arm is a bounded null rather than a measured absence. "Integrin rather than ligand"
  therefore compares one measured decrease against one unbounded null.
- **Preparation independence is unresolved** in this screen
  ([A10 stage 1](../A10_organoid_growth_outcome/reports/STAGE1_IDENTITY_AUDIT.md)), so
  the observation is a within-screen association.
- **The endpoint is a transcript score.** TGF-beta activation was never measured. The
  proposal names a protein event and the evidence behind it is RNA.

## Register readiness: blocked, and honest about why

No public dataset can test the question. The search in
[reports/PUBLIC_DATA_SEARCH.md](reports/PUBLIC_DATA_SEARCH.md) covered GEO, PRIDE, the
Image Data Resource and the BioImage Archive, and found no deposit that pairs an
epithelial integrin perturbation with a measurement of activated TGF-beta. The nearest
candidate, GSE190821, blocks integrin beta6 in vivo with a replicated design but reads
RNA from whole lung and from the epithelium only, so it cannot separate the fibroblast
compartment and has no activation readout. The decisive test is a bench design.

One bounded thing is feasible now and it addresses a rival rather than the hypothesis:
GSE190821's epithelial translatome arm can ask whether blocking integrin beta6 moves the
epithelium's own programme, which is the "changed epithelium, not the integrin" rival.
That is a proposal in [RATIONALE.md](RATIONALE.md), not authorized work.

## Layout

| Path | Contents |
|---|---|
| [reports/REGISTER_DECISION.md](reports/REGISTER_DECISION.md) | the registration argument, the three rejected alternatives, and what the identifier does not buy |
| [RATIONALE.md](RATIONALE.md) | the mechanism with primary sources, the founding observation, seven rivals, the clinical record, the other-layer verdicts |
| [reports/PUBLIC_DATA_SEARCH.md](reports/PUBLIC_DATA_SEARCH.md) | every query run, the six nearest candidates and why each falls short |
| [PLAN.md](PLAN.md) | four stages, the dataset eligibility gate frozen before any candidate is opened, and the stop rules |
| [`config/a15_question_contract.json`](config/a15_question_contract.json) | estimand, required readout, eligibility gate, prohibitions, and the record that nothing is scored |

There are no `scripts/` or `tables/` directories, because nothing has been computed and
nothing may be computed until a dataset clears the stage 1 gate.

## Seven things a later session must not do

1. **Do not write a study note on the source paper.** The screen's source paper is
   roadmap paper 14, DOI 10.1073/pnas.2606113123, and the owner's reading is recorded as
   not started. This workspace uses the deposit and the repository's own records. The
   precedent is DEVELOPMENT decision 21. The sibling deposit GSE307128 carries a series
   summary that states the authors' conclusions; that text was read while searching and
   is deliberately not used, quoted or relied on anywhere here.
2. **Do not treat the A2 Itgb6 observation as this question's result.** A15 opens with no
   result of its own. The observation is A2's, and the wording A2 proposed for it is
   A2's proposal 2 in its stage 5 synthesis.
3. **Do not add a claim row.** The register ends at C168 and grading is the owner's.
4. **Do not read a transcript score as a TGF-beta activation measurement.** The whole
   point of this question is the readout the founding observation did not have.
5. **Do not compute another per-target contrast in GSE307112 and call it a test.** It is
   the same design, at the same fixed positions, with the same confound. A second
   endpoint on the same wells is not a second experiment.
6. **Do not treat wells, guides or organoids as biological replicates**, and do not treat
   an antibody dose arm as a genetic perturbation of the epithelium.
7. **Do not move between the mechanism and the clinic in either direction.** The failure
   of an anti-alphaVbeta6 antibody in an idiopathic pulmonary fibrosis trial does not
   refute the mechanism, and the mechanism does not support a therapeutic reading.
