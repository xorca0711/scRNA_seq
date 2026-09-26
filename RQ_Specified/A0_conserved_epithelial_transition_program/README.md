# A0 — Conserved processes across epithelial transitions

Created: 2026-09-25. Status: **Feasibility audit complete; primary biological analysis blocked by data requirements**.

Read the [current pilot report](reports/PILOT_REPORT.md) for the executed audit,
capture estimates and remaining dependencies. No candidate expression program
has been learned or scored; P1–P4 remain unexecuted.

## Biological question

Do epithelial cells reuse a conserved regulatory process when passing between
differentiation states, despite differences in tissue, starting identity and
destination identity?

The motivating hypothesis concerns a universal process. This pilot tests a
narrower, observable prediction: **a transcriptional program learned from lung
injury and normal lung development also distinguishes independently supported
intermediate states in another epithelium, beyond generic stress, proliferation
and loss of the starting identity.**

## Read in this order

| File | Purpose |
|---|---|
| [PLAN.md](PLAN.md) | Hypotheses, analysis sequence, controls, outputs and expansion criteria |
| [DATASET_CANDIDATES.md](DATASET_CANDIDATES.md) | Three dataset roles, existing resources, eligibility and unresolved selections |
| [Current pilot report](reports/PILOT_REPORT.md) | Extended feasibility findings, figures and explicit limits |
| [Readiness](readiness.json) | Completed work and blocked downstream stages |
| [Reproduction and input contract](REPRODUCING.md) | Commands, required metadata and unblocking conditions |
| [P0 eligibility report](reports/P0_ELIGIBILITY_REPORT.md) | Executed coverage audit, limitations and next selection steps |
| [Dataset audit](tables/dataset_audit.csv) | Candidate-level eligibility decisions |
| [Decisions](decisions.json) | Dated scope and interpretation decisions |
| [Source manifest](source_manifest.json) | Download provenance and hashes |

## Pilot scope

- Three settings: lung repair, normal lung development and differentiation in
  one other epithelium. Prefer one species for the first pass.
- One primary candidate gene program, evaluated against both trajectory
  endpoints and appropriate stressed-cell controls.
- Discovery in the two lung settings; a frozen transfer test in the other
  tissue. Biological samples, rather than individual cells, determine replication.
- A compact report and three figures supporting an **expand, narrow, stop or
  unresolved** decision.

A positive pilot supports a transferable transcriptional association. It does
not establish universality, causal control, actual cell fate or transition rate.

## Relationship to existing work

A0 asks about conservation **across transitions and tissues**. The adjacent
[A1 workspace](../A1_transitional_epithelial_state_distinction/) concerns
transitional epithelial-state distinction. The existing
[ES1 analysis](../../Thesis/epithelial_state_specificity/README.md) supplies
signature provenance, overlap checks and measurement controls; its data and
previously inspected findings are not untouched validation.

The extended audit found nine complete repair triplets in a proposed days 10–15
window, only one complete developmental barcode group (animal identities missing),
and two qualifying mice among four explicitly mapped intestinal controls. The
skin fallback has 34 capture batches from 19 mice without a recovered mapping.
These are coverage findings, not a test of conservation. See the current report
before proceeding to P1.
