# Implementation after the repository review

The owner authorized corrections, claim reclassification and improvement
orders 1–4 on 22 September 2026, then authorized relevant additional analyses
and figures for portfolio use. Decisions below are assistant assessments under
that authorization, not a claim that the owner individually reviewed each row.

## Scope and work ownership

- Main task: scientific wording, claim decisions, structured evidence contract,
  CI, reproducibility documentation, portfolio entry point and final validation.
- Ligand task: corrected epithelial sender/target populations, resource rankings,
  annotation sensitivity and a supporting figure.
- Statistics task: sample-level enrichment sensitivity, confounding, leave-one-out
  stability and reference CAMERA verification if a working R runtime is available.
- Epithelial task: unify A1/A5, freeze traceable state definitions, run feasible
  descriptive contrasts and explicitly assess independent confirmation.

Historical trials and raw inputs are preserved. New corrected analyses have
separate paths and record which data were already inspected. Changes to claim
status will be linked to evidence in decision 32 of DEVELOPMENT.md.

## Completion criteria

1. Current claims and research framing match the quantity actually measured.
2. Corrected measurements have saved outputs, provenance and targeted checks.
3. Retained enrichment directions are reassessed with sample-level uncertainty;
   unidentifiable contrasts remain limited rather than being solved by adjustment.
4. A1/A5 become one specificity project with honest module provenance and a
   confirmation-cohort eligibility decision; missing independent confirmation is
   reported explicitly.
5. Once these deliverables are assessed, read the supplied Notion target map and
   select subsequent datasets by scientific design and portfolio relevance.

No paid compute, remote publication, or outreach is part of this implementation.

## Completed corrections and decisions

| Order | Delivered | Decision and remaining limit |
|---|---|---|
| 1: claims and framing | Revised current claims, decision 32, generated claim index, updated README/RQs/negative results | Computational measurement and biological fate are separate questions. Historical decisions remain traceable. |
| 2: source and ligand measurement | Raw-count C37 depth sensitivity; corrected C12/C14 senders, actual fibroblast floors and five-resource rankings; two new figures | C37 is descriptive and measurement-dependent. C114 resource absence is refuted. RNA rankings do not establish signalling. |
| 3: programme inference | Official CAMERA, voom/TMM, design sensitivities, 139 leave-one-donor-out refits, subtype checks and figure | C158 not established; C161 exploratory/method-sensitive; C163 still not established. Composition and design confounding remain unresolved. |
| 4: epithelial specificity | Unified A1/A5; full published ADI/AT1/AT2 lists, separate short panels, matched-depth pilot, genotype/context contrasts and external eligibility; new figure | Descriptive pilot completed. Only one external animal passes both floors; replicated confirmation is not established. Full DATP/PATS and developmental-maturation signatures remain unfilled. |
| 5: purpose-aware next stage | Read the supplied Notion target map and linked RQ context; fetched four GEO metadata records; versioned candidate manifest and bounded pilot contract | Prioritize outcome-linked GSE307112 after sample-identity audit. GSE242510 is a secondary descriptive candidate. Biological replication is not inferred from library counts. |

### Main numerical consequences

- **C37:** ten donor pairs remain. At the primary 1,000-UMI budget, six have
  higher epithelial detection; paired p = 0.130859. All 16,064 selected cells
  remain at this budget. The 2,000-UMI sensitivity excludes 1,188 cells and gives
  p = 0.019531; it is reported alongside, not substituted for, the primary result.
- **C80–C82/C111–C114:** the corrected CellChatDB pass retains 22 donors and
  312 pairs. AREG's donor-median rank is 10.5. Canonical exact-EGFR ranking is
  distinct from all resource entries; consensus AREG–EGFR_ERBB2 coverage is 7/22.
- **C161:** of 254 frozen legacy candidates, 244 are eligible in validation.
  The primary estimated-correlation CAMERA yields zero replications; fixed
  correlation 0.01 yields 30, of which 28 overlap the original 62. This supports
  method sensitivity, not either robust validation or biological absence.
- **C165–C168:** full ADI holdout contrasts are about 1.04/1.21 detection points
  in P9 controls and 7.61/6.53 in injured controls across technical seeds.
  These seeds are not biological replicates. Published ADI/AT1 holdout lists
  overlap by 119 genes. One of 25 external animals passes both state floors.

Exact before/after rows, including new C165–C168, are in
[claim_decisions.json](claim_decisions.json). The original
[audit](../../audits/2026-09-22/REPOSITORY_REVIEW.md) describes the pre-correction
state and is intentionally not rewritten as a post-correction report.

### Additional research questions

[A6–A10](../../../RESEARCH_QUESTIONS.md) now distinguish composition from
within-state change, genotype changes in the reference population, shared
transition programmes from mature AT1 identity, receptor representation from
receiver specificity, and RNA programmes from measured organoid outcomes.
Each states a test and a decision boundary. They are proposals generated during
the audit; none is promoted to a biological conclusion.

### Structural work

The claim contract generates a manifest and summary from the narrative register.
Selected numeric bindings check explicit artifact fields and filters; coverage
is visible per row, and unnormalized sample metadata is marked as such.
CI checks both source trees and runs contract tests. New run records archive
earlier bytes before replacement, hash inputs and record code identity.
Historical runs are not retroactively called preregistered or fully hashed.

The [portfolio summary](../../PORTFOLIO_SUMMARY.md) provides four figures and a
short account of the contribution. The [next dataset gate](../../NEXT_DATASET_GATE.md)
and [candidate manifest](../../../analysis/claims/dataset_candidates.json) capture
scientific planning without copying private PI outreach rankings into the repo.

### Verification and cost control

The three independent tasks ran bounded, nonoverlapping corrections: ligand
checks (5 tests), statistical checks (9), and epithelial checks (4). The latter
also reconstructed 1,120 primary scores and checked module source hashes.
Root contract tests cover archive preservation, missing/duplicate claims,
family assignment, displaced rows and selected narrative/numeric drift.
Final integrated verification is recorded in [validation.json](validation.json).

Full analysis reruns were not repeated after unchanged checks passed. A portable
local R runtime restored reference-method testing; no paid compute was used.
An R source edit during the long influence loop caused a trailing parse error
after all outputs were written. The statistics README and run record disclose
the recovery, final parsing and complete 139-refit output audit. There is no
claim of an uninterrupted run. Reliable aggregate token or monetary accounting
was unavailable, so no usage totals are invented.

### What remains genuinely open

Independent biological confirmation, outcome-linked validation, harmonized
macrophage states, age/time identifiability and complete developmental-state
references cannot be supplied by relabelling existing cells. The next screen
has useful imaging endpoints, but biological preparation IDs and plate/guide
structure must be resolved. Quantitative claims outside the selected bindings
are not automatically recomputed by CI. These are recorded limits, not silently
closed tasks.
