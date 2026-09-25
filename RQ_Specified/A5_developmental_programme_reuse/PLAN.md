# A5: prospective test of developmental-gene recruitment in adult injury

**Execution completed after this plan was committed.** Read the
[results and biological interpretation](../A5_A11_shared_component_contract/reports/REVISED_TEST_RESULTS.md).
The prospective rules below and their pre-score commits remain the design record.

25 September 2026. Authorized by the owner's request to proceed with the review
revisions. This plan is committed before downloading or scoring the Strunz count
matrix. See the [biological logic](../A5_A11_shared_component_contract/BIOLOGICAL_LOGIC.md)
and machine-readable `config/strunz_test_contract.json`.

## Question and prior evidence

Do adult Krt8-positive transitional cells recruit the P1 developmental AT1/AT2
signature beyond activated AT2 cells from the same injured mouse? A positive
score is evidence for partial transcriptional recruitment, not embryonic identity
or a shared lineage route. Guo's perinatal signature can include differentiation
and birth stress. Strunz already reported poor overall correspondence with that
developmental signature. This is an explicit rival, not evidence to omit.

The old 94-gene and 51-gene variants were filtered with markers derived from the
Strunz test cohort. They remain descriptive. The shared frozen partition is not
rewritten. An additional source-defined module fixes the independence issue.

## Gene definitions, fixed without Strunz marker selection

- Primary `Guo_AT1_AT2_external`: all 100 published AT1/AT2 genes minus the fixed
  operational set Krt8, Cldn4, Sftpc, Cebpa, Slc4a11; 99 remain.
- `Guo_minus_identity`: subtract Guo's own AT1 and AT2 top-100 lists; 57 remain.
  This challenges ordinary alveolar identity using external selection only.
- `Guo_minus_identity_and_controls`: additionally subtract MSigDB 2024.1 mouse
  Hallmarks P53_PATHWAY, HYPOXIA, INFLAMMATORY_RESPONSE, UNFOLDED_PROTEIN_RESPONSE,
  E2F_TARGETS and G2M_CHECKPOINT. These address perinatal stress and cycling as well
  as generic injury, but do not exhaust every possible rival.
- Guo AT1/AT2 identity axes and the six Hallmarks are descriptive controls. Old
  contract modules, their pairwise components and Strunz identity exclusions are
  reported descriptively, never counted as independent corroboration.

Exact published mouse symbols; no alias repair. Require 70% assay coverage of
each source list. A failed primary coverage gate stops the test.

## Cohort, population and time

GSE141259 high-resolution epithelial experiment. Each sample identifies one mouse;
source methods describe two mice per sampled time point. Restrict to injured mice
from **day 2 through day 21 inclusive**, the window established by the metadata
audit, before any score. Exclude NC controls regardless of their conflicting day
labels. This covers the observed transitional period; later recovery and earlier
onset are outside this estimand.

Primary: author `Krt8+ ADI` versus `AT2 activated`. Both are cells from the same
injured environment. Secondary reference: `AT2` (resting), which may include an
activation difference absent from the primary comparison. Do not switch the
primary reference after seeing results. Author clustering is accepted as the
population definition, not independently validated by this signature test.

Metadata alone gives 24 candidate mice for the activated reference and 26 for
resting AT2. Actual eligibility requires at least **30 cells in each arm after
raw-library depth >= 500 UMI**. At least three paired mice are required. Report all
exclusions. If the primary fails, stop rather than lower the floor or substitute
resting AT2. No animal count from another experiment is added.

## Depth-matched measurement

Use the deposited raw integer UMI matrix and exact barcode alignment. Do not use
metadata `n_counts` as raw library size: it contains non-integer values.

For a gene with k counts in a cell of N total UMI, calculate its probability of
detection in a random sample of 500 UMI without replacement:

`P(detected) = 1 - choose(N-k, 500) / choose(N, 500)`.

Average across present module genes and retained cells, then subtract the
reference from ADI within each mouse. Report percentage points. This integrates
over the random depth-matching draw exactly; technical seeds add no independent
evidence. Zero observed counts give zero probability. Assayed but unexpressed
genes remain in the denominator; genes absent from the matrix do not.

This is a detection score, not the A11 logCPM instrument: do not compare their
effect magnitudes numerically. A 500-UMI budget suits a shallow Drop-seq assay and
is fixed before examining raw depths; failure is recorded without retuning.

## Inference and decisions

Primary estimand: the equal-mouse mean difference for eligible mice in this fixed
window. Use a two-sided one-sample t test on paired differences, alpha 0.05, and
its 95% confidence interval. This estimates a mean and assumes independent mouse
contrasts, with normality for exact small-sample t inference. It is conditional on
the sampled time schedule, not a uniformly sampled population of injury days.

Report direction as positive if the interval is above zero, negative if below
zero, otherwise unresolved. There is no externally justified meaningful-effect
margin for this new instrument; do not label a nonsignificant result equivalent
to absence. Show individual mice, time trends, equal-day mean and leave-one-mouse-
out means to expose heterogeneity and dominance. No selected-day significance tests.

One secondary family, Holm correction across three tests using the same mean/t
procedure: identity-excluded module against activated AT2; identity-and-control-
excluded module against activated AT2; full external module against resting AT2.
An ineligible secondary counts as p=1 for family correction and is reported
ineligible. Each reference uses its own complete mice; display overlap and counts.

Full-module positivity alone supports recruitment of the source signature. A
positive identity-excluded result weakens the identity-only explanation; additional
positivity after control exclusions weakens those specific gene-membership rivals.
Neither proves a developmental mechanism, and exclusion failure may reflect a
smaller/weaker module. The mouse-level results do not establish lineage direction,
successful repair, trajectory ordering, or independent validation of the old
Strunz-filtered modules.

## Launch sequence and existing references

1. Freeze additional modules from hashed Guo/MSigDB sources; preserve old modules.
2. Commit this plan, definitions and A11 revisions before new expression scores.
3. Fetch only the deposited processed count matrix and barcodes; verify integer
   counts, dimensions, barcode identities and source coverage. Raw sequencing is
   unnecessary for this test.
4. Apply fixed depth/cell/unit gates. Score eligible paired mice, run declared
   inference and report all outcomes, including failures or uncertainty.
5. Stop at transcriptional interpretation. A1 already supplies the relevant
   context-dependent perturbation/lineage evidence; further pooled chromatin
   analysis would not provide independent confirmation.

Numerical precedents: ES1 depth-matched detection and biological-unit rules;
Strunz's published developmental comparison and Supplementary Data 3 provenance;
Guo Supplementary Data 2; shared membership and coverage tables; A1 regulatory/fate
report. Source paths, hashes, exclusions and software versions accompany the run.
