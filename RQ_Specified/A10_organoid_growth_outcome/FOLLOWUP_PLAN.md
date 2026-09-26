# A10 follow-up: design diagnostics before model extensions

26 September 2026. The owner authorized the follow-up order in the
[logical review](../../docs/LOGICAL_RATIONALE_REVIEW.md). This plan first fixes a
descriptive audit; a separate model amendment will be committed after the audit
and before new fits. Previous numerical outputs and specifications are immutable.

## Existing references checked before launch

- [Stage 1 identity audit](reports/STAGE1_IDENTITY_AUDIT.md): join rules, incomplete
  preparation identification and 30 control-like libraries absent from the layout.
- [Stage 3](reports/STAGE3_FIT_REPORT.md) and [stage 4](reports/STAGE4_REVISED_REPORT.md):
  885 analysed wells, plate heterogeneity, 15 grouping labels, target-removal
  sensitivities, and the within-group outcome-centring limit.
- Scripts 01/03/06, their run records and the original/revised contracts:
  use their library identities and extraction hashes; do not re-extract the
  303 MB workbook or repeat completed fits except a bounded instrument check
  needed to verify a new scorer.
- The deposited imaging schema contains a `segmentation_model` field; the design
  contains three guide columns. Audit these explicitly before assuming comparable
  imaging or independent guides. Guide presence does not establish editing success.

## Phase A: fixed diagnostic scope

Use the three original hashed metadata files, tracked joined wells/module scores
and the complete public GEO family SOFT metadata. Record retrieval URL, bytes and
hash; cache the raw SOFT outside Git. No scientific model is fitted in this phase.

1. Summarize library coverage, missing imaging, targets, existing RNA size/fraction
   proxies and species-assignment uncertainty by plate and deposited group.
2. Tabulate target overlap and counts across plates, with all targets and after
   excluding TDTOMATO plus the most represented target (the previous sensitivity).
   Report whether target and plate effects are separable by design, rather than
   inferring independence from labels.
3. Summarize day-7/day-14 mean-area ranges and changes, organoid counts and area
   proportion. Record segmentation-model allocation by plate/day and test the
   arithmetic consistency of independently recorded imaging summaries where possible.
4. Inventory sample characteristics and processing protocols in all deposited
   sample records. Extract explicit preparation/isolation/animal/donor/lot fields
   if present. Generic protocol descriptions are not sample-level identifiers.
   No paper note or author request is authorized by this metadata retrieval.
5. Inspect guide allocation and control labels. Record the available measurement
   and missing verification; do not substitute a transcript decrease for editing
   efficiency or a matched functional control.

## Adaptive gates after the diagnostic

- If joins fail, stop dependent fitting until the mismatch is resolved.
- If imaging definitions or processing differ, resolve their comparability before
  interpreting plate holdout. Unresolved incompatible endpoint scales stop that
  comparison rather than invite an outcome-dependent correction.
- If plate and target allocation lack overlap, a plate holdout can at most be a
  joint plate/target shift stress test. It cannot isolate plate biology or prove
  transfer between independent preparations. Skip target-adjusted comparisons
  whose effects are not identifiable.
- A simpler proliferation-versus-other-growth comparison is eligible only after
  score provenance and the outcome metric are fixed. Prespecify its programme
  groups, covariates, preprocessing, uncertainty scope and effect-size criterion.
- The missing historical BH procedure will not be fabricated. With unresolved
  biological replication, new comparisons remain descriptive; any retained-model
  comparison must state which baseline is actually used.
- Target-level functional follow-up requires verified controls and editing
  efficiency. If unavailable, prune it and document the missing evidence.

Phase A writes only new files under `tables/followup_v1/`; subsequent phases use
distinct filenames and refuse overwrites. Results, code and deviations from this
sequence will be recorded before updating the current handoff.
