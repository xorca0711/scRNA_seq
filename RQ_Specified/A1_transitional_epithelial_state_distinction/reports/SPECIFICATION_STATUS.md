# Specification status — 25 September 2026

**Update:** the owner subsequently requested a stages 3–4 challenge and actual
execution. The revised [lineage audit](../LINEAGE_AUDIT.md) and
[first-batch report](FIRST_BATCH_REPORT.md) now govern current status. The
planning-stage record below describes what had run before that authorization;
it is retained as history, not the current execution state.

The prospective plan, source map, metadata conflicts, candidate sample/contrast
contracts and six figure groups are ready for review. A1 remains one question
in the root register, with its implementation workspace here. Structure PR #66
is separate; this specification is on `codex/a1-state-distinction-plan`.

Completed planning work:

- Audited 29 GEO series; saved source URLs, response hashes and GSM metadata.
- Recorded three ENA project run inventories, a Zenodo IMC file catalog and
  PRIDE project metadata. These are availability records, not downloaded assays.
- Specified direct histone/methylation, accessibility, time/lineage, protein,
  spatial and perturbation work packages without equating assay capabilities.
- Drafted three processed-count contrasts and 24 candidate sample rows. All
  remain held for identity, pool/replicate, assay-QC and processed-payload checks.
- Added a guarded bulk-count pilot and a figure plan. Other modalities remain
  specified stages rather than implemented or executed full pipelines.

Validation: 13 repository unit tests passed; local Markdown/claim validation
passed 2,044 checks; the paired R model passed a synthetic 500-feature/four-pair
direction check. No biological counts were fitted and no new scientific figures
or findings were produced. The synthetic check is a software test, not a power
analysis or validation of real-data assumptions.

The next scientific step proposed at that planning checkpoint was to resolve source-unit/file mappings and
inspect the small deposited processed tables. Freeze the first eligible
contrast and measured resource budget before inference. Direct histone profiles
can proceed as descriptive work when their controls and normalization are
verified. Raw-read processing, velocity reconstruction and IMC re-segmentation
need separate measured budgets; no defensible completion time exists yet.
