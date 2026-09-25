# A1: Distinguishing transitional epithelial states

Updated 25 September 2026. **Status: adaptive closure analyses completed.**
Start with [the closure report](reports/EVIDENCE_CLOSURE_REPORT.md) and
[the analysis-reference map](reports/ANALYSIS_REFERENCE_MAP.md). SRA original
filenames unlock all eight CD44 mice: paired contrasts and a direct genotype
interaction now run. HPCS biological labels are verified; stringent K12 changes
are confidence abstentions. The closure ledger separates finished work,
pruned analyses and precise external-input requirements. See also the
[robustness report](reports/ROBUSTNESS_REPORT.md),
[second-batch report](reports/SECOND_BATCH_REPORT.md), [figure gallery](figures/README.md)
and [job status](JOBS.md). Original numerical runs and presentation evidence are
preserved. PRs #70/#71 are merged; delivery branch: `codex/a1-evidence-closure`.

The preceding second batch adds checked IRE1α sensitivities, 24 native-assembly histone tracks
at 23 loci, a one-donor methylation-domain reference, and a newly recovered
HPCS descendant source-composition table (5,333 cells, 22 source labels).
H3 and window sensitivity qualify the histone patterns. HPCS source labels
are not verified independent animals, and current mScarlet status is unavailable.
The results support specific follow-ups, not a universal state taxonomy.

First-batch baseline:
Measured-lineage source reconstruction, a ten-mouse IRE1α RNA contrast, two
descriptive sample PCAs and the histone-input audit have run. Direct epigenetic
state distinction remains unresolved; see the [batch report](reports/FIRST_BATCH_REPORT.md).
The scientific question is registered once, under [A1](../../RESEARCH_QUESTIONS.md).

The aim is to determine whether overlapping epithelial RNA programmes mark
the same regulatory state, distinct states, or successive stages with different
developmental origins and behaviours. Chromatin accessibility, histone marks
and DNA methylation are the central molecular measurements. Lineage tracing,
time-resolved trajectories, protein measurements, spatial morphology and
perturbation responses supply additional, separately assessed evidence.

DATP, PATS, Krt8 ADI, aberrant basaloid/ABI and HPCS retain their original
study definitions. They are not assumed to be synonyms or predetermined
classes. HPCS is a cancer-context state; a shared injury signature does not
make injured epithelium malignant.

- [Analysis plan and execution gates](PLAN.md)
- [Lineage evidence and stages 3–4 revision](LINEAGE_AUDIT.md)
- [Completed batch, findings and remaining work](reports/FIRST_BATCH_REPORT.md)
- [Verified regulatory follow-up and descendant reconstruction](reports/SECOND_BATCH_REPORT.md)
- [Further analyses and the gaps they can address](reports/REMAINING_ANALYSIS_OPTIONS.md)
- [Source/annotation robustness and alternate-TSS results](reports/ROBUSTNESS_REPORT.md)
- [Studies, usable assays and interpretation limits](STUDY_MAP.md)
- [Public metadata audit and conflicts](metadata/README.md)
- [Current preflight result](reports/PREFLIGHT.md)
- [Figure gallery and remaining figure plan](figures/README.md)
- `scripts/`: metadata retrieval and analysis entrypoints
- `config/`: prospective sample and contrast contracts
- `tables/`, `figures/`, `reports/`: real numerical outputs and execution evidence

The first batch links a measured PATS endpoint to separately analysed IRE1α
perturbation RNA, while preserving their different assays and experimental units.
PATS histone comparisons still need deposited-track scaling; ATAC/CD44 inferential
contrasts still need source-identity checks. The induced-cell histone comparison
is now quantified descriptively. Existing multiome figures remain
historical observations of RNA and accessibility, not the answer to this
expanded question. Full sequencing reprocessing is conditional on the
metadata gates and a measured storage/runtime estimate.
