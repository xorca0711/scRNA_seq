# A1: Distinguishing transitional epithelial states

Updated 25 September 2026. **Status: revised plan and first analysis batch completed.**
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
- [Studies, usable assays and interpretation limits](STUDY_MAP.md)
- [Public metadata audit and conflicts](metadata/README.md)
- [Current preflight result](reports/PREFLIGHT.md)
- [Figure gallery and remaining figure plan](figures/README.md)
- `scripts/`: metadata retrieval and analysis entrypoints
- `config/`: prospective sample and contrast contracts
- `tables/`, `figures/`, `reports/`: real numerical outputs and execution evidence

The first batch links a measured PATS endpoint to separately analysed IRE1α
perturbation RNA, while preserving their different assays and experimental units.
Histone comparisons need compatible quantification; ATAC/CD44 inferential
contrasts still need source-identity checks. Existing multiome figures remain
historical observations of RNA and accessibility, not the answer to this
expanded question. Full sequencing reprocessing is conditional on the
metadata gates and a measured storage/runtime estimate.
