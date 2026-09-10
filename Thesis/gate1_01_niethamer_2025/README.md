# Gate 1, paper 1: Niethamer et al. 2025 (GSE262927)

This paper's study note and its analysis trial already exist in the
repository. This folder only points to them so that the roadmap order in
[`../README.md`](../README.md) is visible in one place.

**Citation.** Niethamer TK, Planer JD, Morley MP, et al. *Longitudinal
single-cell profiles of lung regeneration after viral infection reveal
persistent injury-associated cell states.* Cell Stem Cell 2025;32(2):302-321.e6.
DOI 10.1016/j.stem.2024.12.002, PMID 39818203. Data: GEO GSE262927.

| What | Where |
|---|---|
| Annotated reference to the published method (parameters, marker tables, open decisions) | [`../../docs/scRNAseq_workflow_Niethamer2025.md`](../../docs/scRNAseq_workflow_Niethamer2025.md) |
| The published workflow, ordering constraints, statistical unit | [`../../WORKFLOW.md`](../../WORKFLOW.md) |
| The executed reanalysis (generated report) | [`../../analysis/GSE262927/README.md`](../../analysis/GSE262927/README.md) |
| Focused analyses: AT2 to Krt8-high transitional to AT1 trajectory; injury-associated capillary state | [`../../analysis/GSE262927/regeneration_focus/`](../../analysis/GSE262927/regeneration_focus/) |
| Lineage-tracing cohort (8 non-atlas samples) | [`../../analysis/GSE262927/lineage_tracing_cohort/`](../../analysis/GSE262927/lineage_tracing_cohort/) |
| Every decision, before and after reading the paper | [`../../docs/ANALYSIS_RATIONALE.md`](../../docs/ANALYSIS_RATIONALE.md) |
| What actually ran, with the divergence table | [`../../docs/PIPELINE_AS_RUN.md`](../../docs/PIPELINE_AS_RUN.md) |
| Results and negative results | [`../../FINDINGS.md`](../../FINDINGS.md) |
| Full citation block | [`../../REFERENCES.md`](../../REFERENCES.md) |
| The analysis in steps: initial run (Stage 0) versus PI-matched follow-ups (Stage 1) and proposals (Stage 2) | [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md) |

**Status.** Study note: done (in the docs above). Analysis trial: done; the
validated claims are listed in the root README claims table. Four PI-matched
follow-ups (per-day phase view, myeloid compartment, alveolar macrophage
origin by trace window, batch sensitivity on infection round) were run on
2026-09-10 with frozen rules; they are Descriptive only and await owner
review; [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md) lays them out
step by step against the initial run. The five reading-workflow questions of
the roadmap (question, evidence type, reusable variables, one limitation, one
bridge) are answered in prose across the rationale document and are not yet
collected in the roadmap format; that is the remaining item for this folder.
