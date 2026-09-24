# Reproducing the review-context analyses

This folder contains new analyses motivated by the review, with public
primary-study data and source-defined gene lists. It is not a reproduction
of an experiment performed by the review authors. The analysis plans,
specifications and run records retain the criteria used before each stage.

## Runtime and invocation

Run from the repository root. Scientific Python entrypoints use the existing
environment through `analysis/scripts/run_with_environment.py`, for example:

```powershell
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u6_audit_pathway_coverage.py
```

The R scripts load the repository's local R library and call official edgeR
and limma implementations. Their session files record exact versions. The
completed run used R 4.6.1, edgeR 4.10.5 and limma 3.68.5. Python versions and
source/model identities are recorded in the stage specifications. Mapping
and native LR computations run with bounded CPU threads; large matrices are
streamed or aggregated without holding the full cohort in memory.

Acquisition and computational entrypoints generally refuse to overwrite an
existing run record. This protects provenance. Reproduce in a fresh checkout
with the specified ignored cache inputs, or use an explicitly implemented
resume mode; do not delete earlier evidence merely to make a script run.
Report/figure entrypoints can regenerate their derived presentation files.

## Stage dependencies

| Stage | Executable entrypoints | Main evidence directory |
|---|---|---|
| Public source/design audit | `u0_geo_design_audit.py`, U2 source-access entrypoints | `trials/u0_geo_design_audit/`, `trials/u2_remaining_context_access/` |
| Mouse annotation inputs | `u3_lineage_inputs.py`, `u4_freeze_resources.py` | `trials/u3_lineage_annotation/`, `trials/u4_resources/` |
| Mouse niches | `u4_run_mouse_pseudobulk.py`, `u4_mouse_liana.py`, `u4_mouse_compatibility.py` | `trials/u4_mouse_niche/` |
| IPF pathways | `u5_ipf_pathways.R` and its launcher | `trials/u5_ipf_pathways/` |
| IPF full-source extraction and native LR | `u5_full_source_panels.py`, `u5_liana_robustness.py` | `trials/u5_full_source_panels/`, `trials/u5_liana_robustness/` |
| IPF RNA compatibility | `u5_ipf_compatibility.py`, `u5_compatibility_normalize.R` | `trials/u5_ipf_compatibility/` |
| IPF ligand targets | `u2_acquire_ligand_prior.py`, `u5_run_ligand_targets.py`, `u5_filter_ligand_targets.py` | `trials/u5_ligand_targets/` |
| Human assay and reference mapping | `u5_acquire_human_full.py`, `u5_human_pilot_input.py`, `u5_human_reference_mapping.py`, `u5_human_full_run.py` | `trials/u5_human_pilot/`, `trials/u5_human_full/` |
| Human paired stages | `u5_run_human_paired.py` freezes and runs aggregation, pathways, compatibility and target-refit scripts in order | `trials/u5_human_niche/`, `trials/u5_human_ligand_targets/` |
| Human LR and all-QC sources | `u5_human_liana.py`, `u5_human_source_profiles.py` | `trials/u5_human_niche/`, `trials/u5_human_sources/` |
| Human interpretation release | `u5_review_human_full.py`, documented marker review, `u5_filter_human_ligand_targets.py`, corresponding report scripts | Same human directories |
| Spatial context | `u5_acquire_spatial_context.py`, `u5_spatial_measurements.py`, `u5_report_spatial_context.py` | `trials/u5_spatial_context/` |
| Source-program specificity | `u6_score_specificity.py`, `u6_freeze_isr_and_human_modules.py`, `u6_score_isr.py` and report scripts | `trials/u6_specificity/`, `trials/u6_isr_extension/` |
| IPF and human program extension | `u6_ipf_epithelial_pseudobulks.py`, `u6_ipf_program_scores.R`, `u6_report_ipf_specificity.py`; human scores come from the paired pathway stage | `trials/u6_ipf_specificity/`, `trials/u6_human_specificity/` |
| Completion | `u6_audit_pathway_coverage.py`, `u6_report_context_coverage.py`, `u6_summarize_evidence.py`, `u6_validate_release.py` | `trials/u6_completion/` |
| Shared RQ figures | Root scripts 19–21; [methods and commands](../../analysis/figures/rq/il1b_context/REPORT.md) | `analysis/figures/rq/il1b_context/` from the repository root |

This dependency table is an entrypoint index, not a claim that scripts with
required cohort arguments can all be invoked without arguments. Each script's
CLI and adjacent specification define its inputs and supported resume mode.

## Sources, versions and data boundaries

GEO acquisition records retain URLs, hashes, sizes and sample crosswalks.
The human `.raw_counts.mtx.txt.gz` inputs are dense gene-row text, despite
their names. Their bounded count-store parser is separately validated.
Full gene totals, gene uniqueness, nonnegative integer counts and pseudobulk
sum parity are checked before inference.

The human atlas uses the existing pinned HLCA reference and a frozen query
adaptation, with the same gene panel and confidence thresholds throughout.
The 13 pathway sets are frozen from MSigDB 2024.1.Hs resources. Mouse mapping
uses strict HCOP one-to-one orthologs; missing genes retain their source-list
denominators. NicheNet v2 priors are checksum-verified from Zenodo 7074291.
Author HPCS and Han ISR definitions have source URLs, pinned code or document
hashes and explicit overlap-removal lists in their module specifications.

Large public matrices, model objects, per-cell files, full private notes and
copyrighted papers remain in ignored caches. Compact biological-unit tables,
figures, scripts, specifications and validation records form the reviewable
release. Historical raw inputs and results in other paper folders are reused
read-only. The repository claim contract is not automatically expanded by
this exploratory review-context analysis.

## Corrections and historical records

- The initial mouse-resource description reversed eligible treatment-arm
  counts in prose. The current denominator is three IgG and two anti-IL1B
  animals at 100 alveolar cells. Numerical fits used the correct labels;
  the original resource file and its hash are retained in `.history`.
- The full human input run resumed after two libraries to fix retention of
  scvi-tools managers in memory. The previous source and run record remain
  archived. Counts, model weights, labels and criteria were not changed.
- Native LR empty-edge handling was corrected after an independently
  verified empty eligible set; absent edges were not assigned zero scores.
  Archived entrypoints identify the earlier executions.
- Earlier W labels were execution bookkeeping. The final U-stage names
  follow the analysis plan; relocation records preserve the old references.

The final manifest hashes the released files. Figure reviews bind visual
checks to PNG hashes, so an edited figure requires a new review. Local links,
frozen human script hashes, stage checks and absence of raw/private cache
files are validated before release.
