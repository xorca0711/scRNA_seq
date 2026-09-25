# Execution scope

## Latest robustness batch

Scripts 22–25 implement the frozen [robustness scope](../config/robustness_batch.json):

| Script | Purpose |
|---|---|
| `22_prepare_robustness.py` | Pinned author-definition eligibility; records the two notebook size holds without downloading them |
| `23_analyze_hpcs_robustness.py` | Library/design rank, source omissions, raw-category compositions and label-invariant partition agreement |
| `24_analyze_histone_tss.py` | All eligible native transcript starts, original-TSS baseline, both windows and matched H3 controls |
| `25_verify_robustness_batch.py` | Independent source/partition checks, baseline reproduction and alternative-window base-resolution checks |
| `a1_robustness.py` | Tested source-weighting, partition-agreement and parent/strand-aware transcript helpers |

The annotation audit is qualified by its [eligibility scope](../config/robustness_annotation_scope.md).
All completed entrypoints refuse to overwrite their output paths. Use the existing
scientific launcher; script 25 additionally uses the installed scikit-learn.
Read the [report](../reports/ROBUSTNESS_REPORT.md) before planning another run.

## Earlier execution scope

`01_audit_geo_metadata.py` retrieves public GEO design metadata. It archives
existing extracts only with explicit `--refresh`; it does not retrieve counts.
`02_preflight.py` checks local catalog integrity and the candidate pair contract.
Those were the only entrypoints run on real sources during the initial planning
stage. The first numerical batch has now run as described below.

`03_run_paired_counts.py` defaults to reporting holds. Its optional `--execute`
mode requires an authorized, frozen contract, verified independent sample/pool
identities, a hashed raw integer-count TSV and an Rscript path. It prepares
sample-aligned inputs and calls `04_fit_paired_counts.R`. The R runner uses a
paired [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
quasi-likelihood model, saves the feature filter, design, normalization, effects,
sample PCA and session information. Counts stay in ignored `cache/`; each run
gets a new output directory and hashes, never overwriting earlier records.

The count adapter must have one unique feature-ID column and one column per
verified, already aggregated library. It must retain provenance back to the
deposited files. Technical replicates cannot masquerade as additional animals.
Source-specific adapters, peak annotation/QC, normalization diagnostics,
leave-one-pair-out sensitivity and cross-contrast FDR reporting follow payload
inspection. The pilot is not the complete A1 workflow. Histone, trajectory,
spatial-protein and proteomic stages retain separate eligibility gates. Direct
histone/domain quantification and a descriptive descendant reconstruction have
now run in the bounded second batch described below.

From the repository root:

```powershell
python RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/02_preflight.py
python RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/03_run_paired_counts.py --contrast A1_ATAC_TIGIT
```

The current configuration deliberately reports identity and payload holds.
After review, freeze the source-specific QC and normalization decisions before
setting a contrast to `frozen`; do not edit these fields to bypass missing data.

## Verified second batch

| Script | Purpose |
|---|---|
| `11_run_ire1_stability.py`, `12_fit_ire1_stability.R` | Completed omission/within-batch fits; reuse outputs, do not refit |
| `13_fetch_direct_mark_inputs.py`, `14_fetch_followup_sources.py` | Bounded original source acquisition and inventories |
| `15_quantify_direct_marks.py` | Completed exact native-assembly signal/domain quantification; refuses overwrite |
| `16_summarize_second_batch.py` | Cached-data summary/render with required fresh `--run-id`; refuses existing outputs |
| `17_verify_second_batch.py` | Input/output hashes, fit/BH/universe checks, base-resolution controls and independent domain unions |
| `18_fetch_identity_evidence.py` | Bounded public source/ENA/notebook retrieval; downloaded code is never executed |
| `19_recover_hpcs_metadata.py` | Remote HDF5 `/obs` only; 64-MiB/15-minute hard ceilings and exact range validation |
| `20_audit_delivery_evidence.py` | Fresh identity tables, native-size and preserved/corrected presentation checks |
| `21_reconstruct_hpcs_source_composition.py` | Frozen descriptive source counts and figure, author-output checks, no biological inference |

Scientific entrypoints use `analysis/scripts/run_with_environment.py` with
`.venv-x64/Lib/site-packages` and the working Python runtime in the handoff.
New figure example: `16_summarize_second_batch.py --run-id another_presentation`.
The existing `second_batch_verified` run must not be reused as an output name.
Reports, source hashes and limitations are in
[SECOND_BATCH_REPORT.md](../reports/SECOND_BATCH_REPORT.md).

Software checks: repository unit tests exercise missing identities, overlapping
pools, broken pairs, altered counts and column alignment. The optional
`scripts/tests/check_paired_counts.py` exercises the installed R/edgeR runner
with 500 synthetic features and four pairs, verifying the coefficient direction.
Its outputs stay in ignored `tmp/a1_synthetic_qll/` and are not scientific evidence.
This local check requires the existing portable R installation; it is not an
undeclared dependency of the standard-library GitHub CI job.

## Adaptive evidence closure (scripts 26–31)

26 fetches two pinned annotation notebooks and 24 exact SRA records; 27 fetches
small primary texts and the 20-experiment PATS catalog. 28 verifies identities,
HPCS mapping and confidence abstention. 29 orchestrates the new CD44 analysis;
30 fits paired genotype-specific contrasts and their direct interaction, with
conditional pair omissions. 31 independently verifies results and renders two
figures. All completed entrypoints refuse output overwrites. Source inventories
record hashes; author notebook code is never executed.

See [the reference map](../reports/ANALYSIS_REFERENCE_MAP.md),
[frozen contract](../config/closure_analysis_contract.json) and
[closure report](../reports/EVIDENCE_CLOSURE_REPORT.md). Script 29 needs the same
`--rscript` path as script 07 and small NCBI gene-ID lookups. Its compressed count
adapter is ignored; deposited inputs and hashes remain authoritative.

## First real analysis batch

| Script | Purpose |
|---|---|
| `05_fetch_processed_inputs.py` | Bounded processed-count/BED downloads, source URLs and SHA-256 inventory |
| `06_lineage_source_reanalysis.py` | Kobayashi ED4 workbook parsing, zero-denominator handling, nested field/mouse summaries and figure |
| `07_ire1_epithelial_analysis.py` | Exact GSE190821 mouse/compartment/treatment crosswalk, count checks, frozen-contract orchestration and figures |
| `08_fit_ire1.R` | Unpaired batch + sex + treatment edgeR model, sensitivity, PCA and eligible estimated-correlation CAMERA tests |
| `09_descriptive_input_audit.py` | TIGIT/CD44 source-alias PCA; H3K4me3 caller/geometry audit without biological tests |
| `10_render_audited_figures.py` | Re-render existing results and verify unchanged scientific-table hashes |

Python analysis scripts need NumPy, pandas and matplotlib in the scientific
environment. R needs edgeR, limma and statmod. On this workstation, the existing
`analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages`
launcher supplies Python packages, and the portable Rscript is at
`analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe`.
Pass that path with `--rscript` to script 07. Public inputs and Ensembl annotation
snapshots are identified in the inventories/contracts; downloaded assay payloads
stay ignored under `cache/`, selected count adapters under `processed/`.

The completed numerical entrypoints refuse to overwrite existing run records.
Archive a batch's reports, tables and figures before a new numerical run.
Script 10 updates presentation only. Exact first-run source bytes are preserved
in [the execution archive](../reports/execution_sources/2026-09-25/manifest.json).
Do not execute archived source files as a workaround for the overwrite guard.
The [report](../reports/FIRST_BATCH_REPORT.md) links all run records and limitations.
