# Execution scope

`01_audit_geo_metadata.py` retrieves public GEO design metadata. It archives
existing extracts only with explicit `--refresh`; it does not retrieve counts.
`02_preflight.py` checks local catalog integrity and the candidate pair contract.
These are the only entrypoints run on real sources during planning.

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
spatial-protein and proteomic stages currently have specifications, not executed
pipelines or fabricated starter results.

From the repository root:

```powershell
python RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/02_preflight.py
python RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/03_run_paired_counts.py --contrast A1_ATAC_TIGIT
```

The current configuration deliberately reports identity and payload holds.
After review, freeze the source-specific QC and normalization decisions before
setting a contrast to `frozen`; do not edit these fields to bypass missing data.

Software checks: repository unit tests exercise missing identities, overlapping
pools, broken pairs, altered counts and column alignment. The optional
`scripts/tests/check_paired_counts.py` exercises the installed R/edgeR runner
with 500 synthetic features and four pairs, verifying the coefficient direction.
Its outputs stay in ignored `tmp/a1_synthetic_qll/` and are not scientific evidence.
This local check requires the existing portable R installation; it is not an
undeclared dependency of the standard-library GitHub CI job.
