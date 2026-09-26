# A0 reproduction and required inputs

Current result: feasibility audit completed; P1–P4 blocked. The authoritative
status is [readiness.json](readiness.json), with interpretation in the
[pilot report](reports/PILOT_REPORT.md).

## Reproduce the executed work

Run from the repository root in PowerShell using the working x64 Python runtime.
The repository environment launcher supplies the numerical packages. These
commands reuse the source cache and do not access the network or score expression.

```powershell
$a0Python = 'C:/Users/dream/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$a0Root = 'RQ_Specified/A0_conserved_epithelial_transition_program'
$a0Scripts = @('audit_local_coverage.py', 'audit_public_metadata.py',
  'plot_p0_coverage.py', 'write_p0_report.py', 'audit_extended_coverage.py',
  'write_feasibility_report.py', 'verify_a0.py')
foreach ($a0Script in $a0Scripts) {
  & $a0Python analysis/scripts/run_with_environment.py `
    --site-packages .venv-x64/Lib/site-packages "$a0Root/scripts/$a0Script"
  if ($LASTEXITCODE -ne 0) { throw "A0 failed at $a0Script" }
}
```

`audit_public_metadata.py` refreshes the GEO metadata inventory, including the
skin study. `write_feasibility_report.py` applies the current reviewed candidate
decisions. The initial P0 report remains historical. After changing a figure,
visually inspect it before updating the figure-review record.

Public cache files are excluded from Git. `source_manifest.json` contains their
URLs, sizes, retrieval modes, dates and SHA-256 hashes. To fetch a missing source,
use `scripts/fetch_source_metadata.py NAME URL --max-bytes LIMIT`, setting LIMIT
above its recorded size. Add `--gzip-first-line` only for the intestinal matrix
header and `--accept application/octet-stream` for `Negretti_obs.bin`. The header
hash covers the saved decompressed header, not the remote expression matrix.
Full PDF downloads are source documents, not expression data.

Metadata downloads require network access in a sandboxed session. Previously
observed failures: NCBI's main www host did not resolve; its public FTP host over
HTTPS worked. Some PMC and publisher pages were unavailable to automated tools;
author code, GEO, the author viewer and accessible primary-study copies were used.
Do not treat an inaccessible document as evidence that metadata do not exist.

## Required cohort metadata

Provide one row per cell, with these fields or an explicit mapping to them:

| Field | Requirement |
|---|---|
| `cell_id` | Unique and matches the count matrix exactly |
| `library_id` | Original capture/sequencing library |
| `biological_unit_id` | Verified mouse/donor or prespecified independent pool |
| `unit_type` | Individual animal, donor, or independent pool |
| `pool_members` | If pooled, membership and evidence of non-overlap across units |
| `technical_replicate_of` | Connects multiple libraries from the same unit |
| `time_or_age` | Sampling time with units and a defined biological stage |
| `condition`, `genotype`, `sort_strategy` | Establish comparable sampling and exclusions |
| `author_state` | Source label, preserving uncertain/unassigned labels |
| `state_role` | Source-supported starting, intermediate or destination state for one branch |
| `label_evidence_source` | Paper/table/code and independent temporal, fate or tissue evidence |
| `qc_or_doublet_status` | Source QC provenance; missing status must be explicit |

Do not fill missing animal identities with library IDs. Combine technical repeats
only when their shared biological origin is documented. A donor with two libraries
is one unit. Independent pools count as pools, not as the number of their members.

## Conditions for resuming P1

1. D2 and V1 each have ≥3 independently identified units with ≥30 cells in all
   three states. Same-day replication is not obligatory if the paired across-time
   estimand and limits are explicitly frozen.
2. State definitions and transition evidence are reviewed without examining the
   candidate program's enrichment. Early embryonic and postnatal branches stay distinct.
3. Matching raw counts, gene IDs, technical-replicate relationships and assay
   coverage are available. The developmental viewer's SCT export is not raw counts.
4. Freeze selected cohorts, exclusions, label-gene exclusions, scoring, effect
   criteria and source hashes before discovery. Record all prior exposure.
5. Run the original P2–P4 sequence only after these requirements pass.

Capture estimates in `tables/developmental_capture_planning.csv` are illustrative
sampling calculations. They neither satisfy these conditions nor estimate power
for detecting or transferring a gene program.
