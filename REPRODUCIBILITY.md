# Reproducibility guide

This repository supports three levels of verification. The first two use only
tracked files; the full analysis requires the public GEO downloads and the
pinned scientific Python environment.

## 1. Review the result without computation

Start with [`FINDINGS.md`](FINDINGS.md), then follow each figure to its tracked
table and generated per-dataset report. [`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md)
is the authoritative record of what ran. It is generated from the decision
logs rather than maintained as a second handwritten method description.

## 2. Validate the tracked portfolio

Python 3.12 is recommended; no third-party package is needed:

```bash
python analysis/scripts/validate_portfolio.py
python -m compileall -q analysis/scripts
```

The validator checks local documentation links, parses every tracked JSON
artefact, and recomputes the headline counts, medians, cluster purity, and
lineage-tracing range from the tracked CSV tables. GitHub Actions runs the same
checks on every push and pull request.

## 3. Re-run from deposited count matrices

### Inputs

Download the supplementary files from the GEO accession pages linked below and
extract them to this layout:

```text
raw_data/
├── GSE262927/
│   ├── GSE262927_RAW/                 # 33 GSM*.h5 matrices
│   └── GSE262927_CellMetaData.csv
└── GSE178360/
    └── GSE178360_RAW/                 # 3 *filtered_feature_bc_matrix.h5 files
```

- [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927):
  `GSE262927_RAW.tar` and `GSE262927_CellMetaData.csv.gz`
- [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360):
  `GSE178360_RAW.tar`

The exact 51-file inventory from the original run is tracked in
[`analysis/raw_data_inventory.txt`](analysis/raw_data_inventory.txt),
[`analysis/raw_data_inventory.csv`](analysis/raw_data_inventory.csv), and
[`analysis/raw_data_inventory.json`](analysis/raw_data_inventory.json). Only
the count matrices and mouse metadata shown above are analysis inputs; the
downloaded author-supplied RDS objects were inspected as deposited-file context
but were not used for model fitting.

`raw_data/` is gitignored and must remain read-only. The original download set
occupied 7.8 GB. Allow additional space for regenerable AnnData checkpoints;
the complete mouse workflow takes multiple hours on a workstation.

### Environment

```bash
python -m venv .venv
# Windows
.venv\Scripts\activate
# macOS/Linux
source .venv/bin/activate

python -m pip install --upgrade pip
python -m pip install -r analysis/requirements.txt
```

The lock file records the exact Python 3.12 environment used for the final
run. The originating Windows ARM64 machine required an emulated x86-64 Python;
that machine-specific detail is not a general requirement and is documented in
[`AI_CONTEXT.md`](AI_CONTEXT.md).

### Execution order

Run commands from the repository root:

```bash
python analysis/scripts/01_scan_raw_data.py
python analysis/scripts/run_scrna_analysis.py --dataset GSE262927
python analysis/scripts/run_scrna_analysis.py --dataset GSE178360 --integration harmony
python analysis/scripts/06_regeneration_focus.py
python analysis/scripts/07_lineage_tracing_cohort.py
python analysis/scripts/03_write_report.py --dataset GSE262927
python analysis/scripts/03_write_report.py --dataset GSE178360
python analysis/scripts/08_reference_aligned_epithelial_umap.py
python analysis/scripts/09_write_portfolio_pdf.py
python analysis/scripts/05_write_pipeline_as_run.py
python analysis/scripts/validate_portfolio.py
```

All random seeds are fixed at 0. The human Harmony selection is an explicit,
recorded override justified in [`docs/ANALYSIS_RATIONALE.md`](docs/ANALYSIS_RATIONALE.md).
For stage-level recovery, use `--stages` with a comma-separated subset shown by
`python analysis/scripts/run_scrna_analysis.py --help`; invalid names fail
fast. `--reuse-markers` reuses a compatible marker table during iteration.

### Output contract

Each dataset directory contains:

- `logs/decisions.json` and `logs/analysis_log.txt` for provenance;
- `qc/` for thresholds and automated review;
- `tables/` for compact, inspectable evidence;
- `figures/` for portfolio-ready results;
- `processed/` for regenerable AnnData checkpoints (gitignored).

The large raw inputs and processed objects are intentionally excluded from
Git. Figures, compact tables, decision logs, and generated reports are tracked
so the scientific claims remain auditable without recomputing the analysis.
