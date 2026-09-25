# Reproducibility guide

This repository supports three levels of verification. The first two use only
tracked files; the full analysis requires the public GEO downloads and the
pinned scientific Python environment.

## 1. Review the result without computation

Start with [`docs/PORTFOLIO.md`](docs/PORTFOLIO.md) and [`RESEARCH_QUESTIONS.md`](RESEARCH_QUESTIONS.md), then follow each figure to its tracked
table and generated per-dataset report. [`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md)
is the record of the original atlas pipeline, not every later trial. It is generated from the decision
logs rather than maintained as a second handwritten method description.

## 2. Validate the tracked repository

Python 3.12 is recommended; no third-party package is needed:

```bash
python analysis/scripts/validate_repository.py
python analysis/scripts/claim_contract.py --check
python "Research Article/gate1_03_nabhan_2018/nb1/verify_outputs.py"
python -m unittest discover -s analysis/tests -q
python -m compileall -q analysis "Research Article" RQ_Specified
```

The validator checks local documentation links, parses every tracked JSON
artefact, and recomputes the headline counts, medians, cluster purity, and
lineage-tracing range from the tracked CSV tables. GitHub Actions runs the same
checks on every push and pull request. The generated claim manifest binds
selected values to explicit artifacts, filters and aggregation rules. Its
coverage is stated per row; CI does not reproduce all biological claims.

Local portable runtimes may contain ignored third-party Python 2 files. In a
workspace containing `analysis/corrections/statistics/.tools/`, compile the
repository sources with `python -m compileall -q -x '[/\\]\.tools[/\\]' analysis "Research Article" RQ_Specified`.
Those ignored runtime files are absent from a clean checkout and CI.

## Nabhan 2018 and Nb1

The [Nabhan analysis](Research%20Article/gate1_03_nabhan_2018/README.md) has three separate
reproduction paths: a 47-cell deposited FPKM panel, an animal-level raw-count
analysis of GSE262927, and acquisition/eligibility checks for an external
cohort. Each report supplies commands and input provenance. These paths do not
pool FPKM with UMI counts. The local Nb1 run produces five PNG/SVG figures and
seven compact tables; its standard-library verifier can check tracked results
without the large input matrix. Full execution requires the scientific
environment and source inputs listed in its run record.

## Scientific correction workflows (September 2026)

Historical trials are preserved. Corrected results and their run specifications
live in [ligand corrections](analysis/corrections/ligand/README.md),
[statistical corrections](analysis/corrections/statistics/README.md), and
[epithelial specificity](Research%20Article/epithelial_state_specificity/README.md). Those
pages give the exact scripts, inputs, outputs and seeds for each pass. Review
them before launching a large data stream. Large caches and the local R runtime
are ignored by Git; compact results, plots, source definitions and provenance
are retained.

### Python environment and recovery

The original atlas requirements are in `analysis/requirements.txt`; the later
Windows x86-64 stack is pinned in
[`analysis/config/requirements-x64.txt`](analysis/config/requirements-x64.txt).
They serve different workflows. Use a stable installed CPython 3.12 interpreter
of the same architecture as the binary packages. The lock file documents the
rebuild commands; a temporary interpreter location should not be used as the
base of a durable environment.

In the reviewed checkout, both old venv launchers pointed to unavailable base
interpreters. A compatible working CPython 3.12 x64 successfully reused the
existing `.venv-x64/Lib/site-packages`, without reinstalling the analysis stack.
The reusable launcher is explicit about this recovery:

```bash
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages --check
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages analysis/scripts/15_claims_ledger_figure.py
```

Here `python` must be the working compatible interpreter, not the broken venv
launcher. This verifies usable installed packages; it is not proof of a fresh
environment rebuild on another machine. On other platforms create a native
compatible environment and follow the workflow-specific dependencies.

### Reference statistics and provenance

The statistical corrections use official R packages rather than naming a
custom approximation CAMERA. Runtime and package versions, commands and output
tables are recorded with that correction. The original W1 output is unchanged.

The shared run-record helper now content-hashes declared inputs, records code
identity, writes atomically, and preserves previous record bytes in a
content-addressed `.history/` directory before another run. This improves future
runs; it cannot retroactively preregister historical analyses. Corrected passes
are labelled as post-audit specifications on already inspected data.

After an authorized claim edit, regenerate the evidence index and displays:

```bash
python analysis/scripts/claim_contract.py
python analysis/scripts/14_write_negative_results.py
python analysis/scripts/15_claims_ledger_figure.py
```

The ledger figure needs matplotlib; the first two commands use the standard
library. Numeric bindings should change only with a documented evidence change,
not merely to silence a failed check.

## 3. Re-run from deposited count matrices

### Inputs

Download the supplementary files from the GEO accession pages linked below and
extract them to this layout:

```text
raw_data/
├── GSE262927/                         # the two original series; analysed under Research Article/, beside their papers
│   ├── GSE262927_RAW/                 #   33 GSM*.h5 matrices
│   └── GSE262927_CellMetaData.csv
├── GSE178360/
│   └── GSE178360_RAW/                 #   3 *filtered_feature_bc_matrix.h5 files
├── GSE145031/  GSE144468/             # Choi 2020: trials D0 to D7 (MatrixMarket triplets, SOFT)
├── GSE316241/  GSE316243/  GSE316244/ # Cardoso 2026: trials C0 to C14 (h5 and SOFT)
├── GSE310335/                         # Cardoso 2026 human organoids
├── GSE247505/                         # England 2025: GSE247505_RAW (mtx triplets), SOFT
├── GSE131907/  GSE136831/             # human LUAD and IPF: trials E1 to E6, C12, C14
├── GSE135893/  GSE132771/             # human IPF; mouse bleomycin: trials E3, E4
├── GSE310539/                         # Lynch 2026 multiome: *_filtered_feature_bc_matrix.h5,
│                                      #   *_atac_peak_annotation.tsv.gz, SOFT
└── GSE247130/                         # Hassan and Chen 2024 multiome: three Aggregate_* h5,
                                       #   three peak annotations, SOFT
```

The two original series need the files named below; every roadmap deposit is
inventoried by its own reality-check trial (C0, D0, M0), whose run record lists
the exact files, sizes and modification times it read.

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
python analysis/scripts/10_phase_timecourse.py
python analysis/scripts/11_myeloid_focus.py
python analysis/scripts/12_amac_trace_by_window.py
python analysis/scripts/13_myeloid_batch_sensitivity.py
python analysis/scripts/03_write_report.py --dataset GSE262927
python analysis/scripts/03_write_report.py --dataset GSE178360
python analysis/scripts/08_reference_aligned_epithelial_umap.py
python analysis/scripts/05_write_pipeline_as_run.py
python analysis/scripts/validate_repository.py
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
- `figures/` for tracked figures;
- `processed/` for regenerable AnnData checkpoints (gitignored).

The large raw inputs and processed objects are intentionally excluded from
Git. Figures, compact tables, decision logs, and generated reports are tracked
so the scientific claims remain auditable without recomputing the analysis.

## Shared research-question figures

[Questions A1–A14](RESEARCH_QUESTIONS.md) are the canonical register. Scripts
16 and 18 produce the earlier shared figures; scripts 19–21 prepare, render
and validate the IL-1 context package. See its
[methods, commands and preserved provenance](analysis/figures/rq/il1b_context/REPORT.md).
Rendering uses tracked tables. Full preparation and input-hash validation
require ignored source caches; relocation did not rerun embeddings or inference.
Paper-specific release validation covers the 17 Yu study figures separately.
