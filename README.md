# Single-cell RNA-seq: lung regeneration after influenza injury

[![Portfolio checks](https://github.com/xorca0711/scRNA_seq/actions/workflows/portfolio-checks.yml/badge.svg)](https://github.com/xorca0711/scRNA_seq/actions/workflows/portfolio-checks.yml)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/)

An independent **Python/scanpy reanalysis of two public lung scRNA-seq
datasets, built from the raw deposited count matrices** — the authors'
processed objects and annotations were never used during model fitting. The
question: **can the published injury-and-regeneration biology be recovered by
an independent pipeline — and where it can't, why not?**

**Headline findings** (all numbers from tracked artefacts; details and
figures in [`FINDINGS.md`](FINDINGS.md)):

- The **AT2 → Krt8⁺ transitional → AT1** trajectory is recovered; held-out
  author labels are ordered correctly by pseudotime (median 0.013 → 0.179 →
  0.237 → 0.327).
- The transitional state behaves as a true intermediate in time: its median
  per-animal proportion is **27.4% at 11 dpi and 0.3% by 366 dpi**.
- The capillary injury state (iCAP) is its mirror image: median per-animal
  proportions are **2.0% → 37.5% at 25 dpi → still 21.7% at one year** — it
  never resolves.
- Lineage tracing supports a **CAP1 origin** for the injury state (traced at
  33–53% per animal in the Kit line); the CAP2 lines are reported as
  **uninformative, not negative**.
- Blind clustering scores **median purity 0.947** against the held-out author
  labels — with 3 of 29 candidate annotations contradicted and flagged, not
  hidden.
- The two datasets get **opposite integration decisions** (mouse: none;
  human: Harmony), each from measured evidence — the experimental design
  decides, not a default.

**Portfolio report:** [thesis-aware lung scRNA-seq synthesis](output/pdf/lung_scrna_portfolio_thesis_context.pdf) — human distal-lung and AT0 results first, with complementary regeneration, lineage-tracing, validation and limitation sections.

![Mouse atlas UMAP](analysis/GSE262927/figures/umap/UMAP_leiden_clusters_sidelegend.png)

## At a glance

| | |
|---|---|
| Primary dataset | [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927) — mouse, H1N1 injury time course, 33 samples, 162,175 cells (Niethamer et al., *Cell Stem Cell* 2025) |
| Second dataset | [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360) — human distal lung, 3 donors, 27,729 cells (Kadur Lakshminarasimha Murthy et al., *Nature* 2022) |
| Stack | Python 3.12, scanpy/AnnData, Scrublet, harmonypy, PAGA + diffusion pseudotime |
| Validation | Deposited author labels held out of all model fitting; used only post hoc as an answer key |
| Provenance | Decisions machine-logged; reports **generated from artefacts**, never hand-entered |
| Reproducibility | Dependency-free artefact validator + pinned full environment + CI |

## What this repository demonstrates

Experimental-design-aware scRNA-seq analysis · scanpy/AnnData and sparse
count workflows · per-sample QC with recorded rationales · per-capture
doublet handling with a prior-based fallback · batch-effect diagnosis (and
the discipline to not correct) · clustering, blind annotation, and external
grading · PAGA/diffusion pseudotime · lineage-trace calling · reproducible,
generated reporting · critical self-audit, including a refuted finding kept
on display.

## Where to go

| Question | Read |
|---|---|
| What was found? | [`FINDINGS.md`](FINDINGS.md) |
| What has been claimed, what stands behind each claim, and which results have potential? | [`CLAIMS.md`](CLAIMS.md) |
| What exactly ran, with parameters? | [`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md) (generated) |
| How can I validate or reproduce it? | [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) |
| Why each analytical decision? | [`docs/ANALYSIS_RATIONALE.md`](docs/ANALYSIS_RATIONALE.md) |
| Who did what — AI-assisted development and scientific ownership? | [`DEVELOPMENT.md`](DEVELOPMENT.md) |
| Machine/session context for AI assistants | [`AI_CONTEXT.md`](AI_CONTEXT.md) |
| Full per-dataset reports, figures, QC | [`analysis/GSE262927/`](analysis/GSE262927/README.md) · [`analysis/GSE178360/`](analysis/GSE178360/README.md) |
| Current state and known issues | [`PROGRESS.md`](PROGRESS.md) |
| The paper roadmap: study notes, extracted decision criteria, per-paper analysis trials | [`Thesis/README.md`](Thesis/README.md) |

## Repository map

```
├── FINDINGS.md                  results, with figures — start here
├── DEVELOPMENT.md               AI-assisted development disclosure + scientific ownership
├── AI_CONTEXT.md                machine-oriented context for AI sessions
├── REPRODUCIBILITY.md           input layout, validation tiers, full re-run guide
├── analysis/
│   ├── GSE262927/               mouse: report, figures, tables, QC, logs
│   │   ├── regeneration_focus/      AT2→AT1 trajectory + iCAP persistence
│   │   ├── lineage_tracing_cohort/  the 8-sample Cre-driver experiment
│   │   ├── phase_timecourse/        per-dpi UMAP, composition and proliferation by lineage (under review)
│   │   ├── myeloid_focus/           the myeloid compartment by dpi (under review)
│   │   └── epithelial_subanalysis/
│   ├── GSE178360/               human: report, figures, tables, QC, logs
│   ├── scripts/                 the pipeline (shared; --dataset selects series)
│   └── LAYOUT.md                what lives where, and why the cohorts differ
├── docs/
│   ├── PIPELINE_AS_RUN.md       what was actually executed (generated)
│   ├── ANALYSIS_RATIONALE.md    each decision, before and after the papers
│   ├── BACKGROUND_FOR_BIOLOGISTS.md · UMAP_AND_FIGURES.md · DOUBLETS_AND_SCRUBLET.md
│   └── SOUPX/SCRUBLET/SCDS/SLINGSHOT/TRADESEQ.md   tool reference notes
├── Thesis/                      paper roadmap: study notes, extracted decision criteria, per-paper analysis trials
├── WORKFLOW.md                  the reference study's published workflow
├── REFERENCES.md                all papers, DOIs, data accessions
└── PROGRESS.md                  session/handoff state, known issues
```

`raw_data/` (7.8 GB of GEO downloads) and all regenerable `.h5ad` objects are
gitignored; figures and small tables are tracked. No sequencing data, count
matrices or paper PDFs are committed.

## Reproducing

The tracked portfolio can be checked without downloading data or installing
scanpy:

```bash
python analysis/scripts/validate_portfolio.py
python -m compileall -q analysis/scripts
```

For a complete re-run, create a Python 3.12 environment, install the pinned
dependencies, and place the GEO downloads under `raw_data/<accession>/`:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r analysis/requirements.txt

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
python analysis/scripts/09_write_portfolio_pdf.py
python analysis/scripts/05_write_pipeline_as_run.py
```

Exact inputs, expected directory layout, validation tiers, resource notes, and
the output contract are in [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md). Seeds
are fixed at 0 throughout; `--stages` resumes from checkpoints. The original
Windows ARM64 execution constraints are documented in [`AI_CONTEXT.md`](AI_CONTEXT.md)
and [`PROGRESS.md`](PROGRESS.md).

## Relationship to the published pipeline

The reference study's own workflow (STARsolo → SoupX → scds+Scrublet →
Seurat/SCTransform → Slingshot → tradeSeq, in R) is documented in
[`WORKFLOW.md`](WORKFLOW.md), per-tool schematics in
[`docs/`](docs/README.md), and an annotated parameter-level reference in
[`docs/scRNAseq_workflow_Niethamer2025.md`](docs/scRNAseq_workflow_Niethamer2025.md).
**The analysis here is deliberately not a port of that pipeline** — it was
built from the deposited data alone, then compared against the papers after
the fact. The tool-by-tool used/not-used table and every divergence are in
[`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md).

## Licence

Code and original written material are © 2026 Xorca; no reuse licence is
granted. See [`LICENSE`](LICENSE). The source papers and public datasets remain
the property of their respective authors and publishers; citations and
open-access links are in [`REFERENCES.md`](REFERENCES.md).
