# Single-cell RNA-seq: lung regeneration after influenza injury

An end-to-end, from-scratch **Python/scanpy reanalysis of two public GEO
datasets** — a 162,175-cell mouse influenza time course and a 27,729-cell
human distal-lung atlas — that recovers the source papers' headline biology
from the raw deposited count matrices, validates every unsupervised step
against the authors' held-out annotations, and documents its own failures as
carefully as its successes.

**→ [`FINDINGS.md`](FINDINGS.md) — the results, with figures.**

![Mouse atlas UMAP](analysis/GSE262927/figures/umap/UMAP_leiden_clusters_sidelegend.png)

## At a glance

| | |
|---|---|
| Primary dataset | [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927) — mouse, H1N1 injury, 33 samples, uninjured → 1 year post-infection (Niethamer et al., *Cell Stem Cell* 2025) |
| Second dataset | [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360) — human distal lung, 3 donors (Kadur Lakshminarasimha Murthy et al., *Nature* 2022) |
| Stack | Python 3.12, scanpy, Scrublet, harmonypy, PAGA + diffusion pseudotime |
| Validation | Deposited author labels held out of all clustering/trajectory steps, used only as an answer key — median cluster purity **0.947** |
| Headline results | AT2 → Krt8⁺ transitional → AT1 trajectory recovered; the persistent capillary injury state (iCAP) recovered; lineage-tracing cohort supports a CAP1 origin |
| Provenance | Every threshold, decision and substitution is machine-logged (`decisions.json`); the reports and [`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md) are **generated from those artefacts**, never hand-entered |

## What this repository demonstrates

- **Adaptive pipeline design.** The pipeline starts by scanning `raw_data/`
  and letting the files determine the workflow — format detection, species
  detection from gene symbols, per-sample QC thresholds derived from each
  sample's own distributions (median ± MAD, with a rationale string attached
  to every threshold).
- **Defensible statistical decisions.** Batch correction is treated as a
  hypothesis test, not a default — and the two datasets resolve it in
  opposite directions, each for a stated, measured reason
  ([`docs/ANALYSIS_RATIONALE.md`](docs/ANALYSIS_RATIONALE.md)).
- **Blind validation.** Clustering and annotation were done from marker
  panels alone; the deposited labels grade the result afterwards. Agreements
  and the three contradicted clusters are reported with equal prominence.
- **Honest negative results.** A suspected doublet-caller bias against a rare
  population was chased down and overturned by a better-designed test;
  thresholds that never bind are labelled as such
  ([`FINDINGS.md` § 6](FINDINGS.md#6--negative-results-and-self-audits)).
- **Reproducibility under constraints.** Fixed seeds throughout, a pinned
  lockfile, checkpointed stages with resume flags — all running on Windows
  ARM64 via an emulated x86-64 interpreter because key wheels don't exist
  natively ([`PROGRESS.md`](PROGRESS.md)).

## Repository map

```
├── FINDINGS.md                  the results — start here
├── analysis/
│   ├── GSE262927/               mouse: report, figures, tables, QC, logs
│   │   ├── regeneration_focus/      AT2→AT1 trajectory + iCAP persistence
│   │   ├── lineage_tracing_cohort/  the 8-sample Cre-driver experiment
│   │   └── epithelial_subanalysis/
│   ├── GSE178360/               human: report, figures, tables, QC, logs
│   ├── scripts/                 the pipeline (shared; --dataset selects series)
│   └── LAYOUT.md                what lives where, and why the cohorts differ
├── docs/
│   ├── PIPELINE_AS_RUN.md       what was actually executed (generated)
│   ├── ANALYSIS_RATIONALE.md    each decision, before and after the papers
│   ├── BACKGROUND_FOR_BIOLOGISTS.md · UMAP_AND_FIGURES.md · DOUBLETS_AND_SCRUBLET.md
│   └── SOUPX/SCRUBLET/SCDS/SLINGSHOT/TRADESEQ.md   tool reference notes
├── WORKFLOW.md                  the reference study's published workflow
├── REFERENCES.md                all papers, DOIs, data accessions
└── PROGRESS.md                  session/handoff state, known issues
```

`raw_data/` (7.8 GB of GEO downloads) and all regenerable `.h5ad` objects are
gitignored; figures and small tables are tracked. No sequencing data, count
matrices or paper PDFs are committed.

## Reproducing

```bash
pip install -r analysis/requirements.txt
python analysis/scripts/01_scan_raw_data.py
python analysis/scripts/run_scrna_analysis.py --dataset GSE262927
python analysis/scripts/run_scrna_analysis.py --dataset GSE178360 --integration harmony
python analysis/scripts/06_regeneration_focus.py
python analysis/scripts/07_lineage_tracing_cohort.py
python analysis/scripts/03_write_report.py --dataset GSE262927
python analysis/scripts/05_write_pipeline_as_run.py
```

Inputs are the two GEO series downloaded into `raw_data/<accession>/`
(accessions and exact files used: `analysis/raw_data_inventory.txt`). Seeds
are fixed at 0 throughout; `--stages` resumes from checkpoints. End-to-end
wall time on the development machine was roughly 3 hours per full mouse run.

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

Written material in this repository is © the author. The papers it describes
are the property of their respective publishers; see
[`REFERENCES.md`](REFERENCES.md) for links to the open-access versions.
