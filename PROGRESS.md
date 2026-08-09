# Progress and handoff state

Living record of what is done, what is pending, and what a future session needs
to know to continue. Update this before stopping.

Last updated: 2026-08-09. **All planned work is complete**, and the repository
has been restructured as a portfolio: `FINDINGS.md` (results with figures) now
leads, `README.md` is a landing page for the executed analysis, and
`scRNAseq_workflow_Niethamer2025.md` moved to `docs/`. Nothing is running.

---

## Status at a glance

| Item | State |
|---|---|
| Raw-data scan + inventory | **DONE** |
| GSE262927 (mouse) full analysis | **DONE** — 162,175 cells, 33 samples, 29 clusters |
| GSE178360 (human) full analysis | **DONE** — 27,729 cells, 3 samples, 31 clusters, Harmony primary |
| Reports (`README.md` per dataset, analysis logs) | **DONE** |
| `docs/PIPELINE_AS_RUN.md` (what actually ran) | **DONE** |
| Source papers read and divergences documented | **DONE** |
| Doc scope labels (README/WORKFLOW/docs) | **DONE** — banners on all 5 tool pages, README/WORKFLOW relabelled |
| `docs/ANALYSIS_RATIONALE.md` (decisions, before vs after papers) | **DONE** |
| `docs/BACKGROUND_FOR_BIOLOGISTS.md` — Harmony | **DONE** |
| `docs/UMAP_AND_FIGURES.md` | **DONE** |
| `docs/DOUBLETS_AND_SCRUBLET.md` | **DONE** |
| **Focused alveolar regeneration re-analysis** | **DONE** — AT2→transitional→AT1 trajectory recovered |
| Capillary endothelial (iCAP) sub-analysis | **DONE** — persistent injury state recovered |
| AT0 doublet-loss follow-up | **DONE** — tested, not lost |
| Lineage-tracing cohort (8 non-atlas samples) | **DONE** — CAP1 origin supported; CAP2 lines uninformative |
| `analysis/` reorganised into two series subdirectories | **DONE** |
| Pipeline speedups (PCA, scan, threads) | **DONE** — verified output-identical |
| Portfolio restructure (`FINDINGS.md`, README rewrite, root tidy-up) | **DONE** |

Repository state: the scientific analysis and portfolio curation are complete;
no analysis process is running. Branch-specific state belongs in Git/GitHub,
not in this durable handoff document.

---

## Environment (needed to run anything)

This machine is **Windows on ARM64**. `numba`, `llvmlite` and `leidenalg` have
no ARM64 wheels and there is **no C compiler**, so scanpy cannot be built
against the native interpreter. The pipeline runs on an **x86-64 CPython 3.12**
interpreter under emulation:

```
.venv-x64/Scripts/python.exe          <- use THIS for everything
```

Created with `uv`. `harmonypy` is pinned to the pure-Python `0.0.10` (newer
releases need a C++ toolchain). **No R, no Seurat, no Bioconductor.**
Lockfile: `analysis/requirements.txt`.

The console is on a legacy codepage; scripts call
`sys.stdout.reconfigure(encoding="utf-8")` because printing `—` otherwise
raises `UnicodeEncodeError`.

---

## How to re-run

```bash
python analysis/scripts/01_scan_raw_data.py
python analysis/scripts/run_scrna_analysis.py --dataset GSE262927
python analysis/scripts/run_scrna_analysis.py --dataset GSE178360 --integration harmony
python analysis/scripts/06_regeneration_focus.py
python analysis/scripts/07_lineage_tracing_cohort.py
python analysis/scripts/03_write_report.py --dataset GSE262927
python analysis/scripts/03_write_report.py --dataset GSE178360
python analysis/scripts/05_write_pipeline_as_run.py
```

Useful flags: `--stages a,b,c` resumes from checkpoints; `--reuse-markers`
re-uses the ~45 min Wilcoxon table; `--stages log` rewrites the analysis log
from persisted decisions with no compute.

Timings (this machine): per-sample stage ~22 min for 33 samples on 3 workers;
merge+normalise ~24 min; PCA ~2 min; neighbours+UMAP ~34 min; resolution scan
~13 min; Wilcoxon ~45 min; figures+epithelial+finalize ~37 min.

---

## Key results

**GSE262927 (mouse, primary).** 212,701 barcodes → 169,807 after per-sample MAD
QC → 162,175 after Scrublet. 2,500 HVGs, 50 PCs, **no batch correction**, Leiden
0.3 → **29 clusters**. Median cluster purity against the deposited author
cell-type labels **0.947** across 107,626 cells (labels held out from the
clustering). Epithelial sub-analysis: 13,333 cells → 16 subclusters.

**GSE178360 (human).** 36,464 → 29,605 after QC → 27,729 after Scrublet.
**Harmony primary** (justified in `qc/celltype_split_by_sample.csv`); 4/31
clusters are >75% one donor and `donor_driven_clustering_check` reports false.
(An earlier session quoted "20/41 → 4/31"; the uncorrected baseline count is
not preserved in current artefacts, so only the after-state is citable.)

**Focused analyses (mouse).** `regeneration_focus/`: the AT2 -> Krt8+
transitional -> AT1 axis, recovered. Pseudotime orders the deposited labels
correctly (AT2 0.013, transitional 0.179, AT1_AT2 0.237, AT1 0.327) with those
labels held out; the median per-animal transitional proportion peaks at 27.4%
at 11 dpi and falls to 0.3% by 366 dpi. The capillary injury state behaves
oppositely - its median per-animal proportion persists (2.0% -> 37.5% at 25
dpi -> 21.7% at 366 dpi).
`lineage_tracing_cohort/`: the trace-call rule reproduces the authors' own
labels at 100.0000% over 107,626 cells, and in the Kit line (labels CAP1) the
injury state is traced at 33-53% per animal, consistent with a CAP1 origin. The
CAP2-specific lines are **uninformative, not negative** - they label only 2-8%
of total endothelium, so a near-zero rate cannot be distinguished from too few
labelled cells.

---

## Known issues carried forward (IMPORTANT)

1. **The mouse object merges two different experiments.** *(ADDRESSED in the
   focused re-analysis, which restricts to the annotated 25-sample cohort; still
   true of `processed/final_clustered.h5ad` itself.)* The 8 samples with no
   metadata (`EEM-scRNA-249/250/251/288/289/290/291/292`) use different Cre
   drivers (Kit-MerCreMer, Car4-CreERT2, Ednrb-CreERT2), a pre-labelling design
   and a single 19 dpi timepoint. Niethamer et al. analysed them **separately**
   from the 25-sample Ki67 atlas. To reproduce the paper, split them out.
   `trace_call` means different things in the two cohorts.
2. **Mouse composition outputs describe a sort ratio, not the lung.** Cells were
   MACS-fractionated and recombined at 85–90% CD45⁻ : 10–15% CD45⁺. Compare
   composition only *within* a compartment.
3. **Scrublet bias — TESTED AND LARGELY EXONERATED.** A crude co-expression
   gate suggested up to 2× over-removal of the human paper's novel populations.
   Re-measured with a gate that can separate real AT0 from AT2+club doublets
   (`SFTPC+ SCGB3A2+ EPCAM+`, lineage-negative), AT0 is flagged at **3.9% vs a
   6.3% baseline** — below background. The cells flagged inside the crude gate
   have 1.8× the UMIs and more cross-lineage co-expression, i.e. they are
   doublets. Lesson: a co-expression gate cannot audit a co-expression
   detector. General caution still applies; the verdict for this dataset does
   not.
4. **No ambient-RNA correction anywhere.** Both papers used SoupX. Matters most
   for AT0 (defined partly by SFTPC, a classic soup contaminant) and for
   GSE178360 DD046Q (heavy haemoglobin).
5. **Doublet calls are a ranking cut, not a detected threshold** in 32/33 mouse
   samples (non-bimodal Scrublet histograms → expected-rate quantile fallback).
6. **MAD upper bounds never bind** — zero cells removed for excess counts/genes.
7. **The papers' headline findings need subsetting.** *(ADDRESSED —
   `analysis/scripts/06_regeneration_focus.py` subsets both the alveolar
   epithelium and the capillary endothelium.)*
8. **The deposited `.RDS` objects contain the authors' labels** (donor IDs and
   published cell types per barcode). Unreadable here for want of R, but they
   are the route to ground-truth labels.
9. Mouse `decisions.json` was partly reconstructed from run logs by
   `04_recover_decisions.py` (decision persistence was added mid-session). Fresh
   runs do not need this.

10. **Candidate annotations are contradicted by the deposited labels for 3 of
   29 mouse clusters** (0, 22, 25) and are flagged `[CONTRADICTED]` in
   `tables/cluster_annotation_proposals.csv` and in the annotation UMAP.
   Cluster 0 is the important one: the marker panel called it "Transitional"
   but it is 89.2% CAP1 capillary endothelium (its top genes are
   interferon-stimulated). 21 of 29 agree. Trust the deposited label where the
   two conflict.

---

## Optional extensions (not completion gaps)

1. Add ambient-RNA correction if compatible unfiltered mouse droplet matrices
   become available; the current GEO deposit does not provide them.
2. Fit a formal trajectory-DE model if gene-level inference along pseudotime
   becomes a project goal; the current analysis intentionally stops at binned
   programme summaries.

---

## Directory layout

`analysis/` is split by series: `analysis/GSE262927/` (mouse) and
`analysis/GSE178360/` (human), with `scripts/`, `requirements.txt` and
`raw_data_inventory.*` shared at the top level. This departs from the canonical
paths in the original brief, which put the mouse outputs at the top level - see
`analysis/LAYOUT.md` for the mapping and for why the three mouse analyses
(whole atlas / regeneration_focus / lineage_tracing_cohort) use different
cohorts.

## Repository conventions

- `raw_data/` is **read-only** and gitignored (7.8 GB). Never modify it.
- Large regenerable outputs are gitignored: `*.h5ad`, `processed/`,
  `sample_shards/`, `cluster_markers_all.csv`, `*.csv.gz`. Figures and small
  tables **are** tracked.
- `docs/PIPELINE_AS_RUN.md` is **generated** — edit
  `analysis/scripts/05_write_pipeline_as_run.py` and re-run, never the `.md`.
- The tool reference pages under `docs/` are deliberately preserved; they
  describe the published method, not what ran here.
