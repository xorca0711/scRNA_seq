# Progress and handoff state

Living record of what is done, what is pending, and what a future session needs
to know to continue. Update this before stopping.

Last updated: 2026-08-08.

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

Branch: `scrna-adaptive-pipeline` → see `git log`. PR #1 merged; **PR #2 open**:
https://github.com/xorca0711/scRNA_seq/pull/2

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
python analysis/scripts/03_write_report.py --dataset GSE262927
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
**Harmony primary** (justified in `qc/celltype_split_by_sample.csv`); donor-private
clusters fell 20/41 → 4/31 and `donor_driven_clustering_check` now reports false.

---

## Known issues carried forward (IMPORTANT)

1. **The mouse object merges two different experiments.** The 8 samples with no
   metadata (`EEM-scRNA-249/250/251/288/289/290/291/292`) use different Cre
   drivers (Kit-MerCreMer, Car4-CreERT2, Ednrb-CreERT2), a pre-labelling design
   and a single 19 dpi timepoint. Niethamer et al. analysed them **separately**
   from the 25-sample Ki67 atlas. To reproduce the paper, split them out.
   `trace_call` means different things in the two cohorts.
2. **Mouse composition outputs describe a sort ratio, not the lung.** Cells were
   MACS-fractionated and recombined at 85–90% CD45⁻ : 10–15% CD45⁺. Compare
   composition only *within* a compartment.
3. **Scrublet is biased against the human paper's discoveries.** Measured:
   distal-BC 12.9%, SCGB3A2-CC 10.5%, AT0 8.0% called doublet vs 6.3% baseline
   (up to 2×). Not erased (87–92% survived) but systematic. Kadur et al. used
   **no** automated doublet caller.
4. **No ambient-RNA correction anywhere.** Both papers used SoupX. Matters most
   for AT0 (defined partly by SFTPC, a classic soup contaminant) and for
   GSE178360 DD046Q (heavy haemoglobin).
5. **Doublet calls are a ranking cut, not a detected threshold** in 32/33 mouse
   samples (non-bimodal Scrublet histograms → expected-rate quantile fallback).
6. **MAD upper bounds never bind** — zero cells removed for excess counts/genes.
7. **The papers' headline findings need subsetting.** Niethamer's iCAP state is
   invisible at top-level clustering; only the epithelium was subset here, not
   the endothelium.
8. **The deposited `.RDS` objects contain the authors' labels** (donor IDs and
   published cell types per barcode). Unreadable here for want of R, but they
   are the route to ground-truth labels.
9. Mouse `decisions.json` was partly reconstructed from run logs by
   `04_recover_decisions.py` (decision persistence was added mid-session). Fresh
   runs do not need this.

---

## Pending work

1. **Optional analysis follow-ups**, in value order:
   - Split the 8 non-atlas mouse samples into their own object (known issue 1).
   - Endothelial sub-analysis to recover Niethamer's iCAP state (known issue 7).
   - Re-run GSE178360 without doublet filtering to check AT0 recovery
     (known issue 3), ideally with ambient correction (known issue 4).
2. `WORKFLOW.md` has a scope banner, but its QC "acceptance order" section still
   reads as an imperative checklist for this repo. Low priority.

---

## Repository conventions

- `raw_data/` is **read-only** and gitignored (7.8 GB). Never modify it.
- Large regenerable outputs are gitignored: `*.h5ad`, `processed/`,
  `sample_shards/`, `cluster_markers_all.csv`, `*.csv.gz`. Figures and small
  tables **are** tracked.
- `docs/PIPELINE_AS_RUN.md` is **generated** — edit
  `analysis/scripts/05_write_pipeline_as_run.py` and re-run, never the `.md`.
- The tool reference pages under `docs/` are deliberately preserved; they
  describe the published method, not what ran here.
