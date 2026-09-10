# What lives where in `analysis/`

Two GEO series were analysed. **Each owns a subdirectory**, so the two are
symmetric and nothing at the top level is series-specific.

```
analysis/
├── GSE262927/        mouse  — 33 samples, H1N1 lineage-tracing time course
├── GSE178360/        human  — 3 samples, distal airway
├── scripts/          shared — one pipeline; --dataset selects the series
├── requirements.txt  shared — one environment for both
├── raw_data_inventory.txt / .csv   shared — the scan covers both series
└── LAYOUT.md         this file
```

> **Note on paths.** The original brief specified canonical output paths at the
> top level (`analysis/figures/umap/...`, `analysis/processed/final_clustered.h5ad`).
> Those now live under `analysis/GSE262927/`. The move was deliberate: with two
> series in one repository, putting one of them at the top level and the other in
> a subdirectory made it impossible to tell shared files from mouse-specific ones
> at a glance. Every required output still exists, one level deeper.

## Inside each series directory

Both have the same structure:

```
<SERIES>/
├── README.md              the report for this series
├── figures/               umap/ dotplots/ featureplots/ composition/ qc/
├── tables/                markers, annotations, composition, cell metadata
├── qc/                    thresholds, before/after, doublets, batch assessment
├── logs/                  analysis_log.txt, decisions.json
├── inventory/             detected format, per-sample results
├── processed/             .h5ad checkpoints  (gitignored — regenerable)
└── epithelial_subanalysis/
```

`GSE262927/` additionally contains four focused analyses (see below).

## The three mouse analyses differ by cohort — this matters

| Directory | Cells | Cohort | Purpose |
|---|---:|---|---|
| `GSE262927/` (top level) | 162,175 | **all 33 samples** | Whole-atlas survey: 29 clusters, marker tables, composition. General-purpose object. |
| `GSE262927/regeneration_focus/` | 5,694 alveolar<br>43,359 capillary | 25-sample annotated atlas | **The biology.** AT2 → Krt8⁺ transitional → AT1 trajectory, and the persistent injury-induced capillary state. |
| `GSE262927/lineage_tracing_cohort/` | 36,958 | **the 8 non-atlas samples** | A separate experiment (Kit / Car4 / Ednrb Cre, pre-labelled, 19 dpi). Where does the injury state come from? |
| `GSE262927/phase_timecourse/` | 107,626 | 25-sample annotated atlas | Per-day UMAP, per-animal lineage composition and Ki67-trace proliferation by lineage: the paper's three-phase structure, computed from the tracked metadata table only (script 10). Owner review pending. |
| `GSE262927/myeloid_focus/` | 9,997 | 25-sample annotated atlas | The myeloid compartment (atlas clusters 5, 17, 24) re-embedded with labels held out, graded against them, and read per animal and day (script 11). Owner review pending. |

**The whole-atlas object merges two experiments.** All 33 samples are in it,
which is fine for surveying cell types but wrong for anything condition- or
trace-related, because the 8 tracing samples use different Cre drivers and an
inverted labelling design. Both focused directories restrict to the correct
cohort. See [`../docs/ANALYSIS_RATIONALE.md`](../docs/ANALYSIS_RATIONALE.md).

## Where to start

| I want… | Read |
|---|---|
| what was actually run, with parameters | [`../docs/PIPELINE_AS_RUN.md`](../docs/PIPELINE_AS_RUN.md) |
| why each decision was made | [`../docs/ANALYSIS_RATIONALE.md`](../docs/ANALYSIS_RATIONALE.md) |
| background (non-computational reader) | [`../docs/UMAP_AND_FIGURES.md`](../docs/UMAP_AND_FIGURES.md), [`../docs/BACKGROUND_FOR_BIOLOGISTS.md`](../docs/BACKGROUND_FOR_BIOLOGISTS.md), [`../docs/DOUBLETS_AND_SCRUBLET.md`](../docs/DOUBLETS_AND_SCRUBLET.md) |
| the mouse report | [`GSE262927/README.md`](GSE262927/README.md) |
| the human report | [`GSE178360/README.md`](GSE178360/README.md) |

## Re-running

```bash
python analysis/scripts/01_scan_raw_data.py
python analysis/scripts/run_scrna_analysis.py --dataset GSE262927
python analysis/scripts/run_scrna_analysis.py --dataset GSE178360 --integration harmony
python analysis/scripts/06_regeneration_focus.py
python analysis/scripts/07_lineage_tracing_cohort.py
python analysis/scripts/10_phase_timecourse.py
python analysis/scripts/11_myeloid_focus.py
python analysis/scripts/03_write_report.py --dataset GSE262927
python analysis/scripts/05_write_pipeline_as_run.py
```
