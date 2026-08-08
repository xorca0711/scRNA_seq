# What lives where in `analysis/`

Two GEO series were analysed. The primary one (mouse) writes to the **top level**
of this directory, because the task specified those canonical paths; the second
(human) writes to a subdirectory. That means the top level mixes
mouse-specific output with a few genuinely shared files, which is easy to
misread. This file removes the ambiguity.

## Ownership

| Path | Belongs to |
|---|---|
| `figures/` | **GSE262927** (mouse) |
| `tables/` | **GSE262927** (mouse) |
| `qc/` | **GSE262927** (mouse) |
| `logs/` | **GSE262927** (mouse) |
| `processed/` | **GSE262927** (mouse) |
| `inventory/` | **GSE262927** (mouse) |
| `epithelial_subanalysis/` | **GSE262927** (mouse) |
| `regeneration_focus/` | **GSE262927** (mouse) |
| `lineage_tracing_cohort/` | **GSE262927** (mouse) |
| `README.md` | **GSE262927** (mouse) — the report for the mouse series |
| `GSE178360/` | **human** — complete parallel structure inside |
| `scripts/` | **shared** — one pipeline, `--dataset` selects the series |
| `requirements.txt` | **shared** — one environment for both |
| `raw_data_inventory.txt` / `.csv` | **shared** — the scan covers both series |

Rule of thumb: **anything at the top level that is not `scripts/`,
`requirements.txt`, `raw_data_inventory.*` or `GSE178360/` is mouse.**

## The three mouse analyses, and how they differ

The mouse series was analysed three times, for three different questions. They
are not redundant.

| Directory | Cells | Question |
|---|---|---|
| top level (`figures/`, `tables/`, …) | 162,175 | The whole-atlas survey: all 33 samples, 29 clusters, marker tables, composition. The general-purpose object. |
| `regeneration_focus/` | 5,694 alveolar + 43,359 capillary | **The biology.** AT2 → Krt8⁺ transitional → AT1 trajectory, and the persistent injury-induced capillary state. Restricted to the 25-sample annotated cohort. |
| `lineage_tracing_cohort/` | 36,958 | The **8 non-atlas samples**, which are a separate experiment (Kit / Car4 / Ednrb Cre drivers, pre-labelled, 19 dpi). Asks where the injury state comes from. |

**Important:** the whole-atlas object at the top level contains all 33 samples,
i.e. it merges the 25-sample Ki67 atlas with the 8-sample tracing experiment.
That is fine for surveying cell types but wrong for anything condition- or
trace-related. The two focused directories each restrict to the correct cohort.
See `docs/ANALYSIS_RATIONALE.md`.

## Where to start

- **What was actually run, with parameters:** [`../docs/PIPELINE_AS_RUN.md`](../docs/PIPELINE_AS_RUN.md)
- **Why each decision was made:** [`../docs/ANALYSIS_RATIONALE.md`](../docs/ANALYSIS_RATIONALE.md)
- **Background if you are not a computational biologist:**
  [`../docs/UMAP_AND_FIGURES.md`](../docs/UMAP_AND_FIGURES.md),
  [`../docs/BACKGROUND_FOR_BIOLOGISTS.md`](../docs/BACKGROUND_FOR_BIOLOGISTS.md),
  [`../docs/DOUBLETS_AND_SCRUBLET.md`](../docs/DOUBLETS_AND_SCRUBLET.md)
- **Mouse report:** [`README.md`](README.md) · **Human report:** [`GSE178360/README.md`](GSE178360/README.md)
