# What lives where in `analysis/`

`analysis/` holds what is shared across deposits: the pipeline, its
configuration, the repository-level figures and the raw-data inventory.
Paper-specific data and inference live under `Thesis/`, beside the study that
motivated them, since 2026-09-21. Cross-study corrections remain in
`analysis/corrections/`. Shared RQ figures and their compact plotted tables
belong in `analysis/figures/rq/`, even when they draw on one study's results.
The [repository structure contract](../docs/REPOSITORY_STRUCTURE.md) defines
canonical paths and identifier scope across the whole repository.

```
analysis/
├── scripts/          the pipeline (run_scrna_analysis.py; --dataset selects the series),
│                     the focused analyses of the mouse series (06 to 13), the generators
│                     (03, 05, 14, 15), shared RQ figures (16, 18–21) and the validator
├── config/           the one validated palette; the x86-64 environment lock
├── figures/          repository-level claims and question figures
│   └── rq/           shared RQ PNG/SVGs; il1b_context/ holds tables and provenance
├── lib/              reusable provenance and validation helpers
├── corrections/      explicitly separated corrective analyses and records
├── requirements.txt  the native environment
├── raw_data_inventory.txt / .csv / .json   the Stage 0 scan of raw_data/ (51 files)
└── LAYOUT.md         this file
```

## Where the two original series went

| Series | Location | Why there |
|---|---|---|
| GSE262927, mouse, 33 samples, H1N1 time course | [`Thesis/gate1_01_niethamer_2025/GSE262927/`](../Thesis/gate1_01_niethamer_2025/GSE262927/README.md) | Niethamer et al. 2025 deposited it; roadmap paper 1 |
| GSE178360, human, 3 donors, distal lung | [`Thesis/ungated_murthy_2022/GSE178360/`](../Thesis/ungated_murthy_2022/GSE178360/README.md) | Kadur Lakshminarasimha Murthy et al. 2022 deposited it; the paper is outside the roadmap, so the folder carries only a pointer note |

`pipeline_utils.SERIES_DIRS` is the one place the locations are written;
every script and trial derives its paths from it or from the equivalent
`REPO / "Thesis" / ...` constant. Run records and logs written before the
move keep the `analysis/GSE...` paths of their day; they are artefacts, not
pointers, and are not edited.

> **History of the paths.** The original brief specified outputs at the top
> level (`analysis/figures/umap/...`, `analysis/processed/final_clustered.h5ad`).
> They moved once to `analysis/GSE262927/`, so that the two series were
> symmetric and nothing at the top level was series-specific, and once more
> to the paper folder, so that a deposit sits with the study that motivated
> it, as every later deposit already did. Every required output still exists.

## Inside each series directory

Both have the same structure:

```
<SERIES>/
├── README.md              the generated report for this series (03_write_report.py)
├── figures/               umap/ dotplots/ featureplots/ composition/ qc/
├── tables/                markers, annotations, composition, cell metadata
├── qc/                    thresholds, before/after, doublets, batch assessment
├── logs/                  analysis_log.txt, decisions.json
├── inventory/             detected format, per-sample results
├── processed/             .h5ad checkpoints  (gitignored, regenerable)
└── epithelial_subanalysis/
```

`GSE262927/` additionally contains the focused analyses below.

## The mouse analyses differ by cohort, and this matters

| Directory (under `Thesis/gate1_01_niethamer_2025/GSE262927/`) | Cells | Cohort | Purpose |
|---|---:|---|---|
| top level | 162,175 | **all 33 samples** | Whole-atlas survey: 29 clusters, marker tables, composition. General-purpose object. |
| `regeneration_focus/` | 5,694 alveolar; 43,359 capillary | 25-sample annotated atlas | **The biology.** AT2 to Krt8-positive transitional to AT1 trajectory, and the persistent injury-induced capillary state. |
| `lineage_tracing_cohort/` | 36,958 | **the 8 non-atlas samples** | A separate experiment (Kit, Car4, Ednrb Cre; pre-labelled, 19 dpi). Where does the injury state come from? |
| `phase_timecourse/` | 107,626 | 25-sample annotated atlas | Per-day UMAP, per-animal lineage composition and Ki67-trace proliferation by lineage: the paper's three-phase structure, computed from the tracked metadata table only (script 10). Owner review pending. |
| `myeloid_focus/` | 9,997 | 25-sample annotated atlas | The myeloid compartment (atlas clusters 5, 17, 24) re-embedded with labels held out, graded against them, and read per animal and day (script 11). Two follow-ups live inside it: `amac_origin/` reads the Ki67 trace by tamoxifen window (script 12, tracked metadata only); `batch_sensitivity/` re-embeds the days that carry both infection rounds with and without Harmony on round (script 13). Owner review pending. |

**The whole-atlas object merges two experiments.** All 33 samples are in it,
which is fine for surveying cell types but wrong for anything condition- or
trace-related, because the 8 tracing samples use different Cre drivers and an
inverted labelling design. The focused directories restrict to the correct
cohort. See [`../docs/ANALYSIS_RATIONALE.md`](../docs/ANALYSIS_RATIONALE.md).

## Where to start

For the current shared scientific questions, use
[RESEARCH_QUESTIONS.md](../RESEARCH_QUESTIONS.md); the IL-1 context figures have
their own [methods and reproduction commands](figures/rq/il1b_context/REPORT.md).

| I want | Read |
|---|---|
| what was actually run, with parameters | [`../docs/PIPELINE_AS_RUN.md`](../docs/PIPELINE_AS_RUN.md) |
| why each decision was made | [`../docs/ANALYSIS_RATIONALE.md`](../docs/ANALYSIS_RATIONALE.md) |
| background (non-computational reader) | [`../docs/UMAP_AND_FIGURES.md`](../docs/UMAP_AND_FIGURES.md), [`../docs/BACKGROUND_FOR_BIOLOGISTS.md`](../docs/BACKGROUND_FOR_BIOLOGISTS.md), [`../docs/DOUBLETS_AND_SCRUBLET.md`](../docs/DOUBLETS_AND_SCRUBLET.md) |
| the mouse report | [`GSE262927/README.md`](../Thesis/gate1_01_niethamer_2025/GSE262927/README.md) |
| the human report | [`GSE178360/README.md`](../Thesis/ungated_murthy_2022/GSE178360/README.md) |

## Re-running

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
python analysis/scripts/05_write_pipeline_as_run.py
```
