# Progress and handoff state

Living record of what is done, what is pending, and what a future session needs
to know to continue. Update this before stopping.

Last updated: 2026-09-09. The scientific analysis of the two series is
complete (state of 2026-08-09 below): `FINDINGS.md` (results with figures)
leads, `README.md` is a landing page for the executed analysis, and
`scRNAseq_workflow_Niethamer2025.md` lives in `docs/`. On 2026-09-09 a
paper-by-paper roadmap directory was added (`Thesis/`, order taken from the
owner's Notion PI Target Map) with the Sikkema 2023 HLCA study note, its
decision criteria as reviewable JSON, a pipeline-framing proposal, and a
first criteria trial on tracked tables. **Owner retain/reject review of that
material is pending.** Nothing is running.

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
| `Thesis/` roadmap index (11 papers, Notion order, PubMed-verified IDs) | **DONE** |
| Sikkema 2023 (HLCA) study note, `integration_benchmark.json`, `PIPELINE_FRAMING.md` | **DONE, owner review pending** |
| Trial S1: HLCA cluster-entropy criteria on the tracked tables | **DONE (2026-09-09); artefacts in Thesis/gate1_04_sikkema_2023_hlca/trials/** |
| Trial S2: scArches mapping of GSE178360 to the HLCA core | **DONE (2026-09-09)**; 23 of 31 clusters agree with the blind proposals; AT0 is a minority and the AT0 candidate subcluster is mostly AT2 or uncertain; our mapping matches the HLCA authors' own transfer of the same cells at 99.2% (level 3) |
| Trial S3: HLCA consensus-marker annotation of GSE178360 | **DONE (2026-09-09)**; 17 to 19 of 31 clusters agree with the blind proposals; AT0 by marker transfer Not established (scheme-dependent) |
| Trial S4: mouse cluster 23 explained | **DONE (2026-09-09)**; low-count, ambient-like; 78% from EEM-scRNA-289 |
| Trial S5: mouse cluster 5 subclustered and re-graded | **DONE (2026-09-09)**; resolved at Leiden 0.5 (94% of labelled cells in pure subclusters), not at 0.2 |

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

11. **`Thesis/` is now a tracked directory** (2026-09-09). Its Markdown is
   link-checked by `validate_portfolio.py` like everything else. PDFs and
   spreadsheets inside it remain gitignored; only notes, JSON and small
   trial tables are tracked. Each new paper gets its own folder in the
   roadmap order (`Thesis/README.md`), one at a time.
12. **Mouse cluster 23 is a low-count, ambient-like barcode population, not a
   cell type** (trial S4, 2026-09-09): median 1,868 counts against 6,790
   atlas-wide, 84% of its cells detect three or more lineage markers, 78% of
   it comes from EEM-scRNA-289 (a Car4-CreERT2 tracing animal) where it is
   55% of the sample, and only 2% of its annotated-cohort cells carry an
   author label. Exclude it from composition and state claims and flag
   sample 289 in the lineage-tracing cohort. Adding this to the generated
   mouse report is pending the owner's decision.
13. **Human cluster 22 is mast cells, not B cells** (trial S3): both HLCA
   marker schemes call Mast cells and the cluster's own top genes are TPSB2,
   TPSAB1, CPA3 and KIT; the blind panel had no mast set. Correction of the
   human annotation table is pending the owner's decision.
14. **AT0 by HLCA marker transfer is scheme-dependent** (trial S3): the flat
   argmax calls AT0 in most distal secretory and AT2 clusters and passes the
   factor-of-two concordance with the strict gate; the hierarchical scheme
   collapses AT0 to a few percent of the gate and relabels the AT0 candidate
   analogue as AT2. Not established by this route; settled by trial S2
   (next item).
15. **The human AT0 headline needs re-wording** (trial S2, 2026-09-09).
   scArches mapping to the HLCA core, which reproduces the HLCA authors'
   own transfer of the same cells at 99.2% (level 3), labels 119 cells AT0
   with low uncertainty (0.4% of the series; 0.25% in the HLCA core) and
   calls the strict SFTPC+ SCGB3A2+ EPCAM+ gate population mostly pre-TB
   secretory and AT2. Epithelial subcluster 4, the "AT0 candidate
   analogue" (328 cells) in `FINDINGS.md` and the portfolio PDF, is 38%
   AT2, 36% uncertain and 23% AT0. Proposed wording: "an AT0-like minority
   exists; the candidate subcluster is not an AT0 population". Nothing has
   been edited; owner decision.
16. **HLCA label transfer confidently mislabels neutrophils** as classical
   monocytes (cluster 14: mode 0.94, mean uncertainty 0.056) because the
   reference has no neutrophil identity. Uncertainty does not catch an
   absent identity that resembles a present one. Apply an independent
   neutrophil check whenever HLCA labels are used.
17. **Human cluster 30 is a ciliated population private to one donor**
   (trial S2: Multiciliated, mode 0.98, uncertainty 0.02). Its identity is
   settled; it stays out of population claims because it is donor-private.

---

## Optional extensions (not completion gaps)

1. Add ambient-RNA correction if compatible unfiltered mouse droplet matrices
   become available; the current GEO deposit does not provide them.
2. Fit a formal trajectory-DE model if gene-level inference along pseudotime
   becomes a project goal; the current analysis intentionally stops at binned
   programme summaries.

---

3. Trial S2 has run (2026-09-09); its consequences (items 15 to 17 above)
   await the owner's retain/reject decisions. The x64 environment now also
   holds torch 2.14.0 (CPU) and scvi-tools 1.5.0.post1, frozen in
   `Thesis/gate1_04_sikkema_2023_hlca/trials/s2_reference_mapping/requirements_s2_env.txt`;
   `analysis/requirements.txt` is unchanged.
4. Add the next roadmap papers (Choi 2020, Nabhan 2018) as folders under
   `Thesis/` in order, each with its five-question note.

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
