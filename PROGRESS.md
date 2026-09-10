# Progress and handoff state

Living record of what is done, what is pending, and what a future session needs
to know to continue. Update this before stopping.

Last updated: 2026-09-10. The scientific analysis of the two series is
complete (state of 2026-08-09 below): `FINDINGS.md` (results with figures)
leads, `README.md` is a landing page for the executed analysis, and
`scRNAseq_workflow_Niethamer2025.md` lives in `docs/`. On 2026-09-09 a
paper-by-paper roadmap directory was added (`Thesis/`, order taken from the
owner's Notion PI Target Map) with the Sikkema 2023 HLCA study note, its
decision criteria as reviewable JSON, a pipeline-framing proposal, and a
first criteria trial on tracked tables. **Owner retain/reject review of that
material is pending.** On 2026-09-10 two focused analyses of the source
paper's phase structure were added under `analysis/GSE262927/`
(`phase_timecourse/` and `myeloid_focus/`, scripts 10 and 11; items 15 and
16 below). **Owner review of those is pending too.** Nothing is running.

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
| Trial S2: scArches mapping of GSE178360 to the HLCA core | **PLANNED**; environment Not established |
| Trial S3: HLCA consensus-marker annotation of GSE178360 | **DONE (2026-09-09)**; 17 to 19 of 31 clusters agree with the blind proposals; AT0 by marker transfer Not established (scheme-dependent) |
| Trial S4: mouse cluster 23 explained | **DONE (2026-09-09)**; low-count, ambient-like; 78% from EEM-scRNA-289 |
| Trial S5: mouse cluster 5 subclustered and re-graded | **DONE (2026-09-09)**; resolved at Leiden 0.5 (94% of labelled cells in pure subclusters), not at 0.2 |
| Phase-wise view of the Ki67 atlas (`analysis/GSE262927/phase_timecourse/`, script 10) | **DONE (2026-09-10), owner review pending**; per-dpi atlas UMAP, per-animal lineage composition, Ki67-trace proliferation by lineage; trace peaks fall in the paper's window for 4 of 5 lineages (Lymphoid peaks at 11 dpi, not 6); Descriptive only, from tracked metadata |
| Myeloid compartment by dpi (`analysis/GSE262927/myeloid_focus/`, script 11) | **DONE (2026-09-10), owner review pending**; 9,997 cells from atlas clusters 5, 17, 24; 16 blind subclusters at Leiden 0.5 with 87% of labelled cells in pure subclusters; aMAC loss and iMON expansion at 6 dpi with reconstitution by 19 to 42 dpi, consistent with the paper's Figure 3; Descriptive only |
| Alveolar macrophage origin by Ki67 trace window (`myeloid_focus/amac_origin/`, script 12) | **DONE (2026-09-10), owner review pending**; the 2 to 3 dpi window labels most of the 42 dpi aMAC pool (median 79.7%); marrow-inheritance and two-source checks Not established under the 30-cell floor |
| Batch sensitivity of the myeloid embedding, Harmony on infection round (`myeloid_focus/batch_sensitivity/`, script 13) | **DONE (2026-09-10), owner review pending**; rounds already mix within every tested day (enrichment 1.11 to 1.46, threshold 2); the 6 dpi iMON state survives correction as its own subcluster (92% of cells from 6 dpi, 77% and 74% of each animal's iMON cells); the frozen survival rule selected the wrong subcluster and the revision is disclosed |

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
python analysis/scripts/10_phase_timecourse.py
python analysis/scripts/11_myeloid_focus.py
python analysis/scripts/12_amac_trace_by_window.py
python analysis/scripts/13_myeloid_batch_sensitivity.py
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
   analogue as AT2. Not established by this route; trial S2 (reference
   mapping) is the route that can settle it.
15. **The paper's three proliferative phases are visible in the trace, not
   in the deposited cell-cycle call** (`phase_timecourse/`, 2026-09-10).
   With the expected peak windows frozen before the metadata table was
   opened, the median per-animal Ki67-traced fraction in each cohort's
   immediate window peaks at 6 dpi for myeloid cells, 11 dpi for epithelium
   and mesenchyme, and 19 dpi for endothelium (4 of 5 lineages in the
   paper's window; lymphoid cells peak at 11 dpi, one harvest later than
   the paper's immune window). The deposited Seurat cell-cycle call, carried
   as a cross-check, agrees for 1 of 5 lineages because it calls most
   lymphocytes and about 40% of all cells cycling; it is reported, not used.
   Two animals per active-repair day: a ranking, no test. Uses the deposited
   lineage labels descriptively. Owner decision pending.
16. **The myeloid compartment reproduces the paper's Figure 3 from a blind
   embedding** (`myeloid_focus/`, 2026-09-10). Atlas clusters 5, 17 and 24
   in the 25-sample cohort (9,997 cells, 13 with a non-myeloid label) were
   re-embedded with the labels held out; at Leiden 0.5 (fixed from trial S5
   before the run) 87% of labelled cells sit in subclusters of purity at
   least 0.75 (0.2: 48%; 1.0: 92%). Median per-animal aMAC share of myeloid
   cells falls from 30.8% at baseline to 4.6% at 6 dpi and returns to 41.7%
   by 19 dpi and 49.0% by 42 dpi; iMON rises from 1.9% to 56.0% at 6 dpi and
   is back to 1.8% by 42 dpi. Four subclusters remain label mixtures (7, 9,
   12, 13; 1,286 cells). The MACS recombination fixes the myeloid share of
   each library, so only within-myeloid fractions are read. Owner decision
   pending.
17. **The rebuilt alveolar macrophage pool is labelled mainly by the 2 to 3
   dpi window** (`myeloid_focus/amac_origin/`, 2026-09-10). At the common
   42 dpi harvest the median per-animal Ki67-traced fraction of
   aMAC-labelled cells is 79.7% for the 2 to 3 dpi window, 60.0% for 7 to 8,
   29.3% for 14 to 15 and 42.0% for 21 to 22, so the early window wins by the
   frozen rule, as the paper's Figure 3 implies. The two pre-registered
   checks fail closed: the marrow-inheritance reference (cMON and
   neutrophils of the same animal) reached the 30-cell floor in only 3 of 8
   animals, and the within-aMAC subcluster split is evaluable in one animal
   per window (where the early windows label subcluster 2 more than
   subcluster 3 by 28 to 33 points and the late windows show the reverse).
   Both are Not established; the per-animal values are in the tables. Owner
   decision pending.
18. **The 6 dpi inflammatory-monocyte state is not a batch island**
   (`myeloid_focus/batch_sensitivity/`, 2026-09-10). On every active-repair
   day and at 42 dpi the two animals come from different infection rounds
   (Table S3; the rounds also differ in Ki67-Cre dosage), so round is the
   technical key that crosses time; 0, 90 and 366 dpi are single-round and
   untestable. Within each tested day the rounds already mix in the
   uncorrected embedding (same-round kNN enrichment 1.11 to 1.46 against
   the pipeline's failure threshold of 2.0; 1.05 to 1.29 after Harmony).
   After Harmony on round the 6 dpi iMON state remains its own subcluster
   (411 cells, 92% from 6 dpi, iMON purity 0.86, holding 77% and 74% of the
   two animals' iMON cells; uncorrected: 530 cells, 91%, 0.89, 90% and 97%).
   Adjusted Rand index between the partitions 0.837. **Rule revision
   disclosed:** the frozen definition of "the iMON state" (subcluster with
   the most iMON-labelled cells) selected the 11 to 19 dpi monocyte state in
   both embeddings; the first-run outcome is kept in the run record and a
   post hoc definition anchored on the 6 dpi cells is reported alongside
   (DEVELOPMENT decision 15). Owner decision pending.

---

## Optional extensions (not completion gaps)

1. Add ambient-RNA correction if compatible unfiltered mouse droplet matrices
   become available; the current GEO deposit does not provide them.
2. Fit a formal trajectory-DE model if gene-level inference along pseudotime
   becomes a project goal; the current analysis intentionally stops at binned
   programme summaries.

---

3. Run trial S2 (map GSE178360 to the HLCA core with scArches, transfer
   labels with the 0.2 uncertainty cutoff, compare with the blind and
   deposited annotations, AT0 check). Blocked on environment work
   (PyTorch and scArches in the emulated interpreter); the plan and the
   pre-registered success criteria are in
   `Thesis/gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md`.
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
