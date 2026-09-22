# Progress and handoff state

Living record of what is done, what is pending, and what a future session needs
to know to continue. Update this before stopping.

Last updated: 2026-09-21 (the two original series moved into Thesis/ beside their source papers, the README left open-ended, the owner's rule on multi-agent cost recorded, the docs pages cleared of em-dashes, a paper-style figure drawn for each research question, a scan of where a gene set enrichment analysis is admissible, and the first three GSEA trials; see items 37 to 40). Previously 2026-09-20 (two branches of roadmap paper 2 on multiome deposits, an eight-adversary audit that corrected 14 register rows, the root README rewritten around a claims ledger drawn from the register, and RESEARCH_QUESTIONS.md added as the question-first entry point; items 35 and 36), 2026-09-17, 2026-09-15, 2026-09-13, 2026-09-12 and 2026-09-10. The scientific analysis of the two series is
complete (state of 2026-08-09 below): `FINDINGS.md` (results with figures)
leads, `README.md` is a landing page for the executed analysis, and
`scRNAseq_workflow_Niethamer2025.md` lives in `docs/`. On 2026-09-09 a
paper-by-paper roadmap directory was added (`Thesis/`, in the owner's
reading order) with the Sikkema 2023 HLCA study note, its
decision criteria as reviewable JSON, a pipeline-framing proposal, and a
first criteria trial on tracked tables. **Owner retain/reject review of that
material is pending.** On 2026-09-10 two focused analyses of the source
paper's phase structure were added under `Thesis/gate1_01_niethamer_2025/GSE262927/`
(`phase_timecourse/` and `myeloid_focus/`, scripts 10 to 13; items 18 to
21 below). **Owner review of those is pending too.** Later on 2026-09-10 the
S2 trial was merged (PR #8) with the items renumbered chronologically, and
the repository was reframed as an analysis log: the August 2026 summary PDF and the
Krt8-high transitional narrative were displaced to `archive/`, the README
now opens with the claims table, and the injury model is named once (item 22
below). Nothing is running.

---

## Status at a glance

| Item | State |
|---|---|
| Raw-data scan + inventory | **DONE** |
| GSE262927 (mouse) full analysis | **DONE**, 162,175 cells, 33 samples, 29 clusters |
| GSE178360 (human) full analysis | **DONE**, 27,729 cells, 3 samples, 31 clusters, Harmony primary |
| Reports (`README.md` per dataset, analysis logs) | **DONE** |
| `docs/PIPELINE_AS_RUN.md` (what actually ran) | **DONE** |
| Source papers read and divergences documented | **DONE** |
| Doc scope labels (README/WORKFLOW/docs) | **DONE**, banners on all 5 tool pages, README/WORKFLOW relabelled |
| `docs/ANALYSIS_RATIONALE.md` (decisions, before vs after papers) | **DONE** |
| `docs/BACKGROUND_FOR_BIOLOGISTS.md`, Harmony | **DONE** |
| `docs/UMAP_AND_FIGURES.md` | **DONE** |
| `docs/DOUBLETS_AND_SCRUBLET.md` | **DONE** |
| **Focused alveolar regeneration re-analysis** | **DONE**, AT2→transitional→AT1 trajectory recovered |
| Capillary endothelial (iCAP) sub-analysis | **DONE**, persistent injury state recovered |
| AT0 doublet-loss follow-up | **DONE**, tested, not lost |
| Lineage-tracing cohort (8 non-atlas samples) | **DONE**, CAP1 origin supported; CAP2 lines uninformative |
| `analysis/` reorganised into two series subdirectories | **DONE** |
| Pipeline speedups (PCA, scan, threads) | **DONE**, verified output-identical |
| Repository restructure (`FINDINGS.md`, README rewrite, root tidy-up) | **DONE** |
| `Thesis/` roadmap index (11 papers, the owner's reading order, PubMed-verified IDs) | **DONE** |
| Sikkema 2023 (HLCA) study note, `integration_benchmark.json`, `PIPELINE_FRAMING.md` | **DONE, owner review pending** |
| Trial S1: HLCA cluster-entropy criteria on the tracked tables | **DONE (2026-09-09); artefacts in Thesis/gate1_04_sikkema_2023_hlca/trials/** |
| Trial S2: scArches mapping of GSE178360 to the HLCA core | **DONE (2026-09-09)**; 23 of 31 clusters agree with the blind proposals; AT0 is a minority and the AT0 candidate subcluster is mostly AT2 or uncertain; our mapping matches the HLCA authors' own transfer of the same cells at 99.2% (level 3) |
| Trial S3: HLCA consensus-marker annotation of GSE178360 | **DONE (2026-09-09)**; 17 to 19 of 31 clusters agree with the blind proposals; AT0 by marker transfer Not established (scheme-dependent) |
| Trial S4: mouse cluster 23 explained | **DONE (2026-09-09)**; low-count, ambient-like; 78% from EEM-scRNA-289 |
| Trial S5: mouse cluster 5 subclustered and re-graded | **DONE (2026-09-09)**; resolved at Leiden 0.5 (94% of labelled cells in pure subclusters), not at 0.2 |
| Phase-wise view of the Ki67 atlas (`Thesis/gate1_01_niethamer_2025/GSE262927/phase_timecourse/`, script 10) | **DONE (2026-09-10), owner review pending**; per-dpi atlas UMAP, per-animal lineage composition, Ki67-trace proliferation by lineage; trace peaks fall in the paper's window for 4 of 5 lineages (Lymphoid peaks at 11 dpi, not 6); Descriptive only, from tracked metadata |
| Myeloid compartment by dpi (`Thesis/gate1_01_niethamer_2025/GSE262927/myeloid_focus/`, script 11) | **DONE (2026-09-10), owner review pending**; 9,997 cells from atlas clusters 5, 17, 24; 16 blind subclusters at Leiden 0.5 with 87% of labelled cells in pure subclusters; aMAC loss and iMON expansion at 6 dpi with reconstitution by 19 to 42 dpi, consistent with the paper's Figure 3; Descriptive only |
| Alveolar macrophage origin by Ki67 trace window (`myeloid_focus/amac_origin/`, script 12) | **DONE (2026-09-10), owner review pending**; the 2 to 3 dpi window labels most of the 42 dpi aMAC pool (median 79.7%); marrow-inheritance and two-source checks Not established under the 30-cell floor |
| Batch sensitivity of the myeloid embedding, Harmony on infection round (`myeloid_focus/batch_sensitivity/`, script 13) | **DONE (2026-09-10), owner review pending**; rounds already mix within every tested day (enrichment 1.11 to 1.46, threshold 2); the 6 dpi iMON state survives correction as its own subcluster (92% of cells from 6 dpi, 77% and 74% of each animal's iMON cells); the frozen survival rule selected the wrong subcluster and the revision is disclosed |
| Trial S2 merged into main; PROGRESS and DEVELOPMENT renumbered chronologically | **DONE (2026-09-10, PR #8)** |
| Repository hygiene and reframing: `archive/` for displaced material, README opens with the claims table, injury model named once, checks renamed | **DONE (2026-09-10), owner review pending**; see item 22 |
| Cardoso 2026 (Gate 2, paper 5): study note, extracts, trials C0 to C12 and E1 to E6 (`Thesis/gate2_05_cardoso_2026/`) | **DONE (2026-09-12 to 2026-09-13), owner review pending**; Gate 1 not recovered; list A exhausted; rows C19 to C57 and C65 to C84; items 23 to 25 and 27 to 29 |
| Choi 2020 (roadmap paper 2): study note, extracts, trials D0 to D7 with corrected passes D2b and D5b (`Thesis/gate1_02_choi_2020/`) | **DONE (2026-09-15), owner review pending**; entered at the owner's direction after reading (the AI-written note of 2026-09-13 was withdrawn first, items 26 and 31); four of five states, time course, ordering, programmes and the organoid shift reproduce as descriptions; primed AT2 never a cluster; rows C58 to C64 and C85 to C104; item 33 |
| Reading order revised: Gate 2 branches 2C, 2N, 2W; papers 12 to 16 and methods references M1 to M7 added (`Thesis/README.md`, `ROADMAP.json`, `REFERENCES.md`) | **DONE (2026-09-15)**; documentation only; item 31 |
| Generated `NEGATIVE_RESULTS.md`, one validated palette, shared deposit readers | **DONE (2026-09-13)** |
| Roadmap re-ranked: Nabhan papers next, Gate 3 paused (`Thesis/README.md`, `ROADMAP.json`) | **DONE (2026-09-15, PR #26)**; documentation only; see item 30 |
| Branch labels (2C, 2N, 2W, 3A, 3B) carried into every document that names a gate or a branch; the two gate namespaces stated | **DONE (2026-09-15, PR #28)**; documentation only; item 32 |

Repository state: the scientific analysis and its curation are complete;
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
`sys.stdout.reconfigure(encoding="utf-8")` because printing `, ` otherwise
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
3. **Scrublet bias, TESTED AND LARGELY EXONERATED.** A crude co-expression
   gate suggested up to 2× over-removal of the human paper's novel populations.
   Re-measured with a gate that can separate real AT0 from AT2+club doublets
   (`SFTPC+ SCGB3A2+ EPCAM+`, lineage-negative), AT0 is flagged at **3.9% vs a
   6.3% baseline**, below background. The cells flagged inside the crude gate
   have 1.8× the UMIs and more cross-lineage co-expression, i.e. they are
   doublets. Lesson: a co-expression gate cannot audit a co-expression
   detector. General caution still applies; the verdict for this dataset does
   not.
4. **No ambient-RNA correction anywhere.** Both papers used SoupX. Matters most
   for AT0 (defined partly by SFTPC, a classic soup contaminant) and for
   GSE178360 DD046Q (heavy haemoglobin).
5. **Doublet calls are a ranking cut, not a detected threshold** in 32/33 mouse
   samples (non-bimodal Scrublet histograms → expected-rate quantile fallback).
6. **MAD upper bounds never bind**, zero cells removed for excess counts/genes.
7. **The papers' headline findings need subsetting.** *(ADDRESSED:
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
   link-checked by `validate_repository.py` like everything else. PDFs and
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
   analogue" (328 cells) in `FINDINGS.md` and the August 2026 summary PDF, is 38%
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
18. **The paper's three proliferative phases are visible in the trace, not
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
19. **The myeloid compartment reproduces the paper's Figure 3 from a blind
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
20. **The rebuilt alveolar macrophage pool is labelled mainly by the 2 to 3
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
21. **The 6 dpi inflammatory-monocyte state is not a batch island**
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
   (DEVELOPMENT decision 16). Owner decision pending.
22. **Established material is displaced, not deleted** (2026-09-10, owner
   instruction). The August 2026 summary PDF and its generator moved to
   `archive/portfolio_2026-08/`; the Krt8-high transitional narrative
   (FINDINGS section 1, the human KRT8 reference-aligned panels) moved to
   `archive/DISPLACED.md` because it is established outside this
   repository. Their artefacts and scripts stay in place
   (then under `analysis/`; the artefacts moved with their series to `Thesis/`
   on 2026-09-21, item 37) because the validator checks their numbers and the
   scripts regenerate them; only the narrative moved. The CI workflow and validator were
   renamed (`repository-checks.yml`, `validate_repository.py`), the paper's
   own workflow document moved to `docs/WORKFLOW_Niethamer2025.md`, and the
   README now opens with the claims table and names H1N1 once as the injury
   model. Open decision for the owner: whether the alveolar-trajectory
   artefacts (`regeneration_focus/` alveolar figures and tables, script 06's
   alveolar branch, script 08 and its reference-aligned figures) should also
   be relocated physically, which would require re-pointing the generated
   human report and the validator's transitional-abundance checks.

23. **Cardoso 2026 (Gate 2, paper 5) entered out of order; all eight trials
   run** (`Thesis/gate2_05_cardoso_2026/`, 2026-09-12, owner instruction).
   Gate 1 returned a negative result and stopped to characterise, as the gate
   design requires; Gate 2 branches (a) and (b) then ran and both returned
   positive results, while branch (d) is closed by the data.
   - **C0, the data reality check.** All five accessions are public and
     readable: GSE316241, GSE316243, GSE316244, GSE310335, and GSE247505
     (England et al. 2025, the companion paper whose epithelial cells the
     paper's CellChat analysis uses and which the Cardoso data-availability
     statement does not name). Thirty libraries, 123,807 barcodes, three
     distinct gene spaces, one non-gene feature (BSD, the reporter
     construct's selection marker, detected in about 39% of RFP-sorted cells
     against 3 to 5% of niche cells).
   - **The binding constraint.** Every deposited mouse library of this paper
     pools three mice and each genotype contributes one library per sort, so
     **no genotype contrast in the deposit has within-group replication**.
     Directions are describable; nothing is testable; no P value is
     admissible. The exception is GSE247505, with two libraries per arm at
     three time points and a within-animal wild-type control (YFP).
   - **Gate 2 branch (d) is closed by the data, not by choice.** Every
     mesenchymal and immune library is a single time point (2 weeks), so the
     paper's imaging claim that fibroblast reprogramming precedes macrophage
     change cannot be tested transcriptomically. Not run; the reason is the
     result.
   - **C1, Gate 1: the fibrotic fibroblast subset is NOT recovered** by the
     pre-registered three-part rule (at least 80% Red2Kras, highest
     reprogrammed score, Tnc detected in at least 40% of cells). Two
     disclosed reasons, both recorded rather than repaired in place: the
     frozen selection step took only clusters with a *confident* fibroblast
     call, which excluded every cluster the caller labelled reprogrammed
     fibroblast and also cluster 14 (99.6% Red2Kras); and the score criterion
     is dominated by Acta2 and Pdgfrb, which mural cells express, so the
     highest-scoring subcluster is smooth muscle. Also recorded: the paper's
     own quality thresholds would remove 29% of the Red2Kras library against
     16% of the Confetti library, and Scrublet's automatic threshold failed
     on both libraries (0.07% and 0.00% called against a 5% prior).
   - **C1b, stop-and-characterise.** No cluster matches the published
     signature under the frozen criteria, but **cluster 14 misses the Tnc
     floor by one thousandth** (0.399 against 0.400) while carrying 937
     cells, 99.6% Red2Kras, Acta2 0.574, Runx1 0.540, Pdgfrb 0.326. The
     threshold was not moved after the fact. Five Red2Kras-private clusters
     (4, 10, 11, 14, 16) are **not** explained by low counts or doublets,
     unlike mouse cluster 23 of the GSE262927 series. Owner decision needed
     on whether cluster 14 is to be called the published population.
   - **Open question raised by C1b, addressed by C1c (written, not run).**
     The Red2Kras-private clusters carry fibrotic markers (Tnc, Runx1, Acta2)
     and inflammatory markers (Lcn2, Saa3) together, whereas the paper makes
     those distinct populations and says the inflammatory cells lack Tnc and
     appear only from 4 weeks. C1c asks the per-cell co-detection question a
     cluster-level fraction cannot answer.
   - **C1c and C1d answered the open question and found a better one.**
     Fibrotic and inflammatory markers mark separate cells at 2 weeks (9.7%
     double-positive against 8.2% expected under independence), which supports
     the paper's distinctness claim. C1d then explained the tumour-private
     clusters by composition rather than quality: 6.5% of the mesenchymal
     library reads as off-target for its own sort, including 184 mutant
     epithelial cells that are 88% Areg-positive. The paper is unaffected, but
     a reanalysis computing signalling inside that library alone would be.
   - **C2, Gate 2b: the paper's epithelial result reproduces blind.**
     DATP-like cells fall from 46.6% to 22.0% of RFP+ cells on Areg deletion
     and AT2 rises from 29.2% to 59.1%, against the paper's 50.1 to 25.9 and
     21.3 to 54.6. Four of five pre-registered directions are met (the fifth
     unscorable), the DATP-like cluster is depleted 15.4-fold and alveolar
     macrophages 9.4-fold. **What fails to collapse:** of the six genes in the
     paper's fibrotic set, Tnc, Acta2, Fst and Runx2 fall while Pdgfrb and
     Runx1 do not. C2b recomputes the two shares the confidence floor hid,
     the same defect C1 disclosed.
   - **C3, Gate 2a: every frozen test holds.** Areg is higher in the DATP-like
     state than in AT2 cells in all four mutant libraries with the state
     defined by our own clustering; the ligand order Areg > Hbegf > Ereg >
     Tgfa is identical in all four, so **Hbegf ranks second**, ahead of the
     ligand the paper followed into culture; the state is nearly absent from
     wild-type clones of the same animals (0.08% and 1.1%). Replicate
     libraries mix within every arm (1.53 to 1.80 against a threshold of 2.0),
     so no batch correction was applied, by the rule rather than by
     preference. **It is deliberately not a CellChat rerun**: CellChat is
     R-only and this machine has no R, so the trial re-derives the expression
     fact the communication claim rests on and states what it is not.
   - **C5, 2026-09-13: three corrections.** Hbegf is second by abundance but third by enrichment over AT2 in three of four libraries; the 4-day replicates disagree 24-fold on the DATP-like share; the fibrotic response is three graded tiers, not a split. Rationale and next steps: `Thesis/gate2_05_cardoso_2026/DIVERGENCES_AND_NEXT.md`.
   - Status of every row: Descriptive only or Exploratory at best, by the
     replication constraint. Owner retain/reject review pending on all of it.
     The two decisions that matter: whether cluster 14 is the published
     population (claim C22), and whether the Hbegf lead (C33) is worth
     pursuing.

24. **The Hbegf lead was extended into four public datasets, and one of the
   repository's own leads was refuted in the process**
   (`Thesis/gate2_05_cardoso_2026/trials/e*.py`, 2026-09-13, owner
   instruction). Six trials, E1 to E4 plus E1b and E2b, all outside the
   Cardoso deposit. Full account in the newest handoff block below.
   - **Why leave the deposit.** Item 23's binding constraint means no Cardoso
     genotype contrast can be tested. GSE131907 has eleven donors with paired
     tumour and normal lung, so the unit becomes the donor and a test becomes
     admissible. AREG detection is higher in epithelium than in myeloid cells
     within donor, 0.336 against 0.215, paired Wilcoxon p = 0.0020. That is
     claim C37, and the only Validated row the Cardoso work has produced.
   - **A frozen rule refuted the repository's best lead.** Trial E4 found
     Runx1 and Pdgfrb both rise with bleomycin alone in sorted Col1a1-GFP
     mesenchyme with no oncogene present. The reading fixed before the
     matrices were opened therefore applies: the Areg-independent tier of
     claim C29 is what an activated lung fibroblast does after injury, not a
     second tumour signal. C29 keeps its numbers and loses its
     interpretation; Hbegf and Egfr are injury-generic as well.
   - **The Hbegf lead does not transfer across species.** In human lung the
     ligand is myeloid-dominant, alveolar macrophages 0.712 and two
     dendritic-cell subsets above 0.55, agreeing with the published human work
     and disagreeing with what trial C6 found in mouse. The receptor behaves
     differently again: mesenchymal in fibrosis as in the mouse, epithelial in
     adenocarcinoma. Placement tracks the disease, the ligand source tracked
     the species, and neither direction is clean enough for a cross-species
     argument.
   - **Three things that could not be computed, reported rather than
     approximated.** E1b's within-tumour comparison has zero donor pairs,
     because the GSE131907 annotation assigns AT2 only in normal lung; the
     Zhao prediction clears no donor floor in either fibrosis cohort and the
     cohorts disagree in direction; and E2's only testable comparison failed
     to detect its difference at seven donors, p = 0.297.
   - **One correction to an earlier output.** E1's top HBEGF compartment is
     recorded as "neutrophil" and that gate is 83 per cent deposited myeloid
     cells in a deposit annotating no neutrophils, so the outcome is the
     pre-named myeloid-dominant one. Both readings are in the trial log.
   - Rows C37 to C48 in the claims register; DEVELOPMENT decision 19. Owner
     retain/reject review pending on all of it.

25. **Three chosen follow-ups run; list A is now empty and every answer is
   negative or a correction** (`Thesis/gate2_05_cardoso_2026/trials/e6_*.py`,
   `c7_*.py`, `c8_*.py`, 2026-09-13, owner selection). Full account in the
   newest handoff block below.
   - **E6, the first direct test of the paper's axis.** Across 22 donors of
     GSE136831, epithelial AREG against fibroblast EGFR gives Spearman rho
     0.348 at p = 0.112, and against a fibroblast activation score rho -0.150.
     Within the 18 fibrosis donors, 0.276 and -0.013. No donor-level coupling
     at any strength this cohort could detect, which is about rho 0.43 and no
     smaller. It does not contradict the paper: pooling all epithelium dilutes
     the state the claim is about, the same limit E1b measured, and a
     state-resolved version is blocked because that state clears the cell
     floor in only seven donors (E2).
   - **The depth control is the result worth reusing.** E6's only significant
     correlation was a control pair, epithelial TGFA against fibroblast
     activation at p = 0.045, and both members track sequencing depth, so the
     frozen rule refused to read it. The same number on AREG would have looked
     like confirmation. Epithelial AREG itself is clean against depth at 0.124.
   - **C7 closes item A1: the sort contaminant is a mixture, not a state.** Of
     its 184 cells, 77 score DATP-like, 70 AT2, 20 cycling, 13 AT1-like and 4
     Cd177-positive, a modal call of 41.8 per cent that sits below this
     repository's own 50 per cent floor. So about four in ten of the
     contaminating cells are the paper's own signalling population and about
     four in ten are ordinary AT2. Claim C26's caution holds in that exact
     form; the contaminant cannot be called the DATP-like state.
   - **C7 also refuted its own primary measure.** Whole-profile rank
     correlation separated the epithelial query from a fibroblast control by
     0.415 of rho and separated epithelial states from each other by 0.0056,
     with the top four matches spanning three different calls inside 0.015. It
     has compartment resolution and no state resolution. The rule named a
     winner because it was written without a margin requirement; the margin it
     did report is what refuses the answer, and the threshold was not moved.
   - **C8 closes item A2: the tiers are amplitudes.** In the deletion arm Runx1
     and Pdgfrb co-occur in 22.46 per cent of gated fibroblasts against 20.44
     expected under independence, ratio 1.099 inside the frozen band, with the
     depth control moving it by 0.043. The mixture test that disagreed is
     uninformative here, because 31.5 per cent of the cells detect neither gene
     so the score has a spike at zero; trial C1 used the same test on a centred
     score where that does not arise, so the defect is the reuse. Recorded, not
     repaired.
   - **The one surviving lead, post hoc and labelled.** The co-organisation
     sits on the falling tier: Fst with Runx2 at ratio 1.735 in Areg-flox/+ and
     0.816 after deletion, where only 24 of 3,700 cells carry both. Runx1 with
     Pdgfrb is inside the independence band in both arms. Read with E4, Areg
     deletion looks like it removes a co-expressing Fst and Runx2 population
     while the injury-generic genes stay spread across the compartment. One
     library per genotype, 24 cells, so this is a target for a future
     pre-registration and not a claim.
   - Rows C49 to C57 in the claims register; DEVELOPMENT decision 20. Owner
     retain/reject review pending on all of it.

26. **Roadmap paper 2 (Choi 2020) entered, study note written, Gate 0 run**
   (2026-09-13; the attribution is corrected and the folder withdrawn on
   2026-09-15, see item 31 and DEVELOPMENT decision 21). This is the paper that defines the DATP state the Cardoso
   rows lean on.
   - **The study note and extract are written**, in the roadmap's five-question
     format, with every marker set read from the Europe PMC full-text XML
     rather than the PMC web rendering. That rendering strips italicised gene
     symbols, so a note built from it would have had empty marker lists. The
     DATP set is Cldn4, Krt8, Ndrg1, Sprr1a and AW112010, with the negative
     condition of low Pdpn, Hopx and Cav1.
   - **D0, the data reality check, found the same ceiling as the Cardoso
     deposit.** One library per condition throughout, so no contrast carries
     within-group replication: the frozen replicate rule was applied to three
     contrasts and none passed.
   - **Six of the eight matrices are raw 10x barcode whitelists**, not called
     cells. Cell calling for the in vivo half is this repository's job, and the
     cell count will not match the paper's, because Cell Ranger 2.0.2's caller
     will not be reproduced. Gate 0 deliberately did not call cells.
   - **The tdTomato reporter is not a counted feature**, so the lineage split
     cannot be checked from the matrix and a reanalysis has to trust the
     library labels. That is the opposite of the Cardoso deposit, where the BSD
     selection marker gave an independent sort check.
   - **One gene space across all eight libraries** (27,998 features), and every
     marker gene including the negative conditions is present in all of them.
   - **The ATAC-seq accession is bigwig coverage only**, so the epigenetic half
     of the Il1r1-subset claim cannot be re-derived from the deposit at all.
   - **A depth warning for Gate 1:** median genes per barcode spans 891 to
     1,825 across the six libraries any cross-library comparison would use.
   - One reporting defect is disclosed: the first run's headline collapsed the
     raw-versus-filtered question into a single boolean across eight
     libraries, which read as "none are raw" when six of them are. The table
     was correct and no number changed; the summary now names both groups.
   - Rows C58 to C64 in the claims register, a new Stage 4 section. Owner
     retain/reject review pending on all of it.

27. **The Fst and Runx2 lead, pre-registered and run** (`c9_the_fst_runx2_population.py`,
   2026-09-13, owner instruction). The trial separated an answerable question
   from an unanswerable one before it ran, and the unanswerable one is stated
   rather than approximated: whether Areg deletion depletes the population
   cannot be tested, because the deposit has one library per genotype and no
   other Areg-flox fibroblast dataset exists.
   - **What replicates.** In both bleomycin animals of GSE132771, with no
     oncogene anywhere, Fst and Runx2 co-occur above chance against a null that
     permutes within deciles of genes per cell, p = 0.005 and 0.010. Depth is
     held fixed by construction, so this is not a depth artefact.
   - **What does not.** The size. Ratios of 2.84 and 2.32 move by 1.10 and 0.45
     across a median depth split against a frozen limit of 0.25, on 32 and 37
     double-positive cells. About thirty cells cannot support a stable ratio,
     and the control cannot separate depth dependence from sampling noise there.
   - **The marker signature replicates in direction**, 26 of the 27 genes that
     could be checked, in both animals independently. Three could not be
     checked at all, because Cemip2, AI506816 and 1110038B12Rik are missing
     from the older annotation GSE132771 was aligned against.
   - **The magnitude criterion was badly designed and its figure should not be
     quoted.** Nine of 27 clear "at least half the reference magnitude", but
     the median gene sits at 0.535 of its reference difference against a 0.5
     threshold, and the bleomycin libraries carry half the genes per cell. The
     criterion measured depth.
   - **What the cells look like.** Piezo2, Ltbp2, P4ha3, Sdc1 and Prrx2, with
     Bmper, Ror2 and Kif26b in the wider set: a matrix and mechanosensor-gene
     signature that points away from the Areg axis rather than into it.
     Exploratory. The wording was tightened on 2026-09-13 from
     "mechanically responsive", which transcript detection cannot support.
   - **The third disclosed rule defect in this folder.** T1 tested a ratio that
     T5 forbade reading, so the rules were not composable; the implementation
     took the conservative branch and its verdict stands, but the reading
     attached to it is wrong. Rows C54, C56 and now C72.
   - Rows C65 to C72 in the claims register. Owner retain/reject review pending.

28. **The lead closes, and every contradiction now has a figure**
   (`c10_published_state_or_not.py`, `c11_figures_for_the_contradictions.py`,
   2026-09-13, owner instruction).
   - **C10 settles what C9 could not.** The Fst and Runx2 double-positive
     fibroblasts are the published Cthrc1-positive pathological fibroblast of
     Tsukui et al. 2020. Standardised difference on that set is 1.81 in the
     reference library and 1.49 and 1.46 in two independent bleomycin animals,
     and it is the largest of the four sets scored everywhere. Alveolar
     identity is strongly negative and smooth muscle is flat.
   - **Cthrc1 was hidden by a ranking cutoff.** It is detected in 44.9 per cent
     of double-positive against 2.0 per cent of double-negative cells, a
     difference of 0.428 that fell below C9's top-thirty cutoff of 0.510.
     Col3a1 is formally unrankable at 1.000 against 0.987.
   - **So the lead stops here.** This repository re-derived a known fibroblast
     state from a tumour deposit by an unusual route, which is a method note
     rather than a finding, and it explains E4: a pathological fibroblast is by
     definition an injury state. Claim C70's reading is superseded.
   - **One caveat that does not rescue it.** The double-positive call is partly
     a depth call, since 342 of 379 such cells sit in the deep half. Depth
     cannot raise the pathological score while lowering the alveolar score, so
     the opposite signs rule it out.
   - **The fourth disclosed rule defect.** T5 covered a computed half that
     moves and said nothing about a half that cannot be computed; the bleomycin
     shallow halves hold one and four double-positive cells. The first run
     treated not computable as failed and suppressed a reading the reference
     library supports. Code corrected to match the written rule, threshold
     unmoved. Rows C54, C56, C72 and now C79.
   - **C11 draws the six contradiction results** from tracked tables only, each
     with a fail-closed assertion. Drawing them corrected three more
     overstatements, including a panel title claiming every epithelial state
     scored the same and a subtitle putting mouse myeloid cells at the bottom
     when neutrophils rank third.
   - Rows C73 to C79 in the claims register. Owner retain/reject review pending.

29. **The CellChatDB scan, and the guard that made it readable by refusing to
   read it** (`c12_cellchatdb_full_resource_scan.py`, 2026-09-13, owner
   instruction, as the substitute for installing R).
   - **What it is.** CellChat's own curated resource and scoring logic through
     liana, on GSE136831, 26 donors and 313 pairs, with the permutation test
     switched off because its unit is the cell. Not a CellChat rerun, and
     CellChat itself has still never run here.
   - **Why not the mouse deposit**, stated before the run: different gene
     spaces between the series the paper integrates, ligand and receptor from
     different libraries at different depths, one library of three pooled mice
     per genotype, and the 184 Areg-high contaminants of C1d would have formed
     the epithelial side.
   - **The guard fired.** Eleven of the top fifteen pairs target CD44, with
     collagen, laminin and fibronectin ligands, so the full-resource ranking
     reports transcript abundance and T1 and T2 are not read. That is the
     methods result the trial produced, and it explains why naming a shortlist
     in advance, as the paper did, is defensible.
   - **The one clean comparison.** Holding the receptor constant removes the
     promiscuity effect. Among the five EGFR ligands, AREG is first at median
     rank 15.5 of 313 across 26 donors, HBEGF second at 48.5, then TGFA, EREG
     and BTC. The top two match the mouse order from C3; EREG and TGFA swap.
   - **A fifth rule defect, disclosed.** The guard's own text scopes it to the
     two readings it names, while the reading gloss in the same docstring says
     nothing else is read if the guard fires. T4 is reported as provisional
     with the conflict stated. Reading it vindicates the paper rather than any
     hypothesis of ours, which is worth noting given the conflict.
   - Rows C80 to C84 in the claims register. Owner retain/reject review pending.
30. **Roadmap re-ranked (2026-09-15).** The owner re-ranked the reading
   order; at the time this repository had trials built for branch 2W and
   none for branch 2N. Consequences recorded in
   `Thesis/README.md` (gate rules and the status column) and
   `Thesis/ROADMAP.json` (a `reranking_2026_09_15` block and per-paper
   statuses): papers 3 and 6 (Nabhan 2018 and 2023) are next, with proposal
   Nb1 as their trial; Gate 3A and 3B are paused, not closed; the analysis
   contract and the KRT8 pilot remain unwritten. Two pieces of document
   drift were fixed in the same change: `AI_CONTEXT.md` still listed Choi
   2020 as next and Cardoso as trials C0 to C3, and the status table above
   had no Cardoso or Choi row. No analysis ran and no claim status changed.
   - The one owner decision left open by the re-ranking concerns private
     planning and is tracked in the owner's private notes; the roadmap keeps
     Gate 3A paused either way.
31. **Reading order revised, the Choi 2020 folder withdrawn, decision 21
   corrected (2026-09-15, owner instruction).** Three things in one change.
   - **Decision 21 corrected.** The record had said the owner directed a
     return to roadmap order on 2026-09-13. The session transcript shows a
     one-word "proceed" against an agent-written list that carried Choi 2020
     as the agent's own recommendation, and later that day an explicit
     instruction to leave the paper until the owner had read it. DEVELOPMENT
     decision 21 and its responsibility row now say so, and the AI-written
     study note is recorded as REJECTED.
   - **The Choi 2020 folder is withdrawn.** `Thesis/gate1_02_choi_2020/` no
     longer exists on main. The study note and the D-series plan stay in git
     history (PR #19, commit 0cd45b6). Trial D0, its artefacts and the marker
     extract moved unchanged to `Thesis/gate2_05_cardoso_2026/` as an
     extension, with the D0 identifier kept so the run record and rows C58 to
     C64 stay true; the Cardoso plan, its trials index, `CLAIMS.md`,
     `README.md`, `REFERENCES.md` and `AI_CONTEXT.md` point there now, and
     `NEGATIVE_RESULTS.md` was regenerated. The owner adds the paper-2 folder
     after reading the paper. Roadmap row 2 says "withdrawn".
   - **The reading order has Gate 2 branches.** 2C is the Choi axis (paper 5
     done; papers 12 England 2025 and 13 the 2026 IL-1beta review added), 2N
     is Nabhan (paper 6; paper 14 the 2026 PNAS platform added), and 2W is
     Wagner (papers 15 Compass and 16 the JCI 2025 lung-fibrosis circuit
     added), because the method behind proposal W1 had no paper in the
     order. Seven methods references (M1 to M7) are listed with the backbone
     step that uses each. Every identifier was verified against PubMed on
     2026-09-15. Gate labels used by folders stay 1, 2, 3A, 3B; the branch
     letter lives in the table and the JSON.
   - No analysis ran. No claim status changed.
32. **Branch labels carried across every document (2026-09-15).** After
   item 31 the branch letters lived only in `Thesis/README.md` and
   `ROADMAP.json`. They are now in the Niethamer plan's theme table and Stage
   2 proposals (W1 is 2W, Nb1 is 2N, S1 is 3A, D1 is 3B), in the Cardoso
   study note and trials index (branch 2C), in the Stage 3 header of
   `CLAIMS.md` and the Cardoso section of `REFERENCES.md`, and in
   `AI_CONTEXT.md`. One thing this exposed: the Cardoso folder uses its own
   trial gates (0, 1, 2a, 2b, 2d), which are the owner's analysis gates for
   that paper and share the word with the reading-order gates. Both
   documents and the machine context now say the two are different things,
   and no trial gate was renamed. The owner's private notes carry the same
   labels as of the same day. Documentation only.

33. **Paper 2 entered and run (2026-09-15).** After reading Choi 2020 the
   owner directed the folder: study note, extract, trials D1 to D7 with
   every rule frozen before its object was opened, and the backbone checked
   trial by trial against the plan. What reproduces from the deposit, as
   descriptions: four of five states, the DATP time course (0.3, 18.2, 6.3
   percent), the hAT2 to DATP to AT1 ordering inside one library, DATP's
   programmes in vivo, and IL-1beta's shift of the organoid epithelium with
   cell counts within 8 percent of the paper's. What does not: primed AT2
   never separates as a cluster, in vivo or in organoids (rows C88, C96).
   The annotation rule missed AT1 and the organoid stromal cluster because
   every cluster of a sorted AT2 library detects Sftpc; both outcomes stand
   and corrected passes D2b and D5b sit beside them, no threshold moved
   (rows C86, C95; decision 23). Rows C85 to C104; regenerable objects under
   `raw_data/GSE145031/choi_trials/` and `raw_data/GSE144468/choi_trials/`.
   Proposals E1 to E6 listed, not run; E6 (primed AT2 as a graded state) is
   the one the results make pressing and needs the owner's pre-registration.

34. **Two more Cardoso trials, a self-correction, and a numbering collision the
   register could not see** (`c13_epcam_transcript_in_the_transitional_state.py`,
   `c14_does_the_ranking_depend_on_the_database.py`, 2026-09-17, owner
   instruction).
   - **C13 kills the simple explanation for the sort contaminant.** The owner
     asked for a surface-EpCAM check, which three-prime counting cannot do, so
     the trial tested the one thing the data type can settle. Across five
     libraries of two deposits, Epcam transcript is equal or higher in the
     transitional state, never lower, against a frozen rule wanting a ten-point
     deficit everywhere. Transcriptional downregulation therefore does not
     explain the escape, and any surface dimming has to be post-transcriptional,
     which is consistent with the shared-sheddase account and is not evidence
     for it. ADAM17 releases amphiregulin (doi:10.1083/jcb.200307137) and also
     cleaves EpCAM with presenilin-2 (doi:10.1038/ncb1824).
   - **Both cautions on that number, pointing opposite ways.** The positive
     difference is not upregulation, because Epcam is near its detection ceiling
     in both groups and transitional cells are deeper in all five libraries. The
     absence of a deficit is robust for the same reason: a ten-point shortfall
     cannot hide in a gene detected in nine cells out of ten, and the largest
     deficit anywhere is 6.6 points, falling to 0.2 at depth.
   - **Exploratory context, caveat attached.** ADAM17 transcript does not rise,
     its only endogenous inhibitor Timp3 falls 40 points in every library, and
     both iRhom partners rise. Two of those run against the depth gradient and
     are conservative; the iRhom pair runs with it and is weakest; Timp3 loss is
     generic to injured states.
   - **C14 refutes a generalisation this repository had asserted.** After C12,
     section 12 of the divergence document said a whole-resource ligand-receptor
     scan on dissociated tissue always surfaces the matrix. Holding cells,
     donors, compartments and scoring fixed and varying only the curated
     resource, the abundance guard fires in one resource of five: CellChatDB at
     0.80, then 0.40, 0.27, 0.20 and 0.20. The domination is CellChatDB's, not
     the tissue's. Section 12 is corrected in place and points at the
     refutation rather than being rewritten.
   - **What survives better than before.** AREG ranks first among the EGFR
     ligands in every resource that contains them, four of five, with the order
     identical in three. Holding the receptor constant removes the curation
     effect along with the promiscuity effect. The heads of the rankings barely
     agree otherwise: CellChatDB shares zero of its top fifteen with italk.
   - **A numbering collision, and the check that now prevents it.** The Choi
     2020 folder took rows C85 to C104 on 2026-09-15; C13 and C14 numbered from
     C85 again and collided across eleven identifiers. The Choi rows came first
     and are cited by item 33 and decision 23, so they keep their numbers and
     these moved to C105 to C115. Nothing caught it because the register had
     only ever been read by people, so `validate_repository.py` now fails the
     build on a duplicate claim identifier.
   - **Three process failures, all mine.** Two runs were killed by the operating
     system during the largest resource with no traceback, and a truncated
     process listing was twice misread as evidence they had died, so duplicates
     were launched that competed for memory; one run was orphaned by using a
     shell background instead of the harness. C14 now checkpoints each resource
     to its own file and skips it on rerun.
   - Rows C105 to C115; sections 13 and 14 of the divergence document; trials
     index updated. Owner retain/reject review pending.

35. **A branch of paper 2 into another genomic layer, five refusals, one
   retraction and one validated deposit error**
   (`Thesis/gate1_02_choi_2020/datp_epigenetics/`, 2026-09-20, owner
   instruction).
   - **Why the branch is not on the Choi deposit.** Choi 2020 supports the
     epigenetic half of its Il1r1 claim with ATAC-seq, and GSE144598 deposits
     two bigwig coverage tracks with no peaks and no reads. Trial D0 called it
     unusable and it stays unusable, so the branch runs on two 10x multiome
     deposits, GSE310539 (Lynch 2026) and GSE247130 (Hassan and Chen 2024).
     Both are from the Jichao Chen laboratory, which matters below.
   - **The question, and why the obvious one was dropped before it ran.**
     Whether transitional cells differ from AT2 cells in accessibility cannot
     be false: the labels come from the same nucleus and both source papers
     defined the state partly from its chromatin. An adversarial check settled
     it empirically, returning 3,393 peaks at p = 0/200 on a depth-matched
     pseudo-contrast containing no transitional cell at all. The replacement
     asks whether the AT2 identity programme is closed at the chromatin level
     or merely silenced at the RNA level, which has two possible outcomes and
     speaks to a disagreement the source papers state and cannot settle.
   - **The one validated result is a deposit error.** GSE247130's barcode
     suffix map is inverted relative to its GEO sample order. Cebpa is 11 to
     14 times higher in the suffix the order calls the Cebpa knockout, in all
     three files, and Cldn4 and Sox9 agree independently. Anyone reusing that
     deposit as documented has the genotypes backwards.
   - **The substantive positive result, Descriptive only.** The
     CLDN4-positive KRT8-positive alveolar group loses the AT2 identity
     programme in RNA by 12.6 to 17.1 detection points against a sham band, in
     both deposits and at all four downsampling seeds, carried by at least five
     genes, with the AT1 arm flat throughout. The deposits share a laboratory,
     the lineage tools, the injury model and this analysis code, so their
     agreement is consistency and not replication.
   - **One of the owner's two questions is answered, by a refusal.** The
     neonatal wells carry more of the labelled group than any injured well, so
     the whole-trial negative control fired. Krt8 and Cldn4 are expressed
     across immature postnatal alveolar epithelium, so this marker set is
     shared with normal development and is not by itself damage-associated.
   - **The chromatin question is Not established, and the reading that said
     otherwise is retracted.** Trial M1e read "silenced but not closed" and
     trial M2 withdrew it: the positive control that licensed reading the null
     clears in five of eight seed-and-well combinations and fails leave-one-out
     on two of six genes. Three adversarial passes sharpened the withdrawal.
     The offset was absent rather than inert, and the global shift it failed to
     remove is positive in all three wells. The AT2 arm is bounded and the
     control is not, so an unbounded statistic was certifying the power of a
     bounded one. What was measured is that no change was detected at a bar
     corresponding to a 45 per cent loss of distal accessibility at the AT2 loci
     in one deposit and 12 per cent in the other. Where the data lean at all
     they lean toward closing, never toward retention.
   - **Five refusals, five different defects, and two of them are results.**
     A presence rule on an ambient transcript that was a depth filter under
     another name; a threshold whose own arithmetic forbade its answer; a
     marker set that is not specific to the condition; two rules that were not
     composable; and a positive control tested once. The last two are the
     eighth and ninth instances of this repository's disclosed defect shape.
   - **The lead worth acting on.** GSE309751, the bulk ATAC accompanying Lynch,
     carries genuine biological replicates, two to three mice per group across
     six arms including 49 days after infection. It is the first chromatin
     contrast anywhere in this project with animal-level replication, so the
     mouse can be the unit and reversibility, which is what "transient" in
     damage-associated transient progenitor actually asserts, becomes testable.
   - Rows C116 to C128; branch documents under
     `Thesis/gate1_02_choi_2020/datp_epigenetics/`. Owner retain/reject review
     pending.

36. **The second branch, an adversarial audit of both, and the repository
   brought current** (2026-09-20, owner instruction).
   - **The Axin2 and Il1r1 branch** (`Thesis/gate1_02_choi_2020/axin2_il1r1/`)
     assessed the question Choi 2020's own Discussion proposes and ran no
     definitive trial on it, because none is possible: both populations are
     tamoxifen-inducible lineage reporters, no deposit carries both readouts in
     the same cells, and both transcripts sit near five per cent detection at
     about one molecule per positive cell. Four routes were costed. Route B
     (locus co-accessibility) ran three passes and closed under a rule frozen
     before the first, because its positive control never cleared; the gate was
     itself mis-specified, so the closure is a decision, not a demonstration.
     Route A rested on data that was never deposited (the England figure is a
     bulk qPCR panel), and its replacement refused because the Il1r1 deletion
     is invisible to three-prime counting. Route C, a bulk array, was withdrawn
     at the owner's correction: bulk cannot answer a co-occurrence question.
     Rows C134 to C150.
   - **The audit.** Eight adversaries read every artefact behind rows C116 to
     C150, each finding sent to an independent verifier that defaulted to
     refuting it: 45 of 47 findings survived, across 14 rows. Four numbers had
     been quoted that no run produced, including the Sox9 axis of the branch's
     one validated result, which was never in the frozen panel. Three more were
     wrong against the tables they cite. One row named the wrong laboratory;
     one was refuted by the branch's own artefacts. The cause was mechanical:
     summary prose typed as string literals while the dataframe sat in scope.
     Trial M4 now logs the corroboration; thirteen rows corrected from CSVs;
     rows C151 to C154 record the audit itself. Frozen docstrings and run
     records were not edited; the summary beside them carries a banner.
   - **Repository hygiene.** The axin2 branch moved into the `trials/` layout
     its siblings use, with imports, links and eight register paths fixed and
     verified. The root README was rewritten: the Leiden UMAP is replaced by a
     claims ledger drawn from the register (`analysis/figures/claims_ledger.png`,
     `15_claims_ledger_figure.py`), the dataset table now lists all fifteen
     opened deposits with their source studies, and the claims table cites row
     ids that were checked against the register after a first draft got seven
     of them wrong. REFERENCES, REPRODUCIBILITY, AI_CONTEXT, the docs index and
     the roadmap table were brought current.
   - **`RESEARCH_QUESTIONS.md`**, new, is the question-first entry point for a
     reader arriving cold: five questions with what the data say
     and do not, nine bioinformatic phenotypes, and three follow-up analyses
     on public data. Before it landed, six checkers read it and the new README
     against the register and the tables and raised 119 findings; 17 reached
     an independent verifier before the session limit stopped the rest (16
     confirmed, 1 refuted because it had already been fixed), and the
     remainder were checked by hand against the rows and CSVs they cite. The
     confirmed defects fell into five classes: register statuses quoted one
     grade too strong (C13, C16, C40, C48, C88, C119, C120, C136, among
     others); the M2 seed range quoted from the single-budget pass (12.6 to
     17.1) instead of the seeds (12.3 to 16.7); a three-prime mechanism that
     row C149 had withdrawn but that still stood on the A2 page; row C116's
     claim text stating the opposite of its finding (the suffix order is
     inverted, not correct); and animal, donor and library counts typed
     loosely. Fifty-eight edits fixed them across CLAIMS, README,
     RESEARCH_QUESTIONS, the chromatin summary generator and the A2 page, and
     the validator passes at 740 checks.
   - Owner retain/reject review pending on rows C116 to C154 and on the
     question-first framing of `RESEARCH_QUESTIONS.md`.

37. **The two original series moved beside their source papers, and the
   README left open** (2026-09-21, owner instruction).
   - `analysis/GSE262927/` is now `Thesis/gate1_01_niethamer_2025/GSE262927/`
     and `analysis/GSE178360/` is `Thesis/ungated_murthy_2022/GSE178360/`. The
     rule is one folder per paper with each deposit beside the paper that
     produced it. Murthy 2022 is outside the roadmap and unread, so its folder
     carries a pointer note, not a study note (decision 23 still applies).
     `analysis/` keeps the shared pipeline, config, repository-level figures
     and the raw-data inventory; `pipeline_utils.SERIES_DIRS` is the one place
     the two locations are written. 373 files renamed; the 25 figure PDFs of
     the reference-aligned panels had to be force-added again at the new path
     because `*.pdf` is ignored.
   - Rewired by script, not by hand: every local Markdown link that crossed
     the boundary re-relativised, every prose path and script constant
     rewritten, `.gitignore` given the same regenerable-output rules under
     `Thesis/`, the validator's JSON sweep extended to `Thesis/`, and the
     report generator's path anchor taught the new root. Run records and logs
     written before the move keep the paths of their day; they are artefacts,
     not pointers. The archived August 2026 PDF script is not maintained and says
     so at the top of `archive/DISPLACED.md`.
   - The generated series reports and `docs/PIPELINE_AS_RUN.md` were
     regenerated. In doing so the report generator's twelve em-dashes, in the
     tracked reports since Stage 0 against the house rule, were replaced; the
     reports differ from their previous versions in paths and those characters
     only.
   - The root README gained a "Still open" section so it reads as a log in
     progress: papers 1, 2, 4 and 5 entered, Nabhan next, the Wagner branch
     unopened, Gate 3 paused, the paper-1 proposals and the Part D
     deliverables not run, rows C116 to C154 under review.
   - Owner instruction recorded in `AI_CONTEXT.md` and in the cross-project
     harness rules: weigh the cost before launching a multi-agent workflow;
     scripts for the deterministic part, agents only where a script cannot
     judge, at most three lenses. This move was done and checked by script.

38. **A paper-style figure for each research question** (2026-09-21, owner
   instruction: figures of the kind papers show, embeddings, feature and
   violin panels, dotplots, not summary bars).
   - `analysis/scripts/16_research_question_figures.py` draws one figure per
     question in `RESEARCH_QUESTIONS.md` from the analysed objects and writes
     the caption block for each between markers it owns, every caption number
     formatted from a CSV beside the figure. A1: RNA-only embedding of the
     GSE310539 wildtype nuclei with the transitional label, the AT2 identity
     programme in RNA and at the vendor-annotated promoter peaks, and the
     detection fraction of each identity gene at one depth budget (trial M1e's
     rule) in RNA and chromatin side by side. A2: a donor-medianed dotplot of
     the EGFR ligands and receptor by compartment and tissue in GSE131907,
     beside the top-15 resource overlap of trial C14. A3: the myeloid
     embedding by phase with the three states of the question coloured, the
     capillary embedding by phase under the injury-induced score (script 06
     recipe and seed, 43,359 cells reproduced), and the per-animal iCAP time
     course. A4: Axin2 and Il1r1 on the A1 embedding, violins with the
     positive nuclei drawn as points in AT2 nuclei of all four wells, and
     co-detection tiles. A5: one embedding per GSE247130 control well with the
     transitional label, Cldn4 and Krt8.
   - What they are not: evidence. Nothing is tested, no threshold is set, the
     multiome embeddings carry no batch correction, and the first draft of A1
     showed a per-nucleus promoter count by group that depth dominates (the
     lesson of rows C127 and C130); it was replaced by the detection-at-budget
     form before landing. Embeddings are cached under an ignored `processed/`
     folder; `run_record.json` holds parameters, seeds, versions and counts.
   - The owner asked for each question's roadmap branch and a further
     mapping, then decided that mapping does not belong in the repository:
     it lives in his private notes, written from this repository's reading
     of fit. Part D's items are numbered 1 to 3 so they no longer collide
     with proposal D1.

39. **Where a gene set enrichment analysis is admissible** (2026-09-21, owner
   question: is any GSEA established, and on which imported data would a
   trial be plausible). None exists; the pipeline record says "not run" and
   every "enrichment" in the scripts is a fold or kNN enrichment. Trial G0
   (`Thesis/gate1_01_niethamer_2025/trials/g0_gsea_feasibility/`) scanned the
   16 imported accessions from their metadata tables and this
   repository's reality-check records, no matrix opened, under rules frozen
   first (unit of animal, donor or library with at least three per arm; raw
   counts on disk; a contrast a question asks; gene sets hashed before
   ranking; a matched-random-set null beside the permutation null). Verdicts:
   4 plausible (GSE262927, GSE136831, GSE135893, GSE131907), GSE309751 marginal,
   11 not admissible. Row C155. Proposals G1 (GSE262927
   by phase within compartment, the gene-set layer of W1) and G2 (GSE136831
   IPF against control per donor, GSE135893 held out) are in the Stage 2
   table, not run; MSigDB is not on disk and no GSEA package is installed,
   so either trial begins by recording those as inputs.

40. **The first gene set enrichment analyses, G1, G1b and G2** (2026-09-21,
   owner instruction to proceed where G0 found it plausible). Inputs recorded
   first: MSigDB 2024.1 hallmark and GO biological process for mouse and human
   (sha256 in each run record, files under `raw_data/msigdb/`) and gseapy 1.3.1
   in `.venv-x64`. Engine: gseapy pre-ranked GSEA on per-unit pseudobulk
   rankings, every hallmark enrichment score reproduced by an in-house running
   sum to 1e-6, and an expression-matched random-set null beside the
   permutation null; a set clears only when both agree.
   - G1 (GSE262927, 8 active-repair, 12 resolution and 3 long-term
     animals per compartment): the G2M positive control failed in both
     compartments, and trial N1's cycling fractions show why (no early
     proliferation peak in either lineage). Recorded, not moved (C156). Its
     first run also died writing its record (tuple keys in a dictionary); the
     script was corrected and rerun with the same seed and rules.
   - G1b sits beside G1 with the gate replaced by repository injury-state
     marker sets (iMON at 6 and 11 dpi, iCAP at 19 and 25 dpi, each against
     366 dpi, four against three animals); both clear (C157). Its rule was
     written after G1's tables were seen, and its record says so.
   - Read through that gate: myeloid DNA-replication programmes still higher at
     42 and 90 dpi than at a year, Exploratory (C158); nothing separates the
     capillary compartment at three animals, Not established (C159).
   - G2 (GSE136831 discovery, GSE135893 held out, donor as the unit): G2 ran the same day: 254 discovery sets cleared both nulls and 62 replicated (rows C160, C161).
   - Under x86-64 emulation a 1000-permutation GO run takes about 35 minutes per
     contrast; the held-out pass is asked only about the sets that cleared
     discovery, which is what R8 reads anyway.
   - **Owner review done 2026-09-22 (DEVELOPMENT decision 28).** Rows C155
     to C161 are all retained at their proposed statuses. C161 stays
     Validated, narrowed to the donor-level direction, with a third limit
     (gene-sampling nulls; replication does not exclude a shared artefact;
     the eight AT2 sets weakest). C156 gained the timing sentence (trace
     window against cycling at sacrifice), C158 a caveat sentence and the
     cross-species lead, C160 a pointer to C161. The paper-1 plan's count of
     single-library deposits was corrected from ten to eight against G0's
     table.

---

## Handoff: session of 2026-09-15, paper 2 entered and run

This supersedes item 1 of the close-of-session handoff below ("nothing on
paper 2"): the owner read Choi 2020 later the same day and directed the
entry. One pull request on top of PR #29: the paper-2 folder with study
note, extract, trials D1 to D7, corrected passes D2b and D5b, rows C85 to
C104, decision 23, and the pointers that PR #27 had moved to the Cardoso
folder moved back. No other analysis changed.

**Owner decisions open:** retain or reject on rows C85 to C104 (item 33), on
top of the rows listed below; whether to pre-register proposal E6; whether to
download the van den Brink 2017 list (M8) for attack A3.

**What the next session should do, in order.**

1. Paper 3, Nabhan 2018, as the previous handoff says: the study note is the
   owner's to write or direct after reading; an agent may run a Gate 0
   deposit check when asked; the accession is not yet in `ROADMAP.json`.
2. Nothing further on paper 2 unless the owner pre-registers E6 or asks for
   A3.
3. Items 3 to 5 of the previous handoff stand.

The run order for paper 2, if anything is re-run: D1, D2, D2b, D3, D4, D5,
D5b, D6, D7 in the x86-64 environment (D1 about 15 minutes, the rest under
5 each), one at a time.

---

## Handoff: close of session 2026-09-15

Main is at PR #28. Three pull requests landed today, all documentation, all
merged with CI green and their branches deleted. No analysis ran, no artefact
content changed, no claim status changed. The working tree is clean.

| Pull request | What it did | Item |
|---|---|---|
| #26 | roadmap re-ranked: Nabhan papers next, Gate 3 paused; `AI_CONTEXT.md` and the status table above brought in line | 30 |
| #27 | decision 21 corrected from the session transcript; the Choi 2020 folder withdrawn and its deposit check relocated as Cardoso trial D0; the reading order given Gate 2 branches 2C, 2N, 2W, papers 12 to 16 and methods references M1 to M7 | 31 |
| #28 | branch labels carried into every document that names a gate or a branch; the reading-order gates and the Cardoso trial gates stated to be different namespaces | 32 |

Outside git, on the same day, the owner's private notes carry the same
labels.

**Owner decisions open**, unchanged in substance: retain or reject on rows C6,
C9 to C18, C19 to C57, C58 to C64 and C65 to C84 (items 12 to 29); one
planning decision tracked in the owner's private notes (item 30).

**What the next session should do, in order.**

1. Nothing on paper 2. The owner reads Choi 2020 and re-enters its folder;
   an agent does not write that note (decision 21).
2. Paper 3, Nabhan 2018. The study note is the owner's to write or direct. An
   agent may run a Gate 0 deposit check when asked; the accession is not yet
   recorded in `ROADMAP.json`, so finding it is the first step.
3. Proposal Nb1 needs a W1-style per-animal pseudobulk for epithelium and
   mesenchyme before any ligand-receptor step (backbone step 5 after step 3).
4. The Wagner packet needs no new analysis: the one-page summary and 3 to 5
   figures from tracked artefacts.
5. Nothing in Gate 3 until S1 or D1 has a result.

Read the `AI_CONTEXT.md` pitfalls before starting, in particular the two added
today: the owner reads a roadmap paper before any note is written, and the two
gate namespaces are never conflated.

---

## Handoff: session of 2026-09-15, reading order and the Choi relocation

Branch `Claude/reading-order-and-choi-relocation`. Documentation, git moves
and one regenerated page; no script ran on data.

| What | Where | State |
|---|---|---|
| D0 trial, artefacts, `choi_utils.py`, `choi_2020_extracts.json` | moved with `git mv` into `Thesis/gate2_05_cardoso_2026/` (trials and folder root) | complete; history preserved |
| Choi study note and D-series plan | removed from main; recoverable from PR #19 | complete |
| D0 section, TOC row, trial-table row | `Thesis/gate2_05_cardoso_2026/ANALYSIS_TRIAL_PLAN.md`; index row in `trials/README.md` | complete |
| Decision 21 and its responsibility row | `DEVELOPMENT.md` | complete |
| Gate 2 branches, papers 12 to 16, methods M1 to M7 | `Thesis/README.md`, `Thesis/ROADMAP.json`, `REFERENCES.md` | complete; validator passes |
| The owner's private reading-order notes | outside git | revised the same day |

**Next.** Paper 3 (Nabhan 2018) remains next, and paper 2 is the owner's to
re-enter. Nothing in Gate 3 until S1 or D1 has a result.

---

## Handoff: session of 2026-09-15, roadmap re-ranking

Branch `Claude/roadmap-reranking-2026-09-15`. Documentation only; no script
ran and no artefact changed.

| What | Where | State |
|---|---|---|
| Re-ranking recorded in the roadmap | `Thesis/README.md` gate rules and status column; `Thesis/ROADMAP.json` | complete; validator passes |
| Machine context brought in line | `AI_CONTEXT.md` status and `thesis_roadmap` block | complete |
| Status table and item 30 | this file | complete |
| The owner's private notes updated to match | outside git | complete |

**What the next session should do, in order.** (1) Nabhan 2018 (roadmap
paper 3): study note in the five-question format, extracts JSON, and a Gate 0
data reality check; the deposit accession is not yet recorded (`datasets` is
empty in `ROADMAP.json`), so finding it is the first step. (2) Proposal Nb1
depends on a W1-style per-animal pseudobulk for epithelium and mesenchyme,
which has not been built; build it before any ligand-receptor step, per
backbone step 5. (3) Nothing in Gate 3 until S1 or D1 has a result.

---

## Handoff: session of 2026-09-13, paper 2 and the repository chores

Branch `Claude/choi-2020-gate0-and-chores`. Four things landed after the
Cardoso follow-ups, in this order.

| What | Where | State |
|---|---|---|
| The two items left open by PR #12 | `analysis/config/palette.json`, `analysis/scripts/08_*` | closed; one validated palette, and the "AT0 candidate" label corrected |
| A generated negative-results page | `NEGATIVE_RESULTS.md`, `analysis/scripts/14_*` | complete; 23 rows, regenerate after any status change |
| Generic deposit readers promoted to the shared module | `Thesis/gate1_04_sikkema_2023_hlca/trials/trial_utils.py` | complete; `cardoso_utils` re-exports, all 19 Cardoso trial modules still import and two were re-run to prove it |
| Roadmap paper 2, study note and Gate 0 | `Thesis/gate1_02_choi_2020/` | note and D0 complete; D1 not started |

**Why the helper module moved.** The Choi folder needed the MatrixMarket
triplet reader, the SOFT parser and the non-Ensembl feature rule, and a gate1
folder importing them from a gate2 folder is the wrong dependency direction.
They now sit in the shared module beside RunRecord, and `cardoso_utils` is a
re-export so no Cardoso trial changed. Verified by importing all 19 trial
modules and by re-running E2b and E5, which reproduced identical results.

**What Gate 1 on paper 2 has to decide before it runs.** Cell calling needs a
frozen threshold, because six libraries are raw droplet matrices. The primed
AT2 state is defined by the loss of Etv5, Abca3 and Cebpa rather than by a
positive marker, so the scoring rule must say what it does about that in
advance. And the trajectory question is the one part of this paper a deposit
without replication can genuinely address, because an ordering inside a library
needs no between-group comparison.

**Owner decisions now open**, in addition to those in items 23 to 25: rows C58
to C64, and whether the post hoc Fst with Runx2 lead from trial C8 deserves a
pre-registered trial of its own.

---

## Handoff: session of 2026-09-13, the three chosen follow-ups

Branch `Claude/cardoso-2026-three-followups`, opened after
[PR #11](https://github.com/xorca0711/scRNA_seq/pull/11) was merged. Nothing in
`analysis/` was touched.

| Trial | Question | State | Re-run |
|---|---|---|---|
| E6 | donor-level coupling of AREG to EGFR in human fibrosis | complete, null | about 30 minutes first time, then seconds from the cached vectors |
| C7 | which epithelial state the mesenchymal-sort contaminant is | complete, mixture | about 25 minutes; reprocesses the GSE316244 RFP libraries |
| C8 | subpopulation or gradient in the Areg-independent tier | complete, gradient | about 12 minutes; reprocesses the GSE316244 niche libraries |

**Where the work stands.** List A of
[`DIVERGENCES_AND_NEXT.md`](Thesis/gate2_05_cardoso_2026/DIVERGENCES_AND_NEXT.md)
is empty. Every question answerable from data already on disk has been asked,
and five of the six answered negatively or forced a correction. What remains on
that page is list B, which needs the bench, and list C, which the data type
cannot support.

**The three results in one line each.** The paper's axis shows no donor-level
coupling that 22 donors could detect, and the smallest effect the test could
have seen is about rho 0.43. The mesenchymal-sort contaminant is a mixture,
roughly four parts DATP-like to four parts AT2, so it cannot be named as a
single state. The Areg-independent tier is a gradient rather than a
subpopulation, confirming the C5 correction from a second direction.

**Two rule defects are on the record and neither threshold was moved.** C7's
profile-correlation test was written without a margin requirement, so it named
a state it separated by 0.0056 while separating compartments by 0.415; the
margin it reported is what refuses the answer. C8's mixture test prefers two
components because 31.5 per cent of the cells detect neither gene in the score,
which guarantees a spike at zero; trial C1 applied the same test to a centred
score where that does not arise, so the defect is the reuse. One implementation
bug was fixed: E6's depth rule covered any pair in its text but only the
primary tests in code, and the code was corrected to match the text.

**The thing to carry into any future correlational trial here.** E6's only
significant correlation was a control pair whose two members both track
sequencing depth, at p = 0.045. Without that pre-registered control it would
have read as a discovery. Carry the depth control by default.

**Owner decisions now open**, in addition to those in items 23 and 24: rows C49
to C57, and whether the post hoc Fst with Runx2 lead from C8 deserves a
pre-registered trial. That lead rests on 24 double-positive cells in one
library, though the two arms are the best depth-matched pair in this deposit at
3,637 against 3,612 median genes per cell.

**Caches added under `raw_data/` this session** (gitignored, regenerable):
`GSE136831/e6_extracted_rows.npz`. The GSE131907 cache from E1 is still there
and still used by E1b.

---

## Handoff: session of 2026-09-13 (the E series)

Branch `Claude/cardoso-2026-public-extensions`, [PR #11](https://github.com/xorca0711/scRNA_seq/pull/11),
merged 2026-09-13 on the owner's instruction, CI green. [PR #10](https://github.com/xorca0711/scRNA_seq/pull/10) was merged on
2026-09-13 before the last four commits of the previous session landed, so this pull request
carries the C5 figures, trial C6 and the whole E series. Nothing in `analysis/` was touched. Six further trials ran,
all outside the Cardoso deposit, because that deposit cannot test anything and
the standing instruction was to extend the Hbegf lead into public data only if
the first extension did not refute the working picture. It did not, so the rest
ran.

| Trial | Dataset | State | Re-run |
|---|---|---|---|
| E1 | GSE131907 human LUAD | complete; refutation rule did not fire | about 25 minutes first time, then seconds from the cached gene vectors |
| E1b | GSE131907, deposited subtypes | complete | seconds; reads E1's cached `e1_extracted_rows.npz` |
| E2 | GSE136831 human IPF | complete | about 35 minutes; streams a 2.0 GB MatrixMarket file |
| E3 | GSE135893 human PF | complete; same script as E2 | about 20 minutes; streams a 1.0 GB file |
| E2b | E2's own tracked table | complete | seconds; supplies the T4 reading E2 declared and did not compute |
| E4 | GSE132771 mouse bleomycin | complete | about 12 minutes |

**The result that changes a previous conclusion.** E4 refuted the reading of
claim C29. Runx1 and Pdgfrb both rise with bleomycin alone in sorted
Col1a1-GFP mesenchyme with no oncogene present, clearing a rule frozen before
the matrices were opened, so the pre-specified reading applies: their survival
of Areg deletion is what an activated lung fibroblast does after injury, not a
second tumour-specific signal. C29 keeps its numbers and loses its
interpretation. Hbegf and Egfr are injury-generic too, so neither is a
tumour-specific feature of the mesenchyme.

**The one row that is now Validated.** E1 is the first trial in this section
whose unit permits a test, because GSE131907 has eleven donors with paired
tumour and normal lung. AREG detection is higher in epithelium than in myeloid
cells within donor, 0.336 against 0.215, paired Wilcoxon p = 0.0020.

**Two self-corrections to read alongside it.** E1's top HBEGF compartment is
recorded as "neutrophil"; that gate is 83 per cent deposited myeloid cells in
a deposit that annotates no neutrophils, so the outcome is the pre-named
myeloid-dominant one. And at subtype resolution CD1c-positive dendritic cells
detect AREG in 0.618 of cells, above both tumour epithelial states, so the
test supports epithelium over the myeloid average and not epithelium as the
source.

**What did not work, on the record.** E2's only comparison to clear its donor
floor failed to detect the difference it was built on, HBEGF myeloid 0.366
against aberrant basaloid 0.390 over seven donors, p = 0.297. The Zhao
prediction (AREG higher in the transitional state than AT2) is not testable in
either fibrosis cohort, since no comparison clears the five-donor floor and the
two cohorts disagree in direction. E1b's T5 is not computable at all, because
the GSE131907 annotation assigns AT2 only in normal lung and the tumour states
only in tumour lung, so zero donors pair; the floor was not lowered, because
the count is zero rather than small.

Owner decisions now open, in addition to those in item 23: the weakening of
C29, rows C37 to C48, and whether the dendritic-cell and monocyte ligand
source deserves a paragraph in the writeup. It replicates across three
datasets and two diseases, and it is established immunology rather than a
finding of this repository (Zaiss et al. 2015,
doi:10.1016/j.immuni.2015.01.020), so its value is that it constrains an
epithelium-to-fibroblast reading rather than that it is new.

Downloads added to `raw_data/` this session (gitignored, about 6.5 GB):
GSE131907, GSE136831, GSE135893, GSE132771. As before, the tracked
`analysis/raw_data_inventory.*` describes the Stage 0 downloads only and the
validator checks that count; do not re-run `01_scan_raw_data.py` without
updating the validator.

**Two things learned about running these.** The human atlases deposit one
merged MatrixMarket file each, too large to load whole on this machine, so
`trials/mtx_stream.py` extracts selected gene rows in a single streaming pass
and memory stays at a few hundred megabytes. And a compartment marker gate
named after a cell type is not that cell type: always crosstab the gate
against the deposited labels before naming an outcome after it, which is what
caught the neutrophil error.

---

## Handoff: session of 2026-09-12

Branch `Claude/cardoso-2026-gate0-gate1`, [PR #10](https://github.com/xorca0711/scRNA_seq/pull/10),
merged 2026-09-13, but only up to trial C3; the later commits moved to PR #11. Nothing in `analysis/` was touched; all new material is under
`Thesis/gate2_05_cardoso_2026/`. **All eight trials have run**; the trial
plan, the claims register (rows C19 to C36), the root README claims table and
DEVELOPMENT decision 18 are written.

| Trial | Gate | State | Re-run |
|---|---|---|---|
| C0 | 0 | complete | idempotent, about 4 minutes |
| C1 | 1 | complete, verdict not recovered | about 20 minutes |
| C1b | 1 | complete; wrote `c1b_mesenchyme.h5ad` (gitignored) that C1c and C1d reuse | reuses the object if present |
| C1c | 1 | complete | about 1 minute, needs C1b's object |
| C1d | 1 | complete | about 1 minute, needs C1b's object |
| C2 | 2b | complete | about 45 minutes |
| C2b | 2b | complete | seconds; reads only C2's tracked tables |
| C3 | 2a | complete | about 50 minutes |

What is left is yours, not a computation: the retain/reject decisions listed
in item 23, chiefly whether cluster 14 is to be called the published
population, and whether the Hbegf lead (claim C33) is worth pursuing.

Two trials that were considered and not run, each for a stated reason. Gate 2
branch (d), the transcriptomic ordering of fibroblast against macrophage
change, is closed because the deposit has no time course in those
compartments. A CellChat rerun is impossible here because CellChat is R-only
and this machine has no R; trial C3 re-derives the expression fact the
communication claim rests on instead and states in its own output what it is
not.

Downloads added to `raw_data/` this session (gitignored, about 1.1 GB):
GSE316241, GSE316243, GSE316244, GSE310335, GSE247505, each with its GEO SOFT
family file. The tracked `analysis/raw_data_inventory.*` still describes the
Stage 0 downloads only, and the validator checks that count; do not re-run
`01_scan_raw_data.py` without updating the validator.

**Memory constraint learned this session:** the machine has 15.6 GB and each
trial peaks near 3 GB under x86-64 emulation, so these trials must be run one
at a time.

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

`analysis/` is split by series: `Thesis/gate1_01_niethamer_2025/GSE262927/` (mouse) and
`Thesis/ungated_murthy_2022/GSE178360/` (human), with `scripts/`, `requirements.txt` and
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
- `docs/PIPELINE_AS_RUN.md` is **generated**, edit
  `analysis/scripts/05_write_pipeline_as_run.py` and re-run, never the `.md`.
- The tool reference pages under `docs/` are deliberately preserved; they
  describe the published method, not what ran here.
