# Progress and handoff state

Living record of what is done, what is pending, and what a future session needs
to know to continue. Update this before stopping.

Last updated: 2026-09-13 (the three chosen follow-ups, E6, C7 and C8; see the newest handoff block). Previously 2026-09-12 and 2026-09-10. The scientific analysis of the two series is
complete (state of 2026-08-09 below): `FINDINGS.md` (results with figures)
leads, `README.md` is a landing page for the executed analysis, and
`scRNAseq_workflow_Niethamer2025.md` lives in `docs/`. On 2026-09-09 a
paper-by-paper roadmap directory was added (`Thesis/`, order taken from the
owner's Notion PI Target Map) with the Sikkema 2023 HLCA study note, its
decision criteria as reviewable JSON, a pipeline-framing proposal, and a
first criteria trial on tracked tables. **Owner retain/reject review of that
material is pending.** On 2026-09-10 two focused analyses of the source
paper's phase structure were added under `analysis/GSE262927/`
(`phase_timecourse/` and `myeloid_focus/`, scripts 10 to 13; items 18 to
21 below). **Owner review of those is pending too.** Later on 2026-09-10 the
S2 trial was merged (PR #8) with the items renumbered chronologically, and
the repository was reframed as an analysis log: portfolio material and the
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
| Portfolio restructure (`FINDINGS.md`, README rewrite, root tidy-up) | **DONE** |
| `Thesis/` roadmap index (11 papers, Notion order, PubMed-verified IDs) | **DONE** |
| Sikkema 2023 (HLCA) study note, `integration_benchmark.json`, `PIPELINE_FRAMING.md` | **DONE, owner review pending** |
| Trial S1: HLCA cluster-entropy criteria on the tracked tables | **DONE (2026-09-09); artefacts in Thesis/gate1_04_sikkema_2023_hlca/trials/** |
| Trial S2: scArches mapping of GSE178360 to the HLCA core | **DONE (2026-09-09)**; 23 of 31 clusters agree with the blind proposals; AT0 is a minority and the AT0 candidate subcluster is mostly AT2 or uncertain; our mapping matches the HLCA authors' own transfer of the same cells at 99.2% (level 3) |
| Trial S3: HLCA consensus-marker annotation of GSE178360 | **DONE (2026-09-09)**; 17 to 19 of 31 clusters agree with the blind proposals; AT0 by marker transfer Not established (scheme-dependent) |
| Trial S4: mouse cluster 23 explained | **DONE (2026-09-09)**; low-count, ambient-like; 78% from EEM-scRNA-289 |
| Trial S5: mouse cluster 5 subclustered and re-graded | **DONE (2026-09-09)**; resolved at Leiden 0.5 (94% of labelled cells in pure subclusters), not at 0.2 |
| Phase-wise view of the Ki67 atlas (`analysis/GSE262927/phase_timecourse/`, script 10) | **DONE (2026-09-10), owner review pending**; per-dpi atlas UMAP, per-animal lineage composition, Ki67-trace proliferation by lineage; trace peaks fall in the paper's window for 4 of 5 lineages (Lymphoid peaks at 11 dpi, not 6); Descriptive only, from tracked metadata |
| Myeloid compartment by dpi (`analysis/GSE262927/myeloid_focus/`, script 11) | **DONE (2026-09-10), owner review pending**; 9,997 cells from atlas clusters 5, 17, 24; 16 blind subclusters at Leiden 0.5 with 87% of labelled cells in pure subclusters; aMAC loss and iMON expansion at 6 dpi with reconstitution by 19 to 42 dpi, consistent with the paper's Figure 3; Descriptive only |
| Alveolar macrophage origin by Ki67 trace window (`myeloid_focus/amac_origin/`, script 12) | **DONE (2026-09-10), owner review pending**; the 2 to 3 dpi window labels most of the 42 dpi aMAC pool (median 79.7%); marrow-inheritance and two-source checks Not established under the 30-cell floor |
| Batch sensitivity of the myeloid embedding, Harmony on infection round (`myeloid_focus/batch_sensitivity/`, script 13) | **DONE (2026-09-10), owner review pending**; rounds already mix within every tested day (enrichment 1.11 to 1.46, threshold 2); the 6 dpi iMON state survives correction as its own subcluster (92% of cells from 6 dpi, 77% and 74% of each animal's iMON cells); the frozen survival rule selected the wrong subcluster and the revision is disclosed |
| Trial S2 merged into main; PROGRESS and DEVELOPMENT renumbered chronologically | **DONE (2026-09-10, PR #8)** |
| Repository hygiene and reframing: `archive/` for displaced material, README opens with the claims table, injury model named once, checks renamed | **DONE (2026-09-10), owner review pending**; see item 22 |

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
   instruction). The portfolio PDF and its generator moved to
   `archive/portfolio_2026-08/`; the Krt8-high transitional narrative
   (FINDINGS section 1, the human KRT8 reference-aligned panels) moved to
   `archive/DISPLACED.md` because it is established elsewhere (the owner's
   G-SURF submission). Their artefacts and scripts stay under `analysis/`
   because the validator checks their numbers and the scripts regenerate
   them; only the narrative moved. The CI workflow and validator were
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
- `docs/PIPELINE_AS_RUN.md` is **generated**, edit
  `analysis/scripts/05_write_pipeline_as_run.py` and re-run, never the `.md`.
- The tool reference pages under `docs/` are deliberately preserved; they
  describe the published method, not what ran here.
