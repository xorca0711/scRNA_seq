# Niethamer 2025: the analysis in steps, from the initial run to the phase and myeloid follow-ups

This page separates two kinds of work on GSE262927. **Stage 0** is the initial
run (August 2026): an independent pipeline built from the raw count matrices,
asked whether the published biology can be recovered. **Stage 1** is the set of
follow-ups (September 2026) commissioned after the owner fixed the framing
(DEVELOPMENT.md decision 12: cell state,
niche, macrophage and monocyte states, annotation robustness, curation
hygiene; no interferon or influenza narrative). **Stage 2** lists what is
proposed and not run. Every Stage 1 step has a run record with its rules
frozen before the data were opened; every status label follows the root
README vocabulary (Validated, Descriptive only, Exploratory,
Retracted-superseded, Not established). Owner retain/reject review of all
Stage 1 material is pending.

The themes referred to below follow the reading-order branches of
[`Thesis/README.md`](../README.md):
2N (AT2 differentiation, stem-niche ligand and receptor candidates),
2W (immune-metabolic state inference, robust cross-cohort integration),
3A (cytokine and receptor programmes of inflammation resolution and
epithelial repair), 3B (Treg stability and plasticity).

## Stage 0. Initial run: can the published biology be recovered from raw counts

| Step | Question | Evidence and rule | Artefact | Status |
|---|---|---|---|---|
| 0.1 | What is in the deposit, and what QC does each library need | Raw-data scan; per-sample MAD thresholds with recorded rationales; Scrublet per capture | `analysis/raw_data_inventory.*`, `Thesis/gate1_01_niethamer_2025/GSE262927/qc/` | Validated (counts checked by `validate_repository.py`) |
| 0.2 | Does blind clustering recover the deposited cell types | 162,175 cells, 29 Leiden clusters, no batch correction (root decision 2); deposited labels held out and used only as an answer key | `Thesis/gate1_01_niethamer_2025/GSE262927/figures/umap/`, `tables/cluster_vs_author_celltype_fraction.csv` | Validated (median purity 0.947) |
| 0.3 | Where does a marker-panel annotation fail | 3 of 29 candidate annotations contradicted by the labels and kept on display | `tables/cluster_annotation_proposals.csv` | Validated, negative result kept |
| 0.4 | Which epithelial states are present | Epithelial sub-analysis, 16 subclusters | `Thesis/gate1_01_niethamer_2025/GSE262927/epithelial_subanalysis/` | Descriptive only |
| 0.5a | Is the AT2 to Krt8-high transitional to AT1 axis a real ordering | PAGA plus diffusion pseudotime rooted in AT2 on the 25-sample cohort; labels held out | `Thesis/gate1_01_niethamer_2025/GSE262927/regeneration_focus/` | displaced to [`archive/DISPLACED.md`](../../archive/DISPLACED.md) (established elsewhere); artefacts unchanged |
| 0.5b | Does the injury capillary state resolve | Capillary subset reclustered within compartment; per-animal abundance by day | `Thesis/gate1_01_niethamer_2025/GSE262927/regeneration_focus/` | Validated (iCAP 2.0% at baseline to 37.5% at 25 dpi to 21.7% at 366 dpi; has not resolved by 366 dpi) |
| 0.6 | Where does the injury capillary state come from | 8-sample Cre cohort; trace rule reproduces the deposited calls at 100% | `Thesis/gate1_01_niethamer_2025/GSE262927/lineage_tracing_cohort/` | Validated for the Kit line (33 to 53% per animal); CAP2 lines uninformative |
| 0.7 | Can the whole thing be checked without data | Generated reports, pipeline record, validator, CI | `docs/PIPELINE_AS_RUN.md`, `REPRODUCIBILITY.md` | done |

Stage 0 answered the paper's epithelial and endothelial claims (the epithelial
trajectory answer is displaced to `archive/DISPLACED.md`, step 0.5a). It did not
touch the paper's phase structure or its myeloid compartment, and it drew no
per-day pictures; those are the gaps Stage 1 fills.

## Stage 1. Follow-ups: the paper's phase and myeloid claims

All four steps use the 25-sample Ki67 atlas only, treat the animal as the unit,
report medians, and compute no P value (two animals per active-repair day).
Deposited labels are held out of every embedding and used afterwards for
grading and for descriptive composition.

### N1. Phase-wise view of the atlas (script 10, tracked metadata only)

- **Question.** Is the paper's three-phase proliferation (immune first,
  epithelium and mesenchyme second, endothelium last; Figure 2F) visible in the
  reanalysed atlas?
- **Serves.** The backbone figure for every reader. Theme 2W: the immune wave
  comes first. Theme 2N: epithelium and mesenchyme co-peak at 11 dpi, which
  fixes the time point for any stem-niche question.
- **Rule, frozen.** For each lineage, the immediate-window harvest (6, 11, 19,
  25 dpi) with the highest median per-animal Ki67-traced fraction is the peak;
  expected windows taken from the paper before the table was opened.
- **Outcome.** Trace peaks agree for 4 of 5 lineages (myeloid 6, epithelium and
  mesenchyme 11, endothelium 19 dpi; lymphoid peaks at 11 rather than 6). The
  deposited cell-cycle call agrees for 1 of 5 and is reported as a cross-check
  only. Per-day atlas UMAP with equal cell numbers per panel; per-animal
  composition within compartment.
- **Artefacts.** [`Thesis/gate1_01_niethamer_2025/GSE262927/phase_timecourse/`](GSE262927/phase_timecourse/README.md)
- **Status.** Descriptive only; owner review pending (PROGRESS item 18).

### N2. Myeloid compartment by day (script 11)

- **Question.** Does the paper's Figure 3 (aMAC loss and iMON expansion at 6
  dpi, reconstitution over three weeks) reproduce from a blind embedding?
- **Serves.** Theme 2W primary (macrophage and monocyte states as the
  immune-metabolic object); theme 3A secondary (the iMON to aMAC rebuild is the
  resolution of the myeloid response).
- **Rule, frozen.** Cells = atlas clusters 5, 17, 24 in the 25-sample cohort;
  Leiden 0.5 fixed from trial S5; resolved if at least 75% of labelled cells
  sit in subclusters of purity at least 0.75; Figure 3 consistent if aMAC
  median falls at 6 dpi and recovers by 25 or 42 dpi while iMON rises at 6 dpi.
- **Outcome.** 9,997 cells; 16 subclusters; 87% of labelled cells in pure
  subclusters (48% at Leiden 0.2, 92% at 1.0). aMAC 30.8% to 4.6% at 6 dpi to
  49.0% at 42 dpi; iMON 1.9% to 56.0% to 1.8%. Consistent with Figure 3.
  Interstitial macrophages keep rising through 90 dpi (candidate only).
- **Artefacts.** [`Thesis/gate1_01_niethamer_2025/GSE262927/myeloid_focus/`](GSE262927/myeloid_focus/README.md)
- **Status.** Descriptive only; owner review pending (PROGRESS item 19).

### N3. Origin of the rebuilt alveolar macrophage pool, trace only (script 12)

- **Question.** Which tamoxifen window's proliferating cells populate the 42
  dpi aMAC pool, and can the labelling be explained by marrow inheritance?
- **Serves.** Theme 2W primary (origin of a reconstituted macrophage pool
  inferred from a design variable, not a fitted trajectory); theme 3A secondary.
  Prior art: Aegerter et al. 2020 [DOI](https://doi.org/10.1038/s41590-019-0568-x)
  (GSE131252 and companions) and Iliakis et al. 2023
  [DOI](https://doi.org/10.1038/s41590-023-01602-1).
- **Rule, frozen.** Reporter-scored aMAC-labelled cells per animal; 30-cell
  floor; main window = highest median traced fraction at the common 42 dpi
  harvest; marrow reference = the same animal's cMON and neutrophils, exceeded
  by at least 10 points in both animals of a window; within-aMAC split = a
  subcluster pair differing by at least 20 points in both animals.
- **Outcome.** The 2 to 3 dpi window labels most of the 42 dpi pool (median
  79.7%, then 60.0, 29.3, 42.0), as the paper implies. The marrow reference
  was evaluable in 3 of 8 animals and the within-aMAC split in one animal per
  window, so both checks are Not established under the frozen floor; the
  per-animal values (early windows label subcluster 2 more than subcluster 3
  by 28 to 33 points, late windows the reverse) are reported descriptively.
- **Artefacts.** [`Thesis/gate1_01_niethamer_2025/GSE262927/myeloid_focus/amac_origin/`](GSE262927/myeloid_focus/amac_origin/README.md)
- **Status.** Window contribution Descriptive only; the two checks Not
  established; owner review pending (PROGRESS item 20).

### N4. Batch sensitivity of the myeloid embedding (script 13)

- **Question.** Is the 6 dpi iMON state biology or a one-day batch island?
- **Serves.** Robustness for every reader, and the robust-integration theme
  (2W). Belongs in the figure's limitation line for the other themes.
- **Rule, frozen.** Correct on infection round from Table S3 (the two rounds
  also differ in Ki67-Cre dosage), never on sample; test only the days that
  carry both rounds (6, 11, 19, 25, 42 dpi); compare uncorrected and Harmony
  embeddings on the pipeline's kNN mixing metric; the iMON state survives if
  its subcluster keeps iMON purity at least 0.75 and holds at least 20% of each
  6 dpi animal's iMON cells.
- **Outcome.** The rounds already mix within every tested day (same-round kNN
  enrichment 1.11 to 1.46 uncorrected, 1.05 to 1.29 after Harmony, threshold
  2.0). After Harmony the 6 dpi iMON state remains its own subcluster (411
  cells, 92% from 6 dpi, purity 0.86, holding 77% and 74% of the two animals'
  iMON cells); adjusted Rand index between the partitions 0.837. The frozen
  survival rule selected the 11 to 19 dpi monocyte state instead of the 6 dpi
  one; the first-run outcome is kept and a post hoc definition anchored on
  the 6 dpi cells is reported alongside (DEVELOPMENT decision 16).
- **Artefacts.** [`Thesis/gate1_01_niethamer_2025/GSE262927/myeloid_focus/batch_sensitivity/`](GSE262927/myeloid_focus/batch_sensitivity/README.md)
- **Status.** Descriptive only (rule revision disclosed); owner review pending
  (PROGRESS item 21).

### The figure each theme rests on

| Theme | Figure | Sentence it supports |
|---|---|---|
| Stem niche (branch 2N) | `phase_timecourse/figures/proliferation_by_lineage.png` | Epithelium and mesenchyme proliferate in the same window (11 dpi), which is when a stem-niche question should be asked. This entry rests on the proliferation trace only; the Krt8-high trajectory (step 0.5a) is displaced and is not used here. |
| Immune-metabolic states (branch 2W) | `myeloid_focus/figures/myeloid_label_composition_by_dpi.png` and `amac_origin/figures/trace_by_window_late_harvests.png` | Monocyte and macrophage states turn over on a per-animal, trace-dated timeline. |
| Resolution programmes (branch 3A, paused) | `myeloid_focus/figures/myeloid_label_composition_by_dpi.png` | The myeloid response resolves between 11 and 42 dpi; a receptor layer is the next step. |
| Treg stability (branch 3B, paused) | none from this series | Tregs are not resolvable here (one 176-cell candidate cluster that the labels contradict); see Stage 2. |

## Stage 2. Proposed and not run

Branch labels follow the reading order of 2026-09-15
([`Thesis/README.md`](../README.md)): 2N and 2W are the active
branches, 3A and 3B are paused until S1 or D1 returns a result. W1 is read
with papers 15 and 16; Nb1 with papers 3, 6 and 14.

| Id | Branch | Proposal | Needs | Stop condition |
|---|---|---|---|---|
| W1 | 2W | Pseudobulk the myeloid subset per animal by phase and define the reconstituting aMAC programme (lipid and surfactant catabolism, Krt79); compare with the ARG1 and ornithine axis of Yadav et al. 2025 [DOI](https://doi.org/10.1172/JCI188734) | script 11 outputs; per-animal counts | fewer than 3 animals per pooled phase |
| S1 | 3A | Receptor expression (Csf2rb, Il10ra, Il10rb, Tgfbr1, Tgfbr2, Il4ra) on the reconstituting subclusters per phase; STAT3-like protective versus STAT1-like inflammatory programme scores in epithelium per animal and phase (gate 3A condition) | pseudobulk from W1 | no programme separates phases at the animal level |
| Nb1 | 2N | Ligand and receptor candidates between AT2 (cluster 10) and the fibroblast niche (clusters 12, 14, 22) at 11 dpi, only after the pseudobulk step (backbone step 5) | W1-style pseudobulk for epithelium and mesenchyme | candidate pairs not expressed in both partners in at least 2 animals |
| D1 | 3B | Gate 3B check on external lung-Treg single-cell series: Loffredo et al. 2025 [DOI](https://doi.org/10.1172/jci.insight.187245) (GSE292440, GSE277256, GSE277226) and McCullough et al. 2026 [DOI](https://doi.org/10.1093/jimmun/vkag119) (GSE300399) | new raw-data inventory entries | Tregs not separable from other T cells in those series |
| V1 | all | Held-out human validation of the myeloid states: Bailey et al. 2024 [DOI](https://doi.org/10.1038/s41590-024-01975-x) (GSE232628 and companions) and Li et al. 2024 [DOI](https://doi.org/10.1126/scitranslmed.adn0136) (GSE263817) | mouse-to-human ortholog panel; donor-aware pseudobulk | states do not map at the donor level |
| G1 | 2W, 3A | The gene-set layer of W1: per-animal pseudobulk of the myeloid and the capillary compartments by phase (42 dpi n = 8, 90 dpi n = 4, 366 dpi n = 3; active repair pooled across 6 to 25 dpi n = 8; baseline a reference band only), pre-ranked GSEA on the per-animal ranking against MSigDB mouse hallmark and GO biological process plus three repository sets (the aMAC reconstitution programme, the iCAP markers, the ARG1 and ornithine circuit of paper 16), all frozen and hashed before the ranking is opened; a matched-size random-gene-set null beside the permutation null | W1 pseudobulk; the MSigDB download recorded as an input with its hash; gseapy or an in-house pre-ranked implementation | fewer than 3 animals in any tested arm after the 50-cell floor; a gene set added after ranking; the positive control fails (the hallmark G2M checkpoint set must clear in active repair against injury resolution, which the Ki67 trace of trial N1 already shows at the animal level) | **Run 2026-09-21.** The G2M gate failed in both compartments (C156); trial G1b replaced it with repository injury-state marker sets at four against three animals, both clear (C157); read through that gate, the myeloid TEST contrast gives 2 DNA-replication GO programmes still higher at 42 and 90 dpi than at a year (C158) and the capillary contrast gives nothing (C159). [`trials/g1_gsea_by_phase/g1_summary.md`](trials/g1_gsea_by_phase/g1_summary.md), [`trials/g1b_corrected_positive_control/g1b_summary.md`](trials/g1b_corrected_positive_control/g1b_summary.md) |
| G2 | 2W, 2N | IPF against control per donor within compartments (AT2 and transitional epithelium, fibroblasts, macrophages) in GSE136831 (32 against 28 donors; COPD out of scope as in E2), the same statistic as G1 with the human hallmark and GO collections plus the AREG-EGFR and niche sets of A2, replicated in GSE135893 (12 against 10 donors) as the held-out cohort | donor pseudobulk streamed from the raw sparse counts; the 50-cell floor per donor and compartment of E6 | a set that clears in GSE136831 and not in GSE135893 is Not established; the positive control fails (the fibrotic fibroblast set of Tsukui, reference M6, must clear in fibroblasts) | **Run 2026-09-21.** Positive control amended before the run to the EMT hallmark (the Tsukui set is not on disk as a list). G2 ran the same day: 254 discovery sets cleared both nulls and 62 replicated (rows C160, C161). [`trials/g2_gsea_ipf/g2_summary.md`](trials/g2_gsea_ipf/g2_summary.md) |

**Trial G0 (run 2026-09-21)** asked, before any of this is attempted, on which of the
16 imported deposits a ranked-list GSEA is admissible at all under the unit rule:
4 (GSE262927, GSE136831, GSE135893, GSE131907), with GSE309751 marginal and 11 not
admissible, eight of them one library per condition. The scan, its frozen rules and the
reading are in [`trials/g0_gsea_feasibility/g0_summary.md`](trials/g0_gsea_feasibility/g0_summary.md)
(register row C155). G1 and G2 are the proposals it licenses; GSE131907, a paired
epithelial contrast already answered at the gene level (C40), is a third option and not proposed.

Papers above were identified through PubMed on 2026-09-10; GEO accessions
through the NCBI link service on the same day.
