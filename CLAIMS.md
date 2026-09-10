# Claims register

One row per claim this repository makes, with the analyses that stand behind
it, the artefact a reader can open, its status in the root vocabulary
(Validated, Descriptive only, Exploratory, Retracted-superseded, Not
established), and a plain judgement of whether the result has potential: a
figure or sentence a target lab would care about, a lead that needs one more
analysis, or a dead end kept for the record. A claim with no logged artefact
is Not established by definition. Status changes only on the owner's
retain/reject decision (DEVELOPMENT.md); rejected rows stay in the table.

Numbers quoted here come from the artefact named in the row. The stepwise
account of how the analyses relate to each other, and to the target labs, is
[`Thesis/gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md`](Thesis/gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md).

## Stage 0. Initial run (August 2026): recovering the published biology from raw counts

| # | Claim | Analyses behind it | Artefact | Status | Potential |
|---|---|---|---|---|---|
| C1 | Blind clustering recovers the deposited mouse cell types (median purity 0.947 over 107,626 labelled cells) | per-sample QC, Scrublet per capture, no batch correction, Leiden 0.3; labels held out and graded afterwards | `analysis/GSE262927/tables/cluster_vs_author_celltype_fraction.csv`; checked by `validate_repository.py` | Validated | Yes. The annotation-robustness credential every lab reads first; the atlas UMAP is already the README figure. |
| C3 | The injury capillary state has not resolved by 366 dpi (2.0% at baseline, 37.5% at 25 dpi, 21.7% at 366 dpi) | capillary subset reclustered within compartment; per-animal abundance by day | `analysis/GSE262927/regeneration_focus/tables/icap_abundance_per_sample.csv` | Validated | Moderate. It carries the "persistent injury state" theme, but no target lab works on endothelium; use it as the template for the interstitial-macrophage lead (C13). |
| C4 | The injury capillary state has a CAP1 origin (Kit line traced at 33 to 53% per animal); the CAP2 lines are uninformative, not negative | 8-sample Cre cohort; trace rule reproduces the deposited calls at 100% | `analysis/GSE262927/lineage_tracing_cohort/` | Validated (Kit); Not established (CAP2) | Moderate. A lineage-trace method credential; low direct pull for the four labs. |
| C5 | Batch correction is decided per dataset: none for the mouse (replicates within a group mix), Harmony for the human series (donor-private fragmentation) | measured within-group replicate enrichment; explicit logged override for GSE178360 | `analysis/GSE262927/qc/batch_assessment.json`, `analysis/GSE178360/qc/` | Validated (as a decision record) | Yes. The design-decides-integration argument is the Wagner-facing method statement. |
| C6 | The human distal-lung series carries an AT0-like minority; the earlier "AT0 candidate subcluster" is mostly AT2 or uncertain | HLCA marker transfer (trial S3) and scArches mapping (trial S2, 99.2% agreement with the HLCA authors' own transfer) | `Thesis/gate1_04_sikkema_2023_hlca/trials/s2_reference_mapping/` (merged 2026-09-10) | Descriptive only; headline re-wording pending (PROGRESS item 15) | Moderate. The reference-mapping credential is sound; the AT0 claim itself is being demoted, which is the honest outcome; the reference-aligned human epithelial panels are displaced to [`archive/DISPLACED.md`](archive/DISPLACED.md). |
| C7 | Scrublet does not over-remove AT0-like cells (3.9% flagged against a 6.3% baseline); the earlier claim that it did is refuted | first co-expression gate, then a stricter lineage-negative gate | `docs/DOUBLETS_AND_SCRUBLET.md` | Retracted-superseded (refuted, kept on display) | Yes, as a curation-hygiene example: a co-expression gate cannot audit a co-expression detector. |
| C8 | A marker-panel annotator is contradicted by the deposited labels in 3 of 29 clusters; cluster 0 is interferon-responding CAP1, not transitional epithelium | blind panel scoring, then grading | `analysis/GSE262927/tables/cluster_annotation_proposals.csv` | Validated (negative result) | Yes, as the "state fooled a type classifier" example for the robustness pitch. |

## Stage 1. PI-matched follow-ups (2026-09-10): the source paper's phase and myeloid claims

Rules frozen before each run; unit is the animal; medians; no P values (two
animals per active-repair day). Owner review pending on every row.

| # | Claim | Analyses behind it | Artefact | Status | Potential |
|---|---|---|---|---|---|
| C9 | Regeneration proceeds in the paper's three proliferative phases: myeloid cells peak at 6 dpi, epithelium and mesenchyme at 11, endothelium at 19 (4 of 5 lineages in the paper's window) | script 10: Ki67 trace per lineage in each cohort's immediate window; peak windows pre-registered | `analysis/GSE262927/phase_timecourse/tables/proliferation_peak_by_lineage.csv` | Descriptive only | Yes. The backbone figure for all four labs; it also fixes 11 dpi as the time point for any stem-niche question (Nabhan). |
| C10 | Lymphoid proliferation peaks at 11 dpi, one harvest later than the paper's immune window | same run as C9 | same table | Descriptive only (deviation from the paper) | Low. Adaptive cells lagging myeloid cells is expected; a footnote, not a lead. |
| C11 | The deposited cell-cycle call does not measure proliferation here (it calls about 40% of cells and most lymphocytes cycling; agrees with the trace for 1 of 5 lineages) | same run as C9, cross-check | `analysis/GSE262927/phase_timecourse/tables/cycling_fraction_per_animal.csv` | Not established as a measure | None as a result; keep as the caveat that stops anyone quoting the classifier. |
| C12 | Alveolar macrophages fall from 30.8% to 4.6% of myeloid cells at 6 dpi and rebuild to 49.0% by 42 dpi while inflammatory monocytes rise from 1.9% to 56.0% and return to 1.8% (Figure 3 of the paper) | script 11: blind re-embedding of atlas clusters 5, 17, 24; labels held out; per-animal composition by day | `analysis/GSE262927/myeloid_focus/tables/myeloid_label_composition_median_by_dpi.csv` | Descriptive only | High. The Wagner figure (monocyte and macrophage states on a per-animal timeline) and the Saxton resolution sentence. |
| C13 | Interstitial macrophages keep rising through 42 to 90 dpi (2.8% at baseline, 10.3% at 42, 14.0% at 90, 7.7% at 366) instead of resolving | same run as C12 | same table | Exploratory (three to four animals) | Yes, the best new lead: a candidate persistent myeloid state paralleling C3. Needs the per-animal pseudobulk programme (Stage 2, W1) and a human check against the profibrotic monocyte-derived macrophages of post-COVID lungs. |
| C14 | The myeloid compartment resolves into 16 label-coherent subclusters at Leiden 0.5 (87% of labelled cells in pure subclusters; 48% at 0.2, 92% at 1.0) | script 11 grading | `analysis/GSE262927/myeloid_focus/tables/myeloid_grading_sub_r0.5.csv` | Descriptive only | Moderate. Extends the robustness credential to the compartment where trial S5 had found the atlas's worst mixture. |
| C15 | The rebuilt alveolar macrophage pool is labelled mainly by the 2 to 3 dpi tamoxifen window (median traced fraction at 42 dpi 79.7%, then 60.0, 29.3, 42.0 for the later windows) | script 12: trace by window at the common harvest; main window pre-registered | `analysis/GSE262927/myeloid_focus/amac_origin/tables/window_contribution_summary.csv` | Descriptive only | Yes. The Wagner pitch sentence: origin of a reconstituted pool read from a design variable rather than a fitted trajectory. |
| C16 | The aMAC labelling exceeds what marrow-progenitor inheritance would give | script 12: same-animal cMON and neutrophil reference, 30-cell floor | `analysis/GSE262927/myeloid_focus/amac_origin/tables/marrow_reference_check.csv` | Not established (references evaluable in 3 of 8 animals) | Low with this series: classical monocytes are too sparse at 42 dpi. Testable only in a monocyte-rich dataset. |
| C17 | The rebuilt pool has two sources visible as a within-aMAC split in traced fraction (subcluster 2 more traced than 3 in early windows, reverse in late windows) | script 12: per-subcluster trace at 42 dpi | `analysis/GSE262927/myeloid_focus/amac_origin/tables/amac_subcluster_traced_fraction_42dpi.csv` | Not established (evaluable in one animal per window) | Moderate. The direction is consistent across the four evaluable animals; W1 markers (Krt79 mature versus inflammatory aMAC) would say whether the split is the paper's aMAC_a and aMAC_b. |
| C18 | The 6 dpi inflammatory-monocyte state is not a batch island: the two infection rounds mix within every tested day (same-round kNN enrichment 1.11 to 1.46, threshold 2.0), and after Harmony on round the state remains one subcluster with 92% of its cells from 6 dpi holding 77% and 74% of the two animals' iMON cells | script 13: uncorrected versus Harmony-on-round embeddings on the days that carry both rounds; survival criteria pre-registered, with a disclosed post hoc definition of the state (the frozen one selected the 11 to 19 dpi monocyte state) | `analysis/GSE262927/myeloid_focus/batch_sensitivity/tables/imon_state_check.csv` | Descriptive only (rule revision disclosed) | Robustness statement for C12 and C15; the Wagner-facing integration argument. Not a pitch item on its own. |

## Stage 2. Proposed, not run

No claims. The proposals, what each needs, and its stop condition are in the
trial plan (W1 Wagner pseudobulk programme, S1 Saxton receptor and STAT
layers, Nb1 Nabhan niche candidates at 11 dpi, D1 DuPage gate check on
external Treg series, V1 human validation of the myeloid states). The lead
with the most potential is C13, because it turns the repository's
persistent-state theme toward the immune compartment the target labs work
on; the one with the least is C16, which this series cannot settle.

## How to keep this register honest

- Add a row only when a run record or decisions log exists for it.
- Quote numbers from the artefact in the row, never from memory.
- When the owner rejects a row, change its status to Retracted-superseded and
  leave it in place; record the rejection in DEVELOPMENT.md.
- A row's Potential is a judgement, and says so; it does not change the status.

## Displaced (established elsewhere)

Rows moved out of the Stage 0 table on 2026-09-10 because the result is
established outside this repository (the owner's G-SURF submission, weighted
toward the Krt8-high transitional state). The artefacts and scripts stay in
place under `analysis/` and remain validated; the narrative is not extended
here. The row keeps its number so that cross-references elsewhere still
resolve.

| # | Claim | Analyses behind it | Artefact | Status | Potential |
|---|---|---|---|---|---|
| C2 | The AT2 to Krt8-high transitional to AT1 ordering is recovered, and the transitional state is a true intermediate (median per-animal 27.4% at 11 dpi, 0.3% at 366 dpi) | regeneration focus: alveolar subset, PAGA plus diffusion pseudotime rooted in AT2, labels held out | [`archive/DISPLACED.md`](archive/DISPLACED.md) (artefacts unchanged under `analysis/GSE262927/regeneration_focus/`) | Validated | Not extended here; established in the owner's G-SURF submission. The Nabhan entry in the trial plan rests on the proliferation figure, not on this trajectory. |
