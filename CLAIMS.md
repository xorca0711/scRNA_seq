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

## Stage 3. Cardoso 2026, the Gate 2 paper (2026-09-12)

A different deposit with a different ceiling. Every mouse library of this
paper pools three mice and each genotype contributes one library per sort, so
**no genotype contrast in it carries within-group replication**: nothing in
this section is tested, no P value is computed, and the highest status any row
can reach is Descriptive only. The companion series GSE247505 is the
exception, with two libraries per arm and three time points. Trials and their
frozen rules: [`Thesis/gate2_05_cardoso_2026/ANALYSIS_TRIAL_PLAN.md`](Thesis/gate2_05_cardoso_2026/ANALYSIS_TRIAL_PLAN.md).
Owner review pending on every row.

| # | Claim | Analyses behind it | Artefact | Status | Potential |
|---|---|---|---|---|---|
| C19 | All five accessions are public and complete: 30 libraries, 123,807 barcodes, three gene spaces, one non-gene feature | trial C0: every library opened, design parsed from the GEO SOFT files | `trials/c0_data_reality_check/c0_library_inventory.csv` | Descriptive only | The precondition for everything else; also the row that names GSE247505, which the paper's data-availability statement omits. |
| C20 | No genotype contrast in the deposit has within-group biological replication | trial C0: deposited design, three mice pooled per library | `trials/c0_data_reality_check/c0_deposited_design.csv` | Descriptive only | High as a methods statement. It is the reason every row below is Descriptive only, and it is the kind of constraint a reviewer checks first. |
| C21 | The fibroblast-before-macrophage ordering cannot be tested from this deposit: every mesenchymal and immune library is one time point | trial C0 | same | Not established, and not establishable from these data | The honest closure of a branch, which is worth more than a weak answer to it. |
| C22 | The published reprogrammed fibroblast subset is not recovered by the pre-registered rule; the population is nonetheless present as clusters 14 and 16 (1,080 cells, about 99% Red2Kras) | trials C1 and C1b: blind clustering, the paper's marker sets as the caller | `trials/c1_fibroblast_compartment/`, `trials/c1b_characterise_red2kras_private/` | Descriptive only; the wording is an owner decision | Moderate. The interesting half is why the rule failed, not that it did. |
| C23 | The paper's reprogrammed-fibroblast marker set is maximised by smooth muscle, because Acta2 and Pdgfrb are mural markers; Runx1 additionally flags myeloid cells | trial C1b Question C | `trials/c1b_characterise_red2kras_private/c1b_question_c_gated_fibroblast_subclusters.csv` | Descriptive only | Yes. A reusable caution about a published marker set, and the kind of finding a marker-based pipeline is built to produce. |
| C24 | Cluster 14 misses the pre-registered Tnc floor by 0.001 and the floor was not moved | trial C1b Question A | `trials/c1b_characterise_red2kras_private/c1b_question_a_published_signature.csv` | Descriptive only (disclosed near-miss) | As a curation-hygiene example, the counterpart of C7. |
| C25 | Fibrotic and inflammatory markers mark separate cells at 2 weeks (9.7% double-positive against 8.2% expected under independence) | trial C1c: per-cell co-detection | `trials/c1c_fibrotic_inflammatory_overlap/c1c_overlap_by_group.csv` | Descriptive only | Moderate. It supports a paper claim by a route the paper did not use. |
| C26 | The mesenchymal sort carries 6.5% off-target cells, including 184 mutant epithelial cells that are 88% Areg-positive, almost entirely in the tumour arm | trial C1d: compartment-marker detection per cluster | `trials/c1d_sort_purity/c1d_sort_purity_by_cluster.csv` | Descriptive only | High, and practical. Anyone computing ligand-receptor signalling inside that library alone would read a contaminant as the source. |
| C27 | The paper's epithelial composition reproduces blind: DATP-like 46.6% to 22.0% and AT2 29.2% to 59.1% on Areg deletion, against the paper's 50.1 to 25.9 and 21.3 to 54.6 | trial C2: blind clustering of the RFP+ sort, labels never used | `trials/c2_areg_deletion_arm/c2_epithelial_composition.csv` | Descriptive only | High. The closest thing to a validation this deposit allows, and the figure a reader will look at first. |
| C28 | Four of five pre-registered directions are met; the fifth is unscorable because the Cd177-positive state is not resolved | trials C2 and C2b | `trials/c2b_composition_without_the_confidence_floor/c2b_summary.md` | Descriptive only (rule revision disclosed) | Moderate, as the record of a metric that failed the same way twice and was disclosed both times. |
| C29 | On Areg deletion, Tnc, Acta2, Fst and Runx2 fall while Pdgfrb and Runx1 do not | trial C2 Part B: detection fractions within the same series | `trials/c2_areg_deletion_arm/c2_fibrotic_programme_survival.csv` | Exploratory | The best new lead in this stage: the reprogrammed state may have an Areg-dependent matrix half and an Areg-independent Pdgfrb/Runx1 half. Needs a depth-matched control before it is more. |
| C31 | Areg is higher in the DATP-like state than in AT2 cells in all four mutant libraries, with the state defined by this repository's clustering | trial C3: blind clustering of 33,217 cells, ten libraries, two replicates per arm | `trials/c3_areg_state_specificity/c3_T1_state_specificity.csv` | Descriptive only | High, and the best-replicated row in this stage. It is the paper's central expression claim, re-derived without its labels, with a within-animal control. |
| C32 | Areg is the top EGFR ligand in that state in every mutant library, and the order Areg > Hbegf > Ereg > Tgfa is identical in all four | trial C3 | `trials/c3_areg_state_specificity/c3_T2_ligand_ranking.csv` | Descriptive only | High. The paper's top hit survives a change in how the state is defined, which is the sensitivity the owner asked for. |
| C33 | Hbegf ranks second, ahead of Ereg, which is the ligand the paper followed into culture | trial C3 | same | Exploratory | The most testable lead of the session: a second EGFR ligand the paper did not pursue. |
| C34 | The DATP-like state is nearly absent from wild-type clones of the same animals (0.08% at 4 days, 1.1% at 2 weeks) | trial C3 T3, the Red2Onco within-animal control | `trials/c3_areg_state_specificity/c3_T3_internal_control.csv` | Descriptive only | Moderate. A cleaner control than any between-genotype comparison this deposit offers. |
| C35 | Replicate libraries mix within every arm (1.53 to 1.80 against a threshold of 2.0), so no batch correction is applied | trial C3, the repository batch rule exercised where replicates exist | `trials/c3_areg_state_specificity/c3_run_record.json` | Descriptive only (decision record) | Moderate, as the counterpart of C5 in a third dataset. |
| C36 | Whether fibroblasts receive the Areg signal | no analysis; CellChat is R-only and unavailable here | none | Not established | None until a communication analysis is possible; stated so that no reader takes C31 to C33 for a communication result. |
| C30 | Mesothelial cells are depleted 10.4-fold on Areg deletion | trial C2b cluster depletion | `trials/c2b_composition_without_the_confidence_floor/c2b_cluster_depletion.csv` | Exploratory | Moderate. The paper reports mesothelial-like cells as tumour-enriched but does not test their Areg dependence. |

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
