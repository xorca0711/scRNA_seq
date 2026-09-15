# Trials index: Cardoso 2026 and its public-data extensions

Every trial of this paper sits in one flat folder, which is too large to
navigate by filename. This page is the index. Each row gives the question the trial asked,
what it returned, and where the numbers are. The full pre-registration for any
trial is the docstring at the top of its script; the narrative is
[`../ANALYSIS_TRIAL_PLAN.md`](../ANALYSIS_TRIAL_PLAN.md) and the reasoning is
[`../DIVERGENCES_AND_NEXT.md`](../DIVERGENCES_AND_NEXT.md).

**Three naming conventions.** The **C series** works on the Cardoso deposit
itself. The **E series** leaves it for public datasets, because that deposit
carries one library per genotype and cannot support a tested comparison. A
trailing letter (C1b, C2b, E1b, E2b) marks a disclosed corrected or completing
pass over an earlier trial, never a replacement of it: the first outcome stays
in the record.

**A warning about four stored readings.** Trials C7, C8, C9 and C10 each wrote
a verdict that its own rules produced and that later analysis showed to be an
unfair description of the data. In every case the threshold was left where it
was and the defect was disclosed rather than repaired. The column below gives
the corrected reading; the run record gives the original. Rows C54, C56, C72
and C79 of the claims register carry the detail.

## The C series, on the Cardoso deposit

| Trial | Question | What it returned |
|---|---|---|
| [C0](c0_data_reality_check.py) | What does the deposit contain, before anything is fitted | Five accessions, thirty libraries, three gene spaces. **Every mouse library pools three mice and each genotype gives one library, so nothing in the deposit is testable.** Found the unnamed companion series the paper's communication analysis depends on |
| [C1](c1_fibroblast_compartment.py) | Does a blind pipeline recover the published fibrotic fibroblast subset | **Not recovered**, and the pre-registered rule was at fault in four separable ways, all disclosed |
| [C1b](c1b_characterise_red2kras_private.py) | What is there instead | No cluster qualifies; one misses the Tnc floor by 0.001 and the floor was not moved. Two of the paper's six markers are mural, so the score peaks on smooth muscle |
| [C1c](c1c_fibrotic_inflammatory_overlap.py) | Are the fibrotic and inflammatory programmes in the same cells | Separate cells: 9.65 per cent double-positive against 8.2 expected, ratio 1.18 |
| [C1d](c1d_sort_purity.py) | What does the mesenchymal sort actually contain | 6.5 per cent off-target, including **184 mutant epithelial cells that are 88 per cent Areg-positive** |
| [C2](c2_areg_deletion_arm.py) | What does deleting the ligand do to the niche | The paper's epithelial collapse reproduces blind: 46.6 to 22.0 per cent and 29.2 to 59.1 |
| [C2b](c2b_composition_without_the_confidence_floor.py) | The two directions C2's confidence floor scored as nothing | Recovered from C2's own tracked tables; one direction met, one unscorable |
| [C3](c3_areg_state_specificity.py) | Is the ligand a property of the state, across the companion series | Yes in all four mutant libraries, with the four-ligand order identical in all four. **Not a communication analysis, and says so in its own output** |
| [C5](c5_figures_for_the_three_findings.py) | Draw the three findings | Three figures, and three corrections that drawing them forced |
| [C6](c6_who_makes_egfr_ligands.py) | Which compartment makes each EGFR ligand | Endothelium and mesenchyme survive the deletion; myeloid compartments do not. Found residual ligand transcript in 59 per cent of deleted cells and an ambient-RNA term |
| [C7](c7_what_the_sort_contaminant_is.py) | Which epithelial state is the sort contaminant | **A mixture**, four parts regenerative-like to four parts AT2. Its own margin refuted its method: 0.415 of rho between compartments, 0.0056 between states |
| [C8](c8_subpopulation_or_gradient.py) | Is the ligand-independent tier a population or a gradient | **A gradient.** The mixture test that disagreed is uninformative, because a third of the cells detect neither gene |
| [C9](c9_the_fst_runx2_population.py) | Pre-registered: is there a co-expressing Fst and Runx2 population | **The existence replicates in two animals, the size does not.** The marker signature replicates in direction for 26 of 27 genes; its magnitude criterion measured library depth |
| [C10](c10_published_state_or_not.py) | Are those cells the published pathological fibroblast | **Yes, and the lead closes.** Cthrc1 at 44.9 against 2.0 per cent; it had missed C9's top thirty by a ranking cutoff |
| [C11](c11_figures_for_the_contradictions.py) | Draw every result that contradicts an earlier trial | Six figures, and three more overstatements corrected in the drawing |
| [C12](c12_cellchatdb_full_resource_scan.py) | Does anything in the whole CellChatDB outrank the axis the paper followed | **The ranking is unreadable by its own guard**: 11 of the top 15 pairs target CD44. Holding the receptor constant, AREG is first of the five EGFR ligands across 26 donors |

**Closed without running, and the reason is the result.** Gate 2 branch (d),
the transcriptomic ordering of fibroblast against macrophage change, is closed
because every mesenchymal and immune library in the deposit is a single time
point. The paper's ordering claim rests on imaging, and a reanalysis cannot
test it.

## The E series, on public datasets

| Trial | Dataset | What it returned |
|---|---|---|
| [E1](e1_human_luad_ligand_sources.py) | GSE131907, human adenocarcinoma | **The only tested claim in this folder**: AREG higher in epithelium than myeloid within donor, 11 donors, paired p = 0.0020 |
| [E1b](e1b_human_luad_by_subtype.py) | the same, by deposited subtype | Qualifies E1: dendritic cells match or exceed the tumour states. The within-tumour comparison has zero donor pairs, and the floor was not lowered |
| [E2, E3](e2_human_fibrosis_ligand_sources.py) | GSE136831, GSE135893 | The myeloid prediction failed to detect at seven donors. The transitional-state prediction clears no floor in either cohort |
| [E2b](e2b_t4_disease_against_control.py) | E2's own tracked table | The reading E2 declared and never computed. Neither ligand rises clearly with disease |
| [E4](e4_bleomycin_fibrotic_genes.py) | GSE132771, bleomycin | **Refuted claim C29's reading** under a rule frozen before the matrices were opened: five of six genes rise with injury alone |
| [E5](e5_figure_for_the_refutation.py) | C5 and E4 tables | One figure for that refutation |
| [E6](e6_donor_level_axis_coupling.py) | GSE136831, 22 donors | **No donor-level coupling detectable**, and the only significant correlation was a control pair whose variables both track depth |

## The D0 trial, relocated from the withdrawn paper-2 folder (2026-09-15)

| Trial | Dataset | What it returned |
|---|---|---|
| [D0](d0_data_reality_check.py) | GSE145031, GSE144468 (Choi 2020) | **One library per condition**, the same ceiling as this deposit; six of eight matrices are raw barcode whitelists; the tdTomato reporter is not a counted feature; the ATAC accession is coverage only. Ran 2026-09-13 as roadmap paper 2's Gate 0 and moved here when that folder was withdrawn pending the owner's reading. Keeps its identifier so its run record and claims rows C58 to C64 stay true |

## Shared code

| File | Role |
|---|---|
| [`cardoso_utils.py`](cardoso_utils.py) | re-export of the repository's shared helpers; the generic readers moved to `gate1_04_sikkema_2023_hlca/trials/trial_utils.py` on 2026-09-13 |
| [`mtx_stream.py`](mtx_stream.py) | one streaming pass over a MatrixMarket file, for matrices too large to load |
| [`choi_utils.py`](choi_utils.py) | the Choi 2020 deposit's library map and its raw-whitelist rule, used by D0 only; moved here with it on 2026-09-15 |
| [`viz_style.py`](viz_style.py) | reads the validated palette from `analysis/config/palette.json` |

## Conventions every trial in this folder follows

- Frozen rules are written to a run record **before** any data is read, and the
  record is completed afterwards with inputs, sizes, modification times,
  package versions and outputs.
- A rule is never moved after its result is seen. Where a rule was wrong, the
  first outcome stays and a corrected pass sits beside it.
- The statistical unit is the mouse or the donor, never the cell, for anything
  group-level. No trial on the Cardoso deposit itself computes a P value.
- Author labels are held out of every fitting step and used only to grade.
- Regenerable intermediates (`.h5ad`, `.npz`) live under `raw_data/`, which is
  gitignored. Every tracked artefact is a table, a figure or a run record.
