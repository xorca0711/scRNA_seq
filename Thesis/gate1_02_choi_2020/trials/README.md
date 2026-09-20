# Trials index: Choi 2020 (roadmap paper 2)

Every trial of this paper sits in one flat folder. This page is the index:
the question each trial asked, what it returned, and where the numbers are.
The full pre-registration for any trial is the docstring at the top of its
script; the narrative and the outcomes are in
[`../ANALYSIS_TRIAL_PLAN.md`](../ANALYSIS_TRIAL_PLAN.md) and the study note
is [`../README.md`](../README.md).

**Two gate namespaces.** The trial gates here (0, 1, 2, 3) are this folder's
analysis gates. The reading-order gates (1, 2C, 2N, 2W, 3A, 3B) live in
[`../../README.md`](../../README.md); this paper is reading-order Gate 1,
paper 2.

**Two branches sit beside this folder** and use the same conventions:
[`../datp_epigenetics/trials/`](../datp_epigenetics/trials/README.md) for the
transitional state in chromatin, and
[`../axin2_il1r1/trials/`](../axin2_il1r1/trials/README.md) for this paper's own
closing Discussion question. Both are on other people's data, because this
paper's ATAC deposit is coverage tracks only.

**The D series** works on the paper's own deposit (GSE145031, GSE144468). A
trailing letter (D2b, D5b) marks a disclosed corrected pass over an earlier
trial, never a replacement of it: the first outcome stays in the record and
no threshold moves. **The E series** is proposed and not run.

**A warning about two stored readings.** Trials D2 and D5 each wrote a
verdict that its own rule produced and that its own table shows to be an
unfair description of the data: every cluster in a lineage-sorted or
organoid AT2 library detects Sftpc, so the rule's Sftpc-low clauses never
fired. The thresholds were left where they were, the defect is disclosed
(rows C86 and C95), and the corrected readings are in D2b and D5b.

| Trial | Question | What it returned |
|---|---|---|
| [D0](d0_data_reality_check.py) | What does the deposit contain, before anything is fitted | One gene space; six of eight matrices are raw whitelists; no counted reporter; **one library per condition, two mice pooled each, so nothing is testable**; the ATAC deposit is coverage only |
| [D1](d1_cells_and_qc.py) | How many cells under the paper's own filter, their quality, their doublets | 11,683 Tomato-positive cells against the paper's 12,086; the mitochondrial bound is the one that binds; Scrublet flags 0 to 3 per lineage library |
| [D2](d2_state_recovery.py) | Are the five states recoverable blind, and does composition follow the time course | **Not recovered** (four states, AT1 missed); the rule's Sftpc clauses could never fire, disclosed, outcome kept |
| [D2b](d2b_corrected_annotation.py) | The same readings with the Sftpc clauses applied only where Sftpc separates clusters | hAT2, cAT2, DATP and AT1 recovered; DATP 0.3, 18.2 and 6.3 percent at PBS, day 14 and day 28; day-14 cycling AT2 7.1 percent; **primed AT2 assigned nowhere**, not at three resolutions, not in DATP sub-clusters, 18 of 9,546 cells at the cell level |
| [D3](d3_ordering.py) | Does AT2 reach AT1 through pAT2 and DATP | O1 to O6 not computable without pAT2; hAT2 < DATP < AT1 in pseudotime at both k, pooled and inside day 14; DATP to AT1 the strongest PAGA edge |
| [D4](d4_composition_and_controls.py) | Do the readings survive the calling rule, a depth window, the lineage label | Directions hold inside 1,000 to 2,000 genes; the calling rule flagged day 14 by diluting its own denominator (disclosed); the day-14 Tomato-negative library has no alveolar cluster, so that reading is unreadable |
| [D5](d5_organoids.py) | Do the organoids carry the states and does IL-1beta move them | Three states; the stromal cluster and a 481-cell AT1 cluster missed by the same defect; outcome kept |
| [D5b](d5b_corrected_annotation.py) | D5 under the corrected rule | AT1 recovered, stromal cluster removed, epithelial counts within 8 percent of the paper's; IL-1beta moves the epithelium from 78 percent hAT2 to 87 percent of a DATP-marker cluster with identity genes intact; the Figure 7 AT1 marker reading holds on two of three sets (three inside the window) |
| [D6](d6_programmes.py) | Are the DATP programmes highest in DATP | In vivo yes for p53, arrest, hypoxia, interferon-gamma and glycolysis; Hif1a detection highest in cycling AT2; the Ndrg1-free hypoxia edge is a depth reading; organoids hypoxia only |
| [D7](d7_attack_claims.py) | Is pAT2 a depth artefact (A1), DATP a doublet cluster (A2), a dissociation signature (A3), a discrete co-expressing population (A4) | A1 not computable; **A2 no** (Scrublet 0.05, co-detection intermediate); A3 not attempted (list not on disk); A4 saturated (Krt8 in 90 percent of cells, ratio 1.05) |

## Shared code

[`choi_utils.py`](choi_utils.py) holds the library manifest, the paper's
filter, the marker sets, the frozen annotation thresholds, the first-pass
rule `annotate_clusters` and the corrected rule `annotate_clusters_b`. The
generic deposit readers come from
[`../../gate1_04_sikkema_2023_hlca/trials/trial_utils.py`](../../gate1_04_sikkema_2023_hlca/trials/trial_utils.py).
Regenerable objects live under `raw_data/GSE145031/choi_trials/` and
`raw_data/GSE144468/choi_trials/` (gitignored). Run order: D1, D2, D2b, D3,
D4, D5, D5b, D6, D7, in the x86-64 environment.
