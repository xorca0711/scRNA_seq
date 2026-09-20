# Analysis trial plan: Choi, Lee et al. 2020 (roadmap paper 2, Gate 1)

Entered on 2026-09-15 at the owner's direction, after the owner read the
paper (DEVELOPMENT decision 23). The withdrawn folder of 2026-09-13 held only
trial D0; that trial is kept unchanged and its identifier reused. Every trial
writes its frozen rules to a run record before it reads any data, and
completes the record afterwards with inputs, package versions and outputs. A
rule is never moved after its result is seen; where a rule turns out to be
wrong, the first outcome stays in the record and a corrected pass sits beside
it (the precedent of trials C1 and C1b in the Cardoso folder); this folder
has two such passes, D2b and D5b. Every rule also says in advance what would
make its own answer unreadable (decision 22).

Trial identifiers are the D series, for the DATP state this paper defines.
The study note is [`README.md`](README.md) and the machine-readable extract is
[`choi_2020_extracts.json`](choi_2020_extracts.json).

**The constraint that governs every trial here.** Both single-cell accessions
carry one library per condition, and each in vivo library pools two mice. So
no contrast between time points, sorts or treatments has within-group
replication: directions are described, nothing between conditions is tested,
and no P value is admissible on such a comparison. Two things escape the
constraint and are the reason this paper can be re-asked at all: an ordering
inside one library (trial D3), and a co-detection test against a null that is
permuted inside one library (trial D7 A4).

| Trial | Gate | Question | State |
|---|---|---|---|
| D0 | 0 | What does the deposit contain, before anything is fitted | complete (2026-09-13) |
| D1 | 0 | How many cells does each library carry under the paper's own filter, what is their quality, which are doublets | run 2026-09-15; 11,683 Tomato-positive cells, comparable to the paper's 12,086 |
| D2 | 1 | Are the five states recoverable blind, and does their composition follow the paper's time course | run 2026-09-15; four of five states, AT1 missed by a disclosed rule defect; outcome kept |
| D2b | 1 | D2's readings under the corrected rule, no threshold moved | run 2026-09-15; hAT2, cAT2, DATP and AT1 recovered, time course in the paper's direction; primed AT2 assigned nowhere |
| D3 | 1 | Does AT2 reach AT1 through pAT2 and DATP (PAGA, diffusion pseudotime), pooled and inside day 14 alone | run 2026-09-15; O1 to O6 not computable without pAT2; hAT2 < DATP < AT1 holds at every setting |
| D4 | 1 | Do the composition readings survive the calling rule, a depth window, and the lineage label | run 2026-09-15; directions hold inside the depth window; the calling rule measured its own denominator; day-14 Tomato-negative unreadable |
| D5 | 1 | Do the organoid libraries carry the five states, and does IL-1beta move them the way the paper reports | run 2026-09-15; three states; the stromal cluster and AT1 missed by the same defect; outcome kept |
| D5b | 1 | D5's readings under the corrected rule | run 2026-09-15; AT1 recovered, stromal cluster removed, epithelial counts within 8 percent of the paper's; the Figure 7 marker reading computable |
| D6 | 2 | Are the paper's DATP programmes highest in the DATP state, in vivo and in organoids | run 2026-09-15; in vivo yes for p53, arrest, hypoxia, interferon-gamma and glycolysis, not for Hif1a; organoids hypoxia only |
| D7 | 2 | Attack: is pAT2 a depth artefact, is DATP a doublet cluster, is DATP a discrete co-expressing population | run 2026-09-15; A1 not computable; A2 not doublet-like; A3 not attempted; A4 saturated |
| E-type | 3 | Cross-dataset checks that need data not on disk or displaced material; proposed, not run | proposed |

Gate 0 reads the deposit; Gate 1 reproduces what the paper reports from its
own deposit; Gate 2 asks what the paper did not (programmes as a whole, and
the artefact hypotheses); Gate 3 would leave the deposit. These are this
folder's trial gates and are not the reading-order gates of
[`../README.md`](../README.md); this paper sits in reading-order Gate 1.

---

## Contents

| Section | Trial | In one line |
|---|---|---|
| [Backbone check](#the-core-logic-backbone-checked-trial-by-trial) | all | which backbone step, method principle and gate pattern each trial serves |
| [D0](#d0-data-reality-check-gate-0) | data reality check | one library per condition; six raw whitelists |
| [D1](#d1-cells-quality-and-doublets-gate-0) | cells, quality, doublets | the paper's filter on every barcode |
| [D2](#d2-state-recovery-gate-1) | state recovery | five states, seven readings |
| [D2b](#d2b-corrected-annotation-gate-1) | corrected annotation | the Sftpc clauses applied only where Sftpc separates clusters |
| [D3](#d3-ordering-gate-1) | ordering | PAGA and DPT, pooled and inside day 14 |
| [D4](#d4-composition-and-its-controls-gate-1) | composition controls | calling, depth window, Tomato-negative |
| [D5](#d5-organoids-gate-1) | organoids | control against IL-1beta, direction only |
| [D5b](#d5b-organoids-under-the-corrected-rule-gate-1) | organoids, corrected | the same pass on the organoid deposit |
| [D6](#d6-programmes-gate-2) | programmes | p53, arrest, hypoxia, interferon-gamma, glycolysis, Il1r1 |
| [D7](#d7-attack-claims-gate-2) | attacks | depth, doublets, discreteness |
| [Proposed](#proposed-not-run-gate-3) | E-type | cross-dataset checks |

- The branches, and why each exists
---

## The branches, and why each exists

Two branches of this paper were opened on 2026-09-20 at the owner's direction.
Both are on other people's data, and in both cases the reason is a property of
this paper's own deposit rather than a preference.

### `datp_epigenetics/` : the transitional state in chromatin

This paper supports the epigenetic half of its Il1r1 claim with ATAC-seq, and
GSE144598 deposits two bigwig coverage tracks with no peaks and no reads, which
trial D0 recorded as unusable. The branch therefore runs on two 10x multiome
deposits and asks whether the AT2 identity programme is closed at the chromatin
level or merely silenced at the RNA level in the CLDN4-positive transitional
state.

Trials M0 to M3. One validated result, and it is a deposit error: GSE247130's
barcode suffix map is inverted relative to its GEO sample order. One
Descriptive-only positive: the labelled group loses the AT2 identity programme
in RNA by 12.6 to 17.1 detection points. One retraction. The chromatin question
is Not established, with a direction. See
[`datp_epigenetics/README.md`](datp_epigenetics/README.md).

**It also answers one of the owner's two questions, by refusing.** A
CLDN4-positive KRT8-positive transcript call cannot separate neonatal
developmental immaturity from injury-induced transition, so that marker set is
shared with normal development and is not by itself damage-associated.

### `axin2_il1r1/` : this paper's own closing Discussion question

Choi 2020 ends its Discussion by proposing a comparison of Il1r1-positive and
Axin2-positive AT2 cells. Six years on nobody has made it, and the branch
establishes why rather than pretending to make it: both populations are defined
by tamoxifen-inducible lineage reporters, no public deposit carries both
readouts in the same cells, and both transcripts sit near five per cent
detection at about one molecule per positive cell.

The branch assessed four routes and ran two of them. Route B, locus
co-accessibility, closed after three passes under a rule frozen before the first.
Route A turned out to rest on data that was never deposited, and its replacement
refused because the Il1r1 genotype labels of that deposit cannot be verified from
its own matrices. See [`axin2_il1r1/README.md`](axin2_il1r1/README.md).

**The premise needed correcting too.** Neither Axin2-positive nor Il1r1-positive
AT2 cells are established as a distinct subset; the field's own data lean toward
an injury-inducible state in both cases, and the Axin2 fraction differs 20 to 30
fold between two knock-in alleles with nobody having explained the gap in eight
years. This paper cites the subset claim as settled.

---

## The core logic backbone, checked trial by trial

The owner asked that every trial be checked against the shared portfolio
backbone (the six steps on the Notion PI Target Map), the repository's method
principles (`AI_CONTEXT.md`), and the gate pattern set by the Cardoso folder
(decisions 18 and 22). The table is the check. "Backbone step" numbers are
the six steps of the Notion page: 1 manifest, 2 QC and annotation, 3
sample-aware pseudobulk, 4 programme scoring, 5 communication only after
evidence, 6 held-out validation and publication.

| Trial | Backbone step | Method principle it enforces | Gate pattern | Admissible statistic | What would make it unreadable |
|---|---|---|---|---|---|
| D0 | 1 (manifest) | raw data determines the workflow; the replicate rule | Gate 0 before anything is fitted | none; description of the deposit | nothing; it only measures |
| D1 | 2 (QC) | per-library MAD thresholds with rationale, Scrublet per capture with the fallback logged, cell calling frozen before the matrix is opened | Gate 0 | none; counts per library | a bound that removes nothing is reported as not binding, never as applied |
| D2 | 2 (annotation) and 4 (state scoring) | labels held out (none deposited; the paper's text is the answer key), batch correction as a hypothesis test (refused here because library equals condition, enrichment reported), decisions machine-logged | Gate 1 with a frozen annotation rule and a stop-and-characterise clause (dispersed presence) | none between libraries; composition is a description of one library of two pooled mice | fewer than four states assigned makes R2 to R5 not computable |
| D2b, D5b | 2 (annotation) | a rule found wrong is not moved; the first outcome stays and a corrected pass with one disclosed added condition sits beside it (the C1 and C1b precedent) | Gate 1, corrected pass | none | the same as D2 and D5 |
| D3 | 4 (state trajectory) | the one within-library statistic the constraint allows; root chosen before pseudotime is seen; two k values pre-declared | Gate 1 | ordering of medians inside one library; no P value | any state below 30 cells in the subset |
| D4 | 2 and 3 (the sample as the unit) | depth control on every cross-library detection comparison (the E6 lesson); the calling rule tested for sensitivity; the lineage label treated as a within-experiment control | Gate 1 controls | none; descriptions with the D2 readings recomputed | the Tomato-negative gate finding no alveolar cluster |
| D5 | 2 and 4 | same rule as D2 applied unchanged to a second deposit; depth window because the two libraries differ in genes per cell (D0) | Gate 1 on the organoid deposit | none; two libraries | fewer than 30 AT1 cells in a treatment makes G5 not computable |
| D6 | 4 (programme scoring) | scores with control gene sets; a marker that sits in two sets is removed from one of them (Ndrg1); a claim that needs a list from memory is not attempted (Il1r1 by phase) | Gate 2 | none; means per state | DATP below 30 cells |
| D7 | 4 and the repository's negative-results rule | attack hypotheses frozen with their readings; depth-matched comparison (A1); the paper's own doublet cut (A2); permutation inside depth deciles (A4, the C9 method) | Gate 2 | A4 only: permutation p inside one library | A1 and A2 below floor; A4 with fewer than 30 double-positives reports existence, not magnitude |
| E-type | 6 (held-out validation) | leave the deposit to get a testable unit (decision 19) | Gate 3 | would be per-donor or per-animal | not run; see the proposals |

Backbone step 5 (communication) is not touched by this folder: the deposit
sorts epithelium only, so no receiving compartment exists here, and the
repository's rule is that communication is tested only after expression and
sample-level evidence, which this paper's deposit cannot supply beyond
description. Step 3 (pseudobulk differential expression) is not possible with
one library per condition; the animal-as-unit rule is honoured by refusing the
comparison rather than by pooling cells into a false replicate.

---

## D0. Data reality check (Gate 0)

Run on 2026-09-13; unchanged. Script
[`trials/d0_data_reality_check.py`](trials/d0_data_reality_check.py);
artefacts in [`trials/d0_data_reality_check/`](trials/d0_data_reality_check/)
and the summary in
[`trials/d0_data_reality_check/d0_summary.md`](trials/d0_data_reality_check/d0_summary.md).
One gene space across all eight libraries; no non-gene features (so the
tdTomato reporter is not counted and the sort cannot be verified from the
matrix); six of eight matrices are raw whitelists; every marker gene present;
no contrast carries within-group replication; the ATAC accession is coverage
only. Its claims are rows C58 to C64 of the register and its reporting
defect (a single boolean that read as "none are raw") is disclosed there.

---

## D1. Cells, quality and doublets (Gate 0)

Script [`trials/d1_cells_and_qc.py`](trials/d1_cells_and_qc.py); artefacts
in [`trials/d1_cells_and_qc/`](trials/d1_cells_and_qc/). Rules are in the
script docstring and the run record: the paper's own post-Cell-Ranger filter
(more than 500 and fewer than 7,000 genes, more than 2,000 UMI) applied to
every raw barcode; a plain floor (500 counts, 200 genes) for sensitivity
only; MAD bounds per library; Scrublet per capture with the fallback logged;
flagged cells kept in the object with their score.

### Outcome

Run 2026-09-15. Tomato-positive libraries: 12,865 barcodes pass the paper's
filter and 11,683 cells remain after the MAD rule and Scrublet (PBS 4,418 to
4,319; day 14 4,226 to 3,790; day 28 4,221 to 3,574), against the paper's
12,514 captured and 12,086 in Figure 1B. Tomato-negative: 6,552, 4,138 and
5,252 cells; organoids 1,868 control and 2,699 IL-1beta (the paper's 1,286
and 2,584 follow its removal of a stromal cluster, which is trial D5's job).
The binding MAD bound was the mitochondrial fraction in every library; the
count and gene bounds removed nothing except a high-count bound in the PBS
Tomato-negative library and a low-gene bound in both organoid libraries, and
the record reports the others as not binding rather than as applied.
Scrublet's automatic threshold sat between 0.36 and 0.57 in the lineage
libraries and flagged 0 to 3 doublets each; the organoid libraries
(thresholds 0.12 and 0.20) flagged 8 and 2. Median genes per cell 1,568,
2,268 and 2,089 across the three Tomato-positive libraries, a 1.4-fold spread
that trials D4 and D6 control for. Summary
[`trials/d1_cells_and_qc/d1_summary.md`](trials/d1_cells_and_qc/d1_summary.md).

---

## D2. State recovery (Gate 1)

Script [`trials/d2_state_recovery.py`](trials/d2_state_recovery.py);
artefacts in [`trials/d2_state_recovery/`](trials/d2_state_recovery/). The
annotation rule and its thresholds are in
[`trials/choi_utils.py`](trials/choi_utils.py), written before the trial ran.

### Outcome

Run 2026-09-15; outcome kept unchanged. Primary resolution Leiden 1.0 (the
resolution with the most states). States assigned: 0.3 cAT2 and hAT2; 0.5
plus pAT2; 1.0 DATP, cAT2, hAT2, pAT2. AT1 was never assigned and the
dispersed check fell under its floor (0.019 of cells with three or more AT1
canonical genes, 222 cells, against 0.02). Readings: R1 not met (four
states); R2 not met (cAT2 pooled 0.027 against 0.03 to 0.09); R3 met (PBS
hAT2 0.985); R4 not met (day-14 hAT2 0.736 against half of 0.985; pAT2
0.019; DATP 0.176); R5 not met (AT1 zero); R6 not computable; R7
same-library kNN enrichment 2.08 (library equals condition, no correction
admissible). Summary
[`trials/d2_state_recovery/d2_summary.md`](trials/d2_state_recovery/d2_summary.md).

Why the rule failed, read from its own table: every one of the fourteen
clusters detects Sftpc in 100 percent of cells (ambient surfactant transcript
on top of the sort), so the two clauses that require Sftpc detection below
0.5 (rule 1, contaminants; rule 3, AT1) could never fire. Cluster 9 (337
cells, AT1 canonical detection 0.70, AT2 identity 0.23) was labelled hAT2;
cluster 13 (40 cells, immune set 0.50) and cluster 8 (205 cells, ciliated set
0.69, club or airway 0.65) stayed in, and cluster 8 took the pAT2 label under
rule 6 because its identity detection is low for the wrong reason. The paper
removed 255 non-epithelial cells; these two clusters hold 245. The threshold
was not moved; D2b is the corrected pass. This is the fifth disclosed rule
defect of the same shape (rows C54, C56, C72 and C79) and the first outside
the Cardoso folder: a threshold was frozen on a variable that does not vary
in this deposit, which trial D0's marker table (Sftpc present in every
library) could have shown before the rule was written.

---

## D2b. Corrected annotation (Gate 1)

Script
[`trials/d2b_corrected_annotation.py`](trials/d2b_corrected_annotation.py);
artefacts in
[`trials/d2b_corrected_annotation/`](trials/d2b_corrected_annotation/). The
rule is `annotate_clusters_b` in
[`trials/choi_utils.py`](trials/choi_utils.py): D2's rules and thresholds
unchanged, with one condition on the condition, the Sftpc clauses apply only
where at least one cluster falls below the Sftpc bound. Two disclosed
additions: R2 is also reported on the day-14 library alone (the paper's
sentence sits in its day-14 paragraph), and if pAT2 is assigned nowhere the
DATP cluster is sub-clustered (Leiden 0.5) and a cell-level primed-like count
is taken among hAT2-labelled cells. The corrected labels go to a second
object; D2's is untouched.

### Outcome

Run 2026-09-15 on D2's embedding. Primary resolution 1.0. States: 0.3 cAT2,
hAT2; 0.5 cAT2, hAT2; 1.0 AT1, DATP, cAT2, hAT2. Three contaminant clusters
removed at resolution 1.0 (ciliated 205, club or airway 362, immune 40 cells;
607 against the paper's 255; the club or airway cluster, 321 of its cells
from day 28, is the one the paper's cut may have kept). Primed AT2 is
assigned nowhere: the cluster D2 called pAT2 is the ciliated cluster;
sub-clustering the 880 DATP cells gives four sub-clusters of which none meets
the pAT2 rule (identity detection 0.27 to 0.55 and inflammatory 0.12 to 0.29
against 0.86 and 0.22 in the reference; sub-cluster 2 would be hAT2 by the
rule and is kept under the DATP label because the sub-clustering looks for
pAT2 only, which the script states); and among 9,546 hAT2-labelled cells 18
(0.19 percent) detect none of Etv5, Abca3 and Cebpa while detecting at least
two of the five inflammatory genes. Readings: R1 not met (four states); R2
pooled 0.028 not met, R2 on day 14 alone 0.071 met; R3 met (0.984); R4 not
met (hAT2 0.717 against the 0.492 needed, pAT2 zero, DATP 0.182 met on its
own); R5 not met on pAT2 alone (AT1 0.029 to 0.057 and hAT2 0.717 to 0.866
rise from day 14 to day 28; cAT2 0.071 to 0.015 and DATP 0.182 to 0.063
fall, every direction the paper reports); R6 met (DATP AT1 canonical
detection 0.287 against 0.703 in the AT1 cluster, so DATP is Pdpn, Hopx and
Cav1-low as the paper defines it); R7 unchanged. Composition as fractions of
alveolar cells at PBS, day 14 and day 28: hAT2 0.984, 0.717, 0.866; cAT2
0.001, 0.071, 0.015; DATP 0.003, 0.182, 0.063; AT1 0.012, 0.029, 0.057.
Summary
[`trials/d2b_corrected_annotation/d2b_summary.md`](trials/d2b_corrected_annotation/d2b_summary.md);
figure
[`d2b_umap_states.png`](trials/d2b_corrected_annotation/d2b_umap_states.png).
Object `raw_data/GSE145031/choi_trials/tomato_annotated_d2b.h5ad`, read by
D3, D4, D6 and D7.

---

## D3. Ordering (Gate 1)

Script [`trials/d3_ordering.py`](trials/d3_ordering.py); artefacts in
[`trials/d3_ordering/`](trials/d3_ordering/).

### Outcome

Run 2026-09-15 on the D2b object. All six pre-registered readings name pAT2
and are therefore not computable, and the record says so rather than
scoring them. What the tables hold: diffusion pseudotime medians rise hAT2 <
DATP < AT1 at every setting (pooled k 15: 0.099, 0.264, 0.385; day 14 k 15:
0.157, 0.475, 0.480; pooled k 30: 0.156, 0.397, 0.587; day 14 k 30: 0.211,
0.589, 0.592). Inside day 14 the DATP and AT1 medians sit 0.005 and 0.003
apart, so that ordering is existence, not separation (106 AT1 cells). PAGA:
DATP to AT1 is the strongest edge in every graph (0.63 pooled, 0.72 day 14
at k 15); hAT2's strongest partner is DATP pooled (0.086) and cAT2 inside
day 14 (0.080, with DATP and AT1 at 0.075 each, a tie); cAT2's strongest
partner is hAT2 everywhere. The paper's chain through primed AT2 cannot be
asked without that state; the three-state ordering is the part the deposit
carries. Tables
[`trials/d3_ordering/d3_readings.csv`](trials/d3_ordering/d3_readings.csv);
figure
[`d3_dpt_by_state.png`](trials/d3_ordering/d3_dpt_by_state.png).

---

## D4. Composition and its controls (Gate 1)

Script [`trials/d4_composition_and_controls.py`](trials/d4_composition_and_controls.py);
artefacts in [`trials/d4_composition_and_controls/`](trials/d4_composition_and_controls/).

### Outcome

Run 2026-09-15 on the D2b object. C1: the rule flags the day-14 DATP fraction
as calling-sensitive (0.183 to 0.150, a move of 0.033 against 0.02) and the
PBS and day-28 fractions as not (0.001 and 0.018). The flag is a dilution:
844 day-14 barcodes pass the plain floor but not the paper's filter (median
1,108 counts), 7 of them carry three or more DATP markers, and the fraction
moves because the denominator grows. The rule measured the denominator; it is
disclosed here and in row C92, and not moved. C2: inside 1,000 to 2,000
genes (2,723; 1,379; 1,389 cells) the fractions are DATP 0.003, 0.118, 0.032
and AT1 0.013, 0.028, 0.081 at PBS, day 14 and day 28; R4 recomputed there
fails on pAT2's absence and on hAT2's day-14 fall (0.835, less than the
halving the reading asks for), R5 fails on pAT2 alone, and every direction
the paper reports for the other four states holds inside the window. C3: the
day-14 Tomato-negative library carries no alveolar cluster under the frozen
gate (4,138 cells, all non-alveolar), which the plan named in advance as
what makes this reading unreadable; PBS carries 21 alveolar cells; day 28
carries 708 (699 hAT2, 9 DATP-like, dispersed fraction 0.007 against 0.02).
Whether an unlabelled DATP-like cell came from an unlabelled AT2 cell is not
a transcriptome question. Summary
[`trials/d4_composition_and_controls/d4_summary.md`](trials/d4_composition_and_controls/d4_summary.md).

---

## D5. Organoids (Gate 1)

Script [`trials/d5_organoids.py`](trials/d5_organoids.py); artefacts in
[`trials/d5_organoids/`](trials/d5_organoids/).

### Outcome

Run 2026-09-15; outcome kept unchanged. Primary resolution 0.3; states 0.3
DATP, cAT2, hAT2; 0.5 DATP, hAT2; 1.0 DATP, cAT2, hAT2. No contaminant
cluster and no AT1 cluster were called, for the reason D2 disclosed: Sftpc is
detected in 100 percent of cells in all six clusters, so cluster 0 (586
cells, mesenchyme set 0.70) and cluster 3 (365 cells, 0.64) were kept as hAT2
and cluster 2 (481 cells, AT1 canonical 0.79) was labelled hAT2. Readings: G1
not met (three states); G2 met (control hAT2 0.969, DATP 0.009); G3 not met
(pAT2 zero; DATP 0.009 to 0.780, the direction the paper reports); G4 met
(median pseudotime 0.151 control, 0.280 IL-1beta); G5 not computable (no AT1
cells); G6 as G3 and G5. Epithelial counts 1,868 and 2,699 against the
paper's 1,286 and 2,584 after its stromal removal; the 582-cell difference in
the control library is cluster 0. D5b is the corrected pass. Summary
[`trials/d5_organoids/d5_summary.md`](trials/d5_organoids/d5_summary.md).

---

## D5b. Organoids under the corrected rule (Gate 1)

Script
[`trials/d5b_corrected_annotation.py`](trials/d5b_corrected_annotation.py);
artefacts in
[`trials/d5b_corrected_annotation/`](trials/d5b_corrected_annotation/). The
same corrected rule as D2b applied to D5's embedding, with D5's readings
recomputed by D5's own functions; the corrected labels go to a second object.

### Outcome

Run 2026-09-15. Primary resolution 0.3; states 0.3 AT1, DATP, cAT2, hAT2;
0.5 AT1, DATP, hAT2; 1.0 AT1, DATP, cAT2, hAT2. Two mesenchyme clusters
removed (675 control and 276 IL-1beta cells), leaving 1,193 control and
2,423 IL-1beta epithelial cells against the paper's 1,286 and 2,584, within 8
percent. Composition: control hAT2 0.782, AT1 0.169, cAT2 0.035, DATP 0.013;
IL-1beta DATP 0.869, AT1 0.115, cAT2 0.015, hAT2 0.001. Readings: G1 not met
(four states; pAT2 assigned nowhere, as in vivo); G2 met; G3 not met on
pAT2, with the DATP direction met (0.013 to 0.869) in the full object and
inside the window (0.012 to 0.883); G4 met (median pseudotime 0.064 control,
0.558 IL-1beta; by state, control hAT2 0.052, DATP 0.516, AT1 0.844,
IL-1beta hAT2 0.203, DATP 0.555, AT1 0.613); G5 on 202 control and 279
IL-1beta AT1 cells: early AT1 set ratio 0.96 (met; Lmo7 0.97, Pdpn 1.00,
Hopx 0.93), late AT1 0.72 (not met against 0.7; Cav1 0.58 and Spock2 0.40
fall, Aqp5 1.00 and Vegfa 0.91 do not), DATP Figure 7 genes 3.66 (met; Cldn4
3.8, AW112010 5.3, Lhfp 1.9); G6 inside 2,500 to 4,500 genes: early 0.94,
late 0.69 (met), DATP genes 3.54. The paper's Figure 7 reading holds on two
of three sets in the full object and on all three inside the window, and the
late-marker set sits on its own threshold in both.

The 2,121-cell IL-1beta cluster the rule labels DATP (identity detection
0.76, inflammatory 0.62, DATP set 0.42) is the cluster the paper calls
primed AT2 at about 77 percent. The frozen rule sees DATP markers with the
identity genes intact; the paper's identity-loss definition of primed does
not describe these cells at the cluster level. Summary
[`trials/d5b_corrected_annotation/d5b_summary.md`](trials/d5b_corrected_annotation/d5b_summary.md);
figure
[`d5b_umap_states.png`](trials/d5b_corrected_annotation/d5b_umap_states.png).
Object `raw_data/GSE144468/choi_trials/organoid_annotated_d5b.h5ad`, read by
D6.

---

## D6. Programmes (Gate 2)

Script [`trials/d6_programmes.py`](trials/d6_programmes.py); artefacts in
[`trials/d6_programmes/`](trials/d6_programmes/).

### Outcome

Run 2026-09-15 on the D2b and D5b objects. In vivo (9,546 hAT2, 313 cAT2, 880
DATP, 337 AT1): P1 met, DATP has the highest mean score for p53 (0.144),
arrest (0.846), hypoxia (0.012) and interferon-gamma response (0.295); P2
half met, glycolysis highest in DATP (0.253) but Hif1a detection highest in
cAT2 (0.476 against 0.386); P3 met in the full object by a margin of 0.004
(hypoxia without Ndrg1 -0.066 in DATP against -0.070 in hAT2); P5 inside
1,000 to 2,000 genes (5,061; 31; 213; 186 cells): P1 and glycolysis hold, P3
does not (-0.178 against -0.055), so the hypoxia edge without Ndrg1 is a
depth reading, and cAT2 is below the floor there. P4: Il1r1 detected in
0.062 of hAT2, 0.118 of cAT2, 0.061 of DATP and 0.104 of AT1; the G2/M claim
not attempted. Organoids (935 hAT2, 79 cAT2, 2,121 DATP, 481 AT1): only
hypoxia is highest in the DATP-labelled cluster (0.155; without Ndrg1 0.209,
P3 met); p53, arrest and interferon-gamma are not; glycolysis is highest in
AT1 (0.178 against 0.113) and Hif1a detection in cAT2; Il1r1 0.027, 0.101,
0.033, 0.195. The organoid DATP-labelled cluster does not carry the in vivo
programme set, a second reason to read it as the paper's primed state under
another label (D5b). Summary
[`trials/d6_programmes/d6_summary.md`](trials/d6_programmes/d6_summary.md);
figure
[`d6_programme_scores.png`](trials/d6_programmes/d6_programme_scores.png).

---

## D7. Attack claims (Gate 2)

Script [`trials/d7_attack_claims.py`](trials/d7_attack_claims.py); artefacts
in [`trials/d7_attack_claims/`](trials/d7_attack_claims/).

### Outcome

Run 2026-09-15 on the D2b object. A1 not computable: no pAT2 state (0 cells
against a floor of 30); the depth table is reported anyway (median genes per
cell AT1 1,618; hAT2 1,855; DATP 2,578; cAT2 3,093). A2: DATP is not
doublet-like; mean Scrublet score 0.052 (hAT2 0.062, AT1 0.040), no cell
above 0.7 in any state, and Sftpc-plus-two-AT1-gene co-detection 0.238 in
DATP against 0.736 in AT1 and 0.012 in hAT2, intermediate rather than
additive, which is what a transitional state looks like and a doublet does
not. A3 not attempted (the M8 list is not on disk). A4 on 3,661 day-14
alveolar cells: Krt8 detected in 0.901 and Cldn4 in 0.167, 578
double-positive against 551 expected, ratio 1.05, permutation p 0.005 within
depth deciles. The frozen reading is no evidence of a discrete co-expressing
population; the reason is that Krt8 is near-saturated, so the test cannot see
discreteness at all, and the p value reports a small excess that the ratio
band calls no evidence. The same saturation lesson as row C75. Summary
[`trials/d7_attack_claims/d7_summary.md`](trials/d7_attack_claims/d7_summary.md).

---

## Proposed, not run (Gate 3)

Each item names what it needs and why it did not run. None is started
without the owner's word; E6 is the one the trials above make pressing.

| Id | Question | Needs | Why not run |
|---|---|---|---|
| E1 | Does the human transitional epithelium of fibrosis (aberrant basaloid or KRT8-high, GSE136831, 26 donors, on disk from the Cardoso E series) carry the DATP programme set (p53, arrest, hypoxia, interferon-gamma) at the donor level | a frozen human ortholog list for the four programmes and for the DATP markers | a donor-level unit exists there, so this would be the first DATP reading with a testable design; it is a cross-species claim and belongs after the owner's review of C85 to C104 |
| E2 | Does the Krt8 transitional state of GSE262927 (Niethamer 2025) score on the DATP set, per animal and time point | nothing new on disk | displaced material: the Krt8 trajectory is established outside this repository (`archive/DISPLACED.md`); reopening it is the owner's decision, not an extension |
| E3 | Which Tomato-negative cells carry Il1b at day 14 (the paper names interstitial macrophages as the IL-1beta source) | an immune annotation rule frozen before the three Tomato-negative libraries are clustered for that purpose | D4 C3 found the day-14 Tomato-negative library almost entirely non-epithelial, which is the right library for the question; the rule was not written in this session |
| E4 | Does Il1r1 rise in G2/M-phase cycling AT2 cells (Figure S5G) | a downloaded and cited phase gene list (Tirosh et al. 2016, as scanpy's tutorial uses it) | the repository does not take gene lists from memory; the list is one download away |
| E5 | Attack A3: is the DATP programme a dissociation-stress signature | the van den Brink 2017 supplementary list (M8) | the same rule; the trial is written to say so |
| E6 | Is primed AT2 a graded state rather than a cluster: the joint distribution of an identity-loss score (Etv5, Abca3, Cebpa) and the inflammatory score, inside the day-14 library, and in the IL-1beta organoid library | a new pre-registration with its own thresholds and its unreadable-if clause | writing it now would be a rule written after the result it responds to; decision 22 says the first outcome stands and a new question is a new pre-registration, which is the owner's to approve |
