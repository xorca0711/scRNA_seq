# Analysis trials motivated by Cardoso, Lee et al. 2026

The owner's plan for this paper is a set of **gates, not a fixed pipeline**:
each stage ends in a decision to continue, branch or stop, and an unexpected
result is to be reported and acted on rather than finished around. This page
is the pre-registered form of that plan and the log of what each gate
returned. Status words follow the root README claims table (Validated,
Descriptive only, Exploratory, Retracted-superseded, Not established). Every
trial writes a run record (`*_run_record.json`) holding the rules frozen
before any data is read, the inputs with sizes and modification times, the
package versions, and the outputs.

**The constraint that governs every trial here**, established by trial C0 and
repeated wherever it bites: every deposited mouse library of this paper pools
three mice, and each genotype contributes one library per sort. No genotype
contrast in the paper's own deposit has within-group biological replication.
Directions can be read and described; they cannot be tested, and no P value is
admissible on them. The single exception is the companion series GSE247505,
which carries two replicate libraries per arm and three time points.

| Trial | Gate | Question | Status | Where |
|---|---|---|---|---|
| C0 | 0 | What is actually deposited, and what statistical unit can it carry? | run; Descriptive only | [`trials/c0_data_reality_check/c0_summary.md`](trials/c0_data_reality_check/c0_summary.md) |
| C1 | 1 | Does the Pdgfrb+Runx1+Tnc+ fibrotic fibroblast subset reproduce, and does it separate cleanly? | run; **not recovered** by the frozen rule | [`trials/c1_fibroblast_compartment/c1_summary.md`](trials/c1_fibroblast_compartment/c1_summary.md) |
| C1b | 1 | Then what is there instead? (stop-and-characterise, plus a disclosed corrected pass) | run | [`trials/c1b_characterise_red2kras_private/c1b_summary.md`](trials/c1b_characterise_red2kras_private/c1b_summary.md) |
| C1c | 1 | Are the fibrotic and inflammatory programmes in the same cells at 2 weeks? | run | [`trials/c1c_fibrotic_inflammatory_overlap/c1c_summary.md`](trials/c1c_fibrotic_inflammatory_overlap/c1c_summary.md) |
| C1d | 1 | What does the mesenchymal sort actually contain? | run | [`trials/c1d_sort_purity/c1d_summary.md`](trials/c1d_sort_purity/c1d_summary.md) |
| C2 | 2b | Does removing Areg collapse the niche, and does anything survive it? | run; 4 of 5 directions met | [`trials/c2_areg_deletion_arm/c2_summary.md`](trials/c2_areg_deletion_arm/c2_summary.md) |
| C2b | 2b | The same composition, without the confidence floor that hid two directions | run | [`trials/c2b_composition_without_the_confidence_floor/c2b_summary.md`](trials/c2b_composition_without_the_confidence_floor/c2b_summary.md) |
| C3 | 2a | Is Areg a property of the DATP-like state under my own clustering, across the time course? | run; every frozen test holds | [`trials/c3_areg_state_specificity/c3_summary.md`](trials/c3_areg_state_specificity/c3_summary.md) |
| C-d | 2d | Does the transcriptome agree that fibroblast reprogramming precedes macrophage change? | **closed by C0, not run** | see below |
| D0 | 0 | What does the Choi 2020 deposit contain, the origin of the DATP state the C series leans on? (roadmap paper 2's Gate 0) | run; Descriptive only | [`../gate1_02_choi_2020/trials/d0_data_reality_check/d0_summary.md`](../gate1_02_choi_2020/trials/d0_data_reality_check/d0_summary.md) |

---

## Contents

Trials run in the order below, not in numerical order, because the E series was
opened when the deposit turned out to carry no testable contrast. The index with
one-line outcomes is [`trials/README.md`](trials/README.md).

| Section | Trial | In one line |
|---|---|---|
| [C0](#c0-data-reality-check-gate-0) | data reality check | the deposit cannot test anything |
| [C1](#c1-the-fibroblast-compartment-gate-1) | fibroblast compartment | not recovered, and the rule was at fault |
| [C1b, C1c, C1d](#c1b-c1c-c1d-stop-and-characterise-gate-1) | stop and characterise | found the sort contaminant |
| [C2, C2b](#c2-and-c2b-the-areg-deletion-arm-gate-2b) | the ligand deletion arm | the paper's epithelial result reproduces blind |
| [C3](#c3-is-areg-a-property-of-the-datp-like-state-gate-2a) | state specificity | holds in all four mutant libraries |
| [C5](#c5-figures-for-the-three-findings-and-three-corrections-2026-09-13) | figures | three findings drawn, three corrections forced |
| [E series](#e-series-public-data-extensions-of-the-hbegf-lead-2026-09-13) | E1 to E5 | the first admissible test, and E4 refuting claim C29's reading |
| [E6, C7, C8](#e6-c7-and-c8-three-questions-asked-properly-and-three-rules-that-showed-their-own-limits-2026-09-13) | three chosen follow-ups | three negatives and two disclosed rule defects |
| [C9](#c9-the-fst-and-runx2-population-pre-registered-2026-09-13) | the pre-registered lead | existence replicates, size does not |
| [C10](#c10-the-published-pathological-fibroblast-or-not-2026-09-13) | published state or not | the lead closes |
| [C11](#c11-figures-for-the-contradictions-2026-09-13) | contradiction figures | six drawn |
| [C12](#c12-the-whole-cellchatdb-ranked-2026-09-13) | the full resource scan | the ranking reports abundance; AREG still first among EGFR ligands |
| [D0](#d0-the-choi-2020-deposit-check-now-in-the-paper-2-folder) | the Choi 2020 deposit | one library per condition again; now with the rest of that paper's trials |

---

## C0. Data reality check (Gate 0)

Run on 2026-09-12. Script [`trials/c0_data_reality_check.py`](trials/c0_data_reality_check.py);
artefacts in [`trials/c0_data_reality_check/c0_summary.md`](trials/c0_data_reality_check/c0_summary.md).

### Pre-registration

Frozen before any matrix was read: the library list; the rule that a feature
whose identifier is not an Ensembl gene ID is not a gene and leaves the
expression matrix; species assignment by the higher-scoring of two fixed probe
panels; QC measured and not applied; gene spaces compared by the exact ordered
(ID, symbol) list; population presence measured as a detection fraction from
the paper's own marker sets; and **the replicate rule** -- the mouse is the
biological replicate, a library pooling several mice carries no within-group
replication, and independent animals behind a contrast are counted as
libraries, never as pooled mice and never as cells.

### Outcome

All five accessions are public and complete. Thirty libraries, 123,807
barcodes before quality control. Every species call agrees with the deposited
organism. Three distinct gene spaces are present.

| Accession | Libraries | Barcodes | Median genes per cell | What it is |
|---|--:|--:|---|---|
| GSE316241 | 2 | 13,226 | 1,824 to 2,111 | mesenchyme, Confetti vs Red2Kras |
| GSE316243 | 2 | 7,836 | 677 to 1,037 | niche: immune plus stroma |
| GSE316244 | 4 | 33,756 | 2,615 to 3,874 | Areg-flox arm, niche and RFP+ epithelium |
| GSE310335 | 2 | 9,408 | 3,757 to 5,576 | human KRASG12D alveolar organoids |
| GSE247505 | 20 | 59,581 | 1,136 to 3,212 | England 2025 lineage-labelled epithelium |

Four findings set the gates.

**1. No genotype contrast in this paper's deposit has biological replication.**
Each mouse library is a pool of three mice and each genotype contributes one
library per sort. This is not a criticism of the paper, whose tested claims
come from imaging and flow cytometry where n is mice; it is a statement about
what a reanalysis of the deposit can claim.

**2. The mesenchymal and immune compartments have no time course.** Every such
library is two weeks post-induction. The paper's ordering claim, that
fibroblast reprogramming precedes macrophage expansion, rests on
immunofluorescence at 1, 2, 4 and 8 weeks. **Gate 2 branch (d) is therefore
closed as a transcriptomic question and was not run.** Reporting it as
untestable is the result.

**3. The companion series is where the replication is.** GSE247505 carries two
libraries per arm at 4 days, 2 weeks and 12 weeks, and its RFP and YFP
libraries are the mutant and wild-type clones of the same animals: a
within-animal control. It is also the epithelial half of the paper's CellChat
object, and the Cardoso data-availability statement does not name it.

**4. One non-gene feature exists in the whole deposit.** `BSD`, the reporter
construct's blasticidin-resistance marker, in GSE316244 only. It is removed
from the matrix by the frozen rule and carried per cell, where it tracks the
sort: detected in 38.9% and 39.8% of cells in the two RFP-sorted libraries
against 5.0% and 2.9% in the two niche libraries. That makes it a usable check
on the sort, and nothing more: it is a transgene contig with low counts.

The paper's own quality thresholds and this repository's differ in a way that
matters for one library: GSE316243's Confetti library has a median of 677
genes per cell, which is why the paper applied a custom 500-gene floor to
exactly that library. Its Red2Kras partner sits at 1,037. The one comparison
the immune claims rest on is therefore between two libraries of visibly
different depth.

### Claims from C0

| Claim | Status |
|---|---|
| All five accessions are public, complete and readable; 30 libraries, 123,807 barcodes | Descriptive only |
| No genotype contrast in the Cardoso deposit carries within-group biological replication | Descriptive only (read from the deposited design) |
| The deposit cannot test the fibroblast-before-macrophage ordering | Not established, and not establishable from these data |
| BSD detection tracks the RFP sort (about 39% versus 3 to 5%) | Descriptive only |
| The two GSE316243 libraries differ about 1.5-fold in median genes per cell | Descriptive only |

---

## C1. The fibroblast compartment (Gate 1)

Run on 2026-09-12. Script [`trials/c1_fibroblast_compartment.py`](trials/c1_fibroblast_compartment.py);
artefacts in [`trials/c1_fibroblast_compartment/`](trials/c1_fibroblast_compartment/).

### Pre-registration

Frozen before any matrix was read. Repository MAD quality thresholds per
library, with the paper's own thresholds applied in parallel for comparison
only; Scrublet per library with the 10x prior and the automatic-threshold
sanity check; no batch correction in the primary embedding because library,
genotype and batch are one variable, with Harmony on library as a disclosed
sensitivity; log1p(CP10K), seurat_v3 2,000 HVGs, 30 principal components,
k = 30, Leiden 0.5 primary (fixed from trial S5) with 0.2 and 1.0 as
sensitivities, seeds 0. Population calls by score_genes on the paper's own
marker sets, assigned by modal argmax, confident only at 50%.

**The test, frozen.** The reprogrammed fibroblast state is *recovered* if some
fibroblast subcluster has at least 80% of its cells from the Red2Kras library,
the highest mean reprogrammed-fibroblast score of all fibroblast subclusters,
and Tnc detected in at least 40% of its cells. Also frozen: a clean-separation
rule (75% of the subcluster's cells with a kNN majority inside it), an
over-clustering guard (the HLCA 30%-in, 20%-out marker rule), a doublet guard,
and a one-state-versus-gradient rule (Gaussian mixture BIC on the score within
Red2Kras fibroblasts).

### Outcome: **not recovered**

11,690 cells analysed (5,751 Confetti, 5,939 Red2Kras). The subcluster with
the highest reprogrammed-fibroblast score failed two of the three criteria:
it is 55.9% Red2Kras against the 80% floor, and Tnc is detected in 7.7% of its
cells against the 40% floor.

Three things are worth more than the verdict itself.

**The quality rules are not interchangeable, and the difference is
asymmetric.** This repository's per-library MAD rule keeps 6,057 of 6,435
Confetti cells and 6,260 of 6,791 Red2Kras cells. The paper's stated
thresholds would keep 5,390 and 4,797, removing 16% of the control library
and 29% of the tumour library. A composition comparison between two libraries
is sensitive to a 13-point asymmetry in what quality control discards.

**Scrublet's automatic threshold is unusable on these libraries.** It called
0.07% and 0.00% of cells doublets against a 10x prior of about 5%, so the
pre-registered sanity check rejected it and the fallback quantile applied.
The doublet calls here are a ranking cut, not a detection, exactly as the
repository's known issue 5 records for the mouse series.

**The frozen selection rule defeated itself, and this is disclosed rather than
quietly repaired.** The rule selected fibroblasts as "clusters whose
*confident* call is alveolar, adventitial or reprogrammed fibroblast". Every
cluster the caller labelled *reprogrammed fibroblast* (clusters 8, 12 and 15)
fell below the 50% confidence floor, as did cluster 14, which is 99.6%
Red2Kras. The selection therefore removed the candidates the test was about,
and the test ran on 6,281 cells that could not contain the answer. C1's
outcome stands as recorded; trial C1b adds the corrected pass beside it. This
follows DEVELOPMENT decision 13: a rule revision is disclosed, and the
first-run outcome stays in the record.

Two further pre-registered readings did return: the Red2Kras fibroblast score
is **two states rather than a gradient** by the BIC rule (one component 8.88,
two components -495.19, on 2,211 cells), and the Harmony sensitivity gives an
adjusted Rand index of 0.685 against the primary partition, 8 subclusters
against 7. The Harmony run is reported only as a sensitivity, because
correcting on library here would correct on genotype.

### Claims from C1

| Claim | Status |
|---|---|
| The published reprogrammed fibroblast subset is not recovered under this repository's clustering by the frozen three-part rule | Descriptive only (negative result) |
| The paper's quality thresholds would remove 29% of the Red2Kras library against 16% of the Confetti library | Descriptive only |
| Scrublet's automatic threshold fails on both libraries and the calls are a ranking cut | Descriptive only |
| The Red2Kras fibroblast score is two states rather than one gradient | Exploratory (one library, no replication) |
| The frozen selection rule excluded the candidate clusters | Disclosed rule failure; corrected pass in C1b |

---

## C1b, C1c, C1d. Stop and characterise (Gate 1)

Gate 1 says that when the published subset does not separate, the trial stops
and characterises rather than pressing on. Three short trials do that. All
three state in their own rules that they were fixed **after** the preceding
trial's cluster tables had been seen, the same disclosure trial S4 made.

### C1b: what is there instead

Script [`trials/c1b_characterise_red2kras_private.py`](trials/c1b_characterise_red2kras_private.py);
artefacts in [`trials/c1b_characterise_red2kras_private/`](trials/c1b_characterise_red2kras_private/).

**Question A, does any cluster match the published signature** (at least 80%
Red2Kras, with Tnc, Acta2 and either Pdgfrb or Runx1 each detected in at least
40% of cells): **none qualifies**. The nearest is cluster 14, and it misses on
one number:

| | cluster 14 |
|---|---|
| cells | 937 |
| fraction Red2Kras | 0.996 |
| Tnc detected | **0.399** (floor 0.400) |
| Acta2 detected | 0.574 |
| Runx1 detected | 0.540 |
| Pdgfrb detected | 0.326 |

The floor was not moved after the fact. A rule that a result misses by one
thousandth is a rule that was arbitrary at the margin, and saying so is more
useful than either quietly relaxing it or reporting "not recovered" as though
the population were absent.

**Question B, the Red2Kras-private clusters** (at least 90% Red2Kras):
clusters 4, 10, 11, 14 and 16. **None** is explained by the quality criteria
that explained mouse cluster 23 of the GSE262927 series: none is low-count and
ambient-like, and none is doublet-enriched. They are real structure.

**Question C, the corrected selection, post hoc.** Replacing C1's
confident-call selection with a compartment gate (Col1a1 detected; Ptprc,
Pecam1 and Epcam not detected) keeps 8,431 of 11,690 cells and **still returns
not recovered**, for a third reason the first pass hid: the subcluster with the
highest reprogrammed-fibroblast score has Acta2 detected in 99.7% of its cells
and Pdgfra in 2.2%. It is smooth muscle. Two of the six genes in the paper's
own reprogrammed-fibroblast set, Acta2 and Pdgfrb, are mural markers, so a
score over that set is maximised by mural cells unless they are excluded
first. That is a property of the marker set, not of this data, and it is worth
carrying into any future use of it.

### C1c: are the fibrotic and inflammatory programmes in the same cells

Script [`trials/c1c_fibrotic_inflammatory_overlap.py`](trials/c1c_fibrotic_inflammatory_overlap.py).
The paper makes fibrotic and inflammatory fibroblasts distinct populations,
with the inflammatory cells lacking Tnc. A cluster-level detection fraction
cannot tell interleaved subsets from co-expression, so this asks per cell,
against what independence predicts.

| Group | cells | Tnc | Lcn2 or Saa3 | both | expected if independent | ratio |
|---|--:|--:|--:|--:|--:|--:|
| Red2Kras, private clusters | 1,854 | 27.8% | 29.5% | 9.7% | 8.2% | 1.18 |
| Red2Kras, other clusters | 4,085 | 17.0% | 4.5% | 0.7% | 0.8% | 0.96 |
| Confetti, all | 5,751 | 11.5% | 1.2% | 0.1% | 0.1% | 1.06 |

All three sit below the frozen 1.25 threshold: **consistent with separate
cells**, which supports the paper's claim by a route the paper did not use.
Per cluster the separation is sharper than the pooled number suggests:
clusters 4 and 11 are inflammatory-dominant with ratios of 0.74 and 0.43, so
in those clusters the two programmes actively avoid each other, while clusters
14 and 16 sit slightly above independence at 1.30 and 1.66.

### C1d: what the mesenchymal sort actually contains

Script [`trials/c1d_sort_purity.py`](trials/c1d_sort_purity.py). C1b's five
private clusters were "not explained by quality", and their top genes
suggested the explanation is composition rather than quality.

**757 of 11,690 cells, 6.5%, sit in clusters that read as off-target for a
CD45-CD31-EpCAM- sort.** Clusters 8, 10 and 15 read as immune; clusters 11 and
12 read as epithelial. The one that matters:

| cluster | cells | fraction Red2Kras | epithelial markers | Col1a1 | **Areg** |
|---|--:|--:|--:|--:|--:|
| 11 | 184 | 0.989 | 98.4% | 21.2% | **88.0%** |
| 12 | 39 | 0.436 | 97.4% | 23.1% | 69.2% |

An Areg-high mutant epithelial population is sitting inside a library sorted
as mesenchyme, almost entirely in the Red2Kras arm. **The paper's own analysis
is immune to this**, because it took its epithelial cells from the separate
lineage-labelled series rather than from this library. A reanalysis that
computed epithelium-to-fibroblast signalling inside this library alone would
not be: it would be reading a sort contaminant as the signalling source. That
is the single most practically useful thing Gate 1 produced.

The same trial clears cluster 14: 95.3% mesenchymal, 89.9% Col1a1, 3.4%
epithelial. It is a genuine fibroblast population, not a contaminant and not a
doublet.

### What Gate 1 concluded

The published population **is present**, as clusters 14 and 16 together, 1,080
cells at about 99% Red2Kras, carrying the fibrotic programme on a retained
alveolar-fibroblast identity (cluster 14: Pdgfra 69.5%, Col13a1 44.4%, Tcf21
58.0%, with Tnc 39.9%, Acta2 57.4%, Runx1 54.0%, Pdgfrb 32.6%, Sfrp1 55.3%;
cluster 16 the same with Mki67 at 39.2%). This is consistent with the paper's
own description of fibroblasts that gain Pdgfrb and Acta2 with *reduced*, not
absent, Pdgfra.

The pre-registered rule nonetheless returned "not recovered", for four
separable reasons, all disclosed and none repaired in place:

1. the selection step took only clusters with a confident fibroblast call, and
   the candidates were not confident;
2. Acta2 and Pdgfrb in the paper's marker set are mural markers, so the score
   is maximised by smooth muscle;
3. Runx1 is also a myeloid transcription factor, so the score additionally
   flags the myeloid contaminant the sort carries;
4. the Tnc floor was missed by one thousandth.

**The owner's decision.** Whether to record this as "recovered, with a
disclosed threshold miss and three disclosed rule defects" or to leave it as
"not recovered" is a retain/reject call, not a computation. The artefacts
support either wording and the numbers do not change.

### Claims from C1b, C1c and C1d

| Claim | Status |
|---|---|
| The published reprogrammed fibroblast population is present as clusters 14 and 16 (1,080 cells, about 99% Red2Kras) with the fibrotic programme on a retained alveolar identity | Descriptive only; owner decision pending on whether this supersedes C1's "not recovered" |
| Cluster 14 misses the pre-registered Tnc floor by 0.001 | Descriptive only (disclosed near-miss; the floor was not moved) |
| The paper's reprogrammed-fibroblast marker set scores highest on smooth muscle, because Acta2 and Pdgfrb are mural markers | Descriptive only, and a caution for any future use of that set |
| Fibrotic and inflammatory markers mark separate cells at 2 weeks (ratio 1.18 against a 1.25 threshold) | Descriptive only; supports the paper's distinctness claim by a different route |
| The GSE316241 mesenchymal sort carries 6.5% off-target cells, including 184 Areg-high mutant epithelial cells almost entirely in the Red2Kras arm | Descriptive only; a practical caution for ligand-receptor reanalysis within that library |
| The five Red2Kras-private clusters are not explained by low counts or doublets | Descriptive only (negative result) |

---

## C2 and C2b. The Areg deletion arm (Gate 2b)

Run on 2026-09-12. Scripts
[`trials/c2_areg_deletion_arm.py`](trials/c2_areg_deletion_arm.py) and
[`trials/c2b_composition_without_the_confidence_floor.py`](trials/c2b_composition_without_the_confidence_floor.py);
artefacts in [`trials/c2_areg_deletion_arm/`](trials/c2_areg_deletion_arm/)
and [`trials/c2b_composition_without_the_confidence_floor/`](trials/c2b_composition_without_the_confidence_floor/).

The owner's branch (b): does removing the ligand reproduce the collapse at
transcriptome level, **and does anything fail to collapse**. The second half
is the part the paper does not ask.

### Pre-registration

Processing identical to C1. The two sorts are embedded separately, because a
joint embedding would confound sort with genotype. Five directions were taken
from the paper *before* these data were opened, and each is recorded as met or
not met: in Areg-flox/flox relative to Areg-flox/+, the reprogrammed
fibroblast share falls, the DATP-like share falls, the Cd177-positive share
falls, the AT2 share rises, and in alveolar macrophages the inflammatory genes
fall while MHC-II genes rise. Part B, frozen separately: a gene of the
fibrotic set is **Areg-independent** if its detection in flox/flox fibroblasts
stays at or above 80% of its flox/+ value while also exceeding its Confetti
value by at least 5 points.

### Outcome: the paper's epithelial result reproduces blind

| Population | Areg-flox/+ | Areg-flox/flox | The paper (Fig. 4m) |
|---|--:|--:|---|
| DATP-like | 46.6% | 22.0% | 50.1% to 25.9% |
| AT2 | 29.2% | 59.1% | 21.3% to 54.6% |

This is an independent pipeline, with the authors' labels never used, landing
within a few points of the published composition on both populations. The
single largest cluster effect is sharper than the share: the confidently
called DATP-like cluster 1 holds 1,712 flox/+ cells against 111 flox/flox, a
15.4-fold depletion.

**Directions met: four of five, with the fifth unscorable.** Three were met in
C2 as frozen (DATP-like falls, AT2 rises, and in 343 alveolar macrophages
Cxcl2 and Ccl9 fall while H2-Ab1 and H2-Eb1 rise). Two returned 0.0% against
0.0%, and for the same reason C1's rule defeated itself: the composition
metric counted only clusters whose modal call holds at least 50% of their
cells, and neither the reprogrammed fibroblast nor the Cd177-positive call
ever clears that floor. That is not a measurement of absence.

C2b recomputes the same shares with the floor removed, carrying the weakest
mode fraction beside every number. The reprogrammed fibroblast share then
falls from **27.2% to 20.1%**, so the direction is met; the Cd177-positive
call is still unscorable, because no cluster carries it as a modal call at
all, which is a statement about resolution rather than about the population.
The cluster-level depletions are the stronger statement, and they run right
through the paper's cascade:

| Compartment | Cluster | Call | flox/+ | flox/flox | ratio |
|---|---|---|--:|--:|--:|
| RFP+ epithelium | 1 | DATP-like | 1,712 | 111 | 15.4 |
| niche | 14 | mesothelium | 353 | 34 | 10.4 |
| niche | 19 | alveolar macrophage | 310 | 33 | 9.4 |
| niche | 17 | reprogrammed fibroblast | 210 | 45 | 4.7 |

### What fails to collapse

Of the six genes in the paper's own reprogrammed-fibroblast marker set, four
fall when Areg is deleted and two do not:

| Gene | flox/+ | flox/flox | verdict |
|---|--:|--:|---|
| Fst | 0.336 | 0.143 | collapses |
| Runx2 | 0.175 | 0.069 | collapses |
| Tnc | 0.183 | 0.138 | collapses |
| Acta2 | 0.279 | 0.215 | collapses |
| **Pdgfrb** | 0.452 | 0.389 | **survives** |
| **Runx1** | 0.560 | 0.543 | **survives** |

The matrix and contractile half of the programme is Areg-dependent; the
Pdgfrb and Runx1 induction is not, or is much less so. The paper treats the
reprogrammed fibroblast as one state that the Areg-EGFR axis initiates, and
this says the state has at least two separable parts. Runx1 is the
transcription factor the paper itself highlights, which makes its survival the
more interesting half.

**Three limits on that result, all of which bite.** It is one library per
genotype, so this is a difference between two libraries. The Confetti baseline
comes from the other series, which is shallower (median 1,824 to 2,111 genes
per cell against 2,615 to 2,750 here), and a detection fraction rises with
depth, so every cross-series comparison in the survival table is
depth-confounded and the Confetti column should not be read as a quantitative
baseline. Only the within-series flox/+ against flox/flox comparison, where
the two libraries are of similar depth, carries the claim. Two of the four
"Areg-independent" calls in the artefact, Pdgfra and Col13a1, are
alveolar-identity genes that were in the comparison set and should not be read
as part of the fibrotic programme at all.

### Claims from C2 and C2b

| Claim | Status |
|---|---|
| The paper's epithelial composition reproduces from a blind pipeline: DATP-like 46.6% to 22.0%, AT2 29.2% to 59.1%, against the paper's 50.1 to 25.9 and 21.3 to 54.6 | Descriptive only (one library per genotype), and the closest thing to a validation this deposit allows |
| Four of the five pre-registered directions are met; the fifth is unscorable because the Cd177-positive state is not resolved at this resolution | Descriptive only |
| The confidently called DATP-like cluster is depleted 15.4-fold in the Areg-deleted library | Descriptive only |
| Alveolar macrophages are depleted 9.4-fold and shift away from the inflammatory and toward the MHC-II profile | Descriptive only; matches Extended Data Fig. 10e in direction |
| Tnc, Acta2, Fst and Runx2 fall on Areg deletion while Pdgfrb and Runx1 do not | Exploratory. The best new lead here, and it needs a depth-matched control before it is more than that |
| Mesothelial cells are depleted 10.4-fold on Areg deletion | Exploratory; the paper reports mesothelial-like cells as Red2Kras-enriched but does not test their Areg dependence |

---

## C3. Is Areg a property of the DATP-like state (Gate 2a)

Run on 2026-09-12. Script
[`trials/c3_areg_state_specificity.py`](trials/c3_areg_state_specificity.py);
artefacts in [`trials/c3_areg_state_specificity/`](trials/c3_areg_state_specificity/).
33,217 cells from the ten Experiment 1 libraries of GSE247505.

### What this trial is, and what it is not

The owner's branch (a) asked whether Areg stays top when the DATP-like cluster
is defined by our own clustering, what ranks below it, and whether the
ranking holds across the time course. **CellChat itself cannot be run here**:
it is an R package and this machine has no R, and a Python reimplementation
would use a different ligand-receptor resource, so calling its output a
CellChat rerun would be false. What the paper's communication claim rests on
is an expression fact, that Areg is induced specifically in the DATP-like
mutant state and is the top EGFR ligand there, and that fact is re-derived
directly here with three things the paper's own analysis did not use: the
state defined by this repository's clustering, the within-animal wild-type
control that the Red2Onco design provides, and the time course. **No
ligand-receptor probability is computed and nothing here shows that a
fibroblast receives the signal.**

### The batch rule, exercised properly for the first time in this deposit

This is the only part of the Cardoso material with replicate libraries, so it
is the only place where this repository's rule, that batch correction is
decided from measured replicate mixing rather than assumed, can actually run.
Same-library enrichment among the 30 nearest neighbours, within each arm:
1.53, 1.62, 1.74 and 1.80 against the pipeline's failure threshold of 2.0.
**No correction applied**, by the rule and not by preference.

### Outcome: every frozen test holds

**T1, state specificity.** Mean log1p(CP10K) Areg is higher in DATP-like cells
than in AT2 cells in **all four** mutant libraries, at both time points, with
the state defined by our own clustering rather than by the authors' labels:

| Library | Areg in DATP-like | Areg in AT2 |
|---|--:|--:|
| 4 days, replicate 1 | 2.55 | 0.54 |
| 4 days, replicate 2 | 2.12 | 0.49 |
| 2 weeks, replicate 1 | 2.92 | 1.13 |
| 2 weeks, replicate 2 | 2.83 | 0.69 |

**T2, the ranking, including what the paper does not show.** In all four
libraries the order is identical: **Areg > Hbegf > Ereg > Tgfa**. Areg is top
everywhere, so the paper's headline survives a change in how the state is
defined. The second place is the new part: **Hbegf**, consistently, and ahead
of Ereg, which is the ligand the paper went on to test in culture alongside
Areg. Hbegf is also an EGFR ligand, so this is a concrete and checkable
prediction rather than a curiosity.

**T3, the within-animal control.** The DATP-like share of cells is far higher
in the mutant RFP clones than in the wild-type YFP clones of the same animals:
14.7% against 0.08% at 4 days, and 25.8% against 1.1% at 2 weeks. The state
is essentially absent from wild-type clones in the same lung, which is a
cleaner control than any between-genotype comparison in this deposit.

**T4, time.** The DATP-like Areg mean rises from 2.34 at 4 days to 2.88 at 2
weeks. Two libraries per arm: a ranking, not a trend, and no test is computed.

### Claims from C3

| Claim | Status |
|---|---|
| Areg is higher in the DATP-like state than in AT2 cells in all four mutant libraries, with the state defined by this repository's clustering | Descriptive only, and the best-replicated result in this stage (two libraries per arm, two time points, within-animal control) |
| Areg is the top EGFR ligand in the DATP-like state in every mutant library, and the ranking Areg > Hbegf > Ereg > Tgfa is identical in all four | Descriptive only. The paper's top hit survives a change of state definition. |
| Hbegf ranks second, ahead of Ereg | Exploratory, and the most testable thing this stage produced: a second EGFR ligand the paper did not follow. |
| The DATP-like state is nearly absent from wild-type clones in the same animals (0.08% and 1.1%) | Descriptive only |
| Replicate libraries mix within every arm (1.53 to 1.80 against a threshold of 2.0), so no batch correction is applied | Descriptive only (a decision record, the counterpart of C5) |
| Anything about whether fibroblasts receive this signal | Not established, and not establishable without a communication analysis this machine cannot run |

---

## C5. Figures for the three findings, and three corrections (2026-09-13)

Script [`trials/c5_figures_for_the_three_findings.py`](trials/c5_figures_for_the_three_findings.py);
artefacts in [`trials/c5_figures_for_the_three_findings/`](trials/c5_figures_for_the_three_findings/).
The rationale, biology and next steps for all three findings are in
[`DIVERGENCES_AND_NEXT.md`](DIVERGENCES_AND_NEXT.md).

None of the earlier trial figures showed the three findings. Drawing them
forced three corrections to wording used above and in the session summary.
The rows above are left as written; these supersede them.

1. **C3, T2 and claim C33.** Hbegf is second by abundance in every mutant
   library, but third by enrichment over AT2 cells of the same library in three
   of four, where Ereg is second. The sentence "Hbegf ranks ahead of the ligand
   the paper followed into culture" was true of abundance only.
2. **C3, T3 and T4.** The 4-day replicates disagree about 24-fold on the
   DATP-like share (1.2%, 44 cells; 28.2%, 541 cells). The 14.7% median
   describes neither library.
3. **C2, Part B and claim C29.** The fibrotic response is graded in three tiers,
   not split in two: Runx1 97% and Pdgfrb 86% retained, Tnc 76% and Acta2 77%,
   Fst 42% and Runx2 39%. The frozen 80% line made it binary.

A discriminating check on cluster 11 also narrows the contaminant: Epcam mRNA
is present in 78% of its cells, so the RNA does not support mutant cells
escaping the sort by losing Epcam; about one in five cells co-detects Krt8 and
Col1a1, a doublet minority.

---

## E series. Public-data extensions of the Hbegf lead (2026-09-13)

Plan and literature grounding in [`DIVERGENCES_AND_NEXT.md`](DIVERGENCES_AND_NEXT.md)
section 7. The reason for leaving this deposit is stated there and is worth
repeating here: the Cardoso deposit pools mice and gives one library per
genotype, so no trial in the C series could test anything. These datasets have
donors or animals as the unit, which is what makes a test admissible.

Order of operations was set by the owner: run E1 first, and continue only if
E1's pre-registered refutation rule did not fire. It did not.

### E1: which human cells make which EGFR ligand in lung adenocarcinoma

Script [`trials/e1_human_luad_ligand_sources.py`](trials/e1_human_luad_ligand_sources.py);
artefacts in [`trials/e1_human_luad_ligand_sources/`](trials/e1_human_luad_ligand_sources/).
Dataset GSE131907 (Kim et al. 2020, doi:10.1038/s41467-020-16164-1), the human
comparison the paper itself used. 88,144 cells in scope across 11 donors with
paired tumour and normal lung.

**T1, the first admissible test in this folder.** In tumour lung, AREG
detection is higher in epithelial than in myeloid cells, paired within donor.
Median 0.3363 against 0.2147 across 11 donors, Wilcoxon signed-rank
p = 0.001953. The frozen refutation rule did not fire, so E2 to E4 proceeded.

**T2, corrected after the fact.** The run record names the "neutrophil" gate
as the top HBEGF compartment. That label is wrong, and the correction matters
because it changes which pre-named outcome was hit. This deposit annotates no
neutrophils at all, dissociation and 10x capture lose them, and the gate's
crosstab against the deposited cell types shows it is 83.3 per cent deposited
myeloid cells and 12.5 per cent T or NK cells. Read properly, T2's outcome is
the pre-named **myeloid-dominant** one. That supports the published human
emphasis (Hult et al. 2022, doi:10.1165/rcmb.2022-0174OC; Van Hiep et al. 2022,
doi:10.3389/fonc.2022.963896) and does not match what trial C6 found in mouse,
where the Areg-independent share of Hbegf sat in endothelium and mesenchyme and
the myeloid compartments were the ones that collapsed.

**T3.** EGFR detection is highest in epithelium, 0.286. In the mouse deposit
Egfr is mesenchymal, 73 per cent, against 4.7 per cent in RFP-positive
epithelium. The receptor is not in the same place in the two datasets.

**A caveat found after the fact, which motivated E1b.** AREG in epithelium is
not tumour-enriched at compartment level, 0.336 in tumour lung against 0.356 in
matched normal lung, and HBEGF is actually higher in normal lung. Broad
epithelial gating pools malignant states with normal AT2, club and ciliated
cells, so it dilutes exactly the transitional state the paper's claim is about.
Compartment gating answers a coarser question than the one asked.

### E1b: the same question at the resolution the claim lives at

Script [`trials/e1b_human_luad_by_subtype.py`](trials/e1b_human_luad_by_subtype.py);
artefacts in [`trials/e1b_human_luad_by_subtype/`](trials/e1b_human_luad_by_subtype/).
Same dataset, regrouped by the deposited `Cell_subtype` annotation, reading the
gene vectors E1 had already cached. Rules fixed after E1's compartment table was
seen and before any subtype value was computed; post hoc and disclosed as such.

**T5 is not computable with these labels, and was not repaired.** The deposited
annotation is tissue-exclusive. AT2 is assigned only in normal lung, 11 donors
and 2,020 cells, and the tumour states tS1 to tS3 only in tumour lung. No donor
has both, so the paired within-tumour comparison has zero pairs. Lowering the
50-cell floor cannot help, because the count of AT2-labelled cells in tumour
samples is zero, not small. The deposit named malignant epithelium by tissue of
origin, which forecloses the within-tumour contrast.

**T6 confirms the myeloid reading with the deposited labels**, independently of
the marker gates that mislabelled it in E1. Median per-donor HBEGF detection in
tumour lung: alveolar macrophages 0.712, CD163+CD14+ dendritic cells 0.655,
CD1c+ dendritic cells 0.552, monocyte-derived macrophages 0.423, then the
tumour state tS2 at 0.420 and tS1 at 0.244. COL13A1-positive matrix
fibroblasts, the alveolar fibroblast counterpart, sit at 0.087.

**T7 falls below its own floor.** Three donors clear 50 cells in both
populations against a floor of five, so no test was run. The direction is
opposite to the mouse: EGFR 0.406 in the tumour states against 0.226 in
COL13A1-positive matrix fibroblasts, with the fibroblast value lower in all
three donors.

**A qualification of E1's T1, recorded because it cuts against my own result.**
At subtype resolution CD1c+ dendritic cells detect AREG in 0.618 of cells,
above tS1 at 0.555 and tS2 at 0.532. T1 stands exactly as written, epithelium
above the myeloid average, because the myeloid compartment pools those
dendritic cells with monocyte-derived macrophages at 0.159 and alveolar
macrophages at 0.298. What T1 does not license is the sentence "epithelium is
the AREG source in human lung adenocarcinoma". Specific dendritic-cell subsets
match or exceed the tumour states.

### E4: is the Areg-independent tier tumour-specific, or does injury alone make it

Script [`trials/e4_bleomycin_fibrotic_genes.py`](trials/e4_bleomycin_fibrotic_genes.py);
artefacts in [`trials/e4_bleomycin_fibrotic_genes/`](trials/e4_bleomycin_fibrotic_genes/).
Dataset GSE132771 (Tsukui et al. 2020, doi:10.1038/s41467-020-15647-5), the
injury comparison the paper itself used, restricted to the four mouse
Col1a1-GFP-positive libraries. Bleomycin 3,631 and 3,387 cells; untreated 4,063
and 3,271 cells after the C1 quality control and doublet removal.

**The frozen rule fired against my own claim.** Runx1 and Pdgfrb are both
injury-generic, higher in both bleomycin libraries than in both untreated
libraries: Runx1 0.225 and 0.237 against 0.078 and 0.088, Pdgfrb 0.225 and
0.208 against 0.117 and 0.126. The reading fixed in advance for that outcome
therefore applies. Their survival of Areg deletion in trial C2 is what an
activated lung fibroblast does after any injury, not evidence of a second
tumour-specific signal. **Claim C29 weakens accordingly**, and the
second-signal hypothesis loses the observation that motivated it.

Tnc, Fst and Runx2 are injury-generic too, so the graded tiers of C2 do not map
onto an injury-generic and a tumour-specific half. Acta2 is not injury-generic,
because untreated library UT2 reaches 0.423 above both bleomycin values, which
is what a smooth-muscle contribution to a Col1a1-GFP sort looks like. Hbegf and
Egfr are also injury-generic, 0.118 and 0.135 against 0.072 and 0.062 for
Hbegf, so neither is tumour-specific either. Pdgfra and Col13a1 fall with
injury in both replicates, consistent with the loss of alveolar fibroblast
identity that Tsukui et al. describe.

Two animals per group. The frozen rule deliberately required both replicates to
agree rather than averaging them, so these are directions with no P value.

### E2 and E3: the same question in human pulmonary fibrosis

One script, [`trials/e2_human_fibrosis_ligand_sources.py`](trials/e2_human_fibrosis_ligand_sources.py),
two cohorts, artefacts in [`trials/e2_human_fibrosis_ligand_sources/`](trials/e2_human_fibrosis_ligand_sources/).
E2 is GSE136831 (Adams et al. 2020, doi:10.1126/sciadv.aba1983), 60 donors in
scope, with the aberrant basaloid population as the transitional state. E3 is
GSE135893 (Habermann et al. 2020, doi:10.1126/sciadv.aba1972), 22 donors, with
the KRT5-negative KRT17-positive population as the transitional state and
transitional AT2 as a second candidate. Both use the deposited cell-type
grouping; neither is reclustered.

**T1, the one comparison that cleared its floor, does not support the
prediction it was built on.** In the E2 cohort, HBEGF detection in pooled
myeloid cells is 0.366 against 0.390 in aberrant basaloid cells across seven
donors, Wilcoxon signed-rank p = 0.297. The frozen expectation from Hult et al.
2022 was myeloid dominance over the transitional state, and that is not what
seven donors show: the two populations are indistinguishable here. The
direction is not even internally consistent, since myeloid cells are higher in
five of seven donors while the median is lower, which is what a null result
with an outlier looks like. Seven donors is also the fewest this comparison
could have had, because the aberrant basaloid state is rare, so this is a
failure to detect a difference rather than a demonstration of equality.

In the E3 cohort the same comparison has three donors against a floor of five,
so no test was run. Its direction favours myeloid cells, 0.709 against 0.418.

**T3, the Zhao prediction, is not testable in either cohort and its directions
disagree.** Zhao et al. 2024 (doi:10.1016/j.stem.2024.07.004) predicts AREG in
the transitional state above AT2. Aberrant basaloid cells sit at 0.485 against
0.767 in ATII cells, which is the opposite direction, on three donors.
KRT5-negative KRT17-positive cells sit at 0.700 against 0.672 in AT2, and
transitional AT2 at 0.781 against 0.574, both in the predicted direction, on
three and two donors. No comparison clears the five-donor floor, so all three
are reported and none is read. The prediction is neither supported nor refuted
by these deposits.

**The rankings agree with each other and with E1b.** HBEGF is topped by
myeloid populations in both cohorts: cDC2 at 0.444, alveolar macrophages at
0.437 and classical monocytes at 0.404 in E2; macrophages at 0.508 and
monocytes at 0.497 in E3. The transitional state is close behind rather than
absent, aberrant basaloid fourth at 0.390. AREG is topped by dendritic cells
and monocytes in both: plasmacytoid dendritic cells at 0.925, cDC2 at 0.874 and
classical monocytes at 0.782 in E2; cDCs at 0.650 and monocytes at 0.528 in E3,
against ATII at 0.755 and AT2 at 0.711.

That last pattern is the most consistent thing the extensions produced, because
E1b found it independently in lung adenocarcinoma, where CD1c-positive
dendritic cells topped AREG above both tumour epithelial states. Three human
datasets, two diseases and two annotation vocabularies agree that dendritic
cells and monocytes are a major source of both ligands in human lung.

It is not a new fact about biology and should not be written up as one.
Leukocyte-derived AREG is established: Zaiss et al. 2015
(doi:10.1016/j.immuni.2015.01.020) reviews mast cells, basophils, ILC2s and
tissue-resident regulatory T cells as AREG sources alongside epithelium and
mesenchyme. The nearest thing to news here is which leukocytes, since the
dendritic-cell and monocyte emphasis is not the one that review foregrounds,
and in E1b the regulatory T cells that review nominates sit near the bottom at
0.07. What it does do is constrain a reading rather than add a finding: an
epithelium-to-fibroblast account of the AREG axis in human lung leaves out the
compartment with the highest detection, and a mouse sort of RFP-positive
epithelium cannot see it.

**EGFR sits in the mesenchyme here, unlike in adenocarcinoma.** Fibroblasts
reach 0.650 and basal cells 0.774 in E2, with myofibroblasts at 0.449 and
PLIN2-positive fibroblasts at 0.512 in E3. That matches the mouse, where Egfr
is mesenchymal at 73 per cent, and differs from E1 and E1b, where the tumour
epithelial states topped the receptor. Receptor placement therefore looks
disease-dependent rather than species-dependent, which is the opposite of what
the ligand comparison suggested.

### E2b: the T4 reading that E2 declared and did not compute

Script [`trials/e2b_t4_disease_against_control.py`](trials/e2b_t4_disease_against_control.py);
artefacts in [`trials/e2b_t4_disease_against_control/`](trials/e2b_t4_disease_against_control/).
E2 froze four readings and wrote three. T4, disease against control for the
same cell types, was in the rules and never reached a table. E2b supplies it by
reading only E2's tracked per-donor table, which is the same remedy trial C2b
applied to C2, and no rule is altered after the fact.

Thirty-seven cohort and cell-type combinations clear the three-donor floor on
both sides. AREG detection is higher in IPF in 24 of them, median difference
0.027; HBEGF in 24 of 37, median difference 0.019; EREG in 12 of 37 and EGFR in
9 of 37. Twenty-four of thirty-seven is a weak majority, the differences are a
couple of percentage points of detection, the donors are not paired and the
libraries are not depth-matched. The honest reading is that neither ligand
shows clear disease enrichment in these cohorts, which is the same answer E1
gave for tumour against matched normal lung. Both ligands look like properties
of the cell types that carry them rather than of the injury.

### Claims from E2, E3 and E2b

| Claim | Class |
|---|---|
| HBEGF detection in myeloid cells exceeds that of the aberrant basaloid transitional state in IPF | Not established; 7 donors, p = 0.297, and the direction is internally inconsistent |
| Dendritic cells and monocytes are a major source of AREG and HBEGF in human lung, at or above the epithelial states | Descriptive only, and consistent across three datasets and two diseases |
| That leukocyte AREG source is a new finding | No; established in the immunology literature, and recorded here only because it constrains an epithelium-centric reading |
| AREG is higher in the transitional state than in AT2 cells in human fibrosis (the Zhao prediction) | Not established; no cohort clears the donor floor and the two cohorts disagree in direction |
| EGFR sits in the mesenchyme in human fibrosis, as in the mouse, and in the epithelium in adenocarcinoma | Descriptive only; receptor placement tracks disease rather than species |
| AREG or HBEGF is enriched in fibrotic lung relative to control | Not established; 24 of 37 cell types with median differences near 0.02, unpaired and not depth-matched |

### E5: the one figure the E series earns

Script [`trials/e5_figure_for_the_refutation.py`](trials/e5_figure_for_the_refutation.py);
figure [`trials/e5_figure_for_the_refutation/e5_retention_against_injury.png`](trials/e5_figure_for_the_refutation/e5_retention_against_injury.png).
Reads only the tracked tables of C5 and E4 and refuses to draw if the retained
tier is no longer Runx1 and Pdgfrb.

The refutation is a two-panel argument that a table hides. On the left, how
much detection each of the six genes keeps after Areg deletion. On the right,
what bleomycin alone does to the same genes in mesenchyme with no oncogene.
Drawing it sharpened the claim in a way worth recording: five of the six genes
are injury-generic, not only the retained pair, so the tiers carry no
tumour-specific information at all rather than merely having an innocent
explanation for their top two. The figure's first headline said the surviving
genes are the ones injury turns on, which implied a correspondence the panels
do not show, and was replaced.

Acta2 is the instructive exception in both panels. It sits just under the
frozen 80 per cent line on the left and is the one gene that fails the
injury-generic rule on the right, because the two untreated animals disagree
across a range of 0.297 to 0.423. That spread is a smooth-muscle contribution
to a Col1a1-GFP sort, the same contamination trial C1b found in the Cardoso
mesenchymal sort, and it is drawn rather than averaged away.

### Claims from E1, E1b and E4

| Claim | Class |
|---|---|
| In human lung adenocarcinoma, AREG detection is higher in epithelial than in myeloid cells within the same donor (11 donors, paired p = 0.0020) | Validated, the only tested claim in this folder |
| HBEGF in human lung adenocarcinoma is myeloid-dominant, by marker gate and by deposited subtype label alike | Descriptive only, and it contradicts the mouse pattern in C6 |
| Epithelium is the AREG source in human lung adenocarcinoma | Not established; CD1c+ dendritic cells match or exceed the tumour states |
| AREG or HBEGF is tumour-enriched in human lung relative to matched normal | Refuted at compartment level; not computable at subtype level |
| Runx1 and Pdgfrb persistence after Areg deletion reflects a second tumour signal | Refuted by E4; both genes rise with bleomycin alone |
| Hbegf or Egfr is specific to the tumour mesenchyme | Refuted by E4; both are injury-generic |
| A within-tumour comparison of tumour epithelial states against AT2 in GSE131907 | Not establishable; the deposited labels are tissue-exclusive |

---

## E6, C7 and C8: three questions asked properly, and three rules that showed their own limits (2026-09-13)

These are the three jobs the owner chose from the next-step list: a donor-level
test of the paper's axis in human fibrosis, the identity of the mesenchymal-sort
contaminant (item A1), and whether the Areg-independent tier is a subpopulation
or a gradient (item A2, re-aimed after E4 removed its original motivation).

All three returned a negative or an unresolved result, and in two of them the
frozen rule turned out to be the thing that failed rather than the biology.
Both defects are disclosed here and neither threshold was moved.

### E6: the axis shows no donor-level coupling that this cohort could detect

Script [`trials/e6_donor_level_axis_coupling.py`](trials/e6_donor_level_axis_coupling.py);
artefacts in [`trials/e6_donor_level_axis_coupling/`](trials/e6_donor_level_axis_coupling/).
GSE136831, 22 donors clearing the 50-cell floor in both compartments, 18 with
fibrosis and 4 control.

The question was the weakest one that still tests the link rather than
describing the parts: if epithelial AREG drives fibroblast EGFR in human
fibrotic lung, donors with more of the first should carry more of the second.

| Test | rho | p | Reading |
|---|--:|--:|---|
| T8, epithelial AREG against fibroblast EGFR | 0.348 | 0.112 | not significant |
| T9, epithelial AREG against fibroblast activation | -0.150 | 0.506 | not significant |
| T12, T8 within fibrosis donors only | 0.276 | 0.268 | not significant |
| T12, T9 within fibrosis donors only | -0.013 | 0.958 | not significant |

The pre-registered reading for that outcome requires reporting what the test
could have seen rather than implying the axis is absent. At 22 donors a
two-sided Spearman calls |rho| of about 0.43 and no less; within the 18
fibrosis donors, about 0.48. So this is a genuine null for strong coupling and
says nothing about weak coupling. It is also not evidence against the paper,
whose claim is about a state inside the epithelium: pooling all epithelium
dilutes it, which is the same limit trial E1b documented.

**The depth control earned its place.** The one correlation that reached
significance was a control, not the primary: epithelial TGFA against fibroblast
activation, rho 0.432, p 0.045. Both of those variables correlate with their
own compartment's sequencing depth, TGFA at 0.507 and activation at 0.412, so
the frozen rule marks the pair confounded and refuses to read it. Had the same
number appeared on AREG it would have looked like the result the trial was
hunting. Epithelial AREG, by contrast, is clean against depth at 0.124, so the
primary null is not a depth artefact either.

**A defect in the implementation, not the rule.** The rule said any pair whose
two members both correlate with depth is not read. The first version of the
script applied that only to the primary tests, so the significant control pair
came out unmarked. The rule was not changed; the code was corrected to apply it
to the control pairs as the text already required, and the trial was re-run
from the cached gene vectors.

### C7: the sort contaminant is a mixture, and whole-profile correlation cannot name a state

Script [`trials/c7_what_the_sort_contaminant_is.py`](trials/c7_what_the_sort_contaminant_is.py);
artefacts in [`trials/c7_what_the_sort_contaminant_is/`](trials/c7_what_the_sort_contaminant_is/).
Reference: the RFP-sorted epithelium of GSE316244, 11,350 cells in 11 clusters
after the C2 pipeline, named by the paper's own Figure 4l marker sets. Query:
the 184 cells of cluster 11 in the cached C1b mesenchyme object.

**T3, the control, passed and is the reason anything here is readable.**
Cluster 11's best correlation against an epithelial cluster is 0.878, while
cluster 14, the cleanest fibroblast cluster in the same object, reaches only
0.463. A gap of 0.415 says the measure separates compartments.

**T1 named a winner it cannot defend, which is a defect in how T1 was
written.** The top match is an AT1-like cluster at 0.8778, and the margin to
second is 0.0056. The top four matches span three different state calls inside
0.015 of rho:

| Reference cluster | Call | rho |
|---|---|--:|
| 5 | AT1_like | 0.8778 |
| 7 | DATP_like | 0.8722 |
| 0 | DATP_like | 0.8657 |
| 8 | AT2 | 0.8631 |

The rule required the margin to be reported and did not require it to be large
before a state could be named, so the script named one. The margin is what
should be read, and it says the measure has compartment resolution and no state
resolution: a rank correlation over 32,163 shared genes is dominated by what
all epithelium shares. The threshold was not moved after the fact, and the
sentence "the top state is AT1-like" is not a result of this trial.

**T2, which does have state resolution, says mixture.** Scoring the 184 query
cells directly with the same marker sets gives DATP-like 77, AT2 70, cycling
20, AT1-like 13 and Cd177-positive 4. The modal call is DATP-like at 41.8 per
cent of cells, below the 50 per cent floor this repository uses elsewhere, with
AT2 at 38 per cent right behind it.

So the answer to item A1 is that the leak is not state-selective. Cluster 11 is
roughly four parts DATP-like to four parts AT2, with a cycling minority. The
caution from finding one survives in a weaker and more precise form: about four
in ten of the contaminating cells score as the paper's own signalling
population, so a ligand-receptor analysis run inside GSE316241 alone would
indeed be reading some of that population as a mesenchymal source, but it would
also be reading ordinary AT2 cells, and the contaminant cannot be described as
the DATP-like state.

### C8: the retained tier is not a subpopulation, and what deletion removes is a co-expressing pair

Script [`trials/c8_subpopulation_or_gradient.py`](trials/c8_subpopulation_or_gradient.py);
artefacts in [`trials/c8_subpopulation_or_gradient/`](trials/c8_subpopulation_or_gradient/).
Fibroblasts defined by the C6 compartment gate rather than by C2's clustering,
a difference disclosed in the script: 2,751 gated cells in Areg-flox/+ and
3,700 in Areg-flox/flox, at 3,637 and 3,612 median genes per cell, so the two
arms are unusually well depth-matched for this deposit.

**T1 says independence.** In the flox/flox fibroblasts, Runx1 and Pdgfrb
co-occur in 22.46 per cent of cells against 20.44 per cent expected if they
were independent, a ratio of 1.099 against a frozen band of 0.80 to 1.25. The
retained genes do not mark a shared set of cells. T4, the depth control, moves
that ratio by 0.043 across the median split, well inside the 0.25 limit, so the
number is not depth-driven.

**T3 favoured two components and should not have been believed.** The mixture
fit prefers two Gaussians over one by 731 units of BIC, which the frozen
reading treats as a discrete subpopulation, so the two tests disagreed and the
pre-registered outcome for disagreement is "unresolved". That is the honest
outcome of the rule as written, but the rule was poorly chosen and the reason is
arithmetic rather than biological. The score is the mean of log1p Runx1 and
log1p Pdgfrb, and 31.5 per cent of these cells detect neither gene, so the
distribution has a spike at exactly zero plus a continuum. A two-component
mixture will win on that shape whatever the biology is. Trial C1 used the same
test on a `score_genes` output, which is centred and continuous and not exposed
to this; reusing it on a raw two-gene mean is the error, and it is recorded
rather than repaired, because repairing it after seeing the result is what the
frozen-rule discipline exists to prevent. **T3 carries no information here**, so
T1 is the only readable answer: a gradient, not a subpopulation, which is how
trial C5 had already corrected the description of claim C29.

**The post hoc observation, which is the interesting part and is labelled as
post hoc.** T2 was included only to give T1 something to compare against, and
it shows the structure sits on the other tier. In the control arm the two genes
that fall hardest, Fst and Runx2, co-occur at 13.78 per cent against 7.94 per
cent expected, a ratio of 1.735: they mark the same cells. In the deletion arm
that co-organisation is gone, ratio 0.816, and only 24 of 3,700 cells carry
both. The retained pair, Runx1 with Pdgfrb, sits at 1.196 in the control arm
and 1.099 after deletion, inside the independence band in both.

Read together with E4, that suggests what Areg deletion removes is a
co-expressing Fst and Runx2 set of fibroblasts, while Runx1 and Pdgfrb are
spread independently across the compartment whether Areg is present or not, and
E4 already showed those two rise with injury alone. None of this was
pre-registered, the deletion-arm ratio rests on 24 cells, and there is one
library per genotype, so it is a lead for a future pre-registration and not a
claim.

### Claims from E6, C7 and C8

| Claim | Class |
|---|---|
| Epithelial AREG and fibroblast EGFR are correlated across donors in human fibrotic lung | Not established; 22 donors, rho 0.348, p = 0.112, and the test could only have detected rho of about 0.43 |
| Epithelial AREG and fibroblast activation are correlated across donors | Not established; rho -0.150, and within fibrosis donors -0.013 |
| The mesenchymal-sort contaminant is a single epithelial state | Refuted; 41.8% score DATP-like and 38% AT2, below this repository's own confidence floor |
| About four in ten contaminating cells score as the paper's DATP-like state | Descriptive only, and it is the surviving form of the finding-one caution |
| Whole-profile rank correlation can name an epithelial state | Refuted by its own margin: 0.415 between compartments, 0.0056 between states |
| The Areg-independent tier is a discrete subpopulation of fibroblasts | Not established; co-detection ratio 1.099 inside the independence band, and the mixture test that disagreed is uninformative on a zero-inflated score |
| Areg deletion removes a co-expressing Fst and Runx2 set rather than reducing one programme evenly | Exploratory, post hoc, and resting on 24 double-positive cells in one library |

---

## C9. The Fst and Runx2 population, pre-registered (2026-09-13)

Script [`trials/c9_the_fst_runx2_population.py`](trials/c9_the_fst_runx2_population.py);
artefacts in [`trials/c9_the_fst_runx2_population/`](trials/c9_the_fst_runx2_population/).
This is the pre-registered version of trial C8's one post hoc observation, run
on the owner's instruction.

The trial separated an answerable question from an unanswerable one up front.
Answerable: does a co-expressing Fst and Runx2 fibroblast population exist in
an independent injury dataset with more than one animal, and what else marks
it. Unanswerable, and stated rather than approximated: whether Areg deletion
depletes it, because the Cardoso deposit has one library per genotype and no
other Areg-flox fibroblast dataset exists.

Six libraries, all gated by the same compartment rule used in C6 and C8: the
two GSE316244 niche libraries and the four GSE132771 Col1a1-GFP libraries.

### The headline verdict, and why it cannot be taken at face value

The script reports T1 as failing, and that verdict stands in the record. It is
not, however, a fair description of what the data did, and the reason is a
defect in how the rules were written rather than anything in the result.

**T1 and T5 were not composable.** T1's criterion, as written in the
pre-registration, is a co-detection ratio of at least 1.25 together with a
permutation p below 0.05, in both bleomycin libraries. T5, written separately,
says a ratio that moves by more than 0.25 across a median depth split is
depth-driven and is not read. Both fired. The ratios are 2.84 and 2.32 with
permutation p of 0.005 and 0.010, so T1's own criterion is met; the depth
shifts are 1.10 and 0.45, so T5 says the ratio those criteria rest on must not
be read. A rule whose test depends on a quantity another rule forbids reading
cannot be evaluated, and the implementation resolved that by bundling the depth
control into T1 and calling it a failure. That is the conservative branch, so
it is the one that stands, but the pre-registered reading attached to it, "the
observation does not replicate", is simply not what happened. The threshold was
not moved and the verdict was not rewritten; what follows is what the numbers
actually support.

### What replicates, and what does not

| Library | Group | Gated cells | Fst | Runx2 | Double positive | Ratio | Permutation p | Depth shift |
|---|---|--:|--:|--:|--:|--:|--:|--:|
| Expt3_Het_niche | Areg-flox/+ | 2,751 | 0.364 | 0.219 | 379 | 1.73 | 0.005 | 0.06 |
| Expt3_Hom_niche | Areg-flox/flox | 3,700 | 0.132 | 0.060 | 24 | 0.82 | 0.945 | 0.17 |
| Bleo1_GFPp | bleomycin | 2,443 | 0.138 | 0.034 | 32 | 2.84 | 0.005 | 1.10 |
| Bleo2_GFPp | bleomycin | 2,056 | 0.172 | 0.045 | 37 | 2.32 | 0.010 | 0.45 |
| UT1_GFPp | untreated | 2,355 | 0.032 | 0.003 | 0 | 0.00 | 1.000 | 0.00 |
| UT2_GFPp | untreated | 1,913 | 0.031 | 0.003 | 0 | 0.00 | 1.000 | 0.00 |

**Above-chance co-occurrence replicates across animals.** The permutation null
permutes the Runx2 indicator within deciles of genes per cell, so it holds
depth structure fixed by construction. It gives p = 0.005 and p = 0.010 in the
two bleomycin animals, both below the frozen line, in a dataset with no
oncogene anywhere. That is the part of the C8 observation that survives, and it
survives in the injury comparison the paper itself used.

**The size of the effect does not replicate and should not be quoted.** The
ratios are 2.84 and 2.32, but they move by 1.10 and 0.45 across the median
depth split, far outside the 0.25 limit. The cause is almost certainly not
depth dependence: there are 32 and 37 double-positive cells in these
libraries, so halving them leaves about sixteen cells per half and the ratio
estimate is dominated by sampling noise. The depth control as specified cannot
tell depth dependence from small-number noise when the double-positive count is
around thirty, which is a limit on the control rather than a finding about the
biology. The Cardoso reference library, with 379 double-positive cells, has a
depth shift of 0.06 and passes comfortably.

**T2 is met, and it says less than it looks.** There are zero double-positive
cells in both untreated libraries against 32 and 37 in the bleomycin ones, so
the frozen both-replicates rule calls the population injury-induced. But
Runx2 detection in an untreated library is 0.003, so the expected number of
double-positive cells under independence is about a quarter of one cell. Zero
is what independence predicts there. T2 therefore restates the marginal rise of
both genes with injury, which trial E4 had already established, and is not
independent evidence that the two genes mark the same cells.

### What marks those cells, and how far that travels

In the Areg-flox/+ library, the 379 double-positive fibroblasts differ from the
910 double-negative ones across a broad set of genes. The strongest by
detection difference are Piezo2 (0.879 against 0.184), Kif26b, Ltbp2, Chst11,
Rnf149, Basp1, Lhfpl2, Dock5, Bmper and Rftn1. This is one library and cells
are not replicates, so it is a description and no P value was computed.

The replication check across the two bleomycin animals is the informative part,
and it has to be read in two pieces.

| Reading | Result |
|---|---|
| Genes evaluable in both bleomycin libraries | 27 of 30 |
| Same direction in both animals | 26 of 27 |
| Same direction and at least half the magnitude | 9 of 27 |

**The direction of the signature replicates almost completely**, 26 of 27
evaluable genes in both animals independently. **The magnitude criterion is
uninformative and should not be quoted as "9 of 30 replicate".** The median
gene's bleomycin difference is 0.535 of its reference difference, and the
threshold was set at 0.5, so the criterion is splitting the distribution
essentially at its own centre. The attenuation has an obvious cause: the
bleomycin libraries carry median 1,799 and 1,938 genes per cell against 3,637
in the reference, so detection differences compress by roughly the same factor
as the depth. A threshold set without reference to that is measuring library
depth, not reproducibility. This is the third rule defect in this folder that
drawing the numbers out exposed, after C7's missing margin requirement and C8's
mixture test.

The nine genes that clear it anyway are Piezo2, Ltbp2, Rnf149, Basp1, Sdc1,
P4ha3, Prrx2, Megf11 and Cotl1. Taken with the broader direction-replicating
set, that reads as a matrix and mechanosensor-gene signature rather than
anything specific to the Areg axis: Piezo2 encodes a mechanosensitive channel,
Ltbp2 sequesters latent TGF-beta, P4ha3 hydroxylates collagen prolines, Sdc1 is
a matrix co-receptor, and Prrx2 is a mesenchymal transcription factor. Calling that a mechanics programme would claim more than the data support: transcript detection is not channel activity, dissociated tissue carries no mechanical context, and dissociation is itself a mechanical and enzymatic insult, so a Piezo2-led signature has a protocol confound that no reanalysis of dissociated data can remove. That is a description of a cell state, offered as
Exploratory and as a target for a future pre-registration, not as a mechanism.

Three of the thirty genes could not be checked at all, because Cemip2, AI506816
and 1110038B12Rik are absent from the older annotation GSE132771 was aligned
against. Gene-symbol drift between deposits of different vintages costs a tenth
of this comparison, and any future cross-deposit marker check should intersect
symbol vocabularies before choosing its top set.

### Claims from C9

| Claim | Class |
|---|---|
| Fst and Runx2 co-occur above chance in fibroblasts of bleomycin-injured lung, with depth held fixed by the permutation null, in both animals (p = 0.005 and 0.010) | Descriptive only, and it replicates across animals in a dataset with no oncogene |
| The magnitude of that co-occurrence (ratios 2.84 and 2.32) | Not established; the ratios move by 1.10 and 0.45 across the depth split, and about thirty double-positive cells cannot support a stable estimate |
| The population is injury-induced | Descriptive only, and weaker than it looks: zero double-positive cells is what independence predicts when Runx2 sits at 0.003 |
| The double-positive marker signature replicates in direction across two bleomycin animals, 26 of 27 evaluable genes | Descriptive only, and the strongest thing this trial produced |
| Nine of thirty markers replicate in magnitude | Not established as a meaningful figure; median attenuation is 0.535 against a threshold of 0.5, and it tracks the twofold depth difference between deposits |
| Areg deletion depletes this population | Not establishable; one library per genotype and no other Areg-flox fibroblast dataset exists |
| The trial's own T1 verdict, read as "the observation does not replicate" | Retracted as a reading; T1 and T5 were not composable and the conservative branch was taken, but the permutation result replicates in both animals |

---

## C10. The published pathological fibroblast, or not (2026-09-13)

Script [`trials/c10_published_state_or_not.py`](trials/c10_published_state_or_not.py);
artefacts in [`trials/c10_published_state_or_not/`](trials/c10_published_state_or_not/).
Run on the owner's instruction, to settle the question trial C9 raised and
could not answer.

### Why C9 could not answer it

None of the canonical pathological-fibroblast markers appeared in C9's top
thirty discriminating genes, which looked like novelty. It was not evidence.
Ranking by detection difference is biased twice over: Col1a1 is a condition of
the fibroblast gate, so it sits at 100 per cent in both groups and cannot rank
at all, and any gene that saturates or is very sparse is suppressed the same
way. This trial scores the published definitions directly instead.

The null was deliberately strong. Tsukui et al. 2020
(doi:10.1038/s41467-020-15647-5) define their pathological fibroblast by Cthrc1
with Postn, Spp1, Fn1 and Tnc, and their own text lists **Fst** among that
cluster's markers. Fst is one half of C9's group definition, which is why it
was excluded from the scored set here, and why the comparison was worth making.
Fang et al. 2025 (doi:10.1038/s41586-024-08542-2) then showed RUNX2 drives the
alveolar-to-pathological transition, so a Runx2-positive fibroblast being a
pathological fibroblast was the published expectation.

### Outcome: the lead closes

| Library | Group | Double positive | pathological | alveolar | adventitial | smooth muscle |
|---|---|--:|--:|--:|--:|--:|
| Expt3_Het_niche | Areg-flox/+ | 379 | **1.81** | -0.76 | -0.59 | 0.06 |
| Bleo1_GFPp | bleomycin | 32 | **1.49** | -1.03 | -0.17 | 0.28 |
| Bleo2_GFPp | bleomycin | 37 | **1.46** | -1.26 | -0.23 | 0.05 |

Standardised differences between double-positive and double-negative gated
fibroblasts. T1 clears its 0.5 floor in all three libraries, and T4 puts the
pathological set first everywhere. The pre-registered reading therefore
applies: **these are the published pathological fibroblast, and trial C9's
signature is that state reached by an unusual route.**

Two things make the pattern hard to explain any other way. Alveolar identity is
strongly *negative* in the double-positives, which is exactly what Tsukui
describe for a state arising from alveolar fibroblasts that lose their
identity. And the smooth-muscle score is flat at 0.06, 0.28 and 0.05, so this
is not the mural contamination that trials C1b and E4 found distorting other
measures.

### Cthrc1 was there the whole time

| Library | Cthrc1, double positive | Cthrc1, double negative |
|---|--:|--:|
| Expt3_Het_niche | 0.449 | 0.020 |
| Bleo1_GFPp | 0.688 | 0.133 |
| Bleo2_GFPp | 0.622 | 0.058 |

A twenty-two-fold difference in the reference library. It missed C9's top
thirty because its detection difference, 0.428, fell below that list's cutoff
of 0.510. The saturation audit also confirms the bias directly: Col3a1 is
formally unrankable, at 1.000 against 0.987.

### Two honest caveats, and the fourth disclosed rule defect

**The double-positive call is partly a depth call.** Detecting two sparse genes
in the same cell needs depth, and the double-positives concentrate in the deep
half: 342 of 379 in the reference, 31 of 32 and 33 of 37 in the bleomycin
libraries. That does not explain the result, though, and the reason is the
opposite signs. Depth cannot push the pathological score up while pushing the
alveolar score down by a similar amount, and `score_genes` subtracts a matched
control set besides. The direction disagreement between the two sets is what
rules depth out.

**T3 says graded, not discrete.** Blind clustering of the reference library
does not put the double-positives in one cluster at the frozen thresholds, so
they are a graded state, consistent with what C8 found for the retained tier.

**The rule defect.** T5 as written triggers on a computed half that changes
sign or falls below half the full value. In the bleomycin libraries the shallow
half holds one and four double-positive cells, so the statistic cannot be
computed there at all, which the written rule does not cover. The first run
treated "not computable" as "failed" and therefore suppressed a reading that
the reference library supports cleanly, with its own depth control passing at
1.345 against 1.813. The code was corrected to match the written rule and the
three outcomes are now recorded separately. The threshold was not moved. This
is the fourth rule defect this folder has disclosed, after C7's missing margin
requirement, C8's mixture test on a zero-inflated score, and C9's
non-composable T1 and T5.

### What this settles elsewhere

It explains trial E4. The pathological fibroblast is by definition an injury
state, so of course Fst, Runx2, Tnc and the rest are injury-generic; that
result and this one are the same fact seen twice.

It also supersedes claim C70. The signature is not a new programme, so the
question of what to call it is moot, and the matrix and mechanosensor genes are
simply part of a known state's transcriptome.

### Claims from C10

| Claim | Class |
|---|---|
| The Fst and Runx2 double-positive fibroblasts are the published Cthrc1-positive pathological fibroblast | Descriptive only, and agreed by three libraries including two independent animals |
| They have lost alveolar identity, matching the published account of that state's origin | Descriptive only |
| Cthrc1 was absent from C9's top thirty because of a ranking cutoff, not because it fails to discriminate | Descriptive only, and the method caution the trial existed to produce |
| Trial C9's signature is a new fibroblast programme | Refuted; it is a known state reached by an unusual route |
| The double-positives form a distinct population | Not established; graded at the frozen thresholds |
| The result is a depth artefact | Refuted by the opposite signs: depth cannot raise one score while lowering another |

---

## C11. Figures for the contradictions (2026-09-13)

Script [`trials/c11_figures_for_the_contradictions.py`](trials/c11_figures_for_the_contradictions.py);
figures in [`trials/c11_figures_for_the_contradictions/`](trials/c11_figures_for_the_contradictions/).
Reads only tracked tables, computes nothing scientific, and each figure carries
a fail-closed assertion that refuses to draw if its source table stops saying
what the title claims.

Trials C5 and E5 drew the findings and the one refutation that existed then.
Six results since then refute an earlier trial, refute one of this repository's
own methods, or close a lead, and each lived only as a number in a table.

| Figure | What it shows |
|---|---|
| `c11_fig1_e6_depth_control.png` | The only significant correlation in trial E6 was a control pair, and both of its variables track sequencing depth |
| `c11_fig2_c7_resolution.png` | Profile correlation separates compartments by 0.415 of rho and states by 0.0056, and the sort contaminant is a mixture |
| `c11_fig3_c8_amplitudes.png` | The retained pair sits inside the independence band, and 31.5 per cent of cells detect neither gene |
| `c11_fig4_species_divergence.png` | HBEGF is topped by dendritic cells and macrophages in human and by epithelium and endothelium in mouse |
| `c11_fig5_c9_existence_vs_size.png` | Co-occurrence replicates in both animals; the ratio is unstable across the depth split, and the magnitude threshold sits at the median attenuation |
| `c11_fig6_c10_lead_closes.png` | The pathological set is the largest effect while alveolar identity is lost, and Cthrc1 sat just under C9's ranking cutoff |

Three wordings were corrected because drawing them exposed an overstatement, a
pattern this folder has now seen at every figure stage. A panel titled "every
epithelial state scores the same" became "the top four matches span three
different state calls", because the lower-scoring clusters genuinely differ. A
subtitle saying myeloid cells "sit at the bottom in mouse" became a statement
naming the two compartments that do and the exception that does not, since
mouse neutrophils rank third. And two animals of one group had been given one
colour while appearing as separate legend entries, which the palette rule
forbids.

---

## C12. The whole CellChatDB, ranked (2026-09-13)

Script [`trials/c12_cellchatdb_full_resource_scan.py`](trials/c12_cellchatdb_full_resource_scan.py);
artefacts in [`trials/c12_cellchatdb_full_resource_scan/`](trials/c12_cellchatdb_full_resource_scan/).
Run on the owner's instruction, as the cheap substitute for installing R.

### What this is, and what it is not

Trial C3 scored four EGFR ligands because those are the ones the paper named.
This trial scores the authors' entire curated resource, so a pair nobody here
thought to look at had a chance to appear. It uses CellChat's own resource and
CellChat's own scoring logic through LIANA's implementation, on human data, and
reads rankings rather than significance. **The permutation test was switched off
deliberately**: its unit is the cell, and a cell-level P value beside a group
question is the number this repository refuses to emit. It is not a CellChat
rerun of the paper's analysis, and CellChat itself has still never run here.

The mouse deposit was not used, for four reasons stated before the run rather
than discovered during it. The paper's own configuration integrates GSE247505
epithelium with the niche, and trial C0 found those sit in different gene
spaces. Ligand and receptor would then come from different libraries at
different depths, which confounds every product-based score the method
computes. Each genotype contributes one library of three pooled mice, so no
statistic has a unit. And running it on the mesenchymal library alone would put
the 184 Areg-high epithelial contaminants of trial C1d on the epithelial side
of the scan, which is the artefact this repository spent three trials
characterising.

GSE136831 avoids all four: epithelium and fibroblasts come from the same donor
at the same depth, and 26 donors cleared the 50-cell floor in both
compartments. 313 pairs survived the donor floor.

### Outcome: the guard fired, and that is the result

The pre-registered guard asked what fraction of the top fifteen pairs are
collagen or laminin ligands against integrin or syndecan receptors, because
those are abundant in every fibroblast rather than informative about
signalling. Eleven of fifteen qualified, and eleven of the top fifteen target
CD44 alone.

| Rank | Pair | Donors | Median rank |
|--:|---|--:|--:|
| 1 | COL4A2 to CD44 | 25 | 3.0 |
| 2 | COL6A2 to CD44 | 21 | 4.0 |
| 3 | FN1 to CD44 | 26 | 6.5 |
| 4 | COL4A1 to CD44 | 25 | 9.0 |
| 5 | LAMB3 to CD44 | 26 | 9.0 |
| 6 | APP to CD74 | 26 | 9.5 |
| 7 | MDK to LRP1 | 26 | 10.5 |

So **the full-resource ranking reports transcript abundance, and T1 and T2 are
not read.** That is the honest primary outcome, and it is a methods result
worth more than the ranking would have been: a whole-resource ligand-receptor
scan on dissociated data puts matrix proteins against a promiscuous receptor at
the top, because those ligands are the most abundant transcripts a fibroblast
has. The guard caught the known failure mode exactly as written.

One note on the guard's own arithmetic. It undercounts slightly, because it
tests for a collagen or laminin ligand or an integrin or syndecan receptor, and
FN1 to CD44 is plainly an abundance pair that matches neither clause. The
undercount is in the conservative direction, since the guard fired anyway.

### The one comparison that survives, with its rule conflict disclosed

T4 asked where the axis sits among the ligands of its own receptor, which is
the comparison trial C3 made in mouse.

| Rank | Pair | Donors | Median rank of 313 | Interquartile range |
|--:|---|--:|--:|---|
| 1 | AREG to EGFR | 26 | 15.5 | 10.2 to 40.0 |
| 2 | HBEGF to EGFR | 26 | 48.5 | 29.2 to 81.5 |
| 3 | TGFA to EGFR | 24 | 115.5 | 65.5 to 175.2 |
| 4 | EREG to EGFR | 22 | 131.5 | 83.5 to 192.2 |
| 5 | BTC to EGFR | 22 | 201.5 | 155.8 to 244.0 |

**AREG is first among the five, across 26 donors, by the authors' own resource
and scoring, with the donor as the unit.** The top two match the mouse order
trial C3 found; EREG and TGFA swap places between the species.

The rule conflict, which is the fifth of its shape in this folder and is
disclosed rather than resolved quietly. T3's own text scopes the guard to T1
and T2. The reading section, written in the same docstring, says that if T3
fails then "nothing else is read", which would include T4. The specific rule
and the general gloss disagree.

T4 is reported here for two reasons, and the reader should weigh them. The
specific rule names T1 and T2 and not T4. And comparing five ligands of one
receptor is structurally immune to what the guard detects, since the
promiscuity of CD44 cannot affect a ranking in which the receptor is held
constant. Against that, the conflict means T4 was not unambiguously
pre-registered as readable, so it should be treated as provisional until a
trial scopes the guard properly. It is worth saying that reading it does not
flatter any hypothesis of this repository's: it vindicates the paper's choice
of shortlist, which is not a conclusion anything here was pushing for.

### What it settles

Nothing in the resource outranks the paper's axis in a way that survives the
guard. Where the comparison is clean, the paper's selection of AREG from among
the EGFR ligands is the right one in human fibrosis too. The discovery
affordance that justified this trial, the chance that an unconsidered pair
would appear, produced matrix-to-CD44 pairs and nothing more.

### Claims from C12

| Claim | Class |
|---|---|
| A whole-resource ligand-receptor ranking on dissociated lung reports transcript abundance: 11 of the top 15 pairs target CD44 | Descriptive only, and the methods result the trial actually produced |
| Some pair in CellChatDB outranks the AREG to EGFR axis | Not established; the ranking that would answer it is unreadable by its own pre-registered guard |
| Among the five EGFR ligands, AREG ranks first in human fibrosis by CellChat's resource and scoring across 26 donors | Descriptive only, and provisional because of a disclosed conflict between the guard's scope and the reading gloss |
| The mouse ligand order transfers to human | Partly: the top two match, and EREG and TGFA swap |
| Any of this is evidence that two cells communicate | Not established, and not establishable without proximity |

---

## D0. The Choi 2020 deposit check, now in the paper-2 folder

Trial D0 ran on 2026-09-13 as roadmap paper 2's Gate 0, sat in this folder
for part of 2026-09-15 while the paper-2 folder was withdrawn, and returned
to [`../gate1_02_choi_2020/`](../gate1_02_choi_2020/ANALYSIS_TRIAL_PLAN.md)
when the owner re-entered the paper the same day. Nothing in it changed; it
keeps its identifier and its claims rows C58 to C64. What it established
still governs this folder's reading of the DATP state: one library per
condition, six raw whitelists, no counted reporter, coverage-only ATAC.
