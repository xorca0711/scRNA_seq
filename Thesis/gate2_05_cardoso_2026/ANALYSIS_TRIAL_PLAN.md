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
