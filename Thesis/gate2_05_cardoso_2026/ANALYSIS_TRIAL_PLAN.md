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
| C0 | 0 | What is actually deposited, and what statistical unit can it carry? | run; Descriptive only | [`trials/c0_data_reality_check/`](trials/c0_data_reality_check/) |
| C1 | 1 | Does the Pdgfrb+Runx1+Tnc+ fibrotic fibroblast subset reproduce, and does it separate cleanly? | run; **not recovered** by the frozen rule | [`trials/c1_fibroblast_compartment/c1_summary.md`](trials/c1_fibroblast_compartment/c1_summary.md) |
| C1b | 1 | Then what is there instead? (stop-and-characterise, plus a disclosed corrected pass) | run | [`trials/c1b_characterise_red2kras_private/c1b_summary.md`](trials/c1b_characterise_red2kras_private/c1b_summary.md) |
| C1c | 1 | Are the fibrotic and inflammatory programmes in the same cells at 2 weeks? | script written, **not yet run** | [`trials/c1c_fibrotic_inflammatory_overlap.py`](trials/c1c_fibrotic_inflammatory_overlap.py) |
| C2 | 2b | Does removing Areg collapse the niche, and does anything survive it? | **running when the session paused** | [`trials/c2_areg_deletion_arm/`](trials/c2_areg_deletion_arm/) |
| C3 | 2a | Is Areg a property of the DATP-like state under my own clustering, across the time course? | script written, **not yet started** | [`trials/c3_areg_state_specificity.py`](trials/c3_areg_state_specificity.py) |
| C-d | 2d | Does the transcriptome agree that fibroblast reprogramming precedes macrophage change? | **closed by C0, not run** | see below |

> **Session state, 2026-09-12.** C0, C1 and C1b are complete and their outcomes
> are written below. C1c and C3 are written and not run; C2 was running when the
> session paused. The handoff block in [`PROGRESS.md`](../../PROGRESS.md) says how
> to resume each one.

---

## C0. Data reality check (Gate 0)

Run on 2026-09-12. Script [`trials/c0_data_reality_check.py`](trials/c0_data_reality_check.py);
artefacts in [`trials/c0_data_reality_check/`](trials/c0_data_reality_check/).

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
