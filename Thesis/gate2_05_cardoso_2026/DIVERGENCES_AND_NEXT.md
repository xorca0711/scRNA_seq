# Why this reanalysis diverges from the paper, what the biology behind each divergence is, and what to do next

Companion to [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md), which records what
each gate did and returned. This page answers a different question: given the
same deposited counts, **why** does an independent pipeline land somewhere
slightly different, which of those differences are biology and which are
bookkeeping, and what would settle each one.

Every number here is read from a trial artefact under [`trials/`](trials/).
Status of everything below: Descriptive only or Exploratory, because no
genotype contrast in this paper's deposit carries within-group replication
(trial C0). Owner review pending.

Figures for the three findings are in
[`trials/c5_figures_for_the_three_findings/`](trials/c5_figures_for_the_three_findings/).

According to PubMed for all external citations, with DOIs given inline.

---

## 1. The seven mechanical reasons a reanalysis moves

None of these is a disagreement with the authors. Each is a place where a
decision had to be made, the paper made one, this repository made another, and
the difference is visible in the output.

| # | The decision | The paper | Here | What it changed |
|---|---|---|---|---|
| 1 | Who annotates | expert annotation of Seurat clusters, with unwanted types removed manually | a score over the paper's own published marker sets, assigned by modal argmax, labels never used to fit anything | the caller has no expert veto, so a cluster whose markers are ambiguous stays ambiguous instead of being named |
| 2 | Clustering and integration | Louvain on Harmony-integrated objects | Leiden (igraph, seed 0) with no correction, because library and genotype are one variable here | cluster boundaries and counts differ; the fibrotic population splits across clusters rather than forming one |
| 3 | Quality thresholds | fixed cuts (>1,000 genes, 2,000 to 50,000 UMIs, <10% mitochondrial) | per-library MAD, the repository's standing rule | the paper's rule would drop 29% of the Red2Kras library against 16% of the Confetti one (trial C1), an asymmetry that lands directly in any composition comparison |
| 4 | Doublets | no dedicated caller named | Scrublet per library, with the automatic threshold rejected when it disagrees with the 10x prior | on these libraries the automatic threshold called 0.00 to 0.07%, so the calls fell back to a rank cut: a ranking, not a detection |
| 5 | What counts as a population | expert judgement | a 50% modal-confidence floor | twice the floor scored a real effect as nothing (C1 dropped the candidate clusters; C2 reported 0.0% against 0.0%), which is why C1b and C2b exist |
| 6 | Whether a marker set is compartment-pure | not an issue when a human annotates | it is the whole caller | two of the six genes in the paper's fibrotic set (Acta2, Pdgfrb) are mural markers and a third (Runx1) is also myeloid, so the score peaks on smooth muscle and flags a myeloid contaminant (trial C1b, Question C) |
| 7 | Sequencing depth across compared libraries | handled inside Seurat integration | detection fractions compared directly | the Areg-arm libraries are deeper (median 2,615 and 2,750 genes per cell) than the Confetti and Red2Kras ones (1,824 and 2,111), so cross-series detection is not comparable and only the within-series contrast carries a claim |

The honest summary of the table: **six of the seven differences are
bookkeeping, and they explain essentially all of the divergence in Gate 1.**
The one that is not bookkeeping is number 6, because it says something about
the published marker set that will bite anyone who reuses it.

---

## 2. Finding one: an Areg-high epithelial population inside the mesenchymal sort

![The sort contaminant](trials/c5_figures_for_the_three_findings/c5_fig1_sort_contaminant.png)

**What the data show.** GSE316241 is deposited as CD45-CD31-EpCAM- mesenchyme.
Of its 11,690 cells, 757 (6.5%) sit in clusters that read as off-target for
that sort (trial C1d). Cluster 11 is the one that matters: 184 cells, 99% from
the Red2Kras library, with epithelial keratins in more than 90% of cells and
**Areg detected in 88%**. It is not a generic epithelial leak. It carries the
paper's own DATP-like markers (Cldn4 60%, Itga2 63%, Sox9 48%) alongside AT2
and AT1 markers, which is what the mutant epithelium looks like in the
lineage-labelled series.

**Why it could be there: three explanations, and what the check says.** Trial
C5 measured the markers that separate them.

| Explanation | What it predicts | What cluster 11 shows | Reading |
|---|---|---|---|
| Mutant cells lose Epcam and slip through the EpCAM- gate | low Epcam | Epcam mRNA in 78% of cells | **not supported at the RNA level**; surface protein can still be lower than mRNA, which 3' counts cannot see |
| Epithelial-fibroblast doublets | epithelial and fibroblast transcripts in the same cell | Krt8 and Col1a1 together in 19%; median Scrublet score 0.062 against 0.035 for the object | **a minority, about one cell in five** |
| Ambient RNA from lysed epithelium | low-level epithelial reads spread across all clusters | Krt8 in 0.3% of fibroblast cluster 14 | **not supported**; ambient contamination would not concentrate in one cluster |

The most parsimonious reading is ordinary sort impurity, with a doublet
minority. FACS purity is never 100%, and the DATP-like mutant population is
expanding fast at two weeks, so a small leak rate produces a visible cluster
almost entirely in the tumour arm.

**Why it matters, and why the paper is not affected.** The paper computed
epithelium-to-fibroblast signalling on epithelium imported from the separate
lineage-labelled series, so this contaminant never entered its CellChat
object. A reanalysis that ran ligand-receptor inference inside GSE316241 alone
would be computing Areg signalling from cells the sort was supposed to
exclude, and would find a tumour-specific Areg source that is an artefact of
the gate. That is the practical caution, and it generalises: a sorted
compartment in a tumour lung will carry the fastest-expanding population of
the compartment it was sorted against.

**What would settle it.** Flow cytometry of EpCAM surface protein on
lineage-labelled DATP-like cells, which the lab could answer from existing
panels. From the data in hand: map cluster 11 onto the RFP-sorted epithelium
of GSE316244, which shares its gene space exactly, and ask which mutant state
it matches.

---

## 3. Finding two: the fibrotic programme does not fall as one unit

![The fibrotic programme under ligand deletion](trials/c5_figures_for_the_three_findings/c5_fig2_fibrotic_programme_split.png)

**What the data show.** In the fibroblasts of GSE316244, comparing
Areg-flox/flox with Areg-flox/+ at similar depth, the six genes of the paper's
own reprogrammed-fibroblast set behave in three tiers rather than two:

| Gene | detection kept after Areg deletion | reading |
|---|--:|---|
| Runx1 | 97% | retained |
| Pdgfrb | 86% | retained |
| Acta2 | 77% | falls by about a quarter |
| Tnc | 76% | falls by about a quarter |
| Fst | 42% | falls by more than half |
| Runx2 | 39% | falls by more than half |

C2's frozen 80% line makes this look binary. It is a gradient, and saying so
is part of the result: the line is where a pre-registered threshold happened to
fall, not where the biology has a joint.

**The biological logic.** The tiers line up with what the genes do rather than
with the label "fibrotic". Runx2 and Fst, the two that mostly disappear, are
the ones tied to the pathological fibroblast programme itself: RUNX2 was shown
to drive an alveolar-to-pathological fibroblast transition and its deletion
reduces pathological fibroblasts and fibrosis (Fang et al., *Nature* 2025,
[DOI](https://doi.org/10.1038/s41586-024-08542-2); note the author correction,
[DOI](https://doi.org/10.1038/s41586-025-09254-x)). Tnc and Acta2, which fall
by about a quarter, are matrix and contractile effectors. Runx1 and Pdgfrb,
which barely move, are the two that other inputs are known to drive: the
alveolar fibroblast lineage passes through an inflammatory state into a
fibrotic one under TGF-beta control, with inflammatory and fibrotic states
separately regulated (Tsukui et al., *Nature* 2024,
[DOI](https://doi.org/10.1038/s41586-024-07660-1)), and Pdgfrb is a receptor rather than a
matrix effector, so its induction reports that the cell has entered an
activated mesenchymal state without saying which ligand put it there.

So the reading that fits both the tiers and the literature is that **Areg-EGFR
is necessary for the matrix and pathological-programme half of fibroblast
reprogramming, and something else drives the Pdgfrb and Runx1 half.** The
paper's own model has candidates for that something else, because its
macrophage arm supplies IL-1beta and its fibroblasts sit next to mutant cells
that make more than one ligand.

**The competing explanation that has to be excluded first.** A retained
detection fraction is also what incomplete deletion looks like. Areg-flox/flox
is a conditional deletion in Sftpc-CreERT2-labelled AT2 cells, so any mutant
cell that escaped recombination still makes Areg, and any fibroblast next to
one still sees it. If escape is common, every gene's retention is inflated and
the tiers could be a dose-response to residual ligand rather than two
pathways. This is checkable in the data: Areg detection in the RFP-sorted
epithelium of the same flox/flox library puts a bound on escape.

**What would settle it.** Within flox/flox fibroblasts, is there a
Runx1-positive Pdgfrb-positive Tnc-negative population, or do all fibroblasts
shift down together? A per-cell co-detection test separates "two programmes"
from "one programme at lower amplitude". That needs the C2 object kept, which
means one re-run.

---

## 4. Finding three: Hbegf is second by abundance, and that is not the same as second in importance

![EGFR ligands, abundance against enrichment](trials/c5_figures_for_the_three_findings/c5_fig3_egfr_ligands_abundance_vs_enrichment.png)

**What the data show, and the correction this forced.** Ranking the four EGFR
ligands by mean expression inside the DATP-like state gives the same order in
all four mutant libraries: Areg, then **Hbegf**, then Ereg, then Tgfa. That was
trial C3's frozen test and it holds. But ranking them by enrichment over the
AT2 cells of the same library gives Areg first and **Ereg** second in three of
the four libraries, with Hbegf third; by the ratio of detection fractions,
Ereg is first in every library and Hbegf is last or second to last. Hbegf is
abundant in the DATP-like state because it is broadly expressed in mutant AT2
cells too.

The earlier session summary said Hbegf "ranks ahead of the ligand the paper
followed into culture". That is true of abundance and misleading as a statement
about the state, and this file is where the correction lives.

**The biological logic, which makes the correction more interesting than the
original claim.** EGFR ligands are not interchangeable at equal abundance,
because they determine what happens to the receptor after binding. HB-EGF
targets essentially all engaged EGFR for lysosomal degradation, whereas
amphiregulin does not target the receptor for degradation and instead drives
recycling, and epiregulin drives complete recycling (Roepstorff et al.,
*Traffic* 2009, [DOI](https://doi.org/10.1111/j.1600-0854.2009.00943.x)).
Epiregulin is also a partial agonist of EGFR dimerisation whose weaker dimers
paradoxically produce **more sustained** signalling and differentiation-type
responses rather than proliferation (Freed et al., *Cell* 2017,
[DOI](https://doi.org/10.1016/j.cell.2017.09.017)).

Read against that, the ranking reverses in meaning. A niche that needs
**sustained** fibroblast EGFR output over days to weeks is better served by a
recycling ligand such as Areg or Ereg than by an abundant degradative one such
as Hbegf, whose own abundance would drive receptor downregulation. The paper's
choice to pursue Areg, and to test Areg together with Ereg in culture, is
consistent with the trafficking biology rather than with the abundance
ranking, and the independent precedent is direct: sustained AREG from
intermediate alveolar stem cells activates fibroblast EGFR and drives
progressive fibrosis, and an anti-AREG antibody blocks it (Zhao et al., *Cell
Stem Cell* 2024, [DOI](https://doi.org/10.1016/j.stem.2024.07.004)).

**So what is Hbegf doing.** Two possibilities, and the data here cannot
separate them. It may be a passenger: a stress-induced EGFR ligand that comes
up with the mutant state and contributes little because it shuts the receptor
down. Or it may be the brake: an abundant degradative ligand that limits how
much fibroblast EGFR signal the niche can sustain, in which case removing
Hbegf would make the Areg axis stronger rather than weaker. The second
possibility is the one worth an experiment, because it predicts the opposite
sign from the obvious one.

**Two further observations from the same table, both thin.** The wild-type YFP
clones of the same 2-week animals contain a few DATP-like cells (21 and 36),
and their ligand order is the same, which hints the ligand profile belongs to
the state rather than to the Kras mutation; 21 cells cannot carry that. And the
two 4-day replicates disagree about 24-fold on how much of the library is
DATP-like (1.2% against 28.2%), so C3's 14.7% median at 4 days describes
neither library, and the 4-day arm of every ligand statement rests on 44 cells
in one of its two replicates.

---

## 5. What to do next, in the order the cost-to-value ratio suggests

### A. Answerable from data already on disk, no new experiment

| # | Question | Method | Cost | What a result would mean |
|---|---|---|---|---|
| A1 | Is the mesenchymal-sort contaminant the same state as the RFP-sorted mutant epithelium? | score cluster 11 against the GSE316244 RFP+ states; the two share a gene space exactly, so no intersection is needed | one short trial | if it matches DATP-like, the contaminant is the paper's own signalling population and the sort caution is sharper; if it does not, the leak is something else |
| A2 | Is the fibrotic gradient one programme at lower amplitude, or two? | per-cell co-detection of Runx1, Pdgfrb and Tnc inside flox/flox fibroblasts | one C2 re-run, keeping the object | decides whether "Areg-independent half" is a real subpopulation or a wording artefact |
| A3 | How much residual Areg is there in the flox/flox arm? | Areg detection in the RFP-sorted flox/flox epithelium | minutes, data in hand | bounds the incomplete-deletion explanation that currently competes with A2 |
| A4 | Does the retention result survive depth matching? | downsample the deeper library's counts to the shallower library's median and recompute retention | one short trial | removes the only within-series confound left on claim C29 |
| A5 | Is the Hbegf abundance pattern present in human LUAD? | the paper's own comparison dataset, Kim et al. 2020, GSE131907, public | one trial, one download | a conserved abundance-versus-enrichment split would raise Hbegf from curiosity to candidate |
| A6 | Does the non-oncogenic injury fibroblast state show the same Runx1 and Pdgfrb retention? | Tsukui et al. bleomycin atlas, GSE132771, public | one trial, one download | separates "Areg-independent" from "injury-generic"; this is the comparison the paper itself used for the fibroblast state |

### B. Needs the bench, and is worth asking Choi about

| # | Question | Experiment | Why it is worth doing |
|---|---|---|---|
| B1 | Do DATP-like cells carry less surface EpCAM than their mRNA suggests? | flow cytometry of EpCAM on lineage-labelled RFP+ cells | settles the sort question, and would explain why mutant epithelium appears in an EpCAM-negative gate in any such experiment |
| B2 | Is Hbegf a passenger or a brake on the Areg axis? | dose-matched recombinant AREG against HB-EGF on alveolar fibroblasts, reading Pdgfrb and Acta2 and EGFR surface level over time; then HB-EGF loss of function in the tri-culture | the trafficking literature predicts HB-EGF gives a transient signal and downregulates the receptor, so its removal could strengthen the niche; that is a falsifiable, counter-intuitive prediction |
| B3 | Which input drives Pdgfrb and Runx1 if not Areg? | the existing tri-culture with IL-1beta or TGF-beta blockade on the Areg-deficient background | would close the one real gap this reanalysis found in the paper's cascade |
| B4 | Does fibroblast reprogramming precede macrophage change transcriptionally? | scRNA-seq of the mesenchymal and immune compartments at 1, 2, 4 and 8 weeks | the deposit has no time course in those compartments, so the paper's ordering claim currently rests on imaging alone (trial C0) |

### C. Out of reach with this data type, and worth stating rather than attempting

- **Tnc isoform usage.** The variable FNIII domains of Tnc cannot be resolved
  from 3' short-read counts. Since Tnc is the ligand in the fibroblast to
  macrophage leg, which isoform is made is a real question, and it needs
  long-read or targeted sequencing. Stated as a limitation, not a plan.
- **Anything about who receives a signal.** No communication analysis was run
  here and CellChat cannot run on this machine. Every ligand statement in
  section 4 is about the sending cell only.
- **Any tested genotype difference in this paper's deposit.** One library per
  genotype, three mice pooled. The companion series is the only place with
  replication, and even there n is two libraries per arm.

---

## 6. Status of the claims on this page

| Claim | Status |
|---|---|
| The seven mechanical reasons, and that six are bookkeeping | Descriptive only, from the trial artefacts |
| The contaminant is sort impurity with a doublet minority, not Epcam loss or ambient RNA | Descriptive only; the RNA cannot exclude lower surface protein |
| The fibrotic response to ligand deletion is graded, in three tiers | Descriptive only (one library per genotype) |
| "Areg-EGFR drives the matrix half, another input drives the Pdgfrb and Runx1 half" | Exploratory, and incomplete deletion is not yet excluded (A3) |
| Hbegf is second by abundance and third by enrichment | Descriptive only |
| "Hbegf may be a brake rather than a passenger" | Exploratory hypothesis, from published trafficking biology, untested here |
| The 4-day replicates disagree about 24-fold | Descriptive only, and a caveat on every 4-day statement |
| Everything in section 5 | Proposed, not run |
