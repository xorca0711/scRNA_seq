# Study note: Choi, Lee et al. 2020, Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors that mediate alveolar regeneration

Roadmap paper 2, Gate 1. Written 2026-09-13, after papers 1, 4 and 5, because
paper 5 was entered out of order on the owner's instruction. Metadata verified
against PubMed; the full text was read from the Europe PMC XML rather than the
PMC web rendering, because that rendering strips italicised gene symbols and
would have left every marker set in this note incomplete.

**Choi J, Park J-E, Tsagkogeorga G, Yanagita M, Koo B-K, Han N, Lee J-H.**
*Inflammatory signals induce AT2 cell-derived damage-associated transient
progenitors that mediate alveolar regeneration.*
**Cell Stem Cell** 2020;27(3):366-382.e7.
DOI: [10.1016/j.stem.2020.06.020](https://doi.org/10.1016/j.stem.2020.06.020) ·
PMID [32750316](https://pubmed.ncbi.nlm.nih.gov/32750316/) ·
PMC [PMC7487779](https://pmc.ncbi.nlm.nih.gov/articles/PMC7487779/)

**Why this paper is in the roadmap.** It is the origin of the DATP state that
three other papers in this repository lean on, and its first author leads one
of the target labs. Paper 5 (Cardoso 2026) carries a DATP-like state into
tumour initiation and cites this work for the state itself; paper 1
(Niethamer 2025) supplies the Krt8-high transitional axis this repository has
already re-derived. Reading this paper is what makes those two comparable,
because it defines the state, names the signal that induces it, and states the
condition under which the state becomes pathological rather than regenerative.

---

## 0. Reading workflow (the five questions the roadmap requires)

**Question tested.** During alveolar regeneration, what is the path from an
AT2 stem cell to a mature AT1 cell, what signal moves a cell along it, and what
happens to regeneration when that signal does not switch off?

**Evidence type.** Mixed, and weighted toward genetics rather than
transcriptomes. Single-cell RNA sequencing supplies the states and their
ordering in two systems, in vivo lineage-labelled epithelium and AT2 organoids.
Lineage tracing supplies the claim that DATPs come from AT2 cells and go on to
become AT1 cells, which no transcriptome can establish on its own.
Immunofluorescence and flow cytometry supply the abundances and the time
course. Organoid culture with recombinant IL-1beta supplies sufficiency, a
HIF1alpha inhibitor and a conditional deletion supply necessity, and an
Il1r1 reporter supplies the claim that only a subset of AT2 cells can respond.
Section 2 classifies each link, because a reanalysis of the deposit can only
re-ask the questions in the scRNA rows.

**Reusable variables.** The DATP marker set (Cldn4, Krt8, Ndrg1, Sprr1a and
AW112010) together with the paper's explicit negative condition, that DATPs
carry much lower canonical AT1 markers (Pdpn, Hopx, Cav1) than mature AT1
cells; the five-state vocabulary the paper resolves in both systems
(homeostatic AT2, cycling AT2, primed AT2, DATP, AT1); the primed-AT2
signature as a loss of AT2 identity genes (Etv5, Abca3, Cebpa down); the
DATP programme groupings the paper reports from Gene Ontology, namely p53
signalling (Trp53, Mdm2, Ccnd1, Gdf15), proliferation arrest (Cdkn1a, Cdkn2a),
hypoxia (Hif1a, Ndrg1) and interferon-gamma response (Ifngr1, Ly6a, Irf7,
Cxcl16); the glycolysis genes the paper nominates downstream of HIF1alpha
(Pgk1, Pkm, Slc16a3); and Il1r1 as the marker that splits AT2 cells into
responder and non-responder subsets.

**One limitation, stated as it bears on reuse.** The single-cell deposit has
**one library per condition**. GSE145031 is six libraries, three time points
(PBS control, day 14, day 28 after bleomycin) by two sorts (lineage-labelled
Tomato-positive and Tomato-negative), with no replicate library at any time
point. GSE144468 is two libraries, control and IL-1beta-treated organoids. So
no transcriptomic contrast in this deposit carries within-group replication,
which is the same ceiling paper 5 turned out to have, and any difference a
reanalysis finds between time points is a difference between two libraries
rather than a tested difference between two groups of animals. The trajectory
claims are a different matter: an ordering inside one library does not need
between-group replication, and that is the part of this paper a reanalysis can
genuinely re-ask.

**One thing to check in the figures rather than assume.** Several
quantifications plot one dot per histological section, with three to six mice
per group, and report mean with standard error. Where the dot is a section, the
plotted unit is not the animal, and this repository's standing rule is that
sections from one animal are technical replicates. This is not a criticism of
the paper's conclusions, which rest on genetics as much as on counts; it is a
note that any figure reused from here must be re-read for what its n refers to
before it is quoted as a per-animal effect.

**One bridge.** This repository has already re-derived the AT2 to Krt8-high
transitional to AT1 ordering from paper 1's influenza series, per-animal and
with replication, and that work is displaced to `archive/DISPLACED.md` as
established elsewhere. This paper gives the same axis a cause (IL-1beta from
interstitial macrophages), a responder subset (Il1r1-positive AT2 cells), a
mechanism (HIF1alpha-driven glycolysis) and a failure mode (chronic IL-1beta
stalls AT1 maturation and DATPs accumulate). The honest comparison is
therefore not one trajectory against another but a cause tested in one injury
model (bleomycin, here) against a time course with replication in another
(influenza, paper 1). Paper 5 then asks what the same state does when the
injury is an oncogene instead.

---

## 1. The cascade as the paper states it

1. Bleomycin injury expands interstitial macrophages at day 7, while alveolar
   macrophages fall, and both return to homeostatic levels by day 28.
2. Interstitial-macrophage-derived IL-1beta acts on the subset of AT2 cells
   that express Il1r1.
3. Those cells lose AT2 identity and become primed AT2 cells.
4. Primed AT2 cells convert into DATPs, a transient state carrying Cldn4,
   Krt8, Ndrg1, Sprr1a and AW112010, with p53 and interferon-gamma programmes
   up, proliferation arrested, and canonical AT1 markers still low.
5. The conversion runs through HIF1alpha-driven glycolysis, and blocking it
   blocks AT1 differentiation.
6. DATPs differentiate into mature AT1 cells, which lineage tracing
   establishes rather than infers.
7. If IL-1beta persists, AT1 maturation does not complete, DATPs accumulate,
   and regeneration is impaired. This is the step that turns a regenerative
   state into a pathological one, and it is the step paper 5 echoes in a
   tumour context.

---

## 2. Which evidence carries which link

A reanalysis of the deposit can only re-ask the questions in the scRNA rows.
Everything else needs the bench, and saying so up front is the point of this
table.

| Link | Evidence in the paper | Reusable from the deposit |
|---|---|---|
| The five states exist and are distinguishable | scRNA-seq of lineage-labelled epithelium, and of organoids | **Yes**, this is the Gate 1 question |
| DATPs sit between primed AT2 and AT1 in pseudotime | PAGA plus diffusion pseudotime on the same data | **Yes**, and it needs no replication, because it is an ordering inside a library |
| DATPs derive from AT2 cells | lineage tracing with an AT2 reporter, about 10% of labelled cells Krt8-positive at day 14 | No. Transcriptomes cannot establish origin |
| DATPs become AT1 cells | lineage tracing from an Ndrg1 driver | No |
| IL-1beta is sufficient to make DATPs | organoids with recombinant IL-1beta, imaging and flow | Partly: the deposit has the organoid transcriptomes, so the state change is re-derivable, the sufficiency claim is not |
| Only Il1r1-positive AT2 cells respond | Il1r1 reporter mice, functional and epigenetic comparison | No, and the ATAC-seq deposit cannot help (section 4) |
| HIF1alpha glycolysis is required | extracellular acidification, glucose uptake, digoxin, conditional deletion | No. The glycolysis gene expression is re-derivable; the requirement is not |
| Interstitial macrophages are the IL-1beta source | flow cytometry across the injury course | No. The deposit sorts epithelium, not macrophages |
| Chronic IL-1beta stalls AT1 maturation | prolonged exposure in vivo and in organoids | No |

---

## 3. The deposited data, and what it can carry

Four accessions are linked to this paper. GSE144553 is the SuperSeries.

| Accession | Content | Libraries | What it can support |
|---|---|---|---|
| GSE145031 | scRNA-seq of AT2 lineage-traced epithelium | 6: PBS, day 14, day 28, each split Tomato-positive and Tomato-negative | State recovery and trajectory ordering; no tested between-time-point contrast |
| GSE144468 | scRNA-seq of AT2 organoids | 2: control and IL-1beta | The state change under IL-1beta as a direction, not a test |
| GSE144598 | ATAC-seq of AT2 subsets | 2 files deposited as bigwig coverage only | **Almost nothing.** Coverage tracks carry no peaks and no reads, so the epigenetic claim cannot be re-derived from the supplementary files. Raw reads would have to come from SRA |
| GSE144553 | SuperSeries over the three above | 12 samples | Bookkeeping only |

The Tomato-negative libraries are worth noting, because they are an unusual and
useful control: the same lung, the same dissociation, and the cells the AT2
lineage did **not** label. Any claim that a state is AT2-derived can be checked
for whether the state also appears unlabelled, which is a within-animal
comparison rather than a between-group one.

## 4. Their pipeline, as stated in Methods

Cell Ranger 2.0.2 for alignment and counting. Seurat v2.0 for clustering and
markers, with doublet handling described as removing cells judged unlikely by
prior knowledge and by higher UMI count rather than by a dedicated caller. The
processed Seurat object was then converted for trajectory work in Scanpy 1.3.6,
with k-nearest neighbours recomputed at k = 15, PAGA for topology and the
diffusion pseudotime function with default parameters. MACS2 for ATAC peak
calling. Imaging quantification in ImageJ on Leica confocal sections, at least
five sections and ten alveolar regions from three mice per group.

Two differences from this repository's standing pipeline matter for a
reanalysis. There is no ambient-RNA correction and no dedicated doublet
detector, and the quality thresholds are not the per-library median absolute
deviation rule used here. Both are ordinary for 2020 and neither is a fault;
they are the reasons a blind re-run will not reproduce cluster boundaries
exactly, which is the same bookkeeping list section 1 of the paper 5 divergence
document sets out.

## 5. Status of the claims in this note

Everything above is a summary of the paper, not a result of this repository. No
trial has run yet. The Gate 0 data reality check is
[`trials/`](trials/), and Gate 1 will ask whether the five states and the
primed-AT2 to DATP to AT1 ordering are recoverable from the deposited counts
with the paper's labels held out of the fitting.
