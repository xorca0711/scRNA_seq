# Study note: Choi, Lee et al. 2020, Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors that mediate alveolar regeneration

Roadmap paper 2, Gate 1. Written on 2026-09-15 at the owner's direction after
the owner read the paper (DEVELOPMENT decision 23). It supersedes an
AI-written note of 2026-09-13 that the owner withdrew (decision 21); that
note stays in git history and nothing in it is retained without re-reading.
The paper's marker sets were read from the Europe PMC XML and checked against
the PDF; the PDF renders one gene symbol ambiguously, and the extract says so.

**Choi J, Park J-E, Tsagkogeorga G, Yanagita M, Koo B-K, Han N, Lee J-H.**
*Inflammatory signals induce AT2 cell-derived damage-associated transient
progenitors that mediate alveolar regeneration.*
**Cell Stem Cell** 2020;27(3):366-382.e7.
DOI: [10.1016/j.stem.2020.06.020](https://doi.org/10.1016/j.stem.2020.06.020) ·
PMID [32750316](https://pubmed.ncbi.nlm.nih.gov/32750316/) ·
PMC [PMC7487779](https://pmc.ncbi.nlm.nih.gov/articles/PMC7487779/)

**Why this paper is in the roadmap.** It defines the DATP state that roadmap
paper 5 (Cardoso 2026) carries into tumour initiation, and its first author
leads the owner's own laboratory. Paper 1 (Niethamer 2025) supplies the
Krt8-high transitional axis this repository has already re-derived; this
paper gives that axis a cause (IL-1beta), a responder subset (Il1r1-positive
AT2 cells), a mechanism (HIF1alpha-driven glycolysis) and a failure mode
(chronic IL-1beta stalls AT1 maturation).

---

## 0. The five questions

**Question tested.** During alveolar regeneration after bleomycin, what is
the path from an AT2 cell to a mature AT1 cell, which niche signal moves a
cell along it, and what happens to regeneration when that signal does not
switch off?

**Evidence type.** Mixed, and weighted toward genetics and imaging. scRNA-seq
of lineage-labelled epithelium at PBS, day 14 and day 28 (12,086 cells after
the paper's contaminant removal; Figure 1) and of AT2 organoids with and
without IL-1beta (Figure 2F to 2H) supplies the states, their composition and
their ordering. Lineage tracing from SPC, Ndrg1 and Krt8 drivers supplies
origin and fate (Figures 3, S4). Organoid co-culture with macrophages,
recombinant IL-1beta, digoxin, 2-deoxyglucose and conditional deletion of
Il1r1 and Hif1a supply sufficiency and necessity (Figures 2, 4, S5, S6). An
Il1r1 reporter and ATAC-seq supply the responder subset (Figures 5, 6). Human
IPF and adenocarcinoma sections supply the disease bridge (Figure 7I to 7K,
S7). Section 2 sorts each link by what a reanalysis of the deposit can re-ask.

**Reusable variables.** Five states (hAT2, cAT2, pAT2, DATP, AT1) with the
paper's own markers: canonical AT2 (Sftpc, Sftpa1, Lyz2); AT2 identity genes
that fall in primed cells (Etv5, Abca3, Cebpa) and lipid genes that fall with
them (Acly, Hmgcr, Hmgcs1); an inflammatory set that rises in primed cells
(Ptges, Orm1, Tmem173, Ifitm2, Ifitm3, and a lipocalin printed as "Lcn1");
cycling genes (Cdk1, Mki67, Cenpa); the DATP set (Cldn4, Krt8, Ndrg1, Sprr1a,
AW112010) with its negative condition (low Pdpn, Hopx, Cav1); early AT1
(Lmo7, Pdpn, Hopx) and late AT1 (Aqp5, Vegfa, Cav1, Spock2); the DATP
programmes from gene ontology (p53: Trp53, Mdm2, Ccnd1, Gdf15; arrest: Cdkn1a,
Cdkn2a; hypoxia: Hif1a, Ndrg1; interferon-gamma: Ifngr1, Ly6a, Irf7, Cxcl16);
glycolysis downstream of HIF1alpha (Pgk1, Pkm, Slc16a3); Il1r1 as the
responder marker. Time points PBS, day 14, day 28. The replicate is two
pooled mice per library, one library per condition.

**One limitation, stated as it bears on reuse.** One library per condition.
Every difference between time points, sorts or treatments in this deposit is
a difference between two libraries, each pooling two mice, and cannot be
tested. The paper's tested claims rest on imaging and flow cytometry, where n
is mice or sections; the transcriptome is the map, not the evidence. Where a
figure plots one dot per section, the plotted unit is not the animal.

**One thing to check in the figures rather than assume.** Figure 1C's
composition per time point is read from the deposit in trial D2, but the
paper gives no per-cluster numbers in its text beyond "approximately 6%
cycling", so the readings are qualitative directions, not numeric targets.

**One bridge.** This repository has already re-derived the AT2 to Krt8-high
transitional to AT1 ordering from paper 1's influenza series, per animal and
with replication (displaced to `archive/`). This paper gives the same axis
its inducer and its failure mode in a second injury model. The honest
comparison is a cause tested in one model against a time course with
replication in another, and the DATP definition here is what paper 5's
"regenerative-like" state leans on.

---

## 1. The cascade as the paper states it

1. Bleomycin injury expands interstitial macrophages at day 7 while alveolar
   macrophages fall; both return to homeostatic levels by day 28 (Figure S2).
2. Interstitial-macrophage-derived IL-1beta acts on the subset of AT2 cells
   that express Il1r1 (Figures 2, 5).
3. Those cells lose AT2 identity and become primed AT2 cells; the priming
   signature rises as cycling AT2 cells pass from S to G2/M, where Il1r1 is
   also up (Figure S5G).
4. Primed AT2 cells convert into DATPs, which carry Cldn4, Krt8, Ndrg1, Sprr1a
   and AW112010, with p53, arrest, hypoxia and interferon-gamma programmes up
   and canonical AT1 markers still low (Figure 1D, S1D).
5. The conversion runs through HIF1alpha-driven glycolysis; digoxin, 2-DG and
   Hif1a deletion block it and block AT1 differentiation (Figure 4, S6).
6. DATPs differentiate into mature AT1 cells, and some revert to AT2, by
   lineage tracing from Ndrg1 and Krt8 (Figure 3, S4).
7. If IL-1beta persists, AT1 maturation does not complete: early AT1 markers
   are reached, late markers are not, DATP genes stay up; withdrawal or 2-DG
   rescues it (Figure 7A to 7H). KRT8-positive CLDN4-positive cells are seen
   in IPF and adenocarcinoma tissue, and the IPF KRT17-positive population
   shares the DATP signature (Figure 7I to 7K, S7H).

---

## 2. Which evidence carries which link, and what the deposit can re-ask

| Link | Evidence in the paper | Re-askable from the deposit | Trial |
|---|---|---|---|
| The five states exist and are distinguishable | scRNA-seq, in vivo and organoid | Yes | D2, D2b, D5, D5b: four of five; primed AT2 not as a cluster |
| DATPs sit between primed AT2 and AT1 in pseudotime | PAGA plus diffusion pseudotime | Yes, and needs no replication (an ordering inside a library) | D3 |
| Composition changes across PBS, day 14, day 28 | scRNA-seq, one library each | Direction only | D2b, D4: the paper's directions for hAT2, cAT2, DATP and AT1 |
| DATPs carry p53, arrest, hypoxia, interferon-gamma and glycolysis programmes | gene ontology and heatmaps | Yes | D6 |
| IL-1beta raises primed and DATP fractions and stalls late AT1 markers in organoids | organoid scRNA-seq, one library per treatment | Direction only | D5b: the DATP-marker fraction rises, the late-marker set sits on its threshold |
| DATPs derive from AT2 cells | SPC lineage tracing | No | none |
| DATPs become AT1 cells and can revert to AT2 | Ndrg1 and Krt8 lineage tracing | No | none |
| IL-1beta is sufficient and required | organoids, Il1r1 deletion, imaging, flow | No; only the state change is re-derivable | D5 |
| Only Il1r1-positive AT2 cells respond | Il1r1 reporter, ATAC-seq | No; the ATAC deposit is coverage only | none |
| HIF1alpha glycolysis is required | ECAR, glucose uptake, digoxin, deletion | No; the gene expression is re-derivable | D6 |
| Interstitial macrophages are the IL-1beta source | flow cytometry; the non-lineage libraries | Partly, from the Tomato-negative libraries; not attempted | proposed |
| Il1r1 rises in G2/M cycling cells | scRNA-seq with phase assignment | Not attempted: phase assignment needs an external gene list | none |

---

## 3. The deposited data, and what it can carry

| Accession | Content | Libraries | What it can support |
|---|---|---|---|
| GSE145031 | scRNA-seq of AT2 lineage-traced epithelium | 6: PBS, day 14, day 28, each split Tomato-positive and Tomato-negative; two mice pooled per library | State recovery, ordering, direction of composition; the Tomato-negative libraries as a within-experiment control |
| GSE144468 | scRNA-seq of AT2 organoids | 2: control and IL-1beta, deposited already filtered | The state change under IL-1beta as a direction |
| GSE144598 | ATAC-seq of AT2 subsets | 2 files, bigwig coverage only | Almost nothing without SRA |
| GSE144553 | SuperSeries | 12 | Bookkeeping |

## 4. Their pipeline, as stated in Methods

Cell Ranger 2.0.2 on GRCm38; cells kept with more than 500 and fewer than
7,000 genes and more than 2,000 UMI; scanpy 1.3.6 defaults for normalisation,
HVG, PCA, neighbours and Louvain; contaminants (214 ciliated, 16 mesenchyme,
25 immune of 12,514 captured) removed after a first unsupervised round; the
EpCAM-negative stromal cluster removed from organoids; Scrublet per sample
(sim_doublet_ratio 2, n_neighbors 30, expected rate 0.1, score above 0.7)
plus over-clustering at resolution 20 with a 0.6 cluster mean, then
prior-knowledge checks; k = 15 neighbours, PAGA and default DPT for the
trajectory, cycling cells excluded. This repository differs in three ways it
records rather than hides: cell calling is the paper's filter on the raw
whitelist instead of Cell Ranger's caller (trial D1); Scrublet uses the 10x
expected rate with scanpy's threshold and a logged fallback; clustering is
Leiden at three pre-declared resolutions.

## 5. Status of the claims in this note

Everything above is a summary of the paper. The trials and their outcomes are
in [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md) and indexed in
[`trials/README.md`](trials/README.md); the claims they produce are rows C85
to C104 of the root register, all awaiting the owner's retain or reject.

What came back on 2026-09-15, in one paragraph. Four of the five states
(hAT2, cycling AT2, DATP, AT1), the DATP time course, the hAT2 to DATP to AT1
ordering inside one library, DATP's p53, arrest, hypoxia, interferon-gamma
and glycolysis programmes in vivo, and IL-1beta's shift of the organoid
epithelium all reproduce from the deposit as descriptions. Primed AT2 does
not separate as a cluster anywhere, and the organoid cluster the paper calls
77 percent primed keeps its AT2 identity genes. DATP is not a doublet
cluster. Two readings are not computable (the within-AT1 Figure 7 comparison
in the first organoid pass; the chain through primed AT2), one is saturated
(Krt8 and Cldn4 co-detection), and the first annotation rule missed AT1 and
the organoid stromal cluster through a disclosed defect, corrected beside it
without moving a threshold. Nothing is tested: one library per condition.
