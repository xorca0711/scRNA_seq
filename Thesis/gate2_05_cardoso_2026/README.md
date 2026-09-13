# Study note: Cardoso, Lee et al. 2026, Early fibrotic niches establish tumour-permissive microenvironments

Roadmap position: **Gate 2, paper 5** in the ordered paper roadmap of the
[UC Berkeley SAP PI Target Map](https://app.notion.com/p/3d1151616b44814c8697ff3af2f8f831),
entered ahead of papers 2, 3 and 6 on the owner's instruction. The paper is a
Choi-lab extension of the Gate 1 spine: the same damage-associated transient
progenitor (DATP) state that mediates alveolar regeneration is co-opted by
KrasG12D-mutant AT2 cells, and the question becomes what that state does to
the cells around it.

**Citation.** Cardoso EC, Lee H, England FJ, Cho H, Lu R, Varankar SS, Park MS,
Rekhtman N, Koo B-K, Simons BD, Choi J, Lee J-H. *Early fibrotic niches
establish tumour-permissive microenvironments.* **Nature**
2026;653(8113):254-264. DOI
[10.1038/s41586-026-10399-6](https://doi.org/10.1038/s41586-026-10399-6),
PMID [42020743](https://pubmed.ncbi.nlm.nih.gov/42020743/), PMC
[PMC13149335](https://pmc.ncbi.nlm.nih.gov/articles/PMC13149335/). Open
access, CC BY 4.0. Received 16 March 2025, accepted 11 March 2026.

**Companion paper the analysis depends on.** England FJ, Bordeu I, Ng M-E,
Bang J, Kim B, Choi J, Cardoso EC, Koo B-K, Simons BD, Lee J-H. *Sustained
NF-kB activation allows mutant alveolar stem cells to co-opt a regeneration
program for tumor initiation.* **Cell Stem Cell** 2025;32(3):375-390.e9. DOI
[10.1016/j.stem.2025.01.011](https://doi.org/10.1016/j.stem.2025.01.011),
PMID [39978341](https://pubmed.ncbi.nlm.nih.gov/39978341/). This is reference 7
of the Nature paper and the source of the RFP+ mutant epithelial cells used in
its cell-communication analysis; its data are GEO
[GSE247505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247505),
which the Nature paper's own data-availability statement does not name.

> Provenance: extracted from the article text and the Extended Data captions
> (PyMuPDF text layer of the owner's local PDF, cross-checked against the
> Europe PMC full text) and from the deposited GEO SOFT family files, on
> 2026-09-12 by an AI session; not yet reviewed by the owner. Numbers in
> sections 1 to 5 are quoted from the paper, not recomputed. What this
> repository measured itself is in
> [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md) and under
> [`trials/`](trials/), and is kept separate from the quotations here.
> According to PubMed, the bibliographic metadata above.

---

## 0. Reading workflow (the five questions the roadmap requires)

**Question tested.** During the window between oncogenic activation and
visible tumour, does the mutant epithelial cell actively build the
microenvironment that permits its own expansion, and if so, through which
signal, in which order, and is the resulting circuit reversible?

**Evidence type.** Mixed, and the mixture is the point: no single claim in the
cascade rests on transcriptomes alone. Single-cell RNA sequencing supplies the
states and the candidate ligands; confocal immunofluorescence supplies the
spatial arrangement and the entire time course; organoid co-culture and
tri-culture supply sufficiency; and four in vivo perturbations (KrasG12D
inhibition, EGFR inhibition, AT2-specific *Areg* deletion, diphtheria-toxin
ablation of Pdgfra+ fibroblasts) supply necessity. Section 2 below classifies
every link of the cascade by the evidence that carries it, because the
deposited transcriptomes can only be used to re-ask the questions in the
"scRNA" rows.

**Reusable variables.** The fibrotic-fibroblast marker set (Pdgfrb, Runx1,
Runx2, Tnc, Fst, Acta2, with reduced Pdgfra); the inflammatory-fibroblast set
(Lcn2, Saa3, Sfrp1, Cxcl12) and the paper's claim that it is *distinct* from
the fibrotic set and lacks Tnc; the reprogrammed alveolar-macrophage set
(Msr1, Cdh1, Ch25h, Ear6, Fstl1, with MHC-II low and Cxcl2/Cxcl16 high); the
DATP-like epithelial set (Sox9, Cldn4, Itga2, Ndrg1, Krt8, with Areg high);
the mutant-state panel of Fig. 4l (AT2, Cd177+, DATP-like, cycling,
AT1-like); the two CellChat legs and their reported top hits (epithelium to
fibroblast: Areg-Egfr; fibroblast to epithelium: Wnt5a, Igf1, Spp1); and the
paper's own statistical convention, that the animal is the unit and
within-mouse fields are averaged to a mouse-level value before testing.

**One limitation, stated as it bears on reuse.** Every mouse single-cell
library in the deposit is a *pool of three mice* and a *single time point*
(two weeks after induction, or four days for part of the companion series).
There is therefore no biological replication inside any genotype contrast in
the deposited data, and no transcriptomic time course of the mesenchymal or
immune compartments at all. The paper's temporal claim, that fibroblast
reprogramming precedes macrophage expansion, is carried entirely by
immunofluorescence at 1, 2, 4 and 8 weeks (Fig. 2a-e, n = 3 to 6 mice per
time point). A reanalysis of the deposit cannot test that ordering, and any
genotype difference it finds is a difference between two libraries, not a
tested difference between two groups of animals.

**One bridge.** The repository's existing mouse series (GSE262927, influenza
injury) and this one share the AT2 to Krt8-high transitional axis and a
macrophage compartment that is remodelled after the epithelium changes. The
honest comparison is not "tumour versus infection" but a narrower one: the
non-oncogenic injury course has a per-animal time course with replication and
no oncogene; this deposit has an oncogene, a fibrotic fibroblast state and a
reprogrammed macrophage state but neither replication nor time. Each covers
the other's blind spot, which is the argument for reading them together
rather than merging them.

---

## 1. The cascade as the paper states it

```
KrasG12D activation in AT2 cells (Sftpc-CreERT2; Red2Kras)
  |
  +-> mutant AT2 adopt a DATP-like regenerative state (Sox9+, Krt8+, Areg-high)
  |        |
  |        +-- Areg --> EGFR on Pdgfra+ alveolar fibroblasts
  |                       |
  |                       +-> fibrotic ("reprogrammed") fibroblasts
  |                             Pdgfrb+ Runx1+ Runx2+ Tnc+ Fst+ Acta2+, Pdgfra-low
  |                             |
  |                             +-- Tnc --> TLR4 on alveolar macrophages
  |                             |             |
  |                             |             +-> reprogrammed AMs: Msr1+, Cdh1+,
  |                             |                 Ch25h+, MHC-II-low, Cxcl2/Cxcl16 high
  |                             |                   |
  |                             |                   +-- Cxcl2/Cxcr2 --> neutrophils
  |                             |                   +-- Cxcl16/Cxcr6 --> gamma-delta T cells
  |                             |                   +-- IL-1beta --> Lcn2+ inflammatory
  |                             |                        fibroblasts at the tumour periphery
  |                             |
  |                             +-- Wnt5a, Igf1, Spp1 --> back to the mutant epithelium,
  |                                   sustaining the DATP-like state
  |
  +-> immunosuppressive shift: Treg and PD-1+ T cells up, CD8 T cells exhausted
```

Spatial grammar the paper insists on: fibrotic fibroblasts occupy tumour
cores and borders; Lcn2+ inflammatory fibroblasts appear only from 4 weeks and
only at the periphery; macrophages accumulate in *inter*-tumour regions while
neutrophils and gamma-delta T cells infiltrate tumour cores.

---

## 2. Which evidence carries which link

This is the table that governs reuse. Only the rows marked **scRNA** can be
re-asked from the deposit; the rest can be quoted but not re-tested here.

| Link | Evidence class | Where | Replication as reported |
|---|---|---|---|
| A reprogrammed fibroblast cluster exists only in Red2Kras lungs | **scRNA** | Fig. 1b-d, ED Fig. 1b-d | 1 library per genotype, 3 mice pooled each |
| Those fibroblasts are Pdgfrb+Acta2+Runx1+ and sit next to RFP+ mutant cells | imaging | Fig. 1e, ED Fig. 1h | n = 3 mice |
| They arise from alveolar fibroblasts | trajectory inference (Monocle 3) | ED Fig. 1e, o | same 1 library |
| ... and lineage-labelled Pdgfra+ cells do upregulate Pdgfrb around tumours | lineage tracing plus ablation | ED Fig. 8a-d | n = 2 mice |
| The same fibroblast state is shared with bleomycin injury | **scRNA integration** with ref. 24 | ED Fig. 1i-o | external dataset |
| A reprogrammed AM population exists only in Red2Kras lungs | **scRNA** | Fig. 1f-h | 1 library per genotype, 3 mice pooled each |
| AMs expand, lose MHC-II, gain Msr1 | flow cytometry, imaging | Fig. 1i-m | n = 3 per group |
| Those AMs are resident, not monocyte-derived | lineage tracing (CCR2-CreERT2;ZsGreen) plus engraftment | ED Fig. 3a-b | n = 3 mice |
| AM depletion shrinks tumours and blocks neutrophil and gamma-delta recruitment | perturbation (clodronate) | Fig. 1n-p, ED Fig. 3c-l | n = 4 vs 3 mice |
| **Fibroblast reprogramming precedes macrophage remodelling** | **imaging time course only** | Fig. 2a-e | n = 3 to 6 mice per time point |
| Lcn2+ inflammatory fibroblasts are late, peripheral and Tnc-negative | **scRNA** for the state, imaging for the timing and position | ED Fig. 4l-o, Fig. 2g-h | 1 library; imaging n = 3 |
| DATP-like cells, not AT2 cells, are the signalling source to fibroblasts | **scRNA (CellChat)**, differential-source comparison | ED Fig. 5c-h | integration of 2 series |
| Areg is the top ligand of that leg, and is DATP-specific | **scRNA (CellChat)** plus imaging | ED Fig. 5e-h, Fig. 3a-b | imaging n = 2 mice |
| Areg induces a fibrotic programme in fibroblasts and only in fibroblasts | organoid and culture perturbation | ED Fig. 6c-n | n = 3 to 6 experiments |
| EGFR inhibition blocks fibroblast reprogramming and AM rewiring | perturbation (gefitinib) | Fig. 3c-g, ED Fig. 8e-r | n = 2 to 4 per arm |
| Mutant-cell-derived Areg is required in vivo | genetic perturbation (Areg-flox) plus **scRNA** | Fig. 4a-m, ED Figs. 9, 10 | 15 vs 15 mice for burden; 2 libraries per genotype for scRNA |
| Tnc from fibrotic fibroblasts reprogrammes AMs through TLR4 | **scRNA** for co-expression, in vitro for causation | ED Fig. 4d-k | n = 3 experiments; no in vivo Tnc perturbation |
| Fibroblasts feed back on the epithelium through Wnt5a, Igf1, Spp1 | **scRNA (CellChat)** plus gefitinib organoids | ED Fig. 5i-j, ED Fig. 7o-v | integration; n = 3 to 4 experiments |
| The circuit is reversible by KrasG12D inhibition | perturbation (MRTX1133), imaging readout | Fig. 3h-p | n = 3 vs 3 mice |
| The programme is conserved in human LUAD | **scRNA reanalysis** of ref. 45 plus imaging | ED Fig. 11 | 5 matched pairs, EGFR wild-type only |
| Human AT2 cells with KRASG12D reproduce it | **scRNA** of organoids plus co-culture | Fig. 5c-i | 2 libraries; n = 3 experiments |
| EGFR-L858R engages the same circuit | engraftment plus imaging | ED Fig. 12 | n = 2 mice |

---

## 3. What the paper's CellChat analysis actually did

Worth stating precisely, because it is the part most easily misread as a
single global run. Extended Data Fig. 5 integrates two datasets: lineage-
labelled RFP+ mutant epithelial cells from the companion paper (England 2025,
GSE247505) and the fibroblasts of this paper's Red2Kras and Confetti lungs
(GSE316241). On that joint object they ran, and reported separately:

1. **an aggregated network** over epithelial and fibroblast states (panel c);
2. **a differential-source comparison**: the same target (alveolar
   fibroblasts) with AT2 cells versus DATP-like cells as the source (panel d),
   which is the comparison that licenses "DATP-like cells, not AT2 cells, are
   the hub";
3. **the ranked outgoing signals** of the DATP-like to alveolar-fibroblast leg
   (panel e), from which the EGF pathway and then the Areg-Egfr pair are taken
   as top by communication probability (panel f);
4. **expression backing** for that pair (panels g, h); and
5. **the reverse leg**, fibroblasts to epithelium, whose top hits are Wnt5a,
   Igf1 and Spp1 (panels i, j).

The paper does not report what ranks second to fifth on either leg, does not
repeat the analysis at more than one time point (it cannot; see section 4),
and does not test the ranking's sensitivity to how the DATP-like cluster is
defined. Those three gaps are what make the owner's Gate 2(a) worth running
rather than re-running.

**The standing caveat applies unchanged.** CellChat infers communication from
ligand-receptor co-expression in dissociated cells and has no proximity
information. The paper supplies the missing proximity separately, by imaging;
a rerun here would not.

---

## 4. The deposited data, and what it can carry

Accessions from the paper's data-availability statement, plus the companion
series, all verified public on 2026-09-12 (measurements in
[`trials/c0_data_reality_check/`](trials/c0_data_reality_check/)):

| Accession | Content | Libraries | Design as deposited |
|---|---|---|---|
| [GSE316241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316241) | mesenchyme (CD45-CD31-EpCAM-), Confetti and Red2Kras | 2 | 2 weeks post-induction; 3 mice pooled per library |
| [GSE316243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316243) | niche: 1:1 immune (CD45+) and stromal (CD45-EpCAM-) | 2 | 2 weeks; 3 mice pooled per library |
| [GSE316244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316244) | Areg-flox arm: niche and RFP+ epithelium, het and hom | 4 | 2 weeks; 3 mice pooled per library |
| [GSE310335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310335) | human alveolar organoids, control and KRASG12D | 2 | day 7 after doxycycline induction |
| [GSE247505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247505) | England 2025 RFP+ and YFP+ lineage-labelled epithelium | 20 | 4 days, 2 weeks, 12 weeks; two replicates per arm |

Three consequences follow, and they set the gates:

- **No genotype contrast in the mouse deposit has biological replication.**
  Each library pools three mice, and each genotype contributes one library per
  sort. A difference between Confetti and Red2Kras, or between Areg-flox/+ and
  Areg-flox/flox, is a difference between two libraries. It can be described;
  it cannot be tested, and no P value is admissible. The paper's own tested
  claims use imaging and flow cytometry, where n is mice, which is consistent.
- **The mesenchymal and immune compartments have no time course.** Every
  deposited library of those compartments is two weeks post-induction. The
  ordering claim (fibroblast before macrophage) is therefore untestable from
  the deposit, which closes the owner's Gate 2(d) as a transcriptomic
  question. Only the epithelium has time points, and only in the companion
  series.
- **The libraries do not share one gene space.** The deposits were aligned
  with four different CellRanger versions against GRCm38, so the reference
  gene lists differ; the companion series and this paper's niche series agree,
  the mesenchyme series does not, and any integration (the CellChat object
  included) needs an explicit intersection on Ensembl gene ID. The exact
  numbers are in the C0 tables.

---

## 5. Their pipeline, as stated in Methods

| Step | What they used |
|---|---|
| Alignment | CellRanger, GRCm38 (mouse) or GRCh38 GENCODE v38 (human); version differs per experiment |
| QC | most libraries: < 10% mitochondrial, > 1,000 genes, > 2,000 and < 50,000 UMIs; the Confetti-versus-Red2Kras immune library and the human organoids instead: 500 to 7,000 genes, > 2,000 UMIs |
| Framework | Seurat, with Scanpy 1.9.1 used for part of the work |
| Normalisation and features | log-normalise, scale, 2,000 variable genes, 30 principal components |
| Clustering | Louvain, UMAP for display |
| Integration | Harmony across Seurat objects; for the injury comparison, `FindIntegrationAnchors` with **10,000 anchor features** then `IntegrateData` |
| Trajectories | Monocle 3 |
| Communication | CellChat, standard protocol |
| Enrichment | g:Profiler on the top differentially expressed genes, P < 0.05, biological process |
| Doublets | not stated beyond the UMI ceiling; no dedicated caller named |
| Ambient RNA | not mentioned |
| Statistics | animal or independent experiment as the unit; within-mouse fields averaged to a mouse-level value; two-tailed unpaired t-test or Mann-Whitney; nested analyses stated to agree with the averaged ones |

Two divergences from this repository's practice are worth naming now, because
they will show up in any comparison. They used no dedicated doublet caller and
no ambient-RNA correction, and they clustered with Louvain on Harmony-corrected
objects. This repository runs Scrublet per capture and decides batch correction
per dataset from measured replicate mixing, which is not available here
because there are no replicates to measure mixing between.

---

## 6. Status of the claims in this note

| Claim | Status here |
|---|---|
| Everything in sections 1 to 5 as a description of what the paper did and reported | Descriptive only (transcribed from the paper; owner review pending) |
| The evidence classification in section 2 | Descriptive only (an AI session's reading of the figure captions; owner review pending) |
| "The deposited mouse libraries carry no within-genotype biological replication" | see [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md), trial C0 |
| "The deposit cannot test the fibroblast-before-macrophage ordering" | see trial C0 |
| Any statement about what this repository's own clustering finds | Not established until the Gate 1 trial has run |
