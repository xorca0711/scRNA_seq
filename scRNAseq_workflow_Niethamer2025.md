# scRNA-seq Workflow — *Longitudinal single-cell profiles of lung regeneration after viral infection reveal persistent injury-associated cell states*

**Target paper.** Niethamer TK, Planer JD, Morley MP, … Vaughan AE, Morrisey EE.
*Cell Stem Cell* 2025;32(2):302–321.e6. DOI: [10.1016/j.stem.2024.12.002](https://doi.org/10.1016/j.stem.2024.12.002)
PMID 39818203 · PMC11805657 · Preprint: bioRxiv 2024.05.24.595801

**Data.** GEO **GSE262927**. Interactive browser and code links are given in the paper as shortened URLs (bit.ly/49kMFJa and bit.ly/3vqeXUZ respectively); resolve them from the published Data & Code Availability section.

> ### ⚠ Filing correction
> The paper is **already in this folder, under the wrong name.**
> `Lv-2024-Alveolar-regeneration-by-airway-sec/Lv-2024-Alveolar-regeneration-by-airway-sec+Suppl.pdf` is not a Lv 2024 paper — it is **Niethamer et al. 2025**, main text plus supplement, 42 pages. The three spreadsheets beside it are that paper's supplementary tables, not Lv's:
>
> | File | Is actually | Contents |
> |---|---|---|
> | `mmc2.xlsx` | **Table S1** | 29 adult human lung donors — age, sex, race/ethnicity, disease (AAT, BOS, COPD, healthy, IPF, LAM) |
> | `mmc3.xlsx` | **Table S2** | 29 pediatric/developmental human samples — age in days, epoch (zero-day → adult) |
> | `mmc4.xlsx` | **Table S3** | **33 mouse scRNA-seq samples** — GEO sample names, tamoxifen dpi, infection date, euthanasia dpi, age, genotype, sex, max weight loss |
>
> Rename the folder and file. If you actually want Lv et al. on airway secretory cells, it is not here and needs downloading. **Table S3 (mmc4.xlsx) is the sample manifest for everything below** — keep it open while you work.

**Why this document exists.** The five method PDFs sitting in this folder — SoupX, Scrublet, scds, Slingshot, tradeSeq — are not a random collection. They are precisely the non-Seurat components of this paper's pipeline. This document places each one in its correct slot, states the parameters the paper actually used, and flags the decisions you have to make yourself because the paper leaves them implicit.

---

## 1. The experimental design you are analysing

You cannot organise the computation without holding the design in your head, because almost every analysis step in this paper is *conditioned on time point and on lineage-trace status*.

**Injury model.** PR8-GP33 H1N1 influenza A (genetically targeted mice) or A/Puerto Rico/8/34 H1N1 (C57BL/6J), intranasal, 50 µL sterile saline. Dose titrated to 1 LD50 for the targeted lines; 1:500,000 dilution for B6 (15–20% weight loss, no death). Injury is *spatially heterogeneous* — zones of severe damage, damage and normal architecture coexist in the same lobe. This heterogeneity is the reason the paper can find a cell state that localises to damaged regions.

**Inclusion criteria (matters for batch interpretation).** Experiments in week 1 post-infection excluded mice losing <5% body weight; experiments after week 1 excluded mice losing <15%. So "6 dpi" and "42 dpi" samples are not drawn from identically-filtered animal populations.

**The longitudinal atlas.** Mki67<sup>ires-CreERT2</sup>;Rosa26<sup>LSL-tdTomato</sup> mice, tamoxifen-pulsed in five staggered cohorts (tamoxifen at 2+3, 7+8, 14+15 or 21+22 dpi; 2 × 200 mg/kg by oral gavage), then sequenced. **25 libraries, one per mouse, 123,189 cells merged:**

| Phase | dpi | n libraries |
|---|---|---|
| Baseline | uninjured | 2 |
| Active repair | 6 · 11 · 19 · 25 | 2 each (8) |
| Injury resolution | 42 | 8 (2 per tamoxifen cohort) |
| Injury resolution | 90 | 4 (1 per cohort) |
| Long-term homeostasis | 366 (1 year) | 3 |

The design is **unbalanced by construction** — 42 dpi carries four times the library count of 6 dpi, because four tamoxifen cohorts converge on that harvest date. Any per-time-point comparison of cell proportions must account for this, and any pseudobulk-style aggregation will be dominated by 42 dpi unless you weight or downsample.

Mice were 6–21 weeks old at first treatment, both sexes, on a **mixed (non-backcrossed) genetic background**, and the atlas spans **three independent influenza infection rounds** (Oct 2021, Jul 2022, Dec 2022). Infection round is therefore a real batch variable, it is recorded in Table S3, and the paper does not model it. If you re-analyse this data, that column is the first thing to test as a covariate.

**Non-infectious control injury.** Hyperoxia (95% O₂, 72 h, then 7 days normoxia) was used to show the iCAP state is not influenza-specific — imaging only, not sequenced.

**Composition control — read this before you touch the QC thresholds.** Cells were MACS-fractionated into CD45<sup>+</sup> and CD45<sup>−</sup> and then *recombined at a fixed 10–15% CD45<sup>+</sup> : 85–90% CD45<sup>−</sup> ratio* before loading. This is deliberate: it prevents the massive post-influenza immune infiltrate from swamping the endothelial and epithelial compartments. Consequence for you: **immune cell proportions in this dataset are engineered, not biological.** Any statement about "% immune cells over time" must come from the flow cytometry, not from cluster proportions. Non-immune *relative* proportions remain interpretable.

**Second, independent dataset.** Lineage-traced mice given tamoxifen **before** infection and all sequenced at **19 dpi**, to assign the iCAP state a cell of origin:

| Line | Traces | Tamoxifen | n | Efficiency / specificity |
|---|---|---|---|---|
| Kit<sup>MerCreMer</sup> | CAP1 + bone-marrow-derived | 3 × 200 mg/kg, −16 d | 3 | — |
| Car4<sup>CreERT2</sup> *(new line)* | CAP2, mature aMAC | 3 × 200 mg/kg, −15 d | 2 | 50.6% / 74.6% |
| Ednrb<sup>CreERT2</sup> *(new line)* | CAP2 | 1 × 50 mg/kg, −14 d | 3 | 50.1% / 91.2% |

Critically, **this cohort was processed and clustered separately from the Mki67 atlas** — they were not merged for primary clustering. Comparison is done by matching cluster identities, not by integration. The differing tamoxifen doses across lines are deliberate (they tune labelling efficiency per allele) but mean labelling rates are not directly comparable between lines without the efficiency/specificity figures above.

**Chemistry.** 10x Chromium, **Next GEM Single Cell 3′ v3.1**. Mouse libraries sequenced by Azenta on a NovaSeq; human on NovaSeq 6000. Cell suspensions loaded at 800–1,200 cells/µL. Digestion: collagenase-I 480 U/mL + dispase 100 µL/mL + DNase 2 µL/mL, 35 min at 37 °C, 100 µm then 40 µm filters, ACK lysis. **Reads per cell is not reported** — check the GEO submission if you need it.

---

## 2. Core pipeline — the order that matters

```
FASTQ
  │
  ├─ [0] STARsolo v2.7.9a → mm39/GRCm39 → raw UMI count matrix
  │
  ├─ [1] SoupX v1.6.0            ambient RNA removal      (per library, on raw+filtered)
  │
  ├─ [2] scds  +  Scrublet       doublet scoring           (per library, pre-merge)
  │
  ├─ [3] merge libraries within cohort
  │      + custom ROSA26 / tdTomato read annotation
  │
  ├─ [4] cell QC filters   %mito ≤15 · nFeature ≥50 · nFeature ≤ MAD-based ceiling
  │
  ├─ [5] CellCycleScoring
  │
  ├─ [6] SCTransform  (regress: percent.mt, nFeature_RNA, nCount_RNA)
  │
  ├─ [7] PCA → Louvain clustering (res 1.0) → UMAP
  │
  ├─ [8] marker-based annotation → *manual cluster culling* → final atlas
  │
  └─ downstream modules (§4): subclustering · trace calling · scoring ·
                              CellChat · gseGO · Slingshot+tradeSeq
```

The three orderings that are easy to get wrong and that this paper gets right:

1. **SoupX before doublet detection.** Ambient contamination inflates the appearance of mixed-lineage expression, which is exactly the signal doublet callers key on. Deconvolve the soup first or you will over-call doublets in the compartments with the most ambient signal (here: the epithelium, because surfactant transcripts dominate lung soup).
2. **Both SoupX and doublet calling run per library, before merging.** The ambient profile is a property of one 10x channel, and the expected doublet rate is a property of one channel's loading density. Pooling first destroys both.
3. **Doublet calling before, not after, the count-based QC filters.** Doublets sit at the *high* end of the nFeature distribution; if you cut the tail first you have thrown away the evidence.

---

## 3. Stage-by-stage detail

### [0] Alignment — STARsolo v2.7.9a

Reference **mm39 / GRCm39**. Note the mouse reference is a recent build; if you compare to public lung datasets aligned to mm10, gene-level identifiers will mostly transfer but coordinates and some gene models will not.

*Why STARsolo rather than CellRanger:* it produces both raw (all barcodes) and filtered matrices and is faster; SoupX needs the raw matrix, so keep it.

**Custom step you must not skip:** reads mapping to the ROSA26 locus were annotated separately, because the lineage-trace readout depends on distinguishing the *recombined* (tdTomato-expressing) allele from the unrecombined one. This requires a modified reference or a supplementary counting pass — see §4.2.

### [1] Ambient RNA — SoupX v1.6.0

**Reference:** `Young-2020-Soupx-removes-ambient-rna-contamina.pdf`

Input: raw (droplet-level, unfiltered) matrix **including empty droplets** + filtered cell matrix + a preliminary clustering. Output: corrected integer count matrix, a drop-in replacement for the original counts.

SoupX is the **only** tool in this pipeline that needs the unfiltered matrix — the ambient profile is estimated from empty droplets. Keep STARsolo's raw output; do not discard it after cell calling.

**Parameters from Young & Behjati 2020:**

| Parameter | Value / guidance |
|---|---|
| `N_emp` (UMI ceiling defining "empty") | Any value **< 100** works; correlation with true background is best at **< 10**. The paper used ≤ 10. |
| Contamination fraction ρ | Assumed **constant within a 10x channel**. Estimate per channel, never pooled. |
| Estimation mode | Automated (`autoEstCont`) takes the **modal** ρ across many single-marker estimates; manual mode uses bimodal soup genes. |
| Manual marker choice for solid tissue | **Haemoglobin genes (Hbb-bs, Hba-a1…)** — the authors name Hb as "a sensible choice for most solid-tissue experiments." Rank the top 500 background genes by bimodality to find others. |
| Cell selection for manual mode | Poisson test at **FDR 0.05**; exclude the *entire cluster* containing any cell with genuine expression. |

The automated mode requires a clustering, but only to find markers, so "consistent estimates will be obtained for any sensible clustering" — a rough Louvain pass on uncorrected counts is enough. This means the pipeline has an unavoidable loop: raw counts → quick cluster → SoupX → re-run the full workflow on corrected counts.

**Why this step is load-bearing here, not hygiene.** The lung is a hard case for ambient RNA: dissociated AT2 cells release surfactant transcripts (*Sftpc*, *Sftpb*, *Scgb1a1*) at extreme abundance, and these appear as low-level "expression" in endothelial and immune cells. This paper's central claim is that endothelial cells acquire non-canonical programs (MHC-II, glycolysis, *Sparcl1*, *Ntrk2*). Uncorrected soup is a direct alternative explanation. Two results from the SoupX paper make the point concretely: widespread *HBB* outside the erythroid lineage in fetal liver turned out to be **contamination, not doublets**, and correcting the soup **increased cross-batch mixing entropy** — i.e. ambient RNA manufactures batch effects, which in a 8-time-point design would look exactly like time-dependent biology.

Useful property: **overcorrection is comparatively safe.** Contamination is preferentially removed from genes closest to background, so an aggressive global ρ is "unlikely to completely remove the expression of genes that are truly markers." When in doubt, correct harder.

### [2] Doublets — scds AND Scrublet

**References:** `Bais-2020-Scds-computational-annotation-of-do.pdf`, `Wolock-2019-Scrublet-computational-identificati.pdf`

The paper runs **both**. The two are genuinely orthogonal: scds's `cxds` is a co-expression/annotation-based score computed directly on binarized counts, while Scrublet simulates doublets and asks how many of a cell's nearest neighbours are simulated. Running both and reviewing the disagreements is more conservative than either alone, and it is what lets the authors later *defend keeping* one suspicious population and discarding another.

**Parameters:**

| Tool | Setting | Value |
|---|---|---|
| **Scrublet** | expected doublet rate `r̂` | From the 10x loading table — roughly 0.8% per 1,000 cells recovered. Not a constant; set per library. |
| | `k` (neighbours) | `round(0.5 · sqrt(n_cells))` |
| | `r` (simulated:observed ratio) | **≥ 2**; the paper's own analyses used 5–10 |
| | HVG defaults | ≥3 counts in ≥3 cells, top 15% by V-score |
| | threshold | Set on the histogram of **simulated** doublet scores, between its two peaks. `skimage.filters.threshold_minimum` automates it, but the authors "still recommend visual inspection." |
| **scds** | `binThresh` | 0 |
| | `ntop` | 500 most variable genes (Binomial variance) |
| | `hybrid` score | min-max normalise `cxds` and `bcds` to [0,1] and **add** — no weighting parameter |
| | threshold | **None produced.** scds ranks cells; it does not call doublets. You must impose a cutoff, normally from the expected 10x rate. |

**Two self-consistency checks Scrublet gives you for free**, both worth reporting: (1) the simulated-score histogram should be **bimodal** — if it isn't, `r̂` is probably wrong; (2) the detected fraction should satisfy `d ≈ r̂ · f_D`, where `f_D` is the detectable (neotypic) doublet fraction — if it doesn't, the threshold is wrong. Scrublet is also worth **running iteratively**: in the paper's bone marrow example, a second pass found 34 more doublets hiding inside a dense doublet cluster that the first pass missed.

**The shared blind spot, which is exactly this paper's problem.** Both tools detect only **heterotypic / neotypic** doublets — doublets of two transcriptionally similar cells are, in Wolock's words, "virtually indistinguishable from singlets." Both also assume every cell state contributing to a doublet is present as a singlet in the data. In scds's benchmark, every dataset except the species-mixing control contained a sizeable set of experimentally annotated doublets that **no method recovered**. So a clean doublet score is weak evidence of a genuine cell state, and this is precisely why the paper could not settle the AT1_AT2 and CAP1/CAP2 questions on doublet score alone.

**This is the step that most affects this paper's biology.** Two cases show why:

- **AT1_AT2 cells (kept).** A cluster co-expressing AT1 and AT2 markers. Retained because UMI and feature counts were only modestly elevated and *in line with other epithelial cells*, doublet scores were unremarkable, and the population has been independently described in mouse and human. It became a real finding — a poised/partially differentiated state.
- **CAP1/CAP2 mixed endothelial cluster (discarded).** Co-expressed CAP1 and CAP2 markers but had *dramatically* higher UMI and feature counts. Removed, with the honest caveat that the authors could not distinguish true doublets from binucleated or syncytial cells.

The lesson to carry into your own analysis: **doublet score alone did not decide either case.** The decision rule was *doublet score + UMI/feature distribution relative to the compartment + prior literature*. Write that rule down explicitly in your methods.

### [3] Merge + lineage annotation

Libraries are merged **within cohort** (Mki67 atlas separate from the Kit/Car4/Ednrb experiment). Variable feature selection and ROSA26 read annotation happen at this point.

### [4] Cell QC filters

A cell is **excluded** if any of:

| Filter | Threshold |
|---|---|
| Mitochondrial reads | > 15% |
| Feature count | < 50 |
| Feature count | > 2 × median absolute deviation |

Two things to notice.

**The lower feature bound of 50 is extremely permissive** by usual standards (200–500 is typical). This is a deliberate choice for an injury time course: dying, stressed and low-RNA cells are part of the biology at 6 dpi, and a hard floor of 500 would preferentially delete the injured state. The cost is admitted later — the paper's final culling step removes "clusters without a strong marker gene signature, which were often near the low end of the feature count cutoff." So the permissive floor is paid for with a *cluster-level* rather than cell-level cleanup.

**The upper bound is MAD-based, not a fixed number.** Written as "> 2 × MAD" in the paper; the standard formulation is *median + n·MAD*. Use `scuttle::isOutlier(..., nmads = 2, type = "higher")` or compute `median(x) + 2*mad(x)` explicitly and state which you used. A 2-MAD ceiling is aggressive — it will trim genuinely high-complexity cells — which is consistent with the authors also relying on it as a doublet backstop.

> **Deviation worth considering for your own run.** Current scverse best practice applies MAD-based cutoffs to *log-transformed* counts and typically uses 3 or 5 MADs, computed **per sample** rather than across the merged object. With 8 time points and a real biological gradient in cell complexity, a single pooled threshold will filter time points unevenly. If you reproduce this, compute the MAD per library and report both versions.

### [5] Cell cycle

`Seurat::CellCycleScoring()` — S and G2M scores assigned. Note the paper **scores but does not regress** cell cycle. Correct choice here: proliferation *is* the phenotype (the entire Ki67-tracing design exists to capture it), so regressing it out would destroy the signal.

### [6] Normalisation — SCTransform

`SCTransform(vars.to.regress = c("percent.mt", "nFeature_RNA", "nCount_RNA"))`.

Regressing `nFeature_RNA` *and* `nCount_RNA` on top of SCTransform's own depth model is more aggressive than the SCTransform default (which already handles sequencing depth). Given the permissive feature floor, this is a reasonable defence against complexity-driven clustering, but be aware it is a non-default choice and can flatten real differences in transcriptional output between quiescent and activated states.

**No batch-correction / integration algorithm was used** — no Harmony, no CCA/RPCA anchors, no scVI. Libraries are merged and SCTransform-normalised only. This is a substantive methodological point to defend or revisit: it means time point and batch are confounded to whatever degree they are confounded. The authors' implicit justification is the observation that "similar cellular contributions to each compartment from each collection day" were seen after subsetting by time point — i.e. no time point collapsed into its own island. If you re-run this, **test integration explicitly** and show that the iCAP cluster is not a batch artefact; that is the single most obvious reviewer challenge to the paper's central claim.

### [7] Clustering and embedding

- Louvain graph-based clustering.
- **Resolution 1.0** for the initial clustering of all cells.
- **Resolution 0.4–1.2** for re-clustering within a lineage (chosen empirically per compartment).
- UMAP for 2D projection.

Number of PCs retained is **not stated** in the paper. You will have to choose and document it (elbow/JackStraw, or a fixed 30–50).

### [8] Annotation and cluster culling

Identities assigned from canonical markers, cross-referenced against published references and **LungMAP**. `FindMarkers()` for inter-cluster differences.

Clusters were then **removed** on three grounds:
1. No strong marker gene signature (often low-feature-count clusters — see §4 above);
2. Markers from multiple cellular compartments in the same cell *despite passing doublet thresholds*;
3. Rare cell types out of scope: platelets, erythrocytes, mast cells, basophils.

This manual culling is a real analytical step with real degrees of freedom. Record the cluster IDs you removed and why.

---

## 4. Downstream modules

These hang off the annotated object. They are largely independent of one another; the dependency is only that each requires a *subset and re-cluster* first.

### 4.1 The subset-and-recluster motif

Nearly every biological result in this paper comes from the same move: **take one compartment, re-cluster it at higher resolution, and look again.** The paper is explicit that this is necessary — the iCAP state "was distinguished only when examining the subset of capillary ECs and not in the larger dataset."

Compartments subset in the paper:

| Subset | Resolution outcome |
|---|---|
| Myeloid | iMON_a/b, aMAC_a/b, iMAC, cDC, pDC, neutrophils |
| Alveolar epithelium (AT1, AT2, transitional, AT1_AT2; airway excluded) | AT1_a/b/c, AT2, transitional |
| T/NK lymphocytes | CD8-TRM identified |
| Capillary endothelium | CAP1_a–d, iCAP_a–c, CAP2 |

Build this as a reusable function, not as copy-pasted scripts. Input: parent object + cell-type vector + resolution. Output: re-normalised, re-clustered child object carrying its parent annotations.

### 4.2 Lineage-trace calling from transcripts

This is the most bespoke piece of the pipeline and is worth reproducing carefully.

Cells are annotated from reads mapping to either the **871 bp region overlapping the 3× SV40 polyA** or the **1 kb region overlapping the bGH polyA** downstream of the tdTomato coding sequence. These two regions were chosen because of higher sequencing depth.

| Call | Rule |
|---|---|
| **Traced** | > 50% of ROSA26-locus transcripts derive from the recombined allele |
| **Not_detected** | No ROSA26-locus transcripts at all |
| **Untraced** | ROSA26 transcripts present, but < 50% from the recombined allele |

The same 50% threshold was used for every Cre line. The three-way call — rather than binary traced/untraced — is what makes the tracing quantitative: `Not_detected` cells are censored, not counted as negative, so tracing percentages are computed on cells with evidence.

Recombination efficiency and specificity are line-dependent and should be reported alongside any tracing percentage: Car4<sup>CreERT2</sup> 50.6% efficient / 74.6% specific; Ednrb<sup>CreERT2</sup> 50.1% / 91.2%.

### 4.3 Gene signature scoring

**Interferon score:** built from **87 interferon-stimulated genes**, computed per cell (`AddModuleScore`), then compared across compartments and time points. This produced the paper's key mechanistic link — the endothelium carries the highest IFN score of any non-immune compartment at 6 dpi, which explains its delayed proliferation.

**iCAP score:** top **50 iCAP-specific genes**, applied to an external organ-specific EC atlas to ask whether iCAPs resemble other-organ endothelium. Note the direction of use — the signature is defined in-dataset and then *applied outward* to public data. Keep the gene list as a versioned file.

Other module scores used: mature AT1 function, collagen IV / laminin ECM, MHC-II.

### 4.4 Cell–cell communication — CellChat v1.6

**Cell types were evenly downsampled before analysis.** This is essential and easy to forget: CellChat's interaction probabilities scale with cluster size, so unequal cluster sizes across a time course produce spurious differences in signalling strength.

Ingoing / outgoing / combined signalling matrices were pulled manually out of `netAnalysis_signalingRole_heatmap()` and re-plotted with **ComplexHeatmap** + **viridis** rather than using CellChat's own plotting.

### 4.5 Enrichment — clusterProfiler

`gseGO()` via **clusterProfiler**, with **enrichplot** and **DOSE**. **10,000 permutations.** Input gene lists were the genes varying most along each pseudotime axis (see 4.6) — so the enrichment is trajectory-driven, not cluster-marker-driven. That coupling is the point: it describes what *changes during a transition*, not what *marks an endpoint*.

### 4.6 Trajectories — Slingshot → tradeSeq

**References:** `Street-2018-Slingshot-cell-lineage-and-pseudoti.pdf`, `Van-den-berge-2020-Trajectory-based-differential-expre.pdf`

**The single most important methodological decision here: the trajectories are *supervised*.** From the paper: analyses were "supervised to generate trajectories corresponding to known or hypothesised axes of differentiation… accomplished by removing other cell types or cellular subtypes from the data set during trajectory identification." Cells not included are rendered grey in the figures.

This is a legitimate and common approach, but it must be declared plainly, because it means **the trajectory does not discover the lineage — it measures progression along a lineage you have already asserted.** Every trajectory in this paper is a hypothesis test of a pre-specified axis, not an unbiased lineage inference.

Trajectories built:

| Trajectory | Start → End | What it established |
|---|---|---|
| Myeloid A | iMON_b → iMON_a → aMAC_b → aMAC_a | Bone-marrow-derived monocytes reconstitute alveolar macrophages; lipid/surfactant catabolism programs turn on |
| Myeloid B | iMON_b → iMAC | Contrasting axis; host-defence, complement, MHC-II programs |
| Alveolar epithelium | AT2 → transitional → AT1 | Identified AT1_c as an immature, persistent AT1 state |
| Capillary EC | CAP1 / CAP2 → iCAP | Bidirectional origin of the iCAP state |

#### Slingshot inputs and the two knobs that decide the answer

Slingshot takes **normalised** expression in a reduced-dimension space, **cluster labels**, and a **required start cluster** (end clusters optional). It builds a minimum spanning tree on cluster centroids using a covariance-scaled, Mahalanobis-like distance — not plain Euclidean — then fits simultaneous principal curves. Lineages are root-to-leaf paths through the MST, so the topology is by construction a single connected tree.

Two choices dominate the result, and neither is stated in the Niethamer methods:

- **Number of clusters K.** In Street's own simulations, too few clusters made Slingshot **miss the branching event entirely** (K = 3 gave pseudotimes matching neither true lineage), while too many produced spurious branches. The clustering *algorithm* mattered far less than K. Since this paper re-clusters each compartment at a resolution chosen "empirically" between 0.4 and 1.2, K is effectively a free parameter — report the resolution and resulting K for every trajectory.
- **Dimensionality reduction and number of dimensions.** Street reports "a potentially large impact on the final result": on identical simulated data, 3-D PCA gave highly variable accuracy while 4-D PCA gave consistently high accuracy. There is no safe default.

Slingshot yields point estimates only — no uncertainty on either lineage assignment or pseudotime. It also produces per-cell **lineage weights** from each cell's projection distance to each curve; those weights are the input tradeSeq consumes.

#### tradeSeq — what to feed it and which test to use

`fitGAM()` takes exactly three things: the **raw count matrix** (not normalised — depth enters as a `log(N_i)` offset inside the negative-binomial GAM), the per-cell pseudotimes, and the per-cell lineage weights from Slingshot. Known covariates (batch, time point) can be added as fixed effects — **this is tradeSeq's only condition-handling mechanism in the 2020 paper**; `conditionTest` was added to the package later and is not described there.

| Parameter | Value |
|---|---|
| Knots `K` | Package default **6**; select with `evaluateK()`, which uses **AIC** (BIC was rejected as favouring overly complex models) |
| `K` values chosen in the paper's own case studies | 3–6; both real 10x datasets landed on **6** |
| Knot placement | Quantiles of pseudotime by default |
| Smoothing parameter λ | Generalized cross-validation; **shared across lineages**, as are the spline basis functions — this is what makes between-lineage contrasts valid |
| M (evaluation points for pattern tests) | 100 |
| FDR | Benjamini–Hochberg; use **`stageR`** for stage-wise control when running several tests per gene |

**Which test answers which question** — the mapping the paper leaves implicit:

| Question | Test |
|---|---|
| Does this gene change at all along the lineage? *(this is the paper's "genes that varied along the pseudotime trajectory" — the input to gseGO)* | `associationTest` |
| What separates progenitors from the differentiated endpoint? | `startVsEndTest` |
| What are the markers of the two lineages' terminal cell types? | `diffEndTest` |
| Do two lineages differ anywhere along their length? *(best screening test; the only one insensitive to the choice of K)* | `patternTest` |
| What changes around a specific branch point or window? | `earlyDETest(knots = c(a, b))` |
| Which genes diverge transiently and then re-converge? | `patternTest` significant **and** `diffEndTest` not significant |

For the CAP1/CAP2 → iCAP question specifically, `patternTest` is the right screen and `earlyDETest` around the divergence is the right follow-up; `diffEndTest` alone would miss the fact that both lineages arrive at a similar iCAP endpoint.

#### One caveat that bites this paper directly

The paper rescales pseudotime to 0–1 ("normalised pseudotime") so that trajectories can be compared on one axis. **A linear 0–1 rescale is not the same as making pseudotimes comparable across lineages**, and Van den Berge et al. warn about exactly this: if one lineage is longer, a gene with a genuinely identical expression pattern can occupy 75% of the short lineage and 25% of the long one, and `patternTest` will call it differentially expressed. Their recommended remedy is **dynamic time warping** before comparison, not linear rescaling. Note also that `patternTest` already length-normalises internally by evaluating at M equally spaced points, so the 0–1 rescale is a *plotting* convenience and should not be relied on as a statistical correction.

Two further caveats to state in your methods: pseudotimes are treated as **fixed and known**, so no uncertainty from the trajectory step propagates into the p-values; and because the same data are used for trajectory inference and for differential expression, the authors themselves recommend treating the p-values "simply as useful numerical summaries for ranking the genes for further inspection."

**Still unspecified in the paper and therefore yours to decide and report:** the reduction and number of dimensions given to `slingshot()`, whether end clusters were fixed or only the start, the `k` used in `fitGAM()`, and which tradeSeq test produced the gene lists sent to gseGO.

### 4.7 Isoform-level check (one-off)

For *Ntrk2*, BAM files from two mice (uninjured control, 366 dpi Ki67-traced) were processed with **deepTools** into bigWig coverage tracks and viewed in the UCSC browser, to show that iCAPs express only kinase-deficient *Ntrk2* isoforms. A useful pattern: when a marker gene's *function* matters, 10x 3′ data can still answer isoform questions at the coverage level even though it cannot quantify isoforms.

---

## 5. Annotation reference — the marker genes that define the atlas

Assembled from the paper's dot plots and text. Use this as the checklist when you re-annotate; if a cluster you generate does not map onto one of these, that is either a new finding or a QC failure, and you should decide which before proceeding.

**Compartment level:** Immune *Ptprc* · Mesenchyme *Col1a1* · Epithelium *Epcam* · Endothelium *Pecam1*

| Compartment | Cell type | Markers |
|---|---|---|
| **Endothelium** | Arterial | *Bmx* |
| | Venous | *Slc6a2* |
| | Lymphatic | *Mmrn1* |
| | **CAP1** | ***Gpihbp1*, *Kit*** |
| | **CAP2** | ***Ednrb*, *Car4*** |
| **Epithelium** | Ciliated | *Foxj1* |
| | Secretory | *Scgb3a2*, *Scgb1a1* |
| | AT2 | *Sftpc*, *Lamp3* |
| | AT1 | *Hopx* |
| | Krt5 (dysplastic/basal-like) | *Krt5* |
| | Alveolar transitional | *Krt8*, *Cldn4* |
| | AT1_AT2 | co-expression only — **no specific genes** |
| **Mesenchyme** | AF1 | *Pdgfra* |
| | AF2 | *Pdgfrb* |
| | Adventitial fibroblast | *Dcn* |
| | Peribronchial fibroblast | *Aspn* |
| | VSMC | *Cnn1* |
| | Mesothelium | *Msln* |
| **Myeloid** | aMAC | *Siglecf*, *Itgax*; mature marker *Krt79* |
| | iMAC | *C1qb*, *C1qa*, *Apoe*, *H2-Aa* |
| | cMON | *Cx3cr1*, *Ly6c2* |
| | pMON | *Adgre4* |
| | iMON | intermediate *Ly6c2* |
| | cDC1 / cDC2 / maDC / pDC | *Xcr1*,*Flt3* / *Cd209a* / *Ccr7* / *Ccr9* |
| | Neutrophil | *S100a8*, *Cxcr2* |
| **Lymphoid** | B / plasma / T / NK / ILC | *Cd79a*,*Ms4a1* / *Jchain* / *Cd3e* / *Nkg7* / *Il7r* |
| | **CD8 T-RM** | *Itgae* (CD103), *Itga1* (CD49a), *Cd69*, *Ifng*-high |

### The three injury-associated states — what to look for and when

| State | Appears | Resolves | Defining signature |
|---|---|---|---|
| **Alveolar transitional** (= Krt8⁺ ADI / DATP / PATS) | peaks **11 dpi** | rare after **25 dpi** | *Krt8*, *Cldn4*; most heavily Ki67-traced epithelial type |
| **AT1_c** — immature AT1 | from **11 dpi** | gone by **366 dpi** | ***Nnat***-high, *Ager* normal, *Fos*/*Jun* high. Contrast AT1_a (terminal): *Gpm6a*, *Igfbp2* |
| **iCAP** — injury-induced capillary EC | from **11 dpi** | **persists ≥ 366 dpi** | ***Sparcl1*, *Ntrk2*** (kinase-deficient isoforms only); high *Gpihbp1*/*Kit*, low *Ednrb*/*Car4*; MHC-II (*H2-Aa*, *Cd74*, *Ciita*), *Ifngr1*; *Ankrd37*, *Fxyd5*, *Bnip3*; glycolysis (*Pfkl*, *Tpi1*, *Ldha*, *Aldoa*, *Pkm*, *Eno1*, *Gpi1*) |

**Subcluster inventories** produced by the subset-and-recluster step (§4.1):

- Alveolar epithelium → AT2_a–d · Alveolar_transitional · AT1_AT2 · AT1_a · AT1_b · **AT1_c**
- Capillary endothelium → CAP1_a–d · **iCAP_a, iCAP_b, iCAP_c** · CAP2
  *CAP1_c/d are 6 dpi interferon-responders (*Atf3*, *Irf7*, *Isg15*, *Plvap*) and were excluded; iCAP_c was excluded for lacking* Ntrk2*. **iCAP_a is Kit(CAP1)-derived, iCAP_b is Car4/Ednrb(CAP2)-derived** — this is the paper's central result.*
- Myeloid → iMON_a/b · aMAC_a/b · iMAC · cDC · pDC · neutrophils
  *iMON_b (6 dpi): *Ifit3*, *Oasl1*, *Spp1*, *Il1b*, *Tnf*. iMON_a (11 dpi): *C1qa*, *Apoe*, *H2-Aa*. aMAC_a (mature, Car4-traced): *Krt79*, TFs *Pparg*/*Rxra*/*Bhlhe40*/*Bhlhe41*. aMAC_b (inflammatory, Kit-traced): low *Siglecf*.*

**Transcription factor analysis** used the **AnimalTFDB** gene list, taking the top 25 TFs differentially expressed over each pseudotime axis. Keep that list versioned alongside the ISG and iCAP signatures.

---

## 6. Environment

| Component | Version in paper |
|---|---|
| STARsolo | 2.7.9a |
| Reference | mm39 / GRCm39 |
| SoupX | 1.6.0 |
| scds | version not stated |
| Scrublet | version not stated (Python) |
| R | 4.3 |
| Seurat | 4.9 (development build between v4 and v5) |
| CellChat | 1.6 |
| clusterProfiler / enrichplot / DOSE | not stated |
| Slingshot / tradeSeq | not stated |
| ComplexHeatmap, viridis, MetBrewer | not stated |
| deepTools *(Ntrk2 coverage only)* | not stated |
| AnimalTFDB *(TF gene list)* | not stated |
| LungDamage MATLAB program *(histology, not scRNA)* | [github.com/WALIII/LungDamage](https://github.com/WALIII/LungDamage) |

**Practical notes for reproduction.**

- **Seurat v4.9 is a pre-release snapshot.** It is not on CRAN. If you install Seurat v5 you will hit the Assay5/BPCells object model, and `SCTransform` + `FindMarkers` behave differently around layers. Either pin `Seurat 4.4.0` (last stable v4) and accept small differences, or use v5 and say so.
- **The pipeline is bilingual.** Scrublet is Python, everything else is R. Plan for one conda environment (Python: scrublet, scanpy for I/O) and one R environment, exchanging `.h5ad`/`.rds` or writing doublet scores out as a CSV keyed by barcode. The CSV route is simpler and more auditable.
- **No R installation was detected on this machine** — only Python 3.12. You will need R ≥4.3 plus Bioconductor (`scds`, `slingshot`, `tradeSeq`, `ComplexHeatmap`, `clusterProfiler`) before any of this runs.

---

## 7. How to read the PDFs in this folder

Read them in pipeline order, not alphabetical order. For each, the question to answer is in the right column.

| # | PDF | Read it to answer |
|---|---|---|
| 0 | `Lv-2024-…+Suppl.pdf` — **actually Niethamer 2025** | The paper itself. Read STAR Methods (near the end) first, then Results. |
| 1 | Young 2020 — **SoupX** | How is the contamination fraction estimated, and what does the automated mode assume about my clustering? |
| 2 | Wolock 2019 — **Scrublet** | What expected doublet rate should I set, and why can't it see homotypic doublets? |
| 3 | Bais 2020 — **scds** | What do cxds, bcds and hybrid scores measure, and why does scds refuse to give me a threshold? |
| 4 | Street 2018 — **Slingshot** | How do cluster number and choice of reduced dimensions determine the lineage graph? |
| 5 | Van den Berge 2020 — **tradeSeq** | Which test answers which biological question, and how do I pick `k`? |

Two notes on what these papers do *not* contain, so you do not go looking:

- **Street 2018 never names an R function.** No `slingshot()`, `getLineages`, `getCurves`, `slingPseudotime`. It describes the two-stage algorithm those functions implement; the API is package documentation, not paper content.
- **Van den Berge 2020 contains no `conditionTest`.** It was added to tradeSeq after publication. The paper's only mechanism for handling condition or batch is the fixed-effect covariate matrix `U` in `fitGAM`. If you want to compare trajectories *between* time points, that is beyond what the cited paper supports and needs its own justification.

---

## 8. Open decisions — what you must specify that the paper does not

Each of these is a place where two people following this paper would produce different numbers. Collect answers before you write your own methods section.

1. **Number of principal components** retained before clustering and UMAP.
2. **UMAP hyperparameters** — `n_neighbors`, `min_dist`, metric.
3. **Number of HVGs**, and SCTransform v1 vs v2.
4. **MAD formulation** — `median + 2·MAD` vs `2 × MAD`, on raw or log counts, per-sample or pooled.
5. **SoupX contamination estimation** — automated modal estimate vs manual marker-based (haemoglobin genes are the recommended solid-tissue choice), and the resulting ρ per library.
6. **Doublet thresholds** — Scrublet's `r̂` and call threshold; the scds cutoff (scds gives none); and the rule for combining the two callers (union, intersection, or manual review).
7. **`FindMarkers` settings** — test, log-fold-change and p-value thresholds. Defaults are Wilcoxon, but the paper does not say.
8. **Slingshot inputs** — reduction type, number of dimensions, and whether end clusters were fixed or only the start.
9. **`fitGAM` knots `k`**, and **which tradeSeq test** produced the gene lists sent to gseGO.
10. **Integration** — whether you accept the paper's no-correction approach or add and justify one. Infection round (three dates, in Table S3) is the obvious covariate to test.
11. **Pseudobulk DE method** for the human disease comparison — no package is named.
12. **The 87 ISG list, the top-50 iCAP gene list, and the AnimalTFDB TF list** — none is reproduced in the methods text. Extract and version them.

**Confirmed absent from the paper, so do not go looking:** RNA velocity (not performed), bleomycin injury (not used, only cited), CellRanger (STARsolo instead), reads per cell, per-time-point cell recovery counts, and any integration/batch-correction step.

---

*Compiled from the published STAR Methods, Results and Key Resources Table of Niethamer et al. 2025 (PMC11805657), retrieved via PubMed Central. Article DOI: [10.1016/j.stem.2024.12.002](https://doi.org/10.1016/j.stem.2024.12.002).*
