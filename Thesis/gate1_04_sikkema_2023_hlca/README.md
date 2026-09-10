# Study note: Sikkema et al. 2023, An integrated cell atlas of the lung in health and disease (HLCA)

Roadmap position: **Gate 1, paper 4 of 4 ("reference framework")** in the
ordered paper roadmap of the
[UC Berkeley SAP PI Target Map](https://app.notion.com/p/3d1151616b44814c8697ff3af2f8f831).
The owner's own reading note is the Notion page
[HLCA study note](https://app.notion.com/p/3d6151616b4480708b3ed6bd7fbb2949);
this file is the repository-side expansion of the two headings that note
leaves sparse: *Integration method / benchmarking pipeline* and *list of
scores / decision criteria*. The machine-readable version of every number
below is [`integration_benchmark.json`](integration_benchmark.json); the
consequences for this repository's pipeline are in
[`PIPELINE_FRAMING.md`](PIPELINE_FRAMING.md); the first analysis trial is in
[`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md).

**Citation.** Sikkema L, Ramirez-Suastegui C, Strobl DC, Gillett TE,
Zappia L, Madissoon E, Markov NS, Zaragosi L-E, ... Luecken MD, Theis FJ,
for the Lung Biological Network Consortium. *An integrated cell atlas of the
lung in health and disease.* **Nature Medicine** 2023;29(6):1563-1577.
DOI [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2),
PMID [37291214](https://pubmed.ncbi.nlm.nih.gov/37291214/). Open access,
CC BY 4.0.

**Companion method paper (the benchmark).** Luecken MD, Buttner M,
Chaichoompu K, et al. *Benchmarking atlas-level data integration in
single-cell genomics.* **Nature Methods** 2022;19:41-50. DOI
[10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8),
PMID 34949812, PMC8748196. The HLCA uses this "scIB" pipeline unchanged.

**Resources named in the paper.**

| Resource | Where |
|---|---|
| HLCA counts, embedding, annotations, metadata | cellxgene collection `6f6d381a-7701-4781-935c-db10d30de293` |
| HLCA core reference model and embedding (for mapping new data) | Zenodo [10.5281/zenodo.7599104](https://doi.org/10.5281/zenodo.7599104) |
| Landing page, reproducibility code, mapping tutorial, count pipeline | github.com/LungCellAtlas/{HLCA, HLCA_reproducibility, mapping_data_to_the_HLCA, scRNAseq_pipelines} |
| Automated mapping portal | FASTGenomics (`beta.fastgenomics.org/p/hlca`); Azimuth and CellTypist models also exist but are not evaluated in the paper |

**Local files (gitignored, never committed).** Main PDF plus supplementary
items MOESM1 to MOESM8. What each contains and how it can be reused is in
[`SUPPLEMENTARY_INVENTORY.md`](SUPPLEMENTARY_INVENTORY.md).

> Provenance of this note: extracted from the paper text (pdftotext), the
> supplementary figure captions, Supplementary Fig. 1 (read from the rendered
> image), and Supplementary Tables 1 to 16 (MOESM3 workbook) on 2026-09-09 by
> an AI session; not yet reviewed by the owner. Numbers are quoted from the
> paper, not recomputed. Where the paper contradicts itself the contradiction
> is recorded, not resolved.

---

## 0. Reading workflow (the five questions the roadmap requires)

- **Question tested.** Can 14 healthy lung scRNA-seq datasets from different
  labs, protocols and platforms be integrated into one reference whose cell
  identities are a cross-study consensus rather than one lab's opinion, and
  can that reference then annotate new data and expose disease states?
- **Evidence type.** Computational only. Every claim is a property of
  integrated transcriptomes: benchmark scores, cluster entropy, expert
  consensus, label-transfer uncertainty, mixed-model coefficients, sc-LDSC
  enrichment, CIBERSORTx deconvolution. No lineage tracing, perturbation or
  spatial data was generated; one spatially annotated external dataset (ref.
  40) was mapped to test whether spatial cell types can be recovered.
- **Reusable variables.** The 5-level annotation hierarchy and its
  label-harmonisation table; the consensus marker sets (Supplementary Table 6,
  local MOESM3 sheet 6); the AT0 and pre-terminal-bronchiole secretory
  definitions; the batch-associated gene list (Supplementary Table 7); the
  covariate encodings; the label-transfer uncertainty cutoff (0.2); the IPF
  alveolar-fibroblast and SPP1-high monocyte-derived macrophage signatures
  (Supplementary Tables 13 to 15).
- **One limitation for GSE262927.** The HLCA is human and adult; the authors
  state that mapping mouse or in vitro data "may require further method
  development". Nothing in this paper licenses mapping the mouse influenza
  time course onto the HLCA. Its *decision criteria* transfer; its *reference
  embedding* does not.
- **One bridge.** GSE178360 (the human distal-lung series already analysed in
  this repository) is an HLCA *extension* dataset ("Tata_unpubl": 3 subjects,
  3 samples, 21,700 cells, 10x 3' v3), not a core dataset, so the HLCA core
  model never saw it; mapping it to the core with scArches and comparing the
  transferred labels with this repository's blind annotations is a legitimate,
  non-circular test of the human analysis, including the AT0 call.

---

## 1. What the paper built, checked against the text

| Object | Content |
|---|---|
| HLCA core | 14 datasets (11 published, 3 unpublished), 107 individuals, 166 samples, 584,444 cells; healthy tissue only (asthma, fibrosis and similar excluded); integrated with scANVI |
| HLCA extended | 37 further datasets (1,797,714 cells, 380 individuals; healthy and diseased, cells and nuclei, 10x 3'/5', Drop-seq, Seq-Well) mapped onto the core with scArches; total 49 datasets, 486 individuals, 2.4 million cells. Fig. 1 says "35 additional datasets"; the Results text says 37. Supplementary Table 1 totals 487 subjects, 743 samples, 2,382,658 cells |
| Consensus annotation | 5-level hierarchy, 61 cell identities, each present in at least 4 of 14 core datasets; 41% of cells kept their original label, 28% were refined, 31% substantially re-annotated |
| New or rare identities | migratory DCs (n=312; CCR7, LAD1, COL19 family), haematopoietic stem cells (n=60; SPINK2, STMN, PRSS57, CD34), hillock-like epithelium (n=4,600; KRT6A, KRT13, KRT14), AT0 (n=1,440; SFTPB+, SCGB3A2+, SFTPC-high, SCGB3A1-low), pre-terminal-bronchiole secretory (n=4,393; SFTPB+, SCGB3A2+, SFTPC-low, SCGB3A1-high), FAM83D+ smooth muscle (n=335); ionocytes, tuft and neuroendocrine cells (0.08, 0.01, 0.02% of the core) each resolved into their own cluster |
| Covariate modules | cell-type-specific gene programmes for sex, age, BMI, smoking and proximal-distal position (CCF score) from a pseudo-bulk mixed model |
| Reference uses | label transfer with uncertainty; GWAS to cell type (sc-LDSC); bulk deconvolution (CIBERSORTx); IPF alveolar-fibroblast state; SPP1-high profibrotic monocyte-derived macrophages shared by IPF, late COVID-19 and carcinoma |

---

## 2. Integration method and benchmarking pipeline (the scheme)

### 2.1 The order of operations

```
raw count matrices from 14 core datasets (shared as counts, not FASTQ; 4 realigned to GRCh38/Ensembl 84)
  |
  +- gene-type filter (Cell Ranger Ensembl 84 gene types, 33,694 IDs)
  +- cells with < 200 genes removed; genes in < 10 cells removed
  +- SCRAN size factors, estimated WITHIN Louvain (r 0.5) clusters built on 50 PCs, k 15
  |     cells with size factor < 0.01 or > 1e6 normalised counts removed (267 cells)
  |
  +- decide the BATCH UNIT: study, or study split by 10x version / processing site?
  |     principal-component regression (PCexpl) vs 10 random shufflings; split if > 1.5 SD above random
  |
  +- consensus HVG ranking across datasets (Cell Ranger flavour per dataset) -> 2,000 (HVG) or 6,000 (FULL)
  |
  +- BENCHMARK on a subset: 10 studies / 13 datasets, 72 donors, 124 samples, 372,111 cells
  |     10 tools x {HVG, FULL} x {scaled, unscaled} x {gene, embed, graph output} = 35 runs
  |     scored by 12 scIB metrics; overall = 0.4 * batch-removal + 0.6 * bio-conservation
  |
  +- WINNER: scANVI (coarse level-3 labels, 2,000 HVGs, unscaled raw counts, batch = dataset)
  |     30 latent dimensions -> kNN k 30 -> Leiden r 0.01 (level 1) -> nested Leiden r 0.2 (levels 2 to 5)
  |
  +- ENTROPY SCREEN of the ~100 clusters: label entropy, donor entropy; 6 doublet clusters dropped -> 94
  +- expert consensus annotation (6 lung biologists) -> 61 identities on 5 levels
  +- pseudo-bulk marker selection, 5 hierarchical rounds
  |
  +- scArches "surgery" on the frozen scANVI model -> map 37 new datasets
        kNN label transfer (k 50) + per-cell uncertainty; uncertainty > 0.2 => "unknown"
```

Two properties of this order matter more than any parameter. The batch
*unit* is decided by a test before integration, not assumed to be the sample
or donor; and the integration method is chosen by a scored benchmark on a
subset before the full atlas is built. Both are decisions this repository
currently makes by argument (see [`PIPELINE_FRAMING.md`](PIPELINE_FRAMING.md)).

### 2.2 Benchmark design

| Element | Value in the paper |
|---|---|
| Framework | scIB (Luecken et al. 2022), scIB v0.1.1; default integration parameters unless stated |
| Data | benchmark subset only: 10 studies split into 13 datasets (Results text says "12 datasets"; Methods enumerate 13 after splits), 72 donors, 124 samples, 372,111 cells |
| Cell and gene filter for the benchmark | cells with >= 500 UMI; genes in >= 5 cells |
| Tools (10) | BBKNN, ComBat, Conos, fastMNN, Harmony, Scanorama, scANVI, scVI, Seurat v3 RPCA, scGen |
| Why "12 methods" | Scanorama and fastMNN were each scored twice, once on their corrected gene matrix and once on their embedding output |
| Feature sets | HVG = 2,000 consensus HVGs; FULL = 6,000 (fastMNN could not run on FULL within memory) |
| Scaling | each non-count method run scaled (mean 0, SD 1) and unscaled; scVI and scANVI raw counts only; scGen 2,000 unscaled HVGs only |
| Input | scVI/scANVI: raw counts of HVGs; all others: SCRAN-normalised log data subset to HVGs |
| Supervision | scANVI and scGen received level-3 harmonised labels (level 2 for blood-vessel endothelium, fibroblast lineage, mesothelium, smooth muscle); 4,499 cells without a label at that level were dropped from the benchmark because scGen cannot take unlabelled cells |
| Batch variable | dataset (not donor, not sample) |
| Memory cap | 376 GB; runs exceeding it were excluded |
| Score | 12 metrics (4 batch, 8 bio); category score = mean of its metrics; overall = 0.4 batch + 0.6 bio; each metric min-max scaled to 0 to 1 across runs |

### 2.3 Ranking outcome (Supplementary Fig. 1, read from the figure; no numeric scores are printed)

| Rank | Method | Features | Scaling | Output |
|---:|---|---|---|---|
| 1 | scANVI (labels) | HVG | unscaled | embedding |
| 2 | scGen (labels) | HVG | unscaled | gene matrix |
| 3 | fastMNN | HVG | unscaled | gene matrix |
| 4 | scANVI (labels) | FULL | unscaled | embedding |
| 5 | BBKNN | HVG | unscaled | graph |
| 6 | Seurat v3 RPCA | HVG | unscaled | gene matrix |
| 7 | Scanorama | HVG | scaled | embedding |
| 8 | BBKNN | FULL | unscaled | graph |
| 9 | fastMNN | HVG | unscaled | embedding |
| 10 | ComBat | FULL | unscaled | gene matrix |
| 11 | ComBat | HVG | unscaled | gene matrix |
| 12 | Seurat v3 RPCA | FULL | unscaled | gene matrix |
| 13 | scVI | HVG | unscaled | embedding |
| 14 | scVI | FULL | unscaled | embedding |
| 18 | Harmony (best Harmony run) | HVG | scaled | embedding |

The full 35-row order is in `integration_benchmark.json`. Three readings
that matter: the two label-aware methods take ranks 1 and 2, so the win is
partly a win for supervision; the unsupervised scVI ranks 13 to 14; Harmony's
best run is 18th. The paper's second, independent argument for scANVI is
rare-cell resolution: after nested clustering, Harmony and Seurat RPCA "failed
to separate" ionocytes, tuft and neuroendocrine cells into their own clusters
(Supplementary Fig. 3b, recall and precision per level-3 cluster).

### 2.4 Final integration parameters (scANVI, scvi-tools 0.8.1)

2,000 HVGs, raw counts, batch = dataset, labels = coarse level-3 set with
"unlabeled" kept as a class (not removed). Two layers; 30 latent dimensions;
encode covariates True; deeply inject covariates False; layer norm both; batch
norm none; gene likelihood negative binomial; 500 unsupervised epochs then 200
semi-supervised; early stopping on ELBO (patience 10, threshold 0, LR
patience 8, LR factor 0.1) for the unsupervised phase and on accuracy on the
full dataset (patience 10, threshold 0.001) for the semi-supervised phase.
Gene-level analyses (DE, covariate models) use uncorrected counts; only
clustering and visualisation use the latent space.

### 2.5 Clustering

kNN k 30 on the 30 latent dimensions (the authors note this k "could impair
the detection of very rare cell types"); Leiden r 0.01 for level 1; for each
level-1 cluster a new kNN graph (k 30) and Leiden r 0.2 for level 2; levels 3,
4, 5 with k 15, 10, 10 and r 0.2; clusters named by lineage (1.2 = third
largest child of cluster 1). UMAP from the k 30 graph.

---

## 3. Scores and decision criteria, with the decision each one drives

| Score or rule | Definition | Threshold or setting | Decision it drives |
|---|---|---|---|
| scIB batch metrics (4) | PCR batch, batch ASW (silhouette), graph iLISI, graph connectivity | mean of the four = batch score | integration-method choice (weight 0.4) |
| scIB bio-conservation metrics (8) | NMI cluster/label, ARI cluster/label, cell-type ASW, isolated-label F1, isolated-label silhouette, graph cLISI, cell-cycle conservation, HVG conservation | mean of the eight = bio score | integration-method choice (weight 0.6) |
| Overall score | 0.4 * batch + 0.6 * bio | rank order | scANVI selected |
| Rare-cell recall and precision | per level-3 cluster, for ionocytes, tuft, neuroendocrine | reported, no cutoff | tiebreaker against Harmony and Seurat RPCA |
| PCexpl (principal-component regression) | fraction of variance in 50 PCs explained by a covariate (10x version or processing site) | > mean + 1.5 SD of 10 random shufflings | split a study into separate batches |
| Consensus HVG | genes highly variable in all datasets, then all-but-one, and so on | 2,000 (default) or 6,000 | feature set |
| Core cell filters | genes detected per cell; cells per gene; SCRAN size factor; normalised total | >= 200 genes; >= 10 cells; size factor >= 0.01; total <= 1e6 | cell and gene inclusion |
| Label entropy | Shannon entropy (natural log) of label fractions in a cluster, labelled cells only; NA if < 20% labelled | high if > 0.56 (entropy of a 75/25 split) | flags cross-study disagreement (33 of 94 clusters) and doublet clusters (6 removed) |
| Donor entropy | same over donor fractions | low if < 0.43 (95/5 split over all other donors) | flags donor-private clusters (14 of 94; used later to drop IPF cluster 5, donor entropy 0.17, 96% one donor) |
| Marker consistency filters | pseudo-bulk per sample and cell type (>= 10 cells; >= 3 for types with < 100 cells), t-test vs rest | expressed in >= 80% of pseudo-bulks; >= 50% of cells per pseudo-bulk (mean); <= 20% of out-group pseudo-bulks; 5 rounds relax these (Supplementary Table 16) | marker set per cell type |
| Minimum sample number | correlation of covariate variance with sample count | 40 samples | below it, covariate attribution is declared spurious and not reported |
| Variance inflation factor | multicollinearity between covariates at donor level | VIF > 5 excludes the cell type | which cell types get gene-level covariate models |
| Covariate encodings | smoking never 0 / former 0.5 / current 1; age 25 to 64 years mapped to 0 to 1; BMI 21.32 to 36.86 mapped to 0 to 1; nose as a binary; CCF 0 (nose) to 1 (alveolus) | fixed | model design |
| Mixed model | log(norm count) ~ 1 + age + sex + BMI + smoking + nose + CCF + (1 \| dataset); edgeR calcNormFactors, voom with mixed effects | BH within each covariate | gene modules per covariate; cameraPR GSEA on GO BP (MSigDB 7.1) |
| Label-transfer uncertainty | u = 1 - weighted fraction of the k = 50 reference neighbours carrying the transferred label | > 0.2 = "unknown" (ROC-calibrated on 12 datasets: TPR 0.879, FPR 0.495); 0.3 used in the two demonstration mappings | annotate vs abstain; per-dataset mean uncertainty decides "mapped well" (27 of 37) |
| Disease-signature contrast | cells with u < 0.2 vs u > 0.4 within one transferred label; Wilcoxon; top 20 up-regulated genes with FDR < 0.05 and mean expression >= 0.1 | fixed | disease signature score (scanpy score_genes style) |
| Cross-dataset clustering resolution | Leiden on joint 30-dim embedding, k 30 | r 0.3, chosen so that no dataset is isolated in one cluster | IPF fibroblast and MDM clusters |
| Cluster-marker DE filters | Wilcoxon on counts normalised to 7,666 (atlas median) | expressed in >= 30% in-cluster and <= 20% out-cluster; BH < 0.05 | Supplementary Tables 14 to 15 |
| Spatial-cell-type cluster search | resolutions 0.1 to 100; k 30 below r 25, k 15 above | min recall 25% and precision 25%; take highest recall unless precision drops > 33%; take the next cluster if precision doubles and recall drops <= 20% | recovering identities absent from the core |
| sc-LDSC | top 1,000 DE genes per cell type; 100 kb window; 1 cM LD window; HapMap 3 SNPs; X, Y and HLA excluded | BH per disease; depression GWAS as negative control | GWAS to cell type |
| Deconvolution inclusion | cell types > 2% of cells in the tissue; 200 cells sampled each; CIBERSORTx defaults | types with > 60% zero proportions excluded from testing; pseudo-bulk benchmark decides inclusion | which proportions are tested (Wilcoxon, BH) |

---

## 4. Annotation hierarchy (the "extract" item the roadmap asks for)

Five levels, built "recursively breaking up coarser labels into finer ones
where finer labels were available" (Supplementary Table 4 maps every
original label of every dataset onto the tree; local MOESM3 sheet 4).
Level 1 composition of the core: epithelial 48%, immune 38%, endothelial 9%,
stromal 4%. Of all cells, 94%, 66% and 7% carry a level-3, level-4 and
level-5 label respectively. Each cluster was annotated at the finest
defensible level and coarser levels were inferred upward. 18 of the 61 final
identities consist mostly of cells that were *mislabeled* in their source
study; the mislabels usually agree with the final call at a coarser level
(Supplementary Fig. 2). Rare identities the core could **not** separate:
regulatory T cells, gamma-delta T cells, ILCs, NKT cells, megakaryocytes;
these were mostly transferred as CD4 T, CD8 T or NK cells when present in
mapped data, which is a known ceiling for this reference on the Gate 3B
(Treg) branch of the roadmap.

Example of the epithelial branch (from Supplementary Table 5):
Epithelial > Airway epithelium > Secretory > Transitional Club-AT2 >
{AT0 (n=1,440), pre-TB secretory (n=4,393)}; and Epithelial > Airway
epithelium > Basal > {Basal resting, Suprabasal, Hillock-like}.

---

## 5. Reference mapping and uncertainty handling

- **Method.** scArches "architectural surgery" on the frozen scANVI model:
  the base model is not retrained; only adaptor weights for the new batch are
  learnt (freeze dropout True; surgery epochs 500; weight decay 0; early
  stopping on ELBO over the full dataset, patience 10, threshold 0.001).
- **Gene space.** The query is subset to the same 2,000 HVGs; HVGs absent
  from the query are zero-filled. Raw counts are used (one dataset that had
  ambient-RNA removal applied upstream was used as provided).
- **Label transfer.** kNN (k 50) in the joint embedding; the label is the
  edge-weighted majority; uncertainty u = 1 - weighted proportion of neighbours
  with that label.
- **Calibration.** Twelve mapped datasets with level-3 or level-4 original
  labels were harmonised to the HLCA tree; the ROC of "transferred label
  correct" against u gave the cutoff 0.2 near the elbow (TPR 0.879, FPR 0.495).
  Level-5 labels were excluded from calibration because their uncertainty
  reflects granularity, not batch effect.
- **What high uncertainty means in practice.** (i) transitions between
  continuous states; (ii) identities absent from the core (cancer cells,
  erythroblasts, Schwann cells, chondrocytes, nerve-associated fibroblasts);
  (iii) disease-altered states (IPF alveolar fibroblasts and macrophages);
  (iv) residual batch effect (nuclei, Drop-seq, COVID-19 and paediatric
  datasets carried higher mean uncertainty).
- **Demonstrations.** Healthy multimodal dataset: 68% correct, 14% incorrect,
  18% unknown at cutoff 0.3. Lung cancer: 77% correct, 1% incorrect, 22%
  unknown, with cancer cells and erythroblasts correctly left unknown; the
  claimed tumour-endothelium marker ACKR1 turned out to be a general venous
  marker in the healthy core.

---

## 6. Donor coverage (and why it sets the power of everything downstream)

Core: 107 individuals, 166 samples; harmonised ethnicity 65% European, 14%
African, 2% admixed American, 2% Asian, 2% mixed, 0.4% Pacific Islander, 14%
not annotated; smoking 52% never, 16% former, 15% active, 17% NA; sex 60%
male; age 10 to 76 years; BMI 20 to 49 (30% NA). Extended atlas: 486
individuals, but paediatric donors and several diseases are only in the
extension. Ethnicity is self-reported and harmonised to 1000 Genomes
super-populations; the authors call SNP-based ancestry "preferable". The
paper's own rule that covariate attribution needs at least 40 samples per
cell type, and that below that number technical and biological covariates
are confounded, is the single most transferable statement in the paper for a
25-sample mouse time course.

---

## 7. Shared profibrotic macrophage states (the last "extract" item)

Monocyte-derived macrophages (MDMs) from the core and all mapped datasets
(datasets or diseases with < 50 MDMs excluded) were clustered on the joint
embedding (k 30, Leiden 0.3). Four biological subtypes plus an intermediate:

| Cluster | Name | Markers | Enriched in |
|---|---|---|---|
| 0 | profibrotic, SPP1-high | SPP1, LPL, CHIT1, CHI3L1, MMP9, FDX1 | IPF; late COVID-19 with post-COVID fibrosis; carcinoma |
| 2 | inflammatory, CCL2-high | CCL2, IL1RN, S100A12 | early COVID-19 bronchoalveolar lavage |
| 4 | inflammatory, C1QA-high | CCL18, IL18, C1QA, TREM2 | COVID-19 pneumonia; carcinoma |
| 3 | MARCO-high (more differentiated) | MARCO, MCEMP1 | non-diseased samples |
| 1 | MDMs (intermediate) | | |

The IPF alveolar-fibroblast state was defined the same way (cluster 0 of the
joint fibroblast clustering: CCL2, COL1A1, CTHRC1, MMP19, SERPINE1, HIF1A;
Supplementary Table 14), after the donor-entropy rule removed a single-donor
cluster. Both signatures are candidate "persistent dysplastic remodelling"
programmes for the portfolio question, in human tissue only.

---

## 8. Biology extracted for the portfolio question

- **AT0 and pre-TB secretory** are defined as two children of "Transitional
  Club-AT2", with SFTPC and SCGB3A1 levels as the discriminators. This is the
  human counterpart of the transitional states the repository studies in
  mouse, and the population whose doublet-loss audit is documented in
  [`../../docs/DOUBLETS_AND_SCRUBLET.md`](../../docs/DOUBLETS_AND_SCRUBLET.md).
- **Interferon and sex.** IFNAR1 is lower in female lymphatic endothelium;
  the authors link it to sex-differential interferon responses.
- **Proximal-distal programmes.** Oxidative phosphorylation, MHC-I antigen
  presentation, IL-1 and TNF signalling and planar cell polarity fall toward
  the distal airway in secretory, multiciliated and basal cells; KRT8 and
  KRT19 (cornification and keratinisation) fall distally in multiciliated and
  secretory cells; NKX2-1, NFIB, GATA6, BMP4, SOX9 rise distally in
  multiciliated cells. KRT8 here is an airway positional marker, not an injury
  marker; that distinction has to travel with any Krt8 claim from GSE262927.
- **BMI.** AT2: lower cellular respiration and differentiation programmes;
  secretory: lower insulin-response pathway; alveolar macrophages: higher
  JAK/STAT-linked inflammatory programmes; plasma cells: lower immune
  response, higher respiration, cell cycle and DNA repair.
- **GWAS.** FVC variants enrich in smooth muscle, alveolar, peribronchial and
  myofibroblast genes (adjusted P 0.07 each); asthma in T cells (0.005); COPD
  in myofibroblasts (0.04); lung adenocarcinoma trends to AT2 (0.18). IPF was
  underpowered; depression served as the negative control.
- **Bulk deconvolution in severe COPD (GOLD 3 to 4, n 83 vs 281 controls).**
  Capillary endothelium up (adjusted P 0.0004); alveolar macrophages,
  interstitial macrophages, AT2 and dendritic cells down; smooth muscle up
  (P 1.85e-6). No compositional change with asthma or inhaled corticosteroids
  in the two airway bulk sets.

---

## 9. Software versions stated

Scanpy 1.9.1; scvi-tools 0.8.1 (scANVI); scArches 0.3.5; scIB 0.1.1;
scikit-learn 0.24.1; R 4.1.1 (covariate models) and 4.0.3 (GSEA); edgeR
3.28.1; lme4 1.1-27.1; limma 3.46.0; LDSC 1.0.1; CIBERSORTx 1.0.

---

## 10. Status of the claims in this note

| Claim | Status here |
|---|---|
| Everything in sections 1 to 9 as a description of what the paper did | Descriptive only (transcribed from the paper; owner review pending) |
| "scANVI is the right integration method for this repository's data" | Not established (no benchmark run on this data) |
| "The entropy thresholds 0.56 and 0.43 apply unchanged to GSE262927 and GSE178360" | Not established (donor-entropy threshold depends on donor count; see `ANALYSIS_TRIAL_PLAN.md`) |
| "GSE178360 can be mapped to the HLCA core without circularity" | Descriptive only (Supplementary Table 1 lists it as an extension dataset, not core) |
| "HLCA labels can be transferred to the mouse series" | Not established, and the authors caution against it |
| "AT0 in GSE178360 can be identified by transferring the HLCA marker sets" | Not established (trial S3: the flat and hierarchical schemes disagree; see `ANALYSIS_TRIAL_PLAN.md`) |
| "AT0 in GSE178360 by HLCA reference mapping" | Descriptive only (trial S2): a minority of 119 confident cells; the AT0 candidate analogue subcluster is mostly AT2 or uncertain; our mapping matches the HLCA authors' own transfer of the same cells |
| "GSE178360 maps well to the HLCA core" | Descriptive only (trial S2): unknown fraction at 0.3 between 12 and 24% per donor, inside the paper's healthy range |
| "Mouse cluster 23 is a cell population" | Descriptive only: it is a low-count, ambient-like barcode set (trial S4) |
| "The label disagreement in mouse cluster 5 is a resolution artefact" | Descriptive only (trial S5: true for 94% of labelled cells at Leiden 0.5, not at 0.2) |
| "Human cluster 22 is mast cells" | Descriptive only (trials S2 and S3 and the model classifier agree; correction of the blind table pending owner decision) |
