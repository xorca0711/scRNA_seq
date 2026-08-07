# Background: how these figures are made, and why the decisions were made

Written for a reader comfortable with lung biology and the bench, but not with
the computational side. Intuition first, maths only where it earns its place.

Companion to [`ANALYSIS_RATIONALE.md`](ANALYSIS_RATIONALE.md), which records
what was decided for *this* dataset. This file explains the concepts behind
those decisions. Sources are listed at the end of each section.

**Contents**

1. [Batch effects and Harmony](#1-batch-effects-and-harmony) — why samples separate, what correction does, and when it destroys your experiment
2. [Doublets and Scrublet](#2-doublets-and-scrublet) — what a doublet is, how it is detected, and why the method can delete real cell types
3. **[UMAP and the standard figures](UMAP_AND_FIGURES.md)** — separate file: how
   the picture is built from counts, what it does and does not mean, how to read
   dot plots and feature plots, and why cluster-marker p-values are a ranking
   device rather than a hypothesis test

---

# 1. Batch effects and Harmony

## 1.1 What a batch effect actually is

Every step between "cells in a tube" and "counts in a matrix" leaves a
fingerprint. A **batch effect** is expression variation that comes from *how* a
sample was processed rather than *what* the cells were doing.

Physical causes:

- **Different days, different hands.** Dissociation time, enzyme lot, time on
  ice. Longer or warmer dissociation induces a stress program (*Fos*, *Jun*,
  heat-shock genes) and preferentially kills fragile cells.
- **Different 10x runs.** Chemistry version (v2 vs v3 recover very different
  gene counts per cell), lane, bead lot, ambient RNA load, sequencing depth.
- **Different donors or animals.** Genotype, age, sex, agonal state. This is
  *biological* variation, but still nuisance variation if your question is
  "what cell types are here".
- **Different reference annotations.** If sample A was aligned against one
  GENCODE build and sample B against another, genes gain or lose exons, are
  renamed, or vanish. Whole genes can look "differentially expressed" purely
  from annotation. **This happened in GSE178360** — sample DD073R used a
  different GRCh38 build (36,601 genes vs 33,538).

**Why it shows up on a UMAP.** UMAP does not look at genes; it looks at the top
principal components. A batch effect nudges *every* cell in a sample slightly in
the same direction. That shared nudge is exactly the kind of large, coordinated
pattern PCA is built to find, so it lands on an early PC — sometimes PC1 — and
UMAP faithfully renders it as physical separation.

The signature to look for is **parallel structure**: instead of one AT2 blob,
one macrophage blob, one endothelial blob, you see *the same set of shapes
repeated once per sample*, side by side. A single cluster containing only one
sample is more ambiguous — it could be a real sample-specific population.

## 1.2 What Harmony does

Harmony (Korsunsky et al. 2019) is deliberately narrow, and that is its
strength.

**It works on PCA coordinates, not on genes.** You hand it the cell-by-PC matrix
(20–50 PCs) plus a batch label per cell. It never touches expression values.

It then loops until stable:

1. **Soft-cluster the cells.** Every cell gets fractional membership — 60%
   cluster 3, 30% cluster 7, 10% cluster 1 — rather than a hard assignment. This
   matters because real tissue is continuous; an AT2-to-AT1 intermediate should
   not be forced into one bin.
2. **Apply a diversity penalty.** This is what "maximum diversity clustering"
   means. Ordinary clustering is perfectly happy to build a cluster containing
   only sample B. Harmony adds a penalty making such clusters expensive, so it
   prefers clusters whose sample composition resembles the dataset overall.
   *Analogy:* seating a dinner party where each table is scored both on "do
   these people have something in common" **and** "is every family represented".
3. **Learn a correction per cluster, per batch** — how far each batch's centroid
   sits from the cluster centroid, via a ridge-regularised linear model.
4. **Move each cell** by a *weighted blend* of the corrections from the clusters
   it partly belongs to. Correction is therefore cell-specific and smooth.

**Output: corrected coordinates only.** In scanpy this lands in
`adata.obsm['X_pca_harmony']`. There are no genes in it.

**Why that matters practically.** You may build the neighbour graph, cluster and
draw UMAPs on the Harmony embedding. You must **not** run differential
expression on it. The Bioconductor OSCA book states that a correction algorithm
"is not obliged to preserve relative differences in per-gene expression" and
that correction "will inevitably introduce artificial agreement across batches".
Harmony's own authors recommend mixed-effects models on the original data
instead. In practice: **cluster on corrected coordinates, test genes on
uncorrected counts** — which is what this pipeline does (markers come from
log-normalised values, and the condition comparison uses per-sample pseudobulk).

## 1.3 The central danger: over-correction

**A correction algorithm cannot tell technical difference from biological
difference.** It only sees the label you hand it. Say "sample is the batch" and
it will make samples indistinguishable — whether or not they genuinely differ.

The failure case is **confounding**. If sample 1 = day 0, sample 2 = day 7,
sample 3 = day 21, then "sample" and "timepoint" are the same variable.
Correcting on sample is mathematically identical to deleting the timepoint
effect. You get a beautiful, well-mixed UMAP that has erased the regeneration
biology you ran the experiment to find — **and nothing in the output warns you.**

> **This is exactly why GSE262927 was not corrected.** Every one of its 33
> samples belongs to exactly one experimental group, so `sample_id` *is* the
> treatment variable. See `ANALYSIS_RATIONALE.md`, Decision 1.

Benchmarks make the trade-off explicit rather than solving it. Luecken et al.
(2022) scored 16 methods on 13 tasks using two *opposing* metric families —
batch removal versus biological conservation — precisely because you can always
buy more of one with less of the other. They note that "retaining batch effects
in a dataset to preserve all nuanced biological variation may be preferable".
Maan et al. (2024) add that when cell-type *proportions* differ between samples,
rare populations are damaged disproportionately.

## 1.4 When should you integrate?

| Situation | Integrate? |
|---|---|
| Biological replicates of one condition, run on different days | **Yes** — textbook |
| Multiple donors; question is "what cell types exist" | **Yes**, batch = donor |
| Your data + a public atlas, different chemistries | **Yes**, and expect to need it |
| Each sample is a different timepoint/treatment, one sample per condition | **No** — this destroys the signal |
| Replicates nested inside conditions (3 mice × 3 timepoints) | Batch = **mouse/run**, never = timepoint |
| Samples already overlap well | **No** — nothing to fix |

The design principle: **you need replication at the level below your variable of
interest.** With ≥2 independent samples per condition you can correct on
sample-within-condition and still measure the condition. With one sample per
condition you cannot, and no algorithm rescues that.

> **GSE178360 is row 2** — three healthy donors, question is "what cell types
> exist in distal airway". Integration is appropriate, and the original authors
> did it too (Seurat CCA anchors).

Always **look at the uncorrected UMAP first** and keep it. This pipeline retains
both embeddings and writes a before/after figure.

## 1.5 How do you tell whether you need it?

The core idea is a **k-nearest-neighbour mixing check**. For each cell, take its
*k* nearest neighbours and ask what fraction come from the same sample. Compare
to the null: if sample *s* holds 12% of all cells, you expect ~12% of any cell's
neighbours to be from *s* under perfect mixing.

Named versions of the same idea:

- **kBET** — χ² test comparing local batch composition to global, reported as a
  rejection rate. Low = well mixed.
- **LISI / iLISI / cLISI** — "how many different labels do I effectively see
  among my neighbours". iLISI on *batch* should be **high**; cLISI on *cell
  type* should be **low**. One number for mixing, one for damage.
- **Silhouette (ASW)** — on *batch*, near zero is good; on *cell type*, high is
  good.

**The interpretation trap, and it caught this analysis.** With many samples, a
big-sounding enrichment can still mean good mixing, because chance level is low.
With 33 equal samples chance is ~3%; observing 12.7% is a "4× enrichment" —
alarming-sounding — yet **87% of every cell's neighbours still come from other
samples**. That is good mixing. Conversely with 2 samples chance is 50%, so 85%
same-sample is only "1.7× enriched" but means near-total separation.

> **Always quote the absolute fraction alongside the baseline, never the ratio
> alone.** GSE262927: 12.7% observed vs 3.2% expected → mixes well.
> GSE178360: 88.6% observed vs 33.7% expected → badly separated.

## 1.6 Alternatives to Harmony

- **Seurat CCA anchors** — mutual nearest-neighbour "anchor" pairs in a shared
  correlation space; outputs a corrected *expression* matrix. Powerful, but
  Seurat's own docs warn it "may also lead to overcorrection". *Used by Kadur et
  al. on this very dataset.*
- **Seurat RPCA** — same framework, projecting each dataset into the other's PC
  space. Faster and explicitly "more conservative".
- **scVI / scANVI** — variational autoencoders modelling raw counts with batch
  as a covariate. Strong on complex tasks; scANVI additionally uses cell-type
  labels, so it is told what to protect — but is only as good as those labels.
- **BBKNN** — skips coordinate correction, building the neighbour graph from
  each cell's top neighbours *within each batch*. Very fast, leans aggressive.
- **Scanorama**, **fastMNN** — mutual-nearest-neighbour families; middle of the
  road in benchmarks.
- **Not correcting** — a real option, and the right one under confounding.

**Where the field disagrees.** Rankings differ between benchmarks: Tran et al.
(2020) put Harmony, LIGER and Seurat 3 on top, while Luecken et al. (2022) found
Harmony competitive on *simple* tasks but scANVI/Scanorama/scVI better on
complex ones — largely because the two studies weighted mixing versus
conservation differently. The practical conclusion is not "method X wins" but
**run two, measure both mixing and cell-type preservation, and choose the
weakest correction that fixes your problem.**

### Sources

- [Korsunsky et al. 2019, *Nature Methods* — Harmony](https://www.nature.com/articles/s41592-019-0619-0) ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6884693/))
- [Luecken et al. 2022, *Nature Methods* — Benchmarking atlas-level data integration](https://www.nature.com/articles/s41592-021-01336-8) ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC8748196/))
- [Tran et al. 2020, *Genome Biology* — A benchmark of batch-effect correction methods](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1850-9)
- [Büttner et al. 2019, *Nature Methods* — kBET](https://www.nature.com/articles/s41592-018-0254-1)
- [Argelaguet et al. 2021, *Nature Biotechnology* — Computational principles and challenges in single-cell data integration](https://www.nature.com/articles/s41587-021-00895-7)
- [Maan et al. 2024, *Nature Biotechnology* — Impacts of dataset imbalance on integration](https://www.nature.com/articles/s41587-023-02097-9)
- [Single-cell Best Practices — Integration](https://www.sc-best-practices.org/cellular_structure/integration.html)
- [OSCA — Using the corrected values](https://bioconductor.org/books/release/OSCA.multisample/using-corrected-values.html)
- [scanpy — `harmony_integrate`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html)
- [Seurat — Fast integration using RPCA](https://satijalab.org/seurat/articles/integration_rpca)
- [Stuart et al. 2019, *Cell* — Seurat v3 anchors](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)

---

# 2. Doublets and Scrublet

*(pending — research in progress)*

---

# 3. UMAP and the standard figures

Moved to its own file: **[`UMAP_AND_FIGURES.md`](UMAP_AND_FIGURES.md)**.
