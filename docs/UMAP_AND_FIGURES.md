# UMAP and the standard figures — how they are made and how to read them

For a reader comfortable with lung biology but not with the computational side.
Part of the background set: see also
[`BACKGROUND_FOR_BIOLOGISTS.md`](BACKGROUND_FOR_BIOLOGISTS.md) (batch effects,
Harmony, doublets) and [`ANALYSIS_RATIONALE.md`](ANALYSIS_RATIONALE.md) (what
was decided for this dataset, and why).

---

## 1. From counts to a matrix

A droplet experiment produces a **cell × gene count matrix**: one row per cell
barcode, one column per gene (~20,000–30,000), each entry the number of distinct
mRNA molecules detected. Counting relies on the **UMI** — a short random barcode
attached to each transcript *before* PCR — so all PCR duplicates sharing a UMI
collapse to one molecule. UMIs are what make the numbers quantitative rather
than a reflection of amplification luck.

The problem: **cells are not sequenced to equal depth.** One cell may yield
2,000 UMIs, its neighbour 20,000 — driven by capture efficiency and cell size,
not biology. A raw count of 10 means different things in each. Two steps fix it:

- **CP10K.** Divide every count in a cell by that cell's total, multiply by
  10,000. Every cell now sums to 10,000, so a value means "molecules per 10,000
  captured" — a proportion, comparable between cells. The 10,000 is arbitrary
  bookkeeping to keep numbers readable.
- **log1p**, i.e. log(x + 1). The +1 avoids log(0). The log matters because
  expression spans orders of magnitude: a housekeeping gene at 5,000 and a
  transcription factor at 5 would otherwise mean the housekeeper alone dominates
  every distance calculation. Logging turns *fold-changes* into equal steps — a
  2-fold change is the same distance at 5 or at 5,000, which is how biologists
  actually think — and stops highly expressed genes drowning everything else.

One correction to a common belief: the many zeros are mostly **not** a special
technical failure mode. Droplet UMI data are consistent with ordinary count
sampling, so a zero usually means "not sampled at this depth", not "broken
measurement". That matters when reading feature plots (§7).

## 2. Highly variable genes

Next, ~2,000–3,000 genes are selected and the rest set aside **for the
structural steps only**. Why discard 90%? Most genes are either off everywhere,
or on at similar levels everywhere. They say nothing about *which cell is which*
and only add noise to the distance between cells.

"Variable" does not mean "highest variance", because in count data variance
rises automatically with mean expression. It means **more variable than expected
for how highly it is expressed** — genes are binned by mean expression and
scored against their own expression-level peers.

What you lose: a gene outside the HVG set contributes nothing to the neighbour
graph, PCA, UMAP or clustering. If a rare population is defined by a gene that
missed the cut, it may never separate. **But feature plots, dot plots and marker
testing use the full matrix**, so any gene can still be plotted and tested. This
pipeline deliberately keeps all genes in the object for that reason.

## 3. Scaling and PCA

**Scaling** z-scores each gene across cells so all enter the next step with
equal weight rather than being ranked by absolute magnitude.

**PCA is a rotation.** Picture each cell as a point in 2,500-dimensional space
(one axis per HVG). That cloud is not spherical — it is stretched along
particular directions, because genes move together in programmes (proliferation,
AT1 identity, interferon response). PCA finds those directions of greatest
stretch. PC1 is the axis along which cells differ most, PC2 the next at right
angles, and so on. **Each PC is a weighted recipe over genes**, so a PC is
effectively a data-derived gene module and a cell's PC score is how strongly it
runs that module.

Typically **30–50 PCs** are kept (50 here). Two reasons:

- **Denoising.** Leading PCs capture coordinated multi-gene variation, which is
  what biology looks like. Trailing ones capture uncorrelated per-gene noise,
  which does not concentrate in any direction. Dropping them removes more noise
  than signal.
- **The curse of dimensionality.** Distances in 2,500 dimensions are dominated
  by noise and become nearly uniform; in 30–50 they behave well, and everything
  downstream runs orders of magnitude faster.

This is why **PCA comes before UMAP**: UMAP runs on PC scores, not on genes.

## 4. The neighbour graph

For each cell, find the *k* closest other cells in PC space (k = 15 here) and
draw an edge to each. The result is a **k-nearest-neighbour graph** — a social
network of cells in which an edge means "transcriptionally similar".

This graph is the real object of interest. Both UMAP *and* clustering are
computed from it, which is why they usually agree: they are two readouts of one
structure, **not independent confirmations of each other**.

## 5. UMAP, and its caveats

UMAP arranges that graph in 2D so **cells that were neighbours in
high-dimensional space stay close on the page** — attractive forces along edges,
repulsion between unconnected points, settled by gradient descent, much like a
spring-and-repulsion network diagram.

The caveats are not pedantry. They are the difference between a defensible
figure legend and a false claim.

- **Distance between well-separated clusters is largely meaningless.** UMAP
  optimises *local* neighbourhoods. Once two clusters share no edges, nothing in
  the objective pins down how far apart they sit. Two blobs at opposite corners
  are not "more different" than two adjacent blobs.
- **Cluster size and area are largely meaningless.** UMAP does not preserve
  density. Area does not indicate cell number; spread does not indicate
  heterogeneity.
- **It is stochastic.** Different random seeds give visibly different pictures,
  as do different `n_neighbors`, `min_dist`, and PC counts. This pipeline fixes
  the seed at 0 — but that buys reproducibility, not correctness.
- **It can invent separations and continuities.** The UMAP documentation itself
  warns it "can also create false tears in clusters". A smooth-looking bridge
  between two populations is **not** evidence of a differentiation trajectory.
- **So it is never evidence on its own.** Claims must rest on marker expression,
  quantification, and ideally staining showing the cells exist in tissue.

**The debate, fairly stated.** Chari & Pachter's *The specious art of single-cell
genomics* is the strongest published critique — they embedded real data into a
literal elephant shape and it scored comparably to standard UMAP on common
metrics. Lause, Berens & Kobak replied in the same journal, agreeing that 2D
embeddings distort distances but arguing the chosen metrics were the wrong ones,
and that under neighbourhood-preservation metrics UMAP does retain real
biological structure. The defensible middle ground, and where most groups have
landed: **UMAP is a legitimate visual index of the neighbour graph, and an
illegitimate measuring instrument.**

## 6. Clustering with Leiden

Leiden is **community detection** on the same kNN graph: find groups of cells
more densely connected to each other than to the rest — cliques in the social
network. It runs on the graph, **not** on the UMAP coordinates.

The **resolution** parameter sets coarseness: higher = more, smaller clusters.
It is a continuous dial, not a discovery.

**There is no objectively correct number of clusters.** At resolution 0.3 this
dataset gives 29; at 1.0 it gives 43. Neither is "true" — they are zoom levels,
and biology genuinely is hierarchical (AT2 cells are one type, and also several
states). What you owe a reader is *justification*: coherent, distinct markers,
and stability. See `ANALYSIS_RATIONALE.md` Decision 4 for how the resolution was
chosen here — and why the obvious criterion turned out to be useless.

**Why Leiden replaced Louvain:** Traag et al. proved Louvain can return badly
connected — even internally disconnected — communities, meaning a "cluster"
could be two unrelated groups of cells sharing one label. Leiden guarantees
connected communities. (Note that Niethamer et al. used Louvain at resolution
1.0; resolution values are not comparable across pipelines.)

## 7. Reading the standard figures

**Cluster UMAP.** Colours are cluster *identifiers*, nothing more. No ordering,
no magnitude, no relationship. Cluster 3 is not between clusters 2 and 4.

**Feature plot.** One gene's normalised expression painted onto the UMAP
coordinates. Two things to check. The **colour scale is usually per-plot**, so a
gene expressed at trivial absolute levels can look brilliantly hot. And **a grey
cell is not a negative cell** — at typical depth a moderately expressed
transcript is missed in many cells simply because it was not sampled. Absence of
colour is weak evidence of absence; a coherent *block* of colour is much
stronger evidence.

**Dot plot** — the workhorse, and the most often misread. It has **two
independent channels**:

- **Dot size = the fraction of cells in that cluster with a non-zero count** (breadth)
- **Dot colour = mean expression** (intensity)

You need both, because they dissociate:

- A **small, dark** dot = a few cells expressing very highly — a rare
  subpopulation hiding inside the cluster, or a doublet artefact.
- A **large, pale** dot = nearly every cell expressing modestly — a genuine
  cluster-wide identity gene.

Both can be biologically real, and they mean opposite things.

**"Scaled" expression** means each gene's colour has been z-scored *across
clusters* before plotting. Without it, a gene at 500 CP10K saturates the scale
and a beautifully specific gene at 5 CP10K is invisible. With it, colour answers
"which clusters express this gene relatively most?" — a within-row comparison.
**Scaled colours are therefore not comparable between rows**, and a dark dot can
correspond to a low absolute level. This pipeline writes both versions and says
which is which in the filename (`..._dotplot` vs `..._dotplot_scaled`).

**Violin plot.** The distribution of one gene per cluster. Its advantage over a
bar of means is that it shows *shape* — bimodality (an on/off subpopulation)
versus a single broad peak. Expect a fat blob at zero.

**Stacked composition bar plot.** One bar per sample, split by cluster
proportion. Two warnings: proportions are **compositional** — if one population
expands, all others shrink by arithmetic necessity, whether or not they changed
in absolute number — and proportions from n = 1 or 2 animals carry no error bar.
**For the mouse data here there is a third**: cells were MACS-sorted and
recombined at roughly 85:15 CD45⁻:CD45⁺, so immune-versus-non-immune proportions
describe the sort, not the lung. See `PIPELINE_AS_RUN.md`.

## 8. Marker genes, and why their p-values are not what they look like

**The Wilcoxon rank-sum test** pools one gene's values from two groups of cells,
ranks them all, and asks whether one group systematically occupies the higher
ranks. It uses order, not magnitude, so it is robust to the skewed, zero-heavy
shape of scRNA-seq data.

**"One cluster vs all other cells"** means exactly that: cluster 5 is group A,
every remaining cell is group B. Consequence worth remembering — a gene shared
by cluster 5 and one similar neighbouring cluster may fail to appear, because
the rest of the dataset dilutes it.

**Why the p-values are not valid in the usual sense.** The clusters were
*defined from the same expression data* used to test them. Leiden's job is
literally to find groups that differ; testing whether they differ then asks a
question whose answer was guaranteed by construction. This is **double dipping**
(post-selection inference). Even pure noise, clustered and then tested, yields
floods of p < 10⁻⁵⁰. Compounding it, cells from one animal are
**pseudoreplicates** — treating 5,000 cells from one mouse as n = 5,000 inflates
false positives dramatically. Scanpy's own `rank_genes_groups` documentation now
carries this warning and points users toward pseudobulk.

**Practical rule: treat marker-gene statistics as a ranking device, not a
hypothesis test.** Report effect sizes (log fold-change, percent expressing) and
use p-values to order candidates. Do not write "significantly enriched,
p = 1e-40" for a cluster marker.

**Why pseudobulk for condition comparisons.** To ask "does gene X change in AT2
cells between injured and control", sum raw counts across AT2 cells **within
each animal** to make one bulk-like profile per animal, then compare across
animals with a bulk framework. Now n = number of animals — the real biological
replication. Squair et al. showed across 18 datasets with matched bulk ground
truth that cell-level methods can return hundreds of "DE genes" where no
biological difference exists. **That is why the condition comparison here is
per-sample pseudobulk — and why it is still labelled exploratory**: with two
control samples, unreplicated within group and confounded with sex, no method
makes that a confirmatory result.

---

## Sources

- [Chari & Pachter 2023, *PLoS Comput Biol* — The specious art of single-cell genomics](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011288)
- [Lause, Berens & Kobak 2024, *PLoS Comput Biol* — The art of seeing the elephant in the room](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1012403) *(the published reply)*
- [Kobak & Berens 2019, *Nat Commun* — The art of using t-SNE for single-cell transcriptomics](https://www.nature.com/articles/s41467-019-13056-x)
- [Kobak & Linderman 2021, *Nat Biotechnol* — Initialization is critical for preserving global data structure](https://www.nature.com/articles/s41587-020-00809-z)
- [Becht et al. 2019, *Nat Biotechnol* — Dimensionality reduction for visualizing single-cell data using UMAP](https://www.nature.com/articles/nbt.4314)
- [UMAP docs — How UMAP Works](https://umap-learn.readthedocs.io/en/latest/how_umap_works.html) · [Clustering caveats](https://umap-learn.readthedocs.io/en/latest/clustering.html) · [FAQ](https://umap-learn.readthedocs.io/en/latest/faq.html)
- [Understanding UMAP — Google PAIR interactive explainer](https://pair-code.github.io/understanding-umap/)
- [Traag, Waltman & van Eck 2019, *Sci Rep* — From Louvain to Leiden](https://www.nature.com/articles/s41598-019-41695-z)
- [Kiselev, Andrews & Hemberg 2019, *Nat Rev Genet* — Challenges in unsupervised clustering of single-cell RNA-seq data](https://www.nature.com/articles/s41576-018-0088-9)
- [Gao, Bien & Witten 2024, *JASA* — Selective inference for hierarchical clustering](https://www.tandfonline.com/doi/full/10.1080/01621459.2022.2116331)
- [Neufeld et al. 2024, *Biostatistics* — Inference after latent variable estimation (count splitting)](https://academic.oup.com/biostatistics/article/25/1/270/6893953)
- [Squair et al. 2021, *Nat Commun* — Confronting false discoveries in single-cell differential expression](https://www.nature.com/articles/s41467-021-25960-2)
- [Zimmerman, Espeland & Langefeld 2021, *Nat Commun* — A practical solution to pseudoreplication bias](https://www.nature.com/articles/s41467-021-21038-1)
- [Crowell et al. 2020, *Nat Commun* — muscat](https://www.nature.com/articles/s41467-020-19894-4)
- [Svensson 2020, *Nat Biotechnol* — Droplet scRNA-seq is not zero-inflated](https://www.nature.com/articles/s41587-019-0379-5)
- [Single-cell best practices](https://www.sc-best-practices.org/) — [normalization](https://www.sc-best-practices.org/preprocessing_visualization/normalization.html) · [feature selection](https://www.sc-best-practices.org/preprocessing_visualization/feature_selection.html) · [dimensionality reduction](https://www.sc-best-practices.org/preprocessing_visualization/dimensionality_reduction.html) · [differential expression](https://www.sc-best-practices.org/conditions/differential_gene_expression.html)
- [scanpy — `rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.rank_genes_groups.html) *(carries the inflated-p-value warning)* · [`dotplot`](https://scanpy.readthedocs.io/en/stable/api/scanpy.pl.dotplot.html)
