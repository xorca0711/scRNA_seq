# Doublets and Scrublet — what they are, how detection works, and how it can delete real cells

For a reader comfortable with lung biology but not with the computational side.
Part of the background set: see also
[`BACKGROUND_FOR_BIOLOGISTS.md`](BACKGROUND_FOR_BIOLOGISTS.md) (batch effects
and Harmony), [`UMAP_AND_FIGURES.md`](UMAP_AND_FIGURES.md), and
[`ANALYSIS_RATIONALE.md`](ANALYSIS_RATIONALE.md) (what was decided here).

**This section matters more than the others for this project**, because the
measured bias described in §3 applies directly to the human dataset analysed in
this repository.

---

## 1. What a doublet physically is

In a 10x run, cells are squeezed one at a time into oil droplets (GEMs), each
containing a bead carrying a unique DNA barcode. Every mRNA captured in that
droplet is stamped with the same barcode, and downstream the software assumes
**one barcode = one cell**.

A **doublet** is when two cells land in the same droplet. Their transcripts get
the same stamp, so the matrix reports the *sum* of two transcriptomes as though
it were one cell. Nothing in the sequencing flags this; the barcode looks
completely ordinary.

Two flavours, and only one is really detectable:

- **Homotypic** — two cells of the same type (AT2 + AT2). The sum still looks
  like an AT2 cell, just with roughly double the UMIs. Nearly impossible to
  detect from transcriptome alone, and mostly harmless: they sit with their
  parent cells and create no fake population.
- **Heterotypic** (or *neotypic*) — two different cell types (AT2 + macrophage).
  The sum expresses `Sftpc` **and** `Cd68`. On a UMAP this often lands in the
  gap between the two parent clusters and **looks exactly like a novel hybrid or
  intermediate cell state.** This is the damaging kind, and the kind algorithms
  can find.

**Why the rate scales with loading.** Encapsulation is essentially random
(Poisson), so loading more cells into the same number of droplets raises the
chance any droplet catches two. 10x's rule of thumb is **~0.8% multiplets per
1,000 cells recovered** (≈8% at 10,000 cells).

Treat that number as a *prior*, not a measurement. It is a manufacturer estimate
from species-mixing ("barnyard") experiments, it is a linear approximation, and
it is **chemistry-specific** — GEM-X chemistry roughly halves it to ~0.4% per
1,000. Species mixing only directly observes *cross-species* doublets, so the
total is inferred. It is a reasonable expectation-setter for your lane; it is
not a property of your particular lung digest.

> This pipeline uses that 0.8%/1,000 figure as the prior against which
> Scrublet's automatic threshold is sanity-checked. See §4.

## 2. How Scrublet works

Scrublet (Wolock, Lopez & Klein 2019) is simple in concept:

1. **Simulate doublets from your own data.** Pick two observed cells at random
   and add their counts gene by gene: `counts(fake) = counts(A) + counts(B)`.
   Repeat tens of thousands of times. You now have synthetic doublets built from
   exactly the cell types present in *your* sample.
2. **Put real and simulated cells in one space** — normalise, take variable
   genes, run PCA on the combined set.
3. **Score each real cell by its neighbourhood.** Build a k-nearest-neighbour
   graph and ask, for each real cell: *what fraction of my nearest neighbours
   are simulated doublets?* A real cell sitting in a region packed with
   simulated doublets is probably one itself. That fraction, rescaled by the
   expected doublet rate, is the **doublet score** (0–1).

**Choosing the threshold — and why bimodality matters.** Scrublet plots a
histogram of the scores *of the simulated doublets*. It is usually **bimodal**,
and the two peaks are meaningful:

- a **low-score peak** = "embedded" doublets — pairs of similar cells whose sum
  lands back inside a real cluster. Undetectable (essentially the homotypic
  case).
- a **high-score peak** = "neotypic" doublets — pairs of dissimilar cells whose
  sum lands where no real cell lives. Detectable.

The threshold goes in the **valley between the two peaks**, and is then applied
to the real cells.

**When the histogram is not bimodal, the automatic threshold is arbitrary.** The
implementation simply takes the minimum between two modes; with no clean valley
it still returns a number, and that number means little. The authors suggest two
self-consistency checks: is the histogram genuinely bimodal, and does the
fraction called ≈ (expected rate) × (detectable fraction)?

> **This is exactly what happened here.** In **32 of 33** mouse samples the
> automatic threshold produced a call rate irreconcilable with the 10x prior, so
> the pipeline substituted the expected-rate quantile and recorded the
> substitution per sample in `qc/doublet_summary.csv`. The consequence, stated
> plainly: **the doublet calls in this analysis are a ranking cut, not a
> detected threshold.**

One structural limitation from the paper: the method needs every cell state
contributing to a doublet to *also be present as a singlet*. If you sorted away
a parent population, its doublets are invisible.

## 3. The critical failure mode — real intermediate cells look like doublets

This is the part that matters most for lung biology, and it is why this section
exists.

Scrublet and its relatives do not detect "two cells in a droplet". They detect
**transcriptomes that look like a blend of two other transcriptomes in the
dataset**. A genuine cell co-expressing two lineage programmes is, to the
algorithm, indistinguishable from an artefactual sum.

The Bioconductor OSCA book states it directly: *"Cells in the middle of a
trajectory are always intermediate between other cells and are liable to be
incorrectly detected as doublets."*

Lung populations at live risk:

- **Krt8⁺ ADI / DATP / PATS** transitional cells bridging AT2 → AT1 — by
  construction they express residual AT2 and emerging AT1 genes.
- **AT0** cells co-expressing `SFTPC` (alveolar) and `SCGB3A2` (airway), the
  bipotent progenitor of Kadur et al. 2022. An `SFTPC⁺SCGB3A2⁺` barcode is
  *precisely* what an AT2 + club-cell doublet looks like.
- **AT1/AT2 double-positive** cells generally.
- **Aberrant basaloid** cells in IPF, mixing basal, mesenchymal and senescence
  programmes.

### Measured on this dataset

Because GSE178360 is the Kadur dataset, this was quantified directly rather than
left as a worry. Cells were classified by the paper's own marker definitions and
cross-tabulated against Scrublet's calls:

| Population | cells | called doublet | vs. baseline |
|---|---:|---:|---:|
| distal-BC (`TP63⁺SFTPB⁺`) | 955 | 12.9% | **2.03×** |
| SCGB3A2-CC (`FOXJ1⁺SCGB3A2⁺SFTPB⁺`) | 1,559 | 10.5% | **1.66×** |
| AT0 (`SFTPC⁺SCGB3A2⁺`) | 5,680 | 8.0% | **1.26×** |
| all cells | 29,605 | 6.3% | 1.00× |

Scrublet removed these populations at up to **twice** the background rate. It
did **not** erase them — 87–92% survived — but the bias is systematic and runs
in the direction that loses the study's discoveries. **Kadur et al. used no
automated doublet caller at all**, relying on manual cluster inspection. That
now reads as a deliberate choice rather than an omission.

*(Caveat: those marker definitions are crude — any non-zero expression — so they
are generous upper-bound sets. The enrichment ratio is the meaningful number,
not the absolute counts.)*

### How the field handles this

1. **Orthogonal validation is the standard of proof.** Strunz et al. confirmed
   Krt8⁺ ADI by immunofluorescence, flow cytometry, mass spectrometry and
   lineage tracing — not transcriptomics alone. AT0 was backed by smFISH and
   spatial transcriptomics. *If a "novel intermediate" exists only in your count
   matrix, it is not yet a cell type.*
2. **Counter-model the specific doublet.** Strunz et al. defended an MHC-II⁺
   club-cell state by showing it did **not** resemble artificially generated
   club + dendritic-cell doublets — they simulated the exact doublet that would
   explain their finding away, and showed it did not match. This is a cheap,
   powerful check.
3. **It is a risk, not an inevitability.** DoubletFinder's authors report their
   method classified known hybrid states as singlets, and OSCA notes real
   trajectories are often non-linear enough that simulated doublets miss the
   real intermediates.
4. **The biology itself is contested** — a 2023 *JCI* commentary argues some
   "transitional" states may be stable intermediates rather than cells in
   transit.

**An honesty note about the benchmarks:** the published doublet benchmarks do
**not** measure this failure mode. They score how many *true* doublets a method
recovers; nothing in that design detects "you deleted a rare real cell type". A
high AUPRC tells you nothing about whether your ADI or AT0 cluster is safe.

## 4. Why detection must run per capture

A cell from lane A and a cell from lane B were never in the same droplet — never
in the same instrument run. Therefore:

- Summing a lane-A and a lane-B cell fabricates an entity that cannot exist,
  contaminating the simulated reference.
- Batch effects between runs dominate the PCA, so "neighbourhood" starts
  reflecting *which lane* a cell came from rather than whether it is a blend.
- The expected doublet rate is a property of **one capture's loading**, not of a
  merged 200,000-cell object.

Both major references say so explicitly. sc-best-practices: *"Doublet detection
methods should not be run on aggregated scRNA-seq data representing multiple
batches."* OSCA: apply *"only to libraries generated in the same experimental
batch."*

> This pipeline runs Scrublet independently on every capture, before merging.

If one biological sample was split across two lanes, the unit is the **lane**.

## 5. Alternatives, and what the benchmarks actually say

| Tool | Approach |
|---|---|
| **Scrublet** | simulated doublets + kNN fraction (Python) |
| **DoubletFinder** | simulated doublets + kNN; tunes *k* by maximising bimodality (R/Seurat) |
| **scDblFinder** | simulated doublets, multi-scale kNN + gradient-boosted tree, iterative (R) |
| **scds** | `cxds` uses gene **co-expression only, no simulation**; `bcds` boosted-tree; `hybrid` combines |
| **DoubletDetection** | simulated doublets + iterative clustering, hypergeometric test |
| **solo** | scVI latent space + neural-net classifier |

**Xi & Li 2021** benchmarked 9 methods on 16 real datasets with experimental
ground truth. DoubletFinder had the best mean AUPRC; solo the best AUROC; cxds
was fastest. Two findings worth internalising:

- **Absolute accuracy is modest** — the best method reached mean AUPRC 0.537.
  The authors call this "the general difficulty in detecting doublets".
- **Scrublet was less accurate but more stable** across data subsets than
  DoubletFinder. A real trade-off, not a strict ranking.

**scDblFinder** was added in a later protocol paper and achieved the highest
mean AUPRC and AUROC on the same datasets while being among the fastest;
sc-best-practices recommends it on that basis. Treat the leaderboard with normal
caution — benchmark and tool share an ecosystem.

**On combining callers — a correction to a common belief.** It is often said
that stacking doublet callers does not help. Xi & Li found the opposite: the
`hybrid` ensemble "largely improved on its both base methods without much
running time increase". Their recommendation is to pair an accurate method with
an **algorithmically distinct** one (e.g. DoubletFinder + cxds), because the top
methods correlate heavily and stacking similar ones adds little. Demuxafy (2024)
operationalises this, showing consensus across methods improves droplet
assignment. Practical reading: **intersection** (flagged by ≥2 methods) is
conservative and protects rare real states; **union** is aggressive.

> Note that Niethamer et al. ran **both scds and Scrublet** on the mouse data.
> This pipeline ran Scrublet alone.

## 6. Practical guidance

**Score early, delete late.** Run detection per lane, store the score and call
in the metadata, then cluster and embed the *full* object. Removing before you
know what you have means you can never audit what you lost.

**Always inspect what got flagged.** Overlay the doublet score on the UMAP and
cross-tabulate flagged cells against clusters.

*Signs a suspicious cluster really is a doublet artefact:*
- a high **proportion** of its cells flagged, not scattered singletons
- elevated UMI and gene counts (~2× its neighbours)
- expression is essentially the **union** of two parent clusters, with no genes
  of its own
- it sits between two clusters in the embedding
- its abundance matches the expected doublet rate

*Signs it is real biology:*
- it appears at **consistent proportions across independent lanes and donors** —
  doublet rate is a per-lane technical property, so a state present at 2% in
  every sample is unlikely to be artefactual
- it has **its own markers** not explained by either parent (e.g. `Krt8`,
  `Cldn4`, `Krt19` for ADI)
- normal UMI counts
- protein, in situ, spatial or lineage-tracing support

**Don't blindly take the top N%** — look at the simulated-score histogram
(`sc.pl.scrublet_score_distribution`). If it is not bimodal, the automatic
threshold is untrustworthy: set it manually and record what you did.

**Do a sensitivity check.** If a conclusion depends on a borderline population,
re-run the key analysis with and without doublet removal and confirm the finding
survives.

**Report it.** Per lane: tool and version, expected rate used, threshold, number
and percentage removed. This pipeline writes all of that to
`qc/doublet_summary.csv`.

---

## Sources

- [Wolock, Lopez & Klein 2019, *Cell Systems* — Scrublet](https://doi.org/10.1016/j.cels.2018.11.005) ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/))
- [Scrublet source and usage notes](https://github.com/swolock/scrublet)
- [`scanpy.pp.scrublet` docs](https://scanpy.readthedocs.io/en/stable/api/generated/scanpy.pp.scrublet.html)
- [McGinnis, Murrow & Gartner 2019, *Cell Systems* — DoubletFinder](https://doi.org/10.1016/j.cels.2019.03.003)
- [Xi & Li 2021, *Cell Systems* — Benchmarking computational doublet-detection methods](https://doi.org/10.1016/j.cels.2020.11.008) ([PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC7897250/))
- [Xi & Li 2021, *STAR Protocols* — protocol extension adding scDblFinder](https://doi.org/10.1016/j.xpro.2021.100699)
- [Germain et al. 2021, *F1000Research* — scDblFinder](https://pmc.ncbi.nlm.nih.gov/articles/PMC9204188/)
- [Neavin et al. 2024, *Genome Biology* — Demuxafy](https://doi.org/10.1186/s13059-024-03224-8)
- [Single-cell best practices — Quality control](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)
- [OSCA — Doublet detection](https://bioconductor.org/books/3.15/OSCA.advanced/doublet-detection.html)
- [10x Genomics — Chromium Next GEM Single Cell 3′ v3.1 User Guide (multiplet rate table)](https://cdn.10xgenomics.com/image/upload/v1668017706/support-documents/CG000315_ChromiumNextGEMSingleCell3-_GeneExpression_v3.1_DualIndex__RevE.pdf)
- [10x Genomics — GEM-X technology (~0.4% per 1,000 cells)](https://www.10xgenomics.com/blog/the-next-generation-of-single-cell-rna-seq-an-introduction-to-gem-x-technology)
- [Strunz et al. 2020, *Nat Commun* — Krt8⁺ transitional stem cell state](https://doi.org/10.1038/s41467-020-17358-3)
- [Choi et al. 2020, *Cell Stem Cell* — DATP](https://doi.org/10.1016/j.stem.2020.06.020)
- [Kobayashi et al. 2020, *Nat Cell Biol* — PATS](https://doi.org/10.1038/s41556-020-0542-8)
- [Kadur Lakshminarasimha Murthy et al. 2022, *Nature* — AT0 bipotent progenitor](https://doi.org/10.1038/s41586-022-04541-3)
- ["Stuck in the Middle with You": intermediate cell states are not always in transition, *JCI* 2023](https://doi.org/10.1172/JCI174633)
