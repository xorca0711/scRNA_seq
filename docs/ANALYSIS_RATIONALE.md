# Analysis rationale — what I decided, why, and what changed after reading the papers

This document exists because the analysis was run in two distinct passes, and
the second pass changed some of the reasoning behind the first.

**Pass 1 — data only.** The brief was explicit: *"Do not begin by assuming what
kind of scRNA-seq data are present… let the structure and properties of the raw
data determine the workflow."* So the pipeline was built by reading the
deposited files and nothing else. No papers, no prior knowledge of the studies.

**Pass 2 — after reading the source publications.** The two papers behind these
accessions live in `Thesis/Primary/`. They were read only after the analysis was
complete. Some decisions were confirmed. One was right for the wrong reason.
Several facts turned out to be unknowable from the deposited files alone.

Both passes are recorded below, honestly, including where pass 1 was wrong.
For the parameters actually used, see [`PIPELINE_AS_RUN.md`](PIPELINE_AS_RUN.md).

---

## Part 1 — The decisions, before and after

### Decision 1: Should the samples be batch-corrected? (GSE262927, mouse)

**Pass 1 reasoning.** Thirty-three samples. The obvious move is to "integrate"
them. But look at what `sample_id` actually is here: every sample belongs to
exactly one `experimental_group` (`H1N1_tam02_sac06`, `H1N1_homeostasis_d03`,
and so on). Sample identity and experimental condition are the *same variable*.
Correcting on sample would therefore delete the influenza time course — the
entire point of the study.

So the question was reframed. Not *"do samples separate?"* — they might separate
for real biological reasons — but *"do samples that share an experimental group,
i.e. genuine replicate animals, fail to mix?"* That is a question about
technical noise only, because replicates within a group differ by animal, not by
treatment.

Measured: within-group replicate enrichment **1.374** (a value of 1.0 means
replicates mix as well as random). Globally, 12.7% of each cell's neighbours
came from its own sample against 3.2% expected by chance. Verdict: **no
correction**.

**Pass 2 verdict: CONFIRMED.** Niethamer et al. performed **no batch correction
and no integration of any kind**. Their Methods contain no Harmony, no anchors,
no CCA — libraries were merged and clustered directly. The only thing regressed
out is technical (mitochondrial %, nFeature, nCount) inside SCTransform. They
also substituted a compositional sanity check across timepoints in place of
correction, which is essentially what the mixing assessment did.

**What pass 2 added.** There *is* a legitimate technical batch variable, and it
is not in the deposited metadata: `rep1`/`rep2` are influenza infection waves
about nine months apart, and they are balanced across the acute timepoints — so
they could be tested or corrected without harming the time signal. That variable
lives in the paper's supplementary table `mmc4.xlsx`. The 90 and 366 dpi
timepoints come exclusively from a third wave, where wave and timepoint are
perfectly confounded and correction is impossible.

---

### Decision 2: Should the samples be batch-corrected? (GSE178360, human)

This one went the other way, and it is the more interesting story.

**Pass 1, first attempt — declined.** The series ships **no metadata file at
all**. The reasoning recorded at the time: *"no condition metadata exists, so
donor, batch and biology are completely confounded and no correction can be
justified as purely technical."* Harmony was computed but kept as a clearly
labelled *supplementary* embedding; the uncorrected version stayed primary.

**Pass 1, second attempt — reversed, on evidence.** The QC review then flagged
that most clusters were nearly private to one donor. Rather than trust an
impression, the specific question was tested: *are single cell types fragmenting
along donor lines?*

| Cell type | split across | one-donor purity of each |
|---|---|---|
| Fibroblast | 4 clusters | 97% / 54% / 81% / 100% |
| Ciliated | 4 clusters | 97% / 86% / 46% / 100% |
| Mono-macrophage | 3 clusters | 54% / 95% / 99% |

Three human lungs do not contain four kinds of fibroblast. That is donor offset
splitting one cell type, not biology. Harmony was made primary, and this test
was written into the pipeline as an automated check
(`donor_driven_clustering_check`). Afterwards: donor-private clusters fell from
20 of 41 to **4 of 31**, and the check reports `donor_driven: false`.

**Pass 2 verdict: RIGHT ANSWER, UNSOUND REASONING.**

The answer holds. Kadur et al. describe three healthy-donor biological
replicates of one anatomical preparation (microdissected distal airways under
1 mm), and **the authors themselves integrated across the same three donors**
using Seurat CCA anchors with 10,000 features. Harmony substitutes for a
correction the original study also performed.

But the *stated reason* was a non-sequitur. "No metadata file exists, therefore
no designed contrast exists" does not follow — missing metadata files are
routine even for designed experiments. The safe conclusion required checking the
source, and it happens to hold for a reason pass 1 could not have known: the
paper's scientific axis is **proximal-to-distal zonation, a gradient *within*
each sample**, which sample-level correction cannot touch. Had the three samples
been three anatomical regions, correcting on sample would have destroyed the
finding, and the recorded justification would not have caught it.

There is also a *better* argument for correcting than the one given: DD073R was
quantified against a different GRCh38 annotation build (36,601 genes vs 33,538),
which is a genuine technical artefact perfectly confounded with sample.
Correction actively helps there.

**Lesson.** Absence of metadata is not evidence of absence of design. Check the
source publication before concluding that nothing is at risk.

---

### Decision 3: How aggressively to remove doublets

**Pass 1 reasoning.** Scrublet, run independently per 10x capture — doublets can
only form between cells that shared a droplet, so simulating them across samples
is meaningless. When Scrublet's automatic threshold disagreed wildly with the
10x multiplet-rate prior (~0.8% per 1,000 cells recovered), it was replaced by
that expected-rate quantile, and the substitution recorded per sample. This
happened in **32 of 33** mouse samples, because the simulated-doublet score
histograms were not bimodal.

**Pass 2 verdict: A REAL LIMITATION, now measured.**

Kadur et al. used **no automated doublet caller at all** — only manual
inspection of clusters. That is not sloppiness. Every novel population in that
paper is defined by co-expression of two lineages' markers:

- **AT0** (the bipotent progenitor of the title) = `SFTPC⁺` (alveolar) + `SCGB3A2⁺` (airway)
- **SCGB3A2-CC** = `FOXJ1⁺` (ciliated) + `SCGB3A2⁺` + `SFTPB⁺` (secretory)
- **distal-BC** = `TP63⁺` (basal) + `SFTPB⁺` (alveolar)

That is precisely the signature a doublet caller is built to remove. So this was
measured directly on the data:

| Population | cells | called doublet | vs. baseline |
|---|---:|---:|---:|
| distal-BC (`TP63⁺SFTPB⁺`) | 955 | 12.9% | **2.03×** |
| SCGB3A2-CC (`FOXJ1⁺SCGB3A2⁺SFTPB⁺`) | 1,559 | 10.5% | **1.66×** |
| AT0 (`SFTPC⁺SCGB3A2⁺`) | 5,680 | 8.0% | **1.26×** |
| all cells | 29,605 | 6.3% | 1.00× |

Scrublet removed these populations at up to twice the background rate. It did
**not** erase them — 87–92% survived, so they remain in the object — but the
bias is systematic and runs in the direction that loses the study's discoveries.

*(Caveat: these marker definitions are crude — any non-zero expression of each
gene — so they are generous upper-bound sets, not the real populations. The
enrichment is the meaningful number, not the absolute counts.)*

### Follow-up: was the AT0 population actually lost? **No.**

The enrichment above was measured with a deliberately crude gate - *any*
non-zero `SFTPC` and `SCGB3A2`. That gate does not distinguish a real AT0 cell
from an actual AT2 + club-cell doublet, which is precisely the thing it was
being used to worry about. Re-measured with a definition that can tell them
apart - `SFTPC+ SCGB3A2+ EPCAM+` **and** negative for `PTPRC`, `PECAM1` and
`COL1A1` - the picture reverses:

| AT0 definition | cells | called doublet | vs. 6.3% baseline |
|---|---:|---:|---:|
| crude (`SFTPC+ SCGB3A2+`) | 5,680 | 8.0% | 1.26x |
| **strict (+EPCAM+, lineage-negative)** | **1,676** | **3.9%** | **0.62x - below baseline** |

And the cells Scrublet flagged inside the crude gate look exactly like
doublets rather than like AT0:

| | flagged | kept |
|---|---:|---:|
| median UMI | 13,705 | 7,691 |
| median genes | 3,860 | 2,610 |
| co-expressing another lineage | 73.5% | 53.3% |

Nearly twice the library size and substantially more cross-lineage
co-expression. **Scrublet was mostly right.** The strict AT0 population is
*under*-flagged relative to background, i.e. it survived the filter intact.

**The correction that matters, and it is a general one:** the apparent 2x bias
was an artefact of the marker gate used to measure it, not of the doublet
caller. A co-expression gate cannot be used to audit a tool whose whole job is
detecting co-expression - you have to add a discriminating criterion (here,
lineage-negativity and library size) before the question is answerable. The
concern was worth testing; the first test was the wrong test.

The general guidance in the sections above still stands - automated doublet
callers *can* remove real intermediate states, the original authors did use
manual curation, and flagged cells should be inspected before deletion. What
changed is the verdict for this specific dataset.

**What should be done differently.** Anyone pursuing AT0 specifically should
re-run without doublet filtering, or inspect the flagged cells against these
marker combinations before deleting them. Skipping ambient-RNA correction
compounds the problem: surfactant transcripts are classic ambient contaminants,
and AT0 is defined partly by `SFTPC`, so some `SFTPC⁺SCGB3A2⁺` cells may be
airway cells wearing borrowed surfactant mRNA. The original authors ran SoupX
before normalisation precisely for this reason.

---

### Decision 4: Which clustering resolution

**Pass 1 reasoning.** Four resolutions were scanned. The obvious criterion —
"does every cluster have marker genes distinguishing it from the rest of the
data?" — turned out to be useless: **all four resolutions scored 100%**. Judged
against all other cells, even 43 clusters looked perfectly justified.

So a second criterion was added: is each cluster distinguishable from its
*nearest neighbouring cluster*? That is what over-clustering actually looks
like — splitting one population in half, where each half still differs from the
rest of the data but not from its own sibling.

| resolution | clusters | marker-supported | min DE genes between nearest clusters |
|---:|---:|---:|---:|
| **0.3** | **29** | 100% | **9** |
| 0.5 | 32 | 100% | 8 |
| 0.8 | 40 | 100% | **0** |
| 1.0 | 43 | 100% | **0** |

At resolutions 0.8 and 1.0 there are cluster pairs that **no gene separates**.
Resolution 0.3 was chosen. All resolutions are retained in the object as
`leiden_res*` columns so nothing is lost.

**Pass 2 note.** Niethamer et al. used **Louvain at resolution 1.0** for the full
atlas. That is a coarser-looking number than it sounds, because it follows
SCTransform normalisation and a different neighbour graph — resolution values
are not comparable across pipelines. More importantly, they then **manually
culled clusters** (removing ones without a marker signature, ones with
multi-compartment markers, and deliberately dropping platelets, erythrocytes,
mast cells and basophils). That hand-curation is not reproducible from the
deposited files, which is why cell counts here differ from theirs.

---

### Decision 5: Whether to treat the 33 mouse samples as one experiment

**Pass 1 reasoning.** Eight samples had no metadata rows. The brief said not to
silently remove samples, so all 33 were merged and the eight were labelled
`unannotated` rather than dropped.

**Pass 2 verdict: WRONG, and this one matters.**

Those eight are a **different experiment**. They use different Cre driver lines
(Kit-MerCreMer, Car4-CreERT2, Ednrb-CreERT2), a *pre*-labelling design with
tamoxifen given 14–16 days **before** infection, and a single 19 dpi timepoint,
with infections run months later. The paper states plainly that the Ki67 line
was "analyzed separately" from those three lines — which is exactly why the
metadata CSV covers 25 of 33 samples.

So the final mouse object blends two experiments. Cell-type clustering is
largely unaffected — a fibroblast is a fibroblast — but `trace_call` means
*proliferation history* in the 25 and *cell-type identity at homeostasis* in the
eight, and the two must never be pooled. To reproduce the paper, split them.

**Lesson.** "Absent metadata" was treated as "incomplete annotation". It
actually meant "not part of this object". Being unwilling to drop data was the
right instinct; assuming the leftovers belonged to the same experiment was not.

---

## Part 2 — Things the deposited files could not tell me

Recorded because they change interpretation, not computation:

- **Cell-type proportions in the mouse data are an artefact.** Cells were
  MACS-fractionated and deliberately recombined at 85–90% CD45⁻ to 10–15%
  CD45⁺. Every composition output describes the sort ratio, not the lung.
  Compare composition only *within* a compartment.
- **The mouse homeostasis controls were never infected.** The `H1N1_` prefix on
  `H1N1_homeostasis_d03` is a naming artefact.
- **`SiteA`/`SiteB` are not two reporters.** They are two sub-regions of the
  single Ai14 ROSA26 allele — the SV40-polyA STOP cassette and the bGH-polyA
  after tdTomato. Removing them from the expression matrix was still correct.
  Pleasingly, the rule derived empirically from the data (`Traced` iff
  SiteB > SiteA) matches the published definition exactly, which independently
  confirms which site marks the recombined allele.
- **Both papers' headline findings require sub-clustering.** Niethamer's
  injury-induced capillary endothelial state is explicitly invisible at
  top-level clustering. Only the epithelium was subset here.
- **The GSE178360 `.RDS` objects contain the authors' own annotations** — donor
  IDs and published cell-type labels per barcode. They could not be read here
  (no R), but they are the route to ground-truth labels.

---

## Part 3 — Background for reading these figures

*(See the sections below for plain-language explanations of UMAP, Harmony and
Scrublet, written for a reader without a computational background.)*

---

## Part 4 - The focused re-analysis (the big goal, finally addressed)

Everything above describes the whole-atlas run. That run produced the required
UMAPs, dot plots and marker tables, but it **did not deliver the biology this
repository is about**. The evidence was unambiguous once looked for: the
deposited `Alveolar_transitional` (233 cells) and `AT1_AT2` (429 cells) both sat
*inside* the AT2 cluster, and no trajectory was fitted - so the
AT2 -> Krt8+ transitional -> AT1 **ordering**, which is the entire biological
claim, was nowhere quantified.

`analysis/scripts/06_regeneration_focus.py` fixes that. It does **not** redo QC,
merging or normalisation - those checkpoints are sound and cost hours. It reuses
`final_clustered.h5ad` and adds the four things that were missing.

### What changed, and why each mattered

**1. Restricted to the 25-sample Ki67 atlas.** The other 8 samples are a
different experiment. Without this filter the AT1 arm nearly doubles with cells
from Kit/Car4/Ednrb-driver mice - cluster 16 is only 41.8% cohort.

**2. A clean alveolar subset: clusters [10, 16] only.** Cluster 16 is *not*
optional - 30 transitional cells sit at the AT1 end, and cluster 10 alone
captures only 85% of them. Clusters 18 (ciliated) and 19 (Krt5+ dysplastic) are
excluded exactly as the source study did; including them would import >1,200
airway cells to gain 3 transitional ones. Result: **5,694 cells** capturing
229/233 transitional (98.3%) and 429/429 AT1_AT2 (100%).

*Implementation detail that matters:* the 180 Secretory cells hiding inside
cluster 10 are removed **by author label, not by dropping their subcluster** -
that subcluster also contains 24 transitional cells, 10% of the entire
transitional population.

**3. A data-derived transitional score, benchmarked rather than assumed.**
The canonical panel performs *worse* on this data than a score derived from it:

| score | genes | AUROC vs AT2 | vs AT1 | vs all |
|---|---:|---:|---:|---:|
| **core (data-derived)** | 8 | 0.988 | **0.984** | **0.987** |
| canonical panel | 11 | 0.980 | 0.883 | 0.947 |

The canonical panel loses 10 AUROC points against AT1 specifically. Measured
per gene, `Tnc` (0.510), `Hbegf` (0.574) and `Cdkn1a` (0.651) do not
discriminate here at all and were dropped; `Cldn4` is the near-binary gate
(70.8% expressing in transitional vs 2.4% in AT2). Final set: **Krt8, Krt18,
Cldn4, Sfn, Clu, S100a14, Gpx2, Tnip3**.

**4. PAGA + diffusion pseudotime**, rooted deliberately in the AT2 subcluster
with the highest AT2 score, because AT2 is the known progenitor.

### The result, and why it is trustworthy

Pseudotime orders the deposited labels correctly - **and those labels were never
used to fit it**, so this is a genuine independent check rather than a
restatement:

| deposited label | median pseudotime | n |
|---|---:|---:|
| AT2 | 0.013 | 3,277 |
| Alveolar_transitional | 0.179 | 229 |
| AT1_AT2 | 0.237 | 429 |
| AT1 | 0.327 | 1,759 |

The transitional state is now **resolved as its own subcluster** (subcluster 1
holds 168 of 229 transitional cells) where the whole-atlas run had it invisible
inside AT2.

And the timecourse independently reproduces the source paper's stated finding
that transitional cells peak at 11 dpi and are rare after 25 dpi:

| dpi | 0 | 6 | **11** | 19 | 25 | 42 | 90 | 366 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| median % transitional | 0.9 | 4.1 | **27.4** | 4.7 | 2.1 | 1.0 | 0.5 | 0.3 |

Three independent lines - pseudotime ordering, marker-score specificity, and the
injury timecourse - agree. None of them was available from the whole-atlas run.

### Caveats that must travel with these figures

- **Pseudotime is a rank ordering with no units.** It is not time, not
  comparable between runs, and depends on the root choice. The root here is an
  assumption (AT2 is the progenitor), stated rather than discovered.
- **PAGA edges are connectivity, not flux.** They do not indicate direction or
  rate of differentiation.
- **A trajectory can be fitted to data that is not a continuum.** The support
  here is that it independently recovers the correct label ordering and the
  correct injury timecourse - not the fit itself.
- The transitional cluster is a useful anchor, but the state is better treated
  as a **score along a continuum**: the cells the cluster misses sit at the AT2
  and AT1 ends, which is what a real continuum forced through a discrete
  boundary looks like.
- Composition percentages are still subject to the MACS sort caveat in
  `PIPELINE_AS_RUN.md`; they are comparable *within* the alveolar compartment,
  which is how they are used here.
