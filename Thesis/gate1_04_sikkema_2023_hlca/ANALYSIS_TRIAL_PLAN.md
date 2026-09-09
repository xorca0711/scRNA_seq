# Analysis trials motivated by Sikkema et al. 2023

Two trials. **S1** applies the paper's cluster-screening criteria to tables
already tracked in this repository (standard library only; no cell-level
computation). **S2** is the reference-mapping trial the roadmap names for
this paper; it is planned and pre-registered here, not run. Status words
follow the root README claims table (Validated, Descriptive only,
Exploratory, Retracted-superseded, Not established).

---

## S1. Cluster-entropy criteria on the tracked cluster tables

Run on 2026-09-09. Script and artefacts:
[`trials/s1_entropy_criteria.py`](trials/s1_entropy_criteria.py),
[`trials/s1_entropy_criteria.json`](trials/s1_entropy_criteria.json),
[`trials/s1_entropy_criteria.md`](trials/s1_entropy_criteria.md).

### Pre-registration (frozen before the tables were read)

| Rule | Value | Origin |
|---|---|---|
| Label entropy is high | > 0.56 (natural log) | HLCA: entropy of a 75/25 two-label cluster |
| Label entropy is reported only if | >= 20% of the cluster carries a deposited label | HLCA |
| Donor entropy is low | < entropy of a cluster with 95% of cells from one donor and 5% spread evenly over the other n minus 1, computed **per stratum** | HLCA rule; the HLCA value 0.43 encodes 107 donors and is not reused |
| Minimum cells per cluster and stratum | 50; below it the pair is NA | HLCA multidisease rule (50 MDMs per dataset) |
| Mouse strata | days post infection for the 25 annotated animals (2, 2, 2, 2, 2, 8, 4, 3 animals at 0, 6, 11, 19, 25, 42, 90, 366 dpi); the 8 lineage-tracing animals by Cre line (Kit 3, Car4 2, Ednrb 3) | the design; the two cohorts are separate experiments (root rationale, decision 5) |
| Human strata | one stratum of 3 healthy donors | the design |
| Negative control | the mouse donor entropy computed unstratified over all 33 animals | pre-registered expectation: it would over-flag time-point-specific states |

Inputs, all tracked: the mouse cluster-by-deposited-label counts, the mouse
and human cluster-by-sample counts, the mouse annotation-proposal table (for
the CONTRADICTED flags), and the two per-sample tables that carry day and
Cre-line for each mouse sample.

### Outcome

**Mouse, label entropy (24 of 29 clusters assessable).**

| Cluster | Cells | Labelled | Entropy | Top deposited label | Reading |
|--:|--:|--:|--:|---|---|
| 5 | 6,462 | 73% | 1.85 | classical monocyte, 22% | mixed myeloid cluster; the HLCA's own worst case was the same disagreement (monocyte, macrophage, DC2). Subcluster before any myeloid claim. |
| 7 | 2,038 | 23% | 1.14 | CAP1, 65% | just above the 20% reporting floor; endothelial mixture |
| 10 | 5,763 | 71% | 0.75 | AT2, 80% | alveolar epithelium where AT2, transitional and AT1/AT2 labels co-occur by design of an injury time course; the HLCA's "continuous transition" case, not a mislabel |
| 22 | 5,152 | 66% | 0.65 | AF2, 70% | fibroblast subtypes AF1 and AF2; the only high-entropy cluster that is also one of the three CONTRADICTED clusters |

Not assessable (fewer than 20% labelled): clusters 2, 11, 23, 27, 28.
Clusters 23 (2,912 cells, 0.5% labelled) and 27 (54 cells, none labelled)
consist almost entirely of cells the authors did not annotate.

Relation to the existing contradiction flags (clusters 0, 22, 25): only 22
overlaps. Clusters 0 and 25 have *low* label entropy (0.43 and 0.09) because
the deposited labels agree with each other (CAP1 89%, peribronchial
fibroblast 99%); the contradiction there was between this repository's
marker panel and the deposition, which label entropy cannot see. The two
checks are complementary, not redundant.

**Mouse, donor entropy.**

- Within stratum: one cluster flagged, **23**, in the 42 dpi stratum (8
  animals, threshold 0.296) and in the Car4 tracing pair (threshold 0.199).
  Together with its 0.5% labelled fraction, cluster 23 reads as an
  animal-private population of cells that the authors removed before
  annotation. It was not flagged by any existing check and should be excluded
  from any composition or state claim.
- Unstratified over 33 animals (threshold 0.372): **no cluster flagged.** The
  pre-registered expectation that the unstratified version would over-flag
  time-point-specific states did not materialise: with 33 animals the 95/5
  threshold is low enough that time-point-specific clusters keep enough
  entropy. Stratification made the criterion *more* sensitive here, not less.
  Recorded as observed.
- Many cluster-by-stratum pairs are NA because two animals per time point
  rarely contribute 50 cells of a small cluster. In this design the 50-cell
  floor, not the entropy threshold, is the binding constraint.

**Human, donor entropy (3 donors, threshold 0.233, Harmony embedding).**

- Low: **cluster 30** (222 cells, entropy 0.072, 98.7% from one donor).
  NA: clusters 17 and 29 (fewer than 50 cells).
- Against the existing check: in the tracked cluster-by-sample table three
  clusters exceed 75% from one donor (24 at 87%, 27 at 80%, 30 at 99%); the
  HLCA-style 95% criterion keeps 24 and 27 and flags only 30. The figure
  "4 of 31 clusters above 75% one donor" in the handoff notes comes from a
  different table (`qc/celltype_split_by_sample.csv`, counted per proposed
  identity) and is not the same quantity; neither number is wrong, they count
  different things.
- Label entropy: **Not established.** The deposited human labels are inside
  RDS objects that cannot be read on this machine.

### Claims from S1

| Claim | Status |
|---|---|
| Mouse clusters 5, 7, 10, 22 show high cross-source label disagreement by the HLCA criterion | Descriptive only |
| Mouse cluster 23 is animal-private within stratum and almost entirely unlabelled by the authors | Descriptive only |
| Human cluster 30 is donor-private by the HLCA criterion after Harmony | Descriptive only |
| Label entropy detects disagreement inside the deposition, not marker-panel error (clusters 0, 25) | Descriptive only |
| Stratified donor entropy is the right form of the criterion for the mouse design | Exploratory (one dataset; the unstratified control did not behave as predicted) |

### What S1 changes

Nothing in the existing results; it adds flags. Proposed consequences,
pending owner decision: record cluster 23 in the mouse report's known
issues; subcluster cluster 5 before any myeloid claim; adopt the entropy
screen as gate G3 in [`PIPELINE_FRAMING.md`](PIPELINE_FRAMING.md).

---

## S2. Reference mapping of GSE178360 to the HLCA core

Status: **PLANNED, not run. Environment: Not established.**

**Question.** Do the HLCA's transferred labels agree with this repository's
blind annotations of the human distal-lung series, and is the AT0 population
recovered by a route that uses neither this repository's marker gates nor
Scrublet? GSE178360 is an HLCA *extension* dataset (Tata_unpubl), never seen
by the core model, so the test is not circular.

**Pre-registered design.**

| Element | Setting | Origin |
|---|---|---|
| Query | the post-QC, post-Scrublet matrix (27,729 cells), raw counts, per donor | this repository |
| Gene space | the HLCA's 2,000 HVGs; symbols harmonised through the Ensembl/HGNC route; HVGs absent from the query zero-filled | paper, Methods |
| Reference | HLCA core scANVI model, Zenodo 10.5281/zenodo.7599104 | paper |
| Mapping | scArches surgery: freeze dropout, 500 surgery epochs, weight decay 0, early stopping on ELBO (patience 10, threshold 0.001) | paper |
| Label transfer | kNN k 50 in the joint embedding; u = 1 minus the weighted fraction of neighbours with the transferred label | paper |
| Unknown cutoff | u > 0.2 (primary); 0.3 reported as sensitivity | paper (calibrated) |
| Unit of comparison | the donor (n = 3); every quantity reported per donor; no P values | owner's rule |

**Comparisons.**

1. Transferred HLCA level-3 and level-4 labels against the blind cluster
   annotations (confusion table per donor).
2. Against the deposited author labels, only if the RDS objects become
   readable.
3. Transferred AT0 and pre-TB secretory counts against the SFTPC+ SCGB3A2+
   EPCAM+ lineage-negative gate count from the doublet audit; concordant if
   within a factor of two per donor, otherwise discordant and reported as
   such.
4. Per-cluster mean uncertainty, with cluster 30 (donor-private in S1) as a
   named check, and the same mapping on the unintegrated embedding's
   clusters to see whether Harmony changed what the reference recognises.

**Success and failure, frozen.** The paper's healthy demonstration (68%
correct, 14% incorrect, 18% unknown at cutoff 0.3) is the comparator for
agreement; the trial is informative if per-donor agreement is reported with
its unknown fraction, whatever the value. A mean uncertainty for the query
above the range the paper reports for healthy adult 10x datasets would be
recorded as "mapped poorly", not explained away. Expected negative results:
the Zenodo model may not load under current scArches or scvi-tools; training
under x86-64 emulation may be too slow to finish; zero-filling may remove a
material fraction of the 2,000 HVGs.

**Environment gate, before anything else.** Attempt a dry-run install of
scvi-tools and scArches into the emulated interpreter with uv; refuse if any
package would need compilation; record the outcome in this file. No install
is attempted on the native ARM64 interpreter.

---

## Not planned

Mapping GSE262927 (mouse) to the HLCA: declined. The reference is human and
adult, and the authors state that cross-species mapping needs further method
development. The HLCA's criteria apply to the mouse series (S1); its
embedding does not.
