# Analysis trials motivated by Sikkema et al. 2023

Five trials. **S1, S3, S4 and S5 have run** and their artefacts are under
[`trials/`](trials/). **S2** is pre-registered here and not run. Status
words follow the root README claims table (Validated, Descriptive only,
Exploratory, Retracted-superseded, Not established). Every trial writes a
run record (`*_run_record.json`) with the rules frozen before any data is
read, the inputs with sizes and modification times, the package versions,
and the outputs.

Portfolio framing, on the owner's instruction of 2026-09-09: results are
described by cell state, repair and niche biology, macrophage and monocyte
states, annotation robustness and curation hygiene, matched to the target
labs (AT2 differentiation and stem-niche communication; donor-aware
cross-cohort annotation; inflammation-resolution and immune-epithelial
programmes; regulatory T cell representation). The mouse series is an
injury time course, and that fact is stated where it constrains a
comparison, but no result below is framed as infection or interferon
biology.

| Trial | Question | Status | Where |
|---|---|---|---|
| S1 | Do the HLCA cluster-entropy criteria flag anything in the tracked cluster tables? | run; Descriptive only | [`trials/s1_entropy_criteria.md`](trials/s1_entropy_criteria.md) |
| S3 | What does the HLCA consensus marker set call the human clusters, and does it recover AT0? | run; Descriptive only, AT0 Not established | [`trials/s3_hlca_marker_annotation/s3_summary.md`](trials/s3_hlca_marker_annotation/s3_summary.md) |
| S4 | What is mouse cluster 23, flagged by S1 as animal-private and unlabelled? | run; Descriptive only | [`trials/s4_cluster23_qc/s4_summary.md`](trials/s4_cluster23_qc/s4_summary.md) |
| S5 | Is the label disagreement in mouse cluster 5 a resolution artefact? | run; Descriptive only | [`trials/s5_cluster5_subclusters/s5_summary.md`](trials/s5_cluster5_subclusters/s5_summary.md) |
| S2 | Does scArches mapping to the HLCA core agree with the blind human annotations and recover AT0? | planned; environment Not established | section S2 below |

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
| Mouse strata | days post injury for the 25 annotated animals (2, 2, 2, 2, 2, 8, 4, 3 animals at 0, 6, 11, 19, 25, 42, 90, 366 days); the 8 lineage-tracing animals by Cre line (Kit 3, Car4 2, Ednrb 3) | the design; the two cohorts are separate experiments (root rationale, decision 5) |
| Human strata | one stratum of 3 healthy donors | the design |
| Negative control | the mouse donor entropy computed unstratified over all 33 animals | pre-registered expectation: it would over-flag time-point-specific states |

### Outcome

**Mouse, label entropy (24 of 29 clusters assessable).** High: clusters 5
(6,462 cells, entropy 1.85, top deposited label classical monocyte at 22%),
7 (2,038 cells, 1.14, CAP1 65%, only 23% labelled), 10 (5,763 cells, 0.75,
AT2 80%: the alveolar epithelium, where AT2, transitional and AT1/AT2 labels
co-occur by design) and 22 (5,152 cells, 0.65, AF2 70%: fibroblast subtypes;
the only high-entropy cluster that is also one of the three CONTRADICTED
clusters). Not assessable: 2, 11, 23, 27, 28. Clusters 0 and 25, the other
two CONTRADICTED clusters, have low label entropy (0.43, 0.09) because the
deposited labels agree with each other; their contradiction was between this
repository's marker panel and the deposition, which label entropy cannot
see. The two checks are complementary.

**Mouse, donor entropy.** Within stratum, one cluster flagged: 23, in the
42-day stratum (threshold 0.296) and in the Car4 tracing pair (threshold
0.199). Unstratified over 33 animals (threshold 0.372): no cluster flagged;
the pre-registered expectation did not materialise and stratification made
the criterion more sensitive, not less. The 50-cell floor, not the entropy
threshold, is the binding constraint of a two-animals-per-time-point design.

**Human, donor entropy (3 donors, threshold 0.233, Harmony embedding).**
Low: cluster 30 (222 cells, entropy 0.072, 98.7% from one donor). NA: 17,
29. In the tracked cluster-by-sample table three clusters exceed 75% from one
donor (24, 27, 30); the 95/5 criterion flags only 30. Label entropy: Not
established (deposited human labels unreadable without R).

### Claims from S1

| Claim | Status |
|---|---|
| Mouse clusters 5, 7, 10, 22 show high cross-source label disagreement by the HLCA criterion | Descriptive only |
| Mouse cluster 23 is animal-private within stratum and almost entirely unlabelled by the authors | Descriptive only (explained by S4) |
| Human cluster 30 is donor-private by the HLCA criterion after Harmony | Descriptive only (S3 adds that its identity is scheme-dependent) |
| Label entropy detects disagreement inside the deposition, not marker-panel error | Descriptive only |
| Stratified donor entropy is the right form of the criterion for the mouse design | Exploratory |

---

## S3. HLCA consensus-marker annotation of the human series

Run on 2026-09-09 (two runs; both recorded). Scripts and artefacts:
[`trials/s3_extract_hlca_markers.py`](trials/s3_extract_hlca_markers.py)
(writes the tracked marker table
[`trials/hlca_supp_table6_markers.csv`](trials/hlca_supp_table6_markers.csv)
from the local Supplementary Table 6, CC BY 4.0),
[`trials/s3_hlca_marker_annotation.py`](trials/s3_hlca_marker_annotation.py),
and the folder [`trials/s3_hlca_marker_annotation/`](trials/s3_hlca_marker_annotation/)
with the run record, the per-cluster assignment table, the cluster-by-type
score matrix, the AT0 tables, a heatmap and a three-panel UMAP.

### Pre-registration

| Rule | Value |
|---|---|
| Scores | scanpy `score_genes` on log1p(CP10K); ctrl_size 50, n_bins 25, random_state 0; three tiers per type as given in the sheet (compartment, intermediate level, own markers) |
| Scorable type | at least one own marker present in the object (see revision 1) |
| Per cell, primary (flat) | compartment = argmax of the four compartment scores; type = argmax of own-marker scores within that compartment |
| Per cluster | mode of the per-cell type; confident if the mode holds >= 50% of cells |
| Versus the blind proposals | fixed synonym table: agree / partial (same compartment) / disagree / reference lacks identity (neutrophil, platelet, erythroid are absent from the 61) |
| AT0 check | HLCA-argmax AT0 versus the strict gate of the doublet audit (SFTPC > 0, SCGB3A2 > 0, EPCAM > 0, PTPRC = PECAM1 = COL1A1 = 0, raw counts), per donor; concordant if within a factor of two for every donor; the epithelial subcluster 4 "AT0 candidate analogue" is also scored |
| Unit | donor; no P values |

**Rule revisions after the first run, disclosed.** The first version
required two own markers per type and treated every "Full atlas" marker as a
compartment marker. The sheet gives a single own marker to B cells, Plasma
cells, Pericytes, Tuft, EC aerocyte capillary and Smooth muscle FAM83D+, and
gives Neuroendocrine only atlas-level own markers, so eight types were
excluded and a fifth, spurious compartment appeared. Corrections: minimum
own markers one; compartments restricted to the four in the sheet; the
renamed gene C8orf4 read as TCIM; Neuroendocrine placed under Airway
epithelium by Supplementary Table 5. Because the flat argmax called AT0 in
most secretory and AT2 clusters, a **hierarchical assignment**
(compartment, then intermediate level, then type) was added as a post hoc
sensitivity analysis; the flat rule stays the primary. First-run outcome
kept in the run record: 53 of 61 types scorable; 16 agree, 12 partial, 1
disagree, 2 reference-lacks; AT0 over gate 1.64, 1.83, 1.10 per donor.

### Outcome (second run; 61 of 61 types scorable)

Agreement with the existing blind proposals over the 31 clusters:

| Scheme | agree | partial, same compartment | partial, compartment-free | disagree | reference lacks identity |
|---|--:|--:|--:|--:|--:|
| flat (primary) | 17 | 10 | 1 | 1 | 2 |
| hierarchical (sensitivity) | 19 | 8 | 1 | 1 | 2 |

Cluster-level findings worth keeping:

- **Cluster 22 is mast cells, not B cells.** The blind panel had no mast set
  and called it B_cell with low confidence; both HLCA schemes call Mast
  cells (mode 0.68), and the cluster's own top genes in the existing table
  are TPSB2, TPSAB1, CPA3 and KIT. Proposed correction to the human
  annotation table, pending owner decision.
- Cluster 20 is B cells (0.74 flat, 0.78 hierarchical); cluster 10 refines
  Endothelial to EC aerocyte capillary; cluster 0 is Interstitial
  macrophages (perivascular) and cluster 1 is DC2, both inside the blind
  Mono_Mac call; cluster 12, blind Transitional with runner-up Basal, is
  Basal resting (0.69, 0.75); cluster 26, blind AT1, is AT1 only in the
  hierarchical scheme (0.70) and pre-TB secretory in the flat one (0.46).
- The two clusters the reference cannot name behave as the paper predicts:
  neutrophils (cluster 14) are forced to Classical monocytes, and the
  low-confidence platelet-like cluster 15 has no compartment majority at all
  (0.41), which is the signature of an identity absent from the reference.
- Cluster 29 (31 cells, blind Neuroendocrine with ASCL1, CALCA, CHGA) is
  not recovered by the HLCA neuroendocrine own markers (CELF3, SLC6A17,
  CDK5R2); the flat scheme says AT0 (0.23) and the hierarchical scheme
  Ionocyte (0.26). Too small to carry a claim either way.
- Cluster 30, donor-private in S1, is AT0 (0.89) in the flat scheme and AT2
  (0.46) in the hierarchical one while its own top genes are ciliary. It
  remains suspicious and should not be cited as a population.

**AT0.** Strict gate per donor: 596, 880, 134 cells.

| Scheme | AT0 per donor | AT0 over gate | Within factor 2 | Epithelial subcluster 4 (n 328) |
|---|---|---|---|---|
| flat (primary) | 1,018; 1,631; 152 | 1.71; 1.85; 1.13 | yes for all three | AT0 280, pre-TB secretory 37, AT2 4 |
| hierarchical (sensitivity) | 106; 167; 8 | 0.18; 0.19; 0.06 | no for all three | AT2 251, Tuft 37, AT1 24, AT0 0 |

The two schemes give opposite answers from the same 61 marker sets. The
reason is structural, not a bug: AT0's own markers (SFTPB, SCGB3A2, SFTA2)
are broadly expressed across distal secretory and AT2 cells, AT2's own
markers (MFSD2A, TCIM, C11orf96) are weak, and the HLCA places AT0 under
Airway epithelium although AT0-like cells score alveolar on the
intermediate markers (SFTA2 appears in both sets), so the hierarchical route
sends them to AT2 before AT0 can be considered. The HLCA itself did not
define AT0 by marker scoring; it defined it by nested clustering of the
integrated embedding plus expert consensus. Marker transfer therefore cannot
settle AT0 here in either direction. The strict gate count sits between the
two schemes. The concordance rule passed for the pre-registered primary and
failed for the sensitivity; both are reported and the AT0 claim by this
route is **Not established**. The reference-mapping trial S2 is the route
that can settle it.

### Claims from S3

| Claim | Status |
|---|---|
| The blind annotations agree with the HLCA consensus at type level for 17 (flat) or 19 (hierarchical) of 31 clusters, with a further 10 or 8 agreeing at compartment level | Descriptive only |
| Cluster 22 is mast cells | Descriptive only; proposed correction pending owner decision |
| Neutrophils and the platelet-like cluster are identities the HLCA cannot name, and behave accordingly | Descriptive only (matches the paper's own account) |
| AT0 identity by HLCA marker transfer | Not established (scheme-dependent; see above) |
| Cluster 30 is a real population | Not established (donor-private; identity scheme-dependent) |

---

## S4. Mouse cluster 23

Run on 2026-09-09. Script and artefacts:
[`trials/s4_cluster23_qc.py`](trials/s4_cluster23_qc.py),
[`trials/s4_cluster23_qc/`](trials/s4_cluster23_qc/) (run record, QC medians
by group, lineage-marker detection, sample composition, boxplots).

### Pre-registration

Rules were set after the cluster-level QC summary (median 1,868 counts, 78%
from one sample) had been seen and before any per-cell data was read; this
is stated in the run record.

| Rule | Value |
|---|---|
| Comparators | all other cells; a random sample of the same size from outside the cluster (seed 0); the dominant sample outside the cluster |
| Low-count call | median total counts and median genes both < 0.5 x the atlas medians |
| Ambient-like call | fraction of cells detecting >= 3 of 5 lineage markers (Sftpc, Scgb1a1, Ptprc, Pecam1, Col1a1) >= 2 x the random comparator |
| Cohort split | tracing samples carry no author labels at all; the labelled fraction is compared within the annotated cohort only |
| Verdict | both calls true: "low-count, ambient-like"; otherwise the single call, or "unexplained" |

### Outcome: **low-count, ambient-like**

| Group | Cells | Median counts | Median genes | Top-20-gene fraction |
|---|--:|--:|--:|--:|
| cluster 23 | 2,912 | 1,868 | 1,096 | 15.4% |
| all other cells | 159,263 | 6,790 | 2,481 | 21.7% |
| random comparator | 2,912 | 6,774 | 2,467 | 21.7% |
| sample EEM-scRNA-289 outside cluster 23 | 1,870 | 10,563 | 3,455 | 17.6% |

- 84.3% of cluster-23 cells detect at least three of the five lineage
  markers, against 4.5% in the random comparator; Sftpc and Scgb1a1 are
  detected in 99% of its cells alongside Pecam1 (61%) and Col1a1 (73%).
  The flat top-20-gene fraction points the same way: no dominant programme,
  a little of everything.
- 78% of the cluster comes from EEM-scRNA-289, a Car4-CreERT2 tracing
  animal, where it makes up 55% of that sample's cells; that sample's median
  counts are 4,141, the lowest of the tracing cohort.
- Within the annotated cohort, 2.2% of cluster-23 cells carry an author
  label, against 86% cohort-wide; the 14 labelled cells are CAP1, plasma
  cell, B cell, AT2 and fibroblast. Consistent with the authors' QC having
  removed these barcodes before annotation.

### Claims from S4

| Claim | Status |
|---|---|
| Cluster 23 is a low-count, ambient-like barcode population, not a cell type | Descriptive only |
| EEM-scRNA-289 has half of its cells in that population and should be flagged in the lineage-tracing cohort analysis | Descriptive only; proposed known issue, pending owner decision |

---

## S5. Mouse cluster 5 subclusters

Run on 2026-09-09. Script and artefacts:
[`trials/s5_cluster5_subclusters.py`](trials/s5_cluster5_subclusters.py),
[`trials/s5_cluster5_subclusters/`](trials/s5_cluster5_subclusters/) (run
record, grading tables at both resolutions, subcluster-by-label counts, top
genes, comparator-panel detection fractions, UMAP).

### Pre-registration

| Rule | Value |
|---|---|
| Cells | all 6,462 cells of atlas cluster 5; no batch correction (root decision 2) |
| Features and embedding | seurat_v3 HVGs on raw counts (2,000), scale (max 10), PCA 30, random_state 0 |
| Graph and clustering | kNN k 30 on 30 PCs; Leiden 0.2 primary and 0.5 sensitivity (the HLCA level-2 rule); igraph flavour, 2 iterations, seed 0 |
| Grading | deposited labels held out; purity and label entropy per subcluster; "resolved by resolution" if >= 75% of labelled cells fall in subclusters with purity >= 0.75 at the primary resolution |
| Interpretation only | Wilcoxon top genes; detection fraction of a fixed comparator panel including the HLCA monocyte-derived macrophage subtypes through mouse orthologs |
| Unit | animal for composition; no P values |

### Outcome: **not resolved at Leiden 0.2; resolved at Leiden 0.5**

| Resolution | Subclusters | Labelled cells in subclusters with purity >= 0.75 | Cell-weighted mean label entropy |
|---|--:|--:|--:|
| atlas cluster 5, before | 1 | | 1.85 |
| 0.2 (primary) | 9 | 59.6% | 0.55 |
| 0.5 (sensitivity) | 13 | 94.3% | 0.39 |

The pre-registered primary fails its own rule and the sensitivity passes it;
both are reported. At 0.5 the states are clean at the animal level:
classical monocytes (purity 0.96; Ly6c2, Ccr2, Plac8), non-classical
monocytes (0.97; Ace, Nr4a1, Ear2), interstitial macrophages (0.95; C1qa,
C1qb, C1qc, Apoe, Trem2, Mrc1), cDC1 (0.89 and 0.95; Xcr1, Clec9a), cDC2
(0.86), migratory DC (0.94; Ccr7, Fscn1), and two monocyte-like subclusters
labelled iMON (0.84, 0.89). Three subclusters stay impure: one small
interstitial-macrophage subcluster (98 cells, 0.51), one carrying T-cell
genes (147 cells; Cd2, Ltb, Trbc2) and one carrying endothelial genes (164
cells; Cdh5, Pecam1, Cldn5), which read as myeloid-lymphoid and
myeloid-endothelial doublets the per-capture Scrublet run did not remove.

Two further observations:

- A 97-cell subcluster with Siglech, Bst2, Tcf4 and Ccr9 detected in nearly
  every cell is plasmacytoid-DC-like; the deposition labels it cDC2 with
  purity 0.99 because the deposited label set has no pDC class. This is a
  refinement over the deposited labels, not a disagreement with them.
- The subcluster labelled iMON at the primary resolution with Isg15, Irf7
  and Ifit genes as its top markers is 75% from one animal at day 6. By the
  S1 donor rule it is animal-private within its stratum and is reported as
  such; no narrative is attached to it.

HLCA comparators through mouse orthologs: the C1QA-high monocyte-derived
macrophage programme (C1qa 0.92, Trem2 0.65, Mrc1 0.73 detected) maps onto
the interstitial-macrophage subcluster; the CCL2-high programme (Ccl2 0.75,
Il1rn 0.73) onto the day-6 monocyte subcluster; the SPP1-high profibrotic
programme is not detected at population level anywhere in cluster 5 (Spp1
at most 0.33, Chil3 the only member above 0.5, and only in monocytes), and
MARCO-high cells are absent because alveolar macrophages sit in another
atlas cluster.

### Claims from S5

| Claim | Status |
|---|---|
| The cross-source label disagreement of atlas cluster 5 is a resolution artefact for 94% of labelled cells at Leiden 0.5, k 30, but not at 0.2 | Descriptive only |
| A plasmacytoid-DC-like population exists inside the deposited cDC2 label | Descriptive only |
| Two doublet-like subclusters (311 cells) survive per-capture Scrublet | Descriptive only |
| The HLCA SPP1-high profibrotic macrophage programme is not detected in cluster 5 | Descriptive only (negative) |

### What S3, S4 and S5 change

Nothing in the tracked results has been edited; the generated reports are
left as they are until the owner decides. Proposed, pending decision:
record cluster 23 and sample EEM-scRNA-289 in the mouse known issues;
correct human cluster 22 to mast cells with a note; add the cluster-5
subcluster table to the mouse epithelial-style focused analyses; adopt the
entropy screen and the resolution rule as gates G3 in
[`PIPELINE_FRAMING.md`](PIPELINE_FRAMING.md).

---

## S2. Reference mapping of GSE178360 to the HLCA core

Status: **PLANNED, not run. Environment: Not established.** After S3, this
is the only route that can settle the AT0 question, because it reproduces
how the HLCA itself defined AT0 (position in the integrated embedding plus
label transfer), rather than marker scoring.

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

**Comparisons.** (1) Transferred level-3 and level-4 labels against the
blind annotations and against the S3 assignments, per donor. (2) Against the
deposited author labels, only if the RDS objects become readable. (3)
Transferred AT0 and pre-TB secretory counts against the strict gate and
against both S3 schemes; concordant if within a factor of two per donor.
(4) Per-cluster mean uncertainty, with cluster 30 as a named check.

**Success and failure, frozen.** The paper's healthy demonstration (68%
correct, 14% incorrect, 18% unknown at cutoff 0.3) is the comparator for
agreement. A mean uncertainty above the range the paper reports for healthy
adult 10x datasets is recorded as "mapped poorly". Expected negative
results: the Zenodo model may not load under current scArches or scvi-tools;
training under x86-64 emulation may be too slow; zero-filling may remove a
material fraction of the 2,000 HVGs.

**Environment gate, before anything else.** Dry-run install of scvi-tools
and scArches into the emulated interpreter with uv; refuse if any package
would need compilation; record the outcome here. No install on the native
ARM64 interpreter.

---

## Not planned

Mapping GSE262927 (mouse) to the HLCA: declined. The reference is human and
adult, and the authors state that cross-species mapping needs further method
development. The HLCA's criteria apply to the mouse series (S1, S4, S5); its
embedding does not.
