# Pipeline as run

**This file is generated** by `analysis/scripts/05_write_pipeline_as_run.py`
from the artefacts the pipeline itself wrote. It records what was *actually
executed*, as opposed to what the reference material in `docs/` describes.
Re-generate it after any analysis run; do not hand-edit it.

> **Scope.** The other pages in `docs/` (SoupX, scds, Slingshot, tradeSeq,
> Scrublet) are **tool reference notes**. They are kept deliberately, and they
> do *not* imply that every tool described there was part of this analysis.
> The table immediately below is the authoritative statement of which ones were.

## Which of the documented tools were actually used

| Tool | Status in this analysis | Notes |
|---|---|---|
| [Scrublet](SCRUBLET.md) | **USED** | run per capture via `sc.pp.scrublet`, before merging |
| [SoupX](SOUPX.md) | **NOT USED** | no ambient-RNA correction was applied to either dataset |
| [scds](SCDS.md) | **NOT USED** | Scrublet was the only doublet caller; scds is R-only and R is unavailable here |
| [Slingshot](SLINGSHOT.md) | **NOT USED** | no trajectory or pseudotime analysis was performed |
| [tradeSeq](TRADESEQ.md) | **NOT USED** | follows from Slingshot not being run; no trajectory-based DE |

## What was NOT done

This section exists so that nobody reading `docs/` concludes that these steps
are represented in the results. They are not.

### Ambient RNA / soup correction

*Reference tools: SoupX, CellBender, decontX*

The raw (unfiltered) droplet matrices needed for SoupX are present for GSE178360 but not for GSE262927, whose deposited files are cell-called only. No correction was applied to either dataset, so ambient contamination is uncorrected - this matters most for GSE178360 sample DD046Q, which carries heavy haemoglobin signal.

### Trajectory / pseudotime

*Reference tools: Slingshot, PAGA, DPT, Monocle*

Not run. The AT2 -> transitional -> AT1 ordering discussed in the literature is NOT quantified anywhere in these outputs; the epithelial sub-analysis resolves the states as clusters only.

### Trajectory-based differential expression

*Reference tools: tradeSeq*

Not run; depends on a trajectory that was not fitted.

### RNA velocity

*Reference tools: scVelo, velocyto*

Not run; requires spliced/unspliced counts that the deposited files do not contain.

### Cell-cycle scoring or regression

*Reference tools: Seurat CellCycleScoring*

Not computed. `cell_cycle_phase` in the metadata was READ from the deposited GSE262927 author table, never recalculated, and was never regressed out.

### SCTransform / covariate regression

*Reference tools: Seurat SCTransform*

Not used; normalisation was normalize_total + log1p.

### Reference mapping / label transfer

*Reference tools: scVI, scANVI, Seurat anchors*

Not run. Author cell-type labels were used only as an independent check of the clustering, never to assign identities.

### Pathway / gene-set enrichment

*Reference tools: GSEA, GO, fgsea*

Not run.

### Reading the deposited Seurat objects

*Reference tools: R, Seurat*

The four GSE178360 `.RDS` objects were never opened. They are valid R 3.6.3 RDS v3 files, but R is not installed and there is no C compiler, so the analysis started from the `filtered_feature_bc_matrix.h5` files.


## How this analysis differs from the source publications

The two papers behind these accessions are in `Thesis/Primary/`. They were
**not** consulted while the pipeline was being built - the analysis was driven
purely by the deposited files, as the brief required - and were read only
afterwards. This section records where the executed analysis agrees with, and
where it departs from, the published work.

### GSE262927 - Niethamer et al. 2025, *Cell Stem Cell* 32:302-321

| | Published | This analysis |
|---|---|---|
| Alignment | STAR-Solo v2.7.9a, mm39 | not re-aligned; deposited matrices used as given |
| Ambient RNA | **SoupX v1.6.0** | **none** |
| Doublets | **scds AND Scrublet** | **Scrublet only** |
| Normalisation | **SCTransform**, regressing mito/nFeature/nCount | normalize_total(1e4) + log1p, no regression |
| Clustering | **Louvain, resolution 1.0** | **Leiden, resolution 0.3** |
| Batch correction | **none - libraries merged directly** | **none** - same conclusion |
| Trajectory | **Slingshot + tradeSeq** | **none** |
| Cell filtering | thresholds **plus manual cluster curation** | thresholds only |

**Agreement on the decision that mattered.** The authors performed no
integration of any kind; they merged libraries and clustered. This pipeline
independently reached the same conclusion from the mixing statistics, on the
grounds that sample is nested inside experimental group.

**Three corrections to what the deposited files alone suggested:**

1. **The 8 samples carrying no metadata are a different experiment**, not an
   unannotated part of the atlas. `EEM-scRNA-249/250/251` use Kit-MerCreMer,
   `-288/289` use Car4-CreERT2 and `-290/291/292` use Ednrb-CreERT2, all with
   tamoxifen given *before* infection and all sacrificed at 19 dpi. The paper
   analysed them separately from the Ki67 atlas, which is exactly why the
   metadata CSV covers 25 of 33 samples. **This pipeline merged all 33 into one
   object, so the final object blends two experiments.** `trace_call` means
   proliferation history in the 25 and cell-type identity at homeostasis in the
   8; the two are not comparable.
2. **`SiteA`/`SiteB` are not two reporter contigs.** They are two sub-regions of
   the single Ai14 ROSA26 allele - the SV40-polyA STOP cassette and the
   bGH-polyA downstream of tdTomato. Removing them from the expression matrix
   remains correct. The rule this pipeline derived empirically from the data
   (`Traced` iff SiteB > SiteA) matches the published definition
   (`SiteB/(SiteA+SiteB) > 0.5`), which independently confirms that SiteB marks
   the recombined allele.
3. **Cell-type proportions are an artefact of the sort, not biology.** Cells
   were MACS-fractionated and deliberately recombined at 85-90% CD45-negative
   to 10-15% CD45-positive. **Every composition output in `figures/composition/`
   and `tables/cluster_percent_within_*.csv` therefore describes the sorting
   ratio, not the lung.** Compare compositions only within a compartment.

Also worth knowing: the homeostasis controls were **never infected** - the
`H1N1_` prefix on `H1N1_homeostasis_d03` is a naming artefact. And the paper's
headline finding, an injury-induced capillary endothelial state, is explicitly
only resolvable after subsetting the capillary endothelium; it does not appear
at top-level clustering, and this analysis subset only the epithelium.

### GSE178360 - Kadur Lakshminarasimha Murthy et al. 2022, *Nature* 604:111-119

| | Published | This analysis |
|---|---|---|
| Alignment | Cell Ranger v5.0.1, GRCh38 GENCODE v32 | not re-aligned; deposited matrices used as given |
| Ambient RNA | **SoupX, before normalisation** | **none** |
| Doublets | **manual cluster inspection only** | **Scrublet, automated** |
| Normalisation | **SCTransform** | normalize_total(1e4) + log1p |
| Integration | **Seurat CCA anchors, 10,000 features** | **Harmony** - same intent, different tool |
| Trajectory | RNA velocity, PAGA, Slingshot, Monocle3, scEpath, MuTrans | **none** |
| Cells retained | ~21,700 | 27,729 |

**The integration decision is correct, but this pipeline's stated reason was
not.** The three samples are healthy-donor biological replicates of one
anatomical preparation (microdissected distal airways under 1 mm), and the
authors themselves integrated across the same three donors. The pipeline's
recorded justification - that no condition metadata exists, therefore no
designed contrast can be destroyed - is a non-sequitur: a missing metadata file
is not evidence of a missing design. The conclusion happens to hold because the
paper's scientific axis, proximal-to-distal zonation, is a gradient *within*
each sample that sample-level correction cannot touch. The divergent gene count
in DD073R is a genuine technical batch effect confounded with sample, which
correction actively helps.

**A measured limitation of the doublet filtering.** Every novel population in
this paper is defined by co-expression of two lineages' markers, which is
precisely what an automated doublet caller is built to remove - and the authors
used none. Measured over the marker definitions directly:

| Population | cells | called doublet | vs. dataset baseline |
|---|---:|---:|---:|
| distal-BC (TP63+ SFTPB+) | 955 | 12.9% | 2.03x |
| SCGB3A2-CC (FOXJ1+ SCGB3A2+ SFTPB+) | 1,559 | 10.5% | 1.66x |
| AT0 (SFTPC+ SCGB3A2+) | 5,680 | 8.0% | 1.26x |
| all cells | 29,605 | 6.3% | 1.00x |

Scrublet removed these populations at up to twice the background rate. It did
not erase them - 87-92% survived - but the bias is systematic and in the
direction that loses the paper's discoveries. Anyone pursuing AT0 should re-run
without doublet filtering, or inspect the flagged cells before removing them.
Skipping ambient-RNA correction compounds this: surfactant transcripts are
classic soup contaminants, and AT0 is defined partly by SFTPC.

**The deposited `.RDS` objects do contain the authors' annotations** - donor
IDs and published cell-type labels per barcode - even though this analysis
could not read them for want of R. They are the route to ground-truth labels
if R becomes available.


## The workflow that was executed

### GSE262927 (mouse, 33 samples)

| Step | What was actually run |
|---|---|
| Input directory | `raw_data/GSE262927/GSE262927_RAW` |
| Input | 10x CellRanger HDF5 (v3 /matrix layout) |
| Species detection | from gene symbols in the matrix, not assumed -> **mouse** |
| Non-gene features removed | ['SiteA', 'SiteB'] |
| Gene space | inner join: 55359-55359 per sample -> 55359 shared |
| Metadata | raw_data/GSE262927/GSE262927_CellMetaData.csv |
| QC thresholds | per sample, MAD-derived (see `analysis/qc/qc_thresholds.csv`) |
| Cells | 212,701 in -> 169,807 after QC -> 162,175 after doublet removal |
| Doublets | Scrublet per capture, 7,632 removed |
| Normalisation | normalize_total(target_sum=1e4) + log1p; the raw counts are preserved in processed/postQC.h5ad and re-attached as layers['counts'] of the final object |
| HVG | seurat_v3 on raw counts, n_top_genes=2500, batch_key=sample_id -> 2500 genes |
| Scaling | unit-variance scaling (z-score, clipped at +/-10) on the HVG subset, then randomized PCA |
| PCA | 50 (variance-ratio elbow at 11, 90% cumulative variance at 51, bracketed to [20,50]) |
| Batch correction | not applied, and the optional supplementary Harmony embedding was skipped: correction is not warranted here (within-group replicate enrichment 1.374 <= 2.0, i.e. replicates inside an experimental group already mix), and at 162,175 cells across 33 batches Harmony would cost hours for an embedding that would not be used. Pass --integration harmony to force it. The mixing statistics are in qc/batch_assessment.json. |
| Neighbours | n_neighbors=15, n_pcs=50, use_rep=X_pca, metric=euclidean |
| UMAP | scanpy defaults min_dist=0.5, spread=1.0, random_state=0 |
| Clustering | Leiden (igraph flavour, 2 iterations) scanned at resolutions [0.3, 0.5, 0.8, 1.0]. Primary resolution: 0.3. No resolution satisfied all three criteria, so the resolution with the fewest under-separated neighbouring cluster pairs (1) and the fewest undersized clusters (2) was taken. Every finer resolution produced neighbouring clusters that no gene separates, which is the signature of over-clustering. All resolutions are retained as leiden_res* columns. |
| Final clusters | **29** at resolution 0.3 |
| Marker test | Wilcoxon rank-sum (scanpy rank_genes_groups, one cluster vs all remaining cells) on log1p(CP10K) values, with expressing fractions; table re-used from the previous run of this same pipeline |
| Annotation | candidate identities scored per cluster from z-scored mean expression of curated panels, down-weighted by the fraction of cells expressing each marker; numeric Leiden labels are preserved and no candidate is promoted to a definitive label |
| Epithelial sub-analysis | performed on 13333 cells from clusters ['10', '16', '18', '19']; 16 epithelial subclusters at Leiden resolution 0.6 |
| Condition DE | EXPLORATORY ONLY. Sample-level pseudobulk with n=2 homeostasis and n=23 infected samples. The homeostasis arm is not replicated within a group and is confounded with sex, so this is not a confirmatory differential-expression result. |

Resolution scan (why that resolution was chosen):

| resolution | n_clusters | pct_clusters_marker_supported | min_DE_genes_between_nearest_clusters | nearest_cluster_pairs_under_threshold |
|---|---|---|---|---|
| 0.3 | 29.0 | 100.0 | 9.0 | 1.0 |
| 0.5 | 32.0 | 100.0 | 8.0 | 1.0 |
| 0.8 | 40.0 | 100.0 | 0.0 | 3.0 |
| 1.0 | 43.0 | 100.0 | 0.0 | 3.0 |

Batch mixing statistics actually measured:

```json
{
  "knn_same_sample_enrichment": 4.181,
  "mean_same_sample_neighbour_fraction": 0.1272,
  "expected_by_chance": 0.0319,
  "within_group_replicate_enrichment": 1.374,
  "replicate_structure_available": true,
  "decision_rule": "Integrate only if biological replicates within the SAME experimental group fail to mix (within-group replicate kNN enrichment > 2.0). Sample identity alone is not evidence of a technical batch effect here, because every sample belongs to exactly one experimental group - correcting on sample would also remove the experimental effect.",
  "integration_recommended": false,
  "supplementary_harmony_warranted": true
}
```

Full record: `analysis/logs/analysis_log.txt`, `analysis/logs/decisions.json`,
`analysis/README.md`.


### GSE178360 (human, 3 samples)

| Step | What was actually run |
|---|---|
| Input directory | `raw_data/GSE178360/GSE178360_RAW` |
| Input | 10x CellRanger HDF5 (v3 /matrix layout) |
| Species detection | from gene symbols in the matrix, not assumed -> **human** |
| Non-gene features removed | none |
| Gene space | inner join: 33538-36601 per sample -> 32732 shared |
| Metadata | none present in this series |
| QC thresholds | per sample, MAD-derived (see `analysis/GSE178360/qc/qc_thresholds.csv`) |
| Cells | 36,464 in -> 29,605 after QC -> 27,729 after doublet removal |
| Doublets | Scrublet per capture, 1,876 removed |
| Normalisation | normalize_total(target_sum=1e4) + log1p; the raw counts are preserved in processed/postQC.h5ad and re-attached as layers['counts'] of the final object |
| HVG | seurat_v3 on raw counts, n_top_genes=2500, batch_key=sample_id -> 2500 genes |
| Scaling | sc.pp.scale(max_value=10) on the HVG subset |
| PCA | 50 (variance-ratio elbow at 6, 90% cumulative variance at 51, bracketed to [20,50]) |
| Batch correction | Harmony (harmonypy) on sample_id used as the PRIMARY embedding, requested explicitly (--integration harmony). The justification recorded for this dataset is that single proposed cell types were fragmenting into sample-private clusters in the uncorrected embedding (see qc/celltype_split_by_sample.csv): one cell type is not several cell types in several donors. This series carries no condition metadata, so there is no designed experimental contrast that correcting on sample could destroy. The unintegrated X_pca and its UMAP are retained alongside it for comparison. |
| Neighbours | n_neighbors=15, n_pcs=50, use_rep=X_pca_harmony, metric=euclidean |
| UMAP | scanpy defaults min_dist=0.5, spread=1.0, random_state=0 |
| Clustering | Leiden (igraph flavour, 2 iterations) scanned at resolutions [0.3, 0.5, 0.8, 1.0]. Primary resolution: 1.0. It was the finest resolution at which no cluster falls below 20 cells, at least 90% of clusters carry >=5 specific marker genes against the rest of the data, and every cluster is also separated from its NEAREST neighbouring cluster by at least 10 genes (adj. p<0.05, |log2FC|>1, >25% expressing) - the last condition is what prevents simply picking whichever resolution yields the most clusters All resolutions are retained as leiden_res* columns. |
| Final clusters | **31** at resolution 1.0 |
| Marker test | Wilcoxon rank-sum (scanpy rank_genes_groups, one cluster vs all remaining cells) on log1p(CP10K) values, with expressing fractions; table re-used from the previous run of this same pipeline |
| Annotation | candidate identities scored per cluster from z-scored mean expression of curated panels, down-weighted by the fraction of cells expressing each marker; numeric Leiden labels are preserved and no candidate is promoted to a definitive label |
| Epithelial sub-analysis | performed on 7446 cells from clusters ['11', '12', '16', '18', '24', '25', '26', '29', '30']; 20 epithelial subclusters at Leiden resolution 0.6 |
| Condition DE | not performed |

Resolution scan (why that resolution was chosen):

| resolution | n_clusters | pct_clusters_marker_supported | min_DE_genes_between_nearest_clusters | nearest_cluster_pairs_under_threshold |
|---|---|---|---|---|
| 0.3 | 21.0 | 100.0 | 27.0 | 0.0 |
| 0.5 | 25.0 | 100.0 | 28.0 | 0.0 |
| 0.8 | 29.0 | 100.0 | 12.0 | 0.0 |
| 1.0 | 31.0 | 100.0 | 12.0 | 0.0 |

Batch mixing statistics actually measured:

```json
{
  "knn_same_sample_enrichment": 2.658,
  "mean_same_sample_neighbour_fraction": 0.8861,
  "expected_by_chance": 0.3367,
  "within_group_replicate_enrichment": null,
  "replicate_structure_available": false,
  "decision_rule": "No condition or group metadata exists for this series, so donor, batch and biology are completely confounded and no correction can be justified as purely technical. The unintegrated embedding is kept as primary; a Harmony embedding is computed alongside it, clearly labelled as supplementary, if samples separate strongly.",
  "integration_recommended": false,
  "supplementary_harmony_warranted": true
}
```

Full record: `analysis/GSE178360/logs/analysis_log.txt`, `analysis/GSE178360/logs/decisions.json`,
`analysis/GSE178360/README.md`.

## Execution order

```
01_scan_raw_data.py          recursive read-only scan of raw_data/ -> inventory
        |
run_scrna_analysis.py --dataset <GSE> --stages ...
        |
        +-- samples    per sample: load -> QC metrics -> MAD thresholds
        |               -> filter -> Scrublet -> .h5ad shard
        +-- merge      concat (barcodes prefixed) -> drop doublets
        |               -> filter genes -> join author metadata
        +-- norm       HVG on counts (seurat_v3) -> save counts checkpoint
        |               -> normalize_total + log1p
        +-- pca        scale HVGs -> randomized PCA -> batch assessment
        |               -> optional Harmony -> neighbours -> UMAP
        +-- cluster    Leiden resolution scan -> pick by size, marker support
        |               and nearest-neighbour separability
        +-- markers    Wilcoxon -> marker tables -> candidate annotation
        |               -> cross-check vs deposited author labels
        +-- figures    UMAPs, dot plots, feature plots, composition, QC review
        +-- epi        epithelial subset -> HVG -> PCA -> UMAP -> Leiden -> markers
        +-- finalize   final object (+counts layer), pseudobulk, analysis log
        |
03_write_report.py           README.md + final terminal summary
05_write_pipeline_as_run.py  this file
```

Stages are resumable: `--stages cluster,markers` re-enters from the saved
checkpoint. `--reuse-markers` re-uses the Wilcoxon table rather than
recomputing it. `--stages log` rewrites the analysis log from the persisted
decisions without recomputing anything.

## Environment actually used

Windows on ARM64. `numba`, `llvmlite` and `leidenalg` publish no ARM64 wheels
and no C compiler is installed, so the pipeline runs on an **x86-64 CPython
3.12** interpreter under emulation (`.venv-x64/`, created with `uv`).
`harmonypy` is pinned to the pure-Python `0.0.10` because newer releases need a
C++ toolchain. **No R, no Seurat, no Bioconductor** - which is why the R-only
tools in `docs/` could not have been used.

| Package | Version |
|---|---|
| `scanpy` | 1.12.3 |
| `anndata` | 0.13.2 |
| `numpy` | 2.4.6 |
| `scipy` | 1.18.0 |
| `pandas` | 3.0.5 |
| `matplotlib` | 3.11.1 |
| `scikit-learn` | 1.9.0 |
| `leidenalg` | 0.12.0 |
| `igraph` | 1.0.0 |
| `umap-learn` | 0.5.12 |
| `numba` | 0.66.0 |
| `harmonypy` | 0.0.10 |
| `scikit-image` | not installed |
| `scikit-misc` | 0.5.2 |
| `h5py` | 3.16.0 |
| `statsmodels` | 0.14.6 |
| `seaborn` | 0.13.2 |

Full lockfile: `analysis/requirements.txt`. Random seed fixed at 0 throughout
(NumPy, PCA, UMAP, Leiden, Scrublet).

## Reproducing

```bash
python analysis/scripts/01_scan_raw_data.py
python analysis/scripts/run_scrna_analysis.py --dataset GSE262927
python analysis/scripts/run_scrna_analysis.py --dataset GSE178360 --integration harmony
python analysis/scripts/03_write_report.py --dataset GSE262927
python analysis/scripts/05_write_pipeline_as_run.py
```
