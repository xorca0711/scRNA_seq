# scRNA_seq

Notes and working documents for a single-cell RNA-seq analysis of lung regeneration after influenza injury.

The reference point is **Niethamer et al., *Cell Stem Cell* 2025** — a longitudinal atlas of 123,189 cells spanning uninjured lung through one year post-infection, which identifies transient injury states in the myeloid and epithelial compartments and one endothelial state (iCAP) that never resolves.

- **[`WORKFLOW.md`](WORKFLOW.md)** — end-to-end operational sequence, ordering
  constraints, subset-and-recluster loop, lineage-trace calling, and QC
  acceptance order.
- **[`docs/`](docs/README.md)** — one workflow schematic per tool: SoupX,
  Scrublet, scds, Slingshot, tradeSeq.
- **[`scRNAseq_workflow_Niethamer2025.md`](scRNAseq_workflow_Niethamer2025.md)** —
  the annotated pipeline reference: study design, stage-by-stage parameters,
  marker-gene annotation tables, and the twelve parameters the reference study
  leaves unspecified.
- **[`REFERENCES.md`](REFERENCES.md)** — all six papers with DOIs, PMC links,
  and software repositories.

---

## What has actually been run

The documents above are **notes and a plan**. An analysis has since been
executed, and it does **not** use every tool listed here.

- **[`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md)** — the authoritative
  record of what was executed: real parameters, cell counts, batch decisions,
  which tools were used, and how it diverges from the source papers. Generated
  from the pipeline's own outputs.
- **[`analysis/README.md`](analysis/README.md)** — full report for the primary
  dataset. **[`analysis/GSE178360/README.md`](analysis/GSE178360/README.md)** —
  the second dataset.
- **[`analysis/scripts/`](analysis/scripts/)** — the pipeline itself.

| | Executed |
|---|---|
| Framework | **Python / scanpy** (no R, no Seurat — R is unavailable in this environment) |
| Datasets | **GSE262927** mouse, 33 samples, 162,175 cells, 29 clusters · **GSE178360** human, 3 samples, 27,729 cells, 31 clusters |
| Of the six tools below | only **Scrublet** was used. SoupX, scds, Slingshot and tradeSeq were **not**. |
| Not done | ambient-RNA correction, trajectory/pseudotime, trajectory DE, cell-cycle regression |

`raw_data/` is never modified and is not committed.

---

## 1. Toolbox

Six tools, each owning one decision the others cannot make.

| Stage | Tool | Owns | Schematic |
|---|---|---|---|
| Alignment | STARsolo 2.7.9a | Barcode, UMI, and gene assignment against mm39 | — |
| Ambient RNA | SoupX 1.6.0 | How much of each count is cell-free background | [`docs/SOUPX.md`](docs/SOUPX.md) |
| Doublets | Scrublet | Whether a barcode looks like a simulated cell pair | [`docs/SCRUBLET.md`](docs/SCRUBLET.md) |
| Doublets | scds | Whether a barcode co-expresses genes that rarely co-occur | [`docs/SCDS.md`](docs/SCDS.md) |
| Cell state | Seurat 4.9 | Normalization, clustering, annotation, marker DE | — |
| Trajectory | Slingshot | Lineage topology and pseudotime ordering | [`docs/SLINGSHOT.md`](docs/SLINGSHOT.md) |
| Trajectory DE | tradeSeq | Which genes change, and in what sense | [`docs/TRADESEQ.md`](docs/TRADESEQ.md) |

Neither doublet caller supersedes the other. They fail differently, which is
the reason for running both — and in the reference study, doublet score alone
decided neither of the two ambiguous populations.

---

## 2. The pipeline in brief

```
FASTQ
  └─ STARsolo 2.7.9a (mm39)          alignment, UMI counting
     └─ SoupX 1.6.0                  ambient RNA removal        │ per library
        └─ scds + Scrublet           doublet scoring            │ pre-merge
           └─ merge within cohort    + ROSA26 trace annotation
              └─ QC filters          %mito ≤15 · nFeature ≥50 · ≤2×MAD
                 └─ CellCycleScoring
                    └─ SCTransform   regress %mt, nFeature, nCount
                       └─ PCA → Louvain (res 1.0) → UMAP
                          └─ marker annotation + cluster culling
                             │
                             ├─ subset & recluster (res 0.4–1.2)
                             ├─ lineage-trace calling
                             ├─ module scoring (ISG, iCAP)
                             ├─ CellChat 1.6
                             ├─ clusterProfiler gseGO
                             └─ Slingshot → tradeSeq
```

Three orderings matter and are easy to get wrong: ambient correction runs **before** doublet calling, both run **per library before merging**, and doublet calling runs **before** count-based filtering — doublets sit in the high-`nFeature` tail you would otherwise cut.

## 3. Data

Public data for the target paper: GEO [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927).

No sequencing data, count matrices or PDFs are stored in this repository.

## 4. Environment

The pipeline is bilingual — Scrublet is Python, everything else is R.

- R ≥ 4.3, Bioconductor: `scds`, `slingshot`, `tradeSeq`, `ComplexHeatmap`, `clusterProfiler`, `enrichplot`, `DOSE`; CRAN: `Seurat`, `SoupX`, `CellChat`, `MetBrewer`
- Python: `scrublet`, `scanpy` (for I/O)

Seurat v4.9 is a pre-release snapshot and is not on CRAN. Pin `Seurat 4.4.0` (last stable v4) or move to v5 and document the change — the Assay5 object model alters `SCTransform` and `FindMarkers` behaviour around layers.

## 5. Licence

Written material in this repository is © the author. The papers it describes are the property of their respective publishers; see [`REFERENCES.md`](REFERENCES.md) for links to the open-access versions.
