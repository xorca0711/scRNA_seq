# scRNA_seq

Notes and working documents for a single-cell RNA-seq analysis of lung regeneration after influenza injury.

The reference point is **Niethamer et al., *Cell Stem Cell* 2025** — a longitudinal atlas of 123,189 cells spanning uninjured lung through one year post-infection, which identifies transient injury states in the myeloid and epithelial compartments and one endothelial state (iCAP) that never resolves.

## Contents

| File | What it is |
|---|---|
| [`scRNAseq_workflow_Niethamer2025.md`](scRNAseq_workflow_Niethamer2025.md) | The analysis pipeline organised stage by stage — parameters, ordering constraints, downstream modules, marker-gene annotation reference, and the parameters the paper leaves unspecified |
| [`REFERENCES.md`](REFERENCES.md) | The target paper and the five method papers, with DOIs, PMC links and software repositories |

## The pipeline in brief

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

## Data

Public data for the target paper: GEO [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927).

No sequencing data, count matrices or PDFs are stored in this repository.

## Environment

The pipeline is bilingual — Scrublet is Python, everything else is R.

- R ≥ 4.3, Bioconductor: `scds`, `slingshot`, `tradeSeq`, `ComplexHeatmap`, `clusterProfiler`, `enrichplot`, `DOSE`; CRAN: `Seurat`, `SoupX`, `CellChat`, `MetBrewer`
- Python: `scrublet`, `scanpy` (for I/O)

Seurat v4.9 is a pre-release snapshot and is not on CRAN. Pin `Seurat 4.4.0` (last stable v4) or move to v5 and document the change — the Assay5 object model alters `SCTransform` and `FindMarkers` behaviour around layers.

## Licence

Written material in this repository is © the author. The papers it describes are the property of their respective publishers; see [`REFERENCES.md`](REFERENCES.md) for links to the open-access versions.
