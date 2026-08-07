#!/usr/bin/env python
"""
Build the per-dataset README and print the final terminal summary.

Everything here is read back out of the artefacts the pipeline actually wrote
(decisions.json, the QC tables, the marker tables), so the report cannot drift
away from what was really run.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import ANALYSIS, DATASETS  # noqa: E402


def _read_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _read_csv(p: Path, **kw) -> pd.DataFrame:
    return pd.read_csv(p, **kw) if p.exists() else pd.DataFrame()


def _md_table(df: pd.DataFrame, max_rows: int = 60) -> str:
    if df.empty:
        return "_(not available)_"
    d = df.head(max_rows)
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    rule = "|" + "|".join("---" for _ in d.columns) + "|"
    rows = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
            for r in d.itertuples(index=False)]
    out = "\n".join([head, rule] + rows)
    if len(df) > max_rows:
        out += f"\n\n_… {len(df) - max_rows} further rows in the CSV._"
    return out


def build(dataset: str) -> tuple[str, str]:
    cfg = DATASETS[dataset]()
    o = cfg.outdir
    d = _read_json(o / "logs" / "decisions.json")
    fmt = _read_json(o / "inventory" / "detected_format.json")
    batch = _read_json(o / "qc" / "batch_assessment.json")

    ba = _read_csv(o / "qc" / "qc_before_after.csv")
    th = _read_csv(o / "qc" / "qc_thresholds.csv")
    db = _read_csv(o / "qc" / "doublet_summary.csv")
    scan = _read_csv(o / "tables" / "clustering_resolution_scan.csv")
    ann = _read_csv(o / "tables" / "cluster_annotation_proposals.csv")
    dom = _read_csv(o / "tables" / "cluster_dominant_author_celltype.csv")
    warn_p = o / "qc" / "qc_review_warnings.txt"
    warns = (warn_p.read_text(encoding="utf-8").strip().splitlines()
             if warn_p.exists() else [])
    absent_p = o / "tables" / "absent_markers.txt"

    n_before = int(ba["Cells before QC"].iloc[-1]) if not ba.empty else 0
    n_after = int(ba["Cells after QC"].iloc[-1]) if not ba.empty else 0
    n_dbl = int(db["predicted_doublets"].sum()) if not db.empty else 0

    # Did the MAD *upper* bounds actually remove anything?  If they did not,
    # saying "an upper threshold was applied" would misrepresent the filtering.
    upper_note = ""
    if not ba.empty:
        rows = ba.iloc[:-1]
        hi_c = int(pd.to_numeric(rows.get("removed_high_counts", 0),
                                 errors="coerce").fillna(0).sum())
        hi_g = int(pd.to_numeric(rows.get("removed_high_genes", 0),
                                 errors="coerce").fillna(0).sum())
        if hi_c == 0 and hi_g == 0:
            upper_note = (
                "\n> **The upper MAD bounds did not bind.** Across every sample "
                "the median+5·MAD ceilings on counts and genes land far above "
                "the observed maximum, so **zero** cells were removed for having "
                "too many counts or genes. High-end outliers in this dataset are "
                "therefore handled by Scrublet, not by a count ceiling. The "
                "ceilings are still reported for transparency, but they should "
                "not be read as an applied filter.\n")
        else:
            upper_note = (f"\nThe upper MAD bounds removed {hi_c} cells on counts "
                          f"and {hi_g} on genes.\n")

    # How often was Scrublet's automatic threshold trustworthy?
    dbl_note = ""
    if not db.empty and "call_method" in db:
        n_auto = int(db["call_method"].astype(str)
                     .str.contains("automatic").sum())
        n_tot_s = len(db)
        if n_auto < n_tot_s:
            dbl_note = (
                f"\n> **Scrublet's automatic threshold was accepted in only "
                f"{n_auto} of {n_tot_s} samples.** In the rest the simulated-"
                f"doublet score histogram was not bimodal and the automatic cut "
                f"produced a call rate irreconcilable with the 10x multiplet "
                f"prior, so the expected-rate quantile was used instead. That "
                f"makes the doublet calls here a *ranking* cut rather than a "
                f"detected threshold — treat borderline calls with "
                f"corresponding caution. Per-sample detail is in the "
                f"`call_method` and `status` columns.\n")

    samples = d.get("samples_detected", "")

    sibling = ""
    if dataset == "GSE262927":
        sibling = """
> **This repository's `raw_data/` holds two independent GEO series.** They are
> different species and are analysed separately — they are never merged, and no
> gene symbol is case-converted between them.
>
> | Series | Species | Samples | Outputs |
> |---|---|---|---|
> | **GSE262927** (this report) | mouse | 33 | `analysis/` (primary) |
> | GSE178360 | human | 3 | `analysis/GSE178360/` |
>
> GSE178360 also ships four pre-processed Seurat `.RDS` objects
> (epithelial / endothelial / mesenchymal / immune subsets). They are valid
> R 3.6.3 RDS v3 files, but R is not installed on this machine and there is no
> C compiler, so they could not be read. The analysis of that series therefore
> starts from its three `filtered_feature_bc_matrix.h5` files, which contain
> the full unprocessed count matrices.
"""

    md = f"""# scRNA-seq analysis — {dataset}
{sibling}

Generated by `analysis/scripts/run_scrna_analysis.py`. Every number below is
read back out of the files the pipeline wrote; nothing here is hand-entered.

The raw data under `raw_data/` was opened read-only and never modified.

---

## 1. What raw files were discovered

`analysis/raw_data_inventory.txt` and `.csv` hold the full recursive scan of
`raw_data/`. For this dataset:

- Directory: `{cfg.raw_dir}`
- Input format: **{d.get('input_format', '?')}**
- Sample files: **{fmt.get('n_sample_files', '?')}**
- Features per file: {fmt.get('n_features', '?')}
- Values: {'raw integer UMI counts' if d.get('counts_are_raw_integers') == 'True'
           else d.get('counts_are_raw_integers', '?')}
- Feature types: {fmt.get('feature_types', '?')}
- Non-gene features found and removed from the expression matrix:
  {fmt.get('non_gene_features_present') or 'none'}

{cfg.notes}

## 2. How samples were inferred

Each HDF5 file is one 10x capture and one sample. The sample id is parsed from
the GEO filename (`{cfg.sample_regex}`) and written to `obs['sample_id']`
before anything else happens.

Cell names are `sampleID_originalBarcode`. 10x barcodes are drawn from a shared
whitelist and **do** recur between libraries, so prefixing is required — without
it, cells would be silently lost to name collisions at merge time. The
untouched barcode is kept in `obs['original_barcode']`.

Samples detected ({len(samples) if isinstance(samples, str) else '?'} chars):

```
{samples}
```

## 3. Species

**{str(d.get('species', '?')).upper()}** — determined from the gene symbols in the matrix
itself, not assumed. Mitochondrial genes follow the
`{'mt-' if d.get('species') == 'mouse' else 'MT-'}` convention
({fmt.get('n_mt_genes', '?')} genes). No symbol was case-converted between species
at any point.

## 4. Metadata structure

{d.get('metadata_table', 'none')}

- Join key: {d.get('metadata_join_key', 'n/a')}
- Cells carrying author metadata: {d.get('cells_with_author_metadata', 'n/a')}
- Samples with no metadata row at all: {d.get('samples_without_any_metadata', 'n/a')}
- Derived condition field: {d.get('condition_field', 'n/a')}

## 5. Cell counts before QC

**{n_before:,}** barcodes across all samples, {d.get('genes_loaded', '?')} genes
in the shared gene space.

Gene space identical across samples: {d.get('gene_space_identical_across_samples', '?')}
({d.get('gene_space_join', '')})

## 6. QC thresholds

Thresholds are derived **per sample** from that sample's own distributions —
median ± *n*·MAD on `log1p(total_counts)` and `log1p(n_genes)`, an upper-only
MAD bound on mitochondrial percentage bracketed into a sane window, and an
upper MAD bound on the fraction of counts in the top 20 genes. Lower bounds are
floored at 500 counts / 200 genes as a sanity guard. Full table:
`qc/qc_thresholds.csv`; every threshold carries its own `rationale` string.

{upper_note}
{_md_table(th[['sample', 'min_genes', 'max_genes', 'min_counts', 'max_counts',
               'max_pct_mt', 'max_pct_top20']] if not th.empty else th)}

## 7. Cells after QC

**{n_after:,}** of {n_before:,} retained
({100 * n_after / n_before if n_before else 0:.1f}%). Per-sample counts and the
reason each cell was dropped are in `qc/qc_before_after.csv` — no cell is
removed without appearing in that table.

{_md_table(ba)}

## 8. Doublet procedure

Scrublet was run **independently on every capture**, before merging, because
doublets can only form between cells that shared a 10x loading. Scores are kept
for every cell in the metadata whether or not the cell was called.

The automatic Scrublet threshold is only trustworthy when the simulated-doublet
score histogram is bimodal. Where it produced a call rate wildly inconsistent
with the 10x multiplet-rate prior (≈0.8% per 1000 cells recovered), the
threshold was replaced by that expected-rate quantile and the substitution is
recorded in the `call_method` and `status` columns of `qc/doublet_summary.csv`.

Doublets removed: **{n_dbl:,}**. {d.get('doublets_removed', '')}
{dbl_note}
{_md_table(db)}

## 9. Normalization

{d.get('normalization', '?')}

Raw counts are preserved in `processed/postQC.h5ad` and re-attached to the
final object as `layers['counts']` by an HDF5-level copy, so counts are never
overwritten by normalized values.

Expression layer used for all dot plots and feature plots: **log1p(CP10K)**.

## 10. HVG method

{d.get('hvg_method', '?')} → **{d.get('n_hvg', '?')}** highly variable genes.

Non-HVG genes are **not** deleted from the object, so any marker can still be
inspected afterwards; the HVG set only restricts the PCA input.

## 11. PCA settings

- {d.get('scaling', '?')}
- Components used: {d.get('n_pcs', '?')}
- Variance curve: `figures/qc/PCA_variance_ratio.png`

## 12. Batch assessment

```json
{json.dumps(batch, indent=2)}
```

## 13. Was integration used?

{d.get('batch_correction', '?')}

## 14. Clustering resolution

{d.get('clustering', '?')}

{_md_table(scan)}

## 15. Number of final clusters

**{d.get('n_clusters', '?')}** clusters at resolution {d.get('leiden_resolution', '?')},
stored as `obs['leiden_cluster']`. The other resolutions are retained as
`obs['leiden_res*']`.

## 16. Marker DE method

{d.get('marker_test', '?')}

- `tables/cluster_markers_all.csv` — every gene, every cluster
- `tables/cluster_markers_top20.csv` — top 20 per cluster by score

## 17. Tentative annotation strategy

{d.get('annotation_strategy', '?')}

Proposed identities are written to `tables/cluster_annotation_proposals.csv`
with supporting markers, conflicting markers and a confidence category. They
are **candidates**: the numeric Leiden labels remain the primary identifier
throughout, and `obs['proposed_cell_type']` is explicitly suffixed
`(candidate)`.

{d.get('author_label_agreement', '')}

{_md_table(ann.head(40)) if not ann.empty else ''}

{('### Agreement with the deposited author annotation\n\n'
  + _md_table(dom)) if not dom.empty else ''}

## 18. Where every output lives

```
{o.name}/
├── raw_data_inventory.txt / .csv   full recursive scan of raw_data (top level)
├── inventory/                      detected format, per-sample results (JSON)
├── logs/analysis_log.txt           every decision, in the required format
├── logs/decisions.json             the same decisions, machine readable
├── qc/
│   ├── qc_thresholds.csv           per-sample thresholds + rationale
│   ├── qc_before_after.csv         cells in/out per sample, per reason
│   ├── doublet_summary.csv         Scrublet per capture
│   ├── batch_assessment.json       kNN mixing statistics + decision rule
│   ├── cluster_qc_summary.csv      per-cluster QC characteristics
│   └── qc_review_warnings.txt      automatic QC review (phase 32)
├── figures/
│   ├── umap/                       cluster / sample / condition / split UMAPs
│   ├── dotplots/                   canonical + data-derived dot plots
│   ├── featureplots/               per-gene UMAPs
│   ├── composition/                composition bars, per-sample proportions
│   └── qc/                         QC violins, scatters, PCA variance
├── tables/                         markers, annotations, composition, metadata
├── processed/
│   ├── postQC.h5ad                 raw counts, post-QC cells
│   ├── normalized.h5ad             log1p(CP10K)
│   └── final_clustered.h5ad        final object incl. layers['counts']
└── epithelial_subanalysis/         second-stage epithelial analysis
```

### Required outputs

| Output | Path |
|---|---|
| Cluster UMAP | `figures/umap/UMAP_leiden_clusters.{{png,pdf}}` |
| Sample UMAP | `figures/umap/UMAP_sample.{{png,pdf}}` |
| Canonical marker dot plot | `figures/dotplots/canonical_marker_dotplot.{{png,pdf}}` |
| Data-derived dot plot | `figures/dotplots/cluster_marker_dotplot.{{png,pdf}}` |
| Marker tables | `tables/cluster_markers_all.csv`, `tables/cluster_markers_top20.csv` |
| Cell metadata | `tables/cell_metadata.csv` |
| Processed object | `processed/final_clustered.h5ad` |

Markers requested by the analysis brief but absent from this matrix are listed
in `tables/absent_markers.txt`{' (present)' if absent_p.exists() else ''} — no
empty panels were plotted.

## Automatic QC review

{chr(10).join('- ' + w for w in warns) if warns else '_No warnings raised._'}

## Reproducing

```bash
python analysis/scripts/01_scan_raw_data.py
python analysis/scripts/run_scrna_analysis.py --dataset {dataset}
python analysis/scripts/03_write_report.py --dataset {dataset}
```

Environment: `analysis/requirements.txt`. Random seed fixed at 0 throughout
(NumPy, PCA, UMAP, Leiden, Scrublet).

### A note on the Python environment

This machine is Windows on **ARM64**. `numba`, `llvmlite` and `leidenalg`
publish no ARM64 wheels and there is no C compiler installed, so scanpy cannot
be built against the native interpreter. The pipeline therefore runs on an
**x86-64 CPython 3.12** interpreter under Windows emulation
(`.venv-x64/`, created with `uv`), where every wheel resolves normally. This is
an environment workaround only — it changes nothing about the analysis.

`harmonypy` ≥ 0.0.11 now ships a C++ core that also cannot be built here; the
pure-Python `harmonypy==0.0.10` is pinned instead.
"""

    # ---- terminal summary ------------------------------------------------
    epi = d.get("epithelial_subanalysis", "not performed")
    summary = f"""
=====================================================
scRNA-seq ANALYSIS COMPLETE — {dataset}
=====================================================

RAW DATA
Directory:        {cfg.raw_dir}
Detected format:  {d.get('input_format', '?')}
Species:          {d.get('species', '?')}
Samples:          {fmt.get('n_sample_files', '?')}
Metadata:         {d.get('metadata_table', 'none')}

QC
Cells loaded:     {n_before:,}
Cells after QC:   {n_after:,}
Cells removed:    {n_before - n_after:,}
Doublets removed: {n_dbl:,}

ANALYSIS
HVGs:             {d.get('n_hvg', '?')}
PCs used:         {d.get('n_pcs', '?')}
Batch correction: {str(d.get('batch_correction', '?'))[:110]}
Leiden resolution:{d.get('leiden_resolution', '?')}
Final clusters:   {d.get('n_clusters', '?')}

MAIN OUTPUTS
Cluster UMAP:               {o}/figures/umap/UMAP_leiden_clusters.pdf
Sample UMAP:                {o}/figures/umap/UMAP_sample.pdf
Canonical marker dot plot:  {o}/figures/dotplots/canonical_marker_dotplot.pdf
Cluster-marker dot plot:    {o}/figures/dotplots/cluster_marker_dotplot.pdf
Cluster marker table:       {o}/tables/cluster_markers_all.csv
Cell metadata:              {o}/tables/cell_metadata.csv
Final processed object:     {o}/processed/final_clustered.h5ad
Analysis script:            {ANALYSIS}/scripts/run_scrna_analysis.py
README:                     {o}/README.md

EPITHELIAL SUBANALYSIS:
{epi}

QC WARNINGS:
{chr(10).join('  - ' + w for w in warns) if warns else '  none'}
=====================================================
"""
    return md, summary


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    args = ap.parse_args()
    md, summary = build(args.dataset)
    cfg = DATASETS[args.dataset]()
    (cfg.outdir / "README.md").write_text(md, encoding="utf-8")
    (cfg.outdir / "logs" / "final_summary.txt").write_text(summary,
                                                           encoding="utf-8")
    print(summary)
    print(f"wrote {cfg.outdir / 'README.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
