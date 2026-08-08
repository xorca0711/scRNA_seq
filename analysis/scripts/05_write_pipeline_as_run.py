#!/usr/bin/env python
"""
Emit docs/PIPELINE_AS_RUN.md - the authoritative record of the workflow and
parameters that were actually executed.

Everything is read back out of the artefacts the pipeline wrote
(logs/decisions.json, the QC tables, requirements.txt), so this document cannot
drift away from what really happened.  Re-run it after any analysis run.

The tool reference pages under docs/ are deliberately untouched: they describe
tools worth knowing about, not necessarily tools this run used.  Which of them
were used is stated explicitly below.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import ANALYSIS, DATASETS, REPO  # noqa: E402

# Tools the repo keeps reference pages for, and whether this run used them.
# Paths are relative to docs/, where the generated file lives.
TOOL_STATUS = [
    ("Scrublet", "SCRUBLET.md", "USED",
     "run per capture via `sc.pp.scrublet`, before merging"),
    ("SoupX", "SOUPX.md", "NOT USED",
     "no ambient-RNA correction was applied to either dataset"),
    ("scds", "SCDS.md", "NOT USED",
     "Scrublet was the only doublet caller; scds is R-only and R is unavailable here"),
    ("Slingshot", "SLINGSHOT.md", "NOT USED",
     "not run in the main pipeline; the focused alveolar re-analysis "
     "(analysis/scripts/06_regeneration_focus.py) instead uses scanpy's "
     "PAGA + diffusion pseudotime for the AT2 -> AT1 trajectory"),
    ("tradeSeq", "TRADESEQ.md", "NOT USED",
     "no trajectory-based DE was performed with any tool"),
]

# Divergences from the two source publications, established by reading them
# (Thesis/Primary/). These are stated here because the deposited files alone
# do not reveal them.
DIVERGENCES = """
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
| Trajectory | **Slingshot + tradeSeq** | **PAGA + diffusion pseudotime**, in the focused alveolar re-analysis only |
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

"""

NOT_PERFORMED = [
    ("Ambient RNA / soup correction", "SoupX, CellBender, decontX",
     "The raw (unfiltered) droplet matrices needed for SoupX are present for "
     "GSE178360 but not for GSE262927, whose deposited files are cell-called "
     "only. No correction was applied to either dataset, so ambient "
     "contamination is uncorrected - this matters most for GSE178360 sample "
     "DD046Q, which carries heavy haemoglobin signal."),
    ("Trajectory / pseudotime", "Slingshot, PAGA, DPT, Monocle",
     "Not run in the MAIN pipeline: the whole-atlas and epithelial "
     "sub-analysis outputs resolve the transitional states as clusters only. "
     "It WAS run in the focused re-analysis "
     "(analysis/scripts/06_regeneration_focus.py): PAGA topology plus "
     "diffusion pseudotime rooted in AT2, on the 25-sample annotated cohort. "
     "See analysis/GSE262927/regeneration_focus/."),
    ("Trajectory-based differential expression", "tradeSeq",
     "Not run; marker programmes are summarised along pseudotime bins in the "
     "focused re-analysis, but no formal trajectory-DE model was fitted."),
    ("RNA velocity", "scVelo, velocyto",
     "Not run; requires spliced/unspliced counts that the deposited files do "
     "not contain."),
    ("Cell-cycle scoring or regression", "Seurat CellCycleScoring",
     "Not computed. `cell_cycle_phase` in the metadata was READ from the "
     "deposited GSE262927 author table, never recalculated, and was never "
     "regressed out."),
    ("SCTransform / covariate regression", "Seurat SCTransform",
     "Not used; normalisation was normalize_total + log1p."),
    ("Reference mapping / label transfer", "scVI, scANVI, Seurat anchors",
     "Not run. Author cell-type labels were used only as an independent "
     "check of the clustering, never to assign identities."),
    ("Pathway / gene-set enrichment", "GSEA, GO, fgsea",
     "Not run."),
    ("Reading the deposited Seurat objects", "R, Seurat",
     "The four GSE178360 `.RDS` objects were never opened. They are valid R "
     "3.6.3 RDS v3 files, but R is not installed and there is no C compiler, "
     "so the analysis started from the `filtered_feature_bc_matrix.h5` files."),
]


def _json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def _csv(p: Path) -> pd.DataFrame:
    return pd.read_csv(p) if p.exists() else pd.DataFrame()


def dataset_block(name: str) -> str:
    cfg = DATASETS[name]()
    o = cfg.outdir
    d = _json(o / "logs" / "decisions.json")
    fmt = _json(o / "inventory" / "detected_format.json")
    batch = _json(o / "qc" / "batch_assessment.json")
    ba = _csv(o / "qc" / "qc_before_after.csv")
    db = _csv(o / "qc" / "doublet_summary.csv")
    scan = _csv(o / "tables" / "clustering_resolution_scan.csv")

    n_before = int(ba["Cells before QC"].iloc[-1]) if not ba.empty else 0
    n_after = int(ba["Cells after QC"].iloc[-1]) if not ba.empty else 0
    n_dbl = int(db["predicted_doublets"].sum()) if not db.empty else 0
    rel = o.relative_to(REPO).as_posix()

    # Paths recorded at runtime are absolute and carry the local username;
    # show them relative to the repo instead.
    def _rel(v) -> str:
        s = str(v)
        try:
            return Path(s).relative_to(REPO).as_posix() if Path(s).is_absolute() else s
        except (ValueError, OSError):
            # Recorded on a machine where the repo lived elsewhere: rebase at
            # the first repo-level directory name found in the path.
            parts = Path(s).parts
            for anchor in ("raw_data", "analysis", "docs"):
                if anchor in parts:
                    return Path(*parts[parts.index(anchor):]).as_posix()
            return s

    meta_txt = _rel(d.get("metadata_table", "none"))
    raw_txt = _rel(cfg.raw_dir)

    scan_md = ""
    if not scan.empty:
        cols = ["resolution", "n_clusters", "pct_clusters_marker_supported",
                "min_DE_genes_between_nearest_clusters",
                "nearest_cluster_pairs_under_threshold"]
        cols = [c for c in cols if c in scan.columns]
        scan_md = ("\n| " + " | ".join(cols) + " |\n|"
                   + "|".join("---" for _ in cols) + "|\n"
                   + "\n".join("| " + " | ".join(str(r[c]) for c in cols) + " |"
                               for _, r in scan.iterrows()) + "\n")

    return f"""
### {name} ({d.get('species', '?')}, {fmt.get('n_sample_files', '?')} samples)

| Step | What was actually run |
|---|---|
| Input directory | `{raw_txt}` |
| Input | {d.get('input_format', '?')} |
| Species detection | from gene symbols in the matrix, not assumed -> **{d.get('species', '?')}** |
| Non-gene features removed | {fmt.get('non_gene_features_present') or 'none'} |
| Gene space | {d.get('gene_space_join', 'n/a')} |
| Metadata | {meta_txt} |
| QC thresholds | per sample, MAD-derived (see `{rel}/qc/qc_thresholds.csv`) |
| Cells | {n_before:,} in -> {n_after:,} after QC -> {n_after - n_dbl:,} after doublet removal |
| Doublets | Scrublet per capture, {n_dbl:,} removed |
| Normalisation | {d.get('normalization', '?')} |
| HVG | {d.get('hvg_method', '?')} -> {d.get('n_hvg', '?')} genes |
| Scaling | {d.get('scaling', '?')} |
| PCA | {d.get('n_pcs', '?')} |
| Batch correction | {d.get('batch_correction', '?')} |
| Neighbours | {d.get('neighbors', '?')} |
| UMAP | {d.get('umap', '?')} |
| Clustering | {d.get('clustering', '?')} |
| Final clusters | **{d.get('n_clusters', '?')}** at resolution {d.get('leiden_resolution', '?')} |
| Marker test | {d.get('marker_test', '?')} |
| Annotation | {d.get('annotation_strategy', '?')} |
| Epithelial sub-analysis | {d.get('epithelial_subanalysis', 'not performed')} |
| Condition DE | {d.get('condition_differential_expression', 'not performed')} |

Resolution scan (why that resolution was chosen):
{scan_md}
Batch mixing statistics actually measured:

```json
{json.dumps(batch, indent=2)}
```

Full record: `{rel}/logs/analysis_log.txt`, `{rel}/logs/decisions.json`,
`{rel}/README.md`.
"""


def main() -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass

    req = (ANALYSIS / "requirements.txt")
    pkgs = {}
    if req.exists():
        for line in req.read_text(encoding="utf-8").splitlines():
            if "==" in line:
                k, v = line.split("==", 1)
                pkgs[k.strip().lower()] = v.strip()
    key_pkgs = ["scanpy", "anndata", "numpy", "scipy", "pandas", "matplotlib",
                "scikit-learn", "leidenalg", "igraph", "umap-learn", "numba",
                "harmonypy", "scikit-image", "scikit-misc", "h5py",
                "statsmodels", "seaborn"]

    tool_rows = "\n".join(
        f"| [{n}]({p}) | **{s}** | {why} |" for n, p, s, why in TOOL_STATUS)
    not_done = "\n".join(
        f"### {what}\n\n*Reference tools: {tools}*\n\n{why}\n"
        for what, tools, why in NOT_PERFORMED)
    env_rows = "\n".join(f"| `{k}` | {pkgs.get(k, 'not installed')} |"
                         for k in key_pkgs)

    md = f"""# Pipeline as run

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
{tool_rows}

## What was NOT done

This section exists so that nobody reading `docs/` concludes that these steps
are represented in the results. They are not.

{not_done}
{DIVERGENCES}
## The workflow that was executed
{dataset_block('GSE262927')}
{dataset_block('GSE178360')}
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
{env_rows}

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
"""

    out = REPO / "docs" / "PIPELINE_AS_RUN.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(md, encoding="utf-8")
    print(f"wrote {out} ({len(md.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
