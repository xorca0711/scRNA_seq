#!/usr/bin/env python
"""
=============================================================================
Adaptive scRNA-seq pipeline - UMAPs, marker dot plots, cluster markers, QC
=============================================================================

The workflow is *not* fixed in advance.  `01_scan_raw_data.py` inspects the
raw-data tree first; the configuration in `pipeline_utils.py` encodes what that
scan actually found, and this script adapts to it (species, gene-space
reconciliation, presence or absence of a metadata table, lineage-reporter
contigs that must leave the expression matrix, and so on).

The raw-data directory is read-only throughout.

Run one dataset end to end:

    python run_scrna_analysis.py --dataset GSE262927

or resume from a checkpoint:

    python run_scrna_analysis.py --dataset GSE262927 --stages cluster,markers

Stages: samples, merge, norm, pca, cluster, markers, figures, epi, finalize
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))

import markers as MK  # noqa: E402
from pipeline_utils import (ANALYSIS, DATASETS, RANDOM_SEED, REPO, DatasetConfig,  # noqa: E402
                            ensure_dirs, figure_of, free_mem, log, mem_report,
                            sample_number, save_fig, setup_matplotlib)

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

plt = setup_matplotlib()
import matplotlib.pyplot as _plt  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402

sc.settings.verbosity = 1
sc.settings.n_jobs = max(1, (os.cpu_count() or 4) - 2)
np.random.seed(RANDOM_SEED)

DECISIONS: dict[str, object] = {}
_DECISIONS_PATH: Path | None = None


def _repo_rel(path) -> str:
    """Render current or previously recorded paths without checkout details."""
    value = str(path).replace("\\", "/")
    repo_prefix = REPO.as_posix().rstrip("/") + "/"
    if value.lower().startswith(repo_prefix.lower()):
        return value[len(repo_prefix):]
    parts = [part for part in value.split("/") if part]
    for anchor in ("raw_data", "analysis", "docs"):
        if anchor in parts:
            return "/".join(parts[parts.index(anchor):])
    return value


def load_decisions(path: Path) -> None:
    """Carry decisions forward across staged runs.

    The pipeline can be resumed one stage at a time, and each stage only
    records its own choices.  Without this, finishing on a late stage would
    overwrite the record with a handful of entries and lose everything the
    earlier stages decided.
    """
    global _DECISIONS_PATH
    _DECISIONS_PATH = path
    if path.exists():
        try:
            prior = json.loads(path.read_text(encoding="utf-8"))
            DECISIONS.update(prior)
            log(f"carried forward {len(prior)} decision(s) from a previous stage")
        except Exception as exc:  # noqa: BLE001
            log(f"could not read prior decisions ({exc}); starting fresh")


def save_decisions() -> None:
    if _DECISIONS_PATH is None:
        return
    _DECISIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
    _DECISIONS_PATH.write_text(
        json.dumps({k: str(v) for k, v in DECISIONS.items()}, indent=2),
        encoding="utf-8")


def record(key: str, value) -> None:
    """Every analytical choice that matters ends up in the log and README."""
    DECISIONS[key] = value
    log(f"  DECISION {key} = {value}")
    save_decisions()          # survive a crash in a later stage


# =============================================================================
# PHASE 1 - discovery
# =============================================================================
def scan_input_directory(cfg: DatasetConfig) -> list[dict]:
    """List the sample-level inputs actually present for this dataset."""
    import re
    files = sorted(cfg.raw_dir.glob(cfg.h5_glob))
    tasks = []
    for p in files:
        m = re.match(cfg.sample_regex, p.name)
        if not m:
            log(f"  skipping unrecognised filename: {p.name}")
            continue
        sample = m.group(1)
        gsm = p.name.split("_")[0]
        tasks.append({"path": str(p), "sample": sample, "gsm": gsm})
    log(f"scan_input_directory: {len(tasks)} sample files under {cfg.raw_dir}")
    return tasks


def detect_input_format(cfg: DatasetConfig, tasks: list[dict]) -> dict:
    """Confirm the on-disk format and the species from the data themselves."""
    import h5py
    info: dict = {"n_sample_files": len(tasks)}
    if not tasks:
        return info
    with h5py.File(tasks[0]["path"], "r") as f:
        top = list(f.keys())
        info["hdf5_top_level"] = top
        info["format"] = ("10x CellRanger HDF5 (v3 /matrix layout)"
                          if "matrix" in top else f"unrecognised: {top}")
        g = f["matrix"]
        info["shape_genes_x_cells"] = [int(x) for x in g["shape"][:]]
        info["data_dtype"] = str(g["data"].dtype)
        probe = np.asarray(g["data"][:200000])
        info["values_are_integers"] = bool(np.all(probe == np.floor(probe)))
        info["explicit_zeros_in_probe"] = int((probe == 0).sum())
        names = np.array([b.decode() for b in g["features"]["name"][:]])
        info["n_features"] = int(names.size)
        sset = set(names)
        n_mouse = sum(1 for s in ["Sftpc", "Krt8", "Trp63", "Pecam1", "Ptprc"]
                      if s in sset)
        n_human = sum(1 for s in ["SFTPC", "KRT8", "TP63", "PECAM1", "PTPRC"]
                      if s in sset)
        info["species_detected"] = "mouse" if n_mouse > n_human else "human"
        info["n_mt_genes"] = int(sum(1 for s in names if s.startswith(
            MK.mito_prefix(info["species_detected"]))))
        if "feature_type" in g["features"]:
            ft, ct = np.unique(
                [b.decode() for b in g["features"]["feature_type"][:]],
                return_counts=True)
            info["feature_types"] = dict(zip(ft.tolist(), ct.tolist()))
        info["non_gene_features_present"] = [
            s for s in cfg.non_gene_features if s in sset]
    if info["species_detected"] != cfg.species:
        raise SystemExit(
            f"Species detected from the matrix ({info['species_detected']}) "
            f"does not match the configured species ({cfg.species}). Refusing "
            f"to continue rather than analyse with the wrong marker panel.")
    record("input_format", info["format"])
    record("species", info["species_detected"])
    record("counts_are_raw_integers", info["values_are_integers"])
    return info


# =============================================================================
# PHASE 4-8 - per-sample load, QC, doublets
# =============================================================================
def load_samples(cfg: DatasetConfig, tasks: list[dict], dirs: dict,
                 workers: int) -> pd.DataFrame:
    from sample_worker import _worker_entry

    payload = [dict(t, species=cfg.species, shard_dir=str(dirs["shards"]),
                    non_gene_features=cfg.non_gene_features,
                    use_gene_ids_as_index=cfg.use_gene_ids_as_index)
               for t in tasks]

    results: list[dict] = []
    log(f"load_samples: {len(payload)} samples on {workers} worker process(es)")
    if workers <= 1:
        for t in payload:
            results.append(json.loads(_worker_entry(t)))
            log(f"  done {results[-1].get('sample')} {mem_report()}")
    else:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_worker_entry, t): t["sample"] for t in payload}
            for fut in as_completed(futs):
                r = json.loads(fut.result())
                results.append(r)
                if "error" in r:
                    log(f"  ERROR {r.get('sample')}: {r['error']}")
                else:
                    log(f"  done {r['sample']}: "
                        f"{r['n_cells_before_qc']} -> {r['n_cells_after_qc']} cells, "
                        f"{r['n_predicted_doublets']} doublets")

    errs = [r for r in results if "error" in r]
    if errs:
        for e in errs:
            log(f"FAILED {e['sample']}: {e.get('traceback', '')[:2000]}")
        raise SystemExit(f"{len(errs)} sample(s) failed; refusing to continue "
                         f"with a silently incomplete dataset.")

    results.sort(key=lambda r: r["sample"])
    (dirs["inventory"] / "per_sample_results.json").write_text(
        json.dumps(results, indent=2, default=str), encoding="utf-8")

    # --- threshold table (phase 7) ---------------------------------------
    th = pd.DataFrame([r["thresholds"] for r in results])
    th.to_csv(dirs["qc"] / "qc_thresholds.csv", index=False)

    # --- before/after table (phase 7) ------------------------------------
    ba = pd.DataFrame([{
        "Sample": r["sample"],
        "GSM": r.get("gsm", ""),
        "Cells before QC": r["n_cells_before_qc"],
        "Cells after QC": r["n_cells_after_qc"],
        "Cells removed": r["n_cells_before_qc"] - r["n_cells_after_qc"],
        "% retained": round(100 * r["n_cells_after_qc"] / r["n_cells_before_qc"], 2),
        **{f"removed_{k}": v for k, v in r["removal_reasons"].items()},
    } for r in results])
    ba.loc["TOTAL"] = ["TOTAL", "", ba["Cells before QC"].sum(),
                       ba["Cells after QC"].sum(), ba["Cells removed"].sum(),
                       round(100 * ba["Cells after QC"].sum()
                             / ba["Cells before QC"].sum(), 2)] + \
                      [ba[c].sum() for c in ba.columns if c.startswith("removed_")]
    ba.to_csv(dirs["qc"] / "qc_before_after.csv", index=False)

    # --- doublet summary (phase 8) ---------------------------------------
    db = pd.DataFrame([{
        "sample": r["sample"],
        "cells_after_qc": r["n_cells_after_qc"],
        "predicted_doublets": r["n_predicted_doublets"],
        "doublet_rate_pct": round(100 * r["doublet_rate"], 3),
        "expected_doublet_rate_pct": round(100 * r.get("expected_doublet_rate",
                                                       float("nan")), 3),
        "scrublet_threshold": r["doublet_threshold"],
        "call_method": r.get("doublet_call_method", ""),
        "status": r["doublet_status"],
    } for r in results])
    db.to_csv(dirs["qc"] / "doublet_summary.csv", index=False)

    log(f"load_samples: {int(ba.loc['TOTAL', 'Cells before QC'])} barcodes in, "
        f"{int(ba.loc['TOTAL', 'Cells after QC'])} pass QC")
    return pd.DataFrame(results)


def calculate_qc(dirs: dict, cfg: DatasetConfig) -> pd.DataFrame:
    """Collect the per-sample pre-QC cell tables and draw the QC figures."""
    frames = []
    for p in sorted(dirs["shards"].glob("*__preQC_obs.*")):
        frames.append(pd.read_parquet(p) if p.suffix == ".parquet"
                      else pd.read_csv(p))
    pre = pd.concat(frames, ignore_index=True)
    out = dirs["qc"] / "merged_preQC_cell_metrics.csv.gz"
    pre.to_csv(out, index=False, compression="gzip")
    log(f"calculate_qc: {len(pre)} pre-QC barcodes -> {out.name}")

    order = sorted(pre["sample_id"].unique())
    _qc_violin(pre, order, dirs)
    _qc_scatter(pre, dirs)
    return pre


def _qc_violin(pre: pd.DataFrame, order: list[str], dirs: dict) -> None:
    metrics = [("total_counts", "UMI counts per cell", True),
               ("n_genes_by_counts", "Genes per cell", True),
               ("pct_counts_mt", "Mitochondrial %", False),
               ("pct_counts_ribo", "Ribosomal %", False)]
    n = len(order)
    fig, axes = _plt.subplots(len(metrics), 1,
                              figsize=(max(8, 0.36 * n + 2), 3.0 * len(metrics)),
                              sharex=True)
    for ax, (col, label, logscale) in zip(axes, metrics):
        data = [pre.loc[pre["sample_id"] == s, col].to_numpy() for s in order]
        parts = ax.violinplot(data, showextrema=False, widths=0.85)
        for pc in parts["bodies"]:
            pc.set_facecolor("#4C72B0")
            pc.set_alpha(0.75)
            pc.set_edgecolor("none")
        med = [np.median(d) for d in data]
        ax.scatter(range(1, n + 1), med, s=6, color="white", zorder=3)
        if logscale:
            ax.set_yscale("log")
        ax.set_ylabel(label)
        ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    axes[-1].set_xticks(range(1, n + 1))
    axes[-1].set_xticklabels(order, rotation=90, fontsize=6)
    axes[0].set_title("Per-sample QC distributions (before filtering)")
    fig.tight_layout()
    save_fig(fig, dirs["fig_qc"], "QC_violin_per_sample")


def _qc_scatter(pre: pd.DataFrame, dirs: dict) -> None:
    sub = pre.sample(min(60000, len(pre)), random_state=RANDOM_SEED)
    fig, axes = _plt.subplots(1, 3, figsize=(15, 4.4))
    sc0 = axes[0].scatter(sub["total_counts"], sub["n_genes_by_counts"],
                          c=sub["pct_counts_mt"], s=1.2, cmap="viridis",
                          vmin=0, vmax=25, rasterized=True)
    axes[0].set_xscale("log"); axes[0].set_yscale("log")
    axes[0].set_xlabel("UMI counts"); axes[0].set_ylabel("Genes detected")
    axes[0].set_title("Counts vs genes (colour = % mito)")
    fig.colorbar(sc0, ax=axes[0], label="% mito")

    axes[1].scatter(sub["total_counts"], sub["pct_counts_mt"], s=1.2,
                    alpha=0.3, color="#C44E52", rasterized=True)
    axes[1].set_xscale("log")
    axes[1].set_xlabel("UMI counts"); axes[1].set_ylabel("% mitochondrial")
    axes[1].set_title("Mitochondrial fraction vs depth")

    ok = sub["passes_qc"].astype(bool)
    axes[2].scatter(sub.loc[~ok, "total_counts"], sub.loc[~ok, "n_genes_by_counts"],
                    s=1.2, alpha=0.5, color="#BBBBBB", label="removed",
                    rasterized=True)
    axes[2].scatter(sub.loc[ok, "total_counts"], sub.loc[ok, "n_genes_by_counts"],
                    s=1.2, alpha=0.5, color="#55A868", label="kept",
                    rasterized=True)
    axes[2].set_xscale("log"); axes[2].set_yscale("log")
    axes[2].set_xlabel("UMI counts"); axes[2].set_ylabel("Genes detected")
    axes[2].set_title("Cells kept vs removed by QC")
    axes[2].legend(markerscale=8)
    fig.tight_layout()
    save_fig(fig, dirs["fig_qc"], "QC_scatter_global")


# =============================================================================
# PHASE 9 - merge
# =============================================================================
def merge_samples(cfg: DatasetConfig, dirs: dict) -> ad.AnnData:
    shards = sorted(dirs["shards"].glob("*.h5ad"))
    if not shards:
        raise SystemExit("no sample shards found - run the 'samples' stage first")
    log(f"merge_samples: concatenating {len(shards)} shards {mem_report()}")

    adatas = {}
    for p in shards:
        a = ad.read_h5ad(p)
        adatas[p.stem] = a
    var_sets = [frozenset(a.var_names) for a in adatas.values()]
    identical = all(v == var_sets[0] for v in var_sets)
    record("gene_space_identical_across_samples", identical)
    join = "inner"
    merged = ad.concat(adatas, axis=0, join=join, label=None,
                       index_unique=None, merge="first")
    n_var_before = [a.n_vars for a in adatas.values()]
    record("gene_space_join", f"{join} join: "
                              f"{min(n_var_before)}-{max(n_var_before)} per sample "
                              f"-> {merged.n_vars} shared")
    del adatas
    free_mem()

    if cfg.use_gene_ids_as_index:
        # Samples were keyed on Ensembl ID so that references with different
        # annotation builds could be intersected safely.  Now that the gene
        # space is shared, switch the index to symbols (made unique) because
        # every marker panel is expressed in symbols.
        from pipeline_utils import make_unique as _mu
        merged.var["gene_id"] = merged.var_names.to_numpy()
        merged.var_names = pd.Index(
            _mu(merged.var["gene_symbol"].astype(str).to_numpy()).astype(str))
        record("var_index",
               "samples were intersected on Ensembl gene ID (the references "
               "differ between samples), then the index was switched to "
               "make.unique'd gene symbols for marker lookup")

    if merged.obs_names.duplicated().any():
        raise SystemExit("duplicate cell names after concatenation - the "
                         "sample-prefixing step failed")
    log(f"merge_samples: {merged.n_obs} cells x {merged.n_vars} genes "
        f"{mem_report()}")
    return merged


def merge_metadata(cfg: DatasetConfig, adata: ad.AnnData, dirs: dict) -> ad.AnnData:
    """Attach the author-supplied per-cell metadata, where one exists."""
    if cfg.metadata_csv is None or not cfg.metadata_csv.exists():
        record("metadata_table", "none present in this series")
        adata.obs["has_author_metadata"] = False
        return adata

    md = pd.read_csv(cfg.metadata_csv)
    log(f"merge_metadata: {len(md)} annotated cells, columns {list(md.columns)}")

    # the leading column is an empty write.csv rownames artefact
    empty = [c for c in md.columns if md[c].isna().all()]
    if empty:
        md = md.drop(columns=empty)
        record("metadata_empty_columns_dropped", empty)
    # celltype and subtype were verified identical during the scan
    if "celltype" in md and "subtype" in md and md["celltype"].equals(md["subtype"]):
        md = md.drop(columns=["subtype"])
        record("metadata_redundant_columns_dropped", ["subtype (identical to celltype)"])

    md["sample_num"] = md["orig.ident"].map(sample_number)
    if md["sample_num"].isna().any():
        raise SystemExit("could not parse a sample number out of every "
                         "orig.ident value")
    md["join_key"] = md["sample_num"] + "|" + md["cb"].astype(str)
    if md["join_key"].duplicated().any():
        raise SystemExit("metadata join key (sample, barcode) is not unique")

    obs = adata.obs
    key = (pd.Series(obs["sample_id"].astype(str).map(sample_number).to_numpy(),
                     index=obs.index)
           + "|" + obs["original_barcode"].astype(str))
    md = md.set_index("join_key")

    rename = {"orig.ident": "orig_ident", "nCount_RNA": "author_nCount_RNA",
              "nFeature_RNA": "author_nFeature_RNA",
              "percent.mito": "author_percent_mito", "phase": "cell_cycle_phase",
              "celltype": "author_celltype", "lineage": "author_lineage"}
    md = md.rename(columns=rename)

    carry = [c for c in ["orig_ident", "experimental_group", "tamoxifen_start_day",
                         "sacrifice_day", "sex", "trace_call", "author_lineage",
                         "author_celltype", "cell_cycle_phase",
                         "author_nCount_RNA", "author_nFeature_RNA",
                         "author_percent_mito"] if c in md.columns]
    for c in carry:
        adata.obs[c] = key.map(md[c]).to_numpy()

    matched = key.isin(md.index)
    adata.obs["has_author_metadata"] = matched.to_numpy()
    n_match = int(matched.sum())
    record("metadata_table", str(cfg.metadata_csv))
    record("metadata_join_key", "(sample number parsed from orig.ident, bare "
                                "16nt barcode) - barcodes are not globally unique")
    record("cells_with_author_metadata", f"{n_match} / {adata.n_obs} "
                                         f"({100 * n_match / adata.n_obs:.1f}%)")

    annotated_samples = sorted(
        adata.obs.loc[matched.to_numpy(), "sample_id"].astype(str).unique())
    all_samples = sorted(adata.obs["sample_id"].astype(str).unique())
    unannotated = [s for s in all_samples if s not in annotated_samples]
    record("samples_without_any_metadata", unannotated)

    # Derived directly from the deposited metadata - nothing invented.
    if "experimental_group" in adata.obs:
        grp = adata.obs["experimental_group"].astype("object")
        cond = np.where(pd.isna(grp), "unannotated",
                        np.where(grp.astype(str).str.contains("homeostasis"),
                                 "homeostasis", "H1N1_infected"))
        adata.obs["condition"] = pd.Categorical(cond)
        record("condition_field",
               "derived from experimental_group: 'homeostasis' when the group "
               "name says homeostasis, otherwise 'H1N1_infected'; samples with "
               "no metadata row are labelled 'unannotated'")
    return adata


def _tidy_obs_for_h5ad(adata: ad.AnnData) -> None:
    # pandas 3 gives text columns the new `str` dtype rather than `object`, so
    # testing for `object` alone silently leaves string columns holding NaN.
    # scanpy then fails sorting categories ("'<' not supported between
    # instances of 'float' and 'str'").  Convert anything that is neither
    # numeric, boolean nor already categorical.
    for c in adata.obs.columns:
        s = adata.obs[c]
        if isinstance(s.dtype, pd.CategoricalDtype):
            if s.isna().any():
                adata.obs[c] = s.cat.add_categories(["NA"]).fillna("NA")
            continue
        if pd.api.types.is_numeric_dtype(s) or pd.api.types.is_bool_dtype(s):
            continue
        # missing values must read as "NA", not the string "nan"
        filled = s.astype("object").where(s.notna(), "NA").astype(str)
        adata.obs[c] = pd.Categorical(filled)
    for c in adata.var.columns:
        v = adata.var[c]
        if not (pd.api.types.is_numeric_dtype(v) or pd.api.types.is_bool_dtype(v)):
            adata.var[c] = v.astype("object").where(v.notna(), "NA").astype(str)


# =============================================================================
# PHASE 10 - normalisation and HVGs
# =============================================================================
def _load_checkpoint(path: Path) -> ad.AnnData:
    """Read a checkpoint and re-normalise its obs dtypes.

    A round trip through HDF5 can bring text columns back as a nullable string
    dtype holding NaN, which scanpy cannot sort into plot categories, so every
    resumed stage starts from the same tidy state a fresh run would have.
    """
    a = ad.read_h5ad(path)
    _tidy_obs_for_h5ad(a)
    return a


def normalize_data(adata: ad.AnnData, dirs: dict, n_hvg: int) -> ad.AnnData:
    """HVGs are chosen on the counts, then the counts are normalised in place.

    Selecting HVGs before normalising means the pipeline never has to hold a
    second full copy of the matrix, which matters at this cell number.
    """
    log(f"normalize_data: HVG selection on raw counts {mem_report()}")
    batch = "sample_id" if adata.obs["sample_id"].nunique() > 1 else None
    try:
        sc.pp.highly_variable_genes(adata, flavor="seurat_v3",
                                    n_top_genes=n_hvg, batch_key=batch)
        hvg_method = f"seurat_v3 on raw counts, n_top_genes={n_hvg}, batch_key={batch}"
    except Exception as exc:  # noqa: BLE001
        log(f"  seurat_v3 unavailable ({exc}); falling back to seurat on lognorm")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor="seurat", batch_key=batch)
        hvg_method = "seurat (fallback) on log-normalised data"
        record("normalization", "normalize_total(target_sum=1e4) + log1p")
        record("hvg_method", hvg_method)
        record("n_hvg", int(adata.var["highly_variable"].sum()))
        return adata

    record("hvg_method", hvg_method)
    record("n_hvg", int(adata.var["highly_variable"].sum()))

    log("normalize_data: writing the counts checkpoint before normalising")
    _tidy_obs_for_h5ad(adata)
    adata.write_h5ad(dirs["processed"] / "postQC.h5ad", compression="gzip")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.uns["expression_layer"] = "log1p(counts per 10,000)"
    record("normalization", "normalize_total(target_sum=1e4) + log1p; the raw "
                            "counts are preserved in processed/postQC.h5ad and "
                            "re-attached as layers['counts'] of the final object")
    log(f"normalize_data: done {mem_report()}")
    return adata


# =============================================================================
# PHASE 11-13 - PCA, batch inspection, integration, neighbours
# =============================================================================
def choose_n_pcs(adata: ad.AnnData, lo: int = 20, hi: int = 50) -> int:
    """Elbow of the variance-ratio curve, bracketed to a sensible range."""
    vr = adata.uns["pca"]["variance_ratio"]
    drops = -np.diff(vr)
    thresh = drops[0] * 0.02
    knee = int(np.argmax(drops < thresh)) + 1 if np.any(drops < thresh) else hi
    cum = np.cumsum(vr)
    n_cum = int(np.searchsorted(cum, 0.90)) + 1
    n = int(np.clip(max(knee, min(n_cum, hi)), lo, hi))
    record("n_pcs", f"{n} (variance-ratio elbow at {knee}, 90% cumulative "
                    f"variance at {n_cum}, bracketed to [{lo},{hi}])")
    return n


def run_pca(adata: ad.AnnData, dirs: dict, n_comps: int = 50) -> int:
    n_hvg = int(adata.var["highly_variable"].sum())
    big = adata.n_obs > 100_000
    if big:
        # Build the dense HVG matrix directly, chunk by chunk, instead of
        # going through scanpy's own paths.  Sparse ARPACK needs hundreds of
        # mat-vec passes (~54 min here) and scanpy's chunked route densifies
        # every gene before applying mask_var (~11 min at only 27k cells).
        # Filling a preallocated n_cells x n_hvg float32 array never holds more
        # than one chunk of surplus memory, and it makes ordinary scaled,
        # randomized PCA affordable at this cell number.
        from sklearn.decomposition import PCA as _PCA

        hv = np.flatnonzero(adata.var["highly_variable"].to_numpy())
        X = np.empty((adata.n_obs, hv.size), dtype=np.float32)
        step = 20_000
        for s in range(0, adata.n_obs, step):
            e = min(s + step, adata.n_obs)
            blk = adata.X[s:e][:, hv]
            X[s:e] = blk.toarray() if sp.issparse(blk) else np.asarray(blk)
        log(f"run_pca: dense HVG matrix {X.shape} "
            f"({X.nbytes / 1024 ** 3:.2f} GB) {mem_report()}")

        mu = X.mean(axis=0)
        X -= mu
        sd = X.std(axis=0)
        sd[sd == 0] = 1.0
        X /= sd
        np.clip(X, -10, 10, out=X)

        pca = _PCA(n_components=n_comps, svd_solver="randomized",
                   random_state=RANDOM_SEED)
        adata.obsm["X_pca"] = np.ascontiguousarray(pca.fit_transform(X),
                                                   dtype=np.float32)
        adata.uns["pca"] = {
            "variance_ratio": pca.explained_variance_ratio_,
            "variance": pca.explained_variance_,
        }
        loadings = np.zeros((adata.n_vars, n_comps), dtype=np.float32)
        loadings[hv] = pca.components_.T.astype(np.float32)
        adata.varm["PCs"] = loadings
        del X
        free_mem()
        record("scaling", "unit-variance scaling (z-score, clipped at +/-10) on "
                          "the HVG subset, then randomized PCA")
    else:
        sub = adata[:, adata.var["highly_variable"].to_numpy()].copy()
        sc.pp.scale(sub, max_value=10)
        sc.pp.pca(sub, n_comps=n_comps, svd_solver="arpack",
                  random_state=RANDOM_SEED)
        adata.obsm["X_pca"] = sub.obsm["X_pca"]
        adata.uns["pca"] = sub.uns["pca"]
        adata.varm["PCs"] = np.zeros((adata.n_vars, n_comps), dtype=np.float32)
        adata.varm["PCs"][adata.var["highly_variable"].to_numpy()] = sub.varm["PCs"]
        del sub
        free_mem()
        record("scaling", "sc.pp.scale(max_value=10) on the HVG subset")

    vr = adata.uns["pca"]["variance_ratio"]
    fig, ax = _plt.subplots(figsize=(5.2, 3.6))
    ax.plot(np.arange(1, len(vr) + 1), vr, "o-", ms=3, lw=1, color="#4C72B0")
    ax.set_xlabel("Principal component"); ax.set_ylabel("Variance ratio")
    ax.set_title("PCA variance explained")
    n_pcs = choose_n_pcs(adata)
    ax.axvline(n_pcs, ls="--", color="#C44E52", lw=1)
    ax.text(n_pcs + 0.6, max(vr) * 0.8, f"{n_pcs} PCs used", fontsize=8,
            color="#C44E52")
    fig.tight_layout()
    save_fig(fig, dirs["fig_qc"], "PCA_variance_ratio")
    log(f"run_pca: done, using {n_pcs} PCs {mem_report()}")
    return n_pcs


def inspect_batch(adata: ad.AnnData, n_pcs: int, dirs: dict) -> dict:
    """Decide from the data whether technical batch correction is warranted.

    Sample identity is confounded with experimental group in this design (each
    sample belongs to exactly one group), so 'do samples separate?' cannot by
    itself justify integration - separation may be the biology.  The question
    that *can* be answered is whether samples that share an experimental group,
    i.e. genuine replicates that differ only by animal, fail to mix.
    """
    from sklearn.neighbors import NearestNeighbors

    rng = np.random.default_rng(RANDOM_SEED)
    n_sub = min(30000, adata.n_obs)
    idx = rng.choice(adata.n_obs, n_sub, replace=False)
    X = adata.obsm["X_pca"][idx, :n_pcs]
    samples = adata.obs["sample_id"].astype(str).to_numpy()[idx]

    k = 30
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    _, ind = nn.kneighbors(X)
    ind = ind[:, 1:]
    same = (samples[ind] == samples[:, None]).mean(axis=1)

    freq = pd.Series(samples).value_counts(normalize=True)
    expected = freq.reindex(samples).to_numpy()
    # 1.0 = neighbours are exactly as mixed as random; large values = separation
    enrich = float(np.mean(same / expected))

    out: dict = {"knn_same_sample_enrichment": round(enrich, 3),
                 "mean_same_sample_neighbour_fraction": round(float(same.mean()), 4),
                 "expected_by_chance": round(float(expected.mean()), 4)}

    # replicate-level mixing: samples inside one experimental group
    rep_enrich = np.nan
    if "experimental_group" in adata.obs:
        groups = adata.obs["experimental_group"].astype(str).to_numpy()[idx]
        vals = []
        for g in pd.unique(groups):
            if g in ("nan", "NA", "unannotated"):
                continue
            m = groups == g
            if m.sum() < 200 or len(np.unique(samples[m])) < 2:
                continue
            sub_same = (samples[ind][m] == samples[m][None].T)
            in_group = groups[ind][m] == g
            denom = in_group.sum(axis=1)
            ok = denom > 0
            if ok.sum() < 50:
                continue
            frac = (sub_same & in_group)[ok].sum(axis=1) / denom[ok]
            fr = pd.Series(samples[m]).value_counts(normalize=True)
            exp = float((fr ** 2).sum())
            vals.append(float(frac.mean()) / exp)
        if vals:
            rep_enrich = float(np.mean(vals))
    out["within_group_replicate_enrichment"] = (round(rep_enrich, 3)
                                                if rep_enrich == rep_enrich
                                                else None)
    has_groups = rep_enrich == rep_enrich
    out["replicate_structure_available"] = bool(has_groups)

    # Decision rule, stated up front so it cannot be tuned after seeing a UMAP.
    if has_groups:
        needs = rep_enrich > 2.0
        out["decision_rule"] = (
            "Integrate only if biological replicates within the SAME "
            "experimental group fail to mix (within-group replicate kNN "
            "enrichment > 2.0). Sample identity alone is not evidence of a "
            "technical batch effect here, because every sample belongs to "
            "exactly one experimental group - correcting on sample would also "
            "remove the experimental effect.")
    else:
        needs = False
        out["decision_rule"] = (
            "No condition or group metadata exists for this series, so donor, "
            "batch and biology are completely confounded and no correction can "
            "be justified as purely technical. The unintegrated embedding is "
            "kept as primary; a Harmony embedding is computed alongside it, "
            "clearly labelled as supplementary, if samples separate strongly.")
    out["integration_recommended"] = bool(needs)
    out["supplementary_harmony_warranted"] = bool(enrich > 2.0)
    (dirs["qc"] / "batch_assessment.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    record("batch_assessment", out)
    return out


def _run_harmony(adata: ad.AnnData, n_pcs: int) -> bool:
    try:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"][:, :n_pcs],
                                   adata.obs, ["sample_id"],
                                   random_state=RANDOM_SEED)
        adata.obsm["X_pca_harmony"] = np.ascontiguousarray(
            np.asarray(ho.Z_corr).T, dtype=np.float32)
        return True
    except Exception as exc:  # noqa: BLE001
        log(f"  harmony failed: {type(exc).__name__}: {exc}")
        return False


def run_integration_if_needed(adata: ad.AnnData, assessment: dict, n_pcs: int,
                              force: str | None) -> tuple[str, str | None]:
    """Return (primary representation, supplementary representation or None)."""
    if force == "none":
        record("batch_correction", "explicitly disabled by --integration none")
        return "X_pca", None

    want_primary = assessment.get("integration_recommended") or force == "harmony"
    want_supp = assessment.get("supplementary_harmony_warranted")

    if want_primary:
        if _run_harmony(adata, n_pcs):
            if force == "harmony":
                why = ("requested explicitly (--integration harmony). The "
                       "justification recorded for this dataset is that single "
                       "proposed cell types were fragmenting into "
                       "sample-private clusters in the uncorrected embedding "
                       "(see qc/celltype_split_by_sample.csv): one cell type is "
                       "not several cell types in several donors. This series "
                       "carries no condition metadata, so there is no designed "
                       "experimental contrast that correcting on sample could "
                       "destroy")
            else:
                why = ("biological replicates within the same experimental "
                       "group failed to mix (within-group enrichment "
                       f"{assessment.get('within_group_replicate_enrichment')} "
                       f"> 2.0)")
            record("batch_correction",
                   f"Harmony (harmonypy) on sample_id used as the PRIMARY "
                   f"embedding, {why}. The unintegrated X_pca and its UMAP are "
                   f"retained alongside it for comparison.")
            return "X_pca_harmony", "X_pca"
        record("batch_correction", "Harmony was warranted but failed; "
                                   "continuing unintegrated")
        return "X_pca", None

    # A supplementary embedding is a nice-to-have, not a result.  Harmony cost
    # grows with cells x batches (roughly 11 min per iteration at 160k cells
    # across 33 batches), so it is not worth hours of compute when the
    # assessment has already concluded that correction is not warranted.
    SUPP_MAX_CELLS = 60_000
    if want_supp and adata.n_obs > SUPP_MAX_CELLS:
        record("batch_correction",
               f"not applied, and the optional supplementary Harmony embedding "
               f"was skipped: correction is not warranted here (within-group "
               f"replicate enrichment "
               f"{assessment.get('within_group_replicate_enrichment')} <= 2.0, "
               f"i.e. replicates inside an experimental group already mix), and "
               f"at {adata.n_obs:,} cells across "
               f"{adata.obs['sample_id'].nunique()} batches Harmony would cost "
               f"hours for an embedding that would not be used. Pass "
               f"--integration harmony to force it. The mixing statistics are "
               f"in qc/batch_assessment.json.")
        return "X_pca", None

    if want_supp and _run_harmony(adata, n_pcs):
        if assessment.get("replicate_structure_available"):
            why = ("Replicates within an experimental group mix acceptably "
                   f"(within-group enrichment "
                   f"{assessment.get('within_group_replicate_enrichment')} "
                   f"<= 2.0), and sample identity is fully nested inside "
                   f"experimental group, so correcting on sample would remove "
                   f"the experimental effect as well.")
        else:
            why = ("This series has no condition metadata, so donor, batch and "
                   "biology cannot be separated and no correction can be "
                   "justified as purely technical.")
        record("batch_correction",
               f"NOT applied to the primary embedding. {why} A Harmony "
               f"embedding on sample_id is nevertheless provided as a clearly "
               f"labelled SUPPLEMENTARY representation (X_pca_harmony, "
               f"X_umap_harmony) because samples do separate in PCA space "
               f"(kNN same-sample enrichment "
               f"{assessment.get('knn_same_sample_enrichment')}). It must not "
               f"be used for condition comparisons.")
        return "X_pca", "X_pca_harmony"

    record("batch_correction",
           "not applied; samples already mix in the unintegrated embedding "
           f"(kNN same-sample enrichment "
           f"{assessment.get('knn_same_sample_enrichment')})")
    return "X_pca", None


def build_embeddings(adata: ad.AnnData, n_pcs: int, primary: str,
                     secondary: str | None, dirs: dict,
                     n_neighbors: int = 15) -> None:
    """Build the primary neighbour graph and UMAP, plus an optional second one.

    The secondary embedding is computed first and stashed, so that the graph
    and X_umap left on the object at the end are unambiguously the primary.
    """
    if secondary:
        log(f"build_embeddings: secondary embedding from {secondary}")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                        use_rep=secondary, random_state=RANDOM_SEED,
                        key_added="neighbors_alt")
        sc.tl.umap(adata, neighbors_key="neighbors_alt",
                   random_state=RANDOM_SEED)
        tag = "harmony" if "harmony" in secondary else "unintegrated"
        adata.obsm[f"X_umap_{tag}"] = adata.obsm["X_umap"].copy()
        adata.uns["secondary_embedding"] = {"rep": secondary, "umap_key":
                                            f"X_umap_{tag}"}

    log(f"build_embeddings: primary embedding from {primary}")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
                    use_rep=primary, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)
    record("neighbors", f"n_neighbors={n_neighbors}, n_pcs={n_pcs}, "
                        f"use_rep={primary}, metric=euclidean")
    record("umap", f"scanpy defaults min_dist=0.5, spread=1.0, "
                   f"random_state={RANDOM_SEED}")

    if secondary:
        _integration_comparison_figure(adata, dirs, primary, secondary)


def _integration_comparison_figure(adata: ad.AnnData, dirs: dict,
                                   primary: str, secondary: str) -> None:
    alt_key = adata.uns["secondary_embedding"]["umap_key"]
    lbl = {"X_pca": "unintegrated PCA", "X_pca_harmony": "Harmony-corrected"}
    samples = adata.obs["sample_id"].astype(str)
    cats = sorted(samples.unique())
    cmap = _plt.get_cmap("tab20")
    colours = {c: cmap(i % 20) for i, c in enumerate(cats)}
    cvec = samples.map(colours).to_numpy()

    fig, axes = _plt.subplots(1, 2, figsize=(13, 5.6))
    for ax, key, name in [(axes[0], "X_umap", primary),
                          (axes[1], alt_key, secondary)]:
        um = adata.obsm[key]
        ax.scatter(um[:, 0], um[:, 1], s=0.6, c=list(cvec), rasterized=True,
                   linewidths=0)
        ax.set_title(f"{lbl.get(name, name)}"
                     f"{'  (PRIMARY)' if key == 'X_umap' else '  (supplementary)'}",
                     fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    handles = [_plt.Line2D([], [], marker="o", ls="", ms=4, color=colours[c],
                           label=c) for c in cats]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5),
               fontsize=6, ncol=1 if len(cats) < 20 else 2, title="sample")
    fig.suptitle("Effect of batch correction on sample mixing", y=1.02)
    fig.tight_layout()
    save_fig(fig, dirs["fig_umap"], "UMAP_integration_before_after")


# =============================================================================
# PHASE 14 - clustering
# =============================================================================
def _marker_support(adata: ad.AnnData, key: str, species: str,
                    min_markers: int = 5) -> tuple[int, int]:
    """How many clusters at this resolution are actually supported by markers?

    Uses the fast t-test rather than Wilcoxon because this runs once per
    candidate resolution; the final marker table is still Wilcoxon.  Genes that
    are only mitochondrial / ribosomal / heat-shock / haemoglobin do not count
    as support, so a cluster of stressed cells cannot look well-supported.
    """
    # Restricted to the HVGs.  This is ~13x less work than the full gene space
    # at no cost to the criterion: a cluster whose only distinguishing genes
    # are non-variable ones is not a cluster the embedding could have found,
    # and restricting the candidate pool makes the test stricter, not looser.
    # It also makes this consistent with _pairwise_separation, which already
    # works in HVG space.
    tmp = f"_scan_{key}"
    mask = ("highly_variable" if "highly_variable" in adata.var else None)
    sc.tl.rank_genes_groups(adata, groupby=key, method="t-test", pts=True,
                            key_added=tmp, mask_var=mask)
    df = sc.get.rank_genes_groups_df(adata, group=None, key=tmp)
    df = df[(df["pvals_adj"] < 0.05) & (df["logfoldchanges"] > 1.0)
            & (df["pct_nz_group"] > 0.25)]
    df = df[~df["names"].map(lambda g: MK.is_uninformative(g, species))]
    counts = df.groupby("group", observed=True).size()
    all_groups = adata.obs[key].astype(str).unique()
    counts = counts.reindex(all_groups, fill_value=0)
    del adata.uns[tmp]
    return int((counts >= min_markers).sum()), int(len(all_groups))


def _group_stats(X: sp.csr_matrix, mask: np.ndarray):
    Xs = X[mask]
    n = Xs.shape[0]
    mean = np.asarray(Xs.mean(axis=0)).ravel()
    sq = np.asarray(Xs.multiply(Xs).mean(axis=0)).ravel()
    var = np.maximum(sq - mean ** 2, 0.0) * (n / max(n - 1, 1))
    frac = np.asarray((Xs > 0).mean(axis=0)).ravel()
    return n, mean, var, frac


def _pairwise_separation(adata: ad.AnnData, key: str, species: str,
                         max_per_cluster: int = 400, min_de: int = 10
                         ) -> tuple[int, int]:
    """Are neighbouring clusters actually distinguishable from each other?

    A cluster can look well-supported against "all other cells" while being
    essentially identical to its nearest neighbour - that is what
    over-clustering looks like.  For every cluster this compares it with its
    nearest neighbouring cluster (PCA centroid distance) and counts genes that
    separate the pair.  Returns (smallest number of separating genes over all
    pairs, number of pairs below the threshold).
    """
    from scipy import stats as sstats
    from statsmodels.stats.multitest import multipletests

    rng = np.random.default_rng(RANDOM_SEED)
    labels = adata.obs[key].astype(str).to_numpy()
    cats = pd.unique(labels)
    if len(cats) < 2:
        return 10 ** 6, 0

    hv = (adata.var["highly_variable"].to_numpy()
          if "highly_variable" in adata.var else np.ones(adata.n_vars, bool))
    informative = np.array([not MK.is_uninformative(g, species)
                            for g in adata.var_names[hv]])

    P = adata.obsm["X_pca"]
    cent = np.vstack([P[labels == c].mean(axis=0) for c in cats])
    d = np.linalg.norm(cent[:, None, :] - cent[None, :, :], axis=2)
    np.fill_diagonal(d, np.inf)

    # Subsample cells BEFORE slicing genes: taking the HVG columns out of the
    # full matrix first would allocate a second near-full-size copy.
    idx_by_cat = {}
    for c in cats:
        ii = np.flatnonzero(labels == c)
        if ii.size > max_per_cluster:
            ii = rng.choice(ii, max_per_cluster, replace=False)
        idx_by_cat[c] = ii
    keep_idx = np.concatenate([idx_by_cat[c] for c in cats])
    Xf = adata.X.tocsr() if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    X = Xf[keep_idx][:, hv].tocsr()
    sub_labels = labels[keep_idx]
    del Xf

    seen: set[tuple[str, str]] = set()
    worst, n_bad = 10 ** 6, 0
    for i, c in enumerate(cats):
        j = int(np.argmin(d[i]))
        pair = tuple(sorted((c, cats[j])))
        if pair in seen:
            continue
        seen.add(pair)
        ma = sub_labels == pair[0]
        mb = sub_labels == pair[1]
        na, mua, va, fa = _group_stats(X, ma)
        nb, mub, vb, fb = _group_stats(X, mb)
        se = np.sqrt(va / max(na, 1) + vb / max(nb, 1))
        ok = se > 0
        t = np.zeros_like(mua); t[ok] = (mua[ok] - mub[ok]) / se[ok]
        df_ = np.full_like(mua, max(na + nb - 2, 1), dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            num = (va / na + vb / nb) ** 2
            den = (va ** 2 / (na ** 2 * max(na - 1, 1))
                   + vb ** 2 / (nb ** 2 * max(nb - 1, 1)))
            df_[den > 0] = num[den > 0] / den[den > 0]
        p = np.ones_like(mua)
        p[ok] = 2 * sstats.t.sf(np.abs(t[ok]), df_[ok])
        padj = np.ones_like(p)
        padj[ok] = multipletests(p[ok], method="fdr_bh")[1]
        lfc = (mua - mub) / np.log(2.0)
        strong = ((padj < 0.05) & (np.abs(lfc) > 1.0)
                  & ((fa > 0.25) | (fb > 0.25)) & informative)
        n_strong = int(strong.sum())
        worst = min(worst, n_strong)
        if n_strong < min_de:
            n_bad += 1
    return worst, n_bad


def run_clustering(adata: ad.AnnData, dirs: dict, species: str,
                   resolutions=(0.3, 0.5, 0.8, 1.0)) -> str:
    min_size = max(20, int(0.0005 * adata.n_obs))
    min_de = 10
    stats = []
    for r in resolutions:
        key = f"leiden_res{r}"
        sc.tl.leiden(adata, resolution=r, key_added=key, flavor="igraph",
                     n_iterations=2, directed=False, random_state=RANDOM_SEED)
        vc = adata.obs[key].value_counts()
        supported, total = _marker_support(adata, key, species)
        worst, n_bad = _pairwise_separation(adata, key, species, min_de=min_de)
        stats.append({
            "resolution": r, "n_clusters": int(vc.size),
            "smallest_cluster": int(vc.min()), "largest_cluster": int(vc.max()),
            "min_cluster_size_threshold": min_size,
            "clusters_below_min_size": int((vc < min_size).sum()),
            "marker_supported_clusters": supported,
            "pct_clusters_marker_supported": round(100 * supported / total, 1),
            "min_DE_genes_between_nearest_clusters": int(worst),
            "nearest_cluster_pairs_under_threshold": int(n_bad),
            "pct_cells_in_largest": round(100 * vc.max() / adata.n_obs, 2)})
        log(f"  leiden res={r}: {vc.size} clusters, smallest {vc.min()}, "
            f"{supported}/{total} marker-supported, "
            f"{n_bad} under-separated neighbour pair(s) (worst={worst} DE genes)")
    df = pd.DataFrame(stats)
    df.to_csv(dirs["tables"] / "clustering_resolution_scan.csv", index=False)

    # A resolution is acceptable when clusters are not fragments, are each
    # distinguishable from the rest of the data, AND are distinguishable from
    # their own nearest neighbour.  The last condition is what stops the scan
    # from simply preferring whichever resolution makes the most clusters.
    ok = df[(df["clusters_below_min_size"] == 0)
            & (df["pct_clusters_marker_supported"] >= 90.0)
            & (df["nearest_cluster_pairs_under_threshold"] == 0)]
    if len(ok):
        chosen = float(ok["resolution"].max())
        why = (f"It was the finest resolution at which no cluster falls below {min_size} "
               f"cells, at least 90% of clusters carry >=5 specific marker "
               f"genes against the rest of the data, and every cluster is also "
               f"separated from its NEAREST neighbouring cluster by at least "
               f"{min_de} genes (adj. p<0.05, |log2FC|>1, >25% expressing) - "
               f"the last condition is what prevents simply picking whichever "
               f"resolution yields the most clusters")
    else:
        best = df.sort_values(
            ["nearest_cluster_pairs_under_threshold", "clusters_below_min_size",
             "pct_clusters_marker_supported"],
            ascending=[True, True, False]).iloc[0]
        chosen = float(best["resolution"])
        why = (f"No resolution satisfied all three criteria, so the resolution "
               f"with the fewest under-separated neighbouring cluster pairs "
               f"({int(best['nearest_cluster_pairs_under_threshold'])}) and the "
               f"fewest undersized clusters "
               f"({int(best['clusters_below_min_size'])}) was taken. Every finer "
               f"resolution produced neighbouring clusters that no gene "
               f"separates, which is the signature of over-clustering.")
    adata.obs["leiden_cluster"] = adata.obs[f"leiden_res{chosen}"]
    n_cl = adata.obs["leiden_cluster"].nunique()
    record("clustering", f"Leiden (igraph flavour, 2 iterations) scanned at "
                         f"resolutions {list(resolutions)}. Primary resolution: "
                         f"{chosen}. {why} All resolutions are retained as "
                         f"leiden_res* columns.")
    record("leiden_resolution", chosen)
    record("n_clusters", int(n_cl))
    return f"leiden_res{chosen}"


# =============================================================================
# PHASE 16 - markers
# =============================================================================
def find_markers(adata: ad.AnnData, dirs: dict, groupby: str = "leiden_cluster",
                 prefix: str = "cluster", reuse: bool = False) -> pd.DataFrame:
    # The Wilcoxon test costs ~45 min at this cell count.  When only a
    # downstream step changed (annotation, plotting), re-using the table this
    # pipeline itself wrote is exact, not an approximation - but it is opt-in,
    # and it refuses to serve a table whose clusters no longer match.
    cached = dirs["tables"] / f"{prefix}_markers_all.csv"
    if reuse and cached.exists():
        df = pd.read_csv(cached)
        want = set(adata.obs[groupby].astype(str).unique())
        have = set(df["cluster"].astype(str).unique())
        if want == have:
            log(f"find_markers: reusing {cached.name} "
                f"({len(df):,} rows, {len(have)} clusters)")
            record("marker_test", "Wilcoxon rank-sum (scanpy rank_genes_groups, "
                                  "one cluster vs all remaining cells) on "
                                  "log1p(CP10K) values, with expressing "
                                  "fractions; table re-used from the previous "
                                  "run of this same pipeline")
            return df
        log(f"find_markers: cached table has clusters {sorted(have)[:5]}... but "
            f"the object has {sorted(want)[:5]}... - recomputing")

    log(f"find_markers: Wilcoxon rank-sum over {adata.obs[groupby].nunique()} "
        f"groups {mem_report()}")
    sc.tl.rank_genes_groups(adata, groupby=groupby, method="wilcoxon", pts=True,
                            key_added=f"rank_{groupby}")
    df = sc.get.rank_genes_groups_df(adata, group=None, key=f"rank_{groupby}")
    df = df.rename(columns={"group": "cluster", "names": "gene",
                            "logfoldchanges": "log2FC", "pvals": "pval",
                            "pvals_adj": "pval_adj",
                            "pct_nz_group": "frac_expressing_in_cluster",
                            "pct_nz_reference": "frac_expressing_other"})
    df.to_csv(dirs["tables"] / f"{prefix}_markers_all.csv", index=False)

    top = (df.sort_values(["cluster", "scores"], ascending=[True, False])
             .groupby("cluster", observed=True).head(20))
    top.to_csv(dirs["tables"] / f"{prefix}_markers_top20.csv", index=False)
    record("marker_test", "Wilcoxon rank-sum (scanpy rank_genes_groups, "
                          "one cluster vs all remaining cells) on "
                          "log1p(CP10K) values, with expressing fractions")
    return df


# =============================================================================
# PHASE 23 - assisted annotation
# =============================================================================
def _cluster_mean_expression(adata: ad.AnnData, genes: list[str],
                             groupby: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    genes = [g for g in genes if g in adata.var_names]
    if not genes:
        return pd.DataFrame(), pd.DataFrame()
    sub = adata[:, genes]
    X = sub.X
    X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
    groups = adata.obs[groupby].astype(str)
    cats = sorted(groups.unique(), key=lambda s: (len(s), s))
    mean = pd.DataFrame(index=cats, columns=genes, dtype=float)
    frac = pd.DataFrame(index=cats, columns=genes, dtype=float)
    codes = groups.to_numpy()
    for c in cats:
        m = codes == c
        Xc = X[m]
        mean.loc[c] = np.asarray(Xc.mean(axis=0)).ravel()
        frac.loc[c] = np.asarray((Xc > 0).mean(axis=0)).ravel()
    return mean, frac


def annotate_clusters(adata: ad.AnnData, markers_df: pd.DataFrame, dirs: dict,
                      species: str, groupby: str = "leiden_cluster",
                      outfile: str = "cluster_annotation_proposals.csv"
                      ) -> pd.DataFrame:
    panel = MK.get_panel(species)
    genes = sorted({g for gs in panel.values() for g in gs})
    mean, frac = _cluster_mean_expression(adata, genes, groupby)
    if mean.empty:
        return pd.DataFrame()

    z = (mean - mean.mean(axis=0)) / (mean.std(axis=0).replace(0, np.nan))
    z = z.fillna(0.0)

    # The cluster's own strongest DE genes, used to corroborate panel scores.
    # Relative enrichment alone is not enough: a small panel can win on
    # z-score while the cluster's actual top genes point somewhere else
    # entirely (ambient haemoglobin scoring as "erythroid" is the classic
    # case). A panel only counts as strong evidence if some of its genes are
    # genuinely among the genes that define the cluster.
    top_de: dict[str, set[str]] = {}
    for cl_, sub_ in markers_df.groupby(markers_df["cluster"].astype(str),
                                        observed=True):
        s = sub_[(sub_["pval_adj"] < 0.05)
                 & (sub_["frac_expressing_in_cluster"] > 0.20)]
        s = s.sort_values("scores", ascending=False)
        top_de[cl_] = set(g for g in s["gene"].head(100)
                          if not MK.is_uninformative(g, species))

    rows = []
    for cl in mean.index:
        de_set = top_de.get(str(cl), set())
        scores, de_frac = {}, {}
        for grp, gs in panel.items():
            gs = [g for g in gs if g in mean.columns]
            if not gs:
                continue
            # a group counts only if the genes are both enriched and detected
            z_term = float(np.mean(
                [z.loc[cl, g] * min(1.0, frac.loc[cl, g] / 0.10) for g in gs]))
            hits = len(set(gs) & de_set)
            d_term = hits / max(1, min(len(gs), 4))
            de_frac[grp] = d_term
            scores[grp] = z_term + 2.5 * d_term
        ranked = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        best, best_s = ranked[0]
        second, second_s = ranked[1] if len(ranked) > 1 else ("", np.nan)

        support = [g for g in panel[best] if g in mean.columns
                   and z.loc[cl, g] > 0.5 and frac.loc[cl, g] > 0.20]
        conflict = []
        for grp, s in ranked[1:]:
            if s > 0.8 * best_s and s > 0.75:
                conflict += [g for g in panel[grp] if g in mean.columns
                             and z.loc[cl, g] > 0.8 and frac.loc[cl, g] > 0.25]
        conflict = [g for g in dict.fromkeys(conflict) if g not in support]

        gap = best_s - (second_s if second_s == second_s else 0.0)
        corroborated = de_frac.get(best, 0.0)
        # High confidence requires the winning panel to be corroborated by the
        # cluster's own DE genes, not just to win on relative enrichment.
        if best_s > 1.2 and gap > 0.5 and len(support) >= 2 and corroborated >= 0.5:
            conf = "High"
        elif best_s > 0.8 and len(support) >= 2 and corroborated > 0:
            conf = "Moderate"
        elif best_s > 0.4:
            conf = "Low"
        else:
            conf = "Uncertain"

        cl_de = (markers_df[markers_df["cluster"].astype(str) == str(cl)]
                 .sort_values("scores", ascending=False))
        top_names = [g for g in cl_de["gene"].head(40).tolist()
                     if not MK.is_uninformative(g, species)][:10]
        only_junk = len(top_names) == 0
        notes = []
        if only_junk:
            notes.append("top DE genes are exclusively mitochondrial / "
                         "ribosomal / heat-shock / haemoglobin - treat as a "
                         "low-quality or stressed population, not a cell type")
            conf = "Uncertain"
        if conflict:
            notes.append("co-expresses markers of more than one lineage - "
                         "check for doublets")
        rows.append({
            "Cluster": cl,
            "n_cells": int((adata.obs[groupby].astype(str) == str(cl)).sum()),
            "Proposed identity": f"{best} (candidate)",
            "Runner-up": second,
            "Panel score": round(best_s, 3),
            "Runner-up score": round(second_s, 3) if second_s == second_s else "",
            "DE corroboration": round(corroborated, 3),
            "Supporting markers": ", ".join(support[:8]),
            "Conflicting markers": ", ".join(conflict[:8]),
            "Top data-derived genes": ", ".join(top_names),
            "Confidence": conf,
            "Notes": "; ".join(notes),
        })

    df = pd.DataFrame(rows)

    # Report the deposited label alongside each proposal where one exists.
    # This is REPORTED ONLY and never enters the scoring, so the agreement
    # statistic stays an independent check - but a reader can immediately see
    # when a marker-panel proposal contradicts the deposited annotation, which
    # is exactly the case a candidate label must not be trusted in.
    if "author_celltype" in adata.obs and "has_author_metadata" in adata.obs:
        m = adata.obs["has_author_metadata"].astype(str).isin(["True", "true"])
        if m.sum():
            sub = adata.obs.loc[m.to_numpy()]
            ct = pd.crosstab(sub[groupby].astype(str),
                             sub["author_celltype"].astype(str))
            frac = ct.div(ct.sum(axis=1).replace(0, np.nan), axis=0)
            dom = frac.idxmax(axis=1)
            pur = frac.max(axis=1)
            key = df["Cluster"].astype(str)
            df["Deposited label (dominant)"] = key.map(dom).fillna("n/a")
            df["Deposited label %"] = (100 * key.map(pur)).round(1)
            df["Annotated cells"] = key.map(ct.sum(axis=1)).fillna(0).astype(int)

            def _flag(r):
                d = str(r["Deposited label (dominant)"])
                if d in ("n/a", "nan"):
                    return ""
                prop = str(r["Proposed identity"]).replace(" (candidate)", "")
                p = prop.lower().replace("_", "")
                dl = d.lower().replace("_", "")
                # Cell-cycle and stress calls describe a STATE, not a lineage,
                # so they are not comparable to a lineage label and must not be
                # reported as a contradiction.
                if p in {"proliferating"}:
                    return "not comparable (state, not lineage)"
                if p in dl or dl in p:
                    return "agrees"
                # families naming the same population differently
                fam = [{"endothelial", "endothelialcap", "cap1", "cap2",
                        "venousendothelium", "arterialendothelium"},
                       {"fibroblast", "fibroactivated", "af1", "af2",
                        "adventitialfibroblast", "peribronchialfibroblast"},
                       {"monomac", "alveolarmac", "amac", "imac", "cmon",
                        "imon", "pmon", "dc", "cdc1", "cdc2", "madc",
                        "neutrophil"},
                       {"tcell", "cd4t", "cd8t", "treg", "tlymphocyte",
                        "nk", "nkcell"},
                       {"bcell", "blymphocyte", "plasma", "plasmacell"},
                       {"smcpericyte", "vsmc"},
                       {"lymphaticec", "lymphaticendothelium"},
                       # Krt5 IS the deposited name for basal cells
                       {"basal", "krt5"},
                       {"club", "secretory"},
                       {"transitional", "alveolartransitional"},
                       {"at1", "at1at2"}, {"at2", "at1at2"}]
                for f in fam:
                    if p in f and dl in f:
                        return "agrees (same lineage)"
                if int(r.get("Annotated cells", 0) or 0) < 50:
                    return "too few annotated cells to judge"
                return "DISAGREES - treat the proposal as unreliable"

            df["Deposition check"] = df.apply(_flag, axis=1)
            n_dis = int(df["Deposition check"].astype(str)
                        .str.startswith("DISAGREES").sum())
            record("annotation_vs_deposited_labels",
                   f"{n_dis} of {len(df)} clusters have a marker-panel proposal "
                   f"that contradicts the dominant deposited label; those rows "
                   f"are flagged in {outfile} and their proposals should not be "
                   f"used")

    df.to_csv(dirs["tables"] / outfile, index=False)
    mean.to_csv(dirs["tables"] / "canonical_marker_mean_expression_by_cluster.csv")
    frac.to_csv(dirs["tables"] / "canonical_marker_fraction_expressing_by_cluster.csv")

    # A proposal the deposited annotation contradicts must not travel as a
    # plain label - it would be read as a finding. Mark it in the value itself,
    # so the flag survives into every figure and export, not only the table.
    labels = df["Proposed identity"].astype(str).copy()
    if "Deposition check" in df.columns:
        bad = df["Deposition check"].astype(str).str.startswith("DISAGREES")
        labels = labels.where(~bad, labels + " [CONTRADICTED]")
    adata.obs["proposed_cell_type"] = (
        adata.obs[groupby].astype(str)
        .map(dict(zip(df["Cluster"].astype(str), labels))))
    record("annotation_strategy",
           "candidate identities scored per cluster from z-scored mean "
           "expression of curated panels, down-weighted by the fraction of "
           "cells expressing each marker; numeric Leiden labels are preserved "
           "and no candidate is promoted to a definitive label")
    return df


def validate_against_author_labels(adata: ad.AnnData, dirs: dict,
                                   groupby: str = "leiden_cluster") -> None:
    if "author_celltype" not in adata.obs:
        return
    m = adata.obs["has_author_metadata"].astype(bool).to_numpy()
    if m.sum() == 0:
        return
    sub = adata.obs.loc[m]
    ct = pd.crosstab(sub[groupby].astype(str), sub["author_celltype"].astype(str))
    ct.to_csv(dirs["tables"] / "cluster_vs_author_celltype_counts.csv")
    frac = ct.div(ct.sum(axis=1), axis=0)
    frac.round(4).to_csv(dirs["tables"] / "cluster_vs_author_celltype_fraction.csv")

    best = frac.idxmax(axis=1)
    purity = frac.max(axis=1)
    pd.DataFrame({"cluster": best.index, "dominant_author_celltype": best.values,
                  "fraction_of_cluster": purity.round(4).values,
                  "n_annotated_cells": ct.sum(axis=1).values}).to_csv(
        dirs["tables"] / "cluster_dominant_author_celltype.csv", index=False)
    record("author_label_agreement",
           f"median cluster purity against the deposited author cell-type "
           f"labels = {purity.median():.3f} over {int(m.sum())} annotated cells "
           f"(independent check only; author labels were not used to build the "
           f"clustering)")

    fig, ax = _plt.subplots(figsize=(max(7, 0.42 * frac.shape[1] + 3),
                                     max(4, 0.32 * frac.shape[0] + 2)))
    im = ax.imshow(frac.to_numpy(), aspect="auto", cmap="magma", vmin=0, vmax=1)
    ax.set_xticks(range(frac.shape[1]))
    ax.set_xticklabels(frac.columns, rotation=90, fontsize=6)
    ax.set_yticks(range(frac.shape[0]))
    ax.set_yticklabels(frac.index, fontsize=6)
    ax.set_xlabel("Deposited author cell type"); ax.set_ylabel("Leiden cluster")
    ax.set_title("Cluster composition vs deposited annotation")
    fig.colorbar(im, ax=ax, label="fraction of cluster")
    fig.tight_layout()
    save_fig(fig, dirs["fig_comp"], "cluster_vs_author_celltype_heatmap")


# =============================================================================
# PHASE 15, 21, 22 - figures
# =============================================================================
def generate_umaps(adata: ad.AnnData, dirs: dict, prefix: str = "") -> None:
    failed: list[str] = []

    def _umap(color, name, **kw):
        if color not in adata.obs and color not in adata.var_names:
            return
        # One awkward column must not cost the whole figure stage, which sits
        # an hour downstream of the last checkpoint.
        try:
            ax = sc.pl.umap(adata, color=color, show=False, frameon=False, **kw)
            save_fig(figure_of(ax), dirs["fig_umap"], f"{prefix}{name}")
        except Exception as exc:  # noqa: BLE001
            failed.append(f"{name} (color={color}): {type(exc).__name__}: {exc}")
            log(f"  WARNING could not plot {name}: {exc}")
            _plt.close("all")

    _umap("leiden_cluster", "UMAP_leiden_clusters", legend_loc="on data",
          legend_fontsize=7, legend_fontoutline=2.0, size=3,
          title="Leiden clusters", palette=sc.pl.palettes.default_102)
    _umap("leiden_cluster", "UMAP_leiden_clusters_sidelegend", size=3,
          title="Leiden clusters", palette=sc.pl.palettes.default_102)
    _umap("sample_id", "UMAP_sample", size=2, title="Sample of origin")
    _umap("condition", "UMAP_condition", size=2, title="Condition")
    _umap("experimental_group", "UMAP_experimental_group", size=2,
          title="Experimental group")
    _umap("author_lineage", "UMAP_author_lineage", size=2,
          title="Deposited lineage annotation")
    _umap("author_celltype", "UMAP_author_celltype", size=2,
          title="Deposited cell-type annotation")
    _umap("trace_call", "UMAP_trace_call", size=2, title="Lineage-trace call")
    _umap("sex", "UMAP_sex", size=2, title="Sex")
    _umap("proposed_cell_type", "UMAP_proposed_cell_type", size=2,
          title="Proposed identity (candidate)")

    for qc, name in [("pct_counts_mt", "UMAP_qc_pct_mito"),
                     ("total_counts", "UMAP_qc_total_counts"),
                     ("n_genes_by_counts", "UMAP_qc_n_genes"),
                     ("doublet_score", "UMAP_qc_doublet_score")]:
        if qc in adata.obs:
            ax = sc.pl.umap(adata, color=qc, show=False, frameon=False, size=2,
                            cmap="viridis")
            save_fig(figure_of(ax), dirs["fig_qc"], f"{prefix}{name}")

    # split views keep identical coordinates so panels stay comparable
    for field, fname in [("condition", "UMAP_split_by_condition"),
                         ("experimental_group", "UMAP_split_by_group")]:
        if field not in adata.obs:
            continue
        cats = [c for c in adata.obs[field].astype(str).unique()]
        cats = sorted(cats)
        if len(cats) < 2 or len(cats) > 20:
            continue
        ncol = min(4, len(cats))
        nrow = int(np.ceil(len(cats) / ncol))
        fig, axes = _plt.subplots(nrow, ncol, figsize=(3.6 * ncol, 3.4 * nrow),
                                  squeeze=False)
        um = adata.obsm["X_umap"]
        vals = adata.obs[field].astype(str).to_numpy()
        for i, c in enumerate(cats):
            ax = axes[i // ncol][i % ncol]
            ax.scatter(um[:, 0], um[:, 1], s=0.8, c="#E5E5E5", rasterized=True,
                       linewidths=0)
            m = vals == c
            ax.scatter(um[m, 0], um[m, 1], s=0.8, c="#C44E52", rasterized=True,
                       linewidths=0)
            ax.set_title(f"{c}  (n={int(m.sum())})", fontsize=8)
            ax.set_xticks([]); ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
        for j in range(len(cats), nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        fig.suptitle(f"UMAP split by {field} (shared coordinates)", y=1.0)
        fig.tight_layout()
        save_fig(fig, dirs["fig_umap"], f"{prefix}{fname}")

    if failed:
        (dirs["fig_umap"] / f"{prefix}FAILED_plots.txt").write_text(
            "\n".join(failed), encoding="utf-8")
        record("umap_plots_failed", failed)
    log(f"generate_umaps: done ({len(failed)} plot(s) skipped)")


def generate_dotplots(adata: ad.AnnData, dirs: dict, species: str,
                      markers_df: pd.DataFrame, groupby: str = "leiden_cluster",
                      panel: dict | None = None, prefix: str = "") -> None:
    panel = panel or MK.get_panel(species)
    available = set(adata.var_names)
    present, absent = MK.filter_panel(panel, available)

    lines = ["Canonical markers requested but ABSENT from the expression matrix",
             "=" * 70,
             f"dataset gene space: {adata.n_vars} genes", ""]
    for grp, gs in absent.items():
        lines.append(f"{grp}: {', '.join(gs)}")
    if not absent:
        lines.append("(none - every requested canonical marker is present)")
    (dirs["tables"] / f"{prefix}absent_markers.txt").write_text(
        "\n".join(lines), encoding="utf-8")

    layer_note = "mean log1p(CP10K); dot size = fraction of cells expressing"

    # ---- Dot plot A: canonical lineage panel ----------------------------
    dp = sc.pl.dotplot(adata, present, groupby=groupby, show=False,
                       return_fig=True, standard_scale=None,
                       colorbar_title="mean log1p(CP10K)",
                       size_title="fraction expressing")
    fig = figure_of(dp)
    fig.suptitle(f"Canonical lineage markers - {layer_note}", fontsize=9, y=1.01)
    save_fig(fig, dirs["fig_dot"], f"{prefix}canonical_marker_dotplot")

    dp = sc.pl.dotplot(adata, present, groupby=groupby, show=False,
                       return_fig=True, standard_scale="var",
                       colorbar_title="scaled mean expression",
                       size_title="fraction expressing")
    save_fig(figure_of(dp), dirs["fig_dot"],
             f"{prefix}canonical_marker_dotplot_scaled")

    # ---- Dot plot B: data-derived markers -------------------------------
    top = {}
    for cl in adata.obs[groupby].cat.categories if hasattr(
            adata.obs[groupby], "cat") else sorted(
            adata.obs[groupby].astype(str).unique()):
        sub = markers_df[markers_df["cluster"].astype(str) == str(cl)]
        sub = sub[(sub["pval_adj"] < 0.05)
                  & (sub["frac_expressing_in_cluster"] > 0.25)]
        sub = sub[~sub["gene"].map(lambda g: MK.is_uninformative(g, species))]
        sub = sub.sort_values("scores", ascending=False)
        picked = []
        for g in sub["gene"]:
            if g not in top and g not in [x for v in top.values() for x in v]:
                picked.append(g)
            if len(picked) == 3:
                break
        if picked:
            top[str(cl)] = picked
    dp = sc.pl.dotplot(adata, top, groupby=groupby, show=False, return_fig=True,
                       standard_scale="var",
                       colorbar_title="scaled mean expression",
                       size_title="fraction expressing")
    fig = figure_of(dp)
    fig.suptitle("Top data-derived cluster markers (Wilcoxon, adj. p<0.05, "
                 ">25% expressing)", fontsize=9, y=1.01)
    save_fig(fig, dirs["fig_dot"], f"{prefix}cluster_marker_dotplot")

    pd.DataFrame([{"cluster": k, "genes": ", ".join(v)} for k, v in top.items()]
                 ).to_csv(dirs["tables"] / f"{prefix}dotplot_selected_genes.csv",
                          index=False)

    # ---- compartment panels so nothing becomes unreadable ---------------
    for comp, groups in [
        ("epithelial", ["AT2", "AT1", "Transitional", "Basal", "Club",
                        "Ciliated", "Neuroendocrine", "Pan_epithelial"]),
        ("immune", ["Alveolar_Mac", "Mono_Mac", "Neutrophil", "DC", "T_cell",
                    "CD4_T", "CD8_T", "Treg", "NK", "B_cell", "Plasma",
                    "Pan_immune"]),
        ("stromal_endothelial", ["Fibroblast", "Fibro_activated", "SMC_pericyte",
                                 "Mesothelium", "Endothelial",
                                 "Endothelial_cap", "Lymphatic_EC"]),
    ]:
        sel = {g: present[g] for g in groups if g in present}
        if not sel:
            continue
        dp = sc.pl.dotplot(adata, sel, groupby=groupby, show=False,
                           return_fig=True, standard_scale="var",
                           colorbar_title="scaled mean expression",
                           size_title="fraction expressing")
        save_fig(figure_of(dp), dirs["fig_dot"],
                 f"{prefix}marker_dotplot_{comp}")
    log("generate_dotplots: done")


def generate_featureplots(adata: ad.AnnData, dirs: dict, species: str,
                          prefix: str = "") -> None:
    wanted = MK.get_feature_genes(species)
    present = [g for g in wanted if g in adata.var_names]
    absent = [g for g in wanted if g not in adata.var_names]
    (dirs["tables"] / f"{prefix}absent_featureplot_genes.txt").write_text(
        ("Requested feature-plot genes NOT present in the matrix:\n"
         + ("\n".join(absent) if absent else "(none)")), encoding="utf-8")
    if not present:
        return
    ax = sc.pl.umap(adata, color=present, show=False, frameon=False, size=2,
                    cmap="viridis", ncols=4, vmax="p99")
    save_fig(figure_of(ax), dirs["fig_feat"], f"{prefix}featureplots_key_markers")
    for g in present:
        a = sc.pl.umap(adata, color=g, show=False, frameon=False, size=3,
                       cmap="viridis", vmax="p99", title=g)
        save_fig(figure_of(a), dirs["fig_feat"], f"{prefix}feature_{g}",
                 formats=("png",))
    log(f"generate_featureplots: {len(present)} genes plotted, "
        f"{len(absent)} absent and recorded")


# =============================================================================
# PHASE 25 - composition
# =============================================================================
def composition_tables(adata: ad.AnnData, dirs: dict,
                       groupby: str = "leiden_cluster") -> None:
    obs = adata.obs
    per_sample = obs["sample_id"].value_counts().rename_axis("sample_id"
                                                            ).reset_index(name="n_cells")
    per_sample.to_csv(dirs["tables"] / "cells_per_sample.csv", index=False)

    ct = pd.crosstab(obs["sample_id"].astype(str), obs[groupby].astype(str))
    ct.to_csv(dirs["tables"] / "cluster_by_sample_counts.csv")
    pct = ct.div(ct.sum(axis=1), axis=0).mul(100).round(3)
    pct.to_csv(dirs["tables"] / "cluster_percent_within_sample.csv")

    if "condition" in obs:
        ctc = pd.crosstab(obs["condition"].astype(str), obs[groupby].astype(str))
        ctc.to_csv(dirs["tables"] / "cluster_by_condition_counts.csv")
        ctc.div(ctc.sum(axis=1), axis=0).mul(100).round(3).to_csv(
            dirs["tables"] / "cluster_percent_within_condition.csv")
    if "experimental_group" in obs:
        ctg = pd.crosstab(obs["experimental_group"].astype(str),
                          obs[groupby].astype(str))
        ctg.to_csv(dirs["tables"] / "cluster_by_experimental_group_counts.csv")

    # stacked composition per sample
    fig, ax = _plt.subplots(figsize=(max(8, 0.34 * ct.shape[0] + 2), 5))
    bottom = np.zeros(ct.shape[0])
    cmap = _plt.get_cmap("tab20")
    for i, cl in enumerate(pct.columns):
        ax.bar(pct.index, pct[cl].to_numpy(), bottom=bottom,
               color=cmap(i % 20), width=0.85, label=cl, linewidth=0)
        bottom += pct[cl].to_numpy()
    ax.set_ylabel("% of cells in sample"); ax.set_ylim(0, 100)
    ax.set_xticks(range(len(pct.index)))
    ax.set_xticklabels(pct.index, rotation=90, fontsize=6)
    ax.legend(ncol=2, fontsize=5.5, bbox_to_anchor=(1.005, 1), loc="upper left",
              title="cluster")
    ax.set_title("Cluster composition per sample")
    fig.tight_layout()
    save_fig(fig, dirs["fig_comp"], "composition_stacked_per_sample")

    # per-sample cluster proportions keep replication visible
    if "condition" in obs:
        smap = obs.groupby("sample_id", observed=True)["condition"].agg(
            lambda s: s.astype(str).mode().iat[0])
        long = pct.reset_index().melt(id_vars=pct.index.name or "index",
                                      var_name="cluster", value_name="pct")
        long = long.rename(columns={long.columns[0]: "sample_id"})
        long["condition"] = long["sample_id"].map(smap)
        long.to_csv(dirs["tables"] / "cluster_proportion_per_sample_long.csv",
                    index=False)

        clusters = list(pct.columns)
        ncol = min(6, len(clusters))
        nrow = int(np.ceil(len(clusters) / ncol))
        fig, axes = _plt.subplots(nrow, ncol, figsize=(2.5 * ncol, 2.3 * nrow),
                                  squeeze=False)
        conds = sorted(long["condition"].dropna().unique())
        for i, cl in enumerate(clusters):
            ax = axes[i // ncol][i % ncol]
            d = long[long["cluster"] == cl]
            for j, c in enumerate(conds):
                v = d.loc[d["condition"] == c, "pct"].to_numpy()
                ax.scatter(np.full(v.shape, j) + np.random.uniform(-.12, .12, v.size),
                           v, s=10, alpha=0.8)
            ax.set_xticks(range(len(conds)))
            ax.set_xticklabels(conds, rotation=45, fontsize=5, ha="right")
            ax.set_title(f"cluster {cl}", fontsize=7)
            ax.tick_params(labelsize=6)
        for j in range(len(clusters), nrow * ncol):
            axes[j // ncol][j % ncol].axis("off")
        fig.suptitle("Cluster proportion per SAMPLE (one point = one sample)",
                     fontsize=9)
        fig.tight_layout()
        save_fig(fig, dirs["fig_comp"], "composition_per_sample_by_condition")
    log("composition_tables: done")


# =============================================================================
# PHASE 32 - automatic QC review
# =============================================================================
def donor_driven_clustering_check(adata: ad.AnnData, dirs: dict,
                                  groupby: str = "leiden_cluster") -> dict:
    """Is the embedding separating cell types, or separating samples?

    The honest test is not "do samples separate" - they may differ
    biologically - but "does one cell type get split into several clusters,
    each belonging to a different sample".  A cell type that fragments along
    sample lines is being driven by donor/technical offset, because the same
    cell type is not four different cell types in four animals or donors.
    """
    out: dict = {"applicable": False}
    if "proposed_cell_type" not in adata.obs or adata.obs["sample_id"].nunique() < 2:
        return out

    ct = pd.crosstab(adata.obs[groupby].astype(str),
                     adata.obs["sample_id"].astype(str))
    purity = ct.div(ct.sum(axis=1), axis=0).max(axis=1)
    ident = (adata.obs.groupby(adata.obs[groupby].astype(str), observed=True)
             ["proposed_cell_type"].agg(lambda s: s.astype(str).mode().iat[0]))

    rows, n_split, n_multi = [], 0, 0
    for name, clusters in ident.groupby(ident).groups.items():
        clusters = list(clusters)
        if len(clusters) < 2:
            continue
        n_multi += 1
        private = [c for c in clusters if purity.get(c, 0) > 0.75]
        if len(private) >= 2:
            n_split += 1
        rows.append({"proposed_identity": name, "n_clusters": len(clusters),
                     "clusters": ", ".join(map(str, clusters)),
                     "n_sample_private_clusters": len(private),
                     "max_donor_purity": round(float(
                         max(purity.get(c, 0) for c in clusters)), 3)})

    df = pd.DataFrame(rows)
    if not df.empty:
        df.sort_values("n_clusters", ascending=False).to_csv(
            dirs["qc"] / "celltype_split_by_sample.csv", index=False)

    frac = n_split / n_multi if n_multi else 0.0
    out = {"applicable": True,
           "identities_spanning_multiple_clusters": n_multi,
           "identities_split_into_sample_private_clusters": n_split,
           "fraction_split_by_sample": round(frac, 3),
           "donor_driven": bool(n_multi >= 3 and frac >= 0.5)}
    (dirs["qc"] / "donor_driven_clustering_check.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    record("donor_driven_clustering_check", out)
    return out


def qc_review(adata: ad.AnnData, dirs: dict, groupby: str = "leiden_cluster"
              ) -> list[str]:
    obs = adata.obs
    g = obs.groupby(obs[groupby].astype(str), observed=True)
    summary = pd.DataFrame({
        "n_cells": g.size(),
        "median_pct_mt": g["pct_counts_mt"].median(),
        "median_total_counts": g["total_counts"].median(),
        "median_n_genes": g["n_genes_by_counts"].median(),
        "mean_doublet_score": g["doublet_score"].mean(),
        "n_samples_present": g["sample_id"].nunique(),
    })
    dom = (pd.crosstab(obs[groupby].astype(str), obs["sample_id"].astype(str))
           .pipe(lambda d: d.div(d.sum(axis=1), axis=0)))
    summary["max_single_sample_fraction"] = dom.max(axis=1)
    summary["dominant_sample"] = dom.idxmax(axis=1)
    summary.to_csv(dirs["qc"] / "cluster_qc_summary.csv")

    warn: list[str] = []
    gm_mt = obs["pct_counts_mt"].median()
    gm_ct = obs["total_counts"].median()
    for cl, row in summary.iterrows():
        if row["median_pct_mt"] > max(2 * gm_mt, 10):
            warn.append(f"cluster {cl}: median mitochondrial fraction "
                        f"{row['median_pct_mt']:.1f}% is far above the dataset "
                        f"median ({gm_mt:.1f}%) - possible dying/low-quality cells")
        if row["median_total_counts"] < 0.4 * gm_ct:
            warn.append(f"cluster {cl}: median depth {row['median_total_counts']:.0f} "
                        f"UMI is well below the dataset median ({gm_ct:.0f}) - "
                        f"check whether this cluster is driven by depth")
        if row["mean_doublet_score"] == row["mean_doublet_score"] and \
                row["mean_doublet_score"] > 0.25:
            warn.append(f"cluster {cl}: mean Scrublet score "
                        f"{row['mean_doublet_score']:.2f} - possible doublet cluster")
        if row["max_single_sample_fraction"] > 0.75 and row["n_cells"] > 100:
            warn.append(f"cluster {cl}: {100 * row['max_single_sample_fraction']:.0f}% "
                        f"of its cells come from a single sample "
                        f"({row['dominant_sample']}) - may be sample-specific "
                        f"rather than a shared population")
    # Aggregate view: a handful of sample-dominated clusters is normal, but if
    # most clusters are private to one sample the embedding is being driven by
    # sample identity rather than by cell type, and that needs saying once and
    # clearly rather than as N separate per-cluster lines.
    dd = donor_driven_clustering_check(adata, dirs, groupby)
    if dd.get("donor_driven"):
        warn.insert(0, (
            f"OVERALL: {dd['identities_split_into_sample_private_clusters']} of "
            f"{dd['identities_spanning_multiple_clusters']} proposed cell types "
            f"that span more than one cluster are split into clusters each "
            f"belonging to a different sample. The same cell type is not "
            f"several different cell types in several donors/animals, so this "
            f"embedding is being driven by sample identity rather than by cell "
            f"type. If no experimental contrast would be destroyed by it, "
            f"batch correction is warranted - see qc/celltype_split_by_sample.csv."))

    n_priv = int((summary["max_single_sample_fraction"] > 0.75).sum())
    n_tot = len(summary)
    frac_priv = n_priv / max(n_tot, 1)
    n_samples = obs["sample_id"].nunique()
    if n_samples > 1 and frac_priv > 0.30:
        headline = (
            f"OVERALL: {n_priv} of {n_tot} clusters ({100 * frac_priv:.0f}%) draw "
            f">75% of their cells from a single sample. Cell types are normally "
            f"shared between samples, so an embedding this sample-private is "
            f"being driven substantially by sample identity. See "
            f"qc/batch_assessment.json for the mixing statistics and the "
            f"integration decision, and figures/umap/"
            f"UMAP_integration_before_after.* if a supplementary corrected "
            f"embedding was produced.")
        warn.insert(0, headline)
    record("sample_private_clusters",
           f"{n_priv}/{n_tot} clusters are >75% one sample")

    (dirs["qc"] / "qc_review_warnings.txt").write_text(
        "\n".join(warn) if warn else "no automatic QC warnings raised",
        encoding="utf-8")
    for w in warn:
        log(f"  QC WARNING: {w}")
    return warn


# =============================================================================
# PHASE 24 - epithelial sub-analysis
# =============================================================================
def epithelial_subanalysis(adata: ad.AnnData, dirs: dict, species: str,
                           n_hvg: int = 2000) -> dict | None:
    ev = MK.EPITHELIAL_EVIDENCE[species]
    pos = [g for g in ev["positive"] if g in adata.var_names]
    neg = [g for g in ev["negative"] if g in adata.var_names]
    if len(pos) < 3:
        return None

    mean_p, frac_p = _cluster_mean_expression(adata, pos, "leiden_cluster")
    mean_n, frac_n = _cluster_mean_expression(adata, neg, "leiden_cluster")
    zp = ((mean_p - mean_p.mean()) / mean_p.std().replace(0, np.nan)).fillna(0)
    zn = ((mean_n - mean_n.mean()) / mean_n.std().replace(0, np.nan)).fillna(0)
    score = zp.mean(axis=1) - zn.mean(axis=1)

    pan = [g for g in (["Epcam", "Cdh1"] if species == "mouse"
                       else ["EPCAM", "CDH1"]) if g in frac_p.columns]
    pan_frac = frac_p[pan].mean(axis=1) if pan else pd.Series(0.0, index=score.index)

    # Coherent evidence, not a single gene: a cluster qualifies if the
    # epithelial panel clearly outweighs the immune/endothelial/fibroblast
    # panels AND either a pan-epithelial marker is broadly detected or the
    # epithelial signal is strong on its own.  The second branch matters
    # because EPCAM/CDH1 capture is genuinely poor in some epithelial states
    # (ciliated cells especially), and requiring it outright drops real
    # epithelium.
    sel = score[(score > 0.5)
                & ((pan_frac > 0.25) | (score > 0.8))].index.tolist()
    tbl = pd.DataFrame({"cluster": score.index, "epithelial_score": score.round(3),
                        "pan_epithelial_fraction": pan_frac.round(3),
                        "selected": [c in sel for c in score.index]})
    tbl.to_csv(dirs["epi_tab"] / "epithelial_cluster_selection.csv", index=False)
    if not sel:
        record("epithelial_subanalysis", "not performed - no cluster met the "
                                         "combined epithelial evidence criterion")
        return None

    record("epithelial_subanalysis_selection",
           f"clusters {sel} selected on combined evidence: mean z of {pos} "
           f"minus mean z of {neg} > 0.5, AND either >25% of cells expressing "
           f"a pan-epithelial marker or an epithelial score > 0.8. Selection "
           f"never rests on a single gene.")

    epi = adata[adata.obs["leiden_cluster"].astype(str).isin(sel)].copy()
    log(f"epithelial_subanalysis: {epi.n_obs} cells in {len(sel)} clusters")
    if epi.n_obs < 300:
        record("epithelial_subanalysis", f"skipped - only {epi.n_obs} cells")
        return None

    sc.pp.highly_variable_genes(epi, flavor="seurat", n_top_genes=n_hvg,
                                batch_key="sample_id" if
                                epi.obs["sample_id"].nunique() > 1 else None)
    sub = epi[:, epi.var["highly_variable"].to_numpy()].copy()
    sc.pp.scale(sub, max_value=10)
    sc.pp.pca(sub, n_comps=min(50, sub.n_obs - 1, sub.n_vars - 1),
              svd_solver="arpack", random_state=RANDOM_SEED)
    epi.obsm["X_pca"] = sub.obsm["X_pca"]
    epi.uns["pca"] = sub.uns["pca"]
    del sub
    free_mem()

    n_pcs = min(30, epi.obsm["X_pca"].shape[1])
    sc.pp.neighbors(epi, n_neighbors=15, n_pcs=n_pcs, random_state=RANDOM_SEED)
    sc.tl.umap(epi, random_state=RANDOM_SEED)
    sc.tl.leiden(epi, resolution=0.6, key_added="leiden_cluster",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=RANDOM_SEED)
    n_cl = epi.obs["leiden_cluster"].nunique()

    edirs = dict(dirs)
    edirs["fig_umap"] = dirs["epi_fig"]
    edirs["fig_dot"] = dirs["epi_fig"]
    edirs["fig_feat"] = dirs["epi_fig"]
    edirs["fig_qc"] = dirs["epi_fig"]
    edirs["fig_comp"] = dirs["epi_fig"]
    edirs["tables"] = dirs["epi_tab"]

    mk = find_markers(epi, edirs, "leiden_cluster", prefix="epithelial_cluster")
    generate_umaps(epi, edirs, prefix="epithelial_")
    generate_dotplots(epi, edirs, species, mk, panel=MK.get_epi_panel(species),
                      prefix="epithelial_")
    generate_featureplots(epi, edirs, species, prefix="epithelial_")
    annotate_clusters(epi, mk, edirs, species,
                      outfile="epithelial_cluster_annotation_proposals.csv")

    epi.obs.to_csv(dirs["epi_tab"] / "epithelial_cell_metadata.csv")
    _tidy_obs_for_h5ad(epi)
    epi.write_h5ad(dirs["epi"] / "epithelial_clustered.h5ad", compression="gzip")
    record("epithelial_subanalysis",
           f"performed on {epi.n_obs} cells from clusters {sel}; "
           f"{n_cl} epithelial subclusters at Leiden resolution 0.6")
    res = {"n_cells": int(epi.n_obs), "clusters": sel, "n_subclusters": int(n_cl)}
    del epi
    free_mem()
    return res


# =============================================================================
# PHASE 26 - exploratory pseudobulk
# =============================================================================
def pseudobulk_condition(dirs: dict, groupby: str = "leiden_cluster") -> None:
    """Aggregate counts per sample so that replication stays at sample level."""
    post = dirs["processed"] / "postQC.h5ad"
    final_obs = dirs["tables"] / "cell_metadata.csv"
    if not post.exists() or not final_obs.exists():
        return
    log("pseudobulk_condition: aggregating counts per sample")
    obs = pd.read_csv(final_obs, index_col=0)
    a = ad.read_h5ad(post)
    common = a.obs_names.intersection(obs.index)
    a = a[common].copy()
    labels = obs.loc[common]

    samples = labels["sample_id"].astype(str)
    cats = sorted(samples.unique())
    M = np.zeros((len(cats), a.n_vars), dtype=np.float64)
    X = a.X.tocsr()
    for i, s in enumerate(cats):
        M[i] = np.asarray(X[(samples == s).to_numpy()].sum(axis=0)).ravel()
    pb = pd.DataFrame(M, index=cats, columns=a.var_names)
    pb.to_csv(dirs["tables"] / "pseudobulk_counts_per_sample.csv.gz",
              compression="gzip")

    if "condition" in labels.columns:
        cond = labels.groupby("sample_id", observed=True)["condition"].agg(
            lambda s: s.astype(str).mode().iat[0]).reindex(cats)
        n_ctrl = int((cond == "homeostasis").sum())
        n_inf = int((cond == "H1N1_infected").sum())
        note = (f"EXPLORATORY ONLY. Sample-level pseudobulk with n={n_ctrl} "
                f"homeostasis and n={n_inf} infected samples. The homeostasis "
                f"arm is not replicated within a group and is confounded with "
                f"sex, so this is not a confirmatory differential-expression "
                f"result.")
        if n_ctrl >= 2 and n_inf >= 2:
            from scipy import stats
            cpm = np.log1p(pb.div(pb.sum(axis=1), axis=0) * 1e6)
            a_ = cpm[cond.to_numpy() == "homeostasis"]
            b_ = cpm[cond.to_numpy() == "H1N1_infected"]
            keep = (pb > 0).sum(axis=0) >= max(3, len(cats) // 4)
            t, p = stats.ttest_ind(a_.loc[:, keep], b_.loc[:, keep],
                                   equal_var=False)
            from statsmodels.stats.multitest import multipletests
            ok = ~np.isnan(p)
            padj = np.full(p.shape, np.nan)
            padj[ok] = multipletests(p[ok], method="fdr_bh")[1]
            res = pd.DataFrame({
                "gene": cpm.columns[keep], "t": t, "pval": p, "pval_adj": padj,
                "mean_logCPM_homeostasis": a_.loc[:, keep].mean().to_numpy(),
                "mean_logCPM_H1N1": b_.loc[:, keep].mean().to_numpy()})
            res["log2FC_H1N1_vs_homeostasis"] = (
                res["mean_logCPM_H1N1"] - res["mean_logCPM_homeostasis"])
            res.sort_values("pval").to_csv(
                dirs["tables"] / "pseudobulk_exploratory_H1N1_vs_homeostasis.csv",
                index=False)
        (dirs["tables"] / "pseudobulk_README.txt").write_text(note,
                                                              encoding="utf-8")
        record("condition_differential_expression", note)
    del a
    free_mem()


# =============================================================================
# save + report
# =============================================================================
def attach_counts_layer_on_disk(final_path: Path, counts_path: Path) -> bool:
    """Copy X from the counts checkpoint into layers/counts of the final file.

    Done at the HDF5 level so that two full matrices are never resident in RAM
    at the same time.
    """
    import h5py

    def _index_of(f) -> np.ndarray:
        """Read the obs index, which anndata may store either as a plain
        dataset or as a nullable-string-array group of values + mask."""
        obs = f["obs"]
        node = obs[obs.attrs["_index"]]
        if isinstance(node, h5py.Group):
            return np.asarray(node["values"][:])
        return np.asarray(node[:])

    try:
        with h5py.File(final_path, "r") as f:
            fin_names = _index_of(f)
        with h5py.File(counts_path, "r") as f:
            cnt_names = _index_of(f)
        if fin_names.shape != cnt_names.shape or not np.array_equal(fin_names,
                                                                    cnt_names):
            log("  cell order differs between final and counts checkpoint - "
                "counts layer not attached (postQC.h5ad still holds the counts)")
            return False
        with h5py.File(counts_path, "r") as src, \
                h5py.File(final_path, "a") as dst:
            if "layers" not in dst:
                grp = dst.create_group("layers")
                grp.attrs["encoding-type"] = "dict"
                grp.attrs["encoding-version"] = "0.1.0"
            if "counts" in dst["layers"]:
                del dst["layers"]["counts"]
            src.copy(src["X"], dst["layers"], name="counts")
        log("  counts layer attached to the final object on disk")
        return True
    except Exception as exc:  # noqa: BLE001
        log(f"  could not attach counts layer ({exc}); postQC.h5ad still holds "
            f"the raw counts")
        return False


def save_outputs(adata: ad.AnnData, dirs: dict) -> None:
    keep = [c for c in [
        "sample_id", "gsm", "original_barcode", "orig_ident",
        "experimental_group", "condition", "tamoxifen_start_day",
        "sacrifice_day", "sex", "trace_call", "cell_cycle_phase",
        "has_author_metadata", "author_lineage", "author_celltype",
        "total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo",
        "pct_counts_hb", "pct_counts_in_top_20_genes", "doublet_score",
        "predicted_doublet", "leiden_cluster", "proposed_cell_type",
        "reporter_SiteA", "reporter_SiteB",
    ] + [c for c in adata.obs.columns if c.startswith("leiden_res")]
        if c in adata.obs.columns]
    meta = adata.obs[keep].copy()
    meta.insert(0, "cell_id", adata.obs_names)
    um = np.asarray(adata.obsm["X_umap"])
    meta["UMAP_1"] = um[:, 0]
    meta["UMAP_2"] = um[:, 1]
    meta.to_csv(dirs["tables"] / "cell_metadata.csv", index=False)
    log(f"save_outputs: cell_metadata.csv with {len(meta)} cells")

    _tidy_obs_for_h5ad(adata)
    final = dirs["processed"] / "final_clustered.h5ad"
    adata.write_h5ad(final, compression="gzip")
    log(f"save_outputs: {final}")
    attach_counts_layer_on_disk(final, dirs["processed"] / "postQC.h5ad")


def write_analysis_log(dirs: dict, cfg: DatasetConfig, fmt: dict,
                       warnings_list: list[str]) -> None:
    d = DECISIONS
    L = [
        "=" * 70, f"ANALYSIS LOG - {cfg.name}", "=" * 70, "",
        "=== DATA DISCOVERY ===",
        f"Raw data directory: {_repo_rel(cfg.raw_dir)}",
        f"Input format: {d.get('input_format')}",
        f"Species: {d.get('species')}",
        f"Number of samples: {fmt.get('n_sample_files')}",
        f"Samples detected: {d.get('samples_detected')}",
        f"Metadata detected: {_repo_rel(d.get('metadata_table'))}",
        f"Counts appear raw/normalized/processed: "
        f"{'raw integer UMI counts' if d.get('counts_are_raw_integers') else 'not integer'}",
        f"Non-gene features removed: {fmt.get('non_gene_features_present')}",
        f"Gene space identical across samples: "
        f"{d.get('gene_space_identical_across_samples')}",
        f"Gene space join: {d.get('gene_space_join')}", "",
        "=== QC ===",
        f"Initial cells: {d.get('cells_loaded')}",
        f"Initial genes: {d.get('genes_loaded')}",
        f"Per-sample cell counts: qc/qc_before_after.csv",
        f"QC thresholds: qc/qc_thresholds.csv (per sample, MAD-derived)",
        f"Cells removed by QC: {d.get('cells_removed_qc')}",
        f"Cells retained after QC: {d.get('cells_after_qc')}",
        f"Doublets removed: {d.get('doublets_removed')}",
        f"Genes retained after min_cells filter: {d.get('genes_after_filter')}", "",
        "=== NORMALIZATION ===",
        f"Method: {d.get('normalization')}",
        f"HVG method: {d.get('hvg_method')}",
        f"Number of HVGs: {d.get('n_hvg')}",
        f"Scaling: {d.get('scaling')}", "",
        "=== DIMENSION REDUCTION ===",
        f"Number of PCs: {d.get('n_pcs')}",
        f"Neighbor parameters: {d.get('neighbors')}",
        f"UMAP parameters: {d.get('umap')}", "",
        "=== BATCH HANDLING ===",
        f"Batch variable considered: sample_id",
        f"Assessment: {json.dumps(d.get('batch_assessment'), indent=2)}",
        f"Correction applied: {d.get('batch_correction')}", "",
        "=== CLUSTERING ===",
        f"Algorithm: {d.get('clustering')}",
        f"Resolution: {d.get('leiden_resolution')}",
        f"Number of clusters: {d.get('n_clusters')}", "",
        "=== ANNOTATION ===",
        f"Strategy: {d.get('annotation_strategy')}",
        f"Author-label agreement: {d.get('author_label_agreement')}",
        f"Epithelial sub-analysis: {d.get('epithelial_subanalysis')}", "",
        "=== OUTPUT ===",
        "Main UMAP: figures/umap/UMAP_leiden_clusters.{png,pdf}",
        "Sample UMAP: figures/umap/UMAP_sample.{png,pdf}",
        "Dot plots: figures/dotplots/canonical_marker_dotplot.{png,pdf}, "
        "figures/dotplots/cluster_marker_dotplot.{png,pdf}",
        "Marker tables: tables/cluster_markers_all.csv, "
        "tables/cluster_markers_top20.csv",
        "Cell metadata: tables/cell_metadata.csv",
        "Processed object: processed/final_clustered.h5ad", "",
        "=== AUTOMATIC QC REVIEW ===",
    ]
    L += [f"  - {w}" for w in warnings_list] or ["  (no warnings raised)"]
    L += ["", "=== ALL RECORDED DECISIONS ===",
          json.dumps({k: str(v) for k, v in d.items()}, indent=2)]
    (dirs["logs"] / "analysis_log.txt").write_text("\n".join(L), encoding="utf-8")
    save_decisions()
    log(f"write_analysis_log: {dirs['logs'] / 'analysis_log.txt'}")


# =============================================================================
# driver
# =============================================================================
def main() -> int:
    valid_stages = {"samples", "merge", "norm", "pca", "cluster", "markers",
                    "figures", "epi", "finalize", "log"}
    ap = argparse.ArgumentParser(
        description="Run the reproducible scanpy workflow for one GEO series.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--dataset", required=True, choices=sorted(DATASETS),
                    help="GEO series to analyse")
    ap.add_argument(
        "--stages", default="all",
        help="comma-separated stages, or 'all'; valid stages: "
             + ", ".join(sorted(valid_stages)),
    )
    ap.add_argument("--workers", type=int, default=3,
                    help="parallel workers for per-sample processing")
    ap.add_argument("--n-hvg", type=int, default=2500,
                    help="number of highly variable genes")
    ap.add_argument("--integration", default=None,
                    choices=["none", "harmony"],
                    help="override the evidence-based integration decision")
    ap.add_argument("--outdir", default=None,
                    help="output directory (defaults to analysis/<dataset>)")
    ap.add_argument("--reuse-markers", action="store_true",
                    help="re-use the cluster marker table written by a previous "
                         "run instead of recomputing the Wilcoxon test; only "
                         "honoured when its clusters match the current object")
    args = ap.parse_args()

    cfg = DATASETS[args.dataset](Path(args.outdir) if args.outdir else None)
    dirs = ensure_dirs(cfg)
    stages = (valid_stages - {"log"} if args.stages == "all"
              else {stage.strip() for stage in args.stages.split(",")
                    if stage.strip()})
    unknown_stages = stages - valid_stages
    if not stages:
        ap.error("--stages must name at least one stage")
    if unknown_stages:
        ap.error("unknown --stages value(s): "
                 + ", ".join(sorted(unknown_stages)))

    load_decisions(dirs["logs"] / "decisions.json")
    log(f"=== {cfg.name} === outdir={cfg.outdir}")
    log(f"scanpy {sc.__version__}, anndata {ad.__version__}, "
        f"numpy {np.__version__} {mem_report()}")

    tasks = scan_input_directory(cfg)
    fmt = detect_input_format(cfg, tasks)
    (dirs["inventory"] / "detected_format.json").write_text(
        json.dumps(fmt, indent=2, default=str), encoding="utf-8")
    record("samples_detected", [t["sample"] for t in tasks])

    if "log" in stages:
        # Rewrite the analysis log from the persisted decisions without
        # recomputing anything.  Useful after a run that was resumed stage by
        # stage, where the log on disk reflects only the final stage.
        wp = dirs["qc"] / "qc_review_warnings.txt"
        prior_warn = [w for w in (wp.read_text(encoding="utf-8").splitlines()
                                  if wp.exists() else [])
                      if w.strip() and not w.startswith("no automatic")]
        write_analysis_log(dirs, cfg, fmt, prior_warn)
        log("rewrote the analysis log from persisted decisions")
        if stages == {"log"}:
            return 0

    if "samples" in stages:
        load_samples(cfg, tasks, dirs, args.workers)
        calculate_qc(dirs, cfg)

    adata = None
    if "merge" in stages:
        adata = merge_samples(cfg, dirs)
        record("cells_loaded", int(pd.read_csv(
            dirs["qc"] / "qc_before_after.csv")["Cells before QC"].iloc[-1]))
        record("genes_loaded", int(adata.n_vars))
        record("cells_after_qc", int(adata.n_obs))
        record("cells_removed_qc",
               int(DECISIONS["cells_loaded"]) - int(adata.n_obs))

        n_before = adata.n_obs
        adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
        record("doublets_removed", f"{n_before - adata.n_obs} cells called by "
                                   f"Scrublet run independently per capture; "
                                   f"scores retained in the metadata for every "
                                   f"cell")
        sc.pp.filter_genes(adata, min_cells=3)
        record("genes_after_filter", int(adata.n_vars))
        adata = merge_metadata(cfg, adata, dirs)
        _tidy_obs_for_h5ad(adata)
        adata.write_h5ad(dirs["processed"] / "merged_postQC.h5ad",
                         compression="gzip")

    if "norm" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "merged_postQC.h5ad")
        adata = normalize_data(adata, dirs, args.n_hvg)
        adata.write_h5ad(dirs["processed"] / "normalized.h5ad", compression="gzip")

    if "pca" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "normalized.h5ad")
        n_pcs = run_pca(adata, dirs)
        assess = inspect_batch(adata, n_pcs, dirs)
        primary, secondary = run_integration_if_needed(
            adata, assess, n_pcs, args.integration)
        build_embeddings(adata, n_pcs, primary, secondary, dirs)
        adata.write_h5ad(dirs["processed"] / "neighbors_umap.h5ad",
                         compression="gzip")

    if "cluster" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "neighbors_umap.h5ad")
        run_clustering(adata, dirs, cfg.species)
        adata.write_h5ad(dirs["processed"] / "clustered.h5ad", compression="gzip")

    mk = None
    if "markers" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "clustered.h5ad")
        mk = find_markers(adata, dirs, reuse=args.reuse_markers)
        annotate_clusters(adata, mk, dirs, cfg.species)
        validate_against_author_labels(adata, dirs)

    warn: list[str] = []
    if "figures" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "clustered.h5ad")
        if mk is None:
            mk = pd.read_csv(dirs["tables"] / "cluster_markers_all.csv")
        generate_umaps(adata, dirs)
        generate_dotplots(adata, dirs, cfg.species, mk)
        generate_featureplots(adata, dirs, cfg.species)
        composition_tables(adata, dirs)
        warn = qc_review(adata, dirs)

    if "finalize" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "clustered.h5ad")
        save_outputs(adata, dirs)

    if "epi" in stages:
        if adata is None:
            adata = _load_checkpoint(dirs["processed"] / "clustered.h5ad")
        epithelial_subanalysis(adata, dirs, cfg.species)

    if "finalize" in stages:
        del adata
        free_mem()
        pseudobulk_condition(dirs)
        write_analysis_log(dirs, cfg, fmt, warn)

    log(f"DONE {cfg.name} {mem_report()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
