"""
Per-sample stage of the pipeline, isolated in its own module so that it can be
executed in a separate process on Windows (spawn start method).

One worker handles exactly one 10x HDF5 file:

  load -> drop non-gene features -> QC metrics -> data-driven thresholds
       -> filter -> Scrublet on that capture alone -> write an .h5ad shard

Doing doublet detection here, before any merging, is deliberate: Scrublet
simulates doublets by combining pairs of observed transcriptomes, and pairs may
only be drawn from cells that could physically have been co-encapsulated, i.e.
from the same 10x capture.
"""

from __future__ import annotations

import os

# Cap BLAS threads BEFORE numpy is imported.  Each worker process otherwise
# grabs every core for its linear algebra, so N workers oversubscribe the
# machine N-fold and spend their time in contention rather than arithmetic.
# Must be set at import time; numpy reads these once.
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "2")

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from pipeline_utils import (RANDOM_SEED, SampleThresholds, apply_thresholds,
                            derive_thresholds, make_unique)
import markers as MK


def _load_10x_h5(path: Path, use_ids_as_index: bool) -> "anndata.AnnData":
    import anndata as ad
    import h5py

    with h5py.File(path, "r") as f:
        g = f["matrix"]
        shape = tuple(int(x) for x in g["shape"][:])          # (n_genes, n_cells)
        data = np.asarray(g["data"][:], dtype=np.float32)
        indices = np.asarray(g["indices"][:], dtype=np.int32)
        indptr = np.asarray(g["indptr"][:], dtype=np.int64)
        barcodes = np.array([b.decode() for b in g["barcodes"][:]], dtype=object)

        feats = g["features"]
        names = np.array([b.decode() for b in feats["name"][:]], dtype=object)
        ids = (np.array([b.decode() for b in feats["id"][:]], dtype=object)
               if "id" in feats else names.copy())
        ftype = (np.array([b.decode() for b in feats["feature_type"][:]],
                          dtype=object) if "feature_type" in feats
                 else np.array(["Gene Expression"] * len(names), dtype=object))
        genome = (np.array([b.decode() for b in feats["genome"][:]], dtype=object)
                  if "genome" in feats else np.array([""] * len(names),
                                                     dtype=object))

    # The file stores the matrix CSC with cells as columns; transposing that
    # CSC gives a CSR with cells as rows without touching the data arrays.
    X = sp.csr_matrix((data, indices, indptr), shape=(shape[1], shape[0]))
    X.eliminate_zeros()          # these files carry explicitly stored zeros

    var = pd.DataFrame({"gene_symbol": names, "gene_id": ids,
                        "feature_type": ftype, "genome": genome})
    if use_ids_as_index:
        var.index = pd.Index(np.asarray(ids, dtype=str), name=None)
    else:
        var.index = pd.Index(make_unique(names).astype(str), name=None)

    obs = pd.DataFrame(index=pd.Index(np.asarray(barcodes, dtype=str)))
    adata = ad.AnnData(X=X, obs=obs, var=var)
    return adata


def process_sample(task: dict) -> dict:
    """Run the whole per-sample stage.  `task` is a plain dict so that it
    pickles cleanly across the process boundary."""
    warnings.filterwarnings("ignore")
    import scanpy as sc

    sc.settings.verbosity = 0
    path = Path(task["path"])
    sample = task["sample"]
    species = task["species"]
    shard_dir = Path(task["shard_dir"])
    non_gene = task["non_gene_features"]
    use_ids = task["use_gene_ids_as_index"]

    result: dict = {"sample": sample, "path": str(path)}

    adata = _load_10x_h5(path, use_ids)
    result["n_barcodes_raw"] = int(adata.n_obs)
    result["n_features_raw"] = int(adata.n_vars)
    result["gsm"] = task.get("gsm", "")

    # ---- lineage-tracing reporter contigs are not genes ------------------
    if non_gene:
        present = [g for g in non_gene if g in adata.var_names]
        if present:
            rep = adata[:, present].X
            rep = np.asarray(rep.todense()) if sp.issparse(rep) else np.asarray(rep)
            for i, g in enumerate(present):
                adata.obs[f"reporter_{g}"] = rep[:, i].astype(np.float32)
            adata = adata[:, ~adata.var_names.isin(present)].copy()
            result["non_gene_features_removed"] = present

    # keep only Gene Expression features if the file mixes modalities
    if "feature_type" in adata.var:
        gex = adata.var["feature_type"].astype(str) == "Gene Expression"
        if not gex.all():
            result["non_gex_features_removed"] = int((~gex).sum())
            adata = adata[:, gex.to_numpy()].copy()

    # ---- QC metrics ------------------------------------------------------
    mt_pref = MK.mito_prefix(species)
    ribo_pref = MK.ribo_prefixes(species)
    hb_pref = MK.hb_prefixes(species)
    sym = adata.var["gene_symbol"].astype(str)

    adata.var["mt"] = sym.str.startswith(mt_pref).to_numpy()
    # exclude ribosomal pseudogenes: they inflate the fraction without meaning
    ribo_hit = sym.str.startswith(ribo_pref) & ~sym.str.contains(r"-ps", case=False,
                                                                regex=True)
    adata.var["ribo"] = ribo_hit.to_numpy()
    adata.var["hb"] = sym.str.startswith(hb_pref).to_numpy()

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20],
        log1p=True, inplace=True)

    result["n_mt_genes"] = int(adata.var["mt"].sum())
    result["n_ribo_genes"] = int(adata.var["ribo"].sum())
    result["n_hb_genes"] = int(adata.var["hb"].sum())

    obs = adata.obs
    pre = pd.DataFrame({
        "sample_id": sample,
        "barcode": adata.obs_names.to_numpy(),
        "total_counts": obs["total_counts"].to_numpy(),
        "n_genes_by_counts": obs["n_genes_by_counts"].to_numpy(),
        "pct_counts_mt": obs["pct_counts_mt"].to_numpy(),
        "pct_counts_ribo": obs["pct_counts_ribo"].to_numpy(),
        "pct_counts_hb": obs["pct_counts_hb"].to_numpy(),
        "pct_counts_in_top_20_genes": obs["pct_counts_in_top_20_genes"].to_numpy(),
    })

    # ---- data-driven thresholds -----------------------------------------
    th = derive_thresholds(obs, sample)
    keep = apply_thresholds(obs, th)
    pre["passes_qc"] = keep.to_numpy()

    # per-reason bookkeeping so nothing is removed silently
    reasons = {
        "low_counts": int((obs["total_counts"] < th.min_counts).sum()),
        "high_counts": int((obs["total_counts"] > th.max_counts).sum()),
        "low_genes": int((obs["n_genes_by_counts"] < th.min_genes).sum()),
        "high_genes": int((obs["n_genes_by_counts"] > th.max_genes).sum()),
        "high_mito": int((obs["pct_counts_mt"] > th.max_pct_mt).sum()),
        "high_top20": int(
            (obs["pct_counts_in_top_20_genes"] > th.max_pct_top20).sum()),
    }
    result["removal_reasons"] = reasons
    result["thresholds"] = {
        "sample": sample, "min_genes": th.min_genes, "max_genes": th.max_genes,
        "min_counts": th.min_counts, "max_counts": th.max_counts,
        "max_pct_mt": round(th.max_pct_mt, 4),
        "max_pct_top20": round(th.max_pct_top20, 4),
        "n_mads_counts": th.n_mads_counts, "n_mads_genes": th.n_mads_genes,
        "n_mads_mt": th.n_mads_mt, "rationale": th.rationale,
    }
    result["n_cells_before_qc"] = int(adata.n_obs)
    result["n_cells_after_qc"] = int(keep.sum())

    adata = adata[keep.to_numpy()].copy()

    # ---- doublets, per capture ------------------------------------------
    dbl_score = np.full(adata.n_obs, np.nan, dtype=np.float32)
    dbl_call = np.zeros(adata.n_obs, dtype=bool)
    dbl_thresh = np.nan
    dbl_status = "not run"
    call_method = "none"

    # 10x multiplet rate is roughly 0.8% per 1000 cells recovered; that is the
    # prior we hold Scrublet's automatic threshold against.
    expected_rate = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    result["expected_doublet_rate"] = round(expected_rate, 4)

    if adata.n_obs >= 100:
        try:
            sc.pp.scrublet(adata, random_state=RANDOM_SEED, verbose=False,
                           expected_doublet_rate=expected_rate)
            dbl_score = adata.obs["doublet_score"].to_numpy().astype(np.float32)
            dbl_call = adata.obs["predicted_doublet"].to_numpy().astype(bool)
            dbl_thresh = float(
                adata.uns.get("scrublet", {}).get("threshold", np.nan))
            auto_rate = float(dbl_call.mean())
            dbl_status = "ok (automatic threshold)"
            call_method = "scrublet automatic threshold"
            # The automatic threshold is only trustworthy when the simulated
            # doublet score histogram is bimodal.  When it is not, it can call
            # essentially nothing (or almost everything); fall back to the
            # expected-rate quantile and say so rather than silently accepting
            # a doublet rate that cannot be right.
            if auto_rate < 0.2 * expected_rate or auto_rate > 3 * expected_rate:
                dbl_thresh = float(np.quantile(dbl_score, 1 - expected_rate))
                dbl_call = dbl_score >= dbl_thresh
                dbl_status = (f"ok (automatic threshold rejected: it called "
                              f"{100 * auto_rate:.2f}% vs an expected "
                              f"{100 * expected_rate:.2f}%)")
                call_method = (f"top {100 * expected_rate:.2f}% of Scrublet "
                               f"scores (10x multiplet-rate prior)")
        except Exception as exc:  # noqa: BLE001
            dbl_status = f"failed: {type(exc).__name__}: {exc}"
    else:
        dbl_status = "skipped (fewer than 100 cells)"
    result["doublet_call_method"] = call_method

    adata.obs["doublet_score"] = dbl_score
    adata.obs["predicted_doublet"] = dbl_call
    result["doublet_status"] = dbl_status
    result["doublet_threshold"] = dbl_thresh
    result["n_predicted_doublets"] = int(np.nansum(dbl_call))
    result["doublet_rate"] = (float(np.nansum(dbl_call)) / adata.n_obs
                              if adata.n_obs else 0.0)

    dmap = dict(zip(adata.obs_names.to_numpy(), dbl_score))
    cmap = dict(zip(adata.obs_names.to_numpy(), dbl_call))
    pre["doublet_score"] = pre["barcode"].map(dmap)
    pre["predicted_doublet"] = pre["barcode"].map(cmap)

    # ---- identity is fixed here and never re-derived downstream ----------
    adata.obs["sample_id"] = sample
    adata.obs["gsm"] = task.get("gsm", "")
    adata.obs["original_barcode"] = adata.obs_names.to_numpy()
    adata.obs_names = [f"{sample}_{bc}" for bc in adata.obs_names]

    # keep only what the merged object needs; var is rebuilt at merge time
    adata.var = adata.var[["gene_symbol", "gene_id", "mt", "ribo", "hb"]]

    shard = shard_dir / f"{sample}.h5ad"
    adata.write_h5ad(shard, compression="gzip")
    result["shard"] = str(shard)
    result["n_cells_written"] = int(adata.n_obs)
    result["nnz"] = int(adata.X.nnz)

    pre_path = shard_dir / f"{sample}__preQC_obs.parquet"
    try:
        pre.to_parquet(pre_path, index=False)
    except Exception:  # noqa: BLE001 - pyarrow may be absent
        pre_path = shard_dir / f"{sample}__preQC_obs.csv.gz"
        pre.to_csv(pre_path, index=False, compression="gzip")
    result["preqc_obs"] = str(pre_path)

    return result


def _worker_entry(task: dict) -> str:
    """Top-level callable for ProcessPoolExecutor; returns JSON."""
    try:
        return json.dumps(process_sample(task), default=str)
    except Exception as exc:  # noqa: BLE001
        import traceback
        return json.dumps({
            "sample": task.get("sample"), "error": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc()}, default=str)
