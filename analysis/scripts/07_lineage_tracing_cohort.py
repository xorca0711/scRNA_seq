#!/usr/bin/env python
"""
=============================================================================
The 8 non-atlas samples: lineage tracing of the injury-induced capillary state
=============================================================================

`GSE262927_CellMetaData.csv` annotates 25 of the 33 deposited samples. The
other 8 are not an unannotated remainder of the atlas - they are a separate
experiment, and the source study analysed them separately:

    EEM-scRNA-249/250/251   Kit-MerCreMer      labels CAP1 capillary endothelium
    EEM-scRNA-288/289       Car4-CreERT2       labels CAP2 (and mature aMAC)
    EEM-scRNA-290/291/292   Ednrb-CreERT2      labels CAP2

All were tamoxifen-labelled *before* infection (the opposite logic to the
Ki67 proliferation-tracing atlas) and all were harvested at 19 dpi. The design
question they answer is one the atlas cannot: **where does the injury-induced
capillary state come from?** If it arises from CAP1 it should be Kit-traced; if
from CAP2, Car4/Ednrb-traced; if from both, traced in all three lines.

Tracing is recomputed here from the reporter contigs rather than taken on
trust. The rule - `Traced` iff SiteB/(SiteA+SiteB) > 0.5, `Not_detected` iff
neither is present - reproduces the authors' own `trace_call` on the annotated
cohort with **100.0000% agreement over 107,626 cells, zero disagreements**,
which is what licenses applying it to samples that carry no metadata.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_utils import (ANALYSIS, RANDOM_SEED, figure_of, free_mem, log,  # noqa: E402
                            mem_report, save_fig, setup_matplotlib)

warnings.filterwarnings("ignore")
plt = setup_matplotlib()
import matplotlib.pyplot as _plt  # noqa: E402
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402

sc.settings.verbosity = 1
np.random.seed(RANDOM_SEED)

OUT = ANALYSIS / "GSE262927" / "lineage_tracing_cohort"
FIG = OUT / "figures"
TAB = OUT / "tables"

# Cre driver per sample, from the study's supplementary sample table.
CRE_LINE = {
    "EEM-scRNA-249": "Kit-MerCreMer", "EEM-scRNA-250": "Kit-MerCreMer",
    "EEM-scRNA-251": "Kit-MerCreMer",
    "EEM-scRNA-288": "Car4-CreERT2", "EEM-scRNA-289": "Car4-CreERT2",
    "EEM-scRNA-290": "Ednrb-CreERT2", "EEM-scRNA-291": "Ednrb-CreERT2",
    "EEM-scRNA-292": "Ednrb-CreERT2",
}
# What each line labels at homeostasis, i.e. the candidate cell of origin.
LINE_LABELS = {"Kit-MerCreMer": "CAP1", "Car4-CreERT2": "CAP2",
               "Ednrb-CreERT2": "CAP2"}

ICAP = ["Sparcl1", "Ntrk2"]
CAP1 = ["Gpihbp1", "Kit"]
CAP2 = ["Ednrb", "Car4"]
PAN_EC = ["Pecam1", "Cdh5", "Cldn5", "Emcn"]
NON_EC = ["Ptprc", "Epcam", "Col1a1"]

DECISIONS: dict[str, object] = {}


def record(k, v):
    DECISIONS[k] = v
    log(f"  DECISION {k} = {v}")


def trace_call(site_a: np.ndarray, site_b: np.ndarray) -> np.ndarray:
    """The study's own definition, validated at 100% on the annotated cohort."""
    tot = site_a + site_b
    return np.where(tot == 0, "Not_detected",
                    np.where(site_b / np.maximum(tot, 1e-9) > 0.5,
                             "Traced", "Untraced"))


def validate_trace_rule(src: Path) -> dict:
    """Re-derive the authors' trace_call and confirm the rule before using it."""
    meta = ANALYSIS / "GSE262927" / "tables" / "cell_metadata.csv"
    if not meta.exists():
        return {"validated": False, "reason": "cell_metadata.csv missing"}
    m = pd.read_csv(meta, low_memory=False)
    k = m[m["has_author_metadata"].astype(str).isin(["True", "true"])]
    k = k[k["trace_call"].astype(str).isin(["Traced", "Untraced", "Not_detected"])]
    if not len(k):
        return {"validated": False, "reason": "no annotated trace_call"}
    pred = trace_call(k["reporter_SiteA"].to_numpy(),
                      k["reporter_SiteB"].to_numpy())
    obs = k["trace_call"].astype(str).to_numpy()
    agree = float((pred == obs).mean())
    out = {"validated": bool(agree > 0.999), "n_cells": int(len(k)),
           "agreement": round(agree, 6)}
    record("trace_rule_validation",
           f"the rule 'Traced iff SiteB/(SiteA+SiteB) > 0.5, Not_detected iff "
           f"neither present' reproduces the deposited trace_call on "
           f"{len(k):,} annotated cells with {100 * agree:.4f}% agreement - "
           f"which is what licenses applying it to the 8 samples that carry no "
           f"metadata")
    if not out["validated"]:
        raise SystemExit("trace-call rule failed validation; refusing to apply "
                         "it to unannotated samples")
    return out


def load_cohort(src: Path) -> ad.AnnData:
    log(f"load_cohort: opening {src.name} backed {mem_report()}")
    a = ad.read_h5ad(src, backed="r")
    keep = np.array(a.obs["sample_id"].astype(str).isin(CRE_LINE), dtype=bool)
    log(f"  {keep.sum():,} cells across {len(CRE_LINE)} lineage-tracing samples")
    sub = a[keep].to_memory()
    a.file.close(); del a; free_mem()

    sub.obs["cre_line"] = pd.Categorical(
        sub.obs["sample_id"].astype(str).map(CRE_LINE))
    sub.obs["labels_at_homeostasis"] = pd.Categorical(
        sub.obs["cre_line"].astype(str).map(LINE_LABELS))
    sub.obs["trace_call"] = pd.Categorical(trace_call(
        sub.obs["reporter_SiteA"].to_numpy(),
        sub.obs["reporter_SiteB"].to_numpy()))
    record("cohort",
           f"{sub.n_obs} cells from {sub.obs['sample_id'].nunique()} samples, "
           f"all harvested at 19 dpi with tamoxifen given BEFORE infection. "
           f"Cre lines: " + ", ".join(
               f"{k}={int(v)}" for k, v in
               sub.obs['cre_line'].value_counts().items()))
    record("trace_call_distribution",
           sub.obs["trace_call"].value_counts().to_dict())
    return sub


def embed(a: ad.AnnData, n_hvg: int, resolution: float) -> None:
    batch = "sample_id" if a.obs["sample_id"].nunique() > 1 else None
    sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=n_hvg,
                                batch_key=batch)
    hv = np.flatnonzero(a.var["highly_variable"].to_numpy())
    X = np.asarray(a.X[:, hv].todense() if sp.issparse(a.X) else a.X[:, hv],
                   dtype=np.float32)
    mu = X.mean(0); X -= mu
    sd = X.std(0); sd[sd == 0] = 1.0; X /= sd
    np.clip(X, -10, 10, out=X)
    from sklearn.decomposition import PCA
    n_comps = int(min(50, X.shape[1] - 1, a.n_obs - 1))
    p = PCA(n_components=n_comps, svd_solver="randomized",
            random_state=RANDOM_SEED)
    a.obsm["X_pca"] = np.ascontiguousarray(p.fit_transform(X), dtype=np.float32)
    del X; free_mem()
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=min(30, n_comps),
                    random_state=RANDOM_SEED)
    sc.tl.umap(a, random_state=RANDOM_SEED)
    sc.tl.leiden(a, resolution=resolution, key_added="subcluster",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=RANDOM_SEED)


def score(a, genes, name):
    g = [x for x in genes if x in a.var_names]
    if g:
        sc.tl.score_genes(a, g, score_name=name, ctrl_size=50,
                          random_state=RANDOM_SEED)
    return g


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(ANALYSIS / "GSE262927" / "processed" / "final_clustered.h5ad"))
    args = ap.parse_args()
    for d in (OUT, FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    src = Path(args.src)

    validate_trace_rule(src)
    a = load_cohort(src)

    # ---- restrict to endothelium, by marker evidence (no author labels here)
    embed(a, n_hvg=2000, resolution=0.4)
    for nm, gs in [("panEC_score", PAN_EC), ("nonEC_score", NON_EC),
                   ("iCAP_score", ICAP), ("CAP1_score", CAP1),
                   ("CAP2_score", CAP2)]:
        score(a, gs, nm)

    per = a.obs.groupby("subcluster", observed=True)[
        ["panEC_score", "nonEC_score"]].mean()
    ec_subs = per.index[(per["panEC_score"] > 0.15)
                        & (per["panEC_score"] > per["nonEC_score"])].tolist()
    record("endothelial_subclusters",
           f"{ec_subs} selected on mean pan-endothelial score exceeding both an "
           f"absolute floor and the non-endothelial score; these samples carry "
           f"no deposited labels, so selection is by marker evidence only")
    per.round(3).to_csv(TAB / "subcluster_compartment_scores.csv")

    ec = a[a.obs["subcluster"].astype(str).isin(ec_subs)].copy()
    log(f"endothelial subset: {ec.n_obs:,} cells")
    if ec.n_obs < 500:
        raise SystemExit("endothelial subset too small to analyse")

    embed(ec, n_hvg=2000, resolution=0.4)
    for nm, gs in [("iCAP_score", ICAP), ("CAP1_score", CAP1),
                   ("CAP2_score", CAP2)]:
        score(ec, gs, nm)

    by = ec.obs.groupby("subcluster", observed=True)["iCAP_score"].mean()
    icap_sub = str(by.idxmax())
    ec.obs["is_iCAP"] = (ec.obs["subcluster"].astype(str) == icap_sub)
    record("icap_subcluster",
           f"subcluster {icap_sub} has the highest mean injury-state score "
           f"({by.max():.3f}) among {ec.n_obs} endothelial cells")

    # ---- the actual question -------------------------------------------
    # If the injury state arises from CAP1 only, it is traced in the Kit line
    # alone. If from CAP2 only, in Car4/Ednrb alone. If from both, in all three.
    # `Not_detected` is dropped: it means no reporter coverage, which is a
    # technical dropout, not an untraced cell.
    d = ec.obs[ec.obs["trace_call"].astype(str) != "Not_detected"].copy()
    d["traced"] = d["trace_call"].astype(str) == "Traced"

    rows = []
    for (line, sample), g in d.groupby(["cre_line", "sample_id"], observed=True):
        icap = g[g["is_iCAP"]]
        rest = g[~g["is_iCAP"]]
        rows.append({
            "cre_line": line, "labels_at_homeostasis": LINE_LABELS.get(str(line), "?"),
            "sample": sample,
            "n_endothelial": len(g),
            "n_iCAP": len(icap),
            "pct_traced_in_iCAP": round(100 * icap["traced"].mean(), 2) if len(icap) else np.nan,
            "pct_traced_in_other_EC": round(100 * rest["traced"].mean(), 2) if len(rest) else np.nan,
        })
    per_sample = pd.DataFrame(rows).sort_values(["cre_line", "sample"])
    per_sample.to_csv(TAB / "icap_tracing_by_cre_line.csv", index=False)
    record("icap_tracing_by_line",
           per_sample.groupby("cre_line")["pct_traced_in_iCAP"].median()
           .round(2).to_dict())

    # ---- figures ---------------------------------------------------------
    for col, name, kw in [
        ("cre_line", "UMAP_cohort_cre_line", dict(size=3)),
        ("trace_call", "UMAP_cohort_trace_call", dict(size=3)),
        ("subcluster", "UMAP_EC_subclusters",
         dict(size=4, legend_loc="on data", legend_fontoutline=2)),
        ("iCAP_score", "UMAP_EC_iCAP_score", dict(size=4, cmap="magma")),
    ]:
        obj = ec if name.startswith("UMAP_EC") else a
        if col not in obj.obs:
            continue
        try:
            ax = sc.pl.umap(obj, color=col, show=False, frameon=False, **kw)
            save_fig(figure_of(ax), FIG, name)
        except Exception as exc:
            log(f"  WARNING {name}: {exc}")
            _plt.close("all")

    fig, ax = _plt.subplots(figsize=(7, 4.4))
    lines = [l for l in ["Kit-MerCreMer", "Car4-CreERT2", "Ednrb-CreERT2"]
             if l in set(per_sample["cre_line"].astype(str))]
    colours = {"Kit-MerCreMer": "#4C72B0", "Car4-CreERT2": "#C44E52",
               "Ednrb-CreERT2": "#55A868"}
    for i, l in enumerate(lines):
        g = per_sample[per_sample["cre_line"].astype(str) == l]
        ax.scatter(np.full(len(g), i - 0.12), g["pct_traced_in_other_EC"],
                   s=44, facecolors="none", edgecolors=colours[l], linewidths=1.4)
        ax.scatter(np.full(len(g), i + 0.12), g["pct_traced_in_iCAP"],
                   s=52, color=colours[l])
    ax.set_xticks(range(len(lines)))
    ax.set_xticklabels([f"{l}\n(labels {LINE_LABELS[l]})" for l in lines],
                       fontsize=8)
    ax.set_ylabel("% of cells traced (tdTomato-recombined)")
    ax.set_title("Origin of the injury-induced capillary state at 19 dpi\n"
                 "open = other endothelium, filled = injury state; "
                 "one point per animal", fontsize=9)
    fig.tight_layout()
    save_fig(fig, FIG, "icap_origin_by_cre_line")

    sc.tl.rank_genes_groups(ec, "subcluster", method="wilcoxon", pts=True,
                            key_added="rk")
    (sc.get.rank_genes_groups_df(ec, group=None, key="rk")
       .sort_values(["group", "scores"], ascending=[True, False])
       .groupby("group", observed=True).head(20)
       .to_csv(TAB / "EC_subcluster_markers_top20.csv", index=False))

    keep = [c for c in ["sample_id", "cre_line", "labels_at_homeostasis",
                        "trace_call", "reporter_SiteA", "reporter_SiteB",
                        "subcluster", "is_iCAP", "iCAP_score", "CAP1_score",
                        "CAP2_score", "total_counts", "n_genes_by_counts"]
            if c in ec.obs.columns]
    meta = ec.obs[keep].copy()
    meta.insert(0, "cell_id", ec.obs_names)
    meta.to_csv(TAB / "endothelial_cell_metadata.csv", index=False)

    (OUT / "decisions.json").write_text(
        json.dumps({k: str(v) for k, v in DECISIONS.items()}, indent=2),
        encoding="utf-8")
    log(f"DONE  endothelial n={ec.n_obs}  iCAP subcluster={icap_sub} "
        f"{mem_report()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
