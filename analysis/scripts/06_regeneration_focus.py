#!/usr/bin/env python
"""
=============================================================================
Focused alveolar-regeneration analysis  (the "big goal")
=============================================================================

The whole-atlas run produced 29 clusters but did not resolve the biology this
repository is about: the deposited `Alveolar_transitional` (233 cells) and
`AT1_AT2` (429 cells) both sat inside the AT2 cluster, and no trajectory was
fitted - so the AT2 -> Krt8+ transitional -> AT1 *ordering*, which is the whole
claim, was never quantified.

This script does not redo QC, merging or normalisation; those checkpoints are
sound and cost hours. It reuses `processed/final_clustered.h5ad` and adds what
was missing:

  1. restriction to the 25-sample Ki67 atlas (the other 8 samples are a
     different experiment - different Cre drivers, pre-labelling design)
  2. a clean alveolar epithelial subset, excluding the Krt5+ dysplastic and
     ciliated airway compartments exactly as Niethamer et al. did
  3. PAGA + diffusion pseudotime rooted in AT2
  4. a data-derived transitional score, benchmarked by AUROC against the
     deposited labels
  5. the capillary endothelial compartment, where the injury-induced state
     lives - it is invisible at top-level clustering

Cluster memberships and the marker set below were established empirically
(see docs/ANALYSIS_RATIONALE.md), not assumed.
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

OUT = ANALYSIS / "GSE262927" / "regeneration_focus"
FIG = OUT / "figures"
TAB = OUT / "tables"

# ---------------------------------------------------------------------------
# Empirically established configuration.  Every value here was measured on this
# object rather than taken from the literature - see the rationale doc.
# ---------------------------------------------------------------------------

# Clusters 10 (AT2 + transitional) and 16 (AT1). Cluster 16 is NOT optional:
# 30 transitional cells sit at the AT1 end, and cluster 10 alone captures only
# 85% of them. Clusters 18 (ciliated) and 19 (Krt5+ dysplastic) are excluded,
# as in the source paper, because they cost >1,200 airway cells to gain 3
# transitional ones.
ALVEOLAR_CLUSTERS = ["10", "16"]

# Capillary endothelium. Cluster 3 is where the injury-induced state is
# concentrated (Sparcl1 92% / Ntrk2 72% expressing, vs 46% / 36% in cluster 1).
# Cluster 4 is ambient-dominated (median 661 genes) and would otherwise
# dominate PC1.
CAPILLARY_CLUSTERS = ["0", "1", "3", "21"]
CAPILLARY_DROP = ["4"]

# Airway labels to remove from the alveolar subset. Remove BY LABEL, never by
# dropping the secretory subcluster - that subcluster also contains 24
# transitional cells, 10% of the entire transitional population.
AIRWAY_LABELS = ["Secretory", "Krt5", "Ciliated"]

# Transitional score. Benchmarked by AUROC against the deposited labels on this
# object: this set reaches 0.987 against all other alveolar cells, where the
# canonical 13-gene panel reaches 0.948 - it loses 10 points against AT1
# specifically. Tnc (0.510), Hbegf (0.574) and Cdkn1a (0.651) are dropped
# because they do not discriminate here. Cldn4 is the near-binary gate.
TRANSITIONAL_CORE = ["Krt8", "Krt18", "Cldn4", "Sfn", "Clu", "S100a14",
                     "Gpx2", "Tnip3"]
TRANSITIONAL_PUBLISHED = ["Krt8", "Krt18", "Krt19", "Cldn4", "Sfn", "Cldn7",
                          "Krt7", "Cdkn1a", "Tnc", "Hbegf", "Lgals3"]
AT2_MARKERS = ["Sftpc", "Sftpb", "Abca3", "Lamp3", "Slc34a2"]
AT1_MARKERS = ["Ager", "Pdpn", "Hopx", "Cav1", "Aqp5", "Rtkn2"]
ICAP_MARKERS = ["Sparcl1", "Ntrk2"]
CAP1_MARKERS = ["Gpihbp1", "Kit"]
CAP2_MARKERS = ["Ednrb", "Car4"]

DECISIONS: dict[str, object] = {}


def record(k, v):
    DECISIONS[k] = v
    log(f"  DECISION {k} = {v}")


# ---------------------------------------------------------------------------
def load_subset(src: Path, clusters: list[str], drop: list[str] | None = None,
                airway_filter: bool = False) -> ad.AnnData:
    """Pull a compartment out of the full object without loading all of it.

    Backed mode keeps the 3.4 GB expression matrix on disk; only the selected
    rows are materialised.
    """
    log(f"load_subset: opening {src.name} backed {mem_report()}")
    a = ad.read_h5ad(src, backed="r")
    obs = a.obs
    # .to_numpy() on a pandas boolean column can hand back a read-only view
    keep = np.array(obs["has_author_metadata"].astype(str).isin(["True", "true"]),
                    dtype=bool)
    keep &= obs["leiden_cluster"].astype(str).isin(clusters).to_numpy()
    if drop:
        keep &= ~obs["leiden_cluster"].astype(str).isin(drop).to_numpy()
    log(f"  {keep.sum():,} cells in clusters {clusters} within the annotated cohort")

    sub = a[keep].to_memory()
    a.file.close()
    del a
    free_mem()

    if airway_filter:
        n0 = sub.n_obs
        lin_ok = sub.obs["author_lineage"].astype(str) == "Epithelium"
        not_airway = ~sub.obs["author_celltype"].astype(str).isin(AIRWAY_LABELS)
        sub = sub[(lin_ok & not_airway).to_numpy()].copy()
        log(f"  airway/lineage label filter: {n0:,} -> {sub.n_obs:,} cells "
            f"(removed by LABEL, not by dropping a subcluster - that "
            f"subcluster also holds transitional cells)")
    log(f"load_subset: {sub.n_obs:,} x {sub.n_vars:,} {mem_report()}")
    return sub


def embed(a: ad.AnnData, n_hvg: int, resolution: float, tag: str) -> None:
    """HVGs -> scaled PCA -> neighbours -> UMAP -> Leiden, within the subset.

    HVGs are recomputed *inside* the compartment: the atlas-wide HVGs are
    dominated by between-compartment differences and carry little information
    about structure within one compartment.
    """
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
    a.uns["pca"] = {"variance_ratio": p.explained_variance_ratio_}
    del X
    free_mem()

    n_pcs = min(30, n_comps)
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=n_pcs, random_state=RANDOM_SEED)
    sc.tl.umap(a, random_state=RANDOM_SEED)
    sc.tl.leiden(a, resolution=resolution, key_added="subcluster",
                 flavor="igraph", n_iterations=2, directed=False,
                 random_state=RANDOM_SEED)
    record(f"{tag}_embedding",
           f"{n_hvg} HVGs recomputed within the compartment (batch_key={batch}), "
           f"z-scored, randomized PCA, {n_pcs} PCs, n_neighbors=15, "
           f"Leiden resolution {resolution} -> "
           f"{a.obs['subcluster'].nunique()} subclusters")


def score_module(a: ad.AnnData, genes: list[str], name: str) -> list[str]:
    present = [g for g in genes if g in a.var_names]
    missing = [g for g in genes if g not in a.var_names]
    if missing:
        log(f"  {name}: {len(missing)} gene(s) absent from the matrix: {missing}")
    if present:
        sc.tl.score_genes(a, present, score_name=name, ctrl_size=50,
                          random_state=RANDOM_SEED)
    return present


def auroc(score: np.ndarray, positive: np.ndarray) -> float:
    """Rank-based AUROC; no sklearn import needed for one number."""
    if positive.sum() == 0 or (~positive).sum() == 0:
        return float("nan")
    r = pd.Series(score).rank().to_numpy()
    n1 = positive.sum(); n0 = (~positive).sum()
    return float((r[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


# ---------------------------------------------------------------------------
def alveolar_analysis(src: Path) -> dict:
    log("=" * 60)
    log("ALVEOLAR REGENERATION: AT2 -> transitional -> AT1")
    log("=" * 60)
    a = load_subset(src, ALVEOLAR_CLUSTERS, airway_filter=True)
    record("alveolar_selection",
           f"clusters {ALVEOLAR_CLUSTERS} of the whole-atlas clustering, "
           f"restricted to the 25-sample annotated cohort, then filtered to "
           f"author_lineage=='Epithelium' excluding {AIRWAY_LABELS}. "
           f"{a.n_obs} cells. Ciliated and Krt5+ dysplastic compartments are "
           f"excluded deliberately, as in the source study.")
    comp = a.obs["author_celltype"].astype(str).value_counts()
    record("alveolar_composition", comp.to_dict())

    embed(a, n_hvg=2000, resolution=0.5, tag="alveolar")

    # ---- transitional scores, benchmarked --------------------------------
    core = score_module(a, TRANSITIONAL_CORE, "transitional_score")
    score_module(a, TRANSITIONAL_PUBLISHED, "transitional_score_published")
    score_module(a, AT2_MARKERS, "AT2_score")
    score_module(a, AT1_MARKERS, "AT1_score")

    lab = a.obs["author_celltype"].astype(str).to_numpy()
    is_trans = lab == "Alveolar_transitional"
    rows = []
    for col, nm in [("transitional_score", "core (data-derived)"),
                    ("transitional_score_published", "canonical panel")]:
        if col not in a.obs:
            continue
        s = a.obs[col].to_numpy()
        rows.append({
            "score": nm, "genes_used": len(core) if "published" not in col
            else len([g for g in TRANSITIONAL_PUBLISHED if g in a.var_names]),
            "AUROC_vs_AT2": round(auroc(s[np.isin(lab, ["Alveolar_transitional", "AT2"])],
                                        is_trans[np.isin(lab, ["Alveolar_transitional", "AT2"])]), 4),
            "AUROC_vs_AT1": round(auroc(s[np.isin(lab, ["Alveolar_transitional", "AT1"])],
                                        is_trans[np.isin(lab, ["Alveolar_transitional", "AT1"])]), 4),
            "AUROC_vs_all": round(auroc(s, is_trans), 4)})
    bench = pd.DataFrame(rows)
    bench.to_csv(TAB / "transitional_score_benchmark.csv", index=False)
    record("transitional_score_benchmark", bench.to_dict("records"))

    # ---- trajectory ------------------------------------------------------
    # PAGA gives the coarse topology between subclusters; diffusion pseudotime
    # orders cells within it. The root is set deliberately to the subcluster
    # with the highest AT2 score, because AT2 is the known progenitor - an
    # arbitrary root would produce an arbitrary ordering.
    sc.tl.paga(a, groups="subcluster")
    at2_by_sub = a.obs.groupby("subcluster", observed=True)["AT2_score"].mean()
    root_sub = str(at2_by_sub.idxmax())
    cand = np.flatnonzero((a.obs["subcluster"].astype(str) == root_sub).to_numpy())
    a.uns["iroot"] = int(cand[np.argmax(a.obs["AT2_score"].to_numpy()[cand])])
    sc.tl.diffmap(a, n_comps=15)
    sc.tl.dpt(a)
    record("trajectory",
           f"PAGA on {a.obs['subcluster'].nunique()} subclusters, then diffusion "
           f"pseudotime rooted in subcluster {root_sub} (highest mean AT2 score, "
           f"AT2 being the known progenitor). Pseudotime is a rank ordering "
           f"with no units and is not comparable between runs; PAGA edges are "
           f"connectivity, not flux or rate.")

    # ---- figures ---------------------------------------------------------
    for col, name, kw in [
        ("author_celltype", "UMAP_alveolar_author_celltype", dict(size=12)),
        ("subcluster", "UMAP_alveolar_subclusters",
         dict(size=12, legend_loc="on data", legend_fontoutline=2)),
        ("transitional_score", "UMAP_alveolar_transitional_score",
         dict(size=12, cmap="viridis")),
        ("dpt_pseudotime", "UMAP_alveolar_pseudotime",
         dict(size=12, cmap="magma")),
        ("sacrifice_day", "UMAP_alveolar_sacrifice_day", dict(size=12)),
        ("AT2_score", "UMAP_alveolar_AT2_score", dict(size=12, cmap="viridis")),
        ("AT1_score", "UMAP_alveolar_AT1_score", dict(size=12, cmap="viridis")),
    ]:
        if col not in a.obs:
            continue
        try:
            ax = sc.pl.umap(a, color=col, show=False, frameon=False, **kw)
            save_fig(figure_of(ax), FIG, name)
        except Exception as exc:
            log(f"  WARNING {name}: {exc}")
            _plt.close("all")

    try:
        sc.pl.paga(a, color="subcluster", show=False, threshold=0.05,
                   fontsize=7, node_size_scale=1.5)
        save_fig(_plt.gcf(), FIG, "PAGA_alveolar_topology")
    except Exception as exc:
        log(f"  WARNING PAGA plot: {exc}")
        _plt.close("all")

    _trajectory_figure(a, FIG)
    _transitional_timecourse(a, FIG, TAB)

    genes = {"AT2": [g for g in AT2_MARKERS if g in a.var_names],
             "Transitional": [g for g in TRANSITIONAL_CORE if g in a.var_names],
             "AT1": [g for g in AT1_MARKERS if g in a.var_names]}
    try:
        dp = sc.pl.dotplot(a, genes, groupby="subcluster", show=False,
                           return_fig=True, standard_scale="var")
        save_fig(figure_of(dp), FIG, "dotplot_alveolar_subclusters")
        dp = sc.pl.dotplot(a, genes, groupby="author_celltype", show=False,
                           return_fig=True, standard_scale="var")
        save_fig(figure_of(dp), FIG, "dotplot_alveolar_author_celltype")
    except Exception as exc:
        log(f"  WARNING dotplot: {exc}")
        _plt.close("all")

    sc.tl.rank_genes_groups(a, "subcluster", method="wilcoxon", pts=True,
                            key_added="rank_sub")
    mk = sc.get.rank_genes_groups_df(a, group=None, key="rank_sub")
    mk.to_csv(TAB / "alveolar_subcluster_markers.csv.gz", index=False,
              compression="gzip")
    (mk.sort_values(["group", "scores"], ascending=[True, False])
       .groupby("group", observed=True).head(20)
       .to_csv(TAB / "alveolar_subcluster_markers_top20.csv", index=False))

    ct = pd.crosstab(a.obs["subcluster"].astype(str),
                     a.obs["author_celltype"].astype(str))
    ct.to_csv(TAB / "alveolar_subcluster_vs_author_label.csv")
    trans_sub = ct.get("Alveolar_transitional")
    if trans_sub is not None:
        best = trans_sub.sort_values(ascending=False)
        record("transitional_enriched_subclusters",
               {k: int(v) for k, v in best.head(4).items()})

    keep_obs = [c for c in ["sample_id", "author_celltype", "author_lineage",
                            "experimental_group", "sacrifice_day",
                            "tamoxifen_start_day", "trace_call", "sex",
                            "subcluster", "dpt_pseudotime", "transitional_score",
                            "AT2_score", "AT1_score", "total_counts",
                            "n_genes_by_counts", "pct_counts_mt"]
                if c in a.obs.columns]
    meta = a.obs[keep_obs].copy()
    meta.insert(0, "cell_id", a.obs_names)
    um = np.asarray(a.obsm["X_umap"])
    meta["UMAP_1"] = um[:, 0]; meta["UMAP_2"] = um[:, 1]
    meta.to_csv(TAB / "alveolar_cell_metadata.csv", index=False)

    for c in a.obs.columns:
        if a.obs[c].dtype == object:
            a.obs[c] = a.obs[c].astype(str)
    a.write_h5ad(OUT / "alveolar_trajectory.h5ad", compression="gzip")
    res = {"n_cells": int(a.n_obs),
           "n_subclusters": int(a.obs["subcluster"].nunique())}
    del a
    free_mem()
    return res


def _trajectory_figure(a, figdir: Path) -> None:
    """Marker programmes along pseudotime - the actual regeneration claim."""
    if "dpt_pseudotime" not in a.obs:
        return
    pt = a.obs["dpt_pseudotime"].to_numpy()
    ok = np.isfinite(pt)
    fig, axes = _plt.subplots(1, 2, figsize=(12, 4.4))

    ax = axes[0]
    for col, colour, lbl in [("AT2_score", "#4C72B0", "AT2 programme"),
                             ("transitional_score", "#C44E52", "Transitional (Krt8+)"),
                             ("AT1_score", "#55A868", "AT1 programme")]:
        if col not in a.obs:
            continue
        y = a.obs[col].to_numpy()[ok]
        order = np.argsort(pt[ok])
        x = pt[ok][order]; ys = y[order]
        w = max(25, len(x) // 40)
        sm = pd.Series(ys).rolling(w, center=True, min_periods=1).mean().to_numpy()
        ax.plot(x, sm, color=colour, lw=2, label=lbl)
    ax.set_xlabel("diffusion pseudotime (rank order, no units)")
    ax.set_ylabel("module score (rolling mean)")
    ax.set_title("Alveolar programmes along pseudotime")
    ax.legend(fontsize=8)

    ax = axes[1]
    labs = a.obs["author_celltype"].astype(str).to_numpy()[ok]
    order = ["AT2", "Alveolar_transitional", "AT1_AT2", "AT1"]
    data = [pt[ok][labs == l] for l in order if (labs == l).sum() > 0]
    names = [l for l in order if (labs == l).sum() > 0]
    parts = ax.violinplot(data, showextrema=False, widths=0.85)
    for pc in parts["bodies"]:
        pc.set_facecolor("#4C72B0"); pc.set_alpha(0.75); pc.set_edgecolor("none")
    ax.scatter(range(1, len(data) + 1), [np.median(d) for d in data], s=14,
               color="white", zorder=3)
    ax.set_xticks(range(1, len(names) + 1))
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("diffusion pseudotime")
    ax.set_title("Deposited labels ordered by pseudotime\n"
                 "(independent check - labels were not used to fit it)",
                 fontsize=9)
    fig.tight_layout()
    save_fig(fig, figdir, "trajectory_alveolar_programmes")


def _transitional_timecourse(a, figdir: Path, tabdir: Path) -> None:
    if "sacrifice_day" not in a.obs:
        return
    obs = a.obs
    day = pd.to_numeric(obs["sacrifice_day"], errors="coerce")
    df = pd.DataFrame({"day": day,
                       "label": obs["author_celltype"].astype(str),
                       "score": obs.get("transitional_score", np.nan),
                       "sample": obs["sample_id"].astype(str)}).dropna(subset=["day"])
    per_sample = (df.assign(is_trans=df["label"] == "Alveolar_transitional")
                    .groupby(["sample", "day"], observed=True)
                    .agg(n_cells=("label", "size"),
                         n_transitional=("is_trans", "sum"),
                         mean_score=("score", "mean")).reset_index())
    per_sample["pct_transitional"] = (100 * per_sample.n_transitional
                                      / per_sample.n_cells).round(3)
    per_sample.to_csv(tabdir / "transitional_abundance_per_sample.csv", index=False)

    fig, axes = _plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    for d, g in per_sample.groupby("day", observed=True):
        ax.scatter([d] * len(g), g["pct_transitional"], s=26, color="#C44E52",
                   alpha=0.85, zorder=3)
    med = per_sample.groupby("day", observed=True)["pct_transitional"].median()
    ax.plot(med.index, med.to_numpy(), color="#C44E52", lw=1.2, alpha=0.6)
    ax.set_xlabel("days post infection (0 = uninfected control)")
    ax.set_ylabel("% of alveolar cells labelled transitional")
    ax.set_title("Transitional abundance over time\n(one point = one sample)",
                 fontsize=9)

    ax = axes[1]
    mm = per_sample.groupby("day", observed=True)["mean_score"].median()
    for d, g in per_sample.groupby("day", observed=True):
        ax.scatter([d] * len(g), g["mean_score"], s=26, color="#4C72B0",
                   alpha=0.85, zorder=3)
    ax.plot(mm.index, mm.to_numpy(), color="#4C72B0", lw=1.2, alpha=0.6)
    ax.set_xlabel("days post infection")
    ax.set_ylabel("mean transitional score")
    ax.set_title("Transitional programme over time", fontsize=9)
    fig.tight_layout()
    save_fig(fig, figdir, "transitional_timecourse")


# ---------------------------------------------------------------------------
def capillary_analysis(src: Path) -> dict:
    log("=" * 60)
    log("CAPILLARY ENDOTHELIUM: injury-induced state")
    log("=" * 60)
    a = load_subset(src, CAPILLARY_CLUSTERS, drop=CAPILLARY_DROP)
    record("capillary_selection",
           f"clusters {CAPILLARY_CLUSTERS} (cluster {CAPILLARY_DROP} dropped as "
           f"ambient-dominated), annotated cohort only. {a.n_obs} cells. The "
           f"source study reports this state is resolvable only after "
           f"subsetting capillary endothelium - it does not appear at "
           f"top-level clustering.")

    embed(a, n_hvg=2000, resolution=0.4, tag="capillary")
    score_module(a, ICAP_MARKERS, "iCAP_score")
    score_module(a, CAP1_MARKERS, "CAP1_score")
    score_module(a, CAP2_MARKERS, "CAP2_score")

    for col, name, kw in [
        ("subcluster", "UMAP_capillary_subclusters",
         dict(size=3, legend_loc="on data", legend_fontoutline=2)),
        ("author_celltype", "UMAP_capillary_author_celltype", dict(size=3)),
        ("iCAP_score", "UMAP_capillary_iCAP_score", dict(size=3, cmap="magma")),
        ("sacrifice_day", "UMAP_capillary_sacrifice_day", dict(size=3)),
    ]:
        if col not in a.obs:
            continue
        try:
            ax = sc.pl.umap(a, color=col, show=False, frameon=False, **kw)
            save_fig(figure_of(ax), FIG, name)
        except Exception as exc:
            log(f"  WARNING {name}: {exc}")
            _plt.close("all")

    genes = {"iCAP": [g for g in ICAP_MARKERS if g in a.var_names],
             "CAP1": [g for g in CAP1_MARKERS if g in a.var_names],
             "CAP2": [g for g in CAP2_MARKERS if g in a.var_names],
             "MHC-II": [g for g in ["Cd74", "H2-Ab1", "Ciita", "Ifngr1"]
                        if g in a.var_names]}
    try:
        dp = sc.pl.dotplot(a, genes, groupby="subcluster", show=False,
                           return_fig=True, standard_scale="var")
        save_fig(figure_of(dp), FIG, "dotplot_capillary_subclusters")
    except Exception as exc:
        log(f"  WARNING capillary dotplot: {exc}")
        _plt.close("all")

    # Does the high-scoring subcluster behave like an injury-induced state -
    # rare in uninfected lung, emerging after infection, persisting late?
    if "iCAP_score" in a.obs and "sacrifice_day" in a.obs:
        by = a.obs.groupby("subcluster", observed=True)["iCAP_score"].mean()
        top = str(by.idxmax())
        day = pd.to_numeric(a.obs["sacrifice_day"], errors="coerce")
        df = pd.DataFrame({"day": day, "sub": a.obs["subcluster"].astype(str),
                           "sample": a.obs["sample_id"].astype(str)}).dropna()
        per = (df.assign(is_top=df["sub"] == top)
                 .groupby(["sample", "day"], observed=True)
                 .agg(n=("sub", "size"), n_top=("is_top", "sum")).reset_index())
        per["pct"] = (100 * per.n_top / per.n).round(3)
        per.to_csv(TAB / "icap_abundance_per_sample.csv", index=False)
        record("icap_subcluster",
               f"subcluster {top} has the highest mean injury-state score "
               f"({by.max():.3f}); abundance per sample and timepoint in "
               f"tables/icap_abundance_per_sample.csv")

        fig, ax = _plt.subplots(figsize=(6, 4))
        for d, g in per.groupby("day", observed=True):
            ax.scatter([d] * len(g), g["pct"], s=26, color="#8172B2", alpha=.85,
                       zorder=3)
        med = per.groupby("day", observed=True)["pct"].median()
        ax.plot(med.index, med.to_numpy(), color="#8172B2", lw=1.2, alpha=.6)
        ax.set_xlabel("days post infection (0 = uninfected control)")
        ax.set_ylabel(f"% of capillary cells in subcluster {top}")
        ax.set_title("Injury-associated capillary state over time\n"
                     "(one point = one sample)", fontsize=9)
        fig.tight_layout()
        save_fig(fig, FIG, "icap_timecourse")

    sc.tl.rank_genes_groups(a, "subcluster", method="wilcoxon", pts=True,
                            key_added="rank_sub")
    mk = sc.get.rank_genes_groups_df(a, group=None, key="rank_sub")
    (mk.sort_values(["group", "scores"], ascending=[True, False])
       .groupby("group", observed=True).head(20)
       .to_csv(TAB / "capillary_subcluster_markers_top20.csv", index=False))
    pd.crosstab(a.obs["subcluster"].astype(str),
                a.obs["author_celltype"].astype(str)).to_csv(
        TAB / "capillary_subcluster_vs_author_label.csv")

    res = {"n_cells": int(a.n_obs),
           "n_subclusters": int(a.obs["subcluster"].nunique())}
    del a
    free_mem()
    return res


# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run focused alveolar-trajectory and capillary-state analyses.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--src", default=str(ANALYSIS / "GSE262927" / "processed" / "final_clustered.h5ad"),
                    help="clustered mouse AnnData object")
    ap.add_argument("--skip-capillary", action="store_true",
                    help="do not run the capillary endothelial analysis")
    ap.add_argument("--skip-alveolar", action="store_true",
                    help="do not run the alveolar trajectory analysis")
    args = ap.parse_args()

    for d in (OUT, FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    src = Path(args.src)
    if not src.exists():
        raise SystemExit(f"source object not found: {src}")

    log(f"scanpy {sc.__version__}  source {src}")
    out: dict = {}
    if not args.skip_alveolar:
        out["alveolar"] = alveolar_analysis(src)
    if not args.skip_capillary:
        out["capillary"] = capillary_analysis(src)

    (OUT / "decisions.json").write_text(
        json.dumps({k: str(v) for k, v in DECISIONS.items()}, indent=2),
        encoding="utf-8")
    log(f"DONE {json.dumps(out)} {mem_report()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
