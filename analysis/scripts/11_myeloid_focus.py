#!/usr/bin/env python
"""Myeloid compartment of the Ki67 atlas: blind subclustering, per-dpi UMAP and
per-animal composition through regeneration.

Niethamer et al. 2025 (Cell Stem Cell, doi:10.1016/j.stem.2024.12.002, Figure 3)
report that alveolar macrophages (aMAC) are lost and inflammatory monocytes
(iMON) expand at 6 dpi, and that the aMAC pool is reconstituted over the
following three weeks from both resident aMACs and an iMON-derived trajectory.

The whole-atlas run put the myeloid lineage into three Leiden clusters
(5: monocyte, interstitial macrophage and dendritic cell mixture; 17: aMAC;
24: neutrophils), and trial S5 (Thesis/gate1_04_sikkema_2023_hlca/trials)
found that cluster 5 on its own resolves at Leiden 0.5, not 0.2. This script
takes the whole compartment from the 25-sample Ki67 atlas, re-embeds it with
the deposited labels held out, grades the subclusters against those labels
afterwards, and reads composition per animal and day.

The object is read row-wise (never loaded fully). Rules are frozen in the run
record before the object is opened. The statistical unit is the animal;
medians are reported and no P values are computed.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import sys  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from focus_utils import RunRecord, df_to_markdown, read_csr_rows, shannon, strip_by_day  # noqa: E402
from pipeline_utils import ANALYSIS, RANDOM_SEED, log, save_fig, setup_matplotlib  # noqa: E402

warnings.filterwarnings("ignore")
plt = setup_matplotlib()
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402

sc.settings.verbosity = 0
np.random.seed(RANDOM_SEED)

SERIES = ANALYSIS / "GSE262927"
OBJ = SERIES / "processed" / "final_clustered.h5ad"
OUT = SERIES / "myeloid_focus"
FIG = OUT / "figures"
TAB = OUT / "tables"

CLUSTERS = ["5", "17", "24"]
DAYS = [0, 6, 11, 19, 25, 42, 90, 366]
PHASE = {0: "baseline", 6: "active repair", 11: "active repair", 19: "active repair",
         25: "active repair", 42: "injury resolution", 90: "injury resolution",
         366: "long-term homeostasis"}
MYELOID_LABELS = ["aMAC", "iMAC", "iMON", "cMON", "pMON", "cDC1", "cDC2", "maDC", "Neutrophil"]
LABEL_COLOURS = {"aMAC": "#1f77b4", "iMAC": "#17becf", "iMON": "#2ca02c", "cMON": "#98df8a",
                 "pMON": "#bcbd22", "cDC1": "#d62728", "cDC2": "#ff9896", "maDC": "#e377c2",
                 "Neutrophil": "#9467bd", "other lineage": "#7f7f7f"}
PANEL = {
    "alveolar macrophage": ["Siglecf", "Ear2", "Plet1", "Fabp4", "Krt79"],
    "interstitial macrophage": ["Cd163", "Mrc1", "Folr2", "Lyve1", "C1qa", "C1qb"],
    "HLCA MDM SPP1-high": ["Spp1", "Lpl", "Chil3", "Mmp9"],
    "HLCA MDM CCL2-high": ["Ccl2", "Il1rn", "S100a8"],
    "HLCA MDM C1QA-high": ["Il18", "Trem2"],
    "HLCA MARCO-high": ["Marco", "Mcemp1"],
    "classical monocyte": ["Ly6c2", "Ccr2", "Plac8"],
    "non-classical monocyte": ["Cx3cr1", "Ace", "Nr4a1"],
    "cDC1": ["Xcr1", "Clec9a"],
    "cDC2": ["Cd209a", "Itgax", "H2-Ab1"],
    "migratory DC": ["Ccr7", "Fscn1"],
    "pDC": ["Siglech", "Bst2"],
    "neutrophil": ["S100a9", "Retnlg", "Csf3r"],
    "interferon": ["Isg15", "Ifit3", "Irf7"],
    "proliferating": ["Mki67", "Top2a"],
}
PRIMARY_RES = 0.5
SENSITIVITY_RES = [0.2, 1.0]

RULES = {
    "object": "analysis/GSE262927/processed/final_clustered.h5ad, read row-wise; never loaded fully",
    "cohort": "has_author_metadata == True (the 25-sample Ki67 atlas)",
    "cells": "atlas Leiden clusters 5, 17 and 24, chosen from the blind clustering and the cluster-level "
             "dominant-label table (tables/cluster_dominant_author_celltype.csv). Cluster 23 (ambient-like, "
             "PROGRESS item 12) is not included. Deposited labels are not used to select cells",
    "hvg": {"flavor": "seurat_v3", "n_top_genes": 2000, "layer": "counts", "batch_key": None},
    "scale_max": 10, "n_pcs": 30, "knn_k": 30,
    "batch_correction": "none (root decision 2: sample and time point are the same variable)",
    "leiden": {"primary_resolution": PRIMARY_RES, "sensitivity_resolutions": SENSITIVITY_RES,
               "flavor": "igraph", "n_iterations": 2, "random_state": 0,
               "why_primary_0.5": "trial S5 (2026-09-09) found atlas cluster 5 unresolved at 0.2 and "
                                  "resolved at 0.5; fixed here before this run"},
    "umap": {"random_state": 0},
    "grading": {"purity_threshold": 0.75, "resolved_if_fraction_of_labelled_cells_in_pure_subclusters_ge": 0.75,
                "label_entropy_high": 0.56, "labels": "deposited author_celltype, held out of the embedding"},
    "days": DAYS,
    "phases": {str(k): v for k, v in PHASE.items()},
    "umap_panels": "shared subset coordinates; every panel draws all subset cells in grey and N cells of "
                   "that day in colour, N = smallest per-day count, seed 0, shuffled draw order",
    "composition": "per animal and day: percent of each deposited myeloid label among the animal's "
                   "myeloid-labelled subset cells (descriptive, uses the answer key) and percent of each "
                   "primary subcluster among the animal's subset cells (blind). Cells whose deposited "
                   "lineage is not Myeloid are counted and excluded from the label composition",
    "paper_expectation": {"aMAC": "lost at 6 dpi, reconstituted over the following three weeks",
                          "iMON": "expanded at 6 dpi, then declining"},
    "check_rule": "median per-animal aMAC percent at 6 dpi below the 0 dpi median AND iMON percent at 6 dpi "
                  "above the 0 dpi median; reconstitution if the aMAC median at 25 or 42 dpi exceeds the "
                  "6 dpi median. Two animals per active-repair day: a description, not a test",
    "markers": "Wilcoxon one-vs-rest per primary subcluster, top 20 by score, a ranking device only",
    "comparator_panel": PANEL,
    "unit": "animal; medians; no P values",
    "seed": RANDOM_SEED,
}


def load_subset(rec: RunRecord) -> ad.AnnData:
    if not OBJ.exists():
        raise SystemExit(f"source object not found: {OBJ}")
    rec.add_input(OBJ)
    log(f"opening {OBJ.name} row-wise")
    with h5py.File(OBJ, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        n_cols = int(f["X"].attrs["shape"][1])
        in_cluster = obs["leiden_cluster"].astype(str).isin(CLUSTERS).to_numpy()
        in_cohort = obs["has_author_metadata"].astype(str).eq("True").to_numpy()
        idx = np.where(in_cluster & in_cohort)[0]
        rec.set("n_cells_in_clusters_all_samples", int(in_cluster.sum()))
        X = read_csr_rows(f["X"], idx, n_cols)
        C = read_csr_rows(f["layers"]["counts"], idx, n_cols)
    sub = ad.AnnData(X=X, obs=obs.iloc[idx].copy(), var=var.copy(), layers={"counts": C})
    sub.obs_names = obs.index[idx]
    sub.obs["day"] = pd.to_numeric(sub.obs["sacrifice_day"], errors="coerce")
    if sub.obs["day"].isna().any():
        raise SystemExit("subset cells without sacrifice_day; refusing to guess")
    sub.obs["day"] = sub.obs["day"].astype(int)
    found = set(sub.obs["day"].unique())
    if found != set(DAYS):
        raise SystemExit(f"day set {sorted(found)} differs from the frozen {DAYS}")
    sub.obs["sample_id"] = sub.obs["sample_id"].astype(str)
    lab = sub.obs["author_celltype"].astype(str)
    lin = sub.obs["author_lineage"].astype(str)
    label = np.where(lin == "Myeloid", lab, "other lineage")
    unknown = sorted(set(label[lin == "Myeloid"]) - set(MYELOID_LABELS))
    if unknown:
        raise SystemExit(f"myeloid labels outside the frozen list: {unknown}")
    cats = [k for k in MYELOID_LABELS + ["other lineage"] if k in set(label)]
    sub.obs["label"] = pd.Categorical(label, categories=cats)
    sub.uns["label_colors"] = [LABEL_COLOURS[c] for c in cats]
    rec.set("n_cells", int(sub.n_obs))
    rec.set("n_animals", int(sub.obs["sample_id"].nunique()))
    rec.set("n_cells_other_lineage_label", int((sub.obs["label"] == "other lineage").sum()))
    rec.set("cells_per_cluster", sub.obs["leiden_cluster"].astype(str).value_counts().to_dict())
    rec.set("cells_per_day", {str(k): int(v) for k, v in sub.obs["day"].value_counts().sort_index().items()})
    log(f"subset: {sub.n_obs:,} cells, {sub.obs['sample_id'].nunique()} animals; "
        f"{rec.record['results']['n_cells_other_lineage_label']} cells carry a non-myeloid label")
    return sub


def embed(sub: ad.AnnData, rec: RunRecord) -> str:
    sc.pp.highly_variable_genes(sub, flavor="seurat_v3", n_top_genes=RULES["hvg"]["n_top_genes"], layer="counts")
    emb = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=RULES["scale_max"])
    sc.tl.pca(emb, n_comps=RULES["n_pcs"], random_state=RANDOM_SEED)
    sub.obsm["X_pca"] = emb.obsm["X_pca"]
    del emb
    sc.pp.neighbors(sub, n_neighbors=RULES["knn_k"], use_rep="X_pca", random_state=RANDOM_SEED)
    keys = []
    for res in [PRIMARY_RES] + SENSITIVITY_RES:
        key = f"sub_r{res}"
        sc.tl.leiden(sub, resolution=res, key_added=key, flavor="igraph", n_iterations=2,
                     random_state=RANDOM_SEED, directed=False)
        keys.append(key)
        log(f"Leiden {res}: {sub.obs[key].nunique()} subclusters")
    sc.tl.umap(sub, random_state=RANDOM_SEED)
    rec.set("n_subclusters", {k: int(sub.obs[k].nunique()) for k in keys})
    return f"sub_r{PRIMARY_RES}"


def grade(sub: ad.AnnData, key: str) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    lab = sub.obs["label"].astype(str)
    labelled = lab != "other lineage"
    tab = pd.crosstab(sub.obs.loc[labelled, key].astype(str), lab[labelled])
    rows = []
    for scl in sorted(sub.obs[key].astype(str).unique(), key=int):
        mask = sub.obs[key].astype(str) == scl
        n_all = int(mask.sum())
        if scl in tab.index:
            counts = tab.loc[scl]
            n_lab = int(counts.sum())
            purity = float(counts.max() / n_lab)
            ent = shannon((counts / n_lab).to_numpy())
            top = str(counts.idxmax())
        else:
            n_lab, purity, ent, top = 0, np.nan, np.nan, ""
        samp = sub.obs.loc[mask, "sample_id"].astype(str).value_counts()
        rows.append({"subcluster": scl, "n_cells": n_all, "n_labelled": n_lab,
                     "top_label": top, "purity": round(purity, 3) if n_lab else None,
                     "label_entropy": round(ent, 3) if n_lab else None,
                     "entropy_flag": ("HIGH" if ent > RULES["grading"]["label_entropy_high"] else "low") if n_lab else "NA",
                     "n_animals": int(len(samp)),
                     "max_animal_fraction": round(float(samp.iloc[0] / n_all), 3),
                     "dominant_animal": samp.index[0]})
    df = pd.DataFrame(rows)
    n_lab_total = int(labelled.sum())
    in_pure = int(df.loc[df["purity"].fillna(0) >= RULES["grading"]["purity_threshold"], "n_labelled"].sum())
    summary = {"fraction_labelled_cells_in_pure_subclusters": round(in_pure / n_lab_total, 4),
               "cell_weighted_mean_entropy": round(float(np.nansum(df["label_entropy"].fillna(0) * df["n_labelled"]) / n_lab_total), 4),
               "n_subclusters": int(len(df))}
    return df, tab, summary


def composition(sub: ad.AnnData, primary: str, rec: RunRecord) -> dict:
    obs = sub.obs
    # deposited-label composition (answer key, descriptive)
    my = obs[obs["label"] != "other lineage"]
    ct = pd.crosstab([my["sample_id"], my["day"]], my["label"]).reindex(columns=MYELOID_LABELS, fill_value=0)
    rows = []
    for (s, d), r in ct.iterrows():
        n = int(r.sum())
        for lab in MYELOID_LABELS:
            rows.append({"sample_id": s, "day": int(d), "phase": PHASE[int(d)], "label": lab,
                         "n_cells": int(r[lab]), "n_myeloid_labelled_cells_animal": n,
                         "pct_of_myeloid": round(100 * r[lab] / n, 3) if n else np.nan})
    lab_long = pd.DataFrame(rows).sort_values(["day", "sample_id", "label"])
    p = TAB / "myeloid_label_composition_per_animal.csv"
    lab_long.to_csv(p, index=False)
    rec.add_output(p)
    lab_med = (lab_long.groupby(["day", "label"])
                       .agg(n_animals=("sample_id", "nunique"), median_pct_of_myeloid=("pct_of_myeloid", "median"))
                       .reset_index().round(2))
    lab_med["phase"] = lab_med["day"].map(PHASE)
    p = TAB / "myeloid_label_composition_median_by_dpi.csv"
    lab_med.to_csv(p, index=False)
    rec.add_output(p)

    # blind subcluster composition
    ct2 = pd.crosstab([obs["sample_id"], obs["day"]], obs[primary].astype(str))
    ct2 = ct2[sorted(ct2.columns, key=int)]
    rows = []
    for (s, d), r in ct2.iterrows():
        n = int(r.sum())
        for scl in ct2.columns:
            rows.append({"sample_id": s, "day": int(d), "phase": PHASE[int(d)], "subcluster": scl,
                         "n_cells": int(r[scl]), "n_subset_cells_animal": n,
                         "pct_of_subset": round(100 * r[scl] / n, 3) if n else np.nan})
    sub_long = pd.DataFrame(rows).sort_values(["day", "sample_id", "subcluster"])
    p = TAB / "myeloid_subcluster_composition_per_animal.csv"
    sub_long.to_csv(p, index=False)
    rec.add_output(p)

    # the paper's Figure 3 expectation, by the frozen rule
    med = lab_med.pivot(index="day", columns="label", values="median_pct_of_myeloid").reindex(DAYS)
    amac_lost = bool(med.loc[6, "aMAC"] < med.loc[0, "aMAC"])
    imon_up = bool(med.loc[6, "iMON"] > med.loc[0, "iMON"])
    recon = bool(max(med.loc[25, "aMAC"], med.loc[42, "aMAC"]) > med.loc[6, "aMAC"])
    check = {"median_pct_aMAC": {str(d): float(med.loc[d, "aMAC"]) for d in DAYS},
             "median_pct_iMON": {str(d): float(med.loc[d, "iMON"]) for d in DAYS},
             "median_pct_cMON": {str(d): float(med.loc[d, "cMON"]) for d in DAYS},
             "median_pct_Neutrophil": {str(d): float(med.loc[d, "Neutrophil"]) for d in DAYS},
             "aMAC_lost_at_6dpi": amac_lost, "iMON_expanded_at_6dpi": imon_up,
             "aMAC_reconstituted_by_25_or_42dpi": recon,
             "verdict": ("consistent with Figure 3" if (amac_lost and imon_up and recon)
                         else "not consistent with Figure 3")}
    rec.set("figure3_check", check)

    # figures
    fig, axes = plt.subplots(3, 3, figsize=(13, 10))
    for ax, lab in zip(axes.ravel(), MYELOID_LABELS):
        strip_by_day(ax, lab_long[lab_long["label"] == lab], "pct_of_myeloid", DAYS, LABEL_COLOURS[lab])
        ax.set_title(lab, fontsize=10)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("days post infection")
        ax.set_ylabel("% of myeloid-labelled cells")
    fig.suptitle("Within-myeloid composition per animal, deposited labels (one point = one animal; line = median)",
                 y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for p in save_fig(fig, FIG, "myeloid_label_composition_by_dpi"):
        rec.add_output(p)

    subs = list(ct2.columns)
    ncol = 4
    nrow = int(np.ceil(len(subs) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow))
    axes = np.atleast_1d(axes).ravel()
    cmap = plt.get_cmap("tab20")
    for i, scl in enumerate(subs):
        ax = axes[i]
        strip_by_day(ax, sub_long[sub_long["subcluster"] == scl], "pct_of_subset", DAYS, cmap(int(scl) % 20))
        top = rec.record["results"]["grading"][primary]["top_label_by_subcluster"].get(scl, "")
        ax.set_title(f"subcluster {scl} (mostly {top})", fontsize=9)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("days post infection")
        ax.set_ylabel("% of myeloid subset")
    for ax in axes[len(subs):]:
        ax.set_visible(False)
    fig.suptitle(f"Blind subcluster composition per animal (Leiden {PRIMARY_RES}; labels shown only in titles)",
                 y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for p in save_fig(fig, FIG, "myeloid_subcluster_composition_by_dpi"):
        rec.add_output(p)
    return {"lab_long": lab_long, "lab_med": lab_med, "sub_long": sub_long, "med": med, "check": check}


def umap_figures(sub: ad.AnnData, primary: str, verdict: str, rec: RunRecord) -> None:
    obs = sub.obs
    um = sub.obsm["X_umap"]
    # 1. subclusters versus held-out labels
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    sc.pl.umap(sub, color=primary, ax=axes[0], show=False, legend_loc="on data", size=8,
               title=f"myeloid subclusters (Leiden {PRIMARY_RES}, k {RULES['knn_k']})")
    sc.pl.umap(sub, color="label", ax=axes[1], show=False, legend_fontsize=7, size=8,
               title="deposited label (held out of the embedding)")
    fig.suptitle(f"Myeloid compartment of the Ki67 atlas, n={sub.n_obs:,}; grading verdict: {verdict}")
    fig.tight_layout()
    for p in save_fig(fig, FIG, "UMAP_myeloid_subclusters_vs_labels"):
        rec.add_output(p)

    # 2 and 3. per-day panels, equal cell number per panel
    n_min = int(obs.groupby("day").size().min())
    rec.set("umap_cells_drawn_per_panel", n_min)
    rng = np.random.default_rng(RANDOM_SEED)
    cmap = plt.get_cmap("tab20")
    for colour_by, name in [("label", "UMAP_myeloid_by_dpi_label"), (primary, "UMAP_myeloid_by_dpi_subcluster")]:
        fig, axes = plt.subplots(2, 4, figsize=(16, 8.6))
        for ax, d in zip(axes.ravel(), DAYS):
            ax.scatter(um[:, 0], um[:, 1], s=1.0, c="#d9d9d9", linewidths=0, rasterized=True)
            where = np.where((obs["day"] == d).to_numpy())[0]
            pick = rng.choice(where, n_min, replace=False)
            pick = pick[rng.permutation(len(pick))]
            vals = obs[colour_by].astype(str).to_numpy()[pick]
            if colour_by == "label":
                cols = [LABEL_COLOURS[v] for v in vals]
            else:
                cols = [cmap(int(v) % 20) for v in vals]
            ax.scatter(um[pick, 0], um[pick, 1], s=4, c=cols, linewidths=0, rasterized=True)
            n_animals = obs.loc[obs["day"] == d, "sample_id"].nunique()
            ax.set_title(f"{d} dpi, {PHASE[d]}\n{n_animals} animals, {len(where):,} cells, {n_min:,} drawn", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ax.spines.values():
                s.set_visible(False)
        if colour_by == "label":
            keys = [k for k in MYELOID_LABELS + ["other lineage"] if k in set(obs["label"].astype(str))]
            handles = [plt.Line2D([], [], marker="o", ls="", color=LABEL_COLOURS[k], label=k, markersize=6) for k in keys]
            what = "deposited label"
        else:
            keys = sorted(obs[primary].astype(str).unique(), key=int)
            handles = [plt.Line2D([], [], marker="o", ls="", color=cmap(int(k) % 20), label=f"sub {k}", markersize=6) for k in keys]
            what = f"blind subcluster (Leiden {PRIMARY_RES})"
        fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 10), frameon=False,
                   bbox_to_anchor=(0.5, 0.0))
        fig.suptitle(f"Myeloid compartment by day post infection ({what}; shared coordinates; "
                     "equal cell number per panel)", y=0.995)
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        for p in save_fig(fig, FIG, name):
            rec.add_output(p)


def markers_and_panel(sub: ad.AnnData, primary: str, rec: RunRecord) -> tuple[pd.DataFrame, pd.DataFrame, list]:
    sc.tl.rank_genes_groups(sub, primary, method="wilcoxon", pts=True, key_added="rank_sub")
    mk = sc.get.rank_genes_groups_df(sub, group=None, key="rank_sub")
    top = (mk.sort_values(["group", "scores"], ascending=[True, False])
             .groupby("group", observed=True).head(20))
    p = TAB / "myeloid_subcluster_markers_top20.csv"
    top.to_csv(p, index=False)
    rec.add_output(p)
    top10 = pd.DataFrame({g: list(top.loc[top["group"] == g, "names"].head(10))
                          for g in sorted(sub.obs[primary].astype(str).unique(), key=int)})
    genes = [g for gl in PANEL.values() for g in gl if g in sub.var_names]
    absent = [g for gl in PANEL.values() for g in gl if g not in sub.var_names]
    det = (sub[:, genes].layers["counts"] > 0).toarray()
    panel = (pd.DataFrame(det, columns=genes, index=sub.obs_names)
               .groupby(sub.obs[primary].astype(str)).mean().round(3))
    panel = panel.loc[sorted(panel.index, key=int)]
    panel.index.name = "subcluster"
    p = TAB / "myeloid_panel_detection_fraction.csv"
    panel.to_csv(p)
    rec.add_output(p)
    rec.set("panel_genes_absent", absent)
    try:
        present = {k: [g for g in v if g in sub.var_names] for k, v in PANEL.items()}
        present = {k: v for k, v in present.items() if v}
        dp = sc.pl.dotplot(sub, var_names=present, groupby=primary, standard_scale="var",
                           show=False, return_fig=True)
        dp.savefig(str(FIG / "dotplot_myeloid_panel_by_subcluster.png"), dpi=200, bbox_inches="tight")
        plt.close("all")
        rec.add_output(FIG / "dotplot_myeloid_panel_by_subcluster.png")
    except Exception as exc:  # noqa: BLE001
        log(f"  WARNING dotplot: {exc}")
        rec.set("dotplot_warning", str(exc))
        plt.close("all")
    return top10, panel, absent


def write_metadata(sub: ad.AnnData, keys: list[str], rec: RunRecord) -> None:
    cols = ["sample_id", "day", "tamoxifen_start_day", "sex", "trace_call", "cell_cycle_phase",
            "leiden_cluster", "author_celltype", "author_lineage", "label"] + keys
    out = sub.obs[cols].copy()
    out.insert(0, "cell_id", sub.obs_names)
    out["UMAP_1"] = sub.obsm["X_umap"][:, 0].round(4)
    out["UMAP_2"] = sub.obsm["X_umap"][:, 1].round(4)
    p = TAB / "myeloid_cell_metadata.csv"
    out.to_csv(p, index=False)
    rec.add_output(p)


def write_readme(sub: ad.AnnData, primary: str, grades: dict, comp: dict, top10: pd.DataFrame,
                 panel: pd.DataFrame, absent: list, rec: RunRecord) -> None:
    r = rec.record["results"]
    per_day = (sub.obs.groupby("day").agg(n_animals=("sample_id", "nunique"), n_cells=("sample_id", "size"))
                 .reindex(DAYS))
    per_day.insert(0, "phase", [PHASE[d] for d in per_day.index])
    per_day.index.name = "dpi"
    gtab = grades[primary]["table"]
    med = comp["med"][MYELOID_LABELS].round(2)
    med.index.name = "dpi"
    chk = comp["check"]
    sens = pd.DataFrame([{"resolution": k.replace("sub_r", ""), **{kk: vv for kk, vv in v.items() if kk != "table" and kk != "top_label_by_subcluster"}}
                         for k, v in grades.items()])
    lines = [
        "# Myeloid compartment of the Ki67 atlas (generated)",
        "",
        "Generated by `analysis/scripts/11_myeloid_focus.py`; the frozen rules, inputs and results",
        "are in [`run_record.json`](run_record.json). Owner retain/reject review pending.",
        "",
        "**What is being reproduced.** Niethamer et al. 2025 (Figure 3) report that alveolar",
        "macrophages (aMAC) are lost and inflammatory monocytes (iMON) expand at 6 dpi, and that the",
        "aMAC pool is reconstituted over the following three weeks from resident aMACs and from an",
        "iMON-derived trajectory. The compartment here is defined from the blind atlas clustering",
        "(clusters 5, 17, 24) within the 25-sample Ki67 atlas, re-embedded without the deposited",
        "labels, and graded against them afterwards.",
        "",
        "## Cohort",
        "",
        f"{r['n_cells']:,} cells from {r['n_animals']} animals; {r['n_cells_other_lineage_label']} of them carry a",
        "non-myeloid deposited label and are excluded from the label composition. Cells per atlas cluster:",
        ", ".join(f"{k}: {v:,}" for k, v in sorted(r["cells_per_cluster"].items(), key=lambda kv: int(kv[0]))) + ".",
        "",
        df_to_markdown(per_day),
        "",
        "## 1. Blind subclusters against the held-out labels",
        "",
        "![Subclusters versus labels](figures/UMAP_myeloid_subclusters_vs_labels.png)",
        "",
        f"Verdict by the frozen rule at the primary resolution ({PRIMARY_RES}): **{r['verdict']}**. Fraction of",
        f"labelled cells in subclusters with purity >= 0.75: {grades[primary]['fraction_labelled_cells_in_pure_subclusters']};",
        f"cell-weighted mean label entropy: {grades[primary]['cell_weighted_mean_entropy']}.",
        "",
        df_to_markdown(gtab, index=False),
        "",
        "Sensitivity to resolution:",
        "",
        df_to_markdown(sens, index=False),
        "",
        "## 2. The compartment by day post infection",
        "",
        "![Myeloid UMAP by dpi, deposited label](figures/UMAP_myeloid_by_dpi_label.png)",
        "",
        "![Myeloid UMAP by dpi, blind subcluster](figures/UMAP_myeloid_by_dpi_subcluster.png)",
        "",
        f"Shared subset coordinates; every panel draws {r['umap_cells_drawn_per_panel']:,} cells (the smallest",
        "per-day count, seed 0) over all subset cells in grey.",
        "",
        "## 3. Within-myeloid composition per animal",
        "",
        "![Label composition by dpi](figures/myeloid_label_composition_by_dpi.png)",
        "",
        "Median per-animal percent of each deposited label among the animal's myeloid-labelled cells",
        "(per-animal values in `tables/myeloid_label_composition_per_animal.csv`; the blind subcluster",
        "version is `tables/myeloid_subcluster_composition_per_animal.csv` and",
        "`figures/myeloid_subcluster_composition_by_dpi.png`):",
        "",
        df_to_markdown(med),
        "",
        "**Figure 3 check by the frozen rule:** aMAC median at 6 dpi below baseline: "
        f"{chk['aMAC_lost_at_6dpi']}; iMON median at 6 dpi above baseline: {chk['iMON_expanded_at_6dpi']}; "
        f"aMAC median at 25 or 42 dpi above the 6 dpi value: {chk['aMAC_reconstituted_by_25_or_42dpi']}. "
        f"Verdict: **{chk['verdict']}**. Two animals per active-repair day; this is a description, not a test.",
        "",
        "![Blind subcluster composition](figures/myeloid_subcluster_composition_by_dpi.png)",
        "",
        "## 4. Markers and comparator panel (ranking devices)",
        "",
        "Top 10 Wilcoxon genes per primary subcluster (full top 20 with scores in",
        "`tables/myeloid_subcluster_markers_top20.csv`):",
        "",
        df_to_markdown(top10, index=False),
        "",
        "Detection fraction of the comparator panel per subcluster (`tables/myeloid_panel_detection_fraction.csv`;",
        "the dot plot is `figures/dotplot_myeloid_panel_by_subcluster.png`).",
        f"Panel genes absent from the object: {', '.join(absent) or 'none'}.",
        "",
        "## Caveats",
        "",
        "- The unit is the animal; active-repair days carry two animals each and no P values are computed.",
        "- Myeloid cells are 10 to 15 percent of every library by the MACS recombination, so only",
        "  within-myeloid fractions are interpretable; the absolute myeloid influx is not measurable here.",
        "- The subset inherits the atlas QC and doublet calls; no ambient-RNA correction was applied",
        "  (raw droplet matrices are absent for this series).",
        "- The iMON-to-aMAC trajectory of the paper is not fitted here; this analysis reads states and",
        "  composition only.",
        "",
    ]
    p = OUT / "README.md"
    p.write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(p)


def main() -> int:
    for d in (OUT, FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    rec = RunRecord(OUT / "run_record.json", "myeloid compartment of the Ki67 atlas: blind subclusters, "
                    "per-day UMAP and per-animal composition", RULES,
                    notes="rules frozen before the object was opened; deposited labels held out of the "
                          "embedding and used afterwards for grading and descriptive composition")
    rec.add_output(OUT / "run_record.json")
    sub = load_subset(rec)
    primary = embed(sub, rec)
    keys = [f"sub_r{res}" for res in [PRIMARY_RES] + SENSITIVITY_RES]
    grades = {}
    for key in keys:
        df, tab, summary = grade(sub, key)
        summary["top_label_by_subcluster"] = dict(zip(df["subcluster"], df["top_label"]))
        grades[key] = {"table": df, **summary}
        p = TAB / f"myeloid_grading_{key}.csv"
        df.to_csv(p, index=False)
        rec.add_output(p)
        p = TAB / f"myeloid_subcluster_by_author_label_{key}.csv"
        tab.to_csv(p)
        rec.add_output(p)
    rec.set("grading", {k: {kk: vv for kk, vv in v.items() if kk != "table"} for k, v in grades.items()})
    frac = grades[primary]["fraction_labelled_cells_in_pure_subclusters"]
    verdict = ("resolved by resolution"
               if frac >= RULES["grading"]["resolved_if_fraction_of_labelled_cells_in_pure_subclusters_ge"]
               else "not resolved by resolution")
    rec.set("verdict", verdict)
    log(f"grading verdict at Leiden {PRIMARY_RES}: {verdict} (pure fraction {frac})")

    comp = composition(sub, primary, rec)
    umap_figures(sub, primary, verdict, rec)
    top10, panel, absent = markers_and_panel(sub, primary, rec)
    write_metadata(sub, keys, rec)
    write_readme(sub, primary, grades, comp, top10, panel, absent, rec)
    rec.finish()
    print()
    print(comp["med"][MYELOID_LABELS].round(1).to_string())
    print(f"\nFigure 3 check: {comp['check']['verdict']}")
    print(f"wrote {OUT.relative_to(ANALYSIS.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
