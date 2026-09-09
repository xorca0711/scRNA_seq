#!/usr/bin/env python
"""Trial S5: subcluster mouse cluster 5 (the highest label-entropy cluster) and re-grade.

Portfolio purpose. Cluster 5 is a myeloid mixture whose deposited labels
disagree (monocytes, macrophages, dendritic cells), the same case the HLCA
reports as its worst label-entropy cluster. Resolving it matters for the
stem-niche communication angle (Nabhan fit) and the inflammation-resolution
and immune-epithelial programme angle (Saxton fit), both of which need
macrophage and monocyte states that are real at the animal level. The
HLCA's monocyte-derived macrophage subtypes (SPP1-high, CCL2-high,
C1QA-high, MARCO-high) are used as human comparators through mouse
orthologs. No infection or interferon framing is used.

Frozen rules (written to the run record before any data is read):

* Cells: all cells of atlas cluster 5. No batch correction (root decision 2).
* Features: seurat_v3 HVGs on raw counts, 2,000 genes, no batch key.
* Embedding: scale (max 10) on HVGs, PCA 30 components, random_state 0.
* Graph and clustering, following the HLCA level-2 rule: kNN k 30 on the
  30 PCs; Leiden resolution 0.2 (primary), 0.5 (sensitivity); igraph
  flavour, 2 iterations, random_state 0.
* Grading against the deposited labels (held out of everything above):
  purity = largest label fraction; label entropy with the 0.56 threshold.
  Verdict "resolved by resolution" if >= 75% of labelled cells fall in
  subclusters with purity >= 0.75 at the primary resolution; otherwise
  "not resolved by resolution".
* Interpretation only: Wilcoxon top genes per subcluster (a ranking
  device), and detection fraction of a fixed comparator panel.
* Unit: animal for composition (max single-sample fraction per subcluster
  reported); no P values.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trial_utils import REPO, RunRecord, df_to_markdown, read_csr_rows, shannon  # noqa: E402

OUT = HERE / "s5_cluster5_subclusters"
OUT.mkdir(exist_ok=True)
OBJ = REPO / "analysis" / "GSE262927" / "processed" / "final_clustered.h5ad"
CLUSTER = "5"
PANEL = {
    "HLCA MDM SPP1-high": ["Spp1", "Lpl", "Chil3", "Mmp9", "Fdx1"],
    "HLCA MDM CCL2-high": ["Ccl2", "Il1rn", "S100a8"],
    "HLCA MDM C1QA-high": ["C1qa", "Il18", "Trem2"],
    "HLCA MARCO-high": ["Marco", "Mcemp1"],
    "alveolar macrophage": ["Siglecf", "Ear2", "Plet1", "Fabp4"],
    "interstitial macrophage": ["Cd163", "Mrc1", "Folr2", "Lyve1", "C1qb"],
    "classical monocyte": ["Ly6c2", "Ccr2", "Plac8"],
    "non-classical monocyte": ["Cx3cr1", "Ace", "Nr4a1", "Ear2"],
    "cDC1": ["Xcr1", "Clec9a"],
    "cDC2": ["Cd209a", "Itgax", "H2-Ab1"],
    "migratory DC": ["Ccr7", "Fscn1"],
    "pDC": ["Siglech", "Bst2"],
    "neutrophil": ["S100a9", "Retnlg", "Csf3r"],
    "proliferating": ["Mki67", "Top2a"],
}
RULES = {
    "cluster": CLUSTER,
    "hvg": {"flavor": "seurat_v3", "n_top_genes": 2000, "layer": "counts", "batch_key": None},
    "scale_max": 10, "n_pcs": 30, "knn_k": 30,
    "leiden": {"primary_resolution": 0.2, "sensitivity_resolution": 0.5, "flavor": "igraph", "n_iterations": 2, "random_state": 0},
    "batch_correction": "none",
    "purity_threshold": 0.75, "resolved_if_fraction_of_labelled_cells_in_pure_subclusters_ge": 0.75,
    "label_entropy_high": 0.56,
    "comparator_panel": PANEL,
    "unit": "animal for composition",
}


def main() -> None:
    rec = RunRecord(OUT / "s5_run_record.json", "S5 subcluster mouse cluster 5 and re-grade", RULES,
                    notes="rules frozen before reading the object; deposited labels held out of clustering")
    import anndata as ad  # noqa: E402
    import scanpy as sc  # noqa: E402
    from anndata.io import read_elem  # noqa: E402
    sc.settings.verbosity = 0
    rec.add_input(OBJ)
    with h5py.File(OBJ, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        n_cols = f["X"].attrs["shape"][1]
        idx = np.where(obs["leiden_cluster"].astype(str).to_numpy() == CLUSTER)[0]
        X = read_csr_rows(f["X"], idx, n_cols)
        C = read_csr_rows(f["layers"]["counts"], idx, n_cols)
    sub = ad.AnnData(X=X, obs=obs.iloc[idx].copy(), var=var.copy(), layers={"counts": C})
    sub.obs_names = obs.index[idx]
    rec.set("n_cells", int(sub.n_obs))

    sc.pp.highly_variable_genes(sub, flavor="seurat_v3", n_top_genes=RULES["hvg"]["n_top_genes"], layer="counts")
    emb = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=RULES["scale_max"])
    sc.tl.pca(emb, n_comps=RULES["n_pcs"], random_state=0)
    sub.obsm["X_pca"] = emb.obsm["X_pca"]
    sc.pp.neighbors(sub, n_neighbors=RULES["knn_k"], n_pcs=RULES["n_pcs"], use_rep="X_pca", random_state=0)
    for res in (RULES["leiden"]["primary_resolution"], RULES["leiden"]["sensitivity_resolution"]):
        sc.tl.leiden(sub, resolution=res, key_added=f"sub_r{res}", flavor="igraph", n_iterations=2, random_state=0, directed=False)
    sc.tl.umap(sub, random_state=0)
    primary = f"sub_r{RULES['leiden']['primary_resolution']}"

    labelled = sub.obs["has_author_metadata"].astype(str).eq("True")
    lab = sub.obs["author_celltype"].astype(str)
    grades = {}
    for key in (primary, f"sub_r{RULES['leiden']['sensitivity_resolution']}"):
        rows = []
        tab = pd.crosstab(sub.obs.loc[labelled, key].astype(str), lab[labelled])
        for scl in sorted(sub.obs[key].astype(str).unique(), key=int):
            n_all = int((sub.obs[key].astype(str) == scl).sum())
            if scl in tab.index:
                counts = tab.loc[scl]
                n_lab = int(counts.sum())
                purity = float(counts.max() / n_lab) if n_lab else np.nan
                ent = shannon((counts / n_lab).to_numpy()) if n_lab else np.nan
                top = counts.idxmax() if n_lab else ""
            else:
                n_lab, purity, ent, top = 0, np.nan, np.nan, ""
            samp = sub.obs.loc[sub.obs[key].astype(str) == scl, "sample_id"].astype(str).value_counts()
            rows.append({"subcluster": scl, "n_cells": n_all, "n_labelled": n_lab, "labelled_fraction": round(n_lab / n_all, 3),
                         "top_label": top, "purity": round(purity, 3) if n_lab else None,
                         "label_entropy": round(ent, 3) if n_lab else None,
                         "entropy_flag": ("HIGH" if ent > RULES["label_entropy_high"] else "low") if n_lab else "NA",
                         "n_samples": int(len(samp)), "max_sample_fraction": round(float(samp.iloc[0] / n_all), 3), "dominant_sample": samp.index[0]})
        df = pd.DataFrame(rows)
        n_lab_total = int(labelled.sum())
        in_pure = int(df.loc[df["purity"].fillna(0) >= RULES["purity_threshold"], "n_labelled"].sum())
        grades[key] = {"table": df, "fraction_labelled_cells_in_pure_subclusters": round(in_pure / n_lab_total, 4),
                       "cell_weighted_mean_entropy": round(float(np.nansum(df["label_entropy"].fillna(0) * df["n_labelled"]) / n_lab_total), 4)}
        df.to_csv(OUT / f"s5_grading_{key}.csv", index=False)
        rec.add_output(OUT / f"s5_grading_{key}.csv")
        tab.to_csv(OUT / f"s5_subcluster_by_author_label_{key}.csv")
        rec.add_output(OUT / f"s5_subcluster_by_author_label_{key}.csv")
    frac = grades[primary]["fraction_labelled_cells_in_pure_subclusters"]
    verdict = "resolved by resolution" if frac >= RULES["resolved_if_fraction_of_labelled_cells_in_pure_subclusters_ge"] else "not resolved by resolution"
    rec.set("verdict", verdict)
    rec.set("grading_summary", {k: {kk: vv for kk, vv in v.items() if kk != "table"} for k, v in grades.items()})
    rec.set("atlas_cluster5_label_entropy_before", round(shannon((lab[labelled].value_counts() / labelled.sum()).to_numpy()), 4))

    # markers (ranking device) and comparator panel
    sc.tl.rank_genes_groups(sub, primary, method="wilcoxon", n_genes=10)
    mk = pd.DataFrame({g: [t for t in sub.uns["rank_genes_groups"]["names"][g]] for g in sub.obs[primary].cat.categories})
    mk.to_csv(OUT / "s5_top10_markers_primary.csv", index=False)
    rec.add_output(OUT / "s5_top10_markers_primary.csv")
    genes = [g for gl in PANEL.values() for g in gl if g in sub.var_names]
    absent = [g for gl in PANEL.values() for g in gl if g not in sub.var_names]
    det = (sub[:, genes].layers["counts"] > 0).toarray()
    panel = pd.DataFrame(det, columns=genes, index=sub.obs_names).groupby(sub.obs[primary].astype(str)).mean().round(3)
    panel.index.name = "subcluster"
    panel.to_csv(OUT / "s5_comparator_panel_detection_fraction.csv")
    rec.add_output(OUT / "s5_comparator_panel_detection_fraction.csv")
    rec.set("panel_genes_absent", absent)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5))
    sc.pl.umap(sub, color=primary, ax=axes[0], show=False, legend_loc="on data", size=8, title=f"cluster 5 subclusters (Leiden {RULES['leiden']['primary_resolution']}, k 30)")
    sub.obs["author_label_or_NA"] = np.where(labelled, lab, "not labelled")
    sc.pl.umap(sub, color="author_label_or_NA", ax=axes[1], show=False, legend_fontsize=6, size=8, title="deposited label (held out)")
    fig.suptitle(f"Mouse atlas cluster 5, n={sub.n_obs}; verdict: {verdict}")
    fig.tight_layout(); fig.savefig(OUT / "s5_umap_subclusters_vs_labels.png", dpi=110); plt.close(fig)
    rec.add_output(OUT / "s5_umap_subclusters_vs_labels.png")

    df = grades[primary]["table"]
    lines = ["# Trial S5 output: mouse cluster 5 subclusters", "",
             f"Verdict by the frozen rules: **{verdict}**. Fraction of labelled cells in subclusters with purity >= 0.75: {frac}; "
             f"cell-weighted mean label entropy after subclustering: {grades[primary]['cell_weighted_mean_entropy']} "
             f"(atlas cluster 5 before: {rec.record['results']['atlas_cluster5_label_entropy_before']}).", "",
             "## Grading at the primary resolution", "", df_to_markdown(df, index=False), "",
             "## Top 10 Wilcoxon genes per subcluster (ranking device)", "", df_to_markdown(mk, index=False), "",
             "## Comparator panel, detection fraction per subcluster", "", df_to_markdown(panel), "",
             f"Panel genes absent from the object: {', '.join(absent) or 'none'}.", ""]
    (OUT / "s5_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "s5_summary.md")
    rec.finish()
    print("\n".join(lines[:40]))


if __name__ == "__main__":
    main()
