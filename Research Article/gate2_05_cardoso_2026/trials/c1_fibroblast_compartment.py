#!/usr/bin/env python
"""Trial C1: recover the fibrotic fibroblast state of GSE316241, or fail to (Gate 1).

The owner's Gate 1: take the mesenchymal series, recover the
Pdgfrb+Runx1+Tnc+ "reprogrammed" fibroblast subset, and ask how cleanly it
separates from homeostatic Pdgfra+ alveolar fibroblasts. If it does not
separate, or splits differently from the published account, that is the
result and the trial stops there rather than pressing on. The open branch the
owner named is asked in the same run: is the fibrotic subset one state or a
gradient out of the alveolar fibroblasts?

Design facts that constrain every reading below, from the deposit itself
(trial C0): two libraries, one per genotype, each a pool of three mice, both
two weeks after induction. Library, genotype and batch are one variable.
Nothing here is a tested group difference and no P value is computed for a
genotype contrast; per-cluster Wilcoxon statistics are used only as a
within-object ranking device, never as a between-genotype test.

Frozen rules (written to the run record before any matrix is read):

* Input: the two GSE316241 triplets. Features whose identifier is not an
  Ensembl gene ID are dropped from the matrix (C0 rule).
* QC: this repository's per-sample MAD rule (pipeline_utils.derive_thresholds,
  the same call the GSE262927 and GSE178360 runs used), applied per library.
  The paper's own stated thresholds (>1,000 genes, 2,000 to 50,000 UMIs,
  <10% mitochondrial) are applied in parallel and reported as a comparison
  only; they do not decide which cells are kept.
* Doublets: scanpy's Scrublet per library before merging, with the 10x
  multiplet-rate prior of 0.8% per 1,000 cells recovered, and the same
  automatic-threshold sanity check the pipeline uses (fall back to the
  expected-rate quantile if the automatic call is under a fifth or over three
  times the prior). Predicted doublets are removed.
* Batch correction: none in the primary embedding, because library, genotype
  and batch are the same variable here and correcting on it would remove the
  effect under study. Harmony on library is run as a disclosed sensitivity
  and reported beside the primary, not in place of it.
* Embedding: log1p(CP10K); seurat_v3 HVGs on raw counts, 2,000 genes; scale
  to max 10; PCA 30 components; kNN k = 30; UMAP; all seeds 0.
* Clustering: Leiden (igraph flavour, 2 iterations, seed 0) at resolution
  0.5 as the primary, with 0.2 and 1.0 as sensitivities. 0.5 is fixed from
  this repository's trial S5, not chosen here.
* Population calls: scanpy score_genes on the paper's own marker sets
  (cardoso_2026_extracts.json), ctrl_size 50, n_bins 25, seed 0; a cluster is
  called by the modal per-cell argmax, and the call is "confident" only if
  the mode holds at least 50% of its cells.
* Fibroblast subset: clusters whose confident call is alveolar, adventitial
  or reprogrammed fibroblast. The subset is re-embedded from raw counts with
  the same parameters, because a subset changes the feature selection and the
  doublet composition and must not inherit the parent object's.
* THE TEST, frozen. "The reprogrammed fibroblast state is recovered" if some
  fibroblast subcluster satisfies all three: (i) at least 80% of its cells
  come from the Red2Kras library; (ii) its mean reprogrammed-fibroblast score
  is the highest of all fibroblast subclusters; (iii) Tnc is detected in at
  least 40% of its cells (the paper names Tnc the top differentially
  expressed gene of this population). Failing any one of the three is
  recorded as a failure to recover, and the trial reports what was found
  instead.
* "It separates cleanly" if at least 75% of that subcluster's cells have a
  majority of their 30 nearest neighbours inside the same subcluster.
* Over-clustering guard: at every resolution, a subcluster is "unsupported"
  if no gene is detected in at least 30% of its cells and at most 20% of the
  cells outside it (the HLCA cluster-marker rule). Unsupported subclusters
  are counted and named; a subset whose extra structure is unsupported is
  reported as over-clustering, not as populations.
* Doublet guard: the candidate fibrotic subcluster is flagged if its median
  Scrublet score exceeds the 90th percentile of the fibroblast subset, or if
  pericyte (Cspg4, Notch3, Postn) or smooth-muscle (Myh11, Tagln) markers are
  detected in more than half its cells; a Pdgfrb+Acta2+ population is exactly
  what a fibroblast-mural doublet would look like.
* One state or a gradient: within Red2Kras fibroblasts, fit one- and
  two-component Gaussian mixtures to the reprogrammed-fibroblast score (seed
  0). Two states if the two-component BIC is lower by at least 10; gradient
  if the one-component BIC is lower; undecided in between. Reported beside a
  diffusion pseudotime rooted at the alveolar-fibroblast centroid.
* Unit: the library. There is no animal-level unit in this deposit. No P
  value is computed for any genotype contrast.
"""

from __future__ import annotations

import gc
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import (ENSEMBL_ID_RE, RAW, REPO, RunRecord,  # noqa: E402
                           df_to_markdown, read_mtx_triplet)

sys.path.insert(0, str(REPO / "analysis" / "scripts"))
from pipeline_utils import derive_thresholds, apply_thresholds, make_unique  # noqa: E402

OUT = HERE / "c1_fibroblast_compartment"
OUT.mkdir(exist_ok=True)
EXTRACTS = json.loads((HERE.parent / "cardoso_2026_extracts.json").read_text(encoding="utf-8"))
SETS = EXTRACTS["marker_sets"]

LIBRARIES = [
    {"library": "Expt1_Confetti_mesenchyme", "gsm": "GSM9447763", "genotype": "Confetti", "condition": "homeostasis"},
    {"library": "Expt1_Red2Kras_mesenchyme", "gsm": "GSM9447764", "genotype": "Red2Kras", "condition": "oncogenesis"},
]

# The paper's nine mesenchymal populations, by its own marker sets.
POPULATION_SETS = {
    "alveolar fibroblast": SETS["alveolar_fibroblast"],
    "adventitial fibroblast": SETS["adventitial_fibroblast"],
    "reprogrammed fibroblast": SETS["reprogrammed_fibrotic_fibroblast"],
    "peri-bronchial fibroblast": SETS["peri_bronchial_fibroblast"],
    "smooth muscle": SETS["smooth_muscle"],
    "pericyte": SETS["pericyte"],
    "mesothelium": SETS["mesothelium"],
    "mesothelial-like": SETS["mesothelial_like_enriched_in_tumour"],
    "proliferating": SETS["proliferating"],
}
FIBROBLAST_CALLS = {"alveolar fibroblast", "adventitial fibroblast", "reprogrammed fibroblast"}
DOUBLET_GUARD = {"pericyte": ["Cspg4", "Notch3", "Postn"], "smooth muscle": ["Myh11", "Tagln"]}

PAPER_QC = {"min_genes": 1000, "min_counts": 2000, "max_counts": 50000, "max_pct_mt": 10}
RULES = {
    "accession": "GSE316241",
    "libraries": LIBRARIES,
    "design_constraint": ("one library per genotype, three mice pooled per library, single time point "
                          "(2 weeks); library, genotype and batch are one variable; no P value for a "
                          "genotype contrast"),
    "qc": {"primary": "repository per-sample MAD rule (pipeline_utils.derive_thresholds)",
           "reported_comparison": PAPER_QC},
    "doublets": {"caller": "scanpy Scrublet per library before merging",
                 "prior": "0.008 per 1000 cells recovered, clipped to [0.01, 0.15]",
                 "automatic_threshold_rejected_if": "called rate < 0.2x or > 3x the prior",
                 "action": "predicted doublets removed"},
    "batch_correction": {"primary": "none (library is genotype)", "sensitivity": "Harmony on library"},
    "embedding": {"normalisation": "log1p(CP10K)", "hvg": "seurat_v3 on counts, 2000",
                  "scale_max": 10, "n_pcs": 30, "knn_k": 30, "seed": 0},
    "clustering": {"primary_resolution": 0.5, "sensitivity_resolutions": [0.2, 1.0],
                   "flavor": "igraph", "n_iterations": 2, "seed": 0,
                   "origin_of_primary": "fixed from repository trial S5, not chosen here"},
    "population_calls": {"method": "scanpy score_genes on the paper's marker sets",
                         "ctrl_size": 50, "n_bins": 25, "seed": 0,
                         "cluster_call": "modal per-cell argmax; confident if the mode holds >= 50% of cells",
                         "marker_sets": POPULATION_SETS},
    "fibroblast_subset": "clusters whose confident call is alveolar, adventitial or reprogrammed fibroblast; re-embedded from raw counts with the same parameters",
    "recovery_test": {"red2kras_fraction_min": 0.80,
                      "must_have_highest_mean_reprogrammed_score": True,
                      "tnc_detection_min": 0.40},
    "clean_separation": {"rule": "fraction of the subcluster's cells whose 30 nearest neighbours are mostly in the same subcluster", "threshold": 0.75},
    "over_clustering_guard": {"marker_rule": "detected in >= 30% in-group and <= 20% out-group",
                              "unsupported_subcluster": "no gene satisfies the rule"},
    "doublet_guard": {"scrublet_rule": "median score above the subset's 90th percentile",
                      "marker_rule": "pericyte or smooth-muscle markers detected in > 50% of cells",
                      "panels": DOUBLET_GUARD},
    "state_or_gradient": {"method": "Gaussian mixture on the reprogrammed-fibroblast score within Red2Kras fibroblasts, seed 0",
                          "two_states_if": "BIC(2) < BIC(1) - 10", "gradient_if": "BIC(1) < BIC(2)",
                          "companion": "diffusion pseudotime rooted at the alveolar-fibroblast centroid"},
    "unit": "the library; there is no animal-level unit in this deposit",
}


def load_library(entry: dict, rec: RunRecord):
    import anndata as ad
    import scanpy as sc

    root = RAW / "GSE316241" / "GSE316241_RAW"
    stem = f"{entry['gsm']}_{entry['library']}_"
    matrix = root / f"{stem}matrix.mtx.gz"
    features = root / f"{stem}features.tsv.gz"
    barcodes = root / f"{stem}barcodes.tsv.gz"
    for path in (matrix, features, barcodes):
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(matrix, features, barcodes)

    is_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    X, var = X[:, is_gene], var.loc[is_gene].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))

    adata = ad.AnnData(X=X.astype(np.float32), var=var,
                       obs=pd.DataFrame(index=pd.Index([f"{entry['library']}_{b}" for b in bc])))
    for key in ("library", "genotype", "condition", "gsm"):
        adata.obs[key] = entry[key]

    symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
    adata.var["ribo"] = (symbols.str.startswith(("Rps", "Rpl"))
                         & ~symbols.str.contains("-ps", case=False)).to_numpy()
    adata.var["hb"] = symbols.str.startswith(("Hba", "Hbb")).to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20],
                               log1p=True, inplace=True)

    th = derive_thresholds(adata.obs, entry["library"])
    keep = apply_thresholds(adata.obs, th).to_numpy()
    paper_keep = ((adata.obs["n_genes_by_counts"] > PAPER_QC["min_genes"])
                  & (adata.obs["total_counts"] > PAPER_QC["min_counts"])
                  & (adata.obs["total_counts"] < PAPER_QC["max_counts"])
                  & (adata.obs["pct_counts_mt"] < PAPER_QC["max_pct_mt"])).to_numpy()
    qc_row = {
        "library": entry["library"], "genotype": entry["genotype"],
        "barcodes": int(adata.n_obs),
        "median_counts": int(np.median(adata.obs["total_counts"])),
        "median_genes": int(np.median(adata.obs["n_genes_by_counts"])),
        "median_pct_mt": round(float(np.median(adata.obs["pct_counts_mt"])), 2),
        "kept_repository_rule": int(keep.sum()),
        "kept_paper_rule": int(paper_keep.sum()),
        "kept_by_both": int((keep & paper_keep).sum()),
        "thresholds": th.rationale,
    }
    adata = adata[keep].copy()

    expected = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=expected)
    auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
    method = "scrublet automatic threshold"
    if auto < 0.2 * expected or auto > 3 * expected:
        cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - expected))
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
        method = f"top {100 * expected:.2f}% of scores (automatic rejected: it called {100 * auto:.2f}%)"
    qc_row.update({"expected_doublet_rate": round(expected, 4),
                   "doublet_call_method": method,
                   "doublets_removed": int(adata.obs["predicted_doublet"].sum())})
    adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
    qc_row["cells_analysed"] = int(adata.n_obs)
    return adata, qc_row


def embed(adata, resolutions, key_prefix="leiden", harmony=False):
    import scanpy as sc

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=30, random_state=0)
    adata.obsm["X_pca"] = emb.obsm["X_pca"]
    use_rep = "X_pca"
    if harmony:
        import harmonypy
        ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["library"], max_iter_harmony=20)
        adata.obsm["X_harmony"] = np.asarray(ho.Z_corr).T
        use_rep = "X_harmony"
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep=use_rep, random_state=0)
    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f"{key_prefix}_{res}",
                     flavor="igraph", n_iterations=2, random_state=0, directed=False)
    sc.tl.umap(adata, random_state=0)
    del emb
    gc.collect()
    return adata


def score_populations(adata, sets: dict) -> pd.DataFrame:
    import scanpy as sc

    scores = {}
    for name, genes in sets.items():
        present = [g for g in genes if g in adata.var_names]
        if not present:
            continue
        sc.tl.score_genes(adata, present, score_name=f"score_{name}", ctrl_size=50,
                          n_bins=25, random_state=0)
        scores[name] = adata.obs[f"score_{name}"].to_numpy()
    frame = pd.DataFrame(scores, index=adata.obs_names)
    adata.obs["population_argmax"] = frame.idxmax(axis=1).to_numpy()
    return frame


def call_clusters(adata, key: str) -> pd.DataFrame:
    rows = []
    for cluster in sorted(adata.obs[key].astype(str).unique(), key=int):
        mask = adata.obs[key].astype(str).to_numpy() == cluster
        calls = adata.obs.loc[mask, "population_argmax"].value_counts()
        mode_fraction = float(calls.iloc[0] / mask.sum())
        libs = adata.obs.loc[mask, "genotype"].value_counts(normalize=True)
        rows.append({
            "cluster": cluster, "n_cells": int(mask.sum()),
            "call": calls.index[0], "mode_fraction": round(mode_fraction, 3),
            "confident": mode_fraction >= 0.50,
            "fraction_Red2Kras": round(float(libs.get("Red2Kras", 0.0)), 3),
            "fraction_Confetti": round(float(libs.get("Confetti", 0.0)), 3),
        })
    return pd.DataFrame(rows)


def detection(adata, genes, mask=None) -> dict:
    out = {}
    X = adata.layers["counts"]
    for gene in genes:
        if gene not in adata.var_names:
            out[gene] = np.nan
            continue
        col = np.asarray(X[:, adata.var_names.get_loc(gene)].todense()).ravel() > 0
        out[gene] = float(col[mask].mean() if mask is not None else col.mean())
    return out


def unsupported_subclusters(adata, key: str) -> dict:
    """HLCA cluster-marker rule: a subcluster with no gene at >=30% in / <=20% out."""
    X = (adata.layers["counts"] > 0)
    groups = adata.obs[key].astype(str).to_numpy()
    result = {}
    for cluster in sorted(set(groups), key=int):
        mask = groups == cluster
        if mask.sum() < 10:
            result[cluster] = {"n_cells": int(mask.sum()), "n_supporting_genes": 0, "supported": False}
            continue
        inside = np.asarray(X[mask].mean(axis=0)).ravel()
        outside = np.asarray(X[~mask].mean(axis=0)).ravel()
        hit = (inside >= 0.30) & (outside <= 0.20)
        genes = adata.var_names[hit].tolist()
        result[cluster] = {"n_cells": int(mask.sum()), "n_supporting_genes": int(hit.sum()),
                           "supported": bool(hit.sum() > 0),
                           "examples": genes[:8]}
    return result


def knn_self_majority(adata, key: str, cluster: str) -> float:
    graph = adata.obsp["connectivities"]
    groups = adata.obs[key].astype(str).to_numpy()
    idx = np.where(groups == cluster)[0]
    same = 0
    for i in idx:
        neighbours = graph[i].indices
        if len(neighbours) == 0:
            continue
        if (groups[neighbours] == cluster).mean() > 0.5:
            same += 1
    return float(same / len(idx)) if len(idx) else float("nan")


def main() -> None:
    import scanpy as sc
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc.settings.verbosity = 0
    rec = RunRecord(OUT / "c1_run_record.json",
                    "C1 fibroblast compartment of GSE316241 (Gate 1)", RULES,
                    notes="rules frozen before any matrix was read")

    parts, qc_rows = [], []
    for entry in LIBRARIES:
        adata, qc_row = load_library(entry, rec)
        parts.append(adata)
        qc_rows.append(qc_row)
        print(f"  {entry['library']}: {qc_row['cells_analysed']} cells analysed")
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT / "c1_qc_per_library.csv", index=False)
    rec.add_output(OUT / "c1_qc_per_library.csv")

    merged = ad.concat(parts, join="inner", label=None)
    merged.var = parts[0].var.loc[merged.var_names].copy()
    rec.set("cells_analysed_total", int(merged.n_obs))
    rec.set("genes_shared", int(merged.n_vars))
    del parts
    gc.collect()

    merged = embed(merged, [0.2, 0.5, 1.0])
    score_populations(merged, POPULATION_SETS)
    calls = call_clusters(merged, "leiden_0.5")
    calls.to_csv(OUT / "c1_mesenchymal_clusters.csv", index=False)
    rec.add_output(OUT / "c1_mesenchymal_clusters.csv")
    print(df_to_markdown(calls, index=False))

    fib_clusters = calls.loc[calls["confident"] & calls["call"].isin(FIBROBLAST_CALLS), "cluster"].tolist()
    rec.set("fibroblast_clusters_selected", fib_clusters)
    mask = merged.obs["leiden_0.5"].astype(str).isin(fib_clusters).to_numpy()
    fib = ad.AnnData(X=merged.layers["counts"][mask].copy(),
                     obs=merged.obs.loc[mask, ["library", "genotype", "condition", "gsm",
                                               "doublet_score", "total_counts",
                                               "n_genes_by_counts", "pct_counts_mt"]].copy(),
                     var=merged.var.copy())
    rec.set("fibroblast_cells", int(fib.n_obs))
    fib = embed(fib, [0.2, 0.5, 1.0], key_prefix="sub")
    score_populations(fib, POPULATION_SETS)
    sub_calls = call_clusters(fib, "sub_0.5")

    tnc = detection(fib, ["Tnc"])["Tnc"]
    rows = []
    for cluster in sub_calls["cluster"]:
        m = fib.obs["sub_0.5"].astype(str).to_numpy() == cluster
        row = sub_calls.loc[sub_calls["cluster"] == cluster].iloc[0].to_dict()
        row["mean_reprogrammed_score"] = round(float(fib.obs.loc[m, "score_reprogrammed fibroblast"].mean()), 4)
        row["mean_alveolar_score"] = round(float(fib.obs.loc[m, "score_alveolar fibroblast"].mean()), 4)
        row["median_doublet_score"] = round(float(np.median(fib.obs.loc[m, "doublet_score"])), 4)
        for gene, frac in detection(fib, ["Tnc", "Fst", "Runx1", "Runx2", "Acta2", "Pdgfrb",
                                          "Pdgfra", "Col13a1", "Lcn2", "Saa3", "Cxcl12"], m).items():
            row[f"det_{gene}"] = round(frac, 3) if frac == frac else np.nan
        for panel, genes in DOUBLET_GUARD.items():
            vals = [v for v in detection(fib, genes, m).values() if v == v]
            row[f"det_{panel.replace(' ', '_')}_max"] = round(max(vals), 3) if vals else np.nan
        rows.append(row)
    sub = pd.DataFrame(rows)
    sub.to_csv(OUT / "c1_fibroblast_subclusters.csv", index=False)
    rec.add_output(OUT / "c1_fibroblast_subclusters.csv")
    print(df_to_markdown(sub, index=False))

    # ---- the frozen test -------------------------------------------------
    top_score = sub.loc[sub["mean_reprogrammed_score"].idxmax(), "cluster"]
    candidates = sub[(sub["fraction_Red2Kras"] >= RULES["recovery_test"]["red2kras_fraction_min"])
                     & (sub["cluster"] == top_score)
                     & (sub["det_Tnc"] >= RULES["recovery_test"]["tnc_detection_min"])]
    recovered = bool(len(candidates))
    verdict = {"recovered": recovered,
               "highest_reprogrammed_score_subcluster": str(top_score),
               "criteria": {}}
    row = sub.loc[sub["cluster"] == top_score].iloc[0]
    verdict["criteria"] = {
        "red2kras_fraction": float(row["fraction_Red2Kras"]),
        "red2kras_fraction_passes": bool(row["fraction_Red2Kras"] >= 0.80),
        "has_highest_reprogrammed_score": True,
        "tnc_detection": float(row["det_Tnc"]) if row["det_Tnc"] == row["det_Tnc"] else None,
        "tnc_detection_passes": bool(row["det_Tnc"] >= 0.40) if row["det_Tnc"] == row["det_Tnc"] else False,
    }
    if recovered:
        verdict["clean_separation_fraction"] = round(knn_self_majority(fib, "sub_0.5", str(top_score)), 4)
        verdict["separates_cleanly"] = bool(verdict["clean_separation_fraction"] >= 0.75)
        p90 = float(np.quantile(fib.obs["doublet_score"], 0.90))
        verdict["doublet_guard"] = {
            "subcluster_median_scrublet": float(row["median_doublet_score"]),
            "subset_p90_scrublet": round(p90, 4),
            "flagged_by_scrublet": bool(row["median_doublet_score"] > p90),
            "flagged_by_mural_markers": bool(max(
                row.get("det_pericyte_max", 0) or 0, row.get("det_smooth_muscle_max", 0) or 0) > 0.50),
        }
    rec.set("recovery_verdict", verdict)

    # ---- over-clustering guard at all three resolutions ------------------
    guard = {res: unsupported_subclusters(fib, f"sub_{res}") for res in (0.2, 0.5, 1.0)}
    rec.set("over_clustering_guard", guard)
    pd.DataFrame([{"resolution": res, "subcluster": k, **v} for res, d in guard.items()
                  for k, v in d.items()]).to_csv(OUT / "c1_subcluster_support.csv", index=False)
    rec.add_output(OUT / "c1_subcluster_support.csv")

    # ---- one state or a gradient ----------------------------------------
    from sklearn.mixture import GaussianMixture
    red = fib.obs["genotype"].to_numpy() == "Red2Kras"
    values = fib.obs.loc[red, "score_reprogrammed fibroblast"].to_numpy().reshape(-1, 1)
    bic = {}
    for k in (1, 2):
        gm = GaussianMixture(n_components=k, random_state=0, n_init=3).fit(values)
        bic[k] = float(gm.bic(values))
    if bic[2] < bic[1] - 10:
        shape = "two states"
    elif bic[1] < bic[2]:
        shape = "gradient"
    else:
        shape = "undecided"
    sc.tl.diffmap(fib, random_state=0)
    alveolar = sub.loc[sub["call"] == "alveolar fibroblast", "cluster"].tolist()
    if alveolar:
        anchor = fib.obs["sub_0.5"].astype(str).isin(alveolar).to_numpy()
        root = int(np.argmin(np.linalg.norm(
            fib.obsm["X_diffmap"][:, 1:4] - fib.obsm["X_diffmap"][anchor, 1:4].mean(axis=0), axis=1)))
        fib.uns["iroot"] = root
        sc.tl.dpt(fib)
        dpt_by_sub = fib.obs.groupby("sub_0.5", observed=True)["dpt_pseudotime"].median().round(4)
        rec.set("dpt_median_by_subcluster", {str(k): float(v) for k, v in dpt_by_sub.items()})
    rec.set("state_or_gradient", {"bic_1": round(bic[1], 2), "bic_2": round(bic[2], 2),
                                  "verdict": shape,
                                  "n_red2kras_fibroblasts": int(red.sum())})

    # ---- Harmony sensitivity --------------------------------------------
    fib_h = ad.AnnData(X=fib.layers["counts"].copy(), obs=fib.obs[["library", "genotype"]].copy(),
                       var=fib.var.copy())
    fib_h = embed(fib_h, [0.5], key_prefix="harm", harmony=True)
    cross = pd.crosstab(fib.obs["sub_0.5"].astype(str).to_numpy(),
                        fib_h.obs["harm_0.5"].astype(str).to_numpy())
    cross.to_csv(OUT / "c1_harmony_sensitivity_crosstab.csv")
    rec.add_output(OUT / "c1_harmony_sensitivity_crosstab.csv")
    from sklearn.metrics import adjusted_rand_score
    rec.set("harmony_sensitivity", {
        "adjusted_rand_index_primary_vs_harmony": round(float(adjusted_rand_score(
            fib.obs["sub_0.5"].astype(str), fib_h.obs["harm_0.5"].astype(str))), 4),
        "n_subclusters_primary": int(fib.obs["sub_0.5"].nunique()),
        "n_subclusters_harmony": int(fib_h.obs["harm_0.5"].nunique()),
        "note": "reported as a sensitivity only; library and genotype are one variable"})

    # ---- figures ---------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    sc.pl.umap(merged, color="leiden_0.5", ax=axes[0][0], show=False, legend_loc="on data",
               size=12, title=f"mesenchyme, Leiden 0.5 (n={merged.n_obs})")
    sc.pl.umap(merged, color="genotype", ax=axes[0][1], show=False, size=12, title="genotype (= library)")
    sc.pl.umap(fib, color="sub_0.5", ax=axes[1][0], show=False, legend_loc="on data", size=14,
               title=f"fibroblasts re-embedded (n={fib.n_obs})")
    sc.pl.umap(fib, color="score_reprogrammed fibroblast", ax=axes[1][1], show=False, size=14,
               title="reprogrammed-fibroblast score")
    fig.suptitle(f"GSE316241 Gate 1: reprogrammed fibroblast {'recovered' if recovered else 'NOT recovered'} "
                 f"by the frozen rule; shape: {shape}")
    fig.tight_layout()
    fig.savefig(OUT / "c1_umap_mesenchyme_and_fibroblasts.png", dpi=110)
    plt.close(fig)
    rec.add_output(OUT / "c1_umap_mesenchyme_and_fibroblasts.png")

    lines = [
        "# Trial C1 output: the fibroblast compartment of GSE316241", "",
        f"**Verdict by the frozen rule: the reprogrammed fibroblast state is "
        f"{'RECOVERED' if recovered else 'NOT recovered'}.** "
        f"Shape of the Red2Kras fibroblast score: **{shape}** "
        f"(BIC 1 component {bic[1]:.1f}, 2 components {bic[2]:.1f}).", "",
        "Design constraint carried from trial C0: one library per genotype, three mice pooled per",
        "library, one time point. Every number below describes two libraries, not two groups of",
        "animals; no P value is computed for the genotype contrast.", "",
        "## QC per library", "", df_to_markdown(qc.drop(columns=["thresholds"]), index=False), "",
        "## Mesenchymal clusters (Leiden 0.5, labels from the paper's marker sets)", "",
        df_to_markdown(calls, index=False), "",
        "## Fibroblast subclusters", "", df_to_markdown(sub, index=False), "",
    ]
    (OUT / "c1_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c1_summary.md")
    rec.finish()
    print(f"\nverdict: recovered={recovered}, shape={shape}")


if __name__ == "__main__":
    main()
