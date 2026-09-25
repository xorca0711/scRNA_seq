#!/usr/bin/env python
"""Trial C1b: characterise what GSE316241 has instead (Gate 1, stop-and-characterise).

Trial C1 applied its frozen rule and returned **not recovered**: no fibroblast
subcluster satisfied all three criteria. The owner's Gate 1 says that when the
published subset does not separate, or splits differently, that is the finding
and the trial stops to characterise it rather than pressing on. This is that
characterisation.

It also discloses a rule that defeated itself, in the manner of this
repository's decision 13. C1's frozen selection took "clusters whose
**confident** call is alveolar, adventitial or reprogrammed fibroblast". In
the event, every cluster the marker-score caller labelled *reprogrammed
fibroblast* (clusters 8, 12 and 15) fell below the 50% confidence floor, and
so did cluster 14, which is 99.6% Red2Kras. The selection therefore removed
exactly the candidates the test was about, and the test was then run on a
subset that could not contain the answer. C1's outcome stands as recorded;
this trial adds the corrected pass beside it and reports both.

**Provenance of these rules.** They were fixed after C1's cluster-level table
(cluster sizes, genotype shares and calls) had been seen, and before any
per-cell data was read in this trial. That is the same disclosure trial S4
made, and it is stated here because the rules below are not blind to C1's
output.

Frozen rules:

* The merged mesenchymal object is rebuilt by the identical code path as C1
  (same QC, same Scrublet handling, same features, same seeds), so the cluster
  numbering matches C1's. It is written to disk so later trials need not
  recompute it.
* Every Leiden 0.5 cluster is characterised, not only the selected ones:
  size, genotype share, median counts, genes, mitochondrial percentage and
  Scrublet score, mean score for each of the paper's population sets, and the
  detection fraction of every gene in those sets.
* QUESTION A, the direct one. Does any cluster match the published
  reprogrammed fibroblast? A cluster qualifies if it is at least 80%
  Red2Kras and Tnc, Acta2 and (Pdgfrb or Runx1) are each detected in at least
  40% of its cells. Every cluster's values are reported whether or not it
  qualifies.
* QUESTION B, the Red2Kras-private structure. Any cluster at least 90%
  Red2Kras is characterised in full, and each is classified by the same
  criteria trial S4 used for an unexplained cluster: "low-count, ambient-like"
  if its median counts and genes are both below half the object's medians;
  "doublet-enriched" if its median Scrublet score exceeds the object's 90th
  percentile; otherwise "not explained by quality".
* QUESTION C, the corrected selection, disclosed as post hoc. Fibroblasts are
  taken by a compartment gate on the raw counts instead of by the argmax call:
  Col1a1 detected, and Ptprc, Pecam1 and Epcam not detected. The subset is
  re-embedded with C1's parameters and C1's recovery test is re-applied
  unchanged. Both outcomes are reported; neither replaces the other.
* Top genes per cluster are a Wilcoxon ranking device for interpretation, as
  everywhere else in this repository, never a test between genotypes.
* Unit: the library. No P value.
"""

from __future__ import annotations

import gc
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import (ENSEMBL_ID_RE, RAW, REPO, RunRecord,  # noqa: E402
                           df_to_markdown, read_mtx_triplet)

sys.path.insert(0, str(REPO / "analysis" / "scripts"))
from pipeline_utils import apply_thresholds, derive_thresholds, make_unique  # noqa: E402

OUT = HERE / "c1b_characterise_red2kras_private"
OUT.mkdir(exist_ok=True)
OBJECT = OUT / "c1b_mesenchyme.h5ad"
EXTRACTS = json.loads((HERE.parent / "cardoso_2026_extracts.json").read_text(encoding="utf-8"))
SETS = EXTRACTS["marker_sets"]

LIBRARIES = [
    {"library": "Expt1_Confetti_mesenchyme", "gsm": "GSM9447763", "genotype": "Confetti"},
    {"library": "Expt1_Red2Kras_mesenchyme", "gsm": "GSM9447764", "genotype": "Red2Kras"},
]
POPULATION_SETS = {
    "alveolar fibroblast": SETS["alveolar_fibroblast"],
    "adventitial fibroblast": SETS["adventitial_fibroblast"],
    "reprogrammed fibroblast": SETS["reprogrammed_fibrotic_fibroblast"],
    "inflammatory fibroblast": SETS["inflammatory_fibroblast"],
    "peri-bronchial fibroblast": SETS["peri_bronchial_fibroblast"],
    "smooth muscle": SETS["smooth_muscle"],
    "pericyte": SETS["pericyte"],
    "mesothelium": SETS["mesothelium"],
    "mesothelial-like": SETS["mesothelial_like_enriched_in_tumour"],
    "proliferating": SETS["proliferating"],
}
QUESTION_A = {"red2kras_fraction_min": 0.80, "required_detection": 0.40,
              "required_genes": ["Tnc", "Acta2"], "either_of": ["Pdgfrb", "Runx1"]}
PRIVATE_MIN = 0.90
GATE = {"positive": ["Col1a1"], "negative": ["Ptprc", "Pecam1", "Epcam"]}

RULES = {
    "accession": "GSE316241",
    "provenance_of_these_rules": ("fixed after C1's cluster-level table was seen and before any "
                                  "per-cell data was read in this trial; the same disclosure trial "
                                  "S4 made"),
    "c1_rule_that_defeated_itself": ("C1 selected fibroblasts by confident call; every cluster called "
                                     "reprogrammed fibroblast fell below the 50% confidence floor, so "
                                     "the selection removed the candidates the test was about"),
    "rebuild": "identical code path to C1 (QC, Scrublet, features, seeds), so cluster numbering matches",
    "question_A": QUESTION_A,
    "question_B": {"private_if_red2kras_fraction_at_least": PRIVATE_MIN,
                   "low_count_ambient_like": "median counts and genes both below half the object medians",
                   "doublet_enriched": "median Scrublet score above the object's 90th percentile"},
    "question_C": {"gate": GATE, "then": "re-embed with C1 parameters and re-apply C1's recovery test unchanged",
                   "status": "post hoc, disclosed; does not replace C1's outcome"},
    "top_genes": "Wilcoxon, a ranking device for interpretation only",
    "unit": "the library; no P value",
}


def build(rec: RunRecord):
    import anndata as ad
    import scanpy as sc

    if OBJECT.exists():
        print(f"  reusing {OBJECT.name}")
        return ad.read_h5ad(OBJECT)

    parts = []
    for entry in LIBRARIES:
        root = RAW / "GSE316241" / "GSE316241_RAW"
        stem = f"{entry['gsm']}_{entry['library']}_"
        paths = {k: root / f"{stem}{k}" for k in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")}
        for path in paths.values():
            rec.add_input(path)
        X, var, bc = read_mtx_triplet(paths["matrix.mtx.gz"], paths["features.tsv.gz"],
                                      paths["barcodes.tsv.gz"])
        is_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
        X, var = X[:, is_gene], var.loc[is_gene].reset_index(drop=True)
        var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
        obs = pd.DataFrame(index=pd.Index([f"{entry['library']}_{b}" for b in bc]))
        obs["library"], obs["genotype"] = entry["library"], entry["genotype"]
        adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)

        symbols = adata.var["gene_symbol"].astype(str)
        adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
        adata.var["ribo"] = (symbols.str.startswith(("Rps", "Rpl"))
                             & ~symbols.str.contains("-ps", case=False)).to_numpy()
        adata.var["hb"] = symbols.str.startswith(("Hba", "Hbb")).to_numpy()
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20],
                                   log1p=True, inplace=True)
        th = derive_thresholds(adata.obs, entry["library"])
        adata = adata[apply_thresholds(adata.obs, th).to_numpy()].copy()
        expected = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
        sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=expected)
        auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
        if auto < 0.2 * expected or auto > 3 * expected:
            cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - expected))
            adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
        adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
        parts.append(adata)
        print(f"  {entry['library']}: {adata.n_obs} cells")

    merged = ad.concat(parts, join="inner")
    merged.var = parts[0].var.loc[merged.var_names].copy()
    del parts
    gc.collect()

    merged.layers["counts"] = merged.X.copy()
    sc.pp.normalize_total(merged, target_sum=1e4)
    sc.pp.log1p(merged)
    sc.pp.highly_variable_genes(merged, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = merged[:, merged.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=30, random_state=0)
    merged.obsm["X_pca"] = emb.obsm["X_pca"]
    del emb
    gc.collect()
    sc.pp.neighbors(merged, n_neighbors=30, n_pcs=30, use_rep="X_pca", random_state=0)
    for res in (0.2, 0.5, 1.0):
        sc.tl.leiden(merged, resolution=res, key_added=f"leiden_{res}", flavor="igraph",
                     n_iterations=2, random_state=0, directed=False)
    sc.tl.umap(merged, random_state=0)
    for name, genes in POPULATION_SETS.items():
        present = [g for g in genes if g in merged.var_names]
        if present:
            sc.tl.score_genes(merged, present, score_name=f"score_{name}", ctrl_size=50,
                              n_bins=25, random_state=0)
    merged.write_h5ad(OBJECT, compression="gzip")
    return merged


def detection_matrix(adata, genes: list[str], key: str) -> pd.DataFrame:
    X = adata.layers["counts"]
    groups = adata.obs[key].astype(str).to_numpy()
    order = sorted(set(groups), key=int)
    data = {}
    for gene in genes:
        if gene not in adata.var_names:
            continue
        col = np.asarray(X[:, adata.var_names.get_loc(gene)].todense()).ravel() > 0
        data[gene] = [round(float(col[groups == g].mean()), 3) for g in order]
    return pd.DataFrame(data, index=pd.Index(order, name="cluster"))


def main() -> None:
    import anndata as ad
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc.settings.verbosity = 0
    rec = RunRecord(OUT / "c1b_run_record.json",
                    "C1b characterise the Red2Kras-private mesenchymal structure (Gate 1)",
                    RULES, notes="rules fixed after C1's cluster table was seen; stated in the rules")
    merged = build(rec)
    key = "leiden_0.5"
    rec.set("cells", int(merged.n_obs))

    # ---- per-cluster characterisation ------------------------------------
    med_counts = float(np.median(merged.obs["total_counts"]))
    med_genes = float(np.median(merged.obs["n_genes_by_counts"]))
    p90_doublet = float(np.quantile(merged.obs["doublet_score"], 0.90))
    rows = []
    for cluster in sorted(merged.obs[key].astype(str).unique(), key=int):
        mask = merged.obs[key].astype(str).to_numpy() == cluster
        obs = merged.obs.loc[mask]
        row = {"cluster": cluster, "n_cells": int(mask.sum()),
               "fraction_Red2Kras": round(float((obs["genotype"] == "Red2Kras").mean()), 3),
               "median_counts": int(np.median(obs["total_counts"])),
               "median_genes": int(np.median(obs["n_genes_by_counts"])),
               "median_pct_mt": round(float(np.median(obs["pct_counts_mt"])), 2),
               "median_doublet_score": round(float(np.median(obs["doublet_score"])), 4)}
        for name in POPULATION_SETS:
            if f"score_{name}" in merged.obs:
                row[f"score_{name}"] = round(float(obs[f"score_{name}"].mean()), 4)
        rows.append(row)
    chars = pd.DataFrame(rows)

    panel_genes = sorted({g for genes in POPULATION_SETS.values() for g in genes})
    det = detection_matrix(merged, panel_genes, key)
    det.to_csv(OUT / "c1b_marker_detection_by_cluster.csv")
    rec.add_output(OUT / "c1b_marker_detection_by_cluster.csv")

    # ---- Question A ------------------------------------------------------
    qa = QUESTION_A
    a_rows = []
    for _, row in chars.iterrows():
        cluster = row["cluster"]
        values = {g: float(det.loc[cluster, g]) for g in qa["required_genes"] + qa["either_of"]
                  if g in det.columns}
        required_ok = all(values.get(g, 0) >= qa["required_detection"] for g in qa["required_genes"])
        either_ok = any(values.get(g, 0) >= qa["required_detection"] for g in qa["either_of"])
        private_ok = row["fraction_Red2Kras"] >= qa["red2kras_fraction_min"]
        a_rows.append({"cluster": cluster, "n_cells": row["n_cells"],
                       "fraction_Red2Kras": row["fraction_Red2Kras"],
                       **{f"det_{g}": values.get(g) for g in qa["required_genes"] + qa["either_of"]},
                       "qualifies": bool(required_ok and either_ok and private_ok)})
    qa_frame = pd.DataFrame(a_rows)
    qa_frame.to_csv(OUT / "c1b_question_a_published_signature.csv", index=False)
    rec.add_output(OUT / "c1b_question_a_published_signature.csv")
    rec.set("question_A", {"any_cluster_matches_published_signature": bool(qa_frame["qualifies"].any()),
                           "qualifying_clusters": qa_frame.loc[qa_frame["qualifies"], "cluster"].tolist()})

    # ---- Question B ------------------------------------------------------
    private = chars[chars["fraction_Red2Kras"] >= PRIVATE_MIN].copy()
    verdicts = []
    for _, row in private.iterrows():
        if row["median_counts"] < 0.5 * med_counts and row["median_genes"] < 0.5 * med_genes:
            verdicts.append("low-count, ambient-like")
        elif row["median_doublet_score"] > p90_doublet:
            verdicts.append("doublet-enriched")
        else:
            verdicts.append("not explained by quality")
    private["quality_verdict"] = verdicts
    private.to_csv(OUT / "c1b_question_b_red2kras_private.csv", index=False)
    rec.add_output(OUT / "c1b_question_b_red2kras_private.csv")
    rec.set("question_B", {"object_median_counts": med_counts, "object_median_genes": med_genes,
                           "doublet_p90": round(p90_doublet, 4),
                           "private_clusters": private[["cluster", "n_cells", "fraction_Red2Kras",
                                                        "quality_verdict"]].to_dict("records")})

    sc.tl.rank_genes_groups(merged, key, method="wilcoxon", n_genes=12)
    top = pd.DataFrame({g: list(merged.uns["rank_genes_groups"]["names"][g])
                        for g in merged.obs[key].cat.categories})
    top.to_csv(OUT / "c1b_top_genes_by_cluster.csv", index=False)
    rec.add_output(OUT / "c1b_top_genes_by_cluster.csv")

    chars.to_csv(OUT / "c1b_cluster_characterisation.csv", index=False)
    rec.add_output(OUT / "c1b_cluster_characterisation.csv")

    # ---- Question C: the corrected selection -----------------------------
    counts = merged.layers["counts"]
    def detected(gene):
        return np.asarray(counts[:, merged.var_names.get_loc(gene)].todense()).ravel() > 0
    gate = detected(GATE["positive"][0])
    for gene in GATE["negative"]:
        if gene in merged.var_names:
            gate &= ~detected(gene)
    rec.set("question_C_gate", {"cells_passing": int(gate.sum()),
                                "fraction_of_object": round(float(gate.mean()), 4)})
    fib = ad.AnnData(X=counts[gate].copy(),
                     obs=merged.obs.loc[gate, ["library", "genotype", "doublet_score"]].copy(),
                     var=merged.var.copy())
    fib.layers["counts"] = fib.X.copy()
    sc.pp.normalize_total(fib, target_sum=1e4)
    sc.pp.log1p(fib)
    sc.pp.highly_variable_genes(fib, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = fib[:, fib.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=30, random_state=0)
    fib.obsm["X_pca"] = emb.obsm["X_pca"]
    del emb
    gc.collect()
    sc.pp.neighbors(fib, n_neighbors=30, n_pcs=30, use_rep="X_pca", random_state=0)
    sc.tl.leiden(fib, resolution=0.5, key_added="sub_0.5", flavor="igraph", n_iterations=2,
                 random_state=0, directed=False)
    sc.tl.umap(fib, random_state=0)
    for name, genes in POPULATION_SETS.items():
        present = [g for g in genes if g in fib.var_names]
        if present:
            sc.tl.score_genes(fib, present, score_name=f"score_{name}", ctrl_size=50,
                              n_bins=25, random_state=0)

    det_sub = detection_matrix(fib, ["Tnc", "Fst", "Runx1", "Runx2", "Acta2", "Pdgfrb", "Pdgfra",
                                     "Col13a1", "Col14a1", "Lcn2", "Saa3"], "sub_0.5")
    rows = []
    for cluster in sorted(fib.obs["sub_0.5"].astype(str).unique(), key=int):
        mask = fib.obs["sub_0.5"].astype(str).to_numpy() == cluster
        row = {"cluster": cluster, "n_cells": int(mask.sum()),
               "fraction_Red2Kras": round(float((fib.obs.loc[mask, "genotype"] == "Red2Kras").mean()), 3),
               "mean_reprogrammed_score": round(float(fib.obs.loc[mask, "score_reprogrammed fibroblast"].mean()), 4),
               "median_doublet_score": round(float(np.median(fib.obs.loc[mask, "doublet_score"])), 4)}
        for gene in det_sub.columns:
            row[f"det_{gene}"] = float(det_sub.loc[cluster, gene])
        rows.append(row)
    sub = pd.DataFrame(rows)
    sub.to_csv(OUT / "c1b_question_c_gated_fibroblast_subclusters.csv", index=False)
    rec.add_output(OUT / "c1b_question_c_gated_fibroblast_subclusters.csv")

    top_score = sub.loc[sub["mean_reprogrammed_score"].idxmax()]
    retest = {
        "subcluster": str(top_score["cluster"]),
        "red2kras_fraction": float(top_score["fraction_Red2Kras"]),
        "red2kras_fraction_passes": bool(top_score["fraction_Red2Kras"] >= 0.80),
        "tnc_detection": float(top_score["det_Tnc"]),
        "tnc_detection_passes": bool(top_score["det_Tnc"] >= 0.40),
        "has_highest_reprogrammed_score": True,
    }
    retest["recovered"] = bool(retest["red2kras_fraction_passes"] and retest["tnc_detection_passes"])
    rec.set("question_C_retest", retest)
    print(df_to_markdown(sub, index=False))

    # ---- figure ----------------------------------------------------------
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    sc.pl.umap(merged, color=key, ax=axes[0][0], show=False, legend_loc="on data", size=10,
               title=f"mesenchyme, Leiden 0.5 (n={merged.n_obs})")
    sc.pl.umap(merged, color="genotype", ax=axes[0][1], show=False, size=10, title="genotype (= library)")
    sc.pl.umap(merged, color="Tnc", ax=axes[1][0], show=False, size=10, title="Tnc")
    sc.pl.umap(merged, color="score_reprogrammed fibroblast", ax=axes[1][1], show=False, size=10,
               title="reprogrammed-fibroblast score")
    fig.tight_layout()
    fig.savefig(OUT / "c1b_umap_characterisation.png", dpi=110)
    plt.close(fig)
    rec.add_output(OUT / "c1b_umap_characterisation.png")

    lines = [
        "# Trial C1b output: what GSE316241 has instead (Gate 1, stop-and-characterise)", "",
        "Trial C1's frozen rule returned **not recovered**. Its selection step took clusters whose",
        "*confident* call was a fibroblast identity, and every cluster the caller labelled",
        "reprogrammed fibroblast fell below the confidence floor, so the selection removed the",
        "candidates the test was about. C1's outcome stands as recorded; the corrected pass is",
        "Question C below and is post hoc.", "",
        "## Question A: does any cluster match the published signature", "",
        f"Qualifying clusters: **{qa_frame.loc[qa_frame['qualifies'], 'cluster'].tolist() or 'none'}** "
        f"(at least 80% Red2Kras, with Tnc, Acta2 and Pdgfrb or Runx1 each detected in at least 40% "
        f"of cells).", "", df_to_markdown(qa_frame, index=False), "",
        "## Question B: the Red2Kras-private clusters", "",
        df_to_markdown(private[["cluster", "n_cells", "fraction_Red2Kras", "median_counts",
                                "median_genes", "median_pct_mt", "median_doublet_score",
                                "quality_verdict"]], index=False), "",
        f"Object medians for reference: {int(med_counts)} counts, {int(med_genes)} genes; "
        f"Scrublet 90th percentile {p90_doublet:.4f}.", "",
        "## Question C: the corrected selection, post hoc", "",
        f"Compartment gate (Col1a1 detected; Ptprc, Pecam1, Epcam not detected) keeps "
        f"{int(gate.sum())} of {merged.n_obs} cells. Re-applying C1's recovery test unchanged: "
        f"**{'recovered' if retest['recovered'] else 'still not recovered'}**.", "",
        df_to_markdown(sub, index=False), "",
        "## All clusters, characterised", "", df_to_markdown(chars, index=False), "",
    ]
    (OUT / "c1b_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c1b_summary.md")
    rec.finish()
    print(f"\nQuestion A qualifying clusters: {qa_frame.loc[qa_frame['qualifies'], 'cluster'].tolist()}")
    print(f"Question C retest recovered: {retest['recovered']}")


if __name__ == "__main__":
    main()
