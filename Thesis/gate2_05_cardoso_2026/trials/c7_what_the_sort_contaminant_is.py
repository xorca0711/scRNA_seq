#!/usr/bin/env python
"""Trial C7: which mutant epithelial state is the contaminant in the mesenchymal sort.

Item A1 of the next-step list in DIVERGENCES_AND_NEXT.md. Trial C1d found 184
cells in the CD45-CD31-EpCAM-negative mesenchymal sort of GSE316241 that read
as epithelial, 99 per cent of them from the Red2Kras library, with Areg
detected in 88 per cent. Trial C5 ruled out Epcam loss at the RNA level and
ruled out ambient RNA, leaving ordinary sort impurity with a doublet minority.
What none of that settled is which epithelial state those cells are.

It matters for how sharp the caution is. If they match the DATP-like state,
the contaminant is the very population the paper's signalling analysis is
about, and any ligand-receptor inference run inside GSE316241 alone would read
the paper's own sending cell as a mesenchymal source. If they match ordinary
AT2 or AT1-like cells, the leak is generic epithelium and the caution is
milder.

GSE316241 and GSE316244 share a gene space exactly (trial C0), so the
comparison needs no intersection and no batch correction between them is
implied: the profiles are compared by rank correlation, not merged.

Frozen rules, set before any correlation was computed:

* Query: cells of cluster 11 at leiden_0.5 in the cached mesenchyme object of
  trial C1b. Cluster 14, the cleanest fibroblast cluster in the same object,
  is carried as a negative control.
* Reference: the RFP-sorted epithelial libraries of GSE316244, both genotypes
  (Expt3_Het_RFP and Expt3_Hom_RFP), processed exactly as trial C2 processed
  them: per-library MAD quality control, Scrublet per capture with the 10x
  prior sanity check, seurat_v3 highly variable genes, 30 principal
  components, k = 30 neighbours, Leiden at resolution 0.5, seed 0.
* Reference states are named by the paper's own Figure 4l marker sets (AT2,
  DATP_like, AT1_like, Cd177_positive, cycling) scored with scanpy's
  score_genes and assigned by modal argmax, labels never used to fit anything.
  The 50 per cent confidence floor is reported and is NOT used to exclude a
  cluster, because trials C1 and C2 both showed that floor discarding the very
  populations under test.
* Comparison genes: genes detected in at least 10 per cent of the query
  cluster or of at least one reference cluster, computed on log-normalised
  values. Mitochondrial, ribosomal and haemoglobin genes are excluded, because
  they track quality rather than identity.
* T1, THE TEST. Spearman correlation of the query cluster's mean profile
  against each reference cluster's mean profile. The reference cluster with
  the highest correlation names the contaminant, and the margin to the second
  highest is reported.
* T2. The query cells are scored with the same Figure 4l marker sets directly,
  and the modal call is reported beside T1. Two independent routes to the same
  question.
* T3, THE CONTROL THAT DECIDES WHETHER T1 MEANS ANYTHING. The same correlation
  for cluster 14. RULE: cluster 11's best correlation must exceed cluster 14's
  best correlation by at least 0.10 of Spearman rho. If it does not, the
  measure is reading depth or ambient signal rather than identity, and T1 and
  T2 are reported without a reading.
* THE READING, fixed in advance:
  - top state DATP_like: the contaminant is the paper's own signalling
    population, and the sort caution of finding one is sharper;
  - top state AT2 or AT1_like: the leak is ordinary alveolar epithelium and
    the caution is milder;
  - top state Cd177_positive or cycling: an unexpected state, reported as
    such with no further reading;
  - T3 failed: no reading at all.
* Caveats: 184 query cells, so this is a profile comparison and not a test of
  a difference between groups; about one cell in five co-detects Krt8 and
  Col1a1 (trial C5), so a doublet minority contributes to the query profile.
"""

from __future__ import annotations

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

OUT = HERE / "c7_what_the_sort_contaminant_is"
OUT.mkdir(exist_ok=True)
MESENCHYME = HERE / "c1b_characterise_red2kras_private" / "c1b_mesenchyme.h5ad"
ROOT = RAW / "GSE316244" / "GSE316244_RAW"
REFERENCE_LIBRARIES = [
    {"gsm": "GSM9447781", "library": "Expt3_Het_RFP", "genotype": "Areg-flox/+"},
    {"gsm": "GSM9447782", "library": "Expt3_Hom_RFP", "genotype": "Areg-flox/flox"},
]
STATES = {
    "AT2": ["Sftpc", "Etv5", "Lamp3"],
    "DATP_like": ["Cldn4", "Itga2", "Ndrg1", "Sox9"],
    "AT1_like": ["Ager", "Hopx", "Clic5"],
    "Cd177_positive": ["Cd177", "Cd38", "Dlk1"],
    "cycling": ["Mki67", "Birc5"],
}
QUERY_CLUSTER = "11"
CONTROL_CLUSTER = "14"
DETECTION_FLOOR = 0.10
MARGIN = 0.10

RULES = {
    "question": "item A1: which mutant epithelial state is the contaminant of the mesenchymal sort",
    "query": f"cluster {QUERY_CLUSTER} at leiden_0.5 in the C1b mesenchyme object",
    "control": f"cluster {CONTROL_CLUSTER}, the cleanest fibroblast cluster in the same object",
    "reference": "RFP-sorted epithelium of GSE316244, both genotypes, processed as in trial C2",
    "gene_space": "GSE316241 and GSE316244 share a gene space exactly (C0), so no intersection is needed",
    "state_markers": STATES,
    "confidence_floor": "reported, never used to exclude a cluster; C1 and C2 showed that floor discarding the populations under test",
    "comparison_genes": f"detected in at least {DETECTION_FLOOR:.0%} of the query or of one reference cluster; mt, ribo and hb genes excluded",
    "T1": "Spearman of the query mean profile against each reference cluster mean profile; highest names the state",
    "T2": "the query cells scored directly with the same marker sets, modal call",
    "T3": f"control: cluster {QUERY_CLUSTER} must beat cluster {CONTROL_CLUSTER} by at least {MARGIN} rho, or T1 and T2 get no reading",
    "reading_fixed_in_advance": {
        "DATP_like": "the contaminant is the paper's own signalling population; the finding-one caution is sharper",
        "AT2_or_AT1_like": "the leak is ordinary alveolar epithelium; the caution is milder",
        "Cd177_positive_or_cycling": "an unexpected state, reported with no further reading",
        "control_failed": "the measure reads depth or ambient signal, not identity; no reading",
    },
    "caveats": ["184 query cells; a profile comparison, not a test between groups",
                "about one query cell in five co-detects Krt8 and Col1a1 (C5), a doublet minority"],
}


def load_reference(rec):
    import anndata as ad
    import scanpy as sc

    parts = []
    for entry in REFERENCE_LIBRARIES:
        stem = entry["gsm"] + "_" + entry["library"] + "_"
        matrix = ROOT / (stem + "matrix.mtx.gz")
        features = ROOT / (stem + "features.tsv.gz")
        if not features.exists():
            features = ROOT / (stem + "genes.tsv.gz")
        barcodes = ROOT / (stem + "barcodes.tsv.gz")
        for path in (matrix, features, barcodes):
            rec.add_input(path)
        X, var, bc = read_mtx_triplet(matrix, features, barcodes)
        keep = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
        if keep.sum() == 0:
            keep = np.ones(len(var), dtype=bool)
        X, var = X[:, keep], var.loc[keep].reset_index(drop=True)
        var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
        obs = pd.DataFrame({"library": entry["library"], "genotype": entry["genotype"]},
                           index=pd.Index([entry["library"] + "_" + b for b in bc]))
        adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)
        symbols = adata.var["gene_symbol"].astype(str)
        adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True)
        th = derive_thresholds(adata.obs, entry["library"])
        adata = adata[apply_thresholds(adata.obs, th).to_numpy()].copy()
        rate = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
        sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=rate)
        auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
        if auto < 0.2 * rate or auto > 3 * rate:
            cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - rate))
            adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
        parts.append(adata[~adata.obs["predicted_doublet"].to_numpy()].copy())
        print(entry["library"], parts[-1].n_obs, "cells kept")

    ref = ad.concat(parts, join="inner", index_unique=None)
    ref.var = parts[0].var.loc[ref.var_names].copy()
    ref.layers["counts"] = ref.X.copy()
    import scanpy as sc
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)
    sc.pp.highly_variable_genes(ref, n_top_genes=2000, flavor="seurat_v3", layer="counts")
    sc.pp.pca(ref, n_comps=30, svd_solver="arpack", random_state=0)
    sc.pp.neighbors(ref, n_neighbors=30, random_state=0)
    sc.tl.leiden(ref, resolution=0.5, key_added="leiden", flavor="igraph",
                 n_iterations=2, directed=False, random_state=0)
    for state, genes in STATES.items():
        present = [g for g in genes if g in ref.var_names]
        sc.tl.score_genes(ref, present, score_name="score_" + state, random_state=0)
    return ref


def mean_profile(adata, mask, genes):
    import scipy.sparse as sp

    sub = adata[mask, genes].X
    if sp.issparse(sub):
        return np.asarray(sub.mean(axis=0)).ravel()
    return np.asarray(sub).mean(axis=0)


def detected(adata, mask, genes):
    import scipy.sparse as sp

    sub = adata[mask, genes].X
    if sp.issparse(sub):
        return np.asarray((sub > 0).mean(axis=0)).ravel()
    return (np.asarray(sub) > 0).mean(axis=0)


def main():
    import anndata as ad
    import scanpy as sc
    from scipy.stats import spearmanr

    rec = RunRecord(OUT / "c7_run_record.json", "C7 identity of the mesenchymal-sort contaminant", RULES)
    rec.add_input(MESENCHYME)
    query_object = ad.read_h5ad(MESENCHYME)
    clusters = query_object.obs["leiden_0.5"].astype(str).to_numpy()
    ref = load_reference(rec)
    rec.set("reference_cells", int(ref.n_obs))
    rec.set("reference_clusters", int(ref.obs["leiden"].nunique()))

    shared = [g for g in ref.var_names if g in set(query_object.var_names)]
    drop = set(ref.var_names[ref.var["gene_symbol"].astype(str).str.startswith(("mt-", "Rps", "Rpl", "Hba", "Hbb"))])
    shared = [g for g in shared if g not in drop]
    rec.set("shared_genes_considered", len(shared))

    labels, calls = {}, []
    score_columns = ["score_" + s for s in STATES]
    for cluster in sorted(ref.obs["leiden"].unique(), key=lambda c: int(c)):
        mask = (ref.obs["leiden"] == cluster).to_numpy()
        means = ref.obs.loc[mask, score_columns].mean()
        per_cell = ref.obs.loc[mask, score_columns].to_numpy().argmax(axis=1)
        mode = np.bincount(per_cell, minlength=len(score_columns)).argmax()
        fraction = float((per_cell == mode).mean())
        name = list(STATES)[mode]
        labels[cluster] = name
        calls.append({"cluster": cluster, "n_cells": int(mask.sum()), "call": name,
                      "mode_fraction": round(fraction, 3), "confident": bool(fraction >= 0.5),
                      "top_score": round(float(means.max()), 4)})
    call_frame = pd.DataFrame(calls)
    call_frame.to_csv(OUT / "c7_reference_clusters.csv", index=False)
    rec.add_output(OUT / "c7_reference_clusters.csv")

    rows = []
    for which, cluster_id in (("query cluster " + QUERY_CLUSTER, QUERY_CLUSTER),
                              ("control cluster " + CONTROL_CLUSTER, CONTROL_CLUSTER)):
        qmask = clusters == cluster_id
        if qmask.sum() == 0:
            raise SystemExit(f"cluster {cluster_id} is absent from the cached object")
        qdet = detected(query_object, qmask, shared)
        qprofile = mean_profile(query_object, qmask, shared)
        for cluster in call_frame["cluster"]:
            rmask = (ref.obs["leiden"] == cluster).to_numpy()
            rdet = detected(ref, rmask, shared)
            keep = (qdet >= DETECTION_FLOOR) | (rdet >= DETECTION_FLOOR)
            if keep.sum() < 200:
                continue
            rho = spearmanr(qprofile[keep], mean_profile(ref, rmask, shared)[keep]).statistic
            rows.append({"query": which, "n_query_cells": int(qmask.sum()),
                         "reference_cluster": cluster, "reference_call": labels[cluster],
                         "n_reference_cells": int(rmask.sum()),
                         "genes_compared": int(keep.sum()), "rho": round(float(rho), 4)})
    corr = pd.DataFrame(rows).sort_values(["query", "rho"], ascending=[True, False])
    corr.to_csv(OUT / "c7_profile_correlations.csv", index=False)
    rec.add_output(OUT / "c7_profile_correlations.csv")

    best = {}
    for which, group in corr.groupby("query"):
        ordered = group.sort_values("rho", ascending=False).reset_index(drop=True)
        best[which] = {"top_call": ordered.loc[0, "reference_call"],
                       "top_cluster": ordered.loc[0, "reference_cluster"],
                       "top_rho": float(ordered.loc[0, "rho"]),
                       "second_call": ordered.loc[1, "reference_call"] if len(ordered) > 1 else None,
                       "margin_to_second": round(float(ordered.loc[0, "rho"] - ordered.loc[1, "rho"]), 4)
                       if len(ordered) > 1 else None}
    rec.set("T1_best_match", best)

    qmask = clusters == QUERY_CLUSTER
    scored = query_object[qmask].copy()
    if "log1p" not in query_object.uns:
        raise SystemExit("the cached object is not log-normalised; scoring it would not be comparable")
    for state, genes in STATES.items():
        present = [g for g in genes if g in scored.var_names]
        sc.tl.score_genes(scored, present, score_name="score_" + state, random_state=0)
    per_cell = scored.obs[score_columns].to_numpy().argmax(axis=1)
    counts = np.bincount(per_cell, minlength=len(score_columns))
    t2 = {list(STATES)[i]: int(counts[i]) for i in range(len(STATES))}
    t2_call = list(STATES)[int(counts.argmax())]
    rec.set("T2_direct_scoring", {"counts": t2, "modal_call": t2_call,
                                  "mode_fraction": round(float(counts.max() / counts.sum()), 3)})

    qbest = best["query cluster " + QUERY_CLUSTER]["top_rho"]
    cbest = best["control cluster " + CONTROL_CLUSTER]["top_rho"]
    control_passes = bool(qbest - cbest >= MARGIN)
    rec.set("T3_control", {"query_best_rho": round(qbest, 4), "control_best_rho": round(cbest, 4),
                           "margin": round(qbest - cbest, 4), "required": MARGIN,
                           "passes": control_passes})

    top = best["query cluster " + QUERY_CLUSTER]["top_call"]
    if not control_passes:
        reading = RULES["reading_fixed_in_advance"]["control_failed"]
    elif top == "DATP_like":
        reading = RULES["reading_fixed_in_advance"]["DATP_like"]
    elif top in ("AT2", "AT1_like"):
        reading = RULES["reading_fixed_in_advance"]["AT2_or_AT1_like"]
    else:
        reading = RULES["reading_fixed_in_advance"]["Cd177_positive_or_cycling"]
    if control_passes and t2_call != top:
        reading += (" The two routes disagree: profile correlation says " + top
                    + " and direct scoring says " + t2_call + ", so neither is read as settled.")
    rec.set("reading", reading)

    lines = ["# Trial C7: what the mesenchymal-sort contaminant is", "",
             "T1 best match: " + top + " (rho " + f"{qbest:.4f}" + ", margin to second "
             + str(best["query cluster " + QUERY_CLUSTER]["margin_to_second"]) + ").",
             "T2 direct scoring modal call: " + t2_call + " " + str(t2) + ".",
             "T3 control: cluster " + QUERY_CLUSTER + " best " + f"{qbest:.4f}" + " against cluster "
             + CONTROL_CLUSTER + " best " + f"{cbest:.4f}" + ", "
             + ("passes" if control_passes else "FAILS") + " the " + str(MARGIN) + " margin.", "",
             "Reading: " + reading, "",
             "## Reference clusters of the RFP-sorted epithelium", "",
             df_to_markdown(call_frame, index=False), "",
             "## Profile correlations", "", df_to_markdown(corr, index=False), ""]
    (OUT / "c7_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c7_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
