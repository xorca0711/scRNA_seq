#!/usr/bin/env python
"""Trial S2: map GSE178360 onto the HLCA core with scArches surgery and transfer labels with uncertainty.

Two stages, selected with --stage:

  surgery   prepare the query in the reference gene space, run scArches
            surgery on the HLCA scANVI model (scvi-tools implementation),
            train the query adaptor, save the query latent space.
  transfer  weighted k-nearest-neighbour label transfer from the HLCA core
            embedding with the paper's uncertainty score, then the
            pre-registered comparisons, figures and summary.

Portfolio purpose. The reference-based annotation of the human distal-lung
series by the route the HLCA itself used (position in the integrated
embedding plus label transfer), which is the only route that can settle the
AT0 question after trial S3 showed marker transfer is scheme-dependent.
Serves the AT2-to-AT0 axis (Nabhan fit) and donor-aware, reference-based
annotation (Wagner fit). No infection or interferon framing.

Frozen rules (from ANALYSIS_TRIAL_PLAN.md, section S2, written before this
script ran):

* Query: the post-QC, post-Scrublet matrix (27,729 cells), raw counts,
  genes matched to the model's 2,000 Ensembl IDs (version suffix stripped),
  missing genes zero-filled and counted; batch = one new dataset category
  for the whole series (the paper's batch unit is the dataset); labels =
  "unlabeled".
* Surgery: scvi-tools SCANVI.load_query_data with freeze_dropout True;
  max 500 epochs; weight decay 0; reduce LR on plateau (patience 8, factor
  0.1); early stopping on the validation ELBO, patience 10, min delta
  0.001, with a 90/10 train/validation split (deviation from the paper's
  "full dataset" monitoring, which scvi-tools no longer offers).
* Label transfer: k = 50 nearest reference cells in the joint 30-dim
  latent space; weights = Gaussian kernel of the neighbour distances with
  per-cell bandwidth (2 / std of the 50 distances)^2, as in scArches
  weighted_knn_transfer; transferred label = weighted majority;
  uncertainty u = 1 - weight of that label.
* Unknown: u > 0.2 primary; 0.3 reported as sensitivity.
* Comparisons per donor: transferred labels versus the blind proposals
  and versus the S3 assignments; transferred AT0 and pre-TB secretory
  counts versus the strict gate and versus both S3 schemes, concordant if
  within a factor of two per donor; per-cluster mean uncertainty with
  cluster 30 as a named check; epithelial subcluster 4 crosstab.
* Unit: donor. No P values.

Deviations recorded here rather than hidden: the scarches package
(0.6.1) cannot be imported with anndata 0.13, so the surgery uses
scvi-tools' own scArches implementation and the weighted kNN transfer is
reimplemented from the scArches formula above.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trial_utils import REPO, RunRecord, df_to_markdown, file_facts  # noqa: E402

OUT = HERE / "s2_reference_mapping"
REF = OUT / "reference"
MODEL_LEGACY = REF / "HLCA_reference_model"
MODEL = REF / "HLCA_reference_model_converted"
EMB = REF / "HLCA_full_v1.1_emb.h5ad"
GENE_ORDER = REF / "HLCA_reference_model_gene_order_ids_and_symbols.csv"
HUMAN = REPO / "analysis" / "GSE178360"
OBJ = HUMAN / "processed" / "final_clustered.h5ad"
EPI_OBJ = HUMAN / "epithelial_subanalysis" / "epithelial_clustered.h5ad"
QUERY_MODEL = OUT / "query_model"
LATENT = OUT / "s2_query_latent.csv.gz"
S3_CELLS = HERE / "s3_hlca_marker_annotation" / "s3_cell_assignment.csv.gz"

RULES = {
    "query": {"source": "analysis/GSE178360/processed/final_clustered.h5ad", "layer": "counts", "gene_match": "Ensembl ID, version stripped",
              "missing_genes": "zero-filled and counted", "batch_key": "dataset", "batch_value": "Tata_unpubl_GSE178360", "labels_key": "scanvi_label", "labels_value": "unlabeled"},
    "surgery": {"api": "scvi.model.SCANVI.load_query_data", "freeze_dropout": True, "max_epochs": 500, "weight_decay": 0.0,
                "reduce_lr_on_plateau": True, "lr_patience": 8, "lr_factor": 0.1,
                "early_stopping": {"monitor": "elbo_validation", "patience": 10, "min_delta": 0.001, "mode": "min"},
                "train_size": 0.9, "batch_size": 128, "seed": 0},
    "transfer": {"k": 50, "metric": "euclidean", "weights": "exp(-d / (2/std(d))^2), normalised per query cell (scArches weighted_knn_transfer)",
                 "uncertainty": "1 - weight of the transferred label", "unknown_cutoff_primary": 0.2, "unknown_cutoff_sensitivity": 0.3,
                 "reference_cells": "HLCA core only", "labels": ["ann_level_1", "ann_level_2", "ann_level_3", "ann_level_4", "ann_level_5", "ann_finest_level"]},
    "at0_strict_gate": "SFTPC>0 & SCGB3A2>0 & EPCAM>0 & PTPRC==0 & PECAM1==0 & COL1A1==0 on raw counts",
    "at0_concordance": "within a factor of 2 per donor",
    "unit": "donor",
    "deviations": ["scarches 0.6.1 not importable with anndata 0.13; scvi-tools scArches implementation used; weighted kNN reimplemented",
                   "early stopping on validation ELBO with a 90/10 split instead of the paper's full-dataset monitoring"],
}


def strip_version(ids: pd.Series) -> pd.Series:
    return ids.astype(str).str.split(".").str[0]


# ------------------------------------------------------------------ surgery
def stage_surgery(rec: RunRecord, timing_test: bool = False) -> None:
    import anndata as ad
    import scanpy as sc
    import scipy.sparse as sp
    import scvi
    import torch

    scvi.settings.seed = RULES["surgery"]["seed"]
    torch.set_num_threads(max(1, torch.get_num_threads()))
    for p in (OBJ, MODEL_LEGACY / "var_names.csv", GENE_ORDER):
        rec.add_input(p)
    adata = ad.read_h5ad(OBJ)
    model_genes = pd.read_csv(MODEL_LEGACY / "var_names.csv", header=None)[0].astype(str).tolist()
    query_ids = strip_version(adata.var["gene_id"])
    id_to_col = {}
    for i, g in enumerate(query_ids):
        id_to_col.setdefault(g, i)  # first occurrence wins; duplicates counted below
    dup = int(query_ids.duplicated().sum())
    cols = [id_to_col.get(g, -1) for g in model_genes]
    missing = [g for g, c in zip(model_genes, cols) if c < 0]
    counts = adata.layers["counts"].tocsc()
    present_idx = [i for i, c in enumerate(cols) if c >= 0]
    sub = counts[:, [cols[i] for i in present_idx]].tocsr().astype(np.float32)  # n_obs x n_present
    # scatter the present columns into the model gene order; missing genes stay all-zero
    P = sp.csr_matrix((np.ones(len(present_idx), dtype=np.float32), (np.arange(len(present_idx)), np.asarray(present_idx))),
                      shape=(len(present_idx), len(model_genes)))
    X = (sub @ P).tocsr()
    query = ad.AnnData(X=X, obs=adata.obs[["sample_id", "leiden_cluster", "proposed_cell_type", "total_counts", "n_genes_by_counts"]].copy())
    query.var_names = model_genes
    query.obs[RULES["query"]["batch_key"]] = RULES["query"]["batch_value"]
    query.obs[RULES["query"]["labels_key"]] = RULES["query"]["labels_value"]
    sym = pd.read_csv(GENE_ORDER).set_index("gene_id")["gene_symbol"]
    rec.set("query_gene_space", {"model_genes": len(model_genes), "present": len(present_idx), "missing": len(missing),
                                 "missing_symbols": [sym.get(g, g) for g in missing], "duplicate_query_ids": dup,
                                 "fraction_of_query_counts_retained": float(X.sum() / counts.sum())})
    print("gene space: present", len(present_idx), "missing", len(missing), missing[:10])

    if not (MODEL / "model.pt").exists():
        scvi.model.SCANVI.convert_legacy_save(str(MODEL_LEGACY), str(MODEL), overwrite=True)
    rec.add_input(MODEL / "model.pt")
    scvi.model.SCANVI.prepare_query_anndata(query, str(MODEL))
    model = scvi.model.SCANVI.load_query_data(query, str(MODEL), freeze_dropout=RULES["surgery"]["freeze_dropout"])
    rec.set("query_model_summary", str(model))
    import time
    t0 = time.time()
    model.train(
        max_epochs=1 if timing_test else RULES["surgery"]["max_epochs"],
        train_size=RULES["surgery"]["train_size"],
        batch_size=RULES["surgery"]["batch_size"],
        plan_kwargs={"weight_decay": RULES["surgery"]["weight_decay"], "reduce_lr_on_plateau": RULES["surgery"]["reduce_lr_on_plateau"],
                     "lr_patience": RULES["surgery"]["lr_patience"], "lr_factor": RULES["surgery"]["lr_factor"]},
        check_val_every_n_epoch=1,
        early_stopping=True,
        early_stopping_monitor=RULES["surgery"]["early_stopping"]["monitor"],
        early_stopping_patience=RULES["surgery"]["early_stopping"]["patience"],
        early_stopping_min_delta=RULES["surgery"]["early_stopping"]["min_delta"],
        early_stopping_mode=RULES["surgery"]["early_stopping"]["mode"],
        n_samples_per_label=None,
        accelerator="cpu",
    )
    elapsed = time.time() - t0
    hist = {k: [float(x) for x in np.asarray(v).ravel()] for k, v in model.history.items() if "elbo" in k}
    rec.set("training", {"epochs_run": len(next(iter(hist.values()))) if hist else None, "seconds": round(elapsed, 1),
                         "elbo_train_last": hist.get("elbo_train", [None])[-1], "elbo_validation_last": hist.get("elbo_validation", [None])[-1],
                         "elbo_validation_min": min(hist["elbo_validation"]) if "elbo_validation" in hist else None})
    if timing_test:
        print(f"timing test: 1 epoch in {elapsed:.1f} s; history keys {list(model.history.keys())}")
        return
    (OUT / "s2_training_history.json").write_text(json.dumps(hist), encoding="utf-8")
    rec.add_output(OUT / "s2_training_history.json")
    model.save(str(QUERY_MODEL), overwrite=True)
    latent = model.get_latent_representation()
    pd.DataFrame(latent, index=query.obs_names, columns=[f"z{i}" for i in range(latent.shape[1])]).to_csv(LATENT, compression="gzip")
    rec.set("latent_shape", list(latent.shape))
    # query-side predictions from the classifier itself (coarse label set), for the record only
    pred = model.predict(soft=False)
    rec.set("classifier_coarse_prediction_counts", dict(Counter(pred)))
    print("surgery done; epochs", rec.record["results"]["training"]["epochs_run"], "latent", latent.shape)


# ----------------------------------------------------------------- transfer
def weighted_knn_transfer(ref_emb: np.ndarray, ref_labels: pd.DataFrame, query_emb: np.ndarray, k: int):
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=k, metric=RULES["transfer"]["metric"], n_jobs=-1).fit(ref_emb)
    dist, idx = nn.kneighbors(query_emb)
    stds = np.std(dist, axis=1)
    stds = (2.0 / np.where(stds > 0, stds, 1e-12)) ** 2
    w = np.exp(-dist / stds.reshape(-1, 1))
    w = w / w.sum(axis=1, keepdims=True)
    out = {}
    for col in ref_labels.columns:
        lab = ref_labels[col].astype(str).to_numpy()
        neigh = lab[idx]  # (n_query, k)
        best, unc = [], []
        for i in range(neigh.shape[0]):
            probs = {}
            for j in range(k):
                probs[neigh[i, j]] = probs.get(neigh[i, j], 0.0) + w[i, j]
            b = max(probs, key=probs.get)
            best.append(b); unc.append(max(1.0 - probs[b], 0.0))
        out[col] = np.array(best, dtype=object); out[col + "_uncertainty"] = np.array(unc)
    return pd.DataFrame(out), dist.mean(axis=1)


def stage_transfer(rec: RunRecord) -> None:
    import anndata as ad
    import h5py
    from anndata.io import read_elem
    import scanpy as sc
    sys.path.insert(0, str(HERE))
    from s3_hlca_marker_annotation import SYNONYMS, COMPARTMENT_OF_BLIND, verdict_for  # noqa: E402

    for p in (EMB, LATENT, OBJ):
        rec.add_input(p)
    with h5py.File(EMB, "r") as f:
        obs = read_elem(f["obs"])
        obsm_keys = list(f["obsm"].keys())
        cand = [k for k in obsm_keys if "scanvi" in k.lower()] or [k for k in obsm_keys if "emb" in k.lower()]
        if len(cand) == 1:
            emb_key, emb_src = cand[0], f["obsm"][cand[0]]
        elif isinstance(f["X"], h5py.Dataset) and f["X"].ndim == 2 and f["X"].shape[1] == 30:
            emb_key, emb_src = "X (30-dim scANVI embedding stored as X in HLCA_full_v1.1_emb.h5ad)", f["X"]
        else:
            raise SystemExit(f"could not identify the embedding unambiguously: obsm {obsm_keys}, X {f['X']}")
        if "core_or_extension" in obs.columns:
            core_col = "core_or_extension"
        else:
            cand_cols = [c for c in obs.columns if "core_or" in c.lower() or c.lower().startswith("core")]
            if not cand_cols:
                raise SystemExit(f"no core/extension column in reference obs: {list(obs.columns)[:60]}")
            core_col = cand_cols[0]
        is_core = obs[core_col].astype(str).str.lower().eq("core").to_numpy()
        core_rows = np.where(is_core)[0]
        if len(core_rows) == 0:
            raise SystemExit(f"column {core_col} selected zero core cells; values: {obs[core_col].astype(str).value_counts().head().to_dict()}")
        ref_emb = emb_src[:][core_rows]
    label_cols = [c for c in RULES["transfer"]["labels"] if c in obs.columns]
    ref_labels = obs.iloc[core_rows][label_cols].reset_index(drop=True)
    rec.set("reference", {"embedding_key": emb_key, "core_column": core_col, "n_core_cells": int(len(core_rows)), "n_all_cells": int(len(obs)),
                          "label_columns": label_cols, "finest_label_classes": int(ref_labels[label_cols[-1]].nunique())})
    print("reference core cells", len(core_rows), "labels", label_cols)

    latent = pd.read_csv(LATENT, index_col=0)
    tr, mean_dist = weighted_knn_transfer(ref_emb, ref_labels, latent.to_numpy(dtype=np.float32), RULES["transfer"]["k"])
    tr.index = latent.index
    tr["mean_neighbour_distance"] = mean_dist
    adata = ad.read_h5ad(OBJ)
    assert (adata.obs_names == tr.index).all()
    finest = label_cols[-1]
    u = tr[finest + "_uncertainty"].to_numpy()
    for cut in (RULES["transfer"]["unknown_cutoff_primary"], RULES["transfer"]["unknown_cutoff_sensitivity"]):
        tr[f"{finest}_called_u{cut}"] = np.where(u > cut, "Unknown", tr[finest])
    for c in ("sample_id", "leiden_cluster", "proposed_cell_type"):
        tr[c] = adata.obs[c].astype(str).to_numpy()
    tr.to_csv(OUT / "s2_transfer_per_cell.csv.gz", compression="gzip")

    # ------------------------------------------------ per donor and cluster
    cut = RULES["transfer"]["unknown_cutoff_primary"]
    called = tr[f"{finest}_called_u{cut}"]
    per_donor = []
    for donor, sub in tr.groupby("sample_id"):
        uu = sub[finest + "_uncertainty"]
        per_donor.append({"donor": donor, "n_cells": len(sub), "mean_uncertainty": round(float(uu.mean()), 4), "median_uncertainty": round(float(uu.median()), 4),
                          "unknown_fraction_u0.2": round(float((uu > 0.2).mean()), 4), "unknown_fraction_u0.3": round(float((uu > 0.3).mean()), 4)})
    pd.DataFrame(per_donor).to_csv(OUT / "s2_uncertainty_per_donor.csv", index=False)
    rec.add_output(OUT / "s2_uncertainty_per_donor.csv")

    rows = []
    blind = tr["proposed_cell_type"].str.replace(" (candidate)", "", regex=False)
    s3 = pd.read_csv(S3_CELLS, index_col=0) if S3_CELLS.exists() else None
    for cl in sorted(tr["leiden_cluster"].unique(), key=int):
        m = (tr["leiden_cluster"] == cl).to_numpy(); n = int(m.sum())
        cnt = Counter(called[m]); mode, mode_n = cnt.most_common(1)[0]
        cnt_all = Counter(tr.loc[m, finest]); mode_all, mode_all_n = cnt_all.most_common(1)[0]
        comp = Counter(tr.loc[m, "ann_level_1"]).most_common(1)[0][0] if "ann_level_1" in tr else "NA"
        b = blind[m].iloc[0]
        row = {"cluster": cl, "n_cells": n, "blind_proposal": b, "hlca_level1": comp,
               "transferred_finest_mode": mode_all, "mode_fraction": round(mode_all_n / n, 3),
               "mode_after_unknown_u0.2": mode, "mode_fraction_after_unknown": round(mode_n / n, 3),
               "unknown_fraction_u0.2": round(float((tr.loc[m, finest + "_uncertainty"] > 0.2).mean()), 3),
               "mean_uncertainty": round(float(tr.loc[m, finest + "_uncertainty"].mean()), 4),
               "verdict_vs_blind": verdict_for(b, mode_all, comp if comp in ("Epithelial", "Immune", "Endothelial", "Stroma") else "NA")}
        if s3 is not None:
            s3m = s3.loc[tr.index[m]]
            row["s3_flat_mode"] = Counter(s3m["hlca_type"]).most_common(1)[0][0]
            row["s3_hier_mode"] = Counter(s3m["hlca_type_hier"]).most_common(1)[0][0]
            row["agrees_with_s3_flat"] = row["s3_flat_mode"] == mode_all
            row["agrees_with_s3_hier"] = row["s3_hier_mode"] == mode_all
        rows.append(row)
    clus = pd.DataFrame(rows); clus.to_csv(OUT / "s2_cluster_transfer.csv", index=False); rec.add_output(OUT / "s2_cluster_transfer.csv")
    vc = Counter(clus["verdict_vs_blind"]); rec.set("verdict_counts_vs_blind", dict(vc))
    rec.set("per_donor", per_donor)

    # ------------------------------------------------------------- AT0 check
    counts = adata.layers["counts"]
    gi = {g: list(adata.var_names).index(g) for g in ("SFTPC", "SCGB3A2", "EPCAM", "PTPRC", "PECAM1", "COL1A1")}
    col = lambda g: np.asarray(counts[:, gi[g]].todense()).ravel()  # noqa: E731
    gate = (col("SFTPC") > 0) & (col("SCGB3A2") > 0) & (col("EPCAM") > 0) & (col("PTPRC") == 0) & (col("PECAM1") == 0) & (col("COL1A1") == 0)
    at0_rows = []
    for donor in sorted(tr["sample_id"].unique()):
        m = (tr["sample_id"] == donor).to_numpy()
        g = int(gate[m].sum())
        a_all = int((tr.loc[m, finest] == "AT0").sum()); a_conf = int((called[m] == "AT0").sum())
        p_all = int((tr.loc[m, finest] == "pre-TB secretory").sum()); p_conf = int((called[m] == "pre-TB secretory").sum())
        both = int((gate[m] & (called[m] == "AT0").to_numpy()).sum())
        r = a_conf / g if g else None
        row = {"donor": donor, "strict_gate": g, "transferred_AT0_all": a_all, "transferred_AT0_u_le_0.2": a_conf,
               "transferred_preTB_all": p_all, "transferred_preTB_u_le_0.2": p_conf, "gate_and_AT0": both,
               "AT0_over_gate": round(r, 3) if r is not None else None, "within_factor_2": (0.5 <= r <= 2.0) if r is not None else None}
        if s3 is not None:
            s3m = s3.loc[tr.index[m]]
            row["s3_flat_AT0"] = int((s3m["hlca_type"] == "AT0").sum()); row["s3_hier_AT0"] = int((s3m["hlca_type_hier"] == "AT0").sum())
        at0_rows.append(row)
    at0 = pd.DataFrame(at0_rows); at0.to_csv(OUT / "s2_at0_check_per_donor.csv", index=False); rec.add_output(OUT / "s2_at0_check_per_donor.csv")
    rec.set("at0_per_donor", at0_rows)
    at0_by_cluster = pd.crosstab(tr["leiden_cluster"][called == "AT0"], tr["sample_id"][called == "AT0"])
    at0_by_cluster.to_csv(OUT / "s2_at0_cells_by_cluster_and_donor.csv"); rec.add_output(OUT / "s2_at0_cells_by_cluster_and_donor.csv")

    epi_result = None
    if EPI_OBJ.exists():
        epi = ad.read_h5ad(EPI_OBJ, backed="r")
        key = "leiden_cluster" if "leiden_cluster" in epi.obs else [c for c in epi.obs.columns if c.startswith("leiden")][0]
        common = epi.obs.index.intersection(tr.index)
        tab = pd.crosstab(epi.obs.loc[common, key].astype(str), called.loc[common])
        tab.to_csv(OUT / "s2_epithelial_subcluster_by_transferred_label.csv"); rec.add_output(OUT / "s2_epithelial_subcluster_by_transferred_label.csv")
        if "4" in tab.index:
            row = tab.loc["4"]
            epi_result = {"n": int(row.sum()), "AT0": int(row.get("AT0", 0)), "pre-TB secretory": int(row.get("pre-TB secretory", 0)),
                          "AT2": int(row.get("AT2", 0)), "Unknown": int(row.get("Unknown", 0)), "top5": row.sort_values(ascending=False).head(5).to_dict()}
        epi.file.close()
    rec.set("epithelial_subcluster_4_transferred", epi_result)

    c30 = clus.loc[clus["cluster"] == "30"].to_dict("records")
    rec.set("cluster_30", c30[0] if c30 else None)

    # ------------------------------------------------------------------
    # Post hoc comparison, not pre-registered: the v1.1 embedding file carries
    # the HLCA authors' own label transfer for every extension dataset,
    # including Tata_unpubl, which is GSE178360. Cells are matched by their
    # 16-mer 10x barcode when that barcode is unique on both sides.
    # ------------------------------------------------------------------
    own = None
    tata_mask = obs["dataset"].astype(str).str.contains("Tata", case=False).to_numpy()
    if tata_mask.any():
        import re
        hl = obs.loc[tata_mask].copy()
        hl_cols = [c for c in hl.columns if c.startswith("transf_ann_level_")] + [c for c in ("sample", "donor_id", "ann_finest_level") if c in hl.columns]
        hl = hl[hl_cols]
        hl["barcode16"] = [m.group(0) if (m := re.search(r"[ACGT]{16}", str(i))) else "" for i in hl.index]
        hl.to_csv(OUT / "s2_hlca_own_transfer_tata_unpubl.csv.gz", compression="gzip")  # regenerable from the Zenodo file; not tracked
        hl_unique = hl[hl["barcode16"] != ""].drop_duplicates("barcode16", keep=False).set_index("barcode16")
        ours = pd.DataFrame({"barcode16": [m.group(0) if (m := re.search(r"[ACGT]{16}", str(b))) else "" for b in adata.obs["original_barcode"].astype(str)]}, index=tr.index)
        ours = ours[ours["barcode16"] != ""].drop_duplicates("barcode16", keep=False)
        common = ours["barcode16"][ours["barcode16"].isin(hl_unique.index)]
        own = {"hlca_tata_cells": int(tata_mask.sum()), "hlca_unique_barcodes": int(len(hl_unique)), "our_unique_barcodes": int(len(ours)),
               "matched_cells": int(len(common)), "levels": {}}
        if len(common):
            ours_idx = common.index
            hl_m = hl_unique.loc[common.values]
            for lvl in (3, 4, 5):
                mine = f"ann_level_{lvl}"
                theirs, theirs_u = f"transf_ann_level_{lvl}_label", f"transf_ann_level_{lvl}_uncert"
                if mine not in tr.columns or theirs not in hl_m.columns:
                    continue
                a = tr.loc[ours_idx, mine].astype(str).to_numpy(); ua = tr.loc[ours_idx, mine + "_uncertainty"].to_numpy()
                b = hl_m[theirs].astype(str).to_numpy(); ub = pd.to_numeric(hl_m[theirs_u], errors="coerce").to_numpy()
                missing = {"nan", "None", "NaN", "", "<NA>"}
                has_label = ~np.isin(a, list(missing)) & ~np.isin(b, list(missing))
                a, b, ua, ub = a[has_label], b[has_label], ua[has_label], ub[has_label]
                both_conf = (ua <= cut) & (ub <= cut)
                own["levels"][f"level_{lvl}"] = {
                    "cells_with_label_on_both_sides": int(has_label.sum()),
                    "agreement_all_matched": round(float((a == b).mean()), 4) if has_label.any() else None,
                    "agreement_when_both_confident_u_le_0.2": round(float((a[both_conf] == b[both_conf]).mean()), 4) if both_conf.any() else None,
                    "fraction_both_confident": round(float(both_conf.mean()), 4),
                    "our_mean_uncertainty": round(float(np.nanmean(ua)), 4), "hlca_mean_uncertainty": round(float(np.nanmean(ub)), 4),
                    "uncertainty_spearman": round(float(pd.Series(ua).corr(pd.Series(ub), method="spearman")), 4),
                }
            ct = (pd.crosstab(pd.Series(tr.loc[ours_idx, "ann_level_4"].astype(str).to_numpy(), name="ours_level_4"),
                              pd.Series(hl_m["transf_ann_level_4_label"].astype(str).to_numpy(), name="hlca_own_level_4"))
                  if "transf_ann_level_4_label" in hl_m.columns else None)
            if ct is not None:
                ct.to_csv(OUT / "s2_ours_vs_hlca_own_transfer_level4_crosstab.csv"); rec.add_output(OUT / "s2_ours_vs_hlca_own_transfer_level4_crosstab.csv")
    rec.set("comparison_with_hlca_own_transfer_post_hoc", own)

    # --------------------------------------------------------------- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    adata.obs["hlca_transferred"] = called.to_numpy()
    adata.obs["hlca_uncertainty"] = u
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sc.pl.umap(adata, color="hlca_transferred", ax=axes[0], show=False, legend_loc="on data", legend_fontsize=5, size=4, title=f"HLCA label transfer (finest level, Unknown if u > {cut})")
    sc.pl.umap(adata, color="hlca_uncertainty", ax=axes[1], show=False, size=4, title="label-transfer uncertainty", color_map="viridis")
    fig.tight_layout(); fig.savefig(OUT / "s2_umap_transfer_and_uncertainty.png", dpi=110); plt.close(fig)
    rec.add_output(OUT / "s2_umap_transfer_and_uncertainty.png")
    fig, ax = plt.subplots(figsize=(7, 4))
    for donor, sub in tr.groupby("sample_id"):
        ax.hist(sub[finest + "_uncertainty"], bins=50, histtype="step", label=donor)
    ax.axvline(0.2, color="k", ls="--", lw=1); ax.set_xlabel("uncertainty"); ax.set_ylabel("cells"); ax.legend(); ax.set_title("Label-transfer uncertainty per donor")
    fig.tight_layout(); fig.savefig(OUT / "s2_uncertainty_hist_per_donor.png", dpi=120); plt.close(fig)
    rec.add_output(OUT / "s2_uncertainty_hist_per_donor.png")

    # --------------------------------------------------------------- summary
    lines = ["# Trial S2 output: scArches mapping of GSE178360 to the HLCA core and label transfer", "",
             f"Reference: {rec.record['results']['reference']['n_core_cells']} HLCA core cells, embedding `{emb_key}`; transferred levels {label_cols}; "
             f"k = {RULES['transfer']['k']}; Unknown if u > {cut}. Deviations are listed in the script docstring.", "",
             "## Uncertainty per donor", "", df_to_markdown(pd.DataFrame(per_donor), index=False), "",
             "## Cluster-level transfer", "",
             "| Cluster | Cells | Blind proposal | HLCA level 1 | Transferred finest (mode) | Mode fraction | Unknown fraction | Mean uncertainty | Verdict vs blind | S3 flat | S3 hier |",
             "|--:|--:|---|---|---|--:|--:|--:|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['cluster']} | {r['n_cells']} | {r['blind_proposal']} | {r['hlca_level1']} | {r['transferred_finest_mode']} | {r['mode_fraction']} | "
                     f"{r['unknown_fraction_u0.2']} | {r['mean_uncertainty']} | {r['verdict_vs_blind']} | {r.get('s3_flat_mode', '')} | {r.get('s3_hier_mode', '')} |")
    lines += ["", "Verdicts versus the blind proposals: " + ", ".join(f"{k}: {v}" for k, v in sorted(vc.items())) + ".", "",
              "## AT0 check per donor", "", df_to_markdown(at0, index=False), ""]
    if epi_result:
        lines.append(f"Epithelial subcluster 4 (existing AT0 candidate analogue, n={epi_result['n']}): transferred AT0 {epi_result['AT0']}, "
                     f"pre-TB secretory {epi_result['pre-TB secretory']}, AT2 {epi_result['AT2']}, Unknown {epi_result['Unknown']}; top {epi_result['top5']}.")
    if own:
        lines += ["", "## Post hoc: our transfer versus the HLCA authors' own transfer for the same dataset (Tata_unpubl = GSE178360)", "",
                  f"HLCA Tata_unpubl cells {own['hlca_tata_cells']}; matched by unique 16-mer barcode: {own['matched_cells']} cells.", ""]
        for lvl, v in own["levels"].items():
            lines.append(f"- {lvl} ({v['cells_with_label_on_both_sides']} cells labelled on both sides): agreement {v['agreement_all_matched']} over those cells; {v['agreement_when_both_confident_u_le_0.2']} when both sides are confident "
                         f"(fraction both confident {v['fraction_both_confident']}); mean uncertainty ours {v['our_mean_uncertainty']} vs HLCA {v['hlca_mean_uncertainty']}; "
                         f"Spearman of uncertainties {v['uncertainty_spearman']}.")
    (OUT / "s2_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "s2_summary.md")
    print("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["surgery", "transfer"], required=True)
    ap.add_argument("--timing-test", action="store_true", help="one epoch only; writes a separate record and no outputs")
    args = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    rec_path = OUT / (f"s2_run_record_{args.stage}" + ("_timingtest" if args.timing_test else "") + ".json")
    rec = RunRecord(rec_path, f"S2 reference mapping, stage {args.stage}" + (" (timing test, not a result)" if args.timing_test else ""),
                    RULES, notes="rules frozen in ANALYSIS_TRIAL_PLAN.md before this script was written")
    if args.stage == "surgery":
        stage_surgery(rec, timing_test=args.timing_test)
    else:
        stage_transfer(rec)
    rec.finish()


if __name__ == "__main__":
    main()
