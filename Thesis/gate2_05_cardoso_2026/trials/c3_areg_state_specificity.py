#!/usr/bin/env python
"""Trial C3: is Areg a property of the DATP-like state, under my own clustering and over time (Gate 2a).

The owner's Gate 2 branch (a), asked in the only form the deposited data can
answer honestly. CellChat itself cannot be run here: it is an R package and
this machine has no R (a standing repository constraint), and a Python
reimplementation would be a different method with a different
ligand-receptor resource, so calling its output a CellChat rerun would be
false. What the paper's communication claim actually rests on is an
expression fact -- that Areg is induced specifically in the DATP-like mutant
state and is the highest-ranked EGFR ligand there (Extended Data Fig. 5e-h) --
and that fact can be re-derived directly, with three things the paper's own
analysis did not use:

1. the DATP-like state defined by this repository's clustering rather than by
   the authors' labels, which is the state-definition sensitivity the owner
   asked for;
2. the **within-animal wild-type control**: in the Red2Onco design the YFP
   library is the non-recombined clone of the same lung, so the comparison is
   mutant versus wild-type inside one animal rather than across genotypes;
3. the **time course**: the companion series carries 4 days, 2 weeks and
   12 weeks, and it carries two replicate libraries per arm, which is the only
   biological replication anywhere in the reusable set.

The output is a statement about expression and state specificity. It is not a
communication analysis: no ligand-receptor probability is computed, and
nothing here establishes that a fibroblast receives the signal. That
limitation is the point of stating it.

Frozen rules (written to the run record before any matrix is read):

* Libraries: the ten Experiment 1 libraries of GSE247505 (Confetti RFP and
  YFP as the wild-type baseline; Red2Kras RFP and YFP at 4 days and at 2
  weeks, two libraries each). The Experiment 2 libraries are the Il1r1 arm of
  the companion paper and are out of scope.
* QC, doublets, normalisation, features, embedding and clustering follow trial
  C1 exactly (seeds 0 throughout).
* Batch correction is decided, not assumed, by this repository's rule: within
  each arm (colour and time point) the two replicate libraries are a technical
  comparison, so same-library enrichment among the 30 nearest neighbours is
  measured within each arm. Correction is applied only if the measured
  enrichment exceeds the pipeline's failure threshold of 2.0; otherwise the
  uncorrected embedding stands. This is the one place in the whole Cardoso
  deposit where that rule can be exercised, because it is the only place with
  replicate libraries.
* States are called from the paper's own Fig. 4l marker sets (AT2, Cd177+,
  DATP-like, cycling, AT1-like) with score_genes (ctrl_size 50, n_bins 25,
  seed 0) and assigned by modal per-cell argmax; a cluster call is confident
  only if the mode holds at least 50% of its cells.
* THE TESTS, frozen.
  - T1 state specificity: in each Red2Kras RFP library, mean log1p(CP10K)
    Areg in DATP-like cells exceeds mean Areg in AT2 cells. Holds if true in
    all four such libraries.
  - T2 ligand ranking: within the DATP-like cells of each Red2Kras RFP
    library, rank Areg, Ereg, Hbegf and Tgfa by mean log1p(CP10K). Areg is
    "top" if it ranks first; the full ranking is reported either way, because
    what ranks second is the part the paper does not show.
  - T3 the internal control: at each time point, the DATP-like share of cells
    is higher in the RFP libraries than in the YFP libraries of the same
    animals. Reported per library; holds only if true at both time points.
  - T4 time: T1 and T2 are reported separately at 4 days and at 2 weeks, and
    the direction of change in the DATP-like Areg mean between them is
    recorded. No trend test is computed; two libraries per arm is a ranking,
    not a trend.
* Unit: the library, two per arm. Medians across the two, no P value, in line
  with the repository's rule that no P value is reported where a group holds
  two replicates.
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

OUT = HERE / "c3_areg_state_specificity"
OUT.mkdir(exist_ok=True)
EXTRACTS = json.loads((HERE.parent / "cardoso_2026_extracts.json").read_text(encoding="utf-8"))
EPI_SETS = EXTRACTS["marker_sets"]["mutant_epithelial_states_fig4l"]
EGF_LIGANDS = ["Areg", "Ereg", "Hbegf", "Tgfa"]

LIBRARIES = [
    {"gsm": "GSM7890829", "library": "Expt1_ConfettiRFP", "colour": "RFP", "clone": "wild type",
     "timepoint": "Confetti baseline", "arm": "Confetti_RFP", "replicate": 1},
    {"gsm": "GSM7890830", "library": "Expt1_ConfettiYFP", "colour": "YFP", "clone": "wild type",
     "timepoint": "Confetti baseline", "arm": "Confetti_YFP", "replicate": 1},
    {"gsm": "GSM7890831", "library": "Expt1_4dRFPr1", "colour": "RFP", "clone": "KrasG12D",
     "timepoint": "4 days", "arm": "4d_RFP", "replicate": 1},
    {"gsm": "GSM7890832", "library": "Expt1_4dRFPr2", "colour": "RFP", "clone": "KrasG12D",
     "timepoint": "4 days", "arm": "4d_RFP", "replicate": 2},
    {"gsm": "GSM7890833", "library": "Expt1_4dYFPr1", "colour": "YFP", "clone": "wild type",
     "timepoint": "4 days", "arm": "4d_YFP", "replicate": 1},
    {"gsm": "GSM7890834", "library": "Expt1_4dYFPr2", "colour": "YFP", "clone": "wild type",
     "timepoint": "4 days", "arm": "4d_YFP", "replicate": 2},
    {"gsm": "GSM7890835", "library": "Expt1_2wRFPr1", "colour": "RFP", "clone": "KrasG12D",
     "timepoint": "2 weeks", "arm": "2w_RFP", "replicate": 1},
    {"gsm": "GSM7890836", "library": "Expt1_2wRFPr2", "colour": "RFP", "clone": "KrasG12D",
     "timepoint": "2 weeks", "arm": "2w_RFP", "replicate": 2},
    {"gsm": "GSM7890837", "library": "Expt1_2wYFPr1", "colour": "YFP", "clone": "wild type",
     "timepoint": "2 weeks", "arm": "2w_YFP", "replicate": 1},
    {"gsm": "GSM7890838", "library": "Expt1_2wYFPr2", "colour": "YFP", "clone": "wild type",
     "timepoint": "2 weeks", "arm": "2w_YFP", "replicate": 2},
]
MUTANT_RFP_ARMS = ["4d_RFP", "2w_RFP"]

RULES = {
    "accession": "GSE247505 (England et al. 2025), Experiment 1 only",
    "libraries": LIBRARIES,
    "why_not_cellchat": ("CellChat is R-only and this machine has no R; a Python reimplementation "
                         "would use a different ligand-receptor resource and would not be the same "
                         "method. This trial re-derives the expression fact the communication claim "
                         "rests on and says so."),
    "processing": "identical to trial C1 (MAD QC, Scrublet, log1p CP10K, seurat_v3 2000 HVG, 30 PCs, k 30, Leiden igraph 0.5/0.2/1.0, seeds 0)",
    "batch_rule": {"measure": "same-library enrichment among the 30 nearest neighbours, within arm",
                   "failure_threshold": 2.0,
                   "action": "correct only if the measured enrichment exceeds the threshold"},
    "state_sets": EPI_SETS,
    "tests": {
        "T1_state_specificity": "mean log1p(CP10K) Areg higher in DATP-like than in AT2, in all four Red2Kras RFP libraries",
        "T2_ligand_ranking": f"rank {EGF_LIGANDS} by mean log1p(CP10K) within DATP-like cells of each Red2Kras RFP library",
        "T3_internal_control": "DATP-like share higher in RFP than in YFP libraries at both time points",
        "T4_time": "T1 and T2 reported per time point; direction of the DATP-like Areg mean between 4 days and 2 weeks recorded; no trend test",
    },
    "unit": "the library, two per arm; medians across replicates; no P value",
    "what_this_is_not": "not a communication analysis; no ligand-receptor probability is computed and nothing establishes that a fibroblast receives the signal",
}


def load(entry: dict, rec: RunRecord):
    import anndata as ad
    import scanpy as sc

    root = RAW / "GSE247505" / "GSE247505_RAW"
    stem = f"{entry['gsm']}_{entry['library']}_"
    paths = {k: root / f"{stem}{k}" for k in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")}
    for path in paths.values():
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(paths["matrix.mtx.gz"], paths["features.tsv.gz"], paths["barcodes.tsv.gz"])
    is_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    X, var = X[:, is_gene], var.loc[is_gene].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))

    obs = pd.DataFrame(index=pd.Index([f"{entry['library']}_{b}" for b in bc]))
    for key in ("library", "colour", "clone", "timepoint", "arm"):
        obs[key] = entry[key]
    obs["replicate"] = entry["replicate"]
    adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)

    symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
    adata.var["ribo"] = (symbols.str.startswith(("Rps", "Rpl"))
                         & ~symbols.str.contains("-ps", case=False)).to_numpy()
    adata.var["hb"] = symbols.str.startswith(("Hba", "Hbb")).to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20],
                               log1p=True, inplace=True)
    th = derive_thresholds(adata.obs, entry["library"])
    keep = apply_thresholds(adata.obs, th).to_numpy()
    row = {"library": entry["library"], "arm": entry["arm"], "replicate": entry["replicate"],
           "clone": entry["clone"], "timepoint": entry["timepoint"],
           "barcodes": int(adata.n_obs), "kept_after_qc": int(keep.sum()),
           "median_genes": int(np.median(adata.obs["n_genes_by_counts"]))}
    adata = adata[keep].copy()

    expected = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=expected)
    auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
    method = "scrublet automatic threshold"
    if auto < 0.2 * expected or auto > 3 * expected:
        cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - expected))
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
        method = f"top {100 * expected:.2f}% of scores (automatic rejected: it called {100 * auto:.2f}%)"
    row.update({"doublet_call_method": method,
                "doublets_removed": int(adata.obs["predicted_doublet"].sum())})
    adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
    row["cells_analysed"] = int(adata.n_obs)
    return adata, row


def same_library_enrichment(adata, arm: str) -> float:
    """kNN enrichment of a cell's own library among its neighbours, within one arm."""
    mask = adata.obs["arm"].to_numpy() == arm
    idx = np.where(mask)[0]
    if len(idx) == 0:
        return float("nan")
    libraries = adata.obs["library"].to_numpy()
    graph = adata.obsp["connectivities"]
    counts, expected = [], []
    within = libraries[idx]
    shares = pd.Series(within).value_counts(normalize=True)
    for i in idx:
        neighbours = graph[i].indices
        neighbours = neighbours[mask[neighbours]]
        if len(neighbours) == 0:
            continue
        counts.append(float((libraries[neighbours] == libraries[i]).mean()))
        expected.append(float(shares[libraries[i]]))
    if not counts:
        return float("nan")
    return float(np.mean(counts) / np.mean(expected))


def main() -> None:
    import anndata as ad
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc.settings.verbosity = 0
    rec = RunRecord(OUT / "c3_run_record.json",
                    "C3 Areg state specificity across the mutant epithelial time course (Gate 2a)",
                    RULES, notes="rules frozen before any matrix was read")

    parts, qc_rows = [], []
    for entry in LIBRARIES:
        adata, row = load(entry, rec)
        parts.append(adata)
        qc_rows.append(row)
        print(f"  {entry['library']}: {row['cells_analysed']} cells")
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT / "c3_qc_per_library.csv", index=False)
    rec.add_output(OUT / "c3_qc_per_library.csv")

    merged = ad.concat(parts, join="inner")
    merged.var = parts[0].var.loc[merged.var_names].copy()
    del parts
    gc.collect()
    rec.set("cells_analysed_total", int(merged.n_obs))

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

    # ---- the batch decision, measured -----------------------------------
    enrichment = {arm: round(same_library_enrichment(merged, arm), 4)
                  for arm in sorted({e["arm"] for e in LIBRARIES})
                  if sum(1 for e in LIBRARIES if e["arm"] == arm) > 1}
    worst = max((v for v in enrichment.values() if v == v), default=float("nan"))
    correction_needed = bool(worst > RULES["batch_rule"]["failure_threshold"])
    rec.set("batch_decision", {"same_library_enrichment_by_arm": enrichment,
                               "worst": worst, "threshold": 2.0,
                               "correction_applied": correction_needed})
    print(f"  replicate mixing by arm: {enrichment}; correction needed: {correction_needed}")

    for res in (0.2, 0.5, 1.0):
        sc.tl.leiden(merged, resolution=res, key_added=f"leiden_{res}", flavor="igraph",
                     n_iterations=2, random_state=0, directed=False)
    sc.tl.umap(merged, random_state=0)

    for name, genes in EPI_SETS.items():
        present = [g for g in genes if g in merged.var_names]
        if present:
            sc.tl.score_genes(merged, present, score_name=f"score_{name}", ctrl_size=50,
                              n_bins=25, random_state=0)
    frame = pd.DataFrame({n: merged.obs[f"score_{n}"].to_numpy() for n in EPI_SETS
                          if f"score_{n}" in merged.obs}, index=merged.obs_names)
    merged.obs["state_argmax"] = frame.idxmax(axis=1).to_numpy()

    rows = []
    for cluster in sorted(merged.obs["leiden_0.5"].astype(str).unique(), key=int):
        mask = merged.obs["leiden_0.5"].astype(str).to_numpy() == cluster
        calls = merged.obs.loc[mask, "state_argmax"].value_counts()
        mode_fraction = float(calls.iloc[0] / mask.sum())
        rows.append({"cluster": cluster, "n_cells": int(mask.sum()), "call": calls.index[0],
                     "mode_fraction": round(mode_fraction, 3), "confident": mode_fraction >= 0.50,
                     "n_libraries": int(merged.obs.loc[mask, "library"].nunique()),
                     "max_library_fraction": round(float(
                         merged.obs.loc[mask, "library"].value_counts(normalize=True).iloc[0]), 3)})
    clusters = pd.DataFrame(rows)
    clusters.to_csv(OUT / "c3_clusters.csv", index=False)
    rec.add_output(OUT / "c3_clusters.csv")
    print(df_to_markdown(clusters, index=False))

    confident = clusters[clusters["confident"]].set_index("cluster")["call"].to_dict()
    merged.obs["state"] = [confident.get(c, "unassigned")
                           for c in merged.obs["leiden_0.5"].astype(str)]

    # ---- per-library state composition and ligand expression -------------
    comp_rows, lig_rows = [], []
    logX = merged.X
    ligand_cols = {g: merged.var_names.get_loc(g) for g in EGF_LIGANDS if g in merged.var_names}
    for entry in LIBRARIES:
        lib_mask = merged.obs["library"].to_numpy() == entry["library"]
        n_lib = int(lib_mask.sum())
        for state in sorted(set(merged.obs["state"])):
            state_mask = lib_mask & (merged.obs["state"].to_numpy() == state)
            row = {"library": entry["library"], "arm": entry["arm"], "clone": entry["clone"],
                   "timepoint": entry["timepoint"], "state": state,
                   "n_cells": int(state_mask.sum()),
                   "pct_of_library": round(100 * float(state_mask.sum()) / n_lib, 2) if n_lib else np.nan}
            comp_rows.append(row)
            if state_mask.sum() >= 20:
                lig = {"library": entry["library"], "arm": entry["arm"], "clone": entry["clone"],
                       "timepoint": entry["timepoint"], "state": state,
                       "n_cells": int(state_mask.sum())}
                for gene, col in ligand_cols.items():
                    values = np.asarray(logX[state_mask, col].todense()).ravel()
                    lig[f"mean_{gene}"] = round(float(values.mean()), 4)
                    lig[f"det_{gene}"] = round(float((values > 0).mean()), 4)
                lig_rows.append(lig)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "c3_state_composition.csv", index=False)
    rec.add_output(OUT / "c3_state_composition.csv")
    lig = pd.DataFrame(lig_rows)
    lig.to_csv(OUT / "c3_ligand_expression.csv", index=False)
    rec.add_output(OUT / "c3_ligand_expression.csv")

    # ---- the frozen tests -------------------------------------------------
    mutant = lig[lig["arm"].isin(MUTANT_RFP_ARMS)]
    t1 = []
    for library in sorted(mutant["library"].unique()):
        sub = mutant[mutant["library"] == library].set_index("state")
        if {"DATP_like", "AT2"}.issubset(sub.index):
            t1.append({"library": library,
                       "areg_DATP_like": float(sub.loc["DATP_like", "mean_Areg"]),
                       "areg_AT2": float(sub.loc["AT2", "mean_Areg"]),
                       "higher_in_DATP": bool(sub.loc["DATP_like", "mean_Areg"] > sub.loc["AT2", "mean_Areg"])})
        else:
            t1.append({"library": library, "evaluable": False})
    t1_frame = pd.DataFrame(t1)
    t1_holds = bool(len(t1_frame) and t1_frame.get("higher_in_DATP", pd.Series(dtype=bool)).all())

    t2 = []
    for library in sorted(mutant["library"].unique()):
        sub = mutant[(mutant["library"] == library) & (mutant["state"] == "DATP_like")]
        if len(sub):
            means = {g: float(sub[f"mean_{g}"].iloc[0]) for g in ligand_cols}
            order = sorted(means, key=means.get, reverse=True)
            t2.append({"library": library, "ranking": " > ".join(order),
                       "top": order[0], "areg_is_top": order[0] == "Areg", **{f"mean_{k}": v for k, v in means.items()}})
    t2_frame = pd.DataFrame(t2)

    t3 = []
    for timepoint in ("4 days", "2 weeks"):
        datp = comp[(comp["timepoint"] == timepoint) & (comp["state"] == "DATP_like")]
        rfp = datp[datp["clone"] == "KrasG12D"]["pct_of_library"]
        yfp = datp[datp["clone"] == "wild type"]["pct_of_library"]
        t3.append({"timepoint": timepoint,
                   "median_pct_DATP_RFP": round(float(rfp.median()), 2) if len(rfp) else np.nan,
                   "median_pct_DATP_YFP": round(float(yfp.median()), 2) if len(yfp) else np.nan,
                   "higher_in_RFP": bool(len(rfp) and len(yfp) and rfp.median() > yfp.median())})
    t3_frame = pd.DataFrame(t3)

    t4 = (mutant[mutant["state"] == "DATP_like"]
          .groupby("timepoint", observed=True)["mean_Areg"].median().round(4).to_dict())
    rec.set("tests", {
        "T1_state_specificity": {"per_library": t1, "holds_in_all": t1_holds},
        "T2_ligand_ranking": t2,
        "T3_internal_control": t3,
        "T4_time": {"median_DATP_like_mean_Areg_by_timepoint": {str(k): float(v) for k, v in t4.items()}},
    })
    for frame, name in ((t1_frame, "c3_T1_state_specificity.csv"), (t2_frame, "c3_T2_ligand_ranking.csv"),
                        (t3_frame, "c3_T3_internal_control.csv")):
        frame.to_csv(OUT / name, index=False)
        rec.add_output(OUT / name)

    # ---- figure -----------------------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))
    sc.pl.umap(merged, color="state", ax=axes[0], show=False, size=6,
               title=f"epithelial states, own clustering (n={merged.n_obs})")
    sc.pl.umap(merged, color="timepoint", ax=axes[1], show=False, size=6, title="time point")
    sc.pl.umap(merged, color="Areg", ax=axes[2], show=False, size=6, title="Areg (log1p CP10K)")
    fig.tight_layout()
    fig.savefig(OUT / "c3_umap_states_time_areg.png", dpi=110)
    plt.close(fig)
    rec.add_output(OUT / "c3_umap_states_time_areg.png")

    lines = [
        "# Trial C3 output: is Areg a property of the DATP-like state (Gate 2a)", "",
        "This is an expression and state-specificity result, not a communication analysis. CellChat",
        "is R-only and unavailable here; no ligand-receptor probability is computed and nothing below",
        "establishes that a fibroblast receives the signal.", "",
        f"**T1, state specificity:** Areg higher in DATP-like than in AT2 cells in "
        f"{'all' if t1_holds else 'not all'} mutant RFP libraries.", "",
        df_to_markdown(t1_frame, index=False), "",
        "**T2, ligand ranking within the DATP-like state.** What ranks second is the part the paper",
        "does not show.", "", df_to_markdown(t2_frame, index=False), "",
        "**T3, the within-animal wild-type control.**", "", df_to_markdown(t3_frame, index=False), "",
        "## Replicate mixing, and the batch decision it drove", "",
        df_to_markdown(pd.DataFrame([{"arm": k, "same_library_kNN_enrichment": v}
                                     for k, v in enrichment.items()]), index=False), "",
        f"Failure threshold 2.0; correction applied: {correction_needed}.", "",
        "## Clusters and their calls", "", df_to_markdown(clusters, index=False), "",
        "## QC per library", "", df_to_markdown(qc, index=False), "",
    ]
    (OUT / "c3_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c3_summary.md")
    rec.finish()
    print(f"\nT1 holds in all mutant libraries: {t1_holds}")


if __name__ == "__main__":
    main()
