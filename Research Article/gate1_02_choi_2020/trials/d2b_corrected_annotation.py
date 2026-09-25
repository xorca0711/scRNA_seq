#!/usr/bin/env python
"""Trial D2b: the corrected annotation pass, beside D2's first outcome.

D2's frozen rule returned "not recovered" (four of five states; AT1 missing)
and the cluster tables show why. In a lineage-sorted library every cluster
detects Sftpc, so the rule's "Sftpc below 0.5" clauses for AT1 and for
contaminants never fire: the cluster with AT1 canonical detection 0.70 was
called hAT2, and a ciliated cluster and an immune cluster (245 cells, against
the 255 the paper removed) stayed in, the ciliated one wearing the pAT2
label. That is a defect in the rule, not in the biology, and it is handled
the way trial C1b handled C1's: D2's outcome stays in the record, no
threshold moves, and this pass adds one condition on the condition (the
Sftpc clauses apply only where Sftpc separates clusters;
choi_utils.annotate_clusters_b).

Frozen for this pass, before the object was re-opened:

* Same embedding (D2's), same three resolutions, same thresholds, same
  primary-resolution rule, same readings R1 to R7 as D2.
* Two disclosed additions D2 did not have:
  (a) R2 is also reported on the day-14 library alone, because the paper's
      sentence ("approximately 6% of lineage-labeled cells expressed cell
      cycle markers") sits in its day-14 paragraph and can be read either
      way; the pooled reading is the one D2 froze and stays primary.
  (b) If pAT2 is assigned to no cluster at the primary resolution, the DATP
      cluster(s) are sub-clustered alone (Leiden 0.5, seed 0) and the rule is
      applied to the sub-clusters, since pAT2 and DATP are neighbours on the
      paper's own trajectory and a coarse cluster can hold both; and a
      cell-level dispersed check counts hAT2-labelled cells that detect none
      of the three identity genes and at least two of the five inflammatory
      genes. Sub-clustering is characterisation, not a new threshold.
* The corrected labels are written to a second object
  (tomato_annotated_d2b.h5ad); D2's object is untouched. Trials D3 to D7
  read the corrected object and say so.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from choi_utils import (ANNOTATION_THRESHOLDS, DERIVED, RESOLUTIONS, STATE_ORDER,  # noqa: E402
                        STATE_SETS, TOMATO_LIBRARIES, RunRecord, annotate_clusters_b,
                        apply_style, cluster_detection_table, detection, df_to_markdown,
                        state_colour)

OUT = HERE / "d2b_corrected_annotation"
OUT.mkdir(exist_ok=True)
RULES = {
    "question": "D2's readings under the corrected annotation rule (Sftpc clauses applied only where Sftpc separates clusters); thresholds unchanged",
    "first_pass": "trials/d2_state_recovery/, outcome kept",
    "rule": "choi_utils.annotate_clusters_b", "thresholds": ANNOTATION_THRESHOLDS,
    "additions": {"R2_day14": "reported beside the pooled R2", "subcluster_DATP_if_no_pAT2": "Leiden 0.5, seed 0, rule re-applied",
                  "cell_level_pAT2_check": "hAT2-labelled cells with 0 of 3 identity genes and >=2 of 5 inflammatory genes"},
    "readings": "R1 to R7 as D2",
}


def main() -> None:
    rec = RunRecord(OUT / "d2b_run_record.json", "D2b Choi-2020 corrected annotation", RULES)
    src = DERIVED / "tomato_annotated.h5ad"; rec.add_input(src)
    adata = ad.read_h5ad(src)
    annotations, states_assigned = {}, {}
    for r in RESOLUTIONS:
        key = f"leiden_{r}"
        table = annotate_clusters_b(cluster_detection_table(adata, key))
        table.to_csv(OUT / f"d2b_cluster_annotation_res{r}.csv", index=False); rec.add_output(OUT / f"d2b_cluster_annotation_res{r}.csv")
        annotations[r] = table
        states_assigned[r] = sorted(set(table["state"]) & set(STATE_ORDER))
        adata.obs[f"state_b_{r}"] = adata.obs[key].astype(str).map(dict(zip(table["cluster"], table["state"])))
    full = [r for r in RESOLUTIONS if len(states_assigned[r]) == 5]
    primary = full[0] if full else max(RESOLUTIONS, key=lambda r: len(states_assigned[r]))
    rec.set("states_assigned_per_resolution", {str(r): v for r, v in states_assigned.items()})
    rec.set("primary_resolution", primary)
    adata.obs["state_first_pass"] = adata.obs["state"]
    adata.obs["state"] = adata.obs[f"state_b_{primary}"]
    adata.obs["cluster"] = adata.obs[f"leiden_{primary}"].astype(str)

    # (b) sub-cluster DATP if pAT2 is missing
    sub_table = None
    if "pAT2" not in states_assigned[primary] and "DATP" in states_assigned[primary]:
        m = (adata.obs["state"] == "DATP").to_numpy()
        sub = adata[m].copy()
        sc.pp.neighbors(sub, n_neighbors=30, use_rep="X_pca", random_state=0)
        sc.tl.leiden(sub, resolution=0.5, key_added="sub", flavor="igraph", n_iterations=2, random_state=0)
        sub_table = annotate_clusters_b(cluster_detection_table(sub, "sub"))
        sub_table.to_csv(OUT / "d2b_datp_subclusters.csv", index=False); rec.add_output(OUT / "d2b_datp_subclusters.csv")
        newmap = dict(zip(sub_table["cluster"], sub_table["state"]))
        relabel = sub.obs["sub"].astype(str).map(newmap)
        adata.obs.loc[sub.obs_names, "state"] = relabel.where(relabel.isin(["pAT2", "DATP"]), "DATP").to_numpy()
        adata.obs.loc[sub.obs_names, "cluster"] = (adata.obs.loc[sub.obs_names, "cluster"].astype(str) + "." + sub.obs["sub"].astype(str)).to_numpy()
        rec.set("datp_subclustering", {"n": int(m.sum()), "sub_states": newmap})
    # cell-level pAT2 dispersed check among hAT2-labelled cells
    h = (adata.obs["state"] == "hAT2").to_numpy()
    X = adata.layers["counts"][h]
    idg = [adata.var_names.get_loc(g) for g in STATE_SETS["AT2_identity"] if g in adata.var_names]
    inf = [adata.var_names.get_loc(g) for g in STATE_SETS["pAT2_inflammatory"] if g in adata.var_names]
    k_id = np.asarray((X[:, idg] > 0).sum(axis=1)).ravel(); k_inf = np.asarray((X[:, inf] > 0).sum(axis=1)).ravel()
    pat2_like = (k_id == 0) & (k_inf >= 2)
    rec.set("cell_level_pAT2_like_among_hAT2", {"n": int(pat2_like.sum()), "fraction": float(pat2_like.mean())})

    # composition and readings, as D2
    comp_n = pd.crosstab(adata.obs["library"], adata.obs["state"])
    order = [e["library"] for e in TOMATO_LIBRARIES]
    alv_cols = [c for c in STATE_ORDER if c in comp_n.columns]
    comp_alv = comp_n[alv_cols].reindex(order)
    comp_frac = comp_alv.div(comp_alv.sum(axis=1), axis=0)
    comp = comp_alv.add_suffix("_n").join(comp_frac.add_suffix("_frac"))
    comp["alveolar_cells"] = comp_alv.sum(axis=1)
    cont_cols = [c for c in comp_n.columns if c not in STATE_ORDER]
    comp["contaminant_cells"] = comp_n[cont_cols].reindex(order).sum(axis=1) if cont_cols else 0
    comp.to_csv(OUT / "d2b_composition.csv"); rec.add_output(OUT / "d2b_composition.csv")
    comp_n[cont_cols].reindex(order).to_csv(OUT / "d2b_contaminants.csv"); rec.add_output(OUT / "d2b_contaminants.csv")
    f = comp_frac.fillna(0.0)
    def g(lib, s):
        return float(f.loc[lib, s]) if s in f.columns else 0.0
    pooled = comp_alv.sum(axis=0)
    cat2_pooled = float(pooled.get("cAT2", 0) / max(pooled.sum(), 1))
    cat2_d14 = g("Day14_AT2_Tomato", "cAT2")
    computable = len(states_assigned[primary]) >= 4
    readings = {
        "R1": {"value": states_assigned[primary], "met": len(states_assigned[primary]) == 5},
        "R2": {"value": cat2_pooled, "met": (0.03 <= cat2_pooled <= 0.09) if computable else None},
        "R2_day14_reported": {"value": cat2_d14, "met": (0.03 <= cat2_d14 <= 0.09) if computable else None},
        "R3": {"value": g("PBS_AT2_Tomato", "hAT2"), "met": (g("PBS_AT2_Tomato", "hAT2") > 0.5 and f.loc["PBS_AT2_Tomato"].idxmax() == "hAT2") if computable else None},
        "R4": {"value": {"hAT2_d14": g("Day14_AT2_Tomato", "hAT2"), "hAT2_pbs": g("PBS_AT2_Tomato", "hAT2"), "pAT2_d14": g("Day14_AT2_Tomato", "pAT2"), "DATP_d14": g("Day14_AT2_Tomato", "DATP")},
               "met": (g("Day14_AT2_Tomato", "hAT2") < 0.5 * g("PBS_AT2_Tomato", "hAT2") and g("Day14_AT2_Tomato", "pAT2") >= 0.05 and g("Day14_AT2_Tomato", "DATP") >= 0.05) if computable else None},
        "R5": {"value": {s: {"d14": g("Day14_AT2_Tomato", s), "d28": g("Day28_AT2_Tomato", s)} for s in STATE_ORDER},
               "met": (g("Day28_AT2_Tomato", "AT1") > g("Day14_AT2_Tomato", "AT1") and g("Day28_AT2_Tomato", "hAT2") > g("Day14_AT2_Tomato", "hAT2")
                       and all(g("Day28_AT2_Tomato", s) < g("Day14_AT2_Tomato", s) for s in ("cAT2", "pAT2", "DATP"))) if computable else None},
        "R7": {"value": json.load(open(HERE / "d2_state_recovery" / "d2_run_record.json", encoding="utf-8"))["results"]["same_library_knn_enrichment"],
               "met": "unchanged from D2 (same embedding)"},
    }
    tab = annotations[primary].set_index("state")
    if "DATP" in tab.index and "AT1" in tab.index:
        datp_at1 = float(tab.loc[["DATP"], "AT1_canonical"].mean()); at1_at1 = float(tab.loc[["AT1"], "AT1_canonical"].mean())
        readings["R6"] = {"value": {"DATP": datp_at1, "AT1": at1_at1}, "met": datp_at1 < 0.5 * at1_at1}
    else:
        readings["R6"] = {"value": None, "met": None}
    rec.set("readings", readings)
    pd.DataFrame([{"reading": k, "value": str(v["value"]), "met": v["met"]} for k, v in readings.items()]).to_csv(OUT / "d2b_readings.csv", index=False)
    rec.add_output(OUT / "d2b_readings.csv")

    apply_style(plt)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    xy = adata.obsm["X_umap"]
    for s in STATE_ORDER + sorted(c for c in adata.obs["state"].unique() if c not in STATE_ORDER):
        m = (adata.obs["state"] == s).to_numpy()
        if not m.any():
            continue
        ax.scatter(xy[m, 0], xy[m, 1], s=2, c=state_colour(s) if s in STATE_ORDER else "#cfcec8", linewidths=0)
        ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), f"{s} ({m.sum()})", fontsize=7, ha="center")
    ax.set_title(f"Corrected annotation (D2b), Leiden {primary}"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    fig.tight_layout(); fig.savefig(OUT / "d2b_umap_states.png", dpi=150); rec.add_output(OUT / "d2b_umap_states.png"); plt.close(fig)

    obj = DERIVED / "tomato_annotated_d2b.h5ad"
    adata.write_h5ad(obj); rec.set("annotated_object", str(obj))
    lines = ["# Trial D2b: corrected annotation, beside D2's first outcome", "",
             f"Primary resolution {primary}; states per resolution: " + "; ".join(f"{r}: {', '.join(states_assigned[r]) or 'none'}" for r in RESOLUTIONS), "",
             "## Readings", "", df_to_markdown(pd.read_csv(OUT / "d2b_readings.csv"), index=False), "",
             "## Composition per library (fractions of alveolar cells)", "", df_to_markdown(comp.round(4)), "",
             "## Contaminant clusters removed (cells per library)", "", df_to_markdown(comp_n[cont_cols].reindex(order)) if cont_cols else "none", "",
             "## Cluster annotation at the primary resolution", "", df_to_markdown(annotations[primary].round(3), index=False), ""]
    if sub_table is not None:
        lines += ["## DATP cluster sub-clustered for a primed AT2 sub-population", "", df_to_markdown(sub_table.round(3), index=False), ""]
    lines += ["## Cell-level primed-like check among hAT2-labelled cells", "", str(rec.record["results"]["cell_level_pAT2_like_among_hAT2"])]
    (OUT / "d2b_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "d2b_summary.md")
    rec.finish()
    print("[D2b] done; readings:", {k: v["met"] for k, v in readings.items()}, flush=True)


if __name__ == "__main__":
    main()
