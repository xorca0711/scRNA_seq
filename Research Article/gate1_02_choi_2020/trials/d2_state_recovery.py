#!/usr/bin/env python
"""Trial D2: are the paper's five AT2-lineage states recoverable blind (Gate 1)?

The three Tomato-positive libraries (PBS, day 14, day 28) are merged without
batch correction, because each library is one condition and the library is
the condition: correcting on library would delete the time course. Same-
library neighbour enrichment is measured and reported, not corrected. The
deposit carries no author labels, so nothing is held out; the answer key is
the paper's text (states, their markers, the negative condition, the
qualitative composition per time point) and the grading is against that.

Frozen rules, set before any object was opened (the annotation rule and its
thresholds live in choi_utils.ANNOTATION_THRESHOLDS and were written before
this trial ran):

* Cells: trial D1 objects, predicted doublets removed.
* Embedding: normalise to 10,000, log1p, 2,000 highly variable genes
  (seurat_v3 on counts), scale to 10, 30 principal components, 30 neighbours,
  Leiden (igraph flavour, seed 0) at 0.3, 0.5 and 1.0, UMAP seed 0.
* Annotation: choi_utils.annotate_clusters at each resolution, cluster level,
  from detection fractions of the paper's own marker sets. Contaminant
  clusters (ciliated, mesenchyme, immune, endothelium, airway) are counted
  and removed from every composition.
* Primary resolution: the lowest of the three at which all five states are
  assigned to at least one cluster. If none, the resolution assigning the
  most states, and the missing state is checked for dispersed presence (at
  least 2 percent of alveolar cells carrying three or more of its markers):
  present-but-unresolved is reported as such, never as absent.
* Readings, each with the paper's statement it answers:
  R1 five states present as clusters at the primary resolution (Figure 1B).
  R2 cycling AT2 between 3 and 9 percent of alveolar Tomato-positive cells
     across the three libraries pooled ("approximately 6%").
  R3 PBS: hAT2 is the largest state and exceeds 50 percent (Figure 1C).
  R4 day 14: hAT2 fraction below half of its PBS fraction, and pAT2 and DATP
     each at least 5 percent ("dramatically reduced", "three additional
     populations").
  R5 day 28: AT1 and hAT2 fractions above their day-14 values, and cAT2, pAT2
     and DATP each below their day-14 values ("return to homeostasis").
  R6 the DATP cluster's AT1-canonical detection below half the AT1 cluster's
     (the negative condition).
  R7 same-library kNN enrichment reported; a value above 2 would have
     triggered the repository's batch rule, but no correction is admissible
     here because library equals condition, so the value is a description.
* No P value. One library per condition; every composition number is a
  description of one library that pools two mice.
* What makes R2 to R5 unreadable: if the primary resolution assigns fewer
  than four of the five states, R2 to R5 are reported as not computable
  rather than scored against a partial vocabulary.
"""

from __future__ import annotations

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
                        STATE_SETS, TOMATO_LIBRARIES, RunRecord, annotate_clusters,
                        apply_style, cluster_detection_table, detection,
                        df_to_markdown, dispersed_presence, same_library_enrichment,
                        state_colour)

OUT = HERE / "d2_state_recovery"
OUT.mkdir(exist_ok=True)
RULES = {
    "question": "Are hAT2, cAT2, pAT2, DATP and AT1 recoverable from the Tomato-positive libraries with the paper's markers, and does their composition follow the paper's time course",
    "cells": "trial D1 objects, predicted doublets removed",
    "embedding": {"normalize_total": 1e4, "log1p": True, "hvg": "seurat_v3, 2000, on counts",
                  "scale_max": 10, "n_pcs": 30, "n_neighbors": 30, "leiden": RESOLUTIONS, "seed": 0},
    "batch": "no correction: library equals condition; same-library kNN enrichment reported",
    "annotation": {"rule": "choi_utils.annotate_clusters", "thresholds": ANNOTATION_THRESHOLDS},
    "primary_resolution_rule": "lowest resolution assigning all five states; else the one assigning most, with a dispersed-presence check for the missing state",
    "readings": {
        "R1": "five states present as clusters at the primary resolution",
        "R2": "cAT2 between 0.03 and 0.09 of alveolar Tomato+ cells pooled",
        "R3": "PBS: hAT2 largest and above 0.50",
        "R4": "day 14: hAT2 below 0.5 x PBS; pAT2 and DATP each at least 0.05",
        "R5": "day 28: AT1 and hAT2 above day 14; cAT2, pAT2, DATP below day 14",
        "R6": "DATP AT1-canonical detection below 0.5 x AT1 cluster's",
        "R7": "same-library kNN enrichment, description only",
    },
    "unreadable_if": "fewer than four states assigned at the primary resolution makes R2 to R5 not computable",
    "no_p_values": "one library per condition, two mice pooled per library",
}


def load() -> ad.AnnData:
    parts = []
    for e in TOMATO_LIBRARIES:
        a = ad.read_h5ad(DERIVED / f"{e['library']}.h5ad")
        a = a[~a.obs["predicted_doublet"].astype(bool)].copy()
        parts.append(a)
    merged = ad.concat(parts, join="inner", label="_src", index_unique="-")
    merged.layers["counts"] = merged.X.copy()
    return merged


def embed(adata: ad.AnnData) -> None:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=30, random_state=0)
    adata.obsm["X_pca"] = emb.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep="X_pca", random_state=0)
    for r in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph", n_iterations=2, random_state=0)
    sc.tl.umap(adata, random_state=0)


def main() -> None:
    rec = RunRecord(OUT / "d2_run_record.json", "D2 Choi-2020 state recovery (Gate 1)", RULES)
    for e in TOMATO_LIBRARIES:
        rec.add_input(DERIVED / f"{e['library']}.h5ad")
    adata = load()
    print(f"[D2] {adata.n_obs} cells, {adata.n_vars} genes", flush=True)
    embed(adata)
    enrich = same_library_enrichment(adata, "library")
    rec.set("same_library_knn_enrichment", enrich)

    annotations = {}
    states_assigned = {}
    for r in RESOLUTIONS:
        key = f"leiden_{r}"
        table = annotate_clusters(cluster_detection_table(adata, key))
        table.to_csv(OUT / f"d2_cluster_annotation_res{r}.csv", index=False)
        rec.add_output(OUT / f"d2_cluster_annotation_res{r}.csv")
        annotations[r] = table
        states_assigned[r] = sorted(set(table["state"]) & set(STATE_ORDER))
        adata.obs[f"state_{r}"] = adata.obs[key].astype(str).map(dict(zip(table["cluster"], table["state"])))
    full = [r for r in RESOLUTIONS if len(states_assigned[r]) == 5]
    primary = full[0] if full else max(RESOLUTIONS, key=lambda r: len(states_assigned[r]))
    rec.set("states_assigned_per_resolution", {str(r): v for r, v in states_assigned.items()})
    rec.set("primary_resolution", primary)
    adata.obs["state"] = adata.obs[f"state_{primary}"]
    adata.obs["cluster"] = adata.obs[f"leiden_{primary}"].astype(str)

    alveolar = adata.obs["state"].isin(STATE_ORDER).to_numpy()
    missing = [s for s in STATE_ORDER if s not in states_assigned[primary]]
    dispersed = {}
    for s in missing:
        genes = STATE_SETS["DATP"] if s == "DATP" else STATE_SETS["cAT2"] if s == "cAT2" else \
            STATE_SETS["AT1_canonical"] if s == "AT1" else STATE_SETS["hAT2_canonical"]
        frac, n = dispersed_presence(adata, genes, ANNOTATION_THRESHOLDS["dispersed_markers"], alveolar)
        dispersed[s] = {"fraction": frac, "n": n,
                        "verdict": "present but not resolved by clustering" if frac >= ANNOTATION_THRESHOLDS["dispersed_present_fraction"] else "not found"}
    rec.set("missing_states_dispersed_check", dispersed)

    # composition per library, alveolar states only
    comp_n = pd.crosstab(adata.obs["library"], adata.obs["state"])
    contaminants = comp_n[[c for c in comp_n.columns if c not in STATE_ORDER]] if any(c not in STATE_ORDER for c in comp_n.columns) else pd.DataFrame(index=comp_n.index)
    comp_alv = comp_n[[c for c in STATE_ORDER if c in comp_n.columns]]
    comp_frac = comp_alv.div(comp_alv.sum(axis=1), axis=0)
    order = [e["library"] for e in TOMATO_LIBRARIES]
    comp_alv = comp_alv.reindex(order)
    comp_frac = comp_frac.reindex(order)
    comp = comp_alv.add_suffix("_n").join(comp_frac.add_suffix("_frac"))
    comp["alveolar_cells"] = comp_alv.sum(axis=1)
    comp["contaminant_cells"] = contaminants.reindex(order).sum(axis=1) if len(contaminants.columns) else 0
    comp.to_csv(OUT / "d2_composition.csv")
    rec.add_output(OUT / "d2_composition.csv")
    contaminants.reindex(order).to_csv(OUT / "d2_contaminants.csv")
    rec.add_output(OUT / "d2_contaminants.csv")

    # readings
    f = comp_frac.fillna(0.0)
    def g(lib, s):
        return float(f.loc[lib, s]) if s in f.columns else 0.0
    pooled = comp_alv.sum(axis=0)
    cat2_pooled = float(pooled.get("cAT2", 0) / max(pooled.sum(), 1))
    computable = len(states_assigned[primary]) >= 4
    readings = {
        "R1": {"value": states_assigned[primary], "met": len(states_assigned[primary]) == 5},
        "R2": {"value": cat2_pooled, "met": (0.03 <= cat2_pooled <= 0.09) if computable else None},
        "R3": {"value": g("PBS_AT2_Tomato", "hAT2"),
               "met": (g("PBS_AT2_Tomato", "hAT2") > 0.5 and f.loc["PBS_AT2_Tomato"].idxmax() == "hAT2") if computable else None},
        "R4": {"value": {"hAT2_d14": g("Day14_AT2_Tomato", "hAT2"), "hAT2_pbs": g("PBS_AT2_Tomato", "hAT2"),
                         "pAT2_d14": g("Day14_AT2_Tomato", "pAT2"), "DATP_d14": g("Day14_AT2_Tomato", "DATP")},
               "met": (g("Day14_AT2_Tomato", "hAT2") < 0.5 * g("PBS_AT2_Tomato", "hAT2")
                       and g("Day14_AT2_Tomato", "pAT2") >= 0.05 and g("Day14_AT2_Tomato", "DATP") >= 0.05) if computable else None},
        "R5": {"value": {s: {"d14": g("Day14_AT2_Tomato", s), "d28": g("Day28_AT2_Tomato", s)} for s in STATE_ORDER},
               "met": (g("Day28_AT2_Tomato", "AT1") > g("Day14_AT2_Tomato", "AT1")
                       and g("Day28_AT2_Tomato", "hAT2") > g("Day14_AT2_Tomato", "hAT2")
                       and all(g("Day28_AT2_Tomato", s) < g("Day14_AT2_Tomato", s) for s in ("cAT2", "pAT2", "DATP"))) if computable else None},
        "R7": {"value": enrich, "met": None},
    }
    tab = annotations[primary].set_index("state")
    if "DATP" in tab.index and "AT1" in tab.index:
        datp_at1 = float(tab.loc[["DATP"], "AT1_canonical"].mean())
        at1_at1 = float(tab.loc[["AT1"], "AT1_canonical"].mean())
        readings["R6"] = {"value": {"DATP": datp_at1, "AT1": at1_at1}, "met": datp_at1 < 0.5 * at1_at1}
    else:
        readings["R6"] = {"value": None, "met": None}
    rec.set("readings", readings)
    pd.DataFrame([{"reading": k, "description": RULES["readings"][k], "value": str(v["value"]), "met": v["met"]}
                  for k, v in readings.items()]).to_csv(OUT / "d2_readings.csv", index=False)
    rec.add_output(OUT / "d2_readings.csv")

    # lipocalin ambiguity, per state
    lip = pd.DataFrame({s: detection(adata, ["Lcn1", "Lcn2"], (adata.obs["state"] == s).to_numpy())
                        for s in STATE_ORDER if (adata.obs["state"] == s).any()}).T
    lip.to_csv(OUT / "d2_lipocalin_detection_by_state.csv")
    rec.add_output(OUT / "d2_lipocalin_detection_by_state.csv")

    # figure: UMAP by state and by library
    apply_style(plt)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    xy = adata.obsm["X_umap"]
    for s in STATE_ORDER + [c for c in adata.obs["state"].unique() if c not in STATE_ORDER]:
        m = (adata.obs["state"] == s).to_numpy()
        if not m.any():
            continue
        axes[0].scatter(xy[m, 0], xy[m, 1], s=2, c=state_colour(s) if s in STATE_ORDER else "#cfcec8", linewidths=0, label=f"{s} ({m.sum()})")
        cx, cy = np.median(xy[m, 0]), np.median(xy[m, 1])
        axes[0].text(cx, cy, s, fontsize=8, ha="center", va="center")
    axes[0].set_title(f"Tomato-positive cells by state, Leiden {primary}")
    P = {"PBS_AT2_Tomato": "#2a78d6", "Day14_AT2_Tomato": "#eb6834", "Day28_AT2_Tomato": "#1baf7a"}
    for lib, col in P.items():
        m = (adata.obs["library"] == lib).to_numpy()
        axes[1].scatter(xy[m, 0], xy[m, 1], s=2, c=col, linewidths=0, label=lib.replace("_AT2_Tomato", ""))
    axes[1].legend(markerscale=5, fontsize=8)
    axes[1].set_title("by library (no batch correction)")
    for axx in axes:
        axx.set_xticks([]); axx.set_yticks([]); axx.set_xlabel("UMAP 1"); axx.set_ylabel("UMAP 2")
    fig.tight_layout()
    fig.savefig(OUT / "d2_umap_states_and_libraries.png", dpi=150)
    rec.add_output(OUT / "d2_umap_states_and_libraries.png")
    plt.close(fig)

    obj = DERIVED / "tomato_annotated.h5ad"
    adata.write_h5ad(obj)
    rec.set("annotated_object", str(obj))

    lines = [
        "# Trial D2: state recovery (Gate 1)",
        "",
        f"{adata.n_obs} Tomato-positive cells after D1; primary resolution Leiden {primary}; states assigned per resolution: "
        + "; ".join(f"{r}: {', '.join(states_assigned[r]) or 'none'}" for r in RESOLUTIONS) + ".",
        f"Same-library kNN enrichment {enrich:.3f} (1 = mixed; no correction admissible, library equals condition).",
        "",
        "## Readings",
        "",
        df_to_markdown(pd.read_csv(OUT / "d2_readings.csv"), index=False),
        "",
        "## Composition per library (alveolar states; fractions of alveolar cells)",
        "",
        df_to_markdown(comp.round(4)),
        "",
        "## Cluster annotation at the primary resolution",
        "",
        df_to_markdown(annotations[primary].round(3), index=False),
        "",
        "## The printed lipocalin: detection of Lcn1 and Lcn2 by state",
        "",
        df_to_markdown(lip.round(3)),
    ]
    if dispersed:
        lines += ["", "## Missing states, dispersed-presence check", "", str(dispersed)]
    (OUT / "d2_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "d2_summary.md")
    rec.finish()
    print("[D2] done; readings:", {k: v["met"] for k, v in readings.items()}, flush=True)


if __name__ == "__main__":
    main()
