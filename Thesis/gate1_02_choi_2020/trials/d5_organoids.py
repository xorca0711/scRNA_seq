#!/usr/bin/env python
"""Trial D5: the organoid deposit, control against IL-1beta.

Two libraries, one per treatment, deposited already filtered. The paper
reports five clusters matching the in vivo states, control organoids mostly
AT2 and AT1 with small pAT2 and DATP fractions, IL-1beta raising the primed
fraction to about 77 percent and the DATP fraction too, pseudotime skewed
toward AT1, and, within AT1-like cells, early AT1 markers comparable, late
AT1 markers lower and DATP genes higher under IL-1beta (Figures 2F to 2H,
7A to 7C). All of that is direction, not test: one library per treatment.

Frozen rules, set before any object was opened:

* Cells: trial D1 organoid objects, predicted doublets removed; the two
  libraries merged without correction (treatment equals library).
* Embedding and annotation exactly as trial D2 (same rule, same thresholds,
  same primary-resolution rule). The stromal cells the paper co-cultured are
  expected to fall into a contaminant cluster (mesenchyme set) and are counted
  against the paper's exclusion.
* Readings:
  G1 five states present at the primary resolution (Figure 2F).
  G2 control: hAT2 plus AT1 above 0.60 of epithelial cells, pAT2 below 0.20,
     DATP below 0.10 (Figure 2G).
  G3 IL-1beta: pAT2 at least 0.50 (paper about 0.77) and DATP fraction above
     control's (Figure 2G).
  G4 DPT (root: the control hAT2 cell with the highest AT2 score) has a higher
     median in IL-1beta cells than in control cells (Figure 2H).
  G5 within AT1 cells, IL-1beta over control detection ratios: early AT1
     between 0.7 and 1.4; late AT1 below 0.7; DATP figure-7 genes above 1.4
     (Figures 7A to 7C). Each gene reported; the reading is on the set mean.
  G6 depth control for G3 and G5: the IL-1beta library has fewer genes per
     cell (D0: 3,193 against 3,598 medians); G3 and G5 are recomputed inside
     a 2,500 to 4,500 gene window.
* Floor: 30 cells per state per library for any reading that names it.
* No P value; two libraries, no replication.
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
from choi_utils import (ANNOTATION_THRESHOLDS, DERIVED_ORGANOID, ORGANOID_LIBRARIES,  # noqa: E402
                        RESOLUTIONS, STATE_ORDER, STATE_SETS, RunRecord, annotate_clusters,
                        apply_style, cluster_detection_table, detection, df_to_markdown,
                        state_colour)

OUT = HERE / "d5_organoids"
OUT.mkdir(exist_ok=True)
WINDOW = (2500, 4500)
FLOOR = 30
RULES = {
    "question": "Do the organoid libraries carry the five states, and does IL-1beta move the composition, the pseudotime and the AT1 markers the way the paper reports",
    "cells": "D1 organoid objects, doublets removed; merged without correction",
    "embedding_and_annotation": "as D2; thresholds " + str(ANNOTATION_THRESHOLDS),
    "readings": {"G1": "five states present", "G2": "control: hAT2+AT1 > 0.60, pAT2 < 0.20, DATP < 0.10",
                 "G3": "IL-1beta: pAT2 >= 0.50 and DATP > control", "G4": "median DPT IL-1beta > control, root control hAT2",
                 "G5": "AT1 cells, IL-1beta/control detection: early 0.7 to 1.4, late < 0.7, DATP figure-7 genes > 1.4",
                 "G6": "G3 and G5 inside a 2,500 to 4,500 gene window"},
    "window_genes": WINDOW, "floor": FLOOR, "no_p_values": True,
}


def embed(adata):
    sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = adata[:, adata.var["highly_variable"]].copy(); sc.pp.scale(emb, max_value=10); sc.tl.pca(emb, n_comps=30, random_state=0)
    adata.obsm["X_pca"] = emb.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep="X_pca", random_state=0)
    for r in RESOLUTIONS:
        sc.tl.leiden(adata, resolution=r, key_added=f"leiden_{r}", flavor="igraph", n_iterations=2, random_state=0)
    sc.tl.umap(adata, random_state=0)


def composition(adata, mask=None):
    obs = adata.obs if mask is None else adata.obs[mask]
    obs = obs[obs["state"].isin(STATE_ORDER)]
    ct = pd.crosstab(obs["treatment"], obs["state"])
    frac = ct.div(ct.sum(axis=1), axis=0)
    out = ct.add_suffix("_n").join(frac.add_suffix("_frac"))
    out["epithelial_cells"] = ct.sum(axis=1)
    return out


def at1_ratios(adata, mask=None):
    rows = []
    m_all = adata.obs["state"].eq("AT1").to_numpy() if mask is None else (adata.obs["state"].eq("AT1").to_numpy() & mask)
    for setname in ("AT1_early", "AT1_late", "DATP_figure7"):
        for g in STATE_SETS[setname]:
            vals = {}
            for t in ("control", "IL-1beta"):
                m = m_all & (adata.obs["treatment"] == t).to_numpy()
                d = detection(adata, [g], m) if m.sum() else pd.Series(dtype=float)
                vals[t] = float(d.get(g, np.nan)); vals[f"n_{t}"] = int(m.sum())
            ratio = vals["IL-1beta"] / vals["control"] if vals.get("control", 0) > 0 else np.nan
            rows.append({"set": setname, "gene": g, **vals, "ratio_il1b_over_control": ratio})
    return pd.DataFrame(rows)


def g5_reading(tab: pd.DataFrame) -> dict:
    out = {}
    means = tab.groupby("set")["ratio_il1b_over_control"].mean()
    out["early"] = bool(0.7 <= means.get("AT1_early", np.nan) <= 1.4) if np.isfinite(means.get("AT1_early", np.nan)) else None
    out["late"] = bool(means.get("AT1_late", np.nan) < 0.7) if np.isfinite(means.get("AT1_late", np.nan)) else None
    out["datp"] = bool(means.get("DATP_figure7", np.nan) > 1.4) if np.isfinite(means.get("DATP_figure7", np.nan)) else None
    out["set_mean_ratios"] = means.round(3).to_dict()
    return out


def main() -> None:
    rec = RunRecord(OUT / "d5_run_record.json", "D5 Choi-2020 organoids, control against IL-1beta", RULES)
    parts = []
    for e in ORGANOID_LIBRARIES:
        p = DERIVED_ORGANOID / f"{e['library']}.h5ad"; rec.add_input(p)
        a = ad.read_h5ad(p); a = a[~a.obs["predicted_doublet"].astype(bool)].copy(); parts.append(a)
    adata = ad.concat(parts, join="inner", label="_src", index_unique="-")
    adata.layers["counts"] = adata.X.copy()
    embed(adata)
    states_assigned, tables = {}, {}
    for r in RESOLUTIONS:
        t = annotate_clusters(cluster_detection_table(adata, f"leiden_{r}"))
        t.to_csv(OUT / f"d5_cluster_annotation_res{r}.csv", index=False); rec.add_output(OUT / f"d5_cluster_annotation_res{r}.csv")
        tables[r] = t; states_assigned[r] = sorted(set(t["state"]) & set(STATE_ORDER))
        adata.obs[f"state_{r}"] = adata.obs[f"leiden_{r}"].astype(str).map(dict(zip(t["cluster"], t["state"])))
    full = [r for r in RESOLUTIONS if len(states_assigned[r]) == 5]
    primary = full[0] if full else max(RESOLUTIONS, key=lambda r: len(states_assigned[r]))
    adata.obs["state"] = adata.obs[f"state_{primary}"]
    rec.set("primary_resolution", primary); rec.set("states_assigned", {str(k): v for k, v in states_assigned.items()})
    contaminants = adata.obs.loc[~adata.obs["state"].isin(STATE_ORDER)].groupby(["treatment", "state"], observed=True).size()
    contaminants.to_csv(OUT / "d5_contaminants.csv"); rec.add_output(OUT / "d5_contaminants.csv")

    comp = composition(adata); comp.to_csv(OUT / "d5_composition.csv"); rec.add_output(OUT / "d5_composition.csv")
    win = ((adata.obs["n_genes"] >= WINDOW[0]) & (adata.obs["n_genes"] <= WINDOW[1])).to_numpy()
    comp_w = composition(adata, win); comp_w.to_csv(OUT / "d5_composition_window.csv"); rec.add_output(OUT / "d5_composition_window.csv")

    def f(tab, t, s):
        col = f"{s}_frac"
        return float(tab.loc[t, col]) if t in tab.index and col in tab.columns and np.isfinite(tab.loc[t, col]) else 0.0

    # G4 pseudotime
    epi = adata[adata.obs["state"].isin(STATE_ORDER) & (adata.obs["state"] != "cAT2")].copy()
    g4 = None
    if (epi.obs["state"] == "hAT2").sum() >= FLOOR:
        sc.tl.score_genes(epi, [g for g in STATE_SETS["hAT2_canonical"] + STATE_SETS["AT2_identity"] if g in epi.var_names],
                          score_name="root_score", ctrl_size=50, n_bins=25, random_state=0)
        cand = np.where(((epi.obs["state"] == "hAT2") & (epi.obs["treatment"] == "control")).to_numpy())[0]
        if len(cand) == 0:
            cand = np.where((epi.obs["state"] == "hAT2").to_numpy())[0]
        epi.uns["iroot"] = int(cand[np.argmax(epi.obs["root_score"].to_numpy()[cand])])
        sc.pp.neighbors(epi, n_neighbors=15, use_rep="X_pca", random_state=0)
        sc.tl.diffmap(epi, n_comps=15); sc.tl.dpt(epi)
        med = epi.obs.groupby("treatment", observed=True)["dpt_pseudotime"].median()
        g4 = bool(med.get("IL-1beta", np.nan) > med.get("control", np.nan)) if all(k in med.index for k in ("control", "IL-1beta")) else None
        med.to_csv(OUT / "d5_dpt_median_by_treatment.csv"); rec.add_output(OUT / "d5_dpt_median_by_treatment.csv")
        epi.obs.groupby(["treatment", "state"], observed=True)["dpt_pseudotime"].median().to_csv(OUT / "d5_dpt_median_by_treatment_state.csv")
        rec.add_output(OUT / "d5_dpt_median_by_treatment_state.csv")

    ratios = at1_ratios(adata); ratios.to_csv(OUT / "d5_at1_marker_ratios.csv", index=False); rec.add_output(OUT / "d5_at1_marker_ratios.csv")
    ratios_w = at1_ratios(adata, win); ratios_w.to_csv(OUT / "d5_at1_marker_ratios_window.csv", index=False); rec.add_output(OUT / "d5_at1_marker_ratios_window.csv")
    n_at1 = adata.obs[adata.obs["state"] == "AT1"].groupby("treatment", observed=True).size()
    g5_ok = all(n_at1.get(t, 0) >= FLOOR for t in ("control", "IL-1beta"))

    readings = {
        "G1": len(states_assigned[primary]) == 5,
        "G2": bool(f(comp, "control", "hAT2") + f(comp, "control", "AT1") > 0.60 and f(comp, "control", "pAT2") < 0.20 and f(comp, "control", "DATP") < 0.10),
        "G3": bool(f(comp, "IL-1beta", "pAT2") >= 0.50 and f(comp, "IL-1beta", "DATP") > f(comp, "control", "DATP")),
        "G4": g4,
        "G5": g5_reading(ratios) if g5_ok else {"not_computable": "fewer than 30 AT1 cells in a treatment", "n_at1": n_at1.to_dict()},
        "G6": {"G3_in_window": bool(f(comp_w, "IL-1beta", "pAT2") >= 0.50 and f(comp_w, "IL-1beta", "DATP") > f(comp_w, "control", "DATP")),
               "G5_in_window": g5_reading(ratios_w) if g5_ok else None},
        "paper_epithelial_counts": {"control": 1286, "IL-1beta": 2584},
        "this_run_epithelial_counts": comp["epithelial_cells"].to_dict() if "epithelial_cells" in comp else None,
    }
    rec.set("readings", readings)

    apply_style(plt)
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    xy = adata.obsm["X_umap"]
    for s in sorted(adata.obs["state"].unique(), key=lambda x: (x not in STATE_ORDER, x)):
        m = (adata.obs["state"] == s).to_numpy()
        axes[0].scatter(xy[m, 0], xy[m, 1], s=3, c=state_colour(s) if s in STATE_ORDER else "#cfcec8", linewidths=0)
        axes[0].text(np.median(xy[m, 0]), np.median(xy[m, 1]), f"{s}", fontsize=8, ha="center")
    axes[0].set_title(f"Organoid cells by state, Leiden {primary}")
    for t, col in (("control", "#5598e7"), ("IL-1beta", "#184f95")):
        m = (adata.obs["treatment"] == t).to_numpy()
        axes[1].scatter(xy[m, 0], xy[m, 1], s=3, c=col, linewidths=0, label=f"{t} ({m.sum()})")
    axes[1].legend(markerscale=4, fontsize=8); axes[1].set_title("by treatment (no correction)")
    for axx in axes:
        axx.set_xticks([]); axx.set_yticks([]); axx.set_xlabel("UMAP 1"); axx.set_ylabel("UMAP 2")
    fig.tight_layout(); fig.savefig(OUT / "d5_umap_states_and_treatment.png", dpi=150); rec.add_output(OUT / "d5_umap_states_and_treatment.png"); plt.close(fig)

    adata.write_h5ad(DERIVED_ORGANOID / "organoid_annotated.h5ad")
    lines = ["# Trial D5: organoids, control against IL-1beta", "", f"Primary resolution {primary}; states per resolution: " +
             "; ".join(f"{r}: {', '.join(states_assigned[r]) or 'none'}" for r in RESOLUTIONS), "", "Readings: " + str(readings), "",
             "## Composition (fractions of epithelial cells)", "", df_to_markdown(comp.round(4)), "",
             "## Composition inside the 2,500 to 4,500 gene window", "", df_to_markdown(comp_w.round(4)), "",
             "## AT1 cells: marker detection, IL-1beta over control", "", df_to_markdown(ratios.round(3), index=False), "",
             "## Cluster annotation at the primary resolution", "", df_to_markdown(tables[primary].round(3), index=False)]
    (OUT / "d5_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "d5_summary.md")
    rec.finish()
    print("[D5] done", {k: readings[k] for k in ("G1", "G2", "G3", "G4")}, flush=True)


if __name__ == "__main__":
    main()
