#!/usr/bin/env python
"""Trial D5b: the organoid deposit under the corrected annotation rule, beside D5.

D5 applied D2's first-pass rule and inherited its defect: every organoid
cluster detects Sftpc in 100 percent of cells, so the Sftpc-low clauses
never fired. A 586-cell control cluster with mesenchyme detection 0.70 (the
size of the stromal cluster the paper removed: 1,868 minus 1,286 is 582) was
kept as hAT2, a 365-cell cluster at 0.64 likewise, and a 481-cell cluster
with AT1 canonical detection 0.79 was labelled hAT2 instead of AT1, which
made reading G5 not computable. D5's outcome stays in the record; this pass
applies choi_utils.annotate_clusters_b (the Sftpc clauses only where Sftpc
separates clusters; every threshold unchanged) to D5's embedding and
recomputes D5's readings G1 to G6 with D5's own functions. The corrected
labels go to a second object (organoid_annotated_d5b.h5ad); trial D6 reads
it and says so.
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
from choi_utils import (ANNOTATION_THRESHOLDS, DERIVED_ORGANOID, RESOLUTIONS, STATE_ORDER,  # noqa: E402
                        STATE_SETS, RunRecord, annotate_clusters_b, apply_style,
                        cluster_detection_table, df_to_markdown, state_colour)
from d5_organoids import FLOOR, WINDOW, at1_ratios, composition, g5_reading  # noqa: E402

OUT = HERE / "d5b_corrected_annotation"
OUT.mkdir(exist_ok=True)
RULES = {
    "question": "D5's readings under the corrected annotation rule (Sftpc clauses applied only where Sftpc separates clusters); thresholds unchanged",
    "first_pass": "trials/d5_organoids/, outcome kept",
    "rule": "choi_utils.annotate_clusters_b", "thresholds": ANNOTATION_THRESHOLDS,
    "readings": "G1 to G6 as D5, computed with D5's functions", "window_genes": WINDOW, "floor": FLOOR, "no_p_values": True,
}


def main() -> None:
    rec = RunRecord(OUT / "d5b_run_record.json", "D5b Choi-2020 organoids, corrected annotation", RULES)
    src = DERIVED_ORGANOID / "organoid_annotated.h5ad"; rec.add_input(src)
    adata = ad.read_h5ad(src)
    adata.obs["state_first_pass"] = adata.obs["state"]
    states_assigned, tables = {}, {}
    for r in RESOLUTIONS:
        t = annotate_clusters_b(cluster_detection_table(adata, f"leiden_{r}"))
        t.to_csv(OUT / f"d5b_cluster_annotation_res{r}.csv", index=False); rec.add_output(OUT / f"d5b_cluster_annotation_res{r}.csv")
        tables[r] = t; states_assigned[r] = sorted(set(t["state"]) & set(STATE_ORDER))
        adata.obs[f"state_b_{r}"] = adata.obs[f"leiden_{r}"].astype(str).map(dict(zip(t["cluster"], t["state"])))
    full = [r for r in RESOLUTIONS if len(states_assigned[r]) == 5]
    primary = full[0] if full else max(RESOLUTIONS, key=lambda r: len(states_assigned[r]))
    adata.obs["state"] = adata.obs[f"state_b_{primary}"]
    rec.set("primary_resolution", primary); rec.set("states_assigned", {str(k): v for k, v in states_assigned.items()})
    contaminants = adata.obs.loc[~adata.obs["state"].isin(STATE_ORDER)].groupby(["treatment", "state"], observed=True).size()
    contaminants.to_csv(OUT / "d5b_contaminants.csv"); rec.add_output(OUT / "d5b_contaminants.csv")

    comp = composition(adata); comp.to_csv(OUT / "d5b_composition.csv"); rec.add_output(OUT / "d5b_composition.csv")
    win = ((adata.obs["n_genes"] >= WINDOW[0]) & (adata.obs["n_genes"] <= WINDOW[1])).to_numpy()
    comp_w = composition(adata, win); comp_w.to_csv(OUT / "d5b_composition_window.csv"); rec.add_output(OUT / "d5b_composition_window.csv")

    def f(tab, t, s):
        col = f"{s}_frac"
        return float(tab.loc[t, col]) if t in tab.index and col in tab.columns and np.isfinite(tab.loc[t, col]) else 0.0

    epi = adata[adata.obs["state"].isin(STATE_ORDER) & (adata.obs["state"] != "cAT2")].copy()
    g4 = None
    if (epi.obs["state"] == "hAT2").sum() >= FLOOR:
        sc.tl.score_genes(epi, [g for g in STATE_SETS["hAT2_canonical"] + STATE_SETS["AT2_identity"] if g in epi.var_names],
                          score_name="root_score", ctrl_size=50, n_bins=25, random_state=0)
        cand = np.where(((epi.obs["state"] == "hAT2") & (epi.obs["treatment"] == "control")).to_numpy())[0]
        epi.uns["iroot"] = int(cand[np.argmax(epi.obs["root_score"].to_numpy()[cand])])
        sc.pp.neighbors(epi, n_neighbors=15, use_rep="X_pca", random_state=0)
        sc.tl.diffmap(epi, n_comps=15); sc.tl.dpt(epi)
        med = epi.obs.groupby("treatment", observed=True)["dpt_pseudotime"].median()
        g4 = bool(med.get("IL-1beta", np.nan) > med.get("control", np.nan))
        med.to_csv(OUT / "d5b_dpt_median_by_treatment.csv"); rec.add_output(OUT / "d5b_dpt_median_by_treatment.csv")
        epi.obs.groupby(["treatment", "state"], observed=True)["dpt_pseudotime"].median().to_csv(OUT / "d5b_dpt_median_by_treatment_state.csv")
        rec.add_output(OUT / "d5b_dpt_median_by_treatment_state.csv")

    ratios = at1_ratios(adata); ratios.to_csv(OUT / "d5b_at1_marker_ratios.csv", index=False); rec.add_output(OUT / "d5b_at1_marker_ratios.csv")
    ratios_w = at1_ratios(adata, win); ratios_w.to_csv(OUT / "d5b_at1_marker_ratios_window.csv", index=False); rec.add_output(OUT / "d5b_at1_marker_ratios_window.csv")
    n_at1 = adata.obs[adata.obs["state"] == "AT1"].groupby("treatment", observed=True).size()
    g5_ok = all(n_at1.get(t, 0) >= FLOOR for t in ("control", "IL-1beta"))
    readings = {
        "G1": len(states_assigned[primary]) == 5,
        "G2": bool(f(comp, "control", "hAT2") + f(comp, "control", "AT1") > 0.60 and f(comp, "control", "pAT2") < 0.20 and f(comp, "control", "DATP") < 0.10),
        "G3": bool(f(comp, "IL-1beta", "pAT2") >= 0.50 and f(comp, "IL-1beta", "DATP") > f(comp, "control", "DATP")),
        "G3_DATP_direction_only": bool(f(comp, "IL-1beta", "DATP") > f(comp, "control", "DATP")),
        "G4": g4,
        "G5": g5_reading(ratios) if g5_ok else {"not_computable": "fewer than 30 AT1 cells in a treatment", "n_at1": n_at1.to_dict()},
        "G5_n_at1": n_at1.to_dict(),
        "G6": {"G3_in_window": bool(f(comp_w, "IL-1beta", "pAT2") >= 0.50 and f(comp_w, "IL-1beta", "DATP") > f(comp_w, "control", "DATP")),
               "G3_DATP_direction_in_window": bool(f(comp_w, "IL-1beta", "DATP") > f(comp_w, "control", "DATP")),
               "G5_in_window": g5_reading(ratios_w) if g5_ok else None},
        "paper_epithelial_counts": {"control": 1286, "IL-1beta": 2584},
        "this_run_epithelial_counts": comp["epithelial_cells"].to_dict(),
    }
    rec.set("readings", readings)

    apply_style(plt)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    xy = adata.obsm["X_umap"]
    for s in sorted(adata.obs["state"].unique(), key=lambda x: (x not in STATE_ORDER, x)):
        m = (adata.obs["state"] == s).to_numpy()
        ax.scatter(xy[m, 0], xy[m, 1], s=3, c=state_colour(s) if s in STATE_ORDER else "#cfcec8", linewidths=0)
        ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), f"{s} ({m.sum()})", fontsize=7, ha="center")
    ax.set_title(f"Organoid cells, corrected annotation (D5b), Leiden {primary}"); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2")
    fig.tight_layout(); fig.savefig(OUT / "d5b_umap_states.png", dpi=150); rec.add_output(OUT / "d5b_umap_states.png"); plt.close(fig)

    obj = DERIVED_ORGANOID / "organoid_annotated_d5b.h5ad"
    adata.write_h5ad(obj); rec.set("annotated_object", str(obj))
    lines = ["# Trial D5b: organoids under the corrected annotation rule, beside D5", "",
             f"Primary resolution {primary}; states per resolution: " + "; ".join(f"{r}: {', '.join(states_assigned[r]) or 'none'}" for r in RESOLUTIONS), "",
             "Readings: " + str(readings), "",
             "## Contaminant clusters removed (cells per treatment)", "", df_to_markdown(contaminants.reset_index().rename(columns={0: "n"}), index=False) if len(contaminants) else "none", "",
             "## Composition (fractions of epithelial cells)", "", df_to_markdown(comp.round(4)), "",
             "## Composition inside the 2,500 to 4,500 gene window", "", df_to_markdown(comp_w.round(4)), "",
             "## AT1 cells: marker detection, IL-1beta over control", "", df_to_markdown(ratios.round(3), index=False), "",
             "## Cluster annotation at the primary resolution", "", df_to_markdown(tables[primary].round(3), index=False)]
    (OUT / "d5b_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "d5b_summary.md")
    rec.finish()
    print("[D5b] done", {k: readings[k] for k in ("G1", "G2", "G3", "G3_DATP_direction_only", "G4", "G5_n_at1")}, flush=True)


if __name__ == "__main__":
    main()
