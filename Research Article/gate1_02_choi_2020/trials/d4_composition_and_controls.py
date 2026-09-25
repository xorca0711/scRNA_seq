#!/usr/bin/env python
"""Trial D4: the time-course composition, and three controls on it.

Trial D2 reads the paper's composition claims from the Tomato-positive
libraries, and trial D2b corrects D2's annotation without moving a
threshold; this trial reads D2b's object and composition. This trial asks whether those readings survive three things the
deposit itself warns about: the cell-calling rule (six raw whitelists), the
twofold spread in genes per cell across libraries (trial D0), and the
lineage label (the Tomato-negative libraries are the same lungs' cells that
the AT2 lineage did not label).

Frozen rules, set before any object was opened:

* Calling sensitivity (C1): for each Tomato-positive library, re-read the raw
  matrix, take the barcodes that pass the plain floor (500 counts, 200
  genes) but not the paper's filter, and count how many carry three or more
  of the five DATP markers. Reading: if adding them would move the library's
  DATP fraction by more than 0.02 in absolute terms, the DATP fraction is
  calling-sensitive; otherwise it is not.
* Depth window (C2): recompute the D2b composition using only cells with 1,000
  to 2,000 detected genes, a window that spans the three libraries' medians
  (1,061; 1,825; 1,689). Reading: D2b's R4 and R5, recomputed inside the
  window, hold or do not.
* Tomato-negative libraries (C3): each non-lineage library is embedded alone
  (same settings as D2; Leiden 0.5), clusters gated to alveolar epithelium
  (Sftpc detection at least 0.5, or AT1-canonical at least 0.5, or DATP at
  least 0.4, and no contaminant set at 0.5 or above), and the corrected
  annotation rule (D2b) applied to those clusters. Readings: the fraction of each library
  that is alveolar epithelium; whether a DATP cluster exists in the day-14
  Tomato-negative library, or DATP markers are dispersed-present there
  (at least 2 percent with three or more markers). A transcriptome cannot say
  whether an unlabelled DATP came from an unlabelled AT2 cell or from
  another origin; that sentence is part of the reading.
* No P value; all libraries are single, two mice pooled each.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from choi_utils import (ANNOTATION_THRESHOLDS, CONTAMINANT_SETS, DERIVED, FLOOR_FILTER,  # noqa: E402
                        NONTOMATO_LIBRARIES, PAPER_FILTER, STATE_ORDER, STATE_SETS,
                        TOMATO_LIBRARIES, RunRecord, annotate_clusters_b,
                        cluster_detection_table, df_to_markdown, dispersed_presence,
                        library_paths, read_mtx_triplet)

OUT = HERE / "d4_composition_and_controls"
OUT.mkdir(exist_ok=True)
WINDOW = (1000, 2000)
RULES = {
    "question": "Do the D2b composition readings survive the calling rule, a depth window, and the lineage label",
    "C1_calling_sensitivity": {"floor": FLOOR_FILTER, "paper": PAPER_FILTER, "datp_markers_min": 3, "move_threshold_abs": 0.02},
    "C2_depth_window_genes": WINDOW,
    "C3_nontomato": {"embedding": "as D2, per library, Leiden 0.5",
                     "alveolar_gate": "Sftpc>=0.5 or AT1_canonical>=0.5 or DATP>=0.4, and no contaminant set >=0.5",
                     "dispersed": {"markers": 3, "fraction": 0.02}},
    "no_p_values": True,
}


def calling_sensitivity(rec: RunRecord, comp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for e in TOMATO_LIBRARIES:
        paths = library_paths("GSE145031", e["library"], e["gsm"])
        for p in paths.values():
            rec.add_input(p)
        X, var, bc = read_mtx_triplet(paths["matrix"], paths["features"], paths["barcodes"])
        counts = np.asarray(X.sum(axis=1)).ravel()
        genes = X.getnnz(axis=1)
        paper = (genes > PAPER_FILTER["min_genes_exclusive"]) & (genes < PAPER_FILTER["max_genes_exclusive"]) & (counts > PAPER_FILTER["min_counts_exclusive"])
        floor = (counts >= FLOOR_FILTER["min_counts"]) & (genes >= FLOOR_FILTER["min_genes"])
        extra = floor & ~paper
        sym = var["gene_symbol"].to_numpy()
        idx = [int(np.where(sym == g)[0][0]) for g in STATE_SETS["DATP"] if (sym == g).any()]
        sub = X[extra][:, idx]
        k = np.asarray((sub > 0).sum(axis=1)).ravel()
        n_extra = int(extra.sum())
        n_datp_like = int((k >= RULES["C1_calling_sensitivity"]["datp_markers_min"]).sum())
        alv = float(comp.loc[e["library"], "alveolar_cells"])
        datp_n = float(comp.loc[e["library"], "DATP_n"]) if "DATP_n" in comp.columns else 0.0
        before = datp_n / alv if alv else np.nan
        after = (datp_n + n_datp_like) / (alv + n_extra) if alv + n_extra else np.nan
        rows.append({"library": e["library"], "floor_only_barcodes": n_extra, "median_counts_extra": float(np.median(counts[extra])) if n_extra else np.nan,
                     "datp_like_among_extra": n_datp_like, "datp_fraction_D2b": before, "datp_fraction_if_added": after,
                     "abs_change": abs(after - before) if np.isfinite(after) and np.isfinite(before) else np.nan})
        del X
    return pd.DataFrame(rows)


def depth_window(adata: ad.AnnData) -> pd.DataFrame:
    lo, hi = WINDOW
    m = (adata.obs["n_genes"] >= lo) & (adata.obs["n_genes"] <= hi) & adata.obs["state"].isin(STATE_ORDER)
    ct = pd.crosstab(adata.obs.loc[m, "library"], adata.obs.loc[m, "state"])
    frac = ct.div(ct.sum(axis=1), axis=0)
    out = ct.add_suffix("_n").join(frac.add_suffix("_frac"))
    out["cells_in_window"] = ct.sum(axis=1)
    return out.reindex([e["library"] for e in TOMATO_LIBRARIES])


def nontomato(rec: RunRecord) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, tables = [], []
    for e in NONTOMATO_LIBRARIES:
        p = DERIVED / f"{e['library']}.h5ad"
        rec.add_input(p)
        a = ad.read_h5ad(p)
        a = a[~a.obs["predicted_doublet"].astype(bool)].copy()
        a.layers["counts"] = a.X.copy()
        sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
        sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000, layer="counts")
        emb = a[:, a.var["highly_variable"]].copy(); sc.pp.scale(emb, max_value=10); sc.tl.pca(emb, n_comps=30, random_state=0)
        a.obsm["X_pca"] = emb.obsm["X_pca"]
        sc.pp.neighbors(a, n_neighbors=30, n_pcs=30, use_rep="X_pca", random_state=0)
        sc.tl.leiden(a, resolution=0.5, key_added="leiden_0.5", flavor="igraph", n_iterations=2, random_state=0)
        table = cluster_detection_table(a, "leiden_0.5")
        contaminated = table[list(CONTAMINANT_SETS)].max(axis=1) >= ANNOTATION_THRESHOLDS["contaminant_min"]
        alveolar = ((table["Sftpc"] >= 0.5) | (table["AT1_canonical"] >= 0.5) | (table["DATP"] >= 0.4)) & ~contaminated
        gated = annotate_clusters_b(table[alveolar].copy()) if alveolar.any() else table[alveolar].copy()
        gated["library"] = e["library"]
        tables.append(gated)
        state_map = dict(zip(gated["cluster"], gated["state"])) if len(gated) else {}
        a.obs["state"] = a.obs["leiden_0.5"].astype(str).map(state_map).fillna("non_alveolar")
        alv_mask = a.obs["state"].isin(STATE_ORDER).to_numpy()
        frac_datp, n_datp = dispersed_presence(a, STATE_SETS["DATP"], 3, alv_mask) if alv_mask.any() else (np.nan, 0)
        comp = a.obs.loc[alv_mask, "state"].value_counts()
        rows.append({"library": e["library"], "cells": int(a.n_obs), "alveolar_cells": int(alv_mask.sum()),
                     "alveolar_fraction": float(alv_mask.mean()),
                     **{f"{s}_n": int(comp.get(s, 0)) for s in STATE_ORDER},
                     "DATP_cluster_present": bool((gated["state"] == "DATP").any()) if len(gated) else False,
                     "DATP_dispersed_fraction": frac_datp, "DATP_dispersed_n": n_datp})
    return pd.DataFrame(rows), pd.concat(tables) if tables else pd.DataFrame()


def main() -> None:
    rec = RunRecord(OUT / "d4_run_record.json", "D4 Choi-2020 composition and its controls", RULES)
    src = DERIVED / "tomato_annotated_d2b.h5ad"
    rec.add_input(src)
    adata = ad.read_h5ad(src)
    comp = pd.read_csv(HERE / "d2b_corrected_annotation" / "d2b_composition.csv", index_col=0)
    rec.add_input(HERE / "d2b_corrected_annotation" / "d2b_composition.csv")

    c1 = calling_sensitivity(rec, comp)
    c1.to_csv(OUT / "d4_calling_sensitivity.csv", index=False); rec.add_output(OUT / "d4_calling_sensitivity.csv")
    c2 = depth_window(adata)
    c2.to_csv(OUT / "d4_depth_window_composition.csv"); rec.add_output(OUT / "d4_depth_window_composition.csv")
    c3, gated = nontomato(rec)
    c3.to_csv(OUT / "d4_nontomato_composition.csv", index=False); rec.add_output(OUT / "d4_nontomato_composition.csv")
    gated.to_csv(OUT / "d4_nontomato_cluster_annotation.csv", index=False); rec.add_output(OUT / "d4_nontomato_cluster_annotation.csv")

    def g(tab, lib, s):
        col = f"{s}_frac"
        return float(tab.loc[lib, col]) if col in tab.columns and lib in tab.index and np.isfinite(tab.loc[lib, col]) else 0.0
    r4_window = (g(c2, "Day14_AT2_Tomato", "hAT2") < 0.5 * g(c2, "PBS_AT2_Tomato", "hAT2")
                 and g(c2, "Day14_AT2_Tomato", "pAT2") >= 0.05 and g(c2, "Day14_AT2_Tomato", "DATP") >= 0.05)
    r5_window = (g(c2, "Day28_AT2_Tomato", "AT1") > g(c2, "Day14_AT2_Tomato", "AT1")
                 and g(c2, "Day28_AT2_Tomato", "hAT2") > g(c2, "Day14_AT2_Tomato", "hAT2")
                 and all(g(c2, "Day28_AT2_Tomato", s) < g(c2, "Day14_AT2_Tomato", s) for s in ("cAT2", "pAT2", "DATP")))
    d14 = c3.set_index("library").loc["Day14_AT2_nonTomato"] if "Day14_AT2_nonTomato" in set(c3["library"]) else None
    readings = {
        "C1_calling_sensitive": {lib: bool(ch > RULES["C1_calling_sensitivity"]["move_threshold_abs"]) if np.isfinite(ch) else None
                                 for lib, ch in zip(c1["library"], c1["abs_change"])},
        "C2_R4_in_window": bool(r4_window), "C2_R5_in_window": bool(r5_window),
        "C3_datp_in_day14_nontomato": (bool(d14["DATP_cluster_present"]) or (np.isfinite(d14["DATP_dispersed_fraction"]) and d14["DATP_dispersed_fraction"] >= 0.02)) if d14 is not None else None,
    }
    rec.set("readings", readings)
    lines = ["# Trial D4: composition and its controls", "", "Readings: " + str(readings), "",
             "## C1 calling sensitivity", "", df_to_markdown(c1.round(4), index=False), "",
             "## C2 composition inside the 1,000 to 2,000 gene window", "", df_to_markdown(c2.round(4)), "",
             "## C3 Tomato-negative libraries", "", df_to_markdown(c3.round(4), index=False), "",
             "A transcriptome cannot say whether an unlabelled DATP-like cell came from an AT2 cell the lineage did not label or from another origin; the count is reported, the origin is not."]
    (OUT / "d4_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "d4_summary.md")
    rec.finish()
    print("[D4] done", readings, flush=True)


if __name__ == "__main__":
    main()
