#!/usr/bin/env python
"""Trial D3: does AT2 reach AT1 through primed AT2 and DATP (PAGA and DPT)?

This is the one part of the paper a deposit with one library per condition
can genuinely re-ask, because an ordering inside a library needs no
between-group replication. It is asked twice: on all Tomato-positive cells
pooled (as the paper did) and on the day-14 library alone (a single
library, two pooled mice, no pooling across conditions).

Frozen rules, set before the annotated object was opened:

* Cells: the corrected annotated Tomato-positive object from trial D2b (D2's
  embedding, D2's thresholds, the Sftpc clauses applied only where Sftpc
  separates clusters), alveolar states only;
  cycling AT2 excluded from pseudotime, as the paper did, but kept in PAGA to
  test the paper's statement that cAT2 sits closest to pAT2.
* Graph: trial D2's PCA; neighbours recomputed at k = 15 (the paper's k) and,
  for robustness, k = 30; PAGA on the state labels; diffusion map with 15
  components; diffusion pseudotime with default parameters.
* Root: the hAT2 cell with the highest score for the paper's AT2 canonical
  and identity genes (score_genes, seed 0), chosen before any pseudotime is
  seen.
* Readings (each with the paper's statement it answers):
  O1 PAGA: hAT2's strongest partner is pAT2 (Figure 1E).
  O2 PAGA: the chain hAT2, pAT2, DATP, AT1 is connected, meaning each of the
     three consecutive pairs is among that node's two strongest edges.
  O3 PAGA: cAT2's strongest partner is pAT2 ("cAT2 assigned closest to pAT2").
  O4 DPT, pooled: median pseudotime rises hAT2 < pAT2 < DATP < AT1 (Figure 1F).
  O5 DPT, day 14 only: the same order.
  O6 robustness: O4 and O5 hold at both k values.
* Floor: a state with fewer than 30 cells in a subset makes any reading that
  names it not computable (the repository's 30-cell floor), reported as such.
* No P value; the ordering is a description of one deposit. The pooled
  ordering mixes three libraries and is reported as the paper reported it,
  with the single-library ordering beside it as the cleaner version.
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
from choi_utils import (DERIVED, STATE_ORDER, STATE_SETS, RunRecord, apply_style,  # noqa: E402
                        df_to_markdown, state_colour)

OUT = HERE / "d3_ordering"
OUT.mkdir(exist_ok=True)
K_VALUES = [15, 30]
FLOOR = 30
RULES = {
    "question": "Is the AT2 to pAT2 to DATP to AT1 ordering recoverable by PAGA and diffusion pseudotime, pooled and within day 14 alone",
    "cells": "D2b corrected annotated Tomato+ object (tomato_annotated_d2b.h5ad), alveolar states; cAT2 excluded from DPT, kept in PAGA",
    "graph": {"pca": "from D2", "k": K_VALUES, "diffmap_comps": 15, "dpt": "default"},
    "root": "hAT2 cell with the highest score_genes(hAT2_canonical + AT2_identity), seed 0",
    "readings": {"O1": "hAT2 strongest PAGA partner is pAT2", "O2": "chain hAT2-pAT2-DATP-AT1 connected (each consecutive pair within the node's top two edges)",
                 "O3": "cAT2 strongest partner is pAT2", "O4": "pooled DPT medians hAT2 < pAT2 < DATP < AT1",
                 "O5": "day-14-only DPT medians in the same order", "O6": "O4 and O5 hold at k=15 and k=30"},
    "floor_cells_per_state": FLOOR,
    "no_p_values": True,
}


def order_ok(med: pd.Series) -> bool:
    seq = ["hAT2", "pAT2", "DATP", "AT1"]
    if any(s not in med.index for s in seq):
        return False
    vals = [med[s] for s in seq]
    return all(vals[i] < vals[i + 1] for i in range(3))


def run_subset(adata: ad.AnnData, label: str, k: int, rec: RunRecord) -> dict:
    sub = adata[adata.obs["state"].isin(STATE_ORDER)].copy()
    sub.obs["state"] = pd.Categorical(sub.obs["state"].astype(str), categories=[s for s in STATE_ORDER if s in set(sub.obs["state"])])
    counts = sub.obs["state"].value_counts()
    sc.pp.neighbors(sub, n_neighbors=k, use_rep="X_pca", random_state=0)
    sc.tl.paga(sub, groups="state")
    conn = pd.DataFrame(sub.uns["paga"]["connectivities"].toarray(),
                        index=sub.obs["state"].cat.categories, columns=sub.obs["state"].cat.categories)
    conn.to_csv(OUT / f"d3_paga_{label}_k{k}.csv")
    rec.add_output(OUT / f"d3_paga_{label}_k{k}.csv")

    def strongest(node, exclude=()):
        row = conn.loc[node].drop(labels=[node, *exclude], errors="ignore")
        return row.idxmax() if len(row) else None

    def top2(node):
        row = conn.loc[node].drop(labels=[node], errors="ignore").sort_values(ascending=False)
        return list(row.index[:2])

    res = {"label": label, "k": k, "n_by_state": counts.to_dict()}
    low = [s for s in STATE_ORDER if s in counts.index and counts[s] < FLOOR]
    res["below_floor"] = low
    res["O1"] = (strongest("hAT2") == "pAT2") if {"hAT2", "pAT2"} <= set(conn.index) and not {"hAT2", "pAT2"} & set(low) else None
    chain = ["hAT2", "pAT2", "DATP", "AT1"]
    if all(s in conn.index for s in chain) and not set(chain) & set(low):
        res["O2"] = all(chain[i + 1] in top2(chain[i]) or chain[i] in top2(chain[i + 1]) for i in range(3))
    else:
        res["O2"] = None
    res["O3"] = (strongest("cAT2") == "pAT2") if {"cAT2", "pAT2"} <= set(conn.index) and not {"cAT2", "pAT2"} & set(low) else None

    # pseudotime without cycling cells
    dpt = sub[sub.obs["state"] != "cAT2"].copy()
    dpt.obs["state"] = pd.Categorical(dpt.obs["state"].astype(str))
    sc.tl.score_genes(dpt, [g for g in STATE_SETS["hAT2_canonical"] + STATE_SETS["AT2_identity"] if g in dpt.var_names],
                      score_name="root_score", ctrl_size=50, n_bins=25, random_state=0)
    h = np.where((dpt.obs["state"] == "hAT2").to_numpy())[0]
    if len(h) == 0:
        res["O4_O5"] = None
        return res
    dpt.uns["iroot"] = int(h[np.argmax(dpt.obs["root_score"].to_numpy()[h])])
    sc.pp.neighbors(dpt, n_neighbors=k, use_rep="X_pca", random_state=0)
    sc.tl.diffmap(dpt, n_comps=15)
    sc.tl.dpt(dpt)
    d = dpt.obs.groupby("state", observed=True)["dpt_pseudotime"]
    table = pd.DataFrame({"median": d.median(), "q25": d.quantile(0.25), "q75": d.quantile(0.75), "n": d.size()})
    table = table.reindex([s for s in STATE_ORDER if s in table.index])
    table.to_csv(OUT / f"d3_dpt_{label}_k{k}.csv")
    rec.add_output(OUT / f"d3_dpt_{label}_k{k}.csv")
    usable = table[table["n"] >= FLOOR]["median"]
    res["dpt_medians"] = table["median"].round(4).to_dict()
    res["order_met"] = order_ok(usable) if all(s in usable.index for s in ("hAT2", "pAT2", "DATP", "AT1")) else None
    res["_dpt_obj"] = dpt
    return res


def main() -> None:
    rec = RunRecord(OUT / "d3_run_record.json", "D3 Choi-2020 ordering by PAGA and DPT", RULES)
    src = DERIVED / "tomato_annotated_d2b.h5ad"
    rec.add_input(src)
    adata = ad.read_h5ad(src)
    results = []
    keep_for_fig = {}
    for k in K_VALUES:
        r_pool = run_subset(adata, "pooled", k, rec)
        r_d14 = run_subset(adata[adata.obs["library"] == "Day14_AT2_Tomato"].copy(), "day14", k, rec)
        keep_for_fig[k] = (r_pool.pop("_dpt_obj", None), r_d14.pop("_dpt_obj", None))
        results += [r_pool, r_d14]
    summary_rows = []
    for r in results:
        summary_rows.append({"subset": r["label"], "k": r["k"], "O1": r.get("O1"), "O2": r.get("O2"), "O3": r.get("O3"),
                             "order_met": r.get("order_met"), "dpt_medians": r.get("dpt_medians"),
                             "below_floor": ",".join(r.get("below_floor", [])) or "none"})
    res_df = pd.DataFrame(summary_rows)
    res_df.to_csv(OUT / "d3_readings.csv", index=False)
    rec.add_output(OUT / "d3_readings.csv")
    pooled15 = next(r for r in results if r["label"] == "pooled" and r["k"] == 15)
    d1415 = next(r for r in results if r["label"] == "day14" and r["k"] == 15)
    readings = {
        "O1": pooled15.get("O1"), "O2": pooled15.get("O2"), "O3": pooled15.get("O3"),
        "O4": pooled15.get("order_met"), "O5": d1415.get("order_met"),
        "O6": all(r.get("order_met") for r in results) if all(r.get("order_met") is not None for r in results) else None,
    }
    rec.set("readings", readings)

    # figure: DPT by state, pooled and day 14, k = 15
    apply_style(plt)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for axx, (label, obj) in zip(axes, [("pooled", keep_for_fig[15][0]), ("day 14 only", keep_for_fig[15][1])]):
        if obj is None:
            continue
        states = [s for s in ["hAT2", "pAT2", "DATP", "AT1"] if s in set(obj.obs["state"])]
        data = [obj.obs.loc[obj.obs["state"] == s, "dpt_pseudotime"].to_numpy() for s in states]
        parts = axx.violinplot(data, showmedians=True, showextrema=False)
        for body, s in zip(parts["bodies"], states):
            body.set_facecolor(state_colour(s)); body.set_alpha(0.8)
        axx.set_xticks(range(1, len(states) + 1)); axx.set_xticklabels([f"{s}\n(n={len(dd)})" for s, dd in zip(states, data)])
        axx.set_title(f"Diffusion pseudotime by state, {label}, k = 15")
        axx.set_ylabel("DPT (root: hAT2)")
    fig.tight_layout()
    fig.savefig(OUT / "d3_dpt_by_state.png", dpi=150)
    rec.add_output(OUT / "d3_dpt_by_state.png")
    plt.close(fig)

    lines = ["# Trial D3: ordering by PAGA and diffusion pseudotime", "",
             "Readings at k = 15 (pooled unless stated): " + ", ".join(f"{k}: {v}" for k, v in readings.items()), "",
             df_to_markdown(res_df, index=False)]
    (OUT / "d3_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "d3_summary.md")
    rec.finish()
    print("[D3] done", readings, flush=True)


if __name__ == "__main__":
    main()
