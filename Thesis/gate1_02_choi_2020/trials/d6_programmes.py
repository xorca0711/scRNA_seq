#!/usr/bin/env python
"""Trial D6: the programmes the paper attaches to DATPs, and the responder gene.

The paper reads DATPs as enriched for p53 signalling, proliferation arrest,
hypoxia and the interferon-gamma response (Figure S1D), and for glycolysis
genes and Hif1a (Figure 4F), in vivo and in organoids. Ndrg1 sits in both the
DATP marker set and the hypoxia programme, so the hypoxia reading is also
taken without it. Il1r1, the responder marker, is reported by state; the
paper's statement that Il1r1 rises in G2/M-phase cycling cells (Figure S5G)
is not attempted here, because assigning cell-cycle phase needs an external
phase gene list this repository would otherwise have to take from memory.

Frozen rules, set before any object was opened:

* Objects: the corrected annotated objects of trials D2b (Tomato-positive)
  and D5b (organoids), alveolar states only.
* Scores: scanpy score_genes per programme (ctrl_size 50, n_bins 25, seed 0)
  on the log-normalised data; and per-gene detection fractions from counts.
* Readings, on the in vivo object (organoid reported beside it):
  P1 DATP has the highest mean score among the five states for each of p53,
     arrest, hypoxia and interferon-gamma response.
  P2 DATP has the highest mean score for glycolysis and the highest Hif1a
     detection.
  P3 P1's hypoxia reading survives removing Ndrg1 (Hif1a alone).
  P4 Il1r1 detection by state, reported; the G2/M claim is not attempted.
  P5 depth control: P1 and P2 recomputed inside the 1,000 to 2,000 gene window.
* Floor: 30 cells per state.
* No P value.
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
from choi_utils import (DERIVED, DERIVED_ORGANOID, STATE_ORDER, STATE_SETS, RunRecord,  # noqa: E402
                        apply_style, detection, df_to_markdown, palette)

OUT = HERE / "d6_programmes"
OUT.mkdir(exist_ok=True)
PROGRAMMES = ["p53", "arrest", "hypoxia", "hypoxia_without_Ndrg1", "ifng_response", "glycolysis"]
WINDOW = (1000, 2000)
FLOOR = 30
RULES = {
    "question": "Are the paper's DATP programmes (p53, arrest, hypoxia, interferon-gamma, glycolysis, Hif1a) highest in the DATP state, in vivo and in organoids",
    "scores": {"method": "score_genes ctrl_size 50, n_bins 25, seed 0", "sets": {k: STATE_SETS[k] for k in PROGRAMMES}},
    "readings": {"P1": "DATP highest mean score for p53, arrest, hypoxia, ifng_response", "P2": "DATP highest glycolysis score and Hif1a detection",
                 "P3": "hypoxia reading holds without Ndrg1", "P4": "Il1r1 detection by state, reported only; G2/M claim not attempted",
                 "P5": "P1 and P2 inside the 1,000 to 2,000 gene window"},
    "floor": FLOOR, "no_p_values": True,
}


def score_table(adata: ad.AnnData, mask: np.ndarray | None = None) -> pd.DataFrame:
    a = adata if mask is None else adata[mask].copy()
    for name in PROGRAMMES:
        genes = [g for g in STATE_SETS[name] if g in a.var_names]
        sc.tl.score_genes(a, genes, score_name=f"score_{name}", ctrl_size=50, n_bins=25, random_state=0)
    rows = []
    for s in STATE_ORDER:
        m = (a.obs["state"] == s).to_numpy()
        if m.sum() == 0:
            continue
        row = {"state": s, "n": int(m.sum())}
        for name in PROGRAMMES:
            row[f"score_{name}"] = float(a.obs.loc[m, f"score_{name}"].mean())
        d = detection(a, ["Hif1a", "Il1r1", "Ndrg1", "Pgk1", "Pkm", "Slc16a3"], m)
        for g, v in d.items():
            row[f"det_{g}"] = float(v)
        rows.append(row)
    return pd.DataFrame(rows)


def readings_from(tab: pd.DataFrame) -> dict:
    t = tab[tab["n"] >= FLOOR].set_index("state")
    if "DATP" not in t.index:
        return {"P1": None, "P2": None, "P3": None, "note": "DATP below floor or absent"}
    p1 = {name: bool(t[f"score_{name}"].idxmax() == "DATP") for name in ("p53", "arrest", "hypoxia", "ifng_response")}
    p2 = {"glycolysis": bool(t["score_glycolysis"].idxmax() == "DATP"), "Hif1a_detection": bool(t["det_Hif1a"].idxmax() == "DATP")}
    p3 = bool(t["score_hypoxia_without_Ndrg1"].idxmax() == "DATP")
    return {"P1": p1, "P1_all": all(p1.values()), "P2": p2, "P2_all": all(p2.values()), "P3": p3}


def main() -> None:
    rec = RunRecord(OUT / "d6_run_record.json", "D6 Choi-2020 programmes and the responder gene", RULES)
    out = {}
    for label, path in (("invivo", DERIVED / "tomato_annotated_d2b.h5ad"), ("organoid", DERIVED_ORGANOID / "organoid_annotated_d5b.h5ad")):
        rec.add_input(path)
        a = ad.read_h5ad(path)
        a = a[a.obs["state"].isin(STATE_ORDER)].copy()
        tab = score_table(a)
        tab.to_csv(OUT / f"d6_scores_{label}.csv", index=False); rec.add_output(OUT / f"d6_scores_{label}.csv")
        out[label] = {"scores": tab, "readings": readings_from(tab)}
        if label == "invivo":
            win = ((a.obs["n_genes"] >= WINDOW[0]) & (a.obs["n_genes"] <= WINDOW[1])).to_numpy()
            tab_w = score_table(a, win)
            tab_w.to_csv(OUT / "d6_scores_invivo_window.csv", index=False); rec.add_output(OUT / "d6_scores_invivo_window.csv")
            out[label]["readings"]["P5"] = readings_from(tab_w)
        il1r1 = detection(a, ["Il1r1"], None)
        out[label]["readings"]["P4_Il1r1_by_state"] = tab.set_index("state")["det_Il1r1"].round(4).to_dict() if "det_Il1r1" in tab else None
        out[label]["readings"]["P4_G2M_claim"] = "not attempted (needs an external phase gene list)"
    rec.set("readings", {k: v["readings"] for k, v in out.items()})

    apply_style(plt)
    P = palette()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for axx, (label, d) in zip(axes, out.items()):
        tab = d["scores"].set_index("state")
        M = tab[[f"score_{p}" for p in PROGRAMMES]]
        M.columns = [p.replace("_", " ") for p in PROGRAMMES]
        M = (M - M.min()) / (M.max() - M.min() + 1e-12)
        im = axx.imshow(M.to_numpy(), aspect="auto", cmap=matplotlib.colors.LinearSegmentedColormap.from_list("ramp", P["sequential_ramp"]))
        axx.set_xticks(range(M.shape[1])); axx.set_xticklabels(M.columns, rotation=40, ha="right", fontsize=8)
        axx.set_yticks(range(M.shape[0])); axx.set_yticklabels([f"{s} (n={int(tab.loc[s, 'n'])})" for s in M.index], fontsize=8)
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                axx.text(j, i, f"{tab.iloc[i][f'score_{PROGRAMMES[j]}']:.2f}", ha="center", va="center", fontsize=7,
                         color=P["ink"] if M.iloc[i, j] < 0.6 else P["surface"])
        axx.set_title(f"Mean programme score by state, {label} (rows scaled 0 to 1; numbers are raw means)", fontsize=9)
    fig.tight_layout(); fig.savefig(OUT / "d6_programme_scores.png", dpi=150); rec.add_output(OUT / "d6_programme_scores.png"); plt.close(fig)

    lines = ["# Trial D6: programmes and the responder gene", ""]
    for label, d in out.items():
        lines += [f"## {label}", "", "Readings: " + str(d["readings"]), "", df_to_markdown(d["scores"].round(4), index=False), ""]
    (OUT / "d6_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "d6_summary.md")
    rec.finish()
    print("[D6] done", {k: {kk: vv for kk, vv in v["readings"].items() if kk in ("P1_all", "P2_all", "P3")} for k, v in out.items()}, flush=True)


if __name__ == "__main__":
    main()
