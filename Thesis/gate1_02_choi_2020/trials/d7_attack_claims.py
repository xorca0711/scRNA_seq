#!/usr/bin/env python
"""Trial D7: the attacks. Four ways the paper's two new states could be
artefacts, each with a rule frozen before the object was opened.

A1 Is primed AT2 a depth artefact? The state is defined by LOSS of Etv5,
   Abca3 and Cebpa. A cell with fewer counts detects fewer genes, so a
   "loss" can be a detection floor. Rule: compare median genes per cell in
   pAT2 and hAT2; then match hAT2 cells to the pAT2 genes-per-cell
   distribution (20 quantile bins, sampling without replacement, seed 0)
   and recompute the identity-set detection in matched hAT2 against pAT2.
   Reading: if the pAT2 over matched-hAT2 identity detection ratio is at
   most 0.7, the loss is not a depth artefact; if it is above 0.7 the state
   is depth-confounded. Not computable if either state is below 30 cells.
A2 Is DATP a doublet cluster? Rule: mean Scrublet score per state; the
   paper's own cluster cut (0.6) and the repository's flagged fraction. And
   the AT2-plus-AT1 co-detection that a doublet would carry: fraction of
   cells detecting Sftpc AND two or more AT1 canonical genes, by state.
   Reading: DATP mean score above 0.6 or a double-marker fraction above
   hAT2's and AT1's both means doublet-like; otherwise not.
A3 Dissociation stress: NOT attempted. It needs the published
   dissociation-induced gene list (van den Brink et al. 2017, DOI
   10.1038/nmeth.4437), which is in that paper's supplement and not on disk;
   the repository does not take gene lists from memory.
A4 Is DATP a discrete co-expressing population or independent noise? Rule:
   in day-14 Tomato-positive alveolar cells, Krt8 and Cldn4 co-detection
   against independence (observed over expected), with a null that permutes
   Cldn4 detection within deciles of total counts (200 permutations, seed 0),
   as trial C9 did. Reading: ratio above 1.25 with permutation p below 0.05
   means co-organised beyond depth; ratio between 0.8 and 1.25 means no
   evidence of a discrete co-expressing population; the ratio is reported
   with the double-positive count, and a count below 30 makes the magnitude
   unreadable (existence only).

All four read D2b's corrected annotated object. No P value crosses
libraries; A4's permutation is inside one library.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from choi_utils import DERIVED, STATE_ORDER, STATE_SETS, RunRecord, detection, df_to_markdown  # noqa: E402

OUT = HERE / "d7_attack_claims"
OUT.mkdir(exist_ok=True)
FLOOR = 30
RULES = {
    "A1": {"identity_set": STATE_SETS["AT2_identity"], "match": "20 quantile bins of n_genes, without replacement, seed 0",
           "reading": "pAT2 / matched hAT2 identity detection ratio <= 0.7 means not a depth artefact"},
    "A2": {"paper_cluster_cut": 0.6, "double_marker": "Sftpc detected AND >=2 of Pdpn, Hopx, Cav1 detected"},
    "A3": "not attempted: needs the published dissociation gene list (van den Brink 2017, 10.1038/nmeth.4437), not on disk",
    "A4": {"genes": ["Krt8", "Cldn4"], "library": "Day14_AT2_Tomato", "permutations": 200, "seed": 0,
           "reading": "ratio > 1.25 and p < 0.05 co-organised; 0.8 to 1.25 no evidence; count < 30 magnitude unreadable"},
    "floor": FLOOR,
}


def genes_matrix(adata, genes, mask):
    X = adata.layers["counts"][mask]
    cols = [adata.var_names.get_loc(g) for g in genes if g in adata.var_names]
    return (X[:, cols] > 0).toarray()


def a1(adata, rec):
    obs = adata.obs
    med = obs[obs["state"].isin(STATE_ORDER)].groupby("state", observed=True)[["total_counts", "n_genes", "pct_mt"]].median()
    med["n"] = obs[obs["state"].isin(STATE_ORDER)].groupby("state", observed=True).size()
    med.to_csv(OUT / "d7_a1_depth_by_state.csv"); rec.add_output(OUT / "d7_a1_depth_by_state.csv")
    p = (obs["state"] == "pAT2").to_numpy(); h = (obs["state"] == "hAT2").to_numpy()
    if p.sum() < FLOOR or h.sum() < FLOOR:
        return {"not_computable": "pAT2 or hAT2 below floor", "n_pAT2": int(p.sum()), "n_hAT2": int(h.sum())}
    rng = np.random.default_rng(0)
    ng = obs["n_genes"].to_numpy()
    edges = np.quantile(ng[p], np.linspace(0, 1, 21))
    h_idx = np.where(h)[0]
    matched = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        want = int(((ng[p] >= lo) & (ng[p] <= hi)).sum())
        pool = h_idx[(ng[h_idx] >= lo) & (ng[h_idx] <= hi)]
        if len(pool) == 0 or want == 0:
            continue
        matched.append(rng.choice(pool, size=min(want, len(pool)), replace=False))
    matched = np.concatenate(matched) if matched else np.array([], dtype=int)
    m_mask = np.zeros(len(obs), dtype=bool); m_mask[matched] = True
    det_p = detection(adata, STATE_SETS["AT2_identity"], p)
    det_h = detection(adata, STATE_SETS["AT2_identity"], h)
    det_hm = detection(adata, STATE_SETS["AT2_identity"], m_mask) if m_mask.sum() else pd.Series(dtype=float)
    tab = pd.DataFrame({"pAT2": det_p, "hAT2_all": det_h, "hAT2_depth_matched": det_hm})
    tab.loc["set_mean"] = tab.mean()
    tab.to_csv(OUT / "d7_a1_identity_detection.csv"); rec.add_output(OUT / "d7_a1_identity_detection.csv")
    ratio_all = float(det_p.mean() / det_h.mean()) if det_h.mean() > 0 else np.nan
    ratio_matched = float(det_p.mean() / det_hm.mean()) if len(det_hm) and det_hm.mean() > 0 else np.nan
    return {"median_genes_pAT2": float(np.median(ng[p])), "median_genes_hAT2": float(np.median(ng[h])),
            "n_matched_hAT2": int(m_mask.sum()), "identity_ratio_pAT2_over_hAT2": ratio_all,
            "identity_ratio_pAT2_over_matched_hAT2": ratio_matched,
            "not_depth_artefact": bool(ratio_matched <= 0.7) if np.isfinite(ratio_matched) else None}


def a2(adata, rec):
    rows = []
    for s in STATE_ORDER:
        m = (adata.obs["state"] == s).to_numpy()
        if m.sum() == 0:
            continue
        D = genes_matrix(adata, ["Sftpc"] + STATE_SETS["AT1_canonical"], m)
        double = float((D[:, 0] & (D[:, 1:].sum(axis=1) >= 2)).mean()) if D.shape[1] >= 4 else np.nan
        rows.append({"state": s, "n": int(m.sum()), "mean_scrublet_score": float(adata.obs.loc[m, "doublet_score"].mean()),
                     "score_gt_0.7_fraction": float((adata.obs.loc[m, "doublet_score"] > 0.7).mean()),
                     "sftpc_and_2plus_AT1_fraction": double})
    tab = pd.DataFrame(rows).set_index("state")
    tab.to_csv(OUT / "d7_a2_doublet_by_state.csv"); rec.add_output(OUT / "d7_a2_doublet_by_state.csv")
    if "DATP" not in tab.index or tab.loc["DATP", "n"] < FLOOR:
        return {"not_computable": "DATP below floor"}
    d = tab.loc["DATP"]
    doublet_like = bool(d["mean_scrublet_score"] > 0.6 or
                        (d["sftpc_and_2plus_AT1_fraction"] > tab.loc["hAT2", "sftpc_and_2plus_AT1_fraction"] if "hAT2" in tab.index else False) and
                        (d["sftpc_and_2plus_AT1_fraction"] > tab.loc["AT1", "sftpc_and_2plus_AT1_fraction"] if "AT1" in tab.index else False))
    return {"DATP_mean_scrublet": float(d["mean_scrublet_score"]), "DATP_double_marker_fraction": float(d["sftpc_and_2plus_AT1_fraction"]),
            "doublet_like": doublet_like}


def a4(adata, rec):
    m = ((adata.obs["library"] == "Day14_AT2_Tomato") & adata.obs["state"].isin(STATE_ORDER)).to_numpy()
    if m.sum() < FLOOR:
        return {"not_computable": "day-14 alveolar cells below floor"}
    D = genes_matrix(adata, ["Krt8", "Cldn4"], m)
    a, b = D[:, 0], D[:, 1]
    both = int((a & b).sum())
    expected = a.mean() * b.mean() * len(a)
    ratio = both / expected if expected > 0 else np.nan
    depth = adata.obs.loc[m, "total_counts"].to_numpy()
    deciles = np.digitize(depth, np.quantile(depth, np.linspace(0, 1, 11)[1:-1]))
    rng = np.random.default_rng(0)
    null = []
    for _ in range(200):
        bp = b.copy()
        for d in np.unique(deciles):
            idx = np.where(deciles == d)[0]
            bp[idx] = rng.permutation(bp[idx])
        null.append(int((a & bp).sum()))
    null = np.array(null)
    p = float((np.sum(null >= both) + 1) / (len(null) + 1))
    res = {"n_cells": int(m.sum()), "krt8_detected": float(a.mean()), "cldn4_detected": float(b.mean()),
           "double_positive": both, "expected_under_independence": float(expected), "ratio": float(ratio),
           "permutation_p_within_depth_deciles": p, "null_mean": float(null.mean())}
    if both < FLOOR:
        res["magnitude"] = "unreadable (fewer than 30 double-positive cells); existence only"
    res["reading"] = ("co-organised beyond depth" if (ratio > 1.25 and p < 0.05) else
                      "no evidence of a discrete co-expressing population" if 0.8 <= ratio <= 1.25 else
                      "outside both bands; reported without a verdict")
    pd.DataFrame([res]).to_csv(OUT / "d7_a4_codetection.csv", index=False); rec.add_output(OUT / "d7_a4_codetection.csv")
    return res


def main() -> None:
    rec = RunRecord(OUT / "d7_run_record.json", "D7 Choi-2020 attack claims", RULES)
    src = DERIVED / "tomato_annotated_d2b.h5ad"; rec.add_input(src)
    adata = ad.read_h5ad(src)
    results = {"A1": a1(adata, rec), "A2": a2(adata, rec), "A3": RULES["A3"], "A4": a4(adata, rec)}
    rec.set("results", results)
    lines = ["# Trial D7: attack claims", ""]
    for k, v in results.items():
        lines += [f"## {k}", "", str(v), ""]
    for f in ("d7_a1_depth_by_state.csv", "d7_a1_identity_detection.csv", "d7_a2_doublet_by_state.csv"):
        if (OUT / f).exists():
            lines += [f"### {f}", "", df_to_markdown(pd.read_csv(OUT / f, index_col=0).round(4)), ""]
    (OUT / "d7_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "d7_summary.md")
    rec.finish()
    print("[D7] done", {k: (v.get("not_depth_artefact") if k == "A1" else v.get("doublet_like") if k == "A2" else v.get("reading") if k == "A4" else "skipped") for k, v in results.items() if isinstance(v, dict)}, flush=True)


if __name__ == "__main__":
    main()
