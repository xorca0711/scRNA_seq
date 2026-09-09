#!/usr/bin/env python
"""Trial S3: annotate the GSE178360 clusters with the HLCA consensus marker sets.

Portfolio purpose. An annotation of the human distal-lung series that comes
from a cross-study consensus (61 identities, Sikkema et al. 2023,
Supplementary Table 6) rather than from this repository's own marker panel.
It serves the AT2-to-AT0 differentiation axis (Nabhan fit) and the
reference-based, donor-aware annotation angle (Wagner fit). No infection or
interferon framing is used; the series is healthy adult distal lung.

Frozen rules (written to the run record before any table is read):

* Scores: scanpy score_genes on log1p(CP10K), ctrl_size 50, n_bins 25,
  random_state 0. Three tiers per HLCA type as given in the sheet:
  tier 1 = compartment markers, tier 2 = intermediate-level markers,
  tier 3 = the type's own markers.
* A type is scorable only if at least one of its own markers is present
  in the object; otherwise it is excluded and listed. Fail closed.
* Per cell (primary, flat): compartment = argmax of the four tier-1 scores;
  type = argmax of the tier-3 scores among scorable types of that
  compartment.
* Per cluster: assigned type = mode of the per-cell type; confident if the
  mode holds >= 50% of the cluster's cells, otherwise "mixed".
* Comparison with the blind proposals uses the synonym table below, fixed
  in advance: agree / partial (same compartment) / disagree / reference
  lacks identity (neutrophil, platelet, erythroid: absent from the 61).
* AT0 check: HLCA-argmax AT0 cells versus the strict gate of the doublet
  audit (SFTPC > 0, SCGB3A2 > 0, EPCAM > 0, PTPRC = PECAM1 = COL1A1 = 0 on
  raw counts), per donor; concordant if the two counts are within a factor
  of two for every donor. The epithelial subcluster called "AT0 candidate
  analogue" (epithelial Leiden 4) is also scored if its object is present.
* Unit of report: the donor. No P values.

Rule revisions after the first run on 2026-09-09 (recorded, not hidden):

1. The first version required two own markers per type. The sheet gives a
   single own marker to B cells, Plasma cells, Pericytes, Tuft, EC aerocyte
   capillary and Smooth muscle FAM83D+, and gives Neuroendocrine only
   atlas-level markers, so that rule excluded eight types and left the B
   and neuroendocrine clusters unassignable. The minimum is now one, so
   the sheet decides which identities exist, not this script.
2. The compartment set is restricted to the four the sheet uses
   (Epithelial, Immune, Endothelial, Stroma); the first run treated the
   neuroendocrine atlas-level markers as a fifth compartment.
3. C8orf4 is read as TCIM (renamed) when the original symbol is absent.
4. A hierarchical assignment (compartment, then intermediate level, then
   type) is added as a post hoc sensitivity analysis because the flat
   argmax called AT0 across most secretory and AT2 clusters (AT0's own
   markers SFTPB, SCGB3A2 and SFTA2 are broadly expressed in distal
   epithelium while AT2's own markers are weak). The flat rule remains the
   primary, pre-registered result and both are reported.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trial_utils import REPO, RunRecord  # noqa: E402

OUT = HERE / "s3_hlca_marker_annotation"
OUT.mkdir(exist_ok=True)
HUMAN = REPO / "analysis" / "GSE178360"
OBJ = HUMAN / "processed" / "final_clustered.h5ad"
EPI_OBJ = HUMAN / "epithelial_subanalysis" / "epithelial_clustered.h5ad"
MARKERS = HERE / "hlca_supp_table6_markers.csv"
COMPARTMENTS = ["Epithelial", "Immune", "Endothelial", "Stroma"]
ALIASES = {"C8orf4": "TCIM"}
# Types whose sheet rows carry no compartment or intermediate marker; placed by Supplementary Table 5.
ORPHAN_PLACEMENT = {"Neuroendocrine": ("Epithelial", "Airway epithelium")}

RULES = {
    "score_genes": {"layer": "X = log1p(CP10K)", "ctrl_size": 50, "n_bins": 25, "random_state": 0},
    "min_own_markers_present": 1,
    "compartments": COMPARTMENTS,
    "gene_aliases": ALIASES,
    "orphan_placement": ORPHAN_PLACEMENT,
    "cell_assignment_primary": "compartment = argmax tier-1 score; type = argmax tier-3 score within compartment",
    "cell_assignment_sensitivity_post_hoc": "compartment, then argmax tier-2 group within compartment, then argmax tier-3 type within that group (types without an intermediate tier compete in every group of their compartment)",
    "cluster_assignment": "mode of per-cell type; confident if mode fraction >= 0.5 else mixed",
    "at0_strict_gate": "SFTPC>0 & SCGB3A2>0 & EPCAM>0 & PTPRC==0 & PECAM1==0 & COL1A1==0 on raw counts",
    "at0_concordance": "HLCA-argmax AT0 count within a factor of 2 of the strict-gate count, per donor",
    "unit": "donor",
    "reference_lacks": ["Neutrophil", "Platelet", "Erythroid", "Regulatory T cells", "gamma-delta T cells", "ILC", "NKT"],
    "rule_revisions_after_first_run": [
        "min own markers 2 -> 1 (the sheet assigns a single own marker to six types and atlas-level markers to Neuroendocrine)",
        "compartments restricted to the four in the sheet",
        "alias C8orf4 -> TCIM",
        "hierarchical assignment added as post hoc sensitivity; flat argmax stays primary",
    ],
    "first_run_outcome_kept_for_the_record": {
        "scorable_types": 53, "excluded": ["B cells", "EC aerocyte capillary", "Neuroendocrine", "Pericytes", "Plasma cells", "SMG serous (nasal)", "Smooth muscle FAM83D+", "Tuft"],
        "verdicts": {"agree": 16, "partial": 12, "disagree": 1, "reference lacks identity": 2},
        "at0_over_gate_per_donor": {"DD046Q": 1.644, "DD047Q": 1.826, "DD073R": 1.097},
        "note": "AT0 was the modal flat call for clusters 11, 13, 16, 18, 29, 30 (AT2, club, pan-epithelial, neuroendocrine and a donor-private ciliated-like cluster by the blind panel)",
    },
}

SYNONYMS = {
    "Mono_Mac": {"Classical monocytes", "Non-classical monocytes", "Monocyte-derived Mph", "Alveolar macrophages",
                 "Alveolar Mph CCL3+", "Alveolar Mph MT-positive", "Alveolar Mph proliferating", "Interstitial Mph perivascular"},
    "Alveolar_Mac": {"Alveolar macrophages", "Alveolar Mph CCL3+", "Alveolar Mph MT-positive", "Alveolar Mph proliferating"},
    "DC": {"DC1", "DC2", "Migratory DCs", "Plasmacytoid DCs"},
    "Endothelial": {"EC arterial", "EC general capillary", "EC aerocyte capillary", "EC venous pulmonary", "EC venous systemic"},
    "Endothelial_cap": {"EC general capillary", "EC aerocyte capillary"},
    "Lymphatic_EC": {"Lymphatic EC mature", "Lymphatic EC differentiating", "Lymphatic EC proliferating"},
    "CD8_T": {"CD8 T cells"},
    "CD4_T": {"CD4 T cells"},
    "T_cell": {"CD4 T cells", "CD8 T cells", "T cells proliferating"},
    "NK": {"NK cells"},
    "B_cell": {"B cells"},
    "Plasma": {"Plasma cells"},
    "Mast": {"Mast cells"},
    "Proliferating": {"T cells proliferating", "AT2 proliferating", "Alveolar Mph proliferating", "Lymphatic EC proliferating"},
    "Fibroblast": {"Adventitial fibroblasts", "Alveolar fibroblasts", "Peribronchial fibroblasts", "Subpleural fibroblasts", "Myofibroblasts"},
    "Fibro_activated": {"Myofibroblasts", "Alveolar fibroblasts"},
    "SMC_pericyte": {"Smooth muscle", "Smooth muscle FAM83D+", "SM activated stress response", "Pericytes"},
    "Mesothelium": {"Mesothelium"},
    "AT2": {"AT2", "AT2 proliferating"},
    "AT1": {"AT1"},
    "Club": {"Club (non-nasal)", "Club (nasal)", "pre-TB secretory"},
    "Goblet": {"Goblet (bronchial)", "Goblet (subsegmental)", "Goblet (nasal)"},
    "Transitional": {"AT0", "pre-TB secretory", "Suprabasal", "Basal resting", "Hillock-like"},
    "Basal": {"Basal resting", "Suprabasal", "Hillock-like"},
    "Ciliated": {"Multiciliated (non-nasal)", "Multiciliated (nasal)", "Deuterosomal"},
    "Neuroendocrine": {"Neuroendocrine"},
    "Ionocyte": {"Ionocyte"},
    "Tuft": {"Tuft"},
    "Pan_epithelial": set(),
    "Pan_immune": set(),
    "Neutrophil": None,
    "Platelet": None,
    "Erythroid": None,
}

COMPARTMENT_OF_BLIND = {
    "Mono_Mac": "Immune", "Alveolar_Mac": "Immune", "DC": "Immune", "CD8_T": "Immune", "CD4_T": "Immune", "T_cell": "Immune",
    "NK": "Immune", "B_cell": "Immune", "Plasma": "Immune", "Mast": "Immune", "Pan_immune": "Immune", "Neutrophil": "Immune",
    "Platelet": "Immune", "Erythroid": "Immune", "Proliferating": None,
    "Endothelial": "Endothelial", "Endothelial_cap": "Endothelial", "Lymphatic_EC": "Endothelial",
    "Fibroblast": "Stroma", "Fibro_activated": "Stroma", "SMC_pericyte": "Stroma", "Mesothelium": "Stroma",
    "AT2": "Epithelial", "AT1": "Epithelial", "Club": "Epithelial", "Goblet": "Epithelial", "Transitional": "Epithelial",
    "Basal": "Epithelial", "Ciliated": "Epithelial", "Neuroendocrine": "Epithelial", "Ionocyte": "Epithelial",
    "Tuft": "Epithelial", "Pan_epithelial": "Epithelial",
}


def verdict_for(blind: str, hlca_type: str, hlca_comp: str) -> str:
    syn = SYNONYMS.get(blind, set())
    if syn is None:
        return "reference lacks identity"
    if hlca_type in syn:
        return "agree"
    if COMPARTMENT_OF_BLIND.get(blind) == hlca_comp:
        return "partial (same compartment)"
    if COMPARTMENT_OF_BLIND.get(blind) is None:
        return "partial (blind label is compartment-free)"
    return "disagree"


def main() -> None:
    rec = RunRecord(OUT / "s3_run_record.json", "S3 HLCA consensus-marker annotation of GSE178360", RULES,
                    notes="rules frozen before reading the object; markers from Supplementary Table 6 (CC BY 4.0); see rule_revisions_after_first_run")
    import anndata as ad  # noqa: E402
    import scanpy as sc  # noqa: E402
    sc.settings.verbosity = 0

    rec.add_input(MARKERS)
    rec.add_input(OBJ)
    markers = pd.read_csv(MARKERS)
    adata = ad.read_h5ad(OBJ)
    present = set(adata.var_names)
    rec.set("object_shape", list(adata.shape))

    def resolve(genes: list[str]) -> tuple[list[str], dict]:
        used, aliased = [], {}
        for g in genes:
            if g in present:
                used.append(g)
            elif ALIASES.get(g) in present:
                used.append(ALIASES[g]); aliased[g] = ALIASES[g]
        return used, aliased

    types = sorted(markers["cell_type"].unique())

    def first_or(series: pd.Series, default: str) -> str:
        return str(series.iloc[0]) if len(series) else default

    comp_genes = {c: sorted(set(markers.loc[(markers.tier == 1) & (markers.marker_for == c), "gene"])) for c in COMPARTMENTS}
    groups = sorted(markers.loc[markers.tier == 2, "marker_for"].unique())
    group_genes = {g: sorted(set(markers.loc[(markers.tier == 2) & (markers.marker_for == g), "gene"])) for g in groups}
    type_genes_raw = {t: sorted(set(markers.loc[(markers.tier == 3) & (markers.cell_type == t), "gene"])) for t in types}
    type_comp = {t: first_or(markers.loc[(markers.tier == 1) & (markers.cell_type == t), "marker_for"], "NA") for t in types}
    type_group = {t: first_or(markers.loc[(markers.tier == 2) & (markers.cell_type == t), "marker_for"], "NA") for t in types}
    for t, (c, g) in ORPHAN_PLACEMENT.items():
        if t in type_comp and type_comp[t] == "NA":
            type_comp[t] = c
        if t in type_group and type_group[t] == "NA":
            type_group[t] = g
    for t in types:
        if type_group[t] == "NA":
            type_group[t] = type_comp[t]  # no intermediate tier: the group is the compartment itself

    type_genes, aliases_used, absent = {}, {}, {}
    for t, gl in type_genes_raw.items():
        used, al = resolve(gl)
        type_genes[t] = used
        aliases_used.update(al)
        absent[t] = [g for g in gl if g not in present and ALIASES.get(g) not in present]
    scorable = [t for t in types if len(type_genes[t]) >= RULES["min_own_markers_present"] and type_comp[t] in COMPARTMENTS]
    excluded = [t for t in types if t not in scorable]
    rec.set("types_total", len(types)); rec.set("types_scorable", len(scorable))
    rec.set("types_excluded", {t: {"own_markers_in_sheet": type_genes_raw[t], "absent": absent[t]} for t in excluded})
    rec.set("absent_markers_by_type", {t: a for t, a in absent.items() if a})
    rec.set("aliases_used", aliases_used)
    rec.set("types_without_intermediate_tier", [t for t in types if type_group[t] == type_comp[t]])
    print("scorable types", len(scorable), "excluded", excluded, "aliases", aliases_used)

    # ---------------------------------------------------------------- scores
    def score(name: str, genes: list[str]) -> str:
        genes = [g for g in genes if g in present]
        col = f"score__{name}"
        sc.tl.score_genes(adata, genes, score_name=col, ctrl_size=RULES["score_genes"]["ctrl_size"],
                          n_bins=RULES["score_genes"]["n_bins"], random_state=RULES["score_genes"]["random_state"])
        return col

    comp_cols = {c: score("comp__" + c, g) for c, g in comp_genes.items()}
    group_cols = {g: score("group__" + g, gl) for g, gl in group_genes.items()}
    type_cols = {t: score("type__" + t, type_genes[t]) for t in scorable}

    comp_mat = adata.obs[[comp_cols[c] for c in COMPARTMENTS]].to_numpy()
    cell_comp = np.array(COMPARTMENTS)[comp_mat.argmax(axis=1)]
    type_arr = np.array(scorable)
    type_mat = adata.obs[[type_cols[t] for t in scorable]].to_numpy()
    comp_of_type = np.array([type_comp[t] for t in scorable])
    group_of_type = np.array([type_group[t] for t in scorable])

    # primary flat assignment
    cell_type = np.full(adata.n_obs, "unassigned", dtype=object)
    for c in COMPARTMENTS:
        cells = np.where(cell_comp == c)[0]
        cols = np.where(comp_of_type == c)[0]
        if len(cells) and len(cols):
            cell_type[cells] = type_arr[cols][type_mat[np.ix_(cells, cols)].argmax(axis=1)]

    # post hoc hierarchical assignment
    cell_group = np.full(adata.n_obs, "NA", dtype=object)
    cell_type_h = np.full(adata.n_obs, "unassigned", dtype=object)
    for c in COMPARTMENTS:
        cells = np.where(cell_comp == c)[0]
        cgroups = sorted({g for g, cc in zip(group_of_type, comp_of_type) if cc == c and g != c and g in group_cols})
        if not len(cells):
            continue
        if cgroups:
            gm = adata.obs[[group_cols[g] for g in cgroups]].to_numpy()[cells]
            chosen = np.array(cgroups)[gm.argmax(axis=1)]
        else:
            chosen = np.full(len(cells), c, dtype=object)
        cell_group[cells] = chosen
        for g in set(chosen):
            sub_cells = cells[chosen == g]
            cols = np.where((comp_of_type == c) & ((group_of_type == g) | (group_of_type == c)))[0]
            if len(cols):
                cell_type_h[sub_cells] = type_arr[cols][type_mat[np.ix_(sub_cells, cols)].argmax(axis=1)]
    adata.obs["hlca_compartment"] = cell_comp
    adata.obs["hlca_type"] = cell_type
    adata.obs["hlca_group_hier"] = cell_group
    adata.obs["hlca_type_hier"] = cell_type_h

    # ------------------------------------------------------- per cluster
    clusters = adata.obs["leiden_cluster"].astype(str)
    blind = adata.obs["proposed_cell_type"].astype(str).str.replace(" (candidate)", "", regex=False)
    rows = []
    for cl in sorted(clusters.unique(), key=int):
        m = (clusters == cl).to_numpy(); n = int(m.sum())
        comp_cnt = Counter(cell_comp[m]); mode_comp, mode_comp_n = comp_cnt.most_common(1)[0]
        cnt = Counter(cell_type[m]); mode, mode_n = cnt.most_common(1)[0]
        cnt_h = Counter(cell_type_h[m]); mode_h, mode_h_n = cnt_h.most_common(1)[0]
        b = blind[m].iloc[0]
        means = {t: float(type_mat[m][:, i].mean()) for i, t in enumerate(scorable) if type_comp[t] == mode_comp}
        top3 = sorted(means.items(), key=lambda kv: -kv[1])[:3]
        per_donor = adata.obs.loc[m].groupby("sample_id", observed=True).size().to_dict()
        rows.append({
            "cluster": cl, "n_cells": n, "blind_proposal": b,
            "hlca_compartment": mode_comp, "compartment_fraction": round(mode_comp_n / n, 3),
            "hlca_type_flat": mode, "flat_mode_fraction": round(mode_n / n, 3),
            "flat_confidence": "confident" if mode_n / n >= 0.5 else "mixed",
            "flat_runner_up": cnt.most_common(2)[1][0] if len(cnt) > 1 else "",
            "flat_runner_up_fraction": round(cnt.most_common(2)[1][1] / n, 3) if len(cnt) > 1 else 0.0,
            "verdict_flat_vs_blind": verdict_for(b, mode, mode_comp),
            "hlca_group_hier": Counter(cell_group[m]).most_common(1)[0][0],
            "hlca_type_hier": mode_h, "hier_mode_fraction": round(mode_h_n / n, 3),
            "hier_confidence": "confident" if mode_h_n / n >= 0.5 else "mixed",
            "verdict_hier_vs_blind": verdict_for(b, mode_h, mode_comp),
            "top3_mean_own_marker_score": "; ".join(f"{t}={v:.3f}" for t, v in top3),
            "cells_per_donor": "; ".join(f"{k}={v}" for k, v in per_donor.items()),
        })
    assign = pd.DataFrame(rows)
    assign.to_csv(OUT / "s3_cluster_assignment.csv", index=False); rec.add_output(OUT / "s3_cluster_assignment.csv")
    mean_scores = pd.DataFrame({t: adata.obs.groupby(clusters, observed=True)[type_cols[t]].mean() for t in scorable})
    mean_scores.index.name = "cluster"
    mean_scores.round(4).to_csv(OUT / "s3_cluster_by_type_mean_score.csv"); rec.add_output(OUT / "s3_cluster_by_type_mean_score.csv")
    adata.obs[["sample_id", "leiden_cluster", "proposed_cell_type", "hlca_compartment", "hlca_type", "hlca_group_hier", "hlca_type_hier"]] \
        .to_csv(OUT / "s3_cell_assignment.csv.gz", index=True, compression="gzip")
    vc_flat = Counter(assign["verdict_flat_vs_blind"]); vc_hier = Counter(assign["verdict_hier_vs_blind"])
    rec.set("verdict_counts_flat", dict(vc_flat)); rec.set("verdict_counts_hier", dict(vc_hier))
    rec.set("clusters_mixed_flat", assign.loc[assign.flat_confidence == "mixed", "cluster"].tolist())
    rec.set("clusters_mixed_hier", assign.loc[assign.hier_confidence == "mixed", "cluster"].tolist())

    # ------------------------------------------------------------- AT0 check
    counts = adata.layers["counts"]
    gi = {g: list(adata.var_names).index(g) for g in ("SFTPC", "SCGB3A2", "EPCAM", "PTPRC", "PECAM1", "COL1A1")}
    col = lambda g: np.asarray(counts[:, gi[g]].todense()).ravel()  # noqa: E731
    gate = (col("SFTPC") > 0) & (col("SCGB3A2") > 0) & (col("EPCAM") > 0) & (col("PTPRC") == 0) & (col("PECAM1") == 0) & (col("COL1A1") == 0)
    adata.obs["at0_strict_gate"] = gate
    at0 = []
    for donor, sub in adata.obs.groupby("sample_id", observed=True):
        idx = adata.obs.index.get_indexer(sub.index)
        g = int(gate[idx].sum())
        a = int((cell_type[idx] == "AT0").sum()); ah = int((cell_type_h[idx] == "AT0").sum())
        p = int((cell_type[idx] == "pre-TB secretory").sum()); ph = int((cell_type_h[idx] == "pre-TB secretory").sum())
        both = int((gate[idx] & (cell_type[idx] == "AT0")).sum()); both_h = int((gate[idx] & (cell_type_h[idx] == "AT0")).sum())
        r = a / g if g else None; rh = ah / g if g else None
        at0.append({"donor": donor, "n_cells": len(idx), "strict_gate": g,
                    "flat_AT0": a, "flat_preTB": p, "flat_gate_and_AT0": both, "flat_AT0_over_gate": round(r, 3) if r is not None else None,
                    "flat_within_factor_2": (0.5 <= r <= 2.0) if r is not None else None,
                    "hier_AT0": ah, "hier_preTB": ph, "hier_gate_and_AT0": both_h, "hier_AT0_over_gate": round(rh, 3) if rh is not None else None,
                    "hier_within_factor_2": (0.5 <= rh <= 2.0) if rh is not None else None})
    at0_df = pd.DataFrame(at0); at0_df.to_csv(OUT / "s3_at0_check_per_donor.csv", index=False); rec.add_output(OUT / "s3_at0_check_per_donor.csv")
    pd.crosstab(clusters[cell_type == "AT0"], adata.obs.loc[cell_type == "AT0", "sample_id"]).to_csv(OUT / "s3_at0_flat_cells_by_cluster_and_donor.csv")
    pd.crosstab(clusters[cell_type_h == "AT0"], adata.obs.loc[cell_type_h == "AT0", "sample_id"]).to_csv(OUT / "s3_at0_hier_cells_by_cluster_and_donor.csv")
    rec.add_output(OUT / "s3_at0_flat_cells_by_cluster_and_donor.csv"); rec.add_output(OUT / "s3_at0_hier_cells_by_cluster_and_donor.csv")

    epi_result = None
    if EPI_OBJ.exists():
        rec.add_input(EPI_OBJ)
        epi = ad.read_h5ad(EPI_OBJ, backed="r")
        key = "leiden_cluster" if "leiden_cluster" in epi.obs else [c for c in epi.obs.columns if c.startswith("leiden")][0]
        common = epi.obs.index.intersection(adata.obs.index)
        epi_lab = epi.obs.loc[common, key].astype(str)
        tab = pd.crosstab(epi_lab, adata.obs.loc[common, "hlca_type"]); tab.to_csv(OUT / "s3_epithelial_subcluster_by_hlca_type_flat.csv")
        tab_h = pd.crosstab(epi_lab, adata.obs.loc[common, "hlca_type_hier"]); tab_h.to_csv(OUT / "s3_epithelial_subcluster_by_hlca_type_hier.csv")
        rec.add_output(OUT / "s3_epithelial_subcluster_by_hlca_type_flat.csv"); rec.add_output(OUT / "s3_epithelial_subcluster_by_hlca_type_hier.csv")
        epi_result = {}
        for name, tb in (("flat", tab), ("hier", tab_h)):
            if "4" in tb.index:
                row = tb.loc["4"]
                epi_result[name] = {"n": int(row.sum()), "AT0": int(row.get("AT0", 0)), "pre-TB secretory": int(row.get("pre-TB secretory", 0)),
                                    "AT2": int(row.get("AT2", 0)), "top4": row.sort_values(ascending=False).head(4).to_dict()}
        epi.file.close()
    rec.set("epithelial_subcluster_4_vs_hlca", epi_result); rec.set("at0_per_donor", at0)

    # --------------------------------------------------------------- figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    top_types = [t for t in pd.unique(pd.concat([assign["hlca_type_flat"], assign["hlca_type_hier"]])) if t in mean_scores.columns]
    fig, axm = plt.subplots(figsize=(max(8, 0.35 * len(top_types)), 9))
    mat = mean_scores[top_types]
    im = axm.imshow(mat.to_numpy(), aspect="auto", cmap="viridis")
    axm.set_yticks(range(len(mat.index))); axm.set_yticklabels(mat.index)
    axm.set_xticks(range(len(top_types))); axm.set_xticklabels(top_types, rotation=90, fontsize=7)
    axm.set_title("GSE178360 clusters x HLCA own-marker score (mean per cluster; assigned types only)")
    fig.colorbar(im, ax=axm, label="score_genes"); fig.tight_layout(); fig.savefig(OUT / "s3_cluster_by_type_heatmap.png", dpi=130); plt.close(fig)
    rec.add_output(OUT / "s3_cluster_by_type_heatmap.png")
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))
    sc.pl.umap(adata, color="hlca_type", ax=axes[0], show=False, legend_loc="on data", legend_fontsize=5, size=4, title="HLCA type, flat argmax (primary)")
    sc.pl.umap(adata, color="hlca_type_hier", ax=axes[1], show=False, legend_loc="on data", legend_fontsize=5, size=4, title="HLCA type, hierarchical (post hoc sensitivity)")
    sc.pl.umap(adata, color="proposed_cell_type", ax=axes[2], show=False, legend_loc="on data", legend_fontsize=5, size=4, title="Blind proposal (existing)")
    fig.tight_layout(); fig.savefig(OUT / "s3_umap_hlca_vs_blind.png", dpi=110); plt.close(fig)
    rec.add_output(OUT / "s3_umap_hlca_vs_blind.png")

    # --------------------------------------------------------------- summary
    lines = ["# Trial S3 output: HLCA consensus-marker annotation of GSE178360", "",
             f"Generated by `s3_hlca_marker_annotation.py`. Types scorable: {len(scorable)} of {len(types)}; excluded: "
             f"{', '.join(excluded) or 'none'}; aliases used: {aliases_used or 'none'}. The flat argmax is the pre-registered primary; "
             "the hierarchical assignment is a post hoc sensitivity (see the script docstring).", "",
             "## Cluster assignment", "",
             "| Cluster | Cells | Blind proposal | Compartment | Flat type | Flat mode | Verdict flat | Hier group | Hier type | Hier mode | Verdict hier |",
             "|--:|--:|---|---|---|--:|---|---|---|--:|---|"]
    for r in rows:
        lines.append(f"| {r['cluster']} | {r['n_cells']} | {r['blind_proposal']} | {r['hlca_compartment']} ({r['compartment_fraction']}) | "
                     f"{r['hlca_type_flat']} | {r['flat_mode_fraction']} | {r['verdict_flat_vs_blind']} | {r['hlca_group_hier']} | "
                     f"{r['hlca_type_hier']} | {r['hier_mode_fraction']} | {r['verdict_hier_vs_blind']} |")
    lines += ["", "Verdicts, flat: " + ", ".join(f"{k}: {v}" for k, v in sorted(vc_flat.items())) + ".",
              "Verdicts, hierarchical: " + ", ".join(f"{k}: {v}" for k, v in sorted(vc_hier.items())) + ".", "",
              "## AT0 check per donor", "",
              "| Donor | Cells | Strict gate | Flat AT0 | Flat pre-TB | Flat gate and AT0 | Flat AT0/gate | Flat within 2x | Hier AT0 | Hier pre-TB | Hier gate and AT0 | Hier AT0/gate | Hier within 2x |",
              "|---|--:|--:|--:|--:|--:|--:|---|--:|--:|--:|--:|---|"]
    for r in at0:
        lines.append(f"| {r['donor']} | {r['n_cells']} | {r['strict_gate']} | {r['flat_AT0']} | {r['flat_preTB']} | {r['flat_gate_and_AT0']} | {r['flat_AT0_over_gate']} | {r['flat_within_factor_2']} | "
                     f"{r['hier_AT0']} | {r['hier_preTB']} | {r['hier_gate_and_AT0']} | {r['hier_AT0_over_gate']} | {r['hier_within_factor_2']} |")
    if epi_result:
        for name, e in epi_result.items():
            lines.append(f"\nEpithelial subcluster 4 (existing AT0 candidate analogue, n={e['n']}), {name}: AT0 {e['AT0']}, pre-TB secretory {e['pre-TB secretory']}, AT2 {e['AT2']}; top {e['top4']}.")
    (OUT / "s3_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8"); rec.add_output(OUT / "s3_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
