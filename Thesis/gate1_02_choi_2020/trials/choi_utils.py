"""Shared helpers for the Choi-2020 trials (roadmap paper 2, Gate 1).

Everything generic comes from the repository's one shared helper module, the
same one the Cardoso trials use: the run record, the MatrixMarket triplet
reader, the GEO SOFT family parser, the non-Ensembl feature rule and the QC
metric helper. Nothing is redefined here.

This module adds only what is specific to this deposit: where its libraries
live, the fact that six of its matrices are unfiltered 10x barcode whitelists
rather than called cells, the paper's own cell filter, the paper's marker
vocabulary for its five states (read from the text, see
choi_2020_extracts.json), and the frozen cluster-annotation rule the trials
share. The rule is written once here so that every trial grades the same way
and no trial can move it after seeing a result.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "raw_data"
HERE = Path(__file__).resolve().parent
EXTRACTS = HERE.parent / "choi_2020_extracts.json"
PALETTE = REPO / "analysis" / "config" / "palette.json"

_SHARED_TRIALS = REPO / "Thesis" / "gate1_04_sikkema_2023_hlca" / "trials"
if str(_SHARED_TRIALS) not in sys.path:
    sys.path.insert(0, str(_SHARED_TRIALS))
from trial_utils import (ENSEMBL_ID_RE, RunRecord, df_to_markdown,  # noqa: E402,F401
                         file_facts, package_versions, parse_soft, qc_metrics,
                         read_mtx_triplet, shannon, utc_now)

# The 10x v2 barcode whitelist. A deposited matrix with this many columns is
# the raw droplet matrix, not called cells, so cell calling is this
# repository's job and the cell count will not match the paper's unless the
# same caller and version are used.
TENX_V2_WHITELIST = 737_280

# Regenerable intermediates live under raw_data/, which is gitignored (the
# convention the Cardoso trials set). Nothing under raw_data/ is tracked.
DERIVED = RAW / "GSE145031" / "choi_trials"
DERIVED_ORGANOID = RAW / "GSE144468" / "choi_trials"

LINEAGE_LIBRARIES = [
    {"gsm": "GSM4304609", "library": "PBS_AT2_Tomato", "timepoint": "PBS", "sort": "Tomato"},
    {"gsm": "GSM4304610", "library": "PBS_AT2_nonTomato", "timepoint": "PBS", "sort": "nonTomato"},
    {"gsm": "GSM4304611", "library": "Day14_AT2_Tomato", "timepoint": "day 14", "sort": "Tomato"},
    {"gsm": "GSM4304612", "library": "Day14_AT2_nonTomato", "timepoint": "day 14", "sort": "nonTomato"},
    {"gsm": "GSM4304613", "library": "Day28_AT2_Tomato", "timepoint": "day 28", "sort": "Tomato"},
    {"gsm": "GSM4304614", "library": "Day28_AT2_nonTomato", "timepoint": "day 28", "sort": "nonTomato"},
]
ORGANOID_LIBRARIES = [
    {"gsm": "GSM4288824", "library": "Control_Organoids", "treatment": "control"},
    {"gsm": "GSM4288825", "library": "Il1b_Organoids", "treatment": "IL-1beta"},
]
ACCESSIONS = {
    "GSE145031": {"role": "scRNA-seq of AT2 lineage-traced epithelium",
                  "libraries": LINEAGE_LIBRARIES},
    "GSE144468": {"role": "scRNA-seq of AT2 organoids",
                  "libraries": ORGANOID_LIBRARIES},
}
TOMATO_LIBRARIES = [e for e in LINEAGE_LIBRARIES if e["sort"] == "Tomato"]
NONTOMATO_LIBRARIES = [e for e in LINEAGE_LIBRARIES if e["sort"] == "nonTomato"]

# The paper's own cell filter (STAR Methods): more than 500 and fewer than
# 7,000 detected genes, more than 2,000 UMI. Applied here to every barcode of
# the raw whitelist in place of Cell Ranger 2.0.2's caller, which is not
# reproduced. The sensitivity rule is the plain floor trial D0 used to size
# the libraries.
PAPER_FILTER = {"min_genes_exclusive": 500, "max_genes_exclusive": 7000, "min_counts_exclusive": 2000}
FLOOR_FILTER = {"min_counts": 500, "min_genes": 200}

# The paper's Scrublet settings, recorded for the divergence table. This
# repository runs scanpy's Scrublet per capture with the 10x expected rate,
# which is its standing rule; the paper's 0.7 score cut and 0.6 cluster mean
# are reported alongside, never substituted.
PAPER_SCRUBLET = {"sim_doublet_ratio": 2, "n_neighbors": 30, "expected_doublet_rate": 0.1,
                  "score_cut": 0.7, "cluster_mean_cut": 0.6}


def library_paths(accession: str, library: str, gsm: str) -> dict[str, Path]:
    """Where this deposit's triplet files sit after the series tar is unpacked."""
    root = RAW / accession / (accession + "_RAW")
    stem = gsm + "_" + library
    return {"matrix": root / (stem + ".mtx.gz"),
            "features": root / (stem + "_gene.tsv.gz"),
            "barcodes": root / (stem + "_barcodes.tsv.gz")}


def extracts() -> dict:
    return json.loads(EXTRACTS.read_text(encoding="utf-8"))


# ---------------------------------------------------------------------------
# The paper's state vocabulary, as marker sets. Every gene is one the paper
# names in its text or figure legends; the source figure is recorded in the
# extract. Sets are detection-fraction sets, not weights.
# ---------------------------------------------------------------------------
STATE_SETS = {
    "hAT2_canonical": ["Sftpc", "Sftpa1", "Lyz2"],
    "AT2_identity": ["Etv5", "Abca3", "Cebpa"],                 # DOWN in primed AT2
    "AT2_lipid": ["Acly", "Hmgcr", "Hmgcs1"],                   # DOWN in primed AT2
    "pAT2_inflammatory": ["Ptges", "Orm1", "Tmem173", "Ifitm2", "Ifitm3"],
    "cAT2": ["Cdk1", "Mki67", "Cenpa"],
    "DATP": ["Cldn4", "Krt8", "Ndrg1", "Sprr1a", "AW112010"],
    "AT1_canonical": ["Pdpn", "Hopx", "Cav1"],                  # LOW in DATP
    "AT1_early": ["Lmo7", "Pdpn", "Hopx"],
    "AT1_late": ["Aqp5", "Vegfa", "Cav1", "Spock2"],
    "DATP_figure7": ["Cldn4", "AW112010", "Lhfp"],
    "p53": ["Trp53", "Mdm2", "Ccnd1", "Gdf15"],
    "arrest": ["Cdkn1a", "Cdkn2a"],
    "hypoxia": ["Hif1a", "Ndrg1"],
    "hypoxia_without_Ndrg1": ["Hif1a"],
    "ifng_response": ["Ifngr1", "Ly6a", "Irf7", "Cxcl16"],
    "glycolysis": ["Pgk1", "Pkm", "Slc16a3"],
    "responder": ["Il1r1"],
}
# The printed text lists a sixth inflammatory gene as "Lcn1". Whether the
# authors meant Lcn1 or Lcn2 cannot be settled from the PDF, so the annotation
# set above leaves the lipocalin out, and trial D2 reports the detection of
# both symbols in every state so the reader can see which one the data carry.
PAPER_TEXT_AMBIGUITY = {"pAT2_inflammatory": {"printed": "Lcn1", "candidates": ["Lcn1", "Lcn2"],
                                              "handling": "excluded from the rule; both reported"}}

# Contaminant panels the paper removed (ciliated, mesenchyme, immune), with the
# genes the paper names, plus the repository's pan-compartment markers.
CONTAMINANT_SETS = {
    "ciliated": ["Foxj1", "Wnt7b", "Cd24a"],
    "mesenchyme": ["Vcam1", "Acta2", "Des", "Pdgfra", "Col1a1", "Col1a2"],
    "immune": ["Ptprc", "Tyrobp", "Il2rg", "Lck"],
    "endothelium": ["Pecam1", "Cdh5", "Cldn5"],
    "club_or_airway": ["Scgb1a1", "Scgb3a2", "Cyp2f2"],
}

STATE_ORDER = ["hAT2", "cAT2", "pAT2", "DATP", "AT1"]


# ---------------------------------------------------------------------------
# Frozen annotation rule (cluster level). Written before any trial ran.
#
#   1. contaminant: a cluster whose mean detection of any CONTAMINANT_SET is
#      at least 0.5 and whose Sftpc detection is below 0.5.
#   2. cAT2: mean detection of the cAT2 set at least 0.35.
#   3. AT1: mean detection of AT1_canonical at least 0.5 and Sftpc detection
#      below 0.5.
#   4. DATP: mean detection of the DATP set at least 0.4 and AT1_canonical
#      mean detection below half of the AT1 clusters' mean (the negative
#      condition); if no AT1 cluster exists, below 0.5 absolute.
#   5. reference hAT2: among clusters with Sftpc detection at least 0.7 that
#      are none of the above, the one with the highest AT2_identity mean
#      detection.
#   6. pAT2: a remaining Sftpc-high cluster whose AT2_identity mean detection
#      is at most 0.7 times the reference AND whose pAT2_inflammatory mean
#      detection is at least 1.5 times the reference.
#   7. hAT2: every other Sftpc-high cluster.
#   8. unassigned: anything left.
#
# What makes the rule's own answer unreadable, stated in advance: if a state
# is absent at every pre-declared resolution (0.3, 0.5, 1.0) but at least
# 2 percent of alveolar cells carry three or more of its markers, the state is
# reported as "present but not resolved by clustering", not as absent.
# ---------------------------------------------------------------------------
ANNOTATION_THRESHOLDS = {
    "contaminant_min": 0.5, "sftpc_low": 0.5, "cAT2_min": 0.35, "AT1_min": 0.5,
    "DATP_min": 0.4, "DATP_negative_ratio": 0.5, "sftpc_high": 0.7,
    "pAT2_identity_ratio": 0.7, "pAT2_inflammatory_ratio": 1.5,
    "dispersed_present_fraction": 0.02, "dispersed_markers": 3,
}
RESOLUTIONS = [0.3, 0.5, 1.0]


def detection(adata, genes, mask=None):
    """Fraction of cells detecting each gene, from the counts layer."""
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    present = [g for g in genes if g in adata.var_names]
    if mask is not None:
        X = X[mask]
    if not present:
        return pd.Series(dtype=float)
    sub = X[:, [adata.var_names.get_loc(g) for g in present]]
    return pd.Series(np.asarray((sub > 0).mean(axis=0)).ravel(), index=present)


def cluster_detection_table(adata, cluster_key):
    """Per-cluster mean detection of every set in STATE_SETS and CONTAMINANT_SETS."""
    rows = []
    for cl in sorted(adata.obs[cluster_key].astype(str).unique(), key=lambda s: int(s) if s.isdigit() else s):
        mask = (adata.obs[cluster_key].astype(str) == cl).to_numpy()
        row = {"cluster": cl, "n_cells": int(mask.sum())}
        sftpc = detection(adata, ["Sftpc"], mask)
        row["Sftpc"] = float(sftpc.get("Sftpc", np.nan))
        for name, genes in {**STATE_SETS, **CONTAMINANT_SETS}.items():
            d = detection(adata, genes, mask)
            row[name] = float(d.mean()) if len(d) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def annotate_clusters(table: pd.DataFrame) -> pd.DataFrame:
    """Apply the frozen rule to a cluster detection table. Returns the table
    with 'state' and 'rule' columns; nothing is fitted."""
    T = ANNOTATION_THRESHOLDS
    table = table.copy()
    state, rule = {}, {}
    for _, r in table.iterrows():
        cl = r["cluster"]
        cont = {k: r[k] for k in CONTAMINANT_SETS if r[k] >= T["contaminant_min"]}
        if cont and r["Sftpc"] < T["sftpc_low"]:
            state[cl], rule[cl] = "contaminant_" + max(cont, key=cont.get), "1"
        elif r["cAT2"] >= T["cAT2_min"]:
            state[cl], rule[cl] = "cAT2", "2"
        elif r["AT1_canonical"] >= T["AT1_min"] and r["Sftpc"] < T["sftpc_low"]:
            state[cl], rule[cl] = "AT1", "3"
    at1 = [cl for cl, s in state.items() if s == "AT1"]
    at1_ref = float(table.set_index("cluster").loc[at1, "AT1_canonical"].mean()) if at1 else np.nan
    for _, r in table.iterrows():
        cl = r["cluster"]
        if cl in state:
            continue
        neg = (r["AT1_canonical"] < T["DATP_negative_ratio"] * at1_ref) if at1 else (r["AT1_canonical"] < T["AT1_min"])
        if r["DATP"] >= T["DATP_min"] and neg:
            state[cl], rule[cl] = "DATP", "4"
    sftpc_high = [r for _, r in table.iterrows() if r["cluster"] not in state and r["Sftpc"] >= T["sftpc_high"]]
    if sftpc_high:
        ref = max(sftpc_high, key=lambda r: r["AT2_identity"])
        ref_id, ref_inf = ref["AT2_identity"], ref["pAT2_inflammatory"]
        state[ref["cluster"]], rule[ref["cluster"]] = "hAT2", "5 (reference)"
        for r in sftpc_high:
            cl = r["cluster"]
            if cl in state:
                continue
            if r["AT2_identity"] <= T["pAT2_identity_ratio"] * ref_id and \
               r["pAT2_inflammatory"] >= T["pAT2_inflammatory_ratio"] * max(ref_inf, 1e-9):
                state[cl], rule[cl] = "pAT2", "6"
            else:
                state[cl], rule[cl] = "hAT2", "7"
    for _, r in table.iterrows():
        state.setdefault(r["cluster"], "unassigned")
        rule.setdefault(r["cluster"], "8")
    table["state"] = table["cluster"].map(state)
    table["rule"] = table["cluster"].map(rule)
    table["AT1_reference_detection"] = at1_ref
    return table


def annotate_clusters_b(table: pd.DataFrame) -> pd.DataFrame:
    """The corrected pass (trial D2b), thresholds unchanged.

    D2's first pass showed that in a lineage-sorted library every cluster
    detects Sftpc (ambient surfactant transcript plus the sort), so the
    rule's "Sftpc below 0.5" clauses for AT1 and for contaminants never fire
    and an AT1 cluster is called hAT2 while ciliated and immune clusters stay
    in. This pass keeps every threshold and adds one condition on the
    condition: the Sftpc clauses are applied only where Sftpc separates
    clusters, meaning at least one cluster falls below the low bound. Where
    no cluster does, Sftpc is uninformative and the clauses are dropped. The
    first pass's outcome stays in the record beside this one.
    """
    T = ANNOTATION_THRESHOLDS
    informative = bool((table["Sftpc"] < T["sftpc_low"]).any())
    if informative:
        out = annotate_clusters(table)
        out["sftpc_informative"] = True
        return out
    # Sftpc is uninformative here: the clauses are dropped by evaluating the
    # same rules on the same table with the Sftpc conditions removed.
    t2 = table.copy()
    state, rule = {}, {}
    for _, r in t2.iterrows():
        cl = r["cluster"]
        cont = {k: r[k] for k in CONTAMINANT_SETS if r[k] >= T["contaminant_min"]}
        if cont:
            state[cl], rule[cl] = "contaminant_" + max(cont, key=cont.get), "1b"
        elif r["cAT2"] >= T["cAT2_min"]:
            state[cl], rule[cl] = "cAT2", "2"
        elif r["AT1_canonical"] >= T["AT1_min"]:
            state[cl], rule[cl] = "AT1", "3b"
    at1 = [cl for cl, s in state.items() if s == "AT1"]
    at1_ref = float(t2.set_index("cluster").loc[at1, "AT1_canonical"].mean()) if at1 else np.nan
    for _, r in t2.iterrows():
        cl = r["cluster"]
        if cl in state:
            continue
        neg = (r["AT1_canonical"] < T["DATP_negative_ratio"] * at1_ref) if at1 else (r["AT1_canonical"] < T["AT1_min"])
        if r["DATP"] >= T["DATP_min"] and neg:
            state[cl], rule[cl] = "DATP", "4"
    rest = [r for _, r in t2.iterrows() if r["cluster"] not in state]
    if rest:
        ref = max(rest, key=lambda r: r["AT2_identity"])
        ref_id, ref_inf = ref["AT2_identity"], ref["pAT2_inflammatory"]
        state[ref["cluster"]], rule[ref["cluster"]] = "hAT2", "5 (reference)"
        for r in rest:
            cl = r["cluster"]
            if cl in state:
                continue
            if r["AT2_identity"] <= T["pAT2_identity_ratio"] * ref_id and \
               r["pAT2_inflammatory"] >= T["pAT2_inflammatory_ratio"] * max(ref_inf, 1e-9):
                state[cl], rule[cl] = "pAT2", "6"
            else:
                state[cl], rule[cl] = "hAT2", "7"
    out = table.copy()
    out["state"] = out["cluster"].map(state).fillna("unassigned")
    out["rule"] = out["cluster"].map(rule).fillna("8")
    out["AT1_reference_detection"] = at1_ref
    out["sftpc_informative"] = False
    return out


def dispersed_presence(adata, genes, min_markers=3, mask=None):
    """Fraction of (masked) cells carrying at least min_markers of the set."""
    X = adata.layers["counts"] if "counts" in adata.layers else adata.X
    present = [g for g in genes if g in adata.var_names]
    if mask is not None:
        X = X[mask]
    sub = X[:, [adata.var_names.get_loc(g) for g in present]]
    k = np.asarray((sub > 0).sum(axis=1)).ravel()
    return float((k >= min_markers).mean()), int((k >= min_markers).sum())


def mad_bounds(x: np.ndarray, nmads: float) -> tuple[float, float]:
    """Median absolute deviation bounds, the repository's QC rule (mirrors
    analysis/scripts/pipeline_utils.mad_bounds without importing across folders)."""
    med = float(np.median(x))
    mad = float(np.median(np.abs(x - med))) * 1.4826
    return med - nmads * mad, med + nmads * mad


def same_library_enrichment(adata, key="library", n_neighbors=30):
    """Mean over cells of (same-library neighbours / expected under mixing)."""
    conn = adata.obsp["connectivities"].tocsr()
    lib = adata.obs[key].astype(str).to_numpy()
    frac = pd.Series(lib).value_counts(normalize=True)
    vals = []
    for i in range(conn.shape[0]):
        nb = conn.indices[conn.indptr[i]:conn.indptr[i + 1]]
        if len(nb) == 0:
            continue
        same = float((lib[nb] == lib[i]).mean())
        vals.append(same / frac[lib[i]])
    return float(np.mean(vals))


def palette() -> dict:
    return json.loads(PALETTE.read_text(encoding="utf-8"))


def apply_style(plt) -> None:
    """The repository's validated palette, applied the way viz_style does it."""
    P = palette()
    plt.rcParams.update({
        "figure.facecolor": P["surface"], "axes.facecolor": P["surface"],
        "savefig.facecolor": P["surface"], "text.color": P["ink"],
        "axes.labelcolor": P["ink"], "xtick.color": P["ink_2"], "ytick.color": P["ink_2"],
        "axes.edgecolor": P["axis"], "axes.spines.top": False, "axes.spines.right": False,
        "grid.color": P["grid"], "font.size": 9, "axes.titlesize": 10, "legend.frameon": False,
    })


STATE_COLOURS_SLOTS = {"hAT2": 1, "pAT2": 2, "DATP": 3, "AT1": 4}   # cAT2 and others take muted


def state_colour(state: str) -> str:
    P = palette()
    slot = STATE_COLOURS_SLOTS.get(state)
    return P["categorical"][str(slot)] if slot else P["muted"]
