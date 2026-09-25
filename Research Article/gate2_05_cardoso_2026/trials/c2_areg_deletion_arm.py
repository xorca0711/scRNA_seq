#!/usr/bin/env python
"""Trial C2: does deleting the ligand collapse the niche, and does anything survive it (Gate 2b).

The owner's Gate 2 branch (b). GSE316244 is the Areg-flox arm: mutant AT2
cells with one or two floxed Areg alleles, profiled as a niche sort
(mesenchyme plus immune) and as an RFP+ epithelial sort. The paper's claim is
that removing the ligand removes the fibrotic fibroblast population and the
reprogrammed epithelial states. Two questions follow, and the second is the
one the paper does not ask: does the collapse reproduce at transcriptome
level, and **is there any part of the programme that does not collapse**.

Design constraint from trial C0, binding on every number here: one library per
genotype per sort, three mice pooled per library, one time point. Directions
are read; nothing is tested; no P value is computed for a genotype contrast.

Frozen rules (written to the run record before any matrix is read):

* Part A, within GSE316244. The two niche libraries and the two RFP libraries
  are processed separately, because they are different sorts and a joint
  embedding would confound sort with genotype.
* QC, doublets, normalisation, features, embedding and clustering follow trial
  C1 exactly (repository MAD rule per library; Scrublet per library with the
  10x prior and the automatic-threshold sanity check; log1p(CP10K); seurat_v3
  2,000 HVGs on counts; scale to 10; 30 PCs; k = 30; Leiden igraph 0.5
  primary with 0.2 and 1.0 as sensitivities; seeds 0). No batch correction:
  library is genotype.
* The BSD feature (the reporter construct's selection marker, the only
  non-gene feature in the deposit) is removed from the matrix by the C0 rule
  and carried per cell, and its detection is reported per cluster as a check
  on the sort. It is a transgene contig with low counts and is used as a
  check, never as a cell-type call.
* Population calls use the paper's own marker sets from
  cardoso_2026_extracts.json, scored with score_genes (ctrl_size 50, n_bins
  25, seed 0) and assigned by modal per-cell argmax; a call is confident only
  if the mode holds at least 50% of the cluster.
* PRE-REGISTERED DIRECTIONS, taken from the paper before these data were
  opened. In Areg-flox/flox relative to Areg-flox/+: (1) the reprogrammed
  fibroblast share of fibroblasts falls (the paper reports 9.2% in flox/+ and
  a share small enough to be omitted at its 5% display floor in flox/flox);
  (2) the DATP-like share of RFP+ cells falls; (3) the Cd177+ share falls;
  (4) the AT2 share rises; (5) in alveolar macrophages Cxcl2, Ccl6 and Ccl9
  fall and MHC-II genes (H2-Ab1, H2-Eb1, Cd74) rise. Each direction is
  recorded as met or not met. A direction that is met is consistent with the
  paper; it is not a test of it, because n is one library per genotype.
* Part B, the question the paper does not ask. The fibrotic programme is
  scored in fibroblasts of four libraries that share a gene space exactly
  (C0 gene space 1): Confetti and Red2Kras mesenchyme from GSE316241, and the
  flox/+ and flox/flox niche libraries here. Scores are compared, and every
  gene of the reprogrammed-fibroblast set is compared as a detection fraction.
  A gene is called **Areg-independent** if its detection in flox/flox
  fibroblasts stays at or above 80% of its detection in flox/+ fibroblasts
  while also exceeding its detection in Confetti fibroblasts by at least 5
  percentage points. This is the frozen definition; it identifies the part of
  the fibrotic programme that survives ligand deletion.
* Part B is SECONDARY and its confounds are named in the output: the two
  series differ in CellRanger version (8.0.0 versus 7.2.0), in sort
  (mesenchyme alone versus mesenchyme plus immune at 1:1) and in tamoxifen
  dosing (three doses versus two). Score comparisons across them carry those
  differences. No joint embedding is built, precisely because a joint
  embedding would present the confound as biology.
* Unit: the library. No P value for any genotype contrast.
"""

from __future__ import annotations

import gc
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import (ENSEMBL_ID_RE, RAW, REPO, RunRecord,  # noqa: E402
                           df_to_markdown, read_mtx_triplet)

sys.path.insert(0, str(REPO / "analysis" / "scripts"))
from pipeline_utils import apply_thresholds, derive_thresholds, make_unique  # noqa: E402

OUT = HERE / "c2_areg_deletion_arm"
OUT.mkdir(exist_ok=True)
EXTRACTS = json.loads((HERE.parent / "cardoso_2026_extracts.json").read_text(encoding="utf-8"))
SETS = EXTRACTS["marker_sets"]

NICHE = [
    {"accession": "GSE316244", "gsm": "GSM9447779", "library": "Expt3_Het_niche",
     "genotype": "Areg-flox/+", "sort": "niche"},
    {"accession": "GSE316244", "gsm": "GSM9447780", "library": "Expt3_Hom_niche",
     "genotype": "Areg-flox/flox", "sort": "niche"},
]
RFP = [
    {"accession": "GSE316244", "gsm": "GSM9447781", "library": "Expt3_Het_RFP",
     "genotype": "Areg-flox/+", "sort": "RFP+ epithelium"},
    {"accession": "GSE316244", "gsm": "GSM9447782", "library": "Expt3_Hom_RFP",
     "genotype": "Areg-flox/flox", "sort": "RFP+ epithelium"},
]
MESENCHYME_REF = [
    {"accession": "GSE316241", "gsm": "GSM9447763", "library": "Expt1_Confetti_mesenchyme",
     "genotype": "Confetti", "sort": "mesenchyme"},
    {"accession": "GSE316241", "gsm": "GSM9447764", "library": "Expt1_Red2Kras_mesenchyme",
     "genotype": "Red2Kras", "sort": "mesenchyme"},
]

NICHE_SETS = {
    "alveolar fibroblast": SETS["alveolar_fibroblast"],
    "adventitial fibroblast": SETS["adventitial_fibroblast"],
    "reprogrammed fibroblast": SETS["reprogrammed_fibrotic_fibroblast"],
    "peri-bronchial fibroblast": SETS["peri_bronchial_fibroblast"],
    "smooth muscle": SETS["smooth_muscle"],
    "pericyte": SETS["pericyte"],
    "mesothelium": SETS["mesothelium"],
    "proliferating": SETS["proliferating"],
    "alveolar macrophage": SETS["alveolar_macrophage"],
    "reprogrammed alveolar macrophage": SETS["reprogrammed_alveolar_macrophage"],
    "monocyte/interstitial macrophage": SETS["monocyte_interstitial_macrophage"],
}
EPI_SETS = {k: v for k, v in SETS["mutant_epithelial_states_fig4l"].items()}
FIBROBLAST_CALLS = {"alveolar fibroblast", "adventitial fibroblast", "reprogrammed fibroblast"}

PRE_REGISTERED = {
    "reprogrammed_fibroblast_share_falls_in_hom": {"paper": "9.2% in flox/+, below the 5% display floor in flox/flox"},
    "DATP_like_share_falls_in_hom": {"paper": "Fig. 4m"},
    "Cd177_positive_share_falls_in_hom": {"paper": "Fig. 4m"},
    "AT2_share_rises_in_hom": {"paper": "Fig. 4m"},
    "AM_inflammatory_falls_and_MHCII_rises_in_hom": {
        "paper": "Extended Data Fig. 10e", "inflammatory": ["Cxcl2", "Ccl6", "Ccl9"],
        "mhc_ii": ["H2-Ab1", "H2-Eb1", "Cd74"]},
}
AREG_INDEPENDENT = {"retains_at_least_fraction_of_het": 0.80,
                    "exceeds_confetti_by_at_least_points": 5.0}

RULES = {
    "accessions": ["GSE316244 (primary)", "GSE316241 (Part B comparison only)"],
    "design_constraint": ("one library per genotype per sort, three mice pooled per library, one time "
                          "point; directions are read, nothing is tested, no P value"),
    "processing": "identical to trial C1 (MAD QC, Scrublet, log1p CP10K, seurat_v3 2000 HVG, 30 PCs, k 30, Leiden igraph 0.5/0.2/1.0, seeds 0)",
    "batch_correction": "none; library is genotype",
    "non_gene_feature": "BSD removed from the matrix and carried per cell as a sort check",
    "sorts_processed_separately": True,
    "population_sets": {"niche": NICHE_SETS, "epithelium": EPI_SETS},
    "pre_registered_directions": PRE_REGISTERED,
    "areg_independent_definition": AREG_INDEPENDENT,
    "part_b_confounds": ["CellRanger 8.0.0 versus 7.2.0", "mesenchyme-only versus 1:1 mesenchyme and immune sort",
                         "three versus two tamoxifen doses"],
    "unit": "the library",
}


def load(entry: dict, rec: RunRecord):
    import anndata as ad
    import scanpy as sc

    root = RAW / entry["accession"] / f"{entry['accession']}_RAW"
    stem = f"{entry['gsm']}_{entry['library']}_"
    paths = {k: root / f"{stem}{k}" for k in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")}
    for path in paths.values():
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(paths["matrix.mtx.gz"], paths["features.tsv.gz"], paths["barcodes.tsv.gz"])

    is_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    bsd = None
    if (~is_gene).any():
        cols = np.where(~is_gene)[0]
        bsd = np.asarray(X[:, cols].sum(axis=1)).ravel()
    X, var = X[:, is_gene], var.loc[is_gene].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))

    obs = pd.DataFrame(index=pd.Index([f"{entry['library']}_{b}" for b in bc]))
    for key in ("library", "genotype", "sort", "accession"):
        obs[key] = entry[key]
    if bsd is not None:
        obs["BSD_counts"] = bsd
        obs["BSD_detected"] = bsd > 0
    adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)

    symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
    adata.var["ribo"] = (symbols.str.startswith(("Rps", "Rpl"))
                         & ~symbols.str.contains("-ps", case=False)).to_numpy()
    adata.var["hb"] = symbols.str.startswith(("Hba", "Hbb")).to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=[20],
                               log1p=True, inplace=True)
    th = derive_thresholds(adata.obs, entry["library"])
    keep = apply_thresholds(adata.obs, th).to_numpy()
    row = {"library": entry["library"], "genotype": entry["genotype"], "sort": entry["sort"],
           "barcodes": int(adata.n_obs), "kept_after_qc": int(keep.sum()),
           "median_genes": int(np.median(adata.obs["n_genes_by_counts"]))}
    adata = adata[keep].copy()

    expected = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=expected)
    auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
    method = "scrublet automatic threshold"
    if auto < 0.2 * expected or auto > 3 * expected:
        cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - expected))
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
        method = f"top {100 * expected:.2f}% of scores (automatic rejected: it called {100 * auto:.2f}%)"
    row.update({"doublet_call_method": method,
                "doublets_removed": int(adata.obs["predicted_doublet"].sum())})
    adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
    row["cells_analysed"] = int(adata.n_obs)
    if "BSD_detected" in adata.obs:
        row["BSD_detected_fraction"] = round(float(adata.obs["BSD_detected"].mean()), 4)
    return adata, row


def embed_and_call(adata, sets: dict, prefix: str):
    import scanpy as sc

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=30, random_state=0)
    adata.obsm["X_pca"] = emb.obsm["X_pca"]
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30, use_rep="X_pca", random_state=0)
    for res in (0.2, 0.5, 1.0):
        sc.tl.leiden(adata, resolution=res, key_added=f"{prefix}_{res}", flavor="igraph",
                     n_iterations=2, random_state=0, directed=False)
    sc.tl.umap(adata, random_state=0)
    del emb
    gc.collect()

    scores = {}
    for name, genes in sets.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"score_{name}", ctrl_size=50,
                              n_bins=25, random_state=0)
            scores[name] = adata.obs[f"score_{name}"].to_numpy()
    frame = pd.DataFrame(scores, index=adata.obs_names)
    adata.obs["call_argmax"] = frame.idxmax(axis=1).to_numpy()
    return adata


def cluster_table(adata, key: str, genotypes: list[str]) -> pd.DataFrame:
    rows = []
    for cluster in sorted(adata.obs[key].astype(str).unique(), key=int):
        mask = adata.obs[key].astype(str).to_numpy() == cluster
        calls = adata.obs.loc[mask, "call_argmax"].value_counts()
        mode_fraction = float(calls.iloc[0] / mask.sum())
        row = {"cluster": cluster, "n_cells": int(mask.sum()), "call": calls.index[0],
               "mode_fraction": round(mode_fraction, 3), "confident": mode_fraction >= 0.50}
        shares = adata.obs.loc[mask, "genotype"].value_counts()
        for genotype in genotypes:
            row[f"n_{genotype}"] = int(shares.get(genotype, 0))
        if "BSD_detected" in adata.obs:
            row["BSD_detected_fraction"] = round(float(adata.obs.loc[mask, "BSD_detected"].mean()), 3)
        rows.append(row)
    return pd.DataFrame(rows)


def composition(table: pd.DataFrame, genotypes: list[str], label: str) -> pd.DataFrame:
    """Per-genotype share of each confident call. Shares within a library only."""
    rows = []
    for genotype in genotypes:
        column = f"n_{genotype}"
        total = table[column].sum()
        for call in sorted(table.loc[table["confident"], "call"].unique()):
            n = int(table.loc[table["confident"] & (table["call"] == call), column].sum())
            rows.append({"compartment": label, "genotype": genotype, "call": call,
                         "n_cells": n, "pct_of_library": round(100 * n / total, 2) if total else np.nan})
    return pd.DataFrame(rows)


def detection_by_group(adata, genes: list[str], group_key: str) -> pd.DataFrame:
    X = adata.layers["counts"]
    groups = adata.obs[group_key].astype(str).to_numpy()
    rows = []
    for gene in genes:
        if gene not in adata.var_names:
            rows.append({"gene": gene, "absent_from_reference": True})
            continue
        col = np.asarray(X[:, adata.var_names.get_loc(gene)].todense()).ravel() > 0
        row = {"gene": gene, "absent_from_reference": False}
        for group in sorted(set(groups)):
            row[group] = round(float(col[groups == group].mean()), 4)
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    import anndata as ad
    import scanpy as sc
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc.settings.verbosity = 0
    rec = RunRecord(OUT / "c2_run_record.json",
                    "C2 Areg deletion arm: does the niche collapse, and what survives (Gate 2b)",
                    RULES, notes="rules frozen before any matrix was read")

    qc_rows, objects = [], {}
    for name, entries in (("niche", NICHE), ("rfp", RFP)):
        parts = []
        for entry in entries:
            adata, row = load(entry, rec)
            parts.append(adata)
            qc_rows.append(row)
            print(f"  {entry['library']}: {row['cells_analysed']} cells")
        merged = ad.concat(parts, join="inner")
        merged.var = parts[0].var.loc[merged.var_names].copy()
        objects[name] = merged
        del parts
        gc.collect()
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(OUT / "c2_qc_per_library.csv", index=False)
    rec.add_output(OUT / "c2_qc_per_library.csv")

    genotypes = ["Areg-flox/+", "Areg-flox/flox"]
    results = {}

    # ---- Part A1: the niche sort ----------------------------------------
    niche = embed_and_call(objects["niche"], NICHE_SETS, "niche")
    niche_table = cluster_table(niche, "niche_0.5", genotypes)
    niche_table.to_csv(OUT / "c2_niche_clusters.csv", index=False)
    rec.add_output(OUT / "c2_niche_clusters.csv")
    print(df_to_markdown(niche_table, index=False))

    fib_clusters = niche_table.loc[niche_table["confident"]
                                   & niche_table["call"].isin(FIBROBLAST_CALLS), "cluster"].tolist()
    fib_mask = niche.obs["niche_0.5"].astype(str).isin(fib_clusters).to_numpy()
    fib = ad.AnnData(X=niche.layers["counts"][fib_mask].copy(),
                     obs=niche.obs.loc[fib_mask, ["library", "genotype", "sort", "accession",
                                                  "doublet_score"]].copy(),
                     var=niche.var.copy())
    fib = embed_and_call(fib, {k: v for k, v in NICHE_SETS.items() if k in FIBROBLAST_CALLS}, "fib")
    fib_table = cluster_table(fib, "fib_0.5", genotypes)
    fib_table.to_csv(OUT / "c2_fibroblast_subclusters.csv", index=False)
    rec.add_output(OUT / "c2_fibroblast_subclusters.csv")
    fib_comp = composition(fib_table, genotypes, "fibroblasts")
    fib_comp.to_csv(OUT / "c2_fibroblast_composition.csv", index=False)
    rec.add_output(OUT / "c2_fibroblast_composition.csv")
    print(df_to_markdown(fib_comp, index=False))

    def share(frame, genotype, call):
        hit = frame[(frame["genotype"] == genotype) & (frame["call"] == call)]
        return float(hit["pct_of_library"].iloc[0]) if len(hit) else 0.0

    het_repro = share(fib_comp, "Areg-flox/+", "reprogrammed fibroblast")
    hom_repro = share(fib_comp, "Areg-flox/flox", "reprogrammed fibroblast")
    results["reprogrammed_fibroblast_share_falls_in_hom"] = {
        "flox_plus_pct": het_repro, "flox_flox_pct": hom_repro,
        "direction_met": bool(hom_repro < het_repro)}

    # ---- Part A2: the RFP+ epithelial sort -------------------------------
    rfp = embed_and_call(objects["rfp"], EPI_SETS, "epi")
    rfp_table = cluster_table(rfp, "epi_0.5", genotypes)
    rfp_table.to_csv(OUT / "c2_epithelial_clusters.csv", index=False)
    rec.add_output(OUT / "c2_epithelial_clusters.csv")
    epi_comp = composition(rfp_table, genotypes, "RFP+ epithelium")
    epi_comp.to_csv(OUT / "c2_epithelial_composition.csv", index=False)
    rec.add_output(OUT / "c2_epithelial_composition.csv")
    print(df_to_markdown(epi_comp, index=False))

    for key, call in (("DATP_like_share_falls_in_hom", "DATP_like"),
                      ("Cd177_positive_share_falls_in_hom", "Cd177_positive"),
                      ("AT2_share_rises_in_hom", "AT2")):
        het = share(epi_comp, "Areg-flox/+", call)
        hom = share(epi_comp, "Areg-flox/flox", call)
        expected_up = key.endswith("rises_in_hom")
        results[key] = {"flox_plus_pct": het, "flox_flox_pct": hom,
                        "direction_met": bool(hom > het) if expected_up else bool(hom < het)}

    # ---- Part A3: alveolar macrophages in the niche sort -----------------
    am_clusters = niche_table.loc[niche_table["confident"]
                                  & niche_table["call"].isin({"alveolar macrophage",
                                                              "reprogrammed alveolar macrophage"}),
                                  "cluster"].tolist()
    am_direction = {"evaluable": bool(am_clusters)}
    if am_clusters:
        am_mask = niche.obs["niche_0.5"].astype(str).isin(am_clusters).to_numpy()
        am = ad.AnnData(X=niche.layers["counts"][am_mask].copy(),
                        obs=niche.obs.loc[am_mask, ["library", "genotype"]].copy(),
                        var=niche.var.copy())
        am.layers["counts"] = am.X.copy()
        infl = PRE_REGISTERED["AM_inflammatory_falls_and_MHCII_rises_in_hom"]["inflammatory"]
        mhc = PRE_REGISTERED["AM_inflammatory_falls_and_MHCII_rises_in_hom"]["mhc_ii"]
        am_det = detection_by_group(am, infl + mhc, "genotype")
        am_det.to_csv(OUT / "c2_alveolar_macrophage_detection.csv", index=False)
        rec.add_output(OUT / "c2_alveolar_macrophage_detection.csv")
        present = am_det[~am_det["absent_from_reference"]].set_index("gene")
        if {"Areg-flox/+", "Areg-flox/flox"}.issubset(present.columns):
            infl_hit = [g for g in infl if g in present.index
                        and present.loc[g, "Areg-flox/flox"] < present.loc[g, "Areg-flox/+"]]
            mhc_hit = [g for g in mhc if g in present.index
                       and present.loc[g, "Areg-flox/flox"] > present.loc[g, "Areg-flox/+"]]
            am_direction.update({"n_cells": int(am.n_obs),
                                 "inflammatory_genes_falling_in_hom": infl_hit,
                                 "mhc_ii_genes_rising_in_hom": mhc_hit,
                                 "direction_met": bool(len(infl_hit) >= 2 and len(mhc_hit) >= 2)})
    results["AM_inflammatory_falls_and_MHCII_rises_in_hom"] = am_direction
    rec.set("pre_registered_directions_outcome", results)

    # ---- Part B: what survives ligand deletion ---------------------------
    ref_parts = []
    for entry in MESENCHYME_REF:
        adata, row = load(entry, rec)
        qc_rows.append(row)
        ref_parts.append(adata)
    ref = ad.concat(ref_parts, join="inner")
    ref.var = ref_parts[0].var.loc[ref.var_names].copy()
    del ref_parts
    gc.collect()
    ref = embed_and_call(ref, NICHE_SETS, "ref")
    ref_table = cluster_table(ref, "ref_0.5", ["Confetti", "Red2Kras"])
    ref_fib = ref_table.loc[ref_table["confident"] & ref_table["call"].isin(FIBROBLAST_CALLS), "cluster"].tolist()
    ref_mask = ref.obs["ref_0.5"].astype(str).isin(ref_fib).to_numpy()

    fibrotic_genes = SETS["reprogrammed_fibrotic_fibroblast"]
    ref_fib_obj = ad.AnnData(X=ref.layers["counts"][ref_mask].copy(),
                             obs=ref.obs.loc[ref_mask, ["library", "genotype"]].copy(),
                             var=ref.var.copy())
    ref_fib_obj.layers["counts"] = ref_fib_obj.X.copy()
    fib.layers["counts"] = fib.layers["counts"]
    det_ref = detection_by_group(ref_fib_obj, fibrotic_genes + SETS["alveolar_fibroblast"], "genotype")
    det_arm = detection_by_group(fib, fibrotic_genes + SETS["alveolar_fibroblast"], "genotype")
    survive = det_ref.merge(det_arm, on=["gene", "absent_from_reference"], how="outer")
    rule = AREG_INDEPENDENT
    calls = []
    for _, row in survive.iterrows():
        het, hom, con = row.get("Areg-flox/+"), row.get("Areg-flox/flox"), row.get("Confetti")
        if any(pd.isna(v) for v in (het, hom, con)) or het == 0:
            calls.append("not evaluable")
            continue
        retained = hom >= rule["retains_at_least_fraction_of_het"] * het
        above = (100 * hom - 100 * con) >= rule["exceeds_confetti_by_at_least_points"]
        calls.append("Areg-independent" if (retained and above) else
                     ("collapses" if hom < het else "unchanged or higher"))
    survive["verdict"] = calls
    survive.to_csv(OUT / "c2_fibrotic_programme_survival.csv", index=False)
    rec.add_output(OUT / "c2_fibrotic_programme_survival.csv")
    rec.set("areg_independent_genes", survive.loc[survive["verdict"] == "Areg-independent", "gene"].tolist())
    rec.set("collapsing_genes", survive.loc[survive["verdict"] == "collapses", "gene"].tolist())
    print(df_to_markdown(survive, index=False))

    # ---- figure and summary ---------------------------------------------
    fig, axes = plt.subplots(1, 3, figsize=(19, 5.6))
    sc.pl.umap(niche, color="niche_0.5", ax=axes[0], show=False, legend_loc="on data", size=8,
               title=f"niche sort, Leiden 0.5 (n={niche.n_obs})")
    sc.pl.umap(fib, color="call_argmax", ax=axes[1], show=False, size=12, title="fibroblasts, call")
    sc.pl.umap(rfp, color="call_argmax", ax=axes[2], show=False, size=8, title="RFP+ epithelium, call")
    fig.tight_layout()
    fig.savefig(OUT / "c2_umaps.png", dpi=110)
    plt.close(fig)
    rec.add_output(OUT / "c2_umaps.png")

    met = {k: v.get("direction_met") for k, v in results.items()}
    lines = [
        "# Trial C2 output: the Areg deletion arm (Gate 2b)", "",
        "Design constraint from trial C0: one library per genotype per sort, three mice pooled per",
        "library, one time point. Directions are read; nothing is tested; no P value is computed.", "",
        f"Pre-registered directions met: **{sum(1 for v in met.values() if v)} of {len(met)}**.", "",
        df_to_markdown(pd.DataFrame([{"pre-registered direction": k, "met": v} for k, v in met.items()]),
                       index=False), "",
        "## QC per library", "", df_to_markdown(pd.DataFrame(qc_rows), index=False), "",
        "## Fibroblast composition by genotype", "", df_to_markdown(fib_comp, index=False), "",
        "## RFP+ epithelial composition by genotype", "", df_to_markdown(epi_comp, index=False), "",
        "## What survives ligand deletion (Part B, secondary)", "",
        "Confounded by CellRanger version, sort composition and tamoxifen dose; see the run record.", "",
        df_to_markdown(survive, index=False), "",
    ]
    (OUT / "c2_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c2_summary.md")
    rec.finish()
    print(f"\ndirections met: {sum(1 for v in met.values() if v)} of {len(met)}")


if __name__ == "__main__":
    main()
