#!/usr/bin/env python
"""Trial G2: gene set enrichment, IPF against control, per donor, within three
compartments of GSE136831, with GSE135893 held out.

The second trial G0 licensed and the strongest design in this repository: the
donor is the unit, 32 against 28 donors in the discovery cohort, 12 against 10
in the held-out one. Every statistic is on per-donor pseudobulks.

Rules, frozen before any count is read
--------------------------------------
R1  Unit: the donor. A donor enters a compartment only with at least 50 cells in
    it (the E6 floor). An arm with fewer than 3 donors stops the contrast.
R2  Cohorts and arms: discovery = GSE136831 (Adams et al. 2020), Disease_Identity
    IPF against Control, COPD out of scope as in E2; held-out = GSE135893
    (Habermann et al. 2020), Diagnosis IPF against Control, the other diagnoses
    out of scope.
R3  Compartments, by the deposited labels (Manuscript_Identity; celltype):
      AT2          = {ATII}                         ; {AT2}
      fibroblasts  = {Fibroblast, Myofibroblast}    ; {Fibroblasts, Myofibroblasts, HAS1 High Fibroblasts, PLIN2+ Fibroblasts}
      macrophages  = {Macrophage, Macrophage_Alveolar}; {Macrophages, Proliferating Macrophages}
    The transitional and aberrant-basaloid populations are NOT pooled with AT2,
    because they are the object of question A5 and would blur the contrast.
R4  Ranking: genes with at least 10 counts in at least 3 donors of the two arms;
    log2 CPM on the compartment pseudobulk; Welch t, IPF minus control; ties
    broken by mean difference then by name. Human symbols from the deposited
    gene tables; a symbol present in both cohorts is required for replication.
R5  Gene sets: MSigDB 2024.1.Hs hallmark and GO biological process, sha256 in
    the run record; 15 to 500 genes after intersection with the tested genes.
R6  Positive control, amended before the run: the plan named the Tsukui
    fibrotic fibroblast set (reference M6), which is not on disk as a list of
    fifteen or more genes; HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION up in IPF
    fibroblasts at FDR < 0.05 in the discovery cohort takes its place, as the
    activated-fibroblast matrix programme every IPF single-cell study reports.
    If it fails, every compartment of the discovery cohort is unreadable.
R7  Statistic: gseapy pre-ranked GSEA, weight 1, 1000 permutations, seed 0; every
    hallmark enrichment score reproduced in-house to 1e-6; every set at FDR < 0.05
    in a discovery compartment is drawn against 500 expression-matched random
    sets and clears only if its matched p is below 0.05.
R8  Replication and reading, fixed now: a hallmark or GO set that clears both
    nulls in a discovery compartment is Exploratory; if it also clears gseapy
    FDR < 0.05 in the same direction in the same compartment of the held-out
    cohort, the proposed status is Validated (held-out replication, donor as the
    unit); a set that clears in discovery and not in the held-out cohort is Not
    established, whatever its discovery FDR. Nothing is read from the held-out
    cohort that did not first clear in discovery.

Outputs (g2_gsea_ipf/): g2_units.csv, g2_hallmark_<cohort>_<compartment>.csv,
g2_positive_control.csv, g2_replication.csv (every discovery set that cleared,
with its held-out result), g2_summary.md, g2_run_record.json; rankings and full
GO tables as csv.gz (ignored, regenerable).
"""
from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(REPO / "Thesis" / "gate1_04_sikkema_2023_hlca" / "trials"))
from trial_utils import RunRecord, df_to_markdown  # noqa: E402
import gsea_utils as gu  # noqa: E402

RAW = REPO / "raw_data"
OUT = HERE / "g2_gsea_ipf"
SEED = 0
CELL_FLOOR = 50
ARM_FLOOR = 3
GMT = {"hallmark": gu.MSIGDB / "h.all.v2024.1.Hs.symbols.gmt", "gobp": gu.MSIGDB / "c5.go.bp.v2024.1.Hs.symbols.gmt"}
COMPARTMENTS = {
    "AT2": ({"ATII"}, {"AT2"}),
    "fibroblasts": ({"Fibroblast", "Myofibroblast"}, {"Fibroblasts", "Myofibroblasts", "HAS1 High Fibroblasts", "PLIN2+ Fibroblasts"}),
    "macrophages": ({"Macrophage", "Macrophage_Alveolar"}, {"Macrophages", "Proliferating Macrophages"}),
}
PC_SET = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
RULES = {
    "R1_unit": "donor; at least 50 cells per donor per compartment; at least 3 donors per arm",
    "R2_cohorts": {"discovery": "GSE136831 IPF vs Control (COPD out of scope)", "held_out": "GSE135893 IPF vs Control (other diagnoses out of scope)"},
    "R3_compartments": {k: {"GSE136831": sorted(v[0]), "GSE135893": sorted(v[1])} for k, v in COMPARTMENTS.items()},
    "R4_ranking": "genes with >= 10 counts in >= 3 donors; log2 CPM; Welch t IPF minus control; deterministic ties",
    "R5_gene_sets": "MSigDB 2024.1.Hs hallmark and GO BP, 15 to 500 genes after intersection",
    "R6_positive_control": f"{PC_SET} up in IPF fibroblasts, discovery cohort, FDR < 0.05 (amendment: the Tsukui set is not on disk as a list)",
    "R7_statistic": "gseapy prerank, weight 1, 1000 permutations, seed 0; in-house ES agreement to 1e-6; matched null 500 draws for discovery sets at FDR < 0.05",
    "R8_reading": "both nulls in discovery: Exploratory; plus same-direction FDR < 0.05 in the held-out compartment: Validated (held-out replication); discovery only: Not established",
}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def cohort_136831(rec):
    p_meta = RAW / "GSE136831" / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
    p_genes = RAW / "GSE136831" / "GSE136831_AllCells.GeneIDs.txt.gz"
    p_mtx = RAW / "GSE136831" / "GSE136831_RawCounts_Sparse.mtx.gz"
    for p in (p_meta, p_genes, p_mtx):
        rec.add_input(p)
    m = pd.read_csv(p_meta, sep="\t", usecols=["CellBarcode_Identity", "Manuscript_Identity", "Disease_Identity", "Subject_Identity"])
    g = pd.read_csv(p_genes, sep="\t")
    genes = g.iloc[:, 1].astype(str).to_numpy()
    return m.rename(columns={"Manuscript_Identity": "label", "Disease_Identity": "disease", "Subject_Identity": "donor"}), genes, p_mtx, 0


def cohort_135893(rec):
    p_meta = RAW / "GSE135893" / "GSE135893_IPF_metadata.csv.gz"
    p_genes = RAW / "GSE135893" / "GSE135893_genes.tsv.gz"
    p_bc = RAW / "GSE135893" / "GSE135893_barcodes.tsv.gz"
    p_mtx = RAW / "GSE135893" / "GSE135893_matrix.mtx.gz"
    for p in (p_meta, p_genes, p_bc, p_mtx):
        rec.add_input(p)
    h = pd.read_csv(p_meta, usecols=["Unnamed: 0", "Diagnosis", "Sample_Name", "celltype"])
    h = h.rename(columns={"Unnamed: 0": "barcode", "Diagnosis": "disease", "Sample_Name": "donor", "celltype": "label"})
    with gzip.open(p_bc, "rt") as fh:
        barcodes = [line.strip() for line in fh]
    with gzip.open(p_genes, "rt") as fh:
        genes = np.array([line.strip().split("\t")[-1] for line in fh])
    # align the metadata (a filtered subset) to the matrix column order
    pos = pd.Series(np.arange(len(barcodes)), index=barcodes)
    h["col"] = pos.reindex(h["barcode"]).to_numpy()
    assert h["col"].notna().all(), "a metadata barcode is missing from the matrix"
    full = pd.DataFrame({"label": None, "disease": None, "donor": None}, index=np.arange(len(barcodes)))
    full.loc[h["col"].astype(int).to_numpy(), ["label", "disease", "donor"]] = h[["label", "disease", "donor"]].to_numpy()
    return full.reset_index(drop=True), genes, p_mtx, 1


def pseudobulk_cohort(name, meta, genes, p_mtx, key, rec):
    """Per donor and compartment, IPF and control only, cell floor applied after counting."""
    meta = meta.copy()
    meta["comp"] = None
    for comp, labels in COMPARTMENTS.items():
        meta.loc[meta["label"].isin(labels[key]), "comp"] = comp
    meta["disease"] = meta["disease"].replace({"Control": "control"})
    sel = meta["comp"].notna() & meta["disease"].isin(["IPF", "control"])
    counts = meta[sel].groupby(["comp", "donor", "disease"]).size().rename("cells").reset_index()
    counts["kept"] = counts["cells"] >= CELL_FLOOR
    kept = counts[counts["kept"]].reset_index(drop=True)
    group_id = {(r.comp, r.donor): i for i, r in kept.iterrows()}
    goc = np.full(len(meta), -1, dtype=np.int64)
    idx = np.flatnonzero(sel.to_numpy())
    keys = list(zip(meta.loc[idx, "comp"], meta.loc[idx, "donor"]))
    goc[idx] = [group_id.get(k, -1) for k in keys]
    log(f"{name}: streaming the matrix into {len(kept)} donor-compartments from {int((goc >= 0).sum()):,} cells")
    pb, nnz = gu.pseudobulk_from_mtx(p_mtx, goc, len(kept))
    rec.set(f"{name}_nonzero_entries_read", int(nnz))
    kept["counts"] = pb.sum(0).astype(int)
    kept["cohort"] = name
    return kept, pb.T, genes  # units by genes


def run_cohort(name, kept, pb, genes, sets, rec, rng, targets=None):
    """GSEA per compartment; returns {compartment: {'hallmark': df, 'gobp': df, 'matched': df}}."""
    out = {}
    for comp in COMPARTMENTS:
        sub = kept[kept.comp == comp]
        ipf = sub.index[sub.disease == "IPF"].to_numpy()
        ctl = sub.index[sub.disease == "control"].to_numpy()
        tag = f"{name}_{comp}"
        if len(ipf) < ARM_FLOOR or len(ctl) < ARM_FLOOR:
            rec.set(f"{tag}_verdict", f"stopped: {len(ipf)} IPF against {len(ctl)} control donors")
            log(f"{tag}: stopped")
            continue
        ranking, facts = gu.ranking_between_arms(pb, genes, ipf, ctl)
        ranking.to_csv(OUT / f"g2_ranking_{tag}.csv.gz", index=False, compression="gzip")
        rec.set(f"{tag}_ranking", {**facts, "donors": [int(len(ipf)), int(len(ctl))]})
        log(f"{tag}: {len(ipf)} IPF against {len(ctl)} control donors, {facts['genes_tested']:,} genes")
        # the held-out cohort is asked only about the sets that cleared discovery (R8),
        # which is also what keeps its permutation cost to minutes
        use_sets = sets
        if targets is not None:
            want = targets.get(comp, {})
            use_sets = {coll: {k: v for k, v in gs.items() if k in want.get(coll, set())} for coll, gs in sets.items()}
            rec.set(f"{tag}_sets_asked", {coll: len(v) for coll, v in use_sets.items()})
        empty = pd.DataFrame(columns=["set", "es", "nes", "p_perm", "fdr", "leading_edge", "size"])
        res = {coll: (gu.run_prerank(ranking, gs, SEED, n_perm=1000) if gs else empty) for coll, gs in use_sets.items()}
        scores, gene_arr = ranking["t"].to_numpy(), ranking["gene"].to_numpy()
        worst = max((abs(gu.enrichment_score(scores, np.isin(gene_arr, list(sets["hallmark"][r["set"]]))) - float(r["es"]))
                     for _, r in res["hallmark"].iterrows()), default=0.0)
        assert worst < 1e-6, f"{tag}: in-house ES disagrees with gseapy by {worst}"
        rec.set(f"{tag}_es_agreement_max_abs_diff", float(worst))
        res["hallmark"].to_csv(OUT / f"g2_hallmark_{tag}.csv", index=False)
        rec.add_output(OUT / f"g2_hallmark_{tag}.csv")
        res["gobp"].to_csv(OUT / f"g2_gobp_{tag}.csv.gz", index=False, compression="gzip")
        matched = []
        if targets is None:  # discovery: matched null for every set at FDR < 0.05
            for coll in ("hallmark", "gobp"):
                for _, r in res[coll][res[coll]["fdr"] < 0.05].iterrows():
                    hit = np.isin(gene_arr, list(sets[coll][r["set"]]))
                    null = gu.matched_random_null(scores, ranking["mean_logcpm"].to_numpy(), hit, 500, rng)
                    summ = gu.matched_null_summary(float(r["es"]), null)
                    matched.append({"compartment": comp, "collection": coll, "set": r["set"], "size": r["size"],
                                    "nes": float(r["nes"]), "fdr": float(r["fdr"]), **summ,
                                    "clears_both": bool(summ["matched_p"] < 0.05)})
        out[comp] = {**res, "matched": pd.DataFrame(matched)}
    return out


def main() -> int:
    OUT.mkdir(exist_ok=True)
    rec = RunRecord(OUT / "g2_run_record.json", "G2 gene set enrichment, IPF against control per donor, held-out replication", RULES)
    rec.set("gene_set_files", {k: gu.gene_set_facts(v) for k, v in GMT.items()})
    sets = {k: gu.read_gmt(v) for k, v in GMT.items()}
    rng = np.random.default_rng(SEED)

    units_all = []
    cohorts = {}
    for name, loader in (("GSE136831", cohort_136831), ("GSE135893", cohort_135893)):
        meta, genes, p_mtx, key = loader(rec)
        kept, pb, genes = pseudobulk_cohort(name, meta, genes, p_mtx, key, rec)
        units_all.append(kept)
        cohorts[name] = (kept, pb, genes)
    units = pd.concat(units_all, ignore_index=True)
    units.to_csv(OUT / "g2_units.csv", index=False)
    rec.add_output(OUT / "g2_units.csv")

    disc = run_cohort("GSE136831", *cohorts["GSE136831"], sets, rec, rng)
    targets = {}
    for comp, d in disc.items():
        m = d["matched"]
        if len(m):
            targets[comp] = {coll: set(m[(m.collection == coll) & m.clears_both]["set"]) for coll in ("hallmark", "gobp")}
    # the positive-control set is always asked of the held-out fibroblasts, for the record
    targets.setdefault("fibroblasts", {"hallmark": set(), "gobp": set()})["hallmark"].add(PC_SET)
    held = run_cohort("GSE135893", *cohorts["GSE135893"], sets, rec, rng, targets=targets)

    # positive control
    pc = disc.get("fibroblasts", {}).get("hallmark", pd.DataFrame())
    row = pc[pc["set"] == PC_SET] if len(pc) else pc
    pc_ok = bool(len(row) and row["nes"].iloc[0] > 0 and row["fdr"].iloc[0] < 0.05)
    pd.DataFrame([{"set": PC_SET, "compartment": "fibroblasts", "cohort": "GSE136831",
                   "nes": float(row["nes"].iloc[0]) if len(row) else np.nan,
                   "fdr": float(row["fdr"].iloc[0]) if len(row) else np.nan, "clears": pc_ok}]).to_csv(OUT / "g2_positive_control.csv", index=False)
    rec.add_output(OUT / "g2_positive_control.csv")
    rec.set("positive_control_clears", pc_ok)
    log(f"positive control {'clears' if pc_ok else 'FAILS'}")

    # replication of every discovery set that cleared both nulls
    rows = []
    for comp, d in disc.items():
        m = d["matched"]
        if not len(m):
            continue
        h = held.get(comp, {})
        for _, r in m[m.clears_both].iterrows():
            ht = h.get(r["collection"], pd.DataFrame())
            hr = ht[ht["set"] == r["set"]] if len(ht) else ht
            same = bool(len(hr) and np.sign(hr["nes"].iloc[0]) == np.sign(r["nes"]) and hr["fdr"].iloc[0] < 0.05)
            rows.append({"compartment": comp, "collection": r["collection"], "set": r["set"], "size": r["size"],
                         "discovery_nes": r["nes"], "discovery_fdr": r["fdr"], "matched_p": r["matched_p"],
                         "held_out_nes": float(hr["nes"].iloc[0]) if len(hr) else np.nan,
                         "held_out_fdr": float(hr["fdr"].iloc[0]) if len(hr) else np.nan,
                         "held_out_tested": bool(len(hr)), "replicates": same,
                         "proposed_status": ("Validated (held-out replication)" if same and pc_ok else
                                             "Not established" if pc_ok else "unreadable")})
    rep = pd.DataFrame(rows)
    rep.to_csv(OUT / "g2_replication.csv", index=False)
    rec.add_output(OUT / "g2_replication.csv")
    rec.set("discovery_sets_clearing_both_nulls", int(len(rep)))
    rec.set("replicated", int(rep["replicates"].sum()) if len(rep) else 0)

    # summary
    md = ["# Trial G2: gene set enrichment, IPF against control per donor, held out", "",
          "Generated by `g2_gsea_ipf.py`; rules R1 to R8 in the docstring and the run record, frozen",
          "before any count was read. Unit: the donor. Discovery GSE136831, held-out GSE135893.", "",
          "## Donors per arm after the 50-cell floor", "",
          df_to_markdown(units.groupby(["cohort", "comp", "disease"]).agg(donors=("donor", "nunique"), cells=("cells", "sum")).reset_index(), index=False), "",
          "## Positive control", "",
          f"{PC_SET} up in IPF fibroblasts, discovery: NES {float(row['nes'].iloc[0]):.2f}, FDR {float(row['fdr'].iloc[0]):.3g}, "
          f"{'clears' if pc_ok else 'FAILS; every discovery reading below is unreadable (R6)'}." if len(row) else "not computed.", "",
          "## Discovery sets clearing both nulls, and their held-out result", ""]
    if len(rep):
        show = rep.copy()
        for c in ("discovery_nes", "held_out_nes"):
            show[c] = show[c].round(2)
        for c in ("discovery_fdr", "matched_p", "held_out_fdr"):
            show[c] = show[c].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
        md.append(df_to_markdown(show.sort_values(["compartment", "replicates", "discovery_fdr"], ascending=[True, False, True]), index=False))
        by = rep.groupby("compartment").agg(cleared=("set", "size"), replicated=("replicates", "sum")).reset_index()
        md += ["", df_to_markdown(by, index=False)]
    else:
        md.append("no discovery set cleared both nulls in any compartment")
    md += ["", "## Reading under R8", "",
           f"{int(rep['replicates'].sum()) if len(rep) else 0} of {len(rep)} discovery sets replicate in the held-out cohort in the same "
           "direction in the same compartment; those carry the proposed status Validated (held-out replication, donor as the unit), "
           "the rest Not established. Proposed statuses await the owner's retain or reject."]
    (OUT / "g2_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    rec.add_output(OUT / "g2_summary.md")
    rec.finish()
    log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
