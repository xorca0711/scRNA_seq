#!/usr/bin/env python
"""Trials E2 and E3: which cells make which EGFR ligand in human pulmonary fibrosis.

Extensions 2 and 3 of the plan in DIVERGENCES_AND_NEXT.md, run as one script
over two independent cohorts so that a pattern seen in one can be checked in
the other:

* E2, the IPF Cell Atlas (Adams et al. 2020, GSE136831), 312,928 cells.
* E3, Habermann et al. 2020 (GSE135893), 220,213 cells.

The question. Trial C6 found that the Areg-independent share of mouse Hbegf
sits in endothelium and in fibroblasts rather than in macrophages, and trial
E1 found that in human lung adenocarcinoma HBEGF is highest in neutrophils and
endothelium, again not in myeloid cells. Published work on fibrosis nominates
the macrophage and the transitional epithelial cell instead (Hult et al. 2022,
doi:10.1165/rcmb.2022-0174OC). Human fibrosis is where that prediction was
made, so it is the fair place to test it.

A second, independent prediction is tested at the same time: sustained AREG
from intermediate alveolar stem cells drives fibroblast EGFR and progressive
fibrosis (Zhao et al. 2024, doi:10.1016/j.stem.2024.07.004), which implies
AREG should be higher in the transitional epithelial population than in AT2
cells of the same donor.

Both cohorts deposit one merged sparse matrix, too large to load here, so
`mtx_stream.extract_gene_rows` makes a single pass and keeps only the genes
below. Nothing is clustered or integrated; the deposited annotations group the
cells, which is the reason to use curated atlases for this question at all.

Frozen rules, set before any expression value was read:

* Cells: the deposited per-cell annotation of each cohort, restricted to the
  disease groups named in its config (control and IPF). Other diagnoses in a
  cohort are out of scope.
* Grouping: the deposited cell-type label. No marker gating, because the
  population this trial is about, the transitional or aberrant basaloid
  epithelial state, cannot be gated with a handful of markers.
* Genes: AREG, EREG, HBEGF, TGFA, EGF, BTC, EPGN and EGFR.
* Measures: per donor, disease group and cell type, the fraction of cells with
  a non-zero count. Counts are raw here, so a detection fraction is the
  comparable measure and mean expression is not reported.
* Floors: a donor contributes a value for a cell type only if it holds at
  least 20 cells of it. A paired test runs only with at least 5 donors.
* T1, the published prediction. In IPF, HBEGF detection is higher in the
  myeloid group than in the transitional epithelial group, paired within
  donor, Wilcoxon signed-rank, two-sided, alpha 0.05.
* T2, the ranking. In IPF, the cell types with the highest median per-donor
  HBEGF detection, reported as a table rather than a single winner.
* T3, the Zhao prediction. In IPF, AREG detection is higher in the
  transitional epithelial group than in AT2 cells, paired within donor, same
  test.
* T4, disease against control. For the same cell types, median per-donor
  detection in IPF and in control, reported side by side, no test, because
  the groups are unpaired and cohort composition differs.
* The unit is the donor throughout. Detection fractions depend on depth, and
  the two cohorts differ in depth and in annotation vocabulary, so numbers are
  compared within a cohort and only the direction is compared across them.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RAW, RunRecord, df_to_markdown  # noqa: E402
from mtx_stream import extract_gene_rows, read_lines  # noqa: E402

OUT = HERE / "e2_human_fibrosis_ligand_sources"
OUT.mkdir(exist_ok=True)

GENES = ["AREG", "EREG", "HBEGF", "TGFA", "EGF", "BTC", "EPGN", "EGFR"]
MIN_CELLS = 20
MIN_DONORS = 5

CONFIGS = {
    "E2_GSE136831": {
        "series": "GSE136831",
        "matrix": "GSE136831_RawCounts_Sparse.mtx.gz",
        "genes": "GSE136831_AllCells.GeneIDs.txt.gz",
        "gene_column": 1, "gene_header": True,
        "barcodes": "GSE136831_AllCells.cellBarcodes.txt.gz",
        "metadata": "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
        "metadata_sep": "\t", "barcode_column": "CellBarcode_Identity",
        "celltype_column": "Manuscript_Identity", "donor_column": "Subject_Identity",
        "disease_column": "Disease_Identity", "groups": ["Control", "IPF"], "disease": "IPF",
        "transitional": "Aberrant_Basaloid", "at2": "ATII",
        "myeloid_keywords": ["Macrophage", "Monocyte", "cDC", "pDC", "Mast", "DC_"],
        "transitional_extra": None,
    },
    "E3_GSE135893": {
        "series": "GSE135893",
        "matrix": "GSE135893_matrix.mtx.gz",
        "genes": "GSE135893_genes.tsv.gz",
        "gene_column": 0, "gene_header": False,
        "barcodes": "GSE135893_barcodes.tsv.gz",
        "metadata": "GSE135893_IPF_metadata.csv.gz",
        "metadata_sep": ",", "barcode_column": None,
        "celltype_column": "celltype", "donor_column": "Sample_Name",
        "disease_column": "Diagnosis", "groups": ["Control", "IPF"], "disease": "IPF",
        "transitional": "KRT5-/KRT17+", "at2": "AT2",
        "myeloid_keywords": ["Macrophage", "Monocyte", "cDC", "pDC", "Mast"],
        "transitional_extra": "Transitional AT2",
    },
}

RULES = {
    "cohorts": {k: v["series"] for k, v in CONFIGS.items()},
    "grouping": "deposited cell-type label; no marker gating",
    "genes": GENES,
    "measure": "fraction of cells with a non-zero raw count",
    "floors": {"cells_per_donor_celltype": MIN_CELLS, "donors_for_a_test": MIN_DONORS},
    "T1": "IPF, HBEGF higher in myeloid than in the transitional epithelial group, paired within donor, Wilcoxon two-sided",
    "T2": "IPF, cell types ranked by median per-donor HBEGF detection",
    "T3": "IPF, AREG higher in the transitional group than in AT2, paired within donor, Wilcoxon two-sided",
    "T4": "median per-donor detection in IPF and control side by side, no test",
    "unit": "the donor",
    "caveats": ["detection depends on depth", "cohorts differ in depth and in annotation vocabulary",
                "numbers compare within a cohort; only direction compares across cohorts"],
}


def load_cohort(cfg, rec):
    series = RAW / cfg["series"]
    matrix = series / cfg["matrix"]
    for key in ("matrix", "genes", "barcodes", "metadata"):
        rec.add_input(series / cfg[key])
    genes = read_lines(series / cfg["genes"], column=cfg["gene_column"])
    if cfg["gene_header"]:
        genes = genes[1:]
    barcodes = read_lines(series / cfg["barcodes"])
    meta = pd.read_csv(series / cfg["metadata"], sep=cfg["metadata_sep"], compression="gzip")
    key = cfg["barcode_column"] or meta.columns[0]
    meta = meta.set_index(key)
    vectors = extract_gene_rows(matrix, genes, set(GENES))
    return np.array(barcodes), meta, vectors


def per_donor(cfg, barcodes, meta, vectors):
    frame = meta.reindex(barcodes)
    frame["celltype"] = frame[cfg["celltype_column"]].astype(str).str.strip(chr(34))
    frame["donor"] = frame[cfg["donor_column"]].astype(str).str.strip(chr(34))
    frame["disease"] = frame[cfg["disease_column"]].astype(str).str.strip(chr(34))
    keep = frame["disease"].isin(cfg["groups"]).to_numpy() & frame["celltype"].notna().to_numpy()
    kept = np.where(keep)[0]
    groups = frame.iloc[kept].groupby(["donor", "disease", "celltype"]).indices
    records = []
    for (donor, disease, celltype), pos in groups.items():
        where = kept[pos]
        row = {"donor": donor, "disease": disease, "celltype": celltype, "n_cells": len(where)}
        for gene, vec in vectors.items():
            row["det_" + gene] = round(float((vec[where] > 0).mean()), 4)
        records.append(row)
    table = pd.DataFrame(records)
    table["evaluable"] = table["n_cells"] >= MIN_CELLS
    return table


def donor_series(table, gene, celltype=None, keywords=None):
    sub = table[table["evaluable"]]
    if celltype is not None:
        sub = sub[sub["celltype"] == celltype]
    else:
        sub = sub[sub["celltype"].apply(lambda c: any(k in c for k in keywords))]
    hits = sub["det_" + gene] * sub["n_cells"]
    pooled = hits.groupby(sub["donor"]).sum() / sub.groupby("donor")["n_cells"].sum()
    return pooled


def compare(a, b, label_a, label_b):
    from scipy.stats import wilcoxon
    joined = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    out = {"group_a": label_a, "group_b": label_b, "n_donors": int(len(joined))}
    if len(joined) == 0:
        return out
    out["median_a"] = round(float(joined["a"].median()), 4)
    out["median_b"] = round(float(joined["b"].median()), 4)
    out["a_higher_in_donors"] = int((joined["a"] > joined["b"]).sum())
    if len(joined) >= MIN_DONORS:
        stat = wilcoxon(joined["a"], joined["b"])
        out["p_value"] = round(float(stat.pvalue), 6)
        out["a_higher"] = bool(out["p_value"] < 0.05 and joined["a"].median() > joined["b"].median())
    else:
        out["p_value"] = None
        out["a_higher"] = None
        out["note"] = "fewer donors than the frozen floor; no test"
    return out


def main():
    rec = RunRecord(OUT / "e2_run_record.json",
                    "E2 and E3 EGFR ligand sources in human pulmonary fibrosis", RULES)
    all_tables, comparisons, rankings = [], [], []
    for name, cfg in CONFIGS.items():
        barcodes, meta, vectors = load_cohort(cfg, rec)
        table = per_donor(cfg, barcodes, meta, vectors)
        table.insert(0, "cohort", name)
        all_tables.append(table)
        ipf = table[table["disease"] == cfg["disease"]]

        myeloid = donor_series(ipf, "HBEGF", keywords=cfg["myeloid_keywords"])
        transitional = donor_series(ipf, "HBEGF", celltype=cfg["transitional"])
        row = compare(myeloid, transitional, "myeloid (pooled)", cfg["transitional"])
        row.update({"cohort": name, "gene": "HBEGF", "test": "T1"})
        comparisons.append(row)

        areg_trans = donor_series(ipf, "AREG", celltype=cfg["transitional"])
        areg_at2 = donor_series(ipf, "AREG", celltype=cfg["at2"])
        row = compare(areg_trans, areg_at2, cfg["transitional"], cfg["at2"])
        row.update({"cohort": name, "gene": "AREG", "test": "T3"})
        comparisons.append(row)

        if cfg["transitional_extra"]:
            extra = donor_series(ipf, "AREG", celltype=cfg["transitional_extra"])
            row = compare(extra, areg_at2, cfg["transitional_extra"], cfg["at2"])
            row.update({"cohort": name, "gene": "AREG", "test": "T3 extra"})
            comparisons.append(row)

        ok = ipf[ipf["evaluable"]]
        rank = ok.groupby("celltype")[["det_HBEGF", "det_AREG", "det_EREG", "det_EGFR"]].median().round(4)
        rank["n_donors"] = ok.groupby("celltype").size()
        rank = rank[rank["n_donors"] >= 3].sort_values("det_HBEGF", ascending=False)
        rank.insert(0, "cohort", name)
        rankings.append(rank.reset_index())
        print(name, "cells with a label in scope:", int(table["n_cells"].sum()))

    table = pd.concat(all_tables, ignore_index=True)
    table.to_csv(OUT / "e2_per_donor_celltype.csv", index=False)
    rec.add_output(OUT / "e2_per_donor_celltype.csv")
    comp = pd.DataFrame(comparisons)
    comp.to_csv(OUT / "e2_paired_comparisons.csv", index=False)
    rec.add_output(OUT / "e2_paired_comparisons.csv")
    rank = pd.concat(rankings, ignore_index=True)
    rank.to_csv(OUT / "e2_ipf_ranking_by_celltype.csv", index=False)
    rec.add_output(OUT / "e2_ipf_ranking_by_celltype.csv")
    rec.set("comparisons", comparisons)
    rec.set("hbegf_top_celltype_per_cohort",
            {c: str(rank[rank["cohort"] == c].iloc[0]["celltype"]) for c in rank["cohort"].unique()})

    lines = ["# Trials E2 and E3: EGFR ligand sources in human pulmonary fibrosis", "",
             "## Pre-registered paired comparisons, IPF donors", "",
             df_to_markdown(comp, index=False), "",
             "## Cell types ranked by median per-donor HBEGF detection, IPF", "",
             df_to_markdown(rank.head(30), index=False), ""]
    (OUT / "e2_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e2_summary.md")
    rec.finish()
    print(df_to_markdown(comp, index=False))


if __name__ == "__main__":
    main()
