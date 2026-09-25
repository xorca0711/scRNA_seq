#!/usr/bin/env python
"""Trial E1: which human cells make which EGFR ligand in early lung adenocarcinoma.

Extension 1 of the plan in DIVERGENCES_AND_NEXT.md, on the dataset the paper
itself used for its human comparison (Kim et al. 2020, Nat Commun,
doi:10.1038/s41467-020-16164-1; GEO GSE131907).

The question. In the mouse deposit, Areg is epithelial and Hbegf is broadly
expressed with its Areg-independent share sitting in endothelium and in the
fibroblasts themselves rather than in macrophages (trial C6). Published work
predicts the opposite emphasis for HBEGF in human disease: myeloid cells and
transitional epithelium in fibrosis (Hult et al. 2022,
doi:10.1165/rcmb.2022-0174OC), and tumour plus myeloid cells in lung
adenocarcinoma (Van Hiep et al. 2022, doi:10.3389/fonc.2022.963896). This
trial asks which pattern human primary lung adenocarcinoma shows.

**Why this dataset changes what is statistically possible.** Every earlier
trial in this folder was limited to one library per genotype, so nothing could
be tested. Here the unit is the donor and there are roughly a dozen per
tissue, so a paired comparison within donor is admissible and one is
pre-registered below. It is the first testable claim in this paper's folder.

Frozen rules, set before any expression value was read:

* Cells: the deposited annotation, restricted to Sample_Origin in {tLung,
  nLung}. Primary tumour lung and normal lung only; lymph node, brain
  metastasis and pleural effusion samples are out of scope.
* Compartment, primary: marker gates on the expression matrix, in the same
  form as trial C6 and in this priority so each cell lands once: immune
  (PTPRC), endothelial (PECAM1 or CDH5), epithelial (EPCAM or KRT8 or KRT18),
  mesenchymal (COL1A1 or COL1A2), otherwise unassigned. Immune cells are split
  by priority into neutrophil (S100A8, S100A9, FCGR3B), myeloid (LYZ, CD68,
  ITGAM, MARCO, CD14), lymphoid (CD3E, CD79A, NKG7) and other immune.
* Compartment, cross-check: the deposited `Cell_type.refined` column, used
  only to report agreement, never to select cells.
* Genes: AREG, EREG, HBEGF, TGFA, EGF, BTC, EPGN and the receptor EGFR.
* Measures per donor, tissue and compartment: the fraction of cells with a
  non-zero value, and the mean of the deposited log2 TPM values. A compartment
  contributes a donor-level value only if it holds at least 50 cells in that
  donor and tissue.
* T1, THE TEST. In tumour lung, the AREG detection fraction of epithelial
  cells exceeds that of myeloid cells, paired within donor. Wilcoxon
  signed-rank, two-sided, alpha 0.05, on donors that clear the 50-cell floor
  in both compartments. The unit is the donor.
* T2, descriptive. In tumour lung, which compartment has the highest median
  per-donor HBEGF detection. Three pre-named outcomes: myeloid-dominant
  supports the published human emphasis; endothelial or mesenchymal-dominant
  matches trial C6 in mouse; epithelial-dominant matches neither and mirrors
  the mouse epithelium instead.
* T3, descriptive. The same table for EGFR, to see whether the receptor sits
  where the mouse data put it.
* REFUTATION RULE, pre-registered because the owner asked for the remaining
  extensions to run only if this one does not refute. E1 refutes the working
  picture if T1 fails, meaning AREG is not higher in epithelial than in
  myeloid cells, or if AREG is detected in under 1% of tumour epithelial
  cells, meaning the dataset cannot address the question. Either outcome
  stops E2 to E4 until the owner decides.
* Caveats carried into every reading: marker gates are not cell-type calls;
  detection fractions depend on depth; tumour epithelium here is malignant
  tissue rather than the earliest mutant state the mouse model captures.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RAW, RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "e1_human_luad_ligand_sources"
OUT.mkdir(exist_ok=True)
SERIES = RAW / "GSE131907"
MATRIX = SERIES / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
ANNOTATION = SERIES / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
CACHE = SERIES / "e1_extracted_rows.npz"   # regenerable, kept out of the tracked tree

LIGANDS = ["AREG", "EREG", "HBEGF", "TGFA", "EGF", "BTC", "EPGN"]
RECEPTOR = ["EGFR"]
GATES = {
    "immune": ["PTPRC"],
    "endothelial": ["PECAM1", "CDH5"],
    "epithelial": ["EPCAM", "KRT8", "KRT18"],
    "mesenchymal": ["COL1A1", "COL1A2"],
}
IMMUNE_SPLIT = {
    "neutrophil": ["S100A8", "S100A9", "FCGR3B"],
    "myeloid": ["LYZ", "CD68", "ITGAM", "MARCO", "CD14"],
    "lymphoid": ["CD3E", "CD79A", "NKG7"],
}
TISSUES = ["tLung", "nLung"]
MIN_CELLS = 50
GENES = sorted({g for v in list(GATES.values()) + list(IMMUNE_SPLIT.values()) for g in v}
               | set(LIGANDS) | set(RECEPTOR))

RULES = {
    "dataset": "GSE131907 (Kim et al. 2020), the human comparison the paper itself used",
    "cells": "deposited annotation, Sample_Origin in tLung or nLung",
    "compartment_primary": GATES,
    "gate_priority": ["immune", "endothelial", "epithelial", "mesenchymal", "unassigned"],
    "immune_split": IMMUNE_SPLIT,
    "compartment_crosscheck": "deposited Cell_type.refined, reported only",
    "genes": GENES,
    "measures": "detection fraction and mean deposited log2 TPM, per donor, tissue and compartment",
    "min_cells_per_donor_compartment": MIN_CELLS,
    "T1": "tumour lung, AREG detection higher in epithelial than myeloid cells, paired within donor, Wilcoxon signed-rank two-sided alpha 0.05",
    "T2": "tumour lung, which compartment has the highest median per-donor HBEGF detection",
    "T3": "the same table for EGFR",
    "refutation_rule": "T1 fails, or AREG detected in under 1% of tumour epithelial cells; either stops E2 to E4",
    "unit": "the donor; this is the first trial in this folder where a test is admissible",
    "caveats": ["marker gates are not cell-type calls", "detection depends on depth",
                "tumour epithelium is malignant tissue, not the earliest mutant state"],
}


def extract_rows(path, wanted):
    """Stream the gzipped matrix once and keep only the wanted gene rows."""
    rows = {}
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            name, sep, rest = line.partition("\t")
            if name not in wanted:
                continue
            values = np.array(rest.rstrip("\n").split("\t"), dtype=np.float32)
            rows[name] = values
            if len(rows) == len(wanted):
                break
    return np.array(header[1:]), rows


def get_matrix(rec):
    if CACHE.exists():
        blob = np.load(CACHE, allow_pickle=False)
        return blob["cells"], {g: blob[g] for g in blob.files if g != "cells"}
    rec.add_input(MATRIX)
    cells, rows = extract_rows(MATRIX, set(GENES))
    np.savez_compressed(CACHE, cells=cells, **rows)
    return cells, rows


def compartments(rows, n):
    def hit(genes):
        m = np.zeros(n, dtype=bool)
        for g in genes:
            if g in rows:
                m |= rows[g] > 0
        return m
    labels = np.full(n, "unassigned", dtype=object)
    for name in ("mesenchymal", "epithelial", "endothelial", "immune"):
        labels[hit(GATES[name])] = name
    immune = labels == "immune"
    labels[immune] = "other immune"
    for name in ("lymphoid", "myeloid", "neutrophil"):
        labels[hit(IMMUNE_SPLIT[name]) & immune] = name
    return labels


def main():
    from scipy.stats import wilcoxon
    rec = RunRecord(OUT / "e1_run_record.json", "E1 human LUAD EGFR ligand sources", RULES)
    rec.add_input(ANNOTATION)
    ann = pd.read_csv(ANNOTATION, sep="\t", compression="gzip").set_index("Index")
    cells, rows = get_matrix(rec)
    n = len(cells)
    frame = ann.reindex(cells)
    frame["compartment"] = compartments(rows, n)
    keep = frame["Sample_Origin"].isin(TISSUES).to_numpy()
    rec.set("cells_in_scope", int(keep.sum()))
    records = []
    groups = frame[keep].groupby(["Sample", "Sample_Origin", "compartment"]).indices
    kept = np.where(keep)[0]
    for (sample, tissue, comp), pos in groups.items():
        where = kept[pos]
        row = {"donor": sample, "tissue": tissue, "compartment": comp, "n_cells": len(where)}
        for gene in LIGANDS + RECEPTOR:
            if gene not in rows:
                continue
            vals = rows[gene][where]
            row["det_" + gene] = round(float((vals > 0).mean()), 4)
            row["mean_" + gene] = round(float(vals.mean()), 4)
        records.append(row)
    table = pd.DataFrame(records)
    table["evaluable"] = table["n_cells"] >= MIN_CELLS
    table.to_csv(OUT / "e1_per_donor_compartment.csv", index=False)
    rec.add_output(OUT / "e1_per_donor_compartment.csv")
    ok = table[table["evaluable"] & (table["tissue"] == "tLung")]
    wide = ok.pivot_table(index="donor", columns="compartment", values="det_AREG")
    paired = wide[["epithelial", "myeloid"]].dropna()
    diff = paired["epithelial"] - paired["myeloid"]
    stat = wilcoxon(paired["epithelial"], paired["myeloid"]) if len(paired) >= 5 else None
    p = float(stat.pvalue) if stat is not None else None
    epi_all = ok[ok["compartment"] == "epithelial"]["det_AREG"]
    t1_pass = bool(p is not None and p < 0.05 and diff.median() > 0)
    floor_pass = bool(len(epi_all) and epi_all.median() >= 0.01)
    rec.set("T1", {"n_donors": int(len(paired)), "median_epithelial": float(paired["epithelial"].median()) if len(paired) else None, "median_myeloid": float(paired["myeloid"].median()) if len(paired) else None, "p_value": p, "passes": t1_pass})
    refuted = not (t1_pass and floor_pass)
    rec.set("refutation_rule_triggered", refuted)
    rec.set("areg_floor_median_epithelial_detection", float(epi_all.median()) if len(epi_all) else None)
    medians = ok.groupby("compartment")[[c for c in ok.columns if c.startswith("det_")]].median().round(4)
    medians["n_donors"] = ok.groupby("compartment").size()
    medians.to_csv(OUT / "e1_tumour_median_by_compartment.csv")
    rec.add_output(OUT / "e1_tumour_median_by_compartment.csv")
    rec.set("T2_hbegf_top_compartment", str(medians["det_HBEGF"].idxmax()))
    rec.set("T3_egfr_top_compartment", str(medians["det_EGFR"].idxmax()))
    cross = pd.crosstab(frame.loc[keep, "compartment"], frame.loc[keep, "Cell_type.refined"])
    cross.to_csv(OUT / "e1_gate_vs_deposited_celltype.csv")
    rec.add_output(OUT / "e1_gate_vs_deposited_celltype.csv")
    lines = ["# Trial E1: EGFR ligand sources in human LUAD", "", "Refutation rule triggered: " + str(refuted), ""]
    lines += ["## Tumour lung, median per-donor detection", "", df_to_markdown(medians), ""]
    (OUT / "e1_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e1_summary.md")
    rec.finish()
    print("\n".join(lines))
    print(rec.record["results"])


if __name__ == "__main__":
    main()
