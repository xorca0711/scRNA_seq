#!/usr/bin/env python
"""Trial E1b: the same human question at the resolution the paper's claim lives at.

Trial E1 gated cells into broad compartments and found AREG epithelial but no
higher in tumour lung than in matched normal lung. That is a dilution artefact
of the gating rather than a contradiction: pooling all epithelium mixes
malignant states with normal AT2, club and ciliated cells, while the paper's
claim is about a specific transitional state inside the epithelium.

GSE131907 deposits a finer annotation that reaches that resolution, including
the tumour epithelial states tS1 to tS3 and, on the receiving side, the
COL13A1-positive matrix fibroblasts that correspond to alveolar fibroblasts.
This trial regroups by that annotation. It reads the gene vectors trial E1
already cached, so it fits nothing and costs seconds.

Provenance: rules fixed after E1's compartment table was seen and before any
subtype-level value was computed. Post hoc and disclosed.

Frozen rules:

* Cells: Sample_Origin in tLung or nLung, grouped by the deposited
  `Cell_subtype`. Subtypes named NA or Undetermined are dropped.
* Tumour epithelial state: tS1, tS2 and tS3 pooled, weighted by cell number.
  The comparison population is AT2.
* Measures and floors as in E1: detection fraction, at least 50 cells for a
  donor to contribute, at least 5 donors for a test.
* T5, the paper's claim at the right resolution. In tumour lung, AREG
  detection is higher in the pooled tumour epithelial states than in AT2
  cells, paired within donor, Wilcoxon signed-rank, two-sided.
* T6. In tumour lung, the subtypes ranked by median per-donor HBEGF detection.
* T7. EGFR detection in COL13A1-positive matrix fibroblasts against the
  tumour epithelial states, paired within donor, same test, to ask whether the
  receiving compartment of the mouse model is also the receiver here.
* The unit is the donor. Detection depends on depth, and the tumour states of
  this cohort are established adenocarcinoma rather than the earliest mutant
  state the mouse model captures.
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

OUT = HERE / "e1b_human_luad_by_subtype"
OUT.mkdir(exist_ok=True)
SERIES = RAW / "GSE131907"
CACHE = SERIES / "e1_extracted_rows.npz"
ANNOTATION = SERIES / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"

GENES = ["AREG", "EREG", "HBEGF", "TGFA", "EGFR"]
TUMOUR_STATES = ["tS1", "tS2", "tS3"]
AT2 = "AT2"
FIBROBLAST = "COL13A1+ matrix FBs"
DROP = {"NA", "Undetermined", "nan"}
MIN_CELLS = 50
MIN_DONORS = 5

RULES = {
    "source": "gene vectors cached by trial E1; deposited Cell_subtype as the grouping",
    "provenance_of_these_rules": "fixed after E1's compartment table was seen, before any subtype value was computed",
    "cells": "Sample_Origin in tLung or nLung; subtypes NA and Undetermined dropped",
    "tumour_epithelial_states": TUMOUR_STATES,
    "comparison_population": AT2,
    "receiving_population": FIBROBLAST,
    "genes": GENES,
    "floors": {"cells_per_donor_subtype": MIN_CELLS, "donors_for_a_test": MIN_DONORS},
    "T5": "tumour lung, AREG higher in pooled tumour epithelial states than in AT2, paired within donor",
    "T6": "tumour lung, subtypes ranked by median per-donor HBEGF detection",
    "T7": "EGFR in COL13A1+ matrix fibroblasts against the tumour epithelial states, paired within donor",
    "unit": "the donor",
    "caveats": ["detection depends on depth",
                "these tumour states are established adenocarcinoma, not the earliest mutant state"],
}


def pooled(table, gene, subtypes):
    sub = table[table["subtype"].isin(subtypes) & table["evaluable"]]
    if sub.empty:
        return pd.Series(dtype=float)
    hits = sub["det_" + gene] * sub["n_cells"]
    return hits.groupby(sub["donor"]).sum() / sub.groupby("donor")["n_cells"].sum()


def compare(a, b, label_a, label_b, gene, test):
    from scipy.stats import wilcoxon

    joined = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    out = {"test": test, "gene": gene, "group_a": label_a, "group_b": label_b,
           "n_donors": int(len(joined))}
    if len(joined):
        out["median_a"] = round(float(joined["a"].median()), 4)
        out["median_b"] = round(float(joined["b"].median()), 4)
        out["a_higher_in_donors"] = int((joined["a"] > joined["b"]).sum())
    if len(joined) >= MIN_DONORS:
        out["p_value"] = round(float(wilcoxon(joined["a"], joined["b"]).pvalue), 6)
        out["a_higher"] = bool(out["p_value"] < 0.05 and out["median_a"] > out["median_b"])
    else:
        out["p_value"] = None
        out["a_higher"] = None
    return out


def main():
    rec = RunRecord(OUT / "e1b_run_record.json",
                    "E1b human LUAD EGFR ligands by deposited subtype", RULES)
    for path in (CACHE, ANNOTATION):
        rec.add_input(path)
    blob = np.load(CACHE, allow_pickle=False)
    cells = blob["cells"]
    vectors = {g: blob[g] for g in GENES if g in blob.files}
    ann = pd.read_csv(ANNOTATION, sep="\t", compression="gzip").set_index("Index")
    frame = ann.reindex(cells)
    frame["subtype"] = frame["Cell_subtype"].astype(str)
    keep = (frame["Sample_Origin"].isin(["tLung", "nLung"]).to_numpy()
            & ~frame["subtype"].isin(DROP).to_numpy())
    kept = np.where(keep)[0]
    rec.set("cells_in_scope", int(len(kept)))

    groups = frame.iloc[kept].groupby(["Sample", "Sample_Origin", "subtype"]).indices
    rows = []
    for (donor, tissue, subtype), pos in groups.items():
        where = kept[pos]
        row = {"donor": donor, "tissue": tissue, "subtype": subtype, "n_cells": len(where)}
        for gene, vec in vectors.items():
            row["det_" + gene] = round(float((vec[where] > 0).mean()), 4)
        rows.append(row)
    table = pd.DataFrame(rows)
    table["evaluable"] = table["n_cells"] >= MIN_CELLS
    table.to_csv(OUT / "e1b_per_donor_subtype.csv", index=False)
    rec.add_output(OUT / "e1b_per_donor_subtype.csv")

    tumour = table[table["tissue"] == "tLung"]
    results = [
        compare(pooled(tumour, "AREG", TUMOUR_STATES), pooled(tumour, "AREG", [AT2]),
                "tumour states tS1-tS3", AT2, "AREG", "T5"),
        compare(pooled(tumour, "EGFR", [FIBROBLAST]), pooled(tumour, "EGFR", TUMOUR_STATES),
                FIBROBLAST, "tumour states tS1-tS3", "EGFR", "T7"),
        compare(pooled(tumour, "HBEGF", TUMOUR_STATES), pooled(tumour, "HBEGF", [AT2]),
                "tumour states tS1-tS3", AT2, "HBEGF", "T5 extra"),
    ]
    comp = pd.DataFrame(results)
    comp.to_csv(OUT / "e1b_paired_comparisons.csv", index=False)
    rec.add_output(OUT / "e1b_paired_comparisons.csv")
    rec.set("comparisons", results)

    ok = tumour[tumour["evaluable"]]
    rank = ok.groupby("subtype")[["det_" + g for g in vectors]].median().round(4)
    rank["n_donors"] = ok.groupby("subtype").size()
    rank = rank[rank["n_donors"] >= 3].sort_values("det_HBEGF", ascending=False)
    rank.to_csv(OUT / "e1b_tumour_ranking_by_subtype.csv")
    rec.add_output(OUT / "e1b_tumour_ranking_by_subtype.csv")
    rec.set("hbegf_top_subtypes", rank.index[:5].tolist())
    rec.set("areg_top_subtypes",
            rank.sort_values("det_AREG", ascending=False).index[:5].tolist())

    lines = ["# Trial E1b: human LUAD EGFR ligands at subtype resolution", "",
             "## Pre-registered paired comparisons, tumour lung", "",
             df_to_markdown(comp, index=False), "",
             "## Subtypes ranked by median per-donor HBEGF detection, tumour lung", "",
             df_to_markdown(rank), ""]
    (OUT / "e1b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e1b_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
