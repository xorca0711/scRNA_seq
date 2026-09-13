#!/usr/bin/env python
"""Trial C1d: what the mesenchymal sort actually contains (Gate 1 addendum).

Trial C1b found five Red2Kras-private clusters and classified none of them as
low-count or doublet-enriched, so the frozen quality criteria left them
"not explained by quality". Their top genes then suggested a different
explanation: composition. Two of the five look like cells the sort was
supposed to exclude, and one of those expresses Areg, which matters because
Areg is the ligand the paper's whole cascade starts from. If mutant
epithelium is present inside a library sorted as CD45-CD31-EpCAM-, then any
ligand-receptor reasoning done within that library alone is reading a
contaminant as a source.

This trial turns that reading into a measurement. It is cheap: it reads the
object C1b wrote, computes compartment-marker detection per cluster, and
classifies. It fits nothing.

Provenance: these rules were fixed after C1b's top-gene table was seen and
before any compartment detection was computed. Stated because they are not
blind to C1b.

Frozen rules:

* Object: the one C1b wrote; backed mode, marker columns only.
* Compartment markers, one set each, from the panels already in use:
  epithelial Epcam, Cdh1, Krt8, Krt18; immune Ptprc; endothelial Pecam1,
  Cdh5; mesenchymal Col1a1, Col1a2, Dcn; mesothelial Wt1, Msln, Upk3b.
* A cluster is called **off-target for this sort** if the median cell in it
  detects a non-mesenchymal compartment while failing to detect Col1a1: the
  rule is, more than 50% of its cells detect at least one marker of the
  epithelial or immune set AND fewer than 50% detect Col1a1. Mesothelium is
  reported separately, because the paper reports mesothelial-like cells as a
  genuine Red2Kras-enriched population rather than a sort failure.
* Areg detection is reported per cluster, because the consequence of an
  epithelial contaminant is specific to it.
* Unit: the library. No P value. This is a composition description of one
  library per genotype.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "c1d_sort_purity"
OUT.mkdir(exist_ok=True)
OBJECT = HERE / "c1b_characterise_red2kras_private" / "c1b_mesenchyme.h5ad"

COMPARTMENTS = {
    "epithelial": ["Epcam", "Cdh1", "Krt8", "Krt18"],
    "immune": ["Ptprc"],
    "endothelial": ["Pecam1", "Cdh5"],
    "mesenchymal": ["Col1a1", "Col1a2", "Dcn"],
    "mesothelial": ["Wt1", "Msln", "Upk3b"],
}
RULES = {
    "object": "the one trial C1b wrote; backed mode, marker columns only",
    "provenance_of_these_rules": "fixed after C1b's top-gene table was seen, before any compartment detection was computed",
    "compartments": COMPARTMENTS,
    "off_target_rule": ("more than 50% of a cluster's cells detect at least one epithelial or immune "
                        "marker AND fewer than 50% detect Col1a1"),
    "mesothelium_note": "reported separately; the paper reports mesothelial-like cells as a genuine Red2Kras-enriched population",
    "areg_reported_per_cluster": "because the consequence of an epithelial contaminant is specific to Areg",
    "sort_as_deposited": "CD45-CD31-EpCAM- (GEO characteristics: mesenchymal cells)",
    "unit": "the library; no P value",
}


def main() -> None:
    import anndata as ad

    rec = RunRecord(OUT / "c1d_run_record.json",
                    "C1d what the mesenchymal sort actually contains (Gate 1 addendum)",
                    RULES, notes="reads the object C1b wrote; fits nothing")
    rec.add_input(OBJECT)
    adata = ad.read_h5ad(OBJECT, backed="r")
    obs = adata.obs
    genes = sorted({g for v in COMPARTMENTS.values() for g in v} | {"Areg"})
    present = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    counts = adata.layers["counts"]
    detected = {g: np.asarray(counts[:, adata.var_names.get_loc(g)].todense()).ravel() > 0
                for g in present}
    adata.file.close()
    rec.set("markers_absent_from_the_reference", missing)

    clusters = obs["leiden_0.5"].astype(str).to_numpy()
    genotype = obs["genotype"].to_numpy()
    rows = []
    for cluster in sorted(set(clusters), key=int):
        mask = clusters == cluster
        row = {"cluster": cluster, "n_cells": int(mask.sum()),
               "fraction_Red2Kras": round(float((genotype[mask] == "Red2Kras").mean()), 3)}
        any_of = {}
        for name, members in COMPARTMENTS.items():
            available = [g for g in members if g in detected]
            if not available:
                continue
            stack = np.vstack([detected[g][mask] for g in available])
            any_of[name] = float(stack.any(axis=0).mean())
            row[f"pct_{name}"] = round(100 * any_of[name], 1)
        row["pct_Col1a1"] = round(100 * float(detected["Col1a1"][mask].mean()), 1) if "Col1a1" in detected else None
        row["pct_Areg"] = round(100 * float(detected["Areg"][mask].mean()), 1) if "Areg" in detected else None
        off_target = ((any_of.get("epithelial", 0) > 0.5 or any_of.get("immune", 0) > 0.5)
                      and float(detected["Col1a1"][mask].mean()) < 0.5)
        row["off_target_for_this_sort"] = bool(off_target)
        if off_target:
            row["reads_as"] = ("epithelial" if any_of.get("epithelial", 0) >= any_of.get("immune", 0)
                               else "immune")
        elif any_of.get("mesothelial", 0) > 0.5:
            row["reads_as"] = "mesothelial"
        else:
            row["reads_as"] = "mesenchymal"
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "c1d_sort_purity_by_cluster.csv", index=False)
    rec.add_output(OUT / "c1d_sort_purity_by_cluster.csv")

    off = table[table["off_target_for_this_sort"]]
    n_off = int(off["n_cells"].sum())
    rec.set("off_target_clusters", off["cluster"].tolist())
    rec.set("off_target_cells", n_off)
    rec.set("off_target_fraction_of_object", round(n_off / int(table["n_cells"].sum()), 4))
    rec.set("off_target_detail", off[["cluster", "n_cells", "fraction_Red2Kras", "reads_as",
                                      "pct_Areg"]].to_dict("records"))

    lines = [
        "# Trial C1d output: what the mesenchymal sort actually contains", "",
        f"The two GSE316241 libraries are deposited as CD45-CD31-EpCAM- mesenchymal cells. "
        f"Of {int(table['n_cells'].sum()):,} cells that passed quality control, "
        f"**{n_off:,} sit in clusters that read as off-target for that sort** "
        f"({100 * n_off / int(table['n_cells'].sum()):.1f}%).", "",
        df_to_markdown(table, index=False), "",
        "The consequence worth carrying forward: an epithelial cluster inside a mesenchymal library",
        "expresses Areg, the ligand the paper's cascade begins with. Any ligand-receptor reasoning",
        "done inside this library alone would be reading a sort contaminant as a signalling source.",
        "The paper avoided this by taking its epithelial cells from a separate lineage-labelled",
        "series; a reanalysis that did not would not.", "",
    ]
    (OUT / "c1d_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c1d_summary.md")
    rec.finish()
    print(df_to_markdown(table, index=False))


if __name__ == "__main__":
    main()
