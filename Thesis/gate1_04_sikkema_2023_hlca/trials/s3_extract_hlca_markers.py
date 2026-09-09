#!/usr/bin/env python
"""Extract Supplementary Table 6 of Sikkema et al. 2023 (consensus markers) to a tracked CSV.

Input: the local, gitignored workbook 41591_2023_2327_MOESM3_ESM.xlsx, sheet
"6 - marker genes". Output: hlca_supp_table6_markers.csv with one row per
(cell_type, tier, gene). Tier 1 = compartment markers ("Full atlas"
reference), tier 2 = level-2 group markers, tier 3 = the type's own markers.
The paper is CC BY 4.0, so the marker lists can be redistributed with
attribution.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import openpyxl

SRC = Path(r"C:\Users\dream\Documents\AC_document\External Thesis\SAP_Thesis study\Gate_1-2_Universal"
           r"\Sikkema-2023-An-integrated-cell-atlas-of-the-lun\Dataset info\41591_2023_2327_MOESM3_ESM.xlsx")
OUT = Path(__file__).resolve().parent / "hlca_supp_table6_markers.csv"


def main() -> None:
    wb = openpyxl.load_workbook(SRC, read_only=True, data_only=True)
    ws = wb["6 - marker genes"]
    rows = list(ws.iter_rows(values_only=True))
    hdr = rows[1]
    types = [h.rsplit("_marker", 1)[0] for h in hdr if h and h.endswith("_marker")]
    out = []
    for t in types:
        i = hdr.index(t + "_marker")
        for r in rows[2:]:
            gene, mfor, ref = r[i], r[i + 1], r[i + 2]
            if not gene:
                continue
            # A marker is the type's own marker when marker_for names the type,
            # even if it was unique at the "Full atlas" level (Neuroendocrine).
            # Tier 1 is reserved for the four compartments.
            if str(mfor).startswith(t):
                tier = 3
            elif ref == "Full atlas" and mfor in ("Epithelial", "Immune", "Endothelial", "Stroma"):
                tier = 1
            else:
                tier = 2
            out.append({"cell_type": t, "tier": tier, "marker_for": str(mfor), "reference": str(ref), "gene": str(gene).strip(),
                        "note": "poss. lowly expressed, non-unique" if "poss." in str(mfor) else ""})
    with OUT.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=["cell_type", "tier", "marker_for", "reference", "gene", "note"])
        w.writeheader()
        w.writerows(out)
    print("types", len(types), "rows", len(out), "->", OUT.name)


if __name__ == "__main__":
    main()
