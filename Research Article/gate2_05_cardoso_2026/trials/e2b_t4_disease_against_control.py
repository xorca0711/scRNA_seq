#!/usr/bin/env python
"""Trial E2b: the T4 comparison trial E2 declared and did not compute.

Trial E2 froze four readings and produced three. T4, median per-donor
detection in disease against control for the same cell types, was written into
the rules and never written to a table. This script supplies it, and it does so
the way trial C2b supplied the missing C2 directions: by reading only the
tracked per-donor table E2 already wrote, so nothing is refitted and no rule is
touched after the fact.

T4 as frozen in E2: "median per-donor detection in IPF and control side by
side, no test". It stays descriptive. The donors are not paired, because a
donor is either a case or a control, so the only honest summary is two medians
and the number of donors behind each. A cell type is reported for a cohort only
if at least three donors clear E2's 50-cell floor on both sides, which keeps
disease-restricted states such as the aberrant basaloid population out of a
comparison they cannot enter.

The question T4 answers is the one E1 answered negatively in tumour lung: does
the ligand rise in the diseased lung, or is it simply present in both. Getting
the same answer in a second disease would say the ligand is a property of the
cell type rather than of the injury.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RunRecord, df_to_markdown  # noqa: E402

SOURCE = HERE / "e2_human_fibrosis_ligand_sources" / "e2_per_donor_celltype.csv"
OUT = HERE / "e2b_t4_disease_against_control"
OUT.mkdir(exist_ok=True)

GENES = ["AREG", "HBEGF", "EREG", "EGFR"]
MIN_DONORS = 3

RULES = {
    "source": "only the tracked per-donor table written by trial E2; nothing is refitted",
    "why_this_exists": "T4 was frozen in E2's rules and never computed; this supplies it without changing any rule",
    "T4": "median per-donor detection in IPF and control side by side, no test",
    "reporting_floor": f"at least {MIN_DONORS} donors on both sides of a cohort and cell type",
    "genes": GENES,
    "unit": "the donor",
    "caveats": ["donors are not paired, because a donor is either a case or a control",
                "disease-restricted states cannot enter this comparison at all",
                "detection depends on depth, and IPF and control libraries were not depth-matched here"],
}


def main():
    rec = RunRecord(OUT / "e2b_run_record.json", "E2b T4 disease against control", RULES)
    rec.add_input(SOURCE)
    table = pd.read_csv(SOURCE)
    if "evaluable" in table.columns:
        table = table[table["evaluable"]]

    rows = []
    for (cohort, celltype), group in table.groupby(["cohort", "celltype"]):
        ipf = group[group["disease"] == "IPF"]
        ctl = group[group["disease"] == "Control"]
        if len(ipf) < MIN_DONORS or len(ctl) < MIN_DONORS:
            continue
        row = {"cohort": cohort, "celltype": celltype,
               "n_donors_ipf": len(ipf), "n_donors_control": len(ctl)}
        for gene in GENES:
            column = "det_" + gene
            if column not in group.columns:
                continue
            a, b = float(ipf[column].median()), float(ctl[column].median())
            row["ipf_" + gene] = round(a, 4)
            row["control_" + gene] = round(b, 4)
            row["delta_" + gene] = round(a - b, 4)
        rows.append(row)

    frame = pd.DataFrame(rows)
    if frame.empty:
        raise SystemExit("no cohort and cell type clears the reporting floor on both sides")
    frame = frame.sort_values(["cohort", "delta_HBEGF"], ascending=[True, False])
    frame.to_csv(OUT / "e2b_t4_disease_against_control.csv", index=False)
    rec.add_output(OUT / "e2b_t4_disease_against_control.csv")

    counts = {}
    for gene in GENES:
        column = "delta_" + gene
        if column not in frame.columns:
            continue
        counts[gene] = {"cell_types_compared": int(frame[column].notna().sum()),
                        "higher_in_ipf": int((frame[column] > 0).sum()),
                        "median_delta": round(float(frame[column].median()), 4)}
    rec.set("direction_counts", counts)

    top = frame.nlargest(min(8, len(frame)), "delta_HBEGF")[
        ["cohort", "celltype", "ipf_HBEGF", "control_HBEGF", "delta_HBEGF",
         "n_donors_ipf", "n_donors_control"]]
    lines = ["# Trial E2b: T4, disease against control", "",
             "One row per cohort and cell type with at least three donors on both sides.",
             "Descriptive by the rule E2 froze; the donors are not paired.", "",
             "## Direction counts per gene", "",
             df_to_markdown(pd.DataFrame(counts).T), "",
             "## Largest HBEGF differences", "", df_to_markdown(top, index=False), ""]
    (OUT / "e2b_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e2b_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
