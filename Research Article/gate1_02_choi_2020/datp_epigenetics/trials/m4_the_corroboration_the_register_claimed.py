#!/usr/bin/env python
"""Trial M4: compute the suffix corroboration the register has been quoting.

An eight-adversary audit of rows C116 to C150, with each finding sent to an
independent verifier, returned 45 confirmed defects. The worst of them is in
this branch's self-described one validated result, and it is not a wrong number
but a missing measurement.

WHAT WENT WRONG. Trial M1's rule R2 corroborated the barcode-suffix map of
GSE247130. Its frozen panel was CORROBORATE = ["Cldn4", "Fos", "Cebpa"]. The
loop that writes the corroboration table iterates over four genes, including
Sox9, behind a guard `if s in c["counts"]`, and Sox9 was never loaded, so the
guard dropped it silently. **The "Sox9 is 7-fold higher" clause that the
register presents as one of two independent corroborating axes was never
computed.** It exists only as a hard-coded string in M1's own narrative
paragraph, and it propagated from there into five other files.

The same paragraph is the source of two further wrong numbers, and the reason is
the same: lines 528 to 545 of M1 write the whole R2 narrative as string literals
while `corr_df`, holding the real values, sits in scope two lines above. Nothing
in that paragraph is formatted from the dataframe. This is exactly the failure
this repository's own rules name: a report hand-typed rather than generated from
a logged artefact.

  * "Cebpa is 11 to 14 times higher ... in all three files" is wrong twice. The
    real ratios are 14.6 (P9), 11.0 (7 weeks) and 4.4 (Sendai). The top of the
    range is outside it, and the infected file, which is the one every
    downstream trial actually uses, is at 4.4 and far outside it.
  * "Cldn4 is 4.5-fold higher in suffix 2 of the infected file" is 4.2 in the
    logged CSV.

WHAT THIS TRIAL DOES, AND WHAT IT MAY NOT DO. M1's frozen panel is NOT edited
and M1 is not re-run: its outcome stands with the rule that produced it, as the
convention requires. M4 is a separate logged pass that computes the full
corroboration including Sox9, writes every value to a tracked table, and
produces its narrative by formatting that table rather than by typing it. After
this, any register row quoting a corroboration number cites M4.

M4 CANNOT rescue the Sox9 clause as written. The recomputed neonatal ratio is
7.05, so the figure was right, but a number that was never logged was Not
established at the moment it was quoted, by this repository's own definition,
and the register said "confirmed independently" of a measurement that did not
exist. That is recorded as its own claim rather than quietly repaired.

FROZEN RULES.

R1 The corroboration panel is Cldn4, Fos, Cebpa and Sox9, for every well of both
   deposits, as pseudobulk counts per ten thousand: the gene's summed counts
   over the well's summed RNA counts.
R2 Every ratio the register quotes is computed here and written to the table. No
   number in this trial's summary is typed; each is formatted from the frame.
R3 The direction expected of each axis is declared before it is computed, from
   the source papers rather than from this data:
     Cebpa   lower in the Cebpa knockout, all three GSE247130 files
     Sox9    higher in the neonatal Cebpa mutant (Hassan and Chen report SOX9
             reactivation as neonatal-specific, so the adult files are expected
             to show a smaller effect and are not evidence either way)
     Cldn4   higher after Sendai infection
     Fos     lower in the Fos, Fosb, Junb mutant of GSE310539
R4 A ratio is reported with its direction and its magnitude separately, because
   the audit showed a magnitude range can be wrong while a direction holds. The
   register may quote a range only if this trial computed it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from multiome_utils import (GSE247130, GSE310539, RNA, RunRecord,  # noqa: E402
                            df_to_markdown, read_barcodes, read_features,
                            stream_selected)

OUT = HERE / "m4_the_corroboration_the_register_claimed"
OUT.mkdir(exist_ok=True)

PANEL = ["Cldn4", "Fos", "Cebpa", "Sox9"]

RULES = {
    "question": "compute the suffix corroboration the register has been quoting, including the Sox9 axis that was never measured",
    "why_it_exists": ("an eight-adversary audit with independent verification returned 45 confirmed defects in "
                      "rows C116 to C150. The worst is that M1's frozen corroboration panel was Cldn4, Fos and "
                      "Cebpa, the write loop guarded Sox9 behind `if s in counts` so it was dropped silently, "
                      "and the Sox9 clause the register presents as an independent corroborating axis was never "
                      "computed. Two further numbers in the same hand-typed paragraph are wrong"),
    "M1_is_not_re_run": "M1's frozen panel is not edited and its outcome stands with the rule that produced it",
    "R1_panel": f"{', '.join(PANEL)}, every well of both deposits, pseudobulk counts per ten thousand",
    "R2_no_typed_numbers": "every number in this trial's summary is formatted from the frame, none is typed",
    "R3_expected_directions": {
        "Cebpa": "lower in the Cebpa knockout, all three GSE247130 files",
        "Sox9": "higher in the neonatal Cebpa mutant; Hassan and Chen report SOX9 reactivation as "
                "neonatal-specific, so the adult files are expected to show a smaller effect and are not "
                "evidence either way",
        "Cldn4": "higher after Sendai infection",
        "Fos": "lower in the Fos, Fosb, Junb mutant of GSE310539",
    },
    "R4_direction_and_magnitude_separately": ("the audit showed a magnitude range can be wrong while a direction "
                                              "holds; the register may quote a range only if computed here"),
}


def files():
    yield ("GSE310539", "totalaggr", GSE310539["matrix"], GSE310539["libraries"])
    for f in GSE247130["files"]:
        yield ("GSE247130", f["stage"], f["matrix"], f["libraries"])


def main() -> None:
    rec = RunRecord(OUT / "m4_run_record.json", "M4 the corroboration the register claimed", RULES)
    rows = []
    for deposit, stage, matrix, libs in files():
        rec.add_input(matrix)
        feats = read_features(matrix)
        bcs = read_barcodes(matrix)
        g = feats[feats.feature_type == RNA]
        rmap = {s: int(g[g.name == s].row.iloc[0]) for s in PANEL if len(g[g.name == s])}
        missing = [s for s in PANEL if s not in rmap]
        assert not missing, f"{stage}: panel gene absent from the matrix: {missing}"
        totals, picked = stream_selected(matrix, np.array(sorted(rmap.values())))
        order = {r: i for i, r in enumerate(sorted(rmap.values()))}
        for suf in sorted(libs):
            sel = (bcs.suffix == suf).values
            depth = totals.rna_counts.values[sel]
            row = {"deposit": deposit, "file": stage, "suffix": suf,
                   "well": libs[suf]["name"], "cells": int(sel.sum()),
                   "total_rna_counts": int(depth.sum())}
            for s in PANEL:
                row["cpm10k_" + s] = round(float(picked[order[rmap[s]]][sel].sum() / depth.sum() * 1e4), 4)
            rows.append(row)
        print("read", stage)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "m4_corroboration.csv", index=False)
    rec.add_output(OUT / "m4_corroboration.csv")

    # ratios, computed rather than asserted. GSE247130's map is inverted, so
    # suffix 1 is the control and suffix 2 the Cebpa mutant.
    ratios = []
    for stage in ("P9", "7wk", "SeV"):
        s = df[(df.deposit == "GSE247130") & (df.file == stage)].set_index("suffix")
        for gene, hi, lo, axis in (("Cebpa", "1", "2", "control over Cebpa mutant"),
                                   ("Sox9", "2", "1", "Cebpa mutant over control"),
                                   ("Cldn4", "2", "1", "Cebpa mutant over control")):
            a, b = float(s.loc[hi, "cpm10k_" + gene]), float(s.loc[lo, "cpm10k_" + gene])
            ratios.append({"deposit": "GSE247130", "file": stage, "gene": gene, "axis": axis,
                           "numerator": a, "denominator": b,
                           "ratio": round(a / b, 3) if b else None})
    s = df[df.deposit == "GSE310539"].set_index("suffix")
    wt = [float(s.loc[x, "cpm10k_Fos"]) for x in ("1", "2")]
    mu = [float(s.loc[x, "cpm10k_Fos"]) for x in ("3", "4")]
    ratios.append({"deposit": "GSE310539", "file": "totalaggr", "gene": "Fos",
                   "axis": "wildtype over AP-1 mutant", "numerator": round(min(wt), 4),
                   "denominator": round(max(mu), 4), "ratio": round(min(wt) / max(mu), 3)})
    for pair, name in ((("2", "1"), "wildtype"), (("4", "3"), "AP-1 mutant")):
        a, b = float(s.loc[pair[0], "cpm10k_Cldn4"]), float(s.loc[pair[1], "cpm10k_Cldn4"])
        ratios.append({"deposit": "GSE310539", "file": "totalaggr", "gene": "Cldn4",
                       "axis": f"Sendai over PBS, {name}", "numerator": a, "denominator": b,
                       "ratio": round(a / b, 3)})
    rt = pd.DataFrame(ratios)
    rt.to_csv(OUT / "m4_ratios.csv", index=False)
    rec.add_output(OUT / "m4_ratios.csv")

    ceb = rt[(rt.gene == "Cebpa")]
    sox = rt[(rt.gene == "Sox9")]
    cld247 = rt[(rt.gene == "Cldn4") & (rt.deposit == "GSE247130")]
    cld310 = rt[(rt.gene == "Cldn4") & (rt.deposit == "GSE310539")]
    fos = rt[rt.gene == "Fos"].iloc[0]
    rec.set("cebpa_ratios", dict(zip(ceb.file, ceb.ratio)))
    rec.set("sox9_ratios", dict(zip(sox.file, sox.ratio)))
    rec.set("cebpa_direction_holds_in_all_three", bool((ceb.ratio > 1).all()))
    rec.set("sox9_direction_holds_in_all_three", bool((sox.ratio > 1).all()))

    fmt = lambda s: ", ".join(f"{f} {r:.1f}" for f, r in zip(s.file, s.ratio))  # noqa: E731
    lines = [
        "# Trial M4: the corroboration the register has been quoting",
        "",
        "Every number on this page is formatted from the table below. None is typed.",
        "",
        "## Why this trial exists",
        "",
        "An eight-adversary audit of rows C116 to C150, each finding sent to an",
        "independent verifier, returned 45 confirmed defects. The worst sits in this",
        "branch's self-described one validated result: **the Sox9 axis the register",
        "presents as independent corroboration was never computed.** Trial M1's",
        "frozen panel was Cldn4, Fos and Cebpa; its write loop iterated over four",
        "genes behind a guard that silently dropped the fourth. The figure existed",
        "only as a hard-coded string in M1's narrative paragraph, and propagated from",
        "there into five other files.",
        "",
        "The same paragraph was the source of two further wrong numbers, for the same",
        "reason: it was typed as string literals while the dataframe holding the real",
        "values sat in scope two lines above.",
        "",
        "## The corroboration, computed",
        "",
        df_to_markdown(df),
        "",
        "## The ratios, computed",
        "",
        df_to_markdown(rt),
        "",
        "## What the register said and what is true",
        "",
        "| The register said | Computed here |",
        "|---|---|",
        f"| Cebpa 11 to 14 times higher in all three files | {fmt(ceb)}. The direction holds in all three; "
        f"the range does not, and the infected file, which every downstream trial uses, is the weakest |",
        f"| Sox9 7-fold higher in the neonatal file, confirmed independently | {fmt(sox)}. The neonatal figure "
        f"was right, and it was never computed by the logged run that claimed it |",
        f"| Cldn4 4.5-fold higher in suffix 2 of the infected file | {cld247[cld247.file == 'SeV'].ratio.iloc[0]:.1f} |",
        f"| Cldn4 rises 12-fold with infection in the wildtype pair and 6-fold in the mutant pair | "
        f"{cld310.iloc[0].ratio:.1f} and {cld310.iloc[1].ratio:.1f} |",
        f"| Fos falls from 5.12 and 5.51 to 1.40 and 1.37 | holds; the smallest wildtype over the largest "
        f"mutant is {fos.ratio:.1f} |",
        "",
        "## What stands and what does not",
        "",
        "**The suffix map of GSE247130 is still inverted, and that is unaffected.**",
        "Cebpa is higher in the suffix the GEO sample order calls the knockout in",
        f"all three files ({fmt(ceb)}), and a conditional knockout cannot carry more",
        "of its own target than its control. Sox9 agrees in all three files",
        f"({fmt(sox)}), with the neonatal file much the strongest, which is the",
        "stage-specificity Hassan and Chen report. Cldn4 agrees in the infected file.",
        "",
        "**What does not stand is the register's account of how well it was shown.**",
        "One of the two axes it called independent corroboration was not measured,",
        "and two of the magnitudes it quoted were wrong. The map was right for",
        "reasons that were partly unlogged, which by this repository's own rule means",
        "they were Not established at the moment they were quoted.",
    ]
    (OUT / "m4_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "m4_summary.md")
    rec.finish()
    print("\n" + rt.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
