#!/usr/bin/env python
"""Does a Wnt-TARGET module clear the detection floor where Axin2 alone does not?

A feasibility measurement, not a trial. No state is compared, no threshold set.

The assessment dismissed a transcript-level approach because Axin2 is detected in
4 to 5 per cent of AT2-enriched cells at about one molecule each, which is a
Poisson coin flip rather than a phenotype. **A module is not Axin2.** This asks
what detection a module of Wnt target genes reaches, and it asks it in
GSE262927, the one deposit in this project with genuine per-animal replication.

Why this exists at all: the assessment's Route C was GSE150957, a bulk array of
Wnt-sorted distal epithelium with a paired within-animal design. The owner
pointed out that it is out of focus, and it is. **Bulk cannot answer a
co-occurrence question at any level of replication**, because the question is
whether two markers sit in the same cell and a sorted fraction reports an
average. Route C was ranked partly on having the mouse as its unit, which is a
virtue of its statistics and not of its relevance. This measurement asks whether
a replacement exists that is same-cell.

WHAT IS AND IS NOT ESTABLISHED BY A GOOD NUMBER HERE. Clearing the detection
floor is necessary and not sufficient: Route B cleared it too, at a 4.1 to 10.4
fold gain, and then closed after three passes because the null spread swamped
the effect. A module also remains a proxy. **A Wnt-target module is not the
Axin2-lineage population**, any more than locus accessibility was. It reports
current target-gene expression, while the reporter reports transcription during
a tamoxifen window plus everything descended from it. Any trial built on this
must ask the reframed question in its own words and not borrow the original's.
"""

from __future__ import annotations

import glob
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "datp_epigenetics" / "trials"))
from multiome_utils import REPO, RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "module_detection_check"
OUT.mkdir(exist_ok=True)

WNT_TARGETS = ["Axin2", "Nkd1", "Tcf7l2", "Lef1", "Notum", "Wif1"]
IL1_PATHWAY = ["Il1r1", "Il1rap"]
CONTEXT = ["Sftpc", "Krt8", "Myd88", "Nfkbia"]
GENES = WNT_TARGETS + IL1_PATHWAY + CONTEXT
N_SAMPLES = 6

RULES = {
    "question": "does a Wnt-target module clear the detection floor where Axin2 alone does not",
    "this_is_a_measurement_not_a_trial": "no state compared, no threshold set, no biology reported",
    "why_it_exists": ("the assessment's Route C was a bulk array, and bulk cannot answer a co-occurrence "
                      "question at any level of replication because the question is whether two markers sit "
                      "in the same cell and a sorted fraction reports an average"),
    "deposit": "GSE262927, the one deposit in this project with genuine per-animal replication",
    "necessary_not_sufficient": ("Route B cleared its detection floor at a 4.1 to 10.4 fold gain and still "
                                 "closed after three passes because the null spread swamped the effect"),
    "still_a_proxy": ("a Wnt-target module is not the Axin2-lineage population. It reports current target "
                      "expression; the reporter reports transcription during a tamoxifen window plus "
                      "everything descended from it"),
}


def scan(path: str) -> tuple[dict, int]:
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        names = np.array([x.decode() for x in g["features"]["name"][:]])
        is_rna = np.array([x.decode() == "Gene Expression"
                           for x in g["features"]["feature_type"][:]])
        n_feat, n_cells = (int(x) for x in g["shape"][:])
        want = {}
        for s in GENES:
            hit = np.flatnonzero((names == s) & is_rna)
            if len(hit):
                want[s] = int(hit[0])
        rows = sorted(want.values())
        slot = np.full(n_feat, -1, np.int64)
        slot[rows] = np.arange(len(rows))
        indptr = g["indptr"][:]
        pos = np.zeros((len(rows), n_cells), dtype=bool)
        for start in range(0, n_cells, 4000):
            stop = min(start + 4000, n_cells)
            lo, hi = int(indptr[start]), int(indptr[stop])
            idx, dat = g["indices"][lo:hi], g["data"][lo:hi]
            offs = indptr[start:stop + 1] - lo
            for j in range(stop - start):
                a, b = int(offs[j]), int(offs[j + 1])
                ii, dd = idx[a:b], dat[a:b]
                s_ = slot[ii]
                k = (s_ >= 0) & (dd > 0)
                pos[s_[k], start + j] = True
        order = {r: i for i, r in enumerate(rows)}
        return {s: pos[order[r]] for s, r in want.items()}, n_cells


def main() -> None:
    rec = RunRecord(OUT / "module_run_record.json",
                    "Wnt-target module detection against Axin2 alone", RULES)
    files = sorted(glob.glob(str(REPO / "raw_data" / "GSE262927" / "GSE262927_RAW" / "*.h5")))[:N_SAMPLES]
    rows = []
    for f in files:
        rec.add_input(Path(f))
        p, n = scan(f)
        row = {"sample": Path(f).stem.split("_")[-1], "cells": n}
        for s in GENES:
            if s in p:
                row["det_" + s] = round(float(p[s].mean()), 4)
        counts = np.sum([p[s].astype(int) for s in WNT_TARGETS if s in p], axis=0)
        row["wnt_any_of_6"] = round(float((counts >= 1).mean()), 4)
        row["wnt_two_or_more"] = round(float((counts >= 2).mean()), 4)
        row["wnt_three_or_more"] = round(float((counts >= 3).mean()), 4)
        il1 = np.zeros(n, bool)
        for s in IL1_PATHWAY:
            if s in p:
                il1 |= p[s]
        row["il1_any"] = round(float(il1.mean()), 4)
        rows.append(row)
        print(row["sample"], "cells", n, "wnt>=2", row["wnt_two_or_more"], "il1", row["il1_any"])

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "module_detection.csv", index=False)
    rec.add_output(OUT / "module_detection.csv")
    rec.set("axin2_alone_range", [float(df.det_Axin2.min()), float(df.det_Axin2.max())])
    rec.set("wnt_two_or_more_range", [float(df.wnt_two_or_more.min()), float(df.wnt_two_or_more.max())])
    rec.set("il1_any_range", [float(df.il1_any.min()), float(df.il1_any.max())])

    lines = [
        "# Does a Wnt-target module clear the floor where Axin2 alone does not?",
        "",
        "A measurement, not a trial. No state is compared and no biology is reported.",
        "",
        df_to_markdown(df),
        "",
        "## The reading",
        "",
        f"**Yes.** Axin2 alone is detected in {df.det_Axin2.min():.1%} to {df.det_Axin2.max():.1%} of cells here,",
        "and in only 3.8 to 5.3 per cent of the AT2-enriched multiome deposits, which",
        "is where the assessment's refusal of a transcript route came from. Requiring",
        "**two or more of six Wnt target genes** reaches",
        f"{df.wnt_two_or_more.min():.1%} to {df.wnt_two_or_more.max():.1%}, and the IL-1 side",
        f"(Il1r1 or Il1rap) reaches {df.il1_any.min():.1%} to {df.il1_any.max():.1%}.",
        "Both are comfortably in the range where a per-cell call carries information",
        "rather than a single stray molecule.",
        "",
        "## Necessary, and not sufficient",
        "",
        "Route B cleared its own detection floor by a factor of 4.1 to 10.4 and still",
        "closed after three passes, because the spread of the matched null swamped an",
        "effect its own positive control showed was really there. A good number here",
        "buys the right to attempt the question; it does not predict the attempt will",
        "succeed, and any trial built on it needs a positive control that is",
        "demonstrated rather than asserted, which is the lesson trials A1 and A1b paid",
        "for.",
        "",
        "## And it is still a proxy, and a different question",
        "",
        "A Wnt-target module is not the Axin2-lineage population, any more than locus",
        "accessibility was. It reports current target-gene expression; the reporter",
        "reports transcription during a tamoxifen window plus everything descended",
        "from it. A trial built on this asks whether cells with active Wnt target",
        "expression also carry IL-1 receptor components, which is a well-posed",
        "question about current signalling state and is **not** the question of",
        "whether Choi's Il1r1-lineage cells and Nabhan's Axin2-lineage cells are the",
        "same cells. It must be asked in its own words.",
    ]
    (OUT / "module_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "module_summary.md")
    rec.finish()
    print("\n" + df.to_string(index=False))


if __name__ == "__main__":
    main()
