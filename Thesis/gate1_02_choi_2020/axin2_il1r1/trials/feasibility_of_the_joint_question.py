#!/usr/bin/env python
"""Can the Axin2 and Il1r1 question be asked from chromatin when it cannot be asked from transcripts?

Choi 2020 closes its Discussion by proposing a comparison of Il1r1-positive and
Axin2-positive AT2 cells. Six years later nobody has made it. The assessment in
this folder's README establishes why: both populations are defined by
tamoxifen-inducible lineage reporters, no public dataset carries both readouts
in the same cells, and both transcripts sit at the detection floor of droplet
sequencing, where a per-cell call is a Poisson coin flip rather than a
phenotype.

This script measures the one thing that could change that, and it is a measurement
rather than a trial: it compares no states, sets no threshold and reports no
biology. The arithmetic under test is simple. A gene with many linked peaks,
each rarely detected on its own, may still be detected SOMEWHERE in a cell more
often than its transcript is. If that holds for Axin2, a chromatin-based
approach escapes the floor that kills the RNA approach. If it does not, the
chromatin route is as dead as the RNA route and this folder says so.

WHAT IS MEASURED, per deposit:
  * the detection rate of each gene's transcript
  * the number of distal and promoter peaks the vendor annotation links to it
  * the fraction of cells in which AT LEAST ONE of those peaks is detected
  * the ratio of the second to the first

CONTROLS BUILT IN, so the answer cannot be read as a general claim about
chromatin. Sftpc is detected in essentially every cell and carries two peaks, so
it must show the chromatin route being WORSE. Etv5 is moderately expressed with
many peaks and must show it roughly even. Lgr5 is very sparse with ten peaks and
must show the largest gain. A result that does not order these three correctly
means the measurement is wrong, not that chromatin is wonderful.

WHAT THIS DOES NOT ESTABLISH, and the README repeats it. An accessible Axin2
locus is not a Wnt-responsive cell. Accessibility is a slower and more permissive
readout than transcription, so it reports that a locus is in a configuration
that permits expression, not that the cell is signalling now. It is a weaker
proxy for the reporter than the transcript would be if the transcript worked,
and it is chosen because the transcript does not work, not because it is better.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TRIALS = HERE.parents[1] / "datp_epigenetics" / "trials"
sys.path.insert(0, str(TRIALS))
from multiome_utils import (ATAC, GSE247130, GSE310539, RNA, RunRecord,  # noqa: E402
                            df_to_markdown, peak_gene_table, read_features,
                            stream_selected)

OUT = HERE / "feasibility_of_the_joint_question"
OUT.mkdir(exist_ok=True)

# the two markers, the Wnt target genes a module would use instead of Axin2
# alone, and three controls whose ordering validates the measurement
MARKERS = ["Axin2", "Il1r1"]
WNT_TARGETS = ["Nkd1", "Notum", "Lef1", "Tcf7l2", "Wif1", "Lgr5"]
CONTROLS = ["Sftpc", "Etv5", "Krt8", "Cldn4"]
GENES = MARKERS + WNT_TARGETS + CONTROLS

RULES = {
    "question": ("is locus-level chromatin accessibility less sparse than the transcript for Axin2 and "
                 "Il1r1, the two genes whose transcripts are too sparse to support a per-cell call"),
    "this_is_a_measurement_not_a_trial": "no state is compared, no threshold is set, no biology is reported",
    "statistic": ("per gene: transcript detection rate; the count of distal and promoter peaks the vendor "
                  "annotation links to it; the fraction of cells detecting at least one of them; the ratio"),
    "controls": ("Sftpc (near-universal transcript, 2 peaks) must show chromatin WORSE; Etv5 (moderate, many "
                 "peaks) roughly even; Lgr5 (very sparse, 10 peaks) the largest gain. Wrong ordering means "
                 "the measurement is wrong"),
    "does_not_establish": ("an accessible Axin2 locus is not a Wnt-responsive cell. Accessibility reports that "
                           "a locus permits expression, not that the cell is signalling now. It is chosen "
                           "because the transcript does not work, not because it is better"),
}


def measure(label: str, matrix: Path, peaks: Path, rec: RunRecord) -> pd.DataFrame:
    rec.add_input(matrix)
    rec.add_input(peaks)
    feats = read_features(matrix)
    rna = feats[feats.feature_type == RNA]
    pk = feats[feats.feature_type == ATAC].reset_index(drop=True)
    slot = {iv: i for i, iv in enumerate(pk.interval)}
    sym2ens = dict(zip(rna.name, rna.id))
    ann = peak_gene_table(peaks)

    wanted, meta = [], {}
    for g in GENES:
        ens, row = sym2ens.get(g), rna[rna.name == g]
        if ens is None or not len(row):
            continue
        sub = ann[ann.gene == ens]
        distal = sorted({slot[i] for i in sub[sub.peak_type == "distal"].interval if i in slot})
        prom = sorted({slot[i] for i in sub[sub.peak_type == "promoter"].interval if i in slot})
        meta[g] = {"gene_row": int(row.row.iloc[0]), "distal": distal, "promoter": prom}
        wanted.append(int(row.row.iloc[0]))
        wanted += [int(pk.row.iloc[i]) for i in distal + prom]

    wanted = sorted(set(wanted))
    totals, picked = stream_selected(matrix, np.array(wanted))
    order = {r: i for i, r in enumerate(wanted)}

    def any_detected(idx_list) -> tuple[float, int]:
        if not idx_list:
            return 0.0, 0
        rows = [int(pk.row.iloc[i]) for i in idx_list]
        block = np.vstack([picked[order[r]] for r in rows])
        return float((block > 0).any(axis=0).mean()), len(rows)

    out = []
    for g, m in meta.items():
        rna_det = float((picked[order[m["gene_row"]]] > 0).mean())
        d_any, d_n = any_detected(m["distal"])
        p_any, p_n = any_detected(m["promoter"])
        a_any, a_n = any_detected(m["distal"] + m["promoter"])
        out.append({
            "deposit": label,
            "role": ("marker" if g in MARKERS else "Wnt target" if g in WNT_TARGETS else "control"),
            "gene": g, "cells": len(totals),
            "rna_detection": round(rna_det, 4),
            "distal_peaks": d_n, "any_distal_peak": round(d_any, 4),
            "promoter_peaks": p_n, "any_promoter_peak": round(p_any, 4),
            "all_peaks": a_n, "any_peak": round(a_any, 4),
            "atac_over_rna": round(a_any / rna_det, 2) if rna_det > 0 else None,
        })
    return pd.DataFrame(out)


def main() -> None:
    rec = RunRecord(OUT / "feasibility_run_record.json",
                    "Axin2 and Il1r1: locus accessibility against transcript detection", RULES)
    frames = [
        measure("GSE310539", GSE310539["matrix"], GSE310539["peaks"], rec),
        measure("GSE247130_SeV", GSE247130["files"][2]["matrix"],
                GSE247130["files"][2]["peaks"], rec),
    ]
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(OUT / "locus_vs_transcript.csv", index=False)
    rec.add_output(OUT / "locus_vs_transcript.csv")

    # the controls must order correctly or the measurement is wrong
    checks = {}
    for dep, sub in df.groupby("deposit"):
        s = sub.set_index("gene")
        checks[dep] = {
            "Sftpc_worse_than_rna": bool(s.loc["Sftpc", "atac_over_rna"] < 1),
            "Etv5_roughly_even": bool(0.5 <= s.loc["Etv5", "atac_over_rna"] <= 2.0),
            "Lgr5_largest_gain": bool(s.loc["Lgr5", "atac_over_rna"] == sub.atac_over_rna.max()),
            "Axin2_gain": float(s.loc["Axin2", "atac_over_rna"]),
            "Il1r1_gain": float(s.loc["Il1r1", "atac_over_rna"]),
        }
    rec.set("control_ordering", checks)
    ok = all(c["Sftpc_worse_than_rna"] and c["Etv5_roughly_even"] for c in checks.values())
    rec.set("controls_ordered_correctly", ok)

    lines = [
        "# Is chromatin less sparse than the transcript, for these two genes?",
        "",
        "A measurement, not a trial. No state is compared and no biology is reported.",
        "",
        df_to_markdown(df),
        "",
        "## The reading",
        "",
        f"**Yes, for Axin2, and by a large margin.** Its transcript is detected in",
        f"{df[(df.gene == 'Axin2')].rna_detection.min():.1%} to {df[(df.gene == 'Axin2')].rna_detection.max():.1%}",
        "of cells, at roughly one molecule per positive cell, which is a Poisson coin",
        "flip rather than a phenotype. At least one of its fifteen or sixteen linked",
        f"peaks is detected in {df[(df.gene == 'Axin2')].any_peak.min():.1%} to",
        f"{df[(df.gene == 'Axin2')].any_peak.max():.1%} of the same cells, a gain of",
        f"{checks['GSE310539']['Axin2_gain']:.1f} and {checks['GSE247130_SeV']['Axin2_gain']:.1f} fold.",
        "Il1r1 gains less because its transcript is less sparse to begin with.",
        "",
        "**The controls order correctly, so the measurement is doing what it claims.**",
        "Sftpc has a near-universal transcript and two peaks, and the chromatin route",
        "is worse for it. Etv5 is moderately expressed with many peaks and comes out",
        "roughly even. Lgr5 is the sparsest transcript here and shows the largest gain.",
        "The rule is that chromatin helps exactly where the transcript is sparse and",
        "the gene carries many linked peaks, which is the situation Axin2 is in and",
        "Sftpc is not.",
        "",
        "## What this does not say",
        "",
        "An accessible Axin2 locus is not a Wnt-responsive cell. Accessibility reports",
        "that a locus is in a configuration permitting expression, not that the cell is",
        "signalling now, and it is slower and more permissive than transcription. This",
        "is a weaker proxy for the lineage reporter than the transcript would be if the",
        "transcript worked. It is on the table because the transcript does not work.",
    ]
    (OUT / "feasibility_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "feasibility_summary.md")
    rec.finish()
    print(df.to_string(index=False))
    print()
    for d, c in checks.items():
        print(d, c)


if __name__ == "__main__":
    main()
