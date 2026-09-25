#!/usr/bin/env python
"""Trial E5: the one figure the E series earns, which is the refutation of C29.

Figure 2 of trial C5 shows what Areg deletion does to the six genes of the
paper's reprogrammed-fibroblast set: Runx1 and Pdgfrb keep almost all their
detection, Fst and Runx2 lose most of theirs. That figure states a measurement
and remains correct. What it invited was an interpretation, that a second
tumour-specific signal sustains the retained tier, and trial E4 refutes it.

The refutation is a two-panel argument and is hard to see in a table, which is
the only reason this figure exists. On the left, how much detection each gene
keeps after Areg deletion, from C2. On the right, what bleomycin alone does to
the same genes in sorted Col1a1-GFP mesenchyme carrying no oncogene, from E4.
A gene high on the left and high on the right is not evidence of a tumour
signal, because injury alone produces it.

Reads only the tracked tables of C5 and E4. Nothing is refitted. The figure
fails closed if the two sources stop agreeing about which genes are retained.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import viz_style as V  # noqa: E402
from cardoso_utils import RunRecord, df_to_markdown  # noqa: E402

RETENTION = HERE / "c5_figures_for_the_three_findings" / "c5_fibrotic_retention.csv"
INJURY = HERE / "e4_bleomycin_fibrotic_genes" / "e4_detection_by_library.csv"
OUT = HERE / "e5_figure_for_the_refutation"
OUT.mkdir(exist_ok=True)

GENES = ["Runx1", "Pdgfrb", "Acta2", "Tnc", "Fst", "Runx2"]
RULES = {
    "sources": "only the tracked tables of trials C5 and E4; nothing is refitted",
    "left_panel": "detection kept after Areg deletion, per cent of Areg-flox/+, from C2 via C5",
    "right_panel": "detection in bleomycin and untreated Col1a1-GFP mesenchyme, per animal, from E4",
    "argument": "a gene retained after Areg deletion that also rises with injury alone is not evidence of a tumour-specific signal",
    "fail_closed": "the figure is refused if the retained tier is no longer Runx1 and Pdgfrb",
    "unit": "one library per genotype on the left; one animal per library, two per group, on the right",
}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rec = RunRecord(OUT / "e5_run_record.json", "E5 figure for the C29 refutation", RULES)
    for path in (RETENTION, INJURY):
        rec.add_input(path)
    keep = pd.read_csv(RETENTION).set_index("gene")
    injury = pd.read_csv(INJURY)

    bleo = injury[injury["group"] == "bleomycin"]
    untr = injury[injury["group"] == "untreated"]
    rows = []
    for gene in GENES:
        column = "det_" + gene
        values_b = bleo[column].astype(float).tolist()
        values_u = untr[column].astype(float).tolist()
        rows.append({"gene": gene,
                     "retention_pct": round(float(keep.loc[gene, "retention_pct"]), 1),
                     "bleomycin_min": min(values_b), "bleomycin_max": max(values_b),
                     "untreated_min": min(values_u), "untreated_max": max(values_u),
                     "injury_generic": bool(min(values_b) > max(values_u))})
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "e5_retention_against_injury.csv", index=False)
    rec.add_output(OUT / "e5_retention_against_injury.csv")

    top_two = set(table.nlargest(2, "retention_pct")["gene"])
    if top_two != {"Runx1", "Pdgfrb"}:
        raise SystemExit("the retained tier is no longer Runx1 and Pdgfrb; refusing to draw this figure")
    both = table[(table["retention_pct"] >= 80) & table["injury_generic"]]["gene"].tolist()
    rec.set("retained_and_injury_generic", both)
    rec.set("injury_generic_count", int(table["injury_generic"].sum()))

    V.apply(plt)
    order = table.sort_values("retention_pct").reset_index(drop=True)
    fig = plt.figure(figsize=(13.6, 5.4), dpi=150)
    generic_n = int(table["injury_generic"].sum())
    fig.text(0.012, 0.945, "Injury alone turns on " + str(generic_n) + " of the 6 genes, whether or not "
             "they survive Areg deletion", fontsize=15, fontweight="semibold", color=V.INK)
    fig.text(0.012, 0.898, "Left: fibroblasts of GSE316244, one library per genotype. Right: sorted "
             "Col1a1-GFP mesenchyme of GSE132771, two bleomycin and two untreated animals, no oncogene "
             "present.", fontsize=10, color=V.INK_2)
    fig.text(0.012, 0.862, "The tiers on the left carry no tumour-specific information: Runx1 and "
             "Pdgfrb persist without Areg and also rise with bleomycin, so their persistence is "
             "fibroblast activation.", fontsize=10, color=V.INK_2)

    left = fig.add_axes([0.075, 0.12, 0.36, 0.68])
    left.set_xlim(0, 118)
    left.set_ylim(-0.6, len(order) - 0.4)
    left.set_yticks(range(len(order)), order["gene"])
    V.recessive(left)
    left.set_xlabel("detection kept after Areg deletion (% of Areg-flox/+)")
    left.set_title("Areg deletion, trials C2 and C5", fontsize=11)
    scale = V.px_scale(left)
    for i, row in order.iterrows():
        colour = V.FLOX_FLOX if row["retention_pct"] >= 80 else V.DEEMPH
        V.hbar(left, i, 0, float(row["retention_pct"]), 13, colour, scale)
        left.text(float(row["retention_pct"]) + 2.5, i, f"{row['retention_pct']:.0f}%",
                  fontsize=9, color=V.INK_2, va="center")
    left.axvline(80, color=V.AXIS, lw=1, ls=(0, (4, 3)), zorder=1)
    left.text(80.5, len(order) - 0.55, "C2 frozen rule: 80%", fontsize=8.5,
              color=V.MUTED, va="bottom")

    right = fig.add_axes([0.545, 0.12, 0.40, 0.68])
    top = float(max(order["bleomycin_max"].max(), order["untreated_max"].max())) * 1.22
    right.set_xlim(0, top)
    right.set_ylim(-0.6, len(order) - 0.4)
    right.set_yticks(range(len(order)), order["gene"])
    V.recessive(right)
    right.set_xlabel("fraction of mesenchymal cells with the gene detected")
    right.set_title("Bleomycin injury alone, trial E4", fontsize=11)
    for i, row in order.iterrows():
        for lo, hi, colour, name in ((row["untreated_min"], row["untreated_max"], V.MUTED, "untreated"),
                                     (row["bleomycin_min"], row["bleomycin_max"], V.SLOT[2], "bleomycin")):
            right.plot([lo, hi], [i, i], color=colour, lw=3.2, solid_capstyle="round", zorder=2)
            right.scatter([lo, hi], [i, i], s=26, color=colour, edgecolors=V.SURFACE,
                          linewidths=1.4, zorder=3,
                          label=name + " (2 animals)" if i == 0 else None)
        if row["injury_generic"]:
            right.text(float(row["bleomycin_max"]) + top * 0.025, i, "rises with injury",
                       fontsize=8.5, color=V.INK_2, va="center")
    right.legend(loc="lower right", handletextpad=0.4)

    path = OUT / "e5_retention_against_injury.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    rec.add_output(path)

    lines = ["# Trial E5: the figure for the C29 refutation", "",
             "Genes retained after Areg deletion at or above 80 per cent that also rise with "
             "bleomycin alone: " + (", ".join(both) or "none") + ".", "",
             df_to_markdown(table, index=False), ""]
    (OUT / "e5_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e5_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
