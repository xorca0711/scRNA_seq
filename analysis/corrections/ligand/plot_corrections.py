"""Figures from corrected evidence only; does not launch analyses."""
from __future__ import annotations

import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import HERE, ROOT, RESOURCES, TARGETS

RESULTS = HERE / "results"
FIG = RESULTS / "figures"
OLD = ROOT / "Research Article" / "gate2_05_cardoso_2026" / "trials"
COLORS = {"historical": "#858d99", "corrected": "#087e8b"}


def save(fig, name):
    fig.savefig(FIG / (name + ".png"), dpi=170, bbox_inches="tight")
    fig.savefig(FIG / (name + ".svg"), bbox_inches="tight")
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    old = pd.read_csv(OLD / "c14_does_the_ranking_depend_on_the_database" / "c14_resource_summary.csv").set_index("resource")
    new = pd.read_csv(RESULTS / "lr" / "resource_summary.csv").set_index("resource")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
    history = pd.read_csv(OLD / "c12_cellchatdb_full_resource_scan" / "c12_per_donor_pairs.csv")
    now = pd.read_csv(RESULTS / "lr" / "cellchatdb_donor_pair_winners.csv")
    for x, (label, table) in enumerate([("Historical\n26 donors", history), ("Corrected\n22 donors", now)]):
        head = table.loc[table["rank"] <= 15]
        mural = head.source.isin(["Pericyte", "SMC"]).mean()
        axes[0].bar(x, mural * 100, color=list(COLORS.values())[x], width=.6)
        axes[0].text(x, mural * 100 + 1.5, f"{int(round(mural*len(head)))}/{len(head)}", ha="center")
    axes[0].set(xticks=[0, 1], xticklabels=["Historical\n26 donors", "Corrected\n22 donors"],
                ylim=(0, 38), ylabel="Mural senders in donor top 15 (%)",
                title="A  Sender provenance, CellChatDB")
    positions = np.arange(len(RESOURCES))
    axes[1].barh(positions - .17, old.loc[RESOURCES].abundance_share_top15 * 100,
                 height=.32, color=COLORS["historical"], label="Historical")
    axes[1].barh(positions + .17, new.loc[RESOURCES].abundance_share_top15 * 100,
                 height=.32, color=COLORS["corrected"], label="Corrected")
    axes[1].set(yticks=positions, yticklabels=RESOURCES, xlim=(0, 100),
                xlabel="Historical abundance flag in top 15 (%)", title="B  Flag composition, not signalling")
    axes[1].invert_yaxis()
    axes[1].legend(frameon=False, fontsize=9, loc="lower right")
    hist_ranks = pd.read_csv(OLD / "c14_does_the_ranking_depend_on_the_database" / "c14_rankings_by_resource.csv")
    for i, resource in enumerate(RESOURCES):
        h = hist_ranks.loc[hist_ranks.resource.eq(resource) & hist_ranks.pair.eq("AREG to EGFR")]
        if h.empty or pd.isna(new.loc[resource, "areg_median_rank"]):
            axes[2].text(.5, i, "No exact AREG–EGFR pair", va="center", ha="center", fontsize=9,
                         transform=axes[2].get_yaxis_transform())
            continue
        a, b = h.iloc[0].median_rank, new.loc[resource, "areg_median_rank"]
        axes[2].plot([a, b], [i, i], color="#b8c0c8", linewidth=2)
        axes[2].scatter(a, i, color=COLORS["historical"], s=45)
        axes[2].scatter(b, i, color=COLORS["corrected"], s=45, marker="D")
    axes[2].set(yticks=positions, yticklabels=RESOURCES, xlabel="AREG–EGFR median within-donor rank",
                title="C  Rank in each resource's pair universe")
    axes[2].invert_yaxis()
    fig.suptitle("Corrected ligand scan: epithelial senders and actual fibroblast cell floors", fontsize=14)
    save(fig, "lr_scope_correction")

    donors = pd.read_csv(RESULTS / "c37" / "per_donor.csv")
    summary = json.loads((RESULTS / "c37" / "summary.json").read_text())
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
    fields = ["original_detection", "detection_at_500", "detection_at_1000", "detection_at_2000"]
    wide = {f: donors.pivot(index="donor", columns="compartment", values=f).dropna() for f in fields}
    all_donors = sorted(set().union(*(set(w.index) for w in wide.values())))
    for donor in all_donors:
        ys = [float(w.loc[donor, "epithelial"] - w.loc[donor, "myeloid"]) if donor in w.index else np.nan for w in wide.values()]
        axes[0].plot(range(4), ys, color="#b3c5cf", linewidth=1, marker="o", markersize=3)
    medians = [summary[f]["median_paired_difference"] for f in fields]
    axes[0].plot(range(4), medians, color="#087e8b", linewidth=3, marker="D", label="Median donor difference")
    axes[0].axhline(0, linestyle="--", color="#555", linewidth=.8)
    axes[0].set(xticks=range(4), xticklabels=[f"Original\nn={summary[fields[0]]['n_donors']}"] +
        [f"{b} UMIs\nn={summary['detection_at_'+str(b)]['n_donors']}" for b in [500, 1000, 2000]],
        ylabel="Epithelial − myeloid AREG detection", title="A  Donor-paired annotation/depth sensitivity")
    axes[0].legend(frameon=False, fontsize=9)
    w = wide["detection_at_1000"]
    axes[1].scatter(w.myeloid, w.epithelial, color="#087e8b", s=38)
    limit = max(w.epithelial.max(), w.myeloid.max()) * 1.1
    axes[1].plot([0, limit], [0, limit], "--", color="#555", linewidth=.8)
    axes[1].set(xlim=(0, limit), ylim=(0, limit), xlabel="Myeloid expected detection at 1000 UMIs",
                ylabel="Epithelial expected detection at 1000 UMIs", title="B  Each point is one donor")
    axes[1].set_aspect("equal", adjustable="box")
    fig.suptitle("C37 post hoc sensitivity: deposited annotation and molecule-matched detection", fontsize=13)
    save(fig, "c37_annotation_depth")


if __name__ == "__main__":
    main()
