#!/usr/bin/env python
"""Trial C5: figures for the three findings, and the corrections drawing them forced.

The session summary named three findings the paper does not report. None of
the earlier trial figures shows any of them; they existed only as tables.
This trial draws them from the tracked tables of C1d, C2 and C3 and from the
object C1b wrote, and fits nothing new.

Drawing them exposed places where the earlier wording overstated the tables.
The corrections are outputs of this trial, not quiet edits elsewhere:

1. Hbegf ranks second by ABUNDANCE inside the DATP-like state, not by
   ENRICHMENT over AT2 cells of the same library, where it is broadly
   expressed too.
2. The two 4-day mutant replicates disagree about 24-fold on the DATP-like
   share (1.2% and 28.2%); the 14.7% median C3 reported describes neither.
3. The fibrotic "split" of C2 is graded, not binary; the frozen 80% line makes
   the binary.

A fourth table is a discriminating check, not a conclusion: for the Areg-high
epithelial cluster inside the mesenchymal sort, it reports the markers that
separate three explanations (mutant cells with low Epcam escaping an EpCAM-
gate, epithelial-fibroblast doublets, or ambient RNA).

Provenance: fixed after the C1d, C2 and C3 tables had been seen. Post hoc and
disclosed.

Frozen rules:

* Enrichment of a ligand = mean log1p(CP10K) in DATP-like cells minus mean in
  AT2 cells of the same library (the log-fold-like quantity). The ratio of
  detection fractions is tabulated as a second notion and not plotted.
* Only the four Red2Kras RFP libraries are plotted for ligands. The wild-type
  YFP DATP-like cells (21 and 36) are tabulated with a flag.
* Retention = detection in Areg-flox/flox fibroblasts over Areg-flox/+, the
  quantity C2's 80% rule used. The Confetti column is not plotted
  (other series, shallower, depth-confounded).
* Cluster 11 check: detection of Epcam, Cdh1, Krt8, Krt18, AT2, DATP-like and
  AT1 markers, Col1a1 and Pdgfra, and co-detection of Krt8 with Col1a1, in
  clusters 11, 12 and 14 and the whole object.
* Palette validated with the dataviz validator on the light surface (see
  viz_style.py). Unit: the library. No P value.
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

OUT = HERE / "c5_figures_for_the_three_findings"
OUT.mkdir(exist_ok=True)
C1D = HERE / "c1d_sort_purity" / "c1d_sort_purity_by_cluster.csv"
C2 = HERE / "c2_areg_deletion_arm" / "c2_fibrotic_programme_survival.csv"
C3_LIG = HERE / "c3_areg_state_specificity" / "c3_ligand_expression.csv"
C3_COMP = HERE / "c3_areg_state_specificity" / "c3_state_composition.csv"
OBJECT = HERE / "c1b_characterise_red2kras_private" / "c1b_mesenchyme.h5ad"

LIGANDS = ["Areg", "Hbegf", "Ereg", "Tgfa"]
LIGAND_COLOUR = {g: V.SLOT[i + 1] for i, g in enumerate(LIGANDS)}
FIBROTIC = ["Runx1", "Pdgfrb", "Acta2", "Tnc", "Fst", "Runx2"]
MUTANT = [("Expt1_4dRFPr1", "4 days, rep 1"), ("Expt1_4dRFPr2", "4 days, rep 2"),
          ("Expt1_2wRFPr1", "2 weeks, rep 1"), ("Expt1_2wRFPr2", "2 weeks, rep 2")]
CLUSTER11_GENES = ["Epcam", "Cdh1", "Krt8", "Krt18", "Sftpc", "Lamp3", "Cldn4", "Sox9", "Itga2",
                   "Ager", "Hopx", "Areg", "Col1a1", "Pdgfra"]

RULES = {
    "provenance_of_these_rules": "fixed after the C1d, C2 and C3 tables had been seen; post hoc and disclosed",
    "enrichment": "mean log1p(CP10K) DATP-like minus AT2, same library; detection ratio tabulated, not plotted",
    "ligand_libraries_plotted": [lib for lib, _ in MUTANT],
    "retention": "Areg-flox/flox over Areg-flox/+ detection; Confetti column not plotted (depth-confounded)",
    "cluster11_check_genes": CLUSTER11_GENES,
    "palette": {"surface": V.SURFACE, "grouped_bars_adjacent": list(V.SLOT.values()),
                "scatter_all_pairs": [V.SLOT[1], V.SLOT[2]], "dumbbell_all_pairs": [V.FLOX_PLUS, V.FLOX_FLOX],
                "validation": "all hard checks pass (light); contrast WARN on aqua, yellow and the lighter blue, relieved by the tables"},
    "unit": "the library; no P value",
}


# ---- derived tables ---------------------------------------------------------
def ligand_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    lig = pd.read_csv(C3_LIG)
    rows = []
    for library in lig["library"].unique():
        sub = lig[lig["library"] == library].set_index("state")
        if not {"DATP_like", "AT2"}.issubset(sub.index):
            continue
        d, a = sub.loc["DATP_like"], sub.loc["AT2"]
        for gene in LIGANDS:
            rows.append({
                "library": library, "clone": d["clone"], "timepoint": d["timepoint"],
                "n_DATP_like": int(d["n_cells"]), "n_AT2": int(a["n_cells"]), "ligand": gene,
                "mean_DATP_like": round(float(d[f"mean_{gene}"]), 4),
                "mean_AT2": round(float(a[f"mean_{gene}"]), 4),
                "enrichment_log1p_diff": round(float(d[f"mean_{gene}"] - a[f"mean_{gene}"]), 4),
                "detection_ratio": (round(float(d[f"det_{gene}"] / a[f"det_{gene}"]), 2)
                                    if a[f"det_{gene}"] > 0 else None),
                "plotted": library in dict(MUTANT),
                "flag": "" if library in dict(MUTANT) else "wild-type clone, fewer than 50 DATP-like cells",
            })
    table = pd.DataFrame(rows)
    ranks = []
    for library, sub in table.groupby("library", sort=False):
        by_abund = sub.sort_values("mean_DATP_like", ascending=False)["ligand"].tolist()
        by_enrich = sub.sort_values("enrichment_log1p_diff", ascending=False)["ligand"].tolist()
        by_ratio = sub.dropna(subset=["detection_ratio"]).sort_values("detection_ratio", ascending=False)["ligand"].tolist()
        ranks.append({"library": library, "n_DATP_like": int(sub["n_DATP_like"].iloc[0]),
                      "rank_by_abundance": " > ".join(by_abund),
                      "rank_by_enrichment_over_AT2": " > ".join(by_enrich),
                      "rank_by_detection_ratio": " > ".join(by_ratio),
                      "Hbegf_place_abundance": by_abund.index("Hbegf") + 1,
                      "Hbegf_place_enrichment": by_enrich.index("Hbegf") + 1})
    return table, pd.DataFrame(ranks)


def replicate_table() -> pd.DataFrame:
    comp = pd.read_csv(C3_COMP)
    datp = comp[(comp["state"] == "DATP_like") & comp["arm"].isin(["4d_RFP", "2w_RFP", "4d_YFP", "2w_YFP"])]
    return datp[["library", "arm", "clone", "timepoint", "n_cells", "pct_of_library"]].reset_index(drop=True)


def retention_table() -> pd.DataFrame:
    t = pd.read_csv(C2).set_index("gene").loc[FIBROTIC]
    frame = pd.DataFrame({"gene": FIBROTIC, "flox_plus": t["Areg-flox/+"].to_numpy(),
                          "flox_flox": t["Areg-flox/flox"].to_numpy(), "c2_verdict": t["verdict"].to_numpy()})
    frame["retention_pct"] = (100 * frame["flox_flox"] / frame["flox_plus"]).round(1)
    frame["reading"] = np.select([frame["retention_pct"] >= 80, frame["retention_pct"] >= 60],
                                 ["retained", "falls by about a quarter"], "falls by more than half")
    return frame.sort_values("retention_pct", ascending=False).reset_index(drop=True)


def cluster11_check(adata) -> pd.DataFrame:
    counts = adata.layers["counts"]
    clusters = adata.obs["leiden_0.5"].astype(str).to_numpy()
    det = {g: np.asarray(counts[:, adata.var_names.get_loc(g)].todense()).ravel() > 0
           for g in CLUSTER11_GENES if g in adata.var_names}
    rows = []
    for label, mask in (("cluster 11", clusters == "11"), ("cluster 12", clusters == "12"),
                        ("cluster 14", clusters == "14"), ("whole object", np.ones(len(clusters), bool))):
        row = {"group": label, "n_cells": int(mask.sum())}
        for gene, hit in det.items():
            row[f"pct_{gene}"] = round(100 * float(hit[mask].mean()), 1)
        row["pct_Krt8_and_Col1a1"] = round(100 * float((det["Krt8"] & det["Col1a1"])[mask].mean()), 1)
        row["median_doublet_score"] = round(float(np.median(adata.obs.loc[mask, "doublet_score"])), 4)
        rows.append(row)
    return pd.DataFrame(rows)


# ---- figure 1: the sort contaminant ---------------------------------------
def figure_contaminant(adata, plt) -> None:
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D

    purity = pd.read_csv(C1D)
    purity["cluster"] = purity["cluster"].astype(str)
    reads = dict(zip(purity["cluster"], purity["reads_as"]))
    umap = adata.obsm["X_umap"]
    clusters = adata.obs["leiden_0.5"].astype(str).to_numpy()
    reading = np.array([reads.get(c, "mesenchymal") for c in clusters])
    areg = np.asarray(adata[:, "Areg"].X.todense()).ravel()

    fig = plt.figure(figsize=(17.5, 6.6), dpi=150)
    fig.text(0.012, 0.955, "A mesenchymal library that is not all mesenchyme", fontsize=15,
             fontweight="semibold", color=V.INK)
    fig.text(0.012, 0.915, "GSE316241, sorted CD45−CD31−EpCAM−, 11,690 cells after QC. One library per genotype, "
             "three mice pooled each: a description of two libraries, not a test.", fontsize=10, color=V.INK_2)

    ax = fig.add_axes([0.02, 0.06, 0.28, 0.76])
    on = ~np.isin(reading, ["epithelial", "immune"])
    ax.scatter(umap[on, 0], umap[on, 1], s=2, c=V.DEEMPH, linewidths=0, rasterized=True)
    for label, colour in (("immune", V.SLOT[2]), ("epithelial", V.SLOT[1])):
        m = reading == label
        ax.scatter(umap[m, 0], umap[m, 1], s=6, c=colour, linewidths=0, rasterized=True)
    c11 = clusters == "11"
    cx, cy = float(np.median(umap[c11, 0])), float(np.median(umap[c11, 1]))
    ax.annotate("cluster 11: 184 cells,\n99% from the Red2Kras library", xy=(cx, cy),
                xytext=(cx - 1.0, cy + 3.2), fontsize=9, color=V.INK, ha="center",
                arrowprops=dict(arrowstyle="-", color=V.INK_2, lw=0.8))
    V.bare(ax)
    ax.set_title("Clusters that are off-target for this sort")
    ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=6, mfc=V.SLOT[1], mec=V.SURFACE, label="reads as epithelial"),
                       Line2D([], [], marker="o", ls="", ms=6, mfc=V.SLOT[2], mec=V.SURFACE, label="reads as immune"),
                       Line2D([], [], marker="o", ls="", ms=6, mfc=V.DEEMPH, mec=V.SURFACE, label="mesenchymal or mesothelial")],
              loc="lower left", handletextpad=0.3)

    ax = fig.add_axes([0.335, 0.06, 0.28, 0.76])
    order = np.argsort(areg, kind="stable")
    xy, val = umap[order], areg[order]
    zero = val == 0
    ax.scatter(xy[zero, 0], xy[zero, 1], s=2, c=V.DEEMPH, linewidths=0, rasterized=True)
    points = ax.scatter(xy[~zero, 0], xy[~zero, 1], s=6, c=val[~zero],
                        cmap=LinearSegmentedColormap.from_list("ramp", V.RAMP), vmin=0,
                        vmax=float(np.quantile(val[~zero], 0.99)), linewidths=0, rasterized=True)
    V.bare(ax)
    ax.set_title("Areg expression on the same embedding")
    cax = fig.add_axes([0.345, 0.10, 0.11, 0.016])
    bar = fig.colorbar(points, cax=cax, orientation="horizontal")
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, colors=V.MUTED, length=0)
    bar.set_label("log1p(CP10K); gray = not detected", fontsize=8, color=V.INK_2)

    ax = fig.add_axes([0.715, 0.10, 0.27, 0.72])
    frame = purity.sort_values("pct_Areg", ascending=True).reset_index(drop=True)
    ax.set_xlim(0, 118)
    ax.set_ylim(-0.7, len(frame) - 0.3)
    ax.set_yticks(range(len(frame)), [f"cluster {c}" for c in frame["cluster"]])
    ax.set_xticks([0, 25, 50, 75, 100])
    V.recessive(ax)
    ax.set_xlabel("cells with Areg detected (%)")
    scale = V.px_scale(ax)
    colour_of = {"epithelial": V.SLOT[1], "immune": V.SLOT[2]}
    for i, row in frame.iterrows():
        V.hbar(ax, i, 0, float(row["pct_Areg"]), 10, colour_of.get(row["reads_as"], V.DEEMPH), scale)
        if row["reads_as"] == "epithelial":
            ax.text(float(row["pct_Areg"]) + 2, i, f"{row['pct_Areg']:.0f}%, {row['n_cells']} cells",
                    va="center", fontsize=8.5, color=V.INK)
    ax.set_title("Areg-positive share by cluster")
    fig.savefig(OUT / "c5_fig1_sort_contaminant.png")
    plt.close(fig)


# ---- figure 2: the fibrotic programme under ligand deletion ----------------
def figure_fibrotic(retention: pd.DataFrame, plt) -> None:
    from matplotlib.lines import Line2D

    frame = retention.iloc[::-1].reset_index(drop=True)
    fig = plt.figure(figsize=(14.5, 5.6), dpi=150)
    fig.text(0.012, 0.945, "Deleting Areg does not remove the fibrotic programme as one unit",
             fontsize=15, fontweight="semibold", color=V.INK)
    fig.text(0.012, 0.9, "Fibroblasts of GSE316244, the six genes of the paper's reprogrammed-fibroblast "
             "set. One library per genotype, similar depth (median 2,615 and 2,750 genes per cell).",
             fontsize=10, color=V.INK_2)
    ax = fig.add_axes([0.075, 0.13, 0.4, 0.66])
    ax.set_xlim(0, 0.62)
    ax.set_ylim(-0.6, len(frame) - 0.4)
    ax.set_yticks(range(len(frame)), frame["gene"])
    V.recessive(ax)
    ax.set_xlabel("fraction of fibroblasts with the gene detected")
    for i, row in frame.iterrows():
        ax.plot([row["flox_plus"], row["flox_flox"]], [i, i], color=V.AXIS, lw=2,
                solid_capstyle="round", zorder=2)
        ax.scatter([row["flox_plus"]], [i], s=70, color=V.FLOX_PLUS, edgecolors=V.SURFACE,
                   linewidths=2, zorder=3)
        ax.scatter([row["flox_flox"]], [i], s=70, color=V.FLOX_FLOX, edgecolors=V.SURFACE,
                   linewidths=2, zorder=4)
    ax.set_title("Detection, Areg-flox/+ to Areg-flox/flox")
    key_plus = Line2D([], [], marker="o", ls="", ms=8, mfc=V.FLOX_PLUS,
                      mec=V.SURFACE, label="Areg-flox/+")
    key_flox = Line2D([], [], marker="o", ls="", ms=8, mfc=V.FLOX_FLOX,
                      mec=V.SURFACE, label="Areg-flox/flox")
    ax.legend(handles=[key_plus, key_flox], loc="lower right")

    ax = fig.add_axes([0.575, 0.13, 0.4, 0.66])
    ax.set_xlim(0, 118)
    ax.set_ylim(-0.6, len(frame) - 0.4)
    ax.set_yticks(range(len(frame)), frame["gene"])
    ax.set_xticks([0, 20, 40, 60, 80, 100])
    V.recessive(ax)
    ax.set_xlabel("detection kept after Areg deletion (% of Areg-flox/+)")
    scale = V.px_scale(ax)
    for i, row in frame.iterrows():
        V.hbar(ax, i, 0, float(row["retention_pct"]), 14, V.SLOT[1], scale)
        label_x = float(row["retention_pct"]) + 2
        label_x = 82.0 if 70 <= label_x <= 84 else label_x
        ax.text(label_x, i, f"{row['retention_pct']:.0f}%",
                va="center", fontsize=9, color=V.INK)
    rule_ink = V.INK_2
    ax.axvline(80, color=rule_ink, lw=1, zorder=4)
    ax.text(81, len(frame) - 0.55, "C2 frozen rule: 80%", fontsize=8.5, color=rule_ink, va="bottom")
    ax.set_title("How much of each gene survives")
    fig.savefig(OUT / "c5_fig2_fibrotic_programme_split.png")
    plt.close(fig)


# ---- figure 3: EGFR ligands in the DATP-like state -------------------------
def ligand_panel(ax, table, column, xlabel, title, label_gene):
    labels = dict(MUTANT)
    libs = [lib for lib, _ in MUTANT]
    step = 1.0
    ax.set_ylim(len(libs) * step - 0.35, -0.65)
    counts = {l: int(table.loc[table["library"] == l, "n_DATP_like"].iloc[0]) for l in libs}
    ax.set_yticks([i * step for i in range(len(libs))],
                  [labels[l] + chr(10) + format(counts[l], ",") + " DATP-like cells" for l in libs])
    V.recessive(ax)
    ax.set_xlabel(xlabel)
    return libs, step


def draw_ligand_bars(ax, table, column, libs, step, label_gene):
    scale = V.px_scale(ax)
    thickness = 9
    offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * 0.2
    for i, lib in enumerate(libs):
        sub = table[table["library"] == lib].set_index("ligand")
        for k, gene in enumerate(LIGANDS):
            value = float(sub.loc[gene, column])
            y = i * step + offsets[k]
            V.hbar(ax, y, 0, max(value, 0), thickness, LIGAND_COLOUR[gene], scale)
            if gene == label_gene:
                ax.text(max(value, 0) + 0.04, y, f"{gene} {value:.2f}",
                        va="center", fontsize=8.5, color=V.INK)


TITLE_3 = "Abundance and enrichment give Hbegf different places"
SUB_3 = ("DATP-like cells of the four mutant libraries of GSE247505. Enrichment is the "
         "difference in mean log1p(CP10K) from AT2 cells of the same library.")


def figure_ligands(table, plt) -> None:
    from matplotlib.patches import Patch

    fig = plt.figure(figsize=(14.5, 6.0), dpi=150)
    fig.text(0.012, 0.95, TITLE_3, fontsize=15, fontweight="semibold", color=V.INK)
    fig.text(0.012, 0.905, SUB_3, fontsize=10, color=V.INK_2)
    panels = (
        ([0.12, 0.15, 0.35, 0.64], "mean_DATP_like", 3.4, "mean log1p(CP10K) in DATP-like cells",
         "Abundance: Hbegf second in every library", "Hbegf"),
        ([0.62, 0.15, 0.35, 0.64], "enrichment_log1p_diff", 2.4, "DATP-like minus AT2, same library (log1p CP10K)",
         "Enrichment over AT2: Hbegf second in one of four", "Ereg"),
    )
    for rect, column, xmax, xlabel, title, label_gene in panels:
        ax = fig.add_axes(rect)
        ax.set_xlim(0, xmax)
        libs, step = ligand_panel(ax, table, column, xlabel, title, label_gene)
        draw_ligand_bars(ax, table, column, libs, step, label_gene)
        ax.set_title(title)
    keys = [Patch(facecolor=LIGAND_COLOUR[g], edgecolor="none", label=g) for g in LIGANDS]
    fig.legend(handles=keys, loc="lower center", ncol=4, bbox_to_anchor=(0.5, 0.0),
               handlelength=1.2, columnspacing=1.6)
    fig.savefig(OUT / "c5_fig3_egfr_ligands_abundance_vs_enrichment.png")
    plt.close(fig)


def main() -> None:
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    V.apply(plt)
    rec = RunRecord(OUT / "c5_run_record.json", "C5 figures for the three findings",
                    RULES, notes="plots and derives from tracked tables; fits nothing")
    for path in (C1D, C2, C3_LIG, C3_COMP, OBJECT):
        rec.add_input(path)

    ligands, ranks = ligand_tables()
    replicates = replicate_table()
    retention = retention_table()
    adata = ad.read_h5ad(OBJECT)
    check11 = cluster11_check(adata)
    plotted = ranks[ranks["library"].isin(dict(MUTANT))]
    hbegf_second_abundance = int((plotted["Hbegf_place_abundance"] == 2).sum())
    hbegf_second_enrichment = int((plotted["Hbegf_place_enrichment"] == 2).sum())
    if (hbegf_second_abundance, hbegf_second_enrichment) != (4, 1):
        raise SystemExit("figure 3 titles no longer match the table; refusing to draw them")

    outputs = {"c5_ligand_abundance_and_enrichment.csv": ligands,
               "c5_ligand_rankings.csv": ranks,
               "c5_datp_like_share_per_library.csv": replicates,
               "c5_fibrotic_retention.csv": retention,
               "c5_cluster11_discriminating_check.csv": check11}
    for name, frame in outputs.items():
        frame.to_csv(OUT / name, index=False)
        rec.add_output(OUT / name)

    figure_contaminant(adata, plt)
    figure_fibrotic(retention, plt)
    figure_ligands(ligands[ligands["plotted"]], plt)
    for name in ("c5_fig1_sort_contaminant.png", "c5_fig2_fibrotic_programme_split.png",
                 "c5_fig3_egfr_ligands_abundance_vs_enrichment.png"):
        rec.add_output(OUT / name)
    rec.set("hbegf_second_by_abundance_in_libraries", hbegf_second_abundance)
    rec.set("hbegf_second_by_enrichment_in_libraries", hbegf_second_enrichment)
    rec.set("datp_like_share_4d_replicates", replicates[replicates["arm"] == "4d_RFP"][["library", "n_cells", "pct_of_library"]].to_dict("records"))
    rec.set("retention", retention[["gene", "retention_pct", "reading"]].to_dict("records"))
    rec.set("cluster11_check", check11.to_dict("records"))

    lines = ["# Trial C5 output: figures for the three findings, and the corrections they forced", "",
             "## Ligand rankings by three notions", "", df_to_markdown(ranks, index=False), "",
             "## DATP-like share per library", "", df_to_markdown(replicates, index=False), "",
             "## Fibrotic gene retention after Areg deletion", "", df_to_markdown(retention, index=False), "",
             "## Cluster 11 discriminating check", "", df_to_markdown(check11, index=False), ""]
    (OUT / "c5_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c5_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
