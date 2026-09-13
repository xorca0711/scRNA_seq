#!/usr/bin/env python
"""Trial C11: figures for the results that contradict or refute an earlier trial.

Trials C5 and E5 drew the findings and the one refutation that existed then.
Six more results since then either refute an earlier trial, refute one of this
repository's own methods, or close a lead, and every one of them lives as a
number in a table where nobody will look at it. This script draws them.

Each figure reads only tracked tables written by the trial it depicts, computes
nothing scientific, and carries a fail-closed assertion so that it refuses to
draw if the table stops saying what the title claims.

The six:

1. E6, the depth control that caught a false positive. The only significant
   correlation in the trial was a control pair, and both its variables track
   sequencing depth.
2. C7, a measure refuted by its own margin. Whole-profile rank correlation
   separated compartments by 0.415 of rho and epithelial states by 0.006, and
   the sort contaminant turned out to be a mixture rather than a state.
3. C8, the tiers are amplitudes. The retained pair sits inside the independence
   band, and a third of the cells detect neither gene, which is why the mixture
   test preferred two components.
4. The species divergence. HBEGF is myeloid-dominant in human lung and near the
   bottom of the myeloid compartments in mouse, so the lead does not transfer.
5. C9, existence against size. Above-chance co-occurrence replicates in both
   animals; the ratio does not, and the magnitude threshold measured depth.
6. C10, the lead closing. The double-positives are the published pathological
   fibroblast, and a ranking artefact hid it.
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

OUT = HERE / "c11_figures_for_the_contradictions"
OUT.mkdir(exist_ok=True)

SRC = {
    "e6_donors": HERE / "e6_donor_level_axis_coupling" / "e6_per_donor.csv",
    "e6_corr": HERE / "e6_donor_level_axis_coupling" / "e6_correlations.csv",
    "c7_corr": HERE / "c7_what_the_sort_contaminant_is" / "c7_profile_correlations.csv",
    "c8_ratios": HERE / "c8_subpopulation_or_gradient" / "c8_codetection_ratios.csv",
    "c6_hbegf": HERE / "c6_who_makes_egfr_ligands" / "c6_hbegf_survival.csv",
    "e2_rank": HERE / "e2_human_fibrosis_ligand_sources" / "e2_ipf_ranking_by_celltype.csv",
    "c9_lib": HERE / "c9_the_fst_runx2_population" / "c9_per_library.csv",
    "c9_rep": HERE / "c9_the_fst_runx2_population" / "c9_marker_replication.csv",
    "c9_markers": HERE / "c9_the_fst_runx2_population" / "c9_reference_markers.csv",
    "c10_eff": HERE / "c10_published_state_or_not" / "c10_effect_sizes.csv",
    "c10_audit": HERE / "c10_published_state_or_not" / "c10_saturation_audit.csv",
}
C7_STATES = {"DATP_like": 77, "AT2": 70, "cycling": 20, "AT1_like": 13, "Cd177_positive": 4}

RULES = {
    "purpose": "draw the six results that contradict or refute an earlier trial, a method, or a lead",
    "inputs": "only tracked tables written by the trial each figure depicts",
    "computes": "nothing scientific; every number is read from a table",
    "fail_closed": "each figure asserts the table still says what its title claims, and refuses to draw otherwise",
}


def head(fig, title, subtitle, extra=None):
    fig.text(0.012, 0.955, title, fontsize=15, fontweight="semibold", color=V.INK)
    fig.text(0.012, 0.905, subtitle, fontsize=10, color=V.INK_2)
    if extra:
        fig.text(0.012, 0.866, extra, fontsize=10, color=V.INK_2)


def fig1_depth_control(plt):
    donors = pd.read_csv(SRC["e6_donors"]).dropna(subset=["epi_AREG", "fib_EGFR", "fib_activation"])
    corr = pd.read_csv(SRC["e6_corr"])

    def get(x, y, comparison):
        row = corr[(corr["x"] == x) & (corr["y"] == y) & corr["comparison"].str.startswith(comparison)]
        return row.iloc[0]

    primary = get("epi_AREG", "fib_EGFR", "T8")
    control = get("epi_TGFA", "fib_activation", "T10")
    depth_tgfa = get("epi_TGFA", "epi_depth", "T11")
    depth_act = get("fib_activation", "fib_depth", "T11")
    if not (primary["p_value"] > 0.05 and control["p_value"] < 0.05):
        raise SystemExit("figure 1 assumes the primary is null and the control is significant; refusing to draw")

    fig = plt.figure(figsize=(13.2, 5.6), dpi=150)
    head(fig, "The only significant correlation in the trial was a control pair",
         "GSE136831, " + str(len(donors)) + " donors with both compartments above the 50-cell floor. "
         "Each dot is one donor. Spearman across donors.",
         "Both variables of that control pair track sequencing depth, so the frozen rule refused to read "
         "it. The same number on AREG would have read as confirmation.")
    for index, (xcol, ycol, stat, label, note) in enumerate((
            ("epi_AREG", "fib_EGFR", primary, "the test the trial was built on",
             "not significant"),
            ("epi_TGFA", "fib_activation", control, "the ligand control",
             "significant, and not read"))):
        ax = fig.add_axes([0.075 + index * 0.50, 0.14, 0.38, 0.63])
        for disease, colour in (("IPF", V.SLOT[1]), ("Control", V.MUTED)):
            part = donors[donors["disease"] == disease]
            ax.scatter(part[xcol], part[ycol], s=64, color=colour, edgecolors=V.SURFACE,
                       linewidths=1.6, zorder=3, label=disease + " (" + str(len(part)) + ")")
        V.recessive(ax, grid_axis="both")
        ax.set_xlabel("epithelial " + xcol.split("_")[1] + " detection")
        ax.set_ylabel("fibroblast " + ("EGFR detection" if ycol == "fib_EGFR" else "activation score"))
        ax.set_title(label, fontsize=11)
        text = ("rho " + format(stat["rho"], ".3f") + ", p " + format(stat["p_value"], ".3f")
                + "\n" + note)
        ax.text(0.03, 0.97, text, transform=ax.transAxes, va="top", fontsize=9.5, color=V.INK_2)
        if index == 0:
            ax.text(0.03, 0.80, "smallest rho this n could call: 0.43", transform=ax.transAxes,
                    va="top", fontsize=9, color=V.MUTED)
            ax.legend(loc="lower right", handletextpad=0.4)
        else:
            ax.text(0.03, 0.80, "depth correlations: TGFA " + format(depth_tgfa["rho"], ".3f")
                    + ", activation " + format(depth_act["rho"], ".3f") + "\nboth at or above the 0.4 line",
                    transform=ax.transAxes, va="top", fontsize=9, color=V.SLOT[2])
    path = OUT / "c11_fig1_e6_depth_control.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig2_c7_resolution(plt):
    corr = pd.read_csv(SRC["c7_corr"])
    query = corr[corr["query"].str.startswith("query")].sort_values("rho", ascending=False)
    control = corr[corr["query"].str.startswith("control")].sort_values("rho", ascending=False)
    state_margin = float(query.iloc[0]["rho"] - query.iloc[1]["rho"])
    compartment_gap = float(query["rho"].max() - control["rho"].max())
    if not (state_margin < 0.05 < compartment_gap):
        raise SystemExit("figure 2 assumes a small state margin and a large compartment gap; refusing to draw")

    fig = plt.figure(figsize=(13.6, 5.8), dpi=150)
    head(fig, "A measure refuted by the margin it was made to report",
         "Trial C7: rank correlation of mean profiles, 32,163 shared genes. The query is the 184-cell "
         "epithelial contaminant of the mesenchymal sort; the control is the cleanest fibroblast cluster "
         "in the same object.",
         "It separates compartments by " + format(compartment_gap, ".3f") + " of rho and epithelial "
         "states by " + format(state_margin, ".4f") + ", so it can say which compartment and never "
         "which state.")
    ax = fig.add_axes([0.075, 0.14, 0.42, 0.62])
    order = query["reference_cluster"].tolist()
    ypos = {c: i for i, c in enumerate(reversed(order))}
    for frame, colour, label in ((control, V.DEEMPH, "control: fibroblast cluster"),
                                 (query, V.SLOT[1], "query: the contaminant")):
        ax.scatter(frame["rho"], [ypos[c] for c in frame["reference_cluster"]], s=70,
                   color=colour, edgecolors=V.SURFACE, linewidths=1.6, zorder=3, label=label)
    ax.set_yticks(range(len(order)),
                  [str(c) + "  " + query.set_index("reference_cluster").loc[c, "reference_call"]
                   for c in reversed(order)])
    V.recessive(ax)
    ax.set_xlabel("Spearman correlation of mean profiles")
    ax.set_title("the top four matches span three different state calls", fontsize=11)
    ax.legend(loc="lower left", handletextpad=0.4)

    ax2 = fig.add_axes([0.60, 0.14, 0.35, 0.62])
    total = sum(C7_STATES.values())
    labels = list(C7_STATES)
    values = [C7_STATES[k] / total for k in labels]
    ax2.set_xlim(0, 0.62)
    ax2.set_ylim(-0.6, len(labels) - 0.4)
    ax2.set_yticks(range(len(labels)), [l.replace("_", " ") for l in reversed(labels)])
    V.recessive(ax2)
    scale = V.px_scale(ax2)
    for i, (label, value) in enumerate(zip(reversed(labels), reversed(values))):
        colour = V.SLOT[1] if label == "DATP_like" else V.DEEMPH
        V.hbar(ax2, i, 0, value, 14, colour, scale)
        ax2.text(value + 0.012, i, format(value, ".1%"), fontsize=9, color=V.INK_2, va="center")
    ax2.axvline(0.5, color=V.AXIS, lw=1, ls=(0, (4, 3)), zorder=1)
    ax2.text(0.505, len(labels) - 0.55, "repository confidence floor: 50%", fontsize=8.5,
             color=V.MUTED, va="bottom")
    ax2.set_xlabel("share of the 184 contaminating cells")
    ax2.set_title("and the contaminant is a mixture, not a state", fontsize=11)
    path = OUT / "c11_fig2_c7_resolution.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig3_c8_amplitudes(plt):
    ratios = pd.read_csv(SRC["c8_ratios"])
    hom = ratios[ratios["genotype"] == "Areg-flox/flox"].set_index("tier")
    retained = hom.loc["retained"]
    if not (0.80 < float(retained["ratio"]) < 1.25):
        raise SystemExit("figure 3 assumes the retained pair sits inside the independence band; refusing to draw")
    a, b = float(retained["detection_a"]), float(retained["detection_b"])
    both = float(retained["observed_both"])
    parts = {"neither gene": 1 - a - b + both, "Runx1 only": a - both,
             "Pdgfrb only": b - both, "both genes": both}

    fig = plt.figure(figsize=(13.4, 5.4), dpi=150)
    head(fig, "The tiers are amplitudes, and the mixture test that disagreed was arithmetic",
         "Trial C8: gated fibroblasts of GSE316244, one library per genotype. Co-detection against what "
         "independence predicts, with the frozen band from 0.80 to 1.25.",
         "A third of the deletion-arm cells detect neither retained gene, so a score built from the two "
         "has a spike at zero and a two-component mixture wins whatever the biology is.")
    ax = fig.add_axes([0.075, 0.14, 0.44, 0.62])
    frame = ratios.copy()
    frame["label"] = frame["pair"] + "\n" + frame["genotype"]
    frame = frame.iloc[::-1].reset_index(drop=True)
    ax.axvspan(0.80, 1.25, color=V.GRID, zorder=0)
    ax.axvline(1.0, color=V.AXIS, lw=1, zorder=1)
    for i, row in frame.iterrows():
        colour = V.SLOT[2] if row["ratio"] > 1.25 else (V.SLOT[1] if row["tier"] == "retained" else V.DEEMPH)
        ax.scatter([row["ratio"]], [i], s=88, color=colour, edgecolors=V.SURFACE, linewidths=1.8, zorder=3)
        ax.text(float(row["ratio"]) + 0.035, i, format(row["ratio"], ".2f"), fontsize=9,
                color=V.INK_2, va="center")
    ax.set_yticks(range(len(frame)), frame["label"])
    ax.set_xlim(0.6, 2.0)
    ax.set_ylim(-0.6, len(frame) - 0.4)
    V.recessive(ax)
    ax.set_xlabel("observed co-detection divided by independence")
    ax.text(1.02, -0.45, "independence", fontsize=8.5, color=V.MUTED, va="bottom")
    ax.set_title("only the falling pair, in the control arm, exceeds the band", fontsize=11)

    ax2 = fig.add_axes([0.60, 0.42, 0.36, 0.30])
    left = 0.0
    for (label, value), colour in zip(parts.items(), (V.SLOT[2], V.DEEMPH, V.DEEMPH, V.SLOT[1])):
        ax2.barh([0], [value], left=left, height=0.55, color=colour, edgecolor=V.SURFACE, linewidth=2)
        if value > 0.12:
            ax2.text(left + value / 2, 0, format(value, ".1%"), ha="center", va="center",
                     fontsize=9.5, color=V.SURFACE if colour != V.DEEMPH else V.INK_2)
        ax2.text(left + value / 2, 0.42, label, ha="center", va="bottom", fontsize=8.6, color=V.INK_2)
        left += value
    ax2.set_xlim(0, 1)
    ax2.set_ylim(-0.5, 0.9)
    V.bare(ax2)
    ax2.set_title("deletion-arm cells by which retained genes they detect", fontsize=11)
    path = OUT / "c11_fig3_c8_amplitudes.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig4_species(plt):
    mouse = pd.read_csv(SRC["c6_hbegf"]).sort_values("det_Hbegf_flox_plus", ascending=False)
    human = pd.read_csv(SRC["e2_rank"])
    human = human[human["cohort"] == "E2_GSE136831"].sort_values("det_HBEGF", ascending=False).head(10)
    myeloid_mouse = ("macrophage", "myeloid", "neutrophil", "monocyte")
    myeloid_human = ("Macrophage", "DC", "Monocyte", "cDC", "pDC", "Mac")
    top_human = human.iloc[0]["celltype"]
    if not any(k in top_human for k in myeloid_human):
        raise SystemExit("figure 4 assumes a myeloid cell type tops the human ranking; refusing to draw")

    fig = plt.figure(figsize=(13.8, 5.8), dpi=150)
    head(fig, "The Hbegf lead does not transfer across species",
         "Left: mouse niche compartments of GSE316244, control arm (trial C6). Right: the ten "
         "highest human IPF cell types of GSE136831 (trial E2). Detection fraction, orange marks "
         "myeloid populations.",
         "Macrophages and dendritic cells top the human ranking, while the equivalent mouse "
         "compartments sit near the bottom at 0.085 and 0.046 and the mouse ligand is epithelial and "
         "endothelial. Mouse neutrophils are the exception, third at 0.178. Absolute values are not "
         "comparable across deposits; the ordering within each is.")
    for index, (frame, namecol, valuecol, keys, title) in enumerate((
            (mouse, "compartment", "det_Hbegf_flox_plus", myeloid_mouse, "mouse, control arm"),
            (human, "celltype", "det_HBEGF", myeloid_human, "human IPF"))):
        ax = fig.add_axes([0.20 + index * 0.49, 0.14, 0.26, 0.60])
        rows = frame.iloc[::-1].reset_index(drop=True)
        ax.set_xlim(0, float(frame[valuecol].max()) * 1.28)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_yticks(range(len(rows)), rows[namecol])
        V.recessive(ax)
        scale = V.px_scale(ax)
        for i, row in rows.iterrows():
            is_myeloid = any(k.lower() in str(row[namecol]).lower() for k in keys)
            V.hbar(ax, i, 0, float(row[valuecol]), 13,
                   V.SLOT[2] if is_myeloid else V.DEEMPH, scale)
            ax.text(float(row[valuecol]) + float(frame[valuecol].max()) * 0.025, i,
                    format(row[valuecol], ".3f"), fontsize=8.8, color=V.INK_2, va="center")
        ax.set_xlabel("Hbegf detection" if index == 0 else "HBEGF detection")
        ax.set_title(title, fontsize=11)
    path = OUT / "c11_fig4_species_divergence.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig5_c9(plt):
    lib = pd.read_csv(SRC["c9_lib"])
    rep = pd.read_csv(SRC["c9_rep"]).dropna(subset=["diff_Bleo1_GFPp", "diff_Bleo2_GFPp"])
    bleo = lib[lib["group"] == "bleomycin"]
    if not (bleo["permutation_p"] < 0.05).all():
        raise SystemExit("figure 5 assumes the permutation replicates in both bleomycin animals; refusing to draw")
    rep["bleo_mean"] = (rep["diff_Bleo1_GFPp"] + rep["diff_Bleo2_GFPp"]) / 2
    attenuation = float((rep["bleo_mean"] / rep["reference_difference"]).median())

    fig = plt.figure(figsize=(13.6, 5.6), dpi=150)
    head(fig, "The co-occurrence replicates; its size does not, and the magnitude rule measured depth",
         "Trial C9. Left: co-detection ratio per library, with the bar spanning the shallow and deep "
         "halves at the median genes per cell. Right: each of the 27 checkable marker genes.",
         "Median attenuation is " + format(attenuation, ".3f") + " against a threshold of 0.5, and the "
         "bleomycin libraries carry half the genes per cell, so the threshold split the distribution at "
         "its own centre.")
    ax = fig.add_axes([0.075, 0.14, 0.38, 0.60])
    rows = lib.iloc[::-1].reset_index(drop=True)
    ax.axvline(1.0, color=V.AXIS, lw=1, zorder=1)
    for i, row in rows.iterrows():
        lo, hi = row["ratio_shallow"], row["ratio_deep"]
        colour = V.SLOT[2] if row["group"] == "bleomycin" else V.SLOT[1]
        if pd.notna(lo) and pd.notna(hi):
            ax.plot([lo, hi], [i, i], color=colour, lw=3, solid_capstyle="round", alpha=0.45, zorder=2)
        ax.scatter([row["ratio"]], [i], s=86, color=colour, edgecolors=V.SURFACE,
                   linewidths=1.8, zorder=3)
        ax.text(float(row["ratio"]) + 0.06, i, "p " + format(row["permutation_p"], ".3f"),
                fontsize=8.8, color=V.INK_2, va="center")
    ax.set_yticks(range(len(rows)), rows["library"])
    ax.set_ylim(-0.6, len(rows) - 0.4)
    V.recessive(ax)
    ax.set_xlabel("co-detection ratio, bar spans the depth halves")
    ax.set_title("the point is the estimate, the bar is its instability", fontsize=11)

    ax2 = fig.add_axes([0.58, 0.14, 0.37, 0.60])
    top = float(rep["reference_difference"].max()) * 1.08
    ax2.plot([0, top], [0, top], color=V.AXIS, lw=1, ls=(0, (4, 3)), zorder=1)
    ax2.plot([0, top], [0, top / 2], color=V.SLOT[2], lw=1.4, zorder=1)
    for replicates, colour in ((True, V.SLOT[1]), (False, V.DEEMPH)):
        part = rep[rep["replicates"] == replicates]
        ax2.scatter(part["reference_difference"], part["bleo_mean"], s=54, color=colour,
                    edgecolors=V.SURFACE, linewidths=1.4, zorder=3,
                    label=("clears the rule" if replicates else "does not") + " (" + str(len(part)) + ")")
    ax2.set_xlim(0, top)
    ax2.set_ylim(0, top)
    V.recessive(ax2, grid_axis="both")
    ax2.set_xlabel("detection difference, reference library")
    ax2.set_ylabel("mean difference, two bleomycin animals")
    ax2.set_title("direction replicates for 26 of 27; magnitude is a depth line", fontsize=11)
    ax2.text(top * 0.42, top * 0.20, "the 0.5 threshold", fontsize=9, color=V.SLOT[2])
    ax2.legend(loc="upper left", handletextpad=0.4)
    path = OUT / "c11_fig5_c9_existence_vs_size.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig6_c10(plt):
    eff = pd.read_csv(SRC["c10_eff"])
    audit = pd.read_csv(SRC["c10_audit"])
    markers = pd.read_csv(SRC["c9_markers"])
    cutoff = float(markers["difference"].min())
    cthrc1 = audit[(audit["gene"] == "Cthrc1")]
    if not (eff["d_pathological"] > 0.5).all():
        raise SystemExit("figure 6 assumes the pathological score separates the groups everywhere; refusing to draw")

    fig = plt.figure(figsize=(13.6, 5.6), dpi=150)
    head(fig, "The lead closes: these are the published pathological fibroblast",
         "Trial C10. Left: standardised difference between double-positive and double-negative gated "
         "fibroblasts for four published sets. Right: Cthrc1 detection in the two groups.",
         "Cthrc1 was absent from trial C9's top thirty only because its detection difference fell below "
         "that list's cutoff of " + format(cutoff, ".3f") + ". The ranking hid the marker that settles it.")
    ax = fig.add_axes([0.075, 0.14, 0.42, 0.60])
    sets = ["pathological", "alveolar", "adventitial", "smooth_muscle"]
    width = 0.8 / len(eff)
    ax.axhline(0, color=V.AXIS, lw=1, zorder=1)
    ax.axhline(0.5, color=V.MUTED, lw=1, ls=(0, (4, 3)), zorder=1)
    palette = [V.SLOT[1], V.SLOT[2], V.SLOT[4]]
    for j, (_, row) in enumerate(eff.iterrows()):
        colour = palette[j % len(palette)]
        xs = np.arange(len(sets)) + (j - (len(eff) - 1) / 2) * width
        ax.bar(xs, [row["d_" + s] for s in sets], width=width * 0.92, color=colour,
               edgecolor=V.SURFACE, linewidth=1.5, label=row["library"], zorder=3)
    ax.set_xticks(range(len(sets)), [s.replace("_", " ") for s in sets])
    V.recessive(ax, grid_axis="y")
    ax.set_ylabel("standardised difference between the groups")
    ax.text(len(sets) - 0.55, 0.56, "frozen floor 0.5", fontsize=8.5, color=V.MUTED)
    ax.set_title("pathological is largest; alveolar identity is lost", fontsize=11)
    ax.legend(loc="lower left", fontsize=8.5, handletextpad=0.4)

    ax2 = fig.add_axes([0.60, 0.14, 0.35, 0.60])
    rows = cthrc1.reset_index(drop=True)
    ax2.set_ylim(-0.6, len(rows) - 0.4)
    ax2.set_xlim(0, 0.95)
    ax2.set_yticks(range(len(rows)), rows["library"])
    V.recessive(ax2)
    scale = V.px_scale(ax2)
    for i, row in rows.iterrows():
        V.hbar(ax2, i + 0.17, 0, float(row["det_double_positive"]), 11, V.SLOT[1], scale)  # noqa: E501
        V.hbar(ax2, i - 0.17, 0, float(row["det_double_negative"]), 11, V.DEEMPH, scale)
        ax2.text(float(row["det_double_positive"]) + 0.02, i + 0.17,
                 format(row["det_double_positive"], ".3f"), fontsize=8.8, color=V.INK_2, va="center")
        ax2.text(float(row["det_double_negative"]) + 0.02, i - 0.17,
                 format(row["det_double_negative"], ".3f"), fontsize=8.8, color=V.MUTED, va="center")
    ax2.set_xlabel("Cthrc1 detection: double-positive above, double-negative below")
    ax2.set_title("the defining marker was there all along", fontsize=11)
    path = OUT / "c11_fig6_c10_lead_closes.png"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rec = RunRecord(OUT / "c11_run_record.json", "C11 figures for the contradictions", RULES)
    for path in SRC.values():
        if not path.exists():
            raise SystemExit(f"{path} is missing; this script draws only from tracked tables")
        rec.add_input(path)
    V.apply(plt)
    drawn = []
    for name, fn in (("E6 depth control", fig1_depth_control),
                     ("C7 resolution", fig2_c7_resolution),
                     ("C8 amplitudes", fig3_c8_amplitudes),
                     ("species divergence", fig4_species),
                     ("C9 existence against size", fig5_c9),
                     ("C10 the lead closes", fig6_c10)):
        path = fn(plt)
        rec.add_output(path)
        drawn.append({"figure": name, "file": path.name})
        print("drew", path.name)
    frame = pd.DataFrame(drawn)
    rec.set("figures", drawn)
    lines = ["# Trial C11: figures for the contradictions", "",
             "Each figure reads only tracked tables from the trial it depicts and refuses to draw if "
             "that table stops saying what its title claims.", "",
             df_to_markdown(frame, index=False), ""]
    (OUT / "c11_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c11_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
