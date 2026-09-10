#!/usr/bin/env python
"""Phase-wise view of the Ki67 atlas: per-dpi UMAP, lineage composition and
lineage-level proliferation, from tracked artefacts only.

Niethamer et al. 2025 (Cell Stem Cell, doi:10.1016/j.stem.2024.12.002) describe
regeneration after H1N1 injury as three proliferative phases: immune cells
first (2 to 6 dpi), epithelium and mesenchyme second (7 to 15 dpi), endothelium
last (14 to 22 dpi), with every lineage back at baseline by 21 to 28 dpi (their
Figures 1 and 2). This script asks whether that phase structure is visible in
the reanalysed atlas using only the tracked per-cell metadata table
(analysis/GSE262927/tables/cell_metadata.csv). No expression matrix is read
and nothing is re-embedded.

Three readouts, all per animal:

1. The atlas UMAP split by day post infection, coloured by deposited lineage,
   with the same number of cells drawn in every panel so that density is
   comparable between days.
2. Lineage composition per animal and day. Immune versus non-immune fractions
   are set by the MACS recombination (10 to 15 percent CD45-positive), so they
   are reported as descriptive only; within-compartment fractions are the
   interpretable ones.
3. Proliferation per lineage: (a) the Ki67 trace, i.e. the fraction of
   reporter-scored cells that are Traced, read in each cohort's immediate
   window (tamoxifen at 2, 7, 14 or 21 dpi, harvest four days later at 6, 11,
   19 or 25 dpi) and at the common 42 dpi harvest; (b) the deposited
   cell-cycle phase call (S or G2M; the authors' Seurat CellCycleScoring
   output, read from their metadata table and never recalculated) per day.

Rules are frozen in the run record before the table is opened. The deposited
lineage labels are used as given: this is a descriptive reproduction of the
paper's phase structure, not a blind validation. The statistical unit is the
animal; medians are reported and no P values are computed, because the
active-repair days carry two animals each.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from focus_utils import RunRecord, df_to_markdown, strip_by_day  # noqa: E402
from pipeline_utils import ANALYSIS, RANDOM_SEED, save_fig, setup_matplotlib  # noqa: E402

plt = setup_matplotlib()

SERIES = ANALYSIS / "GSE262927"
META = SERIES / "tables" / "cell_metadata.csv"
OUT = SERIES / "phase_timecourse"
FIG = OUT / "figures"
TAB = OUT / "tables"

LINEAGES = ["Endothelium", "Mesenchyme", "Epithelium", "Myeloid", "Lymphoid"]
LINEAGE_COLOURS = {"Endothelium": "#1f77b4", "Epithelium": "#ff7f0e", "Lymphoid": "#2ca02c",
                   "Mesenchyme": "#d62728", "Myeloid": "#9467bd"}
COMPARTMENT = {"Endothelium": "CD45neg", "Mesenchyme": "CD45neg", "Epithelium": "CD45neg",
               "Myeloid": "CD45pos", "Lymphoid": "CD45pos"}
DAYS = [0, 6, 11, 19, 25, 42, 90, 366]
PHASE = {0: "baseline", 6: "active repair", 11: "active repair", 19: "active repair",
         25: "active repair", 42: "injury resolution", 90: "injury resolution",
         366: "long-term homeostasis"}
IMMEDIATE_DAYS = [6, 11, 19, 25]
IMMEDIATE_WINDOWS = {2: 6, 7: 11, 14: 19, 21: 25}  # tamoxifen start day -> harvest day
PAPER_PEAK = {"Myeloid": "2 to 6 dpi (immune)", "Lymphoid": "2 to 6 dpi (immune)",
              "Epithelium": "7 to 15 dpi", "Mesenchyme": "7 to 15 dpi",
              "Endothelium": "14 to 22 dpi"}
PAPER_PEAK_DAYS = {"Myeloid": [6], "Lymphoid": [6], "Epithelium": [11, 19],
                   "Mesenchyme": [11, 19], "Endothelium": [19, 25]}

RULES = {
    "input": "analysis/GSE262927/tables/cell_metadata.csv (tracked per-cell metadata and atlas "
             "UMAP coordinates written by run_scrna_analysis.py); no expression data are read",
    "cohort": "has_author_metadata == True: the 25-sample Ki67 atlas. The 8 pre-labelled tracing "
              "samples carry no day post infection and are excluded",
    "lineage": "deposited author_lineage used as given (descriptive; not a blind validation)",
    "compartments": COMPARTMENT,
    "days": DAYS,
    "phases": {str(k): v for k, v in PHASE.items()},
    "umap_panels": "shared atlas coordinates; every panel draws all cohort cells in grey and exactly "
                   "N cells of that day in colour, N = the smallest per-day cell count, sampled with "
                   "seed 0 and drawn in shuffled order",
    "composition": "per animal: percent of each lineage among all cells (sort-engineered, descriptive "
                   "only) and among its compartment (CD45neg = Endothelium + Mesenchyme + Epithelium; "
                   "CD45pos = Myeloid + Lymphoid)",
    "trace": "cells with trace_call in {Traced, Untraced}; Not_detected is reporter dropout and is "
             "excluded, the rule of 07_lineage_tracing_cohort.py; pct_traced = 100 * Traced / "
             "(Traced + Untraced) per animal and lineage",
    "trace_immediate_windows": {f"tamoxifen {k} dpi": f"harvest {v} dpi" for k, v in IMMEDIATE_WINDOWS.items()},
    "trace_common_harvest": "42 dpi: four tamoxifen windows with two animals each; 90 dpi (one animal "
                            "per window) and 366 dpi are shown as single points",
    "cycling": "cell_cycle_phase in {S, G2M} as percent of a lineage's cells per animal and day; the "
               "column is the authors' Seurat CellCycleScoring call read from the deposited metadata "
               "(docs/PIPELINE_AS_RUN.md), never recalculated here",
    "peak_rule": "for each lineage, the immediate-window day (6, 11, 19, 25) with the highest median "
                 "per-animal pct_traced is the observed peak; it agrees with the paper if it falls in "
                 "the expected set. The same rule is applied to pct_cycling over the same four days",
    "paper_expected_peak": PAPER_PEAK,
    "paper_expected_peak_days": PAPER_PEAK_DAYS,
    "unit": "animal; medians reported; no P values (two animals per active-repair day)",
    "seed": RANDOM_SEED,
}


def load_cohort(rec: RunRecord) -> pd.DataFrame:
    if not META.exists():
        raise SystemExit(f"input table missing: {META}")
    rec.add_input(META)
    cols = ["cell_id", "sample_id", "experimental_group", "tamoxifen_start_day", "sacrifice_day",
            "sex", "trace_call", "cell_cycle_phase", "has_author_metadata", "author_lineage",
            "author_celltype", "leiden_cluster", "UMAP_1", "UMAP_2"]
    dtypes = {c: str for c in cols if c not in ("UMAP_1", "UMAP_2")}
    df = pd.read_csv(META, usecols=cols, dtype=dtypes, low_memory=False)
    n_all = len(df)
    df = df[df["has_author_metadata"] == "True"].copy()
    df["day"] = pd.to_numeric(df["sacrifice_day"], errors="coerce")
    if df["day"].isna().any():
        raise SystemExit("cohort cells without sacrifice_day; refusing to guess")
    df["day"] = df["day"].astype(int)
    df["tam"] = pd.to_numeric(df["tamoxifen_start_day"], errors="coerce")
    if df["tam"].isna().any():
        raise SystemExit("cohort cells without tamoxifen_start_day; refusing to guess")
    df["tam"] = df["tam"].astype(int)
    found_days = set(df["day"].unique())
    if found_days != set(DAYS):
        raise SystemExit(f"day set {sorted(found_days)} differs from the frozen {DAYS}")
    if df["author_lineage"].isna().any():
        raise SystemExit("cohort cells without a deposited lineage")
    found_lin = set(df["author_lineage"].unique())
    if found_lin != set(LINEAGES):
        raise SystemExit(f"lineage set {sorted(found_lin)} differs from the frozen {LINEAGES}")
    rec.set("n_cells_in_table", n_all)
    rec.set("n_cells_cohort", int(len(df)))
    rec.set("n_animals_cohort", int(df["sample_id"].nunique()))
    per_day = (df.groupby("day").agg(n_animals=("sample_id", "nunique"), n_cells=("sample_id", "size"))
                 .reindex(DAYS))
    per_day["phase"] = [PHASE[d] for d in per_day.index]
    rec.set("per_day", {str(d): {"n_animals": int(r.n_animals), "n_cells": int(r.n_cells)}
                        for d, r in per_day.iterrows()})
    print(f"cohort: {len(df):,} cells, {df['sample_id'].nunique()} animals")
    print(per_day.to_string())
    return df


def umap_by_day(df: pd.DataFrame, rec: RunRecord) -> None:
    n_min = int(df.groupby("day").size().min())
    rng = np.random.default_rng(RANDOM_SEED)
    fig, axes = plt.subplots(2, 4, figsize=(16, 8.6))
    for ax, d in zip(axes.ravel(), DAYS):
        ax.scatter(df["UMAP_1"], df["UMAP_2"], s=0.3, c="#d9d9d9", linewidths=0, rasterized=True)
        sub = df[df["day"] == d]
        pick = sub.iloc[rng.choice(len(sub), n_min, replace=False)]
        pick = pick.iloc[rng.permutation(len(pick))]
        ax.scatter(pick["UMAP_1"], pick["UMAP_2"], s=1.4,
                   c=pick["author_lineage"].map(LINEAGE_COLOURS), linewidths=0, rasterized=True)
        ax.set_title(f"{d} dpi, {PHASE[d]}\n{sub['sample_id'].nunique()} animals, "
                     f"{len(sub):,} cells, {n_min:,} drawn", fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    handles = [plt.Line2D([], [], marker="o", ls="", color=LINEAGE_COLOURS[k], label=k, markersize=6)
               for k in LINEAGES]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.suptitle("Ki67 atlas by day post infection (shared UMAP coordinates, deposited lineage; "
                 "equal cell number per panel)", y=0.995)
    fig.tight_layout(rect=(0, 0.035, 1, 0.97))
    for p in save_fig(fig, FIG, "UMAP_atlas_by_dpi_lineage"):
        rec.add_output(p)
    rec.set("umap_cells_drawn_per_panel", n_min)


def composition(df: pd.DataFrame, rec: RunRecord) -> tuple[pd.DataFrame, pd.DataFrame]:
    ct = pd.crosstab([df["sample_id"], df["day"]], df["author_lineage"]).reindex(columns=LINEAGES, fill_value=0)
    rows = []
    for (s, d), r in ct.iterrows():
        n_animal = int(r.sum())
        comp_tot = {"CD45neg": int(r[["Endothelium", "Mesenchyme", "Epithelium"]].sum()),
                    "CD45pos": int(r[["Myeloid", "Lymphoid"]].sum())}
        for lin in LINEAGES:
            comp = COMPARTMENT[lin]
            rows.append({"sample_id": s, "day": int(d), "phase": PHASE[int(d)], "lineage": lin,
                         "compartment": comp, "n_cells": int(r[lin]), "n_cells_animal": n_animal,
                         "n_cells_compartment": comp_tot[comp],
                         "pct_of_all_cells": round(100 * r[lin] / n_animal, 3),
                         "pct_within_compartment": (round(100 * r[lin] / comp_tot[comp], 3)
                                                    if comp_tot[comp] else np.nan)})
    long = pd.DataFrame(rows).sort_values(["day", "sample_id", "lineage"])
    p = TAB / "lineage_composition_per_animal.csv"
    long.to_csv(p, index=False)
    rec.add_output(p)
    med = (long.groupby(["day", "lineage"])
               .agg(n_animals=("sample_id", "nunique"),
                    median_pct_of_all_cells=("pct_of_all_cells", "median"),
                    median_pct_within_compartment=("pct_within_compartment", "median"))
               .reset_index().round(2))
    med["phase"] = med["day"].map(PHASE)
    p = TAB / "lineage_composition_median_by_dpi.csv"
    med.to_csv(p, index=False)
    rec.add_output(p)

    fig, axes = plt.subplots(2, 5, figsize=(17, 6.8), sharex=True)
    for j, lin in enumerate(LINEAGES):
        col = LINEAGE_COLOURS[lin]
        sub = long[long["lineage"] == lin]
        ax = axes[0, j]
        strip_by_day(ax, sub, "pct_of_all_cells", DAYS, col)
        ax.set_title(f"{lin}\n% of all cells (sort-engineered)", fontsize=9)
        ax.set_ylim(bottom=0)
        ax = axes[1, j]
        strip_by_day(ax, sub, "pct_within_compartment", DAYS, col)
        ax.set_title(f"% of {'CD45-negative' if COMPARTMENT[lin] == 'CD45neg' else 'CD45-positive'} cells",
                     fontsize=9)
        ax.set_xlabel("days post infection")
        ax.set_ylim(bottom=0)
    axes[0, 0].set_ylabel("% per animal")
    axes[1, 0].set_ylabel("% per animal")
    fig.suptitle("Lineage composition per animal (one point = one animal; line = median). "
                 "Top row is set by the MACS recombination and is descriptive only.", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for p in save_fig(fig, FIG, "lineage_composition_by_dpi"):
        rec.add_output(p)
    return long, med


def proliferation(df: pd.DataFrame, rec: RunRecord) -> dict:
    scored = df[df["trace_call"].isin(["Traced", "Untraced"])]
    nd = df.groupby(["sample_id", "day"])["trace_call"].apply(lambda s: round(100 * (s == "Not_detected").mean(), 2))
    rec.set("pct_not_detected_per_animal", {f"{s} ({d} dpi)": float(v) for (s, d), v in nd.items()})
    rec.set("pct_not_detected_overall", round(100 * (df["trace_call"] == "Not_detected").mean(), 2))

    g = (scored.groupby(["sample_id", "day", "tam", "author_lineage"])
               .agg(n_scored=("trace_call", "size"),
                    n_traced=("trace_call", lambda s: int((s == "Traced").sum())))
               .reset_index().rename(columns={"author_lineage": "lineage"}))
    g["pct_traced"] = (100 * g["n_traced"] / g["n_scored"]).round(3)
    g["phase"] = g["day"].map(PHASE)
    g["window"] = np.where(g["tam"] > 0, "tamoxifen " + g["tam"].astype(str) + " dpi",
                           "tamoxifen " + (-g["tam"]).astype(str) + " d before harvest (uninjured)")
    imm = g["tam"].map(IMMEDIATE_WINDOWS)
    g["series"] = np.where(g["day"] == 0, "baseline",
                           np.where(imm == g["day"], "immediate_window", "later_harvest"))
    g = g.sort_values(["day", "tam", "sample_id", "lineage"])
    p = TAB / "traced_fraction_per_animal.csv"
    g.to_csv(p, index=False)
    rec.add_output(p)

    c = (df.groupby(["sample_id", "day", "author_lineage"])
           .agg(n_cells=("cell_cycle_phase", "size"),
                n_cycling=("cell_cycle_phase", lambda s: int(s.isin(["S", "G2M"]).sum())))
           .reset_index().rename(columns={"author_lineage": "lineage"}))
    c["pct_cycling"] = (100 * c["n_cycling"] / c["n_cells"]).round(3)
    c["phase"] = c["day"].map(PHASE)
    c = c.sort_values(["day", "sample_id", "lineage"])
    p = TAB / "cycling_fraction_per_animal.csv"
    c.to_csv(p, index=False)
    rec.add_output(p)

    # peaks by the frozen rule
    imm_med = (g[g["series"].isin(["immediate_window", "baseline"])]
               .groupby(["day", "lineage"])["pct_traced"].median().unstack("lineage"))
    cyc_med = c.groupby(["day", "lineage"])["pct_cycling"].median().unstack("lineage")
    rows = []
    for lin in LINEAGES:
        s = imm_med.loc[IMMEDIATE_DAYS, lin]
        peak = int(s.idxmax())
        cs = cyc_med.loc[IMMEDIATE_DAYS, lin]
        cpeak = int(cs.idxmax())
        row = {"lineage": lin,
               "median_pct_traced_0dpi": round(float(imm_med.loc[0, lin]), 2)}
        for d in IMMEDIATE_DAYS:
            row[f"median_pct_traced_{d}dpi"] = round(float(s[d]), 2)
        row.update({"observed_trace_peak_day": peak,
                    "observed_cycling_peak_day": cpeak,
                    "paper_expected_peak": PAPER_PEAK[lin],
                    "trace_agrees_with_paper": peak in PAPER_PEAK_DAYS[lin],
                    "cycling_agrees_with_paper": cpeak in PAPER_PEAK_DAYS[lin]})
        rows.append(row)
    peaks = pd.DataFrame(rows)
    p = TAB / "proliferation_peak_by_lineage.csv"
    peaks.to_csv(p, index=False)
    rec.add_output(p)
    rec.set("proliferation_peaks", peaks.set_index("lineage").to_dict(orient="index"))
    rec.set("n_lineages_trace_peak_agrees_with_paper", int(peaks["trace_agrees_with_paper"].sum()))
    rec.set("n_lineages_cycling_peak_agrees_with_paper", int(peaks["cycling_agrees_with_paper"].sum()))

    # figure 1: immediate windows (trace) and cycling by day
    imm_days = [0] + IMMEDIATE_DAYS
    fig, axes = plt.subplots(2, 5, figsize=(17, 6.8))
    for j, lin in enumerate(LINEAGES):
        col = LINEAGE_COLOURS[lin]
        ax = axes[0, j]
        sub = g[(g["lineage"] == lin) & (g["series"].isin(["immediate_window", "baseline"]))]
        strip_by_day(ax, sub, "pct_traced", imm_days, col)
        ax.set_title(f"{lin}\nKi67-traced, immediate window\n(expected peak {PAPER_PEAK[lin]})", fontsize=9)
        ax.set_xlabel("harvest day (tamoxifen 4 d earlier; 0 = uninjured)")
        ax.set_ylim(bottom=0)
        ax = axes[1, j]
        sub = c[c["lineage"] == lin]
        strip_by_day(ax, sub, "pct_cycling", DAYS, col)
        ax.set_title("cells in S or G2M", fontsize=9)
        ax.set_xlabel("days post infection")
        ax.set_ylim(bottom=0)
    axes[0, 0].set_ylabel("% traced of reporter-scored cells")
    axes[1, 0].set_ylabel("% cycling per animal")
    fig.suptitle("Proliferation per lineage (one point = one animal; line = median)", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for p in save_fig(fig, FIG, "proliferation_by_lineage"):
        rec.add_output(p)

    # figure 2: the common 42 dpi harvest by tamoxifen window, 90 dpi as hollow points
    windows = [2, 7, 14, 21]
    fig, axes = plt.subplots(1, 5, figsize=(17, 3.8), sharey=False)
    for j, lin in enumerate(LINEAGES):
        col = LINEAGE_COLOURS[lin]
        ax = axes[j]
        sub42 = g[(g["lineage"] == lin) & (g["day"] == 42)]
        strip_by_day(ax, sub42, "pct_traced", windows, col, day_col="tam")
        sub90 = g[(g["lineage"] == lin) & (g["day"] == 90)]
        strip_by_day(ax, sub90, "pct_traced", windows, col, day_col="tam", median=False, hollow=True)
        ax.set_title(f"{lin}\nfilled: 42 dpi harvest (2 animals per window)\nhollow: 90 dpi (1 animal)",
                     fontsize=9)
        ax.set_xlabel("tamoxifen window start (dpi)")
        ax.set_ylim(bottom=0)
    axes[0].set_ylabel("% traced of reporter-scored cells")
    fig.suptitle("Ki67 trace read at a common harvest: which window labelled each lineage", y=1.02)
    fig.tight_layout()
    for p in save_fig(fig, FIG, "traced_fraction_at_common_harvest_by_window"):
        rec.add_output(p)

    return {"traced": g, "cycling": c, "peaks": peaks, "imm_med": imm_med, "cyc_med": cyc_med}


def write_readme(df: pd.DataFrame, med: pd.DataFrame, prol: dict, rec: RunRecord) -> None:
    per_day = (df.groupby("day").agg(n_animals=("sample_id", "nunique"), n_cells=("sample_id", "size"))
                 .reindex(DAYS))
    per_day.insert(0, "phase", [PHASE[d] for d in per_day.index])
    per_day.index.name = "dpi"
    peaks = prol["peaks"].copy()
    show = peaks[["lineage", "median_pct_traced_0dpi"] + [f"median_pct_traced_{d}dpi" for d in IMMEDIATE_DAYS]
                 + ["observed_trace_peak_day", "observed_cycling_peak_day", "paper_expected_peak",
                    "trace_agrees_with_paper", "cycling_agrees_with_paper"]]
    comp_med = (med.pivot(index="day", columns="lineage", values="median_pct_within_compartment")
                   .reindex(DAYS)[LINEAGES])
    comp_med.index.name = "dpi"
    cyc = prol["cyc_med"].reindex(DAYS)[LINEAGES].round(2)
    cyc.index.name = "dpi"
    n_agree = int(peaks["trace_agrees_with_paper"].sum())
    n_agree_c = int(peaks["cycling_agrees_with_paper"].sum())
    lines = [
        "# Phase-wise view of the Ki67 atlas (generated)",
        "",
        "Generated by `analysis/scripts/10_phase_timecourse.py` from the tracked table",
        "`../tables/cell_metadata.csv`; the frozen rules, inputs and results are in",
        "[`run_record.json`](run_record.json). Owner retain/reject review pending.",
        "",
        "**What is being reproduced.** Niethamer et al. 2025 report that proliferation after",
        "influenza injury is asynchronous across compartments: immune cells proliferate first",
        "(2 to 6 dpi), epithelium and mesenchyme second (7 to 15 dpi), endothelium last",
        "(14 to 22 dpi), and every lineage is back at baseline by 21 to 28 dpi (their Figures 1",
        "and 2, and the three-phase model of Figure 2F). The readouts here use the deposited",
        "lineage labels as given, so this is a descriptive reproduction, not a blind validation.",
        "",
        "## Cohort",
        "",
        f"{rec.record['results']['n_cells_cohort']:,} cells from {rec.record['results']['n_animals_cohort']} animals "
        "(the 25-sample Ki67 atlas; the 8 pre-labelled tracing samples have no day post infection and are excluded).",
        "",
        df_to_markdown(per_day),
        "",
        "## 1. Atlas UMAP by day post infection",
        "",
        "![Atlas UMAP by dpi](figures/UMAP_atlas_by_dpi_lineage.png)",
        "",
        f"Shared atlas coordinates; every panel draws {rec.record['results']['umap_cells_drawn_per_panel']:,} cells "
        "(the smallest per-day count, seed 0) over all cohort cells in grey, so that density is comparable",
        "between days. Panel area and cluster distance carry no meaning (see `docs/UMAP_AND_FIGURES.md`).",
        "",
        "## 2. Lineage composition per animal",
        "",
        "![Lineage composition by dpi](figures/lineage_composition_by_dpi.png)",
        "",
        "Median per-animal percent of each lineage **within its compartment** (CD45-negative:",
        "Endothelium, Mesenchyme, Epithelium; CD45-positive: Myeloid, Lymphoid). The top row of the",
        "figure, percent of all cells, is set by the MACS recombination (10 to 15 percent CD45-positive)",
        "and is descriptive only. Per-animal values: `tables/lineage_composition_per_animal.csv`.",
        "",
        df_to_markdown(comp_med),
        "",
        "## 3. Proliferation per lineage",
        "",
        "![Proliferation by lineage](figures/proliferation_by_lineage.png)",
        "",
        "Top row: Ki67 trace in each cohort's immediate window (tamoxifen at 2, 7, 14 or 21 dpi,",
        "harvest four days later). Cells whose reporter was not detected are excluded as dropout,",
        "the rule of `07_lineage_tracing_cohort.py`; the not-detected share is "
        f"{rec.record['results']['pct_not_detected_overall']}% of cohort cells overall and is listed per animal in the run record.",
        "Bottom row: the deposited cell-cycle call (S or G2M) per day. Per-animal values:",
        "`tables/traced_fraction_per_animal.csv` and `tables/cycling_fraction_per_animal.csv`.",
        "",
        "**Peak window by the frozen rule** (highest median per-animal value over the 6, 11, 19 and 25 dpi",
        "harvests; two animals per day, so this is a ranking, not a test):",
        "",
        df_to_markdown(show, index=False),
        "",
        f"Trace peaks agree with the paper's expected window for {n_agree} of 5 lineages; cycling peaks",
        f"agree for {n_agree_c} of 5.",
        "",
        "Median percent of cells in S or G2M per lineage and day:",
        "",
        df_to_markdown(cyc),
        "",
        "![Trace at a common harvest](figures/traced_fraction_at_common_harvest_by_window.png)",
        "",
        "At the common 42 dpi harvest (two animals per tamoxifen window) the trace reads which window",
        "labelled each lineage, including the labelled cells' progeny; 90 dpi animals (one per window)",
        "are hollow points.",
        "",
        "## Caveats",
        "",
        "- The unit is the animal. Active-repair days carry two animals; 42 dpi carries eight because",
        "  four tamoxifen cohorts converge on it. No P values are computed.",
        "- Immune versus non-immune proportions describe the sort, not the lung. Only within-compartment",
        "  fractions are interpretable.",
        "- The trace marks cells that proliferated during the tamoxifen window and their progeny; a",
        "  high traced fraction at a later harvest can reflect expansion after the window as well as",
        "  proliferation within it.",
        "- Cell-cycle phase is the authors' per-cell Seurat CellCycleScoring call, read from the deposited",
        "  metadata and never recalculated (`docs/PIPELINE_AS_RUN.md`). It assigns S or G2M to any cell",
        "  with a positive phase score, so it over-calls cycling (most lymphocytes score as cycling) and",
        "  is shown only as a trace-independent cross-check, not as a measurement of proliferation.",
        "- The two baseline animals carry different tamoxifen windows (3 and 42 days before harvest) and",
        "  are not comparable with each other on the trace readout.",
        "",
    ]
    p = OUT / "README.md"
    p.write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(p)


def main() -> int:
    for d in (OUT, FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    rec = RunRecord(OUT / "run_record.json", "phase timecourse of the Ki67 atlas from tracked metadata",
                    RULES, notes="rules frozen before the metadata table was opened; deposited lineage "
                                 "labels used descriptively; unit is the animal")
    rec.add_output(OUT / "run_record.json")
    df = load_cohort(rec)
    umap_by_day(df, rec)
    _, med = composition(df, rec)
    prol = proliferation(df, rec)
    write_readme(df, med, prol, rec)
    rec.finish()
    print()
    print(prol["peaks"].to_string(index=False))
    print(f"\nwrote {OUT.relative_to(ANALYSIS.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
