#!/usr/bin/env python
"""Where does the rebuilt alveolar macrophage pool come from? A trace-only read.

Niethamer et al. 2025 (Figure 3) report that alveolar macrophages (aMAC) are
lost at 6 dpi and reconstituted over the following three weeks from two
sources: resident aMACs that resemble the homeostatic population, and cells
arriving along a trajectory from the inflammatory monocytes (iMON) present at
6 dpi. This script reads the Ki67 trace instead of fitting a trajectory. The
atlas carries four tamoxifen windows (2 to 3, 7 to 8, 14 to 15, 21 to 22 dpi)
and a common 42 dpi harvest with two animals per window, so the traced
fraction of aMAC-labelled cells at 42 dpi says which window's proliferating
cells (and their progeny) populate the rebuilt pool.

Two things the trace cannot separate on its own, and how the rules handle
them. First, Ki67-CreERT2 also labels bone-marrow progenitors that divide
during the window, and their progeny inherit the label; the short-lived,
marrow-derived classical monocytes (cMON) and neutrophils of the same animal
are read as a reference for that inherited labelling. Second, resident aMAC
proliferation and iMON-derived reconstitution are both traceable in the early
windows; the blind subclusters from 11_myeloid_focus.py are read for a
within-aMAC split in traced fraction, which is reported descriptively.

Input is the tracked per-cell table written by 11_myeloid_focus.py. Rules are
frozen in the run record before the table is opened. The unit is the animal;
medians are reported and no P values are computed.
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

MYELOID = ANALYSIS / "GSE262927" / "myeloid_focus"
META = MYELOID / "tables" / "myeloid_cell_metadata.csv"
OUT = MYELOID / "amac_origin"
FIG = OUT / "figures"
TAB = OUT / "tables"

LABELS = ["aMAC", "iMON", "iMAC", "cMON", "Neutrophil"]
MARROW_REFERENCE = ["cMON", "Neutrophil"]
LABEL_COLOURS = {"aMAC": "#1f77b4", "iMAC": "#17becf", "iMON": "#2ca02c", "cMON": "#98df8a",
                 "Neutrophil": "#9467bd"}
WINDOWS = [2, 7, 14, 21]
WINDOW_NAME = {2: "2 to 3 dpi", 7: "7 to 8 dpi", 14: "14 to 15 dpi", 21: "21 to 22 dpi"}
IMMEDIATE = {2: 6, 7: 11, 14: 19, 21: 25}
SUBCLUSTER_KEY = "sub_r0.5"
MIN_CELLS = 30
MARROW_MARGIN = 10.0   # percentage points
SPLIT_MARGIN = 20.0    # percentage points

RULES = {
    "input": "analysis/GSE262927/myeloid_focus/tables/myeloid_cell_metadata.csv (tracked; written by "
             "11_myeloid_focus.py from atlas clusters 5, 17, 24 of the 25-sample Ki67 atlas)",
    "cells": "cells with a deposited myeloid label in the frozen list and trace_call in {Traced, Untraced}; "
             "Not_detected is reporter dropout and is excluded (rule of 07_lineage_tracing_cohort.py)",
    "labels_read": LABELS,
    "marrow_reference_labels": MARROW_REFERENCE,
    "windows": {str(k): v for k, v in WINDOW_NAME.items()},
    "harvests": {"common": "42 dpi, two animals per window", "later": "90 dpi (one animal per window) and "
                 "366 dpi (7 dpi window: one animal; 14 dpi window: two)",
                 "immediate": {f"window {k} dpi": f"harvest {v} dpi" for k, v in IMMEDIATE.items()}},
    "pct_traced": "100 * Traced / (Traced + Untraced) per animal and label; reported only when the animal "
                  f"has at least {MIN_CELLS} scored cells of that label, otherwise left empty",
    "question_1_window_contribution": "the window whose two 42 dpi animals give the highest median aMAC "
                                      "traced fraction is the main contributor; paper expectation: the 2 to 3 "
                                      "or 7 to 8 dpi window (immune proliferation peaks at 2 to 7 dpi and the "
                                      "aMAC pool is rebuilt over the following three weeks)",
    "question_2_marrow_inheritance": "per 42 dpi animal, aMAC traced fraction minus the same animal's cMON "
                                     "and neutrophil traced fractions; a window's aMAC labelling 'exceeds the "
                                     f"marrow reference' if aMAC exceeds both references by at least {MARROW_MARGIN} "
                                     "percentage points in both animals of that window",
    "question_3_within_amac_split": f"per 42 dpi animal, traced fraction of aMAC-labelled cells per blind "
                                    f"subcluster ({SUBCLUSTER_KEY} from 11_myeloid_focus.py) with at least "
                                    f"{MIN_CELLS} scored cells; a subcluster pair whose traced fractions differ "
                                    f"by at least {SPLIT_MARGIN} points in the same direction in both animals of a "
                                    "window is reported as a candidate two-source signature. Descriptive only",
    "unit": "animal; medians; no P values (two animals per window at 42 dpi, one at 90 dpi)",
    "seed": RANDOM_SEED,
}


def load(rec: RunRecord) -> pd.DataFrame:
    if not META.exists():
        raise SystemExit(f"input table missing: {META}; run 11_myeloid_focus.py first")
    rec.add_input(META)
    df = pd.read_csv(META, dtype=str)
    need = {"sample_id", "day", "tamoxifen_start_day", "trace_call", "label", SUBCLUSTER_KEY}
    missing = need - set(df.columns)
    if missing:
        raise SystemExit(f"input table lacks columns {sorted(missing)}")
    df["day"] = pd.to_numeric(df["day"], errors="coerce").astype(int)
    df["tam"] = pd.to_numeric(df["tamoxifen_start_day"], errors="coerce").astype(int)
    n_all = len(df)
    df = df[df["label"].isin(LABELS) & df["trace_call"].isin(["Traced", "Untraced"])].copy()
    df["traced"] = df["trace_call"] == "Traced"
    rec.set("n_cells_in_table", n_all)
    rec.set("n_cells_scored_in_frozen_labels", int(len(df)))
    rec.set("n_animals", int(df["sample_id"].nunique()))
    return df


def per_animal(df: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    g = (df.groupby(group_cols).agg(n_scored=("traced", "size"), n_traced=("traced", "sum")).reset_index())
    g["n_traced"] = g["n_traced"].astype(int)
    g["pct_traced"] = np.where(g["n_scored"] >= MIN_CELLS, (100 * g["n_traced"] / g["n_scored"]).round(2), np.nan)
    return g


def series_of(day: int, tam: int) -> str:
    if day == 0:
        return "baseline"
    if IMMEDIATE.get(tam) == day:
        return "immediate"
    if day == 42:
        return "common_42"
    return f"later_{day}"


def main() -> int:
    for d in (OUT, FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    rec = RunRecord(OUT / "run_record.json", "origin of the rebuilt alveolar macrophage pool, read from the "
                    "Ki67 trace by tamoxifen window", RULES,
                    notes="rules frozen before the table was opened; deposited labels used descriptively")
    rec.add_output(OUT / "run_record.json")
    df = load(rec)

    # per animal and label
    lab = per_animal(df, ["sample_id", "day", "tam", "label"])
    lab["window"] = lab["tam"].map(lambda t: WINDOW_NAME.get(t, f"{t} d before harvest (uninjured)"))
    lab["series"] = [series_of(d, t) for d, t in zip(lab["day"], lab["tam"])]
    lab = lab.sort_values(["label", "day", "tam", "sample_id"])
    p = TAB / "traced_fraction_per_animal_label.csv"
    lab.to_csv(p, index=False)
    rec.add_output(p)

    # question 1: window contribution at 42 dpi
    c42 = lab[(lab["day"] == 42)]
    med = (c42.groupby(["label", "tam"])["pct_traced"].median().unstack("tam").reindex(columns=WINDOWS))
    med.index.name = "label"
    later = lab[lab["day"].isin([90, 366])].pivot_table(index="label", columns=["day", "tam"], values="pct_traced", aggfunc="median")
    summary = med.copy()
    summary.columns = [f"median_pct_traced_42dpi_window_{w}" for w in WINDOWS]
    for (d, t) in sorted(later.columns):
        summary[f"pct_traced_{d}dpi_window_{t}"] = later[(d, t)]
    summary = summary.round(2)
    p = TAB / "window_contribution_summary.csv"
    summary.to_csv(p)
    rec.add_output(p)
    amac_row = med.loc["aMAC"]
    if amac_row.isna().any():
        raise SystemExit("aMAC traced fraction missing for a window at 42 dpi; refusing to rank")
    main_window = int(amac_row.idxmax())
    q1 = {"median_pct_traced_42dpi_by_window": {str(w): float(amac_row[w]) for w in WINDOWS},
          "main_contributing_window": main_window,
          "agrees_with_paper": main_window in (2, 7)}
    rec.set("question_1", q1)

    # question 2: marrow reference per animal
    rows = []
    for (s, t), g in c42.groupby(["sample_id", "tam"]):
        vals = g.set_index("label")["pct_traced"]
        a = vals.get("aMAC", np.nan)
        refs = {r: vals.get(r, np.nan) for r in MARROW_REFERENCE}
        rows.append({"sample_id": s, "window_start_dpi": t, "pct_traced_aMAC": a,
                     **{f"pct_traced_{r}": v for r, v in refs.items()},
                     **{f"aMAC_minus_{r}": (round(a - v, 2) if not (np.isnan(a) or np.isnan(v)) else np.nan) for r, v in refs.items()},
                     "exceeds_both_references_by_margin": bool(all((not np.isnan(a)) and (not np.isnan(v)) and (a - v >= MARROW_MARGIN) for v in refs.values()))})
    marrow = pd.DataFrame(rows).sort_values(["window_start_dpi", "sample_id"])
    p = TAB / "marrow_reference_check.csv"
    marrow.to_csv(p, index=False)
    rec.add_output(p)
    q2 = {}
    for t in WINDOWS:
        sub = marrow[marrow["window_start_dpi"] == t]
        q2[str(t)] = {"n_animals": int(len(sub)),
                      "both_animals_exceed": bool(len(sub) == 2 and sub["exceeds_both_references_by_margin"].all())}
    rec.set("question_2", q2)

    # question 3: within-aMAC split by blind subcluster at 42 dpi
    a42 = df[(df["day"] == 42) & (df["label"] == "aMAC")]
    subc = per_animal(a42, ["sample_id", "tam", SUBCLUSTER_KEY]).rename(columns={SUBCLUSTER_KEY: "subcluster"})
    subc = subc.sort_values(["tam", "sample_id", "subcluster"])
    p = TAB / "amac_subcluster_traced_fraction_42dpi.csv"
    subc.to_csv(p, index=False)
    rec.add_output(p)
    q3 = {}
    for t in WINDOWS:
        animals = sorted(subc.loc[subc["tam"] == t, "sample_id"].unique())
        pivot = subc[subc["tam"] == t].pivot(index="subcluster", columns="sample_id", values="pct_traced")
        pivot = pivot.dropna(how="any")
        found = []
        subs = list(pivot.index)
        for i in range(len(subs)):
            for j in range(len(subs)):
                if i == j:
                    continue
                diffs = pivot.loc[subs[i]] - pivot.loc[subs[j]]
                if len(animals) == 2 and (diffs >= SPLIT_MARGIN).all():
                    found.append({"more_traced": subs[i], "less_traced": subs[j],
                                  "differences": {a: round(float(d), 2) for a, d in diffs.items()}})
        q3[str(t)] = {"animals": animals, "subclusters_with_enough_cells_in_both_animals": subs,
                      "candidate_two_source_pairs": found}
    rec.set("question_3", q3)

    # ---------------- figures ----------------
    fig, axes = plt.subplots(1, 5, figsize=(19, 4.2))
    for ax, l in zip(axes, LABELS):
        col = LABEL_COLOURS[l]
        s42 = lab[(lab["label"] == l) & (lab["day"] == 42)]
        strip_by_day(ax, s42, "pct_traced", WINDOWS, col, day_col="tam")
        s90 = lab[(lab["label"] == l) & (lab["day"] == 90)]
        strip_by_day(ax, s90, "pct_traced", WINDOWS, col, day_col="tam", median=False, hollow=True)
        s366 = lab[(lab["label"] == l) & (lab["day"] == 366)].dropna(subset=["pct_traced"])
        if len(s366):
            pos = {d: i for i, d in enumerate(WINDOWS)}
            ax.scatter(s366["tam"].map(pos) + 0.25, s366["pct_traced"], marker="^", s=30, color=col, alpha=0.9, zorder=3)
        if lab[(lab["label"] == l) & (lab["day"].isin([42, 90, 366]))]["pct_traced"].isna().all():
            ax.text(0.5, 0.5, f"fewer than {MIN_CELLS} scored cells\nper animal at these harvests",
                    ha="center", va="center", transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_title(l, fontsize=10)
        ax.set_xlabel("tamoxifen window start (dpi)")
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("% traced of reporter-scored cells")
    fig.suptitle("Ki67 trace at late harvests: filled = 42 dpi (two animals per window), hollow = 90 dpi (one), "
                 "triangle = 366 dpi. cMON and neutrophils are the marrow-inheritance reference.", y=1.03)
    fig.tight_layout()
    for p in save_fig(fig, FIG, "trace_by_window_late_harvests"):
        rec.add_output(p)

    fig, axes = plt.subplots(1, 5, figsize=(19, 4.2))
    imm_days = [0] + [IMMEDIATE[w] for w in WINDOWS]
    for ax, l in zip(axes, LABELS):
        col = LABEL_COLOURS[l]
        s = lab[(lab["label"] == l) & (lab["series"].isin(["immediate", "baseline"]))]
        strip_by_day(ax, s, "pct_traced", imm_days, col)
        if s["pct_traced"].isna().all():
            ax.text(0.5, 0.5, f"fewer than {MIN_CELLS} scored cells\nper animal", ha="center", va="center",
                    transform=ax.transAxes, fontsize=8, color="grey")
        ax.set_title(l, fontsize=10)
        ax.set_xlabel("harvest day (0 = uninjured)")
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("% traced of reporter-scored cells")
    fig.suptitle("Ki67 trace in each cohort's immediate window (harvest four days after tamoxifen; "
                 "one point = one animal; line = median)", y=1.03)
    fig.tight_layout()
    for p in save_fig(fig, FIG, "trace_immediate_windows"):
        rec.add_output(p)

    subs_all = sorted(subc["subcluster"].dropna().unique(), key=int)
    cmap = plt.get_cmap("tab20")
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.2), sharey=True)
    for ax, t in zip(axes, WINDOWS):
        s = subc[(subc["tam"] == t)].dropna(subset=["pct_traced"])
        animals = sorted(s["sample_id"].unique())
        for k, a in enumerate(animals):
            sa = s[s["sample_id"] == a]
            x = [subs_all.index(v) + (k - 0.5) * 0.3 for v in sa["subcluster"]]
            ax.bar(x, sa["pct_traced"], width=0.28, color=[cmap(int(v) % 20) for v in sa["subcluster"]],
                   edgecolor="black" if k else "none", linewidth=0.6, label=a)
        ax.set_xticks(range(len(subs_all)))
        ax.set_xticklabels([f"sub {v}" for v in subs_all], rotation=45, fontsize=8)
        ax.set_title(f"window {WINDOW_NAME[t]}\n42 dpi animals: {', '.join(animals)}", fontsize=9)
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("% traced of aMAC-labelled cells")
    fig.suptitle(f"Within-aMAC split at 42 dpi by blind subcluster (bars: one per animal; "
                 f"only subclusters with >= {MIN_CELLS} scored aMAC cells in that animal)", y=1.03)
    fig.tight_layout()
    for p in save_fig(fig, FIG, "amac_subcluster_trace_42dpi"):
        rec.add_output(p)

    # ---------------- README ----------------
    marrow_show = marrow.copy()
    n_eval = int(marrow[["pct_traced_aMAC"] + [f"pct_traced_{r}" for r in MARROW_REFERENCE]].notna().all(axis=1).sum())
    rec.set("question_2_animals_with_both_references_evaluable", n_eval)
    q2_note = (f"Both references reached the {MIN_CELLS}-cell floor in {n_eval} of {len(marrow)} animals at 42 dpi "
               "(classical monocytes are sparse by then), so the rule is Not established for every window with "
               "fewer than two evaluable animals; the per-animal differences above are descriptive.")
    q3_lines = []
    for t in WINDOWS:
        pairs = q3[str(t)]["candidate_two_source_pairs"]
        if pairs:
            for pr in pairs:
                q3_lines.append(f"- window {WINDOW_NAME[t]}: subcluster {pr['more_traced']} more traced than "
                                f"subcluster {pr['less_traced']} in both animals (differences "
                                + ", ".join(f"{a}: {d}" for a, d in pr["differences"].items()) + " points).")
        else:
            q3_lines.append(f"- window {WINDOW_NAME[t]}: no subcluster pair differs by at least {SPLIT_MARGIN} points "
                            f"in both animals (subclusters with enough cells in both animals: "
                            f"{', '.join(q3[str(t)]['subclusters_with_enough_cells_in_both_animals']) or 'none'}).")
        # per-animal values where at least two subclusters are evaluable, for the reader
        for a, g in subc[(subc["tam"] == t)].dropna(subset=["pct_traced"]).groupby("sample_id"):
            if len(g) >= 2:
                q3_lines.append("  - " + a + ": " + ", ".join(f"sub {r.subcluster} {r.pct_traced}% (n={r.n_scored})"
                                                              for r in g.itertuples()) + ".")
    rec.set("question_3_per_animal_values", {t: subc[(subc["tam"] == t)].dropna(subset=["pct_traced"])
                                             .to_dict(orient="records") for t in WINDOWS})
    lines = [
        "# Origin of the rebuilt alveolar macrophage pool, read from the Ki67 trace (generated)",
        "",
        "Generated by `analysis/scripts/12_amac_trace_by_window.py` from the tracked table",
        "`../tables/myeloid_cell_metadata.csv`; frozen rules, inputs and results in",
        "[`run_record.json`](run_record.json). Owner retain/reject review pending.",
        "",
        "**Question.** Niethamer et al. 2025 (Figure 3) propose two sources for the aMAC pool rebuilt after",
        "its loss at 6 dpi: resident aMACs and an iMON-derived trajectory. Here the Ki67 trace is read",
        "by tamoxifen window at the common 42 dpi harvest (two animals per window), with classical",
        "monocytes and neutrophils of the same animal as the reference for labelling inherited from",
        "bone-marrow progenitors that divided during the window.",
        "",
        f"Cells: {rec.record['results']['n_cells_scored_in_frozen_labels']:,} reporter-scored cells with a deposited label "
        f"in {', '.join(LABELS)} from {rec.record['results']['n_animals']} animals.",
        "",
        "## 1. Which window labelled the 42 dpi aMAC pool",
        "",
        "![Trace by window at late harvests](figures/trace_by_window_late_harvests.png)",
        "",
        "Median per-animal traced fraction at 42 dpi by window (two animals each), with the single 90 dpi",
        "and 366 dpi animals alongside (`tables/window_contribution_summary.csv`; per animal in",
        "`tables/traced_fraction_per_animal_label.csv`):",
        "",
        df_to_markdown(summary),
        "",
        f"By the frozen rule the main contributing window for aMAC is **{WINDOW_NAME[main_window]}** "
        f"(agrees with the paper's expectation of an early window: {q1['agrees_with_paper']}).",
        "",
        "## 2. Is the aMAC labelling explained by marrow inheritance",
        "",
        df_to_markdown(marrow_show, index=False),
        "",
        "Windows where aMAC exceeds both marrow references by at least "
        f"{MARROW_MARGIN} points in both animals: "
        + (", ".join(WINDOW_NAME[t] for t in WINDOWS if q2[str(t)]["both_animals_exceed"]) or "none") + ". " + q2_note,
        "",
        "## 3. Within-aMAC split by blind subcluster at 42 dpi",
        "",
        "![aMAC subcluster trace](figures/amac_subcluster_trace_42dpi.png)",
        "",
        *q3_lines,
        "",
        "Per-animal values: `tables/amac_subcluster_traced_fraction_42dpi.csv`.",
        "",
        "## 4. The same populations in their immediate windows",
        "",
        "![Immediate windows](figures/trace_immediate_windows.png)",
        "",
        "## Caveats",
        "",
        "- Two animals per window at 42 dpi and one at 90 dpi: rankings and margins, no tests.",
        "- The trace marks cells that divided during the window and all their progeny. Labelled marrow",
        "  progenitors keep producing labelled monocytes after the window, which is why the marrow",
        "  reference is read; the reference is itself imperfect because monocyte and neutrophil half-lives",
        "  differ from the interval between window and harvest.",
        "- Resident aMAC proliferation and iMON-derived reconstitution are both traceable in the early",
        "  windows; the within-aMAC split by subcluster is a descriptive signature, not a lineage assignment.",
        "- Deposited labels are used as given; the subclusters come from an embedding that held them out.",
        "",
    ]
    p = OUT / "README.md"
    p.write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(p)
    rec.finish()
    print(summary.to_string())
    print("\nquestion 1:", q1)
    print("question 2:", q2)
    print("question 3:", {k: len(v["candidate_two_source_pairs"]) for k, v in q3.items()})
    print(f"wrote {OUT.relative_to(ANALYSIS.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
