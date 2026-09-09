#!/usr/bin/env python
"""Trial S4: what is mouse cluster 23, the cluster S1 flagged as animal-private and unlabelled?

Portfolio purpose. Curation hygiene: before any cluster of the mouse atlas
is cited as a population, a cluster that is 99.5% unlabelled by the authors
and dominated by one animal must be explained. This is the kind of check a
donor-aware, reproducibility-first portfolio (Wagner fit) has to show.

Frozen rules (written to the run record before any data is read). Note:
the cluster-level QC summary (median counts 1,868, 78% from one sample) had
already been seen when these rules were set; the per-cell data had not.

* Comparators: (A) all cells outside cluster 23; (B) a random sample of the
  same size from outside cluster 23 (seed 0), used for the expression-based
  metrics; (C) cells of the dominant sample outside cluster 23.
* Low-count call: median total counts of cluster 23 < 0.5 x the atlas
  median AND median genes < 0.5 x the atlas median.
* Ambient-like call: the fraction of cells with >= 3 of 5 lineage markers
  detected (Sftpc, Scgb1a1, Ptprc, Pecam1, Col1a1; raw counts > 0) is
  >= 2 x that fraction in comparator B.
* Cohort split: fraction of cluster-23 cells from the 8 lineage-tracing
  samples (which carry no author labels at all) versus from the 25
  annotated samples; within the annotated cohort, the labelled fraction of
  cluster 23 is compared with the cohort-wide labelled fraction.
* Verdict: "low-count, ambient-like" if both calls are true; "low-count
  only", "ambient-like only", or "unexplained" otherwise. Descriptive only.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from trial_utils import REPO, RunRecord, df_to_markdown, read_csr_rows  # noqa: E402

OUT = HERE / "s4_cluster23_qc"
OUT.mkdir(exist_ok=True)
OBJ = REPO / "analysis" / "GSE262927" / "processed" / "final_clustered.h5ad"
CLUSTER = "23"
LINEAGE = ["Sftpc", "Scgb1a1", "Ptprc", "Pecam1", "Col1a1"]
RULES = {
    "cluster": CLUSTER,
    "low_count_rule": "median total_counts < 0.5 * atlas median AND median n_genes < 0.5 * atlas median",
    "ambient_rule": "fraction of cells with >= 3 of 5 lineage markers detected >= 2 x the random comparator",
    "lineage_markers": LINEAGE,
    "random_comparator": {"size": "same as cluster", "seed": 0},
    "tracing_samples": ["EEM-scRNA-249", "EEM-scRNA-250", "EEM-scRNA-251", "EEM-scRNA-288", "EEM-scRNA-289",
                        "EEM-scRNA-290", "EEM-scRNA-291", "EEM-scRNA-292"],
    "unit": "cell-level QC summarised per group; animal composition reported per sample",
}


def main() -> None:
    rec = RunRecord(OUT / "s4_run_record.json", "S4 explain mouse cluster 23", RULES,
                    notes="rules frozen after the cluster-level QC summary was seen, before per-cell data")
    from anndata.io import read_elem  # noqa: E402
    rec.add_input(OBJ)
    with h5py.File(OBJ, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        n_cols = f["X"].attrs["shape"][1]
        cl = obs["leiden_cluster"].astype(str).to_numpy()
        in23 = np.where(cl == CLUSTER)[0]
        rng = np.random.default_rng(RULES["random_comparator"]["seed"])
        others = np.where(cl != CLUSTER)[0]
        rand = np.sort(rng.choice(others, size=len(in23), replace=False))
        gidx = [list(var.index).index(g) for g in LINEAGE]
        counts23 = read_csr_rows(f["layers"]["counts"], in23, n_cols)[:, gidx].toarray()
        countsR = read_csr_rows(f["layers"]["counts"], rand, n_cols)[:, gidx].toarray()

    atlas_med_counts = float(np.median(obs["total_counts"]))
    atlas_med_genes = float(np.median(obs["n_genes_by_counts"]))
    sample = obs["sample_id"].astype(str)
    dominant = sample.iloc[in23].value_counts().index[0]
    groups = {
        "cluster_23": obs.index[in23],
        "all_other_cells": obs.index[others],
        "random_comparator": obs.index[rand],
        f"dominant_sample_{dominant}_outside_23": obs.index[(sample == dominant).to_numpy() & (cl != CLUSTER)],
    }
    qc_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_ribo", "pct_counts_in_top_20_genes", "doublet_score"]
    qc = pd.DataFrame({g: obs.loc[idx, qc_cols].median() for g, idx in groups.items()}).T
    qc["n_cells"] = [len(idx) for idx in groups.values()]
    qc.index.name = "group"
    qc.round(3).to_csv(OUT / "s4_qc_medians_by_group.csv")
    rec.add_output(OUT / "s4_qc_medians_by_group.csv")

    # lineage co-expression
    def coexp(mat):
        det = (mat > 0)
        return {g: float(det[:, i].mean()) for i, g in enumerate(LINEAGE)} | {"ge3_of_5": float((det.sum(axis=1) >= 3).mean()),
                                                                              "ge2_of_5": float((det.sum(axis=1) >= 2).mean())}
    co = pd.DataFrame({"cluster_23": coexp(counts23), "random_comparator": coexp(countsR)}).T
    co.index.name = "group"
    co.round(4).to_csv(OUT / "s4_lineage_marker_detection.csv")
    rec.add_output(OUT / "s4_lineage_marker_detection.csv")

    # cohort split and labelled fraction
    tracing = obs["sample_id"].astype(str).isin(RULES["tracing_samples"]).to_numpy()
    lab = obs["has_author_metadata"].astype(str).eq("True").to_numpy()
    c23 = np.zeros(len(obs), bool); c23[in23] = True
    annotated_cohort = ~tracing
    cohort = {
        "cluster23_cells": int(c23.sum()),
        "cluster23_from_tracing_samples": int((c23 & tracing).sum()),
        "cluster23_from_annotated_samples": int((c23 & annotated_cohort).sum()),
        "cluster23_annotated_cohort_labelled": int((c23 & annotated_cohort & lab).sum()),
        "cluster23_annotated_cohort_labelled_fraction": round(float((c23 & annotated_cohort & lab).sum() / max(1, (c23 & annotated_cohort).sum())), 4),
        "atlas_annotated_cohort_labelled_fraction": round(float((annotated_cohort & lab).sum() / annotated_cohort.sum()), 4),
        "labels_of_the_labelled_cluster23_cells": obs.loc[c23 & lab, "author_celltype"].astype(str).value_counts().to_dict(),
    }
    comp = obs.loc[c23].groupby("sample_id", observed=True).size().sort_values(ascending=False)
    comp_df = pd.DataFrame({"cells_in_cluster23": comp,
                            "fraction_of_cluster23": (comp / comp.sum()).round(4),
                            "fraction_of_that_sample": (comp / sample.value_counts().reindex(comp.index)).round(4),
                            "sample_median_total_counts": obs.groupby("sample_id", observed=True)["total_counts"].median().reindex(comp.index).round(0),
                            "cohort": ["tracing" if s in RULES["tracing_samples"] else "annotated" for s in comp.index]})
    comp_df.index.name = "sample_id"
    comp_df.to_csv(OUT / "s4_sample_composition.csv")
    rec.add_output(OUT / "s4_sample_composition.csv")

    # verdict
    low_count = bool(qc.loc["cluster_23", "total_counts"] < 0.5 * atlas_med_counts and qc.loc["cluster_23", "n_genes_by_counts"] < 0.5 * atlas_med_genes)
    ambient = bool(co.loc["cluster_23", "ge3_of_5"] >= 2 * max(co.loc["random_comparator", "ge3_of_5"], 1e-9))
    verdict = {(True, True): "low-count, ambient-like", (True, False): "low-count only", (False, True): "ambient-like only", (False, False): "unexplained"}[(low_count, ambient)]
    rec.set("atlas_median_total_counts", atlas_med_counts)
    rec.set("atlas_median_n_genes", atlas_med_genes)
    rec.set("dominant_sample", dominant)
    rec.set("qc_medians", qc.round(3).to_dict())
    rec.set("lineage_detection", co.round(4).to_dict())
    rec.set("cohort_split", cohort)
    rec.set("low_count_call", low_count)
    rec.set("ambient_like_call", ambient)
    rec.set("verdict", verdict)

    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    for ax, col, log in zip(axes, ["total_counts", "n_genes_by_counts", "pct_counts_in_top_20_genes"], [True, True, False]):
        data = [obs.loc[idx, col].to_numpy() for idx in groups.values()]
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(groups) + 1)); ax.set_xticklabels([g.replace("_", "\n") for g in groups], fontsize=7)
        ax.set_title(col); ax.set_yscale("log" if log else "linear")
    fig.suptitle(f"Mouse cluster 23 versus comparators; verdict: {verdict}")
    fig.tight_layout(); fig.savefig(OUT / "s4_qc_boxplots.png", dpi=120); plt.close(fig)
    rec.add_output(OUT / "s4_qc_boxplots.png")

    lines = ["# Trial S4 output: mouse cluster 23", "",
             f"Verdict by the frozen rules: **{verdict}**. Low-count call {low_count}; ambient-like call {ambient}.", "",
             "## QC medians", "", df_to_markdown(qc.round(3)), "",
             "## Lineage-marker detection (raw counts > 0)", "", df_to_markdown(co.round(4)), "",
             "## Cohort split", ""] + [f"- {k}: {v}" for k, v in cohort.items()] + \
            ["", "## Sample composition of cluster 23", "", df_to_markdown(comp_df), ""]
    (OUT / "s4_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "s4_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
