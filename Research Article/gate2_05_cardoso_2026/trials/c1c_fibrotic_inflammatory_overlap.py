#!/usr/bin/env python
"""Trial C1c: are fibrotic and inflammatory fibroblasts separate cells here (Gate 1 addendum).

Trial C1b found that the Red2Kras-private fibroblast structure carries the
fibrotic markers (Tnc, Runx1, Acta2, Pdgfrb) and the inflammatory markers
(Lcn2, Saa3) in the same clusters. A cluster-level detection fraction cannot
tell two interleaved subsets from genuine co-expression, and the difference
matters: the paper's Fig. 2 and its discussion make fibrotic and inflammatory
fibroblasts **distinct populations**, with the inflammatory cells appearing
only from 4 weeks, sitting at the tumour periphery, and specifically *lacking*
Tnc. These libraries are 2 weeks old, where the paper expects inflammatory
fibroblasts to be rare.

This trial asks the per-cell question directly. It is cheap: it reads the
object C1b already wrote and computes co-detection; it fits no model and
changes no label.

Provenance: these rules were fixed after C1b's cluster tables were seen and
before any per-cell co-detection was computed. Stated here because they are
not blind to C1b.

Frozen rules:

* Object: the one C1b wrote (same QC, doublet removal, embedding and Leiden
  labels); read in backed mode, three marker columns only.
* Cells considered: the Red2Kras cells of the clusters C1b called
  Red2Kras-private (at least 90% Red2Kras), and, as a comparator, the
  Red2Kras cells of every other cluster.
* Fibrotic-positive: Tnc detected. Inflammatory-positive: Lcn2 or Saa3
  detected. Both on raw counts, detection meaning a non-zero count.
* THE TEST. If the two programmes mark separate cells, the observed
  co-detection rate should not exceed what independence predicts. Report the
  observed fraction of double-positive cells against the product of the two
  marginal fractions, as a ratio. Frozen reading: a ratio at or below 1.25 is
  "consistent with separate cells"; above 1.25 is "co-expressed more often
  than independence predicts". No P value: one library, three pooled mice.
* The same quantities are reported for the 12-week and 4-day comparison
  libraries only if they exist in this object; they do not, and that absence
  is itself reported, because the paper's timing claim cannot be checked here.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import REPO, RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "c1c_fibrotic_inflammatory_overlap"
OUT.mkdir(exist_ok=True)
OBJECT = HERE / "c1b_characterise_red2kras_private" / "c1b_mesenchyme.h5ad"
PRIVATE_TABLE = HERE / "c1b_characterise_red2kras_private" / "c1b_question_b_red2kras_private.csv"

RULES = {
    "object": "the one trial C1b wrote; read in backed mode, marker columns only",
    "provenance_of_these_rules": "fixed after C1b's cluster tables were seen, before any per-cell co-detection was computed",
    "fibrotic_positive": "Tnc detected (non-zero raw count)",
    "inflammatory_positive": "Lcn2 or Saa3 detected (non-zero raw count)",
    "groups": "Red2Kras cells of the C1b Red2Kras-private clusters; Red2Kras cells of all other clusters as comparator",
    "test": {"statistic": "observed double-positive fraction divided by the product of the marginals",
             "consistent_with_separate_cells_if": "<= 1.25", "co_expressed_if": "> 1.25"},
    "no_p_value": "one library, three pooled mice",
    "timing_claim": "not checkable here: this object is a single 2-week time point",
}


def main() -> None:
    import anndata as ad

    rec = RunRecord(OUT / "c1c_run_record.json",
                    "C1c fibrotic and inflammatory marker co-detection (Gate 1 addendum)",
                    RULES, notes="reads the object C1b wrote; fits nothing")
    rec.add_input(OBJECT)
    rec.add_input(PRIVATE_TABLE)
    private = pd.read_csv(PRIVATE_TABLE)["cluster"].astype(str).tolist()
    rec.set("private_clusters_from_c1b", private)

    adata = ad.read_h5ad(OBJECT, backed="r")
    obs = adata.obs
    genes = ["Tnc", "Lcn2", "Saa3"]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        raise SystemExit(f"missing from the object: {missing}")
    cols = [adata.var_names.get_loc(g) for g in genes]
    counts = adata.layers["counts"]
    detected = {}
    for gene, col in zip(genes, cols):
        detected[gene] = np.asarray(counts[:, col].todense()).ravel() > 0
    adata.file.close()

    fibrotic = detected["Tnc"]
    inflammatory = detected["Lcn2"] | detected["Saa3"]
    red2kras = (obs["genotype"] == "Red2Kras").to_numpy()
    in_private = obs["leiden_0.5"].astype(str).isin(private).to_numpy()

    rows = []
    for label, mask in (("Red2Kras, C1b private clusters", red2kras & in_private),
                        ("Red2Kras, all other clusters", red2kras & ~in_private),
                        ("Confetti, all clusters", ~red2kras)):
        n = int(mask.sum())
        if n == 0:
            continue
        p_fib = float(fibrotic[mask].mean())
        p_inf = float(inflammatory[mask].mean())
        p_both = float((fibrotic & inflammatory)[mask].mean())
        expected = p_fib * p_inf
        ratio = float(p_both / expected) if expected > 0 else np.nan
        rows.append({
            "group": label, "n_cells": n,
            "pct_Tnc": round(100 * p_fib, 2),
            "pct_Lcn2_or_Saa3": round(100 * p_inf, 2),
            "pct_double_positive": round(100 * p_both, 2),
            "pct_expected_if_independent": round(100 * expected, 2),
            "observed_over_expected": round(ratio, 3) if ratio == ratio else None,
            "reading": ("consistent with separate cells" if ratio == ratio and ratio <= 1.25
                        else "co-expressed more often than independence predicts"),
        })
    overall = pd.DataFrame(rows)
    overall.to_csv(OUT / "c1c_overlap_by_group.csv", index=False)
    rec.add_output(OUT / "c1c_overlap_by_group.csv")

    per_cluster = []
    clusters = obs["leiden_0.5"].astype(str).to_numpy()
    for cluster in sorted(set(clusters), key=int):
        mask = (clusters == cluster) & red2kras
        n = int(mask.sum())
        if n < 50:
            continue
        p_fib = float(fibrotic[mask].mean())
        p_inf = float(inflammatory[mask].mean())
        p_both = float((fibrotic & inflammatory)[mask].mean())
        expected = p_fib * p_inf
        per_cluster.append({
            "cluster": cluster, "n_Red2Kras_cells": n,
            "is_C1b_private": cluster in private,
            "pct_Tnc": round(100 * p_fib, 2),
            "pct_Lcn2_or_Saa3": round(100 * p_inf, 2),
            "pct_double_positive": round(100 * p_both, 2),
            "observed_over_expected": round(float(p_both / expected), 3) if expected > 0 else None,
        })
    per_cluster_frame = pd.DataFrame(per_cluster)
    per_cluster_frame.to_csv(OUT / "c1c_overlap_by_cluster.csv", index=False)
    rec.add_output(OUT / "c1c_overlap_by_cluster.csv")
    rec.set("overlap_by_group", rows)
    rec.set("timing_claim_checkable", False)

    lines = [
        "# Trial C1c output: are fibrotic and inflammatory fibroblasts separate cells here", "",
        "The paper makes fibrotic and inflammatory fibroblasts distinct populations, with the",
        "inflammatory cells appearing from 4 weeks, sitting at the tumour periphery and lacking Tnc.",
        "These libraries are a single 2-week time point, so the timing half of that claim cannot be",
        "checked here at all. What can be checked is whether the two marker programmes mark the same",
        "cells at 2 weeks.", "",
        df_to_markdown(overall, index=False), "",
        "## By cluster (Red2Kras cells only, clusters with at least 50)", "",
        df_to_markdown(per_cluster_frame, index=False), "",
        "One library, three pooled mice: no P value, and the ratio is a description of one library.", "",
    ]
    (OUT / "c1c_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c1c_summary.md")
    rec.finish()
    print(df_to_markdown(overall, index=False))


if __name__ == "__main__":
    main()
