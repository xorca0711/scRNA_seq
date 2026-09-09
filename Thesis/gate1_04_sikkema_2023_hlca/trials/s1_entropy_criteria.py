#!/usr/bin/env python
"""Trial S1: apply the HLCA cluster-entropy criteria to this repository's tracked tables.

Sikkema et al. 2023 (Nat Med, doi:10.1038/s41591-023-02327-2) screen every
cluster of the integrated atlas with two Shannon entropies:

* label entropy over cell-type labels (natural log; unlabelled cells excluded;
  not available when fewer than 20% of the cluster is labelled). High above
  0.56, the entropy of a 75/25 two-label cluster. High values flag
  cross-study disagreement or doublet clusters.
* donor entropy over donors. Low below the entropy of a cluster with 95% of
  its cells from one donor and 5% spread evenly over the others. In the HLCA
  (107 donors) that is 0.43; the threshold depends on the donor count and is
  recomputed here for every dataset and stratum.

This script reads only tracked CSV tables (standard library only), freezes
the thresholds from the paper before opening them, and writes one JSON
artefact and one Markdown summary next to itself. Deposited author labels are
used as an answer key only; nothing is fitted.

Pre-registered rules (fixed before any table was read):
  label_entropy_high      > 0.56
  min_labelled_fraction   >= 0.20 for label entropy to be reported
  donor_entropy_low       < H(0.95, 0.05 spread over n-1 donors), per stratum
  min_cells_per_stratum   >= 50 cells of the cluster in the stratum, else NA
  mouse strata            days post infection for the 25 annotated animals;
                          the 8 lineage-tracing animals form one stratum
                          per Cre line (separate experiment, all 19 dpi)
  human strata            none (3 healthy donors, one stratum)
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import date
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MOUSE = REPO / "analysis" / "GSE262927"
HUMAN = REPO / "analysis" / "GSE178360"

RULES = {
    "label_entropy_high": 0.56,
    "min_labelled_fraction": 0.20,
    "donor_private_fraction": 0.95,
    "min_cells_per_stratum": 50,
    "log_base": "natural",
    "source": "Sikkema et al. 2023, Methods: Thresholds for high label/donor entropy and doublet clusters; multidisease analysis (50-cell rule)",
}


def shannon(fractions: list[float]) -> float:
    return -sum(p * math.log(p) for p in fractions if p > 0)


def donor_threshold(n_donors: int, private: float = RULES["donor_private_fraction"]) -> float | None:
    """Entropy of a cluster with `private` of its cells from one donor and the rest spread evenly."""
    if n_donors < 2:
        return None
    rest = (1 - private) / (n_donors - 1)
    return shannon([private] + [rest] * (n_donors - 1))


def read_matrix(path: Path, key: str) -> dict[str, dict[str, float]]:
    """Read a CSV whose first column is `key` and remaining columns are numeric."""
    out: dict[str, dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            k = row.pop(key)
            out[k] = {c: float(v) for c, v in row.items() if v not in ("", None)}
    return out


def transpose(m: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = defaultdict(dict)
    for r, cols in m.items():
        for c, v in cols.items():
            out[c][r] = v
    return dict(out)


def label_entropy_table(counts_by_cluster: dict[str, dict[str, float]],
                        total_by_cluster: dict[str, float]) -> list[dict]:
    rows = []
    # Clusters with no labelled cell at all are absent from the counts table;
    # they must still appear, as "not assessable".
    all_clusters = set(counts_by_cluster) | set(total_by_cluster)
    for cluster in sorted(all_clusters, key=lambda c: int(c)):
        labels = counts_by_cluster.get(cluster, {})
        n_labelled = sum(labels.values())
        total = total_by_cluster.get(cluster, 0.0)
        frac_labelled = n_labelled / total if total else 0.0
        top_label, top_n = max(labels.items(), key=lambda kv: kv[1]) if labels else (None, 0.0)
        if frac_labelled < RULES["min_labelled_fraction"] or n_labelled == 0:
            rows.append({"cluster": cluster, "n_cells": int(total), "labelled_fraction": round(frac_labelled, 4),
                         "label_entropy": None, "flag": "NA (< 20% labelled)",
                         "top_label": top_label, "top_label_fraction": None})
            continue
        h = shannon([v / n_labelled for v in labels.values()])
        rows.append({"cluster": cluster, "n_cells": int(total), "labelled_fraction": round(frac_labelled, 4),
                     "label_entropy": round(h, 4),
                     "flag": "HIGH" if h > RULES["label_entropy_high"] else "low",
                     "top_label": top_label, "top_label_fraction": round(top_n / n_labelled, 4)})
    return rows


def donor_entropy_table(cells_by_sample: dict[str, dict[str, float]],
                        strata: dict[str, str]) -> list[dict]:
    """cells_by_sample[sample][cluster] = n; strata[sample] = stratum name."""
    by_stratum: dict[str, list[str]] = defaultdict(list)
    for s in cells_by_sample:
        by_stratum[strata.get(s, "unassigned")].append(s)
    clusters = sorted({c for cols in cells_by_sample.values() for c in cols}, key=lambda c: int(c))
    rows = []
    for stratum, samples in sorted(by_stratum.items(), key=lambda kv: kv[0]):
        thr = donor_threshold(len(samples))
        for cluster in clusters:
            n_per = {s: cells_by_sample[s].get(cluster, 0.0) for s in samples}
            n = sum(n_per.values())
            row = {"stratum": stratum, "n_donors": len(samples), "threshold": round(thr, 4) if thr else None,
                   "cluster": cluster, "n_cells": int(n)}
            if thr is None:
                row.update({"donor_entropy": None, "flag": "NA (single donor in stratum)"})
            elif n < RULES["min_cells_per_stratum"]:
                row.update({"donor_entropy": None, "flag": "NA (< 50 cells)"})
            else:
                h = shannon([v / n for v in n_per.values()])
                top = max(n_per.values()) / n
                row.update({"donor_entropy": round(h, 4), "max_donor_fraction": round(top, 4),
                            "flag": "LOW (donor-private)" if h < thr else "ok"})
            rows.append(row)
    return rows


def mouse_strata() -> dict[str, str]:
    strata: dict[str, str] = {}
    with (MOUSE / "regeneration_focus" / "tables" / "transitional_abundance_per_sample.csv").open(newline="", encoding="utf-8-sig") as h:
        for row in csv.DictReader(h):
            strata[row["sample"]] = f"day_{int(float(row['day'])):03d}"
    with (MOUSE / "lineage_tracing_cohort" / "tables" / "icap_tracing_by_cre_line.csv").open(newline="", encoding="utf-8-sig") as h:
        for row in csv.DictReader(h):
            strata[row["sample"]] = f"tracing_{row['cre_line']}"
    return strata


def contradicted_mouse_clusters() -> list[str]:
    out = []
    with (MOUSE / "tables" / "cluster_annotation_proposals.csv").open(newline="", encoding="utf-8-sig") as h:
        for row in csv.DictReader(h):
            if "DISAGREES" in row.get("Deposition check", ""):
                key = next((k for k in row if "cluster" in k.lower()), None)
                out.append(str(row[key]) if key else "?")
    return out


def main() -> None:
    result: dict = {
        "trial": "S1 entropy criteria from Sikkema et al. 2023 applied to tracked cluster tables",
        "run_date": date.today().isoformat(),
        "rules_frozen_before_reading_tables": RULES,
        "hlca_reference_thresholds": {"label_entropy_high": 0.56, "donor_entropy_low_107_donors": 0.43,
                                       "donor_threshold_recomputed_here": True},
        "inputs": [],
        "datasets": {},
    }

    # ---------------------------------------------------------------- mouse
    m_counts = read_matrix(MOUSE / "tables" / "cluster_vs_author_celltype_counts.csv", "leiden_cluster")
    m_by_sample = read_matrix(MOUSE / "tables" / "cluster_by_sample_counts.csv", "sample_id")
    m_total = {c: sum(cols.values()) for c, cols in transpose(m_by_sample).items()}
    result["inputs"] += [
        "analysis/GSE262927/tables/cluster_vs_author_celltype_counts.csv",
        "analysis/GSE262927/tables/cluster_by_sample_counts.csv",
        "analysis/GSE262927/tables/cluster_annotation_proposals.csv",
        "analysis/GSE262927/regeneration_focus/tables/transitional_abundance_per_sample.csv",
        "analysis/GSE262927/lineage_tracing_cohort/tables/icap_tracing_by_cre_line.csv",
    ]
    m_label = label_entropy_table(m_counts, m_total)
    strata = mouse_strata()
    m_donor_strat = donor_entropy_table(m_by_sample, strata)
    m_donor_global = donor_entropy_table(m_by_sample, {s: "all_33_samples_unstratified" for s in m_by_sample})
    contradicted = contradicted_mouse_clusters()
    high = [r["cluster"] for r in m_label if r["flag"] == "HIGH"]
    low_strat = sorted({r["cluster"] for r in m_donor_strat if r["flag"].startswith("LOW")}, key=int)
    low_global = sorted({r["cluster"] for r in m_donor_global if r["flag"].startswith("LOW")}, key=int)
    result["datasets"]["GSE262927"] = {
        "n_clusters": len(m_label),
        "n_samples": len(m_by_sample),
        "strata": {k: sorted(v for v, s in strata.items() if s == k) for k in sorted(set(strata.values()))},
        "label_entropy": m_label,
        "label_entropy_summary": {
            "clusters_high": high,
            "clusters_na": [r["cluster"] for r in m_label if r["flag"].startswith("NA")],
            "contradicted_by_deposited_labels_in_annotation_table": contradicted,
            "overlap_high_and_contradicted": sorted(set(high) & set(contradicted), key=int),
        },
        "donor_entropy_stratified": m_donor_strat,
        "donor_entropy_unstratified_negative_control": m_donor_global,
        "donor_entropy_summary": {
            "clusters_low_within_stratum": low_strat,
            "clusters_low_unstratified": low_global,
            "note": "unstratified flags are expected to include time-point-specific states because animal identity and time point coincide by design",
        },
    }

    # ---------------------------------------------------------------- human
    h_by_sample = read_matrix(HUMAN / "tables" / "cluster_by_sample_counts.csv", "sample_id")
    result["inputs"].append("analysis/GSE178360/tables/cluster_by_sample_counts.csv")
    h_donor = donor_entropy_table(h_by_sample, {s: "three_healthy_donors" for s in h_by_sample})
    result["datasets"]["GSE178360"] = {
        "n_clusters": len({c for cols in h_by_sample.values() for c in cols}),
        "n_samples": len(h_by_sample),
        "label_entropy": None,
        "label_entropy_status": "Not established: deposited author labels for GSE178360 are inside .RDS objects that cannot be read here (no R)",
        "donor_entropy": h_donor,
        "donor_entropy_summary": {
            "threshold_for_3_donors": round(donor_threshold(3), 4),
            "clusters_low": sorted({r["cluster"] for r in h_donor if r["flag"].startswith("LOW")}, key=int),
            "clusters_na": sorted({r["cluster"] for r in h_donor if r["flag"].startswith("NA")}, key=int),
            "comparison": "the existing check counts clusters with > 75% of cells from one donor (4 of 31 after Harmony)",
        },
    }

    (HERE / "s1_entropy_criteria.json").write_text(json.dumps(result, indent=2), encoding="utf-8")

    # ---------------------------------------------------------- markdown
    lines = ["# Trial S1 output: HLCA entropy criteria on tracked cluster tables", "",
             f"Generated by `s1_entropy_criteria.py` on {result['run_date']} from tracked CSV tables only. "
             "Thresholds were frozen from the paper before the tables were read (see the script docstring). "
             "Full detail in `s1_entropy_criteria.json`.", ""]
    g = result["datasets"]["GSE262927"]
    lines += ["## GSE262927 (mouse), label entropy against deposited labels", "",
              "| Cluster | Cells | Labelled fraction | Label entropy | Flag | Top deposited label | Top fraction |",
              "|--:|--:|--:|--:|---|---|--:|"]
    for r in g["label_entropy"]:
        lines.append(f"| {r['cluster']} | {r['n_cells']} | {r['labelled_fraction']} | "
                     f"{'' if r['label_entropy'] is None else r['label_entropy']} | {r['flag']} | "
                     f"{r['top_label'] or ''} | {'' if r['top_label_fraction'] is None else r['top_label_fraction']} |")
    s = g["label_entropy_summary"]
    lines += ["", f"High label entropy (> 0.56): clusters {', '.join(s['clusters_high']) or 'none'}. "
              f"Not assessable (< 20% labelled): {', '.join(s['clusters_na']) or 'none'}. "
              f"Clusters already marked CONTRADICTED in the annotation table: {', '.join(s['contradicted_by_deposited_labels_in_annotation_table']) or 'none'}; "
              f"overlap with the high-entropy set: {', '.join(s['overlap_high_and_contradicted']) or 'none'}.", ""]
    lines += ["## GSE262927 (mouse), donor entropy within stratum", "",
              "| Stratum | Animals | Threshold | Clusters flagged LOW (donor-private) | Clusters NA (< 50 cells) |",
              "|---|--:|--:|---|---|"]
    by_stratum: dict[str, list[dict]] = defaultdict(list)
    for r in g["donor_entropy_stratified"]:
        by_stratum[r["stratum"]].append(r)
    for stratum, rows in by_stratum.items():
        low = [r["cluster"] for r in rows if r["flag"].startswith("LOW")]
        na = [r["cluster"] for r in rows if r["flag"].startswith("NA")]
        lines.append(f"| {stratum} | {rows[0]['n_donors']} | {rows[0]['threshold']} | {', '.join(low) or 'none'} | {', '.join(na) or 'none'} |")
    d = g["donor_entropy_summary"]
    lines += ["", f"Unstratified negative control (33 animals as if one stratum, threshold "
              f"{g['donor_entropy_unstratified_negative_control'][0]['threshold']}): clusters flagged "
              f"{', '.join(d['clusters_low_unstratified']) or 'none'}. Within-stratum union: "
              f"{', '.join(d['clusters_low_within_stratum']) or 'none'}. {d['note']}.", ""]
    hh = result["datasets"]["GSE178360"]
    lines += ["## GSE178360 (human), donor entropy over three donors", "",
              f"Label entropy: {hh['label_entropy_status']}.", "",
              f"Donor-entropy threshold for 3 donors: {hh['donor_entropy_summary']['threshold_for_3_donors']}. "
              f"Clusters flagged LOW: {', '.join(hh['donor_entropy_summary']['clusters_low']) or 'none'}. "
              f"NA (< 50 cells): {', '.join(hh['donor_entropy_summary']['clusters_na']) or 'none'}. "
              f"{hh['donor_entropy_summary']['comparison']}.", "",
              "| Cluster | Cells | Donor entropy | Max donor fraction | Flag |", "|--:|--:|--:|--:|---|"]
    for r in hh["donor_entropy"]:
        lines.append(f"| {r['cluster']} | {r['n_cells']} | {'' if r.get('donor_entropy') is None else r['donor_entropy']} | "
                     f"{'' if r.get('max_donor_fraction') is None else r['max_donor_fraction']} | {r['flag']} |")
    (HERE / "s1_entropy_criteria.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({"mouse_high_label_entropy": s["clusters_high"], "mouse_contradicted": s["contradicted_by_deposited_labels_in_annotation_table"],
                      "mouse_low_donor_within_stratum": d["clusters_low_within_stratum"], "mouse_low_donor_unstratified": d["clusters_low_unstratified"],
                      "human_low_donor": hh["donor_entropy_summary"]["clusters_low"], "human_threshold": hh["donor_entropy_summary"]["threshold_for_3_donors"]}, indent=1))


if __name__ == "__main__":
    main()
