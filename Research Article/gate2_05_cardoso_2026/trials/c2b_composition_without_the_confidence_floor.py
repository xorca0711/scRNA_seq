#!/usr/bin/env python
"""Trial C2b: the same composition, without the confidence floor that hid two directions.

Trial C2 reported three of its five pre-registered directions met. Two of the
two that failed did so in the same way, and for the same reason trial C1's
rule defeated itself: the composition metric counts only clusters whose modal
population call holds at least 50% of their cells, and the reprogrammed
fibroblast and Cd177-positive calls never clear that floor. Both directions
were therefore scored as 0.0% against 0.0%, which is not a measurement of
absence, it is a measurement of nothing.

This trial recomputes the same shares with the floor removed. It reads only
the cluster tables C2 already wrote; it loads no matrix, fits nothing, and
changes no label. C2's outcome stands as recorded and this sits beside it.

Provenance: fixed after C2's cluster tables were seen. Post hoc and disclosed.

Frozen rules:

* Inputs: `c2_niche_clusters.csv` and `c2_epithelial_clusters.csv`, the
  tracked outputs of trial C2.
* Every cluster is assigned its modal call, whatever the mode's size. The
  mode fraction is carried alongside every number so a reader can see how
  weak each assignment is.
* Shares are within a library and within a compartment, exactly as C2
  computed them; only the floor differs.
* The two directions C2 could not score are re-scored: the reprogrammed
  fibroblast share of fibroblasts, and the Cd177-positive share of RFP+
  cells. A direction is "met" if the flox/flox share is lower than the
  flox/+ share.
* Reported beside each: the single most depleted cluster of each call, since
  a share computed over weakly-called clusters is a weaker statement than one
  cluster's own count.
* Unit: the library. One library per genotype per sort. No P value.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "c2b_composition_without_the_confidence_floor"
OUT.mkdir(exist_ok=True)
SOURCE = HERE / "c2_areg_deletion_arm"
GENOTYPES = ["Areg-flox/+", "Areg-flox/flox"]
FIBROBLAST_CALLS = {"alveolar fibroblast", "adventitial fibroblast", "reprogrammed fibroblast"}

RULES = {
    "inputs": ["c2_areg_deletion_arm/c2_niche_clusters.csv",
               "c2_areg_deletion_arm/c2_epithelial_clusters.csv"],
    "provenance_of_these_rules": "fixed after C2's cluster tables were seen; post hoc and disclosed",
    "change_from_c2": "every cluster is assigned its modal call whatever the mode's size; C2 required the mode to hold at least 50%",
    "what_is_unchanged": "shares are within a library and within a compartment, as C2 computed them",
    "directions_rescored": ["reprogrammed fibroblast share of fibroblasts falls in flox/flox",
                            "Cd177-positive share of RFP+ cells falls in flox/flox"],
    "also_reported": "the single most depleted cluster of each call, because a share over weakly called clusters is a weaker statement than one cluster's count",
    "unit": "the library; one per genotype per sort; no P value",
}


def shares(table: pd.DataFrame, keep: set[str] | None, label: str) -> pd.DataFrame:
    rows = []
    subset = table if keep is None else table[table["call"].isin(keep)]
    for genotype in GENOTYPES:
        column = f"n_{genotype}"
        total = int(subset[column].sum())
        for call in sorted(subset["call"].unique()):
            hit = subset[subset["call"] == call]
            n = int(hit[column].sum())
            rows.append({
                "compartment": label, "genotype": genotype, "call": call,
                "n_cells": n,
                "pct_of_compartment": round(100 * n / total, 2) if total else None,
                "n_clusters": int(len(hit)),
                "weakest_mode_fraction": round(float(hit["mode_fraction"].min()), 3),
            })
    return pd.DataFrame(rows)


def main() -> None:
    rec = RunRecord(OUT / "c2b_run_record.json",
                    "C2b composition without the confidence floor", RULES,
                    notes="reads only C2's tracked cluster tables; loads no matrix")
    niche = pd.read_csv(SOURCE / "c2_niche_clusters.csv")
    epi = pd.read_csv(SOURCE / "c2_epithelial_clusters.csv")
    for path in (SOURCE / "c2_niche_clusters.csv", SOURCE / "c2_epithelial_clusters.csv"):
        rec.add_input(path)

    fib = shares(niche, FIBROBLAST_CALLS, "fibroblasts")
    epi_shares = shares(epi, None, "RFP+ epithelium")
    fib.to_csv(OUT / "c2b_fibroblast_composition.csv", index=False)
    epi_shares.to_csv(OUT / "c2b_epithelial_composition.csv", index=False)
    rec.add_output(OUT / "c2b_fibroblast_composition.csv")
    rec.add_output(OUT / "c2b_epithelial_composition.csv")

    def share_of(frame, genotype, call):
        hit = frame[(frame["genotype"] == genotype) & (frame["call"] == call)]
        return float(hit["pct_of_compartment"].iloc[0]) if len(hit) else 0.0

    outcome = {}
    for frame, call, key in ((fib, "reprogrammed fibroblast", "reprogrammed_fibroblast_share_falls_in_hom"),
                             (epi_shares, "Cd177_positive", "Cd177_positive_share_falls_in_hom")):
        het, hom = share_of(frame, GENOTYPES[0], call), share_of(frame, GENOTYPES[1], call)
        outcome[key] = {"flox_plus_pct": het, "flox_flox_pct": hom,
                        "direction_met": bool(hom < het) if (het or hom) else None,
                        "scorable": bool(het or hom)}
    rec.set("directions_rescored", outcome)

    # the single most depleted cluster of each call, as the stronger statement
    per_cluster = []
    for table, label in ((niche, "niche"), (epi, "RFP+ epithelium")):
        for _, row in table.iterrows():
            het, hom = int(row[f"n_{GENOTYPES[0]}"]), int(row[f"n_{GENOTYPES[1]}"])
            per_cluster.append({
                "compartment": label, "cluster": row["cluster"], "call": row["call"],
                "mode_fraction": row["mode_fraction"], "n_cells": int(row["n_cells"]),
                "n_flox_plus": het, "n_flox_flox": hom,
                "ratio_flox_plus_over_flox_flox": round(het / hom, 2) if hom else None,
            })
    depletion = pd.DataFrame(per_cluster).sort_values(
        "ratio_flox_plus_over_flox_flox", ascending=False, na_position="last")
    depletion.to_csv(OUT / "c2b_cluster_depletion.csv", index=False)
    rec.add_output(OUT / "c2b_cluster_depletion.csv")
    rec.set("most_depleted_clusters",
            depletion.head(6)[["compartment", "cluster", "call", "n_flox_plus", "n_flox_flox",
                               "ratio_flox_plus_over_flox_flox"]].to_dict("records"))

    lines = [
        "# Trial C2b output: the same composition, without the confidence floor", "",
        "Trial C2 scored two of its five pre-registered directions as 0.0% against 0.0%, because the",
        "reprogrammed fibroblast and Cd177-positive calls never reached the 50% mode floor. That is",
        "not a measurement of absence. Removing the floor, and carrying the weakest mode fraction",
        "beside every number so the reader can see how weak each call is:", "",
        "## Fibroblasts", "", df_to_markdown(fib, index=False), "",
        "## RFP+ epithelium", "", df_to_markdown(epi_shares, index=False), "",
        "## The directions C2 could not score", "",
        df_to_markdown(pd.DataFrame([{"direction": k, **v} for k, v in outcome.items()]), index=False), "",
        "## Cluster-level depletion, the stronger statement", "",
        "A share computed over weakly called clusters is weaker than one cluster's own count.",
        "Clusters most depleted in the Areg-flox/flox library:", "",
        df_to_markdown(depletion.head(8)[["compartment", "cluster", "call", "mode_fraction",
                                          "n_flox_plus", "n_flox_flox",
                                          "ratio_flox_plus_over_flox_flox"]], index=False), "",
        "One library per genotype per sort: these are two libraries, not two groups of animals.", "",
    ]
    (OUT / "c2b_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c2b_summary.md")
    rec.finish()
    print(df_to_markdown(pd.DataFrame([{"direction": k, **v} for k, v in outcome.items()]), index=False))
    print()
    print(df_to_markdown(depletion.head(8)[["compartment", "cluster", "call", "n_flox_plus",
                                            "n_flox_flox", "ratio_flox_plus_over_flox_flox"]], index=False))


if __name__ == "__main__":
    main()
