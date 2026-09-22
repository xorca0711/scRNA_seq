"""Corrected epithelial-source resource scan, retaining historical analyses intact."""
from __future__ import annotations

import argparse
import gc
import importlib.metadata
import itertools
import json
import math
import time
import warnings

import numpy as np
import pandas as pd

from common import (HERE, ROOT, LIGANDS, RESOURCES, TARGETS, abundance_guard,
                    allowed_pairs, eligible_donors, sha256, write_json)

OUT = HERE / "results" / "lr"
CACHE = ROOT / "raw_data" / "GSE136831" / "c12_epithelial_stromal.h5ad"
HISTORICAL = ROOT / "Thesis" / "gate2_05_cardoso_2026" / "trials"


def aggregate(winners, donor_count):
    result = (winners.groupby("pair").agg(n_donors=("donor", "nunique"),
        median_rank=("rank", "median"), ligand=("ligand_complex", "first"),
        receptor=("receptor_complex", "first"),
        winning_sources=("source", lambda x: ";".join(sorted(set(x)))),
        winning_targets=("target", lambda x: ";".join(sorted(set(x))))).reset_index())
    result = result.loc[result.n_donors >= max(3, math.ceil(donor_count / 2))]
    result = result.sort_values(["median_rank", "pair"]).reset_index(drop=True)
    result["overall_rank"] = np.arange(1, len(result) + 1)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resources", nargs="+", choices=RESOURCES, default=RESOURCES)
    args = parser.parse_args()
    warnings.filterwarnings("ignore", category=FutureWarning)
    import anndata as ad
    import scanpy as sc
    from liana.method import cellchat
    from liana.resource import select_resource

    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    checkpoints = OUT / "checkpoints"
    checkpoints.mkdir(exist_ok=True)
    input_paths = [CACHE, HERE / "specification.json", HERE / "common.py", HERE / "c12_c14_corrected.py"]
    hashes = {str(p.relative_to(ROOT)): sha256(p) for p in input_paths}
    fingerprint = json.dumps(hashes, sort_keys=True)
    versions = {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy", "anndata", "scanpy", "liana"]}
    fingerprint += json.dumps(versions, sort_keys=True)
    data = ad.read_h5ad(CACHE, backed="r")
    eligibility = eligible_donors(data.obs)
    eligibility.to_csv(OUT / "donor_eligibility.csv")
    donors = eligibility.index[eligibility.eligible].tolist()
    epithelial = set(data.obs.loc[data.obs.category.eq("Epithelial"), "celltype"])
    if epithelial & TARGETS:
        raise ValueError("deposited source and target cell types overlap")
    write_json(OUT / "scope.json", {"senders": sorted(epithelial), "targets": sorted(TARGETS),
        "eligible_donors": donors, "n_donors": len(donors), "fingerprint": fingerprint})
    print(f"LR: {len(donors)} donors with >=50 actual source and target cells", flush=True)
    for resource in args.resources:
        table = select_resource(resource)
        table.to_csv(OUT / f"resource_{resource}.csv", index=False)
    # Read only one donor's full gene matrix at a time. Molecule normalization is
    # identical to C12; the source restriction is applied before LIANA scoring.
    for donor in donors:
        pending = []
        for resource in args.resources:
            stem = checkpoints / f"{resource}_{donor}"
            sidecar = stem.with_suffix(".json")
            scored_path = stem.with_suffix(".csv.gz")
            if (sidecar.exists() and scored_path.exists()
                and json.loads(sidecar.read_text())["fingerprint"] == fingerprint
                and json.loads(sidecar.read_text())["sha256"] == sha256(scored_path)):
                continue
            pending.append(resource)
        if not pending:
            continue
        keep = data.obs.donor.eq(donor) & (data.obs.category.eq("Epithelial") | data.obs.celltype.isin(TARGETS))
        part = data[keep].to_memory()
        present = part.obs.celltype.value_counts()
        part = part[part.obs.celltype.isin(present[present >= 10].index)].copy()
        groups = allowed_pairs(part.obs)
        if groups.empty:
            raise ValueError(f"{donor}: no subtype pairs above the ten-cell floor")
        sc.pp.normalize_total(part, target_sum=1e4)
        sc.pp.log1p(part)
        for resource in pending:
            tick = time.time()
            # groupby_pairs prevents scoring disallowed origins in the first place.
            result = cellchat(part, groupby="celltype", resource_name=resource,
                groupby_pairs=groups, expr_prop=0.1, min_cells=10, use_raw=False,
                n_perms=None, verbose=False, inplace=False)
            if result is None or result.empty:
                raise ValueError(f"{resource}/{donor}: no result; refusing silent exclusion")
            if not result.source.isin(epithelial).all() or not result.target.isin(TARGETS).all():
                raise AssertionError("scorer emitted a forbidden source or target")
            result["donor"] = donor
            result["resource"] = resource
            result["pair"] = result.ligand_complex + " to " + result.receptor_complex
            columns = ["donor", "resource", "pair", "ligand_complex", "receptor_complex",
                       "source", "target", "lr_probs"]
            scored = result[columns].sort_values(["pair", "source", "target"])
            path = checkpoints / f"{resource}_{donor}.csv.gz"
            scored.to_csv(path, index=False)
            write_json(path.with_suffix(".json"), {"fingerprint": fingerprint, "sha256": sha256(path),
                "rows": len(scored), "elapsed_seconds": round(time.time() - tick, 1)})
            # Use the same sidecar name for cache inspection above.
            write_json(checkpoints / f"{resource}_{donor}.json", {"fingerprint": fingerprint,
                "sha256": sha256(path), "rows": len(scored)})
            print(f"LR {resource}/{donor}: {len(scored)} allowed subtype pairs, {time.time()-tick:.1f}s", flush=True)
            del result, scored
        del part
        gc.collect()
    data.file.close()

    summaries, rankings, comparisons, egfr_components = [], [], [], []
    tops = {}
    for resource in args.resources:
        frames = [pd.read_csv(checkpoints / f"{resource}_{donor}.csv.gz") for donor in donors]
        scored = pd.concat(frames, ignore_index=True)
        scored.to_csv(OUT / f"{resource}_all_source_target_scores.csv.gz", index=False)
        # Frozen deterministic tie order makes repeat runs stable. This is the
        # original maximum-then-ordinal-rank procedure with explicit tie breaks.
        winners = scored.sort_values(["donor", "lr_probs", "pair", "source", "target"],
            ascending=[True, False, True, True, True]).drop_duplicates(["donor", "pair"]).copy()
        winners["rank"] = winners.groupby("donor").cumcount() + 1
        winners.to_csv(OUT / f"{resource}_donor_pair_winners.csv", index=False)
        ranked = aggregate(winners, len(donors))
        ranked["resource"] = resource
        rankings.append(ranked)
        top = ranked.head(15)
        tops[resource] = set(top.pair)
        exact = ranked.loc[ranked.receptor.eq("EGFR") & ranked.ligand.isin(LIGANDS)]
        component = ranked.loc[ranked.receptor.str.split("_").apply(lambda s: "EGFR" in s)]
        egfr_components.append(component)
        order = exact.ligand.tolist()
        old_path = HISTORICAL / "c14_does_the_ranking_depend_on_the_database" / f"c14_partial_{resource}.csv"
        old = pd.read_csv(old_path)
        merged = old[["pair", "median_rank", "overall_rank", "n_donors"]].merge(
            ranked[["pair", "median_rank", "overall_rank", "n_donors"]], on="pair", how="outer",
            suffixes=("_historical", "_corrected"))
        merged["resource"] = resource
        comparisons.append(merged)
        share = float(np.mean([abundance_guard(l, r) for l, r in zip(top.ligand, top.receptor)]))
        areg = exact.loc[exact.ligand.eq("AREG")]
        summaries.append({"resource": resource, "n_donors": len(donors), "pairs": len(ranked),
            "top_pair": top.iloc[0].pair, "abundance_share_top15": share,
            "historical_guard_fires": share > .5,
            "exact_egfr_ligand_order": ", ".join(order),
            "areg_rank_among_exact_egfr": order.index("AREG") + 1 if "AREG" in order else None,
            "areg_median_rank": float(areg.iloc[0].median_rank) if len(areg) else None,
            "egfr_component_pairs": len(component),
            "shared_historical_top15": len(set(old.head(15).pair) & tops[resource]),
            "note": "Guard is a composition flag, not a validity test for signalling"})
    pd.concat(rankings, ignore_index=True).to_csv(OUT / "rankings_by_resource.csv", index=False)
    pd.concat(comparisons, ignore_index=True).to_csv(OUT / "historical_rank_comparison.csv", index=False)
    pd.concat(egfr_components, ignore_index=True).to_csv(OUT / "egfr_component_pairs.csv", index=False)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "resource_summary.csv", index=False)
    overlaps = [{"resource_a": a, "resource_b": b, "shared_top15": len(tops[a] & tops[b]),
        "jaccard": len(tops[a] & tops[b]) / len(tops[a] | tops[b])}
        for a, b in itertools.combinations(args.resources, 2)]
    pd.DataFrame(overlaps).to_csv(OUT / "top15_overlap.csv", index=False)
    record = {"inputs_sha256": hashes, "versions": versions, "resources_completed": args.resources,
        "resource_sha256": {r: sha256(OUT / f"resource_{r}.csv") for r in args.resources},
        "elapsed_seconds": round(time.time() - started, 1), "donors": donors,
        "design": "Corrected post-audit sensitivity; same cohort, no independent validation",
        "provenance_checks": {"all_sources_epithelial": True, "all_targets_fibroblast": True,
                              "actual_compartment_floor": 50}}
    write_json(OUT / "run_record.json", record)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
