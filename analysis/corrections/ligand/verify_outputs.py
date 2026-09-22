"""Check completed corrected evidence and record exact artifact provenance."""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from common import HERE, ROOT, LIGANDS, RESOURCES, TARGETS, paired_summary, sha256, write_json


def main():
    base = HERE / "results"
    lr = base / "lr"
    c37 = base / "c37"
    scope = json.loads((lr / "scope.json").read_text())
    allowed = set(scope["senders"])
    donors = set(scope["eligible_donors"])
    eligibility = pd.read_csv(lr / "donor_eligibility.csv", index_col="donor")
    assert (eligibility.loc[list(donors), ["n_epithelial", "n_fibroblast"]] >= 50).all().all()
    checks = {"eligible_donors": len(donors), "resources": {}, "warnings": []}
    # Verify that the subtype floor did not silently reduce either scored
    # compartment below the donor floor after eligibility was assessed.
    import anndata as ad
    cached = ad.read_h5ad(ROOT / "raw_data/GSE136831/c12_epithelial_stromal.h5ad", backed="r")
    post_floor = []
    for donor in sorted(donors):
        obs = cached.obs.loc[cached.obs.donor.eq(donor)
            & (cached.obs.category.eq("Epithelial") | cached.obs.celltype.isin(TARGETS))]
        counts = obs.celltype.value_counts()
        obs = obs.loc[obs.celltype.isin(counts[counts >= 10].index)]
        post_floor.append({"donor": donor, "scored_epithelial": int(obs.category.eq("Epithelial").sum()),
                           "scored_fibroblast": int(obs.celltype.isin(TARGETS).sum())})
    cached.file.close()
    post_floor = pd.DataFrame(post_floor)
    assert (post_floor[["scored_epithelial", "scored_fibroblast"]] >= 50).all().all()
    post_floor.to_csv(lr / "post_subtype_floor_counts.csv", index=False)
    checks["minimum_scored_compartment_cells"] = {k: int(v) for k, v in
        post_floor[["scored_epithelial", "scored_fibroblast"]].min().to_dict().items()}
    historical = ROOT / "Thesis" / "gate2_05_cardoso_2026" / "trials"
    old = pd.read_csv(historical / "c12_cellchatdb_full_resource_scan" / "c12_per_donor_pairs.csv")
    hashes = {}
    receptor_coverage = []
    rankings = pd.read_csv(lr / "rankings_by_resource.csv")
    for resource in RESOURCES:
        winner_path = lr / f"{resource}_donor_pair_winners.csv"
        all_path = lr / f"{resource}_all_source_target_scores.csv.gz"
        winners, scored = pd.read_csv(winner_path), pd.read_csv(all_path)
        assert set(winners.donor) == donors == set(scored.donor)
        assert winners.source.isin(allowed).all() and scored.source.isin(allowed).all()
        assert winners.target.isin(TARGETS).all() and scored.target.isin(TARGETS).all()
        assert not winners.duplicated(["donor", "pair"]).any()
        assert scored.lr_probs.between(0, 1).all()
        maxima = scored.groupby(["donor", "pair"]).lr_probs.max()
        joined = winners.set_index(["donor", "pair"]).lr_probs.reindex(maxima.index)
        np.testing.assert_allclose(joined, maxima, rtol=1e-12, atol=1e-12)
        for _, frame in winners.groupby("donor"):
            assert sorted(frame["rank"]) == list(range(1, len(frame) + 1))
        checks["resources"][resource] = {"scored_subtype_pairs": len(scored),
            "donor_pair_winners": len(winners), "forbidden_senders": 0, "forbidden_targets": 0,
            "donor_count": len(donors)}
        curated = pd.read_csv(lr / f"resource_{resource}.csv")
        curated = curated.loc[curated.ligand.isin(LIGANDS)
            & curated.receptor.str.split("_").apply(lambda parts: "EGFR" in parts)]
        retained = set(rankings.loc[rankings.resource.eq(resource), "pair"])
        for row in curated.itertuples(index=False):
            pair = row.ligand + " to " + row.receptor
            present = winners.loc[winners.pair.eq(pair)]
            receptor_coverage.append({"resource": resource, "ligand": row.ligand,
                "receptor": row.receptor, "n_donors_scored": present.donor.nunique(),
                "eligible_donors": len(donors), "required_donors": int(np.ceil(len(donors) / 2)),
                "median_rank_before_coverage_floor": present["rank"].median() if len(present) else np.nan,
                "retained_in_ranking": pair in retained})
        for path in [winner_path, all_path,
                     historical / "c14_does_the_ranking_depend_on_the_database" / f"c14_partial_{resource}.csv"]:
            hashes[str(path.relative_to(ROOT))] = sha256(path)
        if resource == "cellchatdb":
            shared = old.merge(winners, on=["donor", "pair", "source", "target"], suffixes=("_old", "_new"))
            checks["historical_shared_cellchat_winners"] = {
                "n": len(shared), "max_score_absolute_difference": float(abs(shared.lr_probs_old - shared.lr_probs_new).max()),
                "score_spearman": float(shared.lr_probs_old.corr(shared.lr_probs_new, method="spearman")),
                "scope": "Only historical winners retaining the same allowed source and target; not proof of biological validity"}
            finite = shared.loc[shared.lr_probs_old.between(1e-15, 1-1e-15)
                & shared.lr_probs_new.between(1e-15, 1-1e-15)].copy()
            finite["odds_ratio"] = ((finite.lr_probs_new / (1-finite.lr_probs_new))
                                   / (finite.lr_probs_old / (1-finite.lr_probs_old)))
            scales = finite.groupby("donor").odds_ratio.agg(["min", "max", "median"])
            spread = float(((scales["max"] - scales["min"]) / scales["median"]).max())
            assert spread < 1e-8, "shared scores changed beyond a common donor-level scale"
            checks["historical_shared_cellchat_winners"]["max_within_donor_odds_scale_relative_range"] = spread
            checks["historical_shared_cellchat_winners"]["scale_note"] = (
                "Shared score changes are a common within-donor rescaling in probability-odds space; "
                "LIANA scales by the maximum expression of the scoped cells. This preserves their relative ordering.")
    pd.DataFrame(receptor_coverage).to_csv(lr / "canonical_egfr_receptor_coverage.csv", index=False)
    table = pd.read_csv(c37 / "per_donor.csv")
    summary = json.loads((c37 / "summary.json").read_text())
    for field in ["original_detection", "raw_detection", "detection_at_500", "detection_at_1000", "detection_at_2000"]:
        recomputed = paired_summary(table, field)
        assert recomputed["n_donors"] == summary[field]["n_donors"]
        assert recomputed["epithelial_higher"] == summary[field]["epithelial_higher"]
        for stat in ["median_epithelial", "median_myeloid", "median_paired_difference", "p_two_sided"]:
            np.testing.assert_allclose(recomputed[stat], summary[field][stat], rtol=1e-10, atol=1e-12)
    for budget in [500, 1000, 2000]:
        assert (table[f"n_at_{budget}"] <= table.n_cells).all()
        assert table[f"detection_at_{budget}"].notna().equals(table[f"n_at_{budget}"] >= 50)
        assert table[f"detection_at_{budget}"].dropna().between(0, 1).all()
    checks["c37_summary_recomputed"] = True
    checks["c37_original_raw_detection_disagreements"] = summary["cache_raw_detection_disagreements"]
    checks["c37_cells_excluded_by_budget"] = {str(b): int((table.n_cells - table[f"n_at_{b}"]).sum())
                                            for b in [500, 1000, 2000]}
    for path in list(lr.glob("*.csv")) + [c37 / "per_donor.csv", c37 / "summary.json"]:
        hashes[str(path.relative_to(ROOT))] = sha256(path)
    checks["artifacts_sha256"] = hashes
    write_json(base / "verification.json", checks)
    print(json.dumps({k: v for k, v in checks.items() if k != "artifacts_sha256"}, indent=2))


if __name__ == "__main__":
    main()
