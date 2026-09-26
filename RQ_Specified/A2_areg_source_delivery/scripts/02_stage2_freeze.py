"""A2 stage 2: freeze the test. Reads stage 1 outputs and covariates only.

This stage computes no endpoint. It writes the machine-readable freeze that stage 3
must obey, including the exact null for every declared reference set, so that no
threshold is hand-typed. It refuses to overwrite an existing freeze.

Three changes to the original contract are recorded here with their reasons. All
three follow from stage 1 covariates, none from an endpoint value.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
A10 = ROOT / "RQ_Specified/A10_organoid_growth_outcome"
OUT = HERE / "config"
TABLES = HERE / "tables"
CONTRACT = HERE / "config/a2_delivery_contract.json"
FREEZE = OUT / "a2_stage2_freeze.json"
QC = A10 / "cache/GSE307112_xenome_stats.csv.gz"
TOTALS = A10 / "tables/library_totals.tsv"

DEPTH_FLOOR = 100_000
MATCHED_NEIGHBOURS = 20
ALPHA = 0.05
AXIS = ["AREG", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ITGB6"]
DROPPED = {"ERBB4": "not expressed in the perturbed compartment: 0.367 log2 CPM when not targeted and exactly zero when targeted, so there is no receptor to remove"}
LOW_ABUNDANCE = {"EGFR": "1.967 log2 CPM when not targeted, about 2.9 counts per million, so this contrast is weaker than Erbb2 or Erbb3"}


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda: h.read(1 << 24), b""):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace("\\", "/")


def target_of(library: str) -> str:
    parts = library.split("_", 2)
    return parts[2].upper() if len(parts) > 2 else ""


def rank_sum_null(sizes: list[int], alpha: float) -> dict:
    """Exact one-sided lower null for the sum of within-unit ranks, sizes may differ."""
    dist = {0: 1}
    for n in sizes:
        nxt: dict[int, int] = {}
        for total, count in dist.items():
            for r in range(1, n + 1):
                nxt[total + r] = nxt.get(total + r, 0) + count
        dist = nxt
    space = 1
    for n in sizes:
        space *= n
    cumulative = 0
    tail = {}
    for s in sorted(dist):
        cumulative += dist[s]
        tail[s] = cumulative / space
    critical = max((s for s in sorted(dist) if tail[s] <= alpha), default=None)
    smallest = tail[min(dist)]
    return {
        "unit_sizes": sizes,
        "rank_space": space,
        "smallest_attainable_one_sided_p": smallest,
        "smallest_attainable_p_scientific": "%.3g" % smallest,
        "critical_rank_sum": critical,
        "p_at_critical": round(tail[critical], 5) if critical is not None else None,
        "mean_percentile_required": round(sum(critical / n for n in sizes) / len(sizes) / len(sizes), 4)
        if critical is not None else None,
    }


def joint_null(sizes: list[int], critical: int, min_below: int) -> float:
    """Exact size of the rule: rank sum at or below critical, and below median in at least min_below units."""
    dp = {(0, 0): 1}
    for n in sizes:
        half = n // 2
        nxt: dict[tuple[int, int], int] = {}
        for (below, total), count in dp.items():
            for r in range(1, n + 1):
                key = (below + (1 if r <= half else 0), total + r)
                nxt[key] = nxt.get(key, 0) + count
        dp = nxt
    space = 1
    for n in sizes:
        space *= n
    good = sum(c for (below, total), c in dp.items() if below >= min_below and total <= critical)
    return good / space


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.parse_args()
    if FREEZE.exists():
        raise SystemExit("Refusing to overwrite the existing freeze: %s" % FREEZE.name)
    required = ["stage1_run.json", "stage1_layout.tsv", "stage1_knockout_axis.tsv",
                "stage1_endpoint_coverage.tsv", "stage1_covariates.tsv", "stage1_power.json"]
    missing = [n for n in required if not (TABLES / n).exists()]
    if missing:
        raise SystemExit("stage 1 has not run: missing %s" % missing)
    stage1 = json.loads((TABLES / "stage1_run.json").read_text(encoding="utf-8"))
    if stage1.get("endpoint_scored") is not False:
        raise SystemExit("stage 1 record does not assert that no endpoint was scored")
    if not stage1["item1_layout"]["matches_the_contract"]:
        raise SystemExit("stage 1 layout did not match the contract")
    if stage1["item5_precedent"]["verdict"].startswith("NOT blind"):
        raise SystemExit("stage 1 found the test is not blind; this freeze must be re-specified in the open")

    depth: dict[str, int] = {}
    epi_fraction: dict[str, float] = {}
    with open(TOTALS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            mouse, human = int(r["mouse_total_counts"]), int(r["human_total_counts"])
            depth[r["library"]] = human
            epi_fraction[r["library"]] = mouse / (mouse + human) if (mouse + human) else 0.0
    with gzip.open(QC, "rt", encoding="utf-8-sig", newline="") as fh:
        plate3 = [r for r in csv.DictReader(fh) if r["plate"] == "plate3"]
    by_unit: dict[str, list[str]] = {}
    for r in plate3:
        by_unit.setdefault(r["plate_rep"], []).append(r["library name"])
    units = sorted(by_unit)

    full_sizes = [len(by_unit[u]) for u in units]
    floor_sizes = [sum(1 for n in by_unit[u] if depth[n] >= DEPTH_FLOOR) for u in units]
    matched_sizes = [min(MATCHED_NEIGHBOURS, len(by_unit[u])) for u in units]

    areg_wells = {}
    for u in units:
        hit = [n for n in by_unit[u] if target_of(n) == "AREG"]
        if len(hit) != 1:
            raise SystemExit("unit %s does not hold exactly one Areg well" % u)
        areg_wells[u] = {"library": hit[0], "fibroblast_total_counts": depth[hit[0]],
                         "clears_depth_floor": depth[hit[0]] >= DEPTH_FLOOR,
                         "epithelial_count_fraction": round(epi_fraction[hit[0]], 4)}
    if not all(v["clears_depth_floor"] for v in areg_wells.values()):
        raise SystemExit("an Areg well is below the declared depth floor; re-specify before freezing")

    primary = rank_sum_null(full_sizes, ALPHA)
    primary["exact_size_with_consistency_requirement"] = round(
        joint_null(full_sizes, primary["critical_rank_sum"], len(full_sizes) - 1), 5)
    s1 = rank_sum_null(floor_sizes, ALPHA)
    s2 = rank_sum_null(matched_sizes, ALPHA)

    below_floor = []
    for u in units:
        for n in by_unit[u]:
            t = target_of(n)
            if t in AXIS + ["TDTOMATO"] and depth[n] < DEPTH_FLOOR:
                below_floor.append({"unit": u, "target": t, "library": n, "fibroblast_total_counts": depth[n]})

    freeze = {
        "schema": "a2-stage2-freeze/v1",
        "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "authorizes": "stage 3, leg 1 only; leg 2 keeps its own declarations in the contract",
        "contract": {"path": rel(CONTRACT), "sha256": sha256(CONTRACT)},
        "stage1": {"run_record_sha256": sha256(TABLES / "stage1_run.json"),
                   "endpoint_scored": stage1["endpoint_scored"]},
        "inputs": {rel(TOTALS): sha256(TOTALS), rel(QC): sha256(QC)},
        "endpoints_unchanged_from_the_contract": {
            "primary": "HALLMARK_TGF_BETA_SIGNALING, human, 54 of 54 members present",
            "co_primary": "fibroblast activation score COL1A1, ACTA2, POSTN, CTHRC1, TNC, 5 of 5 present",
            "scoring": "mean over members of log2(count / the well's human total * 1e6 + 1), the same prior and denominator convention A10 used",
            "co_primary_role": "must agree in direction; the primary carries the decision",
        },
        "reference_set": {
            "primary": "all wells of the Areg well's own unit, %s per unit" % full_sizes,
            "no_exclusion_in_the_primary": "excluding wells changes the rank universe, so the primary keeps every well and the depth restrictions are declared sensitivities instead",
        },
        "adjustment": {
            "covariates": ["log2 of the well's fibroblast total counts", "epithelial count fraction"],
            "epithelial_count_fraction_definition": "mouse total counts divided by the sum of mouse and human total counts, from the A10 library totals table",
            "clarification_of_the_contract": "the contract said epithelial read fraction; this freeze fixes the count-based fraction, because it comes from the same counts the endpoint uses and is reproducible from a tracked table without re-deriving the species assignment",
            "method": "within each unit, ordinary least squares of the endpoint on an intercept and the two covariates, across that unit's wells; the residual is the analysed quantity",
            "why": "stage 1 found fibroblast depth spanning four orders of magnitude and all four Areg wells above their unit median, so the depth control is load-bearing rather than decorative; C51 is the precedent",
        },
        "primary_test": {
            "statistic": "the ascending rank of the Areg well's residual among its unit's residuals, summed over the four units",
            "direction": "one-sided lower, declared before any value is read",
            "alpha": ALPHA,
            "null": primary,
            "decision_threshold": "reject only if the rank sum is at or below %d and the Areg well is below its unit median in at least %d of %d units"
                                  % (primary["critical_rank_sum"], len(full_sizes) - 1, len(full_sizes)),
            "exact_size": primary["exact_size_with_consistency_requirement"],
            "what_a_rejection_means": "the Areg well is exceptional among the perturbations in its own unit, in the declared direction; it is not a comparison against an unperturbed population, and it is conservative when other wells in the unit act on the same pathway",
        },
        "declared_sensitivities": {
            "S1_depth_floor": {"definition": "reference set restricted to wells with at least %d fibroblast counts" % DEPTH_FLOOR,
                               "unit_sizes": floor_sizes, "null": s1,
                               "why": "a mean over 54 genes in a shallow well reports detection rather than programme state"},
            "S2_depth_matched": {"definition": "reference set restricted to the %d wells nearest the Areg well in log fibroblast depth within its unit, the Areg well included" % MATCHED_NEIGHBOURS,
                                 "unit_sizes": matched_sizes, "null": s2,
                                 "why": "the Areg wells are systematically deeper than their units, so a depth-matched reference tests the result against that asymmetry directly"},
            "S3_unadjusted": {"definition": "the same statistic on the raw endpoint, with no residualization",
                              "role": "sensitivity only, never the primary"},
            "requirement_for_a_supported_reading": "the primary must clear its threshold and the effect must hold in direction in at least one of S1 and S2, because the depth asymmetry is known in advance",
        },
        "declared_diagnostic": "the within-unit association between the endpoint and log fibroblast depth, reported next to the result; stage 1 reasoned that a depth artefact would push the score upward, which opposes the declared direction, and stage 3 must show whether that holds",
        "discriminating_contrasts": {
            "decision_attached": False,
            "within_plate_3": ["EGFR", "ERBB2", "ERBB3", "ITGB6"],
            "across_plates": ["HBEGF, on plate 4, labelled weaker for that reason"],
            "dropped": DROPPED,
            "low_abundance_flag": LOW_ABUNDANCE,
            "uninterpretable_wells": {"floor": DEPTH_FLOOR, "wells": below_floor,
                                      "rule": "a contrast whose well is below the floor is reported as uninterpretable, not as evidence either way"},
        },
        "control_wells": {
            "role_downgraded": True,
            "reason": "four of the eight plate-3 control wells hold fewer than %d fibroblast counts, so they do not define the unperturbed location on this plate" % DEPTH_FLOOR,
            "new_role": "descriptive context only",
        },
        "areg_wells": areg_wells,
        "decision_rules_unchanged": {
            "supported": "the primary clears its threshold and its consistency requirement, the co-primary agrees in direction, at least one depth sensitivity holds in direction, and the receptor contrasts do not show the same shift",
            "inconclusive": "any other outcome",
            "precise_absence": "unavailable at this design, declared before any test",
            "descriptive_only": "always, while preparation independence is unresolved",
        },
        "changes_from_the_contract": [
            {"change": "Erbb4 dropped from the discriminating set", "reason": DROPPED["ERBB4"], "source": "stage 1 item 2"},
            {"change": "Egfr flagged as low-abundance", "reason": LOW_ABUNDANCE["EGFR"], "source": "stage 1 item 2"},
            {"change": "a fibroblast depth floor for reading single wells, plus two declared depth sensitivities and a supported-reading requirement", "reason": "fibroblast depth spans four orders of magnitude and the Areg wells are systematically deeper than their units", "source": "stage 1 item 4"},
            {"change": "the control wells' anchor role downgraded to descriptive", "reason": "half of the plate-3 control wells are shallow", "source": "stage 1 item 4"},
            {"change": "epithelial read fraction fixed as the count-based epithelial fraction", "reason": "reproducible from a tracked table and derived from the same counts as the endpoint", "source": "this freeze"},
        ],
        "unchanged": ["the hypothesis", "both endpoints", "the declared direction", "the unit", "the prohibitions", "precise absence remaining unavailable"],
        "endpoint_scored": False,
    }
    FREEZE.write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")
    print("froze the test in %s" % rel(FREEZE))
    print("primary: unit sizes %s, critical rank sum %s, exact size %s, smallest attainable p %s"
          % (full_sizes, primary["critical_rank_sum"], freeze["primary_test"]["exact_size"],
             primary["smallest_attainable_p_scientific"]))
    print("S1 floor: unit sizes %s, critical %s, smallest p %s" % (floor_sizes, s1["critical_rank_sum"], s1["smallest_attainable_p_scientific"]))
    print("S2 matched: unit sizes %s, critical %s, smallest p %s" % (matched_sizes, s2["critical_rank_sum"], s2["smallest_attainable_p_scientific"]))
    print("uninterpretable axis or control wells below the floor: %d" % len(below_floor))


if __name__ == "__main__":
    main()
