"""A2 stage 2, second freeze: withdraw the rank test and declare what the design supports.

The first freeze (config/a2_stage2_freeze.json) is preserved unchanged and withdrawn.
Three facts, each verified from the deposit or from this repository's own tables,
remove the inference it declared:

1. The recipient makes the ligand. Stage 1b records human AREG at 9.846 mean log2 CPM
   in 99.2 per cent of plate-3 wells, above the mouse epithelial Areg it would
   replace within that compartment's own denominator. Removing the epithelial source
   does not remove AREG from the culture, so no necessity claim is available.
2. There is no randomization. Fifty of the 53 plate-3 targets occupy one fixed well
   position in all four units, Areg always at F07, so target is confounded with plate
   position and guide pool and the four units are copies of one layout. A
   uniform-rank null has nothing behind it.
3. The reference set is not neutral. Plate 3 carries EGFR, ERBB2, ERBB3, ERBB4, AREG
   and the MAP kinase cascade, so ranking Areg among its plate-mates asks whether it
   is extreme among other perturbations of the same pathway.

A fourth fact improves what remains: the screen's own in-plate control is TIGIT, at
six wells per unit on every plate, which A10's source design check already recorded.
With TDTOMATO that gives eight control wells per unit instead of two.

This freeze therefore declares an effect size against the screen's controls with a
direction count, and no p-value, which is what the register card promised all along.
No endpoint is scored here.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
A10 = ROOT / "RQ_Specified/A10_organoid_growth_outcome"
OUT = HERE / "config"
TABLES = HERE / "tables"
V1 = OUT / "a2_stage2_freeze.json"
V2 = OUT / "a2_stage2_freeze_v2.json"
CONTRACT = OUT / "a2_delivery_contract.json"
QC = A10 / "cache/GSE307112_xenome_stats.csv.gz"
TOTALS = A10 / "tables/library_totals.tsv"

ELIGIBILITY_FLOOR = 100_000
CONTROLS = ["TIGIT", "TDTOMATO"]
DEPTH_BAND = 4.0


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


def main() -> None:
    argparse.ArgumentParser().parse_args()
    if V2.exists():
        raise SystemExit("Refusing to overwrite %s" % V2.name)
    if not V1.exists():
        raise SystemExit("the first freeze is missing; it must be preserved, not replaced")
    for name in ["stage1_run.json", "stage1b_run.json", "stage1b_fibroblast_covariates.tsv"]:
        if not (TABLES / name).exists():
            raise SystemExit("missing %s; stage 1 and its addendum must both have run" % name)

    addendum = json.loads((TABLES / "stage1b_run.json").read_text(encoding="utf-8"))
    if addendum.get("endpoint_scored") is not False:
        raise SystemExit("the addendum record does not assert that no endpoint was scored")

    depth: dict[str, int] = {}
    epi: dict[str, float] = {}
    with open(TOTALS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            m, h = int(r["mouse_total_counts"]), int(r["human_total_counts"])
            depth[r["library"]] = h
            epi[r["library"]] = m / (m + h) if (m + h) else 0.0
    with gzip.open(QC, "rt", encoding="utf-8-sig", newline="") as fh:
        rows = [r for r in csv.DictReader(fh) if r["plate"] == "plate3"]

    by_unit: dict[str, list[str]] = {}
    positions: dict[str, set] = {}
    for r in rows:
        by_unit.setdefault(r["plate_rep"], []).append(r["library name"])
        positions.setdefault(r["crispr_target"].upper(), set()).add(r["well"])
    fixed = sum(1 for t, w in positions.items() if len(w) == 1)

    units = sorted(by_unit)
    design = {}
    for u in units:
        names = by_unit[u]
        areg = [n for n in names if target_of(n) == "AREG"][0]
        ctrl = [n for n in names if target_of(n) in CONTROLS]
        ctrl_ok = [n for n in ctrl if depth[n] >= ELIGIBILITY_FLOOR]
        lo, hi = depth[areg] / DEPTH_BAND, depth[areg] * DEPTH_BAND
        ctrl_band = [n for n in ctrl_ok if lo <= depth[n] <= hi]
        design[u] = {
            "wells": len(names),
            "eligible_wells": sum(1 for n in names if depth[n] >= ELIGIBILITY_FLOOR),
            "areg_library": areg,
            "areg_fibroblast_counts": depth[areg],
            "areg_epithelial_count_fraction": round(epi[areg], 4),
            "control_wells": len(ctrl),
            "control_wells_eligible": len(ctrl_ok),
            "control_fibroblast_counts": sorted(depth[n] for n in ctrl),
            "control_wells_within_depth_band": len(ctrl_band),
        }
    min_ctrl = min(v["control_wells_eligible"] for v in design.values())
    min_band = min(v["control_wells_within_depth_band"] for v in design.values())

    # Which targets are readable at all, per unit, at the eligibility floor.
    clearance = {}
    for t in ["AREG", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ITGB6"]:
        ok = 0
        counts = []
        for u in units:
            hit = [n for n in by_unit[u] if target_of(n) == t]
            if not hit:
                continue
            counts.append(depth[hit[0]])
            if depth[hit[0]] >= ELIGIBILITY_FLOOR:
                ok += 1
        clearance[t] = {"units_clearing_the_floor": ok, "units_present": len(counts),
                        "fibroblast_counts_by_unit": counts}

    freeze = {
        "schema": "a2-stage2-freeze/v2",
        "frozen_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "supersedes": {"file": rel(V1), "sha256": sha256(V1),
                       "status": "withdrawn, preserved unchanged",
                       "report": "reports/STAGE2_WITHDRAWN.md"},
        "contract": {"path": rel(CONTRACT), "sha256": sha256(CONTRACT)},
        "stage1": {"run_record_sha256": sha256(TABLES / "stage1_run.json"),
                   "addendum_sha256": sha256(TABLES / "stage1b_run.json")},
        "why_the_first_freeze_is_withdrawn": [
            {"reason": "the recipient makes the ligand",
             "evidence": "stage 1b: human AREG at 9.846 mean log2 CPM, detected in 99.2 per cent of plate-3 wells, against mouse epithelial Areg at 7.267 in its own compartment",
             "consequence": "removing the epithelial source does not remove AREG from the culture, so neither necessity nor sufficiency is available from this design"},
            {"reason": "there is no randomization",
             "evidence": "%d of %d plate-3 targets occupy one fixed well position in all four units, Areg always at F07" % (fixed, len(positions)),
             "consequence": "target is confounded with plate position and guide pool, the four units are copies of one layout, and a uniform-rank null has nothing behind it"},
            {"reason": "the reference set is not neutral",
             "evidence": "plate 3 carries AREG, EGFR, ERBB2, ERBB3, ERBB4 and the MAP kinase cascade among its 53 targets",
             "consequence": "a rank among plate-mates asks whether Areg is extreme among other perturbations of the same pathway, which is not the declared estimand"},
            {"reason": "the denominator held ineligible wells",
             "evidence": "plate-3 fibroblast totals run from 1 count upward, and the first freeze applied its floor only to sensitivities",
             "consequence": "some ranks were deterministic rather than random draws"},
            {"reason": "the consistency requirement was arithmetically vacuous",
             "evidence": "the joint rule's exact size equalled the rank-sum size, 0.04903",
             "consequence": "it added no stringency while reading as if it did"},
        ],
        "what_replaces_it": {
            "inference": "an effect size with a direction count, descriptive, with no p-value, no alpha and no significance claim",
            "why": "the register card already promised an effect size and a consistency count rather than a significance test, and the design supports exactly that",
        },
        "endpoints": {
            "primary": {"name": "fibroblast activation score",
                        "genes": ["COL1A1", "ACTA2", "POSTN", "CTHRC1", "TNC"],
                        "promoted_from": "co-primary in the first freeze",
                        "why": "it is the readout both cited mechanism papers used, alpha-smooth-muscle actin and collagen, so it is the mechanism-matched endpoint"},
            "secondary": {"name": "HALLMARK_TGF_BETA_SIGNALING",
                          "demoted_from": "primary in the first freeze",
                          "why": "a pathway-membership set whose 54 members include 16 negative regulators and feedback genes, and which contains none of the five activation genes, so its direction is mixed and it is not what the cited papers measured"},
            "scoring": "mean over members of log2(count / the well's human total * 1e6 + 1)",
            "inherited": "the gene list only, chosen before this screen was touched; neither the E6 operator nor the C51 precedent transfers to a per-well log2 CPM score",
        },
        "eligibility": {
            "floor_human_counts": ELIGIBILITY_FLOOR,
            "justification": "at this depth a gene at 50 counts per million is expected about five times, so a five-gene or 54-gene mean estimates a programme rather than reporting detection",
            "applies_to": "every well in every comparison, including the controls",
            "eligible_wells_per_unit": [design[u]["eligible_wells"] for u in units],
        },
        "control_set": {
            "targets": CONTROLS,
            "why": "TIGIT is the screen's own in-plate control, recorded in A10's source design check, at six wells per unit; the first freeze used TDTOMATO alone",
            "wells_per_unit": [design[u]["control_wells"] for u in units],
            "eligible_per_unit": [design[u]["control_wells_eligible"] for u in units],
            "minimum_eligible_controls_in_any_unit": min_ctrl,
        },
        "primary_comparison": {
            "statistic": "within each unit, the Areg well's endpoint minus the mean endpoint of that unit's eligible control wells, in log2 CPM units",
            "reported": "one effect per unit, the four effects, their median, and the count of units in the predicted direction",
            "predicted_direction": "lower in the Areg well, declared before any value is read",
            "no_p_value": True,
            "why_no_p_value": "there is no randomization to support one, and the four units are copies of one layout",
            "depth_matched_variant": {
                "definition": "controls restricted to those within a factor of %g of the Areg well's fibroblast total" % DEPTH_BAND,
                "controls_within_band_per_unit": [design[u]["control_wells_within_depth_band"] for u in units],
                "minimum_in_any_unit": min_band,
                "why": "the Areg wells are deep and several control wells are shallow, so the unmatched contrast confounds depth with target",
                "status": "required alongside the unmatched contrast, not optional",
            },
            "epithelial_fraction": {
                "role": "reported as a two-sided sensitivity, not a primary covariate",
                "why": "epithelial abundance plausibly lies on the causal path from an Areg knockout to the fibroblast read, so adjusting for it removes part of the effect; about a fifth of reads are also unassigned to either species",
                "pre_declared_reading": "an effect that survives only without the epithelial adjustment, or only with it, is reported as mediation-ambiguous rather than as artefact or as signal",
            },
        },
        "discriminating_contrasts": {
            "decision_attached": False,
            "no_threshold": "no numeric threshold is declared, because none is defensible at one well per target per unit; these contrasts are descriptive",
            "areg_receptor": ["EGFR"],
            "heterodimer_partner": ["ERBB2"],
            "not_areg_receptors": {"ERBB3": "binds neuregulins, not AREG",
                                   "ERBB4": "binds neuregulins, HB-EGF, betacellulin and epiregulin, not AREG; also unexpressed here"},
            "reclassified_as": "non-AREG-receptor perturbation controls, not tests of AREG reception",
            "itgb6": "epithelial integrin contrast, descriptive",
            "hbegf": "dropped as an effect-size comparison; it sits on plate 4 and the inherited measurement contract forbids treating cross-plate ranks as commensurable effect sizes",
            "egfr_caveat": "epithelial Egfr is 1.967 log2 CPM when not targeted, near the detection floor, so the autocrine rival may not be testable in this screen at all",
            "readability_at_the_floor": clearance,
            "readability_note": "Erbb2 is the weakest contrast by readability, not Egfr: it clears the fibroblast eligibility floor in one unit of four, while Egfr clears three and Erbb3 four. The first freeze ranked these the other way round.",
            "no_null_supports_absence": "every recorded reduction in this axis is partial, so no arm's null can exclude a role, and the receptor arm contributes only as a consistency check",
        },
        "hard_limits_stated_before_any_score": [
            "The recipient expresses the ligand, so leg 1 can at most bound the incremental contribution of the epithelial source on top of an unremoved autocrine source.",
            "Ligand redundancy is not excluded: stage 1b records fibroblast EREG at 7.659 log2 CPM in 96.3 per cent of wells, and no Ereg, Tgfa or Btc well exists in this screen.",
            "Delivery is not separable from abundance in a single well with one source compartment and no spatial variation; that contrast belongs to the spatial layer.",
            "Neither leg is blind in the strict sense. The primary endpoint's inputs are already cached from A10's extraction, and leg 2's outcome is already computed and graded (C50). What is new in each leg is the contrast, not the measurement.",
            "Preparation independence is unresolved, so every result is a within-screen association.",
            "The cross-species assumption that mouse amphiregulin activates human EGFR is carried and unverified.",
            "The culture medium and matrix are not in the deposit, so exogenous EGF or a TGF-beta receptor inhibitor cannot be ruled out; either would mask or clamp the endpoint, and a null is pharmacologically ambiguous until the methods are checked.",
        ],
        "leg_2_reclassified": {
            "status": "exploratory, not a declared test",
            "why": "the inherited depth gate is likely to refuse it: the logged E6 table already records the activation score against fibroblast depth at rho 0.4116 and fibroblast EGFR at 0.4918, and both members of the new pair are fibroblast detection fractions",
            "receptor_layer_omitted": "the cited mechanism makes recipient EGFR obligatory, and leg 2 does not measure it because A9 owns the receptor layer; the estimand is therefore the post-receptor integrin and latent-complex layer, not recipient licensing as a whole",
            "recipient_in_the_cited_work": "a PDGFRB-positive pericyte, not a lung fibroblast, so extension to fibroblasts is an assumption carried on the same footing as the cross-species one; a mural stratum is reported beside the fibroblast stratum rather than instead of it",
            "controls_required": [
                "the C51 depth rule, unchanged",
                "a co-regulation control: genes matched to the machinery set on detection frequency in the same fibroblasts and not TGF-beta responsive, since TGFB1, THBS1, ITGB8 and LTBP1 are themselves TGF-beta inducible",
                "a per-gene expression eligibility floor, as MC4 requires",
                "the IPF-only stratum beside the pooled result",
            ],
            "named_rival": "co-regulation within one activation programme, which a positive result cannot exclude",
        },
        "design_facts": design,
        "fixed_position_targets": {"targets_with_one_well_position": fixed, "targets_total": len(positions)},
        "endpoint_scored": False,
    }
    V2.write_text(json.dumps(freeze, indent=2) + "\n", encoding="utf-8")
    print("wrote %s" % rel(V2))
    print("eligible wells per unit: %s" % [design[u]["eligible_wells"] for u in units])
    print("control wells per unit %s, eligible %s, within depth band %s"
          % ([design[u]["control_wells"] for u in units],
             [design[u]["control_wells_eligible"] for u in units],
             [design[u]["control_wells_within_depth_band"] for u in units]))
    print("targets at one fixed well position: %d of %d" % (fixed, len(positions)))
    for u in units:
        print("  %s Areg %d counts, controls %s" % (u, design[u]["areg_fibroblast_counts"], design[u]["control_fibroblast_counts"]))


if __name__ == "__main__":
    main()
