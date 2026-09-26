"""A2 stage 3 robustness, post hoc: is the Itgb6 signal a viability or composition effect.

Declared post hoc and labelled as such. The second freeze did not declare this check;
it is a confound inspection prompted by the stage 3 result, and it carries no decision.

The question: the epithelial Itgb6 knockout lowers the human fibroblast activation
score by about 0.94 log2 CPM in every unit where its well is readable. If that well
also holds far fewer fibroblast reads, far less epithelium, or a collapsed organoid
culture, the score change could be a viability or composition artefact rather than a
signalling one. Depth was already matched in the primary, so this adds the deposited
imaging outcome and the compartment totals.

A10 owns the organoid growth outcome. Nothing here fits a growth model, re-reads A10's
models or makes a growth claim; the imaging columns are used only as a confound check.

Standard library only. Refuses to overwrite.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import statistics
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
A10 = ROOT / "RQ_Specified/A10_organoid_growth_outcome"
OUT = HERE / "tables"
IMAGING = A10 / "cache/GSE307112_imaging_outputs.csv.gz"
SCORES = OUT / "stage3_well_scores.tsv"
OUTPUT = "stage3_robustness.tsv"
RECORD = "stage3_robustness_run.json"
CONTROLS = {"TIGIT", "TDTOMATO"}
AXIS = ["AREG", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ITGB6"]
FLOOR = 100_000


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda: h.read(1 << 24), b""):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace("\\", "/")


def main() -> None:
    argparse.ArgumentParser().parse_args()
    if (OUT / OUTPUT).exists():
        raise SystemExit("Refusing to overwrite %s" % OUTPUT)
    if not SCORES.exists():
        raise SystemExit("stage 3 has not run")

    wells = list(csv.DictReader(open(SCORES, encoding="utf-8"), delimiter="\t"))
    by_lib = {w["library"]: w for w in wells}

    imaging: dict[tuple, dict] = {}
    with gzip.open(IMAGING, "rt", encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            if r["plate"] != "plate3" or r["day"] != "day14":
                continue
            imaging[(r["plate_replicate"], r["well"].upper())] = r

    def imaging_for(library: str) -> dict | None:
        """The deposit keys imaging as plate3-1 where a library name says 3-1."""
        parts = library.split("_")
        if len(parts) < 3:
            return None
        unit, well = parts[0], parts[1].upper()
        return imaging.get(("plate" + unit, well))

    rows = []
    for unit in sorted({w["unit"] for w in wells}):
        here = [w for w in wells if w["unit"] == unit]
        ctrl = [w for w in here if w["target"] in CONTROLS and w["eligible"] == "yes"]
        ctrl_fib = [int(w["fibroblast_total_counts"]) for w in ctrl]
        ctrl_img = [imaging_for(w["library"]) for w in ctrl]
        ctrl_area = [float(i["organoids_area_mean"]) for i in ctrl_img if i and i["organoids_area_mean"]]
        ctrl_count = [float(i["organoids_count"]) for i in ctrl_img if i and i["organoids_count"]]
        for target in AXIS:
            hit = [w for w in here if w["target"] == target]
            if not hit:
                continue
            w = hit[0]
            img = imaging_for(w["library"])
            rows.append({
                "unit": unit, "target": target, "library": w["library"],
                "eligible": w["eligible"],
                "activation_score": w["activation_score"],
                "fibroblast_total_counts": w["fibroblast_total_counts"],
                "control_median_fibroblast_counts": int(statistics.median(ctrl_fib)) if ctrl_fib else None,
                "fibroblast_counts_ratio_to_control": round(int(w["fibroblast_total_counts"]) / statistics.median(ctrl_fib), 3) if ctrl_fib else None,
                "epithelial_count_fraction": w["epithelial_count_fraction"],
                "control_median_epithelial_fraction": round(statistics.median([float(c["epithelial_count_fraction"]) for c in ctrl]), 4) if ctrl else None,
                "day14_organoid_count": img["organoids_count"] if img else None,
                "control_median_day14_count": round(statistics.median(ctrl_count), 1) if ctrl_count else None,
                "day14_area_mean": img["organoids_area_mean"] if img else None,
                "control_median_day14_area": round(statistics.median(ctrl_area), 1) if ctrl_area else None,
                "day14_area_minus_control": round(float(img["organoids_area_mean"]) - statistics.median(ctrl_area), 3)
                if img and img["organoids_area_mean"] and ctrl_area else None,
                "day14_count_ratio_to_control": round(float(img["organoids_count"]) / statistics.median(ctrl_count), 3)
                if img and img["organoids_count"] and ctrl_count and statistics.median(ctrl_count) else None,
            })

    with open(OUT / OUTPUT, "w", encoding="utf-8", newline="") as fh:
        fields = list(rows[0].keys())
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("none" if r.get(k) is None else r.get(k)) for k in fields})

    summary = {}
    for target in AXIS:
        sel = [r for r in rows if r["target"] == target and r["eligible"] == "yes"]
        ar = [r["day14_area_minus_control"] for r in sel if r["day14_area_minus_control"] is not None]
        fr = [r["fibroblast_counts_ratio_to_control"] for r in sel if r["fibroblast_counts_ratio_to_control"] is not None]
        summary[target] = {
            "eligible_units": len(sel),
            "median_day14_area_minus_control": round(statistics.median(ar), 3) if ar else None,
            "median_fibroblast_count_ratio_to_control": round(statistics.median(fr), 3) if fr else None,
        }

    record = {
        "stage": "3 robustness, post hoc",
        "declared": "post hoc; not in the second freeze; no decision attached",
        "first_attempt": "tables/stage3_robustness_attempt1_failed_join/, preserved; it joined imaging on the wrong key prefix",
        "question": "whether the Itgb6 signal travels with a viability or composition change",
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": {rel(IMAGING): sha256(IMAGING), rel(SCORES): sha256(SCORES)},
        "a10_boundary": "the imaging columns are a confound check only; no growth model is fitted and no growth claim is made",
        "eligibility_floor": FLOOR,
        "summary": summary,
    }
    (OUT / RECORD).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print("post hoc robustness, eligible wells only, ratios against each unit's eligible controls:")
    for target in AXIS:
        s = summary[target]
        print("  %-6s units %d, day-14 log area minus control %s, fibroblast count ratio %s"
              % (target, s["eligible_units"], s["median_day14_area_minus_control"],
                 s["median_fibroblast_count_ratio_to_control"]))


if __name__ == "__main__":
    main()
