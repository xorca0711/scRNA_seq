"""A2 stage 3, leg 1: the epithelial source contribution as an effect size.

Obeys config/a2_stage2_freeze_v2.json and nothing else. There is no p-value, no
alpha and no significance claim: the freeze withdrew those because the four plate-3
units are copies of one layout and the recipient supplies the ligand.

What this computes, per plate-replicate unit:

* the primary endpoint, the five-gene fibroblast activation score, and the
  secondary, the Hallmark TGF-beta set, on the human reads of every plate-3 well;
* the Areg well's score minus the mean of that unit's eligible control wells
  (TIGIT plus TDTOMATO), unmatched and restricted to controls within a factor of
  four of the Areg well's fibroblast depth;
* the same contrast for the other axis targets, with no threshold and no decision;
* an epithelial-fraction sensitivity, read as mediation-ambiguous either way;
* the declared depth diagnostic.

Standard library only. Hash-verified inputs. Refuses to overwrite.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import re
import statistics
import sys
import zipfile
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
A10 = ROOT / "RQ_Specified/A10_organoid_growth_outcome"
CACHE = A10 / "cache"
OUT = HERE / "tables"
REPORTS = HERE / "reports"
FREEZE = HERE / "config/a2_stage2_freeze_v2.json"
WORKBOOK = CACHE / "GSE307112_gene_counts.xlsx"
QC = CACHE / "GSE307112_xenome_stats.csv.gz"
TOTALS = A10 / "tables/library_totals.tsv"
WORKBOOK_SHA = "60736e5c153a4fb92908f19148ebb5ef2b59f7ade97fb1056deee7d0205acc33"
HUMAN_SHEET = "xl/worksheets/sheet2.xml"
GMT_RELATIVE = "raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt"

ACTIVATION = ["COL1A1", "ACTA2", "POSTN", "CTHRC1", "TNC"]
SECONDARY_SET = "HALLMARK_TGF_BETA_SIGNALING"
PRIOR = 1.0
FLOOR = 100_000
BAND = 4.0
CONTROLS = {"TIGIT", "TDTOMATO"}
AXIS = ["AREG", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ITGB6"]
OUTPUTS = ["stage3_well_scores.tsv", "stage3_contrasts.tsv", "stage3_diagnostics.tsv", "stage3_run.json"]

CELL = re.compile(rb'<c r="([A-Z]+)\d+"[^>]*?>(?:<v>([^<]*)</v>|<is><t>([^<]*)</t></is>)</c>')
ROW = re.compile(rb"<row [^>]*?>(.*?)</row>", re.S)
LIBNAME = re.compile(r"^\d-\d_[A-Z]\d{2}_")


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda: h.read(1 << 24), b""):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    return str(p.relative_to(ROOT)).replace("\\", "/")


def col_index(ref: bytes) -> int:
    n = 0
    for ch in ref:
        n = n * 26 + (ch - 64)
    return n - 1


def target_of(library: str) -> str:
    parts = library.split("_", 2)
    return parts[2].upper() if len(parts) > 2 else ""


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        mean_rank = (i + j) / 2 + 1
        for k in range(i, j + 1):
            out[order[k]] = mean_rank
        i = j + 1
    return out


def pearson(x: list[float], y: list[float]) -> float:
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx <= 0 or syy <= 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def spearman(x: list[float], y: list[float]) -> float:
    return pearson(ranks(x), ranks(y))


def ols_residuals(y: list[float], x: list[float]) -> list[float]:
    """Residuals of y on an intercept and one regressor."""
    n = len(y)
    mx, my = sum(x) / n, sum(y) / n
    sxx = sum((a - mx) ** 2 for a in x)
    if sxx <= 0:
        return [b - my for b in y]
    slope = sum((a - mx) * (b - my) for a, b in zip(x, y)) / sxx
    intercept = my - slope * mx
    return [b - (intercept + slope * a) for a, b in zip(x, y)]


def write_tsv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("none" if r.get(k) in (None, "") else r.get(k)) for k in fields})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT)
    data_root = ap.parse_args().data_root.resolve()

    OUT.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite: %s" % existing)

    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    if freeze["schema"] != "a2-stage2-freeze/v2":
        raise SystemExit("wrong freeze schema: %s" % freeze["schema"])
    if freeze["primary_comparison"]["no_p_value"] is not True:
        raise SystemExit("the freeze must forbid a p-value")
    if freeze["eligibility"]["floor_human_counts"] != FLOOR:
        raise SystemExit("floor disagrees with the freeze")
    if sha256(WORKBOOK) != WORKBOOK_SHA:
        raise SystemExit("workbook hash differs from the recorded download")

    gmt = data_root / GMT_RELATIVE
    hallmark: list[str] = []
    for line in gmt.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if parts and parts[0] == SECONDARY_SET:
            hallmark = parts[2:]
    if len(hallmark) != 54:
        raise SystemExit("%s has %d members, expected 54" % (SECONDARY_SET, len(hallmark)))

    wanted = {g.upper() for g in ACTIVATION} | {g.upper() for g in hallmark}
    recorded_human, recorded_mouse = {}, {}
    with open(TOTALS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            recorded_human[r["library"]] = int(r["human_total_counts"])
            recorded_mouse[r["library"]] = int(r["mouse_total_counts"])
    with gzip.open(QC, "rt", encoding="utf-8-sig", newline="") as fh:
        unit_of = {r["library name"]: r["plate_rep"] for r in csv.DictReader(fh) if r["plate"] == "plate3"}

    print("streaming the human sheet for %d endpoint genes" % len(wanted), flush=True)
    zf = zipfile.ZipFile(WORKBOOK)
    keep: dict[str, dict[int, int]] = {}
    totals: dict[int, int] = {}
    order: dict[int, str] = {}
    header, sym_col, seen = None, None, 0
    with zf.open(HUMAN_SHEET) as handle:
        buffer, tail = b"", False
        while True:
            chunk = handle.read(1 << 23)
            if not chunk:
                tail = True
            buffer += chunk
            cut = buffer.rfind(b"</row>")
            if tail:
                block, buffer = buffer, b""
            elif cut >= 0:
                block, buffer = buffer[: cut + 6], buffer[cut + 6 :]
            else:
                block = b""
            for m in ROW.finditer(block):
                cells = CELL.findall(m.group(1))
                if header is None:
                    header = {col_index(c[0]): (c[2] or c[1]).decode() for c in cells}
                    sym_col = next(i for i, v in header.items() if v == "symbol")
                    order = {i: v for i, v in header.items() if LIBNAME.match(v)}
                    totals = {i: 0 for i in order}
                    continue
                seen += 1
                symbol, values = None, {}
                for ref, num, txt in cells:
                    i = col_index(ref)
                    if i == sym_col:
                        symbol = txt.decode()
                    elif i in totals and num:
                        v = int(num) if b"." not in num else int(float(num))
                        totals[i] += v
                        if v:
                            values[i] = v
                if symbol and symbol.upper() in wanted:
                    keep.setdefault(symbol.upper(), {}).update(values)
                if seen % 20000 == 0:
                    print("   human: %d rows" % seen, flush=True)
            if tail:
                break

    cols = sorted(order)
    names = [order[i] for i in cols]
    mismatch = [n for i, n in zip(cols, names) if totals[i] != recorded_human.get(n)]
    if mismatch:
        raise SystemExit("STOP: recomputed human totals disagree with the A10 record for %d libraries" % len(mismatch))
    print("human totals cross-check against A10: 0 mismatches in %d libraries" % len(names), flush=True)

    found_act = [g for g in ACTIVATION if g.upper() in keep]
    found_hall = [g for g in hallmark if g.upper() in keep]
    if len(found_act) != 5:
        raise SystemExit("activation score incomplete: %s" % found_act)
    print("endpoint members recovered: activation %d/5, hallmark %d/54" % (len(found_act), len(found_hall)), flush=True)

    def score(col: int, genes: list[str]) -> float:
        total = max(totals[col], 1)
        vals = [math.log2(keep.get(g.upper(), {}).get(col, 0) / total * 1e6 + PRIOR) for g in genes]
        return sum(vals) / len(vals)

    wells = []
    for i, name in zip(cols, names):
        if name not in unit_of:
            continue
        fib = totals[i]
        mus = recorded_mouse[name]
        wells.append({
            "library": name, "unit": unit_of[name], "target": target_of(name),
            "fibroblast_total_counts": fib,
            "epithelial_count_fraction": round(mus / (mus + fib), 4) if (mus + fib) else 0.0,
            "activation_score": round(score(i, ACTIVATION), 4),
            "hallmark_tgf_beta_score": round(score(i, found_hall), 4),
            "eligible": "yes" if fib >= FLOOR else "no",
        })
    write_tsv(OUT / "stage3_well_scores.tsv", wells)
    units = sorted({w["unit"] for w in wells})

    diagnostics = []
    for unit in units:
        elig = [w for w in wells if w["unit"] == unit and w["eligible"] == "yes"]
        for endpoint in ["activation_score", "hallmark_tgf_beta_score"]:
            y = [w[endpoint] for w in elig]
            d = [math.log2(w["fibroblast_total_counts"]) for w in elig]
            e = [w["epithelial_count_fraction"] for w in elig]
            diagnostics.append({
                "unit": unit, "endpoint": endpoint, "eligible_wells": len(elig),
                "spearman_endpoint_vs_log_depth": round(spearman(y, d), 4),
                "spearman_endpoint_vs_epithelial_fraction": round(spearman(y, e), 4),
            })
    write_tsv(OUT / "stage3_diagnostics.tsv", diagnostics)

    contrasts = []
    for endpoint in ["activation_score", "hallmark_tgf_beta_score"]:
        for target in AXIS:
            for unit in units:
                here = [w for w in wells if w["unit"] == unit]
                hit = [w for w in here if w["target"] == target]
                if not hit:
                    continue
                w = hit[0]
                ctrl = [c for c in here if c["target"] in CONTROLS and c["eligible"] == "yes"]
                lo, hi = w["fibroblast_total_counts"] / BAND, w["fibroblast_total_counts"] * BAND
                band = [c for c in ctrl if lo <= c["fibroblast_total_counts"] <= hi]
                elig = [c for c in here if c["eligible"] == "yes"]
                adj = None
                if len(elig) >= 5 and w["eligible"] == "yes":
                    y = [c[endpoint] for c in elig]
                    e = [c["epithelial_count_fraction"] for c in elig]
                    res = ols_residuals(y, e)
                    by_lib = {c["library"]: r for c, r in zip(elig, res)}
                    ctrl_res = [by_lib[c["library"]] for c in ctrl if c["library"] in by_lib]
                    if ctrl_res:
                        adj = round(by_lib[w["library"]] - sum(ctrl_res) / len(ctrl_res), 4)
                contrasts.append({
                    "endpoint": endpoint, "target": target, "unit": unit,
                    "library": w["library"],
                    "well_eligible": w["eligible"],
                    "fibroblast_total_counts": w["fibroblast_total_counts"],
                    "score": w[endpoint],
                    "controls_used": len(ctrl),
                    "control_mean": round(sum(c[endpoint] for c in ctrl) / len(ctrl), 4) if ctrl else None,
                    "effect_vs_controls": round(w[endpoint] - sum(c[endpoint] for c in ctrl) / len(ctrl), 4) if ctrl else None,
                    "controls_in_depth_band": len(band),
                    "effect_vs_depth_matched_controls": round(w[endpoint] - sum(c[endpoint] for c in band) / len(band), 4) if band else None,
                    "effect_after_epithelial_adjustment": adj,
                    "interpretable": "yes" if w["eligible"] == "yes" and len(band) >= 2 else "no",
                })
    write_tsv(OUT / "stage3_contrasts.tsv", contrasts)

    summary = {}
    for endpoint in ["activation_score", "hallmark_tgf_beta_score"]:
        for target in AXIS:
            rows = [c for c in contrasts if c["endpoint"] == endpoint and c["target"] == target]
            usable = [c for c in rows if c["interpretable"] == "yes" and c["effect_vs_depth_matched_controls"] is not None]
            raw = [c["effect_vs_controls"] for c in rows if c["effect_vs_controls"] is not None]
            matched = [c["effect_vs_depth_matched_controls"] for c in usable]
            adj = [c["effect_after_epithelial_adjustment"] for c in usable if c["effect_after_epithelial_adjustment"] is not None]
            summary["%s|%s" % (endpoint, target)] = {
                "units_with_a_contrast": len(rows),
                "interpretable_units": len(usable),
                "median_effect_vs_controls": round(statistics.median(raw), 4) if raw else None,
                "median_effect_depth_matched": round(statistics.median(matched), 4) if matched else None,
                "units_lower_than_controls_depth_matched": sum(1 for v in matched if v < 0),
                "median_effect_after_epithelial_adjustment": round(statistics.median(adj), 4) if adj else None,
                "per_unit_depth_matched": matched,
            }

    record = {
        "stage": "3, leg 1",
        "governed_by": {"file": rel(FREEZE), "sha256": sha256(FREEZE), "schema": freeze["schema"]},
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": {rel(WORKBOOK): WORKBOOK_SHA, rel(QC): sha256(QC), rel(TOTALS): sha256(TOTALS)},
        "msigdb_gmt_sha256": sha256(gmt),
        "human_rows_scanned": seen,
        "human_total_mismatches_against_a10": 0,
        "endpoint_members_recovered": {"activation": len(found_act), "hallmark": len(found_hall)},
        "eligibility_floor": FLOOR,
        "depth_band_factor": BAND,
        "control_targets": sorted(CONTROLS),
        "eligible_wells_per_unit": {u: sum(1 for w in wells if w["unit"] == u and w["eligible"] == "yes") for u in units},
        "no_p_value": True,
        "summary": summary,
    }
    (OUT / "stage3_run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print("\n=== primary endpoint, activation score ===")
    for target in AXIS:
        s = summary["activation_score|%s" % target]
        print("  %-6s interpretable %d/%d units, median vs controls %s, depth-matched %s, lower in %s units, after epithelial adj %s"
              % (target, s["interpretable_units"], s["units_with_a_contrast"], s["median_effect_vs_controls"],
                 s["median_effect_depth_matched"], s["units_lower_than_controls_depth_matched"],
                 s["median_effect_after_epithelial_adjustment"]))
    print("\n=== secondary endpoint, Hallmark TGF-beta ===")
    for target in AXIS:
        s = summary["hallmark_tgf_beta_score|%s" % target]
        print("  %-6s median depth-matched %s, lower in %s of %s interpretable units"
              % (target, s["median_effect_depth_matched"], s["units_lower_than_controls_depth_matched"], s["interpretable_units"]))
    print("\n=== depth diagnostic, primary endpoint ===")
    for d in diagnostics:
        if d["endpoint"] == "activation_score":
            print("  %s: rho vs log depth %s, vs epithelial fraction %s (%d eligible wells)"
                  % (d["unit"], d["spearman_endpoint_vs_log_depth"], d["spearman_endpoint_vs_epithelial_fraction"], d["eligible_wells"]))
    print("\nwrote %s" % OUTPUTS)


if __name__ == "__main__":
    main()
