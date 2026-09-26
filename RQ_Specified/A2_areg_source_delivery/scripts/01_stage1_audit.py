"""A2 stage 1: instrument, eligibility and power audit. Computes no endpoint.

The plan authorizes this stage only. It establishes from files whether the design
can carry the declared test, and it must not score either endpoint on any well.

Seven items, in the plan's order:

1. plate-3 layout, recomputed from the deposited per-library QC table;
2. knockout validation for the EGFR-axis targets, reusing A10's recorded metric
   and extending it to Erbb2 and Erbb4, which A10 never tested;
3. endpoint gene coverage against the workbook's human gene index;
4. a covariate audit of depth and composition for the plate-3 wells;
5. a precedent check, that no per-target fibroblast endpoint contrast exists;
6. power for the declared rank statistic, computed exactly;
7. the cross-species assumption, attempted from local resources.

Standard library only. Inputs are hash-verified against the records that created
them. Refuses to overwrite its own outputs.
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
import subprocess
import sys
import time
import zipfile
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
A10 = ROOT / "RQ_Specified/A10_organoid_growth_outcome"
CACHE = A10 / "cache"
OUT = HERE / "tables"
REPORTS = HERE / "reports"
CONTRACT = HERE / "config/a2_delivery_contract.json"

WORKBOOK = CACHE / "GSE307112_gene_counts.xlsx"
QC = CACHE / "GSE307112_xenome_stats.csv.gz"
DESIGN = CACHE / "GSE307112_plate_design.csv.gz"
EXPECTED_SHA = {
    WORKBOOK: "60736e5c153a4fb92908f19148ebb5ef2b59f7ade97fb1056deee7d0205acc33",
    QC: "045b9da623db54a59948f1e8e6b68c1587fa4786a6838a51a890dfe59064edf5",
    DESIGN: "8f6950e6286ed438748bfcb12077e683abe74aa525678096f9ad9bea282b8ad1",
}
SHEETS = {"mouse": "xl/worksheets/sheet1.xml", "human": "xl/worksheets/sheet2.xml"}

AXIS = ["AREG", "EGFR", "ERBB2", "ERBB3", "ERBB4", "ITGB6"]
CONTROL = "TDTOMATO"
SECONDARY = "HBEGF"
PRIOR = 1.0
ACTIVATION = ["COL1A1", "ACTA2", "POSTN", "CTHRC1", "TNC"]
PRIMARY_SET = "HALLMARK_TGF_BETA_SIGNALING"
GMT_RELATIVE = "raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt"
COVERAGE_FLOOR = 0.90

OUTPUTS = [
    "stage1_layout.tsv",
    "stage1_knockout_axis.tsv",
    "stage1_endpoint_coverage.tsv",
    "stage1_covariates.tsv",
    "stage1_precedent.json",
    "stage1_power.json",
    "stage1_run.json",
]

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


def col_letter(index: int) -> str:
    out = ""
    n = index + 1
    while n:
        n, r = divmod(n - 1, 26)
        out = chr(65 + r) + out
    return out


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("none" if r.get(k) in (None, "") else r.get(k)) for k in fieldnames})


def target_of(library: str) -> str:
    parts = library.split("_", 2)
    return parts[2].upper() if len(parts) > 2 else ""


# ---------------------------------------------------------------- item 1
def item1_layout(libraries: list[dict]) -> tuple[list[dict], dict]:
    units: dict[tuple[str, str], dict[str, int]] = {}
    for r in libraries:
        key = (r["plate"], r["plate_rep"])
        units.setdefault(key, {})
        t = r["crispr_target"].upper()
        units[key][t] = units[key].get(t, 0) + 1
    rows = []
    for (plate, rep) in sorted(units):
        counts = units[(plate, rep)]
        row = {"plate": plate, "unit": rep, "wells_total": sum(counts.values()),
               "distinct_targets": len(counts), "control_wells": counts.get(CONTROL, 0),
               "hbegf_wells": counts.get(SECONDARY, 0)}
        for t in AXIS:
            row[t.lower() + "_wells"] = counts.get(t, 0)
        rows.append(row)
    plate3 = [r for r in rows if r["plate"] == "plate3"]
    verdict = {
        "plate3_units": len(plate3),
        "plate3_wells_each": sorted({r["wells_total"] for r in plate3}),
        "axis_one_well_per_plate3_unit": all(r[t.lower() + "_wells"] == 1 for r in plate3 for t in AXIS),
        "control_two_wells_per_plate3_unit": all(r["control_wells"] == 2 for r in plate3),
        "axis_absent_outside_plate3": all(r[t.lower() + "_wells"] == 0 for r in rows if r["plate"] != "plate3" for t in AXIS),
        "hbegf_plate": sorted({r["plate"] for r in rows if r["hbegf_wells"]}),
        "hbegf_wells_per_unit": sorted({r["hbegf_wells"] for r in rows if r["hbegf_wells"]}),
        "control_units": sum(1 for r in rows if r["control_wells"]),
    }
    verdict["matches_the_contract"] = bool(
        verdict["plate3_units"] == 4
        and verdict["plate3_wells_each"] == [60]
        and verdict["axis_one_well_per_plate3_unit"]
        and verdict["control_two_wells_per_plate3_unit"]
        and verdict["axis_absent_outside_plate3"]
        and verdict["hbegf_plate"] == ["plate4"]
    )
    return rows, verdict


# ---------------------------------------------------------------- sheet scans
def scan_mouse(zf: zipfile.ZipFile, wanted: set[str]) -> dict:
    """One pass over the mouse sheet: library totals plus the wanted gene rows."""
    keep: dict[str, dict[int, int]] = {}
    totals: dict[int, int] = {}
    order: dict[int, str] = {}
    header = None
    sym_col = None
    seen = 0
    start = time.time()
    wanted_upper = {w.upper() for w in wanted}
    with zf.open(SHEETS["mouse"]) as handle:
        buffer = b""
        tail = False
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
                symbol = None
                values: dict[int, int] = {}
                for ref, num, txt in cells:
                    i = col_index(ref)
                    if i == sym_col:
                        symbol = txt.decode()
                    elif i in totals and num:
                        v = int(num) if b"." not in num else int(float(num))
                        totals[i] += v
                        if v:
                            values[i] = v
                if symbol and symbol.upper() in wanted_upper:
                    keep.setdefault(symbol, {}).update(values)
                if seen % 20000 == 0:
                    print("   mouse: %d rows, %.0fs" % (seen, time.time() - start), flush=True)
            if tail:
                break
    return {"totals": totals, "order": order, "keep": keep, "rows": seen}


def scan_human_symbols(zf: zipfile.ZipFile) -> dict:
    """Symbol index only. No count is read, so no endpoint can be computed here."""
    with zf.open(SHEETS["human"]) as handle:
        head = handle.read(1 << 22)
    m = ROW.search(head)
    if not m:
        raise SystemExit("could not read the human header row")
    header = {col_index(c[0]): (c[2] or c[1]).decode() for c in CELL.findall(m.group(1))}
    sym_col = next(i for i, v in header.items() if v == "symbol")
    letter = col_letter(sym_col)
    libs = sum(1 for v in header.values() if LIBNAME.match(v))
    pattern = re.compile(('<c r="%s(\\d+)"[^>]*?><is><t>([^<]*)</t></is></c>' % letter).encode())
    symbols: list[str] = []
    with zf.open(SHEETS["human"]) as handle:
        buffer = b""
        tail = False
        while True:
            chunk = handle.read(1 << 23)
            if not chunk:
                tail = True
            buffer += chunk
            cut = buffer.rfind(b"</c>")
            if tail:
                block, buffer = buffer, b""
            elif cut >= 0:
                block, buffer = buffer[: cut + 4], buffer[cut + 4 :]
            else:
                block = b""
            for hit in pattern.finditer(block):
                if hit.group(1) != b"1":
                    symbols.append(hit.group(2).decode())
            if tail:
                break
    return {"symbol_column": letter, "libraries_in_header": libs, "symbols": symbols}


# ---------------------------------------------------------------- item 2
def item2_knockout(scan: dict, recorded_totals: dict[str, int], genes: list[str]) -> tuple[list[dict], dict]:
    order = scan["order"]
    cols = sorted(order)
    names = [order[i] for i in cols]
    mismatches = [n for i, n in zip(cols, names) if scan["totals"][i] != recorded_totals.get(n)]
    rows = []
    for gene in genes:
        found = next((g for g in scan["keep"] if g.upper() == gene.upper()), None)
        if found is None:
            rows.append({"gene": gene, "found_in_sheet": "no", "libraries_targeted": 0,
                         "mean_log2cpm_targeted": None, "mean_log2cpm_other": None,
                         "difference": None, "verdict": "absent from the mouse index"})
            continue
        counts = scan["keep"][found]
        cpm = []
        for i, name in zip(cols, names):
            total = max(scan["totals"][i], 1)
            cpm.append(math.log2(counts.get(i, 0) / total * 1e6 + PRIOR))
        hit = [j for j, name in enumerate(names) if target_of(name) == gene.upper()]
        other = [j for j in range(len(names)) if j not in set(hit)]
        if not hit:
            rows.append({"gene": found, "found_in_sheet": "yes", "libraries_targeted": 0,
                         "mean_log2cpm_targeted": None, "mean_log2cpm_other": None,
                         "difference": None, "verdict": "not a target in this screen"})
            continue
        mt = sum(cpm[j] for j in hit) / len(hit)
        mo = sum(cpm[j] for j in other) / len(other)
        rows.append({"gene": found, "found_in_sheet": "yes", "libraries_targeted": len(hit),
                     "mean_log2cpm_targeted": round(mt, 3), "mean_log2cpm_other": round(mo, 3),
                     "difference": round(mt - mo, 3),
                     "verdict": "reduced when targeted" if mt < mo else "NOT reduced when targeted"})
    return rows, {"library_total_mismatches": len(mismatches), "mismatched_libraries": mismatches[:5],
                  "libraries_compared": len(names)}


# ---------------------------------------------------------------- item 6
def item6_power(units: int, wells: int, alpha: float) -> dict:
    """Exact null for the declared statistic: ranks uniform in each unit."""
    # joint distribution of (below-median count, rank sum) by dynamic programming
    half = wells // 2
    dp = {(0, 0): 1}
    for _ in range(units):
        nxt: dict[tuple[int, int], int] = {}
        for (below, total), n in dp.items():
            for r in range(1, wells + 1):
                key = (below + (1 if r <= half else 0), total + r)
                nxt[key] = nxt.get(key, 0) + n
        dp = nxt
    space = wells ** units
    sums: dict[int, int] = {}
    for (below, total), n in dp.items():
        sums[total] = sums.get(total, 0) + n
    cumulative = 0
    tail: dict[int, float] = {}
    for s in sorted(sums):
        cumulative += sums[s]
        tail[s] = cumulative / space
    critical = max((s for s in sorted(sums) if tail[s] <= alpha), default=None)
    joint = sum(n for (below, total), n in dp.items()
                if below >= units - 1 and critical is not None and total <= critical)
    min_p = tail[min(sums)]
    return {
        "units": units,
        "wells_per_unit": wells,
        "rank_space": space,
        "smallest_attainable_one_sided_p": min_p,
        "smallest_attainable_p_scientific": "%.3g" % min_p,
        "alpha": alpha,
        "critical_rank_sum_at_alpha": critical,
        "p_at_critical_sum": round(tail[critical], 5) if critical is not None else None,
        "mean_rank_needed": round(critical / units, 2) if critical is not None else None,
        "mean_percentile_needed": round(critical / units / wells, 4) if critical is not None else None,
        "exact_size_of_the_joint_rule": round(joint / space, 5),
        "joint_rule": "rank sum at or below the critical value and the well below its unit median in at least %d of %d units" % (units - 1, units),
        "reading": "The joint rule is conservative: its exact size is below the nominal alpha, so a rejection is not an inflated test. The design cannot reach a small p by consistency alone; it needs the well near the bottom of its unit in most units.",
    }


# ---------------------------------------------------------------- main
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT,
                    help="checkout that holds raw_data/msigdb; defaults to this repository root")
    args = ap.parse_args()
    data_root = args.data_root.resolve()

    OUT.mkdir(exist_ok=True)
    REPORTS.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite: %s" % existing)

    for path, expected in EXPECTED_SHA.items():
        if not path.exists():
            raise SystemExit("missing input: %s" % path)
        actual = sha256(path)
        if actual != expected:
            raise SystemExit("hash mismatch for %s: %s" % (path.name, actual))
    print("all three cached inputs match their recorded hashes", flush=True)

    with gzip.open(QC, "rt", encoding="utf-8-sig", newline="") as fh:
        libraries = list(csv.DictReader(fh))
    print("libraries in the QC table: %d" % len(libraries), flush=True)

    layout_rows, layout_verdict = item1_layout(libraries)
    write_tsv(OUT / "stage1_layout.tsv", list(layout_rows[0].keys()), layout_rows)
    print("item 1 layout: matches the contract = %s" % layout_verdict["matches_the_contract"], flush=True)
    if not layout_verdict["matches_the_contract"]:
        raise SystemExit("STOP RULE: the plate-3 layout differs from the contract; re-specify before continuing")

    recorded = {}
    with open(A10 / "tables/library_totals.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            recorded[r["library"]] = int(r["mouse_total_counts"])

    gmt = data_root / GMT_RELATIVE
    if not gmt.exists():
        raise SystemExit("missing %s; pass --data-root at the checkout that holds raw_data" % gmt)
    gmt_sha = sha256(gmt)
    primary_genes: list[str] = []
    for line in gmt.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if parts and parts[0] == PRIMARY_SET:
            primary_genes = parts[2:]
    if not primary_genes:
        raise SystemExit("%s not found in %s" % (PRIMARY_SET, gmt.name))

    zf = zipfile.ZipFile(WORKBOOK)
    print("scanning the human sheet for its symbol index", flush=True)
    human = scan_human_symbols(zf)
    index = {s.upper() for s in human["symbols"]}
    print("human symbols: %d unique of %d rows" % (len(index), len(human["symbols"])), flush=True)

    cov_rows = []
    for label, genes, floor in [(PRIMARY_SET, primary_genes, COVERAGE_FLOOR),
                                ("fibroblast_activation_score", ACTIVATION, 1.0)]:
        present = [g for g in genes if g.upper() in index]
        missing = [g for g in genes if g.upper() not in index]
        fraction = len(present) / len(genes)
        cov_rows.append({"endpoint": label, "set_genes": len(genes), "genes_in_human_index": len(present),
                         "fraction": round(fraction, 4), "declared_floor": floor,
                         "clears_floor": "yes" if fraction >= floor else "no",
                         "missing_symbols": ",".join(missing) if missing else "none"})
    write_tsv(OUT / "stage1_endpoint_coverage.tsv", list(cov_rows[0].keys()), cov_rows)
    for r in cov_rows:
        print("item 3 coverage %s: %d/%d = %s, clears floor %s" %
              (r["endpoint"], r["genes_in_human_index"], r["set_genes"], r["fraction"], r["clears_floor"]), flush=True)
    if any(r["clears_floor"] == "no" for r in cov_rows):
        raise SystemExit("STOP RULE: an endpoint fails its declared coverage floor")

    print("scanning the mouse sheet for totals and the two unvalidated receptors", flush=True)
    mouse = scan_mouse(zf, {"Erbb2", "Erbb4", "Areg", "Egfr", "Erbb3", "Itgb6", "Hbegf"})
    ko_rows, ko_check = item2_knockout(mouse, recorded, ["Areg", "Egfr", "Erbb2", "Erbb3", "Erbb4", "Itgb6", "Hbegf"])
    write_tsv(OUT / "stage1_knockout_axis.tsv", list(ko_rows[0].keys()), ko_rows)
    print("item 2 totals cross-check against A10: %d mismatches in %d libraries" %
          (ko_check["library_total_mismatches"], ko_check["libraries_compared"]), flush=True)
    if ko_check["library_total_mismatches"]:
        raise SystemExit("STOP RULE: recomputed mouse library totals disagree with the A10 extraction record")
    # Self-validation: the five genes A10 already tested must reproduce exactly.
    a10_recorded = {}
    with open(A10 / "tables/knockout_validation.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            a10_recorded[r["gene"].upper()] = float(r["difference"])
    reproduced = []
    for row in ko_rows:
        key = row["gene"].upper()
        if key in a10_recorded and row["difference"] is not None:
            agrees = abs(row["difference"] - a10_recorded[key]) < 0.0015
            reproduced.append({"gene": row["gene"], "recomputed": row["difference"],
                               "a10_recorded": a10_recorded[key], "agrees": agrees})
    disagree = [r for r in reproduced if not r["agrees"]]
    print("item 2 reproduction of A10: %d genes checked, %d disagree" % (len(reproduced), len(disagree)), flush=True)
    if disagree:
        raise SystemExit("STOP RULE: this script does not reproduce A10's recorded validation: %s" % disagree)
    ko_check["reproduced_a10_genes"] = reproduced

    areg = next(r for r in ko_rows if r["gene"].upper() == "AREG")
    if areg["difference"] is None or areg["difference"] >= 0:
        raise SystemExit("STOP RULE: the Areg transcript does not fall when targeted; leg 1 is not run")

    order = mouse["order"]
    cols = sorted(order)
    names = [order[i] for i in cols]
    human_totals = {}
    with open(A10 / "tables/library_totals.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            human_totals[r["library"]] = int(r["human_total_counts"])
    by_unit: dict[str, list[tuple[str, int, float]]] = {}
    for r in libraries:
        if r["plate"] != "plate3":
            continue
        name = r["library name"]
        fib = human_totals.get(name)
        mus = recorded.get(name)
        if fib is None or mus is None:
            raise SystemExit("library %s missing from the totals table" % name)
        by_unit.setdefault(r["plate_rep"], []).append((name, fib, mus / (mus + fib) if (mus + fib) else 0.0))
    cov_audit = []
    for unit in sorted(by_unit):
        entries = sorted(by_unit[unit], key=lambda e: e[1])
        n = len(entries)
        fibs = [e[1] for e in entries]
        median = fibs[n // 2]
        for name, fib, epi_fraction in entries:
            t = target_of(name)
            if t not in AXIS + [CONTROL]:
                continue
            rank = sum(1 for f in fibs if f <= fib)
            cov_audit.append({"unit": unit, "target": t, "library": name,
                              "fibroblast_total_counts": fib,
                              "depth_rank_in_unit": rank, "wells_in_unit": n,
                              "depth_percentile": round(rank / n, 3),
                              "unit_median_fibroblast_total": median,
                              "epithelial_count_fraction": round(epi_fraction, 4),
                              "extreme_depth": "yes" if rank <= 3 or rank > n - 3 else "no"})
    write_tsv(OUT / "stage1_covariates.tsv", list(cov_audit[0].keys()), cov_audit)
    extreme = [r for r in cov_audit if r["extreme_depth"] == "yes" and r["target"] in AXIS]
    print("item 4 covariates: %d axis or control wells audited, %d axis wells extreme in depth" %
          (len(cov_audit), len(extreme)), flush=True)

    tracked = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True,
                             encoding="utf-8").stdout.splitlines()
    inspected = []
    for name in tracked:
        p = ROOT / name
        if p.suffix.lower() not in {".tsv", ".csv", ".json", ".md", ".py"} or not p.exists():
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except (UnicodeDecodeError, OSError):
            continue
        if PRIMARY_SET in text or "CTHRC1" in text:
            has_target_column = bool(re.search(r"\btarget\b", text)) or bool(LIBNAME.search(text))
            inspected.append({"file": name, "sha256": sha256(p),
                              "mentions_primary_set": PRIMARY_SET in text,
                              "mentions_activation_gene": "CTHRC1" in text,
                              "has_per_library_or_per_target_rows": has_target_column})
    scores_header = (A10 / "tables/module_scores.tsv").read_text(encoding="utf-8").splitlines()[0].split("\t")
    fib_columns = [c for c in scores_header if c.startswith("fib__")]
    precedent = {
        "question": "has any per-target fibroblast endpoint contrast already been computed in this repository",
        "files_inspected": inspected,
        "a10_module_scores_fibroblast_columns": fib_columns,
        "a10_module_scores_contains_the_primary_endpoint": any(PRIMARY_SET in c for c in fib_columns),
        "verdict": None,
    }
    risky = [f for f in inspected if f["mentions_primary_set"] and f["has_per_library_or_per_target_rows"]
             and f["file"].endswith((".tsv", ".csv"))]
    precedent["files_with_both_the_endpoint_and_per_row_structure"] = [f["file"] for f in risky]
    precedent["verdict"] = ("blind: no tracked table pairs the declared endpoint with per-target rows"
                            if not risky else
                            "NOT blind: re-specify stage 2 in the open, naming what has already been seen")
    (OUT / "stage1_precedent.json").write_text(json.dumps(precedent, indent=2) + "\n", encoding="utf-8")
    print("item 5 precedent: %s" % precedent["verdict"], flush=True)

    power = item6_power(units=4, wells=60, alpha=0.05)
    (OUT / "stage1_power.json").write_text(json.dumps(power, indent=2) + "\n", encoding="utf-8")
    print("item 6 power: smallest attainable p %s, critical rank sum %s, joint size %s" %
          (power["smallest_attainable_p_scientific"], power["critical_rank_sum_at_alpha"],
           power["exact_size_of_the_joint_rule"]), flush=True)

    fasta = sorted(str(p.relative_to(data_root)) for p in (data_root / "raw_data").rglob("*.fa*")) \
        if (data_root / "raw_data").exists() else []
    cross_species = {
        "assumption": "the mouse amphiregulin EGF-like domain activates human EGFR",
        "attempted_from": "local resources under raw_data",
        "local_sequence_files_found": fasta,
        "status": "carried and unverified" if not fasta else "local sequences available; comparison not yet run",
        "how_to_settle_it": "align the mature EGF-like domains of UniProt P31955 (mouse Areg) and P15514 (human AREG); the plan restricts this stage to local resources, so no fetch was made",
    }

    record = {
        "stage": "1",
        "purpose": "instrument, eligibility and power audit; no endpoint scored",
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "contract_sha256": sha256(CONTRACT),
        "inputs": {rel(p): s for p, s in EXPECTED_SHA.items()},
        "msigdb_gmt": {"path": GMT_RELATIVE, "sha256": gmt_sha, "data_root": str(data_root)},
        "a10_tables_read": {
            "tables/library_totals.tsv": sha256(A10 / "tables/library_totals.tsv"),
            "tables/knockout_validation.tsv": sha256(A10 / "tables/knockout_validation.tsv"),
            "tables/module_scores.tsv": sha256(A10 / "tables/module_scores.tsv"),
        },
        "item1_layout": layout_verdict,
        "item2_knockout": ko_check,
        "item3_coverage": cov_rows,
        "item4_covariates": {"axis_or_control_wells_audited": len(cov_audit),
                             "axis_wells_extreme_in_depth": [r["target"] + " in " + r["unit"] for r in extreme]},
        "item5_precedent": {"verdict": precedent["verdict"], "files_inspected": len(inspected)},
        "item6_power": power,
        "item7_cross_species": cross_species,
        "human_sheet": {"symbol_column": human["symbol_column"],
                        "libraries_in_header": human["libraries_in_header"],
                        "symbol_rows": len(human["symbols"]),
                        "unique_symbols": len(index)},
        "mouse_sheet": {"rows": mouse["rows"], "libraries_in_header": len(mouse["order"])},
        "endpoint_scored": False,
        "prohibition_honoured": "no human count was read; the human pass collected symbols only",
    }
    (OUT / "stage1_run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print("\nstage 1 complete; wrote %d outputs" % len(OUTPUTS))


if __name__ == "__main__":
    main()
