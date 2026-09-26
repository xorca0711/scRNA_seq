"""A2 stage 1 addendum: fibroblast-side ligand and receptor covariates.

A review of the stage 1 documents found three gaps that stage 1 could not answer:

* the human fibroblasts are unedited, so their own AREG is untouched by an
  epithelial Areg knockout, and nothing recorded the size of that unremoved
  source;
* the recipient machinery the mechanism names (integrin alphaV, integrin beta
  subunits, the latent complex) was never shown to be expressed in these
  fibroblasts at all, and neither was the EGFR that the cited mechanism makes
  obligatory on the recipient;
* four other EGFR ligands are never perturbed in this screen, so ligand
  redundancy is a rival whose size was unmeasured.

This addendum reads those genes on the human side. It is a covariate pass, not an
endpoint pass, and it refuses to run if any requested gene belongs to either
declared endpoint set. TGFB1 and THBS1 are deliberately excluded for that reason:
both are members of the Hallmark TGF-beta set, which remains a declared endpoint.

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
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
A10 = ROOT / "RQ_Specified/A10_organoid_growth_outcome"
CACHE = A10 / "cache"
OUT = HERE / "tables"
WORKBOOK = CACHE / "GSE307112_gene_counts.xlsx"
QC = CACHE / "GSE307112_xenome_stats.csv.gz"
TOTALS = A10 / "tables/library_totals.tsv"
WORKBOOK_SHA = "60736e5c153a4fb92908f19148ebb5ef2b59f7ade97fb1056deee7d0205acc33"
HUMAN_SHEET = "xl/worksheets/sheet2.xml"
GMT_RELATIVE = "raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt"
PRIMARY_SET = "HALLMARK_TGF_BETA_SIGNALING"
ACTIVATION = ["COL1A1", "ACTA2", "POSTN", "CTHRC1", "TNC"]
PRIOR = 1.0

LIGANDS = ["AREG", "EREG", "TGFA", "HBEGF", "BTC", "EPGN"]
RECEPTORS = ["EGFR", "ERBB2", "ERBB3", "ERBB4"]
MACHINERY = ["ITGAV", "ITGB1", "ITGB8", "LTBP1"]
WANTED = LIGANDS + RECEPTORS + MACHINERY
OUTPUTS = ["stage1b_fibroblast_covariates.tsv", "stage1b_run.json"]

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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT)
    data_root = ap.parse_args().data_root.resolve()

    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite: %s" % existing)
    if sha256(WORKBOOK) != WORKBOOK_SHA:
        raise SystemExit("workbook hash differs from the recorded download")

    gmt = data_root / GMT_RELATIVE
    hallmark: list[str] = []
    for line in gmt.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if parts and parts[0] == PRIMARY_SET:
            hallmark = parts[2:]
    if not hallmark:
        raise SystemExit("could not read %s" % PRIMARY_SET)
    forbidden = sorted(set(WANTED) & (set(hallmark) | set(ACTIVATION)))
    if forbidden:
        raise SystemExit("BLINDING GUARD: %s belong to a declared endpoint set; refusing" % forbidden)
    print("blinding guard passed: none of the %d requested genes is in either endpoint set" % len(WANTED))

    recorded = {}
    with open(TOTALS, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            recorded[r["library"]] = int(r["human_total_counts"])
    with gzip.open(QC, "rt", encoding="utf-8-sig", newline="") as fh:
        plate3 = {r["library name"]: r["plate_rep"] for r in csv.DictReader(fh) if r["plate"] == "plate3"}

    import zipfile
    zf = zipfile.ZipFile(WORKBOOK)
    keep: dict[str, dict[int, int]] = {}
    totals: dict[int, int] = {}
    order: dict[int, str] = {}
    header = None
    sym_col = None
    seen = 0
    wanted_upper = set(WANTED)
    with zf.open(HUMAN_SHEET) as handle:
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
                    keep.setdefault(symbol.upper(), {}).update(values)
                if seen % 20000 == 0:
                    print("   human: %d rows" % seen, flush=True)
            if tail:
                break

    cols = sorted(order)
    names = [order[i] for i in cols]
    mismatch = [n for i, n in zip(cols, names) if totals[i] != recorded.get(n)]
    if mismatch:
        raise SystemExit("STOP: recomputed human totals disagree with the A10 record for %d libraries" % len(mismatch))
    print("human totals cross-check against A10: 0 mismatches in %d libraries" % len(names))

    rows = []
    for gene in WANTED:
        counts = keep.get(gene, {})
        vals = []
        for i, name in zip(cols, names):
            if name not in plate3:
                continue
            total = max(totals[i], 1)
            vals.append((name, counts.get(i, 0), math.log2(counts.get(i, 0) / total * 1e6 + PRIOR)))
        if not vals:
            continue
        detected = sum(1 for _, c, _ in vals if c)
        mean_all = sum(v for _, _, v in vals) / len(vals)
        control = [v for n, _, v in vals if target_of(n) == "TDTOMATO"]
        areg_wells = [v for n, _, v in vals if target_of(n) == "AREG"]
        rows.append({
            "gene": gene,
            "class": "EGFR ligand" if gene in LIGANDS else ("receptor" if gene in RECEPTORS else "recipient machinery"),
            "plate3_wells": len(vals),
            "wells_with_any_count": detected,
            "detection_fraction": round(detected / len(vals), 3),
            "mean_log2cpm_all_plate3": round(mean_all, 3),
            "mean_log2cpm_control_wells": round(sum(control) / len(control), 3) if control else None,
            "mean_log2cpm_areg_ko_wells": round(sum(areg_wells) / len(areg_wells), 3) if areg_wells else None,
            "perturbed_in_this_screen": "yes" if gene in {"AREG", "EGFR", "ERBB2", "ERBB3", "ERBB4", "HBEGF", "ITGB6"} else "no",
        })
    with open(OUT / OUTPUTS[0], "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("none" if r[k] in (None, "") else r[k]) for k in r})

    areg = next(r for r in rows if r["gene"] == "AREG")
    egfr = next(r for r in rows if r["gene"] == "EGFR")
    record = {
        "stage": "1b, addendum",
        "purpose": "fibroblast-side ligand, receptor and recipient-machinery covariates; no endpoint scored",
        "why": "a review found the unremoved human AREG source, the recipient EGFR and machinery, and ligand redundancy all unmeasured",
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": {rel(WORKBOOK): WORKBOOK_SHA, rel(QC): sha256(QC), rel(TOTALS): sha256(TOTALS)},
        "msigdb_gmt_sha256": sha256(gmt),
        "blinding_guard": "none of the requested genes belongs to the Hallmark TGF-beta set or the activation score; TGFB1 and THBS1 were excluded for exactly that reason",
        "genes": WANTED,
        "excluded_for_blinding": ["TGFB1", "THBS1"],
        "human_total_mismatches_against_a10": 0,
        "human_rows_scanned": seen,
        "headline": {
            "human_AREG_detection_fraction_in_plate3_wells": areg["detection_fraction"],
            "human_AREG_mean_log2cpm": areg["mean_log2cpm_all_plate3"],
            "human_EGFR_detection_fraction": egfr["detection_fraction"],
            "human_EGFR_mean_log2cpm": egfr["mean_log2cpm_all_plate3"],
        },
        "endpoint_scored": False,
    }
    (OUT / OUTPUTS[1]).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print("\nwrote %s" % OUTPUTS)
    for r in rows:
        print("  %-7s %-20s detected in %5.1f%% of wells, mean log2 CPM %7.3f" %
              (r["gene"], r["class"], 100 * r["detection_fraction"], r["mean_log2cpm_all_plate3"]))


if __name__ == "__main__":
    main()
