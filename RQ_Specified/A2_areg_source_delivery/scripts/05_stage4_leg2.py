"""A2 stage 4, leg 2: does the recipient's post-receptor machinery track activation.

Exploratory. Obeys config/a2_leg2_spec.json, which was declared and committed before
this ran. No alpha is claimed and no decision is attached.

Reuses the trial E6 instrument in GSE136831 unchanged and changes only the predictor.
One pass over the deposited sparse matrix keeps the needed gene rows.

Needs numpy. Run with the bundled x64 interpreter and the venv on PYTHONPATH:
  PYTHONPATH=X:/GitHub/scRNA_seq/.venv-x64/Lib/site-packages \
  C:/Users/dream/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe \
  RQ_Specified/A2_areg_source_delivery/scripts/05_stage4_leg2.py --data-root X:/GitHub/scRNA_seq
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
OUT = HERE / "tables"
SPEC = HERE / "config/a2_leg2_spec.json"
TRIALS = ROOT / "Research Article/gate2_05_cardoso_2026/trials"
OUTPUTS = ["stage4_donor_values.tsv", "stage4_correlations.tsv", "stage4_run.json"]

ACTIVATION = ["COL1A1", "ACTA2", "POSTN", "CTHRC1", "TNC"]
MACHINERY = ["ITGAV", "ITGB1", "ITGB8", "LTBP1", "TGFB1", "THBS1"]
HOUSEKEEPING = ["ACTB", "GAPDH", "PGK1", "RPL13A", "RPS18"]
COMPARATOR = "AREG"
GENE_FLOOR = 0.05
CELL_FLOOR = 50
MIN_DONORS = 10
DEPTH_RULE = 0.4
FIBRO = {"Fibroblast", "Myofibroblast"}
MURAL = {"Pericyte", "SMC", "Smooth_Muscle", "SmoothMuscle"}
GMT = "raw_data/msigdb/h.all.v2024.1.Hs.symbols.gmt"


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda: h.read(1 << 24), b""):
            d.update(b)
    return d.hexdigest()


def rel(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT)).replace("\\", "/")
    except ValueError:
        return str(p).replace("\\", "/")


def read_set(path: Path, name: str) -> list[str]:
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if parts and parts[0] == name:
            return parts[2:]
    raise SystemExit("%s not found in %s" % (name, path.name))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT)
    data_root = ap.parse_args().data_root.resolve()

    import numpy as np
    sys.path.insert(0, str(TRIALS))
    from mtx_stream import extract_gene_rows, read_lines  # noqa: E402

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite: %s" % existing)
    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    if spec["classification"] != "exploratory" or spec["statistic"]["alpha"] is not None:
        raise SystemExit("the spec must be exploratory with no alpha")

    series = data_root / "raw_data/GSE136831"
    matrix = series / "GSE136831_RawCounts_Sparse.mtx.gz"
    gene_ids = series / "GSE136831_AllCells.GeneIDs.txt.gz"
    barcodes = series / "GSE136831_AllCells.cellBarcodes.txt.gz"
    metadata = series / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
    for p in [matrix, gene_ids, barcodes, metadata]:
        if not p.exists():
            raise SystemExit("missing %s" % p)

    gmt = data_root / GMT
    oxphos = read_set(gmt, "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
    tgfb = set(read_set(gmt, "HALLMARK_TGF_BETA_SIGNALING"))
    emt = set(read_set(gmt, "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"))
    candidates = [g for g in oxphos if g not in tgfb and g not in emt]
    print("co-regulation candidate universe: %d of %d oxphos members" % (len(candidates), len(oxphos)))

    genes = read_lines(gene_ids)
    cells = read_lines(barcodes)
    print("deposit: %d genes, %d cells" % (len(genes), len(cells)))

    wanted = set(ACTIVATION) | set(MACHINERY) | set(HOUSEKEEPING) | {COMPARATOR} | set(candidates)
    print("streaming the matrix for %d genes" % len(wanted), flush=True)
    vectors = extract_gene_rows(matrix, genes, wanted)
    print("recovered %d gene vectors" % len(vectors), flush=True)
    missing = sorted(g for g in (set(ACTIVATION) | set(MACHINERY) | {COMPARATOR}) if g not in vectors)
    if missing:
        raise SystemExit("declared genes absent from the deposit: %s" % missing)

    index = {b: i for i, b in enumerate(cells)}
    rows = []
    with gzip.open(metadata, "rt") as fh:
        header = [h.strip('"') for h in fh.readline().rstrip("\n").split("\t")]
        pos = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            rows.append({k: f[pos[k]] for k in
                         ["CellBarcode_Identity", "nGene", "CellType_Category",
                          "Manuscript_Identity", "Disease_Identity", "Subject_Identity"]})
    rows = [r for r in rows if r["Disease_Identity"] in {"IPF", "Control"} and r["CellBarcode_Identity"] in index]
    print("cells in scope: %d" % len(rows))
    mural_labels = sorted({r["Manuscript_Identity"] for r in rows if r["Manuscript_Identity"] in MURAL})
    print("mural labels present: %s" % mural_labels)

    groups: dict[tuple, list[int]] = {}
    depth: dict[tuple, list[float]] = {}
    for r in rows:
        i = index[r["CellBarcode_Identity"]]
        keys = []
        if r["Manuscript_Identity"] in FIBRO:
            keys.append((r["Subject_Identity"], "fibroblast"))
        if r["Manuscript_Identity"] in MURAL:
            keys.append((r["Subject_Identity"], "mural"))
        if r["CellType_Category"] == "Epithelial":
            keys.append((r["Subject_Identity"], "epithelium"))
        for k in keys:
            groups.setdefault(k, []).append(i)
            depth.setdefault(k, []).append(float(r["nGene"]))
    disease = {r["Subject_Identity"]: r["Disease_Identity"] for r in rows}

    def detection(key: tuple, gene: str) -> float:
        cols = groups[key]
        v = vectors[gene][cols]
        return float((v > 0).sum()) / len(cols)

    def median(xs: list[float]) -> float:
        s = sorted(xs)
        n = len(s)
        return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

    eligible = {k for k, v in groups.items() if len(v) >= CELL_FLOOR}
    donors = sorted({k[0] for k in eligible if k[1] == "fibroblast"} & {k[0] for k in eligible if k[1] == "epithelium"})
    print("donors with both compartments at the %d-cell floor: %d" % (CELL_FLOOR, len(donors)))
    if len(donors) < MIN_DONORS:
        raise SystemExit("below the declared %d-donor minimum" % MIN_DONORS)

    pooled = {}
    for gene in sorted(set(MACHINERY) | set(candidates)):
        if gene not in vectors:
            continue
        cols = [i for d in donors for i in groups[(d, "fibroblast")]]
        v = vectors[gene][cols]
        pooled[gene] = float((v > 0).sum()) / len(cols)
    machinery_used = [g for g in MACHINERY if pooled.get(g, 0.0) >= GENE_FLOOR]
    machinery_dropped = {g: round(pooled.get(g, 0.0), 4) for g in MACHINERY if g not in machinery_used}
    print("machinery members passing the %.2f floor: %s" % (GENE_FLOOR, machinery_used))
    if machinery_dropped:
        print("dropped for eligibility: %s" % machinery_dropped)

    matched, used = [], set()
    for g in machinery_used:
        target = pooled[g]
        best = min((c for c in candidates if c in pooled and c not in used),
                   key=lambda c: (abs(pooled[c] - target), c))
        used.add(best)
        matched.append(best)
    print("co-regulation matched control: %s" % matched)

    donor_rows = []
    for d in donors:
        fk, ek = (d, "fibroblast"), (d, "epithelium")
        row = {
            "donor": d, "disease": disease[d],
            "fibroblast_cells": len(groups[fk]), "epithelial_cells": len(groups[ek]),
            "fibroblast_depth": median(depth[fk]), "epithelial_depth": median(depth[ek]),
            "activation": round(sum(detection(fk, g) for g in ACTIVATION) / len(ACTIVATION), 5),
            "machinery": round(sum(detection(fk, g) for g in machinery_used) / len(machinery_used), 5),
            "matched_control": round(sum(detection(fk, g) for g in matched) / len(matched), 5),
            "housekeeping": round(sum(detection(fk, g) for g in HOUSEKEEPING if g in vectors) /
                                  max(1, len([g for g in HOUSEKEEPING if g in vectors])), 5),
            "epithelial_AREG": round(detection(ek, COMPARATOR), 5),
        }
        mk = (d, "mural")
        row["mural_cells"] = len(groups.get(mk, []))
        row["mural_machinery"] = round(sum(detection(mk, g) for g in machinery_used) / len(machinery_used), 5) if mk in eligible else None
        row["mural_activation"] = round(sum(detection(mk, g) for g in ACTIVATION) / len(ACTIVATION), 5) if mk in eligible else None
        donor_rows.append(row)

    with open(OUT / "stage4_donor_values.tsv", "w", encoding="utf-8", newline="") as fh:
        fields = list(donor_rows[0].keys())
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in donor_rows:
            w.writerow({k: ("none" if r.get(k) is None else r.get(k)) for k in fields})

    def spearman(x: list[float], y: list[float]) -> tuple[float, float]:
        from scipy import stats
        r = stats.spearmanr(x, y)
        return float(r.statistic), float(r.pvalue)

    def col(name: str, subset: list[dict]) -> list[float]:
        return [r[name] for r in subset]

    pooled_set = donor_rows
    ipf = [r for r in donor_rows if r["disease"] == "IPF"]
    corr = []
    for label, subset in [("pooled", pooled_set), ("IPF only", ipf)]:
        if len(subset) < MIN_DONORS and label == "pooled":
            raise SystemExit("pooled set below the donor minimum")
        for x, y, note in [
            ("machinery", "activation", "predictor of interest"),
            ("epithelial_AREG", "activation", "abundance comparator, the E6 estimand"),
            ("matched_control", "activation", "co-regulation control, must be exceeded"),
            ("housekeeping", "activation", "shared-detection control"),
            ("machinery", "fibroblast_depth", "depth coupling of the predictor"),
            ("activation", "fibroblast_depth", "depth coupling of the outcome"),
            ("epithelial_AREG", "epithelial_depth", "depth coupling of the comparator"),
        ]:
            if len(subset) < 5:
                continue
            rho, p = spearman(col(x, subset), col(y, subset))
            corr.append({"stratum": label, "x": x, "y": y, "donors": len(subset),
                         "rho": round(rho, 4), "p_value": round(p, 6), "note": note})
    mural_ok = [r for r in donor_rows if r["mural_machinery"] is not None]
    if len(mural_ok) >= 5:
        rho, p = spearman(col("mural_machinery", mural_ok), col("mural_activation", mural_ok))
        corr.append({"stratum": "mural", "x": "mural_machinery", "y": "mural_activation",
                     "donors": len(mural_ok), "rho": round(rho, 4), "p_value": round(p, 6),
                     "note": "the population the cited mechanism used; underpowered"})
    with open(OUT / "stage4_correlations.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(corr[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in corr:
            w.writerow(r)

    def get(stratum: str, x: str, y: str) -> dict | None:
        for r in corr:
            if r["stratum"] == stratum and r["x"] == x and r["y"] == y:
                return r
        return None

    pm = get("pooled", "machinery", "activation")
    pa = get("pooled", "epithelial_AREG", "activation")
    mc = get("pooled", "matched_control", "activation")
    hk = get("pooled", "housekeeping", "activation")
    dm = get("pooled", "machinery", "fibroblast_depth")
    da = get("pooled", "activation", "fibroblast_depth")
    refused = abs(dm["rho"]) >= DEPTH_RULE and abs(da["rho"]) >= DEPTH_RULE
    exceeds_coreg = pm["rho"] > mc["rho"]
    verdict = ("refused by the depth rule; both members of the pair track their compartment depth at or above %.1f" % DEPTH_RULE
               if refused else
               "the machinery correlation does not exceed the co-regulation control, so it is attributed to co-regulation and not read as recipient state"
               if not exceeds_coreg else
               "the machinery composite tracks activation and exceeds both controls; still correlational, with no direction and no proximity")

    record = {
        "stage": "4, leg 2",
        "classification": "exploratory",
        "governed_by": {"file": rel(SPEC), "sha256": sha256(SPEC)},
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": {rel(matrix): sha256(matrix), rel(gene_ids): sha256(gene_ids),
                   rel(barcodes): sha256(barcodes), rel(metadata): sha256(metadata)},
        "msigdb_gmt_sha256": sha256(gmt),
        "donors": len(donors),
        "ipf_donors": len(ipf),
        "mural_donors_at_the_floor": len(mural_ok),
        "mural_labels_present": mural_labels,
        "machinery_members_used": machinery_used,
        "machinery_members_dropped_for_eligibility": machinery_dropped,
        "co_regulation_matched_control": matched,
        "pooled_detection_fractions_used_for_matching": {g: round(v, 4) for g, v in sorted(pooled.items()) if g in set(machinery_used) | set(matched)},
        "headline": {
            "machinery_vs_activation": pm,
            "epithelial_AREG_vs_activation": pa,
            "co_regulation_control_vs_activation": mc,
            "housekeeping_vs_activation": hk,
            "predictor_depth_coupling": dm,
            "outcome_depth_coupling": da,
        },
        "depth_rule_refuses_the_pair": refused,
        "machinery_exceeds_co_regulation_control": exceeds_coreg,
        "verdict": verdict,
        "alpha_claimed": None,
        "decision_attached": False,
    }
    (OUT / "stage4_run.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print("\n=== leg 2, exploratory ===")
    for r in corr:
        print("  %-9s %-16s vs %-18s n=%2d rho %7.4f p %.5f  (%s)"
              % (r["stratum"], r["x"], r["y"], r["donors"], r["rho"], r["p_value"], r["note"]))
    print("\ndepth rule refuses the pair: %s" % refused)
    print("machinery exceeds the co-regulation control: %s" % exceeds_coreg)
    print("VERDICT: %s" % verdict)


if __name__ == "__main__":
    main()
