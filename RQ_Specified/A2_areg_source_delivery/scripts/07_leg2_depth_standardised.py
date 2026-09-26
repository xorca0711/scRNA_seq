"""A2 leg 2, second pass: the same pair on a common molecule budget.

Obeys config/a2_leg2_depth_spec.json, which was committed before any standardised value
existed. Exploratory, and explicitly not an independent test: the unstandardised result is
already logged, so the only question here is whether the measure can be made
depth-independent enough for the C51 rule to let the pair be read.

The measure is the C37 treatment. For a cell with T total molecules and k counts of a
gene, the expected detection at budget N is one minus the hypergeometric probability of
drawing none of the k:

    1 - C(T-k, N) / C(T, N),  and 1 when T-k < N

computed through log gamma functions. Cells below the budget are excluded.

Needs numpy and scipy. Run with the bundled x64 interpreter and the venv on PYTHONPATH.
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
CACHE = HERE / "cache"
SPEC = HERE / "config/a2_leg2_depth_spec.json"
TRIALS = ROOT / "Research Article/gate2_05_cardoso_2026/trials"
OUTPUTS = ["leg2_depth_donor_values.tsv", "leg2_depth_correlations.tsv", "leg2_depth_run.json"]

CELL_FLOOR = 50
MIN_DONORS = 10
DEPTH_RULE = 0.4
FIBRO = {"Fibroblast", "Myofibroblast"}
MURAL = {"Pericyte", "SMC"}


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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=ROOT)
    data_root = ap.parse_args().data_root.resolve()

    import numpy as np
    from scipy import stats
    from scipy.special import gammaln
    sys.path.insert(0, str(TRIALS))
    from mtx_stream import extract_gene_rows, read_lines  # noqa: E402

    OUT.mkdir(exist_ok=True)
    CACHE.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite: %s" % existing)

    spec = json.loads(SPEC.read_text(encoding="utf-8"))
    if spec["endpoint_scored_at_declaration"] is not False:
        raise SystemExit("the spec must assert that nothing was scored at declaration")
    frozen = spec["frozen_from_the_first_pass_and_not_re_derived"]
    ACTIVATION = frozen["outcome"]
    MACHINERY = frozen["predictor"]
    HOUSEKEEPING = frozen["shared_detection_control"]
    MATCHED = frozen["co_regulation_control"]
    COMPARATOR = "AREG"
    PRIMARY = spec["budgets"]["primary"]
    BUDGETS = [PRIMARY] + list(spec["budgets"]["sensitivities"])

    series = data_root / "raw_data/GSE136831"
    matrix = series / "GSE136831_RawCounts_Sparse.mtx.gz"
    gene_ids = series / "GSE136831_AllCells.GeneIDs.txt.gz"
    barcodes = series / "GSE136831_AllCells.cellBarcodes.txt.gz"
    metadata = series / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"

    genes = read_lines(gene_ids, column=1)[1:]
    cells = read_lines(barcodes)
    wanted = sorted(set(ACTIVATION) | set(MACHINERY) | set(HOUSEKEEPING) | set(MATCHED) | {COMPARATOR})
    npz = CACHE / "leg2_gene_counts.npz"
    if npz.exists():
        print("reusing the cached gene vectors")
        z = np.load(npz, allow_pickle=False)
        vectors = {g: z[g] for g in wanted}
    else:
        print("streaming the matrix for %d genes" % len(wanted), flush=True)
        vectors = extract_gene_rows(matrix, genes, set(wanted))
        np.savez_compressed(npz, **{g: vectors[g] for g in wanted})
        print("cached %d gene vectors" % len(vectors))
    missing = [g for g in wanted if g not in vectors]
    if missing:
        raise SystemExit("genes absent from the deposit: %s" % missing)

    index = {b: i for i, b in enumerate(cells)}
    rows = []
    with gzip.open(metadata, "rt") as fh:
        header = [h.strip('"') for h in fh.readline().rstrip("\n").split("\t")]
        pos = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            rows.append({k: f[pos[k]] for k in
                         ["CellBarcode_Identity", "nUMI", "nGene", "CellType_Category",
                          "Manuscript_Identity", "Disease_Identity", "Subject_Identity"]})
    rows = [r for r in rows if r["Disease_Identity"] in {"IPF", "Control"} and r["CellBarcode_Identity"] in index]

    groups: dict[tuple, list[int]] = {}
    umi: dict[tuple, list[int]] = {}
    ngene: dict[tuple, list[float]] = {}
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
            umi.setdefault(k, []).append(int(r["nUMI"]))
            ngene.setdefault(k, []).append(float(r["nGene"]))
    disease = {r["Subject_Identity"]: r["Disease_Identity"] for r in rows}
    eligible = {k for k, v in groups.items() if len(v) >= CELL_FLOOR}
    donors = sorted({k[0] for k in eligible if k[1] == "fibroblast"} & {k[0] for k in eligible if k[1] == "epithelium"})
    print("donors with both compartments at the %d-cell floor: %d" % (CELL_FLOOR, len(donors)))
    if len(donors) < MIN_DONORS:
        raise SystemExit("below the declared donor minimum")

    def expected_detection(key: tuple, gene: str, budget: int) -> tuple[float, int, int]:
        """Mean hypergeometric detection expectation over the group's retained cells."""
        cols = np.asarray(groups[key])
        totals = np.asarray(umi[key], dtype=np.float64)
        keep = totals >= budget
        if not keep.any():
            return float("nan"), 0, int(len(cols))
        t = totals[keep]
        k = np.asarray(vectors[gene][cols[keep]], dtype=np.float64)
        rest = t - k
        p_none = np.zeros_like(t)
        ok = rest >= budget
        if ok.any():
            p_none[ok] = np.exp(
                gammaln(rest[ok] + 1) - gammaln(rest[ok] - budget + 1)
                + gammaln(t[ok] - budget + 1) - gammaln(t[ok] + 1)
            )
        return float(np.mean(1.0 - p_none)), int(keep.sum()), int(len(cols))

    def composite(key: tuple, gene_set: list[str], budget: int) -> float:
        vals = [expected_detection(key, g, budget)[0] for g in gene_set]
        return float(np.mean(vals))

    def median(xs) -> float:
        s = sorted(xs)
        n = len(s)
        return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2

    retention = {}
    all_rows: dict[int, list[dict]] = {}
    for budget in BUDGETS:
        kept = sum(expected_detection((d, "fibroblast"), ACTIVATION[0], budget)[1] for d in donors)
        total = sum(len(groups[(d, "fibroblast")]) for d in donors)
        retention[budget] = {"fibroblast_cells_retained": kept, "fibroblast_cells_total": total,
                             "fraction": round(kept / total, 4)}
        donor_rows = []
        for d in donors:
            fk, ek = (d, "fibroblast"), (d, "epithelium")
            mk = (d, "mural")
            row = {
                "donor": d, "disease": disease[d], "budget": budget,
                "fibroblast_cells": len(groups[fk]), "epithelial_cells": len(groups[ek]),
                "fibroblast_depth": median(ngene[fk]), "epithelial_depth": median(ngene[ek]),
                "activation": round(composite(fk, ACTIVATION, budget), 6),
                "machinery": round(composite(fk, MACHINERY, budget), 6),
                "matched_control": round(composite(fk, MATCHED, budget), 6),
                "housekeeping": round(composite(fk, HOUSEKEEPING, budget), 6),
                "epithelial_AREG": round(expected_detection(ek, COMPARATOR, budget)[0], 6),
                "mural_cells": len(groups.get(mk, [])),
            }
            row["mural_machinery"] = round(composite(mk, MACHINERY, budget), 6) if mk in eligible else None
            row["mural_activation"] = round(composite(mk, ACTIVATION, budget), 6) if mk in eligible else None
            donor_rows.append(row)
        all_rows[budget] = donor_rows
        print("  budget %5d: fibroblast cells retained %d of %d" % (budget, kept, total), flush=True)

    flat = [r for b in BUDGETS for r in all_rows[b]]
    with open(OUT / OUTPUTS[0], "w", encoding="utf-8", newline="") as fh:
        fields = list(flat[0].keys())
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in flat:
            w.writerow({k: ("none" if r.get(k) is None else r.get(k)) for k in fields})

    def sp(x, y):
        r = stats.spearmanr(x, y)
        return float(r.statistic), float(r.pvalue)

    corr = []
    for budget in BUDGETS:
        rows_b = all_rows[budget]
        ipf = [r for r in rows_b if r["disease"] == "IPF"]
        for label, subset in [("pooled", rows_b), ("IPF only", ipf)]:
            for x, y, note in [
                ("machinery", "activation", "predictor of interest"),
                ("epithelial_AREG", "activation", "abundance comparator"),
                ("matched_control", "activation", "co-regulation control, must be exceeded"),
                ("housekeeping", "activation", "shared-detection control"),
                ("machinery", "fibroblast_depth", "depth coupling of the predictor"),
                ("activation", "fibroblast_depth", "depth coupling of the outcome"),
                ("epithelial_AREG", "epithelial_depth", "depth coupling of the comparator"),
            ]:
                rho, p = sp([r[x] for r in subset], [r[y] for r in subset])
                corr.append({"budget": budget, "stratum": label, "x": x, "y": y,
                             "donors": len(subset), "rho": round(rho, 4),
                             "p_value": round(p, 6), "note": note})
        mural = [r for r in rows_b if r["mural_machinery"] is not None]
        if len(mural) >= 5:
            rho, p = sp([r["mural_machinery"] for r in mural], [r["mural_activation"] for r in mural])
            corr.append({"budget": budget, "stratum": "mural", "x": "mural_machinery",
                         "y": "mural_activation", "donors": len(mural), "rho": round(rho, 4),
                         "p_value": round(p, 6), "note": "the population the cited mechanism used; underpowered"})
    with open(OUT / OUTPUTS[1], "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(corr[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in corr:
            w.writerow(r)

    def get(budget, stratum, x, y):
        return next(r for r in corr if r["budget"] == budget and r["stratum"] == stratum and r["x"] == x and r["y"] == y)

    pm = get(PRIMARY, "pooled", "machinery", "activation")
    dm = get(PRIMARY, "pooled", "machinery", "fibroblast_depth")
    da = get(PRIMARY, "pooled", "activation", "fibroblast_depth")
    mc = get(PRIMARY, "pooled", "matched_control", "activation")
    hk = get(PRIMARY, "pooled", "housekeeping", "activation")
    pa = get(PRIMARY, "pooled", "epithelial_AREG", "activation")
    still_refused = abs(dm["rho"]) >= DEPTH_RULE and abs(da["rho"]) >= DEPTH_RULE
    exceeds = pm["rho"] > mc["rho"] and pm["rho"] > hk["rho"]
    if still_refused:
        verdict = ("still refused at the primary budget: the predictor tracks depth at %.4f and the outcome at %.4f, "
                   "so this cohort is closed to this question at this measure" % (dm["rho"], da["rho"]))
    elif not exceeds:
        verdict = ("readable, but the correlation does not exceed both controls, so it is attributed to "
                   "co-regulation or shared detection and is not read as recipient state")
    else:
        verdict = ("readable at the primary budget and exceeding both controls; exploratory and correlational, "
                   "with no direction and no proximity, and not promoted to a test because the unstandardised "
                   "value was seen first")

    record = {
        "stage": "leg 2, second pass, depth standardised",
        "classification": "exploratory, not an independent test",
        "governed_by": {"file": rel(SPEC), "sha256": sha256(SPEC)},
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": rel(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "inputs": {rel(gene_ids): sha256(gene_ids), rel(barcodes): sha256(barcodes),
                   rel(metadata): sha256(metadata)},
        "matrix_note": "the 2.1 GB matrix is hashed in the first pass record, tables/stage4_run.json",
        "budgets": BUDGETS,
        "primary_budget": PRIMARY,
        "cell_retention": retention,
        "donors": len(donors),
        "frozen_sets_reused": frozen,
        "primary_results": {"machinery_vs_activation": pm, "predictor_depth": dm, "outcome_depth": da,
                            "co_regulation_control": mc, "shared_detection_control": hk,
                            "abundance_comparator": pa},
        "sensitivity_budgets": {str(b): {"machinery_vs_activation": get(b, "pooled", "machinery", "activation"),
                                         "predictor_depth": get(b, "pooled", "machinery", "fibroblast_depth"),
                                         "outcome_depth": get(b, "pooled", "activation", "fibroblast_depth")}
                                for b in BUDGETS if b != PRIMARY},
        "still_refused_at_primary": still_refused,
        "exceeds_both_controls": exceeds,
        "verdict": verdict,
        "first_pass_for_comparison": {"measure": "raw detection fraction",
                                      "machinery_vs_activation_rho": 0.4331,
                                      "predictor_depth_rho": 0.7696,
                                      "outcome_depth_rho": 0.4116},
        "alpha_claimed": None,
    }
    (OUT / OUTPUTS[2]).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print("\n=== primary budget %d ===" % PRIMARY)
    for r in corr:
        if r["budget"] == PRIMARY:
            print("  %-9s %-16s vs %-18s n=%2d rho %7.4f p %.5f  (%s)"
                  % (r["stratum"], r["x"], r["y"], r["donors"], r["rho"], r["p_value"], r["note"]))
    print("\n=== depth coupling across budgets, pooled ===")
    for b in BUDGETS:
        print("  budget %5d: predictor %7.4f  outcome %7.4f  machinery-vs-activation %7.4f"
              % (b, get(b, "pooled", "machinery", "fibroblast_depth")["rho"],
                 get(b, "pooled", "activation", "fibroblast_depth")["rho"],
                 get(b, "pooled", "machinery", "activation")["rho"]))
    print("\nstill refused at the primary budget: %s" % still_refused)
    print("VERDICT: %s" % verdict)


if __name__ == "__main__":
    main()
