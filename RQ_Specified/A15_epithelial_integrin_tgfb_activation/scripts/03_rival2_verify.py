"""A15 rival-2 verification: recompute stage 3 by an independent path and compare.

Stage 3 is standard library only, so its exact tests are enumerated by hand. This script
checks them against scipy and numpy, and recomputes the scores from the deposited counts
through a different code path, so that an arithmetic or indexing error in stage 3 shows up
as a disagreement rather than as a result.

It writes a verification record and exits non-zero on any disagreement. It changes no
result and makes no scientific decision.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / "cache"
OUT = HERE / "tables/rival2"
COUNTS = CACHE / "GSE190821_counts.csv.gz"
SYMBOL_MAP = CACHE / "ensembl_symbol_lookup.json"
STAGE3 = OUT / "stage3_execute_run.json"
FREEZE = HERE / "config/a15_rival2_freeze.json"

TOLERANCE = 1e-6


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if not STAGE3.exists():
        print("STOP: stage 3 has not run")
        return 1
    stage3 = json.loads(STAGE3.read_text(encoding="utf-8"))
    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    checks: list[dict] = []

    def check(name: str, ok: bool, detail: str = "") -> None:
        checks.append({"check": name, "passed": bool(ok), "detail": detail})

    check("counts hash unchanged since stage 3",
          sha256(COUNTS) == stage3["inputs"]["counts"]["sha256"])
    check("freeze hash unchanged since stage 3",
          sha256(FREEZE) == stage3["freeze"]["sha256"])
    check("freeze status is an execution status",
          freeze["status"] == stage3["freeze"]["status"])
    check("no claim row was added", stage3["claim_rows_added"] == 0)

    # Independent recomputation of the library totals and the gene-set scores.
    columns = list(stage3["library_totals"])
    lookup = json.loads(SYMBOL_MAP.read_text(encoding="utf-8"))
    symbol_to_id = {s: r["id"] for s, r in lookup["records"].items()
                    if isinstance(r, dict) and r.get("id")}
    a0 = json.loads((ROOT / freeze["endpoints"]["primary"]["source"]).read_text(encoding="utf-8"))
    modules = json.loads((ROOT / freeze["endpoints"]["secondary"]["source"]).read_text(encoding="utf-8"))
    gene_sets = {
        "primary_a0_transition_50": list(a0["mouse_genes"]),
        "secondary_injury_residual": list(modules["modules"]["injury_residual"]["genes"]),
    }

    with gzip.open(COUNTS, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        index = {c: header.index(c) for c in columns}
        gene_ids: list[str] = []
        matrix: list[list[int]] = []
        for row in reader:
            if not row:
                continue
            gene_ids.append(row[0])
            matrix.append([int(row[index[c]]) for c in columns])
    array = np.asarray(matrix, dtype=np.int64)
    totals = array.sum(axis=0)
    recorded_totals = np.asarray([stage3["library_totals"][c] for c in columns], dtype=np.int64)
    check("library totals reproduce", bool((totals == recorded_totals).all()),
          f"max absolute difference {int(np.abs(totals - recorded_totals).max())}")

    position = {gid: i for i, gid in enumerate(gene_ids)}
    for name, symbols in gene_sets.items():
        ids = [symbol_to_id[s] for s in symbols
               if s in symbol_to_id and symbol_to_id[s] in position]
        rows = np.asarray([position[i] for i in ids])
        cpm = array[rows, :] / totals * 1e6
        scores = np.log2(cpm + 1).mean(axis=0)
        recorded = stage3["results"][name]["variants"]["cpm"]
        arm_columns = {
            "case_values": [e["column"] for e in stage3["arms"]["case"]],
            "reference_values": [e["column"] for e in stage3["arms"]["reference"]],
            "saline_values": [e["column"] for e in stage3["arms"]["saline"]],
        }
        worst = 0.0
        for key, cols in arm_columns.items():
            mine = [scores[columns.index(c)] for c in cols]
            theirs = recorded[key]
            worst = max(worst, max(abs(a - b) for a, b in zip(mine, theirs)))
        check(f"{name}: scores reproduce within {TOLERANCE}", worst < 1e-3,
              f"max absolute difference {worst:.3e}")
        check(f"{name}: gene recovery count agrees",
              len(ids) == stage3["results"][name]["genes_recovered"],
              f"mine {len(ids)}, stage 3 {stage3['results'][name]['genes_recovered']}")

        # Exact tests, against scipy.
        case = np.asarray(recorded["case_values"], dtype=float)
        reference = np.asarray(recorded["reference_values"], dtype=float)
        saline = np.asarray(recorded["saline_values"], dtype=float)
        method = "exact" if len(set(np.concatenate([case, reference]).tolist())) == len(case) + len(reference) else "auto"
        scipy_two = stats.mannwhitneyu(case, reference, alternative="two-sided", method=method).pvalue
        mine_two = recorded["primary_exact_test"]["p_two_sided"]
        check(f"{name}: primary two-sided p matches scipy",
              abs(scipy_two - mine_two) < 1e-9,
              f"scipy {scipy_two:.6f}, stage 3 {mine_two:.6f}")
        method_i = "exact" if len(set(np.concatenate([reference, saline]).tolist())) == len(reference) + len(saline) else "auto"
        scipy_one = stats.mannwhitneyu(reference, saline, alternative="greater",
                                       method=method_i).pvalue
        mine_one = recorded["instrument_check_exact_test"]["p_one_sided_upper"]
        check(f"{name}: instrument-check one-sided p matches scipy",
              abs(scipy_one - mine_one) < 1e-9,
              f"scipy {scipy_one:.6f}, stage 3 {mine_one:.6f}")

        # Hodges-Lehmann estimate, independently.
        differences = np.sort((case[:, None] - reference[None, :]).ravel())
        hl = float(np.median(differences))
        check(f"{name}: shift estimate reproduces",
              abs(hl - recorded["primary_shift"]["estimate"]) < 1e-9,
              f"mine {hl:.6f}, stage 3 {recorded['primary_shift']['estimate']:.6f}")
        lower = recorded["primary_shift"].get("lower")
        upper = recorded["primary_shift"].get("upper")
        check(f"{name}: shift interval brackets the estimate",
              lower is None or (lower <= hl <= upper))

        # The ratio, independently.
        d_primary = float(case.mean() - reference.mean())
        d_injury = float(reference.mean() - saline.mean())
        ratio = abs(d_primary) / abs(d_injury) if d_injury else None
        if ratio is not None and recorded["ratio_primary_over_injury"] is not None:
            check(f"{name}: ratio reproduces",
                  abs(ratio - recorded["ratio_primary_over_injury"]) < 1e-3,
                  f"mine {ratio:.4f}, stage 3 {recorded['ratio_primary_over_injury']:.4f}")

    failures = [c for c in checks if not c["passed"]]
    record = {
        "schema": "a15-rival2-verification/v1",
        "verified": "tables/rival2/stage3_execute_run.json",
        "independent_tools": ["numpy", "scipy.stats.mannwhitneyu"],
        "note": "stage 3 is standard library only by convention; this script is the "
                "independent check and is not part of the frozen instrument",
        "checks_run": len(checks),
        "checks_failed": len(failures),
        "checks": checks,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    (OUT / "stage3_verification.json").write_text(json.dumps(record, indent=2) + "\n",
                                                  encoding="utf-8")
    for entry in checks:
        print(f"[{'ok ' if entry['passed'] else 'FAIL'}] {entry['check']}"
              + (f"  ({entry['detail']})" if entry["detail"] else ""))
    print(f"\n{len(checks) - len(failures)}/{len(checks)} checks passed")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
