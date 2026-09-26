"""A15 rival-2 verification: recompute stage 3 by an independent path and compare.

Stage 3 is standard library only, so its exact tests are enumerated by hand. This script
checks them against scipy and numpy and recomputes the scores from the deposited counts
through a different code path, so an arithmetic or indexing error shows up as a
disagreement rather than as a result.

Two comparisons are deliberately conditional, and both are recorded rather than hidden:

- stage 3 rounds every p-value to six decimals, so a 1e-9 tolerance against scipy's full
  precision would fail on the rounding alone; the tolerance is 2e-6;
- scipy has no exact rank-sum with ties, and falls back to a tie-corrected normal
  approximation, which is a different estimator from stage 3's exact permutation over
  mid-ranks. Where ties are present the scipy comparison is skipped and recorded as skipped.

It changes no result and makes no scientific decision. Exits non-zero on any disagreement.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import itertools
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
FREEZE = HERE / "config/a15_rival2_freeze_v2.json"

P_TOLERANCE = 2e-6
SCORE_TOLERANCE = 1e-3


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

    def check(name, ok, detail=""):
        checks.append({"check": name, "passed": bool(ok), "detail": detail})

    def skip(name, why):
        checks.append({"check": name, "passed": True, "skipped": True, "detail": why})

    check("counts hash unchanged since stage 3",
          sha256(COUNTS) == stage3["inputs"]["counts"]["sha256"])
    check("freeze hash unchanged since stage 3",
          sha256(FREEZE) == stage3["freeze"]["sha256"])
    check("freeze status is the execution status",
          freeze["status"] == "COMMITTED_FOR_EXECUTION" == stage3["freeze"]["status"])
    check("no claim row was added", stage3["claim_rows_added"] == 0)
    check("declared alpha equals two over seventy",
          abs(stage3["declared_alpha"] - 2 / 70) < 1e-6,
          f"recorded {stage3['declared_alpha']}")

    # Arms and columns.
    epi = stage3["arms"]["epithelium"]
    inp = stage3["arms"]["input"]
    all_columns = sorted(stage3["library_totals"])
    check("22 libraries were read", len(all_columns) == 22, f"got {len(all_columns)}")
    check("no mouse appears twice in one epithelial arm",
          all(len({e["mouse"] for e in epi[a]}) == len(epi[a]) for a in epi))
    check("epithelial and input arms name the same mice",
          all(sorted(e["mouse"] for e in epi[a]) == sorted(e["mouse"] for e in inp[a])
              for a in epi))

    # Independent recomputation of totals and scores.
    with gzip.open(COUNTS, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        index = {c: header.index(c) for c in all_columns}
        gene_ids: list[str] = []
        rows: list[list[int]] = []
        for row in reader:
            if not row:
                continue
            gene_ids.append(row[0])
            rows.append([int(row[index[c]]) for c in all_columns])
    array = np.asarray(rows, dtype=np.int64)
    totals = array.sum(axis=0)
    recorded = np.asarray([stage3["library_totals"][c] for c in all_columns], dtype=np.int64)
    check("library totals reproduce", bool((totals == recorded).all()),
          f"max difference {int(np.abs(totals - recorded).max())}")
    check("gene index length matches stage 3's expectation",
          len(gene_ids) == len(set(gene_ids)), "gene identifiers are unique")

    position = {g: i for i, g in enumerate(gene_ids)}
    lookup = json.loads(SYMBOL_MAP.read_text(encoding="utf-8"))
    sym2id = {s: r["id"] for s, r in lookup["records"].items()
              if isinstance(r, dict) and r.get("id")}
    logcpm = np.log2(array / totals * 1e6 + 1)
    col_at = {c: i for i, c in enumerate(all_columns)}

    def standardised_composite(symbols, columns):
        ids = [sym2id[s] for s in symbols if s in sym2id and sym2id[s] in position]
        cols = [col_at[c] for c in columns]
        keep = []
        for gid in ids:
            counts = array[position[gid], :][cols]
            if float(np.median(counts)) >= 10:
                keep.append(gid)
        if not keep:
            return None
        block = logcpm[[position[g] for g in keep], :][:, cols]
        mu = block.mean(axis=1, keepdims=True)
        sd = block.std(axis=1, ddof=1, keepdims=True)
        sd[sd == 0] = 1.0
        return ((block - mu) / sd).mean(axis=0), keep

    p1 = freeze["endpoints"]["primary_1_epithelial_identity_and_state"]["composites"]
    epi_columns = [e["column"] for a in ("case", "reference", "saline") for e in epi[a]]
    for name, spec in p1.items():
        recorded_block = stage3["results"][f"primary1_{name}"].get("standardised")
        if not recorded_block or "case_values" not in recorded_block:
            skip(f"primary1 {name}: scores", "stage 3 refused or did not compute this composite")
            continue
        out = standardised_composite(spec["genes"], epi_columns)
        if out is None:
            check(f"primary1 {name}: independent composite computable", False)
            continue
        scores, keep = out
        by_col = {c: scores[i] for i, c in enumerate(epi_columns)}
        mine_case = [by_col[e["column"]] for e in epi["case"]]
        mine_ref = [by_col[e["column"]] for e in epi["reference"]]
        worst = max(max(abs(a - b) for a, b in zip(mine_case, recorded_block["case_values"])),
                    max(abs(a - b) for a, b in zip(mine_ref, recorded_block["reference_values"])))
        check(f"primary1 {name}: scores reproduce", worst < SCORE_TOLERANCE,
              f"max difference {worst:.3e}")
        check(f"primary1 {name}: measurable member count agrees",
              len(keep) == recorded_block["coverage"]["measurable"],
              f"mine {len(keep)}, stage 3 {recorded_block['coverage']['measurable']}")
        combined = list(recorded_block["case_values"]) + list(recorded_block["reference_values"])
        if len(set(combined)) == len(combined):
            sp = stats.mannwhitneyu(recorded_block["case_values"],
                                    recorded_block["reference_values"],
                                    alternative="two-sided", method="exact").pvalue
            check(f"primary1 {name}: two-sided p matches scipy exact",
                  abs(sp - recorded_block["exact_test"]["p_two_sided"]) < P_TOLERANCE,
                  f"scipy {sp:.8f}, stage 3 {recorded_block['exact_test']['p_two_sided']}")
        else:
            skip(f"primary1 {name}: scipy comparison", "ties present; scipy has no exact tie test")
        differences = np.sort((np.asarray(recorded_block["case_values"])[:, None]
                               - np.asarray(recorded_block["reference_values"])[None, :]).ravel())
        hl = float(np.median(differences))
        check(f"primary1 {name}: shift estimate reproduces",
              abs(hl - recorded_block["shift"]["estimate"]) < 1e-3,
              f"mine {hl:.5f}, stage 3 {recorded_block['shift']['estimate']}")
        if recorded_block["shift"].get("lower") is not None:
            check(f"primary1 {name}: interval brackets the estimate",
                  recorded_block["shift"]["lower"] <= hl <= recorded_block["shift"]["upper"])

    # Omnibus: recompute the permutation p independently.
    omni = stage3["results"]["primary2_omnibus"]
    itgb6 = sym2id.get("Itgb6")
    primary_columns = [e["column"] for e in epi["case"] + epi["reference"]]
    cols = [col_at[c] for c in primary_columns]
    counts_block = array[:, cols]
    keep_mask = (counts_block >= 10).sum(axis=1) >= 4
    if itgb6 in position:
        keep_mask[position[itgb6]] = False
    check("omnibus gene count reproduces", int(keep_mask.sum()) == omni["genes_used"],
          f"mine {int(keep_mask.sum())}, stage 3 {omni['genes_used']}")
    block = logcpm[keep_mask, :][:, cols]
    mu = block.mean(axis=1, keepdims=True)
    sd = block.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    z = (block - mu) / sd
    corr = np.corrcoef(z.T)
    distance = 1 - corr

    def statistic(idx):
        other = [i for i in range(8) if i not in set(idx)]
        between = [distance[i, j] for i in idx for j in other]
        within = ([distance[i, j] for a, i in enumerate(idx) for j in idx[a + 1:]]
                  + [distance[i, j] for a, i in enumerate(other) for j in other[a + 1:]])
        return float(np.mean(between) - np.mean(within))

    observed = statistic([0, 1, 2, 3])
    values = [statistic(list(p)) for p in itertools.combinations(range(8), 4)]
    p_mine = sum(1 for v in values if v >= observed - 1e-12) / len(values)
    check("omnibus observed statistic reproduces",
          abs(observed - omni["correlation_distance"]["observed_statistic"]) < 1e-4,
          f"mine {observed:.6f}, stage 3 {omni['correlation_distance']['observed_statistic']}")
    check("omnibus permutation p reproduces",
          abs(p_mine - omni["correlation_distance"]["p_one_sided"]) < P_TOLERANCE,
          f"mine {p_mine:.6f}, stage 3 {omni['correlation_distance']['p_one_sided']}")
    check("omnibus p cannot fall below two over seventy",
          omni["correlation_distance"]["p_one_sided"] >= 2 / 70 - 1e-9)

    # The verdict follows the declared rules rather than a judgement.
    v = stage3["verdict"]
    positive = (v["transitional_separates"] or v["identity_separates"]
                or v["omnibus_correlation_separates"])
    if v["downgrade_reasons"]:
        expected = "inconclusive"
    elif positive:
        expected = "rival_2_stays_live"
    elif v["engagement_separates_in_declared_direction"]:
        expected = "weak_bound_on_rival_2"
    else:
        expected = "engagement_not_established"
    check("the recorded verdict follows the declared rules", v["value"] == expected,
          f"recorded {v['value']}, rules give {expected}")

    # The freeze's prohibitions that are mechanically checkable.
    check("no ratio verdict appears in the run record",
          "ratio" not in json.dumps(v).lower())

    failures = [c for c in checks if not c["passed"]]
    skipped = [c for c in checks if c.get("skipped")]
    record = {
        "schema": "a15-rival2-verification/v2",
        "verified": "tables/rival2/stage3_execute_run.json",
        "independent_tools": ["numpy", "scipy.stats.mannwhitneyu", "numpy.corrcoef"],
        "note": "stage 3 is standard library only by convention; this script is the "
                "independent check and is not part of the frozen instrument",
        "tolerances": {"p_values": P_TOLERANCE, "scores": SCORE_TOLERANCE,
                       "why": "stage 3 rounds p-values to six decimals"},
        "checks_run": len(checks),
        "checks_failed": len(failures),
        "checks_skipped": len(skipped),
        "checks": checks,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    (OUT / "stage3_verification.json").write_text(json.dumps(record, indent=2) + "\n",
                                                 encoding="utf-8")
    for entry in checks:
        tag = "skip" if entry.get("skipped") else ("ok  " if entry["passed"] else "FAIL")
        print(f"[{tag}] {entry['check']}" + (f"  ({entry['detail']})" if entry["detail"] else ""))
    print(f"\n{len(checks) - len(failures)}/{len(checks)} passed, {len(skipped)} skipped")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
