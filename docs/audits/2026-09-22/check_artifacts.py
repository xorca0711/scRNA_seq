"""Read-only audit of selected claims; writes only audit_checks.json beside this file.

This is an audit, not a new biological trial or a change to retained decisions.
Standard library checks work without the scientific environment. NumPy enables
an explicitly post hoc annotation sensitivity check using E1's existing cache.
"""
from __future__ import annotations

import ast
import collections
import csv
import gzip
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import statistics
import subprocess
import sys

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).with_name("audit_checks.json")
INPUTS = {}


def record(path):
    p = REPO / path
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    INPUTS[path] = {"bytes": p.stat().st_size, "sha256": h.hexdigest()}
    return p


def rows(path):
    with record(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def signed_rank(differences):
    """Exact two-sided sign enumeration; average ranks for tied magnitudes."""
    ds = [x for x in differences if x != 0]
    ordered = sorted(range(len(ds)), key=lambda i: abs(ds[i]))
    ranks = [0.0] * len(ds)
    j = 0
    while j < len(ds):
        k = j + 1
        while k < len(ds) and abs(ds[ordered[k]]) == abs(ds[ordered[j]]):
            k += 1
        for pos in ordered[j:k]:
            ranks[pos] = (j + 1 + k) / 2
        j = k
    total = sum(ranks)
    observed = abs(sum(r for r, d in zip(ranks, ds) if d > 0) - total / 2)
    extreme = sum(abs(sum(r for r, b in zip(ranks, bits) if b) - total / 2)
                  >= observed - 1e-10
                  for bits in itertools.product((False, True), repeat=len(ds)))
    return extreme / 2 ** len(ds)


def paired_summary(groups):
    pairs = [g for g in groups.values() if "epithelial" in g and "myeloid" in g]
    a = [g["epithelial"] for g in pairs]
    b = [g["myeloid"] for g in pairs]
    return {"n": len(pairs), "median_epithelial": statistics.median(a),
            "median_myeloid": statistics.median(b),
            "epithelial_higher": sum(x > y for x, y in zip(a, b)),
            "exact_p_two_sided": signed_rank([x - y for x, y in zip(a, b)])}


def main():
    result = {"scope": "Repository audit; no claim statuses changed",
              "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=REPO, text=True).strip(),
              "python": sys.version}
    text = record("CLAIMS.md").read_text(encoding="utf-8")
    claims = {}
    for line in text.splitlines():
        if re.match(r"^\| C\d+ \|", line):
            cells = [c.strip() for c in line.split("|")]
            claims[cells[1]] = cells
    retained = {k: v for k, v in claims.items() if "owner retained" in v[5]}
    block = [k for k, v in retained.items() if "as a block" in v[5]]
    missing = []
    for k, v in retained.items():
        for p in re.findall(r"`([^`]+)`", v[4]):
            if not (REPO / p).exists():
                missing.append({"claim": k, "path": p})
    result["claims"] = {"rows": len(claims), "explicit_owner_retained": len(retained),
                        "block_to_revisit": block, "missing_retained_artifacts": missing}
    files = list((REPO / "analysis/scripts").glob("*.py")) + list((REPO / "Thesis").rglob("*.py"))
    for p in files:
        ast.parse(p.read_text(encoding="utf-8-sig"), filename=str(p))
    result["python_sources_parsed"] = len(files)

    card = "Thesis/gate2_05_cardoso_2026/trials/"
    e1 = card + "e1_human_luad_ligand_sources/"
    cross = rows(e1 + "e1_gate_vs_deposited_celltype.csv")
    row = next(r for r in cross if r["compartment"] == "myeloid")
    total = sum(int(v) for k, v in row.items() if k != "compartment")
    result["E1_pooled_normal_and_tumour_myeloid_gate"] = {
        "n": total, "author_myeloid": int(row["Myeloid cells"]),
        "author_myeloid_pct": 100 * int(row["Myeloid cells"]) / total,
        "author_T_NK": int(row["T/NK cells"]),
        "author_T_NK_pct": 100 * int(row["T/NK cells"]) / total}
    groups = collections.defaultdict(dict)
    for r in rows(e1 + "e1_per_donor_compartment.csv"):
        if r["tissue"] == "tLung" and r["evaluable"] == "True" and r["compartment"] in ("epithelial", "myeloid"):
            groups[r["donor"]][r["compartment"]] = float(r["det_AREG"])
    result["C37_original_table_recomputed"] = paired_summary(groups)

    # Optional raw-cache cross-check, not a confirmatory reanalysis.
    try:
        import numpy as np
    except ImportError:
        result["E1_annotation_sensitivity"] = "not run: NumPy unavailable"
    else:
        cache = record("raw_data/GSE131907/e1_extracted_rows.npz")
        meta = record("raw_data/GSE131907/GSE131907_Lung_Cancer_cell_annotation.txt.gz")
        with np.load(cache, allow_pickle=False) as z:
            cells = z["cells"]
            pos = {str(c): i for i, c in enumerate(cells)}
            areg = z["AREG"] > 0
            immune = z["PTPRC"] > 0
            myeloid_hit = np.logical_or.reduce([z[g] > 0 for g in ("LYZ", "CD68", "ITGAM", "MARCO", "CD14")])
            neutro_hit = np.logical_or.reduce([z[g] > 0 for g in ("S100A8", "S100A9", "FCGR3B")])
            is_original_myeloid = immune & myeloid_hit & ~neutro_hit
        counts = collections.defaultdict(lambda: [0, 0])
        gate_types = collections.Counter()
        with gzip.open(meta, "rt", encoding="utf-8-sig", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r["Sample_Origin"] != "tLung" or r["Index"] not in pos:
                    continue
                i = pos[r["Index"]]
                label = r["Cell_type.refined"]
                if is_original_myeloid[i]:
                    gate_types[label] += 1
                comp = {"Epithelial cells": "epithelial", "Myeloid cells": "myeloid"}.get(label)
                if comp:
                    key = (r["Sample"], comp)
                    counts[key][0] += 1
                    counts[key][1] += int(areg[i])
        grouped = collections.defaultdict(dict)
        for (donor, comp), (n, hits) in counts.items():
            if n >= 50:
                grouped[donor][comp] = hits / n
        result["E1_tumour_only_original_myeloid_gate"] = {
            "n": sum(gate_types.values()), "by_deposited_type": dict(gate_types)}
        result["E1_annotation_sensitivity"] = {
            "scope": "Post hoc diagnostic; E1 cached AREG values, deposited major cell types; >=50 cells per donor/compartment; no depth adjustment",
            **paired_summary(grouped),
            "per_donor": dict(grouped)}

    nieth = "Thesis/gate1_01_niethamer_2025/trials/"
    g2 = rows(nieth + "g2_gsea_ipf/g2_replication.csv")
    result["C161"] = {"discovery_passed": len(g2),
        "replicated": sum(r["replicates"] == "True" for r in g2),
        "by_compartment": {comp: {"discovery": sum(r["compartment"] == comp for r in g2),
                            "replicated": sum(r["compartment"] == comp and r["replicates"] == "True" for r in g2)}
                           for comp in sorted({r["compartment"] for r in g2})}}
    w1 = rows(nieth + "w1_amac_pseudobulk_de/w1_units.csv")
    result["W1_kept_animals"] = dict(collections.Counter(r["arm"] for r in w1 if r["kept"] == "True"))
    result["W1_genotype_exclusions"] = sum(r["reason"] == "Ki67Cre/Cre" for r in w1)
    animal_meta = {r["animal"]: r for r in w1}
    g1 = rows(nieth + "g1_gsea_by_phase/g1_units.csv")
    result["G1_myeloid_design"] = {}
    for phase in ("injury resolution", "long-term"):
        used = [animal_meta[r["animal"]] for r in g1
                if r["compartment"] == "myeloid" and r["phase"] == phase and r["kept"] == "True"]
        result["G1_myeloid_design"][phase] = {
            "n": len(used), "sex": dict(collections.Counter(r["sex"] for r in used)),
            "round": dict(collections.Counter(r["round"] for r in used)),
            "Ki67Cre_Cre": sum("Cre/Cre" in r["genotype"] for r in used)}
    sets = rows(nieth + "w1_amac_pseudobulk_de/w1_sets.csv")
    result["W1_sets_passing_stored_fdr"] = sum(r["fdr"] and float(r["fdr"]) < .05 for r in sets if r["tested"] == "True")
    a4 = rows("analysis/figures/rq/rq_a4_detection.csv")
    result["RQ_A4_figure_Il1r1_detection_range_pct"] = [
        min(float(r["pct_detected"]) for r in a4 if r["gene"] == "Il1r1"),
        max(float(r["pct_detected"]) for r in a4 if r["gene"] == "Il1r1")]
    m3 = rows("Thesis/gate1_02_choi_2020/datp_epigenetics/trials/m3_one_instrument_and_the_right_null/m3_arms.csv")
    result["M3_AT2_estimand"] = [{"well": r["well"], "role": r["role"],
        "raw_difference_relative_pct": 100 * float(r["observed_raw"]) / float(r["reference_baseline"]),
        "null_centered_relative_pct": float(r["relative_pct"])} for r in m3 if r["arm"] == "AT2_identity"]
    result["C49_illustrative_not_formal_Spearman_bounds"] = {
        "scope": "Fisher-z/Pearson approximations; not a replacement for donor-level Spearman bootstrap or power simulation",
        "rho_observed": .3484, "n": 22,
        "approx_interval": [math.tanh(math.atanh(.3484) + s * 1.96 / math.sqrt(19)) for s in (-1, 1)],
        "approx_80pct_power_rho": math.tanh((1.96 + statistics.NormalDist().inv_cdf(.8)) / math.sqrt(19))}
    result["inputs"] = INPUTS
    OUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    compact = {k: v for k, v in result.items() if k not in ("inputs", "python")}
    if isinstance(compact.get("E1_annotation_sensitivity"), dict):
        compact["E1_annotation_sensitivity"] = {k: v for k, v in compact["E1_annotation_sensitivity"].items() if k != "per_donor"}
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
