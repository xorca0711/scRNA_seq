"""Rebuild compact metadata evidence; no expression matrices or network required."""
from collections import Counter, defaultdict
from pathlib import Path
import csv
import gzip
import hashlib
import json
import re

HERE = Path(__file__).resolve().parent


def soft_samples(accession):
    rows = []
    text = gzip.decompress((HERE / f"{accession}_family.soft.gz").read_bytes()).decode()
    for block in text.split("^SAMPLE = ")[1:]:
        lines = block.splitlines()
        fields = defaultdict(list)
        for line in lines[1:]:
            if " = " in line:
                k, v = line.split(" = ", 1)
                fields[k].append(v)
        rows.append({"gsm": lines[0], "title": fields["!Sample_title"][0],
                     "characteristics": fields["!Sample_characteristics_ch1"],
                     "files": [v for k, vals in fields.items()
                               if k.startswith("!Sample_supplementary_file") for v in vals]})
    return rows


def write_csv(name, rows):
    with (HERE / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    file_sizes = {}
    with (HERE / "GSE129605_filelist.txt").open() as f:
        for row in csv.reader(f, delimiter="\t"):
            if row[0] == "File":
                file_sizes[row[1]] = int(row[3])
    sample_rows, manifest = [], []
    for sample in soft_samples("GSE129605"):
        arm, animal = sample["title"].split(":")[0].split(".")
        selected = arm in ("saline", "bleomycin")
        sample_rows.append({"gsm": sample["gsm"], "sample_id": f"{arm}.{animal}",
                            "deposited_numeric_id": animal, "arm": arm, "day": 11,
                            "sex": "male", "strain": "C57BL/6", "age_weeks": "12-13",
                            "selected": selected, "unit_status": "separately deposited sample; mouse mapping requires final provenance check"})
        if selected:
            for url in sample["files"]:
                name = url.rsplit("/", 1)[-1]
                manifest.append({"gsm": sample["gsm"], "sample_id": f"{arm}.{animal}",
                                 "name": name, "url": url.replace("ftp://", "https://"),
                                 "expected_bytes": file_sizes[name]})
    assert Counter(r["arm"] for r in sample_rows if r["selected"]) == {"saline": 4, "bleomycin": 4}
    assert len(manifest) == 24
    write_csv("GSE129605_sample_map.csv", sample_rows)
    write_csv("GSE129605_download_manifest.csv", manifest)

    by_id = {}
    for sample in soft_samples("GSE141259"):
        match = re.search(r"\[(muc\d+)\]", sample["title"])
        if match:
            by_id[match.group(1)] = sample["gsm"]
    counts, group = defaultdict(Counter), {}
    with gzip.open(HERE / "GSE141259_WholeLung_cellinfo.csv.gz", "rt") as f:
        for row in csv.DictReader(f):
            sid = row["orig.ident"]
            counts[sid][row["cell.type"]] += 1
            group[sid] = row["identifier"].removeprefix(sid + "_")
    coverage = []
    for sid, ct in sorted(counts.items()):
        at2 = ct["AT2 cells"] + ct["Activated AT2 cells"]
        fb = ct["Fibroblasts"] + ct["Myofibroblasts"]
        coverage.append({"sample_id": sid, "gsm": by_id.get(sid, "UNMAPPED"), "group": group[sid],
                         "all_cells": sum(ct.values()), "AT2": ct["AT2 cells"],
                         "activated_AT2": ct["Activated AT2 cells"],
                         "fibroblasts": ct["Fibroblasts"], "myofibroblasts": ct["Myofibroblasts"],
                         "AT2_including_activated": at2, "fibroblasts_including_myofibroblasts": fb,
                         "both_at_least_50": at2 >= 50 and fb >= 50})
    assert len(coverage) == 28 and sum(r["all_cells"] for r in coverage) == 29297
    write_csv("GSE141259_compartment_coverage.csv", coverage)
    registry = {
        "screen_date": "2026-09-22", "scope": "Independent mouse AT2/fibroblast Wnt expression generalization; no functional autocrine inference",
        "floors": {"biological_units_per_arm": 3, "cells_per_compartment_per_unit": 50,
                   "note": "Feasibility minima, not a power calculation. Must not lower after viewing Wnt effects."},
        "selected_candidate": {"accession": "GSE129605", "decision": "GO_ACQUISITION_HOLD_INFERENCE",
            "contrast": "day11 bleomycin versus saline; exclude nintedanib", "deposited_samples_per_arm": 4,
            "selected_files": len(manifest), "selected_download_bytes": sum(r["expected_bytes"] for r in manifest),
            "complete_archive_bytes": 118149120,
            "counts_downloaded": (HERE / "GSE129605_download_receipt.json").exists(),
            "remaining_gates": ["Confirm each separately indexed sample is one independent mouse, not a pooled or technical replicate",
                                "Freeze annotation independent of Wnt genes; count AT2 and fibroblasts per mouse after QC",
                                "Require at least 3 mice in each arm each meeting both 50-cell floors",
                                "Establish Wnt feature coverage/detection and batch metadata before effects"],
            "source": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129605",
            "paper_doi": "10.1165/rcmb.2018-0313OC", "limitation": "Day11 is outside the preferred acute1-3day window; bleomycin is not viral injury"},
        "rejected_or_incomplete_alternatives": [
            {"accession": "GSE141259", "decision": "HOLD_CELL_COVERAGE", "reason": "Day3 three bleomycin samples have 4/16/3 fibroblasts including myofibroblasts; all seven PBS samples have <50; only one day3 PBS sample"},
            {"accession": "GSE202325", "decision": "HOLD_NO_UNINFECTED_CONTROL", "reason": "3 young and3 aged per day3/day9; all infected PR8. Supports age contrasts, not injury versus uninjured"},
            {"accession": "GSE292515", "decision": "HOLD_CONTROL_REPLICATION", "reason": "Day7 PR8: each age stratum has 3infected versus2PBS; do not pool ages to pass floor"},
            {"accession": "GSE184854", "decision": "HOLD_DEMULTIPLEXING_AND_CONTROL_UNVERIFIED", "reason": "Two GEO records are WT/CCR2KO mixtures of days3/7/21; sample titles alone are not animal units. Hashtag map not inspected; no GO inferred"}],
        "source_files": []}
    for name in sorted([*HERE.glob("*_family.soft.gz"), HERE/"GSE129605_filelist.txt", HERE/"GSE141259_WholeLung_cellinfo.csv.gz"]):
        accession = name.name.split("_")[0]
        kind = "soft" if name.name.endswith("_family.soft.gz") else "suppl"
        remote_name = "filelist.txt" if name.name.endswith("_filelist.txt") else name.name
        source_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/{kind}/{remote_name}"
        registry["source_files"].append({"name": name.name, "bytes": name.stat().st_size,
            "sha256": hashlib.sha256(name.read_bytes()).hexdigest(), "source_url": source_url})
    audit_path = HERE / "GSE129605_acquisition_audit.json"
    if audit_path.exists():
        registry["selected_candidate"]["decision"] = "ACQUIRED_HOLD_ANNOTATION_AND_ANIMAL_PROVENANCE"
        registry["selected_candidate"]["acquisition_audit"] = json.loads(audit_path.read_text())
    (HERE / "feasibility.json").write_text(json.dumps(registry, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"candidate": registry["selected_candidate"]["decision"], "download_bytes": registry["selected_candidate"]["selected_download_bytes"],
                      "coverage_rows": len(coverage), "day3": [r for r in coverage if r["group"] == "Bleo_d3"]}, indent=2))


if __name__ == "__main__":
    main()
