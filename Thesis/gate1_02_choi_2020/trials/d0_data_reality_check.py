#!/usr/bin/env python
"""Trial D0: what the Choi-2020 deposit actually contains (Gate 0).

Gate 0 reads the deposit before anything is fitted, so that every later trial
in this folder inherits facts rather than assumptions. It answers the questions
that decided what was possible in the Cardoso folder: how many libraries,
what design, are the matrices called cells or raw droplets, is the gene space
the same across libraries, are the paper's own marker genes present, and does
any contrast carry within-group replication.

Nothing here fits a model, clusters anything, or computes a statistic on a
biological comparison. The trial deliberately runs before the study note's
claims are tested, and its outputs are the inputs to Gate 1.

Frozen rules, set before any file was opened:

* Scope: the two single-cell accessions, GSE145031 (six lineage-traced
  libraries) and GSE144468 (two organoid libraries). GSE144598 is ATAC-seq
  deposited as bigwig coverage only and is reported as unusable rather than
  parsed; GSE144553 is the SuperSeries and holds no files of its own.
* Design comes from the deposited GEO SOFT family files, parsed, not
  transcribed by hand.
* FEATURE RULE. A feature whose identifier is not an Ensembl gene ID is not a
  gene. Any such feature is counted and reported separately, because a
  lineage-reporter contig counted as a gene would corrupt every later marker
  score. This is the rule trial C0 of the Cardoso folder applied to the BSD
  selection marker.
* CELL-CALLING RULE. A matrix whose column count equals the 10x v2 barcode
  whitelist (737,280) is the raw droplet matrix, not called cells. Where that
  is so, the trial reports it as a fact about the deposit and does NOT call
  cells itself, because cell calling is a Gate 1 decision that belongs with a
  frozen threshold rather than with an inventory.
* Gene spaces are compared as the exact ordered list of (Ensembl ID, symbol)
  pairs. Two libraries share a gene space only if those lists are identical.
* Marker presence is checked for every set in `choi_2020_extracts.json`,
  including the negative conditions, and reported per accession.
* REPLICATE RULE, frozen here so no later trial can quietly assume otherwise:
  a contrast has within-group replication only if at least two libraries share
  every experimental variable except the one being contrasted. This is
  evaluated and written into the run record.
* The trial reports total non-zero entries and per-barcode count sums from the
  MatrixMarket header and a streaming pass, so nothing large is held in memory
  and no library is loaded as a dense object.
"""

from __future__ import annotations

import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from choi_utils import (ACCESSIONS, ENSEMBL_ID_RE, RAW, TENX_V2_WHITELIST,  # noqa: E402
                        RunRecord, df_to_markdown, library_paths, parse_soft)

OUT = HERE / "d0_data_reality_check"
OUT.mkdir(exist_ok=True)
EXTRACTS = HERE.parent / "choi_2020_extracts.json"
UNUSABLE = {
    "GSE144598": "ATAC-seq deposited as bigwig coverage tracks only; no peaks and no reads, "
                 "so the epigenetic claim cannot be re-derived without going to SRA",
    "GSE144553": "SuperSeries over the other three; holds no files of its own",
}

RULES = {
    "question": "Gate 0: what the Choi-2020 deposit contains, before anything is fitted",
    "scope": "GSE145031 (6 lineage-traced libraries) and GSE144468 (2 organoid libraries)",
    "out_of_scope": UNUSABLE,
    "design_source": "the deposited GEO SOFT family files, parsed rather than transcribed",
    "feature_rule": "a feature whose identifier is not an Ensembl gene ID is not a gene; counted and reported separately",
    "cell_calling_rule": f"a matrix with {TENX_V2_WHITELIST} columns is the raw droplet matrix, not called cells; "
                         "reported as a fact, and cell calling is left to Gate 1 where a threshold can be frozen",
    "gene_space_comparison": "exact ordered list of (Ensembl ID, symbol) pairs; identical or not",
    "marker_presence": "every set in choi_2020_extracts.json, negative conditions included",
    "replicate_rule": "a contrast has within-group replication only if at least two libraries share every "
                      "experimental variable except the one being contrasted",
    "memory": "MatrixMarket headers and one streaming pass; no library is loaded densely",
}


def read_features(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            rows.append((parts[0], parts[1] if len(parts) > 1 else ""))
    return pd.DataFrame(rows, columns=["gene_id", "gene_symbol"])


def count_lines(path: Path) -> int:
    n = 0
    with gzip.open(path, "rt") as handle:
        for _ in handle:
            n += 1
    return n


def matrix_facts(path: Path) -> dict:
    """Header plus a streaming pass: shape, non-zeros, and per-barcode totals."""
    with gzip.open(path, "rt") as handle:
        shape = None
        for line in handle:
            if line.startswith("%"):
                continue
            shape = [int(v) for v in line.split()]
            break
        if shape is None:
            raise ValueError(f"{path} has no dimension line")
        n_genes, n_barcodes, n_entries = shape[0], shape[1], shape[2]
        per_barcode = np.zeros(n_barcodes + 1, dtype=np.int64)
        genes_per_barcode = np.zeros(n_barcodes + 1, dtype=np.int32)
        seen = 0
        for line in handle:
            _, col, val = line.split()
            index = int(col)
            per_barcode[index] += int(val)
            genes_per_barcode[index] += 1
            seen += 1
    occupied = int((per_barcode > 0).sum())
    counts = per_barcode[per_barcode > 0]
    return {
        "n_features": n_genes, "n_barcodes": n_barcodes, "n_entries_declared": n_entries,
        "n_entries_read": seen, "barcodes_with_any_count": occupied,
        "median_counts_per_occupied_barcode": float(np.median(counts)) if occupied else 0.0,
        "barcodes_over_500_counts": int((per_barcode >= 500).sum()),
        "barcodes_over_1000_counts": int((per_barcode >= 1000).sum()),
        "median_genes_per_barcode_over_500": float(
            np.median(genes_per_barcode[per_barcode >= 500])) if (per_barcode >= 500).any() else 0.0,
    }


def main():
    rec = RunRecord(OUT / "d0_run_record.json", "D0 Choi-2020 data reality check", RULES)
    rec.add_input(EXTRACTS)
    extracts = json.loads(EXTRACTS.read_text(encoding="utf-8"))

    # --- design, from the SOFT family files -------------------------------
    design_rows = []
    for accession in ACCESSIONS:
        soft = RAW / accession / (accession + "_family.soft.gz")
        if not soft.exists():
            raise SystemExit(f"{soft} is missing; Gate 0 will not guess a design")
        rec.add_input(soft)
        parsed = parse_soft(soft)
        for sample in parsed.get("samples", []):
            design_rows.append({
                "accession": accession,
                "gsm": sample.get("geo_accession", ""),
                "title": sample.get("title", ""),
                "characteristics": " | ".join(sample.get("characteristics", []))[:220],
                "data_processing": " | ".join(sample.get("data_processing", []))[:220],
            })
    design = pd.DataFrame(design_rows)
    design.to_csv(OUT / "d0_deposited_design.csv", index=False)
    rec.add_output(OUT / "d0_deposited_design.csv")

    # --- inventory, one streaming pass per library ------------------------
    inventory, features_by_library = [], {}
    for accession, spec in ACCESSIONS.items():
        for entry in spec["libraries"]:
            paths = library_paths(accession, entry["library"], entry["gsm"])
            for path in paths.values():
                if not path.exists():
                    raise SystemExit(f"{path} is missing; unpack the series tar before running Gate 0")
                rec.add_input(path)
            var = read_features(paths["features"])
            features_by_library[entry["library"]] = var
            n_barcodes = count_lines(paths["barcodes"])
            facts = matrix_facts(paths["matrix"])
            non_gene = var[~var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False)]
            row = {"accession": accession, "library": entry["library"], "gsm": entry["gsm"]}
            row.update({k: v for k, v in entry.items() if k not in ("library", "gsm")})
            row.update(facts)
            row["barcodes_in_tsv"] = n_barcodes
            row["is_raw_whitelist"] = bool(facts["n_barcodes"] == TENX_V2_WHITELIST)
            row["non_gene_features"] = int(len(non_gene))
            row["non_gene_feature_names"] = ", ".join(non_gene["gene_id"].head(6))
            inventory.append(row)
            print(entry["library"], facts["n_features"], "features,", facts["n_barcodes"],
                  "barcodes, raw whitelist:", row["is_raw_whitelist"])
    inv = pd.DataFrame(inventory)
    inv.to_csv(OUT / "d0_library_inventory.csv", index=False)
    rec.add_output(OUT / "d0_library_inventory.csv")
    rec.set("raw_whitelist_libraries", inv.loc[inv["is_raw_whitelist"], "library"].tolist())
    rec.set("already_filtered_libraries", inv.loc[~inv["is_raw_whitelist"], "library"].tolist())
    rec.set("non_gene_features_total", int(inv["non_gene_features"].sum()))

    # --- gene spaces ------------------------------------------------------
    signatures = {}
    for library, var in features_by_library.items():
        key = hash(tuple(zip(var["gene_id"], var["gene_symbol"])))
        signatures.setdefault(key, []).append(library)
    spaces = pd.DataFrame([
        {"gene_space": i + 1, "n_features": len(features_by_library[libs[0]]),
         "n_libraries": len(libs), "libraries": ", ".join(libs)}
        for i, (_, libs) in enumerate(signatures.items())
    ])
    spaces.to_csv(OUT / "d0_gene_spaces.csv", index=False)
    rec.add_output(OUT / "d0_gene_spaces.csv")
    rec.set("distinct_gene_spaces", int(len(signatures)))

    # --- marker presence --------------------------------------------------
    sets = extracts["marker_sets"]
    flat = {}
    for name, value in sets.items():
        if isinstance(value, list):
            flat[name] = value
        elif isinstance(value, dict):
            for sub, genes in value.items():
                if isinstance(genes, list):
                    flat[name + "." + sub] = genes
    presence = []
    for library, var in features_by_library.items():
        symbols = set(var["gene_symbol"])
        for name, genes in flat.items():
            absent = [g for g in genes if g not in symbols]
            presence.append({"library": library, "marker_set": name, "n_genes": len(genes),
                             "n_absent": len(absent), "absent": ", ".join(absent)})
    pres = pd.DataFrame(presence)
    pres.to_csv(OUT / "d0_marker_presence.csv", index=False)
    rec.add_output(OUT / "d0_marker_presence.csv")
    worst = pres[pres["n_absent"] > 0]
    rec.set("marker_sets_with_absent_genes",
            sorted(set(worst["marker_set"])) if len(worst) else [])

    # --- the replicate rule, applied -------------------------------------
    verdicts = []
    lineage = inv[inv["accession"] == "GSE145031"]
    for variable, holding in (("timepoint", "sort"), ("sort", "timepoint")):
        groups = Counter(zip(lineage[holding], lineage[variable]))
        replicated = [k for k, v in groups.items() if v > 1]
        verdicts.append({"accession": "GSE145031", "contrast": variable,
                         "held_constant": holding,
                         "cells_with_more_than_one_library": len(replicated),
                         "has_within_group_replication": bool(replicated)})
    organoid = inv[inv["accession"] == "GSE144468"]
    verdicts.append({"accession": "GSE144468", "contrast": "treatment",
                     "held_constant": "none",
                     "cells_with_more_than_one_library":
                         int(sum(1 for _, v in Counter(organoid["treatment"]).items() if v > 1)),
                     "has_within_group_replication": bool(
                         any(v > 1 for v in Counter(organoid["treatment"]).values()))})
    rep = pd.DataFrame(verdicts)
    rep.to_csv(OUT / "d0_replication.csv", index=False)
    rec.add_output(OUT / "d0_replication.csv")
    rec.set("any_contrast_with_within_group_replication",
            bool(rep["has_within_group_replication"].any()))

    raw_libs = inv.loc[inv["is_raw_whitelist"], "library"].tolist()
    filtered_libs = inv.loc[~inv["is_raw_whitelist"], "library"].tolist()
    lines = ["# Trial D0: what the Choi-2020 deposit contains", "",
             "Raw 10x barcode whitelists rather than called cells: "
             + str(len(raw_libs)) + " of " + str(len(inv)) + " libraries ("
             + ", ".join(raw_libs) + "). Already filtered to called cells: "
             + (", ".join(filtered_libs) if filtered_libs else "none") + ".",
             "Distinct gene spaces across the eight libraries: " + str(len(signatures)) + ".",
             "Non-gene features found: " + str(int(inv["non_gene_features"].sum())) + ".",
             "Any contrast with within-group replication: "
             + str(bool(rep["has_within_group_replication"].any())) + ".", "",
             "## Accessions not parsed, and why", ""]
    for accession, reason in UNUSABLE.items():
        lines.append("- **" + accession + "**: " + reason)
    lines += ["", "## Library inventory", "",
              df_to_markdown(inv[["accession", "library", "n_features", "n_barcodes",
                                  "barcodes_with_any_count", "barcodes_over_500_counts",
                                  "median_genes_per_barcode_over_500", "is_raw_whitelist",
                                  "non_gene_features"]], index=False), "",
              "## Gene spaces", "", df_to_markdown(spaces, index=False), "",
              "## The replicate rule, applied", "", df_to_markdown(rep, index=False), "",
              "## Marker sets with any gene absent from a library", "",
              df_to_markdown(worst, index=False) if len(worst)
              else "None. Every marker set in the extract is fully present in every library.", ""]
    (OUT / "d0_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "d0_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
