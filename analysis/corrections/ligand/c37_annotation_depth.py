"""Formal post hoc C37 annotation/depth sensitivity; historical files are read only."""
from __future__ import annotations

import gzip
import importlib.metadata
import json
import operator
import time

import numpy as np
import pandas as pd

from common import HERE, ROOT, expected_detection_at_budget, paired_summary, sha256, write_json

OUT = HERE / "results" / "c37"
RAW = ROOT / "raw_data" / "GSE131907"


def extract_depth(meta, raw_path, path, fingerprint):
    sidecar = path.with_suffix(".json")
    if path.exists() and sidecar.exists() and json.loads(sidecar.read_text())["fingerprint"] == fingerprint:
        return pd.read_csv(path, index_col="Index")
    with gzip.open(raw_path, "rt", encoding="utf-8") as f:
        header = f.readline().rstrip("\r\n").split("\t")[1:]
        if len(set(header)) != len(header):
            raise ValueError("raw matrix cell barcodes are not unique")
        position = pd.Index(header).get_indexer(meta.index)
        if np.any(position < 0):
            raise ValueError("annotation cells absent from raw matrix")
        select_values = operator.itemgetter(*map(int, position))
        totals = np.zeros(len(position), dtype=np.int64)
        genes = np.zeros(len(position), dtype=np.int32)
        target = None
        for number, line in enumerate(f, 1):
            gene, _, rest = line.partition("\t")
            # The file is dense text, but only 16k of 208k cells are in scope.
            # Split in C and parse only those counts; verified equal to full-row
            # NumPy parsing on a 50-row raw-matrix benchmark before using it.
            fields = rest.rstrip("\r\n").split("\t")
            if len(fields) != len(header):
                raise ValueError(f"malformed raw matrix row {number}")
            counts = np.array(select_values(fields), dtype=np.int64)
            if np.any(counts < 0):
                raise ValueError(f"negative raw matrix count at row {number}")
            totals += counts
            genes += counts > 0
            if gene == "AREG":
                if target is not None:
                    raise ValueError("duplicate AREG row")
                target = counts.copy()
            if number % 5000 == 0:
                print(f"C37 raw stream: {number} genes", flush=True)
    if target is None:
        raise ValueError("AREG missing from raw matrix")
    result = pd.DataFrame({"total_umi": totals, "n_genes": genes, "areg_umi": target}, index=meta.index)
    result.to_csv(path, index_label="Index")
    write_json(sidecar, {"fingerprint": fingerprint, "gene_rows": number, "raw_cells": len(header),
                         "selected_cells": len(position)})
    return result


def main():
    started = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    inputs = [RAW / "e1_extracted_rows.npz", RAW / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
              RAW / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz", HERE / "specification.json",
              HERE / "common.py", __file__]
    hashes = {str(p): sha256(p) for p in inputs}
    fingerprint = hashes[str(inputs[2])] + hashes[str(inputs[1])]
    meta = pd.read_csv(inputs[1], sep="\t", index_col="Index")
    if not meta.index.is_unique:
        raise ValueError("duplicate deposited annotation cell index")
    meta = meta.loc[meta.Sample_Origin.eq("tLung") & meta["Cell_type.refined"].isin(
        ["Epithelial cells", "Myeloid cells"])].copy()
    meta["compartment"] = meta["Cell_type.refined"].map(
        {"Epithelial cells": "epithelial", "Myeloid cells": "myeloid"})
    with np.load(inputs[0], allow_pickle=False) as cached:
        idx = pd.Index(cached["cells"]).get_indexer(meta.index)
        if np.any(idx < 0):
            raise ValueError("annotation cell missing in E1 cache")
        meta["original_detected"] = cached["AREG"][idx] > 0
    raw = extract_depth(meta, inputs[2], OUT / "selected_cell_depth.csv.gz", fingerprint)
    meta = meta.join(raw, validate="one_to_one")
    meta["raw_detected"] = meta.areg_umi > 0
    records = []
    for (donor, compartment), frame in meta.groupby(["Sample", "compartment"]):
        row = {"donor": donor, "compartment": compartment, "n_cells": len(frame),
               "median_umi": float(frame.total_umi.median()), "median_genes": float(frame.n_genes.median()),
               "original_detection": float(frame.original_detected.mean()) if len(frame) >= 50 else np.nan,
               "raw_detection": float(frame.raw_detected.mean()) if len(frame) >= 50 else np.nan}
        for budget in [500, 1000, 2000]:
            prob = expected_detection_at_budget(frame.total_umi, frame.areg_umi, budget)
            n = int(np.isfinite(prob).sum())
            row[f"n_at_{budget}"] = n
            row[f"detection_at_{budget}"] = float(np.nanmean(prob)) if n >= 50 else np.nan
        records.append(row)
    table = pd.DataFrame(records)
    table.to_csv(OUT / "per_donor.csv", index=False)
    results = {column: paired_summary(table, column) for column in
               ["original_detection", "raw_detection", "detection_at_500", "detection_at_1000", "detection_at_2000"]}
    results["cache_raw_detection_disagreements"] = int((meta.original_detected != meta.raw_detected).sum())
    results["selected_cells"] = len(meta)
    write_json(OUT / "summary.json", results)
    write_json(OUT / "run_record.json", {"analysis": "Post hoc annotation/depth sensitivity, not confirmation",
        "inputs_sha256": hashes, "versions": {p: importlib.metadata.version(p) for p in ["numpy", "pandas", "scipy"]},
        "elapsed_seconds": round(time.time() - started, 1), "specification": "../../specification.json"})
    print(json.dumps(results, indent=2), flush=True)


if __name__ == "__main__":
    main()
