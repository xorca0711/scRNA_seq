"""Validate two public sparse-count inputs without estimating a biological effect."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys
import urllib.request

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import archive_existing_record, code_identity, sha256_file, write_json_atomic


def validate_matrix(path, features, barcodes):
    with gzip.open(features, "rt", encoding="utf-8") as f:
        ids = [line.rstrip().split("\t")[0] for line in f if line.strip()]
    with gzip.open(barcodes, "rt", encoding="utf-8") as f:
        bc = [line.strip() for line in f if line.strip()]
    if len(ids) != len(set(ids)) or len(bc) != len(set(bc)):
        raise ValueError("Nonunique feature IDs or barcodes")
    with gzip.open(path, "rt", encoding="utf-8") as f:
        header = f.readline().strip()
        if header not in ("%%MatrixMarket matrix coordinate integer general", "%%MatrixMarket matrix coordinate real general"):
            raise ValueError(f"Unsupported matrix header: {header}")
        line = next(f)
        while line.startswith("%"):
            line = next(f)
        nfeatures, nbarcodes, nnz = map(int, line.split())
        if (nfeatures, nbarcodes) != (len(ids), len(bc)):
            raise ValueError("Matrix dimension mismatch")
        depth = [0] * nbarcodes
        observed = 0
        for line in f:
            if not line.strip():
                continue
            row, col, value = line.split()
            i, j, x = int(row), int(col), float(value)
            if not 1 <= i <= nfeatures or not 1 <= j <= nbarcodes:
                raise ValueError("Out-of-range matrix coordinate")
            if not math.isfinite(x) or x < 0 or x != int(x):
                raise ValueError("Not a nonnegative integer count")
            depth[j-1] += int(x)
            observed += 1
        if observed != nnz:
            raise ValueError("Declared matrix entry count mismatch")
    return {"features": nfeatures, "barcodes": nbarcodes, "stored_entries": nnz,
            "total_counts": sum(depth), "zero_count_barcodes": sum(x == 0 for x in depth),
            "matrix_format": header, "status": "input_checks_passed"}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--download", action="store_true")
    args = p.parse_args()
    spec_path = PAPER / "trials/u1_input_spec.json"
    spec = json.loads(spec_path.read_text(encoding="utf-8"))
    metadata_path = PAPER / "trials/u0_geo_design_audit" / (spec["source_series"] + ".json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    by_id = {s["accession"]: s for s in metadata["samples"]}
    cache = PAPER / "cache/count_inputs"
    cache.mkdir(parents=True, exist_ok=True)
    out = PAPER / "trials/u1_count_input_preflight"
    out.mkdir(parents=True, exist_ok=True)
    files, results, total = [], [], 0
    for gsm in spec["samples"]:
        selected = {}
        for kind, suffix in (("matrix", "_matrix.mtx.gz"), ("features", "_features.tsv.gz"), ("barcodes", "_barcodes.tsv.gz")):
            urls = [u for u in by_id[gsm]["supplementary_file"] if u.endswith(suffix)]
            if len(urls) != 1:
                raise ValueError(f"Expected one {kind} for {gsm}")
            url = urls[0].replace("ftp://", "https://", 1)
            path = cache / url.rsplit("/", 1)[1]
            if args.download and not path.exists():
                partial = path.with_suffix(path.suffix + ".part")
                try:
                    with urllib.request.urlopen(url, timeout=60) as response, partial.open("wb") as dest:
                        size = 0
                        while chunk := response.read(1024 * 1024):
                            size += len(chunk)
                            if size > spec["maximum_compressed_bytes_per_file"] or total + size > spec["maximum_total_compressed_bytes"]:
                                raise ValueError("Frozen download cap exceeded")
                            dest.write(chunk)
                    partial.replace(path)
                finally:
                    if partial.exists():
                        partial.unlink()
            size = path.stat().st_size
            total += size
            if size > spec["maximum_compressed_bytes_per_file"] or total > spec["maximum_total_compressed_bytes"]:
                raise ValueError("Input exceeds frozen cap")
            files.append({"gsm": gsm, "kind": kind, "url": url, "bytes": size, "sha256": sha256_file(path)})
            selected[kind] = path
        result = {"gsm": gsm, **validate_matrix(selected["matrix"], selected["features"], selected["barcodes"])}
        results.append(result)
        print(json.dumps(result), flush=True)
    with (out / "input_qc.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    record = {"run_utc": datetime.now(timezone.utc).isoformat(), "code": code_identity(ROOT, __file__),
              "spec_sha256": sha256_file(spec_path), "metadata_sha256": sha256_file(metadata_path),
              "inputs": files, "total_compressed_bytes": total, "biological_effects_estimated": False,
              "limits": "Two library examples; not paired animals; no between-arm inference; coordinate duplicates are not tested."}
    archive_existing_record(out / "run_record.json")
    write_json_atomic(out / "run_record.json", record)


if __name__ == "__main__":
    main()
