"""Validate acquired dimensions and count semantics, without annotating cells."""
from pathlib import Path
from collections import defaultdict
import csv
import gzip
import hashlib
import json

HERE = Path(__file__).resolve().parent
RAW = HERE.parents[2] / "raw_data" / "GSE129605"


def main():
    receipt = json.loads((HERE / "GSE129605_download_receipt.json").read_text())
    grouped = defaultdict(dict)
    for row in receipt:
        path = RAW / row["name"]
        assert path.stat().st_size == int(row["expected_bytes"])
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["sha256"]
        kind = "matrix" if "matrix.mtx" in row["name"] else "genes" if "genes.tsv" in row["name"] else "barcodes"
        grouped[(row["gsm"], row["sample_id"])][kind] = path
    rows = []
    for (gsm, sid), files in sorted(grouped.items()):
        with gzip.open(files["genes"], "rt") as f:
            genes = list(csv.reader(f, delimiter="\t"))
        with gzip.open(files["barcodes"], "rt") as f:
            barcodes = [line.strip() for line in f]
        assert len(set(barcodes)) == len(barcodes)
        with gzip.open(files["matrix"], "rt") as f:
            header = f.readline().strip()
            assert header == "%%MatrixMarket matrix coordinate integer general", header
            for line in f:
                if not line.startswith("%"):
                    n_genes, n_cells, n_entries = map(int, line.split())
                    break
            assert n_genes == len(genes) and n_cells == len(barcodes)
            seen, total = 0, 0
            for line in f:
                gene, cell, count = map(int, line.split())
                assert 1 <= gene <= n_genes and 1 <= cell <= n_cells and count > 0
                seen += 1
                total += count
            assert seen == n_entries
        rows.append({"gsm": gsm, "sample_id": sid, "arm": sid.split(".")[0], "day": 11,
                     "deposited_genes": n_genes, "deposited_cells": n_cells,
                     "count_entries": n_entries, "UMI_total": total,
                     "AT2_cells": "UNESTABLISHED", "fibroblast_cells": "UNESTABLISHED",
                     "coverage_gate": "HOLD_ANNOTATION_REQUIRED"})
    assert len(rows) == 8
    with (HERE / "GSE129605_acquisition_audit.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    result = {"status": "HOLD_ANNOTATION_AND_ANIMAL_PROVENANCE", "samples": len(rows),
              "deposited_cells": sum(r["deposited_cells"] for r in rows),
              "download_bytes": sum(int(r["observed_bytes"]) for r in receipt),
              "checks": ["24 exact GEO file sizes and local SHA256 checksums", "unique barcodes within each sample",
                         "gene/barcode dimensions match each matrix", "all count entries positive integers in valid bounds",
                         "entry count matches MatrixMarket header"],
              "not_done": ["animal provenance beyond deposited identifiers", "post-QC compartment annotation",
                           "AT2/fibroblast coverage gate", "Wnt expression inference"]}
    (HERE / "GSE129605_acquisition_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
