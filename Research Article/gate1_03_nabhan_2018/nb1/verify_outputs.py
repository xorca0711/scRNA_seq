"""Verify tracked Nb1 evidence without scientific packages or raw matrices."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
from analysis.lib.repository_paths import recorded_file, resolve_repo_path


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(2**20),b""):
            h.update(block)
    return h.hexdigest()


def rows(name):
    with (HERE/"tables"/name).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--with-inputs",action="store_true",help="Also require/hash all source inputs, including the local H5AD")
    args=parser.parse_args()
    record=json.loads((HERE/"run_record.json").read_text(encoding="utf-8"))
    assert record["status"]=="descriptive_only"
    for item in record["outputs"]:
        assert digest(HERE/item["path"])==item["sha256"], f"Output changed: {item['path']}"
    for item in record["inputs"]:
        p=resolve_repo_path(ROOT,item["path"])
        if args.with_inputs or p.name in ["run_nb1.py","PROTOCOL.md","modules.json"]:
            assert digest(recorded_file(ROOT,item['path'],item['sha256']))==item["sha256"], f"Input changed: {item['path']}"
    coverage=rows("animal_compartment_coverage.csv")
    assert len(coverage)==record["dimensions"]["units"]==175
    assert sum(int(r["cells"]) for r in coverage)==record["dimensions"]["selected_cells"]==34503
    assert len({r["sample_id"] for r in coverage})==record["dimensions"]["animals"]==25
    for day in [0,11]:
        assert sum(r["cell_type"]=="AT2" and float(r["day"])==day and r["eligible_50"]=="True" for r in coverage)==1
    paired=rows("paired_fibroblast_ligands.csv")
    for endpoint in ["cpm","expected_detection_1000"]:
        for gene,direction in [("Wnt2","AF1>AF2"),("Wnt4","AF2>AF1")]:
            subset=[r for r in paired if r["endpoint"]==endpoint and r["gene"]==gene]
            assert len(subset)==21 and all(r["direction"]==direction for r in subset)
    expr=rows("gene_expression_by_animal.csv")
    for r in expr:
        assert r["feature_present"]=="True"
        assert 0 <= int(r["positive_cells"]) <= int(r["cells"])
        assert abs(float(r["cpm"])-int(r["raw_count"])/int(r["library_umis"])*1e6)<1e-7
    print(f"Nb1 evidence verified: {len(record['outputs'])} output hashes, code/protocol/panel hashes, units, reported paired directions and expression denominators; raw inputs {'checked' if args.with_inputs else 'not requested'}.")


if __name__=="__main__":
    main()
