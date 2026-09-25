"""Generate reviewable design tables from retrieved metadata, not expression."""
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
import unicodedata

BASE = Path(__file__).resolve().parent / "u0_geo_design_audit"


def write_csv(name, rows):
    with (BASE / name).open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    inventory, design, patient = [], [], []
    for file in sorted(BASE.glob("GSE*.json")):
        data = json.loads(file.read_text(encoding="utf-8"))
        accession = file.stem
        for s in data["samples"]:
            title = " | ".join(s.get("title", []))
            fields = dict(x.split(": ", 1) for x in s.get("characteristics_ch1", []) if ": " in x)
            pid = re.search(r"patient (\d+)", title)
            assay = re.search(r"\[([^\]]+)\]", title)
            urls = s.get("supplementary_file", [])
            inventory.append({"series": accession, "gsm": s["accession"], "title": title,
                              "assay_label": assay.group(1) if assay else "see_processing_metadata",
                              "patient_label": pid.group(1) if pid else "unresolved",
                              "organism": " | ".join(s.get("organism_ch1", [])),
                              "treatment": unicodedata.normalize("NFKC", fields.get("treatment", "")),
                              "endpoint": fields.get("treatment endpoint", ""),
                              "age": fields.get("age", ""), "file_count": len(urls),
                              "unit_status": "deposited_label_only; independence_not_assumed"})
            if accession == "GSE300288":
                proc = " ".join(s.get("data_processing", []))
                design.append({"gsm": s["accession"], "treatment": unicodedata.normalize("NFKC", fields["treatment"]),
                               "endpoint": fields["treatment endpoint"], "age": fields["age"],
                               "replicate_ordinal": re.search(r"replicate (\d+)", title).group(1),
                               "assay_metadata_conflict": "Formalin-fixed" in fields.get("cell line", "") and "Cell Ranger" in proc,
                               "animal_id": "unresolved", "pool_size": "unresolved",
                               "file_triple_present": all(any(u.endswith(suffix) for u in urls) for suffix in ("_matrix.mtx.gz", "_features.tsv.gz", "_barcodes.tsv.gz"))})
            if accession in ("GSE308103", "GSE307534", "GSE307529"):
                if pid is None:
                    raise ValueError(f"No patient label: {title}")
                states = re.findall(r"\b(AAH|AIS|MIA|LUAD)\b", title)
                state = states[0] if states else ("Normal" if "normal" in title.lower() else "unresolved")
                patient.append({"series": accession, "patient_label": pid.group(1), "histology_from_title": state, "gsm": s["accession"], "repeated_tissue": "second" in title.lower()})
    write_csv("sample_inventory.csv", inventory)
    write_csv("GSE300288_design.csv", design)
    write_csv("human_sample_links.csv", patient)
    counts = Counter((s["endpoint"], s["treatment"], s["age"]) for s in design)
    write_csv("GSE300288_group_counts.csv", [{"endpoint": k[0], "treatment": k[1], "age": k[2], "deposited_libraries": v} for k, v in sorted(counts.items())])
    coverage = defaultdict(lambda: defaultdict(set))
    for row in patient:
        coverage[row["series"]][row["patient_label"]].add(row["histology_from_title"])
    overview = []
    for accession, people in sorted(coverage.items()):
        overview.append({"series": accession, "patient_labels": len(people),
                         "patients_with_precursor_and_luad": sum(bool(states & {"AAH", "AIS", "MIA"}) and "LUAD" in states for states in people.values()),
                         "patients_with_normal_precursor_luad": sum(bool(states & {"AAH", "AIS", "MIA"}) and {"Normal", "LUAD"}.issubset(states) for states in people.values()),
                         "interpretation": "title-derived candidate pairing; verify specimen crosswalk and lesion histology"})
    write_csv("human_pairing_coverage.csv", overview)
    print(json.dumps({"series": len({x['series'] for x in inventory}), "sample_records": len(inventory), "treatment_groups": [dict(endpoint=k[0], treatment=k[1], libraries=v) for k,v in sorted(counts.items())], "human_pairing": overview}, ensure_ascii=False))


if __name__ == "__main__":
    main()
