"""Extract compact design metadata from public GEO family SOFT files.

Keep platform probe tables and contact details out of the review artifact.
Download input to raw_data/dataset_gates/{accession}.soft.txt, then run this.
"""
import collections
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FIELDS = {
    "Series": {"title", "summary", "overall_design", "relation", "supplementary_file", "sample_id"},
    "Sample": {"title", "source_name_ch1", "organism_ch1", "characteristics_ch1",
               "library_strategy", "library_source", "library_selection", "data_processing",
               "supplementary_file", "relation"},
}


def extract(accession):
    source = ROOT / "raw_data/dataset_gates" / f"{accession}.soft.txt"
    series, samples = {}, []
    current, kind = None, None
    digest = hashlib.sha256()
    with source.open("rb") as stream:
        for raw in stream:
            digest.update(raw)
            line = raw.decode("utf-8", errors="replace").strip()
            if line.startswith("^"):
                key, value = line[1:].split(" = ", 1)
                kind = key.title()
                if kind == "Series":
                    current = series
                    current["accession"] = value
                elif kind == "Sample":
                    current = {"accession": value}
                    samples.append(current)
                else:
                    current = None
            elif current is not None and line.startswith(f"!{kind}_") and " = " in line:
                key, value = line.split(" = ", 1)
                field = key[len(kind) + 2:]
                if field in FIELDS[kind]:
                    current.setdefault(field, []).append(value)
    processing = {}
    for sample in samples:
        steps = sample.pop("data_processing", [])
        key = hashlib.sha256(json.dumps(steps, ensure_ascii=False).encode()).hexdigest()[:16]
        processing[key] = steps
        sample["data_processing_id"] = key
    result = {"schema_version": 1, "retrieved": "2026-09-22",
              "source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=all&form=text&view=full",
              "raw_sha256": digest.hexdigest(), "raw_path": str(source.relative_to(ROOT)).replace("\\", "/"),
              "series": series, "samples": samples, "data_processing_definitions": processing,
              "scope": "Selected deposited metadata; sample accession counts are not biological replicate counts."}
    destination = ROOT / "docs/remediation/2026-09-22/dataset_metadata" / f"{accession}.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(accession, "samples", len(samples), "strategies", dict(collections.Counter(
        x for s in samples for x in s.get("library_strategy", []))))
    print("series relations:", series.get("relation", []))


if __name__ == "__main__":
    for accession in ("GSE307351", "GSE242510", "GSE307112", "GSE307128"):
        extract(accession)
