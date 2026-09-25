"""Retrieve bounded public GEO metadata; never interpret GSMs as biological n.

No counts, platform tables, contact details, or private annotations are exported.
The complete downloaded response remains in the ignored paper-local cache.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import urllib.request
from xml.etree import ElementTree as ET

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import archive_existing_record, code_identity, write_json_atomic
ALLOWED = {
    "Series": {"title", "summary", "overall_design", "relation", "supplementary_file", "sample_id", "pubmed_id"},
    "Sample": {"title", "source_name_ch1", "organism_ch1", "characteristics_ch1", "library_strategy", "library_source", "library_selection", "data_processing", "supplementary_file", "relation"},
}


def parse_miniml(raw: bytes) -> dict:
    tree = ET.fromstring(raw)
    ns = {"g": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}
    series, samples = {}, []
    mapping = {"Title": "title", "Summary": "summary", "Overall-Design": "overall_design", "Pubmed-ID": "pubmed_id", "Source": "source_name_ch1", "Organism": "organism_ch1", "Characteristics": "characteristics_ch1", "Library-Strategy": "library_strategy", "Library-Source": "library_source", "Library-Selection": "library_selection", "Data-Processing": "data_processing", "Supplementary-Data": "supplementary_file"}
    for kind in ("Series", "Sample"):
        for node in tree.findall("g:" + kind, ns):
            current = {"accession": node.attrib["iid"]}
            for child in node.iter():
                tag = child.tag.split("}")[-1]
                field = mapping.get(tag)
                if field in ALLOWED[kind]:
                    value = (child.text or "").strip()
                    if tag == "Characteristics":
                        value = child.attrib.get("tag", "unspecified") + ": " + value
                    current.setdefault(field, []).append(value)
                elif tag == "Sample-Ref":
                    current.setdefault("sample_id", []).append(child.attrib["ref"])
                elif tag == "Relation":
                    current.setdefault("relation", []).append(child.attrib["type"] + ": " + child.attrib["target"])
            if kind == "Series":
                if series:
                    raise ValueError("Unexpected multiple series")
                series = current
            else:
                samples.append(current)
    if not series or not samples:
        raise ValueError("Response does not contain GEO series and sample records")
    declared = set(series.get("sample_id", []))
    observed = [s["accession"] for s in samples]
    if len(observed) != len(set(observed)) or declared != set(observed):
        raise ValueError("Duplicate or incomplete sample inventory")
    return {"series": series, "samples": samples}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions", nargs="+")
    parser.add_argument("--offline", action="store_true", help="Use cached metadata only")
    args = parser.parse_args()
    cache = PAPER / "cache" / "geo"
    dest = PAPER / "trials" / "u0_geo_design_audit"
    cache.mkdir(parents=True, exist_ok=True)
    dest.mkdir(parents=True, exist_ok=True)
    audit = []
    for acc in args.accessions:
        if not (acc.startswith("GSE") and acc[3:].isdigit()):
            raise ValueError(f"Invalid GEO series: {acc}")
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=xml&view=full"
        local = cache / f"{acc}.xml"
        provenance = cache / f"{acc}.retrieval.json"
        if not args.offline:
            with urllib.request.urlopen(url, timeout=60) as response:
                raw = response.read(20_000_001)
            if len(raw) > 20_000_000:
                raise ValueError("Metadata response exceeds 20 MB cap")
            parsed = parse_miniml(raw)
            if parsed["series"]["accession"] != acc:
                raise ValueError("Returned series differs from requested accession")
            local.write_bytes(raw)
            provenance.write_text(json.dumps({"source_url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "sha256": hashlib.sha256(raw).hexdigest()}, indent=2) + "\n", encoding="utf-8")
        raw = local.read_bytes()
        meta = json.loads(provenance.read_text(encoding="utf-8"))
        if hashlib.sha256(raw).hexdigest() != meta["sha256"]:
            raise ValueError("Cached response hash mismatch")
        parsed = parse_miniml(raw)
        result = {"schema_version": 1, **meta, **parsed, "scope": "Metadata only; sample records are not independent biological units."}
        (dest / f"{acc}.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        summary = {"accession": acc, "sample_records": len(parsed["samples"]), "title": parsed["series"].get("title", []), "assays": dict(Counter(x for s in parsed["samples"] for x in s.get("library_strategy", []))), "source_sha256": meta["sha256"]}
        audit.append(summary)
        print(json.dumps(summary, ensure_ascii=False))
    record = {"run_utc": datetime.now(timezone.utc).isoformat(), "code": code_identity(ROOT, __file__), "mode": "offline" if args.offline else "retrieval", "accessions": audit, "expression_analysis_run": False}
    archive_existing_record(dest / "run_record.json")
    write_json_atomic(dest / "run_record.json", record)


if __name__ == "__main__":
    main()
