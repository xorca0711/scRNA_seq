"""Retrieve public GEO design metadata only; never download sequencing data.

Example: python RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/01_audit_geo_metadata.py GSE141635
Outputs are accession-level JSON extracts in the A1 metadata directory.
Full responses are cached under ignored tmp/a1_epigenetic_metadata/.
Existing extracts are refused unless --refresh is explicitly supplied.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
import urllib.request

STUDY = Path(__file__).resolve().parents[1]
ROOT = STUDY.parents[1]
SERIES_KEYS = {"title", "geo_accession", "status", "last_update_date", "pubmed_id", "overall_design",
               "type", "sample_id", "supplementary_file", "relation"}
SAMPLE_KEYS = {"title", "geo_accession", "source_name_ch1", "organism_ch1",
               "characteristics_ch1", "molecule_ch1", "library_strategy",
               "library_source", "library_selection", "platform_id", "relation",
               "data_processing", "extract_protocol_ch1", "treatment_protocol_ch1",
               "growth_protocol_ch1", "description"}


def extract(raw: bytes, accession: str, url: str) -> dict:
    series, samples, current = {}, [], None
    for line in raw.decode("utf-8-sig").splitlines():
        if line.startswith("^SERIES = "):
            if line.split(" = ", 1)[1] != accession:
                raise ValueError("Unexpected series in response")
            current = series
        elif line.startswith("^SAMPLE = "):
            current = {"accession": line.split(" = ", 1)[1]}
            samples.append(current)
        elif line.startswith("^"):
            current = None
        elif current is not None and " = " in line:
            key, value = line.split(" = ", 1)
            prefix = "!Series_" if current is series else "!Sample_"
            if not key.startswith(prefix):
                continue
            key = key[len(prefix):]
            keep = SERIES_KEYS if current is series else SAMPLE_KEYS
            if key in keep or key.startswith("supplementary_file"):
                current.setdefault(key, []).append(value)
    expected = set(series.get("sample_id", []))
    observed = [s["accession"] for s in samples]
    if series.get("geo_accession") != [accession] or not expected:
        raise ValueError("No valid public GEO series metadata")
    if expected != set(observed) or len(observed) != len(expected):
        raise ValueError("Missing or duplicate sample records")
    return {
        "accession": accession,
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_url": url,
        "response_sha256": hashlib.sha256(raw).hexdigest(),
        "scope": "GEO catalog metadata; file payloads, biological independence and cross-assay pairing not verified",
        "series": series, "sample_count": len(samples), "samples": samples,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("accessions", nargs="+")
    parser.add_argument("--refresh", action="store_true")
    args = parser.parse_args()
    destination = STUDY / "metadata"
    cache = ROOT / "tmp/a1_epigenetic_metadata"
    destination.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    for accession in args.accessions:
        if not re.fullmatch(r"GSE\d+", accession):
            raise ValueError(f"Invalid accession: {accession}")
        output = destination / f"{accession}.json"
        if output.exists() and not args.refresh:
            print(f"{accession}: retained existing extract", flush=True)
            continue
        url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=all&form=text&view=full"
        request = urllib.request.Request(url, headers={"User-Agent": "scRNA_seq-public-metadata-audit/1.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            raw = response.read()
        record = extract(raw, accession, url)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S")
        (cache / f"{accession}_{stamp}.soft.txt").write_bytes(raw)
        if output.exists():
            archive = destination / "history" / stamp
            archive.mkdir(parents=True, exist_ok=True)
            output.replace(archive / output.name)
        output.write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        print(f"{accession}: {record['sample_count']} GSM records; {record['series'].get('type', [])}", flush=True)
        time.sleep(0.4)


if __name__ == "__main__":
    main()
