"""Audit public GEO sample design and author labels; do not score expression."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

BASE = Path(__file__).resolve().parents[1]
CACHE = BASE / "cache/sources"
OUT = BASE / "tables"
FLOOR = 30


def read_soft(path: Path) -> tuple[dict, list[dict]]:
    series: dict[str, list[str]] = {}
    samples = []
    current = None
    for line in gzip.decompress(path.read_bytes()).decode("utf-8").splitlines():
        if line.startswith("^SAMPLE = "):
            if current is not None:
                samples.append(current)
            current = {"gsm": line.split(" = ", 1)[1]}
        elif line.startswith("!") and " = " in line:
            key, value = line[1:].split(" = ", 1)
            if key.startswith("Series_"):
                series.setdefault(key[7:], []).append(value)
            elif current is not None and key.startswith("Sample_"):
                current.setdefault(key[7:], []).append(value)
    if current is not None:
        samples.append(current)
    return series, samples


def characteristics(sample: dict) -> dict:
    return {s.split(": ", 1)[0].lower(): s.split(": ", 1)[1]
            for s in sample.get("characteristics_ch1", []) if ": " in s}


def coverage(frame: pd.DataFrame, sample: str, state: str, states: list[str],
             sample_fields: list[str], filename: str) -> pd.DataFrame:
    counts = pd.crosstab(frame[sample], frame[state]).reindex(columns=states, fill_value=0)
    facts = frame.groupby(sample)[sample_fields].first()
    if frame.groupby(sample)[sample_fields].nunique(dropna=True).gt(1).any().any():
        raise ValueError("Conflicting sample facts")
    result = facts.join(counts)
    result["start_intermediate_eligible"] = result[states[:2]].ge(FLOOR).all(axis=1)
    result["intermediate_destination_eligible"] = result[states[1:]].ge(FLOOR).all(axis=1)
    result["triplet_eligible"] = result[states].ge(FLOOR).all(axis=1)
    result.reset_index().to_csv(OUT / filename, index=False)
    return result


def main() -> None:
    OUT.mkdir(exist_ok=True)
    extracted = {}
    sample_rows = []
    all_samples = {}
    for path in sorted(CACHE.glob("*_family.soft.gz")):
        acc = path.name.split("_")[0]
        series, samples = read_soft(path)
        all_samples[acc] = samples
        extracted[acc] = {"title": series.get("title"), "overall_design": series.get("overall_design"),
                          "n_geo_samples": len(samples), "pubmed_ids": series.get("pubmed_id", []),
                          "supplementary_files": series.get("supplementary_file", [])}
        for s in samples:
            c = characteristics(s)
            sample_rows.append({"accession": acc, "gsm": s["gsm"], "title": s["title"][0],
                "organism": " | ".join(s.get("organism_ch1", [])), "age": c.get("age", ""),
                "time": c.get("injury timecourse", c.get("time point", "")),
                "experiment": c.get("experiment", ""), "cell_type": c.get("cell type", ""),
                "characteristics_json": json.dumps(c, sort_keys=True),
                "supplementary_files_json": json.dumps(s.get("supplementary_file", []))})
    pd.DataFrame(sample_rows).to_csv(OUT / "geo_sample_metadata.csv", index=False)
    (OUT / "geo_series_metadata.json").write_text(json.dumps(extracted, indent=2) + "\n", encoding="utf-8")

    high = pd.read_csv(CACHE / "GSE141259_HighResolution_cellinfo.csv.gz", sep="\t")
    if not high.cell_barcode.is_unique:
        raise ValueError("Duplicate cell barcodes in high-resolution annotation")
    high_out = coverage(high, "identifier", "cell_type", ["AT2", "Krt8+ ADI", "AT1"],
                        ["sample_id", "time_point"], "d1_strunz_high_resolution_coverage.csv")
    whole = pd.read_csv(CACHE / "GSE141259_WholeLung_cellinfo.csv.gz")
    whole_out = coverage(whole, "orig.ident", "cell.type", ["AT2 cells", "Krt8 ADI", "AT1 cells"],
                         ["grouping"], "d1_strunz_whole_lung_coverage.csv")
    high_summary = high_out.groupby("time_point")[["start_intermediate_eligible",
        "intermediate_destination_eligible", "triplet_eligible"]].sum().reset_index()
    high_summary.to_csv(OUT / "d1_strunz_high_resolution_time_coverage.csv", index=False)

    # A saved matrix header carries author cell labels, not expression values.
    header = next(csv.reader((CACHE / "GSE92332_atlas_header.tsv").open(encoding="utf-8"), delimiter="\t"))
    header = [h for h in header if h]
    if len(set(header)) != len(header):
        raise ValueError("Duplicate intestinal header cell IDs")
    parts = [h.split("_", 2) for h in header]
    if any(len(p) != 3 for p in parts):
        raise ValueError("Unexpected intestinal header format")
    intestine = pd.DataFrame(parts, columns=["batch", "barcode", "author_celltype"])
    counts = pd.crosstab(intestine.batch, intestine.author_celltype)
    counts.reset_index().to_csv(OUT / "v1_haber_batch_state_counts.csv", index=False)
    # GEO explicitly maps four infection-control mice into atlas batches.
    # The source title for mouse 1 has the typo "Atla sbatch 3"; matching the
    # terminal word "batch" recovers that explicit number without guessing it.
    mapping = []
    for s in all_samples["GSE92332"]:
        title = s["title"][0]
        match = re.search(r"Control-Mouse(\d+)\).*?batch (\d+)\]", title)
        if match:
            mapping.append({"gsm": s["gsm"], "batch": "B" + match[2],
                            "mouse_id": "Control-Mouse" + match[1], "source_title": title})
    mapping = pd.DataFrame(mapping)
    if len(mapping) != 4 or not mapping.batch.is_unique or not mapping.mouse_id.is_unique:
        raise ValueError("Four one-to-one source-explicit control mappings expected")
    mapping.to_csv(OUT / "v1_haber_verified_mouse_mapping.csv", index=False)
    linked = intestine.merge(mapping, on="batch", how="inner", validate="many_to_one")
    states = ["Stem", "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal"]
    v1 = coverage(linked, "mouse_id", "author_celltype", states,
                  ["batch", "gsm"], "v1_haber_verified_mouse_coverage.csv")

    provenance = []
    for p in sorted(CACHE.glob("*.provenance.json")):
        entry = json.loads(p.read_text())
        entry["cached_file"] = p.name.removesuffix(".provenance.json")
        entry.setdefault("retrieval_mode", "complete_file")
        provenance.append(entry)
    (BASE / "source_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    result = {"stage": "P0", "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "program_scoring_performed": False,
        "strunz_high_resolution": {"cells": len(high), "samples": len(high_out),
            "triplet_units_all_times": int(high_out.triplet_eligible.sum()),
            "start_intermediate_units_all_times": int(high_out.start_intermediate_eligible.sum()),
            "intermediate_destination_units_all_times": int(high_out.intermediate_destination_eligible.sum()),
            "max_triplet_units_same_time": int(high_summary.triplet_eligible.max())},
        "strunz_whole_lung": {"cells": len(whole), "samples": len(whole_out),
            "AT1_cells_total": int(whole_out["AT1 cells"].sum()),
            "triplet_units": int(whole_out.triplet_eligible.sum())},
        "haber_atlas": {"header_cells": len(intestine), "batches": int(intestine.batch.nunique()),
            "source_explicit_control_mice": len(mapping), "verified_mouse_cells": len(linked),
            "candidate_branch": states, "triplet_units": int(v1.triplet_eligible.sum()),
            "triplet_mouse_ids": list(v1.index[v1.triplet_eligible]),
            "unmapped_batches": sorted(set(intestine.batch) - set(mapping.batch)),
            "transition_evidence_gate": "unresolved; author state labels alone are insufficient"},
        "design_note": "Across-time unit counts are not same-time replication; no cohort selected by this script."}
    (BASE / "public_audit_record.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
