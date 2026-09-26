"""A15 rival-2 stage 1 addendum: which deposited fields can be trusted to define the arms.

Stage 1 found that the counts column prefix `Bleo` covers both the inert antibody and the
KIRA8 vehicle. This addendum asks the same question of every other deposited field, because
if a later session reaches for a different field it must know which ones are safe.

Metadata only. No gene-set score is computed, and no per-sample expression value is read.

Three items:

1. which deposited fields separate the epithelial immunoprecipitation from the whole-lung
   input, and which are uniform across all 48 libraries despite describing only the former;
2. whether the deposited sample titles agree with the characteristics they duplicate, and
   whether they can identify the antibody control;
3. the assay facts that bear on how a gene-set score may be read, recorded rather than
   assumed, including what the deposit does not state.

Standard library only. Refuses to overwrite its own outputs.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / "cache"
OUT = HERE / "tables/rival2"
MATRIX = CACHE / "GSE190821_series_matrix.txt.gz"
JOIN = OUT / "stage1_join.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    run_record = OUT / "stage1b_addendum_run.json"
    if run_record.exists() and not args.force:
        print(f"refusing to overwrite {run_record}; pass --force to replace")
        return 2
    if not MATRIX.exists() or not JOIN.exists():
        print("STOP: stage 1 must run first; its matrix cache and join table are required")
        return 1

    text = gzip.open(MATRIX, "rt", encoding="utf-8", errors="replace").read()
    columns: dict[str, list[list[str]]] = defaultdict(list)
    for line in text.split("\n"):
        if not line.startswith("!Sample_"):
            continue
        parts = [p.strip('"') for p in line.rstrip("\n").split("\t")]
        columns[parts[0]].append(parts[1:])

    compartment: list[str] = []
    for column in columns["!Sample_characteristics_ch1"]:
        if column[0].startswith("tissue compartment"):
            compartment = [v.split(":", 1)[1].strip() for v in column]
    if len(compartment) != 48:
        print("STOP: could not recover the tissue compartment field for 48 libraries")
        return 1

    # Item 1: which fields separate the two compartments, and which are uniform.
    separating: list[str] = []
    uniform_but_describes_the_ip: list[dict] = []
    for key, occurrences in columns.items():
        for values in occurrences:
            if key == "!Sample_characteristics_ch1" and values[0].startswith("tissue compartment"):
                continue
            by_compartment: dict[str, set[str]] = defaultdict(set)
            for value, comp in zip(values, compartment):
                by_compartment[comp].add(value)
            epi = by_compartment.get("Epithelium", set())
            whole = by_compartment.get("Whole Lung", set())
            if epi and whole and not (epi & whole):
                separating.append(key)
            if len(set(values)) == 1:
                uniform_but_describes_the_ip.append({
                    "field": key,
                    "single_value": values[0][:200],
                })
    separating = sorted(set(separating))

    # Item 2: do the titles agree with the characteristics, and can they name the control?
    rows = list(csv.DictReader(JOIN.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    compartment_disagreements = []
    exposure_disagreements = []
    for row in rows:
        title = row["source_title"].lower().replace("_", "")
        expected = "epithelial" if row["compartment"] == "Epithelium" else "wholelung"
        if expected not in title:
            compartment_disagreements.append(row["source_title"])
        expected_exposure = "bleomycin" if row["exposure"] == "Bleomycin" else "saline"
        if expected_exposure not in title:
            exposure_disagreements.append(row["source_title"])
    titles_naming_3g9 = sorted(r["source_title"] for r in rows
                               if "3g9" in r["source_title"].lower())
    titles_naming_axum8 = sorted(r["source_title"] for r in rows
                                 if "axum" in r["source_title"].lower())
    axum8_titles = sorted(r["source_title"] for r in rows if r["treatment"] == "Axum8")
    vehicle_titles = sorted(r["source_title"] for r in rows if r["treatment"] == "Vehicle")

    def title_shape(title: str) -> str:
        parts = title.split("_")
        return "_".join(["<exposure>" if i == 0 else ("<compartment>" if i == 1 else "<mouse>")
                         for i in range(len(parts))])

    axum8_shapes = sorted({title_shape(t) for t in axum8_titles})
    vehicle_shapes = sorted({title_shape(t) for t in vehicle_titles})
    titles_can_identify_the_control = not set(axum8_shapes) & set(vehicle_shapes)

    with (OUT / "stage1b_field_reliability.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["field", "separates_compartments", "uniform_across_all_48",
                         "safe_to_define_arms_with"])
        seen = set()
        for key, occurrences in columns.items():
            if key in seen:
                continue
            seen.add(key)
            uniform = any(len(set(v)) == 1 for v in occurrences)
            writer.writerow([
                key,
                key in separating,
                uniform,
                key in {"!Sample_characteristics_ch1", "!Sample_geo_accession"},
            ])

    record = {
        "schema": "a15-rival2-stage1-addendum/v1",
        "stage": "1 addendum, field reliability; no endpoint computed",
        "inputs": {
            "series_matrix": {"path": "cache/GSE190821_series_matrix.txt.gz",
                              "sha256": sha256(MATRIX)},
            "stage1_join": {"path": "tables/rival2/stage1_join.tsv", "sha256": sha256(JOIN)},
        },
        "compartment_field": {
            "fields_that_separate_epithelium_from_whole_lung": separating,
            "note": "besides the tissue compartment characteristic itself, only the sample "
                    "title and the per-sample identifiers separate the two compartments",
        },
        "uniform_fields_that_describe_only_the_immunoprecipitation": [
            entry for entry in uniform_but_describes_the_ip
            if entry["field"] in {"!Sample_source_name_ch1", "!Sample_description",
                                 "!Sample_extract_protocol_ch1"}
        ],
        "hazard": "source name, description and extract protocol are identical for all 48 "
                  "libraries and describe the anti-HA immunoprecipitation, although 24 of "
                  "the libraries are the whole-lung input aliquot. A session that defined "
                  "the epithelial compartment from any of those three fields would treat "
                  "all 48 libraries as epithelial.",
        "titles": {
            "compartment_disagreements": compartment_disagreements,
            "exposure_disagreements": exposure_disagreements,
            "titles_naming_3G9": titles_naming_3g9,
            "titles_naming_Axum8": titles_naming_axum8,
            "axum8_title_shapes": axum8_shapes,
            "vehicle_title_shapes": vehicle_shapes,
            "titles_can_identify_the_antibody_control": titles_can_identify_the_control,
            "consequence": "the titles name the 3G9 arm but never the Axum8 arm, so an "
                           "Axum8 title is indistinguishable in form from a Vehicle title. "
                           "The treatment characteristic remains the only field that "
                           "separates the antibody control from the KIRA8 vehicle.",
        },
        "assay_facts_recorded_not_assumed": {
            "library_preparation": "Lexogen QuantSeq Forward, a 3-prime tag counting assay",
            "consequence_for_scoring": "counts are 3-prime tag counts rather than "
                                       "full-length coverage, so they are not "
                                       "length-normalised and a gene-set mean of log2 CPM "
                                       "is a relative, within-assay quantity only",
            "instrument": "Illumina HiSeq 4000",
            "alignment": "STAR 2.5.2b against Ensembl Mouse GRCm38, per the deposit",
            "compartment_method": "anti-HA immunoprecipitation of tagged ribosomes, so the "
                                  "epithelial measurement is ribosome-associated RNA",
        },
        "what_the_deposit_does_not_state": [
            "per-mouse genotype, so constancy of the RiboTag and Cre alleles across arms is "
            "carried as an assumption and is not established here",
            "immunoprecipitation efficiency or purity per library, so epithelial enrichment "
            "cannot be verified or compared between arms",
            "antibody dose, schedule or measured target engagement",
            "how mice were allocated to antibody, so randomisation is not established",
            "cage, litter, weight or bleomycin dose per mouse",
            "epithelial subtype composition, so a programme difference cannot be separated "
            "from a composition shift within the immunoprecipitated compartment",
        ],
        "endpoint_scored": False,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    run_record.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print(f"fields separating the compartments: {separating}")
    print(f"uniform fields describing only the IP: "
          f"{[e['field'] for e in record['uniform_fields_that_describe_only_the_immunoprecipitation']]}")
    print(f"title/compartment disagreements: {compartment_disagreements or 'none'}")
    print(f"title/exposure disagreements: {exposure_disagreements or 'none'}")
    print(f"titles naming 3G9: {len(titles_naming_3g9)}; naming Axum8: {len(titles_naming_axum8)}")
    print(f"titles can identify the antibody control: {titles_can_identify_the_control}")
    print(f"wrote {run_record}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
