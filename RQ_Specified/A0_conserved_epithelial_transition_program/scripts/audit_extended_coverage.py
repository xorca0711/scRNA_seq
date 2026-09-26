"""Close the metadata audit and estimate illustrative capture requirements.

No expression values are read. Library suffixes are provisional provenance
groups, never inferred animal IDs. The cellxgene decoder implements only the
JSON column subset of the cached official 1.1.1 FlatBuffers schema.
"""
from __future__ import annotations

import hashlib
import json
import struct
from pathlib import Path

import pandas as pd
from scipy.stats import binom

from audit_public_metadata import read_soft, characteristics

BASE = Path(__file__).resolve().parents[1]
CACHE = BASE / "cache/sources"
OUT = BASE / "tables"
FLOOR = 30
STATES = ["AT2", "Transitional", "AT1"]
AGES = ["E12", "E15", "E16", "E18", "P0", "P3", "P5", "P7", "P14"]


def decode_annotations(payload: bytes) -> pd.DataFrame:
    def uint(pos):
        return struct.unpack_from("<I", payload, pos)[0]

    def ushort(pos):
        return struct.unpack_from("<H", payload, pos)[0]

    def field(table, index):
        vtable = table - struct.unpack_from("<i", payload, table)[0]
        offset = ushort(vtable + 4 + 2 * index) if 4 + 2 * index < ushort(vtable) else 0
        if not offset:
            raise ValueError("Required FlatBuffers field absent")
        return table + offset

    def pointer(pos):
        return pos + uint(pos)

    def array(table):
        start = pointer(field(table, 0))
        length = uint(start)
        value = json.loads(payload[start + 4:start + 4 + length].decode("utf-8"))
        if not isinstance(value, list):
            raise ValueError("JSON array expected")
        return value

    root = uint(0)
    rows, columns = uint(field(root, 0)), uint(field(root, 1))
    if payload[field(root, 3)] != 5:
        raise ValueError("Column index must be a JSONEncodedArray")
    names = array(pointer(field(root, 4)))
    vector = pointer(field(root, 2))
    assert uint(vector) == columns == len(names)
    data = {}
    for i, name in enumerate(names):
        column = pointer(vector + 4 + 4 * i)
        if payload[field(column, 0)] != 5:
            raise ValueError("This metadata decoder accepts only JSONEncodedArray columns")
        values = array(pointer(field(column, 1)))
        assert len(values) == rows
        data[name] = values
    return pd.DataFrame(data)


def capture_target(probabilities, required=FLOOR, success=0.95):
    """Smallest N with union-bound P(all three counts >= required) >= success."""
    probabilities = list(probabilities)
    if min(probabilities) <= 0:
        return None
    lower_bound = lambda n: max(0.0, 1.0 - float(binom.cdf(required - 1, n, probabilities).sum()))
    lo, hi = 3 * required - 1, 3 * required
    while lower_bound(hi) < success:
        hi *= 2
        if hi > 10_000_000:
            return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if lower_bound(mid) >= success:
            hi = mid
        else:
            lo = mid
    assert lower_bound(hi) >= success and lower_bound(hi - 1) < success
    return hi


def main():
    schema = json.loads((CACHE / "Negretti_epi_schema.json").read_text())["schema"]
    cells = decode_annotations((CACHE / "Negretti_obs.bin").read_bytes())
    assert list(cells.columns) == ["name_0", "celltype", "timepoint"]
    assert len(cells) == schema["dataframe"]["nObs"] == 10918
    assert cells.name_0.is_unique and cells.notna().all().all()
    assert set(cells.timepoint) == set(AGES)
    cells["barcode_suffix"] = cells.name_0.str.split("_", n=1).str[1]
    assert cells.barcode_suffix.notna().all()
    cells["encoded_library_group"] = cells.timepoint + ":" + cells.barcode_suffix
    groups = cells.groupby(["timepoint", "barcode_suffix", "encoded_library_group"], sort=False)
    lib = groups.celltype.value_counts().unstack(fill_value=0).reset_index()
    lib["epithelial_cells"] = lib[list(cells.celltype.unique())].sum(axis=1)
    lib["complete_at_floor"] = lib[STATES].ge(FLOOR).all(axis=1)
    lib["unit_identity_status"] = "barcode_suffix_group; no public mouse/hash identity"
    lib.to_csv(OUT / "d2_negretti_encoded_group_coverage.csv", index=False)
    age = pd.crosstab(cells.timepoint, cells.celltype).reindex(AGES)
    age["epithelial_cells"] = age.sum(axis=1)
    age["pooled_age_complete_at_floor"] = age[STATES].ge(FLOOR).all(axis=1)
    age["eligible_independent_units"] = "not_established"
    age.to_csv(OUT / "d2_negretti_age_coverage.csv")
    planning = []
    for time in [x for x in AGES if x.startswith("P")]:
        row = age.loc[time]
        p = (row[STATES].astype(float) / float(row.epithelial_cells)).tolist()
        p_half = [p[0], p[1] / 2, p[2]]
        planning.append({"timepoint": time, "observed_epithelial_cells": int(row.epithelial_cells),
            "observed_intermediate_cells": int(row.Transitional), "intermediate_fraction": p[1],
            "retained_epithelial_cells_per_independent_unit": capture_target(p),
            "target_if_intermediate_fraction_halves": capture_target(p_half),
            "assumed_triplet_capture_probability_lower_bound": 0.95,
            "kind": "illustrative_capture_planning_not_gene_effect_power"})
    pd.DataFrame(planning).to_csv(OUT / "developmental_capture_planning.csv", index=False)

    _, skin_samples = read_soft(CACHE / "GSE67602_family.soft.gz")
    skin = pd.DataFrame([{"gsm": s["gsm"], "cell_id": s["title"][0],
        "capture_batch": s["title"][0].split("_")[0],
        "author_celltype": characteristics(s)["cell type level 1"]} for s in skin_samples])
    assert len(skin) == 1422 and skin.cell_id.is_unique
    sc = pd.crosstab(skin.capture_batch, skin.author_celltype)
    sc["all_cells"] = sc.sum(axis=1)
    sc["mouse_identity_status"] = "unmapped; 34 capture batches are not 19 mice"
    sc.to_csv(OUT / "v1_joost_capture_batch_coverage.csv")
    skin.author_celltype.value_counts().rename_axis("author_celltype").to_csv(
        OUT / "v1_joost_published_state_counts.csv", header=["cells"])

    repair = pd.read_csv(OUT / "d1_strunz_high_resolution_coverage.csv")
    repair["day"] = repair.time_point.str.extract(r"(\d+)")[0].astype(int)
    repair["source_label_normal_control"] = repair.sample_id.str.startswith("NC-")
    repair["in_proposed_day10_15_window"] = repair.day.between(10, 15) & ~repair.source_label_normal_control
    repair["coverage_pass_in_proposed_window"] = repair.triplet_eligible & repair.in_proposed_day10_15_window
    repair.to_csv(OUT / "d1_strunz_proposed_window.csv", index=False)

    record = {"stage": "P0_extended_feasibility", "expression_values_read": False,
        "cell_floor": FLOOR, "minimum_independent_units": 3,
        "negretti": {"cells": len(cells), "encoded_library_groups": len(lib),
            "complete_encoded_groups": int(lib.complete_at_floor.sum()),
            "complete_pooled_ages": int(age.pooled_age_complete_at_floor.sum()),
            "biological_replicates_recovered": False,
            "complete_group_ids": lib.loc[lib.complete_at_floor, "encoded_library_group"].tolist(),
            "export_expression_type": "SCT data slot according to author export code; raw counts not obtained"},
        "strunz": {"available_libraries": len(repair),
            "complete_triplets_all_times": int(repair.triplet_eligible.sum()),
            "complete_triplets_proposed_day10_15": int(repair.coverage_pass_in_proposed_window.sum()),
            "paper_high_resolution_mice": 36,
            "available_annotation_mice_assumption": "one available enriched library per source mouse; paper reports two mice per sampled time",
            "version_difference": "paper 34575 cells/36 mice; cached annotation 32559 cells/32 libraries; excluded material not reconstructed",
            "window_status": "candidate design based on repair interval and annotation coverage; no P1 freeze"},
        "joost": {"cells": len(skin), "capture_batches": len(sc), "paper_mice": 19,
            "max_cells_per_capture_batch": int(sc.all_cells.max()),
            "mouse_mapping_recovered": False},
        "sources": {n: hashlib.sha256((CACHE / n).read_bytes()).hexdigest() for n in [
            "Negretti_obs.bin", "Negretti_epi_schema.json", "Negretti_cellxgene.rmd",
            "GSE67602_family.soft.gz", "Joost_fulltext.xml", "Strunz_fulltext.xml",
            "Haber_Supplementary_Table1.xlsx"]}}
    (BASE / "extended_audit_record.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"negretti_cells": len(cells), "negretti_complete_encoded_groups": int(lib.complete_at_floor.sum()),
        "repair_complete_triplets_day10_15": int(repair.coverage_pass_in_proposed_window.sum()),
        "skin_batches": len(sc), "primary_analysis_ready": False}))


if __name__ == "__main__":
    main()
