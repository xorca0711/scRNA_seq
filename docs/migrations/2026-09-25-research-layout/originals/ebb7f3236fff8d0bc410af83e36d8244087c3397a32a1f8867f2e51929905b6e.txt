"""Read-only consistency check of saved source reproduction artifacts."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

import reproduce


def main():
    base = Path(__file__).resolve().parent
    provenance = json.loads((base / "provenance.json").read_text(encoding="utf-8"))
    specification = json.loads((base / "specification.json").read_text(encoding="utf-8"))
    for item in provenance["input_files"]:
        assert reproduce.sha256(reproduce.ROOT / item["path"]) == item["sha256"], item["path"]
    assert reproduce.sha256(base / "reproduce.py") == provenance["code_sha256"]
    assert reproduce.sha256(base / "specification.json") == provenance["specification_sha256"]
    for path, expected in provenance["output_sha256"].items():
        assert reproduce.sha256(base / path) == expected, path
    cells = pd.read_csv(base / "tables/single_cell_fpkm.csv", index_col=0)
    bulk = pd.read_csv(base / "tables/bulk_control_fpkm.csv", index_col=0)
    metadata = pd.read_csv(base / "tables/sample_metadata.csv")
    assert len(cells) == 47 and reproduce.BULK not in cells.index
    assert list(bulk.index) == [reproduce.BULK] and metadata["is_bulk"].sum() == 1
    assert cells.shape[1] == 35 and np.isfinite(cells).all().all()
    counts = pd.read_csv(base / "tables/detection_by_threshold.csv")
    codetection = pd.read_csv(base / "tables/conditional_codetection.csv")
    for threshold in specification["thresholds"]:
        detected = reproduce.detection(cells, threshold["cutoff"])
        rows = counts[counts["threshold"] == threshold["label"]].set_index("gene")
        assert rows["n_detected"].equals(detected.sum().reindex(rows.index).rename("n_detected"))
        np.testing.assert_allclose(rows["fraction"], rows["n_detected"] / 47)
        row = codetection[(codetection["threshold"] == threshold["label"]) &
                          (codetection["estimand"] == "Pdgfra_given_Wnt5a")].iloc[0]
        joint, denominator, fraction = reproduce.conditional_fraction(detected["Pdgfra"], detected["Wnt5a"])
        assert row["joint_n"] == joint and row["denominator_n"] == denominator
        assert np.isclose(row["fraction"], fraction)
    print("PASS: input/code/output hashes; 47-cell cohort; bulk exclusion; 35-gene matrix; all threshold counts; conditional denominator")


if __name__ == "__main__":
    main()
