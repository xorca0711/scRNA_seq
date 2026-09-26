"""Verify source integrity, analytical invariants and report links; record outputs."""
from __future__ import annotations

import ast
import hashlib
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

import pandas as pd
import scipy
from scipy.stats import binom

BASE = Path(__file__).resolve().parents[1]


def digest(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def main():
    sources = []
    for provenance in sorted((BASE / "cache/sources").glob("*.provenance.json")):
        record = json.loads(provenance.read_text())
        path = provenance.with_name(provenance.name.removesuffix(".provenance.json"))
        evidence = digest(path)
        assert record["bytes"] == evidence["bytes"] and record["sha256"] == evidence["sha256"], path.name
        record["cached_file"] = path.name
        sources.append(record)
    (BASE / "source_manifest.json").write_text(json.dumps(sources, indent=2) + "\n", encoding="utf-8")
    files = [p for p in BASE.rglob("*") if p.is_file() and not any(
        x in p.relative_to(BASE).parts for x in ("cache", "processed", "__pycache__"))]
    for path in files:
        if path.suffix == ".py":
            ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        elif path.suffix == ".json":
            json.loads(path.read_text(encoding="utf-8"))

    developmental = pd.read_csv(BASE / "tables/d2_negretti_encoded_group_coverage.csv")
    ages = pd.read_csv(BASE / "tables/d2_negretti_age_coverage.csv")
    states = ["AT2", "Transitional", "AT1"]
    assert developmental.epithelial_cells.sum() == ages.epithelial_cells.sum() == 10918
    assert developmental[states].sum().equals(ages[states].sum())
    assert developmental[states].ge(30).all(axis=1).equals(developmental.complete_at_floor)
    assert int(developmental.complete_at_floor.sum()) == 1
    assert len(developmental) == 14 and int(ages.pooled_age_complete_at_floor.sum()) == 3
    repair = pd.read_csv(BASE / "tables/d1_strunz_proposed_window.csv")
    eligible = repair["coverage_pass_in_proposed_window"]
    assert int(eligible.sum()) == 9
    assert repair.loc[eligible, ["AT2", "Krt8+ ADI", "AT1"]].ge(30).all().all()
    assert not repair.loc[repair.in_proposed_day10_15_window, "sample_id"].str.startswith("NC-").any()
    gut = pd.read_csv(BASE / "tables/v1_haber_verified_mouse_coverage.csv")
    assert gut.mouse_id.is_unique and len(gut) == 4 and gut.triplet_eligible.sum() == 2
    skin = pd.read_csv(BASE / "tables/v1_joost_capture_batch_coverage.csv")
    assert len(skin) == 34 and skin.all_cells.sum() == 1422

    planning = pd.read_csv(BASE / "tables/developmental_capture_planning.csv")
    for row in planning.to_dict("records"):
        source = ages.set_index("timepoint").loc[row["timepoint"]]
        probabilities = (source[states].astype(float) / source.epithelial_cells).to_numpy()
        for key, p in [("retained_epithelial_cells_per_independent_unit", probabilities),
                       ("target_if_intermediate_fraction_halves", probabilities * [1, .5, 1])]:
            n = int(row[key])
            assert 1 - binom.cdf(29, n, p).sum() >= .95
            assert 1 - binom.cdf(29, n - 1, p).sum() < .95

    links = 0
    for path in files:
        if path.suffix != ".md":
            continue
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", path.read_text(encoding="utf-8")):
            if target.startswith(("http:", "https:", "mailto:", "#")):
                continue
            target = unquote(target.split("#", 1)[0].strip("<>"))
            assert (path.parent / target).exists(), (path.name, target)
            links += 1
    readiness = json.loads((BASE / "readiness.json").read_text())
    assert not readiness["primary_analysis_ready"] and not readiness["full_pilot_complete"]
    assert not readiness["expression_program_learned"] and not readiness["transfer_test_performed"]
    report = (BASE / "reports/PILOT_REPORT.md").read_text(encoding="utf-8")
    assert not re.search(r"\$[a-zA-Z_]", report)
    review = json.loads((BASE / "figure_review.json").read_text())
    for entry in review["figures"]:
        assert entry["sha256"] == digest(BASE / entry["path"])["sha256"]
        assert entry["status"] == "visually_checked"

    validation = {"stage": "extended_feasibility", "result": "PASS", "python_syntax": True,
        "json_parse": True, "source_hashes_verified": len(sources), "local_links_verified": links,
        "developmental_annotation_totals_reconciled": 10918, "complete_developmental_encoded_groups": 1,
        "complete_repair_triplets_days10_15": 9, "complete_verified_intestinal_mice": 2,
        "skin_batch_totals_reconciled": 1422, "capture_targets_minimal_under_stated_bound": True,
        "normal_control_exclusion_verified": True, "figures_visually_checked": len(review["figures"]),
        "downstream_status_honest_and_blocked": True,
        "initial_local_audit": "Previously validated; unchanged local source and local audit outputs were reused."}
    (BASE / "validation.json").write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    artifacts = {p.relative_to(BASE).as_posix(): digest(p) for p in sorted(files)
                 if p.name != "execution_record.json"}
    execution = {"completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P0_extended_feasibility", "full_pilot_complete": False,
        "status": readiness["current_status"], "python_executable": sys.executable,
        "python_version": sys.version.split()[0], "pandas_version": pd.__version__,
        "scipy_version": scipy.__version__, "numerical_packages_source": ".venv-x64/Lib/site-packages",
        "validation": "PASS", "program_scoring_performed": False, "artifacts": artifacts}
    (BASE / "execution_record.json").write_text(json.dumps(execution, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"validation": "PASS", "sources": len(sources), "artifacts": len(artifacts),
                      "primary_analysis_ready": False}))


if __name__ == "__main__":
    main()
