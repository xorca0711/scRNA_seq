#!/usr/bin/env python
"""Validate the lightweight, tracked contract of the portfolio repository.

This check deliberately uses only the Python standard library. It verifies the
published headline numbers against tracked tables, checks machine-readable
artefacts, and catches broken local Markdown links without requiring raw data
or the scanpy environment.
"""

from __future__ import annotations

import csv
import json
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


REPO = Path(__file__).resolve().parents[2]
LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
SKIP_DIRS = {".git", ".venv", ".venv-x64", ".claude", "raw_data", "Thesis"}


class Validation:
    def __init__(self) -> None:
        self.failures: list[str] = []
        self.checks = 0

    def require(self, condition: bool, message: str) -> None:
        self.checks += 1
        if not condition:
            self.failures.append(message)

    def equal(self, actual, expected, label: str) -> None:
        self.require(actual == expected, f"{label}: expected {expected!r}, got {actual!r}")


def read_json(path: str) -> dict:
    return json.loads((REPO / path).read_text(encoding="utf-8"))


def read_csv(path: str) -> list[dict[str, str]]:
    with (REPO / path).open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def median_by(rows: list[dict[str, str]], group: str, value: str) -> dict[float, float]:
    grouped: dict[float, list[float]] = defaultdict(list)
    for row in rows:
        grouped[float(row[group])].append(float(row[value]))
    return {key: statistics.median(values) for key, values in grouped.items()}


def markdown_files() -> list[Path]:
    return [
        path
        for path in REPO.rglob("*.md")
        if not any(part in SKIP_DIRS for part in path.relative_to(REPO).parts)
    ]


def check_markdown_links(result: Validation) -> None:
    for document in markdown_files():
        for match in LINK_RE.finditer(document.read_text(encoding="utf-8")):
            raw_target = match.group(1).strip().strip("<>")
            if raw_target.startswith(("http://", "https://", "mailto:", "#")):
                continue
            target = unquote(raw_target.split("#", 1)[0])
            if not target:
                continue
            # A Markdown title after the URL is not used in this repository,
            # but stripping it makes this check friendly to that syntax.
            target = target.split(' "', 1)[0]
            resolved = (document.parent / target).resolve()
            result.require(
                resolved.exists(),
                f"broken link in {document.relative_to(REPO)}: {raw_target}",
            )


def check_machine_readable_files(result: Validation) -> None:
    for path in (REPO / "analysis").rglob("*.json"):
        try:
            json.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeError, json.JSONDecodeError) as exc:
            result.require(False, f"invalid JSON {path.relative_to(REPO)}: {exc}")
        else:
            result.require(True, str(path))

    inventory = read_json("analysis/raw_data_inventory.json")
    result.equal(inventory["raw_root"], "raw_data/", "portable raw-data root")
    result.equal(inventory["n_files"], 51, "raw-data inventory file count")
    result.require(
        all(item["path"].startswith("raw_data/") for item in inventory["files"]),
        "raw-data inventory contains a non-portable path",
    )
    result.require(
        not (REPO / "analysis/GSE262927/inventory/raw_files.json").exists(),
        "obsolete dataset-specific raw-data inventory still exists",
    )


def check_headline_results(result: Validation) -> None:
    mouse = read_json("analysis/GSE262927/logs/decisions.json")
    human = read_json("analysis/GSE178360/logs/decisions.json")
    result.equal(int(mouse["cells_after_qc"]) - int(mouse["doublets_removed"].split()[0]),
                 162_175, "mouse cells analysed")
    result.equal(int(mouse["n_clusters"]), 29, "mouse cluster count")
    result.equal(int(human["cells_after_qc"]) - int(human["doublets_removed"].split()[0]),
                 27_729, "human cells analysed")
    result.equal(int(human["n_clusters"]), 31, "human cluster count")

    proposals = read_csv("analysis/GSE262927/tables/cluster_annotation_proposals.csv")
    result.equal(len(proposals), 29, "mouse annotation proposals")
    contradicted = sum("DISAGREES" in row["Deposition check"] for row in proposals)
    result.equal(contradicted, 3, "contradicted mouse annotations")

    fractions = read_csv("analysis/GSE262927/tables/cluster_vs_author_celltype_fraction.csv")
    purity = statistics.median(
        max(float(value) for key, value in row.items() if key != "leiden_cluster")
        for row in fractions
    )
    result.equal(round(purity, 3), 0.947, "median cluster purity")

    transitional = median_by(
        read_csv("analysis/GSE262927/regeneration_focus/tables/transitional_abundance_per_sample.csv"),
        "day", "pct_transitional",
    )
    result.equal(round(transitional[11.0], 1), 27.4, "11 dpi transitional median")
    result.equal(round(transitional[366.0], 1), 0.3, "366 dpi transitional median")

    icap = median_by(
        read_csv("analysis/GSE262927/regeneration_focus/tables/icap_abundance_per_sample.csv"),
        "day", "pct",
    )
    result.equal(round(icap[0.0], 1), 2.0, "homeostatic iCAP median")
    result.equal(round(icap[25.0], 1), 37.5, "25 dpi iCAP median")
    result.equal(round(icap[366.0], 1), 21.7, "366 dpi iCAP median")

    tracing = read_csv(
        "analysis/GSE262927/lineage_tracing_cohort/tables/icap_tracing_by_cre_line.csv"
    )
    kit_rates = [
        float(row["pct_traced_in_iCAP"])
        for row in tracing
        if row["cre_line"] == "Kit-MerCreMer"
    ]
    result.equal(len(kit_rates), 3, "Kit lineage-tracing animals")
    result.equal(round(min(kit_rates)), 33, "minimum Kit iCAP trace rate")
    result.equal(round(max(kit_rates)), 53, "maximum Kit iCAP trace rate")


def main() -> int:
    result = Validation()
    required = [
        "README.md",
        "FINDINGS.md",
        "REPRODUCIBILITY.md",
        "LICENSE",
        "CITATION.cff",
        ".github/workflows/portfolio-checks.yml",
        "analysis/requirements.txt",
        "analysis/raw_data_inventory.json",
        "docs/PIPELINE_AS_RUN.md",
    ]
    for relative in required:
        result.require((REPO / relative).exists(), f"required file missing: {relative}")

    check_markdown_links(result)
    check_machine_readable_files(result)
    check_headline_results(result)

    if result.failures:
        print(f"Portfolio validation FAILED ({len(result.failures)} of {result.checks} checks):")
        for failure in result.failures:
            print(f"  - {failure}")
        return 1

    print(f"Portfolio validation passed: {result.checks} checks")
    print("  local Markdown links resolve")
    print("  tracked JSON artefacts parse")
    print("  headline counts and biological results match tracked tables")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
