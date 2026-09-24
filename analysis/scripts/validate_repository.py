#!/usr/bin/env python
"""Validate the lightweight, tracked contract of the repository.

This check deliberately uses only the Python standard library. It verifies the
published headline numbers against tracked tables, checks machine-readable
artefacts, and catches broken local Markdown links without requiring raw data
or the scanpy environment.
"""

from __future__ import annotations

import csv
import json
import os
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


REPO = Path(__file__).resolve().parents[2]
LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
SKIP_DIRS = {".git", ".venv", ".venv-x64", ".venv-repro", ".tools", ".claude", "raw_data", "__pycache__", "cache"}


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


def source_files(root: Path, suffix: str) -> list[Path]:
    files = []
    for directory, subdirectories, names in os.walk(root):
        subdirectories[:] = [name for name in subdirectories if name not in SKIP_DIRS]
        files.extend(Path(directory) / name for name in names if name.endswith(suffix))
    return files


def markdown_files() -> list[Path]:
    return source_files(REPO, ".md")


def heading_slugs(text: str) -> set[str]:
    """GitHub-flavoured anchor slugs for every heading in a document.

    Lowercase, drop anything that is not a word character, space or hyphen,
    then hyphenate the spaces. Duplicate headings get a numeric suffix on
    GitHub, so every slug is also accepted with one.
    """
    slugs: set[str] = set()
    for line in text.splitlines():
        if not line.startswith("#"):
            continue
        heading = line.lstrip("#").strip()
        if not heading:
            continue
        slug = re.sub(r"[^\w\s-]", "", heading.lower())
        # Each whitespace character becomes one hyphen; runs are NOT collapsed,
        # which is what GitHub does and why "5 · Negative results" anchors as
        # "5--negative-results".
        slug = re.sub(r"\s", "-", slug).strip("-")
        if slug:
            slugs.add(slug)
            slugs.update(f"{slug}-{n}" for n in range(1, 6))
    return slugs


def check_claim_ids_unique(result: Validation) -> None:
    """Every row of the claims register must own its identifier.

    Two papers numbering from the same point produced eleven duplicated
    identifiers on 2026-09-17, which nothing caught because the register was
    only ever read by humans. A duplicate makes every cross-reference to that
    number ambiguous, so it fails the build.
    """
    claims = REPO / "CLAIMS.md"
    if not claims.exists():
        return
    seen: dict[str, int] = defaultdict(int)
    for line in claims.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\|\s*(C\d+)\s*\|", line)
        if match:
            seen[match.group(1)] += 1
    duplicates = sorted((k for k, v in seen.items() if v > 1), key=lambda s: int(s[1:]))
    result.require(
        not duplicates,
        "duplicate claim identifier in CLAIMS.md: " + ", ".join(duplicates),
    )
    for identifier in sorted(seen, key=lambda s: int(s[1:])):
        result.require(True, f"claim {identifier} is unique")


def check_markdown_links(result: Validation) -> None:
    slug_cache: dict[Path, set[str]] = {}

    def slugs_for(path: Path) -> set[str]:
        if path not in slug_cache:
            try:
                slug_cache[path] = heading_slugs(path.read_text(encoding="utf-8"))
            except OSError:
                slug_cache[path] = set()
        return slug_cache[path]

    for document in markdown_files():
        body = document.read_text(encoding="utf-8")
        for match in LINK_RE.finditer(body):
            raw_target = match.group(1).strip().strip("<>")
            if raw_target.startswith(("http://", "https://", "mailto:")):
                continue
            # A same-document anchor: check it against this document's headings.
            if raw_target.startswith("#"):
                fragment = unquote(raw_target[1:]).split(' "', 1)[0]
                if fragment:
                    result.require(
                        fragment in slugs_for(document),
                        f"broken anchor in {document.relative_to(REPO)}: {raw_target}",
                    )
                continue
            target = unquote(raw_target.split("#", 1)[0])
            fragment = unquote(raw_target.split("#", 1)[1]) if "#" in raw_target else ""
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
            # A fragment on another Markdown file is checked too, so a table of
            # contents cannot drift from the headings it points at.
            if fragment and resolved.exists() and resolved.suffix == ".md":
                result.require(
                    fragment.split(' "', 1)[0] in slugs_for(resolved),
                    f"broken anchor in {document.relative_to(REPO)}: {raw_target}",
                )


def check_machine_readable_files(result: Validation) -> None:
    for path in source_files(REPO / "analysis", ".json") + source_files(REPO / "Thesis", ".json"):
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
        not (REPO / "Thesis/gate1_01_niethamer_2025/GSE262927/inventory/raw_files.json").exists(),
        "obsolete dataset-specific raw-data inventory still exists",
    )


def check_headline_results(result: Validation) -> None:
    mouse = read_json("Thesis/gate1_01_niethamer_2025/GSE262927/logs/decisions.json")
    human = read_json("Thesis/ungated_murthy_2022/GSE178360/logs/decisions.json")
    result.equal(int(mouse["cells_after_qc"]) - int(mouse["doublets_removed"].split()[0]),
                 162_175, "mouse cells analysed")
    result.equal(int(mouse["n_clusters"]), 29, "mouse cluster count")
    result.equal(int(human["cells_after_qc"]) - int(human["doublets_removed"].split()[0]),
                 27_729, "human cells analysed")
    result.equal(int(human["n_clusters"]), 31, "human cluster count")

    proposals = read_csv("Thesis/gate1_01_niethamer_2025/GSE262927/tables/cluster_annotation_proposals.csv")
    result.equal(len(proposals), 29, "mouse annotation proposals")
    contradicted = sum("DISAGREES" in row["Deposition check"] for row in proposals)
    result.equal(contradicted, 3, "contradicted mouse annotations")

    fractions = read_csv("Thesis/gate1_01_niethamer_2025/GSE262927/tables/cluster_vs_author_celltype_fraction.csv")
    purity = statistics.median(
        max(float(value) for key, value in row.items() if key != "leiden_cluster")
        for row in fractions
    )
    result.equal(round(purity, 3), 0.947, "median cluster purity")

    transitional = median_by(
        read_csv("Thesis/gate1_01_niethamer_2025/GSE262927/regeneration_focus/tables/transitional_abundance_per_sample.csv"),
        "day", "pct_transitional",
    )
    result.equal(round(transitional[11.0], 1), 27.4, "11 dpi transitional median")
    result.equal(round(transitional[366.0], 1), 0.3, "366 dpi transitional median")

    icap = median_by(
        read_csv("Thesis/gate1_01_niethamer_2025/GSE262927/regeneration_focus/tables/icap_abundance_per_sample.csv"),
        "day", "pct",
    )
    result.equal(round(icap[0.0], 1), 2.0, "homeostatic iCAP median")
    result.equal(round(icap[25.0], 1), 37.5, "25 dpi iCAP median")
    result.equal(round(icap[366.0], 1), 21.7, "366 dpi iCAP median")

    tracing = read_csv(
        "Thesis/gate1_01_niethamer_2025/GSE262927/lineage_tracing_cohort/tables/icap_tracing_by_cre_line.csv"
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
        ".github/workflows/repository-checks.yml",
        "analysis/requirements.txt",
        "analysis/raw_data_inventory.json",
        "docs/PIPELINE_AS_RUN.md",
    ]
    for relative in required:
        result.require((REPO / relative).exists(), f"required file missing: {relative}")

    check_markdown_links(result)
    check_claim_ids_unique(result)
    check_machine_readable_files(result)
    check_headline_results(result)
    sys.path.insert(0, str(REPO))
    from analysis.lib.research_layout import check_research_layout
    check_research_layout(REPO, result)
    from claim_contract import generated, verify_bindings, MANIFEST, SUMMARY, BINDINGS
    bindings = read_json(str(BINDINGS))["bindings"]
    for error in verify_bindings(bindings, REPO):
        result.require(False, error)
    for relative, expected in zip((MANIFEST, SUMMARY), generated(REPO)):
        path = REPO / relative
        result.require(path.exists() and path.read_text(encoding="utf-8") == expected,
                       f"stale generated evidence contract: {relative}")

    if result.failures:
        print(f"Repository validation FAILED ({len(result.failures)} of {result.checks} checks):")
        for failure in result.failures:
            print(f"  - {failure}")
        return 1

    print(f"Repository validation passed: {result.checks} checks")
    print("  local Markdown links resolve")
    print("  tracked JSON artefacts parse")
    print("  selected headline values and explicit numeric bindings match tracked tables")
    print("  claim index and summary match the register; this is not full scientific reproduction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
