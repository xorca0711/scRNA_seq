"""A15 rival-2 side-branch, stage 1: design and instrument audit. Computes no endpoint.

The side-branch declared in PLAN.md bounds rival 2 only: whether blocking epithelial
integrin beta6 moves the epithelium's own programme. It is not a test of A15, and this
script must not compute any gene-set score.

Nine items, each established from the deposited files rather than from the series
abstract:

1. the mouse-level design, recomputed from the series matrix;
2. the counts-column join, built by rule and verified as a bijection over all 48
   libraries, because the deposited column names use two naming conventions;
3. an independent check that the join rule reproduces A1's ten recorded assignments
   for this same deposit, which A1 built for a different contrast;
4. the naming hazard this audit found, reported whatever it implies;
5. the confounding audit: batch, sex, and whether an injury control shares the
   antibody experiment's batch;
6. coverage of the frozen endpoint gene sets against the counts gene index, through
   an Ensembl symbol lookup, with exact failures reported and no alias repair;
7. a depth and detection covariate audit for the eligible libraries;
8. power, computed exactly for the declared rank statistic at the observed sizes;
9. a precedent check, recording in full that A1 already used this deposit.

Standard library only. The counts file is hash-verified against A1's recorded hash.
Refuses to overwrite its own outputs.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
import sys
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / "cache"
OUT = HERE / "tables/rival2"

ACCESSION = "GSE190821"
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190821/matrix/"
    "GSE190821_series_matrix.txt.gz"
)
COUNTS_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn/GSE190821/suppl/"
    "GSE190821_counts.csv.gz"
)
MATRIX = CACHE / "GSE190821_series_matrix.txt.gz"
COUNTS = CACHE / "GSE190821_counts.csv.gz"
SYMBOL_MAP = CACHE / "ensembl_symbol_lookup.json"

A1 = ROOT / "RQ_Specified/A1_transitional_epithelial_state_distinction"
A1_CONFIG = A1 / "config/ire1_kira8.json"
A1_QC = A1 / "tables/ire1/count_QC.tsv"

A0_FROZEN = ROOT / "RQ_Specified/A0_conserved_epithelial_transition_program/tables/pilot_v1/frozen_programme.json"
A5A11_FROZEN = ROOT / "RQ_Specified/A5_A11_shared_component_contract/tables/frozen_modules.json"

PERTURBED_GENE = "Itgb6"
COMPARTMENT_TOKEN = {"E": "Epithelium", "I": "Whole Lung"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=300) as response, dest.open("wb") as fh:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            fh.write(chunk)


def parse_series_matrix(path: Path) -> tuple[list[dict], dict]:
    text = gzip.open(path, "rt", encoding="utf-8", errors="replace").read()
    titles: list[str] = []
    accessions: list[str] = []
    characteristics: list[list[str]] = []
    series: dict[str, list[str]] = defaultdict(list)
    for line in text.split("\n"):
        if not line.strip():
            continue
        parts = [p.strip('"') for p in line.rstrip("\n").split("\t")]
        key, values = parts[0], parts[1:]
        if key == "!Sample_title":
            titles = values
        elif key == "!Sample_geo_accession":
            accessions = values
        elif key == "!Sample_characteristics_ch1":
            characteristics.append(values)
        elif key.startswith("!Series_"):
            series[key[len("!Series_"):]].extend(values)
    libraries = [{"gsm": accessions[i], "title": titles[i]} for i in range(len(accessions))]
    for column in characteristics:
        field = column[0].split(":")[0].strip()
        for i, value in enumerate(column):
            libraries[i][field] = value.split(":", 1)[1].strip()
    return libraries, dict(series)


def join_column(column: str, mouse_ids: set[str]) -> tuple[str | None, str | None, list[str]]:
    """Resolve a counts column to (mouse, compartment) by rule, with the rules used.

    Two deposited naming conventions appear, and two token variants: batch S061 mice
    are zero-padded to three digits, and one whole-lung column carries an FT suffix.
    Every rule applied is returned so the audit can report which columns needed which.
    """
    tokens = column.split("_")
    compartment = None
    index = None
    for i, token in enumerate(tokens):
        if token in COMPARTMENT_TOKEN:
            compartment = COMPARTMENT_TOKEN[token]
            index = i
            break
    if compartment is None:
        return None, None, ["no compartment token"]
    candidates = [t for j, t in enumerate(tokens) if j != index and j > 0]
    rules: list[str] = []
    for token in candidates:
        if token in mouse_ids:
            return token, compartment, ["exact"]
        stripped = token.lstrip("0")
        if stripped and stripped in mouse_ids:
            rules.append("zero-padded")
            return stripped, compartment, rules
        if token.endswith("FT"):
            base = token[:-2]
            if base in mouse_ids:
                rules.append("FT suffix")
                return base, compartment, rules
            base_stripped = base.lstrip("0")
            if base_stripped and base_stripped in mouse_ids:
                rules.append("FT suffix and zero-padded")
                return base_stripped, compartment, rules
    return None, compartment, ["no mouse token matched"]


def ensembl_symbols(symbols: list[str]) -> dict:
    """Symbol to stable gene ID for mouse, through the same Ensembl REST endpoint A1 used.

    Current stable IDs and labels only. No genomic coordinates are transported to
    GRCm38, and no alias repair is attempted: a symbol that does not resolve exactly
    is reported as a failure.
    """
    if SYMBOL_MAP.exists():
        cached = json.loads(SYMBOL_MAP.read_text(encoding="utf-8"))
        if set(cached.get("requested_symbols", [])) >= set(symbols):
            return cached
    records: dict[str, dict] = {}
    unique = sorted(set(symbols))
    for start in range(0, len(unique), 500):
        chunk = unique[start:start + 500]
        payload = json.dumps({"symbols": chunk}).encode("utf-8")
        request = urllib.request.Request(
            "https://rest.ensembl.org/lookup/symbol/mus_musculus",
            data=payload,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        with urllib.request.urlopen(request, timeout=300) as response:
            records.update(json.loads(response.read().decode("utf-8")))
    result = {
        "source_url": "https://rest.ensembl.org/lookup/symbol/mus_musculus",
        "retrieved_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "annotation_scope": "Current stable gene IDs and labels; no genomic coordinates "
                            "transported to GRCm38; no alias repair",
        "requested_symbols": unique,
        "records": records,
    }
    SYMBOL_MAP.parent.mkdir(parents=True, exist_ok=True)
    SYMBOL_MAP.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def exact_rank_sum_min_p(n1: int, n2: int) -> float:
    return 1.0 / math.comb(n1 + n2, n1)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    run_record = OUT / "stage1_design_run.json"
    if run_record.exists() and not args.force:
        print(f"refusing to overwrite {run_record}; pass --force to replace")
        return 2

    fetch(MATRIX_URL, MATRIX)
    fetch(COUNTS_URL, COUNTS)

    counts_sha = sha256(COUNTS)
    a1_config = json.loads(A1_CONFIG.read_text(encoding="utf-8"))
    if counts_sha != a1_config["count_source_sha256"]:
        print("STOP: counts file does not match A1's recorded hash for this deposit")
        print(f"  mine: {counts_sha}")
        print(f"  A1:   {a1_config['count_source_sha256']}")
        return 1

    libraries, series = parse_series_matrix(MATRIX)
    if len(libraries) != 48:
        print(f"STOP: expected 48 libraries, found {len(libraries)}")
        return 1

    # Item 1: mouse-level design.
    mice: dict[str, dict] = {}
    for record in libraries:
        mouse = record["mouse identifier"]
        entry = mice.setdefault(mouse, {"compartments": []})
        entry["compartments"].append(record["tissue compartment"])
        for field in ("mouse sex", "exposure", "treatment", "treatment type", "batch"):
            previous = entry.get(field)
            if previous is not None and previous != record[field]:
                print(f"STOP: mouse {mouse} disagrees on {field}")
                return 1
            entry[field] = record[field]
    compartment_shapes = Counter(tuple(sorted(e["compartments"])) for e in mice.values())
    one_each = compartment_shapes == Counter({("Epithelium", "Whole Lung"): len(mice)})

    # Item 2: the counts-column join.
    with gzip.open(COUNTS, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        gene_index = [row[0] for row in reader if row]
    columns = header[1:]
    mouse_ids = set(mice)
    resolved: dict[tuple[str, str], list[str]] = {}
    rules_used: dict[str, list[str]] = {}
    unresolved: list[str] = []
    for column in columns:
        mouse, compartment, rules = join_column(column, mouse_ids)
        if mouse is None:
            unresolved.append(column)
            continue
        resolved.setdefault((mouse, compartment), []).append(column)
        if rules != ["exact"]:
            rules_used[column] = rules
    duplicates = {f"{k[0]}|{k[1]}": v for k, v in resolved.items() if len(v) > 1}
    metadata_pairs = {(r["mouse identifier"], r["tissue compartment"]) for r in libraries}
    missing_in_counts = sorted(f"{m}|{c}" for m, c in metadata_pairs - set(resolved))
    extra_in_counts = sorted(f"{m}|{c}" for m, c in set(resolved) - metadata_pairs)
    bijection = (
        not unresolved and not duplicates and not missing_in_counts and not extra_in_counts
        and len(resolved) == 48
    )
    column_for = {k: v[0] for k, v in resolved.items()}

    with (OUT / "stage1_join.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "counts_column", "gsm", "source_title", "mouse", "sex", "batch",
            "exposure", "treatment", "treatment_type", "compartment", "join_rule",
        ])
        for record in sorted(libraries, key=lambda r: (r["batch"], r["treatment"],
                                                       r["tissue compartment"],
                                                       r["mouse identifier"])):
            key = (record["mouse identifier"], record["tissue compartment"])
            column = column_for.get(key, "")
            writer.writerow([
                column, record["gsm"], record["title"], record["mouse identifier"],
                record["mouse sex"], record["batch"], record["exposure"],
                record["treatment"], record["treatment type"],
                record["tissue compartment"],
                ",".join(rules_used.get(column, ["exact"])),
            ])

    # Item 3: independent check of the join rule against A1's recorded assignments.
    a1_rows = list(csv.DictReader(A1_QC.read_text(encoding="utf-8").splitlines(),
                                  delimiter="\t"))
    a1_disagreements = []
    for row in a1_rows:
        key = (row["mouse"], row["compartment"])
        mine = column_for.get(key)
        if mine != row["sample_id"]:
            a1_disagreements.append({
                "a1_sample_id": row["sample_id"], "a1_mouse": row["mouse"],
                "my_column_for_that_mouse": mine,
            })
    a1_check = {
        "a1_rows_checked": len(a1_rows),
        "a1_contrast": f"{a1_config['reference']} against {a1_config['case']}",
        "a1_excluded": a1_config["excluded"],
        "disagreements": a1_disagreements,
        "rule_reproduces_a1": not a1_disagreements,
    }

    # Item 4: the naming hazard.
    prefix_groups: dict[str, set[str]] = defaultdict(set)
    for record in libraries:
        key = (record["mouse identifier"], record["tissue compartment"])
        column = column_for.get(key)
        if not column:
            continue
        tokens = column.split("_")
        prefix = "_".join(
            t for t in tokens
            if t not in COMPARTMENT_TOKEN and not t.rstrip("RFT").isdigit()
            and not t.lstrip("0").isdigit()
        )
        prefix_groups[prefix].add(f"{record['exposure']}/{record['treatment']}")
    ambiguous_prefixes = {p: sorted(g) for p, g in prefix_groups.items() if len(g) > 1}

    # Item 5: the confounding audit.
    def arm(exposure: str, treatment: str) -> list[str]:
        return sorted(m for m, e in mice.items()
                      if e["exposure"] == exposure and e["treatment"] == treatment)

    arms = {
        "bleomycin_3G9": arm("Bleomycin", "3G9"),
        "bleomycin_Axum8": arm("Bleomycin", "Axum8"),
        "saline_Axum8": arm("Saline", "Axum8"),
        "bleomycin_Vehicle": arm("Bleomycin", "Vehicle"),
        "bleomycin_KIRA8": arm("Bleomycin", "KIRA8"),
        "saline_Vehicle": arm("Saline", "Vehicle"),
    }
    arm_detail = {
        name: {
            "n_mice": len(members),
            "mice": members,
            "batches": sorted({mice[m]["batch"] for m in members}),
            "sex": dict(sorted(Counter(mice[m]["mouse sex"] for m in members).items())),
        }
        for name, members in arms.items()
    }
    antibody_batches = sorted(
        set(arm_detail["bleomycin_3G9"]["batches"])
        | set(arm_detail["bleomycin_Axum8"]["batches"])
        | set(arm_detail["saline_Axum8"]["batches"])
    )
    confounding = {
        "primary_arms_share_one_batch": (
            arm_detail["bleomycin_3G9"]["batches"]
            == arm_detail["bleomycin_Axum8"]["batches"]
            and len(arm_detail["bleomycin_3G9"]["batches"]) == 1
        ),
        "antibody_experiment_batches": antibody_batches,
        "sex_matched_between_primary_arms": (
            arm_detail["bleomycin_3G9"]["sex"] == arm_detail["bleomycin_Axum8"]["sex"]
        ),
        "injury_control_in_same_batch": (
            arm_detail["saline_Axum8"]["batches"] == antibody_batches
        ),
        "vehicle_is_not_the_antibody_control": a1_config["excluded"],
    }

    # Item 6: endpoint coverage, through an Ensembl symbol lookup.
    frozen_a0 = json.loads(A0_FROZEN.read_text(encoding="utf-8"))
    a0_genes = list(frozen_a0["mouse_genes"])
    frozen_modules = json.loads(A5A11_FROZEN.read_text(encoding="utf-8"))
    injury_residual = list(frozen_modules["modules"]["injury_residual"]["genes"])
    shared_remodelling = list(frozen_modules["modules"]["shared_remodelling"]["genes"])
    a1_markers = list(a1_config["markers"])

    requested = sorted(set(a0_genes) | set(injury_residual) | set(shared_remodelling)
                       | set(a1_markers) | {PERTURBED_GENE})
    lookup = ensembl_symbols(requested)
    symbol_to_id = {
        symbol: record["id"]
        for symbol, record in lookup["records"].items()
        if isinstance(record, dict) and record.get("id")
    }
    index_set = set(gene_index)

    def coverage(name: str, genes: list[str]) -> dict:
        unresolved_symbols = sorted(set(genes) - set(symbol_to_id))
        mapped = {g: symbol_to_id[g] for g in genes if g in symbol_to_id}
        present = {g: i for g, i in mapped.items() if i in index_set}
            # A symbol may resolve at Ensembl and still be absent from a GRCm38 count
            # matrix; both losses are reported separately rather than pooled.
        return {
            "declared": len(set(genes)),
            "symbol_unresolved": unresolved_symbols,
            "mapped": len(mapped),
            "present_in_counts": len(present),
            "fraction_present": round(len(present) / len(set(genes)), 4) if genes else 0.0,
            "mapped_but_absent_from_counts": sorted(set(mapped) - set(present)),
            "contains_perturbed_gene": PERTURBED_GENE in set(genes),
        }

    endpoint_coverage = {
        "a0_frozen_transition_50": coverage("a0", a0_genes),
        "a5a11_injury_residual": coverage("injury", injury_residual),
        "a5a11_shared_remodelling": coverage("shared", shared_remodelling),
        "a1_frozen_marker_panel": coverage("markers", a1_markers),
    }

    with (OUT / "stage1_endpoint_coverage.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene_set", "symbol", "ensembl_id", "in_counts_index"])
        for name, genes in (("a0_frozen_transition_50", a0_genes),
                            ("a5a11_injury_residual", injury_residual),
                            ("a5a11_shared_remodelling", shared_remodelling),
                            ("a1_frozen_marker_panel", a1_markers)):
            for symbol in sorted(set(genes)):
                gid = symbol_to_id.get(symbol, "")
                writer.writerow([name, symbol, gid, bool(gid) and gid in index_set])

    # Item 7: depth and detection for the eligible libraries. Library totals are a
    # covariate, not the endpoint; no gene set is scored.
    eligible_columns = {}
    for name in ("bleomycin_3G9", "bleomycin_Axum8", "saline_Axum8"):
        for mouse in arms[name]:
            column = column_for[(mouse, "Epithelium")]
            eligible_columns[column] = (name, mouse)
    column_position = {c: i + 1 for i, c in enumerate(columns)}
    totals = {c: 0 for c in eligible_columns}
    detected = {c: 0 for c in eligible_columns}
    with gzip.open(COUNTS, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        next(reader)
        for row in reader:
            if not row:
                continue
            for column in eligible_columns:
                value = int(row[column_position[column]])
                totals[column] += value
                if value > 0:
                    detected[column] += 1

    with (OUT / "stage1_covariates.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["counts_column", "arm", "mouse", "library_total", "genes_detected"])
        for column, (name, mouse) in sorted(eligible_columns.items(),
                                            key=lambda kv: (kv[1][0], kv[1][1])):
            writer.writerow([column, name, mouse, totals[column], detected[column]])

    depth_floor = 100_000  # A1's recorded minimum library count for this deposit
    below_floor = sorted(c for c in eligible_columns if totals[c] < depth_floor)
    depth_span = {
        "minimum": min(totals.values()),
        "maximum": max(totals.values()),
        "ratio_max_over_min": round(max(totals.values()) / min(totals.values()), 3),
        "floor": depth_floor,
        "below_floor": below_floor,
    }

    # Item 8: exact power.
    n_treated = len(arms["bleomycin_3G9"])
    n_control = len(arms["bleomycin_Axum8"])
    n_saline = len(arms["saline_Axum8"])
    power = {
        "primary_contrast": {
            "arms": f"bleomycin 3G9 n={n_treated} against bleomycin Axum8 n={n_control}",
            "assignments": math.comb(n_treated + n_control, n_treated),
            "smallest_attainable_one_sided_p": round(exact_rank_sum_min_p(n_treated, n_control), 6),
            "smallest_attainable_two_sided_p": round(2 * exact_rank_sum_min_p(n_treated, n_control), 6),
            "requires": "complete separation of the two arms on the score",
        },
        "instrument_check": {
            "arms": f"bleomycin Axum8 n={n_control} against saline Axum8 n={n_saline}",
            "assignments": math.comb(n_control + n_saline, n_saline),
            "smallest_attainable_one_sided_p": round(exact_rank_sum_min_p(n_control, n_saline), 6),
            "requires": "complete separation of the two arms on the score",
        },
        "precise_absence": "unavailable at these sizes, declared before any test",
    }

    # Item 9: precedent, recorded in full.
    precedent_hits = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".md", ".json", ".py", ".tsv", ".csv"}:
            continue
        if CACHE in path.parents or ".git" in path.parts:
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        if ACCESSION in text:
            precedent_hits.append(str(path.relative_to(ROOT)).replace("\\", "/"))
    precedent_hits.sort()
    precedent = {
        "files_mentioning_accession": precedent_hits,
        "prior_analysis": {
            "question": "A1",
            "config": str(A1_CONFIG.relative_to(ROOT)).replace("\\", "/"),
            "study_named_by_a1": a1_config["study"],
            "compartment": "Epithelium",
            "exposure": "Bleomycin",
            "arms_scored_by_a1": [a1_config["reference"], a1_config["case"]],
            "arms_explicitly_excluded_by_a1": a1_config["excluded"],
            "counts_hash_matches": True,
            "this_contrast_touched_by_a1": False,
            "instrument_reusable_in_this_session": False,
            "why_not": "A1's instrument is edgeR and limma in R; no R interpreter is "
                       "available in this session, so its TMM normalisation and "
                       "quasi-likelihood fit cannot be reused. A score-level statistic "
                       "is declared instead and is not presented as A1's instrument.",
        },
    }

    stop_rules = {
        "counts_hash_mismatch": {
            "rule": "the counts file differs from A1's recorded hash",
            "triggered": False,
        },
        "join_not_bijective": {
            "rule": "the counts columns do not resolve one-to-one onto the 48 libraries",
            "triggered": not bijection,
        },
        "join_rule_disagrees_with_a1": {
            "rule": "the join rule does not reproduce A1's recorded assignments",
            "triggered": bool(a1_disagreements),
        },
        "arms_too_small": {
            "rule": "fewer than 3 mice in either primary arm",
            "triggered": n_treated < 3 or n_control < 3,
        },
        "treatment_confounded_with_batch": {
            "rule": "the primary arms do not occupy the same single batch",
            "triggered": not confounding["primary_arms_share_one_batch"],
        },
        "no_injury_control_in_batch": {
            "rule": "no saline arm shares the antibody experiment's batch",
            "triggered": not confounding["injury_control_in_same_batch"],
        },
        "endpoint_coverage": {
            "rule": "the primary endpoint recovers fewer than 0.90 of its declared genes",
            "triggered": endpoint_coverage["a0_frozen_transition_50"]["fraction_present"] < 0.90,
        },
        "primary_endpoint_contains_perturbed_gene": {
            "rule": "the primary endpoint contains the perturbed gene",
            "triggered": endpoint_coverage["a0_frozen_transition_50"]["contains_perturbed_gene"],
        },
        "eligible_library_below_depth_floor": {
            "rule": f"an eligible epithelial library falls below {depth_floor} counts",
            "triggered": bool(below_floor),
        },
        "prior_analysis_scored_this_contrast": {
            "rule": "a prior analysis already scored the 3G9 against Axum8 contrast",
            "triggered": False,
        },
    }
    triggered = sorted(k for k, v in stop_rules.items() if v["triggered"])

    record = {
        "schema": "a15-rival2-stage1/v1",
        "stage": "1, design and instrument audit; no gene-set score computed",
        "scope": "bounds A15 rival 2 only; not a test of A15",
        "accession": ACCESSION,
        "series_title": series.get("title", [""])[0],
        "series_overall_design": series.get("overall_design", [""])[0],
        "series_pubmed_id": series.get("pubmed_id", [""])[0] if series.get("pubmed_id") else "",
        "platform": series.get("platform_id", [""])[0],
        "inputs": {
            "series_matrix": {"url": MATRIX_URL, "bytes": MATRIX.stat().st_size,
                              "sha256": sha256(MATRIX)},
            "counts": {"url": COUNTS_URL, "bytes": COUNTS.stat().st_size,
                       "sha256": counts_sha,
                       "verified_against": "RQ_Specified/A1_transitional_epithelial_state_distinction/config/ire1_kira8.json"},
            "ensembl_symbol_lookup": {"path": "cache/ensembl_symbol_lookup.json",
                                      "sha256": sha256(SYMBOL_MAP)},
            "a0_frozen_programme": {"path": str(A0_FROZEN.relative_to(ROOT)).replace("\\", "/"),
                                    "sha256": sha256(A0_FROZEN)},
            "a5a11_frozen_modules": {"path": str(A5A11_FROZEN.relative_to(ROOT)).replace("\\", "/"),
                                     "sha256": sha256(A5A11_FROZEN)},
            "a1_config": {"path": str(A1_CONFIG.relative_to(ROOT)).replace("\\", "/"),
                          "sha256": sha256(A1_CONFIG)},
        },
        "design": {
            "libraries": len(libraries),
            "mice": len(mice),
            "one_epithelial_and_one_whole_lung_per_mouse": one_each,
            "compartment_shapes": {str(k): v for k, v in compartment_shapes.items()},
            "arms": arm_detail,
        },
        "join": {
            "bijective_over_48_libraries": bijection,
            "columns_needing_a_rule_beyond_exact_match": rules_used,
            "unresolved_columns": unresolved,
            "duplicate_pairs": duplicates,
            "metadata_pairs_missing_from_counts": missing_in_counts,
            "counts_pairs_missing_from_metadata": extra_in_counts,
        },
        "independent_check_against_a1": a1_check,
        "naming_hazard": {
            "prefixes_covering_more_than_one_group": ambiguous_prefixes,
            "consequence": "column names alone cannot identify the antibody control arm; "
                           "the treatment field in the series metadata is the only source "
                           "that separates Axum8 from Vehicle",
        },
        "confounding": confounding,
        "counts_structure": {
            "genes_in_index": len(gene_index),
            "unique_gene_ids": len(index_set),
            "sample_columns": len(columns),
            "identifier_type": "Ensembl mouse stable gene IDs",
            "identifier_example": gene_index[:3],
        },
        "endpoint_coverage": endpoint_coverage,
        "perturbed_gene": {
            "symbol": PERTURBED_GENE,
            "ensembl_id": symbol_to_id.get(PERTURBED_GENE, ""),
            "in_counts_index": symbol_to_id.get(PERTURBED_GENE, "") in index_set,
            "in_endpoint_sets": {k: v["contains_perturbed_gene"]
                                 for k, v in endpoint_coverage.items()},
        },
        "covariates": {"eligible_epithelial_libraries": len(eligible_columns),
                       "depth": depth_span},
        "power": power,
        "precedent": precedent,
        "stop_rules": stop_rules,
        "stop_rules_triggered": triggered,
        "endpoint_scored": False,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    run_record.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print(f"counts hash matches A1's record: True")
    print(f"libraries {len(libraries)}, mice {len(mice)}, join bijective: {bijection}")
    print(f"join rules beyond exact match: {rules_used}")
    print(f"rule reproduces A1's 10 assignments: {a1_check['rule_reproduces_a1']}")
    print(f"ambiguous column prefixes: {ambiguous_prefixes}")
    print(f"3G9 n={n_treated} {arm_detail['bleomycin_3G9']['mice']} batches "
          f"{arm_detail['bleomycin_3G9']['batches']} sex {arm_detail['bleomycin_3G9']['sex']}")
    print(f"Axum8 bleo n={n_control} {arm_detail['bleomycin_Axum8']['mice']} batches "
          f"{arm_detail['bleomycin_Axum8']['batches']} sex {arm_detail['bleomycin_Axum8']['sex']}")
    print(f"Axum8 saline n={n_saline} {arm_detail['saline_Axum8']['mice']}")
    print(f"same batch: {confounding['primary_arms_share_one_batch']}, "
          f"sex matched: {confounding['sex_matched_between_primary_arms']}, "
          f"injury control in batch: {confounding['injury_control_in_same_batch']}")
    for name, cov in endpoint_coverage.items():
        print(f"coverage {name}: {cov['present_in_counts']}/{cov['declared']} "
              f"({cov['fraction_present']}) unresolved={cov['symbol_unresolved'][:4]} "
              f"perturbed_gene_inside={cov['contains_perturbed_gene']}")
    print(f"depth: {depth_span['minimum']} to {depth_span['maximum']} "
          f"(ratio {depth_span['ratio_max_over_min']}), below floor: {below_floor or 'none'}")
    print(f"power: primary min one-sided p = "
          f"{power['primary_contrast']['smallest_attainable_one_sided_p']}, "
          f"instrument check = {power['instrument_check']['smallest_attainable_one_sided_p']}")
    print(f"stop rules triggered: {triggered or 'none'}")
    print(f"wrote {run_record}")
    return 1 if triggered else 0


if __name__ == "__main__":
    raise SystemExit(main())
