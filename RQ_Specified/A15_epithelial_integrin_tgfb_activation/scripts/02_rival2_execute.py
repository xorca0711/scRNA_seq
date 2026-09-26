"""A15 rival-2 stage 3: compute every declared contrast at once, under freeze v2.

The freeze forbids a gate between any check and the primary, so this script computes the
covariates, the engagement control, both primary endpoints, the context contrast and every
sensitivity in one non-interactive pass and emits a single run record with the verdict the
declared rules produce. No human decision sits in the middle.

It bounds A15 rival 2 only. It cannot test A15.

Standard library only, so every exact test is enumerated here; script 03 checks them against
scipy and numpy independently. Refuses to overwrite its own outputs.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / "cache"
OUT = HERE / "tables/rival2"
FREEZE = HERE / "config/a15_rival2_freeze_v3.json"
COUNTS = CACHE / "GSE190821_counts.csv.gz"
SYMBOL_MAP = CACHE / "ensembl_symbol_lookup.json"
JOIN = OUT / "stage1_join.tsv"
STAGE1 = OUT / "stage1_design_run.json"
A1_CONFIG = ROOT / "RQ_Specified/A1_transitional_epithelial_state_distinction/config/ire1_kira8.json"

READY = "COMMITTED_FOR_EXECUTION"
ALPHA = 2 / 70  # 0.028571, the exact two-sided floor at four against four
PURITY_THRESHOLD = 1.0
MIN_MEDIAN_COUNT = 10
MEASURABLE_FRACTION = 0.80


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def mean(xs):
    return sum(xs) / len(xs)


def median(xs):
    s = sorted(xs)
    n = len(s)
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def stdev(xs):
    if len(xs) < 2:
        return 0.0
    mu = mean(xs)
    return math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))


def mann_whitney_u(a, b):
    return sum(1 for x in a for y in b if x > y) + 0.5 * sum(1 for x in a for y in b if x == y)


def exact_rank_sum(a, b):
    combined = list(a) + list(b)
    n1, n = len(a), len(a) + len(b)
    observed = mann_whitney_u(a, b)
    below = above = total = 0
    for pick in itertools.combinations(range(n), n1):
        chosen = set(pick)
        g1 = [combined[i] for i in pick]
        g2 = [combined[i] for i in range(n) if i not in chosen]
        u = mann_whitney_u(g1, g2)
        total += 1
        if u <= observed:
            below += 1
        if u >= observed:
            above += 1
    return {
        "u": observed, "assignments": total,
        "p_one_sided_lower": round(below / total, 6),
        "p_one_sided_upper": round(above / total, 6),
        "p_two_sided": round(min(1.0, 2 * min(below, above) / total), 6),
    }


def null_u_counts(n1, n2):
    values = list(range(n1 + n2))
    counts = [0] * (n1 * n2 + 1)
    for pick in itertools.combinations(values, n1):
        chosen = set(pick)
        g1 = [values[i] for i in pick]
        g2 = [values[i] for i in range(len(values)) if i not in chosen]
        counts[int(mann_whitney_u(g1, g2))] += 1
    return counts


def hodges_lehmann(a, b, level=0.95):
    differences = sorted(x - y for x in a for y in b)
    m = len(differences)
    estimate = median(differences)
    counts = null_u_counts(len(a), len(b))
    total = sum(counts)
    # Trim k from each end, where k is the largest value whose interval still attains the
    # level. Coverage of the k-trimmed interval is 1 - 2 P(U <= k), so k = 0 is the full
    # range of the pairwise differences. v2 trimmed one too many and mislabelled the
    # coverage; simulation in this repository showed that interval covered 0.9418.
    cumulative = 0
    trim = None
    attained = None
    for u, count in enumerate(counts):
        cumulative += count
        coverage = 1 - 2 * cumulative / total
        if coverage >= level and 2 * u < m:
            trim, attained = u, round(coverage, 4)
        else:
            break
    if trim is None:
        return {"estimate": round(estimate, 4), "lower": round(differences[0], 4),
                "upper": round(differences[-1], 4), "attained_coverage": None,
                "note": "no interval attains the requested level at these sample sizes"}
    return {"estimate": round(estimate, 4), "lower": round(differences[trim], 4),
            "upper": round(differences[m - 1 - trim], 4), "attained_coverage": attained,
            "trimmed_from_each_end": trim}


def pearson(xs, ys):
    mx, my = mean(xs), mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else 0.0


def contrast(case, reference, label):
    test = exact_rank_sum(case, reference)
    return {
        "label": label,
        "case_values": [round(v, 4) for v in case],
        "reference_values": [round(v, 4) for v in reference],
        "case_mean": round(mean(case), 4), "reference_mean": round(mean(reference), 4),
        "case_sd": round(stdev(case), 4), "reference_sd": round(stdev(reference), 4),
        "case_range": [round(min(case), 4), round(max(case), 4)],
        "reference_range": [round(min(reference), 4), round(max(reference), 4)],
        "mean_difference": round(mean(case) - mean(reference), 4),
        "median_difference": round(median(case) - median(reference), 4),
        "pooled_within_arm_sd": round(math.sqrt((stdev(case) ** 2 + stdev(reference) ** 2) / 2), 4),
        "exact_test": test,
        "shift": hodges_lehmann(case, reference),
        "separates_at_declared_alpha": test["p_two_sided"] <= ALPHA + 1e-12,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    run_record = OUT / "stage3_execute_run.json"
    if run_record.exists() and not args.force:
        print(f"refusing to overwrite {run_record}; pass --force to replace")
        return 2

    freeze = json.loads(FREEZE.read_text(encoding="utf-8"))
    if freeze.get("status") != READY:
        print(f"STOP: freeze status is {freeze.get('status')!r}, not {READY!r}")
        return 1
    stage1 = json.loads(STAGE1.read_text(encoding="utf-8"))
    if stage1["stop_rules_triggered"]:
        print(f"STOP: stage 1 stop rules triggered: {stage1['stop_rules_triggered']}")
        return 1
    counts_sha = sha256(COUNTS)
    if counts_sha != freeze["dataset"]["counts_sha256"]:
        print("STOP: counts hash does not match the freeze")
        return 1
    if sha256(SYMBOL_MAP) != stage1["inputs"]["ensembl_symbol_lookup"]["sha256"]:
        print("STOP: the Ensembl symbol lookup changed since stage 1")
        return 1

    rows = list(csv.DictReader(JOIN.read_text(encoding="utf-8").splitlines(), delimiter="\t"))
    elig = freeze["eligibility"]

    def select(exposure, treatment, compartment):
        return sorted((r for r in rows if r["compartment"] == compartment
                       and r["exposure"] == exposure and r["treatment"] == treatment),
                      key=lambda r: r["mouse"])

    case_mice = elig["arms_primary"]["case"]["mice"]
    ref_mice = elig["arms_primary"]["reference"]["mice"]
    sal_mice = elig["arms_context_only"]["saline"]["mice"]
    epi = {
        "case": select("Bleomycin", "3G9", "Epithelium"),
        "reference": select("Bleomycin", "Axum8", "Epithelium"),
        "saline": select("Saline", "Axum8", "Epithelium"),
    }
    inp = {
        "case": select("Bleomycin", "3G9", "Whole Lung"),
        "reference": select("Bleomycin", "Axum8", "Whole Lung"),
        "saline": select("Saline", "Axum8", "Whole Lung"),
    }
    for name, declared in (("case", case_mice), ("reference", ref_mice), ("saline", sal_mice)):
        for group, what in ((epi, "epithelium"), (inp, "input")):
            got = [r["mouse"] for r in group[name]]
            if got != sorted(declared):
                print(f"STOP: {what} arm {name} is {got}, the freeze declares {sorted(declared)}")
                return 1

    libraries = {r["counts_column"]: r for g in (epi, inp) for arm in g.values() for r in arm}
    epi_primary = [r["counts_column"] for r in epi["case"] + epi["reference"]]

    # Gene sets.
    lookup = json.loads(SYMBOL_MAP.read_text(encoding="utf-8"))
    sym2id = {s: r["id"] for s, r in lookup["records"].items()
              if isinstance(r, dict) and r.get("id")}
    a1 = json.loads(A1_CONFIG.read_text(encoding="utf-8"))
    p1 = freeze["endpoints"]["primary_1_epithelial_identity_and_state"]
    panels = {
        "transitional": p1["composites"]["transitional"]["genes"],
        "identity": p1["composites"]["identity"]["genes"],
    }
    declared_here = p1["declared_here_not_inherited"]["genes"]
    engagement_genes = freeze["engagement_control_declared_before_any_value_is_read"]["genes"]
    purity = freeze["purity_and_composition_covariates_from_the_paired_input"]
    epi_markers = purity["epithelial_enrichment_expected_positive"]
    non_epi_markers = purity["non_epithelial_de_enrichment_expected_negative"]
    macrophage = purity["macrophage_guard"]["genes"]
    a0 = json.loads((ROOT / "RQ_Specified/A0_conserved_epithelial_transition_program/tables/pilot_v1/frozen_programme.json").read_text(encoding="utf-8"))
    a0_genes = list(a0["mouse_genes"])
    modules = json.loads((ROOT / "RQ_Specified/A5_A11_shared_component_contract/tables/frozen_modules.json").read_text(encoding="utf-8"))
    injury_residual = list(modules["modules"]["injury_residual"]["genes"])

    # Fail closed: a declared panel that does not map would otherwise be reported as
    # "refused" or, worse for the covariates, as "no breach", which reads as a result when
    # nothing was computed at all.
    declared_panels = {
        "primary1_transitional": panels["transitional"],
        "primary1_identity": panels["identity"],
        "engagement_control": engagement_genes,
        "epithelial_enrichment": epi_markers,
        "non_epithelial_de_enrichment": non_epi_markers,
        "macrophage_guard": macrophage,
        "declared_here": declared_here,
    }
    unmapped = {name: sorted(set(g) - set(sym2id))
                for name, g in declared_panels.items() if set(g) - set(sym2id)}
    if unmapped:
        print("STOP: declared panels are not present in the Ensembl symbol lookup, so any "
              "covariate result would be vacuous rather than negative:")
        for name, missing in unmapped.items():
            print(f"  {name}: {missing}")
        print("Re-run script 01 with --force so the lookup covers every declared symbol.")
        return 1

    itgb6_id = sym2id.get("Itgb6")
    targeted = set()
    for genes in list(panels.values()) + [declared_here, engagement_genes, epi_markers,
                                          non_epi_markers, macrophage, a0_genes,
                                          injury_residual, ["Itgb6"]]:
        for s in genes:
            if s in sym2id:
                targeted.add(sym2id[s])

    # One pass: whole matrix for the 22 libraries, plus totals.
    with gzip.open(COUNTS, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        pos = {c: header.index(c) for c in libraries}
        totals = {c: 0 for c in libraries}
        matrix: dict[str, dict[str, int]] = {}
        for row in reader:
            if not row:
                continue
            gid = row[0]
            values = {c: int(row[pos[c]]) for c in libraries}
            for c, v in values.items():
                if v < 0:
                    print(f"STOP: negative count {gid} {c}")
                    return 1
                totals[c] += v
            matrix[gid] = values
    if len(matrix) != stage1["counts_structure"]["genes_in_index"]:
        print("STOP: gene count differs from stage 1")
        return 1

    def cpm(gid, column):
        return math.log2(matrix[gid][column] / totals[column] * 1e6 + 1)

    def measurable(gid, columns):
        return median([matrix[gid][c] for c in columns]) >= MIN_MEDIAN_COUNT

    def resolve(genes, columns):
        mapped = [(s, sym2id[s]) for s in genes if s in sym2id and sym2id[s] in matrix]
        keep = [(s, g) for s, g in mapped if measurable(g, columns)]
        return mapped, keep

    def composite(genes, columns, standardise):
        mapped, keep = resolve(genes, columns)
        refused = len(keep) < MEASURABLE_FRACTION * max(1, len(set(genes)))
        if not keep:
            return None, {"declared": len(set(genes)), "mapped": len(mapped),
                          "measurable": 0, "refused": True}
        per_gene = {}
        for _, gid in keep:
            vals = {c: cpm(gid, c) for c in columns}
            if standardise:
                mu, sd = mean(list(vals.values())), stdev(list(vals.values()))
                vals = {c: ((v - mu) / sd if sd else 0.0) for c, v in vals.items()}
            per_gene[gid] = vals
        scores = {c: mean([per_gene[g][c] for _, g in keep]) for c in columns}
        return scores, {"declared": len(set(genes)), "mapped": len(mapped),
                        "measurable": len(keep), "refused": refused,
                        "measurable_symbols": [s for s, _ in keep],
                        "dropped_symbols": sorted(set(s for s, _ in mapped) - set(s for s, _ in keep))}

    epi_all = [r["counts_column"] for arm in epi.values() for r in arm]
    inp_all = [r["counts_column"] for arm in inp.values() for r in arm]
    # v3: every composite is standardised across exactly the libraries entering its own
    # contrast, with no reference to arm labels, so the exact tests stay exact.
    epi_contrast = [r["counts_column"] for r in epi["case"] + epi["reference"]]
    epi_injury = [r["counts_column"] for r in epi["reference"] + epi["saline"]]
    inp_contrast = [r["counts_column"] for r in inp["case"] + inp["reference"]]
    inp_injury = [r["counts_column"] for r in inp["reference"] + inp["saline"]]
    results = {}

    # Primary 1.
    for name, genes in panels.items():
        block = {}
        for variant, std in (("standardised", True), ("unstandardised", False)):
            scores, cov = composite(genes, epi_contrast, std)
            if scores is None:
                block[variant] = {"coverage": cov, "refused": True}
                continue
            c = contrast([scores[r["counts_column"]] for r in epi["case"]],
                         [scores[r["counts_column"]] for r in epi["reference"]],
                         f"primary1 {name} {variant}")
            c["coverage"] = cov
            inj_ctx, _ = composite(genes, epi_injury, std)
            if inj_ctx:
                c["injury_context"] = contrast(
                    [inj_ctx[r["counts_column"]] for r in epi["reference"]],
                    [inj_ctx[r["counts_column"]] for r in epi["saline"]],
                    f"injury context {name} {variant}")
            block[variant] = c
        results[f"primary1_{name}"] = block

    # Primary 1, the two genes declared here rather than inherited.
    extra = {}
    for symbol in declared_here:
        gid = sym2id.get(symbol)
        if not gid or gid not in matrix:
            extra[symbol] = {"status": "not recovered"}
            continue
        vals = {c: cpm(gid, c) for c in epi_contrast}
        extra[symbol] = contrast([vals[r["counts_column"]] for r in epi["case"]],
                                 [vals[r["counts_column"]] for r in epi["reference"]],
                                 f"declared-here {symbol}")
        extra[symbol]["measurable"] = measurable(gid, epi_contrast)
    results["declared_here_no_decision"] = extra

    # Primary 2, the omnibus.
    omnibus_genes = [g for g in matrix
                     if g != itgb6_id
                     and sum(1 for c in epi_primary if matrix[g][c] >= MIN_MEDIAN_COUNT) >= 4]
    profiles = {}
    standardised = {}
    for gid in omnibus_genes:
        vals = {c: cpm(gid, c) for c in epi_primary}
        mu, sd = mean(list(vals.values())), stdev(list(vals.values()))
        standardised[gid] = {c: ((v - mu) / sd if sd else 0.0) for c, v in vals.items()}
    for c in epi_primary:
        profiles[c] = [standardised[g][c] for g in omnibus_genes]

    def corr_distance(a, b):
        return 1 - pearson(profiles[a], profiles[b])

    def euclid_distance(a, b):
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(profiles[a], profiles[b])))

    def permute(statistic_fn):
        observed = statistic_fn([r["counts_column"] for r in epi["case"]])
        values = [statistic_fn(list(pick)) for pick in itertools.combinations(epi_primary, 4)]
        at_or_above = sum(1 for v in values if v >= observed - 1e-12)
        p = at_or_above / len(values)
        return {
            "observed_statistic": round(observed, 6),
            "assignments": len(values),
            "p_one_sided": round(p, 6),
            "null_max": round(max(values), 6),
            "null_median": round(median(values), 6),
            "separates_at_declared_alpha": p <= ALPHA + 1e-12,
        }

    corr_pairs = {(a, b): corr_distance(a, b) for i, a in enumerate(epi_primary)
                  for b in epi_primary[i + 1:]}

    def pair_distance(a, b):
        return corr_pairs[(a, b)] if (a, b) in corr_pairs else corr_pairs[(b, a)]

    def centroid_statistic(group_a):
        group_b = [c for c in epi_primary if c not in set(group_a)]
        ca = [mean([standardised[g][c] for c in group_a]) for g in omnibus_genes]
        cb = [mean([standardised[g][c] for c in group_b]) for g in omnibus_genes]
        return math.sqrt(sum((x - y) ** 2 for x, y in zip(ca, cb)))

    def between_minus_within(group_a):
        group_b = [c for c in epi_primary if c not in set(group_a)]
        between = [pair_distance(a, b) for a in group_a for b in group_b]
        within = ([pair_distance(a, b) for i, a in enumerate(group_a) for b in group_a[i + 1:]]
                  + [pair_distance(a, b) for i, a in enumerate(group_b) for b in group_b[i + 1:]])
        return mean(between) - mean(within)

    def dispersion_ratio(group_a):
        group_b = [c for c in epi_primary if c not in set(group_a)]
        wa = mean([pair_distance(a, b) for i, a in enumerate(group_a) for b in group_a[i + 1:]])
        wb = mean([pair_distance(a, b) for i, a in enumerate(group_b) for b in group_b[i + 1:]])
        lo, hi = min(wa, wb), max(wa, wb)
        return hi / lo if lo > 0 else float("inf")

    observed_case = [r["counts_column"] for r in epi["case"]]
    dispersion = dispersion_ratio(observed_case)
    dispersion_threshold = (freeze["endpoints"]["primary_2_omnibus_epithelial_divergence"]
                            ["dispersion_diagnostic_required_first"]["declared_threshold"])
    dispersion_breached = dispersion >= dispersion_threshold
    centroid = permute(centroid_statistic)
    centroid["readable"] = not dispersion_breached
    centroid["separates_and_readable"] = (centroid["separates_at_declared_alpha"]
                                          and not dispersion_breached)
    results["primary2_omnibus"] = {
        "genes_used": len(omnibus_genes),
        "itgb6_excluded": itgb6_id is not None,
        "dispersion_diagnostic": {
            "within_arm_distance_ratio": round(dispersion, 4),
            "declared_threshold": dispersion_threshold,
            "breached": dispersion_breached,
            "consequence_if_breached": "the omnibus is unreadable as a location change and "
                                       "carries no decision",
        },
        "centroid_distance": centroid,
        "between_minus_within_demoted": {
            **permute(between_minus_within),
            "status": "reported, carries no decision",
            "why": "simulation in this repository gives it a rejection rate of 1.000 under a "
                   "pure dispersion difference with no mean shift, against 0.035 under the "
                   "true null",
        },
    }

    # Engagement control, on the paired input.
    eng_scores, eng_cov = composite(engagement_genes, inp_contrast, True)
    eng_unstd, _ = composite(engagement_genes, inp_contrast, False)
    engagement = {"coverage": eng_cov}
    if eng_scores:
        engagement["standardised"] = contrast(
            [eng_scores[r["counts_column"]] for r in inp["case"]],
            [eng_scores[r["counts_column"]] for r in inp["reference"]],
            "engagement control, whole-lung input, standardised")
        engagement["unstandardised"] = contrast(
            [eng_unstd[r["counts_column"]] for r in inp["case"]],
            [eng_unstd[r["counts_column"]] for r in inp["reference"]],
            "engagement control, whole-lung input, unstandardised")
        eng_inj, _ = composite(engagement_genes, inp_injury, True)
        if eng_inj:
            engagement["injury_context"] = contrast(
                [eng_inj[r["counts_column"]] for r in inp["reference"]],
                [eng_inj[r["counts_column"]] for r in inp["saline"]],
                "engagement control injury context")
        d = engagement["standardised"]["mean_difference"]
        engagement["declared_direction"] = "lower"
        engagement["direction_observed"] = "lower" if d < 0 else "higher"
        engagement["separates_in_declared_direction"] = (
            engagement["standardised"]["separates_at_declared_alpha"] and d < 0)
    results["engagement_control"] = engagement

    # Purity covariates: immunoprecipitation minus input, per mouse.
    mouse_pairs = []
    for arm in ("case", "reference", "saline"):
        for e, i in zip(epi[arm], inp[arm]):
            mouse_pairs.append((arm, e["mouse"], e["counts_column"], i["counts_column"]))
    purity_rows = []
    purity_summary = {}
    for group, genes in (("epithelial_enrichment", epi_markers),
                         ("non_epithelial_de_enrichment", non_epi_markers),
                         ("macrophage_guard", macrophage)):
        per_marker = {}
        for symbol in genes:
            gid = sym2id.get(symbol)
            if not gid or gid not in matrix:
                per_marker[symbol] = {"status": "not recovered",
                                      "exceeds_declared_threshold": None,
                                      "note": "not computed; this is not a negative result"}
                continue
            values = {}
            for arm, mouse, ecol, icol in mouse_pairs:
                enrichment = cpm(gid, ecol) - cpm(gid, icol)
                values[(arm, mouse)] = enrichment
                purity_rows.append([group, symbol, arm, mouse, round(enrichment, 4)])
            case_v = [values[("case", m)] for m in sorted(case_mice)]
            ref_v = [values[("reference", m)] for m in sorted(ref_mice)]
            per_marker[symbol] = {
                "case_mean": round(mean(case_v), 4),
                "reference_mean": round(mean(ref_v), 4),
                "arm_difference": round(mean(case_v) - mean(ref_v), 4),
                "exceeds_declared_threshold": abs(mean(case_v) - mean(ref_v)) >= PURITY_THRESHOLD,
                "exact_test": exact_rank_sum(case_v, ref_v),
            }
        purity_summary[group] = per_marker
    results["purity_covariates"] = purity_summary
    purity_breach = sorted(
        f"{g}:{s}" for g, marks in purity_summary.items() for s, m in marks.items()
        if isinstance(m, dict) and m.get("exceeds_declared_threshold"))
    purity_uncomputed = sorted(
        f"{g}:{s}" for g, marks in purity_summary.items() for s, m in marks.items()
        if isinstance(m, dict) and m.get("status") == "not recovered")

    # A0 handling covariate, on the epithelium.
    a0_scores, a0_cov = composite(a0_genes, epi_contrast, True)
    a0_block = {"coverage": a0_cov, "role": "handling and stress covariate, not an endpoint"}
    if a0_scores:
        a0_block["contrast"] = contrast(
            [a0_scores[r["counts_column"]] for r in epi["case"]],
            [a0_scores[r["counts_column"]] for r in epi["reference"]],
            "A0 handling covariate")
    results["a0_handling_covariate"] = a0_block

    # Secondary, descriptive only.
    sec_scores, inj_cov = composite(injury_residual, epi_contrast, True)
    secondary = {"coverage": inj_cov, "status": "descriptive, no decision weight"}
    if sec_scores:
        secondary["contrast"] = contrast(
            [sec_scores[r["counts_column"]] for r in epi["case"]],
            [sec_scores[r["counts_column"]] for r in epi["reference"]],
            "A5/A11 injury residual, descriptive")
        if a0_scores:
            secondary["correlation_with_a0_score_vector"] = round(
                pearson([a0_scores[c] for c in epi_contrast],
                        [sec_scores[c] for c in epi_contrast]), 4)
    results["secondary_injury_residual"] = secondary

    # Sensitivities.
    def size_factors(columns):
        shared = [g for g in targeted if g in matrix and all(matrix[g][c] > 0 for c in columns)]
        if not shared:
            return {c: 1.0 for c in columns}
        log_means = {g: mean([math.log(matrix[g][c]) for c in columns]) for g in shared}
        return {c: math.exp(median([math.log(matrix[g][c]) - log_means[g] for g in shared]))
                for c in columns}

    sf = size_factors(epi_contrast)
    mor = {}
    for name, genes in panels.items():
        _, keep = resolve(genes, epi_contrast)
        if not keep:
            continue
        scores = {c: mean([math.log2(matrix[g][c] / (sf[c] * totals[c] / 1e6) + 1)
                           for _, g in keep]) for c in epi_contrast}
        mor[name] = contrast([scores[r["counts_column"]] for r in epi["case"]],
                             [scores[r["counts_column"]] for r in epi["reference"]],
                             f"median-of-ratios {name}")
    results["sensitivity_median_of_ratios"] = mor

    loo = []
    for name, genes in panels.items():
        scores, _ = composite(genes, epi_contrast, True)
        if not scores:
            continue
        for arm in ("case", "reference"):
            for r in epi[arm]:
                keep_case = [x["counts_column"] for x in epi["case"]
                             if x["counts_column"] != r["counts_column"]]
                keep_ref = [x["counts_column"] for x in epi["reference"]
                            if x["counts_column"] != r["counts_column"]]
                if not keep_case or not keep_ref:
                    continue
                loo.append({"panel": name, "dropped_mouse": r["mouse"], "dropped_arm": arm,
                            "mean_difference": round(mean([scores[c] for c in keep_case])
                                                     - mean([scores[c] for c in keep_ref]), 4)})
    results["sensitivity_leave_one_mouse_out"] = {
        "status": "descriptive, no decision weight, no robustness language", "rows": loo}

    depth = {}
    for name, genes in panels.items():
        scores, _ = composite(genes, epi_contrast, True)
        if scores:
            depth[name] = round(pearson([math.log2(totals[c]) for c in epi_primary],
                                        [scores[c] for c in epi_primary]), 4)
    results["sensitivity_score_against_log_depth"] = {
        "status": "descriptive, no decision weight", "pearson": depth}

    # The verdict, from the declared rules alone.
    def sep(block, variant="standardised"):
        b = block.get(variant) if isinstance(block, dict) else None
        return bool(b and b.get("separates_at_declared_alpha"))

    trans_sep = sep(results["primary1_transitional"])
    ident_sep = sep(results["primary1_identity"])
    omni_sep = results["primary2_omnibus"]["centroid_distance"]["separates_and_readable"]
    omni_euclid_sep = omni_sep
    eng_sep = bool(engagement.get("separates_in_declared_direction"))
    a0_sep = bool(a0_block.get("contrast", {}).get("separates_at_declared_alpha"))
    mor_disagree = any(
        mor.get(n, {}).get("separates_at_declared_alpha") != sep(results[f"primary1_{n}"])
        for n in panels if n in mor)
    unstd_disagree = any(
        sep(results[f"primary1_{n}"], "unstandardised") != sep(results[f"primary1_{n}"])
        for n in panels)
    refused = [k for k in ("primary1_transitional", "primary1_identity")
               if results[k].get("standardised", {}).get("coverage", {}).get("refused")]

    a0_stronger = False
    a0_comparison_meaningful = False
    if a0_block.get("contrast"):
        a0_p = a0_block["contrast"]["exact_test"]["p_two_sided"]
        panel_ps = [results[f"primary1_{n}"]["standardised"]["exact_test"]["p_two_sided"]
                    for n in panels if "standardised" in results[f"primary1_{n}"]
                    and "exact_test" in results[f"primary1_{n}"]["standardised"]]
        # "Separates more strongly" is only meaningful if something separates. When no
        # panel reaches the declared alpha, comparing which non-significant p-value is
        # smaller compares noise, so the comparison is recorded and not acted on. This is a
        # defect in the frozen rule's wording, reported rather than silently reinterpreted.
        a0_comparison_meaningful = bool(panel_ps) and (
            a0_p <= ALPHA + 1e-12 or min(panel_ps) <= ALPHA + 1e-12)
        a0_stronger = a0_comparison_meaningful and a0_p < min(panel_ps)

    downgrades = []
    if purity_breach:
        downgrades.append(f"purity or composition covariate exceeds the declared threshold: {purity_breach}")
    if purity_uncomputed:
        downgrades.append(f"a purity or composition covariate could not be computed, so its "
                          f"silence is not a negative result: {purity_uncomputed}")
    if a0_sep or a0_stronger:
        downgrades.append("the A0 handling covariate separates, or separates more strongly than the panels")
    if mor_disagree:
        downgrades.append("the median-of-ratios verdict differs from the CPM verdict")
    if unstd_disagree:
        downgrades.append("the unstandardised verdict differs from the standardised verdict")
    if results["primary2_omnibus"]["dispersion_diagnostic"]["breached"]:
        downgrades.append("the omnibus dispersion diagnostic is breached, so the omnibus "
                          "cannot be read as a location change")
    if trans_sep != ident_sep and (trans_sep or ident_sep):
        downgrades.append("the two primary-1 composites disagree")
    if refused:
        downgrades.append(f"a composite was refused under the measurability rule: {refused}")

    if downgrades:
        verdict = "inconclusive"
    elif (trans_sep or ident_sep or omni_sep):
        # Named so it cannot overclaim: both the A15 mechanism and rival 2 predict an
        # epithelial change, so a positive does not discriminate between them, and it does
        # not by itself establish target engagement either.
        verdict = "epithelium_differs_but_does_not_discriminate"
    elif eng_sep:
        verdict = "weak_bound_on_rival_2"
    else:
        verdict = "engagement_not_established"

    # Tables.
    with (OUT / "stage3_purity.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["group", "symbol", "arm", "mouse", "ip_minus_input_log2cpm"])
        w.writerows(purity_rows)
    with (OUT / "stage3_leave_one_out.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["panel", "dropped_mouse", "dropped_arm", "mean_difference"])
        for e in loo:
            w.writerow([e["panel"], e["dropped_mouse"], e["dropped_arm"], e["mean_difference"]])
    with (OUT / "stage3_per_gene.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["panel", "symbol", "ensembl_id", "measurable", "case_mean_log2cpm",
                    "reference_mean_log2cpm", "mean_difference", "p_two_sided"])
        for name, genes in list(panels.items()) + [("declared_here", declared_here),
                                                  ("engagement_control", engagement_genes)]:
            columns = inp_contrast if name == "engagement_control" else epi_contrast
            arms = inp if name == "engagement_control" else epi
            for symbol in genes:
                gid = sym2id.get(symbol)
                if not gid or gid not in matrix:
                    w.writerow([name, symbol, gid or "", False, "", "", "", ""])
                    continue
                vals = {c: cpm(gid, c) for c in columns}
                cv = [vals[r["counts_column"]] for r in arms["case"]]
                rv = [vals[r["counts_column"]] for r in arms["reference"]]
                w.writerow([name, symbol, gid, measurable(gid, columns),
                            round(mean(cv), 4), round(mean(rv), 4),
                            round(mean(cv) - mean(rv), 4),
                            exact_rank_sum(cv, rv)["p_two_sided"]])

    record = {
        "schema": "a15-rival2-stage3/v3",
        "stage": "3, execution under freeze v3",
        "scope": freeze["scope"],
        "freeze": {"path": "config/a15_rival2_freeze_v3.json", "sha256": sha256(FREEZE),
                   "status": freeze["status"]},
        "inputs": {"counts": {"sha256": counts_sha},
                   "ensembl_symbol_lookup": {"sha256": sha256(SYMBOL_MAP)},
                   "stage1_join": {"sha256": sha256(JOIN)},
                   "a1_config": {"sha256": sha256(A1_CONFIG)}},
        "declared_alpha": round(ALPHA, 6),
        "arms": {"epithelium": {k: [{"mouse": r["mouse"], "column": r["counts_column"],
                                     "gsm": r["gsm"]} for r in v] for k, v in epi.items()},
                 "input": {k: [{"mouse": r["mouse"], "column": r["counts_column"],
                                "gsm": r["gsm"]} for r in v] for k, v in inp.items()}},
        "library_totals": {c: totals[c] for c in sorted(libraries)},
        "results": results,
        "verdict": {
            "value": verdict,
            "only_informative_branch": "weak_bound_on_rival_2",
            "positive_branch_does_not_discriminate": (
                "both the A15 mechanism and rival 2 predict an epithelial change"),
            "transitional_separates": trans_sep,
            "identity_separates": ident_sep,
            "omnibus_correlation_separates": omni_sep,
            "omnibus_euclidean_separates": omni_euclid_sep,
            "engagement_separates_in_declared_direction": eng_sep,
            "downgrade_reasons": downgrades,
            "a0_rule_note": {
                "frozen_wording": "if the A0 handling covariate separates the arms more "
                                  "strongly than either primary-1 composite does, the "
                                  "primary reading is downgraded",
                "comparison_was_meaningful": a0_comparison_meaningful,
                "defect": "when nothing separates at the declared alpha, comparing which "
                          "non-significant p-value is smaller compares noise. The rule is "
                          "reported as written and its vacuous case is recorded rather "
                          "than reinterpreted after the fact.",
            },
            "purity_uncomputed": purity_uncomputed,
            "rules": freeze["declared_reading"],
        },
        "endpoint_scored": True,
        "claim_rows_added": 0,
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    run_record.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")

    print(f"omnibus genes used: {len(omnibus_genes)}")
    for key in ("primary1_transitional", "primary1_identity"):
        b = results[key].get("standardised", {})
        if "mean_difference" in b:
            print(f"{key}: diff {b['mean_difference']} p {b['exact_test']['p_two_sided']} "
                  f"separates {b['separates_at_declared_alpha']} "
                  f"(measurable {b['coverage']['measurable']}/{b['coverage']['declared']})")
    o = results["primary2_omnibus"]
    print(f"omnibus centroid: stat {o['centroid_distance']['observed_statistic']} "
          f"p {o['centroid_distance']['p_one_sided']} "
          f"readable {o['centroid_distance']['readable']} "
          f"separates_and_readable {o['centroid_distance']['separates_and_readable']}")
    print(f"dispersion ratio {o['dispersion_diagnostic']['within_arm_distance_ratio']} "
          f"threshold {o['dispersion_diagnostic']['declared_threshold']} "
          f"breached {o['dispersion_diagnostic']['breached']}")
    print(f"between-minus-within (demoted): p "
          f"{o['between_minus_within_demoted']['p_one_sided']}")
    if "standardised" in engagement:
        e = engagement["standardised"]
        print(f"engagement: diff {e['mean_difference']} p {e['exact_test']['p_two_sided']} "
              f"direction {engagement['direction_observed']} "
              f"separates_in_direction {engagement['separates_in_declared_direction']}")
    print(f"purity breaches: {purity_breach or 'none'}; "
          f"uncomputed: {purity_uncomputed or 'none'}")
    if a0_block.get("contrast"):
        print(f"A0 handling covariate: diff {a0_block['contrast']['mean_difference']} "
              f"p {a0_block['contrast']['exact_test']['p_two_sided']}")
    print(f"downgrades: {downgrades or 'none'}")
    print(f"VERDICT: {verdict}")
    print(f"wrote {run_record}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
