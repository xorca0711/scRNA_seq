"""Shared, deterministic estimators for the staged A0 exploratory pilot."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
import pandas as pd
from scipy import sparse

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]
TABLES = BASE / "tables/exploratory"
PROC = BASE / "processed/exploratory"
FIG = BASE / "figures/exploratory"
for folder in (TABLES, PROC, FIG):
    folder.mkdir(parents=True, exist_ok=True)


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def excluded(gene, config):
    return gene in config["label_genes"] or bool(re.match(r"^(mt-|Rpl|Rps)", gene, flags=re.I))


def rank_scores(matrix, genes, module, references):
    """Mean probability of beating a fixed reference gene within each cell.

    Ties contribute 0.5; all-zero cells therefore score exactly 0.5.
    Module/reference overlaps are removed from the reference set.
    """
    lookup = {g: i for i, g in enumerate(genes)}
    query = [g for g in module if g in lookup]
    reference = [g for g in references if g in lookup and g not in set(module)]
    if not query or not reference:
        raise ValueError("Empty module or reference after feature alignment")
    to_dense = lambda x: x.toarray() if sparse.issparse(x) else np.asarray(x)
    bg = np.sort(to_dense(matrix[:, [lookup[g] for g in reference]]), axis=1)
    values = to_dense(matrix[:, [lookup[g] for g in query]])
    result = np.empty(matrix.shape[0], dtype=np.float64)
    for i, (r, q) in enumerate(zip(bg, values)):
        result[i] = (np.searchsorted(r, q, side="left") + np.searchsorted(r, q, side="right")).mean() / (2 * len(reference))
    assert np.isfinite(result).all() and ((result >= 0) & (result <= 1)).all()
    return result, {"requested_genes": len(module), "present_genes": len(query),
        "feature_coverage": len(query) / len(module), "reference_coverage": sum(g in lookup for g in references) / len(references),
        "reference_genes_used": len(reference)}


def paired_effects(obs, scores, start="AT2", intermediate="Krt8+ ADI", end="AT1", dataset="repair"):
    obs = obs.copy()
    obs["score"] = scores
    rows = []
    for unit, frame in obs.groupby("unit", sort=True):
        arrays = {s: frame.loc[frame.state == s, "score"].to_numpy() for s in [start, intermediate, end]}
        if any(len(v) == 0 for v in arrays.values()):
            continue
        a, d = arrays[start], arrays[end]
        variance = ((len(a)-1)*np.var(a, ddof=1) + (len(d)-1)*np.var(d, ddof=1)) / (len(a)+len(d)-2) if min(len(a),len(d)) > 1 else np.nan
        scale = np.sqrt(variance)
        for endpoint in [start, end]:
            delta = float(arrays[intermediate].mean() - arrays[endpoint].mean())
            rows.append({"dataset": dataset, "unit": unit, "time": frame.time.iloc[0], "endpoint": endpoint,
                "n_start": len(a), "n_intermediate": len(arrays[intermediate]), "n_end": len(d),
                "all_states_at_least_30": min(map(len, arrays.values())) >= 30,
                "intermediate_mean": float(arrays[intermediate].mean()), "endpoint_mean": float(arrays[endpoint].mean()),
                "difference": delta, "endpoint_pooled_cell_sd": scale,
                "standardized_difference": delta / scale if scale > 0 else np.nan})
    return pd.DataFrame(rows)


def pseudobulk(matrix, obs, genes, states):
    rows, values, detection = [], [], []
    for (unit, state), indexes in obs[obs.state.isin(states)].groupby(["unit", "state"]).groups.items():
        part = matrix[np.asarray(list(indexes))]
        counts = np.asarray(part.sum(axis=0)).ravel()
        assert counts.sum() > 0
        rows.append({"unit": unit, "state": state, "cells": part.shape[0], "total_counts": float(counts.sum())})
        values.append(np.log2(counts * 1e6 / counts.sum() + 1))
        detection.append(np.asarray((part > 0).mean(axis=0)).ravel())
    return pd.DataFrame(rows), np.asarray(values), np.asarray(detection)


def learn_program(facts, logcpm, detection, genes, config, omit=None):
    units = sorted(set(facts.unit) - ({omit} if omit else set()))
    index = {(r.unit, r.state): i for i, r in facts.iterrows()}
    diffs = {e: np.asarray([logcpm[index[(u,"Krt8+ ADI")]] - logcpm[index[(u,e)]] for u in units]) for e in ["AT2","AT1"]}
    det = np.asarray([detection[index[(u,"Krt8+ ADI")]] for u in units])
    table = pd.DataFrame({"gene": genes, "excluded": [excluded(g,config) for g in genes],
        "adi_detection_consistency": (det >= config["minimum_adi_cell_detection"]).mean(axis=0)})
    for e in diffs:
        table[f"median_log2cpm_difference_vs_{e}"] = np.median(diffs[e], axis=0)
        table[f"positive_fraction_vs_{e}"] = (diffs[e] > 0).mean(axis=0)
    table["minimum_median_effect"] = table[[f"median_log2cpm_difference_vs_{e}" for e in diffs]].min(axis=1)
    table["minimum_positive_fraction"] = table[[f"positive_fraction_vs_{e}" for e in diffs]].min(axis=1)
    table["eligible"] = (~table.excluded & (table.adi_detection_consistency >= 2/3)
        & (table.minimum_positive_fraction >= 2/3) & (table.minimum_median_effect >= config["minimum_log2cpm_difference"]))
    selected = table[table.eligible].sort_values(["minimum_median_effect","minimum_positive_fraction","gene"], ascending=[False,False,True]).head(config["maximum_program_genes"]).gene.tolist()
    table["selected"] = table.gene.isin(selected)
    return selected, table


def load_repair():
    return (sparse.load_npz(PROC / "repair_counts.npz"), pd.read_csv(PROC / "repair_obs.csv"),
            json.loads((PROC / "repair_genes.json").read_text()),
            json.loads((BASE / "exploratory_config.json").read_text()))


def resolve_source_path(name):
    """Resolve frozen Windows/pre-migration paths without altering their hashes."""
    normalized = name.replace('\\', '/')
    path = ROOT / normalized
    if not path.exists() and normalized.startswith('Thesis/'):
        path = ROOT / normalized.replace('Thesis/', 'Research Article/', 1)
    return path
