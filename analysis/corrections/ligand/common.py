"""Shared, side-effect-free definitions for the September 2026 ligand correction."""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import gammaln
from scipy.stats import wilcoxon

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
TARGETS = {"Fibroblast", "Myofibroblast"}
RESOURCES = ["cellchatdb", "cellphonedb", "consensus", "connectomedb2020", "italk"]
LIGANDS = ["AREG", "HBEGF", "TGFA", "EREG", "BTC", "EGF", "EPGN"]


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def eligible_donors(obs, floor=50):
    """Count the promised populations, never all stromal cells as targets."""
    work = obs.copy()
    work["is_sender"] = work["category"].eq("Epithelial")
    work["is_target"] = work["celltype"].isin(TARGETS)
    counts = work.groupby("donor", observed=True)[["is_sender", "is_target"]].sum()
    counts.columns = ["n_epithelial", "n_fibroblast"]
    counts["eligible"] = counts.min(axis=1) >= floor
    return counts


def allowed_pairs(obs):
    sources = sorted(obs.loc[obs["category"].eq("Epithelial"), "celltype"].unique())
    targets = sorted(set(obs["celltype"]) & TARGETS)
    return pd.DataFrame(itertools.product(sources, targets), columns=["source", "target"])


def abundance_guard(ligand, receptor):
    """Historical C14 rule, retained solely as a descriptive flag."""
    ligand, receptor = str(ligand).upper(), str(receptor).upper()
    return (ligand.startswith(("COL", "LAM")) or ligand == "FN1"
            or any(s in receptor for s in ("ITG", "SDC")) or receptor == "CD44")


def expected_detection_at_budget(total, target, budget):
    """Exact expected detection after sampling B molecules without replacement.

    No simulated cells become independent observations. Values remain per-cell
    expectations and are averaged within donor before a donor-paired test.
    """
    total = np.asarray(total, dtype=np.int64)
    target = np.asarray(target, dtype=np.int64)
    if budget < 1 or np.any(target < 0) or np.any(total < target):
        raise ValueError("invalid counts or molecule budget")
    answer = np.full(total.shape, np.nan)
    valid = total >= budget
    answer[valid & (target == 0)] = 0.0
    certain = valid & (total - target < budget)
    answer[certain] = 1.0
    m = valid & (target > 0) & ~certain
    n, k = total[m], target[m]
    lp = (gammaln(n - k + 1) + gammaln(n - budget + 1)
          - gammaln(n - k - budget + 1) - gammaln(n + 1))
    answer[m] = -np.expm1(np.minimum(lp, 0))
    return answer


def paired_summary(table, column):
    wide = table.pivot(index="donor", columns="compartment", values=column)
    paired = wide.reindex(columns=["epithelial", "myeloid"]).dropna()
    diff = paired["epithelial"] - paired["myeloid"]
    return {"n_donors": len(paired),
            "median_epithelial": float(paired.epithelial.median()) if len(paired) else None,
            "median_myeloid": float(paired.myeloid.median()) if len(paired) else None,
            "median_paired_difference": float(diff.median()) if len(paired) else None,
            "epithelial_higher": int((diff > 0).sum()),
            "p_two_sided": float(wilcoxon(diff, method="auto").pvalue) if len(paired) >= 5 else None,
            "donors": paired.index.tolist()}
