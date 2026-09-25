"""Shared helpers for the GSEA trials (G1, G2): per-unit pseudobulk, the ranking,
the gene sets and the two nulls.

Everything here follows the unit rule. A pseudobulk is the sum of raw counts over
the cells of one compartment in one animal or donor; the ranking is a statistic
between arms of those units; the gene sets are read from files whose hash is in
the run record; the engine is gseapy's pre-ranked GSEA (gene-label permutation
null), cross-checked by an in-house enrichment score and supplemented by an
expression-matched random-set null of the kind trial M2 used, because a random
set of the same size drawn without regard to expression level is not a fair
comparator for a set of highly expressed genes (the D2 lesson).
"""
from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
MSIGDB = REPO / "raw_data" / "msigdb"


# ---------------------------------------------------------------------------
# gene sets
# ---------------------------------------------------------------------------
def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_gmt(path: Path) -> dict[str, set[str]]:
    sets = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = set(parts[2:])
    return sets


def gene_set_facts(path: Path) -> dict:
    sets = read_gmt(path)
    return {"file": path.relative_to(REPO).as_posix(), "sha256": sha256(path), "n_sets": len(sets),
            "bytes": path.stat().st_size}


# ---------------------------------------------------------------------------
# pseudobulk
# ---------------------------------------------------------------------------
def pseudobulk_from_h5ad_counts(path: Path, group_of_cell: np.ndarray, n_groups: int,
                                block: int = 4000) -> np.ndarray:
    """Sum the raw `layers/counts` (CSR, cells by genes) of an h5ad file into
    n_groups rows. `group_of_cell` holds the group index per cell, -1 to skip.
    Rows are read in contiguous blocks so that the file is streamed once."""
    with h5py.File(path, "r") as f:
        L = f["layers"]["counts"]
        assert dict(L.attrs).get("encoding-type") == "csr_matrix", "counts layer is not CSR"
        n_cells, n_genes = (int(x) for x in L.attrs["shape"])
        indptr = L["indptr"][:]
        out = np.zeros((n_groups, n_genes), dtype=np.int64)
        for start in range(0, n_cells, block):
            stop = min(start + block, n_cells)
            g = group_of_cell[start:stop]
            if not (g >= 0).any():
                continue
            lo, hi = int(indptr[start]), int(indptr[stop])
            data = L["data"][lo:hi]
            idx = L["indices"][lo:hi]
            offs = indptr[start:stop + 1] - lo
            for j in np.flatnonzero(g >= 0):
                a, b = int(offs[j]), int(offs[j + 1])
                np.add.at(out[g[j]], idx[a:b], data[a:b])
    return out


def pseudobulk_from_mtx(path: Path, group_of_cell: np.ndarray, n_groups: int,
                        chunk_rows: int = 10_000_000) -> tuple[np.ndarray, int]:
    """Sum a MatrixMarket coordinate file (genes by cells, 1-based) into
    n_groups columns, streaming in chunks. Returns (n_genes by n_groups, nnz read)."""
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        assert header.startswith("%%MatrixMarket matrix coordinate"), header
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_genes, n_cells, nnz = (int(x) for x in line.split())
    assert len(group_of_cell) == n_cells, (len(group_of_cell), n_cells)
    out = np.zeros((n_genes, n_groups), dtype=np.int64)
    seen = 0
    reader = pd.read_csv(path, sep=" ", comment="%", header=None, names=["g", "c", "v"],
                         skiprows=2, chunksize=chunk_rows, dtype={"g": np.int32, "c": np.int32, "v": np.int64},
                         engine="c")
    for chunk in reader:
        seen += len(chunk)
        grp = group_of_cell[chunk["c"].to_numpy() - 1]
        keep = grp >= 0
        if not keep.any():
            continue
        flat = (chunk["g"].to_numpy()[keep] - 1).astype(np.int64) * n_groups + grp[keep]
        out += np.bincount(flat, weights=chunk["v"].to_numpy()[keep], minlength=n_genes * n_groups).reshape(n_genes, n_groups).astype(np.int64)
    assert seen == nnz, f"read {seen} entries, header says {nnz}"
    return out, seen


# ---------------------------------------------------------------------------
# ranking
# ---------------------------------------------------------------------------
def log_cpm(counts: np.ndarray) -> np.ndarray:
    """counts: units by genes. Library size is the unit's total over the genes kept."""
    lib = counts.sum(1, keepdims=True).astype(np.float64)
    return np.log2(1e6 * counts / np.maximum(lib, 1) + 1.0)


def welch_t(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Per-gene Welch t between two unit groups (rows of a and b are units)."""
    ma, mb = a.mean(0), b.mean(0)
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    se = np.sqrt(va / a.shape[0] + vb / b.shape[0])
    se[se == 0] = np.nan
    t = (ma - mb) / se
    return np.nan_to_num(t, nan=0.0)


def ranking_between_arms(counts: np.ndarray, genes: np.ndarray, arm_a: np.ndarray, arm_b: np.ndarray,
                         min_count: int = 10, min_units: int = 3) -> tuple[pd.Series, dict]:
    """Genes with at least `min_count` counts in at least `min_units` units of
    the two arms together; log2 CPM; Welch t of arm_a minus arm_b; ties broken
    by the mean difference, then by name, so the ranking is deterministic."""
    sub = counts[np.concatenate([arm_a, arm_b])]
    expressed = (sub >= min_count).sum(0) >= min_units
    lc = log_cpm(counts[:, expressed])
    ia = np.arange(len(arm_a)); ib = np.arange(len(arm_a), len(arm_a) + len(arm_b))
    lc_sub = lc[np.concatenate([arm_a, arm_b])]
    t = welch_t(lc_sub[ia], lc_sub[ib])
    diff = lc_sub[ia].mean(0) - lc_sub[ib].mean(0)
    frame = pd.DataFrame({"gene": genes[expressed], "t": t, "diff": diff, "mean_logcpm": lc_sub.mean(0)})
    frame = frame.sort_values(["t", "diff", "gene"], ascending=[False, False, True]).reset_index(drop=True)
    return frame, {"genes_tested": int(expressed.sum()), "genes_total": int(len(genes)),
                   "filter": f"at least {min_count} counts in at least {min_units} units of the two arms"}


# ---------------------------------------------------------------------------
# enrichment score and the nulls
# ---------------------------------------------------------------------------
def enrichment_score(scores: np.ndarray, hit: np.ndarray, weight: float = 1.0) -> float:
    """Classic weighted Kolmogorov-Smirnov running-sum statistic on a ranking
    sorted from the most positive score down. `hit` is a boolean mask."""
    n = len(scores)
    nh = int(hit.sum())
    if nh == 0 or nh == n:
        return 0.0
    w = np.abs(scores) ** weight
    w_hit = np.where(hit, w, 0.0)
    total = w_hit.sum()
    if total == 0:
        return 0.0
    step_hit = w_hit / total
    step_miss = np.where(hit, 0.0, 1.0 / (n - nh))
    run = np.cumsum(step_hit - step_miss)
    i = int(np.argmax(np.abs(run)))
    return float(run[i])


def matched_random_null(scores: np.ndarray, mean_expr: np.ndarray, hit: np.ndarray,
                        n_draws: int, rng: np.random.Generator, n_bins: int = 20) -> np.ndarray:
    """Enrichment scores of random sets that match the set's size and its
    expression-level profile: genes are binned by mean expression and each
    draw takes, from every bin, as many genes as the set has there."""
    bins = pd.qcut(pd.Series(mean_expr), n_bins, labels=False, duplicates="drop").to_numpy()
    per_bin = {b: np.flatnonzero(bins == b) for b in np.unique(bins)}
    need = {b: int(((bins == b) & hit).sum()) for b in per_bin}
    out = np.empty(n_draws)
    for d in range(n_draws):
        mask = np.zeros(len(scores), dtype=bool)
        for b, k in need.items():
            if k:
                mask[rng.choice(per_bin[b], size=k, replace=False)] = True
        out[d] = enrichment_score(scores, mask)
    return out


def matched_null_summary(observed: float, null: np.ndarray) -> dict:
    same_sign = null[np.sign(null) == np.sign(observed)] if observed != 0 else null
    p = float((np.abs(null) >= abs(observed)).mean())
    nes = float(observed / np.mean(np.abs(same_sign))) if len(same_sign) and np.mean(np.abs(same_sign)) > 0 else float("nan")
    return {"matched_p": p, "matched_nes": nes, "matched_draws": int(len(null))}


def run_prerank(ranking: pd.DataFrame, sets: dict[str, set[str]], seed: int, n_perm: int = 1000,
                min_size: int = 15, max_size: int = 500, threads: int = 8) -> pd.DataFrame:
    """gseapy pre-ranked GSEA on the t ranking; returns one row per set kept."""
    import gseapy as gp
    rnk = ranking.set_index("gene")["t"]
    res = gp.prerank(rnk=rnk, gene_sets={k: sorted(v) for k, v in sets.items()}, permutation_num=n_perm,
                     min_size=min_size, max_size=max_size, seed=seed, threads=threads, outdir=None,
                     weight=1.0, verbose=False)
    r = res.res2d.copy()
    r = r.rename(columns={"Term": "set", "ES": "es", "NES": "nes", "NOM p-val": "p_perm", "FDR q-val": "fdr",
                          "Lead_genes": "leading_edge"})
    keep = ["set", "es", "nes", "p_perm", "fdr", "leading_edge"]
    r["size"] = r["Tag %"].astype(str).str.split("/").str[1].astype(int) if "Tag %" in r else np.nan
    return r[keep + ["size"]].reset_index(drop=True)
