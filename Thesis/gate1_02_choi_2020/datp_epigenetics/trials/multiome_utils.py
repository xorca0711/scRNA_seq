"""Shared helpers for the DATP epigenetics branch of roadmap paper 2.

Everything generic (the run record, the SOFT parser, package versions, the
markdown table writer) comes from the repository's one shared helper module,
the same one the Cardoso and Choi trials use. Nothing is redefined here.

This module adds only what is specific to 10x multiome data, which no earlier
trial in this repository has touched:

* Where the two multiome deposits live and what each barcode suffix means.
* A single streaming pass over a 10x `filtered_feature_bc_matrix.h5` that
  never holds the matrix densely. These files carry 126 to 237 million
  non-zero entries across 180,000 to 241,000 features, so a dense load is not
  possible on this machine and a full sparse load is wasteful when the frozen
  gene set is fifteen genes wide.
* The peak-to-gene join, which reads the vendor's own
  `atac_peak_annotation.tsv.gz` rather than re-annotating peaks here.

THE ONE STRUCTURAL FACT THAT SHAPES EVERY TRIAL IN THIS FOLDER. Both deposits
carry one library per experimental condition. GSE310539 aggregates four
libraries (two genotypes by two treatments) into one file, distinguished by
the barcode suffix; GSE247130 deposits three files of two libraries each
(mutant against control at three stages). No condition is repeated, so no
between-condition contrast in either deposit has within-group replication and
none is testable. Only comparisons between cell states INSIDE one library are
admissible, and the statistical unit for those is the library, with the
reading required to hold in every library rather than pooled across them.

SUFFIX ORDER IS AN ASSUMPTION, NOT A FACT. Neither deposit includes the
aggregation CSV that cellranger-arc used, so the map from barcode suffix to
library below is the GEO sample order, which is the cellranger-arc default but
is not guaranteed. Every trial that uses these names must run
`corroborate_suffix_map` and report the outcome; where the corroboration
fails, the libraries are reported as L1, L2 and so on and no condition name is
attached to any number.
"""

from __future__ import annotations

import gzip
import sys
from collections import defaultdict
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[4]
RAW = REPO / "raw_data"
HERE = Path(__file__).resolve().parent
PALETTE = REPO / "analysis" / "config" / "palette.json"

_SHARED_TRIALS = REPO / "Thesis" / "gate1_04_sikkema_2023_hlca" / "trials"
if str(_SHARED_TRIALS) not in sys.path:
    sys.path.insert(0, str(_SHARED_TRIALS))
from trial_utils import (RunRecord, df_to_markdown, file_facts,  # noqa: E402,F401
                         package_versions, parse_soft, utc_now)

# ---------------------------------------------------------------------------
# The deposits
# ---------------------------------------------------------------------------

# GSE310539: Lynch et al. 2026, AP-1 and the injury-induced transitional state.
# One aggregate file, four libraries by barcode suffix.
GSE310539 = {
    "accession": "GSE310539",
    "paper": "Lynch et al. 2026, Am J Respir Cell Mol Biol, doi:10.1093/ajrcmb/aanag157",
    "matrix": RAW / "GSE310539" / "GSE310539_totalaggr_filtered_feature_bc_matrix.h5",
    "peaks": RAW / "GSE310539" / "GSE310539_totalaggr_atac_peak_annotation.tsv.gz",
    "libraries": {
        "1": {"name": "wildtype_PBS", "genotype": "wildtype", "treatment": "PBS"},
        "2": {"name": "wildtype_SeV", "genotype": "wildtype", "treatment": "SeV"},
        "3": {"name": "AP1mut_PBS", "genotype": "AP-1 mutant", "treatment": "PBS"},
        "4": {"name": "AP1mut_SeV", "genotype": "AP-1 mutant", "treatment": "SeV"},
    },
}

# GSE247130: Hassan and Chen 2024, CEBPA restricts AT2 plasticity.
# Three files, two libraries each. Note that the deposit also carries a second
# seven-week matrix under the near-identical name "Aggregate7wk" (no hyphen);
# the hyphenated file is the one with a matching peak annotation and is the one
# used here. The duplicate is recorded as a deposit quirk, not resolved.
GSE247130 = {
    "accession": "GSE247130",
    "paper": "Hassan and Chen 2024, Nat Commun 15:4148, doi:10.1038/s41467-024-48632-3",
    "files": [
        {
            "stage": "P9",
            "matrix": RAW / "GSE247130" / "GSE247130_Aggregate.P9_filtered_feature_bc_matrix.h5",
            "peaks": RAW / "GSE247130" / "GSE247130_Aggregate.P9_atac_peak_annotation.tsv.gz",
            "libraries": {
                "1": {"name": "P9_mutant", "genotype": "Cebpa mutant", "stage": "P9"},
                "2": {"name": "P9_control", "genotype": "control", "stage": "P9"},
            },
        },
        {
            "stage": "7wk",
            "matrix": RAW / "GSE247130" / "GSE247130_Aggregate7-wk_filtered_feature_bc_matrix.h5",
            "peaks": RAW / "GSE247130" / "GSE247130_Aggregate7-wk_atac_peak_annotation.tsv.gz",
            "libraries": {
                "1": {"name": "7wk_mutant", "genotype": "Cebpa mutant", "stage": "7 weeks"},
                "2": {"name": "7wk_control", "genotype": "control", "stage": "7 weeks"},
            },
        },
        {
            "stage": "SeV",
            "matrix": RAW / "GSE247130" / "GSE247130_Aggregate_SeV_filtered_feature_bc_matrix.h5",
            "peaks": RAW / "GSE247130" / "GSE247130_Aggregate_SeV_atac_peak_annotation.tsv.gz",
            "libraries": {
                "1": {"name": "SeV_mutant", "genotype": "Cebpa mutant", "stage": "SeV infected"},
                "2": {"name": "SeV_control", "genotype": "control", "stage": "SeV infected"},
            },
        },
    ],
}

RNA = "Gene Expression"
ATAC = "Peaks"


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------

def read_features(path: Path) -> pd.DataFrame:
    """Feature table only: identifier, symbol, type and genomic interval."""
    with h5py.File(path, "r") as h:
        f = h["matrix"]["features"]
        out = pd.DataFrame({
            "id": [x.decode() for x in f["id"][:]],
            "name": [x.decode() for x in f["name"][:]],
            "feature_type": [x.decode() for x in f["feature_type"][:]],
            "interval": [x.decode() for x in f["interval"][:]],
        })
    out["row"] = np.arange(len(out))
    return out


def read_barcodes(path: Path) -> pd.DataFrame:
    with h5py.File(path, "r") as h:
        bc = [x.decode() for x in h["matrix"]["barcodes"][:]]
    return pd.DataFrame({"barcode": bc, "suffix": [b.split("-")[-1] for b in bc]})


def peak_gene_table(path: Path) -> pd.DataFrame:
    """The vendor's own peak annotation, as (interval, gene, distance, type).

    One peak can be annotated to several genes, so this table is many to many
    and must never be used as if it were a peak index.
    """
    rows = []
    with gzip.open(path, "rt") as fh:
        header = next(fh)
        assert header.startswith("chrom\tstart\tend\tgene"), header[:40]
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6 or not p[3]:
                continue
            rows.append((p[0] + ":" + p[1] + "-" + p[2], p[3], p[4], p[5].strip()))
    return pd.DataFrame(rows, columns=["interval", "gene", "distance", "peak_type"])


def stream_downsampled_detection(path: Path, cells: np.ndarray, target: int,
                                 seed: int = 0, chunk: int = 2000):
    """Read selected cells, downsample each to `target` in-peak fragments, binarise.

    Detection fraction is monotone and saturating in depth, and its derivative
    with respect to log depth is maximal exactly where most ATAC peaks sit, so
    a residual depth gap of a few tens of per cent moves a peak's detection
    fraction by more than any threshold worth setting. Matching on depth leaves
    a within-caliper residual; downsampling every cell to an identical fragment
    budget removes the gap by construction, which is the only treatment that
    makes two groups' detection fractions comparable at face value.

    Each cell is downsampled ONCE and the same binarised vector is reused for
    the real contrast and for every sham partition, so the sham band measures
    the pipeline rather than a second roll of the dice.

    Returns (detected, kept, n_peaks): a list of index arrays, one per cell in
    `cells` order (empty where the cell fell below `target`), a boolean array
    saying which cells survived, and the number of peak features.
    """
    cells = np.asarray(cells, dtype=np.int64)
    order = {int(c): i for i, c in enumerate(cells)}
    rng = np.random.default_rng(seed)

    with h5py.File(path, "r") as h:
        g = h["matrix"]
        n_feat, n_cells = (int(x) for x in g["shape"][:])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        atac_row = np.flatnonzero(ftype == ATAC)
        peak_slot = np.full(n_feat, -1, dtype=np.int64)
        peak_slot[atac_row] = np.arange(len(atac_row))

        keep_col = np.zeros(n_cells, dtype=bool)
        keep_col[cells] = True
        indptr = g["indptr"][:]
        detected = [np.empty(0, dtype=np.int32)] * len(cells)
        kept = np.zeros(len(cells), dtype=bool)

        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            if not keep_col[start:stop].any():
                continue
            lo, hi = int(indptr[start]), int(indptr[stop])
            idx = g["indices"][lo:hi]
            dat = g["data"][lo:hi]
            offs = indptr[start:stop + 1] - lo
            for j in range(stop - start):
                col = start + j
                if not keep_col[col]:
                    continue
                a, b = int(offs[j]), int(offs[j + 1])
                ii, dd = idx[a:b], dat[a:b]
                s = peak_slot[ii]
                m = s >= 0
                pk, cts = s[m].astype(np.int64), dd[m].astype(np.int64)
                total = int(cts.sum())
                slot_i = order[col]
                if total < target:
                    continue
                if total == target:
                    sub = cts
                else:
                    sub = rng.multivariate_hypergeometric(cts, target, method="marginals")
                hit = sub > 0
                detected[slot_i] = pk[hit].astype(np.int32)
                kept[slot_i] = True

    return detected, kept, len(atac_row)


def detection_fraction(detected: list, members: np.ndarray, n_peaks: int) -> np.ndarray:
    """Per-peak fraction of the member cells in which the peak was detected."""
    members = np.asarray(members, dtype=bool)
    idx = np.flatnonzero(members)
    if len(idx) == 0:
        return np.zeros(n_peaks)
    parts = [detected[i] for i in idx if len(detected[i])]
    if not parts:
        return np.zeros(n_peaks)
    counts = np.bincount(np.concatenate(parts), minlength=n_peaks)
    return counts / len(idx)


def detect_prob_at_depth(count: np.ndarray, total: np.ndarray, target: int) -> np.ndarray:
    """Probability a gene with `count` of `total` UMI is seen in `target` draws.

    Sampling `target` UMI without replacement from a cell of `total`, the gene
    is missed with probability C(total-count, target) / C(total, target). This
    is the exact expectation, so the RNA side needs no random draw and no seed:
    it is the analytic equivalent of downsampling every cell to one depth.
    Peaks are handled by sampling instead, because their expectation would have
    to be accumulated across two hundred thousand features per cell.

    Cells with fewer than `target` UMI return NaN and must be dropped by the
    caller, which is the same rule the peak side applies.
    """
    count = np.asarray(count, dtype=np.float64)
    total = np.asarray(total, dtype=np.float64)
    out = np.full(count.shape, np.nan)
    ok = (total >= target) & (total - count - target + 1 > 0)
    from scipy.special import gammaln
    lp = (gammaln(total[ok] - count[ok] + 1) - gammaln(total[ok] - count[ok] - target + 1)
          - gammaln(total[ok] + 1) + gammaln(total[ok] - target + 1))
    out[ok] = 1.0 - np.exp(lp)
    # a gene whose count leaves fewer than target UMI behind is certainly seen
    certain = (total >= target) & ~ok
    out[certain] = 1.0
    return out


def stream_peak_detection(path: Path, cells: np.ndarray, groups: dict,
                          chunk: int = 2000) -> tuple[np.ndarray, list]:
    """Second pass: per-peak detection counts for each named group of cells.

    `cells` is the subset of column indices to read at all, and `groups` maps a
    name to a boolean mask over that subset. Returns an array of shape
    (n_peaks, n_groups) holding, for each peak, the number of cells of that
    group in which the peak has a non-zero count, and the group order.

    The cost is one pass over the file and one bincount per group. Counting is
    done per group rather than per cell because a cell belongs to the real
    contrast and to several sham partitions at once, and incrementing per cell
    would repeat the same scatter twenty times over.
    """
    cells = np.asarray(cells, dtype=np.int64)
    keep = np.zeros(0, dtype=bool)
    order = {int(c): i for i, c in enumerate(cells)}

    with h5py.File(path, "r") as h:
        g = h["matrix"]
        n_feat, n_cells = (int(x) for x in g["shape"][:])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        atac_row = np.flatnonzero(ftype == ATAC)
        # map a global feature row to its position among peaks, -1 for genes
        peak_slot = np.full(n_feat, -1, dtype=np.int64)
        peak_slot[atac_row] = np.arange(len(atac_row))

        keep = np.zeros(n_cells, dtype=bool)
        keep[cells] = True
        indptr = g["indptr"][:]
        per_cell = [None] * len(cells)

        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            if not keep[start:stop].any():
                continue
            lo, hi = int(indptr[start]), int(indptr[stop])
            idx = g["indices"][lo:hi]
            offs = indptr[start:stop + 1] - lo
            for j in range(stop - start):
                col = start + j
                if not keep[col]:
                    continue
                ii = idx[int(offs[j]):int(offs[j + 1])]
                s = peak_slot[ii]
                per_cell[order[col]] = s[s >= 0].astype(np.int32)

    names = list(groups)
    out = np.zeros((len(atac_row), len(names)), dtype=np.int32)
    for gi, name in enumerate(names):
        mask = np.asarray(groups[name], dtype=bool)
        members = [per_cell[i] for i in np.flatnonzero(mask)]
        if not members:
            continue
        out[:, gi] = np.bincount(np.concatenate(members), minlength=len(atac_row))
    return out, names


def match_pairs(a_idx: np.ndarray, b_idx: np.ndarray, a_key: np.ndarray,
                b_key: np.ndarray, caliper: float) -> tuple[np.ndarray, np.ndarray]:
    """Greedy 1:1 matching on one log-scale key, nearest first, no replacement.

    `caliper` is in the units of the key. Both arrays are sorted once and
    walked, so the cost is n log n rather than n squared. Returns the matched
    subsets of a_idx and b_idx in pair order.
    """
    order_b = np.argsort(b_key)
    b_sorted_key = b_key[order_b]
    used = np.zeros(len(b_idx), dtype=bool)
    out_a, out_b = [], []
    for ai in np.argsort(a_key):
        pos = np.searchsorted(b_sorted_key, a_key[ai])
        best, best_d = -1, np.inf
        for p in range(max(0, pos - 40), min(len(b_sorted_key), pos + 40)):
            if used[order_b[p]]:
                continue
            d = abs(b_sorted_key[p] - a_key[ai])
            if d < best_d:
                best, best_d = p, d
        if best >= 0 and best_d <= caliper:
            used[order_b[best]] = True
            out_a.append(a_idx[ai])
            out_b.append(b_idx[order_b[best]])
    return np.array(out_a, dtype=np.int64), np.array(out_b, dtype=np.int64)


def stream_selected(path: Path, selected_rows: np.ndarray, chunk: int = 2000):
    """One pass over the cells, returning per-cell totals and selected features.

    The 10x h5 matrix is CSC over cells, so a chunk of columns is a contiguous
    slice of `data` and `indices` and costs one read each. Returns

      totals : DataFrame, one row per cell, with RNA and ATAC count sums and
               the number of distinct genes and peaks detected
      picked : dense array of shape (n_selected, n_cells) holding the raw
               counts of the selected feature rows only

    Nothing else is retained, so peak memory is one chunk of the matrix.
    """
    selected_rows = np.asarray(selected_rows, dtype=np.int64)
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        n_feat, n_cells = (int(x) for x in g["shape"][:])
        ftype = np.array([x.decode() for x in g["features"]["feature_type"][:]])
        is_rna = ftype == RNA
        is_atac = ftype == ATAC

        # position of each selected row inside the output block, -1 if not wanted
        slot = np.full(n_feat, -1, dtype=np.int64)
        slot[selected_rows] = np.arange(len(selected_rows))

        indptr = g["indptr"][:]
        picked = np.zeros((len(selected_rows), n_cells), dtype=np.int32)
        rna_sum = np.zeros(n_cells, dtype=np.int64)
        atac_sum = np.zeros(n_cells, dtype=np.int64)
        rna_det = np.zeros(n_cells, dtype=np.int32)
        atac_det = np.zeros(n_cells, dtype=np.int32)

        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            lo, hi = int(indptr[start]), int(indptr[stop])
            data = g["data"][lo:hi]
            idx = g["indices"][lo:hi]
            offs = indptr[start:stop + 1] - lo
            for j in range(stop - start):
                a, b = int(offs[j]), int(offs[j + 1])
                ii, dd = idx[a:b], data[a:b]
                r = is_rna[ii]
                t = is_atac[ii]
                col = start + j
                rna_sum[col] = dd[r].sum()
                atac_sum[col] = dd[t].sum()
                rna_det[col] = int(r.sum())
                atac_det[col] = int(t.sum())
                s = slot[ii]
                keep = s >= 0
                if keep.any():
                    picked[s[keep], col] = dd[keep]

    totals = pd.DataFrame({
        "rna_counts": rna_sum,
        "atac_counts": atac_sum,
        "genes_detected": rna_det,
        "peaks_detected": atac_det,
    })
    return totals, picked


# ---------------------------------------------------------------------------
# The suffix corroboration
# ---------------------------------------------------------------------------

def corroborate_suffix_map(labels: pd.DataFrame, expect: dict, column: str,
                           higher_in: str) -> dict:
    """Check the assumed suffix order against a fact the data can settle.

    `labels` carries one row per suffix with a measured quantity in `column`.
    The assumed map passes if the suffixes whose assigned condition matches
    `higher_in` all score above every suffix that does not. This is a check on
    the LABELS, never on a result: it runs on a quantity that no trial in this
    folder reports as a finding.
    """
    hi = [s for s, v in expect.items() if higher_in in str(v.values())]
    lo = [s for s in expect if s not in hi]
    if not hi or not lo:
        return {"checked": False, "reason": "the assumed map has no contrast on this variable"}
    vals = labels.set_index("suffix")[column].to_dict()
    ok = min(vals[s] for s in hi) > max(vals[s] for s in lo)
    return {
        "checked": True,
        "variable": column,
        "expected_higher": hi,
        "expected_lower": lo,
        "values": {s: float(vals[s]) for s in vals},
        "corroborated": bool(ok),
    }
