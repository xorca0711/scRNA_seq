"""Shared helpers for the Cardoso-2026 trials.

Run records come from the Sikkema trial helpers so that there is one
implementation of the frozen-rules record in the repository; this module adds
what the Cardoso series need and the earlier series did not: MatrixMarket
triplet loading (these deposits are mtx/tsv, not 10x HDF5) and a parser for
the deposited GEO SOFT family files, so that the design table in every trial
is read from an artefact rather than transcribed by hand.
"""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "raw_data"

# One implementation of RunRecord for the whole repository.
_SIKKEMA_TRIALS = REPO / "Thesis" / "gate1_04_sikkema_2023_hlca" / "trials"
if str(_SIKKEMA_TRIALS) not in sys.path:
    sys.path.insert(0, str(_SIKKEMA_TRIALS))
from trial_utils import (RunRecord, df_to_markdown, file_facts,  # noqa: E402
                         package_versions, shannon, utc_now)

__all__ = ["REPO", "RAW", "RunRecord", "df_to_markdown", "file_facts",
           "package_versions", "shannon", "utc_now", "read_mtx_triplet",
           "parse_soft", "ENSEMBL_ID_RE", "qc_metrics"]

# A feature whose identifier is not an Ensembl gene ID is not a gene. The
# earlier mouse series carried SiteA/SiteB lineage contigs; this rule is the
# general form of the same decision and is applied before any QC metric.
ENSEMBL_ID_RE = re.compile(r"^ENS[A-Z]*G\d+")


def read_mtx_triplet(matrix: Path, features: Path, barcodes: Path):
    """Read a CellRanger mtx/tsv triplet as (X cells x genes, var, obs_names).

    The deposited matrix is genes x cells; it is transposed once here so that
    every downstream step sees the scverse orientation.
    """
    from scipy.io import mmread

    with gzip.open(matrix, "rb") as handle:
        M = mmread(handle)
    X = sp.csr_matrix(M.T)
    X.eliminate_zeros()

    var = pd.read_csv(features, sep="\t", header=None, compression="gzip",
                      names=["gene_id", "gene_symbol", "feature_type"],
                      dtype=str).fillna("")
    bc = pd.read_csv(barcodes, sep="\t", header=None, compression="gzip",
                     names=["barcode"], dtype=str)["barcode"].to_numpy()
    if X.shape != (len(bc), len(var)):
        raise ValueError(f"{matrix.name}: matrix {X.shape} does not match "
                         f"{len(bc)} barcodes x {len(var)} features")
    return X, var, bc


def parse_soft(path: Path) -> dict:
    """Parse the fields of a GEO SOFT family file that describe the design."""
    series: dict[str, list[str]] = {}
    samples: list[dict] = []
    current: dict | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                current = {"gsm": line.split("=", 1)[1].strip(),
                           "characteristics": [], "data_processing": [],
                           "description": []}
                samples.append(current)
            elif line.startswith("!Series_") and current is None:
                key, _, value = line[1:].partition(" = ")
                series.setdefault(key, []).append(value)
            elif line.startswith("!Sample_") and current is not None:
                key, _, value = line[1:].partition(" = ")
                key = key.replace("Sample_", "")
                if key.startswith("characteristics"):
                    current["characteristics"].append(value)
                elif key.startswith("data_processing"):
                    current["data_processing"].append(value)
                elif key.startswith("description"):
                    current["description"].append(value)
                else:
                    current.setdefault(key, value)
    return {"series": series, "samples": samples}


def qc_metrics(X: sp.csr_matrix, var: pd.DataFrame, species: str) -> pd.DataFrame:
    """Per-cell QC metrics. No cell is filtered here; Gate 0 only measures."""
    symbols = var["gene_symbol"].astype(str)
    if species == "mouse":
        mt = symbols.str.startswith("mt-")
        ribo = symbols.str.startswith(("Rps", "Rpl")) & ~symbols.str.contains("-ps", case=False)
        hb = symbols.str.startswith(("Hba", "Hbb"))
    else:
        mt = symbols.str.startswith("MT-")
        ribo = symbols.str.startswith(("RPS", "RPL")) & ~symbols.str.contains("-ps", case=False)
        hb = symbols.str.startswith(("HBA", "HBB"))

    total = np.asarray(X.sum(axis=1)).ravel()
    with np.errstate(invalid="ignore", divide="ignore"):
        out = pd.DataFrame({
            "total_counts": total,
            "n_genes": X.getnnz(axis=1),
            "pct_mt": 100 * np.asarray(X[:, mt.to_numpy()].sum(axis=1)).ravel() / np.maximum(total, 1),
            "pct_ribo": 100 * np.asarray(X[:, ribo.to_numpy()].sum(axis=1)).ravel() / np.maximum(total, 1),
            "pct_hb": 100 * np.asarray(X[:, hb.to_numpy()].sum(axis=1)).ravel() / np.maximum(total, 1),
        })
    return out
