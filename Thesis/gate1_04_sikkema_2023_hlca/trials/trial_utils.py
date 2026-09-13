"""Shared helpers for every trial in this repository: run records, row-wise
reads, and the generic deposit readers.

Every trial writes a run record before it reads any data table (the frozen
rules) and completes it afterwards (inputs with sizes and modification
times, package versions, outputs). Nothing here fits a model.
"""

from __future__ import annotations

import datetime as dt
import gzip
import json
import platform
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[3]


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def file_facts(path: Path) -> dict:
    st = path.stat()
    return {"path": str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path),
            "bytes": st.st_size,
            "modified": dt.datetime.fromtimestamp(st.st_mtime, dt.timezone.utc).isoformat(timespec="seconds")}


def package_versions() -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for name in ("numpy", "scipy", "pandas", "anndata", "scanpy", "h5py", "sklearn", "leidenalg", "igraph", "openpyxl"):
        try:
            mod = __import__(name)
            out[name] = getattr(mod, "__version__", "unknown")
        except Exception:  # noqa: BLE001
            out[name] = "not importable"
    return out


class RunRecord:
    """Write the frozen rules first; add inputs, results and outputs later."""

    def __init__(self, path: Path, trial: str, rules: dict, notes: str = "") -> None:
        self.path = path
        self.record = {
            "trial": trial,
            "rules_frozen_at": utc_now(),
            "rules": rules,
            "notes": notes,
            "inputs": [],
            "outputs": [],
            "results": {},
            "software": package_versions(),
        }
        self.flush()

    def add_input(self, path: Path) -> None:
        self.record["inputs"].append(file_facts(Path(path)))
        self.flush()

    def add_output(self, path: Path) -> None:
        self.record["outputs"].append(str(Path(path).relative_to(REPO)))

    def set(self, key: str, value) -> None:
        self.record["results"][key] = value

    def finish(self) -> None:
        self.record["finished_at"] = utc_now()
        self.flush()

    def flush(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.record, indent=2, default=_json_default), encoding="utf-8")


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def read_csr_rows(group, rows: np.ndarray, n_cols: int) -> sp.csr_matrix:
    """Read selected rows of an h5ad CSR group without loading the matrix.

    Consecutive row indices are read as one slice so that a cluster of a few
    thousand cells costs a few thousand small reads at most.
    """
    rows = np.asarray(sorted(set(int(r) for r in rows)))
    indptr = group["indptr"][:]
    data_parts, index_parts, new_indptr = [], [], [0]
    # group consecutive rows into runs
    runs = []
    start = prev = rows[0]
    for r in rows[1:]:
        if r == prev + 1:
            prev = r
        else:
            runs.append((start, prev))
            start = prev = r
    runs.append((start, prev))
    for a, b in runs:
        lo, hi = int(indptr[a]), int(indptr[b + 1])
        d = group["data"][lo:hi]
        i = group["indices"][lo:hi]
        data_parts.append(d)
        index_parts.append(i)
        for r in range(a, b + 1):
            new_indptr.append(new_indptr[-1] + int(indptr[r + 1] - indptr[r]))
    data = np.concatenate(data_parts) if data_parts else np.array([], dtype=np.float32)
    indices = np.concatenate(index_parts) if index_parts else np.array([], dtype=np.int32)
    return sp.csr_matrix((data, indices, np.asarray(new_indptr)), shape=(len(rows), n_cols))


def df_to_markdown(df, index: bool = True) -> str:
    """Render a DataFrame as a GitHub Markdown table without the tabulate dependency."""
    cols = list(df.columns)
    header = ([str(df.index.name or "")] if index else []) + [str(c) for c in cols]
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    for idx, row in df.iterrows():
        cells = ([str(idx)] if index else []) + ["" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in row.tolist()]
        lines.append("| " + " | ".join(c.replace("|", "/") for c in cells) + " |")
    return "\n".join(lines)


def shannon(fractions) -> float:
    fr = np.asarray(list(fractions), dtype=float)
    fr = fr[fr > 0]
    return float(-(fr * np.log(fr)).sum())


def donor_threshold(n: int, private: float = 0.95) -> float | None:
    if n < 2:
        return None
    rest = (1 - private) / (n - 1)
    return shannon([private] + [rest] * (n - 1))


# Generic deposit readers. These lived in the Cardoso trials' own helper module
# until 2026-09-13, when the Choi-2020 folder needed them too and a gate1 folder
# importing from a gate2 folder would have been the wrong dependency direction.
# cardoso_utils re-exports them, so every trial written against it still works.

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
