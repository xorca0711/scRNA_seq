"""Shared helpers for the Sikkema-2023 trials: run records and row-wise reads.

Every trial writes a run record before it reads any data table (the frozen
rules) and completes it afterwards (inputs with sizes and modification
times, package versions, outputs). Nothing here fits a model.
"""

from __future__ import annotations

import datetime as dt
import json
import platform
import sys
from pathlib import Path

import numpy as np
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
