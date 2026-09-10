"""Run records, row-wise h5ad reads, Markdown tables and a per-animal strip plot.

Shared by the focused analyses under analysis/scripts/ (10_phase_timecourse.py,
11_myeloid_focus.py). Mirrors the helpers in
Thesis/gate1_04_sikkema_2023_hlca/trials/trial_utils.py so that the analysis/
tree does not import from Thesis/.

A run record is written with the frozen rules BEFORE any data table is opened
and completed afterwards with the inputs (size, modification time), package
versions, results and outputs. Nothing here fits a model.
"""

from __future__ import annotations

import datetime as dt
import json
import platform
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp

REPO = Path(__file__).resolve().parents[2]


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def file_facts(path: Path) -> dict:
    st = path.stat()
    rel = str(path.relative_to(REPO)) if path.is_relative_to(REPO) else str(path)
    return {"path": rel, "bytes": st.st_size,
            "modified": dt.datetime.fromtimestamp(st.st_mtime, dt.timezone.utc).isoformat(timespec="seconds")}


def package_versions() -> dict:
    out = {"python": sys.version.split()[0], "platform": platform.platform()}
    for name in ("numpy", "scipy", "pandas", "matplotlib", "anndata", "scanpy", "h5py", "igraph"):
        try:
            mod = __import__(name)
            out[name] = getattr(mod, "__version__", "unknown")
        except Exception:  # noqa: BLE001
            out[name] = "not importable"
    return out


class RunRecord:
    """Write the frozen rules first; add inputs, results and outputs later."""

    def __init__(self, path: Path, analysis: str, rules: dict, notes: str = "") -> None:
        self.path = path
        self.record = {
            "analysis": analysis,
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
        rel = str(Path(path).relative_to(REPO))
        if rel not in self.record["outputs"]:
            self.record["outputs"].append(rel)

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
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def read_csr_rows(group, rows: np.ndarray, n_cols: int) -> sp.csr_matrix:
    """Read selected rows of an h5ad CSR group without loading the matrix.

    Consecutive row indices are read as one slice so that a compartment of a
    few thousand cells costs a few thousand small reads at most.
    """
    rows = np.asarray(sorted(set(int(r) for r in rows)))
    indptr = group["indptr"][:]
    data_parts, index_parts, new_indptr = [], [], [0]
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
        data_parts.append(group["data"][lo:hi])
        index_parts.append(group["indices"][lo:hi])
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
        cells = ([str(idx)] if index else []) + [
            "" if (isinstance(v, float) and np.isnan(v)) else str(v) for v in row.tolist()]
        lines.append("| " + " | ".join(c.replace("|", "/") for c in cells) + " |")
    return "\n".join(lines)


def shannon(fractions) -> float:
    fr = np.asarray(list(fractions), dtype=float)
    fr = fr[fr > 0]
    return float(-(fr * np.log(fr)).sum())


def strip_by_day(ax, df, ycol: str, days: list, colour: str, day_col: str = "day",
                 median: bool = True, size: int = 26, hollow: bool = False, seed: int = 0) -> None:
    """One point per animal at categorical day positions, medians joined by a line.

    Days are placed at equal spacing (0, 6, 11 ... 366 are not on a linear
    axis), so the x axis reads as an ordered series of harvests, not as time.
    """
    pos = {d: i for i, d in enumerate(days)}
    x = df[day_col].map(pos)
    keep = x.notna()
    d = df[keep]
    xv = x[keep].to_numpy(dtype=float)
    jitter = np.random.default_rng(seed).uniform(-0.12, 0.12, size=len(d))
    if hollow:
        ax.scatter(xv + jitter, d[ycol], s=size, facecolors="none", edgecolors=colour,
                   alpha=0.9, zorder=3, linewidths=1.0)
    else:
        ax.scatter(xv + jitter, d[ycol], s=size, color=colour, alpha=0.85, zorder=3, linewidths=0)
    if median and len(d):
        med = d.groupby(day_col)[ycol].median()
        ax.plot([pos[k] for k in med.index], med.to_numpy(), color=colour, lw=1.2, alpha=0.6, zorder=2)
    ax.set_xticks(range(len(days)))
    ax.set_xticklabels([str(k) for k in days])
