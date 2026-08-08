"""
Shared helpers for the scRNA-seq pipeline: dataset configuration, memory-aware
IO, MAD-based QC thresholds, and figure saving.

Nothing in here writes to the raw-data directory.
"""

from __future__ import annotations

import gc
import os
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "raw_data"
ANALYSIS = REPO / "analysis"

RANDOM_SEED = 0


# ---------------------------------------------------------------------------
# dataset configuration - derived from the Phase-1 scan, not assumed up front
# ---------------------------------------------------------------------------
@dataclass
class DatasetConfig:
    name: str
    raw_dir: Path
    h5_glob: str
    species: str                      # "mouse" | "human"
    outdir: Path
    sample_regex: str                 # captures the sample id from the filename
    metadata_csv: Path | None = None
    # features that are not genes and must leave the expression matrix
    non_gene_features: list[str] = field(default_factory=list)
    # symbols are not unique in this reference -> key the var index on IDs
    use_gene_ids_as_index: bool = False
    notes: str = ""


def gse262927(outdir: Path | None = None) -> DatasetConfig:
    return DatasetConfig(
        name="GSE262927",
        raw_dir=RAW / "GSE262927" / "GSE262927_RAW",
        h5_glob="GSM*.h5",
        species="mouse",
        # Each series owns a subdirectory, so the two are symmetric and nothing
        # at the top level of analysis/ is series-specific. This departs from
        # the canonical paths in the original brief (analysis/figures/...,
        # analysis/processed/final_clustered.h5ad); see analysis/LAYOUT.md.
        outdir=outdir or (ANALYSIS / "GSE262927"),
        sample_regex=r"^GSM\d+_(.+)\.h5$",
        metadata_csv=RAW / "GSE262927" / "GSE262927_CellMetaData.csv",
        non_gene_features=["SiteA", "SiteB"],
        use_gene_ids_as_index=False,
        notes=(
            "Mouse lung, H1N1 influenza lineage-tracing time course. "
            "10x 3' v3, matrices re-exported from R (rhdf5), counts stored as "
            "float64 with explicit zeros. Gene space identical across all 33 "
            "files. SiteA/SiteB are lineage-tracing reporter contigs, not "
            "genes."
        ),
    )


def gse178360(outdir: Path | None = None) -> DatasetConfig:
    return DatasetConfig(
        name="GSE178360",
        raw_dir=RAW / "GSE178360" / "GSE178360_RAW",
        # note the missing underscore in GSM5388413_DD073Rfiltered_...
        h5_glob="GSM*filtered_feature_bc_matrix.h5",
        species="human",
        outdir=outdir or (ANALYSIS / "GSE178360"),
        sample_regex=r"^GSM\d+_(.+?)_?filtered_feature_bc_matrix\.h5$",
        metadata_csv=None,
        non_gene_features=[],
        use_gene_ids_as_index=True,
        notes=(
            "Human lung/airway, 10x 3' v3, CellRanger 5.0.1. No metadata file "
            "exists in the series. DD073R was quantified against a different "
            "GRCh38 annotation build (36,601 vs 33,538 genes) so samples must "
            "be intersected on Ensembl gene ID. Gene symbols are NOT unique."
        ),
    )


DATASETS = {"GSE262927": gse262927, "GSE178360": gse178360}


# ---------------------------------------------------------------------------
# logging
# ---------------------------------------------------------------------------
class Tee:
    """Write to stdout and to a log file at the same time."""

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = open(path, "a", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, msg):
        self.stdout.write(msg)
        self.fh.write(msg)

    def flush(self):
        self.stdout.flush()
        self.fh.flush()

    def close(self):
        self.fh.close()


_T0 = time.time()


def log(msg: str = "") -> None:
    el = time.time() - _T0
    print(f"[{el:7.1f}s] {msg}", flush=True)


def mem_report(tag: str = "") -> str:
    """Report resident memory of this process, best effort."""
    try:
        import ctypes
        import ctypes.wintypes as wt

        class PMC(ctypes.Structure):
            _fields_ = [
                ("cb", wt.DWORD), ("PageFaultCount", wt.DWORD),
                ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPagedPoolUsage", ctypes.c_size_t),
                ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t),
            ]

        pmc = PMC()
        pmc.cb = ctypes.sizeof(PMC)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(pmc),
            pmc.cb)
        gb = pmc.WorkingSetSize / 1024 ** 3
        peak = pmc.PeakWorkingSetSize / 1024 ** 3
        return f"[mem{' ' + tag if tag else ''}: {gb:.2f} GB now, {peak:.2f} GB peak]"
    except Exception:  # noqa: BLE001
        return ""


def free_mem() -> None:
    gc.collect()


# ---------------------------------------------------------------------------
# MAD-based, per-sample QC thresholds
# ---------------------------------------------------------------------------
def mad(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.median(np.abs(x - np.median(x))))


def mad_bounds(x: np.ndarray, nmads: float) -> tuple[float, float]:
    """Median +/- nmads * MAD (unscaled MAD, the scverse convention)."""
    x = np.asarray(x, dtype=float)
    med = float(np.median(x))
    m = mad(x)
    return med - nmads * m, med + nmads * m


@dataclass
class SampleThresholds:
    sample: str
    min_genes: int
    max_genes: int
    min_counts: int
    max_counts: int
    max_pct_mt: float
    max_pct_top20: float
    n_mads_counts: float
    n_mads_genes: float
    n_mads_mt: float
    rationale: str


def derive_thresholds(obs: pd.DataFrame, sample: str,
                      nmads_counts: float = 5.0,
                      nmads_genes: float = 5.0,
                      nmads_mt: float = 3.0,
                      hard_min_genes: int = 200,
                      hard_min_counts: int = 500,
                      mt_floor: float = 5.0,
                      mt_ceiling: float = 20.0) -> SampleThresholds:
    """Choose per-sample thresholds from that sample's own distributions.

    Counts and gene numbers are bounded on the log1p scale, which is where the
    distributions are approximately symmetric; the mitochondrial cut is an
    upper bound only, taken from the sample's own MAD but bracketed into a
    biologically sane window so that a sample with a pathologically wide or
    narrow mito distribution cannot produce an absurd cut.
    """
    lo_c, hi_c = mad_bounds(np.log1p(obs["total_counts"].to_numpy()), nmads_counts)
    lo_g, hi_g = mad_bounds(np.log1p(obs["n_genes_by_counts"].to_numpy()), nmads_genes)
    _, hi_mt = mad_bounds(obs["pct_counts_mt"].to_numpy(), nmads_mt)
    _, hi_t20 = mad_bounds(obs["pct_counts_in_top_20_genes"].to_numpy(), 5.0)

    min_counts = int(max(hard_min_counts, np.expm1(lo_c)))
    max_counts = int(np.expm1(hi_c))
    min_genes = int(max(hard_min_genes, np.expm1(lo_g)))
    max_genes = int(np.expm1(hi_g))
    max_pct_mt = float(min(max(hi_mt, mt_floor), mt_ceiling))
    max_pct_top20 = float(min(hi_t20, 100.0))

    rationale = (
        f"log1p(counts) median+/-{nmads_counts}MAD -> [{min_counts},{max_counts}]; "
        f"log1p(genes) median+/-{nmads_genes}MAD -> [{min_genes},{max_genes}] "
        f"(lower bounds floored at {hard_min_counts} counts / {hard_min_genes} "
        f"genes); pct_mt median+{nmads_mt}MAD={hi_mt:.2f} bracketed to "
        f"[{mt_floor},{mt_ceiling}] -> {max_pct_mt:.2f}; "
        f"pct_top20 median+5MAD -> {max_pct_top20:.2f}"
    )
    return SampleThresholds(
        sample=sample, min_genes=min_genes, max_genes=max_genes,
        min_counts=min_counts, max_counts=max_counts,
        max_pct_mt=max_pct_mt, max_pct_top20=max_pct_top20,
        n_mads_counts=nmads_counts, n_mads_genes=nmads_genes,
        n_mads_mt=nmads_mt, rationale=rationale)


def apply_thresholds(obs: pd.DataFrame, th: SampleThresholds) -> pd.Series:
    """Boolean mask of cells that pass every threshold."""
    keep = (
        (obs["total_counts"] >= th.min_counts)
        & (obs["total_counts"] <= th.max_counts)
        & (obs["n_genes_by_counts"] >= th.min_genes)
        & (obs["n_genes_by_counts"] <= th.max_genes)
        & (obs["pct_counts_mt"] <= th.max_pct_mt)
        & (obs["pct_counts_in_top_20_genes"] <= th.max_pct_top20)
    )
    return keep


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------
def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": 9,
        "axes.titlesize": 11,
        "axes.labelsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "pdf.fonttype": 42,      # editable text in Illustrator
        "ps.fonttype": 42,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })
    return plt


def save_fig(fig, outdir: Path, name: str, formats=("png", "pdf")) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for ext in formats:
        p = outdir / f"{name}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        paths.append(p)
    import matplotlib.pyplot as plt
    plt.close(fig)
    return paths


def figure_of(obj):
    """Get a Figure out of whatever a scanpy plotting call returned."""
    import matplotlib.pyplot as plt
    if obj is None:
        return plt.gcf()
    # scanpy's BasePlot subclasses (DotPlot, MatrixPlot, ...) expose savefig
    # too, so they must be tested for first.
    if hasattr(obj, "make_figure") and hasattr(obj, "fig"):
        if getattr(obj, "fig", None) is None:
            obj.make_figure()
        return obj.fig
    if hasattr(obj, "savefig"):                      # already a Figure
        return obj
    if hasattr(obj, "fig"):
        return obj.fig
    if hasattr(obj, "get_figure"):                   # Axes
        return obj.get_figure()
    if isinstance(obj, dict):
        for v in obj.values():
            if hasattr(v, "get_figure"):
                return v.get_figure()
    if isinstance(obj, (list, tuple)) and obj:
        return figure_of(obj[0])
    return plt.gcf()


# ---------------------------------------------------------------------------
# misc
# ---------------------------------------------------------------------------
def sample_number(text: str) -> str | None:
    """Extract the numeric sample id shared by the h5 filenames and metadata."""
    m = re.search(r"scRNA[-_](\d+)", str(text))
    return m.group(1) if m else None


def make_unique(names) -> np.ndarray:
    """R's make.unique semantics, used when gene symbols repeat."""
    seen: dict[str, int] = {}
    out = []
    for n in names:
        if n in seen:
            seen[n] += 1
            out.append(f"{n}-{seen[n]}")
        else:
            seen[n] = 0
            out.append(n)
    return np.array(out, dtype=object)


def ensure_dirs(cfg: DatasetConfig) -> dict[str, Path]:
    o = cfg.outdir
    d = {
        "root": o,
        "scripts": ANALYSIS / "scripts",
        "logs": o / "logs",
        "inventory": o / "inventory",
        "qc": o / "qc",
        "processed": o / "processed",
        "tables": o / "tables",
        "fig_umap": o / "figures" / "umap",
        "fig_dot": o / "figures" / "dotplots",
        "fig_feat": o / "figures" / "featureplots",
        "fig_comp": o / "figures" / "composition",
        "fig_qc": o / "figures" / "qc",
        "epi": o / "epithelial_subanalysis",
        "epi_fig": o / "epithelial_subanalysis" / "figures",
        "epi_tab": o / "epithelial_subanalysis" / "tables",
        "shards": o / "processed" / "sample_shards",
    }
    for p in d.values():
        p.mkdir(parents=True, exist_ok=True)
    return d
