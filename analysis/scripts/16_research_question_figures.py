#!/usr/bin/env python
"""Render selected A1/A3/A4/A5 panels without editing authored documents.

Use --figures A3 --replot for a cache-only presentation run. Missing/invalid
caches fail closed. --rebuild-embeddings explicitly prepares fresh display
embeddings in the new run directory; neither mode deletes existing caches.
Every invocation writes images, tables and panel facts to a fresh directory
under analysis/figures/rq/renders/. Review before promoting individual images.
The original run_record.json and curated gallery captions remain unchanged.
A2 is owned by script 18; the retired implementation is preserved in the
2026-09-25 RQ migration archive. No model or hypothesis test is performed here.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "analysis" / "scripts"))
sys.path.insert(0, str(REPO / "Research Article" / "gate2_05_cardoso_2026" / "trials"))
sys.path.insert(0, str(REPO / "Research Article" / "gate1_02_choi_2020" / "datp_epigenetics" / "trials"))

import viz_style as vs  # noqa: E402
from pipeline_utils import RANDOM_SEED, SERIES_DIRS  # noqa: E402
from multiome_utils import (ATAC, GSE310539, RNA, peak_gene_table, read_barcodes,  # noqa: E402
                            read_features)
from m1b_corrected_pass import wells_corrected  # noqa: E402
from m1_closed_or_merely_silenced import AT1, AT2_IDENTITY, TRANSITIONAL  # noqa: E402

OUT = REPO / "analysis" / "figures" / "rq"
CACHE = OUT / "processed"
REBUILD = False
ACTIVE_FIGURES = ("A1", "A3", "A4", "A5")
MOUSE = SERIES_DIRS["GSE262927"]

RAMP = LinearSegmentedColormap.from_list("ramp", vs.RAMP)
DOT = dict(s=2.5, linewidths=0, rasterized=True)
LABEL_GENES = ["Cldn4", "Krt8", "Sftpc", "Axin2", "Il1r1"]
PHASES = [("baseline", [0]), ("active repair", [6, 11, 19, 25]),
          ("injury resolution", [42, 90]), ("long-term homeostasis", [366])]
# Preserve historical table keys; display only observed sampling intervals.
PHASE_DISPLAY = {"baseline": "baseline", "active repair": "6–25 dpi",
                 "injury resolution": "42–90 dpi", "long-term homeostasis": "366 dpi"}

RECORD: dict = {
    "script": "analysis/scripts/16_research_question_figures.py",
    "purpose": "visual aids for RESEARCH_QUESTIONS.md; not evidence; nothing tested",
    "seed": RANDOM_SEED,
    "embedding_recipe": {
        "multiome": "RNA only; genes in >= 3 nuclei; normalize_total 1e4; log1p; 2000 HVGs "
                    "(seurat flavour); scale, clip 10; PCA 30; neighbours 15; UMAP; per "
                    "deposit and per well, no batch correction",
        "capillary": "script 06 recipe (load_subset clusters 0,1,3,21 minus 4 within the "
                     "annotated cohort; embed n_hvg 2000, resolution 0.4; iCAP score from "
                     "Sparcl1, Ntrk2)",
        "myeloid": "tracked coordinates from myeloid_focus/tables/myeloid_cell_metadata.csv",
    },
    "labels": {
        "transitional_at_depth_available": "raw Cldn4 >= 1 and raw Krt8 >= 1 (the R4c form of "
                                           "trial M1c without downsampling); display only",
        "at2_nucleus": "raw Sftpc >= 1 and not transitional (the A1 rule)",
    },
    "panels": {},
    "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
}


def log(msg: str) -> None:
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def note(key: str, **kw) -> None:
    RECORD["panels"].setdefault(key, {}).update(kw)


# ---------------------------------------------------------------------------
# Multiome loading: one chunked pass, RNA kept as a sparse matrix for the
# selected wells, ATAC totals and selected promoter peaks for every nucleus.
# ---------------------------------------------------------------------------
def load_multiome(path: Path, suffixes: list[str], promoter_rows: np.ndarray | None,
                  chunk: int = 2000):
    feats = read_features(path)
    is_rna = (feats["feature_type"] == RNA).to_numpy()
    rna_rows = np.flatnonzero(is_rna)
    rna_pos = np.full(len(feats), -1, dtype=np.int64)
    rna_pos[rna_rows] = np.arange(len(rna_rows))
    prom_slot = np.full(len(feats), -1, dtype=np.int64)
    n_prom = 0 if promoter_rows is None else len(promoter_rows)
    if n_prom:
        prom_slot[promoter_rows] = np.arange(n_prom)
    bcs = read_barcodes(path)
    keep = bcs["suffix"].isin(suffixes).to_numpy()
    n_cells = len(bcs)
    atac_tot = np.zeros(n_cells, dtype=np.int64)
    prom = np.zeros((n_prom, n_cells), dtype=np.int32)
    data_parts, ind_parts, ptr = [], [], [0]
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        indptr = g["indptr"][:]
        for start in range(0, n_cells, chunk):
            stop = min(start + chunk, n_cells)
            lo, hi = int(indptr[start]), int(indptr[stop])
            data = g["data"][lo:hi]
            idx = g["indices"][lo:hi]
            offs = indptr[start:stop + 1] - lo
            for j in range(stop - start):
                col = start + j
                a, b = int(offs[j]), int(offs[j + 1])
                ii, dd = idx[a:b], data[a:b]
                r = is_rna[ii]
                atac_tot[col] = dd[~r].sum()
                if n_prom:
                    s = prom_slot[ii]
                    k = s >= 0
                    if k.any():
                        prom[s[k], col] = dd[k]
                if keep[col]:
                    data_parts.append(dd[r].astype(np.float32))
                    ind_parts.append(rna_pos[ii[r]])
                    ptr.append(ptr[-1] + int(r.sum()))
    X = sp.csr_matrix((np.concatenate(data_parts), np.concatenate(ind_parts),
                       np.asarray(ptr, dtype=np.int64)),
                      shape=(int(keep.sum()), len(rna_rows)))
    import anndata as ad
    obs = bcs[keep].reset_index(drop=True)
    obs.index = obs["barcode"].to_numpy()
    var = feats.iloc[rna_rows][["id", "name"]].set_index("name")
    a = ad.AnnData(X=X, obs=obs, var=var)
    a.var_names_make_unique()
    a.obs["atac_counts"] = atac_tot[keep]
    a.layers["counts"] = a.X.copy()
    return a, atac_tot, prom, keep


def promoter_rows_by_gene(path: Path, peaks: Path, genes: list[str]) -> dict[str, np.ndarray]:
    """Rows of the promoter peaks the vendor annotated to each gene."""
    feats = read_features(path)
    atac = feats[feats["feature_type"] == ATAC].set_index("id")
    # the vendor annotation names genes by Ensembl id; the RNA features map them
    rna = feats[feats["feature_type"] == RNA]
    id_of = dict(zip(rna["name"], rna["id"]))
    table = peak_gene_table(peaks)
    prom = table[table["peak_type"] == "promoter"]
    out = {}
    for g in genes:
        if g not in id_of:
            continue
        iv = prom[prom["gene"] == id_of[g]]["interval"]
        rows = atac.loc[atac.index.isin(iv), "row"].to_numpy()
        if len(rows):
            out[g] = np.unique(rows)
    assert out, "no promoter peak matched the AT2 identity genes"
    return out


def promoter_rows_for(path: Path, peaks: Path, genes: list[str]) -> np.ndarray:
    return np.unique(np.concatenate(list(promoter_rows_by_gene(path, peaks, genes).values())))


def ensure_promoter_score(a) -> None:
    """Recompute the promoter accessibility, in total and per gene, when a cached object lacks it."""
    if "atac_promoter_by_gene" in a.obsm and "atac_at2_promoter" in a.obs:
        return
    if not REBUILD:
        raise ValueError("Cache lacks promoter measurements; --replot never rewrites caches")
    from multiome_utils import stream_selected
    path, peaks = GSE310539["matrix"], GSE310539["peaks"]
    by_gene = promoter_rows_by_gene(path, peaks, AT2_IDENTITY)
    rows = np.unique(np.concatenate(list(by_gene.values())))
    _, picked = stream_selected(path, rows)
    bcs = read_barcodes(path)
    pos = pd.Series(np.arange(len(bcs)), index=bcs["barcode"]).loc[a.obs_names].to_numpy()
    slot = {r: i for i, r in enumerate(rows)}
    per_gene = np.stack([picked[[slot[r] for r in by_gene[g]], :][:, pos].sum(0) for g in by_gene], axis=1)
    a.obsm["atac_promoter_by_gene"] = per_gene.astype(np.int32)
    a.uns["promoter_genes"] = list(by_gene)
    a.uns["promoter_peaks_per_gene"] = {g: int(len(r)) for g, r in by_gene.items()}
    prom_sum = per_gene.sum(1)
    a.obs["atac_at2_promoter"] = np.log1p(1e4 * prom_sum / np.maximum(a.obs["atac_counts"].to_numpy(), 1))
    a.uns["n_promoter_peaks"] = int(len(rows))
    a.write_h5ad(CACHE / "gse310539_wildtype.h5ad")
    log(f"  promoter accessibility recomputed over {len(rows)} peaks for {len(by_gene)} genes")


def detection_at_budget(a, groups: list[tuple[str, np.ndarray]]) -> tuple[pd.DataFrame, dict]:
    """Detection fraction of each AT2 identity gene, RNA and promoter chromatin,
    at one depth budget per modality: trial M1e's rule, the 20th percentile of
    the labelled group's counts, applied here across both wildtype wells."""
    from multiome_utils import detect_prob_at_depth
    promoter_genes = [str(g) for g in a.uns["promoter_genes"]]  # an h5ad round trip turns the list into an array
    genes = [g for g in promoter_genes if g in a.var_names]
    rna_counts = a.layers["counts"]
    rna_total = np.asarray(rna_counts.sum(1)).ravel()
    atac_total = a.obs["atac_counts"].to_numpy()
    trans = groups[-1][1]
    budget = {"rna": int(np.percentile(rna_total[trans], 20)), "atac": int(np.percentile(atac_total[trans], 20))}
    rows = []
    for name, mask in groups:
        for j, g in enumerate(genes):
            rc = np.asarray(rna_counts[:, a.var_names.get_loc(g)].todense()).ravel()
            ac = a.obsm["atac_promoter_by_gene"][:, promoter_genes.index(g)]
            pr = detect_prob_at_depth(rc[mask], rna_total[mask], budget["rna"])
            pa = detect_prob_at_depth(ac[mask], atac_total[mask], budget["atac"])
            rows.append({"group": name.replace("\n", " "), "gene": g,
                         "rna_detection_pct": 100 * np.nanmean(pr), "rna_nuclei": int(np.isfinite(pr).sum()),
                         "atac_detection_pct": 100 * np.nanmean(pa), "atac_nuclei": int(np.isfinite(pa).sum())})
    return pd.DataFrame(rows), budget


def heat(ax, M: np.ndarray, rows: list[str], cols: list[str], title: str, vmax: float) -> None:
    ax.imshow(M, cmap=RAMP, vmin=0, vmax=vmax, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.0f}", ha="center", va="center", fontsize=8,
                    color=vs.INK if M[i, j] < 0.6 * vmax else vs.SURFACE)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, fontsize=8)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=8, style="italic")
    ax.set_title(title, fontsize=10)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)


def embed_rna(a, tag: str):
    import scanpy as sc
    sc.pp.filter_genes(a, min_cells=3)
    for g in LABEL_GENES:
        a.obs[f"raw_{g}"] = np.asarray(a.layers["counts"][:, a.var_names.get_loc(g)].todense()).ravel() if g in a.var_names else 0
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=2000)
    sub = a[:, a.var["highly_variable"]].copy()
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=30, random_state=RANDOM_SEED)
    a.obsm["X_pca"] = sub.obsm["X_pca"]
    del sub
    sc.pp.neighbors(a, n_neighbors=15, n_pcs=30, random_state=RANDOM_SEED)
    sc.tl.umap(a, random_state=RANDOM_SEED)
    for name, genes in (("AT2_identity", AT2_IDENTITY), ("transitional", TRANSITIONAL), ("AT1", AT1)):
        present = [g for g in genes if g in a.var_names]
        sc.tl.score_genes(a, present, score_name=f"score_{name}", random_state=RANDOM_SEED)
    a.obs["transitional_label"] = (a.obs["raw_Cldn4"] >= 1) & (a.obs["raw_Krt8"] >= 1)
    a.obs["at2_nucleus"] = (a.obs["raw_Sftpc"] >= 1) & ~a.obs["transitional_label"]
    log(f"  {tag}: {a.n_obs:,} nuclei, {a.n_vars:,} genes embedded")
    return a


def sha256(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def cached(tag: str, build, *, backed=None):
    import anndata as ad
    p = CACHE / f"{tag}.h5ad"
    if REBUILD:
        if p.exists():
            raise FileExistsError(f"Refusing to overwrite display cache: {p}")
        ad.settings.allow_write_nullable_strings = True
        a = build()
        CACHE.mkdir(parents=True, exist_ok=True)
        a.write_h5ad(p)
    else:
        if not p.is_file():
            raise FileNotFoundError(f"Missing cache {p}; no embedding will be rebuilt in --replot mode")
        a = ad.read_h5ad(p, backed=backed)
        if "X_umap" not in a.obsm:
            if getattr(a, "isbacked", False):
                a.file.close()
            raise ValueError(f"Cache lacks X_umap: {p}; no automatic rebuild")
        log(f"  {tag}: cached embedding")
    RECORD.setdefault("input_sha256", {})[p.relative_to(REPO).as_posix()] = sha256(p)
    return a


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------
def umap_axes(ax, title: str) -> None:
    vs.bare(ax)
    ax.set_title(title)
    ax.set_aspect("equal")


def feature(ax, xy, values, title, vmax=None, label=""):
    order = np.argsort(values)
    m = ax.scatter(xy[order, 0], xy[order, 1], c=values[order], cmap=RAMP,
                   vmin=0, vmax=vmax if vmax is not None else np.percentile(values, 99), **DOT)
    umap_axes(ax, title)
    cb = plt.colorbar(m, ax=ax, fraction=0.04, pad=0.02)
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, labelsize=8, colors=vs.MUTED)
    if label:
        cb.set_label(label, color=vs.INK_2, fontsize=8)


def highlight(ax, xy, mask, colour, title, label):
    ax.scatter(xy[~mask, 0], xy[~mask, 1], c=vs.DEEMPH, **DOT)
    ax.scatter(xy[mask, 0], xy[mask, 1], c=colour, **DOT)
    umap_axes(ax, title)
    ax.text(0.02, 0.98, label, transform=ax.transAxes, va="top", ha="left",
            fontsize=8.5, color=vs.INK_2)


def violins(ax, groups: list[tuple[str, np.ndarray, str]], ylabel: str, annotate=None,
            points: bool = False):
    """Violins of every value; with points=True the positive values are also
    drawn as jittered dots, which is what makes a sparse transcript readable."""
    pos = np.arange(len(groups)) + 1
    if points:
        rng = np.random.default_rng(RANDOM_SEED)
        for x, (_, vals, colour) in zip(pos, groups):
            nz = vals[vals > 0]
            ax.scatter(x + rng.uniform(-0.28, 0.28, len(nz)), nz, s=3, color=colour, alpha=0.35,
                       linewidths=0, rasterized=True, zorder=3)
    parts = ax.violinplot([g[1] for g in groups], positions=pos, widths=0.8,
                          showmedians=False, showextrema=False)
    for body, (_, _, colour) in zip(parts["bodies"], groups):
        body.set_facecolor(colour)
        body.set_edgecolor(colour)
        body.set_alpha(0.55)
        body.set_linewidth(0.8)
    for x, (_, vals, colour) in zip(pos, groups):
        ax.plot([x - 0.25, x + 0.25], [np.median(vals)] * 2, color=vs.INK, lw=1.4)
        if annotate is not None:
            ax.text(x, ax.get_ylim()[1] if False else np.max(vals) * 1.03 + 1e-9,
                    annotate(vals), ha="center", va="bottom", fontsize=8, color=vs.INK_2)
    ax.set_xticks(pos)
    ax.set_xticklabels([g[0] for g in groups], fontsize=8.5)
    ax.set_ylabel(ylabel)
    vs.recessive(ax, grid_axis="y")


def fraction_label(vals: np.ndarray) -> str:
    return f"{100 * np.mean(vals > 0):.1f}% detected"


def save(fig, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / name
    fig.savefig(p, dpi=200, bbox_inches="tight")
    plt.close(fig)
    log(f"  wrote {p.relative_to(REPO).as_posix()}")
    return p


# ---------------------------------------------------------------------------
# A1 and A4: GSE310539
# ---------------------------------------------------------------------------
def build_gse310539():
    path, peaks = GSE310539["matrix"], GSE310539["peaks"]
    prom_rows = promoter_rows_for(path, peaks, AT2_IDENTITY)
    log(f"  GSE310539: {len(prom_rows)} promoter peaks annotated to the AT2 identity genes")
    wt = [s for s, m in GSE310539["libraries"].items() if m["genotype"] == "wildtype"]
    a, atac_tot, prom, keep = load_multiome(path, wt, prom_rows)
    a.obs["well"] = a.obs["suffix"].map({s: m["name"] for s, m in GSE310539["libraries"].items()})
    # promoter accessibility of the AT2 identity genes, per 10,000 ATAC counts
    prom_sum = prom.sum(0)[keep]
    a.obs["atac_at2_promoter"] = np.log1p(1e4 * prom_sum / np.maximum(a.obs["atac_counts"].to_numpy(), 1))
    a.uns["n_promoter_peaks"] = int(len(prom_rows))
    return embed_rna(a, "GSE310539 wildtype wells")


def all_well_counts_gse310539() -> pd.DataFrame:
    """Raw counts of the label genes in every nucleus of the four wells."""
    from multiome_utils import stream_selected
    path = GSE310539["matrix"]
    feats = read_features(path)
    rows = {g: int(feats[(feats["feature_type"] == RNA) & (feats["name"] == g)]["row"].iloc[0]) for g in LABEL_GENES}
    order = sorted(rows.items(), key=lambda kv: kv[1])
    totals, picked = stream_selected(path, np.array([r for _, r in order]))
    bcs = read_barcodes(path)
    df = pd.DataFrame({g: picked[i] for i, (g, _) in enumerate(order)})
    df["well"] = bcs["suffix"].map({s: m["name"] for s, m in GSE310539["libraries"].items()}).to_numpy()
    df["rna_counts"] = totals["rna_counts"].to_numpy()
    df["transitional_label"] = (df["Cldn4"] >= 1) & (df["Krt8"] >= 1)
    df["at2_nucleus"] = (df["Sftpc"] >= 1) & ~df["transitional_label"]
    return df


def figure_a1(a) -> dict:
    xy = a.obsm["X_umap"]
    well = a.obs["well"].to_numpy()
    trans = a.obs["transitional_label"].to_numpy()
    colours = {"wildtype_PBS": vs.SLOT[1], "wildtype_SeV": vs.SLOT[2]}
    fig, axes = plt.subplots(3, 3, figsize=(15.0, 12.2), gridspec_kw={"width_ratios": [1, 1, 1.3]})
    ax = axes[0, 0]
    for w in ("wildtype_PBS", "wildtype_SeV"):
        m = well == w
        ax.scatter(xy[m, 0], xy[m, 1], c=colours[w], label=f"{w.replace('_', ' ')} (n = {m.sum():,})", **DOT)
    umap_axes(ax, "a  Wildtype nuclei by well")
    ax.legend(loc="lower left", markerscale=6, handletextpad=0.3)
    frac = {w: 100 * trans[well == w].mean() for w in colours}
    highlight(axes[0, 1], xy, trans, vs.SLOT[2], "b  Transitional label (Cldn4 and Krt8 detected)",
              f"PBS {frac['wildtype_PBS']:.1f}%, SeV {frac['wildtype_SeV']:.1f}% of nuclei")
    feature(axes[0, 2], xy, a.obs["score_AT2_identity"].to_numpy().clip(min=0),
            "c  AT2 identity programme, RNA score", label="score (clipped at 0)")
    feature(axes[1, 0], xy, a.obs["atac_at2_promoter"].to_numpy(),
            f"d  AT2 identity promoters, chromatin ({a.uns['n_promoter_peaks']} peaks)",
            label="log1p counts per 10,000 ATAC counts")
    groups = [("PBS\nreference", (well == "wildtype_PBS") & ~trans, vs.SLOT[1]),
              ("SeV\nreference", (well == "wildtype_SeV") & ~trans, vs.SLOT[2]),
              ("SeV\ntransitional", (well == "wildtype_SeV") & trans, vs.SLOT[2])]
    rna = a.obs["score_AT2_identity"].to_numpy()
    atac = a.obs["atac_at2_promoter"].to_numpy()
    violins(axes[1, 1], [(n, rna[m], c) for n, m, c in groups], "AT2 identity RNA score")
    axes[1, 1].set_title("e  RNA score, by group")
    # panel f: the form of the registered statistic, detection at one depth budget
    det, budget = detection_at_budget(a, [(n, m) for n, m, _ in groups])
    spec_f = axes[1, 2].get_subplotspec()
    axes[1, 2].remove()
    gsf = spec_f.subgridspec(1, 2, wspace=0.45)
    axf1, axf2 = fig.add_subplot(gsf[0, 0]), fig.add_subplot(gsf[0, 1])
    genes = list(dict.fromkeys(det["gene"]))
    cols = [n.replace("\n", " ").replace(" ", "\n") for n, _, _ in groups]
    M_r = det.pivot(index="gene", columns="group", values="rna_detection_pct").loc[genes, [c.replace("\n", " ") for c in cols]].to_numpy()
    M_a = det.pivot(index="gene", columns="group", values="atac_detection_pct").loc[genes, [c.replace("\n", " ") for c in cols]].to_numpy()
    heat(axf1, M_r, genes, cols, f"f  RNA detected,\n% at {budget['rna']:,} UMI", vmax=100)
    heat(axf2, M_a, genes, cols, f"promoter peak detected,\n% at {budget['atac']:,} fragments", vmax=100)
    # panels g to i: the per-nucleus promoter reading the heatmap replaced, put
    # back beside it on the owner's decision of 2026-09-22 (DEVELOPMENT decision
    # 28), with the ATAC depth that dominates it. Not a test; the lesson of C127.
    depth = a.obs["atac_counts"].to_numpy()
    violins(axes[2, 0], [(n, atac[m], c) for n, m, c in groups],
            "log1p promoter counts per 10,000 ATAC counts")
    axes[2, 0].set_title("g  Promoter chromatin, per nucleus, by group")
    violins(axes[2, 1], [(n, np.log10(np.maximum(depth[m], 1)), c) for n, m, c in groups],
            "log10 ATAC counts per nucleus")
    axes[2, 1].set_title("h  ATAC depth, by group")
    violins(axes[2, 2], [(n, atac[m][atac[m] > 0], c) for n, m, c in groups],
            "log1p promoter counts per 10,000 ATAC counts")
    axes[2, 2].set_title("i  Promoter chromatin, nuclei with any signal")
    fig.suptitle("A1  Closed or silenced: the AT2 identity programme in RNA and in chromatin, "
                 "GSE310539 wildtype wells", x=0.01, ha="left", fontsize=12.5, fontweight="semibold")
    fig.tight_layout(w_pad=2.5, rect=(0, 0, 1, 1))
    save(fig, "rq_a1_chromatin.png")
    table = pd.DataFrame([{
        "group": n.replace("\n", " "), "n_nuclei": int(m.sum()),
        "median_rna_score": float(np.median(rna[m])), "median_atac_promoter": float(np.median(atac[m])),
        "median_atac_counts": float(np.median(depth[m])),
        "pct_no_promoter_fragment": float(100 * np.mean(atac[m] == 0)),
        "median_atac_promoter_nonzero": float(np.median(atac[m][atac[m] > 0])),
    } for n, m, _ in groups])
    table.to_csv(OUT / "rq_a1_groups.csv", index=False)
    det.to_csv(OUT / "rq_a1_detection_at_budget.csv", index=False)
    mean_det = det.groupby("group")[["rna_detection_pct", "atac_detection_pct"]].mean()
    note("A1", wells={w: int((well == w).sum()) for w in colours}, transitional_pct=frac,
         promoter_peaks=int(a.uns["n_promoter_peaks"]), promoter_genes=int(len(genes)),
         budget=budget, groups=table.to_dict("records"), mean_detection=mean_det.to_dict("index"))
    return {"table": table, "frac": frac, "n": {w: int((well == w).sum()) for w in colours},
            "budget": budget, "mean_det": mean_det, "n_genes": len(genes)}


def figure_a4(a, counts: pd.DataFrame) -> dict:
    xy = a.obsm["X_umap"]
    fig = plt.figure(figsize=(13.5, 8.2))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[0, 1])
    axv = fig.add_subplot(gs[0, 2:]); axc = [fig.add_subplot(gs[1, i]) for i in range(4)]
    for ax, g, letter in ((ax1, "Axin2", "a"), (ax2, "Il1r1", "b")):
        v = np.asarray(a[:, g].X.todense()).ravel()
        feature(ax, xy, v, f"{letter}  {g}, log-normalised", vmax=max(float(v.max()), 0.1),
                label="log1p CP10k")
        ax.text(0.02, 0.98, f"{100 * np.mean(v > 0):.1f}% of nuclei", transform=ax.transAxes,
                va="top", fontsize=8.5, color=vs.INK_2)
    wells = list(GSE310539["libraries"].values())
    names = [w["name"] for w in wells]
    well_colour = {w: vs.SLOT[i + 1] for i, w in enumerate(names)}
    at2 = counts[counts["at2_nucleus"]]
    rows = []
    groups = []
    for g in ("Axin2", "Il1r1"):
        for w in names:
            sub = at2[at2["well"] == w]
            v = np.log1p(1e4 * sub[g].to_numpy() / np.maximum(sub["rna_counts"].to_numpy(), 1))
            groups.append((w.replace("_", "\n"), v, well_colour[w]))
            rows.append({"gene": g, "well": w, "n_at2_nuclei": len(sub), "pct_detected": 100 * np.mean(sub[g] > 0),
                         "mean_count_in_positive": float(sub.loc[sub[g] > 0, g].mean()) if (sub[g] > 0).any() else 0.0})
    violins(axv, groups, "log1p CP10k in AT2 nuclei", annotate=fraction_label, points=True)
    axv.axvline(4.5, color=vs.GRID, lw=1)
    axv.text(0.25, -0.2, "Axin2", ha="center", va="top", fontsize=9.5, color=vs.INK, style="italic",
             transform=axv.transAxes)
    axv.text(0.75, -0.2, "Il1r1", ha="center", va="top", fontsize=9.5, color=vs.INK, style="italic",
             transform=axv.transAxes)
    axv.set_title("c  Both transcripts in AT2 nuclei (Sftpc-positive, not transitional): violins of all "
                  "nuclei, points are the positive ones", fontsize=10)
    axv.tick_params(axis="x", labelsize=7.5)
    tiles = []
    for k, (ax, w) in enumerate(zip(axc, names), start=1):
        sub = at2[at2["well"] == w]
        ax_det = sub["Axin2"] >= 1
        il_det = sub["Il1r1"] >= 1
        mat = np.array([[(~ax_det & ~il_det).sum(), (~ax_det & il_det).sum()],
                        [(ax_det & ~il_det).sum(), (ax_det & il_det).sum()]], dtype=float)
        pct = 100 * mat / mat.sum()
        ax.imshow(pct, cmap=RAMP, vmin=0, vmax=100)
        for r in range(2):
            for c in range(2):
                ax.text(c, r, f"{pct[r, c]:.1f}%\n({int(mat[r, c]):,})", ha="center", va="center",
                        fontsize=8.5, color=vs.INK if pct[r, c] < 60 else vs.SURFACE)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Il1r1 absent", "Il1r1 detected"], fontsize=8)
        ax.set_yticks([0, 1]); ax.set_yticklabels(["Axin2 absent", "Axin2 detected"], fontsize=8)
        ax.set_title(f"d{k}  {w.replace('_', ' ')}, n = {len(sub):,}", fontsize=9.5)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        tiles.append({"well": w, "n_at2_nuclei": len(sub), "both_pct": pct[1, 1], "axin2_only_pct": pct[1, 0],
                      "il1r1_only_pct": pct[0, 1], "neither_pct": pct[0, 0]})
    fig.suptitle("A4  Axin2 and Il1r1 in the same AT2 nuclei, GSE310539", x=0.01, ha="left",
                 fontsize=12.5, fontweight="semibold")
    fig.tight_layout()
    save(fig, "rq_a4_axin2_il1r1.png")
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "rq_a4_detection.csv", index=False)
    pd.DataFrame(tiles).to_csv(OUT / "rq_a4_codetection.csv", index=False)
    note("A4", detection=rows, codetection=tiles)
    return {"table": table, "tiles": pd.DataFrame(tiles)}


# ---------------------------------------------------------------------------
# A5: GSE247130 control wells
# ---------------------------------------------------------------------------
def build_gse247130(name: str):
    w = next(w for w in wells_corrected() if w["well"] == name)
    a, _, _, _ = load_multiome(w["matrix"], [w["suffix"]], None)
    a.obs["well"] = name
    return embed_rna(a, name)


def figure_a5(objs: dict) -> dict:
    names = ["P9_control", "7wk_control", "SeV_control"]
    titles = {"P9_control": "P9 (uninjured neonatal)", "7wk_control": "seven weeks (uninjured adult)",
              "SeV_control": "SeV infected adult"}
    fig = plt.figure(figsize=(14, 11))
    gs = fig.add_gridspec(3, 4, width_ratios=[1, 1, 1, 0.9])
    rows = []
    vmax = {g: max(float(np.asarray(objs[n][:, g].X.todense()).max()) for n in names) for g in ("Cldn4", "Krt8")}
    for r, n in enumerate(names):
        a = objs[n]
        xy = a.obsm["X_umap"]
        trans = a.obs["transitional_label"].to_numpy()
        highlight(fig.add_subplot(gs[r, 0]), xy, trans, vs.SLOT[2],
                  f"{'abc'[r]}  {titles[n]}, n = {a.n_obs:,}", f"{100 * trans.mean():.2f}% labelled")
        for c, g in enumerate(("Cldn4", "Krt8"), start=1):
            v = np.asarray(a[:, g].X.todense()).ravel()
            feature(fig.add_subplot(gs[r, c]), xy, v, f"{g}", vmax=vmax[g], label="log1p CP10k")
        rows.append({"well": n, "n_nuclei": int(a.n_obs), "transitional_pct": 100 * trans.mean(),
                     "cldn4_pct": 100 * float(np.mean(a.obs["raw_Cldn4"] >= 1)),
                     "krt8_pct": 100 * float(np.mean(a.obs["raw_Krt8"] >= 1))})
    for c, g in enumerate(("Cldn4", "Krt8")):
        ax = fig.add_subplot(gs[c, 3])
        groups = [(titles[n].split(" (")[0], np.asarray(objs[n][:, g].X.todense()).ravel(), vs.SLOT[1 + i])
                  for i, n in enumerate(names)]
        violins(ax, groups, f"{g}, log1p CP10k", annotate=fraction_label, points=True)
        ax.set_title(f"{'de'[c]}  {g} across the three wells")
        ax.tick_params(axis="x", labelsize=8)
    ax = fig.add_subplot(gs[2, 3])
    ax.axis("off")
    ax.text(0, 0.95, "Label: Cldn4 and Krt8 both detected\nat the depth available (display rule).\n"
                     "The registered comparison at one\ndepth budget is trial M1c, row C119.",
            va="top", fontsize=9, color=vs.INK_2)
    fig.suptitle("A5  A transitional marker set in development and after injury, GSE247130 control wells",
                 x=0.01, ha="left", fontsize=12.5, fontweight="semibold")
    fig.tight_layout()
    save(fig, "rq_a5_development.png")
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "rq_a5_wells.csv", index=False)
    note("A5", wells=rows)
    return {"table": table}


# ---------------------------------------------------------------------------
# A3: GSE262927
# ---------------------------------------------------------------------------
def build_capillary():
    import anndata as ad
    reg = importlib.import_module("06_regeneration_focus")
    src = MOUSE / "processed" / "final_clustered.h5ad"
    a = ad.read_h5ad(src, backed="r")
    obs = a.obs
    keep = np.array(obs["has_author_metadata"].astype(str).isin(["True", "true"]), dtype=bool)
    keep &= obs["leiden_cluster"].astype(str).isin(reg.CAPILLARY_CLUSTERS).to_numpy()
    keep &= ~obs["leiden_cluster"].astype(str).isin(reg.CAPILLARY_DROP).to_numpy()
    idx = np.flatnonzero(keep)
    X = a[idx].X
    X = X.to_memory() if hasattr(X, "to_memory") else sp.csr_matrix(X)
    sub = ad.AnnData(X=sp.csr_matrix(X), obs=obs.iloc[idx].copy(), var=a.var.copy())
    a.file.close()
    del a
    # a manual row take keeps every category of the full object, including the
    # eight unannotated samples with no cell here; an empty batch breaks the
    # per-batch HVG call, so drop unused categories as anndata's own subsetting does
    for col in sub.obs.columns:
        if isinstance(sub.obs[col].dtype, pd.CategoricalDtype):
            sub.obs[col] = sub.obs[col].cat.remove_unused_categories()
    log(f"  capillary: {sub.n_obs:,} cells from clusters {reg.CAPILLARY_CLUSTERS} minus {reg.CAPILLARY_DROP}")
    reg.embed(sub, n_hvg=2000, resolution=0.4, tag="capillary")
    reg.score_module(sub, reg.ICAP_MARKERS, "iCAP_score")
    sub.obs["day"] = sub.obs["sacrifice_day"].astype(float).astype(int)
    # the backed read hands back pandas nullable strings, which the h5ad writer
    # refuses by default; plain object columns and index write everywhere
    sub.obs.index = pd.Index(np.asarray(sub.obs.index.astype(str), dtype=object))
    for k in list(sub.obs.columns):
        if not pd.api.types.is_numeric_dtype(sub.obs[k]) and not pd.api.types.is_bool_dtype(sub.obs[k]):
            sub.obs[k] = pd.Series(np.asarray(sub.obs[k].astype(str), dtype=object), index=sub.obs.index)
    sub.var.index = pd.Index(np.asarray(sub.var.index.astype(str), dtype=object))
    return sub


def figure_a3(cap) -> dict:
    my = pd.read_csv(MOUSE / "myeloid_focus" / "tables" / "myeloid_cell_metadata.csv")
    icap = pd.read_csv(MOUSE / "regeneration_focus" / "tables" / "icap_abundance_per_sample.csv")
    fig = plt.figure(figsize=(15, 10.5))
    gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 0.8])
    states = {"aMAC": (vs.SLOT[1], "alveolar macrophage"), "iMAC": (vs.SLOT[2], "interstitial macrophage"),
              "iMON": (vs.SLOT[3], "inflammatory monocyte")}
    xy = my[["UMAP_1", "UMAP_2"]].to_numpy()
    rows = []
    for c, (phase, days) in enumerate(PHASES):
        ax = fig.add_subplot(gs[0, c])
        m = my["day"].isin(days).to_numpy()
        ax.scatter(xy[:, 0], xy[:, 1], c=vs.DEEMPH, **DOT)
        for lab, (colour, _) in states.items():
            mm = m & (my["label"] == lab).to_numpy()
            ax.scatter(xy[mm, 0], xy[mm, 1], c=colour, **DOT)
            rows.append({"phase": phase, "label": lab, "cells": int(mm.sum()), "of_phase_pct": 100 * mm.sum() / max(m.sum(), 1)})
        umap_axes(ax, f"{'abcd'[c]}  {PHASE_DISPLAY[phase]} (n = {m.sum():,})")
    handles = [plt.Line2D([], [], marker="o", ls="", color=colour, label=f"{lab}, {desc}")
               for lab, (colour, desc) in states.items()]
    fig.axes[0].legend(handles=handles, loc="lower left", fontsize=8, handletextpad=0.3)
    cxy = cap.obsm["X_umap"]
    score = cap.obs["iCAP_score"].to_numpy()
    vmax = float(np.percentile(score, 99))
    for c, (phase, days) in enumerate(PHASES):
        ax = fig.add_subplot(gs[1, c])
        m = cap.obs["day"].isin(days).to_numpy()
        ax.scatter(cxy[:, 0], cxy[:, 1], c=vs.DEEMPH, **DOT)
        order = np.flatnonzero(m)[np.argsort(score[m])]
        sm = ax.scatter(cxy[order, 0], cxy[order, 1], c=score[order], cmap=RAMP, vmin=0, vmax=vmax, **DOT)
        umap_axes(ax, f"{'efgh'[c]}  capillary, {PHASE_DISPLAY[phase]} (n = {m.sum():,})")
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.outline.set_visible(False); cb.ax.tick_params(length=0, labelsize=8, colors=vs.MUTED)
    cb.set_label("iCAP score (Sparcl1, Ntrk2)", color=vs.INK_2, fontsize=8)
    ax = fig.add_subplot(gs[2, :])
    days = sorted(icap["day"].unique())
    med = icap.groupby("day")["pct"].median()
    ax.plot(range(len(days)), [med[d] for d in days], color=vs.SLOT[1], lw=2, zorder=1)
    for i, d in enumerate(days):
        v = icap.loc[icap["day"] == d, "pct"]
        ax.scatter(np.full(len(v), i) + np.linspace(-0.12, 0.12, len(v)), v, s=28, color=vs.SLOT[1],
                   edgecolor=vs.SURFACE, linewidth=0.8, zorder=2)
        ax.text(i, med[d] + 2.5, f"{med[d]:.1f}", ha="center", fontsize=8, color=vs.INK_2)
    ax.set_xticks(range(len(days))); ax.set_xticklabels([f"{int(d)} dpi" for d in days])
    ax.set_ylabel("iCAP, % of capillary cells")
    ax.set_title("i  Injury-induced capillary state per animal (points) and the median per day (line)")
    vs.recessive(ax, grid_axis="y")
    fig.suptitle("A3  Myeloid and capillary states across sampling times after infection, GSE262927",
                 x=0.01, ha="left", fontsize=12.5, fontweight="semibold")
    fig.tight_layout()
    save(fig, "rq_a3_persistence.png")
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "rq_a3_myeloid_by_phase.csv", index=False)
    pd.DataFrame({"day": days, "median_icap_pct": [med[d] for d in days],
                  "animals": [int((icap["day"] == d).sum()) for d in days]}).to_csv(OUT / "rq_a3_icap_by_day.csv", index=False)
    note("A3", myeloid_cells=int(len(my)), capillary_cells=int(cap.n_obs),
         icap_median_by_day={str(int(d)): float(med[d]) for d in days},
         animals_by_day={str(int(d)): int((icap["day"] == d).sum()) for d in days})
    return {"med": med, "days": days, "n_my": len(my), "n_cap": int(cap.n_obs)}


# ---------------------------------------------------------------------------
def main() -> int:
    global OUT, CACHE, REBUILD
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--figures", nargs="+", choices=ACTIVE_FIGURES, required=True)
    mode = ap.add_mutually_exclusive_group(required=True)
    mode.add_argument("--replot", action="store_true", help="require existing display caches; never refit embeddings")
    mode.add_argument("--rebuild-embeddings", action="store_true", help="explicitly prepare fresh display embeddings")
    args = ap.parse_args()
    selected = list(dict.fromkeys(args.figures))
    REBUILD = args.rebuild_embeddings
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    OUT = REPO / "analysis" / "figures" / "rq" / "renders" / run_id
    OUT.mkdir(parents=True, exist_ok=False)
    if REBUILD:
        CACHE = OUT / "processed"
    RECORD["selected_figures"] = selected
    RECORD["render_mode"] = "fresh_embeddings" if REBUILD else "existing_embeddings"
    RECORD["script_sha256"] = sha256(Path(__file__))
    RECORD["palette_sha256"] = sha256(REPO / "analysis/config/palette.json")
    RECORD["caption_owner"] = "analysis/figures/rq/README.md (curated); panel facts in this record"
    vs.apply(plt)
    if {"A1", "A4"}.intersection(selected):
        wt = cached("gse310539_wildtype", build_gse310539)
        if "A1" in selected:
            ensure_promoter_score(wt)
            cache_path = CACHE / "gse310539_wildtype.h5ad"
            RECORD["input_sha256"][cache_path.relative_to(REPO).as_posix()] = sha256(cache_path)
            figure_a1(wt)
        if "A4" in selected:
            figure_a4(wt, all_well_counts_gse310539())
        del wt
    if "A5" in selected:
        objs = {n: cached(n, lambda n=n: build_gse247130(n))
                for n in ("P9_control", "7wk_control", "SeV_control")}
        figure_a5(objs)
        del objs
    if "A3" in selected:
        for name in ("myeloid_focus/tables/myeloid_cell_metadata.csv",
                     "regeneration_focus/tables/icap_abundance_per_sample.csv"):
            p = MOUSE / name
            RECORD.setdefault("input_sha256", {})[p.relative_to(REPO).as_posix()] = sha256(p)
        cap = cached("capillary", build_capillary, backed="r")
        try:
            figure_a3(cap)
        finally:
            if getattr(cap, "isbacked", False):
                cap.file.close()
        del cap
    from importlib.metadata import version
    RECORD["versions"] = {p: version(p) for p in ("scanpy", "anndata", "numpy", "pandas", "matplotlib", "umap-learn")}
    RECORD["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    RECORD["output_sha256"] = {p.name: sha256(p) for p in sorted(OUT.iterdir()) if p.is_file()}
    (OUT / "run_record.json").write_text(json.dumps(RECORD, indent=2) + "\n", encoding="utf-8")
    log(f"done: {OUT.relative_to(REPO).as_posix()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
