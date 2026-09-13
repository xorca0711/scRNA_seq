#!/usr/bin/env python
"""Create a Fig. 1c-oriented view of the GSE178360 epithelial UMAP.

The saved epithelial embedding is never recomputed.  This script removes
off-compartment carryover cells and applies one rigid display rotation so the
major epithelial landmarks read in the same direction as Fig. 1c of Murthy
et al. (Nature 2022, doi:10.1038/s41586-022-04541-3).

The labels remain explicitly marked as candidates because they were inferred
independently and do not reproduce the authors' 18 Seurat clusters one-to-one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd


REFERENCE_DOI = "10.1038/s41586-022-04541-3"
ROTATION_DEGREES_CCW = 165.510271
PRIMARY_GENE_COLORS = {
    "KRT8": "#31688e",
    "CLDN4": "#1f9e89",
    "KRT17": "#d95f0e",
    "SFN": "#7b3294",
}
PRIMARY_GENES = ("KRT8", "CLDN4", "KRT17", "SFN")
REFERENCE_MARKER_GENES = (
    "SFTPC", "SCGB3A2", "SFTPB", "SCGB1A1", "KRT5", "TP63",
    "FOXJ1", "AGER", "KRT19",
)
OFF_COMPARTMENT_CONTROL_GENES = (
    "CD3D", "CD4", "CD8A", "FOXP3", "LYZ", "COL1A1", "CTHRC1", "PECAM1",
)

IDENTITIES = {
    "Ciliated (candidate)": {
        "short": "Ciliated",
        "markers": "FOXJ1+ TPPP3+ CAPS+",
        "color": "#f06f6a",
        "label_xy": (-11.8, 2.1),
        "ha": "left",
    },
    "Neuroendocrine (candidate)": {
        "short": "Neuroendocrine",
        "markers": "ASCL1+ CALCA+ CHGA+",
        "color": "#ec5aa2",
        "label_xy": (-2.0, 7.9),
        "ha": "center",
    },
    "Transitional (candidate)": {
        "short": "KRT8 transitional",
        "markers": "KRT8+ CLDN4+ KRT17+ SFN+",
        "color": "#d34e9f",
        "label_xy": (5.2, 10.4),
        "ha": "center",
    },
    "Basal (candidate)": {
        "short": "Basal / distal-basal-like",
        "markers": "KRT5+ TP63+ KRT17+",
        "color": "#00a99d",
        "label_xy": (13.2, 4.3),
        "ha": "left",
    },
    "Club (candidate)": {
        "short": "Secretory / club",
        "markers": "SCGB1A1+ SCGB3A2+ SFTPB+",
        "color": "#ef8a24",
        "label_xy": (12.1, -1.7),
        "ha": "left",
    },
    "AT2 (candidate)": {
        "short": "AT2",
        "markers": "SFTPC+ SFTPB+ NAPSA+",
        "color": "#91bd2c",
        "label_xy": (2.1, -11.9),
        "ha": "left",
    },
    "AT1 (candidate)": {
        "short": "AT1",
        "markers": "AGER+ CAV1+ HOPX+",
        "color": "#16a34a",
        "label_xy": (7.2, -12.5),
        "ha": "left",
    },
}

# Cluster-level analogues are deliberately more granular than the broad
# independent identities above.  A reference label is used verbatim only when
# its defining marker pattern is present; otherwise the data-supported label
# says "-like" or describes what is actually resolved.
CLUSTER_ANNOTATIONS = {
    "6": {
        "label": "Ciliated cells (CC)\nFOXJ1+ CAPS+ (candidate; L6)",
        "color": "#f8766d", "label_xy": (-14.2, 8.0), "ha": "left",
    },
    "9": {
        "label": "Ciliated-low state\nCAPS+ TPPP3+; FOXJ1low (candidate; L9)",
        "color": "#e76f51", "label_xy": (-13.9, -6.6), "ha": "left",
    },
    "11": {
        "label": "Motile ciliated\nFOXJ1+ DNAH+ (candidate; L11)",
        "color": "#ff8a80", "label_xy": (-14.2, 4.2), "ha": "left",
    },
    "12": {
        "label": "FOXJ1+SCGB3A2+ (SCGB3A2-CC)\n(candidate analogue; L12)",
        "color": "#c96f45", "label_xy": (-14.2, 0.0), "ha": "left",
    },
    "13": {
        "label": "Neuroendocrine cells\nASCL1+ CALCA+ CHGA+ (candidate; L13)",
        "color": "#ec5aa2", "label_xy": (-2.2, 8.6), "ha": "center",
    },
    "1": {
        "label": "SFTPB+KRT5+ basal cells (BC)\n(candidate analogue; L1)",
        "color": "#20bfa9", "label_xy": (3.8, 11.3), "ha": "center",
    },
    "5": {
        "label": "SFTPB+KRT5low (Distal-BC-1)\n(candidate analogue; L5)",
        "color": "#00a6b2", "label_xy": (13.6, 9.1), "ha": "left",
    },
    "15": {
        "label": "Differentiating basal cells\nKRT17+ TP63+ KRT5low (candidate; L15)",
        "color": "#6c9eeb", "label_xy": (14.2, 4.2), "ha": "left",
    },
    "3": {
        "label": "SFTPB+SCGB3A2+SCGB1A1low\n(TRB-SC-like candidate; L3)",
        "color": "#f28e2b", "label_xy": (14.0, -1.3), "ha": "left",
    },
    "16": {
        "label": "SFTPB+SCGB3A2+SCGB1A1+\n(pre-TRB-SC analogue; L16)",
        "color": "#eea43b", "label_xy": (-3.2, -3.9), "ha": "center",
    },
    "17": {
        "label": "SCGB1A1+ secretory cells\nSCGB3A2+ SFTPB+ (candidate; L17)",
        "color": "#d9a300", "label_xy": (4.1, -4.1), "ha": "center",
    },
    "4": {
        "label": "SFTPC+SCGB3A2+ (AT0)\n(candidate analogue; L4)",
        "color": "#00bfc4", "label_xy": (6.8, -6.9), "ha": "left",
    },
    "0": {
        "label": "AT2\nSFTPC+ SFTPB+ NAPSA+ (candidate; L0)",
        "color": "#94b72b", "label_xy": (-3.6, -11.3), "ha": "center",
    },
    "18": {
        "label": "AT1\nAGER+ CAV1+ HOPX+ (candidate; L18)",
        "color": "#00a65a", "label_xy": (8.2, -12.1), "ha": "left",
    },
}

UNRESOLVED_REFERENCE_STATES = (
    "Deuterosomal cells",
    "Proliferating airway cells",
    "SFTPB+KRT5- (Distal-BC-2)",
    "SCGB1A1+MUC5B+",
    "MUC5AC+MUC5B+",
    "SFTPB+SCGB3A2+SCGB1A1- (TRB-SC)",
    "Proliferating AT2",
    "Immature AT1",
)

COMPACT_LABELS = {
    "Ciliated (candidate)": "Ciliated",
    "Neuroendocrine (candidate)": "NE",
    "Transitional (candidate)": "KRT8-trans.",
    "Basal (candidate)": "Basal",
    "Club (candidate)": "Secretory",
    "AT1 (candidate)": "AT1",
}

COMPACT_OFFSETS = {
    "Ciliated (candidate)": (-0.2, 0.9),
    "Neuroendocrine (candidate)": (0.0, 0.9),
    "Transitional (candidate)": (0.0, 1.0),
    "Basal (candidate)": (0.0, 1.0),
    "Club (candidate)": (0.0, -1.1),
    "AT1 (candidate)": (0.3, -0.9),
}

CLUSTER_COMPACT_LABELS = {
    "4": ("AT0", (0.9, 0.4)),
    "0": ("AT2", (-0.3, -1.0)),
}

EXPRESSION_CMAP = LinearSegmentedColormap.from_list(
    "gray_viridis", ["#dedede", "#39568c", "#1f9e89", "#fde725"]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("analysis/GSE178360/epithelial_subanalysis/epithelial_clustered.h5ad"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "analysis/GSE178360/epithelial_subanalysis/figures/reference_aligned"
        ),
    )
    return parser.parse_args()


def reference_oriented_coordinates(umap: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mean-center and rotate row-vector coordinates; return xy and center."""
    center = np.asarray(umap, dtype=float).mean(axis=0)
    theta = np.deg2rad(ROTATION_DEGREES_CCW)
    rotation = np.array(
        [[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]]
    )
    return (np.asarray(umap, dtype=float) - center) @ rotation, center


def clean_axis(ax: plt.Axes) -> None:
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="box")


def add_display_axes(
    ax: plt.Axes,
    xy: np.ndarray,
    origin: tuple[float, float] | None = None,
) -> None:
    xmin, ymin = xy.min(axis=0)
    xspan, yspan = np.ptp(xy, axis=0)
    if origin is None:
        x0 = xmin + 0.035 * xspan
        y0 = ymin + 0.035 * yspan
    else:
        x0, y0 = origin
    ax.annotate("", (x0 + 0.085 * xspan, y0), (x0, y0),
                arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "#333333"})
    ax.annotate("", (x0, y0 + 0.085 * yspan), (x0, y0),
                arrowprops={"arrowstyle": "->", "lw": 1.0, "color": "#333333"})
    ax.text(x0 + 0.09 * xspan, y0, "display UMAP2", fontsize=7, va="center")
    ax.text(x0, y0 + 0.095 * yspan, "display UMAP1", fontsize=7,
            ha="center", va="bottom", rotation=90)


def identity_centers(epi: ad.AnnData, xy: np.ndarray) -> dict[str, np.ndarray]:
    values = epi.obs["proposed_cell_type"].astype(str).to_numpy()
    return {identity: np.median(xy[values == identity], axis=0) for identity in IDENTITIES}


def add_cluster_numbers(ax: plt.Axes, epi: ad.AnnData, xy: np.ndarray) -> None:
    clusters = epi.obs["leiden_cluster"].astype(str).to_numpy()
    for cluster in sorted(np.unique(clusters), key=lambda value: int(value)):
        center = np.median(xy[clusters == cluster], axis=0)
        ax.text(
            *center,
            cluster,
            ha="center",
            va="center",
            fontsize=7,
            weight="bold",
            color="#222222",
            bbox={"boxstyle": "circle,pad=0.22", "fc": "white", "ec": "#444444",
                  "lw": 0.65, "alpha": 0.92},
            zorder=5,
        )


def save_figure(fig: plt.Figure, output: Path, stem: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    fig.savefig(output / f"{stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(output / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_composition(epi: ad.AnnData, xy: np.ndarray, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(13.8, 9.5))
    clusters = epi.obs["leiden_cluster"].astype(str).to_numpy()
    for cluster, spec in CLUSTER_ANNOTATIONS.items():
        mask = clusters == cluster
        ax.scatter(
            xy[mask, 0], xy[mask, 1], s=3.8, c=spec["color"], linewidths=0,
            alpha=0.9, rasterized=True,
        )

    for cluster, spec in CLUSTER_ANNOTATIONS.items():
        center = np.median(xy[clusters == cluster], axis=0)
        ax.annotate(
            spec["label"],
            xy=center,
            xytext=spec["label_xy"],
            ha=spec["ha"],
            va="center",
            fontsize=8.5,
            color="#222222",
            linespacing=1.3,
            arrowprops={"arrowstyle": "-", "lw": 0.85, "color": spec["color"],
                        "shrinkA": 3, "shrinkB": 3},
            bbox={"boxstyle": "round,pad=0.25", "fc": "white", "ec": "none",
                  "alpha": 0.86},
            zorder=6,
        )

    # Keep the display-axis key in the empty far-left margin rather than over
    # the ciliated morphology.
    add_display_axes(ax, xy, origin=(-13.9, -12.2))
    clean_axis(ax)
    ax.set_xlim(-14.8, 22.2)
    ax.set_ylim(-13.4, 12.2)
    ax.set_title(
        "GSE178360 epithelial states - Fig. 1c-oriented display",
        fontsize=16, weight="bold", pad=33,
    )
    ax.text(
        0.5, 1.012,
        "6,386 epithelial cells | all 14 retained Leiden regions annotated | literal Fig. 1c notation used only where marker-supported",
        transform=ax.transAxes, ha="center", va="bottom", fontsize=9, color="#4a4a4a",
    )
    ax.text(
        21.6, -4.5,
        "Fig. 1c states not separately resolved here:\n- "
        + "\n- ".join(UNRESOLVED_REFERENCE_STATES),
        ha="right", va="top", fontsize=7.1,
        color="#4a4a4a", linespacing=1.35,
        bbox={"boxstyle": "round,pad=0.45", "fc": "#f7f7f7", "ec": "#bdbdbd",
              "lw": 0.7, "alpha": 0.95},
    )
    ax.text(
        0.5, -0.018,
        "Rigid display rotation only (165.51 degrees CCW); neighborhood graph and UMAP geometry unchanged. "
        "Analogue means closest marker-supported correspondence, not reproduction of the authors' cluster. "
        "Reference: Murthy et al., Nature 2022, Fig. 1c.",
        transform=ax.transAxes, ha="center", va="top", fontsize=7.5, color="#555555",
    )
    save_figure(fig, output, "epithelial_UMAP_proposed_reference_aligned")


def expression_vector(epi: ad.AnnData, gene: str) -> np.ndarray:
    matrix = epi[:, gene].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix).reshape(-1)


def add_compact_identity_labels(
    ax: plt.Axes,
    centers: dict[str, np.ndarray],
    cluster_centers: dict[str, np.ndarray],
) -> None:
    for identity, center in centers.items():
        if identity == "AT2 (candidate)":
            continue
        offset = np.asarray(COMPACT_OFFSETS[identity])
        ax.annotate(
            COMPACT_LABELS[identity],
            xy=center,
            xytext=center + offset,
            ha="center",
            va="center",
            fontsize=6.7,
            color="#2f2f2f",
            arrowprops={"arrowstyle": "-", "lw": 0.55, "color": "#666666"},
            bbox={"boxstyle": "round,pad=0.16", "fc": "white", "ec": "none",
                  "alpha": 0.78},
            zorder=6,
        )
    for cluster, (label, offset) in CLUSTER_COMPACT_LABELS.items():
        center = cluster_centers[cluster]
        ax.annotate(
            label,
            xy=center,
            xytext=center + np.asarray(offset),
            ha="center", va="center", fontsize=6.7, color="#2f2f2f",
            arrowprops={"arrowstyle": "-", "lw": 0.55, "color": "#666666"},
            bbox={"boxstyle": "round,pad=0.16", "fc": "white", "ec": "none",
                  "alpha": 0.78}, zorder=6,
        )


def draw_expression_panel(
    ax: plt.Axes,
    epi: ad.AnnData,
    xy: np.ndarray,
    gene: str,
    centers: dict[str, np.ndarray],
    cluster_centers: dict[str, np.ndarray],
) -> None:
    values = expression_vector(epi, gene)
    positive = values > 0
    vmax = float(np.percentile(values[positive], 99)) if positive.any() else 1.0
    vmax = max(vmax, 1e-9)
    ax.scatter(xy[:, 0], xy[:, 1], s=2.2, c="#dedede", linewidths=0,
               rasterized=True)
    order = np.argsort(values)
    shown = order[values[order] > 0]
    points = ax.scatter(
        xy[shown, 0], xy[shown, 1], s=2.5, c=values[shown], cmap=EXPRESSION_CMAP,
        vmin=0, vmax=vmax, linewidths=0, rasterized=True,
    )
    add_compact_identity_labels(ax, centers, cluster_centers)
    clean_axis(ax)
    ax.set_xlim(-10.8, 12.3)
    ax.set_ylim(-12.2, 9.2)
    ax.set_title(gene, fontsize=13, weight="bold", pad=6)
    cbar = ax.figure.colorbar(points, ax=ax, fraction=0.035, pad=0.018)
    cbar.set_label("log1p(CP10K)", fontsize=7)
    cbar.ax.tick_params(labelsize=6, length=2)


def plot_feature_panel(epi: ad.AnnData, xy: np.ndarray, output: Path) -> None:
    centers = identity_centers(epi, xy)
    clusters = epi.obs["leiden_cluster"].astype(str).to_numpy()
    cluster_centers = {
        cluster: np.median(xy[clusters == cluster], axis=0)
        for cluster in CLUSTER_COMPACT_LABELS
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 10.2))
    for ax, gene in zip(axes.ravel(), PRIMARY_GENES, strict=True):
        draw_expression_panel(ax, epi, xy, gene, centers, cluster_centers)
    fig.suptitle(
        "Reference-oriented epithelial feature projections",
        fontsize=16, weight="bold", y=0.995,
    )
    fig.text(
        0.5, 0.967,
        "KRT8, CLDN4, KRT17 and SFN on the identical epithelial subset and display transform",
        ha="center", fontsize=9, color="#4a4a4a",
    )
    fig.text(
        0.5, 0.012,
        "Labels denote independently inferred candidate regions; expression is clipped at each gene's 99th percentile among positive cells.",
        ha="center", fontsize=7.5, color="#555555",
    )
    fig.subplots_adjust(left=0.03, right=0.97, bottom=0.04, top=0.94,
                        wspace=0.08, hspace=0.12)
    save_figure(fig, output, "epithelial_featureplots_KRT8_CLDN4_KRT17_SFN_reference_aligned")

    for gene in PRIMARY_GENES:
        fig, ax = plt.subplots(figsize=(8.2, 7.6))
        draw_expression_panel(ax, epi, xy, gene, centers, cluster_centers)
        ax.text(
            0.5, -0.01,
            "Fig. 1c-oriented display; candidate regions; rigid rotation only",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
            color="#555555",
        )
        save_figure(fig, output, f"epithelial_feature_{gene}_reference_aligned")

    plot_marker_collection(
        epi, xy, output, REFERENCE_MARKER_GENES,
        "Epithelial reference-marker projections",
        "epithelial_reference_marker_featureplots_reference_aligned",
        centers, cluster_centers,
    )
    plot_marker_collection(
        epi, xy, output, OFF_COMPARTMENT_CONTROL_GENES,
        "Off-compartment control projections",
        "epithelial_off_compartment_control_featureplots_reference_aligned",
        centers, cluster_centers,
    )


def plot_marker_collection(
    epi: ad.AnnData,
    xy: np.ndarray,
    output: Path,
    genes: tuple[str, ...],
    title: str,
    stem: str,
    centers: dict[str, np.ndarray],
    cluster_centers: dict[str, np.ndarray],
) -> None:
    present = [gene for gene in genes if gene in epi.var_names]
    ncols = 3
    nrows = int(np.ceil(len(present) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 4.4 * nrows), squeeze=False)
    for ax, gene in zip(axes.ravel(), present, strict=False):
        draw_expression_panel(ax, epi, xy, gene, centers, cluster_centers)
    for ax in axes.ravel()[len(present):]:
        ax.axis("off")
    fig.suptitle(title, fontsize=16, weight="bold", y=0.995)
    fig.text(
        0.5, 0.008,
        "Identical Fig. 1c-oriented epithelial subset and display transform; per-gene positive-cell p99 color limit.",
        ha="center", fontsize=7.5, color="#555555",
    )
    fig.subplots_adjust(left=0.025, right=0.98, bottom=0.025, top=0.955,
                        wspace=0.08, hspace=0.14)
    save_figure(fig, output, stem)

    for gene in present:
        if gene in PRIMARY_GENES:
            continue
        fig, ax = plt.subplots(figsize=(8.2, 7.6))
        draw_expression_panel(ax, epi, xy, gene, centers, cluster_centers)
        ax.text(
            0.5, -0.01,
            "Fig. 1c-oriented epithelial subset; candidate regions; rigid rotation only",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5,
            color="#555555",
        )
        save_figure(fig, output, f"epithelial_feature_{gene}_reference_aligned")


def primary_expression_summary(epi: ad.AnnData) -> pd.DataFrame:
    clusters = epi.obs['leiden_cluster'].astype(str).to_numpy()
    order = ('6', '11', '12', '9', '13', '1', '5', '15', '16', '17', '3', '4', '0', '18')
    labels = {
        '6': 'L6  Ciliated (CC)', '11': 'L11  Motile ciliated',
        '12': 'L12  SCGB3A2-CC-like', '9': 'L9  Ciliated-low',
        '13': 'L13  Neuroendocrine', '1': 'L1  Basal (BC-like)',
        '5': 'L5  Distal-BC-1-like', '15': 'L15  Differentiating basal',
        '16': 'L16  pre-TRB-SC-like', '17': 'L17  Secretory',
        '3': 'L3  TRB-SC-like', '4': 'L4  AT0 candidate',
        '0': 'L0  AT2 candidate', '18': 'L18  AT1 candidate',
    }
    rows = []
    for cluster in order:
        mask = clusters == cluster
        for gene in PRIMARY_GENES:
            values = expression_vector(epi, gene)[mask]
            rows.append({
                'leiden_cluster': cluster,
                'cluster_label': labels[cluster],
                'n_cells': int(mask.sum()),
                'gene': gene,
                'mean_log1p_cp10k': float(values.mean()),
                'median_log1p_cp10k': float(np.median(values)),
                'q25_log1p_cp10k': float(np.quantile(values, 0.25)),
                'q75_log1p_cp10k': float(np.quantile(values, 0.75)),
                'q95_log1p_cp10k': float(np.quantile(values, 0.95)),
                'pct_expressing': float(100 * np.mean(values > 0)),
            })
    return pd.DataFrame(rows)


def plot_primary_dotplot(summary: pd.DataFrame, output: Path) -> None:
    order = list(dict.fromkeys(summary['cluster_label']))
    fig, ax = plt.subplots(figsize=(8.8, 8.2))
    for gene_index, gene in enumerate(PRIMARY_GENES):
        gene_data = summary[summary['gene'] == gene].set_index('cluster_label').loc[order]
        scatter = ax.scatter(
            np.full(len(order), gene_index), np.arange(len(order)),
            s=18 + 4.6 * gene_data['pct_expressing'].to_numpy(),
            c=gene_data['mean_log1p_cp10k'].to_numpy(), cmap='viridis',
            vmin=0, vmax=float(summary['mean_log1p_cp10k'].max()),
            edgecolors='#333333', linewidths=0.35,
        )
    ax.set_xticks(range(len(PRIMARY_GENES)), PRIMARY_GENES, fontsize=11, weight='bold')
    ax.set_yticks(range(len(order)), order, fontsize=8.7)
    ax.invert_yaxis()
    ax.grid(color='#e4e8eb', linewidth=0.7)
    ax.set_axisbelow(True)
    ax.set_title('Primary marker expression by retained epithelial region', fontsize=15, weight='bold', pad=30)
    ax.text(
        0.5, 1.025, 'Dot colour: mean log1p(CP10K) | dot size: cells with expression > 0',
        transform=ax.transAxes, ha='center', fontsize=9, color='#555555',
    )
    colorbar = fig.colorbar(scatter, ax=ax, pad=0.03, fraction=0.04)
    colorbar.set_label('Mean log1p(CP10K)', fontsize=9)
    handles = [
        ax.scatter([], [], s=18 + 4.6 * pct, color='#6a6a6a', edgecolors='#333333', linewidths=0.35)
        for pct in (25, 50, 75, 100)
    ]
    ax.legend(
        handles, ['25%', '50%', '75%', '100%'], title='Cells expressing',
        frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.12),
        fontsize=8, title_fontsize=8,
    )
    fig.subplots_adjust(left=0.31, right=0.91, top=0.88, bottom=0.15)
    save_figure(fig, output, 'epithelial_primary_markers_dotplot_reference_aligned')


def plot_primary_violins(epi: ad.AnnData, output: Path) -> None:
    order = ('6', '11', '12', '9', '13', '1', '5', '15', '16', '17', '3', '4', '0', '18')
    labels = ('L6 Ciliated', 'L11 Motile cil.', 'L12 SCGB3A2-CC', 'L9 Ciliated-low',
              'L13 Neuroendo.', 'L1 Basal', 'L5 Distal-BC-1', 'L15 Diff. basal',
              'L16 pre-TRB-SC', 'L17 Secretory', 'L3 TRB-SC', 'L4 AT0',
              'L0 AT2', 'L18 AT1')
    clusters = epi.obs['leiden_cluster'].astype(str).to_numpy()
    fig, axes = plt.subplots(2, 2, figsize=(15.2, 10.4), sharex=True)
    for ax, gene in zip(axes.ravel(), PRIMARY_GENES, strict=True):
        vector = expression_vector(epi, gene)
        distributions = [vector[clusters == cluster] for cluster in order]
        parts = ax.violinplot(
            distributions, positions=np.arange(len(order)), widths=0.82,
            showmeans=False, showmedians=True, showextrema=False, points=80,
        )
        for body in parts['bodies']:
            body.set_facecolor(PRIMARY_GENE_COLORS[gene])
            body.set_edgecolor('#333333')
            body.set_alpha(0.8)
            body.set_linewidth(0.5)
        parts['cmedians'].set_color('#111111')
        parts['cmedians'].set_linewidth(0.9)
        ax.set_title(gene, color=PRIMARY_GENE_COLORS[gene], fontsize=14, weight='bold')
        ax.set_ylabel('log1p(CP10K)', fontsize=9)
        ax.set_ylim(0, max(0.25, float(np.quantile(vector, 0.995)) * 1.05))
        ax.grid(axis='y', color='#e5e7eb', linewidth=0.7)
        ax.set_axisbelow(True)
    for ax in axes[1]:
        ax.set_xticks(np.arange(len(order)), labels, rotation=52, ha='right', fontsize=7.4)
    fig.suptitle('Primary marker distributions across epithelial Leiden regions', fontsize=16, weight='bold', y=0.985)
    fig.text(
        0.5, 0.942,
        'All cells are included, including zeros; violin width is density and the dark bar is the median.',
        ha='center', fontsize=9, color='#555555',
    )
    fig.subplots_adjust(left=0.06, right=0.985, top=0.90, bottom=0.17, wspace=0.14, hspace=0.17)
    save_figure(fig, output, 'epithelial_primary_markers_violin_reference_aligned')


def plot_primary_expression_summaries(epi: ad.AnnData, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    summary = primary_expression_summary(epi)
    summary.to_csv(output / 'epithelial_primary_markers_by_cluster.csv', index=False)
    plot_primary_dotplot(summary, output)
    plot_primary_violins(epi, output)


def write_metadata(
    source: ad.AnnData,
    epi: ad.AnnData,
    center: np.ndarray,
    output: Path,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    identities = epi.obs["proposed_cell_type"].astype(str)
    clusters = epi.obs["leiden_cluster"].astype(str)
    rows = []
    for identity, spec in IDENTITIES.items():
        mask = identities == identity
        rows.append(
            {
                "candidate_identity": identity,
                "n_cells": int(mask.sum()),
                "retained_leiden_clusters": ",".join(
                    sorted(clusters[mask].unique(), key=int)
                ),
                "marker_notation": spec["markers"],
            }
        )
    pd.DataFrame(rows).to_csv(output / "reference_aligned_identity_summary.csv", index=False)

    metadata = {
        "dataset": "GSE178360",
        "source_object": str(Path("analysis/GSE178360/epithelial_subanalysis/epithelial_clustered.h5ad")),
        "source_n_cells": int(source.n_obs),
        "retained_n_cells": int(epi.n_obs),
        "excluded_n_cells": int(source.n_obs - epi.n_obs),
        "subset_rule": "retain cells whose proposed_cell_type is one of the seven epithelial candidate identities",
        "excluded_carryover_compartments": ["immune", "endothelial", "plasma", "mesothelial"],
        "reference": {
            "citation": "Murthy et al., Human distal lung maps and lineage hierarchies reveal a bipotent progenitor, Nature (2022)",
            "doi": REFERENCE_DOI,
            "figure": "Fig. 1c",
            "pdf_page": 2,
        },
        "display_transform": {
            "operation_order": ["mean-center retained subset", "rotate"],
            "centering_vector_original_umap": center.tolist(),
            "rotation_degrees_counterclockwise": ROTATION_DEGREES_CCW,
            "reflection": False,
            "scaling": False,
            "recomputed_umap": False,
            "recomputed_neighbors": False,
        },
        "primary_feature_genes": list(PRIMARY_GENES),
        "epithelial_reference_marker_genes": list(REFERENCE_MARKER_GENES),
        "off_compartment_control_genes": list(OFF_COMPARTMENT_CONTROL_GENES),
        "cluster_level_candidate_analogues": {
            cluster: spec["label"].replace("\n", " | ")
            for cluster, spec in CLUSTER_ANNOTATIONS.items()
        },
        "reference_states_not_separately_resolved": list(UNRESOLVED_REFERENCE_STATES),
        "interpretation_guardrail": "Candidate identities are independently inferred and are not asserted to reproduce the authors' 18 clusters one-to-one.",
    }
    (output / "reference_alignment_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def main() -> None:
    args = parse_args()
    source = ad.read_h5ad(args.input)
    missing_obs = {"proposed_cell_type", "leiden_cluster"} - set(source.obs)
    if missing_obs:
        raise KeyError(f"Missing required obs columns: {sorted(missing_obs)}")
    missing_genes = sorted(set(PRIMARY_GENES) - set(source.var_names))
    if missing_genes:
        raise KeyError(f"Missing requested genes: {missing_genes}")

    keep = source.obs["proposed_cell_type"].astype(str).isin(IDENTITIES)
    epi = source[keep].copy()
    xy, center = reference_oriented_coordinates(epi.obsm["X_umap"])

    plot_primary_expression_summaries(epi, args.output)
    plot_composition(epi, xy, args.output)
    plot_feature_panel(epi, xy, args.output)
    write_metadata(source, epi, center, args.output)
    print(
        f"Wrote reference-aligned epithelial figures for {epi.n_obs:,} cells "
        f"({source.n_obs - epi.n_obs:,} off-compartment cells excluded) to {args.output}"
    )


if __name__ == "__main__":
    main()
