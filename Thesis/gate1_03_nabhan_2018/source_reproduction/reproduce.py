"""Bounded, deterministic GSE109444 processed-FPKM reproduction.

Reads the immutable GEO tar directly (no extraction) and separates its bulk control.
Run from any directory; raw inputs are resolved relative to the repository.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import tarfile

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[2]
RAW = ROOT / "raw_data" / "GSE109444"
BULK = "GSM2943054"
FPKM_COLUMNS = ["tracking_id", "class_code", "nearest_ref_id", "gene_id",
                "gene_short_name", "tss_id", "locus", "length", "coverage",
                "FPKM", "FPKM_lo", "FPKM_hi", "status"]
HOUSEKEEPING = ["Actb", "Gapdh", "Ubc", "Ppia"]
SOURCES = {
    "GSE109444_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109444/soft/GSE109444_family.soft.gz",
    "GSE109444_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109444/suppl/GSE109444_RAW.tar",
    "Nabhan_2018_paper_and_supplement.pdf": "https://desailab.stanford.edu/sites/g/files/sbiybj24296/files/media/file/sc1.pdf",
}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_soft(path):
    """Keep repeated SOFT fields rather than losing their earlier values."""
    series, samples = {}, {}
    current = series
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if line.startswith("^SAMPLE = "):
                accession = line.split(" = ", 1)[1]
                current = samples.setdefault(accession, {})
            elif line.startswith("!") and " = " in line:
                key, value = line.split(" = ", 1)
                current.setdefault(key.lstrip("!"), []).append(value)
    return series, samples


def read_tracking(compressed_bytes):
    # These deposited files have NO header; the first gene is a real observation.
    frame = pd.read_csv(io.BytesIO(gzip.decompress(compressed_bytes)), sep="\t",
                        names=FPKM_COLUMNS, header=None, dtype=str,
                        keep_default_na=False)
    if frame.shape[1] != 13 or frame.iloc[0]["tracking_id"] == "tracking_id":
        raise ValueError("Unexpected FPKM tracking schema")
    for name in ["FPKM", "FPKM_lo", "FPKM_hi"]:
        frame[name] = pd.to_numeric(frame[name], errors="raise")
        if not np.isfinite(frame[name]).all() or (frame[name] < 0).any():
            raise ValueError(f"Invalid values in {name}")
    return frame


def detection(values, threshold):
    return values > 0 if threshold == 0 else values >= threshold


def conditional_fraction(numerator, denominator):
    denominator_n = int(np.asarray(denominator, dtype=bool).sum())
    joint_n = int((np.asarray(numerator, dtype=bool) & denominator).sum())
    return joint_n, denominator_n, joint_n / denominator_n if denominator_n else None


def split_samples(matrix, metadata):
    bulk_ids = metadata.loc[metadata["is_bulk"], "accession"].tolist()
    if bulk_ids != [BULK]:
        raise ValueError(f"Unexpected bulk controls: {bulk_ids}")
    return matrix.drop(index=bulk_ids), matrix.loc[bulk_ids]


def save_csv(frame, name, index=False):
    frame.to_csv(BASE / "tables" / name, index=index, float_format="%.10g",
                 lineterminator="\n")


def make_figures(cells, counts, coexpression, spec, metadata):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "svg.hashsalt": "nabhan2018-source-reproduction-v1",
                         "axes.spines.top": False,
                         "axes.spines.right": False})
    order = cells.sort_values(["Wnt5a", "Pdgfra"], ascending=False,
                              kind="stable").index.tolist()
    order_table = metadata.set_index("accession").loc[order].reset_index()
    order_table.insert(0, "plot_position", np.arange(1, len(order) + 1))
    save_csv(order_table, "plot_cell_order.csv")
    figure_genes = spec["figure2d_genes"]
    context_genes = spec["additional_context_genes"]
    fig = plt.figure(figsize=(15, 13), facecolor="white")
    grid = fig.add_gridspec(2, 2, height_ratios=[25, 10], width_ratios=[47, 10],
                            hspace=0.33, wspace=0.15)
    cap = float(np.ceil(np.log2(1 + cells[figure_genes + context_genes].to_numpy().max())))
    norm = Normalize(0, cap)
    for row, genes in enumerate([figure_genes, context_genes]):
        ax, bar = fig.add_subplot(grid[row, 0]), fig.add_subplot(grid[row, 1])
        im = ax.imshow(np.log2(1 + cells.loc[order, genes].to_numpy().T),
                       aspect="auto", cmap="magma", norm=norm,
                       interpolation="nearest")
        ax.set_yticks(range(len(genes)), labels=genes)
        positions = [0, 9, 19, 29, 39, 46]
        ax.set_xticks(positions, labels=[str(x + 1) for x in positions])
        ax.tick_params(axis="both", length=0)
        ax.set_xlabel("Cell order (Wnt5a, then Pdgfra, descending; accession mapping in table)")
        ax.set_title("A  Published Figure 2D gene panel" if row == 0 else
                     "B  Additional QC and Wnt-processing context", loc="left", pad=12,
                     fontsize=12, fontweight="bold")
        selected = counts.loc[counts["threshold"] == ">=1"].set_index("gene").loc[genes]
        bar.barh(range(len(genes)), selected["n_detected"], color="#277DA8", height=0.72)
        for y, value in enumerate(selected["n_detected"]):
            bar.text(value + 0.8, y, str(int(value)), va="center", fontsize=8)
        bar.set_ylim(len(genes) - 0.5, -0.5)
        bar.set_xlim(0, 53)
        bar.set_yticks([])
        bar.set_xticks([0, 20, 40])
        bar.set_xlabel("Cells with FPKM >= 1")
        bar.set_title("Detection / 47", fontsize=11)
    fig.suptitle("Adult lung mesenchyme: deposited single-cell expression",
                 x=0.095, y=0.975, ha="left", fontsize=18, fontweight="bold")
    fig.text(0.095, 0.945, "GSE109444 | 47 cells; 200-cell bulk excluded | Descriptive reproduction, no animal-level inference",
             fontsize=11, color="#444444")
    fig.subplots_adjust(left=0.095, right=0.89, top=0.905, bottom=0.075)
    cax = fig.add_axes([0.915, 0.54, 0.016, 0.3])
    fig.colorbar(im, cax=cax, label="log2(1 + FPKM)")
    fig.text(0.095, 0.026, "FPKM >= 1 is a declared display convention. Original positivity cutoff and clustering order are not specified.\n"
             "Marker expression does not demonstrate cell identity, spatial proximity, active secretion or functional Wnt signaling.",
             fontsize=9, color="#444444")
    for ext in ["png", "svg"]:
        fig.savefig(BASE / "figures" / f"source_expression_panel.{ext}", dpi=180,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)

    focus = coexpression.loc[coexpression["estimand"] == "Pdgfra_given_Wnt5a"].copy()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.7), gridspec_kw={"width_ratios": [1.25, 1]})
    x = np.arange(len(focus))
    percentages = focus["fraction"] * 100
    axes[0].plot(x, percentages, color="#277DA8", marker="o", markersize=8, linewidth=2)
    axes[0].axhline(74, linestyle="--", color="#999999", linewidth=1)
    axes[0].text(3.12, 68, "Paper: 74%", va="top", fontsize=9, color="#666666")
    for pos, (_, row) in enumerate(focus.iterrows()):
        axes[0].text(pos, 100 * row["fraction"] + 4,
                     f'{int(row["joint_n"])}/{int(row["denominator_n"])}\n{100 * row["fraction"]:.1f}%',
                     ha="center", va="bottom", fontsize=10)
    axes[0].set(xticks=x, xticklabels=focus["threshold"], xlabel="FPKM positivity cutoff (both genes)",
                ylabel="Pdgfra-positive among Wnt5a-positive (%)", ylim=(0, 110), xlim=(-0.3, 4.05))
    axes[0].set_title("A  Conditional co-detection", loc="left", fontweight="bold")
    core_wnts = ["Wnt5a", "Wnt2", "Wnt2b", "Wnt4", "Wnt9a"]
    width = 0.18
    colors = ["#c3dce8", "#80b8d1", "#277DA8", "#143e5e"]
    for j, (label, color) in enumerate(zip([t["label"] for t in spec["thresholds"]], colors)):
        data = counts.loc[counts["threshold"] == label].set_index("gene").loc[core_wnts]
        axes[1].bar(np.arange(5) + (j - 1.5) * width, data["n_detected"], width,
                    color=color, label=label)
    axes[1].set(xticks=range(5), xticklabels=core_wnts, ylabel="Detected cells / 47", ylim=(0, 49))
    axes[1].legend(title="FPKM", frameon=False, ncol=2, loc="upper right", fontsize=8)
    axes[1].set_title("B  Wnt detection sensitivity", loc="left", fontweight="bold")
    fig.suptitle("Positivity thresholds change the numerical reproduction", x=0.07,
                 ha="left", fontsize=16, fontweight="bold")
    fig.text(0.07, 0.02, "Same 47 deposited cells at every threshold; no threshold selected to match 74%.\n"
             "The paper does not specify a Wnt/Pdgfra positivity cutoff. Cell fractions are descriptive, not independent-animal estimates.",
             fontsize=9, color="#444444")
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.2, top=0.83, wspace=0.33)
    for ext in ["png", "svg"]:
        fig.savefig(BASE / "figures" / f"detection_sensitivity.{ext}", dpi=180,
                    metadata={"Date": None} if ext == "svg" else None)
    plt.close(fig)


def main():
    spec = json.loads((BASE / "specification.json").read_text())
    for name in ["tables", "figures"]:
        (BASE / name).mkdir(exist_ok=True)
    genes = spec["figure2d_genes"] + spec["additional_context_genes"]
    series, samples = read_soft(RAW / "GSE109444_family.soft.gz")
    metadata_rows, gene_rows, member_records, schema_rows = [], [], [], []
    for accession, fields in sorted(samples.items()):
        source = fields["Sample_source_name_ch1"][0]
        metadata_rows.append({"accession": accession, "title": fields["Sample_title"][0],
                              "source_name": source, "is_bulk": "bulk" in source.lower(),
                              "animal_id": "not supplied",
                              "characteristics": "; ".join(fields.get("Sample_characteristics_ch1", []))})
    metadata = pd.DataFrame(metadata_rows)
    with tarfile.open(RAW / "GSE109444_RAW.tar", "r") as archive:
        for member in sorted(archive.getmembers(), key=lambda m: m.name):
            if not member.isfile():
                continue
            accession = member.name.split("_", 1)[0]
            if accession not in samples:
                raise ValueError(f"Archive sample without metadata: {accession}")
            compressed = archive.extractfile(member).read()
            frame = read_tracking(compressed)
            targets = frame[frame["gene_short_name"].isin(genes)].copy()
            present = targets["gene_short_name"].value_counts()
            if set(present.index) != set(genes) or (present != 1).any():
                raise ValueError(f"Missing or ambiguous target gene in {accession}: {present.to_dict()}")
            targets.insert(0, "accession", accession)
            gene_rows.append(targets[["accession", "gene_short_name", "tracking_id", "gene_id", "locus",
                                      "FPKM", "FPKM_lo", "FPKM_hi", "status"]])
            schema_rows.append({"accession": accession, "rows": len(frame),
                                "first_tracking_id": frame.iloc[0]["tracking_id"],
                                "target_genes": len(targets),
                                "target_status_not_OK": int((targets["status"] != "OK").sum())})
            member_records.append({"accession": accession, "member": member.name,
                                   "bytes": len(compressed), "sha256": hashlib.sha256(compressed).hexdigest()})
    long = pd.concat(gene_rows, ignore_index=True)
    matrix = long.pivot(index="accession", columns="gene_short_name", values="FPKM")[genes]
    cells, bulk = split_samples(matrix, metadata)
    assert len(cells) == 47 and len(bulk) == 1 and len(samples) == 48
    assert BULK not in cells.index and set(matrix.index) == set(samples)
    save_csv(metadata, "sample_metadata.csv")
    save_csv(pd.DataFrame(schema_rows), "parse_validation.csv")
    save_csv(long, "panel_gene_records.csv")
    save_csv(cells, "single_cell_fpkm.csv", index=True)
    save_csv(bulk, "bulk_control_fpkm.csv", index=True)
    count_rows, co_rows = [], []
    wnt_genes = [g for g in genes if g.startswith("Wnt")]
    for threshold in spec["thresholds"]:
        detected = detection(cells, threshold["cutoff"])
        for gene in genes:
            count_rows.append({"threshold": threshold["label"], "gene": gene,
                               "n_detected": int(detected[gene].sum()), "n_cells": len(cells),
                               "fraction": float(detected[gene].mean()),
                               "median_FPKM_all_cells": float(cells[gene].median()),
                               "max_FPKM": float(cells[gene].max())})
        conditions = {"Wnt5a": detected["Wnt5a"], "Pdgfra": detected["Pdgfra"],
                      "any_Wnt": detected[wnt_genes].any(axis=1)}
        for numerator, denominator in [("Pdgfra", "Wnt5a"), ("Wnt5a", "Pdgfra"),
                                       ("Axin2", "any_Wnt"), ("Porcn", "any_Wnt"),
                                       ("Wls", "any_Wnt")]:
            joint_n, denominator_n, fraction = conditional_fraction(detected[numerator], conditions[denominator])
            co_rows.append({"threshold": threshold["label"], "estimand": f"{numerator}_given_{denominator}",
                            "numerator_gene": numerator, "conditioning_expression": denominator,
                            "joint_n": joint_n, "denominator_n": denominator_n, "fraction": fraction,
                            "n_cells": len(cells)})
    counts, coexpression = pd.DataFrame(count_rows), pd.DataFrame(co_rows)
    save_csv(counts, "detection_by_threshold.csv")
    save_csv(coexpression, "conditional_codetection.csv")
    qc = cells[HOUSEKEEPING].copy()
    qc["n_housekeeping_ge1"] = (qc >= 1).sum(axis=1)
    qc["passes_partial_housekeeping_rule"] = qc["n_housekeeping_ge1"] >= 3
    save_csv(qc, "housekeeping_qc_diagnostic.csv", index=True)
    make_figures(cells, counts, coexpression, spec, metadata)
    summary = {
        "status": "Descriptive source-data reproduction; exact published percentage not assumed",
        "n_single_cells": len(cells), "n_bulk_controls_excluded": len(bulk),
        "partial_housekeeping_qc_pass": int(qc["passes_partial_housekeeping_rule"].sum()),
        "target_gene_status_not_OK": int((long["status"] != "OK").sum()),
        "pdgfra_given_wnt5a": coexpression[coexpression["estimand"] == "Pdgfra_given_Wnt5a"].to_dict("records"),
        "primary_threshold_counts": counts[counts["threshold"] == ">=1"][["gene", "n_detected", "n_cells"]].to_dict("records"),
        "independent_animal_mapping": "Not supplied in deposited sample metadata; no animal-level inference",
    }
    (BASE / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    import matplotlib
    provenance = {
        "downloaded_on_utc_date": "2026-09-22", "accession": "GSE109444",
        "paper_doi": "10.1126/science.aam6603", "paper_pmid": "29420258",
        "series_metadata": {k: v for k, v in series.items() if k.startswith("Series_") and k != "Series_sample_id"},
        "input_files": [{"path": str((RAW / name).relative_to(ROOT)).replace("\\", "/"),
                         "source_url": url, "bytes": (RAW / name).stat().st_size,
                         "sha256": sha256(RAW / name)} for name, url in SOURCES.items()],
        "archive_members": member_records,
        "code_sha256": sha256(__file__), "specification_sha256": sha256(BASE / "specification.json"),
        "runtime": {"python": platform.python_version(), "pandas": pd.__version__,
                    "numpy": np.__version__, "matplotlib": matplotlib.__version__},
        "raw_processing": "Deposited gene-level Cufflinks 2.0.2 FPKM (mm10); no raw-read remapping or re-normalization",
        "gene_resolution": "Exact gene_short_name; one record required per target per sample; no duplicate summation",
        "output_sha256": {str(path.relative_to(BASE)).replace("\\", "/"): sha256(path)
                          for folder in ["tables", "figures"] for path in sorted((BASE / folder).iterdir()) if path.is_file()},
    }
    (BASE / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ["n_single_cells", "n_bulk_controls_excluded", "partial_housekeeping_qc_pass", "pdgfra_given_wnt5a"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    main()
