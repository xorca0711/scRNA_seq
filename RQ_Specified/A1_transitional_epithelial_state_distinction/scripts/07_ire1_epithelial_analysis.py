"""Run the frozen GSE190821 epithelial RiboTag contrast, with source crosswalk."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_mouse(value):
    return str(int(value)) if str(value).isdigit() else str(value)


def select_samples(metadata, columns, contract):
    """Match explicit compartment, mouse ID and treatment, never column order."""
    rows = []
    for record in metadata["samples"]:
        f = record.get("fields", record)
        char = dict(v.split(": ", 1) for v in f["characteristics_ch1"])
        sel = contract["selection"]
        if any(char[k] != sel[k] for k in ("tissue compartment", "exposure")) or char["treatment"] not in sel["treatment"]:
            continue
        mouse = canonical_mouse(char["mouse identifier"])
        matches = [col for col in columns if "E" in col.split("_")
                   and mouse in [canonical_mouse(t) for t in col.split("_")]]
        if len(matches) != 1:
            raise ValueError(f"Ambiguous/missing column for mouse {mouse}: {matches}")
        col = matches[0]
        if ("KIRA8" in col) != (char["treatment"] == "KIRA8") or "3G9" in col:
            raise ValueError("Column treatment disagrees with GSM metadata")
        rows.append(dict(sample_id=col, gsm=f["geo_accession"][0], source_title=f["title"][0],
                         mouse=mouse, sex=char["mouse sex"], batch=char["batch"], group=char["treatment"],
                         assay=f["source_name_ch1"][0], exposure=char["exposure"], compartment=char["tissue compartment"]))
    frame = pd.DataFrame(rows)
    expected = contract["expected_mice_per_group"]
    if frame.groupby("group").size().to_dict() != {"KIRA8": expected, "Vehicle": expected}:
        raise ValueError("Frozen five-versus-five selection did not match")
    if frame.mouse.duplicated().any() or frame.sample_id.duplicated().any() or frame.gsm.duplicated().any():
        raise ValueError("Repeated animal or library")
    if sorted(frame.batch.unique()) != sorted(contract["expected_batches"]):
        raise ValueError("Unexpected batches")
    if not frame.assay.str.contains("Ribotag immunoprecipitation", case=False).all():
        raise ValueError("Wrong assay")
    return frame


def draw_figures(out, contract, marker_map):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none",
                         "axes.spines.top": False, "axes.spines.right": False})
    color = {"Vehicle": "#626B76", "KIRA8": "#167C8F"}
    shape = {"S061": "o", "S135": "^"}
    pc = pd.read_csv(out / "sample_PCA.tsv", sep="\t")
    variance = pd.read_csv(out / "PCA_variance.tsv", sep="\t")
    effects = pd.read_csv(out / "gene_effects.tsv", sep="\t").set_index("feature_id")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    ax = axes[0]
    for (group, batch), block in pc.groupby(["group", "batch"]):
        ax.scatter(block.PC1, block.PC2, c=color[group], marker=shape[batch], label=f"{group} · {batch}", s=52)
    ax.set(xlabel=f"PC1 ({variance.variance_fraction.iloc[0]:.1%})", ylabel=f"PC2 ({variance.variance_fraction.iloc[1]:.1%})", title="A  Sample structure before batch adjustment")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax = axes[1]
    hit = effects.FDR < .05
    ax.scatter(effects.loc[~hit, "logFC"], -np.log10(effects.loc[~hit, "FDR"].clip(lower=1e-300)), c="#B9C0C5", s=5, rasterized=True)
    ax.scatter(effects.loc[hit, "logFC"], -np.log10(effects.loc[hit, "FDR"].clip(lower=1e-300)), c="#167C8F", s=7, label=f"BH FDR < 0.05: {hit.sum():,}", rasterized=True)
    ax.axhline(-np.log10(.05), color="#888888", lw=.7, ls="--")
    ax.axvline(0, color="#888888", lw=.7)
    ax.set(xlabel="log2 fold change (KIRA8 − vehicle)", ylabel="−log10 BH FDR", title="B  Treatment effect, adjusted for batch and sex")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("IRE1α inhibition: epithelial RiboTag RNA after bleomycin", fontsize=13)
    fig.text(.5, .02, "GSE190821 · 5 mice per arm · day 7 RNA · bulk epithelial mixture; no single-cell fate measurement", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .06, 1, .94))
    for ext in ("png", "svg"):
        fig.savefig(BASE / f"figures/a1_ire1_rna_overview.{ext}", dpi=220)
    plt.close(fig)
    values = pd.read_csv(out / "marker_logCPM.tsv", sep="\t")
    meta = pd.read_csv(out / "sample_manifest.tsv", sep="\t")
    fig, axes = plt.subplots(2, 4, figsize=(11, 6.7))
    marker_results = []
    for ax, symbol in zip(axes.flat, contract["markers"]):
        fid = marker_map["records"].get(symbol, {}).get("id")
        val = values[values.feature_id == fid]
        effect = effects.loc[fid] if fid in effects.index else None
        if val.empty:
            ax.text(.5, .5, "Not retained", ha="center", transform=ax.transAxes)
        else:
            for i, group in enumerate(["Vehicle", "KIRA8"]):
                block = meta[meta.group == group].sort_values("mouse")
                offsets = np.linspace(-.15, .15, len(block))
                yy = []
                for dx, row in zip(offsets, block.itertuples()):
                    y = float(val[row.sample_id].iloc[0]); yy.append(y)
                    ax.scatter(i + dx, y, c=color[group], marker=shape[row.batch], s=40, zorder=3)
                ax.plot([i - .24, i + .24], [np.mean(yy)] * 2, c=color[group], lw=1.3)
            ax.grid(axis="y", alpha=.15)
        subtitle = f"log2FC {effect.logFC:+.2f}; q={effect.FDR:.3g}" if effect is not None else "Not tested"
        ax.set(xticks=[0, 1], xticklabels=["Vehicle", "KIRA8"], title=f"{symbol}\n{subtitle}", ylabel="TMM log2CPM")
        marker_results.append(dict(symbol=symbol, feature_id=fid, retained=effect is not None,
                                   logFC=float(effect.logFC) if effect is not None else None,
                                   FDR=float(effect.FDR) if effect is not None else None))
    fig.suptitle("Predefined epithelial and injury-response markers", fontsize=14)
    fig.text(.5, .025, "Each point is one mouse; lines are group means. Circles: S061; triangles: S135.\nExpression display is unadjusted; effects and q-values use batch + sex + treatment across all tested genes.", ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .95))
    for ext in ("png", "svg"):
        fig.savefig(BASE / f"figures/a1_ire1_marker_panel.{ext}", dpi=220)
    plt.close(fig)
    pd.DataFrame(marker_results).to_csv(out / "predefined_marker_effects.tsv", sep="\t", index=False)


def main():
    if (BASE / "reports/ire1_run.json").exists():
        raise SystemExit("Existing IRE1 run: archive its outputs before a new numerical run; use script 10 for figure-only rendering.")
    ap = argparse.ArgumentParser()
    ap.add_argument("--rscript", type=Path, required=True)
    args = ap.parse_args()
    contract_path = BASE / "config/ire1_kira8.json"
    c = json.loads(contract_path.read_text(encoding="utf-8"))
    count_path = BASE / "cache/inputs/GSE190821/GSE190821_counts.csv.gz"
    meta_path = BASE / "metadata/GSE190821.json"
    assert c["status"] == "frozen_before_expression_results"
    assert sha(count_path) == c["count_source_sha256"] and sha(meta_path) == c["metadata_sha256"]
    assert sha(BASE / "metadata/IRE1_published_gene_signatures.csv") == c["signature_sha256"]
    x = pd.read_csv(count_path, index_col=0)
    meta = select_samples(json.loads(meta_path.read_text(encoding="utf-8")), x.columns, c)
    x = x.loc[:, meta.sample_id]
    a = x.to_numpy()
    assert not x.index.duplicated().any() and np.isfinite(a).all() and (a >= 0).all() and (a == np.floor(a)).all()
    assert (x.sum() >= c["qc"]["minimum_library_count"]).all()
    out = BASE / "tables/ire1"; out.mkdir(parents=True, exist_ok=True)
    cache = BASE / "processed/ire1"; cache.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out / "sample_manifest.tsv", sep="\t", index=False)
    x.to_csv(cache / "counts.tsv", sep="\t", index_label="feature_id")
    qc = meta.copy()
    qc["raw_library_count"] = x.sum().to_numpy()
    qc["genes_detected"] = (x > 0).sum().to_numpy()
    qc.to_csv(out / "count_QC.tsv", sep="\t", index=False)
    marker_map = json.loads((BASE / "metadata/ensembl_marker_lookup.json").read_text())
    pd.DataFrame([{"symbol": s, "feature_id": marker_map["records"].get(s, {}).get("id", "unmapped")} for s in c["markers"]]).to_csv(out / "marker_mapping.tsv", sep="\t", index=False)
    signatures = pd.read_csv(BASE / "metadata/IRE1_published_gene_signatures.csv")
    set_map = json.loads((BASE / "metadata/ensembl_signature_lookup.json").read_text())
    rows = []
    for name in c["gene_sets"]["selected"]:
        for symbol in sorted(set(signatures[name].dropna())):
            rows.append({"gene_set": name, "symbol": symbol, "feature_id": (set_map["records"].get(symbol) or {}).get("id", "unmapped")})
    pd.DataFrame(rows).to_csv(out / "gene_set_mapping.tsv", sep="\t", index=False)
    parameters = dict(counts_path=(cache / "counts.tsv").as_posix(), samples_path=(out / "sample_manifest.tsv").as_posix(),
                      counts_md5=hashlib.md5((cache / "counts.tsv").read_bytes()).hexdigest(),
                      samples_md5=hashlib.md5((out / "sample_manifest.tsv").read_bytes()).hexdigest(),
                      primary_model=c["model"], sensitivity_model=c["sensitivity_model"],
                      min_count=c["filter"]["min_count"], min_total_count=c["filter"]["min_total_count"],
                      minimum_library_count=c["qc"]["minimum_library_count"],
                      gene_set_min_fraction=c["gene_sets"]["minimum_fraction_of_source_symbols_in_tested_universe"],
                      gene_set_min_genes=c["gene_sets"]["minimum_genes"], pca_features=c["pca"]["top_variable_features"])
    pd.DataFrame(parameters.items(), columns=["key", "value"]).to_csv(out / "input_contract.tsv", sep="\t", index=False)
    rpath = Path(__file__).with_name("08_fit_ire1.R")
    record = {"status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
              "contract_sha256": sha(contract_path), "source_sha256": sha(count_path), "metadata_sha256": sha(meta_path),
              "python_sha256": sha(Path(__file__)), "R_sha256": sha(rpath), "independent_mice": 10,
              "annotation_sha256": {name: sha(BASE / "metadata" / name) for name in ("ensembl_marker_lookup.json", "ensembl_signature_lookup.json")},
              "claim_ceiling": c["inference"]}
    runpath = BASE / "reports/ire1_run.json"
    runpath.write_text(json.dumps(record, indent=2) + "\n")
    try:
        run = subprocess.run([str(args.rscript.resolve()), str(rpath), str(out)], cwd=ROOT, text=True, capture_output=True)
        (BASE / "reports/ire1_R.log").write_text(run.stdout + run.stderr)
        if run.returncode:
            raise RuntimeError(f"R failed ({run.returncode}): {run.stderr[-2000:]}")
        draw_figures(out, c, marker_map)
        record["status"] = "completed"
        record["fit"] = pd.read_csv(out / "fit_summary.tsv", sep="\t").to_dict("records")[0]
        record["output_sha256"] = {p.relative_to(BASE).as_posix(): sha(p) for p in sorted(out.glob("*.tsv"))}
    except Exception as e:
        record.update(status="failed", error=str(e)); raise
    finally:
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        runpath.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["fit"], indent=2))


if __name__ == "__main__":
    main()
