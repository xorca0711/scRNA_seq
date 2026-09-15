#!/usr/bin/env python
"""Trial D1: call cells from the raw whitelists, measure quality, flag doublets.

Six of the eight Choi-2020 matrices are unfiltered 10x barcode whitelists
(trial D0), so cell calling is this repository's job. Cell Ranger 2.0.2's
caller is not reproduced. Instead the paper's own filter, which it applied
after Cell Ranger to remove empty droplets and doublets, is applied directly
to every barcode: more than 500 and fewer than 7,000 detected genes, more
than 2,000 UMI. Any barcode Cell Ranger would have called and the paper kept
passes this filter; barcodes the caller would have dropped but that carry
more than 2,000 UMI and 500 genes are the only difference, and they are
counted.

Frozen rules, set before any file was opened:

* Calling: the paper's filter, exclusive bounds, on every barcode. A second,
  plain floor (500 counts, 200 genes) is applied for sensitivity only and its
  extra barcodes are counted, never used downstream.
* Quality: the repository's per-library MAD rule on called cells. Lower and
  upper bounds at 5 MAD on log1p(total counts) and log1p(genes); an upper
  bound on mitochondrial fraction at median + 3 MAD, capped at 20 percent.
  Every bound is recorded with the number of cells it removes; a bound that
  removes nothing is reported as not binding, not as an applied filter.
* Doublets: scanpy's Scrublet per capture, seed 0, expected rate 0.8 percent
  per 1,000 called cells (the 10x table), automatic threshold; if the
  automatic threshold is undefined the top expected-rate fraction of scores
  is flagged and the fallback is logged. Flagged cells are kept in the saved
  object with their score and removed by later trials, so the flag can be
  audited. The paper's own Scrublet settings (score above 0.7, cluster mean
  above 0.6) are recorded for the divergence table and not substituted.
* Organoid libraries are already filtered; the same paper filter, MAD rule
  and Scrublet run on them unchanged so every library is treated alike.
* Nothing is normalised, clustered or scored here.
* The library is the unit. Each in vivo library pools two mice (STAR Methods),
  so no count in this trial is a per-animal number.

Outputs: one .h5ad per library under raw_data/<accession>/choi_trials/
(gitignored, regenerable) with a counts layer and the QC and doublet columns;
tracked tables and a run record under trials/d1_cells_and_qc/.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from choi_utils import (ACCESSIONS, DERIVED, DERIVED_ORGANOID, FLOOR_FILTER,  # noqa: E402
                        PAPER_FILTER, PAPER_SCRUBLET, RunRecord, df_to_markdown,
                        library_paths, mad_bounds, qc_metrics, read_mtx_triplet)

OUT = HERE / "d1_cells_and_qc"
OUT.mkdir(exist_ok=True)
RULES = {
    "question": "How many cells does each library carry under the paper's own filter, what is their quality, and which are doublets",
    "calling": PAPER_FILTER,
    "calling_note": "the paper's post-Cell-Ranger filter applied to every raw barcode; Cell Ranger 2.0.2's caller is not reproduced",
    "sensitivity_floor": FLOOR_FILTER,
    "qc": {"nmads_counts_genes": 5, "nmads_mt": 3, "mt_cap_percent": 20,
           "note": "per library; a bound that removes nothing is reported as not binding"},
    "doublets": {"method": "scanpy scrublet per capture, seed 0",
                 "expected_rate": "0.008 per 1,000 called cells (10x table)",
                 "threshold": "automatic; fallback to the top expected-rate fraction if undefined, logged",
                 "paper_settings_recorded_not_used": PAPER_SCRUBLET},
    "unit": "the library; each in vivo library pools two mice (STAR Methods)",
    "paper_reported_counts": {
        "tomato_positive_captured_before_contaminant_removal": 12514,
        "tomato_positive_in_figure_1B": 12086,
        "contaminants_removed_by_paper": {"ciliated": 214, "mesenchyme": 16, "immune": 25},
        "non_lineage_labelled_captured": 14017,
        "non_lineage_labelled_doublet_cluster_removed": 1125,
        "organoid_epithelial_control": 1286,
        "organoid_epithelial_il1b": 2584,
    },
}


def call_and_qc(X: sp.csr_matrix, var: pd.DataFrame, barcodes: np.ndarray, meta: dict,
                rec: RunRecord) -> tuple[ad.AnnData, dict]:
    counts = np.asarray(X.sum(axis=1)).ravel()
    genes = X.getnnz(axis=1)
    paper = (genes > PAPER_FILTER["min_genes_exclusive"]) & (genes < PAPER_FILTER["max_genes_exclusive"]) \
        & (counts > PAPER_FILTER["min_counts_exclusive"])
    floor = (counts >= FLOOR_FILTER["min_counts"]) & (genes >= FLOOR_FILTER["min_genes"])
    row = {**meta, "n_barcodes": int(X.shape[0]),
           "paper_rule_cells": int(paper.sum()), "floor_rule_cells": int(floor.sum()),
           "floor_only_barcodes": int((floor & ~paper).sum()), "paper_only_barcodes": int((paper & ~floor).sum())}

    Xc = X[paper]
    obs = qc_metrics(Xc, var, "mouse")
    obs.index = pd.Index(barcodes[paper], name="barcode")
    keep = np.ones(len(obs), dtype=bool)
    bounds = {}
    for col in ("total_counts", "n_genes"):
        lo, hi = mad_bounds(np.log1p(obs[col].to_numpy()), RULES["qc"]["nmads_counts_genes"])
        m = (np.log1p(obs[col]) >= lo) & (np.log1p(obs[col]) <= hi)
        bounds[col] = {"low": float(np.expm1(lo)), "high": float(np.expm1(hi)),
                       "removed_low": int((np.log1p(obs[col]) < lo).sum()),
                       "removed_high": int((np.log1p(obs[col]) > hi).sum())}
        keep &= m.to_numpy()
    mt = obs["pct_mt"].to_numpy()
    _, hi_mt = mad_bounds(mt, RULES["qc"]["nmads_mt"])
    hi_mt = min(hi_mt, RULES["qc"]["mt_cap_percent"])
    bounds["pct_mt"] = {"high": float(hi_mt), "removed_high": int((mt > hi_mt).sum())}
    keep &= mt <= hi_mt
    row["mad_removed"] = int((~keep).sum())
    row["binding_bounds"] = ";".join(k + ":" + ",".join(s for s, v in b.items() if s.startswith("removed") and v > 0)
                                     for k, b in bounds.items() if any(s.startswith("removed") and v > 0 for s, v in b.items())) or "none"
    Xq = Xc[keep]
    obs = obs[keep].copy()
    adata = ad.AnnData(X=Xq.astype(np.float32), obs=obs,
                       var=var.set_index(var["gene_symbol"].where(~var["gene_symbol"].duplicated(), var["gene_id"])))
    adata.var_names_make_unique()
    adata.var.index.name = "feature"   # the index is the symbol, or the ID where a symbol is duplicated
    adata.layers["counts"] = adata.X.copy()
    for k, v in meta.items():
        adata.obs[k] = v

    n = adata.n_obs
    expected = 0.008 * n / 1000
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=expected)
    thr = float(adata.uns.get("scrublet", {}).get("threshold", np.nan))
    method = "scrublet automatic threshold"
    if not np.isfinite(thr):
        cut = np.quantile(adata.obs["doublet_score"], 1 - expected)
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"] >= cut
        thr = float(cut)
        method = f"fallback: top {100 * expected:.2f}% of scores"
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    row.update({"cells_after_qc": int(n), "scrublet_expected_rate": float(expected),
                "scrublet_threshold": thr, "scrublet_method": method,
                "doublets_flagged": int(adata.obs["predicted_doublet"].sum()),
                "cells_final": int((~adata.obs["predicted_doublet"]).sum()),
                "median_counts": float(np.median(adata.obs["total_counts"])),
                "median_genes": float(np.median(adata.obs["n_genes"])),
                "median_pct_mt": float(np.median(adata.obs["pct_mt"])),
                "paper_scrublet_score_gt_0.7": int((adata.obs["doublet_score"] > 0.7).sum())})
    rec.set(f"bounds_{meta['library']}", bounds)
    return adata, row


def main() -> None:
    rec = RunRecord(OUT / "d1_run_record.json", "D1 Choi-2020 cell calling, quality and doublets", RULES)
    rows = []
    for accession, info in ACCESSIONS.items():
        target = DERIVED if accession == "GSE145031" else DERIVED_ORGANOID
        target.mkdir(parents=True, exist_ok=True)
        for entry in info["libraries"]:
            paths = library_paths(accession, entry["library"], entry["gsm"])
            for p in paths.values():
                rec.add_input(p)
            print(f"[D1] reading {entry['library']}", flush=True)
            X, var, bc = read_mtx_triplet(paths["matrix"], paths["features"], paths["barcodes"])
            meta = {"accession": accession, **entry}
            adata, row = call_and_qc(X, var, bc, meta, rec)
            del X
            out = target / f"{entry['library']}.h5ad"
            adata.write_h5ad(out)
            row["object"] = str(out.relative_to(DERIVED.parents[1]))
            rows.append(row)
            print(f"[D1] {entry['library']}: {row['paper_rule_cells']} called, {row['cells_final']} final", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "d1_cell_calling.csv", index=False)
    rec.add_output(OUT / "d1_cell_calling.csv")

    tom = table[table.get("sort", pd.Series(index=table.index, dtype=object)) == "Tomato"]
    tomato_total = int(tom["cells_final"].sum())
    tomato_called = int(tom["paper_rule_cells"].sum())
    paper = RULES["paper_reported_counts"]
    summary = [
        "# Trial D1: cells, quality and doublets",
        "",
        f"Tomato-positive libraries: {tomato_called} barcodes pass the paper's filter and {tomato_total} remain after the MAD rule "
        f"and Scrublet, against the paper's {paper['tomato_positive_captured_before_contaminant_removal']} captured and "
        f"{paper['tomato_positive_in_figure_1B']} in Figure 1B (their counts precede this repository's doublet removal and "
        "follow Cell Ranger's caller, so the numbers are comparable in size, not identical by construction).",
        "",
        "## Per library",
        "",
        df_to_markdown(table[["library", "n_barcodes", "paper_rule_cells", "floor_rule_cells", "floor_only_barcodes",
                              "mad_removed", "binding_bounds", "scrublet_method", "scrublet_threshold",
                              "doublets_flagged", "cells_final", "median_counts", "median_genes", "median_pct_mt"]],
                       index=False),
        "",
        "The organoid libraries were deposited already filtered; the paper's filter, the MAD rule and Scrublet were run on "
        "them unchanged. The paper reports 1,286 control and 2,584 IL-1beta epithelial cells after removing an EpCAM-negative "
        "stromal cluster, which is trial D5's job, not this one's.",
    ]
    (OUT / "d1_summary.md").write_text("\n".join(summary) + "\n", encoding="utf-8")
    rec.add_output(OUT / "d1_summary.md")
    rec.set("tomato_positive_called", tomato_called)
    rec.set("tomato_positive_final", tomato_total)
    rec.finish()
    print("[D1] done", flush=True)


if __name__ == "__main__":
    main()
