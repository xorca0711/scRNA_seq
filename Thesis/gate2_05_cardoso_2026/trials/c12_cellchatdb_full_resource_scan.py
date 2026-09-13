#!/usr/bin/env python
"""Trial C12: does anything in the whole CellChatDB outrank the axis the paper followed.

Trial C3 scored four EGFR ligands because those are the ones the paper named,
and found Areg top of the four in every mutant library. That is a statement
about a shortlist somebody else drew. This trial scores the authors' entire
curated ligand-receptor resource instead, so that a pair nobody in this project
thought to look at has a chance to appear.

**What this is not.** It is not a CellChat rerun of the paper's analysis. It
uses CellChat's own resource and CellChat's own scoring logic through LIANA's
implementation of them, on different data, and it reads rankings rather than
significance. The permutation test is switched off deliberately: its unit is
the cell, and a cell-level P value presented beside a group question is the
exact number this repository refuses to emit. Nothing here says two cells
communicate, because co-expression carries no proximity.

**Why not the mouse deposit, stated rather than worked around.** Three reasons,
each fatal on its own. The paper's own configuration integrates epithelium from
GSE247505 with the niche from GSE316241 and GSE316243, and trial C0 found those
sit in different gene spaces needing an Ensembl intersection. Ligand and
receptor would then come from different libraries at different sequencing
depths, which confounds every product-based score this method computes. And
each mouse genotype contributes one library of three pooled mice, so no
statistic has a unit. Running it on the mesenchymal library alone would be
worse still: trial C1d found 184 Areg-high epithelial contaminants in that
sort, so the epithelial side of the scan would be the artefact this repository
spent three trials characterising.

GSE136831 avoids all four problems. Epithelium and fibroblasts come from the
same library and the same donor at the same depth, and there are enough donors
to make the donor the unit.

Frozen rules, set before any score was computed:

* Dataset: GSE136831 (Adams et al. 2020, doi:10.1126/sciadv.aba1983), deposited
  raw counts and metadata. Disease_Identity in IPF or Control; COPD out of
  scope, as in trials E2 and E6.
* Cells: CellType_Category Epithelial or Stromal, with the fine
  Manuscript_Identity as the grouping LIANA receives.
* Source populations: every epithelial type. Target populations: Fibroblast and
  Myofibroblast. Pericytes and smooth muscle are excluded as mural, as in E6.
* Normalisation: counts to 10,000 per cell then log1p, which is what the method
  expects.
* Method: LIANA's implementation of CellChat, resource `cellchatdb`, with
  `n_perms=None` so no permutation P value is produced at all.
* Unit: the donor. The method is run separately per donor, and a donor
  contributes only if it holds at least 50 cells in both compartments, the same
  floor as trials E1 to E3 and E6. Ranks are aggregated across donors by the
  median, with the interquartile range reported.

* T1, THE PRIMARY READING. Where AREG to EGFR sits in the ranking of all
  epithelium-to-fibroblast pairs, by median rank across donors, and which pairs
  outrank it.
* T2. The fifteen highest-ranked pairs, with the number of donors each appears
  in.
* T3, THE GUARD THAT DECIDES WHETHER T1 AND T2 CAN BE READ AT ALL. Ligand-
  receptor rankings are known to be dominated by collagen and laminin ligands
  against integrin receptors, because those genes are abundant in every
  fibroblast rather than informative about signalling. The fraction of the top
  fifteen whose ligand starts with COL or LAM or whose receptor contains ITG or
  SDC is computed. RULE: if more than half the top fifteen are such pairs, the
  ranking is reporting transcript abundance, and neither T1 nor T2 is read.
* T4. AREG to EGFR is also located within the subset of pairs that share its
  receptor, so that its rank is reported once against everything and once
  against the other EGFR ligands, which is the comparison trial C3 made.

* THE READING, fixed in advance:
  - T3 fails: the ranking reports abundance; nothing else is read;
  - AREG to EGFR in the top decile and T3 passes: the axis the paper followed
    is prominent when the whole resource is scored in human fibrosis;
  - AREG to EGFR outside the top decile and T3 passes: the axis is not
    prominent against the full resource, which qualifies the paper's choice of
    shortlist without contradicting any of its measurements;
  - either way, the pairs that outrank it are reported as leads and nothing
    more, because this method cannot distinguish a real interaction from
    co-expression.
* NOT TESTED: whether any pair represents an actual interaction, which needs
  proximity; and anything about the mouse model, for the four reasons above.
"""

from __future__ import annotations

import gc
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RAW, RunRecord, df_to_markdown  # noqa: E402
from mtx_stream import open_text, read_lines  # noqa: E402

OUT = HERE / "c12_cellchatdb_full_resource_scan"
OUT.mkdir(exist_ok=True)
SERIES = RAW / "GSE136831"
MATRIX = SERIES / "GSE136831_RawCounts_Sparse.mtx.gz"
GENE_IDS = SERIES / "GSE136831_AllCells.GeneIDs.txt.gz"
BARCODES = SERIES / "GSE136831_AllCells.cellBarcodes.txt.gz"
METADATA = SERIES / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
CACHE = SERIES / "c12_epithelial_stromal.h5ad"   # regenerable, outside the tracked tree

TARGETS = ["Fibroblast", "Myofibroblast"]
DISEASES = ["IPF", "Control"]
MIN_CELLS = 50
TOP_N = 15
ABUNDANCE_SHARE = 0.5
AXIS = ("AREG", "EGFR")

RULES = {
    "question": "does any pair in the whole CellChatDB outrank the axis the paper followed",
    "what_this_is_not": "not a CellChat rerun of the paper's analysis; CellChat's resource and scoring "
                        "logic through LIANA, on different data, read as rankings not significance",
    "permutations": "switched off (n_perms=None); a cell-level P value beside a group question is the "
                    "number this repository refuses to emit",
    "why_not_the_mouse_deposit": [
        "the paper's configuration integrates GSE247505 epithelium with the niche, and C0 found those "
        "sit in different gene spaces",
        "ligand and receptor would come from different libraries at different depths, confounding every "
        "product-based score",
        "one library of three pooled mice per genotype, so no statistic has a unit",
        "the mesenchymal library alone carries 184 Areg-high epithelial contaminants (C1d), which would "
        "make the epithelial side of the scan the artefact this repository characterised",
    ],
    "dataset": "GSE136831, IPF and Control; COPD out of scope as in E2 and E6",
    "cells": "CellType_Category Epithelial or Stromal, grouped by Manuscript_Identity",
    "sources": "every epithelial type", "targets": TARGETS,
    "excluded_targets": "Pericyte and SMC, as mural, following E6",
    "normalisation": "counts to 10,000 per cell then log1p",
    "method": "LIANA implementation of CellChat, resource cellchatdb",
    "unit": "the donor; run per donor, floor of " + str(MIN_CELLS) + " cells in both compartments, "
            "ranks aggregated by median with the interquartile range reported",
    "T1": "where AREG to EGFR sits among all epithelium-to-fibroblast pairs, and what outranks it",
    "T2": "the top " + str(TOP_N) + " pairs with the number of donors each appears in",
    "T3": "GUARD: if more than half the top " + str(TOP_N) + " are collagen or laminin ligands against "
          "integrin or syndecan receptors, the ranking reports abundance and nothing is read",
    "T4": "the same rank reported against the other ligands of the same receptor, the comparison C3 made",
    "not_tested": ["whether any pair is an actual interaction, which needs proximity",
                   "anything about the mouse model, for the four reasons above"],
}


def build_object(rec):
    """One streaming pass, keeping every gene row for the cells in scope only."""
    import anndata as ad

    if CACHE.exists():
        return ad.read_h5ad(CACHE)
    for path in (MATRIX, GENE_IDS, BARCODES, METADATA):
        rec.add_input(path)
    meta = pd.read_csv(METADATA, sep="\t")
    barcodes = read_lines(BARCODES)
    meta = meta.set_index("CellBarcode_Identity").reindex(barcodes)
    keep = (meta["Disease_Identity"].isin(DISEASES).to_numpy()
            & meta["CellType_Category"].isin(["Epithelial", "Stromal"]).to_numpy())
    wanted = np.where(keep)[0]
    position = -np.ones(len(barcodes), dtype=np.int64)
    position[wanted] = np.arange(len(wanted))
    print("cells in scope:", len(wanted), "of", len(barcodes))

    genes = read_lines(GENE_IDS, column=1)[1:]
    rows_chunks, cols_chunks, data_chunks = [], [], []
    rows, cols, data = [], [], []
    seen = 0
    with open_text(MATRIX) as handle:
        shape = None
        for line in handle:
            if line.startswith("%"):
                continue
            shape = [int(v) for v in line.split()]
            break
        if shape[0] != len(genes):
            raise SystemExit(f"{shape[0]} matrix rows against {len(genes)} gene names")
        for line in handle:
            g, c, v = line.split()
            target = position[int(c) - 1]
            if target < 0:
                continue
            rows.append(target)
            cols.append(int(g) - 1)
            data.append(int(v))
            seen += 1
            if len(rows) >= 8_000_000:
                rows_chunks.append(np.asarray(rows, dtype=np.int32))
                cols_chunks.append(np.asarray(cols, dtype=np.int32))
                data_chunks.append(np.asarray(data, dtype=np.float32))
                rows, cols, data = [], [], []
    if rows:
        rows_chunks.append(np.asarray(rows, dtype=np.int32))
        cols_chunks.append(np.asarray(cols, dtype=np.int32))
        data_chunks.append(np.asarray(data, dtype=np.float32))
    del rows, cols, data
    gc.collect()
    X = sp.coo_matrix((np.concatenate(data_chunks),
                       (np.concatenate(rows_chunks), np.concatenate(cols_chunks))),
                      shape=(len(wanted), len(genes))).tocsr()
    del rows_chunks, cols_chunks, data_chunks
    gc.collect()
    print("kept", X.shape, "with", X.nnz, "non-zero entries")
    obs = meta.iloc[wanted][["Subject_Identity", "Disease_Identity", "CellType_Category",
                             "Manuscript_Identity", "nGene"]].copy()
    obs.columns = ["donor", "disease", "category", "celltype", "n_genes"]
    obs.index = pd.Index([barcodes[i] for i in wanted])
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes)))
    adata.var_names_make_unique()
    adata.write_h5ad(CACHE)
    return adata


def main():
    import scanpy as sc
    from liana.method import cellchat

    rec = RunRecord(OUT / "c12_run_record.json", "C12 CellChatDB full-resource scan", RULES)
    adata = build_object(rec)
    rec.set("cells_in_scope", int(adata.n_obs))
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    counts = adata.obs.groupby(["donor", "category"], observed=True).size().unstack(fill_value=0)
    usable = counts[(counts.get("Epithelial", 0) >= MIN_CELLS)
                    & (counts.get("Stromal", 0) >= MIN_CELLS)].index.tolist()
    rec.set("donors_usable", len(usable))
    print("donors clearing the floor in both compartments:", len(usable))
    if len(usable) < 5:
        raise SystemExit("fewer than five usable donors; refusing to aggregate")

    frames = []
    for donor in usable:
        part = adata[adata.obs["donor"] == donor].copy()
        present = part.obs["celltype"].value_counts()
        part = part[part.obs["celltype"].isin(present[present >= 10].index)].copy()
        if not set(TARGETS) & set(part.obs["celltype"]):
            continue
        try:
            cellchat(part, groupby="celltype", resource_name="cellchatdb",
                     expr_prop=0.1, use_raw=False, n_perms=None, verbose=False)
        except Exception as error:                     # noqa: BLE001
            print("donor", donor, "skipped:", type(error).__name__, str(error)[:90])
            del part
            gc.collect()
            continue
        res = part.uns["liana_res"]
        res = res[res["target"].isin(TARGETS)
                  & ~res["source"].isin(TARGETS)].copy()
        if res.empty:
            del part
            gc.collect()
            continue
        res["pair"] = res["ligand_complex"] + " to " + res["receptor_complex"]
        best = res.sort_values("lr_probs", ascending=False).drop_duplicates("pair")
        best["rank"] = np.arange(1, len(best) + 1)
        best["donor"] = donor
        frames.append(best[["donor", "pair", "ligand_complex", "receptor_complex",
                            "source", "target", "lr_probs", "rank"]])
        del part, res
        gc.collect()
        print("donor", donor, "ranked", len(best), "pairs")

    if not frames:
        raise SystemExit("no donor produced a ranking; refusing to report")
    per_donor = pd.concat(frames, ignore_index=True)
    per_donor.to_csv(OUT / "c12_per_donor_pairs.csv", index=False)
    rec.add_output(OUT / "c12_per_donor_pairs.csv")

    agg = (per_donor.groupby("pair")
           .agg(n_donors=("donor", "nunique"), median_rank=("rank", "median"),
                q25=("rank", lambda s: s.quantile(0.25)),
                q75=("rank", lambda s: s.quantile(0.75)),
                median_lr_prob=("lr_probs", "median"),
                ligand=("ligand_complex", "first"), receptor=("receptor_complex", "first"))
           .reset_index())
    agg = agg[agg["n_donors"] >= max(3, len(usable) // 2)]
    agg = agg.sort_values("median_rank").reset_index(drop=True)
    agg["overall_rank"] = np.arange(1, len(agg) + 1)
    agg.to_csv(OUT / "c12_aggregated_ranking.csv", index=False)
    rec.add_output(OUT / "c12_aggregated_ranking.csv")

    top = agg.head(TOP_N)
    def is_abundance(row):
        ligand = str(row["ligand"]).upper()
        receptor = str(row["receptor"]).upper()
        return (ligand.startswith(("COL", "LAM")) or "ITG" in receptor or "SDC" in receptor)
    abundance = int(top.apply(is_abundance, axis=1).sum())
    t3_fails = bool(abundance > ABUNDANCE_SHARE * len(top))
    rec.set("T3", {"top_n": int(len(top)), "abundance_pairs": abundance,
                   "share": round(abundance / max(len(top), 1), 3), "ranking_readable": not t3_fails})

    axis_pair = AXIS[0] + " to " + AXIS[1]
    row = agg[agg["pair"] == axis_pair]
    if row.empty:
        axis_rank, axis_note = None, "the pair does not clear the donor floor in this cohort"
    else:
        axis_rank = int(row.iloc[0]["overall_rank"])
        axis_note = ("median rank " + str(row.iloc[0]["median_rank"]) + " across "
                     + str(int(row.iloc[0]["n_donors"])) + " donors")
    decile = max(1, len(agg) // 10)
    rec.set("T1", {"pair": axis_pair, "overall_rank": axis_rank, "of_pairs": int(len(agg)),
                   "top_decile_boundary": decile, "note": axis_note})

    same_receptor = agg[agg["receptor"] == AXIS[1]].copy()
    same_receptor["rank_within_receptor"] = np.arange(1, len(same_receptor) + 1)
    same_receptor.to_csv(OUT / "c12_same_receptor.csv", index=False)
    rec.add_output(OUT / "c12_same_receptor.csv")
    within = same_receptor[same_receptor["pair"] == axis_pair]
    rec.set("T4", {"rank_among_ligands_of_" + AXIS[1]:
                   int(within.iloc[0]["rank_within_receptor"]) if len(within) else None,
                   "n_ligands_of_that_receptor": int(len(same_receptor))})

    if t3_fails:
        reading = ("the ranking reports transcript abundance, because " + str(abundance) + " of the top "
                   + str(len(top)) + " pairs are collagen or laminin against integrin or syndecan; "
                   "neither T1 nor T2 is read")
    elif axis_rank is not None and axis_rank <= decile:
        reading = ("the axis the paper followed is prominent when the whole resource is scored in human "
                   "fibrosis, at rank " + str(axis_rank) + " of " + str(len(agg)))
    elif axis_rank is not None:
        reading = ("the axis is not prominent against the full resource, at rank " + str(axis_rank)
                   + " of " + str(len(agg)) + ", which qualifies the choice of shortlist without "
                   "contradicting any measurement the paper made")
    else:
        reading = "the axis pair does not clear the donor floor in this cohort, so T1 has no answer"
    rec.set("reading", reading)

    lines = ["# Trial C12: the whole CellChatDB, ranked", "",
             "Donors: " + str(len(usable)) + " clearing the floor in both compartments. Pairs kept: "
             + str(len(agg)) + ".",
             "T3 guard: " + str(abundance) + " of the top " + str(len(top))
             + " pairs are abundance-type; ranking readable: " + str(not t3_fails) + ".",
             "T1: " + axis_pair + " " + axis_note + ("" if axis_rank is None else
             (", overall rank " + str(axis_rank) + " of " + str(len(agg))
              + " (top decile is rank " + str(decile) + " or better)")), "",
             "Reading: " + reading, "",
             "**Not tested:** whether any pair is an actual interaction, and anything about the mouse "
             "model. This is co-expression with no proximity, and no permutation P value was computed.",
             "", "## Top " + str(TOP_N) + " epithelium-to-fibroblast pairs", "",
             df_to_markdown(top[["overall_rank", "pair", "n_donors", "median_rank", "q25", "q75",
                                 "median_lr_prob"]], index=False), "",
             "## Ligands of " + AXIS[1] + ", ranked", "",
             df_to_markdown(same_receptor[["rank_within_receptor", "pair", "n_donors", "median_rank",
                                           "median_lr_prob"]].head(12), index=False), ""]
    (OUT / "c12_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c12_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
