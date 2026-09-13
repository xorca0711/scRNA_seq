#!/usr/bin/env python
"""Trial E4: is the Areg-independent tier tumour-specific, or does injury alone produce it.

Extension 4 of the plan in DIVERGENCES_AND_NEXT.md, on the dataset the paper
itself used to argue that its tumour fibroblast state is shared with injury
(Tsukui et al. 2020, Nat Commun, doi:10.1038/s41467-020-15647-5; GEO
GSE132771).

The question. Trial C2 found that deleting Areg leaves Runx1 and Pdgfrb almost
untouched while Fst and Runx2 lose most of their detection. Two readings of
that survive so far: a second signal drives the Runx1 and Pdgfrb half, or
those two genes are simply what any activated lung fibroblast expresses, in
which case their persistence says nothing about the tumour niche at all. This
trial separates them by asking what bleomycin injury alone does to the same
six genes, in sorted Col1a1-GFP mesenchyme with no oncogene anywhere.

Frozen rules, set before any matrix was read:

* Libraries: the four mouse Col1a1-GFP-positive samples, Bleo1, Bleo2, UT1 and
  UT2. The GFP-negative samples and the human samples of the same series are
  out of scope here.
* Quality control and doublet removal follow trial C1 exactly.
* Genes: the paper's six reprogrammed-fibroblast markers (Tnc, Fst, Runx1,
  Runx2, Acta2, Pdgfrb), the four EGFR ligands, Egfr, and Pdgfra and Col13a1
  as the alveolar-identity comparison.
* Measure: the fraction of cells with a non-zero count, per library.
* THE TEST, frozen. A gene is "injury-generic" if its detection is higher in
  **both** bleomycin libraries than in **both** untreated libraries. A gene
  where that fails is not injury-generic by this rule.
* The reading, fixed in advance so it cannot be chosen after the fact:
  - if Runx1 and Pdgfrb are injury-generic, then their survival of Areg
    deletion in trial C2 is what an activated fibroblast does, not evidence of
    a second tumour signal, and claim C29 weakens accordingly;
  - if they are not injury-generic while Fst and Runx2 are, the tumour reading
    survives and the second-signal hypothesis stays open.
* Two animals per group, so directions only. No P value, and the frozen rule
  deliberately requires both replicates to agree rather than averaging them.
"""

from __future__ import annotations

import gc
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import (ENSEMBL_ID_RE, RAW, REPO, RunRecord,  # noqa: E402
                           df_to_markdown, read_mtx_triplet)

sys.path.insert(0, str(REPO / "analysis" / "scripts"))
from pipeline_utils import apply_thresholds, derive_thresholds, make_unique  # noqa: E402

OUT = HERE / "e4_bleomycin_fibrotic_genes"
OUT.mkdir(exist_ok=True)
ROOT = RAW / "GSE132771" / "GSE132771_RAW"

FIBROTIC = ["Tnc", "Fst", "Runx1", "Runx2", "Acta2", "Pdgfrb"]
OTHER = ["Areg", "Ereg", "Hbegf", "Tgfa", "Egfr", "Pdgfra", "Col13a1"]
LIBRARIES = [
    {"gsm": "GSM3891612", "library": "Bleo1_GFPp", "group": "bleomycin"},
    {"gsm": "GSM3891613", "library": "Bleo2_GFPp", "group": "bleomycin"},
    {"gsm": "GSM3891616", "library": "UT1_GFPp", "group": "untreated"},
    {"gsm": "GSM3891617", "library": "UT2_GFPp", "group": "untreated"},
]
RULES = {
    "dataset": "GSE132771 (Tsukui et al. 2020), the injury comparison the paper itself used",
    "libraries": LIBRARIES,
    "scope": "mouse Col1a1-GFP-positive samples only; GFP-negative and human samples out of scope",
    "qc": "identical to trial C1",
    "genes": FIBROTIC + OTHER,
    "measure": "fraction of cells with a non-zero count, per library",
    "test": "a gene is injury-generic if detection is higher in both bleomycin libraries than in both untreated libraries",
    "reading_fixed_in_advance": {
        "runx1_pdgfrb_injury_generic": "their survival of Areg deletion in C2 is what an activated fibroblast does; claim C29 weakens",
        "runx1_pdgfrb_not_generic_while_fst_runx2_are": "the tumour reading survives and the second-signal hypothesis stays open",
    },
    "unit": "the animal; two per group; directions only, no P value",
}


def load(entry, rec):
    import anndata as ad
    import scanpy as sc

    stem = entry["gsm"] + "_" + entry["library"] + "_"
    matrix = ROOT / (stem + "matrix.mtx.gz")
    features = ROOT / (stem + "genes.tsv.gz")
    if not features.exists():
        features = ROOT / (stem + "features.tsv.gz")
    barcodes = ROOT / (stem + "barcodes.tsv.gz")
    for path in (matrix, features, barcodes):
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(matrix, features, barcodes)
    if "feature_type" not in var or var["feature_type"].isna().all():
        var["feature_type"] = "Gene Expression"
    keep_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    if keep_gene.sum() == 0:
        keep_gene = np.ones(len(var), dtype=bool)
    X, var = X[:, keep_gene], var.loc[keep_gene].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
    obs = pd.DataFrame(index=pd.Index([entry["library"] + "_" + b for b in bc]))
    adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)
    symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True)
    th = derive_thresholds(adata.obs, entry["library"])
    adata = adata[apply_thresholds(adata.obs, th).to_numpy()].copy()
    rate = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=rate)
    auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
    if auto < 0.2 * rate or auto > 3 * rate:
        cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - rate))
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
    return adata[~adata.obs["predicted_doublet"].to_numpy()].copy()


def main():
    rec = RunRecord(OUT / "e4_run_record.json",
                    "E4 bleomycin injury and the fibrotic gene set", RULES)
    rows = []
    for entry in LIBRARIES:
        adata = load(entry, rec)
        row = {"library": entry["library"], "group": entry["group"], "n_cells": int(adata.n_obs)}
        for gene in FIBROTIC + OTHER:
            if gene not in adata.var_names:
                row["det_" + gene] = None
                continue
            values = np.asarray(adata[:, gene].X.todense()).ravel() > 0
            row["det_" + gene] = round(float(values.mean()), 4)
        rows.append(row)
        print(entry["library"], adata.n_obs)
        del adata
        gc.collect()

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "e4_detection_by_library.csv", index=False)
    rec.add_output(OUT / "e4_detection_by_library.csv")

    bleo = table[table["group"] == "bleomycin"]
    untr = table[table["group"] == "untreated"]
    verdicts = []
    for gene in FIBROTIC + OTHER:
        column = "det_" + gene
        values_b = bleo[column].dropna().tolist()
        values_u = untr[column].dropna().tolist()
        if len(values_b) < 2 or len(values_u) < 2:
            verdicts.append({"gene": gene, "verdict": "not evaluable"})
            continue
        generic = min(values_b) > max(values_u)
        verdicts.append({"gene": gene, "bleomycin": values_b, "untreated": values_u,
                         "injury_generic": bool(generic)})
    frame = pd.DataFrame(verdicts)
    frame.to_csv(OUT / "e4_injury_generic_verdicts.csv", index=False)
    rec.add_output(OUT / "e4_injury_generic_verdicts.csv")
    generic_set = [v["gene"] for v in verdicts if v.get("injury_generic")]
    rec.set("injury_generic_genes", generic_set)
    runx1_pdgfrb = all(g in generic_set for g in ("Runx1", "Pdgfrb"))
    rec.set("runx1_and_pdgfrb_injury_generic", bool(runx1_pdgfrb))
    fst_runx2 = all(g in generic_set for g in ("Fst", "Runx2"))
    if runx1_pdgfrb:
        reading = RULES["reading_fixed_in_advance"]["runx1_pdgfrb_injury_generic"]
    elif fst_runx2:
        reading = RULES["reading_fixed_in_advance"]["runx1_pdgfrb_not_generic_while_fst_runx2_are"]
    else:
        reading = "neither pre-specified pattern; no reading was fixed for this outcome"
    rec.set("reading", reading)

    lines = ["# Trial E4: bleomycin injury and the fibrotic gene set", "",
             "Injury-generic genes: " + (", ".join(generic_set) or "none"), "",
             "## Detection per library", "", df_to_markdown(table, index=False), "",
             "## Verdicts", "", df_to_markdown(frame, index=False), ""]
    (OUT / "e4_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e4_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
