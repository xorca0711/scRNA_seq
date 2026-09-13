#!/usr/bin/env python
"""Trial C6: which compartment makes each EGFR ligand, and does any source survive Areg deletion.

Why this trial exists. Trial C5 found Hbegf to be the second most abundant
EGFR ligand in the mutant DATP-like state, and the session then floated the
idea that Hbegf might be a brake on the Areg axis, because HB-EGF sends EGFR
to degradation while Areg recycles it. A literature check refuted that idea
for lung: HB-EGF is pro-fibrotic there, it is made by lung macrophages and by
transitional alveolar epithelial cells, and deleting it from the myeloid
compartment protects mice from bleomycin fibrosis (Hult et al., Am J Respir
Cell Mol Biol 2022, doi:10.1165/rcmb.2022-0174OC). The better hypothesis is
the opposite of a brake: **Hbegf may be a second, Areg-independent driver of
the same fibroblast programme**, which would explain the tier of the fibrotic
signature that trial C2 found surviving Areg deletion (Runx1 97%, Pdgfrb 86%).

That hypothesis makes a prediction this deposit can test without any new data:
if a non-epithelial compartment makes Hbegf, then deleting Areg from AT2 cells
cannot remove it, and an EGFR ligand source survives the genetic perturbation.

This trial therefore asks a plain question of the niche libraries: who makes
which EGFR ligand, and who carries the receptor. It does not cluster, fit or
integrate anything. Compartments are assigned by marker gates on raw counts,
which is a coarse but transparent assignment that cannot be tuned to the
answer.

Provenance: rules fixed after the C5 tables were seen and after the literature
check, and before any matrix was read in this trial. Post hoc and disclosed.

Frozen rules:

* Libraries: the two niche libraries of GSE316243 (Confetti, Red2Kras), the
  two niche and two RFP+ libraries of GSE316244 (Areg-flox/+, Areg-flox/flox),
  and the two mesenchymal libraries of GSE316241, so that every sorted
  compartment of the deposit is represented.
* Quality control and doublet removal follow trial C1 exactly. Non-Ensembl
  features are dropped by the C0 rule.
* Compartment gates on raw counts, applied in this fixed priority so that
  every cell lands in exactly one compartment: immune (Ptprc), then
  endothelial (Pecam1 or Cdh5), then epithelial (Epcam or Krt8 or Krt18),
  then mesenchymal (Col1a1 or Col1a2), then unassigned. Immune cells are
  further split, again by priority: neutrophil (S100a8, S100a9 or Retnlg),
  alveolar-macrophage-like (a myeloid marker with Siglecf, Car4 or Itgax),
  other myeloid (Lyz2, Cd68, Adgre1 or Itgam), lymphoid (Cd3e, Cd79a or
  Ncr1), other immune.
* A compartment is evaluable in a library only if it holds at least 100
  cells. Smaller ones are reported with their count and flagged.
* Reported per compartment and library: cells, and for Areg, Ereg, Hbegf,
  Tgfa and Egfr the detection fraction and the mean log1p(CP10K).
* THE TEST, frozen. Hbegf has a non-epithelial source if some non-epithelial
  compartment with at least 100 cells detects Hbegf in at least 10% of its
  cells. That source survives Areg deletion if the same compartment in the
  Areg-flox/flox niche library holds at least 70% of its flox/+ detection.
* Marker gates are not cell-type calls, sequencing depth differs between
  libraries, and the unit is the library. No P value anywhere.
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

OUT = HERE / "c6_who_makes_egfr_ligands"
OUT.mkdir(exist_ok=True)

LIGANDS = ["Areg", "Ereg", "Hbegf", "Tgfa"]
RECEPTOR = ["Egfr"]
LIBRARIES = [
    {"accession": "GSE316243", "gsm": "GSM9447777", "library": "Expt2_Confetti_niche",
     "arm": "Confetti", "sort": "niche (immune, stroma, endothelium)"},
    {"accession": "GSE316243", "gsm": "GSM9447778", "library": "Expt2_Red2Kras_niche",
     "arm": "Red2Kras", "sort": "niche (immune, stroma, endothelium)"},
    {"accession": "GSE316244", "gsm": "GSM9447779", "library": "Expt3_Het_niche",
     "arm": "Areg-flox/+", "sort": "niche (mesenchyme, immune)"},
    {"accession": "GSE316244", "gsm": "GSM9447780", "library": "Expt3_Hom_niche",
     "arm": "Areg-flox/flox", "sort": "niche (mesenchyme, immune)"},
    {"accession": "GSE316244", "gsm": "GSM9447781", "library": "Expt3_Het_RFP",
     "arm": "Areg-flox/+", "sort": "RFP+ epithelium"},
    {"accession": "GSE316244", "gsm": "GSM9447782", "library": "Expt3_Hom_RFP",
     "arm": "Areg-flox/flox", "sort": "RFP+ epithelium"},
    {"accession": "GSE316241", "gsm": "GSM9447763", "library": "Expt1_Confetti_mesenchyme",
     "arm": "Confetti", "sort": "mesenchyme"},
    {"accession": "GSE316241", "gsm": "GSM9447764", "library": "Expt1_Red2Kras_mesenchyme",
     "arm": "Red2Kras", "sort": "mesenchyme"},
]
GATES = {
    "immune": ["Ptprc"],
    "endothelial": ["Pecam1", "Cdh5"],
    "epithelial": ["Epcam", "Krt8", "Krt18"],
    "mesenchymal": ["Col1a1", "Col1a2"],
}
IMMUNE_SPLIT = {
    "neutrophil": ["S100a8", "S100a9", "Retnlg"],
    "alveolar macrophage-like": ["Siglecf", "Car4", "Itgax"],
    "other myeloid": ["Lyz2", "Cd68", "Adgre1", "Itgam"],
    "lymphoid": ["Cd3e", "Cd79a", "Ncr1"],
}
MYELOID_ANY = ["Lyz2", "Cd68", "Adgre1", "Itgam", "Fcgr3", "Csf1r"]
MIN_CELLS = 100
NON_EPITHELIAL_SOURCE_MIN = 0.10
SURVIVES_MIN = 0.70

RULES = {
    "provenance_of_these_rules": "fixed after the C5 tables and the literature check, before any matrix was read here",
    "hypothesis_tested": ("Hbegf is a second, Areg-independent driver rather than a brake; the refuted "
                          "brake idea and the literature that refutes it are in the module docstring"),
    "libraries": LIBRARIES,
    "qc": "identical to trial C1 (per-library MAD, Scrublet with the 10x prior and the sanity check)",
    "compartment_gates": GATES,
    "gate_priority": ["immune", "endothelial", "epithelial", "mesenchymal", "unassigned"],
    "immune_split": IMMUNE_SPLIT,
    "immune_split_priority": ["neutrophil", "alveolar macrophage-like", "other myeloid", "lymphoid", "other immune"],
    "myeloid_markers_for_the_am_gate": MYELOID_ANY,
    "min_cells_to_evaluate": MIN_CELLS,
    "non_epithelial_source_if_detection_at_least": NON_EPITHELIAL_SOURCE_MIN,
    "survives_areg_deletion_if_retains_at_least": SURVIVES_MIN,
    "genes": LIGANDS + RECEPTOR,
    "caveats": ["marker gates are not cell-type calls", "depth differs between libraries",
                "the unit is the library; no P value"],
}


def load(entry, rec):
    import anndata as ad
    import scanpy as sc

    root = RAW / entry["accession"] / (entry["accession"] + "_RAW")
    stem = entry["gsm"] + "_" + entry["library"] + "_"
    paths = [root / (stem + s) for s in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz")]
    for path in paths:
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(*paths)
    keep_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    X, var = X[:, keep_gene], var.loc[keep_gene].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
    obs = pd.DataFrame(index=pd.Index([entry["library"] + "_" + b for b in bc]))
    adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)
    symbols = adata.var["gene_symbol"].astype(str)
    adata.var["mt"] = symbols.str.startswith("mt-").to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True)
    adata.obs["pct_counts_ribo"] = 0.0
    th = derive_thresholds(adata.obs, entry["library"])
    adata = adata[apply_thresholds(adata.obs, th).to_numpy()].copy()
    rate = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=rate)
    auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
    if auto < 0.2 * rate or auto > 3 * rate:
        cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - rate))
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
    adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()
    return adata


def any_detected(adata, genes):
    """Boolean mask: at least one of these genes has a non-zero count."""
    present = [g for g in genes if g in adata.var_names]
    if not present:
        return np.zeros(adata.n_obs, dtype=bool)
    block = adata[:, present].X
    return np.asarray((block > 0).sum(axis=1)).ravel() > 0


def assign_compartment(adata):
    labels = np.full(adata.n_obs, "unassigned", dtype=object)
    for name in ("mesenchymal", "epithelial", "endothelial", "immune"):
        hit = any_detected(adata, GATES[name])
        labels[hit] = name
    return labels


def split_immune(adata, labels):
    myeloid = any_detected(adata, MYELOID_ANY)
    out = labels.copy()
    immune = labels == "immune"
    out[immune] = "other immune"
    for name in ("lymphoid", "other myeloid", "alveolar macrophage-like", "neutrophil"):
        hit = any_detected(adata, IMMUNE_SPLIT[name]) & immune
        if name == "alveolar macrophage-like":
            hit = hit & myeloid
        out[hit] = name
    return out


def gene_stats(adata, mask):
    total = np.asarray(adata.X.sum(axis=1)).ravel()
    out = {}
    for gene in LIGANDS + RECEPTOR:
        if gene not in adata.var_names:
            continue
        col = np.asarray(adata[:, gene].X.todense()).ravel()
        scale = total[mask]
        vals = np.log1p(col[mask] * 10000.0 / np.maximum(scale, 1.0))
        out["det_" + gene] = round(float((vals > 0).mean()), 4)
        out["mean_" + gene] = round(float(vals.mean()), 4)
    return out


def main():
    rec = RunRecord(OUT / "c6_run_record.json", "C6 who makes EGFR ligands", RULES)
    rows = []
    for entry in LIBRARIES:
        adata = load(entry, rec)
        labels = split_immune(adata, assign_compartment(adata))
        for name in sorted(set(labels)):
            mask = labels == name
            row = {"library": entry["library"], "arm": entry["arm"], "sort": entry["sort"]}
            row["compartment"] = name
            row["n_cells"] = int(mask.sum())
            row["evaluable"] = bool(mask.sum() >= MIN_CELLS)
            row.update(gene_stats(adata, mask))
            rows.append(row)
        print(entry["library"], adata.n_obs)
        del adata
        gc.collect()
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "c6_ligands_by_compartment.csv", index=False)
    rec.add_output(OUT / "c6_ligands_by_compartment.csv")
    ok = table[table["evaluable"]]
    sources = ok[(ok["compartment"] != "epithelial") & (ok["det_Hbegf"] >= NON_EPITHELIAL_SOURCE_MIN)]
    rec.set("non_epithelial_hbegf_sources", sources[["library", "compartment", "n_cells", "det_Hbegf"]].to_dict("records"))
    het = ok[ok["library"] == "Expt3_Het_niche"].set_index("compartment")
    hom = ok[ok["library"] == "Expt3_Hom_niche"].set_index("compartment")
    survival = []
    for name in het.index.intersection(hom.index):
        a = float(het.loc[name, "det_Hbegf"])
        b = float(hom.loc[name, "det_Hbegf"])
        keep = round(b / a, 3) if a > 0 else None
        row = {"compartment": name, "det_Hbegf_flox_plus": a, "det_Hbegf_flox_flox": b}
        row["retained_fraction"] = keep
        row["survives"] = bool(keep is not None and keep >= SURVIVES_MIN)
        survival.append(row)
    surv = pd.DataFrame(survival)
    surv.to_csv(OUT / "c6_hbegf_survival.csv", index=False)
    rec.add_output(OUT / "c6_hbegf_survival.csv")
    rec.set("hbegf_survival", survival)
    lines = ["# Trial C6: who makes EGFR ligands", ""]
    lines += ["## By compartment", "", df_to_markdown(table, index=False), ""]
    lines += ["## Hbegf after Areg deletion", "", df_to_markdown(surv, index=False), ""]
    (OUT / "c6_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c6_summary.md")
    rec.finish()
    print(df_to_markdown(surv, index=False))


if __name__ == "__main__":
    main()
