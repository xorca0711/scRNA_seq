#!/usr/bin/env python
"""Trial C8: is the Areg-independent tier a subpopulation, or one programme at lower amplitude.

Item A2 of the next-step list in DIVERGENCES_AND_NEXT.md, with its meaning
changed by trial E4.

Trial C2 found that Areg deletion leaves Runx1 and Pdgfrb almost untouched
while Fst and Runx2 lose most of their detection, and trial C5 showed the
response is a three-tier gradient rather than a split. Trial E4 then removed
the interpretation that made the retained tier interesting: both genes rise
with bleomycin alone, so their persistence is fibroblast activation rather
than a tumour-specific second signal.

What survives that refutation is a structural question the earlier trials
never asked. Do Runx1 and Pdgfrb mark a distinct set of fibroblasts that
outlives Areg deletion, or does every fibroblast keep a little of the whole
programme so that the tiers are amplitudes rather than populations? The answer
changes how the C29 numbers should be described, and it is answerable from
data on disk.

A note on how the fibroblasts are defined here. Trial C2 defined them by
clustering. This trial defines them by a compartment marker gate instead, the
same gate trials C1b Question C and C6 used: Col1a1-positive and
Ptprc-negative and Pecam1-negative and Epcam-negative. That is a disclosed
difference, not an accident. The gate is transparent, it does not depend on a
clustering that would have to be reproduced, and the per-cell question here
does not need cluster boundaries. It also means the cell set is not identical
to C2's, so the detection fractions are not expected to match C2's exactly and
are recomputed here rather than quoted.

Frozen rules, set before any co-detection was computed:

* Libraries: the niche sorts of GSE316244, Expt3_Het_niche (Areg-flox/+) and
  Expt3_Hom_niche (Areg-flox/flox).
* Quality control and doublet removal follow trial C1 exactly.
* Fibroblast gate: Col1a1 detected, and Ptprc, Pecam1 and Epcam not detected.
* Genes: the paper's six reprogrammed-fibroblast markers. The retained tier is
  Runx1 and Pdgfrb, the falling tier is Fst and Runx2, as C5 measured them.
* T1, THE TEST. Within the Areg-flox/flox fibroblasts, the observed fraction
  of cells co-detecting Runx1 and Pdgfrb against the fraction expected if the
  two were independent. RULE, the same threshold trial C1c used: a ratio at or
  above 1.25 means the retained genes mark the same cells more than chance,
  which is what a surviving subpopulation looks like; a ratio at or below 0.80
  means they mark different cells; between 0.80 and 1.25 is consistent with
  independence, which is what a graded programme looks like.
* T2. The same ratio for the falling pair (Fst with Runx2) and for one mixed
  pair (Runx1 with Fst), in both genotypes, so the T1 number has something to
  be compared against.
* T3, STATE OR GRADIENT. A retained-tier score, the mean of the log-normalised
  Runx1 and Pdgfrb values, is fitted with a one-component and a two-component
  Gaussian mixture in the flox/flox fibroblasts, and the Bayesian information
  criterion decides. This is the same state-versus-gradient test trial C1 used.
  Two components winning means a discrete subpopulation; one component means a
  gradient.
* T4, THE CONTROL THAT DECIDES WHETHER T1 MEANS ANYTHING. Co-detection rises
  with sequencing depth for any pair of genes, so the same three ratios are
  recomputed after splitting the flox/flox fibroblasts at their median genes
  per cell. RULE: if a ratio moves by more than 0.25 between the shallow and
  deep halves, that ratio is depth-driven and is not read.
* THE READING, fixed in advance:
  - T1 at or above 1.25 and T3 favouring two components: the Areg-independent
    tier is a distinct subpopulation of fibroblasts, and claim C29 should be
    described as a population difference;
  - T1 between 0.80 and 1.25 and T3 favouring one component: the tiers are
    amplitudes of one programme, and C29 should be described as a gradient,
    which is how C5 already corrected it;
  - the two disagreeing: reported as unresolved, with both numbers;
  - T4 failing for the Runx1 and Pdgfrb ratio: no reading.
* One library per genotype and three mice pooled in each, so nothing here is
  tested between genotypes. Every number is within a single library.
"""

from __future__ import annotations

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

OUT = HERE / "c8_subpopulation_or_gradient"
OUT.mkdir(exist_ok=True)
ROOT = RAW / "GSE316244" / "GSE316244_RAW"
LIBRARIES = [
    {"gsm": "GSM9447779", "library": "Expt3_Het_niche", "genotype": "Areg-flox/+"},
    {"gsm": "GSM9447780", "library": "Expt3_Hom_niche", "genotype": "Areg-flox/flox"},
]
GATE_POSITIVE = ["Col1a1"]
GATE_NEGATIVE = ["Ptprc", "Pecam1", "Epcam"]
RETAINED = ("Runx1", "Pdgfrb")
FALLING = ("Fst", "Runx2")
MIXED = ("Runx1", "Fst")
ALL_GENES = ["Tnc", "Fst", "Runx1", "Runx2", "Acta2", "Pdgfrb"]
UPPER, LOWER = 1.25, 0.80
DEPTH_SHIFT = 0.25

RULES = {
    "question": "item A2, re-aimed after E4: is the Areg-independent tier a subpopulation or a gradient",
    "libraries": LIBRARIES,
    "qc": "identical to trial C1",
    "fibroblast_definition": "marker gate, Col1a1 positive and Ptprc, Pecam1, Epcam negative, as in C1b Question C and C6",
    "disclosed_difference": "C2 defined fibroblasts by clustering; the cell set here is not identical and detection fractions are recomputed, not quoted from C2",
    "genes": ALL_GENES,
    "tiers": {"retained": list(RETAINED), "falling": list(FALLING), "mixed": list(MIXED)},
    "T1": f"flox/flox co-detection of {RETAINED[0]} and {RETAINED[1]} against independence; "
          f">= {UPPER} same cells, <= {LOWER} different cells, between is independence",
    "T2": "the same ratio for the falling pair and one mixed pair, both genotypes",
    "T3": "one against two Gaussian components on the retained-tier score by BIC, as in C1",
    "T4": f"depth control: split flox/flox fibroblasts at the median genes per cell; a ratio moving "
          f"more than {DEPTH_SHIFT} between halves is depth-driven and is not read",
    "reading_fixed_in_advance": {
        "subpopulation": "the Areg-independent tier is a distinct subpopulation; C29 is a population difference",
        "gradient": "the tiers are amplitudes of one programme; C29 is a gradient, as C5 already corrected",
        "disagreement": "unresolved; both numbers reported",
        "depth_failed": "the ratio is depth-driven; no reading",
    },
    "unit": "the library; one per genotype with three mice pooled, so nothing is tested between genotypes",
}


def load(entry, rec):
    import anndata as ad
    import scanpy as sc

    stem = entry["gsm"] + "_" + entry["library"] + "_"
    matrix = ROOT / (stem + "matrix.mtx.gz")
    features = ROOT / (stem + "features.tsv.gz")
    if not features.exists():
        features = ROOT / (stem + "genes.tsv.gz")
    barcodes = ROOT / (stem + "barcodes.tsv.gz")
    for path in (matrix, features, barcodes):
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(matrix, features, barcodes)
    keep = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    if keep.sum() == 0:
        keep = np.ones(len(var), dtype=bool)
    X, var = X[:, keep], var.loc[keep].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
    obs = pd.DataFrame({"library": entry["library"], "genotype": entry["genotype"]},
                       index=pd.Index([entry["library"] + "_" + b for b in bc]))
    adata = ad.AnnData(X=X.astype(np.float32), var=var, obs=obs)
    adata.var["mt"] = adata.var["gene_symbol"].astype(str).str.startswith("mt-").to_numpy()
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=[20], log1p=True, inplace=True)
    th = derive_thresholds(adata.obs, entry["library"])
    adata = adata[apply_thresholds(adata.obs, th).to_numpy()].copy()
    rate = float(np.clip(0.008 * adata.n_obs / 1000, 0.01, 0.15))
    sc.pp.scrublet(adata, random_state=0, verbose=False, expected_doublet_rate=rate)
    auto = float(adata.obs["predicted_doublet"].to_numpy().mean())
    if auto < 0.2 * rate or auto > 3 * rate:
        cut = float(np.quantile(adata.obs["doublet_score"].to_numpy(), 1 - rate))
        adata.obs["predicted_doublet"] = adata.obs["doublet_score"].to_numpy() >= cut
    adata = adata[~adata.obs["predicted_doublet"].to_numpy()].copy()

    counts = {}
    for gene in set(GATE_POSITIVE + GATE_NEGATIVE + ALL_GENES):
        if gene in adata.var_names:
            counts[gene] = np.asarray(adata[:, gene].X.todense()).ravel()
        else:
            counts[gene] = np.zeros(adata.n_obs, dtype=np.float32)
    gate = np.ones(adata.n_obs, dtype=bool)
    for gene in GATE_POSITIVE:
        gate &= counts[gene] > 0
    for gene in GATE_NEGATIVE:
        gate &= counts[gene] == 0
    frame = pd.DataFrame({g: counts[g][gate] for g in ALL_GENES})
    frame["genes_per_cell"] = adata.obs["n_genes_by_counts"].to_numpy()[gate]
    frame["genotype"] = entry["genotype"]
    frame["library"] = entry["library"]
    print(entry["library"], adata.n_obs, "cells,", int(gate.sum()), "pass the fibroblast gate")
    return frame


def ratio(frame, pair):
    a = (frame[pair[0]] > 0).to_numpy()
    b = (frame[pair[1]] > 0).to_numpy()
    observed = float((a & b).mean())
    expected = float(a.mean() * b.mean())
    return {"pair": pair[0] + " with " + pair[1], "n_cells": int(len(frame)),
            "detection_a": round(float(a.mean()), 4), "detection_b": round(float(b.mean()), 4),
            "observed_both": round(observed, 4), "expected_if_independent": round(expected, 4),
            "ratio": round(observed / expected, 4) if expected > 0 else None}


def verdict(value):
    if value is None:
        return "not evaluable"
    if value >= UPPER:
        return "same cells more than chance"
    if value <= LOWER:
        return "different cells"
    return "consistent with independence"


def main():
    from sklearn.mixture import GaussianMixture

    rec = RunRecord(OUT / "c8_run_record.json",
                    "C8 is the Areg-independent tier a subpopulation or a gradient", RULES)
    frames = [load(entry, rec) for entry in LIBRARIES]
    cells = pd.concat(frames, ignore_index=True)
    cells.groupby(["genotype", "library"]).size().to_frame("n_gated_cells").to_csv(
        OUT / "c8_gated_cells.csv")
    rec.add_output(OUT / "c8_gated_cells.csv")
    rec.set("gated_cells", cells.groupby("genotype").size().to_dict())
    rec.set("median_genes_per_cell", cells.groupby("genotype")["genes_per_cell"].median().to_dict())

    rows = []
    for genotype, group in cells.groupby("genotype"):
        for pair, name in ((RETAINED, "retained"), (FALLING, "falling"), (MIXED, "mixed")):
            row = ratio(group, pair)
            row.update({"genotype": genotype, "tier": name, "verdict": verdict(row["ratio"])})
            rows.append(row)
    pairs = pd.DataFrame(rows)
    pairs.to_csv(OUT / "c8_codetection_ratios.csv", index=False)
    rec.add_output(OUT / "c8_codetection_ratios.csv")

    hom = cells[cells["genotype"] == "Areg-flox/flox"]
    t1 = pairs[(pairs["genotype"] == "Areg-flox/flox") & (pairs["tier"] == "retained")].iloc[0]
    rec.set("T1", {"ratio": t1["ratio"], "verdict": t1["verdict"],
                   "observed": t1["observed_both"], "expected": t1["expected_if_independent"]})

    values = np.log1p(hom[list(RETAINED)].to_numpy()).mean(axis=1).reshape(-1, 1)
    bics = {}
    for k in (1, 2):
        model = GaussianMixture(n_components=k, random_state=0, n_init=5).fit(values)
        bics[k] = round(float(model.bic(values)), 2)
    components = 1 if bics[1] <= bics[2] else 2
    rec.set("T3", {"bic": bics, "components_favoured": components,
                   "delta_bic": round(bics[1] - bics[2], 2)})

    median_depth = float(hom["genes_per_cell"].median())
    halves = []
    for name, part in (("shallow", hom[hom["genes_per_cell"] <= median_depth]),
                       ("deep", hom[hom["genes_per_cell"] > median_depth])):
        for pair, tier in ((RETAINED, "retained"), (FALLING, "falling"), (MIXED, "mixed")):
            row = ratio(part, pair)
            row.update({"half": name, "tier": tier,
                        "median_genes_per_cell": round(float(part["genes_per_cell"].median()), 1)})
            halves.append(row)
    depth = pd.DataFrame(halves)
    depth.to_csv(OUT / "c8_depth_control.csv", index=False)
    rec.add_output(OUT / "c8_depth_control.csv")
    shifts = {}
    for tier in ("retained", "falling", "mixed"):
        sub = depth[depth["tier"] == tier].set_index("half")["ratio"]
        if sub.isna().any() or len(sub) < 2:
            shifts[tier] = None
            continue
        shifts[tier] = round(float(abs(sub["deep"] - sub["shallow"])), 4)
    rec.set("T4_depth_shifts", shifts)
    depth_ok = shifts.get("retained") is not None and shifts["retained"] <= DEPTH_SHIFT

    if not depth_ok:
        reading = RULES["reading_fixed_in_advance"]["depth_failed"]
    elif t1["ratio"] is not None and t1["ratio"] >= UPPER and components == 2:
        reading = RULES["reading_fixed_in_advance"]["subpopulation"]
    elif t1["ratio"] is not None and LOWER < t1["ratio"] < UPPER and components == 1:
        reading = RULES["reading_fixed_in_advance"]["gradient"]
    else:
        reading = RULES["reading_fixed_in_advance"]["disagreement"]
    rec.set("reading", reading)

    lines = ["# Trial C8: subpopulation or gradient", "",
             "T1, flox/flox " + t1["pair"] + ": observed " + f"{t1['observed_both']:.4f}"
             + " against " + f"{t1['expected_if_independent']:.4f}" + " expected, ratio "
             + str(t1["ratio"]) + ", " + str(t1["verdict"]) + ".",
             "T3, mixture BIC: one component " + str(bics[1]) + ", two components " + str(bics[2])
             + ", favouring " + str(components) + ".",
             "T4, depth shift of that ratio: " + str(shifts.get("retained"))
             + " against a limit of " + str(DEPTH_SHIFT) + ", "
             + ("passes" if depth_ok else "FAILS") + ".", "",
             "Reading: " + reading, "",
             "## Co-detection against independence", "", df_to_markdown(pairs, index=False), "",
             "## Depth control, flox/flox split at the median genes per cell", "",
             df_to_markdown(depth, index=False), ""]
    (OUT / "c8_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c8_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
