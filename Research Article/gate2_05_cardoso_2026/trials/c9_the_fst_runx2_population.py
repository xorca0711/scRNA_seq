#!/usr/bin/env python
"""Trial C9: is there a co-expressing Fst and Runx2 fibroblast population, and what marks it.

Trial C8 produced one observation post hoc, and this trial is the
pre-registered version of it. In the Areg-flox/+ niche library, Fst and Runx2
co-occur in 13.78 per cent of gated fibroblasts against 7.94 per cent expected
under independence, a ratio of 1.735, so they mark the same cells. In the
Areg-flox/flox library that structure is gone, ratio 0.816 on 24
double-positive cells of 3,700, while the retained pair Runx1 and Pdgfrb sits
inside the independence band in both arms. C8 labelled that post hoc and
refused to claim it.

Two things have to be separated, and only one of them is answerable.

**Answerable.** Does a co-expressing Fst and Runx2 fibroblast population exist
at all, in an independent dataset with more than one animal per group, and if
so what else marks it. GSE132771 (Tsukui et al. 2020) has two bleomycin and
two untreated animals of sorted Col1a1-GFP mesenchyme with no oncogene, which
is the dataset trial E4 used, so the existence question can be asked with the
animal as the unit and a both-replicates-agree rule.

**Not answerable, and stated rather than approximated.** Whether Areg deletion
depletes that population. The Cardoso deposit has one library per genotype and
there is no other Areg-flox fibroblast dataset, so the depletion that motivated
this trial cannot be tested by anything on disk. Everything this trial says
about the deletion arm is descriptive.

Frozen rules, set before any value was computed:

* Fibroblasts: the compartment gate used by trials C6 and C8, Col1a1 detected
  and Ptprc, Pecam1 and Epcam not detected. Quality control and doublet removal
  follow trial C1.
* Libraries: the two GSE316244 niche libraries (Areg-flox/+ and
  Areg-flox/flox), and the four GSE132771 Col1a1-GFP-positive libraries
  (Bleo1, Bleo2, UT1, UT2).
* Double-positive means Fst and Runx2 both detected; double-negative means
  neither detected.

* T1, EXISTENCE, with two controls. In each library, the co-detection ratio of
  Fst with Runx2 against independence, and a permutation null: cells are binned
  into deciles of genes per cell, the Runx2 indicator is permuted within each
  bin, and the observed co-detection is compared with 200 such permutations.
  The permutation holds depth structure fixed, so it answers whether
  co-occurrence exceeds chance for cells of like depth. ITS UNIT IS THE CELL
  and it is therefore a statement inside one animal, not a group claim.
  THE GROUP RULE: the population is called co-organised in injury only if the
  ratio is at least 1.25 AND the permutation p is below 0.05 in BOTH bleomycin
  libraries.
* T2, INJURY RESPONSE. The double-positive share of gated fibroblasts, per
  library. RULE, as in trial E4: injury-induced only if the share is higher in
  both bleomycin libraries than in both untreated libraries.
* T3, WHAT ELSE MARKS THEM, DESCRIPTIVE. In the Areg-flox/+ library, genes
  ranked by the difference in detection fraction between double-positive and
  double-negative fibroblasts, restricted to genes detected in at least 10 per
  cent of either group, with Fst and Runx2 themselves excluded. The top 30 are
  reported. This is one library and cells are not replicates, so it is a
  description and no P value is computed.
* T4, DO THOSE MARKERS REPLICATE. Each T3 gene is checked in both bleomycin
  libraries. A gene replicates if its detection difference has the same sign
  and at least half the magnitude in both. The reported quantity is how many of
  the 30 replicate, which is a statement across animals rather than across
  cells.
* T5, DEPTH CONTROL. In every library, the co-detection ratio is recomputed on
  the shallow and deep halves at the median genes per cell. A ratio that moves
  by more than 0.25 is depth-driven and is not read.

* THE READING, fixed in advance:
  - T1 passes in both bleomycin libraries: a co-expressing Fst and Runx2
    fibroblast population exists in injury without any oncogene, so the C8
    observation is not a one-library artefact, and T3 with T4 says what marks
    it;
  - T1 fails: the C8 observation does not replicate in an independent injury
    dataset and the lead closes, whatever T3 shows;
  - T2 met as well: the population is injury-induced rather than resident;
  - T5 failing for a library: that library's ratio is not read.
* Nothing here tests the Areg-deletion depletion. Two animals per group in
  GSE132771 and one library per genotype in GSE316244, so no P value is
  computed on any group comparison.
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

OUT = HERE / "c9_the_fst_runx2_population"
OUT.mkdir(exist_ok=True)

PAIR = ("Fst", "Runx2")
GATE_POSITIVE = ["Col1a1"]
GATE_NEGATIVE = ["Ptprc", "Pecam1", "Epcam"]
TOP_N = 30
DETECTION_FLOOR = 0.10
RATIO_FLOOR = 1.25
DEPTH_SHIFT = 0.25
N_PERMUTATIONS = 200
SEED = 0

LIBRARIES = [
    {"series": "GSE316244", "gsm": "GSM9447779", "library": "Expt3_Het_niche",
     "group": "Areg-flox/+", "root": RAW / "GSE316244" / "GSE316244_RAW", "suffix": "features"},
    {"series": "GSE316244", "gsm": "GSM9447780", "library": "Expt3_Hom_niche",
     "group": "Areg-flox/flox", "root": RAW / "GSE316244" / "GSE316244_RAW", "suffix": "features"},
    {"series": "GSE132771", "gsm": "GSM3891612", "library": "Bleo1_GFPp",
     "group": "bleomycin", "root": RAW / "GSE132771" / "GSE132771_RAW", "suffix": "genes"},
    {"series": "GSE132771", "gsm": "GSM3891613", "library": "Bleo2_GFPp",
     "group": "bleomycin", "root": RAW / "GSE132771" / "GSE132771_RAW", "suffix": "genes"},
    {"series": "GSE132771", "gsm": "GSM3891616", "library": "UT1_GFPp",
     "group": "untreated", "root": RAW / "GSE132771" / "GSE132771_RAW", "suffix": "genes"},
    {"series": "GSE132771", "gsm": "GSM3891617", "library": "UT2_GFPp",
     "group": "untreated", "root": RAW / "GSE132771" / "GSE132771_RAW", "suffix": "genes"},
]
REFERENCE = "Areg-flox/+"
BLEOMYCIN = "bleomycin"
UNTREATED = "untreated"

RULES = {
    "question": "the pre-registered version of trial C8's post hoc observation",
    "answerable": "does a co-expressing Fst and Runx2 fibroblast population exist in an independent "
                  "injury dataset with animals, and what else marks it",
    "not_answerable": "whether Areg deletion depletes it; one library per genotype and no other "
                      "Areg-flox fibroblast dataset exists, so this is stated rather than approximated",
    "fibroblast_gate": "Col1a1 detected, Ptprc and Pecam1 and Epcam not detected, as in C6 and C8",
    "qc": "identical to trial C1",
    "pair": list(PAIR),
    "T1": f"co-detection ratio against independence plus a permutation null within deciles of genes per "
          f"cell ({N_PERMUTATIONS} permutations, seed {SEED}); the permutation unit is the CELL, so it is a "
          f"statement inside one animal. GROUP RULE: co-organised only if ratio >= {RATIO_FLOOR} and "
          f"permutation p < 0.05 in BOTH bleomycin libraries",
    "T2": "double-positive share higher in both bleomycin than both untreated libraries, as in E4",
    "T3": f"descriptive: top {TOP_N} genes by detection difference between double-positive and "
          f"double-negative fibroblasts of the {REFERENCE} library, floor {DETECTION_FLOOR}, the pair itself "
          f"excluded; one library and cells are not replicates, so no P value",
    "T4": "a T3 gene replicates if its detection difference has the same sign and at least half the "
          "magnitude in both bleomycin libraries; the reported quantity is how many of the top genes replicate",
    "T5": f"depth control: the ratio must not move by more than {DEPTH_SHIFT} between the shallow and deep "
          f"halves at the median genes per cell, or that library's ratio is not read",
    "reading_fixed_in_advance": {
        "T1_passes": "a co-expressing population exists in injury without an oncogene; the C8 observation is "
                     "not a one-library artefact, and T3 with T4 says what marks it",
        "T1_fails": "the C8 observation does not replicate in an independent injury dataset and the lead "
                    "closes, whatever T3 shows",
        "T2_met": "the population is injury-induced rather than resident",
        "T5_fails": "that library's ratio is not read",
    },
    "unit": "the animal; two per group in GSE132771 and one library per genotype in GSE316244, so no P value "
            "on any group comparison",
}


def load(entry, rec):
    """QC, doublet removal and the fibroblast gate; returns the gated object."""
    import anndata as ad
    import scanpy as sc

    stem = entry["gsm"] + "_" + entry["library"] + "_"
    matrix = entry["root"] / (stem + "matrix.mtx.gz")
    features = entry["root"] / (stem + entry["suffix"] + ".tsv.gz")
    if not features.exists():
        features = entry["root"] / (stem + "features.tsv.gz")
    barcodes = entry["root"] / (stem + "barcodes.tsv.gz")
    for path in (matrix, features, barcodes):
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(matrix, features, barcodes)
    if "feature_type" not in var or var["feature_type"].isna().all():
        var["feature_type"] = "Gene Expression"
    keep = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    if keep.sum() == 0:
        keep = np.ones(len(var), dtype=bool)
    X, var = X[:, keep], var.loc[keep].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
    obs = pd.DataFrame({"library": entry["library"], "group": entry["group"]},
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

    def detected(gene):
        if gene not in adata.var_names:
            return np.zeros(adata.n_obs, dtype=bool)
        return np.asarray(adata[:, gene].X.todense()).ravel() > 0

    gate = np.ones(adata.n_obs, dtype=bool)
    for gene in GATE_POSITIVE:
        gate &= detected(gene)
    for gene in GATE_NEGATIVE:
        gate &= ~detected(gene)
    gated = adata[gate].copy()
    print(entry["library"], adata.n_obs, "cells,", gated.n_obs, "pass the fibroblast gate")
    del adata
    gc.collect()
    return gated


def indicators(gated):
    out = {}
    for gene in PAIR:
        if gene in gated.var_names:
            out[gene] = np.asarray(gated[:, gene].X.todense()).ravel() > 0
        else:
            out[gene] = np.zeros(gated.n_obs, dtype=bool)
    return out


def co_detection(a, b):
    observed = float((a & b).mean())
    expected = float(a.mean() * b.mean())
    return observed, expected, (observed / expected if expected > 0 else np.nan)


def permutation_p(a, b, depth, rng):
    """Permute b within deciles of depth, so co-occurrence is judged at like depth."""
    observed = float((a & b).mean())
    edges = np.quantile(depth, np.linspace(0, 1, 11))
    bins = np.clip(np.digitize(depth, edges[1:-1]), 0, 9)
    hits = 0
    for _ in range(N_PERMUTATIONS):
        shuffled = np.empty_like(b)
        for k in range(10):
            where = np.where(bins == k)[0]
            if len(where) == 0:
                continue
            shuffled[where] = rng.permutation(b[where])
        if float((a & shuffled).mean()) >= observed:
            hits += 1
    return (hits + 1) / (N_PERMUTATIONS + 1)


def detection_frame(gated, mask, floor_reference=None):
    import scipy.sparse as sp

    sub = gated.X[mask]
    values = np.asarray((sub > 0).mean(axis=0)).ravel() if sp.issparse(sub) else (np.asarray(sub) > 0).mean(axis=0)
    return pd.Series(values, index=gated.var_names)


def main():
    rec = RunRecord(OUT / "c9_run_record.json",
                    "C9 the co-expressing Fst and Runx2 fibroblast population", RULES)
    rng = np.random.default_rng(SEED)

    per_library, reference_ranking, bleo_detection = [], None, {}
    for entry in LIBRARIES:
        gated = load(entry, rec)
        ind = indicators(gated)
        a, b = ind[PAIR[0]], ind[PAIR[1]]
        depth = gated.obs["n_genes_by_counts"].to_numpy().astype(float)
        observed, expected, ratio = co_detection(a, b)
        p = permutation_p(a, b, depth, rng)

        median = float(np.median(depth))
        halves = {}
        for name, mask in (("shallow", depth <= median), ("deep", depth > median)):
            if mask.sum() < 50:
                halves[name] = np.nan
                continue
            halves[name] = co_detection(a[mask], b[mask])[2]
        shift = (abs(halves["deep"] - halves["shallow"])
                 if not (np.isnan(halves["deep"]) or np.isnan(halves["shallow"])) else np.nan)

        per_library.append({
            "series": entry["series"], "library": entry["library"], "group": entry["group"],
            "n_gated": int(gated.n_obs), "median_genes_per_cell": round(median, 1),
            "det_Fst": round(float(a.mean()), 4), "det_Runx2": round(float(b.mean()), 4),
            "double_positive_share": round(observed, 4),
            "expected_if_independent": round(expected, 4),
            "ratio": round(float(ratio), 4) if not np.isnan(ratio) else None,
            "permutation_p": round(float(p), 4),
            "n_double_positive": int((a & b).sum()),
            "ratio_shallow": round(float(halves["shallow"]), 4) if not np.isnan(halves["shallow"]) else None,
            "ratio_deep": round(float(halves["deep"]), 4) if not np.isnan(halves["deep"]) else None,
            "depth_shift": round(float(shift), 4) if not np.isnan(shift) else None,
            "depth_control_passes": bool(not np.isnan(shift) and shift <= DEPTH_SHIFT),
        })

        double_positive = a & b
        double_negative = (~a) & (~b)
        if entry["group"] == REFERENCE and double_positive.sum() >= 30:
            up = detection_frame(gated, double_positive)
            down = detection_frame(gated, double_negative)
            keep = (up >= DETECTION_FLOOR) | (down >= DETECTION_FLOOR)
            diff = (up[keep] - down[keep]).drop(labels=[g for g in PAIR if g in up.index],
                                                errors="ignore")
            reference_ranking = pd.DataFrame({
                "gene": diff.abs().sort_values(ascending=False).head(TOP_N).index})
            reference_ranking["det_double_positive"] = up.reindex(reference_ranking["gene"]).round(4).to_numpy()
            reference_ranking["det_double_negative"] = down.reindex(reference_ranking["gene"]).round(4).to_numpy()
            reference_ranking["difference"] = diff.reindex(reference_ranking["gene"]).round(4).to_numpy()
        if entry["group"] == BLEOMYCIN:
            bleo_detection[entry["library"]] = {
                "up": detection_frame(gated, double_positive),
                "down": detection_frame(gated, double_negative),
                "n_double_positive": int(double_positive.sum()),
            }
        del gated
        gc.collect()

    table = pd.DataFrame(per_library)
    table.to_csv(OUT / "c9_per_library.csv", index=False)
    rec.add_output(OUT / "c9_per_library.csv")

    bleo = table[table["group"] == BLEOMYCIN]
    untr = table[table["group"] == UNTREATED]
    t1_pass = bool(len(bleo) == 2
                   and (bleo["ratio"] >= RATIO_FLOOR).all()
                   and (bleo["permutation_p"] < 0.05).all()
                   and bleo["depth_control_passes"].all())
    rec.set("T1", {"bleomycin_ratios": bleo["ratio"].tolist(),
                   "bleomycin_permutation_p": bleo["permutation_p"].tolist(),
                   "depth_controls_pass": bleo["depth_control_passes"].tolist(),
                   "passes": t1_pass})
    t2_pass = bool(len(bleo) == 2 and len(untr) == 2
                   and bleo["double_positive_share"].min() > untr["double_positive_share"].max())
    rec.set("T2", {"bleomycin_shares": bleo["double_positive_share"].tolist(),
                   "untreated_shares": untr["double_positive_share"].tolist(),
                   "injury_induced": t2_pass})

    replication = []
    if reference_ranking is not None and len(bleo_detection) == 2:
        for _, row in reference_ranking.iterrows():
            gene = row["gene"]
            target = float(row["difference"])
            checks = []
            for library, blob in bleo_detection.items():
                if gene not in blob["up"].index:
                    checks.append(None)
                    continue
                value = float(blob["up"][gene] - blob["down"][gene])
                checks.append(value)
            ok = all(v is not None and np.sign(v) == np.sign(target)
                     and abs(v) >= abs(target) / 2 for v in checks)
            entry = {"gene": gene, "reference_difference": round(target, 4), "replicates": bool(ok)}
            for (library, _), value in zip(bleo_detection.items(), checks):
                entry["diff_" + library] = round(value, 4) if value is not None else None
            replication.append(entry)
    rep = pd.DataFrame(replication)
    if len(rep):
        rep.to_csv(OUT / "c9_marker_replication.csv", index=False)
        rec.add_output(OUT / "c9_marker_replication.csv")
        rec.set("T4", {"genes_checked": int(len(rep)),
                       "genes_replicating": int(rep["replicates"].sum())})
    if reference_ranking is not None:
        reference_ranking.to_csv(OUT / "c9_reference_markers.csv", index=False)
        rec.add_output(OUT / "c9_reference_markers.csv")

    if t1_pass:
        reading = RULES["reading_fixed_in_advance"]["T1_passes"]
    else:
        reading = RULES["reading_fixed_in_advance"]["T1_fails"]
    if t2_pass:
        reading += " " + RULES["reading_fixed_in_advance"]["T2_met"]
    failed_depth = table.loc[~table["depth_control_passes"], "library"].tolist()
    if failed_depth:
        reading += (" Depth control failed for " + ", ".join(failed_depth)
                    + ", whose ratio is not read.")
    rec.set("reading", reading)

    lines = ["# Trial C9: the co-expressing Fst and Runx2 fibroblast population", "",
             "T1, existence in injury with both controls: " + ("PASSES" if t1_pass else "FAILS") + ".",
             "T2, injury-induced by the both-replicates rule: " + str(t2_pass) + ".",
             ("T4, markers from the reference arm that replicate in both bleomycin animals: "
              + str(int(rep["replicates"].sum())) + " of " + str(len(rep)) + ".") if len(rep)
             else "T4 not evaluable.", "",
             "Reading: " + reading, "",
             "**Not tested here:** whether Areg deletion depletes this population. One library per "
             "genotype, and no other Areg-flox fibroblast dataset exists.", "",
             "## Per library", "",
             df_to_markdown(table[["series", "library", "group", "n_gated", "median_genes_per_cell",
                                   "det_Fst", "det_Runx2", "double_positive_share",
                                   "expected_if_independent", "ratio", "permutation_p",
                                   "n_double_positive", "depth_shift", "depth_control_passes"]],
                            index=False), ""]
    if reference_ranking is not None:
        lines += ["## What else marks them in the Areg-flox/+ arm (descriptive, one library)", "",
                  df_to_markdown(reference_ranking.head(20), index=False), ""]
    if len(rep):
        lines += ["## Whether those markers replicate across the two bleomycin animals", "",
                  df_to_markdown(rep.head(20), index=False), ""]
    (OUT / "c9_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c9_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
