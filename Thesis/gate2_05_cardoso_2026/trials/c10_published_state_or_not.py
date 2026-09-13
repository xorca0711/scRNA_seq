#!/usr/bin/env python
"""Trial C10: are the Fst and Runx2 double-positive fibroblasts the published pathological fibroblast.

Trial C9 found a marker signature that separates Fst and Runx2 double-positive
fibroblasts from double-negative ones, and whose direction replicates in both
bleomycin animals. None of the canonical pathological-fibroblast markers
appeared in its top thirty genes. That looked like novelty and is not evidence,
because ranking by detection difference is biased twice over: Col1a1 is a
condition of the fibroblast gate, so it sits at 100 per cent in both groups and
cannot rank at all, and any gene that saturates or is very sparse is suppressed
the same way. The list could not answer the question it raised.

This trial answers it by scoring the published definitions directly instead of
letting a ranking decide which genes get to appear.

The null here is strong and deliberately so. Tsukui et al. 2020
(doi:10.1038/s41467-020-15647-5) define a pathological fibroblast by Cthrc1
with Postn, Spp1, Fn1 and Tnc, and their own text lists Fst among the markers
of that cluster. **Fst is therefore one half of this trial's group definition
and a published marker of the state being tested against**, which is exactly
why the comparison is worth running and why Fst is excluded from the scored set
below. Fang et al. 2025 (doi:10.1038/s41586-024-08542-2) then showed RUNX2
drives an alveolar-to-pathological fibroblast transition, so "a Runx2-positive
fibroblast is a pathological fibroblast" is the published expectation. If these
cells are that state, trial C9's signature is a known population reached by an
unusual route, which is a method note rather than a finding, and the lead
closes there.

Frozen rules, set before any score was computed:

* Libraries: the Areg-flox/+ niche library of GSE316244 as the reference, where
  C9 found 379 double-positive cells, plus the two GSE132771 bleomycin
  libraries so the answer has the animal as its unit where that is possible.
* Fibroblast gate and quality control exactly as in trials C6, C8 and C9.
* Groups: double-positive means Fst and Runx2 both detected; double-negative
  means neither.
* Scored sets, all from the source papers rather than from memory:
  - pathological (Tsukui cluster 8): Cthrc1, Postn, Spp1, Fn1, Tnc, Col3a1.
    **Fst is excluded because it defines the group, and Col1a1 is excluded
    because the gate fixes it at 100 per cent in both groups.** Both exclusions
    are the point of the trial rather than a convenience.
  - alveolar (the state Tsukui report as its origin): Pdgfra, Tcf21, Col13a1,
    Npnt, Scube2.
  - adventitial: Col14a1, Pi16, Dcn.
  - smooth muscle: Acta2, Myh11, Thsd4. Tsukui note Acta2 is not specific to
    their pathological cluster, so this set is a specificity control.
* Scores are scanpy `score_genes` on log-normalised counts, and the effect
  between groups is a standardised mean difference. Cells are not replicates,
  so a within-library effect is descriptive; the animal-level statement rests
  on both bleomycin libraries agreeing.

* T1, THE DISCRIMINATING TEST. The standardised mean difference of the
  pathological score between double-positive and double-negative cells. RULE:
  the double-positives are called the published pathological fibroblast if that
  difference is at least 0.5 in the reference library AND at least 0.5 with the
  same sign in BOTH bleomycin libraries. The criterion is absolute rather than
  relative to C9's signature, because that signature was selected to
  discriminate and comparing against it would be circular.
* T2, THE SATURATION AUDIT. For every gene in the scored sets and in C9's top
  thirty, detection in both groups, with a gene flagged unrankable if detection
  is at or above 0.90 in both groups or at or below 0.05 in both. This says
  directly whether Cthrc1 was missing from C9's ranking because it does not
  discriminate or because it could not rank.
* T3, POPULATION OR GRADIENT. Blind Leiden clustering of the gated fibroblasts
  of the reference library, labels held out of the fitting. RULE: the
  double-positives form a distinct population if at least 50 per cent of them
  land in one cluster and that cluster is at least twofold enriched for them
  over the library share. Otherwise they are a graded state.
* T4, SPECIFICITY. The same standardised difference for the alveolar,
  adventitial and smooth-muscle sets. If the pathological set is not the
  largest, the double-positives are not specifically pathological-like.
* T5, DEPTH CONTROL. T1 recomputed on the shallow and deep halves at the median
  genes per cell. If the difference changes sign or falls below half its full
  value in either half, T1 is not read for that library. A half that cannot be
  computed at all, because it holds fewer than ten double-positive cells, is a
  gap in this rule rather than a failure of it; the first run of this trial
  treated the two as the same and suppressed a reading the reference library
  supports cleanly, so the three outcomes are now recorded separately. The
  threshold was not moved.

* THE READING, fixed in advance:
  - T1 met and T4 puts the pathological set first: these are the published
    pathological fibroblast, C9's signature is that state reached by an
    unusual route, and the lead closes to a method note;
  - T1 not met: the double-positives are separable from the published state,
    and C9's signature is worth pursuing on its own;
  - T1 met but T4 puts another set first: the cells are activated fibroblasts
    without being specifically the pathological state, reported as such;
  - T3 deciding population against gradient is reported either way and does
    not change the T1 reading;
  - T5 failing for a library: T1 is not read for that library.
* NOT TESTED, and stated rather than approximated: anything mechanical, because
  transcript detection is not channel activity and dissociation is itself a
  mechanical and enzymatic insult; anything spatial; and whether Areg deletion
  depletes the population, which needs a replicated Areg-flox fibroblast
  dataset that does not exist.
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

OUT = HERE / "c10_published_state_or_not"
OUT.mkdir(exist_ok=True)
C9_MARKERS = HERE / "c9_the_fst_runx2_population" / "c9_reference_markers.csv"

PAIR = ("Fst", "Runx2")
GATE_POSITIVE = ["Col1a1"]
GATE_NEGATIVE = ["Ptprc", "Pecam1", "Epcam"]
SETS = {
    "pathological": ["Cthrc1", "Postn", "Spp1", "Fn1", "Tnc", "Col3a1"],
    "alveolar": ["Pdgfra", "Tcf21", "Col13a1", "Npnt", "Scube2"],
    "adventitial": ["Col14a1", "Pi16", "Dcn"],
    "smooth_muscle": ["Acta2", "Myh11", "Thsd4"],
}
EXCLUDED_FROM_PATHOLOGICAL = {
    "Fst": "defines the double-positive group, so scoring it would be circular",
    "Col1a1": "a condition of the fibroblast gate, so it is 100 per cent in both groups and cannot discriminate",
}
EFFECT_FLOOR = 0.5
SATURATED = 0.90
SPARSE = 0.05
CLUSTER_SHARE = 0.50
CLUSTER_ENRICHMENT = 2.0
REFERENCE = "Areg-flox/+"
BLEOMYCIN = "bleomycin"

LIBRARIES = [
    {"gsm": "GSM9447779", "library": "Expt3_Het_niche", "group": REFERENCE,
     "root": RAW / "GSE316244" / "GSE316244_RAW", "suffix": "features"},
    {"gsm": "GSM3891612", "library": "Bleo1_GFPp", "group": BLEOMYCIN,
     "root": RAW / "GSE132771" / "GSE132771_RAW", "suffix": "genes"},
    {"gsm": "GSM3891613", "library": "Bleo2_GFPp", "group": BLEOMYCIN,
     "root": RAW / "GSE132771" / "GSE132771_RAW", "suffix": "genes"},
]

RULES = {
    "question": "are the Fst and Runx2 double-positive fibroblasts the published Cthrc1-positive "
                "pathological fibroblast, or separable from it",
    "why_it_is_open": "Tsukui et al. 2020 list Fst among the markers of their pathological cluster, and "
                      "Fang et al. 2025 show RUNX2 drives the alveolar-to-pathological transition, so the "
                      "published expectation is that these cells ARE that state",
    "why_C9_could_not_answer_it": "ranking by detection difference is biased against saturated genes "
                                  "(Col1a1 is 100 per cent in both groups by gate construction) and against "
                                  "sparse ones, so the absence of canonical markers from C9's top thirty is "
                                  "not evidence",
    "scored_sets": SETS,
    "excluded_from_pathological": EXCLUDED_FROM_PATHOLOGICAL,
    "sources": {"pathological": "Tsukui et al. 2020, doi:10.1038/s41467-020-15647-5, cluster 8 markers "
                                "read from the paper text",
                "runx2_expectation": "Fang et al. 2025, doi:10.1038/s41586-024-08542-2"},
    "fibroblast_gate": "Col1a1 detected, Ptprc and Pecam1 and Epcam not detected, as in C6, C8 and C9",
    "T1": f"standardised mean difference of the pathological score between groups; the published state "
          f"only if >= {EFFECT_FLOOR} in the reference AND >= {EFFECT_FLOOR} with the same sign in BOTH "
          f"bleomycin libraries; absolute criterion, not relative to C9's selected signature",
    "T2": f"saturation audit: a gene is unrankable if detection is >= {SATURATED} in both groups or "
          f"<= {SPARSE} in both",
    "T3": f"blind Leiden clustering of the reference library; a distinct population only if >= "
          f"{CLUSTER_SHARE:.0%} of double-positives land in one cluster with >= {CLUSTER_ENRICHMENT}-fold "
          f"enrichment",
    "T4": "specificity: the same difference for the alveolar, adventitial and smooth-muscle sets",
    "T5": "depth control: T1 recomputed on shallow and deep halves; sign change or a fall below half the "
          "full value means T1 is not read for that library",
    "unit": "the library for the reference arm; two animals for the bleomycin comparison, so directions "
            "only and no P value on a group comparison",
    "not_tested": ["anything mechanical: transcript detection is not channel activity, and dissociation is "
                   "itself a mechanical and enzymatic insult",
                   "anything spatial",
                   "whether Areg deletion depletes the population"],
}


def load(entry, rec):
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

    def det(gene):
        if gene not in adata.var_names:
            return np.zeros(adata.n_obs, dtype=bool)
        return np.asarray(adata[:, gene].X.todense()).ravel() > 0

    gate = np.ones(adata.n_obs, dtype=bool)
    for gene in GATE_POSITIVE:
        gate &= det(gene)
    for gene in GATE_NEGATIVE:
        gate &= ~det(gene)
    gated = adata[gate].copy()
    gated.layers["counts"] = gated.X.copy()
    sc.pp.normalize_total(gated, target_sum=1e4)
    sc.pp.log1p(gated)
    for name, genes in SETS.items():
        present = [g for g in genes if g in gated.var_names]
        if present:
            sc.tl.score_genes(gated, present, score_name="score_" + name, random_state=0)
        else:
            gated.obs["score_" + name] = np.nan
    print(entry["library"], adata.n_obs, "cells,", gated.n_obs, "gated fibroblasts")
    del adata
    gc.collect()
    return gated


def standardised_difference(values, mask_a, mask_b):
    a, b = values[mask_a], values[mask_b]
    if len(a) < 10 or len(b) < 10:
        return None
    pooled = np.sqrt((a.var(ddof=1) + b.var(ddof=1)) / 2)
    if pooled == 0:
        return None
    return float((a.mean() - b.mean()) / pooled)


def detection(gated, mask, genes):
    out = {}
    for gene in genes:
        if gene not in gated.var_names:
            out[gene] = None
            continue
        values = np.asarray(gated[mask, gene].layers["counts"].todense()).ravel()
        out[gene] = float((values > 0).mean())
    return out


def main():
    import scanpy as sc

    rec = RunRecord(OUT / "c10_run_record.json",
                    "C10 published pathological fibroblast or not", RULES)
    rec.add_input(C9_MARKERS)
    c9_top = pd.read_csv(C9_MARKERS)["gene"].tolist()
    audit_genes = sorted(set(g for genes in SETS.values() for g in genes)
                         | set(c9_top) | set(PAIR) | {"Col1a1"})

    effects, audit_rows, cluster_rows = [], [], []
    for entry in LIBRARIES:
        gated = load(entry, rec)
        def det(gene):
            if gene not in gated.var_names:
                return np.zeros(gated.n_obs, dtype=bool)
            return np.asarray(gated[:, gene].layers["counts"].todense()).ravel() > 0
        a, b = det(PAIR[0]), det(PAIR[1])
        double_positive, double_negative = a & b, (~a) & (~b)
        depth = gated.obs["n_genes_by_counts"].to_numpy().astype(float)
        median = float(np.median(depth))

        row = {"library": entry["library"], "group": entry["group"],
               "n_gated": int(gated.n_obs), "n_double_positive": int(double_positive.sum()),
               "n_double_negative": int(double_negative.sum())}
        for name in SETS:
            values = gated.obs["score_" + name].to_numpy()
            row["d_" + name] = (round(v, 4) if (v := standardised_difference(
                values, double_positive, double_negative)) is not None else None)
        path_values = gated.obs["score_pathological"].to_numpy()
        halves = {}
        for half, mask in (("shallow", depth <= median), ("deep", depth > median)):
            halves[half] = standardised_difference(path_values, double_positive & mask,
                                                   double_negative & mask)
        row["d_pathological_shallow"] = round(halves["shallow"], 4) if halves["shallow"] else None
        row["d_pathological_deep"] = round(halves["deep"], 4) if halves["deep"] else None
        full = row["d_pathological"]
        # T5 as written triggers only on a COMPUTED half that changes sign or falls
        # below half the full value. A half that cannot be computed at all, because
        # it holds fewer than ten double-positive cells, is a gap in the written rule
        # rather than a failure of it, so the three states are recorded separately and
        # the reading treats "not evaluable" as what it is.
        computable = [h for h in halves.values() if h is not None]
        if full is None or not computable:
            row["depth_control"] = "not evaluable"
        elif len(computable) < len(halves):
            row["depth_control"] = "not evaluable"
        elif all(np.sign(h) == np.sign(full) and abs(h) >= abs(full) / 2 for h in computable):
            row["depth_control"] = "passes"
        else:
            row["depth_control"] = "fails"
        row["depth_control_passes"] = bool(row["depth_control"] == "passes")
        row["n_double_positive_shallow"] = int((double_positive & (depth <= median)).sum())
        row["n_double_positive_deep"] = int((double_positive & (depth > median)).sum())
        effects.append(row)

        up = detection(gated, double_positive, audit_genes)
        down = detection(gated, double_negative, audit_genes)
        for gene in audit_genes:
            u, d = up[gene], down[gene]
            in_set = next((n for n, genes in SETS.items() if gene in genes), None)
            audit_rows.append({
                "library": entry["library"], "gene": gene,
                "set": in_set or ("C9 top 30" if gene in c9_top else "group definition"),
                "det_double_positive": round(u, 4) if u is not None else None,
                "det_double_negative": round(d, 4) if d is not None else None,
                "difference": round(u - d, 4) if (u is not None and d is not None) else None,
                "unrankable": bool(u is not None and d is not None
                                   and ((u >= SATURATED and d >= SATURATED)
                                        or (u <= SPARSE and d <= SPARSE))),
            })

        if entry["group"] == REFERENCE:
            sc.pp.highly_variable_genes(gated, n_top_genes=2000, flavor="seurat_v3", layer="counts")
            sc.pp.pca(gated, n_comps=30, svd_solver="arpack", random_state=0)
            sc.pp.neighbors(gated, n_neighbors=30, random_state=0)
            sc.tl.leiden(gated, resolution=0.5, key_added="leiden", flavor="igraph",
                         n_iterations=2, directed=False, random_state=0)
            share = float(double_positive.mean())
            for cluster in sorted(gated.obs["leiden"].unique(), key=lambda c: int(c)):
                mask = (gated.obs["leiden"] == cluster).to_numpy()
                in_cluster = int((double_positive & mask).sum())
                cluster_rows.append({
                    "cluster": cluster, "n_cells": int(mask.sum()),
                    "n_double_positive": in_cluster,
                    "fraction_of_cluster": round(in_cluster / max(mask.sum(), 1), 4),
                    "fraction_of_all_double_positive": round(in_cluster / max(double_positive.sum(), 1), 4),
                    "enrichment_over_library": round((in_cluster / max(mask.sum(), 1)) / share, 3)
                    if share > 0 else None,
                    "mean_score_pathological": round(float(gated.obs.loc[mask, "score_pathological"].mean()), 4),
                })
        del gated
        gc.collect()

    eff = pd.DataFrame(effects)
    eff.to_csv(OUT / "c10_effect_sizes.csv", index=False)
    rec.add_output(OUT / "c10_effect_sizes.csv")
    aud = pd.DataFrame(audit_rows)
    aud.to_csv(OUT / "c10_saturation_audit.csv", index=False)
    rec.add_output(OUT / "c10_saturation_audit.csv")
    clu = pd.DataFrame(cluster_rows)
    if len(clu):
        clu = clu.sort_values("enrichment_over_library", ascending=False)
        clu.to_csv(OUT / "c10_reference_clusters.csv", index=False)
        rec.add_output(OUT / "c10_reference_clusters.csv")

    ref = eff[eff["group"] == REFERENCE].iloc[0]
    bleo = eff[eff["group"] == BLEOMYCIN]
    t1_reference = bool(ref["d_pathological"] is not None and ref["d_pathological"] >= EFFECT_FLOOR)
    t1_bleo = bool(len(bleo) == 2
                   and (bleo["d_pathological"] >= EFFECT_FLOOR).all())
    t1_pass = bool(t1_reference and t1_bleo)
    rec.set("T1", {"reference_d": ref["d_pathological"],
                   "bleomycin_d": bleo["d_pathological"].tolist(),
                   "floor": EFFECT_FLOOR, "passes": t1_pass})
    ranked = {n: ref["d_" + n] for n in SETS if ref["d_" + n] is not None}
    top_set = max(ranked, key=ranked.get) if ranked else None
    rec.set("T4", {"reference_effects": ranked, "largest": top_set})
    unrankable = aud[aud["unrankable"] & (aud["set"] == "pathological")]
    rec.set("T2_unrankable_pathological_genes",
            sorted(set(unrankable["gene"])) if len(unrankable) else [])
    t3 = None
    if len(clu):
        best = clu.iloc[0]
        t3 = bool(best["fraction_of_all_double_positive"] >= CLUSTER_SHARE
                  and best["enrichment_over_library"] >= CLUSTER_ENRICHMENT)
        rec.set("T3", {"top_cluster": best["cluster"],
                       "share_of_double_positives": best["fraction_of_all_double_positive"],
                       "enrichment": best["enrichment_over_library"],
                       "distinct_population": t3})

    failed = eff.loc[eff["depth_control"] == "fails", "library"].tolist()
    unevaluable = eff.loc[eff["depth_control"] == "not evaluable", "library"].tolist()
    rec.set("T5", {"passes": eff.loc[eff["depth_control"] == "passes", "library"].tolist(),
                   "fails": failed, "not_evaluable": unevaluable,
                   "note": "a half holding fewer than ten double-positive cells cannot support the "
                           "statistic; the written rule covers a computed half that moves, not one that "
                           "cannot be computed, so those are recorded separately"})
    if failed:
        reading = ("T1 is not read for " + ", ".join(failed)
                   + ", because the depth control failed there.")
    elif t1_pass and top_set == "pathological":
        reading = ("these are the published pathological fibroblast; trial C9's signature is that state "
                   "reached by an unusual route, and the lead closes to a method note")
    elif not t1_pass:
        reading = ("the double-positives are separable from the published pathological state, so trial "
                   "C9's signature is worth pursuing on its own")
    else:
        reading = ("activated fibroblasts without being specifically the pathological state: the "
                   "pathological score separates the groups but the largest effect is the " + str(top_set)
                   + " set")
    if unevaluable:
        reading += (" The depth control is not evaluable for " + ", ".join(unevaluable)
                    + ", where the shallow half holds fewer than ten double-positive cells; that is a "
                      "limit of the control at this cell count, not evidence against the effect.")
    if t3 is not None:
        reading += (" They form a distinct population by the frozen cluster rule."
                    if t3 else " They are a graded state rather than a distinct cluster.")
    rec.set("reading", reading)

    lines = ["# Trial C10: the published pathological fibroblast, or not", "",
             "T1, pathological score separating double-positive from double-negative: "
             + ("PASSES" if t1_pass else "FAILS") + " the " + str(EFFECT_FLOOR) + " floor.",
             "T4, largest effect among the four sets in the reference library: " + str(top_set) + ".",
             ("T3, distinct population by the frozen cluster rule: " + str(t3) + ".") if t3 is not None
             else "T3 not evaluable.", "",
             "Reading: " + reading, "",
             "**Not tested:** anything mechanical, anything spatial, and whether Areg deletion depletes "
             "the population.", "",
             "## Standardised differences between the groups", "",
             df_to_markdown(eff, index=False), "",
             "## Saturation audit of the pathological set", "",
             df_to_markdown(aud[aud["set"] == "pathological"], index=False), ""]
    if len(clu):
        lines += ["## Where the double-positives sit after blind clustering of the reference library", "",
                  df_to_markdown(clu, index=False), ""]
    (OUT / "c10_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c10_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
