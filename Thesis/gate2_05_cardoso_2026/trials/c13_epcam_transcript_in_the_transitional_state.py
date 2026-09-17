#!/usr/bin/env python
"""Trial C13: does Epcam transcript explain why the transitional state escapes an EpCAM gate.

Trial C1d found 184 epithelial cells inside a CD45-CD31-EpCAM-negative
mesenchymal sort, 88 per cent of them Areg-positive, and trial C7 showed they
are a mixture of roughly four parts transitional state to four parts AT2. Two
explanations for the escape compete, and they differ in where the loss of
surface EpCAM happens.

**Transcriptional.** The transitional state simply expresses less Epcam, so it
stains dimly and falls out of the gate. Simple, and testable here.

**Post-transcriptional.** The cells carry Epcam message but shed the protein
from the surface. This is mechanistically attractive because the sheddase is
shared: ADAM17 releases amphiregulin (Sahin et al. 2004,
doi:10.1083/jcb.200307137) and also cleaves the EpCAM ectodomain with
presenilin-2 (Maetzel et al. 2009, doi:10.1038/ncb1824). A cell becoming an
amphiregulin source would then be the same cell going EpCAM-dim.

**What this trial can and cannot do, stated before it runs.** It cannot measure
surface protein, because three-prime counting sees transcript only, and it
cannot measure shedding at all. It can measure the transcript, and that is
enough to separate the two explanations in one direction: if Epcam message is
clearly lower in the transitional state, the transcriptional account suffices
and nothing post-transcriptional needs invoking. If the message is not lower,
then any surface dimming has to be post-transcriptional, which is consistent
with the shedding account without being evidence for it.

The sheddase genes are carried as exploratory context only, and their negative
would be close to uninformative. ADAM17 activity is controlled overwhelmingly
after translation, through iRhom-dependent trafficking, a protein-disulphide-
isomerase conformational switch and phosphorylation, so its mRNA is a poor
proxy for how much enzyme is working. That is a reason to report the numbers
and not to lean on them.

Frozen rules, set before any value was computed:

* Libraries. Primary: the four KrasG12D RFP libraries of GSE247505, which is
  the companion series with two replicates at each of two time points where the
  transitional state is abundant. Secondary: the two RFP-sorted epithelial
  libraries of GSE316244. The wild-type YFP libraries are reported for counts
  only, because trial C3 found the transitional state is nearly absent there.
* Quality control and doublet removal follow trial C1 exactly.
* States: the paper's own Figure 4l marker sets (AT2, DATP_like, AT1_like,
  Cd177_positive, cycling), scored with `score_genes` and assigned per cell by
  argmax. This differs from trial C3, which clustered first and then called
  clusters; per-cell assignment is used here because the comparison is between
  two states inside one library rather than between libraries, and the
  difference is disclosed rather than hidden.
* Groups: cells called DATP_like against cells called AT2, within a library. A
  library contributes only if both groups hold at least 50 cells.
* Genes: Epcam is the primary. Adam17, Rhbdf1, Rhbdf2 and Timp3 are the
  exploratory sheddase context, Timp3 being the endogenous ADAM17 inhibitor so
  its predicted direction is down. Areg, Krt8 and Cldn4 are sanity checks.
* Measures: detection fraction and mean log-normalised value, per library and
  group.

* T1, THE PRIMARY TEST. Epcam detection in DATP_like against AT2, within each
  library. RULE: Epcam is called transcriptionally lower only if it is lower in
  EVERY evaluable library AND the median drop is at least 10 percentage points.
* T2, EXPLORATORY CONTEXT. The same comparison for Adam17, Rhbdf1, Rhbdf2 and
  Timp3, called only if the direction agrees in every evaluable library.
  Reported with the post-translational caveat attached, never without it.
* T3, DEPTH CONTROL. Median genes per cell for each group in each library, and
  T1 recomputed on the shallow and deep halves of the library. If the two
  groups differ in depth and the T1 direction does not survive both halves, T1
  is not read for that library.
* T4, SANITY. Areg should be higher in DATP_like than in AT2, which trial C3
  established across four libraries. If it is not, the state calls are wrong
  and nothing else in this trial is read.

* THE READING, fixed in advance:
  - T1 met: transcriptional downregulation is sufficient to explain the escape,
    and the shedding account is not needed to explain trial C1d;
  - T1 not met: the transcript does not explain the escape, so any surface
    dimming would have to be post-transcriptional. That is consistent with the
    shedding account and is NOT evidence for it, because this trial cannot see
    protein;
  - T4 failed: nothing is read;
  - T3 failed for a library: T1 is not read for that library.
* NOT TESTED: surface EpCAM protein, shedding activity, ADAM17 enzyme activity,
  and whether any of this causes the escape. Those need flow cytometry on
  lineage-labelled cells, which is a bench experiment.
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

OUT = HERE / "c13_epcam_transcript_in_the_transitional_state"
OUT.mkdir(exist_ok=True)

STATES = {
    "AT2": ["Sftpc", "Etv5", "Lamp3"],
    "DATP_like": ["Cldn4", "Itga2", "Ndrg1", "Sox9"],
    "AT1_like": ["Ager", "Hopx", "Clic5"],
    "Cd177_positive": ["Cd177", "Cd38", "Dlk1"],
    "cycling": ["Mki67", "Birc5"],
}
PRIMARY = "Epcam"
SHEDDASE = ["Adam17", "Rhbdf1", "Rhbdf2", "Timp3"]
SANITY = ["Areg", "Krt8", "Cldn4"]
GENES = [PRIMARY] + SHEDDASE + SANITY
MIN_CELLS = 50
DROP_POINTS = 0.10

LIBRARIES = [
    {"gsm": "GSM7890831", "library": "Expt1_4dRFPr1", "series": "GSE247505", "arm": "KrasG12D 4d", "primary": True},
    {"gsm": "GSM7890832", "library": "Expt1_4dRFPr2", "series": "GSE247505", "arm": "KrasG12D 4d", "primary": True},
    {"gsm": "GSM7890835", "library": "Expt1_2wRFPr1", "series": "GSE247505", "arm": "KrasG12D 2w", "primary": True},
    {"gsm": "GSM7890836", "library": "Expt1_2wRFPr2", "series": "GSE247505", "arm": "KrasG12D 2w", "primary": True},
    {"gsm": "GSM9447781", "library": "Expt3_Het_RFP", "series": "GSE316244", "arm": "Areg-flox/+", "primary": False},
    {"gsm": "GSM9447782", "library": "Expt3_Hom_RFP", "series": "GSE316244", "arm": "Areg-flox/flox", "primary": False},
]
ROOTS = {"GSE247505": RAW / "GSE247505" / "GSE247505_RAW",
         "GSE316244": RAW / "GSE316244" / "GSE316244_RAW"}

RULES = {
    "question": "does Epcam transcript explain why the transitional state escapes an EpCAM-negative gate",
    "competing_explanations": {
        "transcriptional": "the state expresses less Epcam, stains dimly and falls out of the gate",
        "post_transcriptional": "the message is present but the protein is shed; ADAM17 releases amphiregulin "
                                "(doi:10.1083/jcb.200307137) and also cleaves EpCAM with presenilin-2 "
                                "(doi:10.1038/ncb1824)",
    },
    "what_this_cannot_do": ["measure surface protein", "measure shedding", "measure ADAM17 enzyme activity",
                            "establish causation for the escape"],
    "why_the_sheddase_genes_are_exploratory_only":
        "ADAM17 is controlled after translation through iRhom trafficking, a PDI conformational switch and "
        "phosphorylation, so its mRNA is a poor proxy for working enzyme and a negative is near-uninformative",
    "libraries": LIBRARIES,
    "qc": "identical to trial C1",
    "states": STATES,
    "state_assignment": "score_genes then per-cell argmax; differs from C3, which clustered first, and the "
                        "difference is disclosed",
    "groups": "DATP_like against AT2 within a library, both needing at least " + str(MIN_CELLS) + " cells",
    "genes": GENES,
    "T1": f"PRIMARY: Epcam detection lower in DATP_like in EVERY evaluable library AND median drop at least "
          f"{DROP_POINTS:.0%} points",
    "T2": "exploratory: Adam17, Rhbdf1, Rhbdf2, Timp3, called only if the direction agrees in every library, "
          "and never reported without the post-translational caveat",
    "T3": "depth control: median genes per cell per group, and T1 recomputed on shallow and deep halves",
    "T4": "sanity: Areg higher in DATP_like than AT2, as trial C3 established; if not, nothing is read",
    "reading_fixed_in_advance": {
        "T1_met": "transcriptional downregulation suffices to explain the escape; the shedding account is not "
                  "needed to explain trial C1d",
        "T1_not_met": "the transcript does not explain the escape, so any surface dimming would have to be "
                      "post-transcriptional; consistent with the shedding account and NOT evidence for it",
        "T4_failed": "the state calls are wrong and nothing is read",
    },
    "unit": "the library; direction consistency across libraries, no P value on a group comparison",
}


def load(entry, rec):
    import anndata as ad
    import scanpy as sc

    root = ROOTS[entry["series"]]
    stem = entry["gsm"] + "_" + entry["library"] + "_"
    matrix = root / (stem + "matrix.mtx.gz")
    features = root / (stem + "features.tsv.gz")
    if not features.exists():
        features = root / (stem + "genes.tsv.gz")
    barcodes = root / (stem + "barcodes.tsv.gz")
    for path in (matrix, features, barcodes):
        if not path.exists():
            print("missing:", path.name)
            return None
        rec.add_input(path)
    X, var, bc = read_mtx_triplet(matrix, features, barcodes)
    keep = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
    if keep.sum() == 0:
        keep = np.ones(len(var), dtype=bool)
    X, var = X[:, keep], var.loc[keep].reset_index(drop=True)
    var.index = pd.Index(make_unique(var["gene_symbol"].astype(str).to_numpy()).astype(str))
    obs = pd.DataFrame(index=pd.Index([entry["library"] + "_" + b for b in bc]))
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
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    for state, genes in STATES.items():
        present = [g for g in genes if g in adata.var_names]
        if present:
            sc.tl.score_genes(adata, present, score_name="score_" + state, random_state=0)
        else:
            adata.obs["score_" + state] = np.nan
    columns = ["score_" + s for s in STATES]
    adata.obs["state"] = np.array(list(STATES))[adata.obs[columns].to_numpy().argmax(axis=1)]
    print(entry["library"], adata.n_obs, "cells;", dict(adata.obs["state"].value_counts()))
    return adata


def detection(adata, mask, gene):
    if gene not in adata.var_names:
        return None, None
    counts = np.asarray(adata[mask, gene].layers["counts"].todense()).ravel()
    values = np.asarray(adata[mask, gene].X.todense()).ravel()
    return float((counts > 0).mean()), float(values.mean())


def main():
    rec = RunRecord(OUT / "c13_run_record.json",
                    "C13 Epcam transcript in the transitional state", RULES)
    rows, depth_rows = [], []
    for entry in LIBRARIES:
        adata = load(entry, rec)
        if adata is None:
            continue
        state = adata.obs["state"].to_numpy()
        depth = adata.obs["n_genes_by_counts"].to_numpy().astype(float)
        median_depth = float(np.median(depth))
        groups = {"DATP_like": state == "DATP_like", "AT2": state == "AT2"}
        evaluable = all(m.sum() >= MIN_CELLS for m in groups.values())
        for name, mask in groups.items():
            row = {"series": entry["series"], "library": entry["library"], "arm": entry["arm"],
                   "primary_set": entry["primary"], "state": name, "n_cells": int(mask.sum()),
                   "evaluable": evaluable,
                   "median_genes_per_cell": round(float(np.median(depth[mask])), 1) if mask.sum() else None}
            for gene in GENES:
                det, mean = detection(adata, mask, gene)
                row["det_" + gene] = round(det, 4) if det is not None else None
                row["mean_" + gene] = round(mean, 4) if mean is not None else None
            rows.append(row)
        if evaluable:
            for half, hmask in (("shallow", depth <= median_depth), ("deep", depth > median_depth)):
                entry_row = {"library": entry["library"], "half": half}
                for name, mask in groups.items():
                    both = mask & hmask
                    det, _ = detection(adata, both, PRIMARY)
                    entry_row["n_" + name] = int(both.sum())
                    entry_row["det_" + name] = round(det, 4) if det is not None and both.sum() else None
                depth_rows.append(entry_row)
        del adata
        gc.collect()

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "c13_by_library_and_state.csv", index=False)
    rec.add_output(OUT / "c13_by_library_and_state.csv")
    depth = pd.DataFrame(depth_rows)
    if len(depth):
        depth["drop"] = depth["det_AT2"] - depth["det_DATP_like"]
        depth.to_csv(OUT / "c13_depth_split.csv", index=False)
        rec.add_output(OUT / "c13_depth_split.csv")

    ok = table[table["evaluable"]]
    wide = ok.pivot_table(index=["series", "library", "arm"], columns="state",
                          values=["det_" + g for g in GENES] + ["median_genes_per_cell"])
    contrasts = []
    for gene in GENES:
        column = "det_" + gene
        if column not in wide:
            continue
        sub = wide[column].dropna()
        if sub.empty:
            continue
        for index, row in sub.iterrows():
            contrasts.append({"series": index[0], "library": index[1], "arm": index[2], "gene": gene,
                              "det_DATP_like": round(float(row["DATP_like"]), 4),
                              "det_AT2": round(float(row["AT2"]), 4),
                              "difference": round(float(row["DATP_like"] - row["AT2"]), 4)})
    frame = pd.DataFrame(contrasts)
    frame.to_csv(OUT / "c13_contrasts.csv", index=False)
    rec.add_output(OUT / "c13_contrasts.csv")

    def consistent(gene, direction):
        sub = frame[frame["gene"] == gene]
        if sub.empty:
            return None, sub
        signs = np.sign(sub["difference"].to_numpy())
        return bool((signs == direction).all()), sub

    areg_ok, areg = consistent("Areg", 1)
    rec.set("T4_sanity", {"areg_higher_in_DATP_like_everywhere": areg_ok,
                          "per_library": areg["difference"].tolist() if len(areg) else []})

    epcam_lower, epcam = consistent(PRIMARY, -1)
    median_drop = float(-epcam["difference"].median()) if len(epcam) else None
    t1_met = bool(epcam_lower and median_drop is not None and median_drop >= DROP_POINTS)
    rec.set("T1", {"lower_in_every_library": epcam_lower,
                   "median_drop_points": round(median_drop, 4) if median_drop is not None else None,
                   "required": DROP_POINTS, "met": t1_met,
                   "per_library": epcam[["library", "det_DATP_like", "det_AT2", "difference"]].to_dict("records")
                   if len(epcam) else []})

    shed = {}
    for gene, predicted in (("Adam17", 1), ("Rhbdf1", 1), ("Rhbdf2", 1), ("Timp3", -1)):
        agrees, sub = consistent(gene, predicted)
        shed[gene] = {"predicted_direction": "up" if predicted > 0 else "down",
                      "agrees_in_every_library": agrees,
                      "median_difference": round(float(sub["difference"].median()), 4) if len(sub) else None}
    rec.set("T2_exploratory_sheddase", shed)
    rec.set("T2_caveat", RULES["why_the_sheddase_genes_are_exploratory_only"])

    depth_ok = True
    if len(depth):
        for library, part in depth.groupby("library"):
            if part["drop"].notna().all() and not (np.sign(part["drop"]) == np.sign(part["drop"].iloc[0])).all():
                depth_ok = False
    rec.set("T3_depth_direction_survives_both_halves", depth_ok)

    if areg_ok is False:
        reading = RULES["reading_fixed_in_advance"]["T4_failed"]
    elif t1_met:
        reading = RULES["reading_fixed_in_advance"]["T1_met"]
    else:
        reading = RULES["reading_fixed_in_advance"]["T1_not_met"]
    if not depth_ok:
        reading += " The depth control failed in at least one library, whose T1 is not read."
    rec.set("reading", reading)

    lines = ["# Trial C13: Epcam transcript in the transitional state", "",
             "T4 sanity, Areg higher in DATP-like everywhere: " + str(areg_ok) + ".",
             "T1, Epcam lower in DATP-like in every library: " + str(epcam_lower)
             + "; median drop " + (f"{median_drop:.4f}" if median_drop is not None else "n/a")
             + " against a required " + str(DROP_POINTS) + ". Met: " + str(t1_met) + ".", "",
             "Reading: " + reading, "",
             "**Not measured:** surface EpCAM protein, shedding, ADAM17 activity. Those need flow cytometry.",
             "", "## Epcam and context, per library", "",
             df_to_markdown(frame[frame["gene"].isin([PRIMARY, "Areg"])], index=False), "",
             "## Exploratory sheddase context (mRNA is a poor proxy for ADAM17 activity)", "",
             df_to_markdown(frame[frame["gene"].isin(SHEDDASE)], index=False), ""]
    (OUT / "c13_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c13_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
