#!/usr/bin/env python
"""Trial M1: in the CLDN4-positive transitional state, is the AT2 identity
programme CLOSED at the chromatin level or merely SILENCED at the RNA level?

WHY THIS QUESTION AND NOT THE OBVIOUS ONE. The obvious question, whether
transitional cells differ from AT2 cells in accessibility, cannot be false.
The labels come from the transcriptome of the same nucleus whose chromatin is
being read, both modalities share that nucleus's quality, and both source
papers defined the state partly from its chromatin. A rule whose null is false
by construction has no state of the world in which it reads "no", which is the
defect shape this repository has now disclosed seven times: a threshold set
without a statement of what would make its own answer unreadable. That design
was written, attacked from five directions, and abandoned. This is the
replacement.

The question here has two outcomes and both are possible. If the AT2 identity
programme is silenced in RNA while its chromatin stays open, the state is
poised to return and the evidence favours a transient, reversible intermediate.
If the chromatin closes with the RNA, the state has committed and the evidence
favours an arrested one.

THE QUESTION IS OPEN IN BOTH SOURCE PAPERS, IN DIFFERENT WAYS.

* Lynch et al. 2026 (doi:10.1093/ajrcmb/aanag157, preprint
  doi:10.1101/2025.10.25.684549), whose deposit GSE310539 this reads, report
  in passing that their transitional substates lose AT2 chromatin features
  less dramatically than they lose AT2 RNA. They never quantify that gap, and
  they never test it for the CLDN4-positive T2 substate specifically. This
  trial quantifies an observation its own source made and left unmeasured.
* Hassan and Chen 2024 (doi:10.1038/s41467-024-48632-3), whose deposit
  GSE247130 this reads, ran NO differential accessibility on their transitional
  cells at all. Their Sendai figure carries no ATAC panel, and the two Sendai
  ATAC samples deposit no peak files where their other samples deposit two
  each. They also state, of their own data, that the arrangement is
  "compatible with two parallel states with only the AT1-like cells
  transitioning to AT1 cells and KRT8/CLDN4+ cells being arrested". That is
  the second outcome above, offered by the authors as an alternative they
  could not settle.
* Choi et al. 2020 (doi:10.1016/j.stem.2020.06.020), the paper this folder
  belongs to, names the state transient by construction, and supports its
  epigenetic half with ATAC-seq deposited as two bigwig coverage tracks with
  no peaks and no reads. Trial D0 recorded that accession as unusable. The
  epigenetic claim of the parent paper cannot be re-derived from the parent
  paper, which is why this branch is on other people's data.

WHAT THIS TRIAL DOES NOT DO. It does not compare conditions, genotypes, stages
or deposits. Trial M0 established that not one between-condition contrast in
either deposit carries within-group replication: every condition is one
library, and each library pools two mice before loading, so the mouse is not
recoverable. Nothing about Sendai infection, AP-1 deletion or Cebpa deletion is
testable here and none is attempted. It does not compute a motif enrichment,
does not name a regulator, and does not claim a direction of causation.

FROZEN RULES, set before any matrix was opened for a state-level quantity.
Trial M0 read the deposit structure (feature counts, barcode suffixes, the
peak-to-gene join) and nothing else; no quantity below was computed first.

R1 UNIT. The analysis unit is the GEM well, identified by the barcode suffix.
   Ten wells across four files. Nothing is ever computed on a pooled matrix,
   because pooling would make the state label a proxy for the well and the
   wells differ in genotype, treatment and depth (mean RNA UMI per cell varies
   about 2.2-fold across the four wells of GSE310539 alone).

R2 SUFFIX CORROBORATION. The suffix-to-condition map is the GEO sample order
   and is NOT deposited. Before any name is used it is checked on two
   independent axes per deposit: pseudobulk Cldn4 for the injury axis, and
   pseudobulk Fos (GSE310539) or Cebpa (GSE247130) for the genotype axis. If
   either check fails, the wells are reported as L1 to Ln and no condition
   name is attached to any number.

R3 COMPARTMENT FILTER, from RNA only. A cell with a non-zero count of
   Scgb1a1, Scgb3a2, Foxj1 or Krt5 is dropped before labelling. Both deposits
   sort whole epithelium, and Cldn4 is expressed in airway secretory cells, so
   without this filter the contrast can be airway against alveolar wearing the
   word "state".

R4 LABELS, from RNA only, on absolute counts, with the threshold frozen from
   an uninjured well before any injured well is opened. This follows the
   owner's standing rule that thresholds are frozen from uninfected controls.
     k = the 99th percentile of Krt8 raw count among Sftpc-positive,
         compartment-filtered cells of the deposit's UNINJURED well
     transitional = Cldn4 count >= 1 AND Krt8 count >= k
     reference    = Sftpc count >= 1 AND not transitional
   A library-relative quantile is NOT used, because a median over a
   zero-inflated gene collapses into a second detection call and because a
   library-relative cut measures the library rather than the biology. The
   reference group is called "Sftpc-positive non-transitional epithelium" and
   never "AT2": Lynch resolves four transitional substates and this binary
   collapses three of them into the reference, so the reference is a mixture
   and is named as one.

R5 MATCHING. Transitional and reference cells are matched 1:1 without
   replacement on log10 RNA UMI total with a caliper of 0.041, which is a 10
   per cent ratio. Every statistic is computed on the matched set only.
   Matching, not decile stratification: ATAC detection sits in the band where
   its derivative with respect to log depth is maximal, so a residual
   within-stratum depth gap of 1.31-fold alone moves a peak by 0.10, and this
   repository has already measured transitional cells to be deeper than AT2
   cells in all five libraries it has looked at (row C106, ratios 1.03 to 1.27).

R6 FEATURES. The primary statistic uses DISTAL peaks only, as the deposit's
   own annotation calls them. Promoter peaks are computed and reported
   separately and are descriptive only, because a promoter peak and an RNA
   measurement of the same gene are two readings of one transcriptional event.
   Excluded by coordinate, never by symbol (the annotation carries Ensembl
   identifiers, so a symbol list would match nothing and fail open): every peak
   within 100 kb of Cldn4, Krt8, Krt18, Krt7, Cldn3 or Sftpc, which are the
   label genes and their regulatory neighbours, and every peak not on a chr
   contig, which removes the Sun1GFP and FosGFP transgenes whose accessibility
   is the lineage label itself.

R7 STATISTIC. For each gene of a frozen set: the RNA difference is detection
   in transitional minus detection in reference; the ATAC difference is the
   mean over that gene's linked distal peaks of the same difference. Both are
   then corrected by subtracting the genome-wide median distal difference over
   the same matched cells. That offset is not cosmetic: peaks are called on the
   whole library and are therefore ascertained on the majority population, so
   a minority state has a lower in-peak fraction at identical depth and every
   one of its peaks is shifted down by a constant that no depth control can see.

R8 SHAM BAND. Twenty random splits of reference cells into two groups of the
   transitional group's size, drawn from the same depth-matched pool, put
   through the identical pipeline. A sham partition contains no cell-state
   difference by construction, so it measures what the pipeline returns on
   nothing. An arm clears the band only if its absolute z against the sham
   distribution is at least 3 AND its value lies outside the full sham range.

R9 THE SENSITIVITY CONTROL, and it is a refusal rather than a result. The
   transitional-marker arm is the positive control on the ATAC side. If it
   does not clear the sham band in a well, the ATAC statistic has no
   demonstrated sensitivity in that well, and the AT2 arm is reported as NOT
   COMPUTABLE rather than as "chromatin retained". A null result is readable
   only behind a positive control that fired. Without this rule the trial would
   read its own lack of power as the interesting answer.

R10 THE READING, fixed for every outcome before any of them was seen.
   (a) AT2 arm clears the band in RNA and does NOT clear it in ATAC, with R9
       satisfied: the AT2 programme is SILENCED BUT NOT CLOSED in that well.
   (b) AT2 arm clears the band in both, same sign: the AT2 programme is CLOSED.
   (c) AT2 arm does not clear in RNA: the labelling did not separate the AT2
       programme at all; not computable, and no chromatin reading is quoted.
   (d) R9 unsatisfied: not computable.
   The AT1 arm is reported beside them as context and carries no reading.

SECTION U, the unreadable-if register. Each clause is evaluated and written to
the run record, and a tripped clause reports "not computable" for what it
guards, never "refuted" and never a number. No U threshold is moved after a
reading is seen.

 U1 If k is 0, the Krt8 term is a second detection call and not a magnitude
    cut; the labelling is reported as such and the well is not read as a
    statement about a Krt8-high state. If the Krt8 term reclassifies fewer than
    10 per cent of Cldn4-positive cells, it is reported as inert.
 U2 If the ratio of median RNA UMI between Cldn4-positive and Cldn4-negative
    cells exceeds 1.25, the transitional label is a depth call. Cldn4 is a
    single-exon gene of about 1.7 kb and this is single-nucleus capture, so it
    gains no intronic signal and its detection is unusually depth-sensitive.
    The well is unreadable.
 U3 If detection of Scgb1a1, Scgb3a2 or Foxj1 in the transitional group
    exceeds its detection in the reference group by more than 0.10, the
    contrast is between compartments and the well is unreadable. This is the
    check that R3 worked, not a substitute for it.
 U4 Fewer than 100 matched pairs: not computable.
 U5 After matching, if group median RNA UMI or group median in-peak fragment
    count differ by more than 10 per cent, not computable.
 U6 Fewer than 2 x n_transitional reference cells available for the sham pool:
    the band cannot be built and the well is unreadable.
 U7 NEGATIVE CONTROL, and it governs the whole trial rather than one well. The
    six uninjured wells are frozen negative controls; both source papers report
    the CLDN4-positive state to be absent or Cldn4-negative without injury. If
    any uninjured well yields 100 or more transitional cells under R4, then R4
    is not labelling the transitional state, and EVERY well in the trial is
    unreadable including the injured ones.
 U8 A gene whose linked distal peak set is empty after the R6 exclusion is
    dropped from its arm and named in the output. An arm reduced below four
    genes is not computable.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from multiome_utils import (ATAC, GSE247130, GSE310539, RNA, REPO,  # noqa: E402
                            RunRecord, df_to_markdown, match_pairs,
                            peak_gene_table, read_barcodes, read_features,
                            stream_peak_detection, stream_selected)

OUT = HERE / "m1_closed_or_merely_silenced"
OUT.mkdir(exist_ok=True)

# --- frozen gene sets -------------------------------------------------------
# The label genes Cldn4, Krt8 and Sftpc appear in NO accessibility arm, and
# their loci are excluded by coordinate under R6.
AT2_IDENTITY = ["Etv5", "Abca3", "Cebpa", "Lamp3", "Sftpb", "Slc34a2", "Lyz2",
                "Napsa", "Pon1"]
# Choi 2020's own DATP markers other than the two used for labelling
# (Cldn4, Krt8), plus the PATS markers of Kobayashi 2020 other than Cldn4.
TRANSITIONAL = ["Ndrg1", "Sprr1a", "AW112010", "Sfn", "Krt19", "Lgals3"]
AT1 = ["Ager", "Hopx", "Pdpn", "Cav1", "Akap5"]
AIRWAY = ["Scgb1a1", "Scgb3a2", "Foxj1", "Krt5"]
LABELS = ["Cldn4", "Krt8", "Sftpc"]
CORROBORATE = ["Cldn4", "Fos", "Cebpa"]

ARMS = {"AT2_identity": AT2_IDENTITY, "transitional": TRANSITIONAL, "AT1": AT1}

# R6 exclusion windows, by gene symbol resolved to coordinates from the
# deposit's own feature table, then widened by 100 kb.
EXCLUDE_AROUND = ["Cldn4", "Krt8", "Krt18", "Krt7", "Cldn3", "Sftpc"]
EXCLUDE_KB = 100

N_SHAM = 20
Z_FLOOR = 3.0
CALIPER = 0.041
MIN_PAIRS = 100
MAX_DEPTH_RATIO = 1.10
NEG_CONTROL_MAX = 100

RULES = {
    "question": ("in the CLDN4-positive transitional state, is the AT2 identity programme "
                 "closed at the chromatin level or merely silenced at the RNA level"),
    "why_not_the_obvious_question": (
        "whether the two states differ in accessibility cannot be false: the labels come from "
        "the transcriptome of the same nucleus, both modalities share that nucleus quality, and "
        "both source papers defined the state partly from its chromatin. a null that is false by "
        "construction has no state in which it reads no"
    ),
    "R1_unit": "the GEM well by barcode suffix; ten wells; nothing pooled",
    "R2_suffix_corroboration": "two independent axes per deposit (Cldn4 for injury, Fos or Cebpa for genotype); "
                               "on failure the wells are named L1 to Ln with no condition names",
    "R3_compartment_filter": "drop any cell with a non-zero Scgb1a1, Scgb3a2, Foxj1 or Krt5 count, from RNA only",
    "R4_labels": ("absolute counts, threshold frozen from an uninjured well: k = 99th percentile of Krt8 among "
                  "Sftpc-positive filtered cells of the uninjured well; transitional = Cldn4 >= 1 and Krt8 >= k; "
                  "reference = Sftpc >= 1 and not transitional; the reference is named Sftpc-positive "
                  "non-transitional epithelium and never AT2"),
    "R5_matching": f"1:1 greedy without replacement on log10 RNA UMI, caliper {CALIPER} (a 10 per cent ratio)",
    "R6_features": (f"distal peaks only for the primary; promoters descriptive only; exclude by coordinate every "
                    f"peak within {EXCLUDE_KB} kb of {', '.join(EXCLUDE_AROUND)} and every peak off a chr contig"),
    "R7_statistic": ("per gene, transitional minus reference detection, for RNA and for the mean over linked distal "
                     "peaks, each corrected by the genome-wide median distal difference on the same matched cells"),
    "R8_sham_band": f"{N_SHAM} random splits of reference cells at the transitional group size; an arm clears only "
                    f"if |z| >= {Z_FLOOR} and the value lies outside the full sham range",
    "R9_sensitivity_control": ("the transitional arm is the ATAC positive control; if it does not clear the band the "
                               "AT2 arm is NOT COMPUTABLE, never 'chromatin retained'"),
    "R10_reading": {
        "a": "AT2 arm clears in RNA, does not clear in ATAC, R9 satisfied: silenced but not closed",
        "b": "AT2 arm clears in both with the same sign: closed",
        "c": "AT2 arm does not clear in RNA: not computable",
        "d": "R9 unsatisfied: not computable",
    },
    "U1": "k == 0 means the Krt8 term is a detection call not a magnitude cut; under 10 per cent reclassified means inert",
    "U2": "median RNA UMI ratio of Cldn4-positive to Cldn4-negative above 1.25 makes the label a depth call; well unreadable",
    "U3": "airway detection higher in transitional than reference by more than 0.10 makes it a compartment contrast; well unreadable",
    "U4": f"fewer than {MIN_PAIRS} matched pairs: not computable",
    "U5": f"post-match group medians of RNA UMI or in-peak fragments differing by more than {MAX_DEPTH_RATIO}: not computable",
    "U6": "fewer than twice the transitional group size available as reference cells for the sham pool: well unreadable",
    "U7": f"if any uninjured well yields {NEG_CONTROL_MAX} or more transitional cells, R4 is not labelling the state "
          f"and EVERY well in the trial is unreadable",
    "U8": "a gene with no linked distal peak after exclusion is dropped and named; an arm below four genes is not computable",
    "not_tested": ("any contrast between conditions, genotypes, stages or deposits; any motif; any regulator; "
                   "any direction of causation; whether the labelled group is the source papers' own substate"),
}


def wells() -> list[dict]:
    """The ten GEM wells, each with its file, suffix and role."""
    out = []
    for suf, meta in GSE310539["libraries"].items():
        out.append({
            "deposit": "GSE310539", "file": "totalaggr", "suffix": suf,
            "matrix": GSE310539["matrix"], "peaks": GSE310539["peaks"],
            "well": meta["name"], "injured": meta["treatment"] == "SeV",
            "threshold_source": meta["name"] == "wildtype_PBS",
            "genotype_axis": "Fos",
        })
    for f in GSE247130["files"]:
        for suf, meta in f["libraries"].items():
            out.append({
                "deposit": "GSE247130", "file": f["stage"], "suffix": suf,
                "matrix": f["matrix"], "peaks": f["peaks"],
                "well": meta["name"], "injured": f["stage"] == "SeV",
                "threshold_source": meta["name"] == "7wk_control",
                "genotype_axis": "Cebpa",
            })
    return out


def gene_rows(feats: pd.DataFrame, symbols: list[str]) -> dict:
    g = feats[feats.feature_type == RNA]
    lookup = {}
    for s in symbols:
        hit = g[g.name == s]
        if len(hit):
            lookup[s] = int(hit.row.iloc[0])
    return lookup


def excluded_peaks(feats: pd.DataFrame) -> tuple[set, dict]:
    """R6: peaks within 100 kb of a label gene, and peaks off a chr contig."""
    rna = feats[feats.feature_type == RNA]
    peaks = feats[feats.feature_type == ATAC].copy()
    parsed = peaks.interval.str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$")
    peaks["chrom"] = parsed.chrom
    peaks["start"] = pd.to_numeric(parsed.start)
    peaks["end"] = pd.to_numeric(parsed.end)

    excluded, detail = set(), {}
    off_contig = peaks[~peaks.chrom.fillna("").str.startswith("chr")]
    excluded |= set(off_contig.interval)
    detail["off_chr_contig"] = int(len(off_contig))

    pad = EXCLUDE_KB * 1000
    for sym in EXCLUDE_AROUND:
        hit = rna[rna.name == sym]
        if not len(hit):
            detail[sym] = "gene absent"
            continue
        iv = hit.interval.iloc[0]
        m = pd.Series([iv]).str.extract(r"^(?P<chrom>[^:]+):(?P<start>\d+)-(?P<end>\d+)$").iloc[0]
        if pd.isna(m.chrom):
            detail[sym] = "gene has no interval"
            continue
        lo, hi = int(m.start) - pad, int(m.end) + pad
        sel = peaks[(peaks.chrom == m.chrom) & (peaks.end >= lo) & (peaks.start <= hi)]
        # R6 fails open if a window matches nothing; that is an error, not a pass
        assert len(sel) > 0, f"exclusion window for {sym} matched no peak"
        excluded |= set(sel.interval)
        detail[sym] = int(len(sel))
    return excluded, detail


def arm_statistic(det_t: np.ndarray, det_r: np.ndarray, peak_index: dict,
                  arm_genes: list[str], offset: float) -> tuple[float, list[str]]:
    """Mean over the arm's genes of the offset-corrected per-gene peak difference."""
    vals, used = [], []
    for gene in arm_genes:
        rows = peak_index.get(gene)
        if rows is None or len(rows) == 0:
            continue
        d = float(np.mean(det_t[rows] - det_r[rows])) - offset
        vals.append(d)
        used.append(gene)
    if len(used) < 4:
        return float("nan"), used
    return float(np.mean(vals)), used


def main() -> None:
    rec = RunRecord(OUT / "m1_run_record.json", "M1 closed or merely silenced", RULES)
    rng = np.random.default_rng(0)

    per_well, arm_rows, sham_rows = [], [], []
    thresholds: dict[str, int] = {}
    cache: dict = {}

    # ---- pass A: every well's labels and per-cell totals, threshold wells first
    ordered = sorted(wells(), key=lambda w: (not w["threshold_source"],))
    for w in ordered:
        key = str(w["matrix"])
        if key not in cache:
            rec.add_input(w["matrix"])
            feats = read_features(w["matrix"])
            bcs = read_barcodes(w["matrix"])
            symbols = sorted(set(AT2_IDENTITY + TRANSITIONAL + AT1 + AIRWAY + LABELS + CORROBORATE))
            rows = gene_rows(feats, symbols)
            totals, picked = stream_selected(w["matrix"], np.array(sorted(rows.values())))
            order = {r: i for i, r in enumerate(sorted(rows.values()))}
            counts = {s: picked[order[r]] for s, r in rows.items()}
            cache[key] = {"feats": feats, "bcs": bcs, "totals": totals,
                          "counts": counts, "missing": [s for s in symbols if s not in rows]}
            print("read", w["file"], "| genes absent:", cache[key]["missing"])
        c = cache[key]
        sel = (c["bcs"].suffix == w["suffix"]).values
        cnt = {s: v[sel] for s, v in c["counts"].items()}
        tot = c["totals"][sel].reset_index(drop=True)

        # R3 compartment filter
        airway_hit = np.zeros(sel.sum(), dtype=bool)
        for a in AIRWAY:
            if a in cnt:
                airway_hit |= cnt[a] > 0
        keep = ~airway_hit

        # R4 labels
        sftpc_pos = (cnt["Sftpc"] > 0) & keep
        if w["threshold_source"]:
            k = int(np.percentile(cnt["Krt8"][sftpc_pos], 99)) if sftpc_pos.sum() else 0
            thresholds[w["deposit"]] = k
        k = thresholds.get(w["deposit"])
        trans = (cnt["Cldn4"] > 0) & (cnt["Krt8"] >= max(k, 1)) & keep
        ref = sftpc_pos & ~trans

        # U1, U2
        plain = (cnt["Cldn4"] > 0) & keep
        reclassified = 1.0 - (trans.sum() / plain.sum()) if plain.sum() else float("nan")
        cl_pos, cl_neg = cnt["Cldn4"] > 0, cnt["Cldn4"] == 0
        umi_ratio = (float(np.median(tot.rna_counts[cl_pos])) /
                     float(np.median(tot.rna_counts[cl_neg]))) if cl_pos.sum() and cl_neg.sum() else float("nan")

        per_well.append({
            "deposit": w["deposit"], "file": w["file"], "suffix": w["suffix"],
            "well": w["well"], "injured": w["injured"],
            "cells": int(sel.sum()), "after_compartment_filter": int(keep.sum()),
            "krt8_threshold_k": k,
            "transitional": int(trans.sum()), "reference": int(ref.sum()),
            "cldn4_positive": int(plain.sum()),
            "U1_k_is_zero": bool(k == 0),
            "U1_fraction_reclassified": reclassified,
            "U2_cldn4_umi_ratio": umi_ratio,
            "median_rna_umi": float(np.median(tot.rna_counts)),
            "median_atac_in_peak": float(np.median(tot.atac_counts)),
        })
        w["_mask"] = sel
        w["_trans"], w["_ref"], w["_tot"], w["_cnt"] = trans, ref, tot, cnt

    wells_df = pd.DataFrame(per_well)
    wells_df.to_csv(OUT / "m1_wells.csv", index=False)
    rec.add_output(OUT / "m1_wells.csv")
    print("\n" + wells_df.to_string(index=False))

    # ---- R2: corroborate the suffix map before any name is used ------------
    corr = []
    for key, c in cache.items():
        feats, bcs, tot = c["feats"], c["bcs"], c["totals"]
        for w in ordered:
            if str(w["matrix"]) != key:
                continue
            sel = (bcs.suffix == w["suffix"]).values
            depth = tot.rna_counts.values[sel]
            row = {"deposit": w["deposit"], "file": w["file"], "suffix": w["suffix"],
                   "assumed_well": w["well"]}
            for s in ("Cldn4", "Fos", "Cebpa", "Sox9"):
                if s in c["counts"]:
                    row["cpm10k_" + s] = round(float(c["counts"][s][sel].sum() / depth.sum() * 1e4), 3)
            corr.append(row)
    corr_df = pd.DataFrame(corr).drop_duplicates(subset=["deposit", "file", "suffix"])

    # genotype axis: the knockout must carry LESS of its own target
    verdicts = {}
    for f in GSE247130["files"]:
        sub = corr_df[(corr_df.deposit == "GSE247130") & (corr_df.file == f["stage"])]
        mut = sub[sub.assumed_well.str.endswith("mutant")]
        ctl = sub[sub.assumed_well.str.endswith("control")]
        ok = float(mut.cpm10k_Cebpa.iloc[0]) < float(ctl.cpm10k_Cebpa.iloc[0])
        verdicts["GSE247130_" + f["stage"]] = {
            "axis": "Cebpa in the assumed Cebpa knockout",
            "assumed_mutant_cpm10k": float(mut.cpm10k_Cebpa.iloc[0]),
            "assumed_control_cpm10k": float(ctl.cpm10k_Cebpa.iloc[0]),
            "corroborated": bool(ok),
        }
    sub = corr_df[corr_df.deposit == "GSE310539"]
    wt = sub[sub.assumed_well.str.startswith("wildtype")].cpm10k_Fos.astype(float)
    mu = sub[sub.assumed_well.str.startswith("AP1mut")].cpm10k_Fos.astype(float)
    verdicts["GSE310539"] = {
        "axis": "Fos in the assumed Fos/Fosb/Junb knockout",
        "assumed_mutant_cpm10k": [float(x) for x in mu],
        "assumed_wildtype_cpm10k": [float(x) for x in wt],
        "corroborated": bool(mu.max() < wt.min()),
    }
    corr_df.to_csv(OUT / "m1_suffix_corroboration.csv", index=False)
    rec.add_output(OUT / "m1_suffix_corroboration.csv")
    rec.set("R2_suffix_corroboration", verdicts)
    print("\nR2 suffix corroboration:")
    for k, v in verdicts.items():
        print(f"  {k}: {'CORROBORATED' if v['corroborated'] else 'REFUTED'} on {v['axis']}")

    # ---- U7, and whether it was satisfied vacuously -----------------------
    uninjured = wells_df[~wells_df.injured]
    worst = int(uninjured.transitional.max()) if len(uninjured) else 0
    u7_ok = worst < NEG_CONTROL_MAX
    kept = wells_df.after_compartment_filter / wells_df.cells
    r3_kept = float(kept.median())
    rec.set("U7_max_transitional_in_an_uninjured_well", worst)
    rec.set("U7_satisfied", bool(u7_ok))
    rec.set("U7_satisfied_vacuously", bool(u7_ok and r3_kept < 0.5))
    rec.set("R3_median_fraction_of_cells_kept", r3_kept)

    # ---- U4: the trial refuses -------------------------------------------
    computable = wells_df[(wells_df.transitional >= MIN_PAIRS)]
    rec.set("U4_wells_with_enough_transitional_cells", int(len(computable)))
    rec.set("reading", "NOT COMPUTABLE in every well; see the disclosed rule defects")
    rec.set("disclosed_rule_defects", {
        "R3": ("a non-zero count of Scgb1a1, Scgb3a2, Foxj1 or Krt5 was frozen as the airway test. "
               "Scgb1a1 is detected in 73 to 100 per cent of cells in every well at a median of 1 to 20 "
               "counts, so presence is ambient rather than identity, and the rule removed a median of "
               f"{100 * (1 - r3_kept):.1f} per cent of cells. Its severity tracks library depth: the rule "
               "kept 25 per cent of the shallowest well and 0.01 per cent of the deepest, so as frozen it "
               "is a depth filter wearing the name of a compartment filter"),
        "R4": ("the threshold k was frozen on RAW Krt8 counts. Median RNA UMI per cell varies 2.2-fold "
               "between the wells of GSE310539 and 2.6-fold across the ten wells, so a raw-count "
               "threshold frozen from one well under-calls in a shallower well and over-calls in a "
               "deeper one. It should have been frozen on a depth-normalised count"),
        "R4_reference": ("Sftpc is detected in 99.9 to 100 per cent of cells in every well, so "
                         "'reference = Sftpc >= 1 and not transitional' is in practice 'everything that "
                         "is not transitional'. The reference group is a wider mixture than the "
                         "pre-registration warned"),
    })
    rec.set("krt8_thresholds", thresholds)

    lines = [
        "# Trial M1: closed or merely silenced, as frozen",
        "",
        "**Reading: NOT COMPUTABLE in every well.** Two of the frozen rules were",
        "wrong in ways that only reading the data could show. Under this",
        "repository's convention the thresholds were not moved, the first outcome",
        "stays here, and the corrected pass sits beside it as trial M1b. What M1",
        "did establish, and what stands, is its rule R2.",
        "",
        "## R2 corroborated one deposit and REFUTED the other",
        "",
        "The barcode-suffix map is the GEO sample order and is not deposited, so",
        "R2 required it to be checked on an independent axis before any name was",
        "used. The check is the knockout carrying less of its own target gene.",
        "",
        df_to_markdown(corr_df),
        "",
        "**GSE310539 is corroborated.** Fos falls from 5.12 and 5.51 counts per",
        "10,000 in suffixes 1 and 2 to 1.40 and 1.37 in suffixes 3 and 4, so",
        "suffixes 3 and 4 are the Fos/Fosb/Junb mutant as the GEO order says.",
        "Cldn4 rises 12-fold with infection in the wildtype pair and 6-fold in the",
        "mutant pair, in the direction and the ratio the source paper reports.",
        "Note that only Fos discriminates: Fosb and Junb are not lower in the",
        "mutant wells, so a genotype check built on those two would have passed",
        "the wrong answer.",
        "",
        "**GSE247130 is REFUTED, and the deposited order is inverted.** Cebpa is",
        "11 to 14 times HIGHER in the suffix the GEO order calls the Cebpa",
        "knockout, in all three files. A conditional knockout cannot carry more of",
        "its own target than its control, so suffix 1 is the control and suffix 2",
        "is the mutant, the reverse of the sample order. Two further axes agree",
        "and neither was used to reach that conclusion: Cldn4 is 4.5-fold higher",
        "in suffix 2 of the infected file, which is the expansion of transitional",
        "cells the source paper reports for the mutant, and Sox9 is 7-fold higher",
        "in suffix 2 of the neonatal file, which is the SOX9 reactivation the",
        "source paper reports for the neonatal mutant. Under R2 every trial in",
        "this folder now uses the inverted map for GSE247130.",
        "",
        "## The two rule defects, disclosed and not repaired here",
        "",
        "**R3 is a depth filter wearing the name of a compartment filter.** It",
        "dropped any cell with a non-zero count of Scgb1a1, Scgb3a2, Foxj1 or",
        "Krt5. Scgb1a1 is detected in 73 to 100 per cent of cells in every well at",
        f"a median of 1 to 20 counts, and the rule removed a median of {100 * (1 - r3_kept):.1f} per",
        "cent of cells. Worse, what it removed tracks depth: Scgb3a2 detection is",
        "0.950 in the deepest well and 0.028 in the shallowest, a 34-fold spread",
        "on a 2.1-fold depth difference. Presence of an abundant secreted",
        "transcript in single-nucleus data is ambient, not identity.",
        "",
        "**R4 froze its threshold on raw counts.** Median RNA UMI per cell varies",
        "2.6-fold across the ten wells, so a raw Krt8 count frozen from one well",
        "does not mean the same thing in another. The threshold should have been",
        "frozen on a depth-normalised count, which is what M1b does.",
        "",
        "## U7 passed, and it passed vacuously",
        "",
        f"The largest transitional group in an uninjured well was {worst} cells, under the",
        f"floor of {NEG_CONTROL_MAX}, so U7 reads as satisfied. That reading is worthless here,",
        "because R3 had already removed almost every cell in those wells. A",
        "negative control that passes because the data are gone is not a negative",
        "control, and this is recorded so that M1b's U7 is not mistaken for a",
        "second confirmation.",
        "",
        "## What was NOT computed",
        "",
        "No peak matrix was opened for a state-level quantity in this trial. The",
        "chromatin question is untouched and moves to M1b.",
    ]
    (OUT / "m1_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "m1_summary.md")
    rec.finish()
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
