#!/usr/bin/env python
"""Trial A1: are the Axin2 locus and the Il1r1 locus open in the same nuclei?

Choi 2020 proposed comparing Il1r1-positive and Axin2-positive AT2 cells and
nobody has done it. This folder's assessment establishes why: both populations
are defined by tamoxifen-inducible lineage reporters, no deposit carries both
readouts in the same cells, and both transcripts sit near five per cent
detection at about one molecule per positive cell, which is a Poisson coin flip
rather than a phenotype.

The one route that escapes that floor was measured before this trial was
written. At least one of Axin2's fifteen or sixteen linked peaks is detected in
22.1 and 39.4 per cent of cells against a transcript detected in 5.3 and 3.8 per
cent, a gain of 4.1 and 10.4 fold, with three controls ordering correctly in
both deposits (row C139). This trial uses that route.

WHAT IS ASKED, AND WHAT IS NOT. The question is whether the two loci are open
together in the same nucleus MORE OR LESS OFTEN THAN COMPARABLE PAIRS OF LOCI
ARE. It is not whether they are correlated, because at a fixed fragment budget
every pair of loci is positively correlated through the breadth of that
nucleus's open chromatin, so a correlation test would return yes for any pair
and could not be false. That is the defect shape this repository has now
disclosed nine times and it is refused here by construction: the null is a
distribution over MATCHED RANDOM GENE PAIRS, not independence.

THE SCOPE LIMIT, STATED FIRST BECAUSE IT GOVERNS EVERY SENTENCE BELOW. An
accessible Axin2 locus is not a Wnt-responsive cell. Accessibility reports that
a locus is in a configuration permitting expression; it is slower and more
permissive than transcription, and it is a weaker proxy for the lineage reporter
than the transcript would be if the transcript worked. Nothing this trial
produces may be phrased as a statement about Wnt-responsive or IL-1-responsive
cells, about AEPs, or about Choi's Il1r1-lineage population. The object is the
pair of loci, and it is named as such.

FROZEN RULES.

R1 UNIT. The GEM well, by barcode suffix, as trial M0 established and trial M1
   corrected for GSE247130, whose deposited order is inverted. Nothing pooled.
   Every well reported separately. Descriptive only: each well is one library
   pooling two mice.

R2 CELLS. Non-airway cells, by the within-cell ratio rule R3b of the
   neighbouring branch (airway if summed Scgb1a1, Scgb3a2, Foxj1 and Krt5
   exceeds summed Sftpc, Sftpb and Lamp3). The primary analysis set is the
   Sftpc-positive NON-TRANSITIONAL cells, which is this repository's AT2 proxy;
   the CLDN4-positive KRT8-positive group is reported beside it as secondary
   and descriptive, because Choi's claim is that DATPs arise from Il1r1-positive
   cells and the comparison is worth seeing even though it is not tested here.

R3 DEPTH. Every cell downsampled to exactly B in-peak fragments, B being the
   20th percentile over the analysis set of that well; cells below B are
   dropped and the number is reported. This is rule R14 of the neighbouring
   branch, and it exists because detection fraction is monotone and saturating
   in depth, so an undownsampled co-detection statistic measures which nuclei
   were sequenced deeply rather than which loci are open.

R4 LOCUS DETECTION. A locus is detected in a cell if AT LEAST ONE peak the
   vendor annotation links to its gene, promoter or distal, is detected after
   downsampling. Peak identity is never compared across files, because peaks are
   called per aggregate.

R5 THE STATISTIC. For a pair of loci (A, B) in one well:
       observed   = fraction of cells detecting both
       expected   = fraction detecting A times fraction detecting B
       ratio      = observed / expected
   The ratio, not the difference, because the marginals differ between wells and
   between pairs and a difference would not be comparable across either.

R6 THE NULL, AND IT IS OVER GENE PAIRS RATHER THAN OVER CELLS. 300 random gene
   pairs, each matched to (Axin2, Il1r1) on the number of linked peaks within a
   factor of 1.5, drawn from different chromosomes from each other so that
   physical linkage cannot inflate the null, and excluding every gene named in
   this trial. Their ratios give the distribution against which the observed
   ratio is read: z, and the percentile within the 300.

R7 TWO CONTROLS, AND BOTH MUST BEHAVE OR THE WELL IS NOT COMPUTABLE.
   POSITIVE, PHYSICAL: Krt8 and Krt18 sit 24 kilobases apart, so their linked
   peak sets overlap and the ratio MUST come out high. If it does not, the
   measurement is broken rather than the biology surprising.
   POSITIVE, BIOLOGICAL: Etv5 and Abca3, two AT2 identity genes on different
   chromosomes that this repository has already measured moving together. This
   is the control that shows the statistic can detect co-regulation that is not
   physical proximity. If it does not clear, the trial has no demonstrated
   sensitivity and the Axin2-Il1r1 reading is NOT COMPUTABLE, never "no
   relationship". This is rule R9 of the neighbouring branch, which that branch
   learned the hard way.

R8 A SECOND, MORE ROBUST FORM. Axin2 alone is one gene. The trial also reports
   a Wnt-target locus set (Axin2, Nkd1, Notum, Lef1, Wif1, Lgr5) against Il1r1,
   because a module is less hostage to one gene's peak annotation. Both are
   reported; neither is privileged; the single-gene form is the literal question
   and the module form is the one to believe if they disagree.

R9 THE READING, fixed before it was seen.
   (a) ratio above the null at |z| >= 3 and outside the full null range: the two
       loci are open together MORE often than comparable pairs.
   (b) below by the same criterion: LESS often.
   (c) neither: not distinguishable from comparable pairs. Reported as that, and
       not as evidence of independence, because the null's own spread is the
       limit of what the design can see and that spread is reported with it.
   (d) either control failing: NOT COMPUTABLE for that well.

SECTION U, the unreadable-if register.
 U1 Fewer than 500 cells in the analysis set after R3: not computable.
 U2 Either locus detected in fewer than 2 per cent or more than 98 per cent of
    cells: the ratio is unstable at the margins and the well is not computable
    for that pair.
 U3 Fewer than 200 usable matched pairs in the null: not computable.
 U4 The physical control Krt8-Krt18 not clearing: the measurement is broken and
    the whole well is unreadable.
 U5 The biological control Etv5-Abca3 not clearing: no demonstrated sensitivity,
    so every negative in that well is not computable rather than a null result.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
TRIALS = HERE.parent / "datp_epigenetics" / "trials"
sys.path.insert(0, str(TRIALS))
from multiome_utils import (ATAC, RNA, RunRecord, df_to_markdown,  # noqa: E402
                            peak_gene_table, read_barcodes, read_features,
                            stream_downsampled_detection, stream_selected)
from m1_closed_or_merely_silenced import AIRWAY, gene_rows  # noqa: E402
from m1b_corrected_pass import SEED, wells_corrected  # noqa: E402
from m1d_mature_wells_only import NEONATAL  # noqa: E402

OUT = HERE / "a1_locus_co_accessibility"
OUT.mkdir(exist_ok=True)

WNT_SET = ["Axin2", "Nkd1", "Notum", "Lef1", "Wif1", "Lgr5"]
PAIRS = [
    ("Axin2", "Il1r1", "primary"),
    ("WNT_MODULE", "Il1r1", "primary, module form"),
    ("Krt8", "Krt18", "control, physical (24 kb apart)"),
    ("Etv5", "Abca3", "control, biological (AT2 identity, different chromosomes)"),
]
NAMED = set(WNT_SET + ["Il1r1", "Krt8", "Krt18", "Etv5", "Abca3", "Sftpc",
                       "Sftpb", "Lamp3", "Cldn4"] + AIRWAY)
N_NULL = 300
Z_FLOOR = 3.0
DROP_Q = 20
MIN_CELLS = 500
MIN_PAIRS = 200

RULES = {
    "question": "are the Axin2 and Il1r1 loci open in the same nuclei more or less often than comparable pairs of loci",
    "not_asked": ("whether they are correlated; at a fixed fragment budget every pair is positively correlated "
                  "through the breadth of that nucleus's open chromatin, so that question cannot be false"),
    "scope_limit": ("an accessible Axin2 locus is not a Wnt-responsive cell. Nothing here may be phrased as a "
                    "statement about Wnt-responsive or IL-1-responsive cells, about AEPs, or about Choi's "
                    "Il1r1-lineage population. The object is the pair of loci"),
    "R1_unit": "the GEM well, with GSE247130's inverted suffix map applied; Descriptive only",
    "R2_cells": "non-airway by within-cell ratio; primary set is Sftpc-positive non-transitional; the transitional group is secondary",
    "R3_depth": f"every cell downsampled to the {DROP_Q}th percentile of in-peak fragments over the analysis set",
    "R4_locus": "a locus is detected if at least one linked peak, promoter or distal, is detected after downsampling",
    "R5_statistic": "observed co-detection over the product of the marginals, as a ratio",
    "R6_null": f"{N_NULL} random gene pairs matched on linked peak count within a factor of 1.5, on different chromosomes from each other",
    "R7_controls": ("Krt8-Krt18 physical control must clear or the measurement is broken; Etv5-Abca3 biological "
                    "control must clear or there is no demonstrated sensitivity and negatives are not computable"),
    "R8_module_form": f"a Wnt-target set ({', '.join(WNT_SET)}) against Il1r1, reported beside the single-gene form",
    "R9_reading": "above the null, below the null, or not distinguishable from comparable pairs; never 'independent'",
    "U": {"U1": f"fewer than {MIN_CELLS} cells", "U2": "a locus outside 2 to 98 per cent detection",
          "U3": f"fewer than {MIN_PAIRS} usable null pairs", "U4": "physical control fails, well unreadable",
          "U5": "biological control fails, negatives not computable"},
}


def locus_index(matrix: Path, peaks: Path):
    """Every gene's linked peak slots, plus its chromosome."""
    feats = read_features(matrix)
    rna = feats[feats.feature_type == RNA]
    pk = feats[feats.feature_type == ATAC].reset_index(drop=True)
    slot = {iv: i for i, iv in enumerate(pk.interval)}
    sym2ens = dict(zip(rna.name, rna.id))
    ens2sym = {v: k for k, v in sym2ens.items()}
    ann = peak_gene_table(peaks)
    loci, chrom = {}, {}
    for ens, grp in ann.groupby("gene"):
        sym = ens2sym.get(ens)
        if sym is None:
            continue
        rows = sorted({slot[i] for i in grp.interval if i in slot})
        if rows:
            loci[sym] = np.array(rows, dtype=np.int64)
            chrom[sym] = grp.chrom.iloc[0] if "chrom" in grp else grp.interval.iloc[0].split(":")[0]
    return loci, chrom, len(pk)


def null_pairs(loci: dict, chrom: dict, ka: int, kb: int, rng) -> list:
    cand_a = [g for g, r in loci.items() if g not in NAMED and ka / 1.5 <= len(r) <= ka * 1.5]
    cand_b = [g for g, r in loci.items() if g not in NAMED and kb / 1.5 <= len(r) <= kb * 1.5]
    out = []
    for _ in range(N_NULL * 3):
        if len(out) >= N_NULL or not cand_a or not cand_b:
            break
        a = cand_a[int(rng.integers(len(cand_a)))]
        b = cand_b[int(rng.integers(len(cand_b)))]
        if a == b or chrom.get(a) == chrom.get(b):
            continue
        out.append((a, b))
    return out


def main() -> None:
    rec = RunRecord(OUT / "a1_run_record.json", "A1 Axin2 and Il1r1 locus co-accessibility", RULES)
    rng = np.random.default_rng(SEED)
    rows, detail = [], []

    for w in wells_corrected():
        if w["file"] == NEONATAL:
            continue
        rec.add_input(w["matrix"])
        loci, chrom, n_peaks = locus_index(w["matrix"], w["peaks"])
        missing = [g for g in NAMED if g in ("Axin2", "Il1r1", "Krt8", "Krt18", "Etv5", "Abca3")
                   and g not in loci]
        if missing:
            print(w["well"], "missing loci:", missing)
            continue

        feats = read_features(w["matrix"])
        bcs = read_barcodes(w["matrix"])
        symbols = sorted(set(list(NAMED) + ["Krt8", "Cldn4"]))
        rmap = gene_rows(feats, symbols)
        totals, picked = stream_selected(w["matrix"], np.array(sorted(rmap.values())))
        order = {r: i for i, r in enumerate(sorted(rmap.values()))}
        cnt = {s: picked[order[r]] for s, r in rmap.items()}

        sel = (bcs.suffix == w["suffix"]).values
        air = sum(cnt[a] for a in AIRWAY if a in cnt)
        alv = sum(cnt[a] for a in ("Sftpc", "Sftpb", "Lamp3") if a in cnt)
        keep = sel & (air <= alv)
        trans = keep & (cnt["Cldn4"] >= 1) & (cnt["Krt8"] >= 1)
        at2 = keep & (cnt["Sftpc"] >= 1) & ~trans

        for setname, mask in (("AT2", at2), ("transitional", trans)):
            idx = np.flatnonzero(mask)
            if len(idx) < MIN_CELLS:
                rows.append({"well": w["well"], "cell_set": setname, "cells": int(len(idx)),
                             "reading": "not computable", "why": f"U1, fewer than {MIN_CELLS} cells"})
                continue
            atac = totals.atac_counts.values
            B = int(np.percentile(atac[idx], DROP_Q))
            use = idx[atac[idx] >= B]
            print(f"\n{w['well']} [{setname}]: {len(use)} of {len(idx)} cells at B={B}")
            det, kept, npk = stream_downsampled_detection(w["matrix"], use, B, seed=SEED)
            use = use[kept]
            det = [d for d, k in zip(det, kept) if k]
            n = len(use)

            # peak-by-locus incidence, so one sparse product gives every locus at once
            names = sorted(loci)
            col = {g: i for i, g in enumerate(names)}
            r_, c_ = [], []
            for g in names:
                for p in loci[g]:
                    r_.append(p); c_.append(col[g])
            P = sp.csr_matrix((np.ones(len(r_), np.int8), (r_, c_)), shape=(npk, len(names)))
            hits = np.zeros((n, len(names)), dtype=bool)
            for i, d in enumerate(det):
                if len(d):
                    hits[i] = np.asarray(P[d].sum(axis=0)).ravel() > 0

            def frac(g):
                if g == "WNT_MODULE":
                    return hits[:, [col[x] for x in WNT_SET if x in col]].any(axis=1)
                return hits[:, col[g]]

            def ratio(a, b):
                A, Bv = frac(a), frac(b)
                pa, pb = A.mean(), Bv.mean()
                if not (0.02 <= pa <= 0.98 and 0.02 <= pb <= 0.98):
                    return None, pa, pb, float("nan")
                obs = float((A & Bv).mean())
                exp = float(pa * pb)
                return (obs / exp if exp > 0 else float("nan")), float(pa), float(pb), obs

            np_pairs = null_pairs(loci, chrom, len(loci["Axin2"]), len(loci["Il1r1"]), rng)
            nulls = []
            for a, b in np_pairs:
                r, _, _, _ = ratio(a, b)
                if r is not None and np.isfinite(r):
                    nulls.append(r)
            nulls = np.array(nulls)
            if len(nulls) < MIN_PAIRS:
                rows.append({"well": w["well"], "cell_set": setname, "cells": n,
                             "reading": "not computable",
                             "why": f"U3, only {len(nulls)} usable null pairs"})
                continue
            nm, nsd = float(nulls.mean()), float(nulls.std(ddof=1))

            res = {}
            for a, b, role in PAIRS:
                r, pa, pb, obs = ratio(a, b)
                if r is None:
                    res[(a, b)] = None
                    detail.append({"well": w["well"], "cell_set": setname, "pair": f"{a}|{b}",
                                   "role": role, "p_a": pa, "p_b": pb, "reading": "U2"})
                    continue
                z = (r - nm) / nsd if nsd > 0 else float("nan")
                outside = bool(r > nulls.max() or r < nulls.min())
                clears = bool(abs(z) >= Z_FLOOR and outside)
                res[(a, b)] = (r, z, clears)
                detail.append({"well": w["well"], "cell_set": setname, "pair": f"{a}|{b}",
                               "role": role, "p_a": round(pa, 4), "p_b": round(pb, 4),
                               "observed_both": round(obs, 4), "expected": round(pa * pb, 4),
                               "ratio": round(r, 4), "null_mean": round(nm, 4),
                               "null_sd": round(nsd, 4), "null_pairs": len(nulls),
                               "percentile": round(float((nulls < r).mean()), 4),
                               "z": round(z, 3), "clears": clears})
                print(f"   {a:11s} x {b:6s}  p {pa:.3f}/{pb:.3f}  ratio {r:.3f}  "
                      f"null {nm:.3f}+-{nsd:.3f}  z {z:+.2f}  {'CLEARS' if clears else ''}")

            phys = res.get(("Krt8", "Krt18"))
            bio = res.get(("Etv5", "Abca3"))
            row = {"well": w["well"], "cell_set": setname, "cells": n, "B_atac": B,
                   "null_mean": round(nm, 4), "null_sd": round(nsd, 4), "null_pairs": len(nulls)}
            if not phys or not phys[2] or phys[0] < nm:
                row.update(reading="unreadable", why="U4, the physical control Krt8-Krt18 did not clear")
            elif not bio or not bio[2] or bio[0] < nm:
                row.update(reading="not computable",
                           why="U5, the biological control Etv5-Abca3 did not clear, so there is no demonstrated sensitivity")
            else:
                for key, label in ((("Axin2", "Il1r1"), "axin2"), (("WNT_MODULE", "Il1r1"), "wnt_module")):
                    v = res.get(key)
                    if v is None:
                        row[label] = "not computable (U2)"
                    elif v[2]:
                        row[label] = ("open together MORE often than comparable pairs" if v[1] > 0
                                      else "open together LESS often than comparable pairs")
                    else:
                        row[label] = "not distinguishable from comparable pairs"
                    row[label + "_ratio"] = round(v[0], 3) if v else None
                    row[label + "_z"] = round(v[1], 2) if v else None
                row["reading"] = row.get("axin2")
                row["why"] = "R9, read against 300 matched random gene pairs, both controls cleared"
            rows.append(row)
            print("   ->", row.get("reading"))

    res_df, det_df = pd.DataFrame(rows), pd.DataFrame(detail)
    res_df.to_csv(OUT / "a1_readings.csv", index=False)
    det_df.to_csv(OUT / "a1_pairs.csv", index=False)
    rec.add_output(OUT / "a1_readings.csv")
    rec.add_output(OUT / "a1_pairs.csv")
    rec.set("readings", res_df.to_dict("records"))
    rec.finish()
    print("\n" + res_df.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
