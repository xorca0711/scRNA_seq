#!/usr/bin/env python
"""Trial A1c: the same question with a statistic that uses more than one bit per locus.

A1 refused because its physical control was not the arithmetic identity it was
designed to be. A1b replaced that control with a split-locus control, two random
halves of ONE gene's peak set, which must be co-accessible because they sit in
one regulatory domain. A1b then refused too, and the reason is worth more than
either trial: the split-locus control does not clear. Median ratios of 1.05 to
1.31 against null spreads of 0.083 to 0.255, at z of +0.51 to +2.05, never
reaching the floor. **If the design cannot detect co-accessibility between two
halves of the same locus, it cannot detect it between two different loci.**

Neither trial's thresholds are moved and both outcomes stand.

WHY A THIRD PASS IS NOT THRESHOLD SHOPPING, and the distinction is the whole
justification for this file. A1 and A1b reduced each locus to ONE BIT: detected
or not. A locus with fifteen or forty linked peaks carries far more information
than that, and the binary form discards essentially all of it, which is why the
null spread swamps an effect the control shows is really there. This trial
changes the INSTRUMENT, not the threshold, and the arbiter of whether the new
instrument is admissible is the split-locus control, not the Axin2 result. The
decision rule is fixed here, before running:

    If the split-locus control clears under the graded statistic, the instrument
    can see co-accessibility and the Axin2 reading may be taken. If it does not
    clear, ROUTE B IS CLOSED on these deposits and this folder says so, with no
    fourth pass and no further tuning.

The Axin2 against Il1r1 result plays no part in that decision. It has already
been seen twice, at 0.917, 0.921 and 0.960 in A1 and 0.917, 0.921 and 0.960 in
A1b, below the null mean in every computable AT2 well and never close to the
floor. Nothing about the graded statistic was chosen to move it.

R5c THE GRADED STATISTIC. For each cell and each locus, the score is the
    FRACTION of that locus's linked peaks detected after downsampling, a number
    in [0, 1] rather than a bit. The statistic for a pair is the Pearson
    correlation of those two scores across cells. The null is the same
    correlation for the same 300 matched random gene pairs, which absorbs the
    fact that at a fixed fragment budget every pair of loci is positively
    correlated through the breadth of that nucleus's open chromatin.

    Everything else is A1's: the wells, the cells, the depth budget, the
    matched-pair null construction, the biological control, section U and the
    reading. The budget in particular is not raised.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
TRIALS = HERE.parents[1] / "datp_epigenetics" / "trials"
sys.path.insert(0, str(TRIALS))
from multiome_utils import (RunRecord, read_barcodes, read_features,  # noqa: E402
                            stream_downsampled_detection, stream_selected)
from m1_closed_or_merely_silenced import AIRWAY, gene_rows  # noqa: E402
from m1b_corrected_pass import SEED, wells_corrected  # noqa: E402
from m1d_mature_wells_only import NEONATAL  # noqa: E402
from a1_locus_co_accessibility import (DROP_Q, MIN_CELLS, MIN_PAIRS, NAMED,  # noqa: E402
                                       PAIRS, WNT_SET, Z_FLOOR, locus_index,
                                       null_pairs)
from a1b_split_locus_control import N_SPLIT  # noqa: E402

OUT = HERE / "a1c_graded_statistic"
OUT.mkdir(exist_ok=True)

RULES = {
    "question": "A1's question with a graded statistic rather than one bit per locus",
    "why_A1b_refused": ("the split-locus control, two halves of one gene's peak set, did not clear: median "
                        "ratios 1.05 to 1.31 against null spreads 0.083 to 0.255. A design that cannot see "
                        "co-accessibility within one locus cannot see it between two"),
    "R5c_graded": ("per cell and locus, the FRACTION of that locus's linked peaks detected; the statistic is "
                   "the Pearson correlation of two such scores across cells; the null is the same for 300 "
                   "matched random gene pairs"),
    "decision_rule_fixed_before_running": ("if the split-locus control clears under the graded statistic the "
                                           "Axin2 reading may be taken; if it does not, Route B is CLOSED on "
                                           "these deposits, with no fourth pass and no further tuning"),
    "the_axin2_result_plays_no_part": ("it has been seen twice at 0.917, 0.921 and 0.960, below the null mean "
                                       "in every computable AT2 well and never close to the floor"),
    "unchanged": "wells, cells, depth budget, null construction, biological control, section U, reading",
}


def main() -> None:
    rec = RunRecord(OUT / "a1c_run_record.json", "A1c graded statistic", RULES)
    rng = np.random.default_rng(SEED)
    rows, detail = [], []

    for w in wells_corrected():
        if w["file"] == NEONATAL:
            continue
        rec.add_input(w["matrix"])
        loci, chrom, n_peaks = locus_index(w["matrix"], w["peaks"])
        feats = read_features(w["matrix"])
        bcs = read_barcodes(w["matrix"])
        rmap = gene_rows(feats, sorted(set(list(NAMED) + ["Krt8", "Cldn4"])))
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
            det, kept, npk = stream_downsampled_detection(w["matrix"], use, B, seed=SEED)
            det = [d for d, k in zip(det, kept) if k]
            n = len(det)
            print(f"\n{w['well']} [{setname}]: {n} cells at B={B}")

            k_ax = len(loci["Axin2"])
            cand = [g for g, r in loci.items()
                    if g not in NAMED and 2 * k_ax / 1.5 <= len(r) <= 2 * k_ax * 1.5]
            splits, split_names = {}, []
            for g in cand[:N_SPLIT * 3]:
                if len(split_names) >= N_SPLIT:
                    break
                pks = loci[g].copy(); rng.shuffle(pks); h = len(pks) // 2
                splits[g + "__A"], splits[g + "__B"] = pks[:h], pks[h:2 * h]
                split_names.append(g)

            all_loci = dict(loci); all_loci.update(splits)
            names = sorted(all_loci)
            col = {g: i for i, g in enumerate(names)}
            sizes = np.array([len(all_loci[g]) for g in names], dtype=float)
            r_, c_ = [], []
            for g in names:
                for p in all_loci[g]:
                    r_.append(int(p)); c_.append(col[g])
            P = sp.csr_matrix((np.ones(len(r_), np.int8), (r_, c_)), shape=(npk, len(names)))

            # R5c: the graded score, fraction of each locus's peaks detected
            score = np.zeros((n, len(names)), dtype=np.float32)
            for i, d in enumerate(det):
                if len(d):
                    score[i] = np.asarray(P[d].sum(axis=0)).ravel() / sizes
            sc = score - score.mean(axis=0, keepdims=True)
            sd = sc.std(axis=0)

            def corr(a, b):
                if a == "WNT_MODULE":
                    ia = [col[x] for x in WNT_SET if x in col]
                    va = score[:, ia].mean(axis=1)
                    va = va - va.mean()
                    sa = va.std()
                else:
                    ia = col[a]; va = sc[:, ia]; sa = sd[ia]
                ib = col[b]; vb = sc[:, ib]; sb = sd[ib]
                if sa <= 0 or sb <= 0:
                    return None
                return float((va * vb).mean() / (sa * sb))

            nulls = []
            for a, b in null_pairs(loci, chrom, k_ax, len(loci["Il1r1"]), rng):
                c = corr(a, b)
                if c is not None and np.isfinite(c):
                    nulls.append(c)
            nulls = np.array(nulls)
            if len(nulls) < MIN_PAIRS:
                rows.append({"well": w["well"], "cell_set": setname, "cells": n,
                             "reading": "not computable", "why": f"U3, only {len(nulls)} null pairs"})
                continue
            nm, nsd = float(nulls.mean()), float(nulls.std(ddof=1))

            sp_c = [corr(g + "__A", g + "__B") for g in split_names]
            sp_c = [c for c in sp_c if c is not None and np.isfinite(c)]
            sr = float(np.median(sp_c)) if sp_c else float("nan")
            sz = (sr - nm) / nsd if nsd > 0 else float("nan")
            split_clears = bool(sz >= Z_FLOOR)
            print(f"   split-locus control: median r {sr:.4f} over {len(sp_c)} genes, "
                  f"null {nm:.4f}+-{nsd:.4f}, z {sz:+.2f} {'CLEARS' if split_clears else 'FAILS'}")

            res = {}
            for a, b, role in PAIRS:
                c = corr(a, b)
                if c is None:
                    res[(a, b)] = None
                    continue
                z = (c - nm) / nsd if nsd > 0 else float("nan")
                outside = bool(c > nulls.max() or c < nulls.min())
                res[(a, b)] = (c, z, bool(abs(z) >= Z_FLOOR and outside))
                detail.append({"well": w["well"], "cell_set": setname, "pair": f"{a}|{b}",
                               "role": role, "correlation": round(c, 5),
                               "null_mean": round(nm, 5), "null_sd": round(nsd, 5),
                               "percentile": round(float((nulls < c).mean()), 4),
                               "z": round(z, 3), "clears": res[(a, b)][2]})
                print(f"   {a:11s} x {b:6s}  r {c:+.4f}  z {z:+.2f}  "
                      f"{'CLEARS' if res[(a, b)][2] else ''}")
            detail.append({"well": w["well"], "cell_set": setname, "pair": "split-locus median",
                           "role": "control, split locus", "correlation": round(sr, 5),
                           "null_mean": round(nm, 5), "null_sd": round(nsd, 5),
                           "z": round(sz, 3), "clears": split_clears})

            bio = res.get(("Etv5", "Abca3"))
            row = {"well": w["well"], "cell_set": setname, "cells": n, "B_atac": B,
                   "null_mean": round(nm, 5), "null_sd": round(nsd, 5), "null_pairs": len(nulls),
                   "split_control_r": round(sr, 4), "split_control_z": round(sz, 2)}
            if not split_clears:
                row.update(reading="not computable",
                           why="R5c decision rule: the split-locus control did not clear, so Route B is closed on these deposits")
            elif not bio or not bio[2] or bio[0] < nm:
                row.update(reading="not computable", why="U5, the biological control Etv5-Abca3 did not clear")
            else:
                for key, label in ((("Axin2", "Il1r1"), "axin2"), (("WNT_MODULE", "Il1r1"), "wnt_module")):
                    v = res.get(key)
                    row[label + "_r"] = round(v[0], 5) if v else None
                    row[label + "_z"] = round(v[1], 2) if v else None
                    row[label] = ("not computable" if v is None else
                                  ("MORE co-accessible than comparable pairs" if v[2] and v[1] > 0 else
                                   "LESS co-accessible than comparable pairs" if v[2] else
                                   "not distinguishable from comparable pairs"))
                row["reading"] = row.get("axin2")
                row["why"] = "R9, against 300 matched random gene pairs, both controls cleared"
            rows.append(row)
            print("   ->", row.get("reading"))

    res_df, det_df = pd.DataFrame(rows), pd.DataFrame(detail)
    res_df.to_csv(OUT / "a1c_readings.csv", index=False)
    det_df.to_csv(OUT / "a1c_pairs.csv", index=False)
    rec.add_output(OUT / "a1c_readings.csv"); rec.add_output(OUT / "a1c_pairs.csv")
    rec.set("readings", res_df.to_dict("records"))
    rec.set("route_b_closed", not any(r.get("split_control_z", 0) and r["split_control_z"] >= Z_FLOOR
                                      for r in rows))
    rec.finish()
    print("\n" + res_df.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
