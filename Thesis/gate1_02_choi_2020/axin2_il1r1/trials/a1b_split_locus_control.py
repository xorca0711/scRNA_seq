#!/usr/bin/env python
"""Trial A1b: A1's question with a positive control that cannot fail for biological reasons.

A1 returned UNREADABLE in every well because its rule U4 fired: the physical
positive control, Krt8 against Krt18, did not clear the null. Those two genes sit
24 kilobases apart and A1 assumed their linked peak sets would therefore overlap,
making the control trivially true. They do not overlap much, because the vendor
annotation assigns each peak to its NEAREST gene, so the peaks between them are
split between them rather than shared. The control was a real biological question
wearing the costume of an arithmetic identity, and it was the wrong instrument.

A1's thresholds are not moved and its outcome stands beside this one. Across the
eight wells Krt8 against Krt18 returned ratios of 1.11 to 1.25 against null means
near 1.03, consistently elevated and consistently short of the floor, which is
the behaviour of a real but modest effect rather than of a broken measurement.
That is informative and it is why this trial keeps reporting it, as a diagnostic
rather than as a gate.

THE ONE CHANGE. Rule R7's physical control is replaced by a SPLIT-LOCUS control.
Take one gene's linked peak set, split it at random into two halves, and treat
the halves as two loci. Those halves sit in one regulatory domain by
construction, so if the measurement can detect co-accessibility at all it must
detect it here. Twenty such controls are built, on twenty different genes chosen
to give halves matched in size to Axin2's peak count, and the MEDIAN of their
ratios is the gate. A median over twenty removes the risk of one locus that
happens to span two independent domains.

NOTHING ELSE MOVES. The cells, the depth budget, the statistic, the null, the
biological control, the reading and section U are A1's, unchanged. In particular
the budget stays at the twentieth percentile even though A1 showed the null
spread is wide in the shallower deposit, because raising it after seeing the
direction of the Axin2 result would be choosing a parameter to suit an answer.
If the split-locus control also fails at this budget, that is the answer: the
approach is underpowered at this depth, and this folder will say so rather than
tune until something clears.

WHAT A1 ALREADY SHOWED, recorded here so the corrected pass cannot be read as
the first look. In the AT2 sets the Axin2 against Il1r1 ratio came out 0.917,
0.921 and 0.960, below the null mean in every computable well, at z of -1.01,
-1.33 and -0.70. None of that was readable, because the control had failed. It
is repeated here so that a reader can see the corrected pass did not discover a
direction it had not already been shown.
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
from multiome_utils import (RunRecord, df_to_markdown, read_barcodes,  # noqa: E402
                            read_features, stream_downsampled_detection,
                            stream_selected)
from m1_closed_or_merely_silenced import AIRWAY, gene_rows  # noqa: E402
from m1b_corrected_pass import SEED, wells_corrected  # noqa: E402
from m1d_mature_wells_only import NEONATAL  # noqa: E402
from a1_locus_co_accessibility import (DROP_Q, MIN_CELLS, MIN_PAIRS, N_NULL,  # noqa: E402
                                       NAMED, PAIRS, WNT_SET, Z_FLOOR,
                                       locus_index, null_pairs)

OUT = HERE / "a1b_split_locus_control"
OUT.mkdir(exist_ok=True)
N_SPLIT = 20

RULES = {
    "question": "A1's question, with a positive control that cannot fail for biological reasons",
    "why_A1_refused": ("its physical control Krt8 against Krt18 assumed the two loci share peaks; the vendor "
                       "annotation assigns each peak to its nearest gene, so they largely do not, and the "
                       "control was a real biological question rather than an arithmetic identity"),
    "R7b_split_locus_control": (f"{N_SPLIT} genes, each with its linked peak set split at random into two halves "
                                "treated as two loci; the halves sit in one regulatory domain by construction, so "
                                "the measurement must detect them; the median of their ratios is the gate"),
    "unchanged": "cells, depth budget, statistic, null, biological control, reading and section U are A1's",
    "budget_not_raised": ("A1 showed the null spread is wide in the shallower deposit, but raising the budget "
                          "after seeing the direction of the Axin2 result would be choosing a parameter to suit "
                          "an answer. If the split-locus control also fails, the approach is underpowered at this "
                          "depth and this folder says so"),
    "what_A1_already_showed": ("Axin2 against Il1r1 came out 0.917, 0.921 and 0.960, below the null mean in every "
                               "computable AT2 well, at z of -1.01, -1.33 and -0.70, none of it readable"),
}


def main() -> None:
    rec = RunRecord(OUT / "a1b_run_record.json", "A1b split-locus control", RULES)
    rng = np.random.default_rng(SEED)
    rows, detail = [], []

    for w in wells_corrected():
        if w["file"] == NEONATAL:
            continue
        rec.add_input(w["matrix"])
        loci, chrom, n_peaks = locus_index(w["matrix"], w["peaks"])
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
            det, kept, npk = stream_downsampled_detection(w["matrix"], use, B, seed=SEED)
            det = [d for d, k in zip(det, kept) if k]
            n = len(det)
            print(f"\n{w['well']} [{setname}]: {n} cells at B={B}")

            # R7b: build the split-locus controls before anything else is scored
            k_axin2 = len(loci["Axin2"])
            cand = [g for g, r in loci.items()
                    if g not in NAMED and 2 * k_axin2 / 1.5 <= len(r) <= 2 * k_axin2 * 1.5]
            split_names = []
            splits = {}
            for g in cand[:N_SPLIT * 3]:
                if len(split_names) >= N_SPLIT:
                    break
                pks = loci[g].copy()
                rng.shuffle(pks)
                h = len(pks) // 2
                splits[g + "__A"] = pks[:h]
                splits[g + "__B"] = pks[h:2 * h]
                split_names.append(g)

            all_loci = dict(loci)
            all_loci.update(splits)
            names = sorted(all_loci)
            col = {g: i for i, g in enumerate(names)}
            r_, c_ = [], []
            for g in names:
                for p in all_loci[g]:
                    r_.append(int(p)); c_.append(col[g])
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
                pa, pb = float(A.mean()), float(Bv.mean())
                if not (0.02 <= pa <= 0.98 and 0.02 <= pb <= 0.98):
                    return None, pa, pb, float("nan")
                obs, exp = float((A & Bv).mean()), pa * pb
                return (obs / exp if exp > 0 else float("nan")), pa, pb, obs

            nulls = []
            for a, b in null_pairs(loci, chrom, k_axin2, len(loci["Il1r1"]), rng):
                r, _, _, _ = ratio(a, b)
                if r is not None and np.isfinite(r):
                    nulls.append(r)
            nulls = np.array(nulls)
            if len(nulls) < MIN_PAIRS:
                rows.append({"well": w["well"], "cell_set": setname, "cells": n,
                             "reading": "not computable", "why": f"U3, only {len(nulls)} null pairs"})
                continue
            nm, nsd = float(nulls.mean()), float(nulls.std(ddof=1))

            split_ratios = []
            for g in split_names:
                r, _, _, _ = ratio(g + "__A", g + "__B")
                if r is not None and np.isfinite(r):
                    split_ratios.append(r)
            sr = float(np.median(split_ratios)) if split_ratios else float("nan")
            sz = (sr - nm) / nsd if nsd > 0 else float("nan")
            split_clears = bool(abs(sz) >= Z_FLOOR and sz > 0)
            print(f"   split-locus control: median ratio {sr:.3f} over {len(split_ratios)} genes, "
                  f"null {nm:.3f}+-{nsd:.3f}, z {sz:+.2f} {'CLEARS' if split_clears else 'FAILS'}")

            res = {}
            for a, b, role in PAIRS:
                r, pa, pb, obs = ratio(a, b)
                if r is None:
                    res[(a, b)] = None
                    continue
                z = (r - nm) / nsd if nsd > 0 else float("nan")
                outside = bool(r > nulls.max() or r < nulls.min())
                res[(a, b)] = (r, z, bool(abs(z) >= Z_FLOOR and outside))
                detail.append({"well": w["well"], "cell_set": setname, "pair": f"{a}|{b}",
                               "role": role, "p_a": round(pa, 4), "p_b": round(pb, 4),
                               "observed_both": round(obs, 4), "ratio": round(r, 4),
                               "null_mean": round(nm, 4), "null_sd": round(nsd, 4),
                               "percentile": round(float((nulls < r).mean()), 4),
                               "z": round(z, 3), "clears": res[(a, b)][2]})
                print(f"   {a:11s} x {b:6s}  ratio {r:.3f}  z {z:+.2f}  "
                      f"{'CLEARS' if res[(a, b)][2] else ''}")
            detail.append({"well": w["well"], "cell_set": setname, "pair": "split-locus median",
                           "role": "control, split locus", "ratio": round(sr, 4),
                           "null_mean": round(nm, 4), "null_sd": round(nsd, 4),
                           "z": round(sz, 3), "clears": split_clears})

            bio = res.get(("Etv5", "Abca3"))
            row = {"well": w["well"], "cell_set": setname, "cells": n, "B_atac": B,
                   "null_mean": round(nm, 4), "null_sd": round(nsd, 4), "null_pairs": len(nulls),
                   "split_control_ratio": round(sr, 3), "split_control_z": round(sz, 2),
                   "krt8_krt18_ratio": round(res[("Krt8", "Krt18")][0], 3) if res.get(("Krt8", "Krt18")) else None}
            if not split_clears:
                row.update(reading="not computable",
                           why="R7b, the split-locus control did not clear, so the design has no demonstrated "
                               "sensitivity to co-accessibility at this depth")
            elif not bio or not bio[2] or bio[0] < nm:
                row.update(reading="not computable",
                           why="U5, the biological control Etv5-Abca3 did not clear")
            else:
                for key, label in ((("Axin2", "Il1r1"), "axin2"), (("WNT_MODULE", "Il1r1"), "wnt_module")):
                    v = res.get(key)
                    row[label + "_ratio"] = round(v[0], 3) if v else None
                    row[label + "_z"] = round(v[1], 2) if v else None
                    row[label] = ("not computable (U2)" if v is None else
                                  ("open together MORE often than comparable pairs" if v[2] and v[1] > 0 else
                                   "open together LESS often than comparable pairs" if v[2] else
                                   "not distinguishable from comparable pairs"))
                row["reading"] = row.get("axin2")
                row["why"] = "R9, against 300 matched random gene pairs, both controls cleared"
            rows.append(row)
            print("   ->", row.get("reading"))

    res_df, det_df = pd.DataFrame(rows), pd.DataFrame(detail)
    res_df.to_csv(OUT / "a1b_readings.csv", index=False)
    det_df.to_csv(OUT / "a1b_pairs.csv", index=False)
    rec.add_output(OUT / "a1b_readings.csv")
    rec.add_output(OUT / "a1b_pairs.csv")
    rec.set("readings", res_df.to_dict("records"))
    rec.finish()
    print("\n" + res_df.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
