#!/usr/bin/env python
"""Trial M2: does M1e's reading survive its own arithmetic?

M1e read, in two wells of two deposits, that the AT2 identity programme is
silenced in RNA while its distal chromatin does not move. A null on one arm is
only worth anything if it is an absence of effect rather than an absence of
power, and an arm mean is only worth anything if it is not one gene wearing a
gene set. M2 attacks both, and it attacks nothing else: it changes no rule, adds
no well and asks no new biological question. Every number here is a diagnostic
on a reading that already exists.

WHAT IT COMPUTES, all on M1e's own pipeline, settings and seeds unless stated.

 1. PER GENE. The RNA and ATAC difference for every gene of every arm, with the
    number of distal peaks the gene contributed. An arm mean carried by one gene
    is not an arm.

 2. LEAVE ONE OUT. Each arm recomputed with each of its genes removed in turn,
    against the sham band rebuilt the same way. If dropping any single gene
    flips whether an arm clears, the arm is that gene.

 3. THE OFFSET, WHICH WAS INERT. M1e subtracted the genome-wide MEDIAN distal
    difference and that median came out exactly 0.0 in every well, because more
    than half of all distal peaks are detected in neither group at this budget.
    The intended correction therefore did nothing. M2 reports the MEAN global
    distal difference, the trimmed mean, and the arm values recomputed against
    each, so the reading can be checked against a correction that is not
    degenerate. This matters in one direction specifically: the sham band cannot
    absorb a labelled-versus-reference global shift, because both sham groups
    are reference cells.

 4. SEEDS. The whole pipeline rerun at three further downsampling seeds. The
    ATAC downsample is one random draw per cell, so a reading that moves with
    the seed is a reading about the draw.

 5. AN EQUIVALENCE STATEMENT FOR THE NULL. Rather than reporting that the AT2
    ATAC arm failed to clear, M2 reports the largest effect the data would have
    hidden: the sham standard deviation times three, which is the smallest value
    the frozen rule would have called a clearance. Saying "the AT2 arm moved by
    less than X" is a claim; saying "the AT2 arm did not move" is not one the
    data can support and this trial does not make it.

NOTHING HERE MAY CHANGE M1e's READING. If a diagnostic contradicts it, the
contradiction is reported and M1e's reading stands in the record beside it, as
the C1/C1b precedent requires.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from multiome_utils import (RunRecord, detect_prob_at_depth,  # noqa: E402
                            detection_fraction, df_to_markdown, read_barcodes,
                            read_features, stream_downsampled_detection,
                            stream_selected)
from m1_closed_or_merely_silenced import (AIRWAY, ARMS, AT1, AT2_IDENTITY,  # noqa: E402
                                          LABELS, N_SHAM, TRANSITIONAL,
                                          Z_FLOOR, gene_rows)
from m1b_corrected_pass import SEED, arm_value, build_links, wells_corrected
from m1d_mature_wells_only import MIN_CELLS, NEONATAL
from m1e_per_well_budget import DROP_Q

OUT = HERE / "m2_robustness_of_the_m1e_reading"
OUT.mkdir(exist_ok=True)
SEEDS = [SEED, 7, 13, 29]

RULES = {
    "question": "does M1e's reading survive per-gene decomposition, leave-one-out, a non-degenerate offset and a change of seed",
    "changes_no_rule": True,
    "may_not_change_M1e": "if a diagnostic contradicts M1e, the contradiction is reported and M1e's reading stands beside it",
    "equivalence": ("the null on the AT2 ATAC arm is reported as the largest effect the frozen rule would have "
                    "hidden, three sham standard deviations, not as an absence of movement"),
    "seeds": SEEDS,
}


def prepare():
    all_wells = wells_corrected()
    symbols = sorted(set(AT2_IDENTITY + TRANSITIONAL + AT1 + AIRWAY + LABELS))
    files: dict = {}
    rng = np.random.default_rng(SEED)
    for w in all_wells:
        key = str(w["matrix"])
        if key in files:
            continue
        feats = read_features(w["matrix"])
        bcs = read_barcodes(w["matrix"])
        rmap = gene_rows(feats, symbols)
        totals, picked = stream_selected(w["matrix"], np.array(sorted(rmap.values())))
        order = {r: i for i, r in enumerate(sorted(rmap.values()))}
        files[key] = {"bcs": bcs, "totals": totals, "deposit": w["deposit"],
                      "raw": {s: picked[order[r]].astype(np.int64) for s, r in rmap.items()}}
        print("read", w["file"])
    budgets = {}
    for dep in ("GSE310539", "GSE247130"):
        rq = []
        for w in all_wells:
            if w["deposit"] != dep:
                continue
            f = files[str(w["matrix"])]
            sel = (f["bcs"].suffix == w["suffix"]).values
            rq.append(np.percentile(f["totals"].rna_counts.values[sel], 40))
        budgets[dep] = int(min(rq))
    for f in files.values():
        B = budgets[f["deposit"]]
        tot = f["totals"].rna_counts.values
        ok = tot >= B
        f["rna_ok"] = ok
        f["ds"] = {}
        for s, c in f["raw"].items():
            out = np.zeros(len(c), dtype=np.int64)
            good = np.clip(c[ok], 0, tot[ok])
            out[ok] = rng.hypergeometric(np.maximum(good, 0), np.maximum(tot[ok] - good, 0), B)
            f["ds"][s] = out
        air = sum(f["raw"][a] for a in AIRWAY if a in f["raw"])
        alv = sum(f["raw"][a] for a in ("Sftpc", "Sftpb", "Lamp3") if a in f["raw"])
        f["airway"] = air > alv
    for w in all_wells:
        f = files[str(w["matrix"])]
        sel = (f["bcs"].suffix == w["suffix"]).values
        keep = sel & f["rna_ok"] & ~f["airway"]
        w["_trans"] = keep & (f["ds"]["Cldn4"] >= 1) & (f["ds"]["Krt8"] >= 1)
        w["_ref"] = keep & (f["ds"]["Sftpc"] >= 1) & ~w["_trans"]
        w["_neonatal"] = w["file"] == NEONATAL
        w["_mutant"] = "Cebpa_mutant" in w["well"]
    return all_wells, files, budgets


def main() -> None:
    rec = RunRecord(OUT / "m2_run_record.json", "M2 robustness of the M1e reading", RULES)
    all_wells, files, budgets = prepare()

    gene_rows_out, loo_rows, offset_rows, seed_rows = [], [], [], []
    for w in all_wells:
        if w["_neonatal"] or not w["injured"] and "7wk_Cebpa" not in w["well"]:
            continue
        f = files[str(w["matrix"])]
        trans, ref = w["_trans"], w["_ref"]
        if int(trans.sum()) < MIN_CELLS:
            continue
        atac = f["totals"].atac_counts.values
        B = int(np.percentile(atac[trans], DROP_Q))
        t_ok, r_ok = trans & (atac >= B), ref & (atac >= B)
        m, nref = int(t_ok.sum()), int(r_ok.sum())
        if m < MIN_CELLS or nref < 2 * m:
            continue
        q = min(4, nref // m - 1)
        links, retained, _ = build_links(w["matrix"], w["peaks"])
        arms = {a: [g for g in gs if not (w["_mutant"] and g == "Cebpa")] for a, gs in ARMS.items()}
        Brna = budgets[w["deposit"]]
        tot = f["totals"].rna_counts.values
        print(f"\n=== {w['well']}: {m} labelled, B_atac={B}, q={q}")

        for seed in SEEDS:
            rng = np.random.default_rng(seed)
            ti = np.flatnonzero(t_ok)
            ri = rng.permutation(np.flatnonzero(r_ok))[:(q + 1) * m]
            cells = np.concatenate([ti, ri])
            det_cells, kept, n_peaks = stream_downsampled_detection(w["matrix"], cells, B, seed=seed)
            idx = np.arange(len(cells))
            tmask = idx < m
            pool = idx[idx >= m]
            det_t = detection_fraction(det_cells, tmask, n_peaks)
            det_r = detection_fraction(det_cells, np.isin(idx, pool[:q * m]), n_peaks)
            d = (det_t - det_r)[retained]
            off_med, off_mean = float(np.median(d)), float(np.mean(d))
            off_trim = float(np.mean(np.sort(d)[int(.05 * len(d)):int(.95 * len(d))]))

            def rna_arm(tc, rc, genes):
                vals = []
                for g in genes:
                    pt = detect_prob_at_depth(f["raw"][g][tc], tot[tc], Brna)
                    pr = detect_prob_at_depth(f["raw"][g][rc], tot[rc], Brna)
                    vals.append(float(np.nanmean(pt) - np.nanmean(pr)))
                return float(np.mean(vals)) if vals else float("nan")

            sham = {a: {"atac": [], "rna": []} for a in arms}
            sham_loo = {(a, g): [] for a, gs in arms.items() for g in gs}
            for s in range(N_SHAM):
                perm = np.random.default_rng(1000 + s).permutation(pool)
                pa, pb = perm[:m], perm[m:m + q * m]
                da = detection_fraction(det_cells, np.isin(idx, pa), n_peaks)
                db = detection_fraction(det_cells, np.isin(idx, pb), n_peaks)
                offs = float(np.median((da - db)[retained]))
                for a, gs in arms.items():
                    v, used = arm_value(da, db, links, gs, offs)
                    sham[a]["atac"].append(v)
                    sham[a]["rna"].append(rna_arm(cells[pa], cells[pb], used))
                    for g in gs:
                        vv, _ = arm_value(da, db, links, [x for x in gs if x != g], offs)
                        sham_loo[(a, g)].append(vv)

            def zed(obs, band):
                band = np.asarray(band, float)
                sd = float(np.nanstd(band, ddof=1))
                z = (obs - float(np.nanmean(band))) / sd if sd > 0 else float("nan")
                outside = bool(obs > np.nanmax(band) or obs < np.nanmin(band))
                return z, sd, bool(abs(z) >= Z_FLOOR and outside)

            for a, gs in arms.items():
                for off_name, off in (("median", off_med), ("mean", off_mean), ("trimmed", off_trim)):
                    obs, used = arm_value(det_t, det_r, links, gs, off)
                    z, sd, clears = zed(obs, sham[a]["atac"])
                    offset_rows.append({"well": w["well"], "seed": seed, "arm": a,
                                        "offset_kind": off_name, "offset": off,
                                        "observed": obs, "z": z, "clears": clears})
                obs_a, used = arm_value(det_t, det_r, links, gs, off_med)
                obs_r = rna_arm(cells[idx[tmask]], cells[pool[:q * m]], used)
                za, sda, ca = zed(obs_a, sham[a]["atac"])
                zr, sdr, cr = zed(obs_r, sham[a]["rna"])
                seed_rows.append({"well": w["well"], "seed": seed, "arm": a,
                                  "atac_obs": obs_a, "atac_z": za, "atac_clears": ca,
                                  "atac_detectable_floor_3sd": 3 * sda,
                                  "rna_obs": obs_r, "rna_z": zr, "rna_clears": cr})

            if seed != SEED:
                continue
            # per gene and leave one out, on the frozen seed only
            for a, gs in arms.items():
                for g in gs:
                    rowsg = links.get(g)
                    n_pk = 0 if rowsg is None else len(rowsg)
                    da = float(np.mean(det_t[rowsg] - det_r[rowsg])) - off_med if n_pk else float("nan")
                    pt = detect_prob_at_depth(f["raw"][g][cells[idx[tmask]]], tot[cells[idx[tmask]]], Brna)
                    pr = detect_prob_at_depth(f["raw"][g][cells[pool[:q * m]]], tot[cells[pool[:q * m]]], Brna)
                    gene_rows_out.append({"well": w["well"], "arm": a, "gene": g,
                                          "distal_peaks": n_pk, "atac_delta": da,
                                          "rna_delta": float(np.nanmean(pt) - np.nanmean(pr))})
                    obs, used = arm_value(det_t, det_r, links, [x for x in gs if x != g], off_med)
                    z, sd, clears = zed(obs, sham_loo[(a, g)])
                    loo_rows.append({"well": w["well"], "arm": a, "dropped_gene": g,
                                     "genes_left": len(used), "atac_obs": obs, "atac_z": z,
                                     "atac_clears": clears})

    genes_df = pd.DataFrame(gene_rows_out)
    loo_df = pd.DataFrame(loo_rows)
    off_df = pd.DataFrame(offset_rows)
    seed_df = pd.DataFrame(seed_rows)
    for name, df in [("m2_per_gene.csv", genes_df), ("m2_leave_one_out.csv", loo_df),
                     ("m2_offsets.csv", off_df), ("m2_seeds.csv", seed_df)]:
        df.to_csv(OUT / name, index=False)
        rec.add_output(OUT / name)

    print("\n--- per gene, frozen seed")
    print(genes_df.to_string(index=False))
    print("\n--- leave one out, ATAC")
    print(loo_df.to_string(index=False))
    print("\n--- offsets")
    print(off_df[off_df.seed == SEED].to_string(index=False))
    print("\n--- seeds")
    print(seed_df.to_string(index=False))

    rec.set("n_seeds", len(SEEDS))
    rec.finish()
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
