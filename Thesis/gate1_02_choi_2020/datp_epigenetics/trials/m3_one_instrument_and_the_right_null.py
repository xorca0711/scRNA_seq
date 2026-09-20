#!/usr/bin/env python
"""Trial M3: one instrument for both deposits, and the null the sham band could not be.

M2 withdrew M1e's chromatin reading and left the question Not established. A
five-lens adversarial review of that withdrawal then found three defects in M2
itself, and this trial implements their corrections. It does not revisit the
withdrawal, which stands; it asks whether "not established" is the right
description or whether a better instrument and a better null say something
sharper.

THE THREE DEFECTS, and this trial is honest that it did not find them itself.

D1 THE TWO DEPOSITS WERE NOT RUNNING THE SAME INSTRUMENT. Rule R12 removed
   Cebpa from the AT2 arm in the Cebpa-mutant wells, because it is genetically
   deleted there, and kept it everywhere else. Cebpa contributes only 3 distal
   peaks but the largest positive per-gene value in wildtype_SeV, and the review
   reports that it is the only reason that well's AT2 arm is positive at all.
   M2's own leave-one-out agrees: dropping Cebpa gives the most negative of the
   nine drops. Keeping a gene in one deposit and not the other means the two
   arms were never comparable.

D2 THE SHAM BAND IS THE WRONG NULL FOR THE QUESTION BEING ASKED. It permutes
   CELLS, so it asks whether this SPLIT is special. It can never ask whether
   these GENES are special, which is what a gene-set arm claims. The review
   drew peak-count-matched random gene sets and evaluated them on the real
   split: their mean is not zero. That also explains why M1e's genome-wide
   median offset came out exactly 0.0 and M2 called it absent rather than
   inert. It is 0 by construction, because a fifth to a third of arm peaks are
   detected in no cell of either group and the modal peak sits near 0.004.

D3 THE EQUIVALENCE STATISTIC WAS THE WRONG ONE AND ERRED UNSAFELY. M2 reported
   "three sham standard deviations" as the largest effect the rule would have
   hidden. Three sigma is the 50-per-cent-power detection floor, not a bound
   the data support, and it understates the hidden effect. A confidence
   interval is the right statistic and this trial reports one.

FROZEN RULES.

R17 ONE INSTRUMENT. Cebpa leaves the AT2 identity arm in EVERY well, not only
    in the Cebpa-mutant ones, so both deposits run the same eight genes: Etv5,
    Abca3, Lamp3, Sftpb, Slc34a2, Lyz2, Napsa, Pon1.

R18 THE GENE-SET NULL. For each arm and well, 300 random gene sets are drawn,
    matched to that arm gene by gene on the number of distal peaks the gene
    carries (nearest available count, sampling without replacement, excluding
    every gene named in any arm or in the label set), and evaluated on the REAL
    labelled-versus-reference split. The MEAN of that distribution replaces the
    genome-wide median as the ascertainment offset.

R19 AN ARM MUST NOW CLEAR TWO NULLS, NOT ONE.
      z_sham    = (observed - gene-set null mean) / sd of the cell-permutation shams
      z_geneset = (observed - gene-set null mean) / sd of the gene-set null
    plus the arm's percentile within the gene-set null. An arm clears only if
    both absolute z values reach 3. The two nulls answer different questions and
    an arm that beats one and not the other is reported as such.

R20 FIVE HUNDRED SHAMS, NOT TWENTY, and the equivalence statement is a
    confidence interval on the relative scale: the arm value and its interval
    divided by the reference group's own mean detection of that arm's peaks,
    which makes it a per cent of reference accessibility rather than a
    detection-fraction difference nobody can interpret.

R21 NO NORMAL-TAIL P VALUES. A z here is a standardised distance. The exact
    permutation resolution with N shams is 1/(N+1), and the percentile within
    each null is reported beside every z. A z of 12 is not a p of 1e-5.

R22 THE PREDICTIONS ARE PRE-REGISTERED, BECAUSE THE REVIEW MADE THEM. This
    trial implements corrections proposed by an adversarial review that also
    predicted specific values, so agreement confirms that review's arithmetic
    rather than discovering anything independently, and the run record says so.
    Predicted: the AT2 arm without Cebpa in wildtype_SeV goes to about -0.0034
    at z_sham about -1.21; gene-set null means about +0.0009 and +0.0021; the
    AT2 arm in SeV_Cebpa_mutant reaches z about -2.80 at percentile 0 of 300;
    re-centred AT2 z values about -2.66, +0.14 and -1.11; and the transitional
    arm survives re-centring at about +3.93 and +4.46.

WHAT THIS TRIAL MAY NOT DO. It may not resurrect M1e's reading. If the AT2 arm
clears both nulls in the closing direction the reading is that the programme
closes, which is the opposite of what was withdrawn. If it clears neither the
question stays Not established. Nothing here can return "silenced but not
closed", because that reading required a null to be read as retention and the
review has shown the arm is consistently negative under every weighting.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from multiome_utils import (ATAC, RNA, RunRecord, detect_prob_at_depth,  # noqa: E402
                            detection_fraction, df_to_markdown, peak_gene_table,
                            read_features, stream_downsampled_detection)
from m1_closed_or_merely_silenced import (AIRWAY, AT1, AT2_IDENTITY, LABELS,  # noqa: E402
                                          TRANSITIONAL, Z_FLOOR, excluded_peaks)
from m1b_corrected_pass import SEED, build_links, wells_corrected
from m1d_mature_wells_only import MIN_CELLS, NEONATAL
from m1e_per_well_budget import DROP_Q
from m2_robustness_of_the_m1e_reading import prepare

OUT = HERE / "m3_one_instrument_and_the_right_null"
OUT.mkdir(exist_ok=True)

# R17: one instrument, Cebpa out everywhere
AT2_COMMON = [g for g in AT2_IDENTITY if g != "Cebpa"]
ARMS3 = {"AT2_identity": AT2_COMMON, "transitional": TRANSITIONAL, "AT1": AT1}
N_SHAM3 = 500
N_GENESET = 300
RESERVED = set(AT2_IDENTITY + TRANSITIONAL + AT1 + AIRWAY + LABELS
               + ["Krt18", "Krt7", "Cldn3", "Sftpc"])

RULES = {
    "question": "does M2's Not established survive one common instrument and a gene-set null",
    "found_by": "a five-lens adversarial review of M2, not by this repository",
    "R17_one_instrument": f"Cebpa leaves the AT2 arm in every well: {', '.join(AT2_COMMON)}",
    "R18_geneset_null": f"{N_GENESET} random gene sets matched gene by gene on distal peak count, evaluated on the real split; their mean replaces the genome-wide median offset",
    "R19_two_nulls": "an arm must clear both the cell-permutation sham and the gene-set null at |z| >= 3",
    "R20_shams": f"{N_SHAM3}; the equivalence statement is a confidence interval on the relative scale",
    "R21_no_p_values": "a z is a standardised distance; the percentile within each null is reported beside it",
    "R22_predictions_preregistered": {
        "at2_wildtype_no_cebpa": -0.0034, "at2_wildtype_z_sham": -1.21,
        "geneset_null_means": [0.0009, 0.0021],
        "at2_sev_mutant_z_geneset": -2.80,
        "recentred_at2_z": [-2.66, 0.14, -1.11],
        "transitional_recentred_z": [3.93, 4.46],
    },
    "may_not_do": "resurrect M1e's reading; the arm is consistently negative under every weighting",
}


def gene_peak_pool(matrix: Path, peaks: Path) -> dict:
    """Every gene's distal peak slots, for the matched random gene sets."""
    feats = read_features(matrix)
    excl, _ = excluded_peaks(feats)
    pk = feats[feats.feature_type == ATAC].reset_index(drop=True)
    slot = {iv: i for i, iv in enumerate(pk.interval)}
    rna = feats[feats.feature_type == RNA]
    ens2sym = dict(zip(rna.id, rna.name))
    ann = peak_gene_table(peaks)
    ann = ann[(ann.peak_type == "distal") & (~ann.interval.isin(excl))]
    pool: dict[str, np.ndarray] = {}
    for ens, grp in ann.groupby("gene"):
        sym = ens2sym.get(ens)
        if sym is None or sym in RESERVED:
            continue
        rows = sorted({slot[i] for i in grp.interval if i in slot})
        if rows:
            pool[sym] = np.array(rows, dtype=np.int64)
    return pool


def matched_sets(pool: dict, arm_counts: list[int], n: int, seed: int) -> list[list[np.ndarray]]:
    """n random gene sets, each matched gene by gene on distal peak count."""
    rng = np.random.default_rng(seed)
    names = np.array(list(pool))
    counts = np.array([len(pool[g]) for g in names])
    order = np.argsort(counts)
    names, counts = names[order], counts[order]
    out = []
    for _ in range(n):
        chosen, used = [], set()
        for k in arm_counts:
            lo, hi = np.searchsorted(counts, k * 0.7), np.searchsorted(counts, k * 1.45, side="right")
            if hi <= lo:
                j = int(np.argmin(np.abs(counts - k)))
                lo, hi = max(0, j - 25), min(len(counts), j + 25)
            cand = [i for i in range(lo, hi) if names[i] not in used]
            if not cand:
                cand = [i for i in range(len(names)) if names[i] not in used]
            i = int(rng.choice(cand))
            used.add(names[i])
            chosen.append(pool[names[i]])
        out.append(chosen)
    return out


def raw_arm(det_t: np.ndarray, det_r: np.ndarray, rowsets: list[np.ndarray]) -> float:
    vals = [float(np.mean(det_t[r] - det_r[r])) for r in rowsets if len(r)]
    return float(np.mean(vals)) if vals else float("nan")


def main() -> None:
    rec = RunRecord(OUT / "m3_run_record.json", "M3 one instrument and the right null", RULES)
    all_wells, files, budgets = prepare()
    rows, arm_rows = [], []

    for w in all_wells:
        if w["_neonatal"]:
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
        rng = np.random.default_rng(SEED)
        ti = np.flatnonzero(t_ok)
        ri = rng.permutation(np.flatnonzero(r_ok))[:(q + 1) * m]
        cells = np.concatenate([ti, ri])
        role = "test" if w["injured"] else "negative control"
        print(f"\n=== {w['well']} ({role}): {m} labelled, B={B}, q={q}")

        det_cells, kept, n_peaks = stream_downsampled_detection(w["matrix"], cells, B, seed=SEED)
        links, retained, _ = build_links(w["matrix"], w["peaks"])
        pool = gene_peak_pool(w["matrix"], w["peaks"])
        idx = np.arange(len(cells))
        tmask = idx < m
        poolidx = idx[idx >= m]
        det_t = detection_fraction(det_cells, tmask, n_peaks)
        det_r = detection_fraction(det_cells, np.isin(idx, poolidx[:q * m]), n_peaks)

        # one binary cell-by-peak block over only the peaks any arm needs
        for arm, genes in ARMS3.items():
            rowsets = [links[g] for g in genes if g in links and len(links[g])]
            used = [g for g in genes if g in links and len(links[g])]
            if len(used) < 4:
                continue
            counts = [len(r) for r in rowsets]
            obs = raw_arm(det_t, det_r, rowsets)

            # R18 gene-set null, on the real split
            gs = matched_sets(pool, counts, N_GENESET, seed=SEED + 11)
            gvals = np.array([raw_arm(det_t, det_r, s) for s in gs], dtype=float)
            g_mean, g_sd = float(np.nanmean(gvals)), float(np.nanstd(gvals, ddof=1))
            pct = float((gvals < obs).mean())

            # R20 cell-permutation shams at the same group sizes
            svals = np.empty(N_SHAM3, dtype=float)
            for s in range(N_SHAM3):
                perm = np.random.default_rng(5000 + s).permutation(poolidx)
                da = detection_fraction(det_cells, np.isin(idx, perm[:m]), n_peaks)
                db = detection_fraction(det_cells, np.isin(idx, perm[m:m + q * m]), n_peaks)
                svals[s] = raw_arm(da, db, rowsets)
            s_sd = float(np.nanstd(svals, ddof=1))

            z_sham = (obs - g_mean) / s_sd if s_sd > 0 else float("nan")
            z_gs = (obs - g_mean) / g_sd if g_sd > 0 else float("nan")
            clears = bool(abs(z_sham) >= Z_FLOOR and abs(z_gs) >= Z_FLOOR)

            # R20 relative-scale interval, against the reference group's own baseline
            baseline = float(np.mean([float(np.mean(det_r[r])) for r in rowsets]))
            half = 1.96 * s_sd * np.sqrt(1 + 1.0 / N_SHAM3)
            arm_rows.append({
                "well": w["well"], "role": role, "arm": arm, "genes": len(used),
                "gene_list": ",".join(used), "distal_peaks": int(sum(counts)),
                "observed_raw": obs, "geneset_null_mean": g_mean, "geneset_null_sd": g_sd,
                "geneset_percentile": pct, "sham_sd": s_sd,
                "z_sham": z_sham, "z_geneset": z_gs, "clears_both": clears,
                "reference_baseline": baseline,
                "relative_pct": 100 * (obs - g_mean) / baseline if baseline else float("nan"),
                "relative_lo_pct": 100 * (obs - g_mean - half) / baseline if baseline else float("nan"),
                "relative_hi_pct": 100 * (obs - g_mean + half) / baseline if baseline else float("nan"),
            })
            print(f"  {arm:13s} obs {obs:+.5f}  gsnull {g_mean:+.5f}  "
                  f"z_sham {z_sham:+.2f}  z_gs {z_gs:+.2f}  pct {pct:.3f}  "
                  f"rel {100*(obs-g_mean)/baseline:+.1f}% [{100*(obs-g_mean-half)/baseline:+.1f},"
                  f"{100*(obs-g_mean+half)/baseline:+.1f}]  {'CLEARS' if clears else ''}")

        a = [r for r in arm_rows if r["well"] == w["well"]]
        at2 = next((r for r in a if r["arm"] == "AT2_identity"), None)
        tr = next((r for r in a if r["arm"] == "transitional"), None)
        if at2 and tr:
            if not tr["clears_both"]:
                reading, why = "not computable", "R9: the transitional arm does not clear both nulls"
            elif at2["clears_both"]:
                reading = ("AT2 programme CLOSES" if at2["z_sham"] < 0 else "AT2 programme OPENS")
                why = "R19: the AT2 arm clears both nulls"
            else:
                reading, why = "not established", "R19: the AT2 arm clears at most one null"
            rows.append({"well": w["well"], "role": role, "n_labelled": m, "reading": reading,
                         "why": why, "at2_relative_pct": at2["relative_pct"],
                         "at2_ci_lo": at2["relative_lo_pct"], "at2_ci_hi": at2["relative_hi_pct"],
                         "at2_geneset_percentile": at2["geneset_percentile"]})
            print("   ->", reading, "|", why)

    res, arms = pd.DataFrame(rows), pd.DataFrame(arm_rows)
    res.to_csv(OUT / "m3_readings.csv", index=False)
    arms.to_csv(OUT / "m3_arms.csv", index=False)
    rec.add_output(OUT / "m3_readings.csv"); rec.add_output(OUT / "m3_arms.csv")
    rec.set("readings", res.to_dict("records"))
    rec.finish()
    print("\n" + res.to_string(index=False))
    print("\n" + arms.drop(columns=["gene_list"]).to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
