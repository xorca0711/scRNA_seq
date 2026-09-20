#!/usr/bin/env python
"""Trial M1b: the corrected pass on M1's question.

M1 asked whether the AT2 identity programme is CLOSED at the chromatin level
or merely SILENCED at the RNA level in the CLDN4-positive transitional state,
and returned NOT COMPUTABLE in every well because two of its frozen rules were
wrong. Those rules were not moved; M1's outcome stands beside this one. What
M1 established and what this trial inherits is its rule R2: the barcode-suffix
map of GSE247130 is INVERTED relative to the GEO sample order, refuted on
Cebpa in all three files and confirmed independently on Cldn4 and on Sox9.

THE THREE CORRECTIONS, each naming the defect it repairs.

C1 repairs M1's R3. M1 dropped a cell if it carried a non-zero count of
   Scgb1a1, Scgb3a2, Foxj1 or Krt5. In single-nucleus data an abundant
   secreted transcript is ambient in nearly every droplet: Scgb1a1 is detected
   in 73 to 100 per cent of cells in every well, so the rule removed a median
   of 99.7 per cent of cells, and what it removed tracked depth (Scgb3a2
   detection 0.950 in the deepest well against 0.028 in the shallowest). Here
   a cell is airway only if its airway markers outweigh its alveolar markers
   inside the same cell, so the rule removes cells whose identity is airway
   rather than cells that merely carry airway ambient.

C2 repairs M1's R4. M1 froze its Krt8 threshold on RAW counts while median RNA
   UMI per cell varies 2.6-fold across the ten wells, so the same number meant
   different things in different wells. Here every label is called on RNA that
   has first been downsampled to one common depth, so a detection call is not
   a depth call and the thresholds transfer between wells by construction.

C3 repairs a defect M1 did not have but its ancestor design did, and this one
   was measured rather than argued. An adversarial check built a pseudo-contrast
   from reference cells only, depth-matched, containing no transitional cell at
   all, and the original peak-counting statistic returned 3,393 peaks over a
   0.10 detection difference against a permuted null mean of 2,617 at p = 0/200
   in GSE310539, and the same in two further files, with 100 per cent of
   crossing peaks pointing at the deeper group. A statistic that returns
   certainty on data containing no biology is not a statistic. Depth is
   therefore removed by DOWNSAMPLING every cell to an identical fragment
   budget rather than by stratifying or matching, and the pseudo-contrast is
   promoted from a diagnostic to a permanent part of the reading: the sham band.

FROZEN RULES. Rules R1, R6, R7, R9 and R10 of M1 are carried over unchanged;
R2 is carried over with its own refutation applied; R3, R4, R5 and R8 are
replaced as below. Section U is carried over with U2 and U6 restated.

R3b COMPARTMENT FILTER, from RNA only, and RELATIVE rather than absolute. A
    cell is airway and is dropped if its summed Scgb1a1, Scgb3a2, Foxj1 and
    Krt5 count exceeds its summed Sftpc, Sftpb and Lamp3 count. A ratio inside
    one cell is free of that cell's depth and free of a global ambient shift,
    which an absolute cut is not: Sendai virus is an airway-tropic pathogen and
    Scgb1a1 runs about twice as high per 10,000 UMI in the infected wells as in
    their uninjured controls, so a cut frozen from an uninjured well and applied
    absolutely removes half of an infected well. U3 remains the check that this
    worked.

    CALIBRATION, DISCLOSED. This rule and the R5b budget below were both set
    after a first run of this trial that computed LABELS ONLY. That run used an
    absolute airway cut at the uninjured well's 99th percentile and a depth
    budget at the deposit's 10th percentile; it removed 49 per cent of one
    infected well as airway and left 0 to 56 transitional cells per well, so
    the labelling did not survive its own technical settings. Both settings are
    functions of the depth distribution and the compartment markers alone. NO
    ACCESSIBILITY QUANTITY WAS COMPUTED UNDER ANY SETTING before these were
    fixed, and the chromatin statistic below is computed exactly once. The
    first run's numbers are in `m1b_calibration_labels.csv` so that the record
    shows both.

R4b LABELS, from RNA only, on downsampled counts, thresholds frozen from the
    uninjured well before any injured well is opened:
      k = the 99th percentile of downsampled Krt8 among Sftpc-positive,
          non-airway cells of the uninjured well
      transitional = downsampled Cldn4 >= 1 AND downsampled Krt8 >= k
      reference    = downsampled Sftpc >= 1 AND not transitional
    The reference is named "Sftpc-positive non-transitional epithelium" and
    never "AT2". Lynch resolves four transitional substates T1 to T4 and this
    binary collapses T1, T3 and T4 into the reference, so the reference is a
    mixture and is named as one.

R5b DEPTH EQUALISATION REPLACES MATCHING, at the largest budget that keeps
    most of the data. Within a deposit, B_rna and B_atac are the LARGEST values
    such that at least 60 per cent of the cells of EVERY well of that deposit
    reach them, which is the minimum over wells of that well's 40th percentile.
    The budget is therefore a function of the depth distribution alone and
    carries no state-level information. Every cell is downsampled to exactly
    B_rna UMI and exactly B_atac in-peak fragments; a cell below either budget
    is dropped and the number dropped is reported per group. Every detection
    fraction in this trial is computed on the downsampled data and on nothing
    else.
      The deposit's 10th percentile, which the calibration run used, leaves
    560 UMI for GSE310539 against a median of about 4,700, and at that depth
    Cldn4 is essentially undetectable, so the label does not survive. Depth
    equalisation always costs sensitivity; the question is only how much, and
    60 per cent retention is the point chosen in advance of any accessibility
    quantity being computed.
      RNA is downsampled two ways for two purposes, and the reason is stated
    here rather than left to the implementation. LABELS use one seeded
    hypergeometric draw per cell, because a label thresholds an integer
    magnitude. The ARM STATISTIC uses the exact hypergeometric expectation of
    detection, because there the quantity wanted is a detection fraction and
    its expectation is available in closed form with no sampling noise. Peaks
    are downsampled by one seeded draw per cell, reused unchanged for the real
    contrast and for every sham, so the sham band measures the pipeline rather
    than a second roll of the dice.

R8b SHAM BAND. A pool of reference cells of size (q+1) x n_transitional is
    drawn, with q = min(4, floor(n_reference / n_transitional) - 1). The real
    contrast is the transitional group against q x n_transitional reference
    cells. Each of 20 shams draws n_transitional reference cells from the pool
    as a pseudo-transitional group and q x n_transitional of the remainder as
    its comparator, so the sham group sizes equal the real ones. An arm clears
    the band only if its absolute z against the 20 sham values is at least 3
    AND its value lies outside the full sham range.

U2b Because labels are called on downsampled RNA, the Cldn4 detection call is
    no longer a function of depth and M1's U2 no longer applies as written. It
    is replaced: if the transitional and reference groups differ by more than
    10 per cent in median ORIGINAL in-peak fragments AFTER the R5b drop, the
    well is reported with that ratio beside every number, because a surviving
    imbalance in the cells that cleared the budget is still an imbalance in
    which cells were read even though it is not one in how deeply they were read.

U6b Fewer than 2 x n_transitional reference cells: the sham band cannot be
    built and the well is unreadable.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from multiome_utils import (ATAC, GSE247130, GSE310539, RNA, REPO,  # noqa: E402
                            RunRecord, detect_prob_at_depth, detection_fraction,
                            df_to_markdown, peak_gene_table, read_barcodes,
                            read_features, stream_downsampled_detection,
                            stream_selected)
from m1_closed_or_merely_silenced import (AIRWAY, ARMS, AT1, AT2_IDENTITY,  # noqa: E402
                                          EXCLUDE_AROUND, EXCLUDE_KB, LABELS,
                                          N_SHAM, TRANSITIONAL, Z_FLOOR,
                                          excluded_peaks, gene_rows)

OUT = HERE / "m1b_corrected_pass"
OUT.mkdir(exist_ok=True)

MIN_CELLS = 100
MAX_DROP = 0.20
SEED = 0

# M1's rule R2 refuted the GEO sample order for GSE247130 on Cebpa in all three
# files, and Cldn4 and Sox9 agree. Suffix 1 is the control and suffix 2 is the
# Cebpa mutant, the reverse of the deposited order.
GSE247130_INVERTED = {"1": "control", "2": "Cebpa mutant"}

RULES = {
    "question": ("M1's question, on corrected rules: in the CLDN4-positive transitional state, is the "
                 "AT2 identity programme closed at the chromatin level or merely silenced at the RNA level"),
    "relationship_to_M1": "M1's thresholds were not moved; its outcome stands beside this one",
    "inherited_R2_refutation": ("GSE247130's barcode suffix map is inverted relative to the GEO sample order: "
                                "suffix 1 is the control and suffix 2 is the Cebpa mutant, refuted on Cebpa in "
                                "all three files and confirmed on Cldn4 and Sox9"),
    "C1_repairs_R3": ("M1 dropped any cell with a non-zero airway marker count, which is ambient in single-nucleus "
                      "data and removed a median of 99.7 per cent of cells with a severity that tracked depth"),
    "C2_repairs_R4": "M1 froze its threshold on raw counts while median UMI varies 2.6-fold across the ten wells",
    "C3_depth": ("the ancestor peak-counting design was measured to return p = 0/200 on a depth-matched "
                 "pseudo-contrast containing no transitional cell, so depth is removed by downsampling rather "
                 "than by stratifying, and the pseudo-contrast becomes the permanent sham band"),
    "R3b_compartment": "airway if a downsampled marker count reaches its 99th percentile among Sftpc-positive cells of the uninjured well",
    "R4b_labels": "on downsampled RNA; k = 99th percentile of downsampled Krt8 among Sftpc-positive non-airway cells of the uninjured well",
    "R5b_depth": "every cell downsampled to the deposit's 10th-percentile RNA UMI and in-peak fragment budget; cells below are dropped",
    "R6_features": f"distal peaks only for the primary; exclude every peak within {EXCLUDE_KB} kb of "
                   f"{', '.join(EXCLUDE_AROUND)} and every peak off a chr contig",
    "R7_statistic": "per gene, transitional minus reference detection, corrected by the genome-wide median distal difference",
    "R8b_sham": f"{N_SHAM} shams at matched group sizes; an arm clears only if |z| >= {Z_FLOOR} and the value is outside the full sham range",
    "R9_sensitivity": "the transitional arm is the ATAC positive control; if it does not clear, the AT2 arm is NOT COMPUTABLE",
    "R10_reading": "a=silenced not closed, b=closed, c/d=not computable",
    "not_tested": "any contrast between conditions, genotypes, stages or deposits; any motif; any causation",
}


def wells_corrected() -> list[dict]:
    out = []
    for suf, meta in GSE310539["libraries"].items():
        out.append({"deposit": "GSE310539", "file": "totalaggr", "suffix": suf,
                    "matrix": GSE310539["matrix"], "peaks": GSE310539["peaks"],
                    "well": meta["name"], "injured": meta["treatment"] == "SeV",
                    "threshold_source": meta["name"] == "wildtype_PBS"})
    for f in GSE247130["files"]:
        for suf in ("1", "2"):
            geno = GSE247130_INVERTED[suf]
            name = f"{f['stage']}_{'control' if geno == 'control' else 'Cebpa_mutant'}"
            out.append({"deposit": "GSE247130", "file": f["stage"], "suffix": suf,
                        "matrix": f["matrix"], "peaks": f["peaks"],
                        "well": name, "injured": f["stage"] == "SeV",
                        "threshold_source": name == "7wk_control"})
    return out


def arm_value(det_t, det_r, links, genes, offset):
    vals, used = [], []
    for g in genes:
        rows = links.get(g)
        if rows is None or len(rows) == 0:
            continue
        vals.append(float(np.mean(det_t[rows] - det_r[rows])) - offset)
        used.append(g)
    return (float(np.mean(vals)) if len(used) >= 4 else float("nan")), used


def build_links(matrix: Path, peaks: Path) -> tuple[dict, np.ndarray, dict]:
    """Gene symbol -> its DISTAL peak slots, plus the retained distal peak mask."""
    feats = read_features(matrix)
    excl, detail = excluded_peaks(feats)
    peak_feats = feats[feats.feature_type == ATAC].reset_index(drop=True)
    slot = {iv: i for i, iv in enumerate(peak_feats.interval)}
    rna = feats[feats.feature_type == RNA]
    sym2ens = dict(zip(rna.name, rna.id))

    ann = peak_gene_table(peaks)
    ann = ann[(ann.peak_type == "distal") & (~ann.interval.isin(excl))]
    ens2sym = {sym2ens[s]: s for arm in ARMS.values() for s in arm if s in sym2ens}
    links: dict[str, np.ndarray] = {}
    for ens, sym in ens2sym.items():
        rows = ann.loc[ann.gene == ens, "interval"]
        links[sym] = np.array(sorted({slot[i] for i in rows if i in slot}), dtype=np.int64)

    keep_iv = set(ann.interval)
    retained = np.zeros(len(peak_feats), dtype=bool)
    for iv, i in slot.items():
        if iv in keep_iv:
            retained[i] = True
    return links, retained, detail


def main() -> None:
    rec = RunRecord(OUT / "m1b_run_record.json", "M1b corrected pass", RULES)
    rng = np.random.default_rng(SEED)
    all_wells = wells_corrected()
    symbols = sorted(set(AT2_IDENTITY + TRANSITIONAL + AT1 + AIRWAY + LABELS))

    files: dict = {}
    for w in all_wells:
        key = str(w["matrix"])
        if key in files:
            continue
        rec.add_input(w["matrix"])
        rec.add_input(w["peaks"])
        feats = read_features(w["matrix"])
        bcs = read_barcodes(w["matrix"])
        rmap = gene_rows(feats, symbols)
        totals, picked = stream_selected(w["matrix"], np.array(sorted(rmap.values())))
        order = {r: i for i, r in enumerate(sorted(rmap.values()))}
        raw = {s: picked[order[r]].astype(np.int64) for s, r in rmap.items()}
        files[key] = {"bcs": bcs, "totals": totals, "raw": raw, "deposit": w["deposit"]}
        print("read", w["file"], "cells", len(bcs))

    # ---- R5b budgets: the largest that keeps 60 per cent of every well ------
    budgets = {}
    for dep in ("GSE310539", "GSE247130"):
        rna_q, atac_q = [], []
        for w in all_wells:
            if w["deposit"] != dep:
                continue
            f = files[str(w["matrix"])]
            sel = (f["bcs"].suffix == w["suffix"]).values
            rna_q.append(np.percentile(f["totals"].rna_counts.values[sel], 40))
            atac_q.append(np.percentile(f["totals"].atac_counts.values[sel], 40))
        budgets[dep] = {"B_rna": int(min(rna_q)), "B_atac": int(min(atac_q))}
    rec.set("R5b_budgets", budgets)
    print("\nR5b budgets (largest keeping 60 per cent of every well):", budgets)

    # ---- downsampled RNA counts for labelling, and the R3b airway rule -----
    for f in files.values():
        B = budgets[f["deposit"]]["B_rna"]
        tot = f["totals"].rna_counts.values
        ok = tot >= B
        f["rna_ok"] = ok
        f["ds"] = {}
        for s, c in f["raw"].items():
            out = np.zeros(len(c), dtype=np.int64)
            good = np.clip(c[ok], 0, tot[ok])
            out[ok] = rng.hypergeometric(np.maximum(good, 0),
                                         np.maximum(tot[ok] - good, 0), B)
            f["ds"][s] = out
        air = sum(f["raw"][a] for a in AIRWAY if a in f["raw"])
        alv = sum(f["raw"][a] for a in ("Sftpc", "Sftpb", "Lamp3") if a in f["raw"])
        f["airway"] = air > alv

    # ---- R4b: k frozen from the uninjured wells -----------------------------
    thr = {}
    for w in all_wells:
        if not w["threshold_source"]:
            continue
        f = files[str(w["matrix"])]
        sel = (f["bcs"].suffix == w["suffix"]).values & f["rna_ok"] & ~f["airway"]
        sftpc = (f["ds"]["Sftpc"] > 0) & sel
        thr[w["deposit"]] = {"krt8_k": int(np.percentile(f["ds"]["Krt8"][sftpc], 99)),
                             "n_clean": int(sftpc.sum()), "source_well": w["well"]}
    rec.set("thresholds_from_uninjured_wells", thr)
    print("thresholds:", {d: (t["source_well"], t["krt8_k"]) for d, t in thr.items()})

    rows = []
    for w in all_wells:
        f = files[str(w["matrix"])]
        sel = (f["bcs"].suffix == w["suffix"]).values
        ok = sel & f["rna_ok"]
        keep = ok & ~f["airway"]
        k = max(thr[w["deposit"]]["krt8_k"], 1)
        trans = keep & (f["ds"]["Cldn4"] >= 1) & (f["ds"]["Krt8"] >= k)
        ref = keep & (f["ds"]["Sftpc"] >= 1) & ~trans
        w["_trans"], w["_ref"] = trans, ref
        u3 = {}
        for a in ("Scgb1a1", "Scgb3a2", "Foxj1"):
            if a in f["ds"] and trans.sum() and ref.sum():
                u3[a] = round(float((f["ds"][a][trans] > 0).mean()
                                    - (f["ds"][a][ref] > 0).mean()), 3)
        rows.append({
            "deposit": w["deposit"], "file": w["file"], "suffix": w["suffix"],
            "well": w["well"], "injured": w["injured"], "cells": int(sel.sum()),
            "kept_by_budget": round(ok.sum() / max(sel.sum(), 1), 3),
            "airway_dropped": int((ok & f["airway"]).sum()),
            "after_filter": int(keep.sum()),
            "transitional": int(trans.sum()), "reference": int(ref.sum()),
            "transitional_pct": round(100 * trans.sum() / max(keep.sum(), 1), 2),
            "U3_max_airway_excess": max(u3.values()) if u3 else float("nan"),
        })
    labels = pd.DataFrame(rows)
    labels.to_csv(OUT / "m1b_labels.csv", index=False)
    rec.add_output(OUT / "m1b_labels.csv")
    print("\n" + labels.to_string(index=False))

    unin = labels[~labels.injured]
    worst = int(unin.transitional.max())
    u7 = worst < MIN_CELLS
    rec.set("U7_max_transitional_uninjured", worst)
    rec.set("U7_satisfied", bool(u7))
    print(f"\nU7: uninjured wells reach at most {worst} transitional cells "
          f"(floor {MIN_CELLS}) -> {'satisfied' if u7 else 'TRIAL UNREADABLE'}")
    if not u7:
        rec.set("reading", "TRIAL UNREADABLE under U7")
        rec.finish()
        return

    # ---- the chromatin statistic, computed once -----------------------------
    results, arm_rows = [], []
    for w in all_wells:
        trans, ref = w["_trans"], w["_ref"]
        m, nref = int(trans.sum()), int(ref.sum())
        base = {"deposit": w["deposit"], "well": w["well"], "injured": w["injured"],
                "n_transitional": m, "n_reference": nref}
        if m < MIN_CELLS:
            results.append({**base, "reading": "not computable",
                            "why": f"fewer than {MIN_CELLS} transitional cells"})
            continue
        if nref < 2 * m:
            results.append({**base, "reading": "not computable",
                            "why": "U6b, the sham pool is smaller than twice the transitional group"})
            continue

        q = min(4, nref // m - 1)
        f = files[str(w["matrix"])]
        ti = np.flatnonzero(trans)
        ri = rng.permutation(np.flatnonzero(ref))[:(q + 1) * m]
        cells = np.concatenate([ti, ri])
        B = budgets[w["deposit"]]["B_atac"]
        print(f"\n{w['well']}: {m} transitional, q={q}, {len(cells)} cells to the peak pass")
        det_cells, kept, n_peaks = stream_downsampled_detection(w["matrix"], cells, B, seed=SEED)

        drop_t = 1.0 - float(kept[:m].mean())
        drop_r = 1.0 - float(kept[m:].mean())
        if max(drop_t, drop_r) > MAX_DROP:
            results.append({**base, "reading": "not computable",
                            "why": f"R5b dropped {max(drop_t, drop_r):.0%} of a group at the ATAC budget"})
            continue

        links, retained, excl_detail = build_links(w["matrix"], w["peaks"])
        idx = np.arange(len(cells))
        tmask = (idx < m) & kept
        pool = idx[(idx >= m) & kept]
        rmask = np.isin(idx, pool[:q * m])
        det_t = detection_fraction(det_cells, tmask, n_peaks)
        det_r = detection_fraction(det_cells, rmask, n_peaks)
        offset = float(np.median((det_t - det_r)[retained]))

        Brna = budgets[w["deposit"]]["B_rna"]
        tot = f["totals"].rna_counts.values

        def rna_arm(tcells, rcells, genes):
            vals = []
            for g in genes:
                pt = detect_prob_at_depth(f["raw"][g][tcells], tot[tcells], Brna)
                pr = detect_prob_at_depth(f["raw"][g][rcells], tot[rcells], Brna)
                vals.append(float(np.nanmean(pt) - np.nanmean(pr)))
            return float(np.mean(vals)) if vals else float("nan")

        sham = {a: {"atac": [], "rna": []} for a in ARMS}
        for s in range(N_SHAM):
            perm = np.random.default_rng(1000 + s).permutation(pool)
            pa, pb = perm[:m], perm[m:m + q * m]
            da = detection_fraction(det_cells, np.isin(idx, pa), n_peaks)
            db = detection_fraction(det_cells, np.isin(idx, pb), n_peaks)
            off = float(np.median((da - db)[retained]))
            for arm, genes in ARMS.items():
                v, used = arm_value(da, db, links, genes, off)
                sham[arm]["atac"].append(v)
                sham[arm]["rna"].append(rna_arm(cells[pa], cells[pb], used))

        row = {**base, "q": q, "atac_offset": offset,
               "dropped_transitional": round(drop_t, 3),
               "dropped_reference": round(drop_r, 3)}
        clears = {}
        for arm, genes in ARMS.items():
            obs_a, used = arm_value(det_t, det_r, links, genes, offset)
            obs_r = rna_arm(cells[idx[tmask]], cells[pool[:q * m]], used)
            for mod, obs in (("atac", obs_a), ("rna", obs_r)):
                band = np.array(sham[arm][mod], dtype=float)
                sd = float(np.nanstd(band, ddof=1))
                z = (obs - float(np.nanmean(band))) / sd if sd > 0 else float("nan")
                outside = bool(obs > np.nanmax(band) or obs < np.nanmin(band))
                ok_arm = bool(abs(z) >= Z_FLOOR and outside)
                clears[(arm, mod)] = ok_arm
                arm_rows.append({**base, "arm": arm, "modality": mod,
                                 "genes_used": len(used), "observed": obs,
                                 "sham_mean": float(np.nanmean(band)), "sham_sd": sd,
                                 "z": z, "outside_sham_range": outside,
                                 "clears_band": ok_arm})

        if not clears[("transitional", "atac")]:
            row["reading"] = "not computable"
            row["why"] = ("R9: the transitional arm did not clear the band, so the ATAC side "
                          "has no demonstrated sensitivity in this well")
        elif not clears[("AT2_identity", "rna")]:
            row["reading"] = "not computable"
            row["why"] = "R10c: the labelling did not separate the AT2 programme in RNA"
        elif clears[("AT2_identity", "atac")]:
            row["reading"] = "AT2 programme CLOSED"
            row["why"] = "R10b: the AT2 arm clears the band in both modalities"
        else:
            row["reading"] = "AT2 programme SILENCED BUT NOT CLOSED"
            row["why"] = ("R10a: the AT2 arm clears in RNA and not in ATAC, behind a positive "
                          "control that fired")
        results.append(row)
        print("   ->", row["reading"], "|", row["why"])

    res = pd.DataFrame(results)
    arms_df = pd.DataFrame(arm_rows)
    res.to_csv(OUT / "m1b_readings.csv", index=False)
    arms_df.to_csv(OUT / "m1b_arms.csv", index=False)
    rec.add_output(OUT / "m1b_readings.csv")
    rec.add_output(OUT / "m1b_arms.csv")
    rec.set("readings", res.to_dict("records"))
    rec.finish()
    print("\n" + res.to_string(index=False))
    if len(arms_df):
        print("\n" + arms_df.to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
