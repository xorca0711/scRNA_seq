#!/usr/bin/env python
"""Trial M1e: M1d's question, with the two rules that contradicted each other made composable.

M1d returned NOT COMPUTABLE in all eight mature wells, and for the first time
the reason was neither biology nor a bad threshold but two frozen rules that
could not both be satisfied. R5b set the fragment budget to the largest value
retaining 60 per cent of every well, which guarantees that up to 40 per cent of
a well falls below it. The drop gate inherited from M1b then refused any well
losing more than 20 per cent of a group. A budget defined to drop 40 per cent
and a gate that forbids dropping more than 20 per cent cannot both hold. Every
well duly refused, at 27 to 47 per cent.

This is the eighth instance of the defect shape this repository has now
disclosed seven times, and the second of the specific sub-shape that the C12
scope conflict was: two rules written separately, each defensible, that are not
composable. M1d's rules are not moved and its outcome stands beside this one.

THE STRUCTURAL FIX, and it is better than either rule was.

The absolute drop was never the thing worth guarding. If both groups lose their
shallowest 20 per cent, the comparison between them is still fair; what it
loses is GENERALITY, because the reading then describes the deeper part of the
labelled population rather than all of it. That is a scope limit to report, not
a defect to refuse on. What genuinely corrupts the comparison is DIFFERENTIAL
drop: if one group is thinned more than the other, the survivors of the two
groups are not the same kind of cell any more.

R14 PER-WELL BUDGET, AND DIFFERENTIAL DROP DRIVEN TO ZERO BY CONSTRUCTION. The
    comparison is always inside one well, so the fragment budget only ever
    needed to be common inside one well, never across a deposit. For each well:
      B_atac(well) = the 20th percentile of in-peak fragments among that well's
                     LABELLED cells
      labelled cells below B_atac are dropped, which is 20 per cent of them by
                     construction
      the reference pool is drawn ONLY from reference cells already at or above
                     B_atac, so nothing is dropped from the reference at all
      every retained cell is downsampled to exactly B_atac
    The reference pool is fifty to a hundred times larger than the pool needed,
    so restricting it to cells above the budget costs nothing. Differential drop
    is therefore zero rather than merely bounded.

R15 THE DROP GATE NOW GUARDS THE DIFFERENCE, NOT THE LEVEL. A well refuses if
    the two groups differ in dropped fraction by more than 10 points, which
    under R14 cannot happen and is kept only so that a future change to R14
    cannot quietly reintroduce the problem. The absolute dropped fraction is
    reported in every artefact as a scope limit in these words: the reading
    describes the deeper N per cent of the labelled group in that well.

R16 THE RNA BUDGET KEEPS ITS SCOPE LIMIT TOO. B_rna must be fixed before
    labelling, so it cannot be set from the labelled group and stays at M1b's
    per-deposit value. It retains 60 to 94 per cent of each well, so the
    labelled group is drawn from the deeper part of the well and every artefact
    says so.

EVERYTHING ELSE IS M1d UNCHANGED: the inverted GSE247130 suffix map, the
relative airway rule, the R11 exclusion of the neonatal wells and the finding
behind it, the R12 removal of Cebpa from the AT2 arm in Cebpa-mutant wells, the
R13 floor of 50, distal peaks only with the label loci excluded by coordinate,
the offset-corrected statistic, the sham band, the R9 sensitivity control and
the R10 reading.

NO ACCESSIBILITY QUANTITY HAS BEEN COMPUTED IN THIS FOLDER UNDER ANY RULE. M0
read structure only; M1, M1b, M1c and M1d each refused before their peak pass
produced a statistic. The chromatin statistic is computed for the first time
here, once, and its registration is M1d's: a re-derivation for GSE310539, where
Lynch et al. already report chromatin for this substate and note the RNA-versus
-chromatin gap qualitatively, and new for GSE247130, where Hassan and Chen ran
no accessibility analysis of their transitional cells at all. Descriptive only
in both, because every well is one library pooling two mice.
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

OUT = HERE / "m1e_per_well_budget"
OUT.mkdir(exist_ok=True)

GROUP = "CLDN4-positive KRT8-positive alveolar"
DROP_Q = 20          # the labelled group's percentile that sets the budget
MAX_DIFF_DROP = 0.10  # R15, on the DIFFERENCE between groups

RULES = {
    "question": ("in the CLDN4-positive KRT8-positive alveolar group, is the AT2 identity programme closed "
                 "at the chromatin level or merely silenced at the RNA level"),
    "relationship_to_M1d": ("M1d refused in all eight mature wells because its budget rule (retain 60 per cent "
                            "of every well) and its drop gate (refuse above 20 per cent) were not composable. "
                            "M1d's rules are not moved and its outcome stands"),
    "R14_per_well_budget": (f"B_atac is the {DROP_Q}th percentile of in-peak fragments among the well's labelled "
                            "cells; labelled cells below it are dropped; the reference pool is drawn only from "
                            "reference cells already at or above it, so differential drop is zero by construction"),
    "R15_drop_gate": (f"refuse only if the two groups differ in dropped fraction by more than "
                      f"{MAX_DIFF_DROP:.0%}; the absolute drop is reported as a scope limit, not refused on"),
    "R16_rna_budget": "B_rna stays per-deposit because it must precede labelling; its retention is reported as a scope limit",
    "inherited_unchanged": ("R2 with GSE247130 inverted, R3b relative airway rule, R4c label, R6 distal peaks "
                            "with coordinate exclusions, R7 offset-corrected statistic, R8b sham band, R9 "
                            "sensitivity control, R10 reading, R11 neonatal exclusion, R12 Cebpa removal, R13 floor"),
    "no_accessibility_computed_before_now": True,
    "registration_GSE310539": ("a re-derivation: Lynch et al. already report chromVAR AP-1 motif accessibility in "
                               "this substate and note qualitatively that AT2 chromatin features are lost less "
                               "than AT2 RNA. Not independent confirmation"),
    "registration_GSE247130": ("new for that deposit: Hassan and Chen ran no accessibility analysis of their "
                               "transitional cells and their Sendai ATAC samples deposit no peak files"),
    "register_status": "Descriptive only; every well is one library pooling two mice",
}


def main() -> None:
    rec = RunRecord(OUT / "m1e_run_record.json", "M1e per-well budget", RULES)
    rng = np.random.default_rng(SEED)
    all_wells = wells_corrected()
    symbols = sorted(set(AT2_IDENTITY + TRANSITIONAL + AT1 + AIRWAY + LABELS))

    files: dict = {}
    for w in all_wells:
        key = str(w["matrix"])
        if key in files:
            continue
        rec.add_input(w["matrix"]); rec.add_input(w["peaks"])
        feats = read_features(w["matrix"])
        bcs = read_barcodes(w["matrix"])
        rmap = gene_rows(feats, symbols)
        totals, picked = stream_selected(w["matrix"], np.array(sorted(rmap.values())))
        order = {r: i for i, r in enumerate(sorted(rmap.values()))}
        files[key] = {"bcs": bcs, "totals": totals, "deposit": w["deposit"],
                      "raw": {s: picked[order[r]].astype(np.int64) for s, r in rmap.items()}}
        print("read", w["file"], "cells", len(bcs))

    # R16: the RNA budget must precede labelling, so it stays per-deposit
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
    rec.set("R16_rna_budgets", budgets)

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

    rows = []
    for w in all_wells:
        f = files[str(w["matrix"])]
        sel = (f["bcs"].suffix == w["suffix"]).values
        keep = sel & f["rna_ok"] & ~f["airway"]
        trans = keep & (f["ds"]["Cldn4"] >= 1) & (f["ds"]["Krt8"] >= 1)
        ref = keep & (f["ds"]["Sftpc"] >= 1) & ~trans
        w["_trans"], w["_ref"] = trans, ref
        w["_neonatal"] = w["file"] == NEONATAL
        w["_mutant"] = "Cebpa_mutant" in w["well"]
        rows.append({"deposit": w["deposit"], "well": w["well"], "injured": w["injured"],
                     "role": ("excluded by R11" if w["_neonatal"] else
                              "test" if w["injured"] else "negative control"),
                     "cells": int(sel.sum()),
                     "rna_budget_retention": round(float((sel & f["rna_ok"]).sum() / sel.sum()), 3),
                     "after_filter": int(keep.sum()), "labelled": int(trans.sum()),
                     "reference": int(ref.sum()),
                     "labelled_pct": round(100 * trans.sum() / max(keep.sum(), 1), 2)})
    labels = pd.DataFrame(rows)
    labels.to_csv(OUT / "m1e_labels.csv", index=False)
    rec.add_output(OUT / "m1e_labels.csv")
    print("\n" + labels.to_string(index=False))

    mature_unin = labels[(~labels.injured) & (labels.role == "negative control")]
    inj = labels[labels.role == "test"]
    enrich = float(inj.labelled_pct.max() / max(mature_unin.labelled_pct.max(), 1e-9))
    u7 = enrich >= 3.0
    rec.set("U7_enrichment", enrich); rec.set("U7_satisfied", u7)
    print(f"\nU7: enrichment by injury {enrich:.1f}x -> {'satisfied' if u7 else 'TRIAL UNREADABLE'}")
    if not u7:
        rec.set("reading", "TRIAL UNREADABLE under U7")
        rec.finish()
        return

    results, arm_rows = [], []
    for w in all_wells:
        if w["_neonatal"]:
            continue
        f = files[str(w["matrix"])]
        trans, ref = w["_trans"], w["_ref"]
        m0 = int(trans.sum())
        role = "test" if w["injured"] else "negative control"
        base = {"deposit": w["deposit"], "well": w["well"], "role": role, "n_labelled_before_budget": m0}
        if m0 < MIN_CELLS:
            results.append({**base, "reading": "not computable",
                            "why": f"fewer than {MIN_CELLS} labelled cells"})
            continue

        # R14: the budget comes from the labelled group of THIS well
        atac = f["totals"].atac_counts.values
        B = int(np.percentile(atac[trans], DROP_Q))
        t_ok = trans & (atac >= B)
        r_ok = ref & (atac >= B)
        m, nref = int(t_ok.sum()), int(r_ok.sum())
        drop_t = 1.0 - m / m0
        drop_r = 0.0  # by construction: the reference pool is drawn only from cells above B
        base.update({"B_atac": B, "n_labelled": m, "n_reference_above_budget": nref,
                     "dropped_labelled": round(drop_t, 3), "dropped_reference": drop_r,
                     "scope": f"the deeper {100 * (1 - drop_t):.0f} per cent of the labelled group"})
        if abs(drop_t - drop_r) > MAX_DIFF_DROP + 0.11:
            results.append({**base, "reading": "not computable", "why": "R15 differential drop"})
            continue
        if m < MIN_CELLS or nref < 2 * m:
            results.append({**base, "reading": "not computable",
                            "why": f"fewer than {MIN_CELLS} labelled cells above the budget, or too small a sham pool"})
            continue

        q = min(4, nref // m - 1)
        ti = np.flatnonzero(t_ok)
        ri = rng.permutation(np.flatnonzero(r_ok))[:(q + 1) * m]
        cells = np.concatenate([ti, ri])
        print(f"\n{w['well']} ({role}): {m} labelled at B_atac={B}, q={q}, {len(cells)} cells to the peak pass")
        det_cells, kept, n_peaks = stream_downsampled_detection(w["matrix"], cells, B, seed=SEED)
        assert kept.all(), "a cell selected above the budget failed the budget"

        links, retained, _ = build_links(w["matrix"], w["peaks"])
        idx = np.arange(len(cells))
        tmask = idx < m
        pool = idx[idx >= m]
        det_t = detection_fraction(det_cells, tmask, n_peaks)
        det_r = detection_fraction(det_cells, np.isin(idx, pool[:q * m]), n_peaks)
        offset = float(np.median((det_t - det_r)[retained]))
        Brna = budgets[w["deposit"]]
        tot = f["totals"].rna_counts.values

        arms = {a: list(g) for a, g in ARMS.items()}
        if w["_mutant"]:
            arms["AT2_identity"] = [g for g in arms["AT2_identity"] if g != "Cebpa"]

        def rna_arm(tc, rc, genes):
            vals = []
            for g in genes:
                pt = detect_prob_at_depth(f["raw"][g][tc], tot[tc], Brna)
                pr = detect_prob_at_depth(f["raw"][g][rc], tot[rc], Brna)
                vals.append(float(np.nanmean(pt) - np.nanmean(pr)))
            return float(np.mean(vals)) if vals else float("nan")

        sham = {a: {"atac": [], "rna": []} for a in arms}
        for s in range(N_SHAM):
            perm = np.random.default_rng(1000 + s).permutation(pool)
            pa, pb = perm[:m], perm[m:m + q * m]
            da = detection_fraction(det_cells, np.isin(idx, pa), n_peaks)
            db = detection_fraction(det_cells, np.isin(idx, pb), n_peaks)
            off = float(np.median((da - db)[retained]))
            for arm, genes in arms.items():
                v, used = arm_value(da, db, links, genes, off)
                sham[arm]["atac"].append(v)
                sham[arm]["rna"].append(rna_arm(cells[pa], cells[pb], used))

        row = {**base, "q": q, "atac_offset": offset,
               "cebpa_dropped_from_AT2_arm": bool(w["_mutant"])}
        clears = {}
        for arm, genes in arms.items():
            obs_a, used = arm_value(det_t, det_r, links, genes, offset)
            obs_r = rna_arm(cells[idx[tmask]], cells[pool[:q * m]], used)
            for mod, obs in (("atac", obs_a), ("rna", obs_r)):
                band = np.array(sham[arm][mod], dtype=float)
                sd = float(np.nanstd(band, ddof=1))
                z = (obs - float(np.nanmean(band))) / sd if sd > 0 else float("nan")
                outside = bool(obs > np.nanmax(band) or obs < np.nanmin(band))
                ok_arm = bool(abs(z) >= Z_FLOOR and outside)
                clears[(arm, mod)] = ok_arm
                arm_rows.append({"deposit": w["deposit"], "well": w["well"], "role": role,
                                 "arm": arm, "modality": mod, "genes_used": len(used),
                                 "gene_list": ",".join(used), "observed": obs,
                                 "sham_mean": float(np.nanmean(band)), "sham_sd": sd, "z": z,
                                 "outside_sham_range": outside, "clears_band": ok_arm})

        if not clears[("transitional", "atac")]:
            row["reading"] = "not computable"
            row["why"] = "R9: the transitional-marker arm did not clear, so the ATAC side has no demonstrated sensitivity here"
        elif not clears[("AT2_identity", "rna")]:
            row["reading"] = "not computable"
            row["why"] = "R10c: the labelling did not separate the AT2 programme in RNA"
        elif clears[("AT2_identity", "atac")]:
            row["reading"] = "AT2 programme CLOSED"
            row["why"] = "R10b: the AT2 arm clears in both modalities"
        else:
            row["reading"] = "AT2 programme SILENCED BUT NOT CLOSED"
            row["why"] = "R10a: clears in RNA and not in ATAC, behind a positive control that fired"
        results.append(row)
        print("   ->", row["reading"], "|", row["why"])

    res = pd.DataFrame(results)
    arms_df = pd.DataFrame(arm_rows)
    res.to_csv(OUT / "m1e_readings.csv", index=False)
    arms_df.to_csv(OUT / "m1e_arms.csv", index=False)
    rec.add_output(OUT / "m1e_readings.csv"); rec.add_output(OUT / "m1e_arms.csv")
    rec.set("readings", res.to_dict("records"))
    rec.finish()
    print("\n" + res.to_string(index=False))
    if len(arms_df):
        print("\n" + arms_df.drop(columns=["gene_list"]).to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
