#!/usr/bin/env python
"""Trial M1c: M1b's question, with the label set to what this data type supports.

M1b returned NOT COMPUTABLE in all ten wells, and the reason is arithmetic
rather than biological. Its rule R4b required Cldn4 detection AND Krt8 at or
above the 99th percentile of the uninjured well. A 99th-percentile magnitude
cut passes one per cent of the control well BY CONSTRUCTION, and Cldn4 is
detected in five to nine per cent of cells in an infected well, so their
conjunction cannot exceed a few tenths of a per cent anywhere. It returned 0.03
to 0.66 per cent, against the twelve per cent transitional fraction the source
paper reports for its infected wildtype library. M1b's threshold is not moved
and its outcome stands beside this one.

WHAT CHANGED, AND WHY THIS IS NOT THRESHOLD SHOPPING. No accessibility quantity
has been computed in this folder under any rule: M1b refused at the cell-count
gate, before its peak pass ran, in every well. The chromatin statistic is still
computed exactly once, here. What the three refusals established is a fact
about the data type, and it is the honest scope limit of this branch:

  THE PAPERS DEFINE THIS STATE AS CLDN4-POSITIVE AND KRT8-HIGH, AND THEY DO IT
  BY IMMUNOSTAINING. Three-prime single-nucleus counting cannot make a
  "KRT8-high" call of that kind. Krt8 is detected in 30 to 59 per cent of cells
  with a median count of 0 or 1, so a magnitude cut on it is a cut on a
  zero-inflated quantity, which is the defect this repository disclosed at row
  C56. Hassan and Chen state that baseline KRT8 is present in ALL AT2 cells, so
  the discriminating marker in their hands is ectopic CLDN4, not KRT8; and
  Lynch report KRT8-high CLDN4-negative cells in compensatory-growth regions,
  so KRT8-high is not sufficient either.

R4c LABELS, from RNA only, on RNA downsampled to a common depth:
      transitional = downsampled Cldn4 >= 1 AND downsampled Krt8 >= 1
      reference    = downsampled Sftpc >= 1 AND not transitional
    Cldn4 carries the discrimination and Krt8 is required only to be present,
    because that is the strongest form of the papers' definition this assay can
    express. The group is therefore named, in every artefact, the
    CLDN4-POSITIVE KRT8-POSITIVE ALVEOLAR GROUP. It is the transcript-level
    shadow of the CLDN4-positive KRT8-high state of Lynch et al. and of the
    KRT8/CLDN4-positive state of Hassan and Chen, and it is not the same
    object as either. No sentence produced by this trial may call it DATP, T2,
    PATS or ADI.

    The floor it has to clear is its own false-positive rate, and that is what
    U7 measures: Cldn4 is detected in about one per cent of uninjured cells and
    five to nine per cent of infected ones, so the label is enriched about
    eightfold by injury and the uninjured wells remain the negative control.

EVERYTHING ELSE IS M1b UNCHANGED: the inverted GSE247130 suffix map that M1's
rule R2 established, the relative airway rule R3b, the depth budgets R5b, the
distal-only feature rule R6 with its coordinate exclusions, the offset-corrected
statistic R7, the sham band R8b, the sensitivity control R9 and the reading R10.
Section U is unchanged. The labelling genes Cldn4 and Krt8 appear in no
accessibility arm and their loci are excluded by coordinate, so the
transitional arm remains independent of the labelling and R9 keeps its force.
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
                                          EXCLUDE_AROUND, EXCLUDE_KB, LABELS,
                                          N_SHAM, TRANSITIONAL, Z_FLOOR, gene_rows)
from m1b_corrected_pass import (MAX_DROP, MIN_CELLS, SEED, arm_value,  # noqa: E402
                                build_links, wells_corrected)

OUT = HERE / "m1c_label_at_the_depth_available"
OUT.mkdir(exist_ok=True)

GROUP_NAME = "CLDN4-positive KRT8-positive alveolar"

RULES = {
    "question": ("M1b's question with the label set to what three-prime single-nucleus counting supports: "
                 "in the CLDN4-positive KRT8-positive alveolar group, is the AT2 identity programme closed "
                 "at the chromatin level or merely silenced at the RNA level"),
    "relationship_to_M1b": ("M1b refused in all ten wells at the cell-count gate, before its peak pass ran, "
                            "so no accessibility quantity has been computed under any rule and the chromatin "
                            "statistic is still computed exactly once, here. M1b's threshold is not moved"),
    "why_R4b_could_not_work": ("a 99th-percentile magnitude cut passes one per cent of the control well by "
                               "construction, and Cldn4 is detected in five to nine per cent of an infected "
                               "well, so their conjunction cannot exceed a few tenths of a per cent"),
    "R4c_labels": ("transitional = downsampled Cldn4 >= 1 and downsampled Krt8 >= 1; reference = downsampled "
                   "Sftpc >= 1 and not transitional"),
    "scope_limit": ("the group is named the CLDN4-positive KRT8-positive alveolar group in every artefact. It is "
                    "the transcript-level shadow of the CLDN4-positive KRT8-high state of Lynch and the "
                    "KRT8/CLDN4-positive state of Hassan and Chen, both of which are immunostaining definitions, "
                    "and it is not the same object as either. It is never called DATP, T2, PATS or ADI"),
    "inherited_unchanged": ("R2 with GSE247130 inverted, R3b relative airway rule, R5b depth budgets, R6 distal "
                            "peaks with coordinate exclusions, R7 offset-corrected statistic, R8b sham band, "
                            "R9 sensitivity control, R10 reading, section U"),
    "not_tested": "any contrast between conditions, genotypes, stages or deposits; any motif; any causation",
}


def main() -> None:
    rec = RunRecord(OUT / "m1c_run_record.json", "M1c label at the depth available", RULES)
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

    budgets = {}
    for dep in ("GSE310539", "GSE247130"):
        rq, aq = [], []
        for w in all_wells:
            if w["deposit"] != dep:
                continue
            f = files[str(w["matrix"])]
            sel = (f["bcs"].suffix == w["suffix"]).values
            rq.append(np.percentile(f["totals"].rna_counts.values[sel], 40))
            aq.append(np.percentile(f["totals"].atac_counts.values[sel], 40))
        budgets[dep] = {"B_rna": int(min(rq)), "B_atac": int(min(aq))}
    rec.set("R5b_budgets", budgets)
    print("\nbudgets:", budgets)

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

    rows = []
    for w in all_wells:
        f = files[str(w["matrix"])]
        sel = (f["bcs"].suffix == w["suffix"]).values
        keep = sel & f["rna_ok"] & ~f["airway"]
        trans = keep & (f["ds"]["Cldn4"] >= 1) & (f["ds"]["Krt8"] >= 1)
        ref = keep & (f["ds"]["Sftpc"] >= 1) & ~trans
        w["_trans"], w["_ref"] = trans, ref
        u3 = {a: round(float((f["ds"][a][trans] > 0).mean() - (f["ds"][a][ref] > 0).mean()), 3)
              for a in ("Scgb1a1", "Scgb3a2", "Foxj1") if a in f["ds"] and trans.sum() and ref.sum()}
        rows.append({"deposit": w["deposit"], "file": w["file"], "suffix": w["suffix"],
                     "well": w["well"], "injured": w["injured"], "cells": int(sel.sum()),
                     "after_filter": int(keep.sum()),
                     "transitional": int(trans.sum()), "reference": int(ref.sum()),
                     "transitional_pct": round(100 * trans.sum() / max(keep.sum(), 1), 2),
                     "U3_max_airway_excess": max(u3.values()) if u3 else float("nan")})
    labels = pd.DataFrame(rows)
    labels.to_csv(OUT / "m1c_labels.csv", index=False)
    rec.add_output(OUT / "m1c_labels.csv")
    print("\n" + labels.to_string(index=False))

    unin, inj = labels[~labels.injured], labels[labels.injured]
    worst = int(unin.transitional_pct.max())
    enrich = float(inj.transitional_pct.max() / max(unin.transitional_pct.max(), 1e-9))
    u7 = bool(int(unin.transitional.max()) < MIN_CELLS) or enrich >= 3.0
    rec.set("U7_uninjured_max_pct", float(unin.transitional_pct.max()))
    rec.set("U7_injured_max_pct", float(inj.transitional_pct.max()))
    rec.set("U7_enrichment_by_injury", enrich)
    rec.set("U7_satisfied", u7)
    print(f"\nU7: uninjured at most {unin.transitional_pct.max():.2f} per cent, injured up to "
          f"{inj.transitional_pct.max():.2f} per cent, enrichment {enrich:.1f}x -> "
          f"{'satisfied' if u7 else 'TRIAL UNREADABLE'}")
    if not u7:
        rec.set("reading", "TRIAL UNREADABLE under U7")
        rec.finish()
        return

    results, arm_rows = [], []
    for w in all_wells:
        trans, ref = w["_trans"], w["_ref"]
        m, nref = int(trans.sum()), int(ref.sum())
        base = {"deposit": w["deposit"], "well": w["well"], "injured": w["injured"],
                "n_transitional": m, "n_reference": nref}
        if m < MIN_CELLS:
            results.append({**base, "reading": "not computable",
                            "why": f"fewer than {MIN_CELLS} cells in the {GROUP_NAME} group"})
            continue
        if nref < 2 * m:
            results.append({**base, "reading": "not computable", "why": "U6b, sham pool too small"})
            continue

        q = min(4, nref // m - 1)
        f = files[str(w["matrix"])]
        ti = np.flatnonzero(trans)
        ri = rng.permutation(np.flatnonzero(ref))[:(q + 1) * m]
        cells = np.concatenate([ti, ri])
        B = budgets[w["deposit"]]["B_atac"]
        print(f"\n{w['well']}: {m} labelled, q={q}, {len(cells)} cells to the peak pass")
        det_cells, kept, n_peaks = stream_downsampled_detection(w["matrix"], cells, B, seed=SEED)

        drop_t = 1.0 - float(kept[:m].mean())
        drop_r = 1.0 - float(kept[m:].mean())
        if max(drop_t, drop_r) > MAX_DROP:
            results.append({**base, "reading": "not computable",
                            "why": f"R5b dropped {max(drop_t, drop_r):.0%} of a group at the ATAC budget"})
            continue

        links, retained, _ = build_links(w["matrix"], w["peaks"])
        idx = np.arange(len(cells))
        tmask = (idx < m) & kept
        pool = idx[(idx >= m) & kept]
        rmask = np.isin(idx, pool[:q * m])
        det_t = detection_fraction(det_cells, tmask, n_peaks)
        det_r = detection_fraction(det_cells, rmask, n_peaks)
        offset = float(np.median((det_t - det_r)[retained]))
        Brna = budgets[w["deposit"]]["B_rna"]
        tot = f["totals"].rna_counts.values

        def rna_arm(tc, rc, genes):
            vals = []
            for g in genes:
                pt = detect_prob_at_depth(f["raw"][g][tc], tot[tc], Brna)
                pr = detect_prob_at_depth(f["raw"][g][rc], tot[rc], Brna)
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
               "dropped_transitional": round(drop_t, 3), "dropped_reference": round(drop_r, 3)}
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
                arm_rows.append({**base, "arm": arm, "modality": mod, "genes_used": len(used),
                                 "gene_list": ",".join(used), "observed": obs,
                                 "sham_mean": float(np.nanmean(band)), "sham_sd": sd, "z": z,
                                 "outside_sham_range": outside, "clears_band": ok_arm})

        if not clears[("transitional", "atac")]:
            row["reading"] = "not computable"
            row["why"] = ("R9: the transitional-marker arm did not clear the band, so the ATAC side has "
                          "no demonstrated sensitivity in this well")
        elif not clears[("AT2_identity", "rna")]:
            row["reading"] = "not computable"
            row["why"] = "R10c: the labelling did not separate the AT2 programme in RNA"
        elif clears[("AT2_identity", "atac")]:
            row["reading"] = "AT2 programme CLOSED"
            row["why"] = "R10b: the AT2 arm clears the band in both modalities"
        else:
            row["reading"] = "AT2 programme SILENCED BUT NOT CLOSED"
            row["why"] = "R10a: clears in RNA and not in ATAC, behind a positive control that fired"
        results.append(row)
        print("   ->", row["reading"], "|", row["why"])

    res = pd.DataFrame(results)
    arms_df = pd.DataFrame(arm_rows)
    res.to_csv(OUT / "m1c_readings.csv", index=False)
    arms_df.to_csv(OUT / "m1c_arms.csv", index=False)
    rec.add_output(OUT / "m1c_readings.csv")
    rec.add_output(OUT / "m1c_arms.csv")
    rec.set("readings", res.to_dict("records"))
    rec.finish()
    print("\n" + res.to_string(index=False))
    if len(arms_df):
        print("\n" + arms_df.drop(columns=["gene_list"]).to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
