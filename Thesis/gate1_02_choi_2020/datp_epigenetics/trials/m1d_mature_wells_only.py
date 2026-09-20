#!/usr/bin/env python
"""Trial M1d: M1c's question, with the neonatal wells taken out of the design.

M1c returned TRIAL UNREADABLE under its own rule U7, and the reason is worth
more than the trial was. U7 required the uninjured wells to hold less of the
labelled group than the injured ones. They did not: the neonatal P9 wells
labelled 3.69 and 8.07 per cent, above every injured well, so the whole-trial
negative control fired and M1c refused. M1c's threshold is not moved and its
outcome stands beside this one.

THE REASON IS A RESULT, AND IT IS ONE OF THE TWO THINGS THIS BRANCH WAS BUILT
TO FIND OUT. The owner's first question was whether the transition state has a
character specific to regeneration or disease. At the transcript level, in
these deposits, the answer is no for the markers that define it: Krt8 and Cldn4
are expressed across immature postnatal alveolar epithelium, so a
CLDN4-positive KRT8-positive call cannot separate neonatal developmental
immaturity from injury-induced transition. Hassan and Chen's own framing
predicts this and never states it as a limitation of the marker: their whole
argument is that neonatal AT2 cells are plastic and that Cebpa deletion returns
mature cells toward that neonatal state. A marker set shared with normal
development is not a damage-associated marker set.

FROZEN RULES, all of M1c unchanged except the three below.

R11 STAGE RULE. The two P9 neonatal wells are removed from the test set AND
    from the negative-control set, declared here in advance with the reason
    above. They are still labelled and still reported, because their label
    fraction is the finding of the preceding paragraph, but no reading is taken
    from them and they neither support nor refute anything. The negative
    control becomes the four MATURE uninjured wells: wildtype_PBS and AP1mut_PBS
    of GSE310539, and 7wk_control and 7wk_Cebpa_mutant of GSE247130.

R12 GENOTYPE RULE. Cebpa is a member of the AT2 identity arm and is genetically
    deleted in the Cebpa-mutant wells of GSE247130, where Hassan and Chen report
    5,287 normally AT2-specific peaks closing on that deletion. In those wells
    Cebpa is dropped from the AT2 arm and named in the output, and the arm must
    still reach four genes or it is not computable. Without this the AT2 arm in
    a mutant well would measure the knockout rather than the transition.

R13 CELL FLOOR. The floor moves from 100 cells to 50, and the reason is that
    the old floor was inherited from a rule this design no longer contains.
    M1's U4 set 100 because a detection fraction estimated on n cells has
    granularity 1/n and M1's statistic thresholded a single peak at 0.10. This
    design has no per-peak threshold: its statistic is an arm mean over many
    peaks, read against a sham band computed AT THE SAME GROUP SIZES. The sham
    band is therefore the correct small-n guard and it subsumes the old floor,
    because a group too small to measure widens the band until nothing clears
    it. The floor that remains is only what the sham pool needs to exist.

NO ACCESSIBILITY QUANTITY HAS BEEN COMPUTED IN THIS FOLDER UNDER ANY RULE. M0
read structure only; M1, M1b and M1c each refused before their peak pass ran.
The chromatin statistic is computed for the first time here, once.

REGISTRATION OF WHAT A POSITIVE RESULT WOULD MEAN, fixed before it is seen, and
it is not the same for the two deposits.

  GSE310539. Lynch et al. already report chromatin for this substate: chromVAR
  AP-1 motif accessibility spikes in their CLDN4-positive T2, and they note in
  passing that their transitional substates lose AT2 chromatin features less
  dramatically than they lose AT2 RNA. A result here that agrees is a
  QUANTIFICATION OF A PUBLISHED QUALITATIVE OBSERVATION and is registered as a
  re-derivation, not as independent confirmation. What is new is only that the
  observation is gene-anchored, measured against a sham band, and made for the
  CLDN4-positive group specifically rather than for T1, T3 and T4.

  GSE247130. Hassan and Chen ran no accessibility analysis of their transitional
  cells at all: their Sendai figure carries no ATAC panel, and their two Sendai
  ATAC samples deposit no peak files where their other four deposit two each.
  A result here is new for that deposit.

  EITHER WAY the trial can only describe. Every well is one library pooling two
  mice, so the statistical unit is the library and nothing here is a test.
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
from m1b_corrected_pass import MAX_DROP, SEED, arm_value, build_links, wells_corrected

OUT = HERE / "m1d_mature_wells_only"
OUT.mkdir(exist_ok=True)

MIN_CELLS = 50
GROUP = "CLDN4-positive KRT8-positive alveolar"
NEONATAL = "P9"

RULES = {
    "question": ("in the CLDN4-positive KRT8-positive alveolar group, is the AT2 identity programme closed "
                 "at the chromatin level or merely silenced at the RNA level"),
    "relationship_to_M1c": ("M1c refused under U7 because the neonatal P9 wells held more of the labelled group "
                            "than any injured well; its threshold is not moved and its outcome stands"),
    "the_refusal_was_a_result": ("Krt8 and Cldn4 are expressed across immature postnatal alveolar epithelium, so "
                                 "a CLDN4-positive KRT8-positive call cannot separate neonatal developmental "
                                 "immaturity from injury-induced transition. A marker set shared with normal "
                                 "development is not a damage-associated marker set"),
    "R11_stage_rule": "the two P9 neonatal wells leave both the test set and the negative-control set, declared in advance",
    "R12_genotype_rule": "Cebpa is dropped from the AT2 arm in Cebpa-mutant wells and named; the arm must still reach four genes",
    "R13_cell_floor": ("50, not 100: the old floor was set for a per-peak 0.10 threshold this design does not "
                       "contain, and the sham band computed at the same group sizes is the correct small-n guard"),
    "no_accessibility_computed_before_now": True,
    "registration_GSE310539": ("Lynch et al. already report chromVAR AP-1 motif accessibility in this substate and "
                               "note qualitatively that AT2 chromatin features are lost less than AT2 RNA. "
                               "Agreement here is a quantification of a published observation and is registered "
                               "as a re-derivation, not independent confirmation"),
    "registration_GSE247130": ("Hassan and Chen ran no accessibility analysis of their transitional cells; their "
                               "Sendai figure has no ATAC panel and their Sendai ATAC samples deposit no peak "
                               "files. A result here is new for that deposit"),
    "register_status": "Descriptive only. Every well is one library pooling two mice; the unit is the library and nothing here is a test",
    "not_tested": "any contrast between conditions, genotypes, stages or deposits; any motif; any causation",
}


def main() -> None:
    rec = RunRecord(OUT / "m1d_run_record.json", "M1d mature wells only", RULES)
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

    for f in files.values():
        B = budgets[f["deposit"]]["B_rna"]
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
                     "neonatal": w["_neonatal"], "role": ("excluded by R11" if w["_neonatal"]
                                                          else "test" if w["injured"] else "negative control"),
                     "cells": int(sel.sum()), "after_filter": int(keep.sum()),
                     "labelled": int(trans.sum()), "reference": int(ref.sum()),
                     "labelled_pct": round(100 * trans.sum() / max(keep.sum(), 1), 2)})
    labels = pd.DataFrame(rows)
    labels.to_csv(OUT / "m1d_labels.csv", index=False)
    rec.add_output(OUT / "m1d_labels.csv")
    print("\n" + labels.to_string(index=False))

    mature_unin = labels[(~labels.injured) & (~labels.neonatal)]
    inj = labels[labels.injured]
    enrich = float(inj.labelled_pct.max() / max(mature_unin.labelled_pct.max(), 1e-9))
    u7 = enrich >= 3.0
    rec.set("U7_mature_uninjured_max_pct", float(mature_unin.labelled_pct.max()))
    rec.set("U7_injured_max_pct", float(inj.labelled_pct.max()))
    rec.set("U7_enrichment", enrich)
    rec.set("U7_satisfied", u7)
    rec.set("R11_neonatal_label_pct", [float(x) for x in labels[labels.neonatal].labelled_pct])
    print(f"\nU7 on mature wells: uninjured at most {mature_unin.labelled_pct.max():.2f} per cent, "
          f"injured up to {inj.labelled_pct.max():.2f} per cent, enrichment {enrich:.1f}x -> "
          f"{'satisfied' if u7 else 'TRIAL UNREADABLE'}")
    if not u7:
        rec.set("reading", "TRIAL UNREADABLE under U7")
        rec.finish()
        return

    results, arm_rows = [], []
    for w in all_wells:
        if w["_neonatal"]:
            continue
        trans, ref = w["_trans"], w["_ref"]
        m, nref = int(trans.sum()), int(ref.sum())
        role = "test" if w["injured"] else "negative control"
        base = {"deposit": w["deposit"], "well": w["well"], "role": role,
                "n_labelled": m, "n_reference": nref}
        if m < MIN_CELLS or nref < 2 * m:
            results.append({**base, "reading": "not computable",
                            "why": f"fewer than {MIN_CELLS} labelled cells or too small a sham pool"})
            continue

        q = min(4, nref // m - 1)
        f = files[str(w["matrix"])]
        ti = np.flatnonzero(trans)
        ri = rng.permutation(np.flatnonzero(ref))[:(q + 1) * m]
        cells = np.concatenate([ti, ri])
        B = budgets[w["deposit"]]["B_atac"]
        print(f"\n{w['well']} ({role}): {m} labelled, q={q}, {len(cells)} cells to the peak pass")
        det_cells, kept, n_peaks = stream_downsampled_detection(w["matrix"], cells, B, seed=SEED)
        drop_t = 1.0 - float(kept[:m].mean()); drop_r = 1.0 - float(kept[m:].mean())
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

        # R12: Cebpa is deleted in a Cebpa-mutant well and leaves the AT2 arm there
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
                arm_rows.append({**base, "arm": arm, "modality": mod, "genes_used": len(used),
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
    res.to_csv(OUT / "m1d_readings.csv", index=False)
    arms_df.to_csv(OUT / "m1d_arms.csv", index=False)
    rec.add_output(OUT / "m1d_readings.csv"); rec.add_output(OUT / "m1d_arms.csv")
    rec.set("readings", res.to_dict("records"))
    rec.finish()
    print("\n" + res.to_string(index=False))
    if len(arms_df):
        print("\n" + arms_df.drop(columns=["gene_list", "role"]).to_string(index=False))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
