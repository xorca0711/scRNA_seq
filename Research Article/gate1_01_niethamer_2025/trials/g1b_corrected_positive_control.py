#!/usr/bin/env python
"""Trial G1b: the corrected positive-control gate for trial G1.

G1's gate was HALLMARK_G2M_CHECKPOINT up in pooled active repair (6 to 25 dpi)
against injury resolution (42 and 90 dpi). It failed in both compartments, and
trial N1's own cycling fractions say why: in this deposit neither the myeloid nor
the endothelial lineage has an early proliferation peak in the per-animal
cell-cycle fraction (myeloid 34 per cent at 6 dpi against 43 at 90; endothelium
rising through 366), so the expectation, not the ranking, was wrong. This trial
replaces the gate with expectations the register already holds and does not
touch the TEST contrast, whose tables it reads from G1 unchanged.

Rules, frozen before this trial reads any count, with one honesty statement
------------------------------------------------------------------------
R0  Contamination stated: the author of this rule had seen G1's positive-control
    and TEST tables. The new gates below are justified by register rows that
    predate any GSEA (C3, C12, C18) and by marker tables trials 06 and 11 wrote
    from cluster contrasts, not by any G1 enrichment result; but the reader
    should weigh the gates knowing they were written second. Nothing this trial
    licenses can exceed the statuses R8 of G1 allows.
R1  The gates are repository marker sets, each the top-20 markers of the
    subcluster that carries a registered injury state, from tables written by
    the clustering trials (a within-cell contrast, independent of the phase
    contrast tested here, on the same deposit; partially circular and said so):
      myeloid:   the subcluster of trial 11 holding most of the inflammatory
                 monocyte (iMON) cells, whose 6 dpi peak is rows C12 and C18;
                 gate contrast 6 and 11 dpi (4 animals) against 366 dpi (3).
      capillary: subcluster 5 of script 06, the injury-induced capillary state
                 (iCAP, row C3: 37.5 per cent of capillary cells at 25 dpi, 21.7
                 at 366); gate contrast 19 and 25 dpi (4 animals) against 366 (3).
    The gate clears when the set is up in the early arm at gseapy FDR < 0.05
    among the 50 hallmark sets plus itself, and its expression-matched p is
    below 0.05.
R2  Everything else is G1's: unit, floors, compartments, pseudobulk, ranking,
    engine, seed. GO biological process is not run here; the gate needs only the
    hallmark collection and the repository set.
R3  Reading: where the gate clears, G1's TEST tables for that compartment become
    readable and are read under G1's R8; where it fails, the compartment stays
    unreadable and the honest conclusion is that this design cannot detect a
    registered composition shift of that size with four against three animals.

Outputs (g1b_corrected_positive_control/): g1b_units.csv, g1b_pc_sets.csv,
g1b_hallmark_PC_<compartment>.csv, g1b_positive_control.csv, g1b_summary.md,
g1b_run_record.json.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import g1_gsea_by_phase as g1  # noqa: E402
import gsea_utils as gu  # noqa: E402
from trial_utils import RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "g1b_corrected_positive_control"
GATES = {
    "myeloid": {"early_days": [6, 11], "late_days": [366], "source": "trial 11 iMON subcluster top-20 markers", "rows": "C12, C18"},
    "capillary": {"early_days": [19, 25], "late_days": [366], "source": "script 06 capillary subcluster 5 (iCAP) top-20 markers", "rows": "C3"},
}
RULES = {
    "R0_contamination": "rule written after seeing G1's tables; gates justified by rows C3, C12, C18 and by clustering-trial marker tables, not by G1 enrichment results",
    "R1_gates": {k: {**v} for k, v in GATES.items()},
    "R1_clearing": "set up in the early arm at gseapy FDR < 0.05 among hallmark plus itself, and matched p < 0.05",
    "R2_inherited": "G1 rules R1, R2, R5, R7 unchanged; GO BP not run",
    "R3_reading": "gate clears: G1 TEST tables readable under G1 R8; gate fails: compartment stays unreadable",
}


def repository_sets(rec: RunRecord) -> dict[str, set[str]]:
    my_meta = g1.MOUSE / "myeloid_focus" / "tables" / "myeloid_cell_metadata.csv"
    my_mk = g1.MOUSE / "myeloid_focus" / "tables" / "myeloid_subcluster_markers_top20.csv"
    cap_mk = g1.MOUSE / "regeneration_focus" / "tables" / "capillary_subcluster_markers_top20.csv"
    for p in (my_meta, my_mk, cap_mk):
        rec.add_input(p)
    meta = pd.read_csv(my_meta)
    sub_col = "sub_r0.5"
    imon_sub = meta[meta["label"] == "iMON"][sub_col].value_counts().idxmax()
    my = pd.read_csv(my_mk)
    my_group_col = "group" if "group" in my.columns else my.columns[0]
    my_set = my[my[my_group_col].astype(str) == str(imon_sub)]["names"].astype(str).tolist()
    cap = pd.read_csv(cap_mk)
    cap_set = cap[cap["group"].astype(str) == "5"]["names"].astype(str).tolist()
    assert len(my_set) >= 15 and len(cap_set) >= 15, (len(my_set), len(cap_set))
    rec.set("imon_subcluster", str(imon_sub))
    rec.set("imon_cells_in_that_subcluster", int((meta[sub_col].astype(str) == str(imon_sub)) [meta["label"] == "iMON"].sum()))
    sets = {"REPO_MYELOID_IMON_SUBCLUSTER_TOP20": set(my_set), "REPO_CAPILLARY_ICAP_SUBCLUSTER5_TOP20": set(cap_set)}
    pd.DataFrame([{"set": k, "gene": g} for k, v in sets.items() for g in sorted(v)]).to_csv(OUT / "g1b_pc_sets.csv", index=False)
    rec.add_output(OUT / "g1b_pc_sets.csv")
    return sets


def main() -> int:
    OUT.mkdir(exist_ok=True)
    rec = RunRecord(OUT / "g1b_run_record.json", "G1b corrected positive-control gate for G1", RULES)
    rec.set("gene_set_files", {"hallmark": gu.gene_set_facts(g1.GMT["hallmark"])})
    repo = repository_sets(rec)
    hallmark = gu.read_gmt(g1.GMT["hallmark"])
    units, groups, pb, genes, day_of = g1.pseudobulk_compartments(rec, OUT)
    rec.add_output(OUT / "g1b_units.csv")
    rng = np.random.default_rng(g1.SEED)
    gate_set = {"myeloid": "REPO_MYELOID_IMON_SUBCLUSTER_TOP20", "capillary": "REPO_CAPILLARY_ICAP_SUBCLUSTER5_TOP20"}
    pc_rows, gate_ok = [], {}
    for comp, gate in GATES.items():
        name = gate_set[comp]
        sets = {"hallmark": {**hallmark, name: repo[name]}}
        tag = f"PC_{comp}"
        got = g1.gsea_contrast(tag, pb, genes, g1.arm_indices(groups, day_of, comp, gate["early_days"]),
                               g1.arm_indices(groups, day_of, comp, gate["late_days"]), sets, rec, OUT, "g1b")
        if got is None:
            gate_ok[comp] = False
            continue
        ranking, res = got
        hall = res["hallmark"]
        row = hall[hall["set"] == name]
        m = g1.matched_for(ranking, sets["hallmark"], row, rng)
        ok = bool(m and m[0]["clears_both"] and m[0]["nes"] > 0)
        gate_ok[comp] = ok
        pc_rows.append({"compartment": comp, "set": name, "early_days": str(gate["early_days"]), "late_days": str(gate["late_days"]),
                        "animals": f"{len(g1.arm_indices(groups, day_of, comp, gate['early_days']))} against {len(g1.arm_indices(groups, day_of, comp, gate['late_days']))}",
                        **({k: m[0][k] for k in ("size", "nes", "fdr", "matched_p", "matched_nes")} if m else {}), "clears": ok})
        for ctx in ("HALLMARK_G2M_CHECKPOINT", "HALLMARK_INFLAMMATORY_RESPONSE"):
            r = hall[hall["set"] == ctx]
            if len(r):
                pc_rows.append({"compartment": comp, "set": ctx + " (context only)", "early_days": str(gate["early_days"]),
                                "late_days": str(gate["late_days"]), "animals": "", "size": r["size"].iloc[0], "nes": float(r["nes"].iloc[0]),
                                "fdr": float(r["fdr"].iloc[0]), "matched_p": np.nan, "matched_nes": np.nan, "clears": ""})
        g1.log(f"{tag}: gate {'clears' if ok else 'FAILS'}")
    pc = pd.DataFrame(pc_rows)
    pc.to_csv(OUT / "g1b_positive_control.csv", index=False)
    rec.add_output(OUT / "g1b_positive_control.csv")
    rec.set("gate", gate_ok)

    # read G1's TEST tables where the gate clears
    g1_matched = pd.read_csv(g1.OUT / "g1_matched_null.csv") if (g1.OUT / "g1_matched_null.csv").exists() else pd.DataFrame()
    g1_go = pd.read_csv(g1.OUT / "g1_gobp_cleared.csv") if (g1.OUT / "g1_gobp_cleared.csv").exists() else pd.DataFrame()
    for p in ("g1_matched_null.csv", "g1_gobp_cleared.csv"):
        if (g1.OUT / p).exists():
            rec.add_input(g1.OUT / p)
    readable = {}
    for comp, ok in gate_ok.items():
        if not ok:
            readable[comp] = "unreadable"
            continue
        hm = g1_matched[(g1_matched.contrast == "TEST") & (g1_matched.compartment == comp)] if len(g1_matched) else g1_matched
        go = g1_go[g1_go.compartment == comp] if len(g1_go) else g1_go
        readable[comp] = {"hallmark_clearing_both": hm[hm.clears_both]["set"].tolist() if len(hm) else [],
                          "gobp_clearing_both": go[go.clears_both]["set"].tolist() if len(go) else [],
                          "gobp_fdr05": int(len(go))}
    rec.set("test_reading", readable)

    md = ["# Trial G1b: the corrected positive-control gate for G1", "",
          "Generated by `g1b_corrected_positive_control.py`; rules R0 to R3 in the docstring and the run",
          "record. R0 states that this rule was written after G1's tables were seen; the gates rest on",
          "register rows C3, C12 and C18 and on the clustering trials' marker tables.", "",
          "## Gates", "", df_to_markdown(pc.round(4), index=False) if len(pc) else "none computed", ""]
    for comp, ok in gate_ok.items():
        md.append(f"- **{comp}**: gate {'clears' if ok else 'fails'}; G1's TEST contrast in this compartment is "
                  + ("readable under G1 R8." if ok else "still unreadable: with four against three animals this design does not detect a composition shift the register holds (R3)."))
    md += ["", "## G1 TEST tables read under this gate", ""]
    for comp, r in readable.items():
        if r == "unreadable":
            md.append(f"- {comp}: unreadable.")
        else:
            md.append(f"- {comp}: hallmark sets clearing both nulls: {', '.join(r['hallmark_clearing_both']) or 'none'}; "
                      f"GO BP sets at FDR < 0.05: {r['gobp_fdr05']}, clearing both nulls: {', '.join(r['gobp_clearing_both']) or 'none'} "
                      "(direction: positive NES is higher in injury resolution, 42 and 90 dpi, than at 366 dpi).")
    md += ["", "Proposed statuses follow G1 R8 (Exploratory in one compartment, Descriptive only in both; nothing",
           "Validated with three animals at 366 dpi) and await the owner's retain or reject."]
    (OUT / "g1b_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    rec.add_output(OUT / "g1b_summary.md")
    rec.finish()
    g1.log("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
