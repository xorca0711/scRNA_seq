#!/usr/bin/env python
"""Trial C14: is trial C12's result a property of the data or of CellChat's database.

Trial C12 scored CellChat's own resource with CellChat's own logic and found
two things. The full-resource ranking of epithelium-to-fibroblast pairs was
dominated by matrix ligands against CD44, so its own guard refused to read it.
And holding the receptor constant, AREG ranked first of the five EGFR ligands
across 26 donors.

Both could be properties of CellChatDB rather than of the tissue. A curated
resource decides which pairs exist to be ranked, and different curations
disagree substantially: the five resources compared here hold between 1,223 and
4,620 pairs. This trial holds the method, the object, the donors and the
compartments fixed, varies only the resource, and asks whether either
conclusion survives.

**A disclosed change to the guard.** C12's abundance guard tested for a
collagen or laminin ligand, or an integrin or syndecan receptor. Applied to
C12's own output it missed FN1 to CD44, which is plainly an abundance pair, and
it did not name CD44 at all even though eleven of that trial's top fifteen
targeted it. The guard used here is wider: fibronectin is added to the ligand
side and CD44 to the receptor side. This is a new guard defined for this trial
before it ran, applied uniformly to all five resources. **C12's verdict is not
restated under it**, because rescoring a finished trial under a later rule is
exactly what this repository does not do.

Frozen rules, set before any resource was scored:

* Object: the cached epithelial and stromal cells of GSE136831 built by trial
  C12, so the cells, donors and normalisation are identical to that trial.
* Method held constant: LIANA's implementation of CellChat, `n_perms=None`, so
  no permutation P value is produced. Its unit would be the cell, and this
  repository does not emit one beside a group question.
* Resources varied: cellchatdb, cellphonedb, consensus, connectomedb2020,
  italk.
* Sources: every epithelial type. Targets: Fibroblast and Myofibroblast.
  Pericytes and smooth muscle stay excluded as mural, as in E6 and C12.
* Unit: the donor. Run per donor, a donor contributing only with at least 50
  cells in both compartments, ranks aggregated by median, and a pair kept only
  if it appears in at least half the contributing donors.

* T1, THE GUARD, per resource. The share of the top fifteen pairs whose ligand
  matches COL, LAM or FN1, or whose receptor matches ITG, SDC or CD44. RULE: a
  resource whose share exceeds one half is reporting transcript abundance and
  its ranking is not read.
* T2, THE PRIMARY QUESTION. Does the guard fire in every resource?
  - fires in all five: the domination is a property of dissociated tissue
    rather than of any one curation, which generalises C12's methods result;
  - fires in some: it is partly a curation property, and which resources stay
    readable is itself the finding.
* T3. Among the EGFR ligands each resource contains, where does AREG rank?
  Consistency across resources tests whether C12's within-receptor result
  depended on CellChatDB.
* T4. Pairwise overlap of the top-fifteen lists across resources, as a Jaccard
  index. Low overlap means the head of the ranking is curation-dependent even
  where the guard does not fire.
* NOT TESTED: whether any pair is a real interaction. Co-expression carries no
  proximity, in any resource.
"""

from __future__ import annotations

import gc
import itertools
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RAW, RunRecord, df_to_markdown  # noqa: E402

OUT = HERE / "c14_does_the_ranking_depend_on_the_database"
OUT.mkdir(exist_ok=True)
CACHE = RAW / "GSE136831" / "c12_epithelial_stromal.h5ad"

RESOURCES = ["cellchatdb", "cellphonedb", "consensus", "connectomedb2020", "italk"]
TARGETS = ["Fibroblast", "Myofibroblast"]
MIN_CELLS = 50
TOP_N = 15
ABUNDANCE_SHARE = 0.5
EGFR_LIGANDS = ["AREG", "HBEGF", "TGFA", "EREG", "BTC", "EGF", "EPGN"]
LIGAND_PREFIXES = ("COL", "LAM")
LIGAND_EXACT = {"FN1"}
RECEPTOR_MARKS = ("ITG", "SDC")
RECEPTOR_EXACT = {"CD44"}

RULES = {
    "question": "is C12's result a property of the data or of CellChat's database",
    "held_constant": "object, donors, compartments, method (LIANA CellChat, n_perms=None)",
    "varied": RESOURCES,
    "resource_sizes_differ": "1,223 to 4,620 pairs, so curation decides what exists to be ranked",
    "disclosed_guard_change": "C12's guard tested COL or LAM ligands and ITG or SDC receptors; applied to "
                              "C12's own output it missed FN1 to CD44 and never named CD44, which eleven of "
                              "its top fifteen targeted. This trial's guard adds FN1 and CD44, was defined "
                              "before running, and is applied uniformly. C12's verdict is NOT restated "
                              "under it.",
    "sources": "every epithelial type", "targets": TARGETS,
    "unit": "the donor; at least " + str(MIN_CELLS) + " cells in both compartments; median rank; a pair "
            "kept only if it appears in at least half the contributing donors",
    "T1": f"guard per resource: share of the top {TOP_N} that are abundance-type; above one half the "
          f"ranking is not read",
    "T2": "does the guard fire in every resource",
    "T3": "where AREG ranks among the EGFR ligands each resource contains",
    "T4": "pairwise Jaccard overlap of the top-" + str(TOP_N) + " lists across resources",
    "not_tested": "whether any pair is a real interaction; co-expression carries no proximity",
    "checkpointing": "each resource writes its own partial file and is skipped when rerun; the first run "
                     "was killed by the operating system during the consensus resource, with no traceback",
}


def is_abundance(ligand, receptor):
    ligand, receptor = str(ligand).upper(), str(receptor).upper()
    if ligand.startswith(LIGAND_PREFIXES) or ligand in LIGAND_EXACT:
        return True
    return any(m in receptor for m in RECEPTOR_MARKS) or receptor in RECEPTOR_EXACT


def rank_one_resource(adata, donors, resource, rec):
    from liana.method import cellchat

    frames = []
    for donor in donors:
        part = adata[adata.obs["donor"] == donor].copy()
        present = part.obs["celltype"].value_counts()
        part = part[part.obs["celltype"].isin(present[present >= 10].index)].copy()
        if not set(TARGETS) & set(part.obs["celltype"]):
            del part
            gc.collect()
            continue
        try:
            cellchat(part, groupby="celltype", resource_name=resource, expr_prop=0.1,
                     use_raw=False, n_perms=None, verbose=False)
        except Exception as error:                      # noqa: BLE001
            print(" ", resource, donor, "skipped:", type(error).__name__, str(error)[:70])
            del part
            gc.collect()
            continue
        res = part.uns["liana_res"]
        res = res[res["target"].isin(TARGETS) & ~res["source"].isin(TARGETS)].copy()
        if not res.empty:
            res["pair"] = res["ligand_complex"] + " to " + res["receptor_complex"]
            best = res.sort_values("lr_probs", ascending=False).drop_duplicates("pair")
            best["rank"] = np.arange(1, len(best) + 1)
            best["donor"] = donor
            frames.append(best[["donor", "pair", "ligand_complex", "receptor_complex", "lr_probs", "rank"]])
        del part, res
        gc.collect()
    if not frames:
        return None
    per_donor = pd.concat(frames, ignore_index=True)
    agg = (per_donor.groupby("pair")
           .agg(n_donors=("donor", "nunique"), median_rank=("rank", "median"),
                ligand=("ligand_complex", "first"), receptor=("receptor_complex", "first"))
           .reset_index())
    agg = agg[agg["n_donors"] >= max(3, len(donors) // 2)].sort_values("median_rank").reset_index(drop=True)
    agg["overall_rank"] = np.arange(1, len(agg) + 1)
    agg["resource"] = resource
    return agg


def main():
    import anndata as ad

    rec = RunRecord(OUT / "c14_run_record.json", "C14 does the ranking depend on the database", RULES)
    if not CACHE.exists():
        raise SystemExit(f"{CACHE} is missing; run trial C12 first to build it")
    rec.add_input(CACHE)
    adata = ad.read_h5ad(CACHE)
    import scanpy as sc
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    counts = adata.obs.groupby(["donor", "category"], observed=True).size().unstack(fill_value=0)
    donors = counts[(counts.get("Epithelial", 0) >= MIN_CELLS)
                    & (counts.get("Stromal", 0) >= MIN_CELLS)].index.tolist()
    rec.set("donors", len(donors))
    print("donors:", len(donors))

    tables, summary, tops = [], [], {}
    for resource in RESOURCES:
        # Checkpoint per resource. The consensus resource is the largest at 4,620
        # pairs and the first run of this trial was killed by the operating
        # system partway through it, with no Python traceback, losing the two
        # resources already scored. Each resource now writes its own file as soon
        # as it finishes and is skipped on a later run, so a kill costs one
        # resource rather than the whole scan.
        partial = OUT / ("c14_partial_" + resource + ".csv")
        if partial.exists():
            print("reusing checkpoint for", resource)
            agg = pd.read_csv(partial)
        else:
            print("scoring", resource)
            agg = rank_one_resource(adata, donors, resource, rec)
            if agg is not None and not agg.empty:
                agg.to_csv(partial, index=False)
        if agg is None or agg.empty:
            summary.append({"resource": resource, "pairs": 0, "note": "no ranking produced"})
            continue
        tables.append(agg)
        top = agg.head(TOP_N).copy()
        top["abundance_type"] = [is_abundance(l, r) for l, r in zip(top["ligand"], top["receptor"])]
        share = float(top["abundance_type"].mean())
        tops[resource] = set(top["pair"])
        egfr = agg[agg["receptor"] == "EGFR"].copy()
        egfr = egfr[egfr["ligand"].isin(EGFR_LIGANDS)].sort_values("median_rank").reset_index(drop=True)
        areg_rank = None
        if len(egfr):
            match = egfr.index[egfr["ligand"] == "AREG"]
            areg_rank = int(match[0]) + 1 if len(match) else None
        summary.append({
            "resource": resource, "pairs": int(len(agg)),
            "abundance_share_top15": round(share, 3),
            "guard_fires": bool(share > ABUNDANCE_SHARE),
            "n_egfr_ligands_present": int(len(egfr)),
            "areg_rank_among_egfr_ligands": areg_rank,
            "egfr_ligand_order": ", ".join(egfr["ligand"].tolist()) if len(egfr) else None,
            "top_pair": top.iloc[0]["pair"],
        })
    if not tables:
        raise SystemExit("no resource produced a ranking; refusing to report")

    full = pd.concat(tables, ignore_index=True)
    full.to_csv(OUT / "c14_rankings_by_resource.csv", index=False)
    rec.add_output(OUT / "c14_rankings_by_resource.csv")
    frame = pd.DataFrame(summary)
    frame.to_csv(OUT / "c14_resource_summary.csv", index=False)
    rec.add_output(OUT / "c14_resource_summary.csv")

    overlaps = []
    for a, b in itertools.combinations(sorted(tops), 2):
        union = tops[a] | tops[b]
        overlaps.append({"resource_a": a, "resource_b": b,
                         "shared_of_top15": len(tops[a] & tops[b]),
                         "jaccard": round(len(tops[a] & tops[b]) / len(union), 3) if union else None})
    overlap = pd.DataFrame(overlaps)
    if len(overlap):
        overlap.to_csv(OUT / "c14_top15_overlap.csv", index=False)
        rec.add_output(OUT / "c14_top15_overlap.csv")

    fires = frame["guard_fires"].dropna()
    all_fire = bool(len(fires) and fires.all())
    rec.set("T2_guard_fires_everywhere", all_fire)
    rec.set("T1_per_resource", frame.to_dict("records"))
    ranks = frame["areg_rank_among_egfr_ligands"].dropna().unique().tolist()
    rec.set("T3_areg_rank_values", ranks)
    rec.set("T4_median_jaccard",
            round(float(overlap["jaccard"].median()), 3) if len(overlap) else None)

    if all_fire:
        reading = ("the abundance domination is a property of dissociated tissue rather than of any one "
                   "curation: the guard fires in all " + str(len(fires)) + " resources, so C12's methods "
                   "result generalises")
    elif fires.any():
        reading = ("the domination is partly a curation property: the guard fires in "
                   + str(int(fires.sum())) + " of " + str(len(fires)) + " resources, and the ones that "
                   "stay readable are the finding")
    else:
        reading = "the guard fires in no resource, which would mean C12's guard was specific to CellChatDB"
    if ranks and len(ranks) == 1 and ranks[0] == 1:
        reading += ". AREG is first among the EGFR ligands in every resource that contains them, so that "
        reading += "result did not depend on CellChatDB."
    elif ranks:
        reading += ". AREG's rank among the EGFR ligands varies across resources: " + str(ranks) + "."

    rec.set("reading", reading)

    lines = ["# Trial C14: does the ranking depend on the database", "",
             "Method, object, donors and compartments held constant; only the resource varies.", "",
             "Reading: " + reading, "", "## Per resource", "",
             df_to_markdown(frame, index=False), ""]
    if len(overlap):
        lines += ["## Overlap of the top " + str(TOP_N) + " lists", "",
                  df_to_markdown(overlap, index=False), ""]
    lines += ["**Not tested:** whether any pair is a real interaction. Co-expression carries no proximity, "
              "in any resource.", ""]
    (OUT / "c14_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "c14_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
