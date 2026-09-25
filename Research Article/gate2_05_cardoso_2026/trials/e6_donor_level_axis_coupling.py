#!/usr/bin/env python
"""Trial E6: does the paper's axis show donor-level coupling in human fibrosis.

The paper's cascade runs from an Areg-producing epithelial state to an
Egfr-bearing alveolar fibroblast. Nothing in this repository has tested that
link, and nothing in the Cardoso deposit can: one library per genotype, the
unit is the library. Trials E1 to E4 established where each molecule sits, which
is a statement about sending cells only.

This trial asks the weakest question that is still a test of the link itself.
If the axis operates in human fibrotic lung, then donors whose epithelium
carries more AREG should carry more receptor-bearing or more activated
fibroblasts. GSE136831 (Adams et al. 2020, doi:10.1126/sciadv.aba1983) has
enough donors with both compartments sampled for a correlation across donors,
with the donor as the unit.

What this can and cannot be. A correlation across donors carries no direction
and no proximity: it cannot say the epithelium drives the fibroblast, and it
cannot say the two cells ever touch. It can be absent, which would be
informative, and it can be shown to be an artefact of depth or of disease
severity, which is why three controls below are pre-registered rather than
optional. It is a coupling test, not a mechanism test, and it is reported as
such whatever it returns.

Frozen rules, set before any correlation was computed:

* Dataset: GSE136831, deposited raw counts, deposited metadata.
* Donors: Subject_Identity, restricted to Disease_Identity in {IPF, Control}.
  COPD donors are out of scope, matching trial E2.
* Epithelium: CellType_Category == "Epithelial", pooled across its types.
* Fibroblasts: Manuscript_Identity in {Fibroblast, Myofibroblast}. Pericytes
  and smooth muscle are excluded, because they are mural rather than
  fibroblast and trials C1b and E4 both showed a mural contribution distorting
  exactly this kind of measure.
* Measures per donor and compartment: the fraction of cells with a non-zero
  count. A compartment contributes a donor-level value only if it holds at
  least 50 cells in that donor, the same floor as E1 to E3.
* Depth per donor and compartment: the median deposited nGene per cell.
* Activation score: the mean of the detection fractions of COL1A1, ACTA2,
  POSTN, CTHRC1 and TNC in the pooled fibroblasts. CTHRC1 is included because
  Tsukui et al. 2020 define the pathological fibroblast by it.
* At least 10 donors are required for any correlation to be computed.

* T8, THE PRIMARY TEST. Spearman correlation across donors between epithelial
  AREG detection and fibroblast EGFR detection. Two-sided, alpha 0.05.
* T9. The same against the fibroblast activation score.
* T10, LIGAND CONTROL. T8 and T9 recomputed with epithelial TGFA in place of
  AREG. TGFA is an EGFR ligand the paper did not follow and which trial E4
  found is not injury-generic. If TGFA correlates as strongly as AREG, the
  result is not specific to the ligand the paper nominated.
* T11, DEPTH CONTROL. Spearman correlation of every variable against the
  median nGene per cell of its own compartment. RULE: if both members of a
  pair correlate with depth at |rho| >= 0.4, that pair is reported as
  confounded and is not read, whatever its own correlation was.
* T12, DISEASE-SEVERITY CONTROL. T8 and T9 recomputed within IPF donors only.
  A correlation present in the pooled set and absent within IPF is consistent
  with both variables tracking disease rather than each other.

* THE READING, fixed in advance so it cannot be chosen afterwards:
  - T8 or T9 significant, not depth-confounded by T11, and still significant
    in the IPF-only stratum of T12: the axis shows donor-level coupling in
    human fibrotic lung. Still correlational, still no proximity.
  - significant in the pooled set but not within IPF: consistent with both
    variables tracking disease severity; no coupling claim.
  - not significant anywhere: no donor-level coupling is detectable at this
    sample size, and the trial reports the correlation it could have detected
    rather than implying the axis is absent.
  - TGFA matching AREG in T10: any coupling found is a general EGFR-ligand
    pattern, not specific to AREG.
* Caveats carried into every reading: the epithelial pool mixes states, so a
  null result could be dilution, which is the same limit trial E1b found;
  detection fractions depend on depth; and 22 IPF donors bound the smallest
  effect this can see.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import RAW, RunRecord, df_to_markdown  # noqa: E402
from mtx_stream import extract_gene_rows, read_lines  # noqa: E402

OUT = HERE / "e6_donor_level_axis_coupling"
OUT.mkdir(exist_ok=True)
SERIES = RAW / "GSE136831"
MATRIX = SERIES / "GSE136831_RawCounts_Sparse.mtx.gz"
GENE_IDS = SERIES / "GSE136831_AllCells.GeneIDs.txt.gz"
BARCODES = SERIES / "GSE136831_AllCells.cellBarcodes.txt.gz"
METADATA = SERIES / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
CACHE = SERIES / "e6_extracted_rows.npz"   # regenerable, outside the tracked tree

LIGANDS = ["AREG", "TGFA", "HBEGF"]
RECEPTOR = "EGFR"
ACTIVATION = ["COL1A1", "ACTA2", "POSTN", "CTHRC1", "TNC"]
GENES = sorted(set(LIGANDS) | {RECEPTOR} | set(ACTIVATION))

FIBROBLAST = {"Fibroblast", "Myofibroblast"}
DISEASES = ["IPF", "Control"]
MIN_CELLS = 50
MIN_DONORS = 10
DEPTH_RHO = 0.4

RULES = {
    "dataset": "GSE136831 (Adams et al. 2020), deposited raw counts and metadata",
    "donors": "Subject_Identity, Disease_Identity in IPF or Control; COPD out of scope as in E2",
    "epithelium": 'CellType_Category == "Epithelial", pooled',
    "fibroblasts": "Manuscript_Identity in Fibroblast or Myofibroblast; pericytes and smooth muscle excluded as mural",
    "measure": "fraction of cells with a non-zero count",
    "depth": "median deposited nGene per cell, per donor and compartment",
    "activation_score": "mean detection of " + ", ".join(ACTIVATION),
    "floors": {"cells_per_donor_compartment": MIN_CELLS, "donors_per_correlation": MIN_DONORS},
    "T8": "primary: Spearman, epithelial AREG detection against fibroblast EGFR detection, two-sided alpha 0.05",
    "T9": "Spearman, epithelial AREG detection against the fibroblast activation score",
    "T10": "ligand control: T8 and T9 recomputed with epithelial TGFA in place of AREG",
    "T11": "depth control: if both members of a pair correlate with depth at |rho| >= "
           + str(DEPTH_RHO) + ", the pair is confounded and is not read",
    "T12": "disease-severity control: T8 and T9 within IPF donors only",
    "reading_fixed_in_advance": {
        "significant_not_confounded_and_holds_within_ipf":
            "the axis shows donor-level coupling in human fibrotic lung; correlational, no proximity",
        "significant_pooled_but_not_within_ipf":
            "consistent with both variables tracking disease severity; no coupling claim",
        "not_significant":
            "no donor-level coupling detectable at this sample size; report the detectable effect size",
        "tgfa_matches_areg":
            "any coupling is a general EGFR-ligand pattern, not specific to AREG",
    },
    "unit": "the donor",
    "caveats": ["a correlation carries no direction and no proximity",
                "the epithelial pool mixes states, so a null could be dilution as in E1b",
                "detection fractions depend on depth",
                "the IPF donor count bounds the smallest detectable effect"],
}


def load_vectors(rec):
    if CACHE.exists():
        blob = np.load(CACHE, allow_pickle=False)
        return {g: blob[g] for g in blob.files}
    for path in (MATRIX, GENE_IDS, BARCODES):
        rec.add_input(path)
    genes = read_lines(GENE_IDS, column=1)[1:]
    vectors = extract_gene_rows(MATRIX, genes, set(GENES))
    np.savez_compressed(CACHE, **vectors)
    return vectors


def detection(vectors, gene, where):
    if gene not in vectors:
        return np.nan
    return float((vectors[gene][where] > 0).mean())


def spearman(frame, a, b, label, note=""):
    from scipy.stats import spearmanr

    joined = frame[[a, b]].dropna()
    out = {"comparison": label, "x": a, "y": b, "n_donors": int(len(joined)), "note": note}
    if len(joined) >= MIN_DONORS:
        res = spearmanr(joined[a], joined[b])
        out["rho"] = round(float(res.statistic), 4)
        out["p_value"] = round(float(res.pvalue), 6)
        out["significant"] = bool(out["p_value"] < 0.05)
    else:
        out["rho"] = None
        out["p_value"] = None
        out["significant"] = None
        out["note"] = (note + "; " if note else "") + "fewer donors than the frozen floor; no test"
    return out


def detectable_rho(n):
    """The smallest |rho| a two-sided Spearman at alpha 0.05 would call, approximately."""
    if n < 4:
        return None
    return round(float(1.96 / np.sqrt(n - 1)), 3)


def main():
    rec = RunRecord(OUT / "e6_run_record.json",
                    "E6 donor-level coupling of the AREG to EGFR axis in human fibrosis", RULES)
    rec.add_input(METADATA)
    meta = pd.read_csv(METADATA, sep="\t")
    barcodes = read_lines(BARCODES)
    if len(barcodes) != len(meta):
        meta = meta.set_index("CellBarcode_Identity").reindex(barcodes).reset_index()
    vectors = load_vectors(rec)
    n_cells = len(next(iter(vectors.values())))
    if n_cells != len(meta):
        raise SystemExit(f"metadata has {len(meta)} rows but the matrix has {n_cells} cells")

    in_scope = meta["Disease_Identity"].isin(DISEASES).to_numpy()
    is_epi = in_scope & (meta["CellType_Category"] == "Epithelial").to_numpy()
    is_fib = in_scope & meta["Manuscript_Identity"].isin(FIBROBLAST).to_numpy()
    rec.set("cells_in_scope", int(in_scope.sum()))
    rec.set("epithelial_cells", int(is_epi.sum()))
    rec.set("fibroblast_cells", int(is_fib.sum()))

    rows = []
    for donor, group in meta[in_scope].groupby("Subject_Identity"):
        epi = np.where(is_epi & (meta["Subject_Identity"] == donor).to_numpy())[0]
        fib = np.where(is_fib & (meta["Subject_Identity"] == donor).to_numpy())[0]
        row = {"donor": donor, "disease": group["Disease_Identity"].iloc[0],
               "n_epithelial": len(epi), "n_fibroblast": len(fib)}
        if len(epi) >= MIN_CELLS:
            for gene in LIGANDS:
                row["epi_" + gene] = round(detection(vectors, gene, epi), 4)
            row["epi_depth"] = float(np.median(meta["nGene"].to_numpy()[epi]))
        if len(fib) >= MIN_CELLS:
            row["fib_EGFR"] = round(detection(vectors, RECEPTOR, fib), 4)
            scores = []
            for gene in ACTIVATION:
                value = detection(vectors, gene, fib)
                row["fib_" + gene] = round(value, 4)
                scores.append(value)
            row["fib_activation"] = round(float(np.nanmean(scores)), 4)
            row["fib_depth"] = float(np.median(meta["nGene"].to_numpy()[fib]))
        rows.append(row)

    table = pd.DataFrame(rows)
    table.to_csv(OUT / "e6_per_donor.csv", index=False)
    rec.add_output(OUT / "e6_per_donor.csv")
    usable = table.dropna(subset=["epi_AREG", "fib_EGFR", "fib_activation"])
    rec.set("donors_with_both_compartments", int(len(usable)))
    rec.set("donors_by_disease", usable["disease"].value_counts().to_dict())
    rec.set("smallest_detectable_rho", {"pooled": detectable_rho(len(usable)),
                                        "ipf_only": detectable_rho(int((usable["disease"] == "IPF").sum()))})

    ipf = usable[usable["disease"] == "IPF"]
    results = [
        spearman(usable, "epi_AREG", "fib_EGFR", "T8 primary"),
        spearman(usable, "epi_AREG", "fib_activation", "T9"),
        spearman(usable, "epi_TGFA", "fib_EGFR", "T10 ligand control"),
        spearman(usable, "epi_TGFA", "fib_activation", "T10 ligand control"),
        spearman(usable, "epi_AREG", "epi_depth", "T11 depth", "epithelial depth"),
        spearman(usable, "epi_TGFA", "epi_depth", "T11 depth", "epithelial depth"),
        spearman(usable, "fib_EGFR", "fib_depth", "T11 depth", "fibroblast depth"),
        spearman(usable, "fib_activation", "fib_depth", "T11 depth", "fibroblast depth"),
        spearman(ipf, "epi_AREG", "fib_EGFR", "T12 IPF only"),
        spearman(ipf, "epi_AREG", "fib_activation", "T12 IPF only"),
    ]
    frame = pd.DataFrame(results)
    frame.to_csv(OUT / "e6_correlations.csv", index=False)
    rec.add_output(OUT / "e6_correlations.csv")

    def get(label, x, y):
        for r in results:
            if r["comparison"].startswith(label) and r["x"] == x and r["y"] == y:
                return r
        return {}

    depth = {(r["x"]): abs(r["rho"]) if r["rho"] is not None else 0.0
             for r in results if r["comparison"].startswith("T11")}
    verdicts = []
    for label, x, y, stratum in (("T8 primary", "epi_AREG", "fib_EGFR", "pooled"),
                                 ("T9", "epi_AREG", "fib_activation", "pooled"),
                                 ("T10 ligand control", "epi_TGFA", "fib_EGFR", "pooled"),
                                 ("T10 ligand control", "epi_TGFA", "fib_activation", "pooled")):
        primary = get(label, x, y)
        confounded = depth.get(x, 0.0) >= DEPTH_RHO and depth.get(y, 0.0) >= DEPTH_RHO
        within = get("T12", x, y) if x == "epi_AREG" else {}
        control = get("T10", "epi_TGFA", y) if x == "epi_AREG" else {}
        if confounded:
            reading = RULES["reading_fixed_in_advance"]["significant_pooled_but_not_within_ipf"]
            reading = "both variables correlate with depth at or above the frozen line; not read"
        elif x != "epi_AREG":
            reading = ("control pair: " + ("significant" if primary.get("significant") else "not significant")
                       + ", and " + ("depth-confounded" if confounded else "not depth-confounded"))
        elif primary.get("significant") and within.get("significant"):
            reading = RULES["reading_fixed_in_advance"]["significant_not_confounded_and_holds_within_ipf"]
        elif primary.get("significant"):
            reading = RULES["reading_fixed_in_advance"]["significant_pooled_but_not_within_ipf"]
        else:
            reading = RULES["reading_fixed_in_advance"]["not_significant"]
        if (primary.get("significant") and control.get("significant")
                and control.get("rho") is not None and primary.get("rho") is not None
                and abs(control["rho"]) >= abs(primary["rho"])):
            reading += " " + RULES["reading_fixed_in_advance"]["tgfa_matches_areg"]
        verdicts.append({"test": label, "stratum": stratum, "rho": primary.get("rho"),
                         "p_value": primary.get("p_value"),
                         "depth_confounded": bool(confounded), "reading": reading})
    rec.set("verdicts", verdicts)

    lines = ["# Trial E6: donor-level coupling of the AREG to EGFR axis in human fibrosis", "",
             "Donors with both compartments above the 50-cell floor: " + str(len(usable))
             + " (" + ", ".join(f"{k} {v}" for k, v in usable["disease"].value_counts().items()) + ").",
             "Smallest |rho| a two-sided test could call here: pooled "
             + str(detectable_rho(len(usable))) + ", IPF only "
             + str(detectable_rho(len(ipf))) + ".", "",
             "## Pre-registered readings", "",
             df_to_markdown(pd.DataFrame(verdicts), index=False), "",
             "## All correlations, including the three control families", "",
             df_to_markdown(frame, index=False), ""]
    (OUT / "e6_summary.md").write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(OUT / "e6_summary.md")
    rec.finish()
    print("\n".join(lines))


if __name__ == "__main__":
    main()
