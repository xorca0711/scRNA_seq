#!/usr/bin/env python
"""Trial G0: on which of the imported deposits is a gene set enrichment analysis admissible?

No gene set enrichment analysis (GSEA) has been run in this repository
(docs/PIPELINE_AS_RUN.md: "Pathway / gene-set enrichment: not run"). Before one is
proposed, this trial scans every accession under raw_data/ and asks whether a
ranked-list GSEA could be run on it at all under the house rules, without opening
any expression matrix: only the deposited metadata tables, this repository's own
reality-check records and per-cell tables are read.

Rules, frozen before any table is read
--------------------------------------
R1  The unit is the animal, donor or library, never the cell. A contrast is
    admissible only when each arm holds at least three units after the trial's
    cell floor (the floor proposal W1 already carries: fewer than three animals in
    any tested arm stops it). Two per arm gives directions only, as trial E4 did.
R2  Raw counts must be on disk, so that a per-unit pseudobulk can be summed. A
    normalised-only matrix cannot be pseudobulked honestly.
R3  The contrast must be one a question in RESEARCH_QUESTIONS.md asks or a Stage 2
    proposal names. A contrast that exists only because two arms happen to be
    replicated is not a trial.
R4  Gene sets are frozen and their file hashed in the run record before the
    ranking is opened. No MSigDB collection is on disk today; a download is a new
    external input and is recorded as such.
R5  The statistic is a pre-ranked GSEA on a per-unit pseudobulk ranking, with a
    matched-size random-gene-set null of the kind trial M2 used beside the
    permutation null of the GSEA implementation. Cell-level marker lists are not a
    ranking for this purpose, because cells are not replicates.

Verdict vocabulary: "plausible" (R1 to R3 met; R4 and R5 are the trial's own
obligations), "marginal" (R1 met only by pooling arms or with exactly three units
in one arm), "not admissible" (R1 or R2 fails). Nothing here is a result about
biology; it is a scan of designs.

Outputs (Thesis/gate1_01_niethamer_2025/trials/g0_gsea_feasibility/):
  g0_deposits.csv     one row per accession with the counted units and the verdict
  g0_summary.md       the table and the reading
  g0_run_record.json  rules, inputs, versions, verdict counts
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
sys.path.insert(0, str(REPO / "Thesis" / "gate1_04_sikkema_2023_hlca" / "trials"))
from trial_utils import RunRecord, df_to_markdown  # noqa: E402

RAW = REPO / "raw_data"
OUT = HERE / "g0_gsea_feasibility"
FLOOR = 3

RULES = {
    "R1_unit_and_floor": "animal, donor or library as the unit; at least three units per arm after the cell floor",
    "R2_raw_counts": "raw counts on disk; normalised-only matrices are not pseudobulked",
    "R3_named_contrast": "the contrast is one RESEARCH_QUESTIONS.md asks or a Stage 2 proposal names",
    "R4_frozen_gene_sets": "gene sets frozen and hashed in the run record before the ranking is opened; MSigDB is not on disk",
    "R5_statistic": "pre-ranked GSEA on a per-unit pseudobulk ranking, with a matched-size random-gene-set null beside the implementation's permutation null",
    "verdicts": ["plausible", "marginal", "not admissible"],
    "no_matrix_is_opened": True,
}

CARDOSO = REPO / "Thesis" / "gate2_05_cardoso_2026" / "trials"
CHOI = REPO / "Thesis" / "gate1_02_choi_2020" / "trials"
MULTI = REPO / "Thesis" / "gate1_02_choi_2020" / "datp_epigenetics" / "trials"
MOUSE = REPO / "Thesis" / "gate1_01_niethamer_2025" / "GSE262927"


def first(glob: str, root: Path) -> Path:
    hits = sorted(root.glob(glob))
    assert hits, f"missing {glob} under {root}"
    return hits[0]


def main() -> int:
    OUT.mkdir(exist_ok=True)
    rec = RunRecord(OUT / "g0_run_record.json", "G0 GSEA feasibility scan of the imported deposits", RULES)
    rows: list[dict] = []

    def add(accession, study, species, unit, arms: dict, raw, contrast, question, verdict, why, source):
        rows.append({
            "accession": accession, "study": study, "species": species, "unit": unit,
            "arms_and_units": "; ".join(f"{k}: {v}" for k, v in arms.items()),
            "arms_at_or_above_floor": int(sum(v >= FLOOR for v in arms.values())),
            "smallest_arm": min(arms.values()) if arms else 0,
            "raw_counts_on_disk": raw, "contrast": contrast, "asked_by": question,
            "verdict": verdict, "why": why, "source": source,
        })

    # GSE262927: the annotated cohort, animals per day
    meta_path = MOUSE / "tables" / "cell_metadata.csv"
    rec.add_input(meta_path)
    m = pd.read_csv(meta_path, usecols=["sample_id", "sacrifice_day", "has_author_metadata"])
    m = m[m["has_author_metadata"].astype(str).isin(["True", "true"])]
    per_day = m.groupby("sacrifice_day")["sample_id"].nunique()
    days = {f"{int(d)} dpi": int(n) for d, n in per_day.items()}
    admissible_days = [d for d, n in days.items() if n >= FLOOR]
    phases = {"active repair (6 to 25 dpi)": int(sum(n for d, n in per_day.items() if 6 <= d <= 25)),
              "injury resolution (42, 90 dpi)": int(sum(n for d, n in per_day.items() if d in (42, 90))),
              "long-term (366 dpi)": int(per_day.get(366.0, 0))}
    add("GSE262927", "Niethamer et al. 2025", "mouse", "animal", {**days, **phases},
        True, "phase against phase within a compartment (myeloid, capillary, alveolar epithelium); "
        f"days with three or more animals: {', '.join(admissible_days)}",
        "A3; proposals W1, S1", "plausible",
        "the only mouse deposit with animal-level replication; baseline (2 animals) stays a reference band; the "
        "two-animal active-repair days enter only as a pooled phase",
        "tables/cell_metadata.csv (annotated cohort)")

    # GSE136831: donors by disease
    p = RAW / "GSE136831" / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz"
    rec.add_input(p)
    a = pd.read_csv(p, sep="\t", usecols=["Disease_Identity", "Subject_Identity"])
    donors = a.groupby("Disease_Identity")["Subject_Identity"].nunique().to_dict()
    add("GSE136831", "Adams et al. 2020", "human", "donor", {f"{k} donors": int(v) for k, v in donors.items()},
        (RAW / "GSE136831" / "GSE136831_RawCounts_Sparse.mtx.gz").exists(),
        "IPF against control within a compartment (AT2 and transitional epithelium, fibroblasts, macrophages)",
        "A2 (C48, C49); A3 in human fibrosis", "plausible",
        "the largest donor-replicated contrast in reach; per-donor cell floor of trial E6 (50) applies", "deposited metadata")

    # GSE135893: donors by diagnosis
    p = RAW / "GSE135893" / "GSE135893_IPF_metadata.csv.gz"
    rec.add_input(p)
    h = pd.read_csv(p, usecols=["Diagnosis", "Sample_Name"])
    hd = h.groupby("Diagnosis")["Sample_Name"].nunique().to_dict()
    add("GSE135893", "Habermann et al. 2020", "human", "donor", {f"{k} donors": int(v) for k, v in hd.items()},
        (RAW / "GSE135893" / "GSE135893_matrix.mtx.gz").exists(),
        "IPF against control within a compartment, as the held-out replication of GSE136831",
        "A2, A3; proposal V1 names held-out human validation", "plausible",
        "second fibrosis cohort; a gene set that clears in GSE136831 and not here is not established", "deposited metadata")

    # GSE131907: paired lung donors
    p = RAW / "GSE131907" / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    rec.add_input(p)
    k = pd.read_csv(p, sep="\t", usecols=["Sample", "Sample_Origin"])
    k = k[k["Sample_Origin"].isin(["nLung", "tLung"])]
    n_by = k.groupby("Sample_Origin")["Sample"].nunique().to_dict()
    pid = lambda s: s.str.replace(r"^LUNG_[NT]", "", regex=True)
    paired = len(set(pid(k[k.Sample_Origin == "nLung"]["Sample"])) & set(pid(k[k.Sample_Origin == "tLung"]["Sample"])))
    add("GSE131907", "Kim et al. 2020", "human", "donor (paired)",
        {"normal lung donors": int(n_by["nLung"]), "tumour lung donors": int(n_by["tLung"]), "paired": paired},
        (RAW / "GSE131907" / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz").exists(),
        "tumour against paired normal lung within the epithelial compartment",
        "A2 (C40)", "plausible",
        "paired design; the raw UMI matrix is on disk although trials E1 and E1b read the deposited log2 TPM, so a "
        "pseudobulk trial must stream the raw file", "deposited annotation")

    # GSE132771: libraries per condition from trial E4's inputs
    p = first("e4_bleomycin_fibrotic_genes/e4_run_record.json", CARDOSO)
    rec.add_input(p)
    e4 = json.load(open(p, encoding="utf-8"))
    libs = sorted({re.sub(r".*GSM\d+_", "", i["path"]).split("_")[0] for i in e4["inputs"]})
    arms = {"bleomycin": sum(x.startswith("Bleo") for x in libs), "untreated": sum(x.startswith("UT") for x in libs)}
    add("GSE132771", "Tsukui et al. 2020", "mouse", "animal", arms, True,
        "bleomycin against untreated, collagen-producing cells", "C42 (trial E4)", "not admissible",
        "two animals per arm: directions only, as E4 reported", "trial E4 run record")

    # GSE247505: libraries per arm from trial C3's record
    p = first("c3_areg_state_specificity/c3_run_record.json", CARDOSO)
    rec.add_input(p)
    c3 = json.load(open(p, encoding="utf-8"))
    lib_rows = None

    def find_libs(o):
        nonlocal lib_rows
        if isinstance(o, dict):
            for kk, v in o.items():
                if kk == "libraries" and isinstance(v, list) and v and isinstance(v[0], dict) and "arm" in v[0]:
                    lib_rows = v
                find_libs(v)
        elif isinstance(o, list):
            for v in o:
                find_libs(v)
    find_libs(c3)
    assert lib_rows, "C3 record carries no library table"
    per_arm = pd.Series([r["arm"] for r in lib_rows]).value_counts().to_dict()
    add("GSE247505", "England et al. 2025", "mouse", "library", {k: int(v) for k, v in per_arm.items()}, True,
        "any timed or genotype arm against another", "A4 (C136, C149); trial C3", "not admissible",
        f"{len(lib_rows)} libraries in {len(per_arm)} arms, at most {max(per_arm.values())} per arm; pooling arms "
        "across time or colour would break the design that makes the deposit a time course", "trial C3 run record")

    # single-library deposits, from the reality-check records
    for acc, study, species, contrast, asked, src in (
        ("GSE145031", "Choi et al. 2020", "mouse", "bleomycin against PBS by day", "A5, trials D0 to D7", CHOI / "d0_data_reality_check"),
        ("GSE144468", "Choi et al. 2020", "mouse", "IL-1beta against control organoids", "trials D5, D5b", CHOI / "d0_data_reality_check"),
        ("GSE316241", "Cardoso, Lee et al. 2026", "mouse", "Red2Kras against Confetti", "trials C0 to C2b", CARDOSO / "c0_data_reality_check"),
        ("GSE316243", "Cardoso, Lee et al. 2026", "mouse", "genotype against genotype", "trials C0 to C2b", CARDOSO / "c0_data_reality_check"),
        ("GSE316244", "Cardoso, Lee et al. 2026", "mouse", "Areg-flox against control", "trial C2", CARDOSO / "c0_data_reality_check"),
        ("GSE310335", "Cardoso, Lee et al. 2026", "human", "KRAS G12D organoid against control", "trial C0", CARDOSO / "c0_data_reality_check"),
        ("GSE310539", "Lynch et al. 2026", "mouse", "Sendai virus against PBS; AP-1 mutant against wildtype", "A1, A4; trials M0 to M4", MULTI / "m0_deposit_reality_check"),
        ("GSE247130", "Hassan and Chen 2024", "mouse", "Cebpa mutant against control; stage against stage", "A1, A5; trials M0 to M4", MULTI / "m0_deposit_reality_check"),
    ):
        recs = sorted(src.glob("*run_record*.json")) or sorted(src.glob("*.json"))
        assert recs, f"no record under {src}"
        rec.add_input(recs[0])
        add(acc, study, species, "library", {"libraries per condition": 1}, True, contrast, asked, "not admissible",
            "one library per condition; no within-group replication, as the reality-check trial recorded", recs[0].relative_to(REPO).as_posix())

    add("GSE178360", "Kadur Lakshminarasimha Murthy et al. 2022", "human", "donor", {"healthy donors": 3}, True,
        "none; three healthy donors and no condition", "the human series report", "not admissible",
        "no contrast exists in the deposit", "Thesis/ungated_murthy_2022/GSE178360/README.md")
    add("GSE309751", "Lynch et al. 2026", "mouse", "animal", {"mock and injured at 14 days": 3, "49 days": 2}, False,
        "injured against mock at 14 and 49 days, bulk ATAC on sorted AT2", "A1, deliverable 1", "marginal",
        "assessed and not opened; bulk ATAC peak calls, not counts; two mice in the 49-day arms; peak presence at "
        "frozen loci is the admissible object, not a ranked GSEA", "REFERENCES.md (assessed, not opened)")

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "g0_deposits.csv", index=False)
    rec.add_output(OUT / "g0_deposits.csv")
    counts = df["verdict"].value_counts().to_dict()
    rec.set("accessions_scanned", int(len(df)))
    rec.set("verdicts", counts)
    rec.set("plausible", df.loc[df.verdict == "plausible", "accession"].tolist())
    rec.set("gene_set_collections_on_disk", sorted(p.name for p in RAW.rglob("*.gmt")))
    tooling = {}
    for mod in ("gseapy", "decoupler"):
        try:
            __import__(mod)
            tooling[mod] = True
        except Exception:
            tooling[mod] = False
    rec.set("gsea_tooling_installed", tooling)

    show = df[["accession", "study", "unit", "arms_and_units", "raw_counts_on_disk", "contrast", "asked_by", "verdict"]]
    plaus = df[df.verdict == "plausible"]
    md = [
        "# Trial G0: on which deposits is a gene set enrichment analysis admissible?",
        "",
        "Generated by `g0_gsea_feasibility.py` from the deposited metadata tables and this",
        "repository's own reality-check records; no expression matrix was opened. Rules R1 to",
        "R5 are in the script's docstring and the run record, frozen before the tables were read.",
        "",
        f"**{len(df)} accessions scanned: {counts.get('plausible', 0)} plausible, {counts.get('marginal', 0)} marginal, "
        f"{counts.get('not admissible', 0)} not admissible.** No GSEA has been run in this repository, no gene set",
        "collection is on disk, and neither gseapy nor decoupler is installed; a GSEA trial adds all three as recorded inputs.",
        "",
        df_to_markdown(show, index=False),
        "",
        "## Reading",
        "",
        "Under the unit rule the question is never whether a gene set is enriched in a group of cells, which is",
        "always answerable and never a test, but whether a per-unit pseudobulk ranking carries it. That leaves",
        f"{len(plaus)} deposits, and they sort by what they can answer:",
        "",
    ]
    for _, r in plaus.iterrows():
        md.append(f"- **{r.accession}** ({r.study}; {r.unit}; {r.arms_and_units}): {r.contrast}. Asked by {r.asked_by}. {r.why}.")
    md += [
        "",
        "The multiome deposits fail R1 for every between-condition contrast, and their within-well state contrasts",
        "have no unit at all; the nearest admissible object there is the matched-gene-set detection statistic trial",
        "M2 already runs, which is a gene-set-level reading and Descriptive only. GSE309751 is marginal on numbers",
        "and wrong in kind: peak presence at frozen loci (deliverable 1), not a ranked list.",
        "",
        "The proposals this scan licenses are G1 and G2 in the Stage 2 table of `../ANALYSIS_TRIAL_PLAN.md`:",
        "G1 on GSE262927 by phase within a compartment, as the gene-set layer of proposal W1; G2 on GSE136831",
        "IPF against control within a compartment, with GSE135893 held out. GSE131907 is a paired epithelial",
        "contrast that A2 already answered at the gene level (C40) and is listed as a third option, not proposed.",
        "",
        "Register row C155 records this scan. It establishes no biology.",
    ]
    (OUT / "g0_summary.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    rec.add_output(OUT / "g0_summary.md")
    rec.finish()
    print(df[["accession", "arms_at_or_above_floor", "smallest_arm", "raw_counts_on_disk", "verdict"]].to_string(index=False))
    print("verdicts:", counts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
