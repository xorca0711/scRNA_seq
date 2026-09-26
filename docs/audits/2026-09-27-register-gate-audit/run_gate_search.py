"""Search GEO once per register gate, and record every query and candidate verbatim.

One pass, for the questions whose cards state a data gate. A15 is excluded on purpose:
its gate was already audited in
RQ_Specified/A15_epithelial_integrin_tgfb_activation/reports/PUBLIC_DATA_SEARCH.md by the
session that owns that question, and duplicating it would waste effort and risk a
contradictory verdict.

This script only retrieves. It writes candidates to a table and makes no judgement; the
verdicts are written by hand in REPORT.md against what the table shows, and each verdict
names the accessions it rests on.

Standard library only. Refuses to overwrite.
"""
from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
TOP = 6
OUTPUTS = ["candidates.tsv", "search_run.json"]

# Each entry: question, the gate quoted from its register card, and the queries that
# express that gate as a GEO search. Queries are declared here so they are reviewable.
GATES = [
    ("A1", "matched replicated regulation-to-fate linkage remains missing", [
        "lung[Title] AND (ATAC[All Fields] OR multiome[All Fields]) AND lineage[All Fields] AND mouse[Organism]",
        "lung regeneration[All Fields] AND chromatin[All Fields] AND (tdTomato[All Fields] OR reporter[All Fields]) AND mouse[Organism]",
    ]),
    ("A3", "age-matched controls and comparable sampling needed", [
        "lung[Title] AND macrophage[All Fields] AND aging[All Fields] AND mouse[Organism] AND single cell[All Fields]",
        "lung[All Fields] AND (influenza[All Fields] OR injury[All Fields]) AND aged[All Fields] AND macrophage[All Fields] AND mouse[Organism]",
    ]),
    ("A4", "measured activity/history and lineage-linked response needed", [
        "lung[All Fields] AND Axin2[All Fields] AND mouse[Organism]",
        "lung[All Fields] AND Il1r1[All Fields] AND (lineage[All Fields] OR reporter[All Fields]) AND mouse[Organism]",
    ]),
    ("A5", "needs another eligible adult-injury cohort for generality", [
        "lung[Title] AND (bleomycin[All Fields] OR influenza[All Fields]) AND time course[All Fields] AND single cell[All Fields] AND mouse[Organism]",
        "alveolar type 2[All Fields] AND injury[All Fields] AND single cell RNA[All Fields] AND mouse[Organism]",
    ]),
    ("A7", "replicated age-matched animals; current data is one well per condition", [
        "Cebpa[All Fields] AND lung[All Fields]",
        "Cebpa[All Fields] AND (conditional[All Fields] OR knockout[All Fields]) AND single cell[All Fields]",
    ]),
    ("A8", "independent mature endpoints needed", [
        "alveolar type 1[All Fields] AND maturation[All Fields] AND lung[All Fields]",
        "lung[Title] AND alveolar[All Fields] AND (morphometry[All Fields] OR proteomics[All Fields]) AND differentiation[All Fields]",
    ]),
    ("A9", "RNA screen possible; protein/function data needed", [
        "lung fibroblast[All Fields] AND (CITE-seq[All Fields] OR surface protein[All Fields] OR antibody derived[All Fields])",
        "lung[All Fields] AND fibroblast[All Fields] AND EGFR[All Fields] AND (phosphorylation[All Fields] OR activation[All Fields])",
    ]),
    ("A11", "lack of a non-neoplastic injury comparator limits specificity claims", [
        "human[Organism] AND lung[Title] AND (pneumonia[All Fields] OR ARDS[All Fields] OR acute lung injury[All Fields]) AND single cell[All Fields] AND epithelial[All Fields]",
        "human[Organism] AND lung[All Fields] AND single cell[All Fields] AND paired[All Fields] AND (injury[All Fields] OR repair[All Fields]) AND epithelium[All Fields]",
    ]),
    ("A12", "conditional component model; activation unmeasured", [
        "IL1B[All Fields] AND lung[All Fields] AND (reporter[All Fields] OR phospho[All Fields] OR activation[All Fields])",
        "interleukin 1[All Fields] AND lung[All Fields] AND perturbation[All Fields] AND single cell[All Fields]",
    ]),
    ("A13", "IPF has six and three complete triads, below the ten-unit joint-model floor", [
        "human[Organism] AND lung[Title] AND fibrosis[All Fields] AND single cell[All Fields] AND (macrophage[All Fields] AND fibroblast[All Fields])",
        "human[Organism] AND idiopathic pulmonary fibrosis[All Fields] AND single cell[All Fields] AND atlas[All Fields]",
    ]),
    ("A14", "withdrawal and recipient-specific data needed", [
        "lung[All Fields] AND fibrosis[All Fields] AND (resolution[All Fields] OR reversible[All Fields] OR withdrawal[All Fields]) AND single cell[All Fields]",
        "lung[All Fields] AND doxycycline[All Fields] AND (off[All Fields] OR withdrawal[All Fields]) AND mouse[Organism]",
    ]),
]


def sha256(path: Path) -> str:
    d = hashlib.sha256()
    with path.open("rb") as h:
        for b in iter(lambda: h.read(1 << 24), b""):
            d.update(b)
    return d.hexdigest()


def fetch(endpoint: str, params: dict) -> dict:
    url = EUTILS + endpoint + "?" + urllib.parse.urlencode(params)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.load(r)
        except Exception as exc:  # noqa: BLE001
            if attempt == 2:
                raise
            print("   retry after %s" % type(exc).__name__, flush=True)
            time.sleep(3)
    raise RuntimeError("unreachable")


def main() -> None:
    existing = [n for n in OUTPUTS if (HERE / n).exists()]
    if existing:
        raise SystemExit("Refusing to overwrite: %s" % existing)

    rows = []
    log = []
    for question, gate, queries in GATES:
        for qi, query in enumerate(queries, start=1):
            term = query + " AND gse[Entry Type]"
            res = fetch("esearch.fcgi", {"db": "gds", "term": term, "retmax": TOP, "retmode": "json"})
            ids = res["esearchresult"]["idlist"]
            count = int(res["esearchresult"]["count"])
            log.append({"question": question, "query_index": qi, "term": term, "hits": count,
                        "ids_returned": ids})
            print("%s q%d: %d hits" % (question, qi, count), flush=True)
            if ids:
                summ = fetch("esummary.fcgi", {"db": "gds", "id": ",".join(ids), "retmode": "json"})
                for uid in ids:
                    d = summ.get("result", {}).get(uid, {})
                    rows.append({
                        "question": question, "gate": gate, "query_index": qi,
                        "accession": d.get("accession", ""),
                        "n_samples": d.get("n_samples", ""),
                        "taxon": (d.get("taxon", "") or "").replace("\t", " "),
                        "gdstype": (d.get("gdstype", "") or "").replace("\t", " "),
                        "pdat": d.get("PDAT", ""),
                        "title": (d.get("title", "") or "").replace("\t", " ").replace("\n", " "),
                    })
            time.sleep(0.4)

    with open(HERE / OUTPUTS[0], "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    record = {
        "audit": "register data-gate search",
        "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "script": "docs/audits/2026-09-27-register-gate-audit/run_gate_search.py",
        "script_sha256": sha256(Path(__file__).resolve()),
        "archive": "NCBI GEO DataSets, via eutils esearch and esummary",
        "entry_type_filter": "gse",
        "top_summaries_per_query": TOP,
        "questions_searched": [g[0] for g in GATES],
        "a15_excluded": "its gate was audited by the session that owns A15; see RQ_Specified/A15_epithelial_integrin_tgfb_activation/reports/PUBLIC_DATA_SEARCH.md",
        "queries": log,
        "candidates_recorded": len(rows),
        "judgement": "none; this script retrieves only, and REPORT.md carries the verdicts",
    }
    (HERE / OUTPUTS[1]).write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print("\nwrote %d candidates across %d queries" % (len(rows), len(log)))


if __name__ == "__main__":
    main()
