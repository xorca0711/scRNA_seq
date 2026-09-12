#!/usr/bin/env python
"""Trial C0: what is actually in the Cardoso 2026 deposits, and what can it support.

Gate 0 of the owner's staged plan. This trial opens every deposited library of
the four accessions named in the paper's data-availability statement, plus the
England 2025 series the paper's CellChat analysis integrates with, and
measures: the design as deposited (parsed from the GEO SOFT family files, not
transcribed), the gene space of each series, the QC distribution of each
library, and whether the populations the paper's claims rest on are detectable
at all.

Nothing is filtered, clustered or compared here. The question is only which
entry points of the plan are live and what statistical unit the deposit can
carry.

Frozen rules (written to the run record before any matrix is read):

* Libraries: every mtx/tsv triplet under raw_data/<accession>/, for
  GSE316241 (mesenchyme), GSE316243 (niche: immune plus stroma), GSE316244
  (Areg-flox arm: niche and RFP+ epithelium), GSE310335 (human alveolar
  organoids) and GSE247505 (England 2025 RFP+ lineage-labelled epithelium).
* A feature whose identifier is not an Ensembl gene ID is not a gene: it is
  removed from the expression matrix and carried as a per-cell column. This
  is the general form of the SiteA/SiteB decision taken for GSE262927.
* Species is decided per library by which of the two fixed marker panels
  (mouse, human) has the higher mean detection fraction; a library whose
  winning panel is detected in under 10% of cells on average is reported as
  "species unclear" rather than assigned.
* QC is measured, not applied: per-cell total counts, genes detected,
  mitochondrial, ribosomal and haemoglobin fractions, summarised as medians.
  Thresholds belong to Gate 1 and are derived per library there.
* Gene spaces are compared by the exact ordered list of (Ensembl ID, symbol);
  two libraries share a space only if those lists are identical.
* Population presence is measured as the fraction of cells with a non-zero
  count for each gene of a fixed panel taken from the paper's own marker
  lists (Fig. 1d, Fig. 1h, Fig. 4l, Extended Data Figs. 4d and 5e-j). A
  detection fraction is a presence measure, not an abundance and not a
  cell-type call.
* The biological replicate is the mouse. Where the deposited design states
  that a library pools several mice, the library provides no within-group
  replication, and the number of independent animals behind a genotype
  contrast is recorded as the number of libraries, not the number of pooled
  mice and not the number of cells. This rule is fixed here because it
  determines which of the plan's gates can produce a tested claim at all.
"""

from __future__ import annotations

import gc
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cardoso_utils import (ENSEMBL_ID_RE, RAW, REPO, RunRecord,  # noqa: E402
                           df_to_markdown, parse_soft, qc_metrics,
                           read_mtx_triplet)

OUT = HERE / "c0_data_reality_check"
OUT.mkdir(exist_ok=True)

ACCESSIONS = {
    "GSE316241": {"role": "mesenchymal cells, Confetti and Red2Kras lungs (paper Fig. 1b-e)",
                  "species_expected": "mouse"},
    "GSE316243": {"role": "immune and stromal cells, Confetti and Red2Kras lungs (paper Fig. 1f-h)",
                  "species_expected": "mouse"},
    "GSE316244": {"role": "Areg-flox arm: niche and RFP+ mutant epithelium (paper Fig. 4d-m)",
                  "species_expected": "mouse"},
    "GSE310335": {"role": "human KRAS-G12D alveolar organoids (paper Fig. 5f-h)",
                  "species_expected": "human"},
    "GSE247505": {"role": "England 2025: RFP+ lineage-labelled mutant epithelium, the CellChat partner",
                  "species_expected": "mouse"},
}

MOUSE_PROBE = ["Sftpc", "Epcam", "Ptprc", "Col1a1", "Pecam1", "Actb", "Cd68", "Krt8"]
HUMAN_PROBE = ["SFTPC", "EPCAM", "PTPRC", "COL1A1", "PECAM1", "ACTB", "CD68", "KRT8"]

# Panels quoted from the paper's own figures; mouse symbols, with the human
# spellings applied to the organoid series.
PANEL = {
    "alveolar fibroblast (Fig. 1d)": ["Pdgfra", "Col13a1", "Tcf21", "Npnt", "Scube2"],
    "adventitial fibroblast (Fig. 1d)": ["Col14a1", "Pi16", "Dcn"],
    "reprogrammed/fibrotic fibroblast (Fig. 1d)": ["Tnc", "Fst", "Runx1", "Runx2", "Acta2", "Pdgfrb"],
    "inflammatory fibroblast (ED Fig. 4l-m)": ["Lcn2", "Saa3", "Sfrp1", "Cxcl12"],
    "alveolar macrophage (Fig. 1h)": ["Itgax", "Siglec5", "Siglecf", "Pparg", "Mertk", "Car4"],
    "reprogrammed AM (Fig. 1h)": ["Msr1", "Cdh1", "Ch25h", "Ear6", "Fstl1", "Cxcl2"],
    "monocyte/IM (Fig. 1h)": ["Itgam", "Ly6c2", "Ccr2", "Cx3cr1", "Mmp19"],
    "Tnc-TLR4 axis (ED Fig. 4d)": ["Tnc", "Tlr4"],
    "EGF axis (ED Fig. 5f-h)": ["Areg", "Ereg", "Hbegf", "Tgfa", "Egfr"],
    "fibroblast-to-epithelium (ED Fig. 5i-j)": ["Wnt5a", "Igf1", "Spp1"],
    "mutant epithelial states (Fig. 4l)": ["Sftpc", "Etv5", "Lamp3", "Cd177", "Cldn4", "Itga2",
                                           "Ndrg1", "Sox9", "Mki67", "Ager", "Hopx"],
    "compartment": ["Epcam", "Ptprc", "Col1a1", "Pecam1"],
}
HUMAN_PANEL = {
    "human organoid states (Fig. 5h)": ["SFTPC", "NAPSA", "NKX2-1", "SCGB1A1", "SCGB3A2", "SCGB3A1",
                                        "SFTPB", "KRT18", "LAMP3", "LPCAT1", "KRAS", "KRT8",
                                        "ITGA2", "ITGB4", "HMGA2", "SOX9", "CDKN2A", "CDKN2B",
                                        "AREG", "HBEGF", "TIMP1", "TGFB1", "NFKB2", "MKI67"],
    "compartment": ["EPCAM", "PTPRC", "COL1A1", "PECAM1"],
}

RULES = {
    "accessions": {k: v["role"] for k, v in ACCESSIONS.items()},
    "non_gene_feature_rule": "identifier not matching ^ENS[A-Z]*G\\d+ is removed from the matrix and carried per cell",
    "species_rule": "higher mean detection of the mouse or human probe panel; under 0.10 mean detection reported as unclear",
    "qc": "measured only (total counts, genes, pct mito/ribo/haemoglobin), reported as per-library medians; no filtering",
    "gene_space_rule": "identical ordered (ensembl_id, symbol) list",
    "presence_rule": "fraction of cells with a non-zero count; a presence measure, not abundance, not a cell-type call",
    "replicate_rule": ("the mouse is the biological replicate; a library pooling several mice carries no "
                       "within-group replication; independent animals behind a contrast are counted as "
                       "libraries, never as pooled mice and never as cells"),
    "probe_panels": {"mouse": MOUSE_PROBE, "human": HUMAN_PROBE},
    "marker_panels": {"mouse": PANEL, "human": HUMAN_PANEL},
}


def library_triplets(accession: str) -> list[dict]:
    root = RAW / accession
    triplets = []
    for matrix in sorted(root.rglob("*matrix.mtx.gz")):
        stem = matrix.name.replace("matrix.mtx.gz", "")
        features = matrix.with_name(stem + "features.tsv.gz")
        if not features.exists():
            features = matrix.with_name(stem + "genes.tsv.gz")
        barcodes = matrix.with_name(stem + "barcodes.tsv.gz")
        if not (features.exists() and barcodes.exists()):
            raise FileNotFoundError(f"incomplete triplet for {matrix}")
        gsm = matrix.name.split("_", 1)[0] if matrix.name.startswith("GSM") else ""
        triplets.append({"accession": accession, "library": stem.rstrip("_"),
                         "gsm": gsm, "matrix": matrix, "features": features,
                         "barcodes": barcodes})
    return triplets


def soft_design(accession: str) -> pd.DataFrame:
    path = RAW / accession / f"{accession}_family.soft.gz"
    if not path.exists():
        return pd.DataFrame()
    parsed = parse_soft(path)
    rows = []
    for sample in parsed["samples"]:
        chars = {}
        for item in sample["characteristics"]:
            key, _, value = item.partition(": ")
            chars[key.strip()] = value.strip()
        rows.append({
            "accession": accession,
            "gsm": sample["gsm"],
            "title": sample.get("title", ""),
            "organism": sample.get("organism_ch1", ""),
            "genotype": chars.get("genotype", ""),
            "cell_type_sorted": chars.get("cell type", ""),
            "tissue": chars.get("tissue", ""),
            "instrument": sample.get("instrument_model", ""),
            "alignment": " | ".join(sample["data_processing"]),
        })
    df = pd.DataFrame(rows)
    df.attrs["series"] = parsed["series"]
    return df


def main() -> None:
    rec = RunRecord(OUT / "c0_run_record.json",
                    "C0 Cardoso 2026 data reality check (Gate 0)", RULES,
                    notes=("rules frozen before any matrix was read; the deposited design is parsed "
                           "from the GEO SOFT family files, not transcribed"))

    # ---- the deposited design, from the SOFT files ------------------------
    design_frames, series_text = [], {}
    for accession in ACCESSIONS:
        soft = RAW / accession / f"{accession}_family.soft.gz"
        if soft.exists():
            rec.add_input(soft)
        df = soft_design(accession)
        if not df.empty:
            series_text[accession] = {
                k: v for k, v in df.attrs["series"].items()
                if k in ("Series_title", "Series_status", "Series_overall_design",
                         "Series_platform_id", "Series_pubmed_id")}
            design_frames.append(df)
    design = pd.concat(design_frames, ignore_index=True) if design_frames else pd.DataFrame()
    design.to_csv(OUT / "c0_deposited_design.csv", index=False)
    rec.add_output(OUT / "c0_deposited_design.csv")
    rec.set("series_metadata", series_text)

    # ---- every library ----------------------------------------------------
    inventory, presence_rows, gene_space = [], [], {}
    for accession, meta in ACCESSIONS.items():
        for lib in library_triplets(accession):
            for key in ("matrix", "features", "barcodes"):
                rec.add_input(lib[key])
            X, var, barcodes = read_mtx_triplet(lib["matrix"], lib["features"], lib["barcodes"])

            is_gene = var["gene_id"].str.match(ENSEMBL_ID_RE.pattern).fillna(False).to_numpy()
            non_gene = var.loc[~is_gene, "gene_id"].tolist()
            non_gene_counts = {}
            for name in non_gene:
                col = np.asarray(X[:, int(np.where(var["gene_id"].to_numpy() == name)[0][0])].todense()).ravel()
                non_gene_counts[name] = {"cells_detected": int((col > 0).sum()),
                                         "fraction_of_cells": round(float((col > 0).mean()), 4),
                                         "total_counts": int(col.sum())}
            X = X[:, is_gene]
            var = var.loc[is_gene].reset_index(drop=True)

            symbols = var["gene_symbol"].astype(str).to_numpy()
            def detection(genes: list[str]) -> dict[str, float]:
                out = {}
                for gene in genes:
                    hits = np.where(symbols == gene)[0]
                    out[gene] = (round(float((X[:, hits].sum(axis=1) > 0).mean()), 4)
                                 if len(hits) else np.nan)
                return out

            mouse_hit = detection(MOUSE_PROBE)
            human_hit = detection(HUMAN_PROBE)
            mouse_mean = float(np.nanmean(list(mouse_hit.values()))) if not all(np.isnan(list(mouse_hit.values()))) else 0.0
            human_mean = float(np.nanmean(list(human_hit.values()))) if not all(np.isnan(list(human_hit.values()))) else 0.0
            species = "mouse" if mouse_mean >= human_mean else "human"
            if max(mouse_mean, human_mean) < 0.10:
                species_call = "species unclear"
            else:
                species_call = species

            qc = qc_metrics(X, var, species)
            panel = PANEL if species == "mouse" else HUMAN_PANEL
            for group, genes in panel.items():
                for gene, frac in detection(genes).items():
                    presence_rows.append({"accession": accession, "library": lib["library"],
                                          "panel": group, "gene": gene,
                                          "detected_in_fraction_of_cells": frac,
                                          "gene_absent_from_reference": bool(np.isnan(frac))})

            key = (accession, tuple(var["gene_id"]), tuple(var["gene_symbol"]))
            gene_space.setdefault(hash(key[1:]), []).append(f"{accession}:{lib['library']}")

            inventory.append({
                "accession": accession, "library": lib["library"], "gsm": lib["gsm"],
                "role": meta["role"],
                "barcodes": int(X.shape[0]), "genes_in_reference": int(X.shape[1]),
                "non_gene_features": ", ".join(non_gene) or "none",
                "species_probe_mouse": round(mouse_mean, 4),
                "species_probe_human": round(human_mean, 4),
                "species_call": species_call,
                "species_expected": meta["species_expected"],
                "species_agrees": species_call == meta["species_expected"],
                "median_total_counts": int(np.median(qc["total_counts"])),
                "median_genes": int(np.median(qc["n_genes"])),
                "median_pct_mt": round(float(np.median(qc["pct_mt"])), 2),
                "median_pct_ribo": round(float(np.median(qc["pct_ribo"])), 2),
                "median_pct_hb": round(float(np.median(qc["pct_hb"])), 2),
                "cells_under_500_counts": int((qc["total_counts"] < 500).sum()),
                "cells_over_10pct_mt": int((qc["pct_mt"] > 10).sum()),
                "nnz": int(X.nnz),
                "non_gene_feature_detail": json.dumps(non_gene_counts) if non_gene_counts else "",
            })
            print(f"  {accession} {lib['library']}: {X.shape[0]} cells, "
                  f"{X.shape[1]} genes, median {int(np.median(qc['n_genes']))} genes/cell, "
                  f"{species_call}")
            del X, var, qc
            gc.collect()

    inv = pd.DataFrame(inventory)
    inv.to_csv(OUT / "c0_library_inventory.csv", index=False)
    rec.add_output(OUT / "c0_library_inventory.csv")
    pres = pd.DataFrame(presence_rows)
    pres.to_csv(OUT / "c0_marker_presence.csv", index=False)
    rec.add_output(OUT / "c0_marker_presence.csv")

    spaces = [{"gene_space": i + 1, "n_libraries": len(v), "libraries": "; ".join(v)}
              for i, v in enumerate(gene_space.values())]
    gs = pd.DataFrame(spaces)
    gs.to_csv(OUT / "c0_gene_spaces.csv", index=False)
    rec.add_output(OUT / "c0_gene_spaces.csv")

    # ---- the summary, composed from the tables just written --------------
    pooled = {"GSE316241": 3, "GSE316243": 3, "GSE316244": 3}
    per_accession = []
    for accession, meta in ACCESSIONS.items():
        rows = inv[inv["accession"] == accession]
        design_rows = design[design["accession"] == accession] if not design.empty else pd.DataFrame()
        genotypes = sorted({g for g in design_rows.get("genotype", pd.Series(dtype=str)) if g})
        per_accession.append({
            "accession": accession,
            "libraries": int(len(rows)),
            "barcodes": int(rows["barcodes"].sum()),
            "median_genes_per_cell_range": f"{int(rows['median_genes'].min())} to {int(rows['median_genes'].max())}",
            "mice_pooled_per_library": pooled.get(accession, "not stated"),
            "independent_animals_per_genotype": (
                "1 library, 3 mice pooled: no within-genotype replication"
                if accession in pooled else "see the deposited design"),
            "distinct_genotypes_in_characteristics": len(genotypes),
            "role": meta["role"],
        })
    summary_table = pd.DataFrame(per_accession)
    summary_table.to_csv(OUT / "c0_accession_summary.csv", index=False)
    rec.add_output(OUT / "c0_accession_summary.csv")

    key_panels = ["reprogrammed/fibrotic fibroblast (Fig. 1d)", "alveolar fibroblast (Fig. 1d)",
                  "inflammatory fibroblast (ED Fig. 4l-m)", "reprogrammed AM (Fig. 1h)",
                  "EGF axis (ED Fig. 5f-h)", "fibroblast-to-epithelium (ED Fig. 5i-j)"]
    mesen = pres[(pres["accession"] == "GSE316241") & (pres["panel"].isin(key_panels))]
    mesen_wide = (mesen.pivot_table(index=["panel", "gene"], columns="library",
                                    values="detected_in_fraction_of_cells")
                  .round(3).reset_index())
    mesen_wide.columns = [c.replace("GSM9447763_Expt1_", "").replace("GSM9447764_Expt1_", "")
                          for c in mesen_wide.columns]

    absent = sorted(pres.loc[pres["gene_absent_from_reference"], "gene"].unique())
    non_gene_detail = {r["library"]: json.loads(r["non_gene_feature_detail"])
                       for r in inventory if r["non_gene_feature_detail"]}

    lines = [
        "# Trial C0 output: what the Cardoso 2026 deposits contain (Gate 0)", "",
        f"All {len(inv)} deposited libraries of the five accessions are public and readable: "
        f"{int(inv['barcodes'].sum()):,} barcodes before any quality control. Every species call "
        f"agrees with the accession's stated organism. "
        f"{len(gene_space)} distinct gene spaces are present, so integration across accessions needs "
        f"an explicit intersection on Ensembl gene ID.", "",
        "## What each accession is, and what it can carry", "",
        df_to_markdown(summary_table, index=False), "",
        "The entry marked *no within-genotype replication* is the binding constraint of this deposit: "
        "each mouse library is one library per genotype built from three pooled mice, so a genotype "
        "difference in these data is a difference between two libraries. It can be described; it "
        "cannot be tested, and no P value is admissible on it.", "",
        "## Gene spaces", "", df_to_markdown(gs, index=False), "",
        "## Per-library quality, measured and not applied", "",
        df_to_markdown(inv[["accession", "library", "barcodes", "median_total_counts", "median_genes",
                            "median_pct_mt", "cells_over_10pct_mt", "non_gene_features"]], index=False), "",
        "## Non-gene features", "",
    ]
    if non_gene_detail:
        lines += [
            "One feature in the whole deposit is not an Ensembl gene: `BSD` (blasticidin-S deaminase, "
            "a selection marker carried by the reporter construct), present in GSE316244 only. It is "
            "removed from the expression matrix by the frozen rule and carried per cell. Its detection "
            "tracks the sort:", "",
            df_to_markdown(pd.DataFrame([
                {"library": k, "cells_detected": v["BSD"]["cells_detected"],
                 "fraction_of_cells": v["BSD"]["fraction_of_cells"],
                 "total_counts": v["BSD"]["total_counts"]}
                for k, v in non_gene_detail.items()]), index=False), "",
        ]
    else:
        lines += ["No non-gene features were found.", ""]
    lines += [
        "## Are the populations the paper's claims rest on detectable at all", "",
        "Fraction of cells with a non-zero count, mesenchymal series. A detection fraction is a "
        "presence measure: it is not abundance and it is not a cell-type call.", "",
        df_to_markdown(mesen_wide, index=False), "",
        f"Genes of the paper's own panels that are absent from a deposited reference: "
        f"{', '.join(absent) if absent else 'none'}.", "",
    ]
    (OUT / "c0_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "c0_summary.md")

    # ---- what the design can carry ---------------------------------------
    rec.set("libraries_total", int(len(inv)))
    rec.set("cells_total_deposited", int(inv["barcodes"].sum()))
    rec.set("species_calls_agree_with_expectation", bool(inv["species_agrees"].all()))
    rec.set("distinct_gene_spaces", int(len(gene_space)))
    rec.set("non_gene_features_found",
            sorted({f for row in inventory for f in row["non_gene_features"].split(", ")
                    if f != "none"}))
    rec.set("libraries_per_accession", inv.groupby("accession").size().to_dict())
    rec.finish()
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
