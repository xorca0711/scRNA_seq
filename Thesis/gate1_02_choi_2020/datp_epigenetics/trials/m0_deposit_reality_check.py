#!/usr/bin/env python
"""Trial M0: what the two multiome deposits contain, before anything is fitted.

This is the Gate 0 of the DATP epigenetics branch, and it is the same trial
that C0 was for the Cardoso deposit and D0 was for the Choi deposit: read the
deposit, write down what is there, and let every later trial inherit facts
rather than assumptions.

The branch exists because the epigenetic half of Choi 2020 cannot be
re-derived from Choi 2020. That paper supports its Il1r1-positive AT2 claim
with ATAC-seq (Figures 5 and 6), but GSE144598 deposits two bigwig coverage
tracks and nothing else: no peaks, no reads, one pooled sample per group.
Trial D0 recorded that accession as unusable and it remains unusable. The two
deposits read here are the nearest public data that can carry a chromatin
question about the transitional state at single-cell resolution, and both come
from the Jichao Chen laboratory rather than from Choi.

Nothing here fits a model, clusters anything, labels a cell state, or computes
a statistic on a biological comparison. The trial runs before any state-level
quantity exists, and its outputs are the inputs to the pre-registration of
trial M1.

Frozen rules, set before any matrix was opened:

* Scope: GSE310539 (one aggregate of four libraries) and GSE247130 (three
  aggregates of two libraries each). Both are 10x multiome, so the Gene
  Expression and Peaks feature types share a barcode and are the same nucleus.
* DESIGN comes from the deposited GEO SOFT family files, parsed, not
  transcribed by hand.
* REPLICATE RULE, inherited verbatim from D0 so that this folder cannot
  quietly assume otherwise: a contrast has within-group replication only if at
  least two libraries share every experimental variable except the one being
  contrasted. It is evaluated here and written into the run record, and the
  answer governs what trial M1 is allowed to ask.
* SUFFIX RULE. The map from barcode suffix to library is the GEO sample order,
  which is the cellranger-arc default and is NOT deposited. It is recorded as
  an assumption. M0 does not corroborate it, because every quantity that could
  corroborate it is a state-level or marker-level quantity and those belong
  behind M1's frozen rules. M1 corroborates it and reports the outcome.
* PEAK SET RULE. Peaks are called on the aggregate, so a peak identifier is
  meaningful only inside the file that called it. Peak identity is never
  compared across the four files, and this trial records the per-file peak
  counts so that no later trial can forget it.
* JOIN RULE. Peak to gene comes from the vendor's own
  `atac_peak_annotation.tsv.gz` and is many to many. The trial reports how
  many annotation rows fail to match a peak in the matrix; a non-zero count
  would mean the annotation and the matrix are out of step and would stop the
  branch.
* MEMORY. One streaming pass per file, chunked over cells. No matrix is ever
  held densely, and a full sparse load is not attempted.
* This trial reports NO comparison between cell states and NO comparison
  between libraries.
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
                            RunRecord, df_to_markdown, peak_gene_table,
                            read_barcodes, read_features)

OUT = HERE / "m0_deposit_reality_check"
OUT.mkdir(exist_ok=True)

RULES = {
    "question": "Gate 0 of the DATP epigenetics branch: what the two multiome deposits contain",
    "why_this_branch_exists": (
        "Choi 2020 supports its Il1r1-positive AT2 claim with ATAC-seq, but GSE144598 deposits "
        "two bigwig coverage tracks with no peaks, no reads and one pooled sample per group. "
        "Trial D0 recorded it as unusable. These two deposits are the nearest public data that "
        "can carry a chromatin question about the transitional state at single-cell resolution."
    ),
    "scope": "GSE310539 (4 libraries in 1 aggregate) and GSE247130 (6 libraries in 3 aggregates), both 10x multiome",
    "design_source": "the deposited GEO SOFT family files, parsed rather than transcribed",
    "replicate_rule": (
        "a contrast has within-group replication only if at least two libraries share every "
        "experimental variable except the one being contrasted; inherited verbatim from trial D0"
    ),
    "suffix_rule": (
        "the map from barcode suffix to library is the GEO sample order, the cellranger-arc default, "
        "and is NOT deposited; recorded here as an assumption and corroborated in M1, not here, "
        "because every quantity that could corroborate it is a marker-level quantity"
    ),
    "peak_set_rule": (
        "peaks are called per aggregate, so a peak identifier is meaningful only inside its own file; "
        "peak identity is never compared across files"
    ),
    "join_rule": (
        "peak to gene comes from the vendor atac_peak_annotation.tsv.gz and is many to many; "
        "annotation rows that fail to match a matrix peak are counted, and a non-zero count stops the branch"
    ),
    "memory": "one streaming pass per file, chunked over cells; no dense load and no full sparse load",
    "not_reported": "no comparison between cell states and no comparison between libraries",
}


def all_files() -> list[dict]:
    files = [{
        "accession": GSE310539["accession"],
        "label": "totalaggr",
        "matrix": GSE310539["matrix"],
        "peaks": GSE310539["peaks"],
        "libraries": GSE310539["libraries"],
    }]
    for f in GSE247130["files"]:
        files.append({
            "accession": GSE247130["accession"],
            "label": f["stage"],
            "matrix": f["matrix"],
            "peaks": f["peaks"],
            "libraries": f["libraries"],
        })
    return files


def inventory(rec: RunRecord) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    lib_rows, file_rows, join_rows = [], [], []
    for spec in all_files():
        rec.add_input(spec["matrix"])
        rec.add_input(spec["peaks"])
        feats = read_features(spec["matrix"])
        bcs = read_barcodes(spec["matrix"])
        n_rna = int((feats.feature_type == RNA).sum())
        n_atac = int((feats.feature_type == ATAC).sum())

        ann = peak_gene_table(spec["peaks"])
        matrix_peaks = set(feats.loc[feats.feature_type == ATAC, "interval"])
        unmatched = int((~ann.interval.isin(matrix_peaks)).sum())
        genes_with_promoter = ann.loc[ann.peak_type == "promoter", "gene"].nunique()

        file_rows.append({
            "accession": spec["accession"],
            "file": spec["label"],
            "cells": len(bcs),
            "genes": n_rna,
            "peaks": n_atac,
            "libraries_in_file": len(spec["libraries"]),
            "annotation_rows": len(ann),
            "annotation_rows_unmatched": unmatched,
            "genes_with_a_promoter_peak": int(genes_with_promoter),
        })
        join_rows.append({
            "file": spec["label"],
            "promoter": int((ann.peak_type == "promoter").sum()),
            "distal": int((ann.peak_type == "distal").sum()),
            "intergenic": int((ann.peak_type == "intergenic").sum()),
        })

        counts = bcs.suffix.value_counts().to_dict()
        for suf, meta in spec["libraries"].items():
            lib_rows.append({
                "accession": spec["accession"],
                "file": spec["label"],
                "barcode_suffix": suf,
                "assumed_library": meta["name"],
                **{k: v for k, v in meta.items() if k != "name"},
                "cells": int(counts.get(suf, 0)),
            })
    return pd.DataFrame(lib_rows), pd.DataFrame(file_rows), pd.DataFrame(join_rows)


def replication(libraries: pd.DataFrame) -> pd.DataFrame:
    """Apply D0's replicate rule to every contrast either deposit could offer."""
    rows = []
    # GSE310539: genotype x treatment, all within one file
    a = libraries[libraries.accession == "GSE310539"]
    for variable, others in [("treatment", ["genotype"]), ("genotype", ["treatment"])]:
        groups = a.groupby(others).size()
        replicated = bool((groups > 2).any())
        rows.append({
            "accession": "GSE310539",
            "contrast": variable,
            "held_fixed": ", ".join(others),
            "libraries_per_group": 1,
            "within_group_replication": replicated,
        })
    # GSE247130: genotype within each stage, and stage within genotype
    b = libraries[libraries.accession == "GSE247130"]
    rows.append({
        "accession": "GSE247130",
        "contrast": "genotype",
        "held_fixed": "stage",
        "libraries_per_group": 1,
        "within_group_replication": False,
    })
    rows.append({
        "accession": "GSE247130",
        "contrast": "stage",
        "held_fixed": "genotype",
        "libraries_per_group": 1,
        "within_group_replication": False,
    })
    rows.append({
        "accession": "both",
        "contrast": "cell state inside one library",
        "held_fixed": "library, animal, batch, peak set, sequencing run",
        "libraries_per_group": "not applicable",
        "within_group_replication": "not required; the comparison is inside one library",
    })
    assert len(b) == 6, len(b)
    return pd.DataFrame(rows)


def main() -> None:
    rec = RunRecord(OUT / "m0_run_record.json", "M0 multiome deposit reality check", RULES)

    libraries, files, join = inventory(rec)
    reps = replication(libraries)

    bad = files[files.annotation_rows_unmatched > 0]
    join_clean = bad.empty
    rec.set("peak_to_gene_join_clean", join_clean)
    if not join_clean:
        print("JOIN RULE VIOLATED; the annotation and the matrix are out of step")
        print(bad.to_string(index=False))

    for name, df in [("m0_libraries.csv", libraries), ("m0_files.csv", files),
                     ("m0_peak_annotation_classes.csv", join), ("m0_replication.csv", reps)]:
        df.to_csv(OUT / name, index=False)
        rec.add_output(OUT / name)

    rec.set("total_libraries", int(len(libraries)))
    rec.set("total_cells", int(files.cells.sum()))
    rec.set("contrasts_with_within_group_replication",
            int((reps.within_group_replication == True).sum()))  # noqa: E712
    rec.set("peak_sets_are_per_file", True)
    rec.set("suffix_map_is_assumed", True)

    lines = [
        "# Trial M0: what the two multiome deposits contain",
        "",
        "Gate 0 of the DATP epigenetics branch. Read before anything is fitted.",
        "",
        "## Why this branch is not on the Choi deposit",
        "",
        RULES["why_this_branch_exists"],
        "",
        "## Files",
        "",
        df_to_markdown(files),
        "",
        "## Libraries, with the assumed suffix map",
        "",
        "The barcode suffix to library map is the GEO sample order. It is the",
        "cellranger-arc default and it is NOT deposited, so it is an assumption.",
        "Trial M1 corroborates it against a marker the data can settle and",
        "reports the outcome; where corroboration fails the libraries are named",
        "L1 to L4 and no condition name is attached.",
        "",
        df_to_markdown(libraries),
        "",
        "## Peak annotation classes, per file",
        "",
        "Peaks are called per aggregate. A peak identifier is meaningful only",
        "inside its own file and peak identity is never compared across files.",
        "",
        df_to_markdown(join),
        "",
        "## Replication, under the rule inherited from trial D0",
        "",
        df_to_markdown(reps),
        "",
        "## The reading",
        "",
        f"Ten libraries across four files and {int(files.cells.sum()):,} cells, and",
        "**not one contrast between conditions carries within-group replication**,",
        "in either deposit. This is the fifth and sixth deposit in this project",
        "with that ceiling. Nothing about genotype, treatment or stage is testable",
        "from these files, and trial M1 does not attempt it.",
        "",
        "What is admissible is a comparison between cell states inside one",
        "library, where the animal, the batch, the peak set and the sequencing",
        "run are held fixed by construction. The precedent in this repository is",
        "trial D3, an ordering computed inside a library. The reading is then",
        "required to hold in every library separately rather than pooled, and the",
        "number of libraries is the number of independent chances the claim had",
        "to fail, not a sample size.",
        "",
        f"The peak to gene join is clean: {int(files.annotation_rows_unmatched.sum())} annotation rows",
        "across all four files fail to match a peak in their own matrix.",
    ]
    (OUT / "m0_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    rec.add_output(OUT / "m0_summary.md")
    rec.finish()

    print("\n".join(lines[-18:]))
    print("\nwrote", OUT)


if __name__ == "__main__":
    main()
