"""Shared helpers for the Choi-2020 trials.

Everything generic comes from the repository's one shared helper module, the
same one the Cardoso trials use: the run record, the MatrixMarket triplet
reader, the GEO SOFT family parser, the non-Ensembl feature rule and the QC
metric helper. Nothing is redefined here.

This module adds only what is specific to this deposit: where its libraries
live, and the fact that its matrices are unfiltered 10x barcode whitelists
rather than called cells, which every trial in this folder has to handle.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RAW = REPO / "raw_data"

_SHARED_TRIALS = REPO / "Thesis" / "gate1_04_sikkema_2023_hlca" / "trials"
if str(_SHARED_TRIALS) not in sys.path:
    sys.path.insert(0, str(_SHARED_TRIALS))
from trial_utils import (ENSEMBL_ID_RE, RunRecord, df_to_markdown,  # noqa: E402,F401
                         file_facts, package_versions, parse_soft, qc_metrics,
                         read_mtx_triplet, shannon, utc_now)

# The 10x v2 barcode whitelist. A deposited matrix with this many columns is
# the raw droplet matrix, not called cells, so cell calling is this
# repository's job and the cell count will not match the paper's unless the
# same caller and version are used.
TENX_V2_WHITELIST = 737_280

LINEAGE_LIBRARIES = [
    {"gsm": "GSM4304609", "library": "PBS_AT2_Tomato", "timepoint": "PBS", "sort": "Tomato"},
    {"gsm": "GSM4304610", "library": "PBS_AT2_nonTomato", "timepoint": "PBS", "sort": "nonTomato"},
    {"gsm": "GSM4304611", "library": "Day14_AT2_Tomato", "timepoint": "day 14", "sort": "Tomato"},
    {"gsm": "GSM4304612", "library": "Day14_AT2_nonTomato", "timepoint": "day 14", "sort": "nonTomato"},
    {"gsm": "GSM4304613", "library": "Day28_AT2_Tomato", "timepoint": "day 28", "sort": "Tomato"},
    {"gsm": "GSM4304614", "library": "Day28_AT2_nonTomato", "timepoint": "day 28", "sort": "nonTomato"},
]
ORGANOID_LIBRARIES = [
    {"gsm": "GSM4288824", "library": "Control_Organoids", "treatment": "control"},
    {"gsm": "GSM4288825", "library": "Il1b_Organoids", "treatment": "IL-1beta"},
]
ACCESSIONS = {
    "GSE145031": {"role": "scRNA-seq of AT2 lineage-traced epithelium",
                  "libraries": LINEAGE_LIBRARIES},
    "GSE144468": {"role": "scRNA-seq of AT2 organoids",
                  "libraries": ORGANOID_LIBRARIES},
}


def library_paths(accession: str, library: str, gsm: str) -> dict[str, Path]:
    """Where this deposit's triplet files sit after the series tar is unpacked."""
    root = RAW / accession / (accession + "_RAW")
    stem = gsm + "_" + library
    return {"matrix": root / (stem + ".mtx.gz"),
            "features": root / (stem + "_gene.tsv.gz"),
            "barcodes": root / (stem + "_barcodes.tsv.gz")}
