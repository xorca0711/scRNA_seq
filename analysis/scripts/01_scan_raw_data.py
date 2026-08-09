#!/usr/bin/env python
"""
PHASE 1 - Recursive scan and structural inventory of the raw-data directory.

This script makes no assumption about what kind of single-cell data are
present.  It walks the tree, classifies every file by extension *and* by
actually opening it, and writes a structured inventory that the rest of the
pipeline reads to decide how to load the data.

The raw-data directory is opened read-only and is never modified.

Outputs
-------
analysis/raw_data_inventory.txt   human-readable report + directory tree
analysis/raw_data_inventory.csv   one row per file
analysis/raw_data_inventory.json  machine-readable version of the same scan
"""

from __future__ import annotations

import gzip
import json
import os
import sys
import tarfile
from dataclasses import dataclass, field, asdict
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RAW = REPO / "raw_data"
ANALYSIS = REPO / "analysis"

# Marker panels used only to *test* species, never to convert symbols.
MOUSE_PROBES = ["Sftpc", "Krt8", "Trp63", "Ager", "Scgb1a1", "Foxj1",
                "Cd3d", "Lyz2", "Col1a1", "Pecam1", "Actb", "Epcam"]
HUMAN_PROBES = ["SFTPC", "KRT8", "TP63", "AGER", "SCGB1A1", "FOXJ1",
                "CD3D", "LYZ", "COL1A1", "PECAM1", "ACTB", "EPCAM"]


# --------------------------------------------------------------------------
# record type
# --------------------------------------------------------------------------
@dataclass
class FileRecord:
    path: str
    filename: str
    parent: str
    ext: str
    size_bytes: int
    size_human: str
    file_type: str = "unknown"
    inferred_role: str = "unknown"
    dataset: str = ""
    sample: str = ""
    dimensions: str = ""
    contains_raw_counts: str = "unknown"
    contains_metadata: str = "no"
    appears_processed: str = "unknown"
    species_evidence: str = ""
    notes: str = ""
    extra: dict = field(default_factory=dict)


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n:.1f} TB"


# --------------------------------------------------------------------------
# format probes
# --------------------------------------------------------------------------
def probe_hdf5(p: Path, rec: FileRecord) -> None:
    """Identify 10x CellRanger HDF5 (v2 or v3), h5ad, loom, or generic HDF5."""
    try:
        import h5py
    except ImportError:
        rec.notes = "h5py unavailable"
        return

    try:
        with h5py.File(p, "r") as f:
            top = list(f.keys())
            rec.extra["hdf5_top_level_keys"] = top

            # AnnData (.h5ad) written through HDF5
            if {"X", "obs", "var"}.issubset(set(top)):
                rec.file_type = "h5ad (AnnData)"
                rec.inferred_role = "processed AnnData object"
                rec.appears_processed = "yes"
                return

            # loom
            if "matrix" in top and "row_attrs" in top and "col_attrs" in top:
                rec.file_type = "loom"
                rec.inferred_role = "loom expression object"
                return

            # 10x CellRanger v3: /matrix/{data,indices,indptr,shape,barcodes,features}
            group_name = None
            if "matrix" in top and isinstance(f["matrix"], h5py.Group) \
                    and "shape" in f["matrix"]:
                group_name = "matrix"
                rec.file_type = "10x CellRanger HDF5 (v3 layout)"
            else:
                # 10x CellRanger v2: /<genome>/{data,indices,indptr,shape,barcodes,genes}
                for k in top:
                    if isinstance(f[k], h5py.Group) and "shape" in f[k]:
                        group_name = k
                        rec.file_type = "10x CellRanger HDF5 (v2 genome-group layout)"
                        rec.extra["genome_group"] = k
                        break

            if group_name is None:
                rec.file_type = "HDF5 (unrecognised layout)"
                return

            g = f[group_name]
            shape = tuple(int(x) for x in g["shape"][:])
            n_genes, n_cells = shape
            rec.dimensions = f"{n_genes} genes x {n_cells} barcodes"
            rec.extra["n_genes"] = n_genes
            rec.extra["n_barcodes"] = n_cells

            # dtype / integrality -> raw counts?
            data = g["data"]
            rec.extra["data_dtype"] = str(data.dtype)
            probe_n = int(min(len(data), 100_000))
            if probe_n:
                vals = data[:probe_n]
                is_int = bool(np.all(np.equal(np.mod(vals, 1), 0)))
                rec.contains_raw_counts = "yes (integer values)" if is_int \
                    else "no (non-integer values)"
                rec.extra["min_value"] = float(np.min(vals))
                rec.extra["max_value_in_probe"] = float(np.max(vals))
                rec.appears_processed = "no (looks like raw counts)" if is_int else "yes"

            # gene / feature names -> species evidence
            names = None
            if "features" in g and "name" in g["features"]:
                names = g["features"]["name"][:]
                if "genome" in g["features"]:
                    genomes = np.unique(g["features"]["genome"][:])
                    rec.extra["genome"] = [x.decode() for x in genomes]
                if "feature_type" in g["features"]:
                    ft, cnt = np.unique(g["features"]["feature_type"][:],
                                        return_counts=True)
                    rec.extra["feature_types"] = {
                        k.decode(): int(v) for k, v in zip(ft, cnt)}
            elif "gene_names" in g:
                names = g["gene_names"][:]

            if names is not None:
                syms = [x.decode() if isinstance(x, bytes) else str(x)
                        for x in names]
                sset = set(syms)
                mouse_hits = [m for m in MOUSE_PROBES if m in sset]
                human_hits = [h for h in HUMAN_PROBES if h in sset]
                mt_mouse = [s for s in syms if s.startswith("mt-")]
                mt_human = [s for s in syms if s.startswith("MT-")]
                rec.extra["mouse_probe_hits"] = mouse_hits
                rec.extra["human_probe_hits"] = human_hits
                rec.extra["n_mt_mouse_style"] = len(mt_mouse)
                rec.extra["n_mt_human_style"] = len(mt_human)
                if len(mouse_hits) > len(human_hits):
                    rec.species_evidence = (
                        f"MOUSE ({len(mouse_hits)}/{len(MOUSE_PROBES)} "
                        f"title-case markers, {len(mt_mouse)} 'mt-' genes)")
                elif len(human_hits) > len(mouse_hits):
                    rec.species_evidence = (
                        f"HUMAN ({len(human_hits)}/{len(HUMAN_PROBES)} "
                        f"upper-case markers, {len(mt_human)} 'MT-' genes)")
                else:
                    rec.species_evidence = "ambiguous"
                rec.extra["n_duplicate_gene_symbols"] = int(len(syms) - len(sset))

            # raw droplet matrix vs filtered cell matrix
            fname = p.name.lower()
            if "raw_feature_bc" in fname or n_cells > 200_000:
                rec.inferred_role = "10x RAW (unfiltered droplet) count matrix"
            elif "filtered_feature_bc" in fname:
                rec.inferred_role = "10x FILTERED (cell-called) count matrix"
            else:
                rec.inferred_role = ("10x count matrix, cell-called by barcode "
                                     "count" if n_cells < 100_000 else
                                     "10x count matrix (unfiltered?)")
    except Exception as exc:  # noqa: BLE001
        rec.notes = f"probe failed: {type(exc).__name__}: {exc}"


def probe_table(p: Path, rec: FileRecord) -> None:
    """Decide whether a delimited text file is an expression matrix or metadata."""
    opener = gzip.open if p.suffix == ".gz" else open
    try:
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            head = [fh.readline() for _ in range(6)]
        head = [h for h in head if h]
        if not head:
            rec.notes = "empty file"
            return
        sep = "\t" if head[0].count("\t") > head[0].count(",") else ","
        cols = head[0].rstrip("\n").split(sep)
        rec.extra["n_columns"] = len(cols)
        rec.extra["columns_preview"] = cols[:40]
        rec.extra["delimiter"] = "tab" if sep == "\t" else "comma"

        # count rows cheaply
        with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
            n_rows = sum(1 for _ in fh) - 1
        rec.extra["n_data_rows"] = n_rows
        rec.dimensions = f"{n_rows} rows x {len(cols)} columns"

        # An expression matrix has many numeric columns whose names look like
        # barcodes/samples; metadata has few columns with descriptive names.
        body = [ln.rstrip("\n").split(sep) for ln in head[1:]]
        numeric_frac = 0.0
        if body:
            row = body[0]
            n_num = 0
            for v in row[1:]:
                try:
                    float(v)
                    n_num += 1
                except ValueError:
                    pass
            numeric_frac = n_num / max(1, len(row) - 1)
        rec.extra["fraction_numeric_in_first_row"] = round(numeric_frac, 3)

        looks_matrix = len(cols) > 200 and numeric_frac > 0.9
        if looks_matrix:
            rec.file_type = f"delimited text ({rec.extra['delimiter']})"
            rec.inferred_role = "possible dense expression matrix"
            rec.contains_raw_counts = "candidate - verify orientation"
        else:
            rec.file_type = f"delimited text ({rec.extra['delimiter']})"
            rec.inferred_role = "metadata / annotation table"
            rec.contains_metadata = "yes"
            rec.contains_raw_counts = "no"
            rec.appears_processed = "yes (derived annotation)"
    except Exception as exc:  # noqa: BLE001
        rec.notes = f"probe failed: {type(exc).__name__}: {exc}"


def probe_rds(p: Path, rec: FileRecord) -> None:
    """Read the RDS magic header without needing R."""
    try:
        with open(p, "rb") as fh:
            magic = fh.read(8)
        rec.extra["magic_bytes"] = magic.hex()
        if magic[:2] == b"\x1f\x8b":
            rec.file_type = "RDS (gzip-compressed R serialization)"
        elif magic[:2] in (b"XX", b"X\n") or magic[:1] == b"X":
            rec.file_type = "RDS (uncompressed XDR R serialization)"
        elif magic[:3] == b"\xfd7z":
            rec.file_type = "RDS (xz-compressed R serialization)"
        elif magic[:3] == b"BZh":
            rec.file_type = "RDS (bzip2-compressed R serialization)"
        else:
            rec.file_type = "RDS (unrecognised compression)"
        rec.inferred_role = "pre-processed Seurat object (author-supplied)"
        rec.appears_processed = "yes"
        rec.contains_raw_counts = "unknown without R"
        rec.contains_metadata = "likely (Seurat meta.data)"
        rec.notes = ("requires R/Seurat to read; R is not installed in this "
                     "environment")
    except Exception as exc:  # noqa: BLE001
        rec.notes = f"probe failed: {type(exc).__name__}: {exc}"


def probe_tar(p: Path, rec: FileRecord) -> None:
    try:
        with tarfile.open(p, "r:*") as tf:
            members = tf.getnames()
        rec.file_type = "tar archive"
        rec.inferred_role = "GEO supplementary archive (already extracted)"
        rec.extra["n_members"] = len(members)
        rec.extra["members"] = members[:60]
        rec.dimensions = f"{len(members)} members"
        rec.contains_raw_counts = "see extracted members"
    except Exception as exc:  # noqa: BLE001
        rec.notes = f"probe failed: {type(exc).__name__}: {exc}"


# --------------------------------------------------------------------------
# scan
# --------------------------------------------------------------------------
def classify(p: Path) -> FileRecord:
    st = p.stat()
    rec = FileRecord(
        # Persist repository-relative paths so the tracked inventory is
        # portable and does not disclose the analyst's local checkout path.
        path=p.relative_to(REPO).as_posix(),
        filename=p.name,
        parent=p.parent.relative_to(RAW).as_posix() if p.parent != RAW else ".",
        ext="".join(p.suffixes[-2:]) if p.name.endswith(".gz") else p.suffix,
        size_bytes=st.st_size,
        size_human=human(st.st_size),
    )
    # dataset = first path component under raw_data
    rel = p.relative_to(RAW)
    rec.dataset = rel.parts[0] if len(rel.parts) > 1 else ""

    name = p.name
    low = name.lower()

    if low.endswith(".h5") or low.endswith(".hdf5"):
        probe_hdf5(p, rec)
    elif low.endswith(".h5ad"):
        probe_hdf5(p, rec)
    elif low.endswith(".loom"):
        probe_hdf5(p, rec)
    elif low.endswith(".rds") or low.endswith(".rds.gz"):
        if low.endswith(".rds.gz"):
            rec.file_type = "gzip-compressed RDS"
            rec.inferred_role = ("compressed copy of a pre-processed Seurat "
                                 "object")
            rec.appears_processed = "yes"
            rec.notes = "duplicate of the extracted .RDS in the parent folder"
        else:
            probe_rds(p, rec)
    elif low.endswith(".tar"):
        probe_tar(p, rec)
    elif low.endswith((".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz")):
        probe_table(p, rec)
    elif low.endswith((".mtx", ".mtx.gz")):
        rec.file_type = "MatrixMarket"
        rec.inferred_role = "10x sparse count matrix"
        rec.contains_raw_counts = "likely"
    elif low.endswith((".xlsx", ".xls")):
        rec.file_type = "Excel workbook"
        rec.inferred_role = "possible sample sheet / metadata"
        rec.contains_metadata = "likely"
    else:
        rec.file_type = f"other ({rec.ext or 'no extension'})"

    # sample id from the GEO GSM naming convention
    if name.startswith("GSM"):
        stem = name.split(".")[0]
        parts = stem.split("_", 1)
        rec.extra["gsm"] = parts[0]
        if len(parts) > 1:
            tail = parts[1]
            for suffix in ("_filtered_feature_bc_matrix",
                           "filtered_feature_bc_matrix",
                           "_raw_feature_bc_matrix",
                           "raw_feature_bc_matrix"):
                if tail.endswith(suffix):
                    tail = tail[: -len(suffix)]
                    break
            rec.sample = tail.rstrip("_")
    return rec


def build_tree(root: Path, max_entries_per_dir: int = 8) -> str:
    """Render a compact directory tree, collapsing long runs of sibling files."""
    lines = [f"{root.name}/"]

    def walk(d: Path, prefix: str) -> None:
        entries = sorted(d.iterdir(), key=lambda x: (x.is_file(), x.name))
        dirs = [e for e in entries if e.is_dir()]
        files = [e for e in entries if e.is_file()]
        shown_files = files[:max_entries_per_dir]
        hidden = len(files) - len(shown_files)
        items = dirs + shown_files
        for i, e in enumerate(items):
            last = (i == len(items) - 1) and hidden == 0
            conn = "`-- " if last else "|-- "
            if e.is_dir():
                lines.append(f"{prefix}{conn}{e.name}/")
                walk(e, prefix + ("    " if last else "|   "))
            else:
                lines.append(f"{prefix}{conn}{e.name}  [{human(e.stat().st_size)}]")
        if hidden:
            lines.append(f"{prefix}`-- ... and {hidden} more file(s) of the "
                         f"same kind")

    walk(root, "")
    return "\n".join(lines)


def main() -> int:
    if not RAW.is_dir():
        print(f"ERROR: raw data directory not found: {RAW}", file=sys.stderr)
        return 1

    all_files = sorted(p for p in RAW.rglob("*") if p.is_file())
    print(f"Scanning {len(all_files)} files under {RAW} ...")

    records: list[FileRecord] = []
    for p in all_files:
        rec = classify(p)
        records.append(rec)
        print(f"  [{rec.file_type:45s}] {p.relative_to(RAW)}")

    df = pd.DataFrame([
        {k: v for k, v in asdict(r).items() if k != "extra"} | {
            "extra_json": json.dumps(r.extra, default=str)}
        for r in records
    ])
    csv_path = ANALYSIS / "raw_data_inventory.csv"
    df.to_csv(csv_path, index=False)

    # ---- text report -----------------------------------------------------
    out = []
    A = out.append
    A("=" * 78)
    A("RAW DATA INVENTORY")
    A("=" * 78)
    A("Root scanned : raw_data/")
    A(f"Files found  : {len(records)}")
    A(f"Total size   : {human(sum(r.size_bytes for r in records))}")
    A("")
    A("-" * 78)
    A("DIRECTORY STRUCTURE")
    A("-" * 78)
    A(build_tree(RAW))
    A("")

    for ds in sorted({r.dataset for r in records if r.dataset}):
        sub = [r for r in records if r.dataset == ds]
        A("-" * 78)
        A(f"DATASET: {ds}")
        A("-" * 78)
        A(f"  files: {len(sub)}   size: {human(sum(r.size_bytes for r in sub))}")
        kinds = pd.Series([r.file_type for r in sub]).value_counts()
        for k, v in kinds.items():
            A(f"    {v:3d} x {k}")
        spec = {r.species_evidence for r in sub if r.species_evidence}
        if spec:
            A(f"  species evidence: {'; '.join(sorted(spec))}")
        samples = sorted({r.sample for r in sub if r.sample})
        if samples:
            A(f"  samples inferred from filenames ({len(samples)}): "
              f"{', '.join(samples)}")
        A("")

    A("-" * 78)
    A("PER-FILE DETAIL")
    A("-" * 78)
    for r in records:
        A(f"path              : {r.path}")
        A(f"  file type       : {r.file_type}")
        A(f"  inferred role   : {r.inferred_role}")
        A(f"  dataset/sample  : {r.dataset} / {r.sample or '-'}")
        A(f"  size            : {r.size_human}")
        A(f"  dimensions      : {r.dimensions or '-'}")
        A(f"  raw counts?     : {r.contains_raw_counts}")
        A(f"  metadata?       : {r.contains_metadata}")
        A(f"  processed?      : {r.appears_processed}")
        if r.species_evidence:
            A(f"  species         : {r.species_evidence}")
        if r.notes:
            A(f"  notes           : {r.notes}")
        if r.extra:
            for k, v in r.extra.items():
                s = json.dumps(v, default=str)
                if len(s) > 400:
                    s = s[:400] + " ...(truncated)"
                A(f"    - {k}: {s}")
        A("")

    txt_path = ANALYSIS / "raw_data_inventory.txt"
    txt_path.write_text("\n".join(out), encoding="utf-8")

    # ---- machine readable ------------------------------------------------
    payload = {
        "raw_root": "raw_data/",
        "n_files": len(records),
        "files": [asdict(r) for r in records],
    }
    json_path = ANALYSIS / "raw_data_inventory.json"
    json_path.write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8")

    print(f"\nWrote {txt_path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
