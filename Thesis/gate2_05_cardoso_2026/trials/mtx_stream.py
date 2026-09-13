"""Stream a MatrixMarket coordinate file and keep only selected gene rows.

The human fibrosis atlases deposit one merged sparse matrix each, too large to
load whole on this machine. Every trial that needs a handful of genes out of
such a matrix can make a single pass instead: read the triplets, keep the
entries whose gene row is wanted, and build dense vectors over cells. Memory
is then a few gene vectors, not the matrix.

The caller supplies the gene and barcode lists, because the deposits name
those files differently.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt", encoding="utf-8")


def read_lines(path: Path, column: int = 0, sep: str = "\t") -> list[str]:
    """Read one column of a plain or gzipped text file, ignoring a quoted header."""
    out = []
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split(sep)
            out.append(fields[column].strip('"'))
    return out


def extract_gene_rows(matrix: Path, gene_names: list[str], wanted: set[str]) -> dict:
    """Return {gene: dense float32 vector over cells} from a coordinate matrix.

    The MatrixMarket header is read to get the shape, and the file is assumed
    to be genes by cells, which is how these deposits store it. Row indices are
    1-based in the format.
    """
    index_of = {}
    for position, name in enumerate(gene_names):
        if name in wanted:
            index_of[position + 1] = name
    vectors: dict[str, np.ndarray] = {}
    with open_text(matrix) as handle:
        shape = None
        for line in handle:
            if line.startswith("%"):
                continue
            shape = [int(v) for v in line.split()]
            break
        if shape is None:
            raise ValueError(f"{matrix} has no dimension line")
        n_rows, n_cols = shape[0], shape[1]
        if n_rows != len(gene_names):
            raise ValueError(f"{matrix}: {n_rows} rows but {len(gene_names)} gene names")
        for name in index_of.values():
            vectors[name] = np.zeros(n_cols, dtype=np.float32)
        for line in handle:
            first, second, third = line.split()
            row = int(first)
            name = index_of.get(row)
            if name is None:
                continue
            vectors[name][int(second) - 1] = float(third)
    return vectors
