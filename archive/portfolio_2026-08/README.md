# Displaced portfolio material (August 2026)

This directory holds two files that were moved out of the analysis tree on
2026-09-10 with `git mv`, so their history is intact.

| File | What it is |
|---|---|
| [`lung_scrna_portfolio_thesis_context.pdf`](lung_scrna_portfolio_thesis_context.pdf) | A thesis-aware portfolio PDF generated on 2026-08-13 from the tracked analysis artefacts of that date (mouse GSE262927 regeneration and lineage-tracing results, human GSE178360 epithelial results, limitations). |
| [`09_write_portfolio_pdf.py`](09_write_portfolio_pdf.py) | The reportlab script that generated the PDF. It previously lived at `analysis/scripts/09_write_portfolio_pdf.py`. |

## Why they were displaced

The PDF and its generator are portfolio curation material. That material is
established elsewhere and is not part of the ongoing analysis record; the
repository now tracks analysis progress only. Keeping the files under
`output/` and `analysis/scripts/` made them look like current pipeline
outputs, which they are not.

## Status of the script

The script is kept for the record only. It is not maintained, and it is not
listed in any re-run list ([`../../REPRODUCIBILITY.md`](../../REPRODUCIBILITY.md),
the CI workflow, or the per-dataset reports). If it is run from this
directory it still resolves the repository root from its own location,
recreates `output/pdf/` and writes
`output/pdf/lung_scrna_portfolio_thesis_context.pdf` there; that directory is
not part of the tree any more and `*.pdf` is gitignored, so a fresh run leaves
an untracked file behind. The text the script embeds still refers to the
validator and to itself by their pre-2026-09-10 names, and its headline
numbers describe the artefacts as they stood on 2026-08-13.

## Figures

Every figure the PDF embedded was read from [`../../analysis/`](../../analysis/)
and remains there unchanged: the GSE178360 integration UMAP and the
reference-aligned epithelial UMAP and feature plots, the GSE262927
regeneration-focus pseudotime UMAP and time courses, and the lineage-tracing
cohort figure. None of them was moved or edited as part of this displacement.
