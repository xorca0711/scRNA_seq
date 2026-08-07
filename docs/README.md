# Documentation Index

> **Reference material, not a record of what was run.** The pages below
> describe tools and the reference study's design. They do **not** imply that
> every tool described here was used in the analysis in `analysis/` — most were
> not. For what was actually executed, with real parameters and cell counts,
> see **[`PIPELINE_AS_RUN.md`](PIPELINE_AS_RUN.md)**, which is generated from
> the pipeline's own outputs and states tool-by-tool which were used.
>
> Short version: of the five tools documented here, only **Scrublet** was used.
> SoupX, scds, Slingshot and tradeSeq were **not**. `PIPELINE_AS_RUN.md` also
> records where the executed analysis departs from the two source publications.

One schematic per tool. Each file states what the tool consumes and returns,
gives a `flowchart` of its internal decisions, lists its parameters and their
sources, and records the failure modes that are silent rather than loud.

- [`SOUPX.md`](SOUPX.md): ambient RNA estimation from empty droplets, automated
  and marker-based contamination fractions, and why correction must precede
  doublet calling in lung tissue.
- [`SCRUBLET.md`](SCRUBLET.md): simulated-doublet neighbourhood scoring,
  threshold selection on the simulated-score histogram, the two free
  consistency checks, and the limits of neotypic detection.
- [`SCDS.md`](SCDS.md): co-expression (`cxds`) and classifier (`bcds`) scoring,
  their `hybrid` combination, benchmark position against other callers, and the
  fact that scds deliberately returns no threshold.
- [`SLINGSHOT.md`](SLINGSHOT.md): minimum spanning tree topology, simultaneous
  principal curves, the two unreported knobs that decide the answer, and what
  supervision does to the interpretation of a trajectory.
- [`TRADESEQ.md`](TRADESEQ.md): the negative binomial GAM, knot selection by
  AIC, the mapping from biological question to statistical test, and the
  pseudotime-comparability caveat that applies directly to the reference study.

Rationale and background (written after the analysis ran):

- [`ANALYSIS_RATIONALE.md`](ANALYSIS_RATIONALE.md): every major decision in two
  passes — what was decided from the data alone, and what changed after reading
  the two source papers, including where the first pass was wrong.
- [`BACKGROUND_FOR_BIOLOGISTS.md`](BACKGROUND_FOR_BIOLOGISTS.md): batch effects
  and Harmony from first principles, for a reader without a computational
  background — including when correction destroys the experiment.
- [`DOUBLETS_AND_SCRUBLET.md`](DOUBLETS_AND_SCRUBLET.md): what a doublet is, how
  Scrublet works, and the measured bias by which it removed this project's own
  populations of interest at up to twice the background rate.
- [`UMAP_AND_FIGURES.md`](UMAP_AND_FIGURES.md): how the UMAP is built from
  counts, what it does and does not mean, how to read dot plots and feature
  plots, and why cluster-marker p-values are a ranking device, not a test.

Supporting documents:

- [`PIPELINE_AS_RUN.md`](PIPELINE_AS_RUN.md): **what was actually executed** —
  per-dataset parameters, cell counts, batch decisions, which of the tools
  above were used, and how the analysis diverges from the source papers.
  Generated from the pipeline's artefacts; regenerate with
  `analysis/scripts/05_write_pipeline_as_run.py`.
- [`../analysis/README.md`](../analysis/README.md): the full analysis report
  for the primary (mouse) dataset, including QC tables and output locations.
- [`../WORKFLOW.md`](../WORKFLOW.md): the end-to-end sequence, ordering
  constraints, subset-and-recluster loop, lineage-trace calling, and quality
  control acceptance order.
- [`../scRNAseq_workflow_Niethamer2025.md`](../scRNAseq_workflow_Niethamer2025.md):
  the annotated pipeline reference — study design, stage-by-stage parameters,
  marker-gene annotation tables, and the twelve parameters the reference study
  leaves unspecified.
- [`../REFERENCES.md`](../REFERENCES.md): all six papers with DOIs, PMC links,
  and software repositories.

Alignment (STARsolo) and the Seurat stages have no separate schematic. Their
steps are linear and fully covered in [`../WORKFLOW.md`](../WORKFLOW.md); a
per-tool diagram would only restate it.

No PDFs are stored in this repository. Three of the five method papers are
CC BY 4.0 and three are not, so all six are linked to their open-access PMC
versions in [`../REFERENCES.md`](../REFERENCES.md) instead.
