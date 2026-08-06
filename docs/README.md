# Documentation Index

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

Supporting documents:

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
