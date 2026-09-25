# Nb1: Wnt source, response and AT2 state after viral injury

Specification frozen on 2026-09-22 before examining Wnt expression in this run.
The owner has completed Nabhan 2018 and explicitly authorized this sequence.

## Questions and scope

1. Which deposited fibroblast states express the Wnt ligands implicated by
   Nabhan 2018, and is the pattern shared across animals?
2. Are epithelial Wnt transcripts and secretion-machinery transcripts detectable
   after viral injury, alongside the small Wnt-response panel?
3. How do Wnt-response transcripts relate descriptively to AT2 identity and
   proliferation, without selecting cells on Axin2 expression?

This is a transcriptional assessment. Neither expression nor coexpression
establishes secreted ligand, spatial contact, autocrine signaling, stem-cell
function or a productive repair outcome. Nabhan's perturbation experiments
provide mechanistic motivation, not validation of these observational results.

## Input and units

- GSE262927 `processed/final_clustered.h5ad`, raw `layers/counts` only.
- Match the deposited cell metadata by cell ID; preserve animal IDs, day, sex,
  genotype and infection round. Check matrix/metadata alignment explicitly.
- Use deposited author labels AT2, AT1_AT2, AT1, AF1, AF2,
  Adventitial_fibroblast and Peribronchial_fibroblast. AF1/AF2 are the primary
  alveolar fibroblast compartments; the other fibroblasts are context, not
  interchangeable alveolar niche cells. Exclude VSMC and mesothelium.
- Include cells with an animal ID, author annotation and finite sacrifice day.
  Do not define an Axin2-positive subset or select clusters after seeing Wnts.
- The analysis unit is animal × deposited cell type. Primary eligibility is
  >=50 cells per unit. Show 20- and 100-cell sensitivity as coverage checks;
  they cannot create animal replication. Do not pool days or genotypes to
  manufacture a contrast. All-genotype results are descriptive; show the
  heterozygous-only eligibility sensitivity.
- Metadata already show first injury sampling at day 6, two baseline animals
  (one AT2 unit >=50 cells), and two day-11 animals (one AT2 unit >=50 cells).
  Thus neither the <=24-hour source switch nor a replicated baseline-versus-
  day-11 test is supported locally. No local differential-expression p-values.

## Fixed panels and endpoints

- Survey all 19 mouse Wnt ligands: Wnt1, Wnt2, Wnt2b, Wnt3, Wnt3a, Wnt4,
  Wnt5a, Wnt5b, Wnt6, Wnt7a, Wnt7b, Wnt8a, Wnt8b, Wnt9a, Wnt9b, Wnt10a,
  Wnt10b, Wnt11 and Wnt16. Highlight the paper's fibroblast panel
  Wnt5a/Wnt2/Wnt2b/Wnt4/Wnt9a and injury-associated Wnt7b.
- Secretion machinery: Porcn and Wls; RNA alone does not show secretion.
- Wnt-response panel: Axin2 and Lef1, reported in Nabhan 2018. Treat them as
  two measured transcripts, not a comprehensive or validated activity score.
- Identity: the existing published AT2 holdout panel in
  `Research Article/epithelial_state_specificity/modules.json`, preserving its provenance.
- Proliferation: Mki67 and Top2a, displayed individually; their mean normalized
  expression is only a compact two-gene marker summary, not a cell-cycle model.
- Sum raw counts by unit, collapse duplicate gene symbols before calculating
  CPM, and report each gene's positive-cell fraction with its numerator and
  denominator. Report zeros and missing features distinctly.
- For depth sensitivity, calculate the expected positive-cell fraction after
  sampling exactly 1,000 UMIs without replacement from each eligible cell.
  Cells below that budget are excluded from this sensitivity only; report
  their counts. Use the exact hypergeometric probability, not random draws.
  This addresses library depth, not ambient RNA, batch or capture efficiency.
- Identity summary is the mean log2(CPM+1) across present holdout genes.
  Report panel coverage; no inference that this score measures stemness.
- Show per-animal source/response plots and compartment heatmaps. Do not
  interpret pooled-cell correlations, UMAP separation or a few detected
  transcripts as independent biological replication. No significance testing
  or correlation p-values for the local study.

## Decision gates and outputs

Source-data reproduction (GSE109444) is a separate FPKM analysis and must not be
numerically pooled with UMI counts. Its bulk control is excluded from cell
fractions; lack of established animal replication limits inference.

An external cohort needs identifiable independent animals, compatible injury
and control arms, relevant time points, and >=3 eligible animals per arm and
required compartment before a replicated test is designed. This minimum is a
feasibility threshold, not a power calculation. Freeze a separate contrast and
covariates before examining expression if a cohort passes; otherwise stop at a
documented HOLD rather than relabeling the descriptive results as confirmation.

Deliver input hashes, a reproducible script, eligibility and expression tables,
PNG/SVG figures, integrity checks and a report separating observations from
unsupported mechanisms. Raw matrices and large caches remain untracked.

Primary paper: https://doi.org/10.1126/science.aam6603.
