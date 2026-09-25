# Figures for the derived research questions

25 September 2026. Five measured-data figures and one explicitly labeled
experimental-design schematic are attached to the
[RQ proposal](../../../../RESEARCH_QUESTIONS.md#a11-shared-plasticity-versus-neoplasia-associated-context) and
[paper gallery](../../../../Research%20Article/gate2_C3_yu_lee_choi_min_2026/README.md#figure-gallery).
They supplement the original 17 figures. No original annotations, thresholds,
DE fits, pathway results or ligand-target fits were changed.

| Figure | RQ and purpose | Main plotted data |
|---|---|---|
| [A12-S1: source context](../rq_a12_source_context.png) | Enabling source-identity question: reference UMAP, uncertainty, IL1B expression and paired source-allocation violins | [Coordinates and display values](umap_display_values.csv); [all/display coverage](embedding_sampling_coverage.csv); [paired source fractions](D0_paired_unassigned_IL1B.csv) |
| [A11: epithelial programs](../rq_a11_epithelial_programs.png) | A11: patient PCA, paired HPCS violins, all seven HPCS contrasts and a six-program heatmap | [PCA points](AT2_pseudobulk_PCA_values.csv); [PCA features](AT2_PCA_features.csv); [paired scores](D1_HPCS_paired_scores.csv); [all contrasts](D1_HPCS_all_contrasts.csv); [heatmap](D1_program_heatmap_values.csv) |
| [A12a: components](../rq_a12_source_recipient_components.png) | A12 and context for A13: donor-weighted component dot plot and paired component/edge changes | [Dot summaries](D2_component_dot_values.csv); [donors](D2_component_donor_values.csv); [paired components](D2_paired_components.csv); [paired edges](D2_paired_edge_values.csv) |
| [A12b: enrichment](../rq_a12_recipient_enrichment.png) | A12: all 21 broad-receiver LUAD-normal tests under each correlation setting | [All displayed rows](D2_pathway_values.csv) |
| [A12c: ligand targets](../rq_a12_ligand_target_candidates.png) | A12: all 32 eligible focused-triad up-target candidates across AT2-like, fibroblast and macrophage receivers | [Fits, ranks and omission coverage](D2_ligand_target_values.csv) |
| [A14: withdrawal experiment](../rq_a14_withdrawal_design.png) | A14: proposed experimental design; no observed or simulated outcomes | [Experimental question](../../../../RESEARCH_QUESTIONS.md#a14-resolution-versus-persistence-after-signal-withdrawal) |

Every figure has a sibling SVG. UMAP point clouds are rasterized within the SVG
to keep files practical; labels and other plot elements remain vector objects.

## New diagnostic computations

**UMAP.** The frozen human reference-mapping run already saved 30-dimensional
latent vectors and corresponding cell annotations for all 75 libraries.
The figure preparation checks the latent/annotation row order and uses the
existing fine-label confidence cutoff of 0.2. It samples at most 75 cells per
patient, histology and existing source-profile broad label, after pooling the
available library rows into those strata. This yields 34,178 display cells
from 555,480 QC nuclei and retains all 23 patients. Full population counts by
source label agree exactly with the completed source-profile tables.

UMAP uses Euclidean distance, 30 neighbors, min_dist 0.35, 300 epochs and seed
20260924. No reference model or cell labels were refit. All three overlays use
the same coordinates. Display balancing means density cannot estimate captured
cell abundance. The reference-conditioned geometry is not independent
annotation validation or a lineage trajectory. IL1B expression is raw counts
normalized by each cell's full-assay total to 10,000, then log1p-transformed.
Its color scale retains the full observed display range. Source fractions in
the paired violin use all QC cells, not the display subsample.

**PCA.** Seventy broad AT2-like patient-histology pseudobulks from 23 patients
pass the original confidence and 50-cell floors. Raw counts must sum exactly
to recorded full-assay library totals. Diagnostic log2(CPM+1) expression is
filtered to genes with CPM >=1 in at least five aggregates; 2,000 genes are
selected by variance without using histology labels, then centered without
variance scaling. PC1 and PC2 explain 20.25% and 10.45% of variance. No batch
correction or histology-separation optimization is applied. This normalization
is a display diagnostic; the released program and pathway results retain
their original TMM/voom models. Patient pairing and public patient labels remain
available in the PCA table. No PCA separation test was performed.

## Existing quantitative results, newly displayed

The reduced-HPCS violin and program heatmap use existing primary TMM scores;
their paired differences are checked against the released patient table.
The violin has 23 independent patients with paired observations, and the
sparser histology contrasts retain individual points without a density curve.
Black bars are medians; diamonds in difference panels are means. There are no
cell-level significance stars.

Component dot plots weight each eligible donor equally within histology;
they display every specified ligand, receptor and regulatory component.
Normal and LUAD donor coverage differs by compartment, and the dot plot is a
descriptive group summary. The adjacent paired plots use only complete
patients for each specific component or edge. Recovered RNA and compatibility
scores do not assay mature cytokine release or communication.

The enrichment figure retains every eligible declared broad-compartment test
for LUAD-normal, with global q values from the full 279-test family. Filled
triangles meet q<0.05 under the indicated model, and direction is indicated by
both shape and color. Missing marks identify sets not assigned to that receiver.
Primary estimated correlation and fixed-0.01 sensitivity remain distinct.
CAMERA direction/q values are not normalized enrichment scores or effect sizes.

The ligand-target figure includes all eligible candidates for the selected
contrast, direction and source scope, not only IL1B or successful top hits.
Horizontal ranges come from eligible omitted-patient fits and are not
confidence intervals. Candidate-specific omission counts and rank ranges are
in the plotted table. Different receiver backgrounds prevent comparing these
correlations as causal effect sizes. Other contrasts and down-target fits
remain in the original analysis, not silently discarded.

## Remaining scope and reproduction

A13 joint associations and regional spatial inference are not newly fit.
Current spatial maps remain linked in the proposal; missing independent
regions and complete-compartment coverage still limit inference. A14 is a
design schematic only. The source diagnostics do not resolve unassigned cell
identity; additional source review remains a separate scientific task.

Run the scripts from the repository root with the existing scientific runtime:

```powershell
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages analysis/scripts/19_prepare_il1b_rq_figures.py
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages analysis/scripts/20_render_il1b_rq_figures.py
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages analysis/scripts/21_validate_il1b_rq_figures.py
```

Preparation refuses to overwrite an existing run record. Its saved coordinates
can be reused by rendering without repeating UMAP. Rendering may regenerate
presentation artifacts; changed images need renewed visual review before
validation and release.

[Specification](specification.json) · [preparation record](preparation_run_record.json) ·
[input hashes](preparation_inputs.csv) · [render record](render_run_record.json) ·
[render inputs](render_inputs.csv) · [donor QC summaries](donor_QC_summary.csv) ·
[visual review](visual_review.json) · [validation](validation.json).

## Repository layout migration, 25 September 2026

Shared RQ assets now live here; study-specific inference stays under `Research Article/`.
CSV basenames D0/D1/D2 retain their historical source-table IDs. Presentation
labels are A11, A12a–c, A12-S1 and A14, matching the root question register.
The UMAP/PCA preparation was not rerun. Its original entrypoint bytes and run
identity are preserved in `.history/layout_2026-09-25/`;
[the relocation manifest](relocation_manifest.json) records every original
path and hash. Rerendering changes labels and destinations, not inferential results.
Preparation inputs use paths relative to the Yu study; current render inputs
and output paths are relative to the repository root. The validator handles
these explicit namespaces. A clean checkout can replot from tracked tables;
full input-hash validation additionally requires the original ignored caches.
