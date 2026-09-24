# Figures for the derived research questions

25 September 2026. Five measured-data figures and one explicitly labeled
experimental-design schematic are attached to the
[RQ proposal](../../DERIVED_RESEARCH_QUESTIONS.md) and
[paper gallery](../../README.md#d0-source-identity-and-annotation-context).
They supplement the original 17 figures. No original annotations, thresholds,
DE fits, pathway results or ligand-target fits were changed.

| Figure | RQ and purpose | Main plotted data |
|---|---|---|
| [D0: source context](../../figures/derived_D0_source_annotation_context.png) | Enabling source-identity question: reference UMAP, uncertainty, IL1B expression and paired source-allocation violins | [Coordinates and display values](umap_display_values.csv); [all/display coverage](embedding_sampling_coverage.csv); [paired source fractions](D0_paired_unassigned_IL1B.csv) |
| [D1: epithelial programs](../../figures/derived_D1_epithelial_programs.png) | RQ1: patient PCA, paired HPCS violins, all seven HPCS contrasts and a six-program heatmap | [PCA points](AT2_pseudobulk_PCA_values.csv); [PCA features](AT2_PCA_features.csv); [paired scores](D1_HPCS_paired_scores.csv); [all contrasts](D1_HPCS_all_contrasts.csv); [heatmap](D1_program_heatmap_values.csv) |
| [D2a: components](../../figures/derived_D2_source_recipient_components.png) | RQ2 and context for RQ3: donor-weighted component dot plot and paired component/edge changes | [Dot summaries](D2_component_dot_values.csv); [donors](D2_component_donor_values.csv); [paired components](D2_paired_components.csv); [paired edges](D2_paired_edge_values.csv) |
| [D2b: enrichment](../../figures/derived_D2_recipient_enrichment.png) | RQ2: all 21 broad-receiver LUAD-normal tests under each correlation setting | [All displayed rows](D2_pathway_values.csv) |
| [D2c: ligand targets](../../figures/derived_D2_ligand_target_candidates.png) | RQ2: all 32 eligible focused-triad up-target candidates across AT2-like, fibroblast and macrophage receivers | [Fits, ranks and omission coverage](D2_ligand_target_values.csv) |
| [D4: withdrawal experiment](../../figures/derived_D4_withdrawal_experiment_proposal.png) | RQ4: proposed experimental design; no observed or simulated outcomes | [Experimental question](../../DERIVED_RESEARCH_QUESTIONS.md#rq4-resolution-versus-persistence-after-signal-withdrawal) |

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

RQ3 joint associations and regional spatial inference are not newly fit.
Current spatial maps remain linked in the proposal; missing independent
regions and complete-compartment coverage still limit inference. D4 is a
design schematic only. The D0 diagnostics do not resolve unassigned cell
identity; additional source review remains a separate scientific task.

Run the scripts from the repository root with the existing scientific runtime:

```powershell
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u7_prepare_proposal_figures.py
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u7_render_proposal_figures.py
python analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u7_validate_proposal_figures.py
```

Preparation refuses to overwrite an existing run record. Its saved coordinates
can be reused by rendering without repeating UMAP. Rendering may regenerate
presentation artifacts; changed images need renewed visual review before
validation and release.

[Specification](specification.json) · [preparation record](preparation_run_record.json) ·
[input hashes](preparation_inputs.csv) · [render record](render_run_record.json) ·
[render inputs](render_inputs.csv) · [donor QC summaries](donor_QC_summary.csv) ·
[visual review](visual_review.json) · [validation](validation.json).
