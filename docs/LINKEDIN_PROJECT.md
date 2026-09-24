# LinkedIn project update

Draft prepared 25 September 2026 from the repository's current evidence.
This text has not been posted to LinkedIn.

## Title

scRNA_seq: Public Lung Single-Cell and Multiome Reanalysis of Injury, Repair and Pathological Remodelling

## Description — ready to paste

An ongoing computational research project reanalysing public mouse and human
lung single-cell RNA-seq and multiome (RNA + ATAC) data to develop testable
questions about injury, repair and persistent pathological remodelling.

Organising question: which epithelial, macrophage and fibroblast programmes
accompany productive repair, fibrosis and neoplasia-associated plasticity?

Current outputs:
• 14 research questions with supporting observations, alternative explanations
  and explicit requirements for further evidence.
• 168 evidence-register entries, alongside reproducible scripts, result tables,
  provenance records and figure galleries.
• Donor/animal-level analyses spanning injury time courses, pulmonary fibrosis,
  lung lesions, organoids and RNA–chromatin comparisons. The latest human
  lesion analysis includes 75 libraries from 23 patients and 555,480 QC-retained
  cells/nuclei, with paired molecular contrasts and spatial context.

Methods include Scanpy QC and annotation, dataset-specific batch assessment,
HLCA reference mapping with scvi-tools, biological-replicate pseudobulk,
limma/voom and CAMERA pathway analysis, ligand–receptor resource comparisons
with LIANA, conditional ligand–target analysis and programme-specificity checks.
PCA, UMAP, paired distributions and enrichment figures retain sample coverage
and model sensitivity.

A central result is that elevated plasticity-related RNA programmes are not
necessarily specific to neoplasia. Inferred communication and programme scores
generate hypotheses; they do not establish cell fate or causal mechanisms.
Null results, uncertain annotations and corrections remain visible.

Stack: Python, R/limma, Scanpy, scvi-tools, LIANA, pandas, matplotlib, Git and
GitHub Actions. AI-assisted development and scientific review are documented.

Repository: https://github.com/xorca0711/scRNA_seq

## Selected media order and captions

Reassessed 25 September 2026 against the displayed panels, sample units,
underlying tables and recorded figure hashes. Use these four figures plus the
GitHub repository link. Each contributes a different completed analysis;
the selection is based on relevance, readable evidence and coverage of the
project, not on statistical significance. Export the original full-resolution
PNGs without cropping away labels or limitations.

| Order and media title | Caption to paste | Asset |
|---|---|---|
| 1. Epithelial programmes across human lung-lesion contexts | GSE308103: PCA of 70 patient–histology AT2-like pseudobulks, programme scores in 23 paired normal/LUAD patients, and additional paired contrasts. These descriptive RNA patterns motivate specificity questions; histology and programme scores do not identify malignant cells or establish fate. | [PNG](../analysis/figures/rq/rq_a11_epithelial_programs.png) |
| 2. Different Wnt profiles in paired alveolar fibroblast populations | GSE262927: five source-defined Wnt ligands in AF1 and AF2, paired within 21 animals with at least 50 cells in each compartment. Wnt2 is higher in AF1 and Wnt4 in AF2 in every pair. This describes compartment-specific RNA, without establishing secretion, cell contact or an injury effect. | [PNG](../Thesis/gate1_03_nabhan_2018/nb1/figures/05_fibroblast_pairs.png) |
| 3. IL-1 source and recipient RNA in the human lung niche | GSE308103: donor-weighted ligand/receptor detection and paired LUAD–normal changes in reference-compatible macrophage, fibroblast and AT2-like populations. Pair counts vary by compartment. RNA compatibility does not establish functional signalling or the dominant tissue cytokine source. | [PNG](../analysis/figures/rq/rq_a12_source_recipient_components.png) |
| 4. How sequencing depth and reference resources affect AREG inference | Ten tumour donors show how the epithelial–myeloid AREG detection contrast changes after matching RNA depth. A separate 22-donor IPF/control analysis shows resource-dependent AREG–EGFR ranks; rank universes differ. These are measurements of detection and ranking, not signalling strength. | [PNG](../analysis/figures/rq/rq_a2_source_rank.png) |

### Why the earlier selection changed

- The CAMERA correlation-sensitivity panel is useful supporting evidence, but
  a weak introductory image: the primary tests do not meet the declared FDR
  threshold, and its filled symbols are sensitivity findings. It remains in
  the repository with the primary result alongside it.
- The full epithelial-specificity audit is dense and partly a coverage result
  (only one external animal evaluable). It does not add an independently
  validated epithelial result to the lead figure.
- The AREG three-panel summary replaces the more detailed C37 correction plot.
  It retains donor pairing and depth sensitivity while also showing resource
  dependence. Its caption explicitly distinguishes the two cohorts.
- The fibroblast Wnt panel adds a mouse study and a direct, paired biological
  comparison. It is not described as a temporal injury effect or a proven niche circuit.
- The older phase UMAP title implies a year “after repair”; age/harvest
  confounding and absent recovery outcomes limit that wording. The large
  multiome composite needs more explanatory space than this selection provides.
  Both remain available with their full repository context.
- Omit the claims-count graphic and the proposed withdrawal schematic from
  project media: neither is a measured biological result. This does not hide
  null findings or replace them with positive-only selections.

The accompanying local export contains numbered PNGs, copy-ready captions,
the project text and a source/hash manifest. Repository source figures remain
in their canonical locations; export copies are byte-identical.

## Verification and wording changes

- The [canonical register](../RESEARCH_QUESTIONS.md) now has A1–A14; the
  [generated claims index](../analysis/claims/manifest.json) contains 168 entries.
  Entries include methods, limitations and negative results, not 168 discoveries.
- Human counts and inferential limits are documented in the
  [completed review](../Thesis/gate2_C3_yu_lee_choi_min_2026/EVIDENCE_REVIEW.md)
  and [annotation review](../Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u5_human_full/ANNOTATION_REVIEW.md).
- Removed the stale accession/study, total-cell, script, run-record and figure
  totals. Repeated analyses, superseries, companion assays and regenerated
  figures make an unqualified portfolio-wide total misleading.
- Removed the blanket “pre-registered” label: dated specifications and
  post-analysis exploratory questions do not establish unseen-data preregistration.
- Added fibroblast/recipient context, R/limma, correlation sensitivity and the
  newer human/spatial analysis. The project does not claim that every method
  is newer than every source publication or that every programme replicated.
