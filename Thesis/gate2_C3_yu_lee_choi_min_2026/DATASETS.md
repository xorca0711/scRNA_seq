# Public-data shortlist and acquisition order

Verified on 24 September 2026. Ten new series have sample-level GEO MINiML
metadata saved under [U0](trials/u0_geo_design_audit/run_record.json), with
source URLs, retrieval times and response hashes. **GSM counts are records,
not automatically biological replicates.** Priority means suitability for the
specified question, not a declaration that all data are ready to analyse.

## Current execution

The [completion register](WORK_PACKAGES.md) records the completed public-data
analyses. All 75 human RNA libraries, 56 human spatial sections and nine
post-viral matrices were processed. Both IPF cohorts, eligible early mouse
niches and specificity extensions have measured outputs. The shortlist below
preserves source-design questions; endpoint-level resolutions and remaining
limits are in the [evidence review](EVIDENCE_REVIEW.md).

## New verified deposits

| Priority / series | Design actually visible in the deposit | Role and decision |
|---|---|---|
| 1: [GSE300288](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300288) | 31 mouse scRNA-labelled libraries; Gprc5a-null, NNK exposure; IgG, anti-IL-1beta, anti-PD-1 and combination; endpoints 3 and 7 months | Best direct perturbation candidate. Early anti-IL-1beta has 3 libraries; other groups have 4. Resolve animal/pool IDs, missing early replicate 1, treatment history and FFPE-versus-single-cell metadata before inferential release |
| 2: [GSE308103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308103) | 75 human fixed-tissue snRNA sample records; 23 title-derived patient labels, 22 with candidate normal/precursor/LUAD sets | Human within-patient context. Repeated tissue pieces do not increase n. Raw reads are not publicly deposited; processed matrices are listed. Preserve AAH/AIS/MIA separately |
| 3: [GSE307534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307534) | 56 human Visium CytAssist records; 25 patient labels, 24 candidate precursor/LUAD pairs | Spatial companion to the same study; validate lesion annotations, coordinates and cross-assay specimen links. Not independent replication of GSE308103 |
| 4: [GSE267226](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267226) | Human Visium, 3 PASC-PF patients and 2 controls; H5 matrices and tissue images listed | Different pathological repair context, descriptive because control n=2. Tissue-position files were not listed in sampled metadata; find them before spatial neighbourhood analysis |
| 4: [GSE267228](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE267228) | Four aged mice at 60 days post-injury; IgG (2) versus anti-CD8 (2), Visium | Immune perturbation context, descriptive. **No deposited anti-IL-1beta arm.** Do not substitute the paper's functional IL-1beta experiment for this actual design |
| 5: [GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259) | 60 sample records across whole-lung/high-resolution experiments and bleomycin time points | Injury time reference; resolve animal/sort/experiment links before combining records. Existing repository ADI definitions derive from this study, so it is not an independent signature-validation set |
| 6: [GSE277777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277777) | 32 records, including GEX and hashtag libraries, lineage reporters, treatments and a transplant arm | HPCS malignant specificity extension. Deconvolve hash/animal IDs, technical sequencing and model strata; 32 is not mouse n |
| Companion: [GSE222901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222901) | 35 mouse records; Gprc5a-null NNK/saline time course with separate lineage-sorted and spatial experiments | Untreated/stage context for the Peng/Han work. Overall design reports four mice per main group/time point; not every record belongs to that experiment. Check reuse/overlap across publications |
| Companion: [GSE300293](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300293) | 19 mouse spatial records | Spatial context for the same study. Some generic tumour-bearing sample titles coexist with saline characteristics; reconcile before classifying histology |
| Crosswalk only: [GSE307529](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307529) | 54 human WES records, 25 title-derived patient labels | Potential genotype/lesion crosswalk. This is DNA, not an expression cohort. No initial variant download or independent-cohort claim; title/file histology disagreements need source-table resolution |

The six Peng/Han-related deposits are linked by the primary study's
[data availability statement](https://doi.org/10.1016/j.ccell.2025.10.004).
Its public code and processed data are also linked through
[LungPCA_Code](https://github.com/FuduanPeng/LungPCA_Code) and
[Zenodo](https://doi.org/10.5281/zenodo.17172149). Metadata visibility does
not mean that every large spatial download or useful annotation has been
retrieved. Request processed subsets before downloading entire archives.

The injury series are grounded in [Narasimhan et al.](https://doi.org/10.1038/s41586-024-07926-8)
and [Strunz et al.](https://doi.org/10.1038/s41467-020-17358-3); the malignant
comparator in [Chan et al.](https://doi.org/10.1038/s41586-025-09985-x).

## Existing deposits: reuse for a declared purpose

| Existing data | Role in this plan | Limit carried forward |
|---|---|---|
| GSE145031 / GSE144468, Choi 2020 | DATP definition, exposure-direction positive control | One library per condition/arm; no replicated treatment test or deposited withdrawal series established |
| GSE262927, Niethamer 2025 | Non-oncogenic repair reference, within-compartment programmes | Time/age and other design confounds; sparse late epithelial states; already inspected |
| GSE136831 / GSE135893, Adams / Habermann | Required candidate fibrotic niche arm: macrophage/fibroblast LR and pathway analyses, alongside transitional definitions | Extend full-source extracts to include macrophages; old epithelial-stromal caches are incomplete for the triad. Previously inspected, cross-sectional cohorts; donor/subtype coverage and dissociation bias remain |
| GSE247505 and its relevant subseries, England 2025 | KRAS/Il1r1 context | Two replicate libraries per relevant arm; genotype metadata is essential; current evidence already inspected |
| GSE316241 / GSE316243 / GSE316244, Cardoso 2026 | Stromal and epithelial niche hypotheses and known resource sensitivity | One pooled library per genotype/sort; no new treatment inference by pooling cells |
| GSE310539 / GSE247130 | Injury/development/genotype specificity and chromatin-context controls | One well per condition, known suffix-map correction; no reversibility or animal-level chromatin confirmation |

Reuse [existing result tables and limitations](../../RESEARCH_QUESTIONS.md)
first. A fresh question may justify a bounded targeted extraction, but these
are not a queue for rerunning the entire repository.

## Context and compartment gates for the revised niche questions

The [niche plan](NICHE_ANALYSIS_PLAN.md) makes macrophage/fibroblast
ligand-receptor analysis and pathway enrichment required proposed arms.
For GSE300288 and GSE308103, actual eligible cell coverage is now measured
in the [completed context figure](README.md#figure-gallery) and biological-unit
tables. Comparisons require same-sample
sender/receiver coverage for each edge, all three compartments for a joint
triad comparison, and verified animals or patient/lesion IDs. An assay that
captures only epithelial cells cannot answer a missing stromal question.

Keep four roles distinct: direct blockade (GSE300288), fibrotic niche
associations (GSE136831/GSE135893), repair references (GSE141259/GSE262927),
and human precursor/cancer context (GSE308103/GSE307534). CRC/PDAC/breast
examples in the notes generate labelled lung hypotheses, not extra lung
replication cohorts. Age, mechanics and exposure duration are testable only
where independent metadata or measurements support them. No new cross-organ
download programme is included in this revision.

## Leads with incomplete readiness

- **Han 2023 mitochondrial ISR:** the primary paper lists PRJNA865889,
  PRJNA940730, PRJNA940746, PRJNA940973, PRJNA940986 and PRJNA940992 for
  sequencing reads, plus [author code](https://github.com/MinhoLee-DGU/2023.Han.et.al.Nature).
  [Source](https://doi.org/10.1038/s41586-023-06423-8). A reusable processed
  matrix and its biological-unit crosswalk remain unverified. Keep this as a
  conditional developmental specificity comparator; no FASTQ launch.
- **Ciminieri 2023 IL-1beta/fibroblast organoids:** relevant functional
  experiment, [primary paper](https://doi.org/10.1165/rcmb.2022-0209OC).
  No accession was verified from the paper inspected here. Keep the entry
  unresolved rather than assigning an unrelated GEO result.
- **Choi 2021 airway plasticity:** review ref. 44 is relevant to lineage
  alternatives. Public accession and reusable sample design remain to be
  verified; it is not part of the launch-ready count queue.
- **GSE307112/GSE307128:** the repository's
  [separate outcome-linked screen gate](../../docs/NEXT_DATASET_GATE.md) already
  audits these. Organoid growth/morphology provides a valuable orthogonal
  endpoint if well/preparation identities are resolved. It is not a direct
  IL-1beta withdrawal or mature-AT1 fate assay.

## Acquisition and independence rules

Use the [sample inventory](trials/u0_geo_design_audit/sample_inventory.csv),
[GSE300288 design](trials/u0_geo_design_audit/GSE300288_design.csv) and
[candidate human links](trials/u0_geo_design_audit/human_sample_links.csv)
before requesting count files. Keep one-to-many capture/lesion relationships.
Assay companions and reanalysed public cohorts are not extra independent
replications. Use exact metadata-listed file URLs, local checksums and sparse
streaming. Only two GSE300288 count examples have been downloaded so far.
The owner authorized staged acquisition and analysis after final review.
The first seven-library mouse batch is acquired; later-stage gates remain.
The linked Zenodo files are restricted as of launch, and the older DOI in
the author README returns 404; annotated objects are not assumed public.
See [final source findings](FINAL_REVIEW.md).
