# A1 analysis plan: regulatory state, lineage and phenotype

25 September 2026. Prospective reanalysis plan following a targeted literature
and public catalog search; not an unseen-data preregistration. Published
results and earlier repository analyses informed the hypotheses. Execution
has not started beyond metadata retrieval. Sources and assay availability
are in [STUDY_MAP.md](STUDY_MAP.md) and [metadata](metadata/README.md).

## Question and possible answers

Which chromatin programmes distinguish transitional alveolar epithelial
states beyond overlapping RNA markers, and which distinctions agree with
lineage history, protein phenotype, tissue position or perturbation response?

Test three alternatives without deciding the answer through annotation:

1. A common regulatory programme spans contexts; the published names mainly
   capture timing, sampling or experimental differences.
2. A common transition programme contains distinct regulatory branches with
   different protein phenotypes, destinations or responses.
3. Apparently similar RNA states arise through different regulatory routes;
   their similarity reflects generic stress, proliferation or identity loss.

An unresolved or continuous distinction is a valid answer. Keep differentiation,
senescence/stress, proliferation, basal conversion and malignancy on separate
axes. Do not equate RNA entropy with plasticity or loss of AT2 identity with arrest.

| Subquestion within A1 | Primary measurement | Evidence that would strengthen it |
|---|---|---|
| Regulatory distinction | Accessible regions, histone-mark occupancy and, where measured, CpG methylation | Repeated within-study contrasts and a separate assay/cohort |
| Dynamic distinction | Sampling time, transcriptomic branch and recorded lineage label | Traced descendants, mature AT1 contribution or persistent alternative fate |
| Phenotypic distinction | Protein abundance/localization, morphology and tissue neighbourhood | Independent donor replication and functional response |
| Context specificity | Repair, development, fibrosis and cancer analysed separately | Transport of a frozen feature set without using context to define the answer |

A5 owns developmental specificity, A11 shared RNA programmes and A14
withdrawal/future-fate experiments. A1 links to these questions instead of
creating another global question register.

## What each assay can establish

- ATAC/FAIRE measures accessibility, not histone modifications or DNA methylation.
- ChIP/CUT&Tag directly measures a specified mark or occupancy in a population.
  Co-occurring marks in bulk do not establish single-cell bivalency.
- WGBS measures methylation at covered cytosines. Standard bisulfite data do
  not distinguish 5mC from 5hmC. One donor cannot establish population DMRs.
- CTCF binding, motif similarity and co-accessibility do not measure 3D contacts.
  No state-matched Hi-C/Micro-C dataset passed this search's availability audit.
- Genetic lineage tracing measures labelled ancestry/descendants under a defined
  pulse/chase. It is stronger fate evidence than a computational trajectory;
  reporter specificity, label persistence, proliferation and death still matter.
- Protein-sorted RNA is evidence about a protein-enriched population, not
  unbiased proteomics. Imaging mass cytometry measures its antibody panel;
  laser-capture proteomics measures a regional mixture rather than a single cell.
- Perturbation can test dependency within the model. Drug effects require
  viability/composition checks; disappearance of a state is not demonstrated repair.

## Dataset roles and first execution order

**First, finish the sample crosswalk; then prioritise small processed tables.**
The initial shortlist is deliberately broader than the executable set.

| Order | Dataset/work package | First practical deliverable | Current constraint |
|---|---|---|---|
| 0 | Every selected cohort | Sample/animal/donor/pool/assay crosswalk, conflicts and file inventory | Catalog availability is verified, payload integrity and independence are not |
| 1a | GSE154966 HPCS bulk ATAC | Within-source paired accessibility contrast from deposited counts | Four apparent pairing blocks; verify pooling and IDs before inference |
| 1b | GSE289683 + GSE291333; GSE141635; GSE150527 | Direct histone-mark profiles and normal-differentiation reference | Mostly one or two donors/preparations: descriptive, not replicated population inference |
| 2 | GSE273123 CD44-sorted RNA | Protein-enriched population contrast within genotype | Four apparent pairing blocks per genotype; verify animal and column IDs |
| 3 | Zenodo 10930946 / processed IMC resource | Protein-state and neighbourhood comparison at donor level | Obtain processed cell/ROI/donor table, masks and channel dictionary first |
| 4 | Existing GSE310539/GSE247130; GSE290014; GSE327686 | ATAC-derived state map with RNA labels held out | Pooled/single libraries, depth and conflicting metadata limit interpretation |
| 5 | GSE141259; GSE277777; GSE223302 | Time/lineage/perturbation comparison of frozen state axes | Recover sample identities; counts alone do not contain velocity layers or lineage events |
| 6 | PXD058626 and DDBJ perturbation deposits | Regional protein corroboration and treatment-sensitive chromatin | Processed matrix/cross-assay mapping unresolved; sequencing reprocessing may be expensive |

The source cohorts differ in species, genotype and experimental setting.
These work packages provide complementary evidence, not paired modalities from the same
individuals. A programme can be shared without its function being shared.

## Stage 0: metadata and identity contract

The runnable scaffold implements the metadata audit and a guarded paired
bulk-count pilot; see [scripts](scripts/README.md) and
[prospective contracts](config/README.md). Default execution reports holds.
The current 24 title-derived candidate sample rows have no verified biological
unit IDs; all three proposed contrasts remain unfrozen. This is an explicit
boundary between specification and analysis, not a failed biological result.

Use `metadata/<accession>.json` as dated source evidence, not as a curated
sample sheet. Build `config/samples.json` with one row per biological unit ×
assay × time × sort gate; connect technical libraries through explicit IDs.
Required fields are accession, sample ID, biological-unit ID, unit kind,
pool composition, donor/animal, genotype, sex/age if reported, injury/drug,
dose, harvest time, sorting/reporter gate, assay/mark, batch, genome build,
file URLs/checksums and evidence for any cross-assay pairing. Unknown stays null.

Use separate fields for days after injury, culture day and lineage-chase day.
An RNA and ATAC GSM from one multiome well are one unit, not two replicates.
An HTO library identifies samples; it is not a phenotypic protein panel.
Retain superseries/subseries overlap and shared-lab flags.

Resolve the explicit conflicts in the metadata README before affected
comparisons. Do not silently pick titles over characteristics, swap filenames
based on expected biology, or infer a mouse from a sequencing run name.
If sample identities cannot be resolved from public source material, keep
that contrast on hold and proceed with unaffected datasets.

Download processed counts, labels and small peak tables before FASTQs or
large fragment files. Produce a measured bytes/CPU/RAM plan after the first
file-format pilot. The 22 IMC MCD files alone total approximately 15.55 GB;
processed cell tables may avoid re-segmenting the complete archive.

## Stage 1: define states without circular validation

Freeze the complete published gene lists, direction, species, comparator,
selection rule and source supplement for DATP, PATS, ADI and HPCS separately.
The existing ES1 ADI list is available; its short DATP/PATS panels do not
become complete signatures. Recover source lists before promising a full
DATP-versus-PATS classification. Record overlaps and missing genes.

Use native author labels as one annotation view; build a second view from
modality-specific features. ATAC clustering/LSI and regulatory modules are
fit without RNA transition labels. RNA-to-ATAC label transfer may help
interpret a cluster, but it is not independent validation of that cluster.
Within paired multiome, compare RNA and ATAC assignments nucleus by nucleus.
Across unpaired assays, use aggregate correspondence, never invented cell pairs.

Exclude label-defining genes from validation panels; repeat ATAC comparisons
without promoters and nearby peaks of KRT8, CLDN4 and other defining markers.
Keep a prespecified stress/cell-cycle control and mature AT1, AT2, airway/basal
and immune negative-control populations. In cancer, account for CNV/genotype
and tumour purity before interpreting altered accessibility as epigenetic state.

Do not train a classifier to distinguish studies or species and call it a
state classifier. Where external biological labels exist, compare RNA-only,
chromatin-only and combined models using whole-donor/animal holdouts and nested
feature selection. Report balanced accuracy, macro-F1, calibration and an
unassigned category only if the number of independent units supports validation.
Otherwise publish descriptive correspondence and leave classification unrun.

## Stage 2: chromatin and direct modifications

1. **Accessibility:** use genome-matched peak coordinates and fragment-derived
   QC (TSS enrichment, FRiP where computable, nucleosome pattern, depth,
   doublets, blacklist). Freeze assay-specific QC before target contrasts.
   Construct a common peak universe within each study; recount fragments where
   possible. A filtered peak matrix cannot supply all missing QC metrics.
2. **Sample-level contrasts:** aggregate counts within animal/donor and state.
   Fit count-based models to genuine biological replicates, using within-unit
   pairing when verified. Adjust batch/genotype/time only when estimable.
   Do not use integrated/normalized embeddings as differential-test input.
3. **Regulatory programmes:** estimate GC/accessibility-matched motif deviations
   with [chromVAR](https://pubmed.ncbi.nlm.nih.gov/28825706/), and peak-to-gene
   candidates with [Signac](https://doi.org/10.1038/s41592-021-01282-5) where
   paired data allow. AP-1, TP53, NF-kB, TEAD, FOXA/NKX2-1 and CEBP families
   are literature-nominated candidates, not discoveries from this reanalysis.
   Motif families do not uniquely identify the active TF or prove its binding.
4. **Histone marks:** separate narrow promoter/enhancer marks from broad
   H3K27me3 domains. Use matching H3/input controls and documented library or
   spike-in normalization. Audit the p300-inhibitor study's spike-in factors
   before interpreting global H3K27ac changes; library normalization alone can
   conceal global changes. BigWig/bedGraph overlays remain descriptive if
   replicate-level counts cannot be recovered.
5. **DNA methylation:** GSE150527 provides a normal-differentiation reference
   with one WGBS donor at D0/D4/D6. Plot covered CpGs/regions and methylation
   fractions alongside histone/RNA changes. Do not run a donor-level DMR test
   or call D4 a purified DATP/PATS state solely from culture time.
6. **Enhancer interpretation:** distinguish promoter identity loss, candidate
   enhancer activation and measured repression. H3K4me1/H3K27ac combinations
   can support regulatory annotations where both are measured. Do not infer
   H3K27me3 or poised enhancers from absent ATAC signal.

GSE154966's initial paired ATAC contrast is TIGIT-positive versus negative
within the same verified source block. GSE273123's protein-sorted RNA contrast
is CD44-positive versus negative within genotype and verified animal, followed
by a genotype-by-sort interaction only if the design supports it. They test
different populations; neither is automatically an isolated DATP/PATS assay.

For the count pilot, the proposed model is `~ biological_unit_id + group`,
with positive-minus-negative as the fixed direction. The starting filter uses
edgeR `filterByExpr` (min.count 10, min.total.count 15), TMM and robust
quasi-likelihood fitting; review composition/background diagnostics before
freezing these settings. The [pilot implementation](scripts/04_fit_paired_counts.R)
saves every tested feature and its BH correction within the declared contrast.
RNA WT and mutant strata remain separate; testing a genotype interaction or
making a joint claim requires its own family/design, not comparing significance
labels. Leave-one-pair-out checks and assay-specific QC follow payload inspection.

## Stage 3: time, trajectory and measured lineage

GSE141259 is the time-course reference. Keep its whole-lung and enriched
epithelial preparations distinguishable and audit animal overlap. Use native
sampling times and broad lineage compartments to orient an initial PAGA/
diffusion or principal-graph trajectory. Compare an alternative topology and
test stability under animal exclusion, root choice, cell-cycle exclusion and
downsampling. A single attractive UMAP is not a topology test.

Use [CellRank 2](https://doi.org/10.1038/s41592-024-02303-9) time/pseudotime
kernels only as model-based transition hypotheses. RNA velocity requires
spliced/unspliced layers or an explicit raw-read reconstruction and kinetic
diagnostics. Do not manufacture velocity from ordinary expression counts.
Gene trends should be estimated with sample-aware models; cell bootstrap
bands do not express between-animal uncertainty.

In GSE277777, use trace origin, chase duration, sample hashes and treatment
history rather than an HPCS score alone. Compare observed reporter-positive
descendant distributions and the paper's measured ablation/growth endpoints.
Do not count pooled mice as individually observed replicates. A genetic fate
experiment reported in the paper is not automatically available per cell in
the deposited matrix. Normal repair and cancer histories remain separate.

Late tissue can reflect replacement, death, migration or different cells.
“Persistent” requires repeated population evidence; “irreversible” or
“epigenetic memory” needs an appropriate withdrawal/chase experiment.

## Stage 4: proteins, location and function

**Spatial proteins (IMC):** obtain antibody panel, segmentation masks, cell
features, donor/ROI and pathology labels. Distinguish intact AT2, AT1,
intermediate/basaloid and airway epithelium using multiple measured channels.
Reserve independent channels/features for validation so a KRT-defined gate
is not validated by KRT abundance itself. Test morphology and immune/fibroblast
neighbourhoods conditional on epithelial state, local density and pathology.
Use within-ROI spatial nulls preserving tissue compartments; summarise effects
within donor before inference. Regions labelled early/intermediate/advanced
are cross-sectional pathology, not a longitudinal trajectory.

**LCM-MS (PXD058626):** recover the protein abundance table and donor-region
map, including MUC5B genotype. Compare epithelium overlying fibroblastic foci
with prespecified epithelial reference regions using donor blocking when
matched. Audit epithelial/stromal admixture, missingness and detection frequency.
Report an observed-data analysis and an imputation sensitivity; do not
interpret imputed low abundance as protein absence. Regional enrichment is
supporting protein evidence, not single-cell identity or secretion rate.

**Functional anchors:** connect CD44-sorted expression with the source study's
conditioned-medium/coculture tests; Krt8 perturbation with GSE223302's time
course; and HPCS with lineage tracing/ablation. Analyse actual accessible
endpoint tables if their biological units are recoverable. Keep reported
source-paper experiments separate from computations performed here.
The SAGE Perturb-seq preprint is particularly relevant to bifurcating repair
and pathological states, but its public accession remains unverified.

No directly relevant epithelial CITE-seq/ADT dataset passed this targeted
scan. An alveolar macrophage CITE-seq study does not fill that gap. Whole-lung
proteomics or plasma biomarkers cannot establish a transitional-cell proteome.

## Stage 5: cross-study synthesis and decision rules

Within-study effects come first. Compare orthologous gene/regulatory modules
across species; report mapping loss. Genomic liftover is a sensitivity analysis,
not evidence that the same enhancer is conserved. Do not merge all matrices
and remove study effects when study and biological context are confounded.

Three verified independent biological units per arm/pair is the existing
minimum gate for a proposed inferential contrast, not proof of adequate power.
Two-replicate histone data remain usable for descriptive profiles and
concordance, with no promoted population claim. One-library scATAC can map
heterogeneity but cannot yield a replicated condition effect. No threshold
is relaxed to rescue a favoured dataset.

Freeze contrast families, feature universes and QC before estimation.
For eligible count models report effect sizes, uncertainty and BH FDR across
the tested family. Pathway tests use the assayed background; estimated gene
correlation remains primary where CAMERA is used. Depth sensitivities,
leave-one-unit-out results and alternative roots/annotations are labelled.

A regulatory distinction is supported only when it survives relevant QC and
has reproducible direction at the independent-unit level. Distinct protein
or lineage evidence strengthens interpretation, but unrelated cohorts do not
constitute same-cell multiomics. Failure to distinguish states may indicate
continuity, inadequate power, missing assays or technical incompatibility;
those outcomes must be separated.

**Stress-test of the intended conclusion:** none of the currently audited
designs measures every named transitional population side by side with all
modalities. The first analyses can establish within-context contrasts and test
transport of frozen programmes; they cannot by themselves settle a universal
DATP/PATS/HPCS taxonomy. Histone comparisons with two preparations must not be
used as independently validated classification targets. Additional matched,
replicated assays would be needed for that stronger claim.

## Figure plan

| Figure | Panels and scientific purpose | Gate |
|---|---|---|
| A1-1 | Sample/assay design; RNA and ATAC embeddings shown separately; state correspondence with unassigned cells | Correct sample mapping; UMAP descriptive |
| A1-2 | Sample-level ATAC PCA; paired accessibility effects; motif/module heatmap; selected locus tracks | Genuine units for effect inference; loci frozen or labelled exploratory |
| A1-3 | H3K27ac/H3K4me3/H3K27me3 tracks; mark-by-state heatmap; methylation coverage and fractions | Direct measured marks; low replication shown explicitly |
| A1-4 | Chronological state fractions; trajectory/topology sensitivity; observed lineage-descendant matrix | Temporal/lineage metadata; arrows labelled inferred or measured |
| A1-5 | Tissue images; protein distributions; morphology and neighbourhood effects per donor | Antibody/segmentation/ROI audit |
| A1-6 | Regional proteomic heatmap and paired effects; cross-assay evidence matrix | Protein table, donor map and detection audit |

Violin plots may display cells, but inference and error bars use biological
units. Illustrative tissue fields are chosen by a declared rule, not maximum
effect. Avoid an omnibus cross-species UMAP or a diagram that implies causal
arrows unsupported by perturbation data. Figures are produced only after
analysis, with source tables, code hashes and visual review.

## Deliverables and stopping points

Stage 0 ends with a verified sample manifest, explicit holds, input hashes,
estimated resource cost and a frozen primary contrast. Subsequent stages
write `tables/`, `figures/` and `reports/` here, preserving prior run records
on rerun. The canonical A1 summary links to this gallery; source-paper
outputs remain under `Research Article/`.

First review point: resolve the HPCS/IMC/CD44 identifiers and select the
executable contrasts. Later batches can proceed independently where gates
pass. Downloading all raw reads, claiming a complete DATP/PATS/HPCS taxonomy,
or launching every method listed above is not the present plan.
