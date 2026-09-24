# IL-1beta, context-specific niche circuits and epithelial plasticity

**Version 3, 24 September 2026. Owner-authorized staged execution.**
The owner requested final self-review, improvements and launch. The
[adversarial review](FINAL_REVIEW.md) records the fixes and remaining scientific
gates. Acquisition/QC has started on the seven early two-arm libraries;
annotation and design checks control release of biological comparisons. U0
metadata and U1 input checks were preparation, not hypothesis tests. This is a
prospective specification for newly acquired expression data, informed by
published findings and previously inspected repository data; it is not an
unseen-data preregistration.

## 1. Research questions and discriminating observations

The main question is which macrophage-fibroblast-epithelial circuits and
recipient programmes accompany supportive versus persistent pathological
niches, and how their responses to IL-1beta blockade relate to epithelial
state balance. The [required niche analysis](NICHE_ANALYSIS_PLAN.md) supplies
the explicit ligand-receptor, compartment-specific pathway and ligand-target
methods missing from version 1. Its context map follows the Body annotations:
molecular regulation, source/recipient identity, regenerative route, niche
condition, age/history, organ/disease stage and spatial/tissue outcomes.
These dimensions are recorded separately; unmeasured ones remain untested.

| Priority | Question | Testable contrast | Evidence that would weaken the proposed link |
|---|---|---|---|
| Central niche aim N1 | Which source states and recipients carry candidate IL-1 signals? | Fixed macrophage/monocyte-to-fibroblast and epithelial ligand-receptor comparisons; epithelial sources and IL1A comparator measured separately | Cell mixture explains the edge; recipient response is absent; inhibitory context differs |
| Central niche aims N2/N3 | How do fibroblast inflammatory, trophic and matrix programmes relate to reciprocal myeloid/epithelial signalling? | Actual LR inference plus fibroblast and macrophage pseudobulk pathway enrichment, within each repair, IPF or cancer context | Matrix and inflammatory programmes separate, an edge lacks receiver support, or a cross-organ hypothesis is not detected in eligible lung data |
| Perturbation aim N4 / primary phenotypic endpoint | Does direct IL-1beta blockade alter niche responses and epithelial state balance together? | Anti-IL-1beta versus IgG in GSE300288; per-animal KAC fraction plus prespecified LR and recipient pathway families | Niche responses and KAC fraction change independently; annotation/QC explains the epithelial shift. Unchanged IL1B RNA alone does not weaken protein neutralization |
| Secondary | Is a plastic-state inflammatory niche reproducible in human precursor lesions? | Within-patient normal/precursor/LUAD contrasts in GSE308103; paired lesion context in GSE307534 | Association is explained by immune abundance, sequencing depth or shared score genes; lacks within-patient consistency |
| Context aim N5 | Which circuit/response associations differ between repair, fibrosis and neoplasia? | Within-study contrasts and subtype correspondence; age, genotype, injury history and assay retained as context | Similar patterns occur during ordinary injury or apparent differences follow sampling; fibrosis and cancer cannot be placed on one inferred progression axis |
| Specificity | Do DATP/ADI/PATS, KAC and HPCS definitions distinguish contexts beyond generic stress? | Independent source-defined programmes across repair, fibrosis and cancer; developmental ISR comparator if reusable data exist | Correspondence disappears after shared-gene removal or is equally strong in development/mitochondrial stress |
| Experimental follow-up | Does a persistent state become independent of the initiating IL-1beta cue? | Requires exposure, withdrawal/blockade and fate/chromatin measurements with biological replication | The state resolves with cue withdrawal, or persistence is explained by survival/selection rather than memory |

The last question is central to the review but is **not currently answerable
as a definitive public-RNA analysis**. There is no verified common dose or
exposure-duration series here from which to estimate a universal IL-1beta
threshold. Lack of a treatment response will not prove autonomy or memory.

## 2. Data sequence and role assignment

1. **GSE300288: perturbation-first core.** Public 3-month groups contain four
   IgG libraries, three anti-IL-1beta libraries, four anti-PD-1 libraries and
   four combination libraries. At 7 months there are four in each group.
   These are deposited library counts pending animal/pool verification.
2. **GSE308103 and GSE307534: human context and spatial extension.** Metadata
   provide 23 snRNA patient labels, with 22 candidate normal/precursor/LUAD
   triplets, and 25 spatial patient labels, with 24 candidate precursor/LUAD
   pairs. Verify tissue IDs before joining. These assays share a study/cohort.
3. **GSE267226/GSE267228: pathological post-injury context.** Three human
   PASC-PF versus two controls; mouse anti-CD8 versus IgG, two mice per arm.
   Use descriptive donor/animal contrasts, not spot-level replication.
4. **GSE141259 plus existing Choi/Niethamer data: injury reference.** Check
   animal identity, time, experiment and sort strata. Do not rerun a generic
   Krt8 trajectory already established elsewhere. Its published ADI signature
   is already used by this repository, so it is a definition/reference cohort.
5. **GSE136831/GSE135893: fibrotic niche arm.** Use existing human cohorts for
   macrophage/fibroblast LR and pathway questions after donor/subtype coverage
   checks. Extract missing compartments from full sources; prior epithelial-
   stromal correction caches omit macrophages. This is an extended question
   in previously inspected data, not unseen validation or a treatment substitute.
6. **GSE277777: malignant HPCS specificity extension.** Resolve GEX/HTO and
   animal demultiplexing before use. This is a comparator, not independent
   evidence that IL-1beta controls the HPCS.

The [dataset inventory](DATASETS.md) includes additional companions, existing
IPF cohorts and outstanding source retrievals. No single integrated matrix
will mix disease, species, treatment and assay into one latent ordering.

## 3. U2: establish the units and freeze the analysis inputs

**Required before biological scoring.** Create a sample manifest with series,
GSM, file hash, species, assay, animal/donor, pool size, tissue/lesion, treatment,
start/end, age, sex, preparation batch, capture and technical-replicate IDs.
Carry raw deposited values beside every interpreted field and its source.

For GSE300288 specifically:

- Reconcile the FFPE `cell line` characteristic with the single-cell library
  source and Cell Ranger processing description. A passing MatrixMarket check
  verifies format, not whether that metadata field is erroneous.
- Resolve library-to-animal and pooling from methods/supplementary metadata.
  An ordinal such as "replicate 2" is not evidence of a paired animal design.
- Explain the missing early anti-IL-1beta replicate-1 library; do not impute
  the missing library or silently rebalance groups.
- Verify treatment duration and start. A 7-month endpoint following earlier
  treatment does not establish efficacy of starting treatment in established
  cancer. The age difference between endpoints also prevents treating their
  contrast as pure progression.
- Join lesion burden, histology or function to the same animals only where
  identifiers support it. Paper-level rescue results are contextual evidence,
  not automatically matched outcomes for these libraries.

For the human companions, preserve normal, AAH, AIS, MIA and LUAD separately.
Repeated lesions/sections are nested in patients. Resolve conflicting lesion
labels using pathology/source tables; title-derived links are provisional.
The 23 RNA and 25 spatial patient labels cannot be added as independent n.
GSE307529 is a public WES companion for potential lesion-genotype crosswalks;
it is not an RNA dataset and is not in the initial download queue.

**Release gate:** a design is eligible for inferential testing only when its
biological units are resolved, the contrast is estimable, and at least three
independent units per group meet state coverage. This floor is a project
eligibility rule, not a power calculation. Below it, report directions and
coverage without inferential P values. If the core fails, proceed to an
eligible human within-patient association analysis and explicitly drop the
causal treatment claim.

## 4. U3: establish cell states and assay comparability

Begin with deposited processed counts and metadata. Keep their filtering
history explicit; raw integer counts need not be raw unfiltered droplets.
Use sparse/streaming reads and paper-local caches. Generate per-library QC:
counts, detected genes, mitochondrial fraction, barcode calling and doublet
evidence. Reuse the repository's established QC helpers; fresh nuclei and
FFPE targeted assays need assay-specific checks rather than one shared filter.
Ambient correction requires the appropriate empty-droplet data. Record when
that is unavailable. Do not install or alter the existing analysis stack
without an actual implementation need.

Keep unsupervised fits separate from deposited labels. Use source annotations
and HLCA/marker evidence deliberately in the later compartment analysis,
with an explicit label source, uncertainty and sensitivity table. Maintain
AT2, AT1, alveolar transitional, secretory, ciliated and KRT5/TP63 basal
populations separately. Keep malignant, nonmalignant and ambiguous epithelial
labels; CNV inference is corroboration, not sufficient proof of malignancy.
Never call a marker gate a cell type without validation.

Build a versioned programme manifest before comparisons:

| Measurement | Source/definition | Important separation |
|---|---|---|
| AT2, AT1 and ADI identity | Reuse hashed full reported Strunz lists in [existing modules](../epithelial_state_specificity/modules.json) | Source-capped lists, not exhaustive programmes; same-study recovery is not external validation |
| DATP/PATS | Existing short source panels as explicitly labelled proxies; acquire full published lists if available | Do not relabel short panels as full signatures |
| KAC | Obtain the exact published Peng/Han programme and annotation implementation from their supplement/code | No outcome-driven gene selection; independent marker sensitivity and uncertain/unassigned cells retained |
| HPCS | Chan supplementary Table 5 and [author code](https://github.com/dbetel/HPCS_LUAD) | Malignant plasticity, not an interchangeable name for DATP |
| NF-kB, hypoxia, ISR, cell cycle, p53, EMT | Pinned published gene sets; ISR source from Han supplementary Table 1 | NF-kB is a shared response; ISR and glycolysis RNA do not measure flux |
| IL-1 components | IL1B/Il1b, IL1A/Il1a; IL1R1/Il1r1 and IL1RAP/Il1rap; IL1RN, IL1R2 and SIGIRR; processing machinery separately | A source-defined component panel, not a validated activity signature |
| Niche programmes | Source-defined macrophage/monocyte and fibroblast subtypes; separate inflammatory, metabolic, matrix and trophic programmes in [the niche specification](NICHE_ANALYSIS_PLAN.md) | Required sample-level LR and pathway analyses; gene expression is not ligand secretion, receptor activation or matrix stiffness |

Record source, species, orthologue mapping version, measured-gene coverage
and every overlap. Use one-to-one orthologues for explicitly cross-species
comparisons. Compare per-study standardized effects, not raw score magnitudes.
Withhold a full-signature result below a proposed 70% assayed-gene coverage;
report the missing genes and a 50%/80% coverage sensitivity. Freeze this rule
and any state-call cutoff before expression comparisons. Never silently shrink
a gene set to the few genes present on a spatial panel.

Remove assignment genes from within-state outcome scores. Report full and
overlap-removed programme scores separately. Use sample-stratified matched
expression/detection gene sets, depth matching and cell-cycle/stress controls
to distinguish specificity from library complexity and common genes.
Also test KAC assignment after excluding IL-1/NF-kB response genes: a treatment
can change a classifier's inputs without eliminating a cell population. If
classification depends on the perturbed response itself, interpret the primary
fraction as a change in assigned transcriptional phenotype, not cell depletion
or differentiation. Inspect model-matrix rank and condition/batch overlap
before any covariate-adjusted pseudobulk fit.

## 5. U4: direct-perturbation analysis

**Primary contrast:** anti-IL-1beta alone versus IgG at the 3-month endpoint.
**Primary phenotypic endpoint:** per-animal KAC fraction among eligible alveolar
epithelial cells, with frozen KAC definitions from U3. If an independently
validated KAC call cannot be recovered, this endpoint is not evaluable; a
continuous score is a separately labelled secondary analysis, not a replacement.

The two prespecified contrasts are the 3-month and 7-month comparisons, with
the former leading interpretation. Report the mean fraction difference in
percentage points, each animal's value, denominator, and leave-one-animal-out
range. Test labels at the biological-unit level only if exchangeability and
design support it; enumerate permissible allocations for these small groups.
The leave-one-out range is a robustness diagnostic, not a confidence interval.
If random allocation is undocumented, any permutation inference is exploratory
and conditional on exchangeability. Control multiplicity across the two
endpoint contrasts with Holm adjustment. Effect size and uncertainty lead;
small-n significance is not required to keep an informative null result.
The final review establishes a stronger limit: at the deposited 3-versus-4
and 4-versus-4 sample sizes, the two-sided absolute-difference permutation
tests cannot pass Holm across these two endpoints at 0.05. Retain the family,
report attainable P values, and interpret effects/robustness without claiming
equivalence or switching tests to obtain significance.

Proposed initial coverage rules: at least 100 alveolar epithelial cells per
animal for the state fraction; at least 30 cells of a given state for its
pseudobulk. Fractions with zero KAC cells remain valid zeros when denominator
coverage passes; they must not be excluded by the 30-cell state floor. Show
50/100/200 denominator and 20/30/50 state-floor sensitivity without choosing
the most favourable threshold. Do not infer whole-lung abundance from sorted
or assay-biased composition.

Parallel niche arms and additional phenotypic endpoints distinguish mechanisms:

1. AT1/AT2 maturity and KAC programme expression within comparable epithelial
   states, pseudobulked by animal. Fit a count model with library normalization
   and estimable covariates; use unintegrated counts, not corrected embeddings.
2. Required parallel niche arms: **U4-LR**, sample-specific LIANA inference
   and biological-unit RNA-compatibility contrasts; **U4-PW**, macrophage and
   fibroblast pseudobulk enrichment with TMM/voom/CAMERA; **U4-LT**, conditional
   NicheNet receiver-target prioritization. Methods, context-specific candidate
   edges, eligibility and multiplicity are fixed in [NICHE_ANALYSIS_PLAN.md](NICHE_ANALYSIS_PLAN.md).
   Compare abundance separately from within-subtype changes. The proposed
   niche floor is 50 cells/subtype/sample (30/50/100 sensitivity), distinct
   from the 30-cell epithelial-state floor above. No absent compartment is zero.
3. Cell death, proliferation and sampling sensitivity. A lower KAC fraction
   with no mature-cell gain is compatible with depletion, not proof of repair.
4. The four treatment arms allow a secondary anti-IL-1beta-by-anti-PD-1
   interaction only if sample structure supports it. Combination versus IgG
   alone cannot identify an IL-1beta-specific effect.
5. Endpoint-specific treatment responses may be compared descriptively.
   Because ages and exposure histories differ, do not claim the contrast
   identifies an irreversible therapeutic window.

For gene-level secondary DE use BH FDR within the declared contrast family;
for pathway testing use official CAMERA with residual inter-gene correlation
estimated, following repository corrections. Declare separate global BH
families for core LR edges and primary pathways across subtypes and both
endpoints; human and exploratory screens have separately labelled families.
No cell-level P values, no
gene-permutation null that ignores correlation, and no hidden covariate
adjustment when age/treatment/batch cannot be separated.

## 6. U5: human and spatial association analyses

Use GSE308103 first for paired per-patient programme contrasts. Within each
patient and histology, aggregate repeated tissue samples without counting
them as extra patients; preserve lesion-specific estimates as sensitivity.
Analyse AAH, AIS and MIA separately where coverage permits. A pooled precursor
estimate is a declared secondary analysis. Patient-blocked contrasts address
between-person differences; cross-sectional lesions still do not prove a
normal-to-cancer trajectory or common ancestry.

Run **U5-LR/PW/LT** with patient/lesion-aware units and the same context-aware
methods. Normal-to-AAH, normal-to-AIS, normal-to-MIA and normal-to-LUAD paired
contrasts remain separate, with multiplicity across eligible contrasts.
Freeze one source-defined KAC measure and each receiver's NF-kB response
measure for the joint association; other pathway associations are exploratory.
Report marginal paired contrasts and depth/mixture/shared-gene sensitivities.
Associations use patient-level lesion-minus-normal differences and require
at least ten complete patients with both measurements; below that display
points without a correlation test. Use Spearman association of the paired
differences, with patient-level permutation and BH across the frozen receiver
subtypes and eligible histology contrasts; report leave-one-patient-out stability.
Do not pool histologies to increase n. Shared tissue pieces do not add patients.
Avoid a many-covariate model in roughly two dozen patients.
If prediction is explored, use nested patient-level
holdout and compare against histology/mixture-only baselines; do not promise a
classifier as a required output.

For GSE307534, obtain count matrices, tissue masks, spatial coordinates and
histology/lesion annotations. Define anatomical regions independently of the
tested gene scores to avoid selecting a KRT8-high ROI and then reporting KRT8
enrichment as a finding. Visium spots contain mixtures; describe region-level
co-localization. If verified Xenium data are added later, audit panel coverage,
segmentation and cell identities before measuring cell distances. Spatial
permutations stay within section/lesion and preserve tissue structure; animal
or patient is the replication unit for the final contrast.

GSE267226/267228 then asks whether the direction recurs in pathological repair.
Three diseased human donors and two controls, and two mouse samples per arm,
limit conclusions to descriptive context. Missing spatial coordinates block
neighbourhood inference even if expression matrices pass. Anti-CD8 effects
cannot be reported as an anti-IL-1beta transcriptomic experiment.

## 7. U6: specificity, synthesis and next experiments

Use injury and IPF references to ask whether the chosen KAC/HPCS programme
contains information beyond shared stress/transitional genes. GSE222901 and
GSE300293 contextualize the core mouse study; they are not automatically
independent confirmations, because specimens and prior analyses can overlap.
GSE277777 provides a later malignant comparator after hashed-animal recovery.
The Han developmental ISR study is conditional on obtaining reusable processed
data; its listed BioProjects are not permission to start an unbounded FASTQ run.

Report evidence by context, species, independent cohort and intervention.
Use within-study effect estimates; consider a meta-analysis only if estimands,
assays and biological units actually match. Do not treat all transitional
labels as the same state or interpret a shared score as fibrosis progressing
to cancer.

A decisive future experiment would vary IL-1beta exposure duration and
withdrawal, epithelial versus fibroblast recipient context, and relevant
genotype; include mature fate/lineage, viability, protein release and chromatin
readouts from independently prepared cultures/animals. This is a proposed
experimental design, not a result recoverable from current RNA snapshots.

## 8. Deliverables, resource limits and completion

Each released stage writes its own specification, run record, input/code hashes,
sample-exclusion table, result tables and a short interpretation with alternatives.
Plots use the repository palette and show individual biological units.
Keep new figure galleries with this paper's analysis outputs and link them
from the paper README. The owner chose to avoid a partial gallery on the root
README; keep that landing page as navigation to the paper-level galleries.

Planned figures: (1) context/sample/subtype coverage; (2) per-animal niche and
epithelial perturbation contrasts; (3) directional LR matrices; (4) separate
macrophage/fibroblast pathway-enrichment heatmaps; (5) within-patient human
contrasts and eligible ligand-target support; (6) spatial region associations;
(7) programme overlap and context specificity. UMAPs are
annotation diagnostics, not required evidence figures.

Acquire only the first eligible contrast after review, then expand in stages.
Estimate compressed sizes, expanded storage and peak memory from the two
existing examples and file headers before full downloads. Stream one library
at a time on this workstation. Avoid downloading whole superseries TAR files,
raw sequencing or the very large spatial Zenodo collection by default. No
runtime estimate is asserted before the pilot measures it.

Completion means every prespecified contrast is either evaluated with valid
units and its sensitivities, or explicitly classified as unidentifiable with
the exact missing information. Null results, absent state coverage and failed
specificity tests remain outputs. Only then propose claim-register additions
for review; existing historical claims remain unchanged.

The review decisions are the context-specific triad questions, the perturbation
and fibrotic-niche dataset roles, the proposed coverage/multiplicity rules,
and whether the spatial
and HPCS extensions belong in the initial run. The owner authorized launch after final review; scientific release gates
remain in force. See FINAL_REVIEW.md and the stage run records.
