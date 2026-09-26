# A0 pilot analysis plan

Date: 2026-09-25. Status: **draft design; P0 initial audit executed; P1 not frozen**.

Execution status: [P0 eligibility report](reports/P0_ELIGIBILITY_REPORT.md).
The initial audit retains the original floors and reports unresolved selection
gaps. No expression-based discovery or transfer has been performed.

This is an exploratory plan informed by published studies and previous local
analyses, not an unseen-data preregistration. Record subsequent decisions and
amendments with their dates and whether relevant results had already been seen.

## 1. Question, hypotheses and limits

**Primary question:** Does an intermediate-associated transcriptional program
transfer from lung repair and normal lung development to differentiation in
another epithelium?

| Explanation | Discriminating observation |
|---|---|
| Conserved transition-associated program | A frozen program is enriched in intermediates relative to both endpoints across the three settings, with specificity beyond stress and proliferation |
| Conserved repair response | Enrichment follows injury or stressed cells, with weak evidence in normal differentiation |
| Tissue- or trajectory-specific programs | Reproducible within-context effects fail to transfer despite adequate measurements |
| Apparent intermediate caused by endpoint mixtures or technical effects | Signal depends on mixed identities, doublets, labeling genes, depth or integration |

The primary analysis tests a shared **transcriptional state program**. A shared
process could instead involve different genes, protein regulation or conserved
temporal ordering. Failure here does not reject every version of the biological
hypothesis. Intermediate enrichment is not a necessary property of all causal
regulators; monotonic changes are outside the primary endpoint.

## 2. Three-setting design

| Role | Setting | Required comparison |
|---|---|---|
| D1: discovery | Lung injury and regeneration | Starting state, supported intermediate and differentiated destination within an interpretable injury/time stratum |
| D2: discovery | Normal lung development | Starting progenitor, supported intermediate and destination along one specified developmental branch |
| V1: transfer | Normal differentiation in another epithelium | Independently supported intermediate and both endpoints; an intestinal or skin lineage is a candidate |

Use one clearly defined branch per dataset. A developmental progenitor need not
be the same cell type as the injury starting state. Declare each trajectory
explicitly; do not equate fetal, immature, proliferating and transitional cells.

Prefer the same species across all three settings to reduce orthology and
assay confounding. If unavoidable, freeze one-to-one ortholog mapping and report
gene loss before the transfer test. A second injury-only tissue would test
repair conservation, not replace the normal-differentiation question silently.

Keep expression matrices and preprocessing separate initially. Cross-dataset
integration is not required for the primary test.

## 3. Dataset eligibility

Before selecting cohorts, audit metadata and coverage without inspecting the
candidate program's association with the transfer labels.

- Accessible expression matrix, documented normalization or raw counts, gene
  identifiers, biological-unit IDs, library IDs and state annotations.
- Documented evidence for the proposed transition from the source study:
  lineage tracing, time-resolved experiments or functional evidence is preferred.
  Pseudotime placement alone is insufficient to establish biological transition.
- At least **three independent biological units for each primary contrast**.
  This is a feasibility floor, not a power calculation. Pooled animals in one
  library contribute one observable unit; technical replicates add no biological n.
- Default minimum of **30 retained cells per state per biological unit**.
  Prefer units containing all three states; audit the two pairwise contrasts
  separately if complete triplets are unavailable.
- At least **80% of the frozen program's genes represented in the assay**;
  also report actual expression/detection. Feature presence alone is insufficient.
- State contrasts must be distinguishable from batch, genotype, injury and
  developmental stage. If all intermediates occur in one batch or stage and all
  endpoints in another, the proposed state effect is not independently identifiable.
- A source-supported stressed epithelial comparison group, either within the
  selected repair study or a justified control cohort. Its absence limits the
  claim of specificity; do not invent a non-transition label from a low module score.

Unpaired comparisons require replicated groups and an identifiable design; they
are secondary to paired comparisons in this pilot. Small groups and failed
coverage remain visible in the audit. Do not relax floors after seeing effects
to obtain a preferred answer. Numerical defaults may be revised for documented
design reasons before analysis and must be retained in the decision record.

## 4. Execution sequence

### P0 — Audit and select datasets

1. Complete [DATASET_CANDIDATES.md](DATASET_CANDIDATES.md), including accessions,
   source links, biological units, previous exposure and transition evidence.
2. Audit available processed matrices before considering new raw-read processing.
   Reuse local data only when their design fits the question.
3. Tabulate retained cells by sample, state, time/stage, condition and library.
   Check the eligibility requirements independently of program scores.
4. Select one dataset per role. If no eligible V1 is available, report feasibility
   unresolved; another lung cohort cannot establish cross-tissue transfer.

### P1 — Freeze state definitions and analysis settings

1. Define starting, intermediate, destination and stressed-control groups from
   source evidence and annotations. Record uncertainty and all genes used to
   assign labels, including known marker panels behind author annotations.
2. Use author labels explicitly as comparison definitions. This differs from an
   unsupervised atlas-recovery exercise; do not describe label use as blind.
3. Freeze epithelial selection, QC, normalization, gene filtering, scoring,
   sample-level contrasts and sensitivity analyses. Apply per-sample QC and record
   retention; inspect likely doublets and ambient contamination.
4. Select V1 using design/coverage only, then freeze its role. Its outcomes must
   not influence gene selection, weights or score cutoffs. Record any prior
   inspection; a previously analyzed cohort is a transfer check, not pristine validation.
5. Save the configuration, input hashes, software versions and decision record.
   Use the repository's Python environment and seed 0 unless documented otherwise.

### P2 — Discover one candidate program in D1 and D2

1. Form count pseudobulks by biological unit and state when counts are available.
   Estimate intermediate-minus-start and intermediate-minus-destination effects
   separately within the declared time/stage/condition strata.
2. Rank genes by consistent positive effects across **both endpoints in both
   discovery settings**. Use equal biological-unit weighting and sample-level
   estimates; do not rank by cell-level significance.
3. Use a bounded consensus-ranking procedure to select at most 50 genes. Require
   at least 20 eligible genes for the primary module; otherwise report that no
   sufficiently supported common program was identified. Freeze exact expression
   filters and ranking/tie rules in P1, before this discovery step.
4. Exclude state-labeling genes and any known V1 labeling genes from the primary
   instrument before selection. Record overlap among candidate genes, endpoint
   identities and source marker panels. Removing a few label genes does not fully
   remove annotation circularity; independent biological evidence still matters.
5. Keep one primary program. Published ADI/DATP/PATS panels and pathway sets are
   benchmarks with documented provenance, not additional independent discoveries.
   Short panels and full reported marker lists must remain distinguishable.
6. Freeze gene membership and scoring before examining V1 outcomes. If an
   alternative module is developed after V1 inspection, it is exploratory and
   needs a new transfer dataset.

### P3 — Score and transfer

1. Apply the same frozen gene set and score definition separately in each dataset.
   A within-cell rank-based enrichment score is the proposed primary instrument;
   finalize its method and gene-universe handling during P1. Report coverage.
2. Aggregate scores within each sample/state, then calculate the two paired
   differences. Report every unit, median effects, dispersion and eligible n.
   Never use cell resampling as biological uncertainty.
3. Standardize within each dataset using an explicitly frozen reference rule
   when comparing effect magnitudes. Do not interpret raw scores from different
   assays as identical biological units.
4. Evaluate V1 once using the frozen primary analysis. Report both endpoint
   contrasts; enrichment over only the starting state can reflect destination
   maturation rather than a distinct intermediate-associated program.
5. Evaluate uncertainty at the biological-unit level. With three or four units,
   show the individual estimates and leave-one-unit-out stability without a strong
   inferential claim. For larger cohorts, prespecify a sample/block bootstrap or
   appropriate sample-level model. Report multiplicity for secondary modules.

### P4 — Specificity and technical sensitivity

| Challenge | Planned check |
|---|---|
| Generic injury/stress | Score independently sourced inflammatory, p53, hypoxia and immediate-early/AP-1 programs; compare with source-supported stressed epithelial controls |
| Proliferation | Report cell-cycle activity and repeat contrasts in comparable cycling strata where supported |
| Simple loss of starting identity | Retain both endpoint contrasts and score starting/destination identity separately |
| Annotation circularity | Use label-excluded primary genes and assess alternative source-supported state definitions |
| Sequencing depth and quality | Inspect score-versus-depth/mitochondrial relationships; repeat at comparable depth or count budgets where feasible |
| Endpoint mixtures/doublets | Check mixed-marker cells and doublet flags; compare with simulated endpoint mixtures as a technical sensitivity, not as biological replicates |
| Generic gene-set behavior | Compare with fixed size- and expression-matched gene sets; these assess score specificity, not biological replication |
| Dominance of one sample | Display each biological unit and leave-one-unit-out effects |

Do not automatically regress away stress, cell cycle or YAP-related activity:
these could participate in transition. Present unadjusted results and targeted
sensitivity analyses side by side. A stress-free gene list is not the hypothesis.

Incremental sample-level modeling against stress, proliferation and depth is
secondary and only justified when independent n and design support the covariates.
Otherwise use stratified/descriptive comparisons and leave specificity unresolved
where necessary. Do not fit an overparameterized classifier to a few donors.

## 5. Decision criteria

These are pilot investment criteria, not proof thresholds. Freeze any numerical
effect-size target during P1; do not choose one after seeing V1 results.

| Decision | Evidence required |
|---|---|
| **Expand** | Adequate coverage and replication in all roles; positive V1 effects against both endpoints in at least two-thirds of eligible units for each contrast, positive median effects, and no sign reversal of either median when leaving out one unit; substantive specificity controls support more than generic stress/proliferation |
| **Narrow** | Reproducible association restricted to injury, lung, or one trajectory; rewrite the hypothesis around that supported scope |
| **Stop this transcriptional-program pilot** | Technically adequate, identifiable comparisons show no reproducible common program or no meaningful transfer under the frozen analysis |
| **Unresolved** | Low coverage, confounding, uncertain labels, missing controls or wide uncertainty prevent a biological decision |

Directional consistency alone is insufficient for expansion: retain effect
magnitudes, uncertainty and control results. Statistical non-significance alone
does not establish absence. Failure to find one shared gene module cannot rule
out shared regulatory dynamics implemented by different genes.

## 6. Outputs and work boundaries

Create these artifacts only as their corresponding stages run:

- `tables/dataset_audit.csv`: sources, biological units, contrasts, eligibility,
  coverage and exclusion reasons.
- `tables/program_genes.csv`: selected genes, discovery evidence, exclusions and
  benchmark overlap; accompanied by frozen configuration and input hashes.
- `tables/sample_effects.csv`: unit-level scores and both endpoint contrasts.
- `tables/specificity_checks.csv`: control and sensitivity results.
- Three figures: (1) design/coverage audit, (2) sample-level endpoint contrasts
  across contexts, (3) specificity and sensitivity results.
- `reports/PILOT_REPORT.md`: observations, uncertainty, decision and the next
  discriminating experiment; no causal or universal claim from RNA alone.
- `decisions.json`: dated selections, amendments, prior data exposure and reasons.

Cap the first pass at three core datasets and one candidate program. Reuse ES1
provenance and measurement code where appropriate. Start with processed data;
do not expand into atlas integration, raw-read reprocessing, velocity, GRN
reconstruction or perturbation screening without a concrete unresolved question.

A positive pilot motivates another independent tissue and temporal or perturbation
evidence. Distinguishing transition entry from exit requires those additional
measurements; intermediate abundance alone cannot identify either rate.

## 7. Evidence and provenance

- [Ke et al., lung morphogenesis and regeneration](https://pmc.ncbi.nlm.nih.gov/articles/PMC11945641/): precedent for developmental reuse within lung; not evidence of universality across epithelia.
- [Strunz et al., Krt8 transitional state](https://www.nature.com/articles/s41467-020-17358-3): lung regeneration context and reported marker lists.
- [Yui et al., YAP/TAZ-dependent colonic reprogramming](https://pmc.ncbi.nlm.nih.gov/articles/5766831/): mechanistic context for intestinal repair; dataset eligibility remains a separate question.
- [Ge et al., skin lineage infidelity](https://pmc.ncbi.nlm.nih.gov/articles/PMC5510746/): functional plasticity during skin repair; not a substitute for a verified normal-differentiation cohort.
- [ES1 plan](../../Research%20Article/epithelial_state_specificity/PLAN.md) and [results](../../Research%20Article/epithelial_state_specificity/results/SUMMARY.md): reusable local resources and known measurement/replication limits.
