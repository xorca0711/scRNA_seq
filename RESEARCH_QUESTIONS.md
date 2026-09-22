# Research questions and the measurements that can answer them

This repository uses public lung single-cell and multiome data to generate
hypotheses and test whether computational observations survive changes in
annotation, measurement and sample composition. Deposited counts are used where
available; E1 instead uses deposited normalized expression. Author labels were
held out of the original unsupervised atlas fits, but are explicitly used in
later compartment pseudobulks and annotation sensitivity analyses.

The current computational question is:

> Which epithelial and macrophage programme changes repeat across independent
> samples after accounting for cell-state composition, genotype and measurement quality?

The motivating biological follow-up is:

> Which of those changes associate with mature AT1 contribution or persistent
> pathological remodelling in cohorts with independently measured outcomes?

These questions have different evidence requirements. The current collections
have no shared, independently measured repair outcome. A late time point is
not proof of recovery, and a persistent transcriptional state is not proof of
pathology. Fibrosis, tumour initiation and viral injury are separate contexts.

The [claim register](CLAIMS.md) records the detailed evidence. Its
[generated summary](docs/CLAIM_SUMMARY.md) distinguishes scientific status,
review authority and numeric-check coverage. A run record specifies the next
pass; it does not establish that the data or question were previously unseen.
The [September audit](docs/audits/2026-09-22/REPOSITORY_REVIEW.md) and
[implementation record](docs/remediation/2026-09-22/IMPLEMENTATION_STATUS.md)
explain the corrections. Historical trial outputs remain available.

## Part A. Questions, evidence and decision limits

### A1. Which RNA and chromatin changes accompany transitional epithelial states?

A1 and A5 now form one epithelial-state specificity project. The first task is
to measure identity and state programmes while distinguishing injury,
development and genotype. The second is to ask whether those measurements
associate with an independently observed fate. Accessibility cannot substitute
for that fate measurement: productive AT2-to-AT1 differentiation can also lose
AT2 identity, while disappearance of a transitional state can reflect
maturation, death or replacement rather than reversal into AT2.

Across GSE310539 and GSE247130, the operational Cldn4/Krt8-labelled group loses
AT2 RNA detection under the original depth-budget rules (C118). These deposits
share a laboratory, so their agreement is consistency rather than independent
confirmation. The M3 background-centered distal-accessibility contrasts are
approximately -37.8%, -18.3% and -15.5% of reference accessibility; raw contrasts
are -31.8%, -16.0% and -8.5%. These are different estimands (C131).
The uninjured comparator is a **7-week Cebpa mutant**, not healthy wild type.
Its response can reflect genuine genetic plasticity. No well passes all the
original gates, and reference-cell split intervals do not provide animal-level
uncertainty (C133). Chromatin closure, cellular arrest and reversibility remain
unestablished here.

The [unified specificity project](Thesis/epithelial_state_specificity/README.md)
freezes source definitions, separates full signatures from short marker panels,
compares age and genotype explicitly, and records eligibility for independent
confirmation. One well per age/genotype condition supports descriptive
contrasts; thousands of cells do not supply missing biological replication.

GSE309751 bulk ATAC peak calls could test time-associated peak detection at
frozen loci in different mice. They cannot directly test reopening in the same
cells, and binary peak absence is sensitive to sequencing and peak calling.
Quantitative differential accessibility needs suitable quantitative data;
cellular reversibility additionally needs fate information.

<!-- rq-figure:A1 -->
![A1: the AT2 identity programme in RNA and in chromatin](analysis/figures/rq/rq_a1_chromatin.png)

*Figure A1. Wildtype nuclei of GSE310539 (PBS n = 7,340, SeV n = 8,093), RNA-only embedding for display. (a) by well; (b) nuclei with both Cldn4 and Krt8 detected at the depth available (0.8% of PBS, 7.2% of SeV nuclei); (c) the AT2 identity score in RNA; (d) accessibility of the 10 promoter peaks the vendor annotated to the same nine genes, per 10,000 ATAC counts, a per-nucleus view that depth dominates; (e) the RNA score by group, medians as bars (1.12 in PBS reference, 0.31 in SeV transitional); (f) the form the registered statistic takes: for each of the 9 genes with an annotated promoter peak, the fraction of nuclei in which the transcript, and separately the promoter peak, is detected once every nucleus is held to one depth budget (RNA 7,013 UMI, chromatin 3,132 fragments; the 20th percentile of the transitional group, trial M1e's rule). Averaged over the genes, RNA detection is 73% in PBS reference and 57% in SeV transitional nuclei; promoter detection 7% and 6%. (g to i) the per-nucleus reading the heatmap replaced, kept beside it because it shows how depth inverts the answer (row C127): by group, the promoter score is highest in transitional nuclei (median 0.91, against 0.00 in PBS and 0.38 in SeV reference), because those nuclei carry about twice the fragments (median 9,435 against 4,799 and 4,842) and so fewer of them have no promoter fragment at all (39.1% against 51.1% and 50.0%); among nuclei with any signal the transitional group is the lowest (1.31 against 1.54 and 1.60). Read on its own, panel g would support the retracted "silenced but not closed" reading (C120). A visual aid for rows C118 and C121: the registered test adds the matched-gene-set null and the per-well budgets of trials M1e and M2, and these panels do not replace it. Drawn by `analysis/scripts/16_research_question_figures.py`; numbers in `analysis/figures/rq/rq_a1_groups.csv` and `rq_a1_detection_at_budget.csv`.*
<!-- /rq-figure:A1 -->

### A2. Which cells express AREG, and how sensitive are candidate rankings to the resource?

Source expression and database sensitivity are separate questions. C37's
historical permissive marker gate contained many deposited T/NK cells.
Reclassification and corrected comparisons therefore use deposited cell types,
per-donor sample units and explicit molecule-depth sensitivity. The original
0.336 versus 0.215 contrast is a historical gate result, not a clean
annotation-defined epithelial-versus-myeloid estimate. With deposited labels,
the unadjusted direction is epithelial-higher in all 10 donors (0.4591 versus
0.2715, p = 0.001953). At the prespecified primary 1,000-UMI expected-detection
budget only 6/10 retain that direction (0.1224 versus 0.0999, p = 0.130859).
The 500/2,000-UMI sensitivities give p = 0.322266/0.019531. Thus the source
contrast is measurement-dependent; the larger-budget result is not selected
as the answer. All 10 donors and all 16,064 selected cells remain at 500 and
1,000 UMIs; only the 2,000-UMI sensitivity excludes cells (1,188).

C12/C14 originally admitted pericytes and smooth-muscle senders and applied a
stromal rather than fibroblast target floor. The [corrected ligand analysis](analysis/corrections/ligand/README.md)
uses an explicit epithelial allowlist, the actual target population, saved
source/target provenance and a shared eligible donor set. Historical ranking
numbers must not be treated as corrected epithelial results. In the corrected
22-donor CellChatDB pass, AREG's donor-median rank is 10.5 among 312 retained
pairs, and it remains first among the prespecified canonical exact-EGFR ligands.
Consensus contains AREG–EGFR_ERBB2 but it is scored in only 7/22 donors, below
the 11-donor retention rule. Resource encoding and sample coverage both matter.

Expression-derived ranks nominate hypotheses. They do not establish secretion,
ligand processing, receptor activation or functional necessity. Agreement among
resources is not biological replication because the expression data, score and
underlying curation overlap. Coverage and ranking effects should be distinguished.

E6's pooled correlation (rho 0.348, p 0.112, 22 donors) did not meet its
significance criterion. The approximately 0.43 calculation was a significance
threshold, not a bound on undetected coupling or an 80%-power calculation.
The state-resolved question has fewer eligible donors. Future spatial work
should test proximity and activation in replicated samples; proximity alone
still does not establish signalling. Literature novelty requires a dated search.

### A3. Which macrophage programmes vary with phase, and what explains the differences?

The time course and per-animal fractions remain valuable descriptions of
population composition. Relative fractions do not demonstrate absolute
expansion, and the same state at a later time need not consist of the same
persisting cells. Broad macrophage pseudobulks mix subtype proportions with
within-subtype changes.

G1's resolution-versus-long-term ranking is confounded by age, infection round,
sex and genotype as well as composition (C158). W1 excludes some genotype
variation and narrows to deposited aMACs, but its late contrast retains age and
processing limitations. An adjustment cannot identify an injury-time effect
when age-matched uninjured controls are absent.

The [statistical correction work](analysis/corrections/statistics/README.md)
separates reference-method verification, sample-level sensitivity, leave-one-out
stability and subtype contributions. W1's original CAMERA-style approximation
must be distinguished from reference limma CAMERA. Its seven-gene ornithine
set was not tested under its minimum-size rule: this is an unfilled part of
the ARG1/ornithine circuit question, not evidence against that circuit. Official
reference CAMERA preserves the W1 no-hit result. G1's DNA-replication lead does
not pass confound-aware sensitivities. G2 supports 30 frozen candidates across
cohorts with fixed correlation 0.01, but none with estimated correlation; C161
is exploratory and method-sensitive. Neither
metabolite flux nor receiver-side functional response was measured.

<!-- rq-figure:A3 -->
![A3: myeloid and capillary states by phase](analysis/figures/rq/rq_a3_persistence.png)

*Figure A3. GSE262927 annotated cohort. (a to d) the myeloid embedding of trial 11 (9,997 cells, tracked coordinates) by phase, with alveolar macrophages, interstitial macrophages and inflammatory monocytes coloured and every other label in grey; (e to h) the capillary endothelium (43,359 cells, script 06 recipe) by phase, coloured by the injury-induced capillary score; (i) the iCAP fraction per animal with the median per day: 2.0% at baseline, 37.5% at 25 dpi, 21.7% at 366 dpi. A visual aid for rows C3, C12, C13 and C15; the per-animal numbers are the registered ones. Numbers in `rq_a3_myeloid_by_phase.csv` and `rq_a3_icap_by_day.csv`.*
<!-- /rq-figure:A3 -->

### A4. How do current Wnt activity and IL-1 responsiveness overlap in AT2 cells?

Historical lineage marking, current transcript detection, pathway activity and
future fate are different measurements. Wnt-associated maintenance and
IL-1-associated transition can be sequential or context dependent within the
same cells; the question does not require two stable opposing subsets.

The bulk qPCR comparison of sorted cells supports an association, not a
same-cell overlap fraction. Detection near 4–7% applies to GSE145031; the
GSE310539 figure below reports Il1r1 detection of 19.3–29.8%. Neither makes a
single sparse transcript a validated substitute for a lineage reporter.

A useful experiment would specify present activity versus lineage history,
reporter washout, animal-level four-quadrant proportions and a functional
response endpoint. Available counts and accessibility can guide feasibility;
they do not establish reporter overlap or responsiveness. Claims of being the
first or only such comparison require a scoped literature review.

<!-- rq-figure:A4 -->
![A4: Axin2 and Il1r1 in AT2 nuclei](analysis/figures/rq/rq_a4_axin2_il1r1.png)

*Figure A4. GSE310539. (a, b) Axin2 and Il1r1 on the wildtype embedding of Figure A1; (c) both transcripts in AT2 nuclei (Sftpc detected, not transitional) of all four wells, with the fraction detected above each violin (Axin2 4.4 to 7.5%, Il1r1 19.3 to 29.8%); (d) co-detection tiles per well: both detected in 1.0 to 3.1% of AT2 nuclei. A visual aid for rows C134 to C137 and C142: the question is whether the two mark distinct subsets, and at this detection depth the count matrices cannot say. Numbers in `rq_a4_detection.csv` and `rq_a4_codetection.csv`.*
<!-- /rq-figure:A4 -->

### A5. Which transitional signatures are specific to injury rather than development or genotype?

This is the specificity arm of A1, not a separate discovery project. At the
original common RNA budget, the two-transcript label identifies 3.69% of P9
control cells and 8.07% of P9 **Cebpa-mutant** cells. The control result already
refutes injury exclusivity of the classifier (C119); the larger mutant value
must not be described as ordinary development.

Shared Krt8/Cldn4 expression does not establish that full DATP, PATS and ADI
programmes, regulatory mechanisms or fates are equivalent. Those definitions
must retain their original source and gene universe. Analyses should include
label-free versions excluding Krt8 and Cldn4 to expose circular enrichment.
A maturation score does not hold development fixed: the existing age-by-injury
design is not a replicated factorial experiment.

The unified project records which complete source signatures are available,
which comparisons use marker panels, and whether an independent cohort has
adequate animals and epithelial coverage. A failed eligibility check is a
reason to narrow the conclusion, not to substitute cell-level significance.

<!-- rq-figure:A5 -->
![A5: the transitional marker set in development and after injury](analysis/figures/rq/rq_a5_development.png)

*Figure A5. GSE247130 control wells, one RNA-only embedding each: P9 (n = 12,186), seven weeks (n = 7,589) and SeV infected (n = 11,773). (a to c) nuclei with both Cldn4 and Krt8 detected at the depth available (12.47%, 2.08%, 1.27%), then Cldn4 and Krt8 on the same embeddings; (d, e) the two transcripts across wells. A visual aid for row C119; the registered comparison is trial M1c at one depth budget, where the P9 wells labelled more than any injured adult well. Numbers in `rq_a5_wells.csv`.*
<!-- /rq-figure:A5 -->

## Corrected epithelial specificity result

![Full source-defined epithelial panels and external animal coverage](Thesis/epithelial_state_specificity/results/es1_specificity.png)

*ES1 uses 2,000 UMI and label-excluded source definitions. ADI enrichment is
+1.04/+1.21 detection points in neonatal controls and +7.61/+6.53 in injured
adult controls. Only one of 25 external-study animals passes both group floors.
Seven-week control eligibility changes with seed. These are within-well
measurements and coverage checks, not replicated fate inference (C165–C168).*

## Additional RQs emerging from the corrected analyses

These are post hoc candidate questions. They are linked to observed limitations
or patterns, not presented as new established mechanisms or automatic next runs.

### A6. How much of the IPF macrophage proliferation signal is composition, and what remains within a shared noncycling state?

**Motivation:** corrected G2 support depends strongly on correlation assumptions.
The validation cohort has more deposited proliferating macrophages on average in
IPF (3.77% versus 2.11%), but its proliferating subtype has only three IPF donors
and one control above the 50-cell floor. Within-subtype findings in discovery
are not confirmed in the available validation labels.

**Test:** harmonize resident, recruited and cycling states across cohorts before
constructing donor pseudobulks. Decompose cell-fraction and raw transcript
contributions, then compare the same noncycling state under a fixed composition
standard. Require at least three donors per arm and justify power beyond that
minimum. A composition-only explanation is challenged by a reproducible
within-state effect after standardization. Failure to reach significance does
not establish composition-only; excluding a meaningful within-state effect
requires a prespecified equivalence margin and adequate precision.

Evidence: [statistical corrections](analysis/corrections/statistics/README.md).
The current data leave both mixture and within-state explanations plausible.

### A7. Does Cebpa genotype shift the reference AT2 population and compress the apparent transitional contrast?

**Motivation:** ES1's published AT2 holdout contrast is smaller in mutant wells:
P9 control -4.59/-4.79 points versus mutant -0.70/-0.31, and injured adult
control -8.88/-7.34 versus mutant -2.28/-2.43 across technical seeds. Reference
scores also decrease in mutants (C167). A smaller difference need not mean a
more preserved transitional state; it can arise from a changed reference.

**Test:** estimate a genotype-by-state interaction from independent animals at
matched ages, reporting changes in both reference and labelled states. Repeat
with state definitions independent of Sftpc and the scored genes, because
reference selection can itself change with genotype. Comparable genotype shifts
in both states, bounded by a prespecified interaction margin, would argue against
the differential-state explanation. The current one-well comparisons cannot
establish that interaction.

### A8. Does a broad AT1 score capture shared transition programmes rather than late maturation?

**Motivation:** the published ADI and AT1 holdout lists share **119 genes**.
Broad AT1 scores rise in several injured labelled groups while a four-gene
late-AT1 panel is seed-sensitive (C168). Neither the broad score nor the small
late-marker panel alone establishes fate.

**Test:** freeze shared and AT1-unique components before validation; relate both
to mature AT1 protein, morphology or lineage contribution in independent animals.
If the unique component predicts mature endpoints with adequate precision
beyond the shared component, a purely shared-transition explanation is weakened.
If no precise difference can be estimated, retain the question rather than
calling one component biologically absent.

Evidence for A7/A8: [ES1 source definitions and results](Thesis/epithelial_state_specificity/README.md).
Neither requires treating the neonatal label as the same biological state as an
adult injury intermediate.

### A9. Does apparent EGFR ligand specificity reflect receiver biology or receptor representation and coverage?

**Motivation:** C114's resource-absence claim is refuted: consensus contains
AREG–EGFR_ERBB2, but only 7 of 22 donors pass its expression score requirements,
below the 11-donor retention rule. Meanwhile, AREG leads a restricted canonical
exact-EGFR list in four resources. These are different estimands on shared data,
not independent confirmations of one mechanism.

**Test:** freeze the receptor-complex definitions, compare donor-level EGFR and
ERBB2 coexpression within the same fibroblast states, and distinguish missing
resource entries from failed expression eligibility. Require matched biological
replicates and independent protein or receptor-activation measurements before
inferring receiver specificity. If the ranking difference disappears under
matched definitions and donor coverage, the representation explanation gains
support; persistent RNA differences alone still do not prove activation.

This is a post hoc question generated by the correction, not an established
answer. Evidence: [ligand correction](analysis/corrections/ligand/RESULTS.md).

### A10. Do epithelial perturbation responses predict organoid growth and fibroblast responses across independent preparations?

**Motivation:** the next-dataset metadata gate identifies GSE307112 as a source
of species-separated epithelial/fibroblast RNA and well-linked imaging outcomes.
This supplies a possible outcome link absent from the current atlas comparisons.
It is a new prospective question, not a result of the existing reanalyses.

**Test:** freeze epithelial programmes and an imaging endpoint, resolve well,
guide, plate and biological-preparation identities, and test prediction on held-out
biological batches. Compare programme predictions with models using plate and
baseline imaging alone. A precise absence of incremental predictive performance
would weaken the proposed RNA–growth association. If independent batches cannot
be identified, restrict the analysis to descriptive screen associations.
Species-separated bulk RNA does not resolve within-compartment composition;
organoid area is not mature AT1 fate or in vivo repair.

See the [dataset gate and pilot contract](docs/NEXT_DATASET_GATE.md). This takes
priority over another unrestricted ligand-ranking pass because it can test a
different endpoint; the design still has to pass before inference.

## Part B. What a reader can reproduce

- Initial mouse and human atlas recovery, annotation disagreements and integration
  choices: [FINDINGS.md](FINDINGS.md) and [pipeline record](docs/PIPELINE_AS_RUN.md).
- Current observations, corrections and source-specific limitations:
  [CLAIMS.md](CLAIMS.md), [generated evidence index](analysis/claims/manifest.json),
  and [negative results](NEGATIVE_RESULTS.md).
- Data preparation, environment recovery and correction commands:
  [REPRODUCIBILITY.md](REPRODUCIBILITY.md).

## Part C. How evidence is assessed

A donor or animal is the inferential unit where replication exists. Cell
resampling measures technical or within-sample sensitivity, not biological
replication. Gene-set permutations and donor-label permutations test different
nulls. Discovery selection and confirmation are separated; overlapping gene
sets are not counted as independent mechanisms. Leave-one-out direction checks
measure influence, not external replication.

The corrected work preserves historical trials and writes separate results.
Machine checks bind selected numbers to explicit files and filters; unbound
claims remain visibly outside that check's coverage. Status counts combine
observations, methods and decisions and are not a measure of scientific merit.

## Part D. Next decisions

1. Consolidate corrected source, enrichment and epithelial-specificity results
   before promoting a mechanistic claim.
2. Select confirmation datasets using independent samples, explicit outcomes,
   cell-state coverage, age/genotype controls and quantitative assay availability.
3. Start the [outcome-linked dataset pilot](docs/NEXT_DATASET_GATE.md), selected
   after reviewing the supplied project-purpose document. Its sample identities
   and replication must pass before an inferential model is fitted.
4. Undertake spatial or functional follow-up only for surviving candidates with
   a defined discriminating endpoint.

## Part E. Limits stated first

Several deposits contain one pooled library per condition. Same-laboratory
agreement, alternate databases on the same expression matrix, and repeated
resampling do not create independent cohorts. There is no common functional
repair outcome across this collection. The project supplies reproducible
observations and candidate mechanisms; causal claims require additional design.
