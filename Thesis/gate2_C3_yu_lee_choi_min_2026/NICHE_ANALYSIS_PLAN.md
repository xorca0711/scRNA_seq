# Context-specific macrophage, fibroblast and epithelial signalling

**Version 3, 24 September 2026; owner-authorized staged execution.**
The [final review](FINAL_REVIEW.md) records fixes, source-access limits and
scientific release gates. The initial mouse QC and IPF pathway stages are released.
This is a required part of the [main analysis plan](ANALYSIS_TRIAL_PLAN.md).
It operationalizes the distinctions in the owner's [Body annotations](https://app.notion.com/p/3e3151616b44806881baf2d8ffeae89e),
re-fetched on 24 September (page last edited 06:54:04 UTC), and the review's
primary-source map in [SOURCE_SYNTHESIS.md](SOURCE_SYNTHESIS.md). The earlier
plan named niche programmes but did not specify ligand-receptor inference or
macrophage/fibroblast pathway enrichment. This revision fills that omission.
Everything below is a proposed test, not an observed result.

## 1. Context determines the question

The organizing question is: **Which source-to-recipient circuits and recipient
programmes accompany supportive, suppressive or persistent niche states, and
which of these change under IL-1beta blockade?** Epithelial state balance is
the phenotypic anchor. Macrophage and fibroblast responses are central aims.

| Level in the Body notes | Explicit distinction | Measurement and limit |
|---|---|---|
| Molecular regulation | Ligand transcription, maturation/release, receptor competence and endogenous inhibition | Keep IL1B, processing genes, IL1R1/IL1RAP and IL1RN/IL1R2/SIGIRR separate; RNA cannot establish secretion or inflammasome activation |
| Source and recipient cells | Macrophages/monocytes/neutrophils versus epithelial sources; epithelial versus fibroblast recipients | Directional, subtype-resolved ligand-receptor analysis plus recipient-specific enrichment; cell abundance and within-state changes reported separately |
| Regenerative route | Resident AT2 transition versus airway-derived repair | Separate alveolar and airway annotations; Notch/FOSL2 is an airway extension, not a DATP definition |
| Niche condition | Transient stromal remodelling versus inflammatory/matrix persistence; direct epithelial stimulation versus fibroblast-mediated suppression | Fibroblast inflammatory, matrix and trophic programmes tested separately; organoid suppression and mechanical stiffness require functional measurements |
| History and host | Signal onset/withdrawal, duration, age, senescence and genotype | Use verified metadata; neither a late time point nor a SASP-like RNA panel establishes chronic exposure, senescence or memory |
| Tissue and disease | Non-oncogenic repair, established fibrosis, precursor lesions, established malignancy | Fit contrasts within each study/context first; no single fibrosis-to-cancer trajectory |
| Spatial and tissue outcomes | Molecular compatibility, neighbourhood organization, organ architecture/function | RNA edges, regional co-localization and measured function are distinct evidence; Visium cannot establish cell-to-cell contact |

The notes' cross-organ examples inform hypotheses without making those
mechanisms established in lung. In particular, CRC inflammatory CAF chemokines,
PDAC stellate-cell circuits and breast-cancer myeloid cascades remain labelled
by their source context. Endothelial, neutrophil and T-cell branches below
are conditional extensions when those cells are represented.

## 2. Research questions and competing outcomes

| ID | Question | Planned test | Discriminating interpretation |
|---|---|---|---|
| N1: source and recipient | Which myeloid states supply candidate IL-1 signals, and do fibroblasts and epithelial states show different recipient responses? | IL1B/IL1A-to-IL1R1/IL1RAP compatibility by fixed source/target subtype; recipient NF-kB/inflammatory enrichment; abundance versus within-state decomposition | More macrophages alone does not establish stronger within-macrophage output; receptor expression without a recipient response weakens a transcriptional circuit interpretation |
| N2: fibroblast context | Does a fibroblast niche retain trophic support, adopt inflammatory recruitment, or show matrix remodelling, and how does that relate to epithelial plasticity? | Fibroblast pathway tests and outgoing ligand families, separately in injury, IPF and precursor/cancer cohorts | Inflammatory and matrix programmes may diverge; their co-occurrence does not prove that inflammation causes fibrosis or that fibroblasts suppress progenitors |
| N3: reciprocal niche | Which macrophage-to-fibroblast, fibroblast-to-myeloid and fibroblast-to-epithelial edges accompany persistent states? | Fixed directional ligand-receptor tables, recipient target prioritization and available spatial context | A one-sided programme change or an edge driven by cell mixture weakens the proposed reciprocal circuit; a feedback loop is not proved by correlated RNA |
| N4: perturbation | Does IL-1beta blockade change macrophage/fibroblast recipient programmes and epithelial state balance together or separately? | Animal-level anti-IL-1beta versus IgG at each endpoint; parallel niche and KAC contrasts, optional four-arm interaction | Recipient inflammation may decrease despite unchanged IL1B RNA; reduced KAC fraction alone does not show restored supportive stroma or mature fate |
| N5: context specificity | Are these associations shared with repair or specific to fibrosis, precursor lesions or established cancer? | Within-study estimates and subtype/pathway correspondence across eligible cohorts | Shared generic stress or different captured cell types may explain agreement; age, genotype and assay remain possible effect modifiers |

**Neutralizing IL-1beta protein need not lower IL1B RNA or an RNA-based
ligand-receptor score.** An unchanged score cannot refute successful blockade.
Recipient transcriptional changes and matched protein/functional evidence,
where available, determine the strength of the treatment interpretation.
Conversely, NF-kB and STAT3 responses are not specific to IL-1beta.

## 3. Dataset roles and eligibility

| Context | Candidate data | Required compartment/design check |
|---|---|---|
| Direct blockade in a cancer-prone lung model | GSE300288 | Verify animals, assay and treatment history; establish macrophage, fibroblast and epithelial coverage in the same biological samples before a triad analysis |
| Human precursor and LUAD context | GSE308103; GSE307534 spatial companion | Patient/lesion links, fixed-RNA coverage and subtype representation; companions are not independent replication |
| Fibrotic recipient niche | Existing GSE136831 and GSE135893 | Extend full source matrices to macrophages and fibroblasts using donor metadata; these are previously inspected cohorts, not unseen validation |
| Non-oncogenic injury | GSE141259; existing GSE262927 | Resolve animal/sort/time structure and matched compartments; early and late samples cannot silently substitute for dose or withdrawal |
| Pathological post-injury spatial context | GSE267226 and GSE267228 | Descriptive small-n regional evidence; mouse treatment is anti-CD8, not anti-IL-1beta |
| Organoid recipient-context mechanism | Ciminieri 2023; existing Choi/Cardoso evidence | Functional context only until a suitable matched, replicated expression accession is verified; epithelial-only inputs cannot supply missing stromal or myeloid cells |

The existing [ligand correction](../../analysis/corrections/ligand/README.md)
used epithelial senders and fibroblast/myofibroblast receivers. Its subset
caches omit macrophages and cannot be reused as complete triad inputs.
Reuse audited methods and metadata, then create bounded paper-local extracts
from the unchanged full sources. Do not rerun every historical analysis.

Freeze a sample-by-subtype coverage table before effect calculations. Separate
resident alveolar, interstitial/inflammatory macrophage states and monocytes
where source labels support them; do not infer lineage origin from RNA or
force M1/M2 labels. Preserve source-defined alveolar/supportive, inflammatory
and matrix/myofibroblast states with uncertain labels retained. Homologous
labels across species require marker evidence, not matching names alone.

Proposed niche floors: 50 cells of each analysed subtype per biological sample,
sensitivity 30/50/100; at least three independent animals/donors per arm,
or three complete patients for a paired contrast. This is an eligibility
floor, not a power claim. Broad-compartment analyses remain separate from
subtype tests. Pairwise edges need both compartments in the same sample;
joint triad associations need all three. A missing fibroblast state blocks
that comparison, not every epithelial analysis. No borrowing missing cells
from other donors or zero-imputing unobserved compartments. If GSE300288 lacks
eligible fibroblasts, report the treatment niche question as partly unmeasured
and use IPF/human cohorts for association only; they do not replace a blockade test.

## 4. U4-LR / U5-LR: actual ligand-receptor analysis

Freeze a versioned resource export, species mapping and explicit complex
subunits before comparisons. Proposed primary resource: LIANA consensus;
CellChatDB is a separately reported sensitivity. Audit coverage first: a
missing curated interaction is not biological absence. Do not select the
resource that yields the strongest result. The following are hypothesis
families, not a claim that every exact complex has already been mapped.

| Family | Prespecified direction | Source status and intended interpretation |
|---|---|---|
| IL1B and IL1A to IL1R1 with IL1RAP | Each eligible macrophage/monocyte subtype to fibroblast and epithelial subtypes; epithelial sources assessed separately | Core review-motivated lung question; report essential receptor subunits even if a resource encodes a simpler edge; IL1A is a specificity comparator |
| AREG/HBEGF to resource-defined EGFR complexes | Epithelial and myeloid senders to fibroblasts, kept separate | Lung/repository prior and alternative trophic/remodelling route; do not collapse EGFR and EGFR-containing heteromers |
| TGFB1 to TGFBR1/TGFBR2 | Myeloid to fibroblast states | Alternative matrix-associated circuit; precursor RNA is not active TGF-beta |
| CCL2 to CCR2; CXCL12 to CXCR4 | Fibroblast states to eligible myeloid states | Cross-organ stromal-recruitment hypotheses to test in lung; do not assign the ligand to fibroblasts in advance of source measurements |
| CXCL1/2 to CXCR2 | Observed epithelial, stromal or myeloid sources to neutrophils | Lung and cross-organ cancer-context leads; conditional on neutrophil capture; transcript identity does not prove suppressive PMN-MDSC function |
| FGF7/FGF10 to compatible FGFR2 complexes; resource-supported WNT receptor complexes | Fibroblast to epithelial states | Exploratory supportive-niche alternatives; source/resource verification required; common 3-prime RNA does not establish FGFR2b isoform specificity, and WNT ligands are not all canonical |
| VEGFA to KDR/FLT1 resource entries | Myeloid/stromal/epithelial sources to endothelium | Conditional vascular extension motivated by the notes; RNA is not angiogenesis, VEGF isoform usage or barrier function |

IL1RN, IL1R2 and SIGIRR are inhibitory-context measurements, not positive
activation edges. Processing genes and extracellular matrix programmes have
their own panels. A source-agnostic descriptive scan of all observed major
compartments guards against assuming that every IL1B molecule is macrophage-derived.
New edges from that scan remain exploratory.

The core LR family comprises the IL1B, IL1A, AREG, HBEGF and TGFB1 rows.
Chemokine, trophic WNT/FGF and vascular rows are separately labelled exploratory
or conditional families. Freeze exact source/target subtype pairs and resource
complexes in U2 before effects; do not promote a passing exploratory edge.

Execution after approval:

1. On each eligible biological sample separately, run the LIANA CellChat-like
   magnitude implementation already used in the repository, on log1p counts
   normalized to 10,000. Keep fixed source-target labels and molecular pairs.
   Proposed detection fraction is 0.10, with 0.05/0.20 sensitivity. Disable
   cell-permutation P values as evidence of animal-level replication. This
   implementation is not native R CellChat.
2. Report resource, ligand and receptor subunits, mean expression, detection,
   cell counts and within-sample magnitude/rank. Do not pool treatment arms or
   take the highest-scoring subtype per donor. Sample-specific LIANA scores
   are descriptive; their scaling is not assumed comparable across runs.
3. For a common-scale differential measurement, make raw-count pseudobulks
   per sample and subtype, TMM-normalize within each subtype/study using all
   groups together, and calculate log2 CPM with edgeR prior.count=1. For a
   fixed edge, define S = (minimum ligand-subunit logCPM in sender + minimum
   receptor-subunit logCPM in receiver) / 2. This declared RNA-compatibility
   statistic is separate from the LIANA algorithm and from signalling activity.
   Include every assay-covered fixed edge in eligible samples, including
   observed zero counts stabilized by the fixed prior. Missing assays,
   missing subunits and absent cell compartments are not numerical zeros.
4. Estimate treatment mean differences in S, or within-patient lesion
   differences. Display the individual sample values and leave-one-unit-out
   range. Use unit-level permutation tests only when the design supports
   exchangeability; paired sign flips require exchangeable/symmetric paired
   differences under the null. Unsupported designs remain descriptive.
   Inspect ligand and every receptor subunit's pseudobulk changes alongside S
   so an average cannot conceal opposite component changes.
5. Correct the prespecified differential edge family with BH across all
   evaluable core edges, fixed subtype pairs and both treatment endpoints.
   Human histology contrasts have their own declared family; each IPF cohort
   uses an IPF-versus-control family. Screen-derived edges form
   a separate exploratory family. Below-detection LIANA rows, off-panel genes,
   resource absence and insufficient cells have distinct status fields; no
   missing LIANA row is converted to zero. Sensitivities do not replace the
   primary resource, threshold or correction family.

The [LIANA+ primary methods paper](https://www.nature.com/articles/s41556-024-01469-w)
supports separating communication inference from its assumptions. The fixed
RNA-compatibility contrast above is our proposed analysis, not a published
LIANA probability or a direct estimate of cytokine delivery.

## 5. U4-PW / U5-PW: macrophage and fibroblast pathway enrichment

Use **biological-sample pseudobulks, edgeR TMM, limma voom and official
limma::camera**, consistent with the [corrected statistics protocol](../../analysis/corrections/statistics/PROTOCOL.md).
Estimate residual inter-gene correlation (`inter.gene.cor=NA`,
`allow.neg.cor=FALSE`); fixed correlation 0.01 is a labelled sensitivity.
Use raw counts and estimable covariates, patient blocking for paired human
contrasts, and separate endpoint fits for mice. Never use integrated values
or cell-level replication. Record the installed versions used; do not upgrade
the working R stack just for this plan. [Official limma documentation](https://bioconductor.org/packages/release/bioc/html/limma.html).

| Compartment | Proposed gene-set collection | Biological question and limits |
|---|---|---|
| Macrophage/monocyte subtypes | Hallmark TNFA_SIGNALING_VIA_NFKB, INFLAMMATORY_RESPONSE, IL6_JAK_STAT3_SIGNALING, GLYCOLYSIS, OXIDATIVE_PHOSPHORYLATION, HYPOXIA, INTERFERON_GAMMA_RESPONSE | Inflammatory and metabolic response versus cell-number shifts; none proves secretion, metabolic flux or immune suppression |
| Fibroblast subtypes | Hallmark TNFA_SIGNALING_VIA_NFKB, IL6_JAK_STAT3_SIGNALING, TGF_BETA_SIGNALING, WNT_BETA_CATENIN_SIGNALING; frozen Reactome extracellular-matrix organization set | Inflammatory recruitment and matrix programmes may move independently; matrix RNA is not stiffness and fibroblast EMT enrichment is not evidence of an epithelial transition |
| Alveolar epithelial subtypes | NF-kB/inflammatory, hypoxia/glycolysis, p53, cell-cycle and source-defined AT2/AT1/ADI/KAC programmes | Recipient response and identity measured separately; remove assignment/overlap genes in sensitivity |
| Conditional airway/endothelial/T-cell branches | Source-defined Notch/FOSL2, vascular response and cytotoxic/IFN response sets respectively | Retain the distinct route or recipient; these are extensions only if data coverage supports them |

Freeze exact collection releases, set IDs, genes, species conversions and
hashes in U2. Names above are proposed families, not downloaded gene lists.
Require at least 70% assay coverage and at least 10 expressed/testable members
for competitive enrichment, with 50%/80% assay-coverage sensitivity. The gene
universe is all genes passing the frozen, design-aware expression filter in
that subtype, not the entire genome. Short IL-1/inhibitor, processing and
trophic panels remain descriptive if they do not qualify for enrichment.
Acquire a source-backed YAP/TEAD or senescence set only as a labelled extension;
neither becomes a mechanical or senescence assay.

Report CAMERA direction, raw P, global BH q, tested-member count, sample n
and gene-level effect distributions. BH covers all primary pathway sets,
compartments/subtypes and both endpoints in the treatment analysis; human
histology contrasts form a separately declared family. Do not select a
correlation setting or subtype after seeing which yields significance.
Per-cell module scores are visualization aids, not substitutes for enrichment.
Each IPF cohort has its own declared IPF-versus-control family across all
tested compartments and sets; directional agreement is cross-cohort consistency
in previously inspected data. Earlier human CAMERA findings were sensitive
to the correlation setting, so they are not assumed robust support for this
new question. Injury time contrasts require a frozen estimable animal design
before any time-specific inference.

## 6. U4-LT / U5-LT: connect candidate ligands to recipient programmes

For each eligible receiver subtype, use a pinned species-matched
[NicheNet ligand-target prior](https://github.com/saeyslab/nichenetr) to rank
expressed candidate ligands against that receiver's animal/donor-level DE
targets. Keep sender-focused triad candidates and the source-agnostic ranking
separate. Background genes are all tested receiver genes. Define up/down
targets separately at BH q < 0.05; do not loosen the threshold to obtain a
network. With fewer than 10 mapped targets, report insufficient target
evidence and retain the LR and pathway results without a target claim.

This connects a proposed edge to a measured downstream response; it is not
independent validation of the same RNA data or proof of mediation. Report
target-set size/coverage and leave-one-animal/patient-out ranking stability.
No causal mediation model is promised for three versus four animals.
Spatial analysis asks whether candidate sources, receivers and responses
share independently defined regions, with patients as replicates; inferred
spot mixtures cannot establish cell contact. Functional ligand neutralization,
recipient-specific perturbation and organoid output would test the mechanism.

## 7. Outputs and completion rules

### Initial execution boundary (24 September 2026)

The first two-cohort IPF pathway analysis is complete; its resources and
statistical choices were frozen in `trials/u5_ipf_spec/specification.json`.
The first GSE136831 LIANA pass freezes donor/subtype eligibility, resources,
core/extended ligand families and a seeded cap of 500 cells per donor/subtype
in `trials/u5_ipf_liana/specification.json`. This bounds initial computation.
It is a descriptive pass, not completion of the all-cell, detection-threshold,
cell-floor, differential RNA-compatibility or source-agnostic analyses.
Resource-supported WNT and ligand-target extensions remain pending.
No cell-level P values from this pass are used as donor-level evidence.

### Required deliverables

Required niche deliverables, evaluated or explicitly marked ineligible:

- Context/sample/subtype coverage matrix, including unsupported context axes.
- Fixed directional LR table and dot matrix with source/receiver coverage,
  resource status, sample effects and multiplicity labels.
- Separate macrophage and fibroblast pathway-enrichment tables/heatmaps,
  alongside epithelial phenotypic contrasts and individual sample values.
- Ligand-to-recipient-target evidence table when the target gate passes;
  distinguish RNA compatibility, recipient response, spatial support and
  intervention evidence instead of merging them into a causal network score.
- Context comparison: repair, IPF and precursor/cancer estimates displayed
  separately; no pooled progression axis or claim of fibroblast suppression
  from expression alone.

Figures belong in this paper's gallery and use the repository palette.
The main README remains navigation. Freeze implementation specifications,
resource hashes and exclusion tables before biological comparison; record
failure/null results and preserve unchanged historical claims. The owner authorized launch after the final self-review. Missing annotations,
coverage or valid sample design still prevent the affected comparisons.
