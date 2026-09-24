# Figure gallery plan: context-specific IL-1beta niches

Version 1, 24 September 2026. This presentation plan was written after the
initial IPF pathway results were inspected. It does not change the frozen
analysis families, eligibility rules or primary statistical methods. Figure
selection follows the research questions, including null and conflicting
results, rather than a significance ranking.

The gallery belongs in this paper's [README](README.md). The repository's
main README remains navigation. Published figures must link to their plotted
tables, generating script and run record. Planned panels appear in a status
table; missing images are not embedded as placeholders.

## Reading order and current readiness

The completed gallery and evidence register are in [README](README.md#figure-gallery)
and [EVIDENCE_REVIEW](EVIDENCE_REVIEW.md). The panel specifications below
preserve the original presentation plan; unavailable endpoints are explicitly
represented by eligibility tables or narrower measured panels.

The post-analysis [derived questions and figure plan](DERIVED_RESEARCH_QUESTIONS.md)
maps four follow-up RQs and an annotation audit to proposed UMAP, pseudobulk
PCA, violin/paired-point, enrichment and spatial/perturbation panels. Its D0-D4
identifiers are proposals, distinct from the completed F01-F07 figures below.

| ID | Figure and question | Current readiness | Analysis dependency |
|---|---|---|---|
| F01 | Context, experimental design and observed cell coverage: what can each cohort answer? | Completed design and biological-unit coverage figure | U0/U2/U3; current eligibility tables |
| F02 | IL-1 sources, recipients and inhibitory context: who expresses the components? | Completed IPF and all-QC human sources; receptor/inhibitor tables retain assay limits | N1; validated labels and donor-level expression/detection tables |
| F03 | Directional macrophage-fibroblast-epithelial RNA compatibility: which fixed edges are supported? | Completed IPF, mouse and paired human RNA compatibility with native LR sensitivities | N1/N3; U4-LR/U5-LR |
| F04 | Fibroblast, macrophage and epithelial recipient pathways: do response programmes agree? | Completed IPF, mouse and paired human pathways; eligible IPF/human ligand targets | N2/N3; U5-PW; treatment panels conditional on U4-PW |
| F05 | Direct IL-1beta blockade: do epithelial state balance and niche responses change together? | Early niche measurements complete; KAC identity and primary phenotype remain unidentifiable | N4; U3/U4 |
| F06 | Human lesion and spatial context: do associations persist within patients and regions? | All human libraries and spatial/context matrices processed; independent ROI and post-viral coordinates unavailable | N5; verified GSE308103/GSE307534 crosswalk and U5 |
| F07 | Evidence across repair, fibrosis and neoplasia: what agrees, conflicts or remains unmeasured? | Completed specificity panels and evidence matrix; interpretation review ready | N5/U6; completed or explicitly ineligible analyses |

F01-F04 are the first gallery-building batch. F05-F07 are scientific goals,
not promised positive findings. A failed eligibility gate produces an explicit
coverage/ineligibility panel or table, not an invented biological result.
The seven-part outline is an upper-level organization, not a requirement to
produce seven dense composite images; split panels when labels become small.

## F01. Context, design and coverage

**Panel A — study/context map.** Separate rows for direct mouse blockade,
non-oncogenic injury, human IPF, precursor/LUAD, and spatial companions.
Columns show species, assay, intervention or observational contrast, verified
biological units, tissue sampling, and measured context axes. Mark unknown
age/exposure history rather than inferring chronicity from an endpoint.
RNA and spatial companions share a study and cannot be added as independent n.

**Panel B — sample-by-subtype coverage.** Show captured cell counts and the
50-cell niche threshold with explicit eligible/ineligible/annotation-pending
labels. Distinguish pairwise sender-receiver availability from availability
of all three compartments in a single sample. The current 25-donor LR subset
is selected for eligible pairs; it is not necessarily 25 complete triads.
Display donor counts by condition beside each tested contrast.

**Panel C — mouse acquisition/QC.** Show deposited versus retained cells for
each of the seven early libraries, with the exclusion rules. Separate QC
retention from annotation confidence. Captured cell fractions are not whole-
lung abundance because the preparation enriched epithelial cells.

Inputs: `trials/u0_geo_design_audit/sample_inventory.csv`,
`trials/u3_acquire_qc/sample_qc.csv`,
`trials/u5_ipf_liana/donor_subtype_coverage.csv`,
`trials/u5_ipf_liana/fixed_donor_pairs.csv`, and pathway eligibility tables.
Panel labels must distinguish verified animal/donor counts from library counts.

## F02. Sources, recipients and inhibitory context

**Panel A — source expression.** Summarize IL1B and IL1A separately across
observed major compartments, using donor-level expression and detection
fractions. Each donor contributes equally to group summaries; display donor
points or a companion donor table. Include non-macrophage sources when the
full-source extraction supports them. The current restricted triad cannot
establish the dominant source across the entire tissue.

**Panel B — recipient components.** Display IL1R1 and IL1RAP separately in
fibroblast and epithelial subtypes, with donor n and detection. Preserve
source labels; do not infer macrophage origin or force M1/M2 categories.

**Panel C — regulation.** Keep IL1RN, IL1R2 and SIGIRR in a distinct inhibitory-
context panel. Any processing-machinery extension uses a source-frozen gene
list and is labelled RNA expression, not cytokine maturation or secretion.

Keep cohorts/species in separate facets and distinguish cell abundance from
within-subtype expression. Assayed non-detection, off-panel genes and absent
compartments have different symbols and explanatory labels.

## F03. Directional ligand-receptor results

**Panel A — fixed-edge coverage.** List the core IL1B, IL1A, AREG, HBEGF and
TGFB1 families with sender, receiver, curated ligand/receptor subunits and
eligible donor counts. Essential IL1RAP and TGFBR subunits remain visible
even where a resource encodes a simpler edge. Report resource omissions.

**Panel B — within-donor descriptive evidence.** Show LIANA magnitude/rank
within donor, resource and fixed candidate set. Label the initial seeded
500-cell/subtype cap prominently. Do not compare pooled raw LIANA magnitudes
as a calibrated between-donor treatment effect or use cell-permutation P
values as donor evidence. Omitted LIANA rows are not zero-valued edges.

**Panel C — common-scale contrasts, when completed.** Plot the declared
pseudobulk RNA-compatibility score difference for fixed subtype pairs beside
individual donor values and ligand/receptor-component changes. Display the
prespecified adjusted P value only where inference is eligible; otherwise
show descriptive effects. Label leave-one-unit-out ranges as robustness
ranges, not confidence intervals. Keep each cohort/contrast separate.

Core and exploratory chemokine/FGF/WNT/vascular families occupy separate
panels. Neutrophil and endothelial targets need their own coverage checks;
including a ligand in the current resource subset does not evaluate every
recipient branch. Consensus versus CellChatDB comparisons and all-cell,
seed, detection and cell-floor sensitivities accompany, rather than replace,
the first result. A network diagram is optional only after this evidence
table exists; its arrows encode proposed sender/receiver orientation, not
demonstrated signalling or feedback. Edge width cannot imply cytokine flux.

## F04. Recipient pathways and robustness

The current [IPF figure](figures/ipf_pathway_primary_and_sensitivity.png)
is the initial broad-compartment panel. Preserve its filename and links.
Neither cohort has a pathway passing primary global FDR < 0.05. The fixed-
correlation findings remain sensitivity results, including disagreement
between cohorts; they are not promoted to primary discoveries.

**Panel A — fibroblasts.** Separate inflammatory/NF-kB/STAT3, matrix/TGF-beta,
and trophic WNT programmes. **Panel B — macrophages.** Show inflammatory,
interferon and metabolic/hypoxic programmes. **Panel C — epithelial recipients.**
Show the declared AT2 response sets without conflating identity with response.
Include all declared evaluable sets and mark ineligible subtypes explicitly.

**Panel D — assumptions and eligibility.** Place estimated-correlation primary
results beside fixed-0.01 sensitivity, donor n, tested-gene count, assay coverage
and eligible 30/50/100-cell-floor results. The 30-cell result is conditional
on the original cache's broad-compartment >=50-cell donor selection. It does
not recover donors excluded before caching. AT2 broad and sole-subtype views
are duplicates, not biological replication.

CAMERA supplies competitive direction and P/q values, not a GSEA normalized
enrichment score. Do not invent an NES or label a signed significance score
as an effect size. Supporting gene-level log-fold changes can be shown
separately, using the frozen gene lists rather than selecting successful genes.
Treatment pathways become separately labelled panels only after U4 eligibility.

## F05. Perturbation and epithelial state balance

**Panel A — annotation evidence.** Show source-defined state markers and
coverage, AT2/AT1/airway separation, ambiguous assignments and classification
sensitivity after removing IL-1/NF-kB response genes. KRT8-high alone is not KAC.
UMAP is an annotation diagnostic, not evidence of lineage or therapeutic rescue.

**Panel B — phenotype.** Plot every eligible animal's KAC fraction among
alveolar epithelial cells, denominator and mean contrast in percentage points.
Add the source paper's all-epithelial denominator as a labelled sensitivity.
Do not connect unpaired mice with lines. Keep early and later endpoints apart.
If defensible KAC identity is unavailable, retain an ineligibility statement;
a continuous proxy score cannot replace the primary fraction endpoint.

**Panel C — niche response.** Align eligible macrophage/fibroblast pathway
and RNA-compatibility contrasts with the phenotype, with per-arm coverage
and attrition. An unchanged IL1B transcript or compatibility score does not
refute protein neutralization. No matched histology or function panel is
added without verified animal identifiers.

The fixed two-endpoint exact-test/Holm family cannot attain adjusted
significance at the deposited small n under the reviewed setup. Keep effects,
individual animals and uncertainty visible; do not switch tests for stars.

## F06. Human lesions and spatial extension

**Panel A — patient/lesion map.** Preserve normal, AAH, AIS, MIA and LUAD
labels and repeated lesions; connect measurements only for verified patients.
**Panel B — within-patient contrasts.** Display eligible cell-state, fixed-edge
and recipient-pathway contrasts with patient-level values. **Panel C — spatial
support.** Show candidate source/recipient/response distributions in source-
defined anatomical regions, together with patient-level region summaries.

Do not count spots as patients, interpret spot mixtures as cell contact, or
present RNA and spatial companions as independent replication. Tissue images
and external annotations require usable source files and provenance; do not
substitute an illustrative tissue drawing for measured spatial data.

## F07. Context and evidence synthesis

Use an evidence matrix, not a single fibrosis-to-cancer trajectory. Rows are
fixed circuit/recipient-programme questions; columns separate injury repair,
IPF, precursor lesions, malignancy and direct blockade. Each entry states
tested context, direction/effect where comparable, eligible biological n,
and evidence type: expression compatibility, recipient response, spatial
association or intervention. Distinguish unmeasured, ineligible, inconclusive,
discordant and supported-by-the-specified-test outcomes in plain labels.

Keep source-paper functional evidence separate from this repository's
reanalysis. Cross-organ CRC/PDAC/breast mechanisms remain hypotheses in lung.
Display ligand-target support only if its mapped-target gate passes and label
it as dependent on the same expression data, not independent validation.

## Supporting figures

| ID | Content | Role |
|---|---|---|
| S01 | Per-library QC, treatment-blind clusters, marker support, uncertain labels, sample contribution and preparation limitations | Annotation audit; no treatment-coloured discovery used to define KAC |
| S02 | Resource/subunit coverage, missing-row reasons, LR cap/seed/threshold checks, all declared pathway sensitivity tables and leave-one-unit-out diagnostics | Show which conclusions survive assumptions; include failures |
| S03 | Source-defined DATP/ADI/PATS/KAC/HPCS programme overlap, assignment-gene removal, matched-gene/depth/stress controls and optional target-ranking stability | Specificity/target extensions only where their source and sample gates pass |

## Rendering, captions and release checks

- Read colours from `analysis/config/palette.json`; never hard-code a new
  palette. Within a figure, colour has one named meaning. Primary versus
  sensitivity uses position/shape or line style, so colour can retain biological
  identity. Label categorical changes between figures explicitly. Add direct
  labels/shapes and a table view; do not rely on green/amber contrast alone.
- Keep the same subtype and pathway order across matching panels. Use the
  validated sequential ramp for magnitude only. Use aligned zero-centred axes
  and explicit signs for contrasts; no ad hoc red/green heatmap.
- Export a readable PNG for the README plus a vector PDF/SVG for final plots
  where practical. Use standard scientific plotting tools. Preserve the
  existing initial PNG; a later vector export must use its source tables and
  reproduce its meaning, not change the analysis.
- Each caption states the question, cohort/assay, replicate unit, n by arm,
  displayed measure, normalization, statistical family, sensitivity status,
  missing-data encoding and principal limitation. A hypothesis schematic is
  explicitly labelled as a hypothesis rather than a result.
- Rendering scripts read saved outputs; they do not rerun clustering or
  statistical fits. Save plotted data, input hashes, script/version identity
  and run-record links. Inspect each exported figure for readable labels,
  correct legends, clipping and agreement with its table before publishing.
- Update the gallery status only when outputs exist and have been checked.
  No inference is required to publish honest null results or an eligibility
  failure. A figure's completion is distinct from completion of the study.

## Delivery sequence

1. Retain F04's available IPF result and its limitations in the gallery.
2. Assemble F01 and S01 from existing QC, design and coverage outputs.
3. The initial LR run and its output checks have passed. Assemble the
   initial F03 descriptive panels; add the donor-level F02 triad view with
   explicit source-scan scope. Add differential/sensitivity panels as those
   analyses actually finish.
4. Release F05 and F06 only when their annotation/design/access gates are
   satisfied, or document why a question remains unevaluable.
5. Build F07 and the final gallery index from the complete evidence register,
   retaining negative and contradictory results.

Authority: [analysis plan](ANALYSIS_TRIAL_PLAN.md),
[context-specific niche plan](NICHE_ANALYSIS_PLAN.md),
[final methodological review](FINAL_REVIEW.md), and
[initial results](INITIAL_RUN_REPORT.md).
