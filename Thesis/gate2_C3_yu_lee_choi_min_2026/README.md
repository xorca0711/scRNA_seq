# Yu, Lee, Choi_Min and Choi 2026: IL-1beta and stem-cell plasticity

**IL-1beta signaling as a molecular arbiter of stem cell plasticity:
orchestrating the niches of repair, fibrosis, and cancer.** Sua Yu, Seo Hyeon
Lee, Min Seo Choi and Jinwook Choi. *Seminars in Immunology* 83, 102050.
[DOI](https://doi.org/10.1016/j.smim.2026.102050).

This is branch **2C, item 3**, corresponding to stable roadmap paper **13**.
The folder name follows the owner's instruction of 24 September 2026; the
stable paper number remains 13 in the existing roadmap. The owner has read
the review. The synthesis and proposed analyses below are assistant-authored
and remain open to revision.

**Execution status: staged analysis in progress; first results available.**
Read the [initial analysis report](INITIAL_RUN_REPORT.md) for completed mouse
QC/clustering and two-cohort IPF pathway results. No pathway passes the primary
global FDR threshold; fixed-correlation sensitivity results are not substituted
for the primary analysis. The initial donor-level ligand-receptor run completed
for all 25 donors with both resources; output-integrity checks passed. Differential
contrasts, sensitivity analyses and biological interpretation remain pending.
See [the review and launch limits](FINAL_REVIEW.md),
[mouse QC](trials/u3_acquire_qc/run_record.json),
[cluster diagnostics](trials/u3_cluster_review/run_record.json), and
[IPF pathway run](trials/u5_ipf_pathways/run_record.json), and
[ligand-receptor run](trials/u5_ipf_liana/run_record.json).
Ten GEO metadata records and two count-file format checks were completed
during preparation, before the owner requested a planning checkpoint. No
biological effect was computed during that preparation. The current launch
adds mouse QC/annotation diagnostics and donor-level IPF pathway analysis;
KAC and niche inference have their own scientific release gates.

## Read and review

| Artifact | Purpose |
|---|---|
| [Analysis pipeline plan](ANALYSIS_TRIAL_PLAN.md) | Research questions, estimands, stages, controls, decision gates and outputs |
| [Context-specific niche analyses](NICHE_ANALYSIS_PLAN.md) | Required macrophage/fibroblast ligand-receptor and pathway analyses; Body-note context distinctions and conditional ligand-target extension |
| [Public-data shortlist](DATASETS.md) | Ranked new deposits, existing-data roles, unresolved access and design issues |
| [Source and annotation synthesis](SOURCE_SYNTHESIS.md) | What the annotations emphasize, what primary evidence supports, and what remains a hypothesis |
| [Machine-readable dataset manifest](dataset_candidates.json) | Source links, roles, dependencies and release gates |
| [Proposed analysis contract](analysis_contract.json) | Review status, primary estimand, proposed floors and execution boundary |
| [Niche analysis contract](niche_analysis_contract.json) | Context axes, LR/pathway/target methods, eligibility and multiplicity families |
| [Preparation report](trials/PREPARATION_REPORT.md) | What actually ran and what it established |
| [Initial analysis report](INITIAL_RUN_REPORT.md) | Completed modules, null/sensitivity results, remaining gates and execution evidence |
| [Figure gallery plan](FIGURE_GALLERY_PLAN.md) | Seven figure groups, supporting diagnostics, panel specifications and release criteria |
| [Sample inventory](trials/u0_geo_design_audit/sample_inventory.csv) | All 371 deposited sample records; these are not 371 independent animals or donors |

## Proposed research focus

> Which macrophage-fibroblast-epithelial circuits and recipient programmes
> distinguish supportive from persistent pathological niches, and how does
> their response to IL-1beta blockade relate to epithelial state balance?

Version 3 retains ligand-receptor inference and macrophage/fibroblast pathway
enrichment required analysis arms. It follows the Body notes' distinctions
between source/recipient cells, repair routes, niche mechanics, age/history,
fibrosis and distinct cancer contexts. RNA can test selected associations;
stiffness, secretion, cellular origin and functional suppression require
additional measurements. Cross-organ examples remain labelled hypotheses.

The review supplies the organizing hypothesis. Primary studies supply the
experiments, signatures and data. Persistent injury, dysplastic epithelial
repair, fibrosis and malignancy remain distinct outcomes. KRT8 expression
alone cannot assign any of them.

The first proposed biological dataset is **GSE300288**, which includes direct
anti-IL-1beta treatment. The strongest human extension is **GSE308103**, with
candidate within-patient normal/precursor/LUAD contrasts, followed by its
spatial companion **GSE307534**. These are complementary assays from the same
study and partly the same patients, not independent replication cohorts.
The post-injury spatial series **GSE267226/GSE267228** supplies a different
context, but its deposited mouse intervention is anti-CD8, not anti-IL-1beta.
See [the shortlist](DATASETS.md) for the evidence and limitations.

## Relationship to the existing repository

Existing Choi, Cardoso and England trials remain in their source-paper folders.
This branch references their outputs rather than moving or rerunning them.
The [epithelial specificity project](../epithelial_state_specificity/README.md)
provides traceable state definitions and known limitations. The
[outcome-linked dataset gate](../../docs/NEXT_DATASET_GATE.md) remains an
independent route for organoid growth outcomes; it is not replaced by this plan.

The review does not supply a new experimental dataset or a result to reproduce
as a single analysis. This branch evaluates its proposed links across public
primary studies, including situations that could contradict the unified model.

## Figure gallery

The [gallery plan](FIGURE_GALLERY_PLAN.md) orders figures by the research
questions. Only generated, checked figures are embedded below; planned
figures remain listed with their analysis dependencies.

| Figure group | Current gallery status |
|---|---|
| F01 Context, design and cell coverage | Inputs available for initial batch; rendering planned |
| F02 Sources, recipients and inhibitory context | Donor expression panels pending; full-source scan required |
| F03 Directional ligand-receptor results | Initial 25-donor descriptive run and output checks complete; rendering, contrasts and sensitivities pending |
| F04 Recipient pathways | Initial two-cohort IPF figure available below |
| F05 Mouse blockade and epithelial state balance | Annotation and coverage gates remain |
| F06 Human lesion and spatial context | Patient/tissue crosswalk and analysis pending |
| F07 Cross-context evidence synthesis | Planned after component analyses |

### F04: initial IPF pathway results

![IPF pathway primary and sensitivity results](figures/ipf_pathway_primary_and_sensitivity.png)

Primary estimated-correlation and fixed-correlation sensitivity results for
the declared broad-compartment pathways. Filled points indicate global
q < 0.05 within the corresponding cohort and analysis family. See the
[initial report](INITIAL_RUN_REPORT.md) for eligibility and interpretation;
neither cohort has a primary pathway passing this threshold.

## Layout and execution boundary

`trials/` contains specifications, scripts and compact preparation outputs.
`cache/` is ignored and holds downloaded metadata, count examples and the
private annotation snapshot. Copyrighted papers, full private notes, large
matrices and personal planning are not tracked. Raw inputs elsewhere in the
repository were not changed.

U0/U1 are preparation utilities. U2/U3 and the initial IPF pathway arm now
have executable stages and run records. Later biological stages are released
only when their design, annotation and resource specifications are satisfied;
launch authorization is recorded in the analysis contracts.
