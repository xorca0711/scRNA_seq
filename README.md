# Integrative reanalysis of public lung single-cell and multiome data

[![Repository checks](https://github.com/xorca0711/scRNA_seq/actions/workflows/repository-checks.yml/badge.svg)](https://github.com/xorca0711/scRNA_seq/actions/workflows/repository-checks.yml)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/)

An ongoing research project that reanalyses public lung single-cell RNA-seq
and multiome (RNA + ATAC) datasets, applying and evaluating analysis frameworks
to characterize molecular phenotypes, cell-state programmes and data
distributions. Its purpose is to turn these observations into testable research
questions about lung injury, repair and tissue remodelling.

Paper studies live in [Research Article](Research%20Article/README.md); question-specific
plans and workflows live in [RQ_Specified](RQ_Specified/README.md). The shared
[research-question register](RESEARCH_QUESTIONS.md) remains the canonical index.

The organising biological question is:

> Which epithelial and immune-state programmes distinguish productive lung
> repair from persistent remodelling after injury?

The practical work combines atlas reconstruction, programme scoring,
sample-level comparisons and cross-study analysis. It examines which patterns
repeat across independent samples, how distributions vary with cell-state
composition and measurement quality, and which observations motivate a
discriminating experiment. Distinguishing repair from pathology requires
independently measured outcomes; molecular states alone do not establish it.

Most analyses start from deposited counts; E1 uses deposited normalized
expression. Author labels are held out of the original unsupervised atlas fits
and deliberately used in later compartment and annotation-sensitivity analyses.
The claim register links observations to their evidence and limitations.
Exploratory findings and independent confirmation are distinguished explicitly.

## Start here

| To assess | Open |
|---|---|
| Biological hypotheses, evidence and next tests | [Research questions](RESEARCH_QUESTIONS.md); [measurement contracts](docs/RQ_MEASUREMENT_CONTRACTS.md) |
| Current merged analyses and remaining work | [Current handoff](PROGRESS.md); [question-specific status](RQ_Specified/README.md) |
| Biological logic, interpretation limits and follow-up order | [Logical rationale review](docs/LOGICAL_RATIONALE_REVIEW.md) |
| Paper-specific analyses and figure galleries | [Paper roadmap](Research%20Article/README.md) |
| A short portfolio entry and four figures | [Portfolio summary](docs/PORTFOLIO_SUMMARY.md) |
| Three portfolio case studies and what was actually demonstrated | [Portfolio guide](docs/PORTFOLIO.md) |
| Every claim, decision authority and explicit numeric-check coverage | [Claim register](CLAIMS.md) and [generated summary](docs/CLAIM_SUMMARY.md) |
| Corrected analyses and remaining limitations | [Implementation record](docs/remediation/2026-09-22/IMPLEMENTATION_STATUS.md) |
| Commands, dependencies and data requirements | [Reproducibility guide](REPRODUCIBILITY.md) |
| A0: conserved transition programme pilot | [Scientific result](RQ_Specified/A0_conserved_epithelial_transition_program/reports/PILOT_V1_RESULTS.md); [source recovery](RQ_Specified/A0_conserved_epithelial_transition_program/reports/SOURCE_RECOVERY.md) |
| Completed organoid growth analysis and original eligibility gates | [A10 follow-up results](RQ_Specified/A10_organoid_growth_outcome/reports/FOLLOWUP_RESULTS.md); [historical dataset gate](docs/NEXT_DATASET_GATE.md) |
| Nabhan 2018 source reproduction and animal-level Wnt analysis | [Wnt niche analysis](Research%20Article/gate1_03_nabhan_2018/README.md) |

## What the portfolio demonstrates

- **Recovering and auditing an atlas.** The original mouse analysis recovers
  deposited cell types with median cluster purity 0.947, while retaining three
  contradicted marker-based annotations as evidence of state/type confusion.
- **Correcting a biological comparison.** Ligand-source and resource analyses
  explicitly test annotation, molecule depth, sender populations, target-cell
  floors and curation. Expression rankings remain candidate-generating evidence.
- **Testing the statistical interpretation.** Sample-level pseudobulks, reference
  gene-set methods and sensitivity analyses separate reproducible directions
  from significance, mixture effects and design confounding.
- **Defining epithelial-state specificity.** A unified injury/development/genotype
  analysis distinguishes traceable gene programmes from two-marker labels and
  keeps chromatin measurements separate from fate claims.

The current evidence and caveats live in [CLAIMS.md](CLAIMS.md). Earlier owner
retentions are preserved in [DEVELOPMENT.md](DEVELOPMENT.md); the subsequent
reassessment was explicitly delegated by the owner and is identified as such.

## Figure galleries by analysis branch

Existing figures are collected beside their analysis reports. Captions state
the experimental unit and the current interpretation; historical plots do not
override the claim register.

| Branch | Figures and evidence |
|---|---|
| [Shared question gallery](analysis/figures/rq/README.md) | Measured panels, diagnostic figures and proposed designs supporting [A1–A14](RESEARCH_QUESTIONS.md), with captions and sources |
| [A0: conserved epithelial transition programme](RQ_Specified/A0_conserved_epithelial_transition_program/figures/pilot_v1/README.md) | Biological-unit coverage, discovery and failed frozen intestinal transfer; [result and interpretation](RQ_Specified/A0_conserved_epithelial_transition_program/reports/PILOT_V1_RESULTS.md) |
| [A1: regulatory, lineage and functional state distinction](RQ_Specified/A1_transitional_epithelial_state_distinction/figures/README.md) | Histone/H3, lineage, IRE1/CD44 evidence; resolved HPCS mice and region-dependent AP-1 outcomes; [latest report](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/REGULATORY_FATE_REPORT.md) and [reference map](RQ_Specified/A1_transitional_epithelial_state_distinction/reports/ANALYSIS_REFERENCE_MAP.md); causal regulation-to-fate linkage remains unresolved |
| [Niethamer: viral injury and repair](Research%20Article/gate1_01_niethamer_2025/README.md#figure-gallery) | Animal-level tracing/cycling, myeloid composition and programme-inference sensitivity; links to phase, lineage and batch galleries |
| [Choi: epithelial states](Research%20Article/gate1_02_choi_2020/README.md#figure-gallery) | Deposited in vivo and organoid state maps, with links to specificity and chromatin analyses |
| [Choi: chromatin and transitional states](Research%20Article/gate1_02_choi_2020/datp_epigenetics/README.md#figure-gallery) | RNA/promoter displays and developmental specificity |
| [Choi: Axin2 and Il1r1](Research%20Article/gate1_02_choi_2020/axin2_il1r1/README.md#figure-gallery) | Transcript detection and co-detection, with reporter/activity limitations |
| [Nabhan: Wnt niches](Research%20Article/gate1_03_nabhan_2018/README.md#figure-gallery) | Source-expression heatmaps, threshold sensitivity, paired fibroblast and AT2 plots |
| [Sikkema: reference annotation](Research%20Article/gate1_04_sikkema_2023_hlca/README.md#figure-gallery) | HLCA label transfer, uncertainty, marker-score heatmap and annotation/QC checks |
| [Cardoso: ligand sources and niche candidates](Research%20Article/gate2_05_cardoso_2026/README.md#figure-gallery) | Source/depth sensitivity and resource/receptor coverage |
| [Yu, Lee, Choi_Min: IL-1beta niches](Research%20Article/gate2_C3_yu_lee_choi_min_2026/README.md#figure-gallery) | IPF, mouse and paired human niches; spatial context, specificity and evidence review |
| [Murthy: human distal-lung atlas](Research%20Article/ungated_murthy_2022/README.md#figure-gallery) | Donor UMAP and epithelial marker dot plot |
| [Epithelial-state specificity across studies](Research%20Article/epithelial_state_specificity/README.md) | Matched-depth published signatures and cohort eligibility |

## Datasets

The [IL-1beta review branch](Research%20Article/gate2_C3_yu_lee_choi_min_2026/README.md)
has completed its feasible public-data analyses across IPF, early mouse
blockade, human lesions, spatial context and epithelial program specificity.
Its [evidence review](Research%20Article/gate2_C3_yu_lee_choi_min_2026/EVIDENCE_REVIEW.md)
distinguishes measured findings from inconclusive and unidentifiable endpoints.
The [public-data shortlist](Research%20Article/gate2_C3_yu_lee_choi_min_2026/DATASETS.md)
records dataset roles; the table below describes the earlier analyses.

Deposits used in the existing biological analyses, with their source studies and where in
this repository it was read. Deposits are analysed under the roadmap paper
that motivated them. The two series that predate the roadmap sit beside their
source papers too: GSE262927 under paper 1, and GSE178360 under a folder for
its source paper, which is outside the reading list.

| Accession | Species and design | Source study | Read under | Unit and ceiling |
|---|---|---|---|---|
| [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927) | mouse, respiratory-virus (H1N1) injury time course, uninjured to 366 dpi, 33 samples, 162,175 cells | Niethamer et al., *Cell Stem Cell* 2025 | [`Research Article/gate1_01_niethamer_2025/GSE262927/`](Research%20Article/gate1_01_niethamer_2025/GSE262927/README.md), [Stage 1 follow-ups](Research%20Article/gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md), [HLCA trials S4, S5](Research%20Article/gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md) | animal; 2 per active-repair day, 8 at 42 dpi |
| [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360) | human, healthy distal lung, 3 donors, 27,729 cells | Kadur Lakshminarasimha Murthy et al., *Nature* 2022 | [`Research Article/ungated_murthy_2022/GSE178360/`](Research%20Article/ungated_murthy_2022/GSE178360/README.md), [HLCA trials S1 to S3](Research%20Article/gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md) | donor; 3 |
| [GSE145031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145031) | mouse, AT2 lineage-traced epithelium, PBS and bleomycin day 14 and 28, 6 libraries | Choi et al., *Cell Stem Cell* 2020 | [trials D0 to D7](Research%20Article/gate1_02_choi_2020/ANALYSIS_TRIAL_PLAN.md) | one library per condition; six matrices are raw barcode whitelists |
| [GSE144468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144468) | mouse, AT2 organoids with and without IL-1beta, 2 libraries | Choi et al., *Cell Stem Cell* 2020 | [trials D5, D5b](Research%20Article/gate1_02_choi_2020/ANALYSIS_TRIAL_PLAN.md) | one library per arm |
| [GSE316241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316241), [GSE316243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316243), [GSE316244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316244) | mouse, Confetti against Red2Kras mesenchyme and niche, and the Areg-flox arm; 8 libraries, 3 mice pooled each | Cardoso, Lee et al., *Nature* 2026 | [trials C0 to C2b, C5 to C11, C13](Research%20Article/gate2_05_cardoso_2026/trials/README.md) | one library per genotype; nothing between genotypes is testable |
| [GSE310335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310335) | human, KRAS G12D alveolar organoids, 2 libraries | Cardoso, Lee et al., *Nature* 2026 | [trial C0](Research%20Article/gate2_05_cardoso_2026/trials/README.md) | one library per arm |
| [GSE247505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247505) | mouse, Red2Kras clones by time and Il1r1 dosage, 20 libraries, 2 per timed or genotype group | England et al., *Cell Stem Cell* 2025 | [trial C3](Research%20Article/gate2_05_cardoso_2026/trials/README.md), [trial A2, on the GSE247504 sub-series](Research%20Article/gate1_02_choi_2020/axin2_il1r1/trials/README.md) | library; the Il1r1 deletion is not visible in the deposited counts, so genotype rests on metadata (C149) |
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | human, lung adenocarcinoma and normal lung, 11 paired donors | Kim et al., *Nature Communications* 2020 | [trials E1, E1b](Research%20Article/gate2_05_cardoso_2026/trials/README.md); [A11 held-out test](RQ_Specified/A11_lesion_programme_addition/README.md) | donor; the one tested claim in E1; 10 verified tumour-normal pairs, 8 eligible for A11 |
| [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831) | human, idiopathic pulmonary fibrosis and control, 312,928 cells | Adams et al., *Science Advances* 2020 | [trials E2, E6, C12, C14](Research%20Article/gate2_05_cardoso_2026/trials/README.md) | donor; 26 in the resource scans, 22 in the coupling test, 7 and 3 in the rare-state comparisons |
| [GSE135893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893) | human, pulmonary fibrosis and control, 114,396 cells | Habermann et al., *Science Advances* 2020 | [trial E3](Research%20Article/gate2_05_cardoso_2026/trials/README.md) | donor |
| [GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771) | mouse, bleomycin against uninjured, collagen-producing cells | Tsukui et al., *Nature Communications* 2020 | [trial E4](Research%20Article/gate2_05_cardoso_2026/trials/README.md) | animal |
| [GSE310539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310539) | mouse, 10x multiome, sorted epithelium 14 days after Sendai virus or PBS, wildtype and AP-1 mutant, 4 wells, 39,849 nuclei | Lynch et al., *Am J Respir Cell Mol Biol* 2026 | [trials M0 to M4](Research%20Article/gate1_02_choi_2020/datp_epigenetics/trials/README.md), [A1 to A1c](Research%20Article/gate1_02_choi_2020/axin2_il1r1/trials/README.md) | one well per condition, two mice pooled |
| [GSE247130](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247130) | mouse, 10x multiome, AT2 lineage, Cebpa mutant against control at P9, 7 weeks and after Sendai virus, 6 wells, 64,294 nuclei | Hassan and Chen, *Nature Communications* 2024 | [trials M0 to M4](Research%20Article/gate1_02_choi_2020/datp_epigenetics/trials/README.md), [A1 to A1c](Research%20Article/gate1_02_choi_2020/axin2_il1r1/trials/README.md) | one well per condition; **deposited suffix order inverted** |
| [GSE308103](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308103) | human, precursor lesions and lung adenocarcinoma with paired normal lung, 75 libraries | Peng et al., *Cancer Cell* 2026 | [IL-1beta review trials U5, U6](Research%20Article/gate2_C3_yu_lee_choi_min_2026/README.md); A11 discovery | patient; 23, with histology contrasts reusing patients |
| [GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259) | mouse, bleomycin epithelial time course, EpCAM-sorted, 32 deposited samples | Strunz et al., *Nature Communications* 2020 | [A5 test](RQ_Specified/A5_developmental_programme_reuse/README.md); [A0 discovery arm](RQ_Specified/A0_conserved_epithelial_transition_program/README.md) | mouse, one per sample by the published methods; 24 in the A5 primary, 9 in A0 |
| [GSE92332](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92332) | mouse, small intestinal epithelium atlas | Haber et al., *Nature* 2017 | [A0 transfer arm](RQ_Specified/A0_conserved_epithelial_transition_program/README.md) | mouse; 3 |
| [GSE307112](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307112) | mouse type 2 cells with human lung fibroblasts, alveolosphere knockout screen, 886 well libraries, 203 targets, imaging at days 7 and 14 | roadmap paper 14, unread by the owner; the deposit is used, not the paper | [A10](RQ_Specified/A10_organoid_growth_outcome/README.md) | well; 15 RNA batches; independent preparations unresolved |

The last three rows are read only by question-specific workflows under
[`RQ_Specified/`](RQ_Specified/README.md), not under a paper folder. GSE308103 and
GSE131907 are read under paper folders and also serve A11, as discovery and held-out
test respectively. A0's developmental arm uses the Sountoulidis 2023 human airway
atlas; its data source and identities are recorded in that question's
[source recovery](RQ_Specified/A0_conserved_epithelial_transition_program/reports/SOURCE_RECOVERY.md).

Assessed and not opened: [GSE144598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144598)
(the DATP paper's ATAC-seq, coverage tracks only, unusable);
[GSE309751](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309751)
(bulk ATAC-seq companion to GSE310539, 2 to 3 mice per group at 14 and 49
days, the only replicated chromatin contrast in reach). GSE150957 was never
opened and is recorded only in row C146, withdrawn as out of focus because bulk
cannot answer a co-occurrence question. Every opened deposit's citation is in
[`REFERENCES.md`](REFERENCES.md).

| | |
|---|---|
| Stack | Python 3.12, scanpy and AnnData, Scrublet, harmonypy, PAGA and diffusion pseudotime, scvi-tools and scArches for reference mapping, LIANA+ for ligand-receptor inference, h5py streaming for multiome matrices |
| Unit and statistics | The animal or donor is the unit; medians per group; no P value where a group holds two animals; sham bands and matched-gene-set nulls where cells are the only replicate |
| Provenance | Rules frozen in run records before data are opened; every summary generated from a tracked table; an audit of the register against its own artefacts is itself a register row |
| Checks | Dependency-free artefact validator, pinned environments for both the native and the emulated x86-64 interpreter, CI on every push |

## Where to go

| Question | Read |
|---|---|
| What are the questions, what do the data say, what would move them? | [`RESEARCH_QUESTIONS.md`](RESEARCH_QUESTIONS.md) |
| What has been claimed, what stands behind each claim, and what is its status? | [`CLAIMS.md`](CLAIMS.md) |
| What was refuted, retracted or could not be established? | [`NEGATIVE_RESULTS.md`](NEGATIVE_RESULTS.md) (generated) |
| The original two-series analysis, with figures | [`FINDINGS.md`](FINDINGS.md) |
| The paper roadmap: study notes, extracted criteria, per-paper trials | [`Research Article/README.md`](Research%20Article/README.md) |
| What exactly ran, with parameters? | [`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md) (generated) |
| How can I validate or reproduce it? | [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) |
| Why each analytical decision? | [`docs/ANALYSIS_RATIONALE.md`](docs/ANALYSIS_RATIONALE.md) |
| Who decided what, and what was rejected? | [`DEVELOPMENT.md`](DEVELOPMENT.md) |
| Current state, known issues, pending decisions | [`PROGRESS.md`](PROGRESS.md) |
| Machine context for AI sessions | [`AI_CONTEXT.md`](AI_CONTEXT.md) |
| Per-dataset reports, figures, QC for the two original series | [`Research Article/gate1_01_niethamer_2025/GSE262927/`](Research%20Article/gate1_01_niethamer_2025/GSE262927/README.md), [`Research Article/ungated_murthy_2022/GSE178360/`](Research%20Article/ungated_murthy_2022/GSE178360/README.md) |
| What was displaced and why | [`archive/`](archive/DISPLACED.md) |

## Repository map

The [structure and label contract](docs/REPOSITORY_STRUCTURE.md) defines
canonical question, figure, script and paper-record locations.

```
RESEARCH_QUESTIONS.md        biological hypotheses, evidence, tests and execution readiness
CLAIMS.md                    claims register: evidence, status, potential; counts in analysis/claims/
NEGATIVE_RESULTS.md          generated from the register: refuted, retracted, unestablished
FINDINGS.md                  the original two-series analysis with figures
DEVELOPMENT.md               who decided what; rejected output stays visible
PROGRESS.md                  living handoff: state, known issues, pending decisions
AI_CONTEXT.md                machine-oriented context for AI sessions
REPRODUCIBILITY.md           input layout, validation tiers, re-run guide
REFERENCES.md                every source and roadmap paper, DOIs, data accessions
analysis/
  scripts/                   the pipeline, the focused analyses, the generators, the validator
  config/                    the one validated palette; the x86-64 environment lock
  figures/                   repository-level figures: the claims ledger from the register, and rq/,
                             question-linked data figures and explicitly labelled study designs
  LAYOUT.md                  what lives where, and why the mouse cohorts differ
Research Article/                      paper roadmap, one folder per paper, in reading order; each deposit beside its paper
  gate1_01_niethamer_2025/     phase/myeloid follow-ups; Stage 2 designs and executed follow-ups
    GSE262927/                   the mouse series: report, figures, tables, QC, focused analyses
  gate1_02_choi_2020/          trials D0 to D7; two branches:
    datp_epigenetics/            the transitional state in chromatin (M0 to M4)
    axin2_il1r1/                 the paper's own closing question (A1 to A2)
  gate1_03_nabhan_2018/        source reproduction, descriptive Nb1, external eligibility
  gate1_04_sikkema_2023_hlca/  reference mapping to the HLCA (S1 to S5)
  gate2_05_cardoso_2026/       early fibrotic niches (C0 to C14, E1 to E6)
  gate2_C3_yu_lee_choi_min_2026/ IL-1 context analyses and paper-specific F01–F07 gallery
  epithelial_state_specificity/ cross-study signature, genotype and coverage analysis
  ungated_murthy_2022/         a source paper outside the roadmap; pointer note only
    GSE178360/                   the human series: report, figures, tables, QC
docs/                        rationale, background, generated pipeline record, source-study notes
archive/                     displaced material: what moved, when, and why
```

`raw_data/` (GEO downloads) and all regenerable `.h5ad` and
`.npz` objects are gitignored; figures, small tables and run records are
tracked. No sequencing data, count matrices or paper PDFs are committed.

## Reproducing

The evidence-contract checks use only Python's standard library:

```bash
python -m unittest discover -s analysis/tests -q
python analysis/scripts/claim_contract.py --check
python analysis/scripts/validate_repository.py
python -m compileall -q analysis "Research Article" RQ_Specified
```

Scientific reruns require the inputs and dependencies listed in
[REPRODUCIBILITY.md](REPRODUCIBILITY.md), including the separate correction
workflows and their saved numerical outputs. The original atlas pipeline and
historical paper-by-paper trials remain available; seeds and eligibility rules
are specified per analysis rather than assumed to be universal.

## Still open

The collection lacks a common functional repair outcome and several comparisons
lack biological replication. Age, genotype, processing and cell mixtures remain
important alternative explanations. Independent confirmation and future dataset
selection are governed by the [research questions](RESEARCH_QUESTIONS.md) and
[implementation record](docs/remediation/2026-09-22/IMPLEMENTATION_STATUS.md).
A polished figure does not close a design gap.

## Source studies

The mouse series' own published workflow (STARsolo, SoupX, scds and Scrublet,
Seurat, Slingshot, tradeSeq, in R) is documented as reference material in
[`docs/WORKFLOW_Niethamer2025.md`](docs/WORKFLOW_Niethamer2025.md) and
[`docs/scRNAseq_workflow_Niethamer2025.md`](docs/scRNAseq_workflow_Niethamer2025.md).
The analysis here is not a port of that workflow: it was built from the
deposited data alone and compared with the papers afterwards. The
tool-by-tool record of what was used and every divergence is
[`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md). Each roadmap paper's
study note records what that paper did and which of its claims this
repository could and could not re-ask.

## Licence

Code and original written material are copyright 2026 Xorca; no reuse licence
is granted. See [`LICENSE`](LICENSE). The source papers and public datasets
remain the property of their respective authors and publishers; citations and
open-access links are in [`REFERENCES.md`](REFERENCES.md).
