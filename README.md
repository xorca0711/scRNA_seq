# Integrative reanalysis of public lung single-cell and multiome data

[![Repository checks](https://github.com/xorca0711/scRNA_seq/actions/workflows/repository-checks.yml/badge.svg)](https://github.com/xorca0711/scRNA_seq/actions/workflows/repository-checks.yml)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](https://www.python.org/)

An ongoing research project that reanalyses public lung single-cell RNA-seq
and multiome (RNA + ATAC) datasets, applying and evaluating analysis frameworks
to characterize molecular phenotypes, cell-state programmes and data
distributions. Its purpose is to turn these observations into testable research
questions about lung injury, repair and tissue remodelling.

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
| Scientific questions and measurement limits | [Research questions](RESEARCH_QUESTIONS.md) |
| A short portfolio entry and four figures | [Portfolio summary](docs/PORTFOLIO_SUMMARY.md) |
| Three portfolio case studies and what was actually demonstrated | [Portfolio guide](docs/PORTFOLIO.md) |
| Every claim, decision authority and explicit numeric-check coverage | [Claim register](CLAIMS.md) and [generated summary](docs/CLAIM_SUMMARY.md) |
| Corrected analyses and remaining limitations | [Implementation record](docs/remediation/2026-09-22/IMPLEMENTATION_STATUS.md) |
| Commands, dependencies and data requirements | [Reproducibility guide](REPRODUCIBILITY.md) |
| Next outcome-linked analysis and eligibility gates | [Next dataset gate](docs/NEXT_DATASET_GATE.md) |

![Claim status by analysis family](analysis/figures/claims_ledger.png)

*Generated from the register using explicit analysis-family assignments.
Rows include biological observations, method checks and decision records;
status counts are not independent discoveries or a measure of calibration.
Current counts and review states are in the [generated summary](docs/CLAIM_SUMMARY.md).*

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

## Datasets

Every deposit the log has opened, with the study it came from and where in
this repository it was read. Deposits are analysed under the roadmap paper
that motivated them. The two series that predate the roadmap sit beside their
source papers too: GSE262927 under paper 1, and GSE178360 under a folder for
its source paper, which is outside the reading list.

| Accession | Species and design | Source study | Read under | Unit and ceiling |
|---|---|---|---|---|
| [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927) | mouse, respiratory-virus (H1N1) injury time course, uninjured to 366 dpi, 33 samples, 162,175 cells | Niethamer et al., *Cell Stem Cell* 2025 | [`Thesis/gate1_01_niethamer_2025/GSE262927/`](Thesis/gate1_01_niethamer_2025/GSE262927/README.md), [Stage 1 follow-ups](Thesis/gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md), [HLCA trials S4, S5](Thesis/gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md) | animal; 2 per active-repair day, 8 at 42 dpi |
| [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360) | human, healthy distal lung, 3 donors, 27,729 cells | Kadur Lakshminarasimha Murthy et al., *Nature* 2022 | [`Thesis/ungated_murthy_2022/GSE178360/`](Thesis/ungated_murthy_2022/GSE178360/README.md), [HLCA trials S1 to S3](Thesis/gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md) | donor; 3 |
| [GSE145031](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145031) | mouse, AT2 lineage-traced epithelium, PBS and bleomycin day 14 and 28, 6 libraries | Choi et al., *Cell Stem Cell* 2020 | [trials D0 to D7](Thesis/gate1_02_choi_2020/ANALYSIS_TRIAL_PLAN.md) | one library per condition; six matrices are raw barcode whitelists |
| [GSE144468](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144468) | mouse, AT2 organoids with and without IL-1beta, 2 libraries | Choi et al., *Cell Stem Cell* 2020 | [trials D5, D5b](Thesis/gate1_02_choi_2020/ANALYSIS_TRIAL_PLAN.md) | one library per arm |
| [GSE316241](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316241), [GSE316243](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316243), [GSE316244](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE316244) | mouse, Confetti against Red2Kras mesenchyme and niche, and the Areg-flox arm; 8 libraries, 3 mice pooled each | Cardoso, Lee et al., *Nature* 2026 | [trials C0 to C2b, C5 to C11, C13](Thesis/gate2_05_cardoso_2026/trials/README.md) | one library per genotype; nothing between genotypes is testable |
| [GSE310335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310335) | human, KRAS G12D alveolar organoids, 2 libraries | Cardoso, Lee et al., *Nature* 2026 | [trial C0](Thesis/gate2_05_cardoso_2026/trials/README.md) | one library per arm |
| [GSE247505](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247505) | mouse, Red2Kras clones by time and Il1r1 dosage, 20 libraries, 2 per timed or genotype group | England et al., *Cell Stem Cell* 2025 | [trial C3](Thesis/gate2_05_cardoso_2026/trials/README.md), [trial A2, on the GSE247504 sub-series](Thesis/gate1_02_choi_2020/axin2_il1r1/trials/README.md) | library; the Il1r1 deletion is not visible in the deposited counts, so genotype rests on metadata (C149) |
| [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907) | human, lung adenocarcinoma and normal lung, 11 paired donors | Kim et al., *Nature Communications* 2020 | [trials E1, E1b](Thesis/gate2_05_cardoso_2026/trials/README.md) | donor; the one tested claim |
| [GSE136831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136831) | human, idiopathic pulmonary fibrosis and control, 312,928 cells | Adams et al., *Science Advances* 2020 | [trials E2, E6, C12, C14](Thesis/gate2_05_cardoso_2026/trials/README.md) | donor; 26 in the resource scans, 22 in the coupling test, 7 and 3 in the rare-state comparisons |
| [GSE135893](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135893) | human, pulmonary fibrosis and control, 114,396 cells | Habermann et al., *Science Advances* 2020 | [trial E3](Thesis/gate2_05_cardoso_2026/trials/README.md) | donor |
| [GSE132771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771) | mouse, bleomycin against uninjured, collagen-producing cells | Tsukui et al., *Nature Communications* 2020 | [trial E4](Thesis/gate2_05_cardoso_2026/trials/README.md) | animal |
| [GSE310539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310539) | mouse, 10x multiome, sorted epithelium 14 days after Sendai virus or PBS, wildtype and AP-1 mutant, 4 wells, 39,849 nuclei | Lynch et al., *Am J Respir Cell Mol Biol* 2026 | [trials M0 to M4](Thesis/gate1_02_choi_2020/datp_epigenetics/trials/README.md), [A1 to A1c](Thesis/gate1_02_choi_2020/axin2_il1r1/trials/README.md) | one well per condition, two mice pooled |
| [GSE247130](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247130) | mouse, 10x multiome, AT2 lineage, Cebpa mutant against control at P9, 7 weeks and after Sendai virus, 6 wells, 64,294 nuclei | Hassan and Chen, *Nature Communications* 2024 | [trials M0 to M4](Thesis/gate1_02_choi_2020/datp_epigenetics/trials/README.md), [A1 to A1c](Thesis/gate1_02_choi_2020/axin2_il1r1/trials/README.md) | one well per condition; **deposited suffix order inverted** |

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
| The paper roadmap: study notes, extracted criteria, per-paper trials | [`Thesis/README.md`](Thesis/README.md) |
| What exactly ran, with parameters? | [`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md) (generated) |
| How can I validate or reproduce it? | [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) |
| Why each analytical decision? | [`docs/ANALYSIS_RATIONALE.md`](docs/ANALYSIS_RATIONALE.md) |
| Who decided what, and what was rejected? | [`DEVELOPMENT.md`](DEVELOPMENT.md) |
| Current state, known issues, pending decisions | [`PROGRESS.md`](PROGRESS.md) |
| Machine context for AI sessions | [`AI_CONTEXT.md`](AI_CONTEXT.md) |
| Per-dataset reports, figures, QC for the two original series | [`Thesis/gate1_01_niethamer_2025/GSE262927/`](Thesis/gate1_01_niethamer_2025/GSE262927/README.md), [`Thesis/ungated_murthy_2022/GSE178360/`](Thesis/ungated_murthy_2022/GSE178360/README.md) |
| What was displaced and why | [`archive/`](archive/DISPLACED.md) |

## Repository map

```
RESEARCH_QUESTIONS.md        the hypotheses by question; phenotypes; specified follow-up analyses
CLAIMS.md                    claims register: evidence, status, potential (164 rows)
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
                             one paper-style figure per research question from the analysed objects
  LAYOUT.md                  what lives where, and why the mouse cohorts differ
Thesis/                      paper roadmap, one folder per paper, in reading order; each deposit beside its paper
  gate1_01_niethamer_2025/     phase and myeloid follow-ups (N1 to N4); Stage 2 proposals
    GSE262927/                   the mouse series: report, figures, tables, QC, focused analyses
  gate1_02_choi_2020/          trials D0 to D7; two branches:
    datp_epigenetics/            the transitional state in chromatin (M0 to M4)
    axin2_il1r1/                 the paper's own closing question (A1 to A2)
  gate1_04_sikkema_2023_hlca/  reference mapping to the HLCA (S1 to S5)
  gate2_05_cardoso_2026/       early fibrotic niches (C0 to C14, E1 to E6)
  ungated_murthy_2022/         a source paper outside the roadmap; pointer note only
    GSE178360/                   the human series: report, figures, tables, QC
docs/                        rationale, background, generated pipeline record, source-study notes
archive/                     displaced material: what moved, when, and why
```

`raw_data/` (GEO downloads, about 11 GB) and all regenerable `.h5ad` and
`.npz` objects are gitignored; figures, small tables and run records are
tracked. No sequencing data, count matrices or paper PDFs are committed.

## Reproducing

The evidence-contract checks use only Python's standard library:

```bash
python -m unittest discover -s analysis/tests -q
python analysis/scripts/claim_contract.py --check
python analysis/scripts/validate_repository.py
python -m compileall -q analysis Thesis
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
