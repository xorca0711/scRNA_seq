# Thesis roadmap: one folder per paper, in the order the project reads them

This directory builds the project up paper by paper. Each paper gets one
folder holding (1) a study note in the roadmap's five-question format
(question tested, evidence type, reusable variables, one limitation, one
bridge), (2) any parameters or decision criteria extracted from the paper as
reviewable JSON, and (3) the analysis trial that paper motivates, with its
pre-registered plan, its logged outcome, and its negative results. Folders
are added one at a time in the order below, which is the ordered paper
roadmap of the owner's Notion page
[UC Berkeley SAP PI Target Map](https://app.notion.com/p/3d1151616b44814c8697ff3af2f8f831)
(gates 1, 2, 3A, 3B). The same table is machine-readable in
[`ROADMAP.json`](ROADMAP.json).

Paper PDFs and supplementary spreadsheets stay on disk and are gitignored
(`*.pdf`, `*.xlsx`); only text, JSON and small tables are tracked. The
Notion pages hold the owner's reading notes; this directory holds the
reproducible parts.

**Working question (from the roadmap).** Which epithelial and immune-state
programmes distinguish productive lung repair from persistent dysplastic
remodelling after inflammatory injury?

## Order and status

| # | Gate | Paper | DOI | PMID | Role in the roadmap | Folder | Study note | Analysis trial |
|--:|---|---|---|---|---|---|---|---|
| 1 | 1 | Niethamer et al. 2025, *Cell Stem Cell* | [10.1016/j.stem.2024.12.002](https://doi.org/10.1016/j.stem.2024.12.002) | 39818203 | source paper for GSE262927: animals, time points, annotations, known findings, limits | [`gate1_01_niethamer_2025/`](gate1_01_niethamer_2025/README.md) | done (in `docs/`) | done (`analysis/GSE262927/`); PI-matched follow-ups N1 to N4 run 2026-09-10, Descriptive only, owner review pending: [`ANALYSIS_TRIAL_PLAN.md`](gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md) |
| 2 | 1 | Choi et al. 2020, *Cell Stem Cell* | [10.1016/j.stem.2020.06.020](https://doi.org/10.1016/j.stem.2020.06.020) | 32750316 | biological spine: IL-1beta/HIF1alpha-driven AT2 to DATP to AT1 transition | not started | not started | not started |
| 3 | 1 | Nabhan et al. 2018, *Science* | [10.1126/science.aam6603](https://doi.org/10.1126/science.aam6603) | 29420258 | fibroblast Wnt niches maintain AT2 stemness; niche exit permits AT1 differentiation | not started | not started | not started |
| 4 | 1 | Sikkema et al. 2023, *Nature Medicine* (HLCA) | [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) | 37291214 | reference framework: annotation hierarchy, reference mapping, uncertainty handling, donor coverage, shared profibrotic macrophage states | [`gate1_04_sikkema_2023_hlca/`](gate1_04_sikkema_2023_hlca/README.md) | done, owner review pending | S1, S3, S4, S5 run (Descriptive only); S2 planned; see [`ANALYSIS_TRIAL_PLAN.md`](gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md) |
| 5 | 2 | Cardoso, Lee et al. 2026, *Nature* | [10.1038/s41586-026-10399-6](https://doi.org/10.1038/s41586-026-10399-6) | 42020743 | early fibrotic niches; regenerative-like mutant AT2 states coordinate fibroblast and immune remodelling through AREG-EGFR | not started | not started | not started |
| 6 | 2 | Nabhan et al. 2023, *Cell* | [10.1016/j.cell.2023.05.022](https://doi.org/10.1016/j.cell.2023.05.022) | 37321220 | Frizzled-specific Wnt agonists separate regeneration from fibrotic risk; receptor-specific Wnt modules | not started | not started | not started |
| 7 | 3A | Saxton et al. 2021, *Science* | [10.1126/science.abc8433](https://doi.org/10.1126/science.abc8433) | 33737461 | structure-based decoupling of IL-10 pro- and anti-inflammatory functions | not started | not started | not started |
| 8 | 3A | Saxton et al. 2021, *Immunity* | [10.1016/j.immuni.2021.03.008](https://doi.org/10.1016/j.immuni.2021.03.008) | 33852830 | IL-22 tissue-protective vs pro-inflammatory functions decoupled | not started | not started | not started |
| 9 | 3B | DuPage et al. 2015, *Immunity* | [10.1016/j.immuni.2015.01.007](https://doi.org/10.1016/j.immuni.2015.01.007) | 25680271 | Ezh2 maintains regulatory T cell identity after activation | not started | not started | not started |
| 10 | 3B | Wang et al. 2018, *Cell Reports* | [10.1016/j.celrep.2018.05.050](https://doi.org/10.1016/j.celrep.2018.05.050) | 29898397 | targeting EZH2 reprograms intratumoral Tregs | not started | not started | not started |
| 11 | 3B | Zhang et al. 2026, *Science Immunology* | [10.1126/sciimmunol.adx4411](https://doi.org/10.1126/sciimmunol.adx4411) | 41961946 | intratumoral Treg ablation elicits NK-mediated control | not started | not started | not started |

Outside the roadmap but already in the repository: Kadur Lakshminarasimha
Murthy et al. 2022, *Nature* (DOI
[10.1038/s41586-022-04541-3](https://doi.org/10.1038/s41586-022-04541-3),
PMID 35355018), the source paper for GSE178360, analysed in
[`../analysis/GSE178360/`](../analysis/GSE178360/README.md). The HLCA lists
that series as one of its extension datasets (Tata_unpubl), which is why the
Sikkema trial S2 targets it.

DOIs and PMIDs were verified against PubMed on 2026-09-09.

## Gate rules carried over from the roadmap

- **Analysis gate.** After papers 1 to 4, write a one-page analysis contract
  and begin with deposited count matrices plus metadata, not FASTQ. The
  initial question stays narrow: do early versus late KRT8 transitional
  epithelial states associate with different macrophage and fibroblast niche
  programmes in GSE262927?
- **Stop rule.** Finish papers 1 to 4, write the contract, start the pilot.
  Read papers 5 to 11 only when a specific pilot result makes their branch
  relevant.
- **Gate 3A** opens only if the cytokine or resolution signal is strong.
  **Gate 3B** opens only if regulatory T cells are sufficiently represented;
  note that the HLCA core could not separate Tregs from other T cells, which
  is a known ceiling for any reference-based route to that branch.

## Folder contract

```
Thesis/
  README.md                       this index
  ROADMAP.json                    the same table, machine-readable
  gateG_NN_firstauthor_year/
    README.md                     study note: five questions, plus the extracts the roadmap names for that paper
    *.json                        parameters and decision criteria extracted from the paper, reviewable
    PIPELINE_FRAMING.md           (when a paper is a methods reference) what this repository adopts, adapts, declines
    ANALYSIS_TRIAL_PLAN.md        pre-registered trial(s): rule, threshold, dataset, success and failure criteria, outcome
    trials/                       trial scripts and their logged artefacts (small tables and JSON are tracked)
```

Status vocabulary is the one used in the root README claims table:
Validated, Descriptive only, Exploratory, Retracted-superseded, Not
established. A trial result with no artefact under `trials/` is Not
established. Thresholds are frozen in the plan before the trial reads data.

## Local, untracked material

The PDFs and supplementary files for the roadmap papers live in the owner's
external thesis-study folder (organised by gate) and, for the earlier
papers, directly under this directory; both locations are ignored by git.
`docs/scRNAseq_workflow_Niethamer2025.md` records a filing correction for
one mislabelled local PDF.
