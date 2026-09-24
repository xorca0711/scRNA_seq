# Thesis roadmap: one folder per paper, in the order the project reads them

This directory builds the project up paper by paper. Each paper gets one
folder holding (1) a study note in the roadmap's five-question format
(question tested, evidence type, reusable variables, one limitation, one
bridge), (2) any parameters or decision criteria extracted from the paper as
reviewable JSON, and (3) the analysis trial that paper motivates, with its
pre-registered plan, its logged outcome, and its negative results. Folders
are added one at a time in the order below, with owner-directed exceptions recorded in
the table: paper 5 (Cardoso 2026) was entered on the owner's instruction of
2026-09-12, ahead of papers 2, 3 and 6; and paper 2 (Choi 2020) was entered by
the agent on 2026-09-13 on a general "proceed", then withdrawn on 2026-09-15
until the owner had read it, with its deposit check kept under paper 5 as an
extension (DEVELOPMENT decision 21), then re-entered later the same day at
the owner's direction after reading (decision 23). The order below is the owner's
reading order
(gates 1, 2C, 2N, 2W, 3A, 3B; the Gate 2 branches were added on 2026-09-15). The same table is machine-readable in
[`ROADMAP.json`](ROADMAP.json).

On 24 September 2026, after reading paper 13, the owner opened its planning
branch as [`gate2_C3_yu_lee_choi_min_2026/`](gate2_C3_yu_lee_choi_min_2026/README.md).
`C3` denotes the third item in branch 2C; stable roadmap order 13 is unchanged.
The owner subsequently authorized staged execution. Mouse QC/clustering,
two-cohort IPF pathways and an initial 25-donor ligand-receptor analysis are
complete; treatment inference, sensitivities and human/spatial work remain
pending. See the [initial report](gate2_C3_yu_lee_choi_min_2026/INITIAL_RUN_REPORT.md).

Paper PDFs and supplementary spreadsheets stay on disk and are gitignored
(`*.pdf`, `*.xlsx`); only text, JSON and small tables are tracked. The
owner's reading notes are kept privately; this directory holds the
reproducible parts.

**Working question (from the roadmap).** Which epithelial and immune-state
programmes distinguish productive lung repair from persistent dysplastic
remodelling after inflammatory injury?

## Order and status

| # | Gate | Paper | DOI | PMID | Role in the roadmap | Folder | Study note | Analysis trial |
|--:|---|---|---|---|---|---|---|---|
| 1 | 1 | Niethamer et al. 2025, *Cell Stem Cell* | [10.1016/j.stem.2024.12.002](https://doi.org/10.1016/j.stem.2024.12.002) | 39818203 | source paper for GSE262927: animals, time points, annotations, known findings, limits | [`gate1_01_niethamer_2025/`](gate1_01_niethamer_2025/README.md) | done (in `docs/`) | done (`Thesis/gate1_01_niethamer_2025/GSE262927/`); phase and myeloid follow-ups N1 to N4 run 2026-09-10, Descriptive only, owner review pending: [`ANALYSIS_TRIAL_PLAN.md`](gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md) |
| 2 | 1 | Choi et al. 2020, *Cell Stem Cell* | [10.1016/j.stem.2020.06.020](https://doi.org/10.1016/j.stem.2020.06.020) | 32750316 | biological spine: IL-1beta/HIF1alpha-driven AT2 to DATP to AT1 transition | [`gate1_02_choi_2020/`](gate1_02_choi_2020/README.md) (re-entered 2026-09-15 at the owner's direction after reading, decision 23; the AI-written note of 2026-09-13 was withdrawn first and stays in git history, PR #19) | done (2026-09-15, owner-directed) | trials D0 to D7 run 2026-09-15 with a corrected annotation pass D2b beside D2; Gate 1 not fully recovered from the deposit (four of five states; primed AT2 never assigned as a cluster); Descriptive only, owner review pending: [`ANALYSIS_TRIAL_PLAN.md`](gate1_02_choi_2020/ANALYSIS_TRIAL_PLAN.md). Two branches opened 2026-09-20 on other laboratories' multiome deposits, because this paper's ATAC deposit is coverage tracks only: [`datp_epigenetics/`](gate1_02_choi_2020/datp_epigenetics/README.md) (M0 to M4; one validated deposit error, the chromatin question Not established with a direction) and [`axin2_il1r1/`](gate1_02_choi_2020/axin2_il1r1/README.md) (this paper's own closing Discussion question; assessed, two routes run and refused). Rows C116 to C154 |
| 3 | 1 | Nabhan et al. 2018, *Science* | [10.1126/science.aam6603](https://doi.org/10.1126/science.aam6603) | 29420258 | fibroblast Wnt niches maintain AT2 stemness; niche exit permits AT1 differentiation | [`gate1_03_nabhan_2018/`](gate1_03_nabhan_2018/README.md) | owner completed reading 2026-09-22; notes remain private | Source reproduction, descriptive Nb1 and external cohort eligibility: [analysis report](gate1_03_nabhan_2018/README.md) |
| 4 | 1 | Sikkema et al. 2023, *Nature Medicine* (HLCA) | [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) | 37291214 | reference framework: annotation hierarchy, reference mapping, uncertainty handling, donor coverage, shared profibrotic macrophage states | [`gate1_04_sikkema_2023_hlca/`](gate1_04_sikkema_2023_hlca/README.md) | done, owner review pending | S1 to S5 run (Descriptive only); see [`ANALYSIS_TRIAL_PLAN.md`](gate1_04_sikkema_2023_hlca/ANALYSIS_TRIAL_PLAN.md) |
| 5 | 2C | Cardoso, Lee et al. 2026, *Nature* | [10.1038/s41586-026-10399-6](https://doi.org/10.1038/s41586-026-10399-6) | 42020743 | early fibrotic niches; regenerative-like mutant AT2 states coordinate fibroblast and immune remodelling through AREG-EGFR | [`gate2_05_cardoso_2026/`](gate2_05_cardoso_2026/README.md) | done, owner review pending | C0 to C6 plus the E series run 2026-09-12 and 2026-09-13; Gate 1 returned **not recovered**, and six extension trials left this deposit for public data because it carries no within-group replication. E1 on GSE131907 is the only tested claim (p = 0.0020); E4 on GSE132771 refuted the second-signal reading of claim C29: [`ANALYSIS_TRIAL_PLAN.md`](gate2_05_cardoso_2026/ANALYSIS_TRIAL_PLAN.md) |
| 12 | 2C | England et al. 2025, *Cell Stem Cell* | [10.1016/j.stem.2025.01.011](https://doi.org/10.1016/j.stem.2025.01.011) | 39978341 | companion of paper 5: NF-kappaB separates tumour initiation from regeneration; its deposit GSE247505 is the only replicated time course in the Cardoso set and is already used by trial C3 | not started | not started | none planned; its data are analysed under paper 5 |
| 13 | 2C | Yu, Lee, Choi_Min and Choi 2026, *Seminars in Immunology* (review) | [10.1016/j.smim.2026.102050](https://doi.org/10.1016/j.smim.2026.102050) | 42497497 | review-motivated IL-1beta perturbation, state specificity and niche-context questions | [`gate2_C3_yu_lee_choi_min_2026/`](gate2_C3_yu_lee_choi_min_2026/README.md) | owner read; synthesis prepared; staged execution authorized 2026-09-24 | [initial batch complete](gate2_C3_yu_lee_choi_min_2026/INITIAL_RUN_REPORT.md): mouse QC/clustering, IPF pathways and descriptive LR; remaining inference and extensions pending |
| 6 | 2N | Nabhan et al. 2023, *Cell* | [10.1016/j.cell.2023.05.022](https://doi.org/10.1016/j.cell.2023.05.022) | 37321220 | Frizzled-specific Wnt agonists separate regeneration from fibrotic risk; receptor-specific Wnt modules | queued after paper 3 (re-ranking 2026-09-15) | not started | not started |
| 14 | 2N | Nabhan et al. 2026, *PNAS* | [10.1073/pnas.2606113123](https://doi.org/10.1073/pnas.2606113123) | 42418498 | an alveolosphere screen of 201 genes with chimeric RNA-seq of stem-cell effects on the fibroblast niche | not started | not started | not started |
| 15 | 2W | Wagner et al. 2021, *Cell* | [10.1016/j.cell.2021.05.045](https://doi.org/10.1016/j.cell.2021.05.045) | 34216539 | Compass: immune-metabolic state inference from single-cell RNA, the method behind proposal W1 | not started | not started | proposal W1 (Stage 2 of paper 1's plan) |
| 16 | 2W | Yadav et al. 2025, *JCI* | [10.1172/JCI188734](https://doi.org/10.1172/JCI188734) | 40875483 | the lung-fibrosis myeloid-to-mesenchymal ARG1 and ornithine circuit, with Wagner as co-author; the comparison proposal W1 names | not started | not started | proposal W1 (Stage 2 of paper 1's plan) |
| 7 | 3A | Saxton et al. 2021, *Science* | [10.1126/science.abc8433](https://doi.org/10.1126/science.abc8433) | 33737461 | structure-based decoupling of IL-10 pro- and anti-inflammatory functions | paused (Gate 3A paused 2026-09-15; opens if proposal S1 separates repair phases at the animal level) | not started | not started |
| 8 | 3A | Saxton et al. 2021, *Immunity* | [10.1016/j.immuni.2021.03.008](https://doi.org/10.1016/j.immuni.2021.03.008) | 33852830 | IL-22 tissue-protective vs pro-inflammatory functions decoupled | paused (Gate 3A paused 2026-09-15; opens if proposal S1 separates repair phases at the animal level) | not started | not started |
| 9 | 3B | DuPage et al. 2015, *Immunity* | [10.1016/j.immuni.2015.01.007](https://doi.org/10.1016/j.immuni.2015.01.007) | 25680271 | Ezh2 maintains regulatory T cell identity after activation | paused (Gate 3B paused 2026-09-15; opens if proposal D1 finds Tregs separable in an external series) | not started | not started |
| 10 | 3B | Wang et al. 2018, *Cell Reports* | [10.1016/j.celrep.2018.05.050](https://doi.org/10.1016/j.celrep.2018.05.050) | 29898397 | targeting EZH2 reprograms intratumoral Tregs | paused (Gate 3B paused 2026-09-15; opens if proposal D1 finds Tregs separable in an external series) | not started | not started |
| 11 | 3B | Zhang et al. 2026, *Science Immunology* | [10.1126/sciimmunol.adx4411](https://doi.org/10.1126/sciimmunol.adx4411) | 41961946 | intratumoral Treg ablation elicits NK-mediated control | paused (Gate 3B paused 2026-09-15; opens if proposal D1 finds Tregs separable in an external series) | not started | not started |

## Figure galleries

The paper pages now bring existing analysis figures together with their
interpretive limits: [Niethamer](gate1_01_niethamer_2025/README.md#figure-gallery),
[Choi](gate1_02_choi_2020/README.md#figure-gallery),
[Nabhan](gate1_03_nabhan_2018/README.md#figure-gallery),
[Sikkema](gate1_04_sikkema_2023_hlca/README.md#figure-gallery),
[Cardoso](gate2_05_cardoso_2026/README.md#figure-gallery),
[Yu, Lee, Choi_Min](gate2_C3_yu_lee_choi_min_2026/README.md#figure-gallery) and
[Murthy](ungated_murthy_2022/README.md#figure-gallery).
Choi's [chromatin](gate1_02_choi_2020/datp_epigenetics/README.md#figure-gallery)
and [Axin2/Il1r1](gate1_02_choi_2020/axin2_il1r1/README.md#figure-gallery)
branches have their own galleries; the
[cross-study specificity report](epithelial_state_specificity/README.md)
already embeds its evidence figure. Papers without an executed analysis have
no result gallery.

The `#` is a stable identifier assigned when a paper enters the roadmap, and
folder names carry it; the row order above is the reading order. Papers 12 to
16 were added on 2026-09-15 (see the re-ranking bullet below).

### Methods references, read at the backbone step that uses them

Not gated papers: each is read when the backbone step it supports is built.
Identifiers verified against PubMed on 2026-09-15.

| Ref | Paper | DOI | PMID | Backbone step | Used in this repository |
|---|---|---|---|---|---|
| M1 | Squair et al. 2021, *Nature Communications*: confronting false discoveries in single-cell differential expression | [10.1038/s41467-021-25960-2](https://doi.org/10.1038/s41467-021-25960-2) | 34584091 | 3, sample-aware pseudobulk | the animal-as-unit rule throughout; no pseudobulk differential expression yet (proposal W1) |
| M2 | Lotfollahi et al. 2022, *Nature Biotechnology*: scArches reference mapping | [10.1038/s41587-021-01001-7](https://doi.org/10.1038/s41587-021-01001-7) | 34462589 | 2, reference mapping | trial S2 |
| M3a | Dimitrov et al. 2024, *Nature Cell Biology*: LIANA+ | [10.1038/s41556-024-01469-w](https://doi.org/10.1038/s41556-024-01469-w) | 39223377 | 5, communication | trial C12 (liana) |
| M3b | Jin et al. 2021, *Nature Communications*: CellChat | [10.1038/s41467-021-21246-9](https://doi.org/10.1038/s41467-021-21246-9) | 33597522 | 5, communication | the Cardoso paper's tool; R-only, never run here |
| M4 | Zaiss et al. 2015, *Immunity*: amphiregulin in immunity, inflammation and repair | [10.1016/j.immuni.2015.01.020](https://doi.org/10.1016/j.immuni.2015.01.020) | 25692699 | 5, the constraint on epithelium-centric AREG readings | claim C45 |
| M5a | Kobayashi et al. 2020, *Nature Cell Biology*: the transitional (PATS) state | [10.1038/s41556-020-0542-8](https://doi.org/10.1038/s41556-020-0542-8) | 32661339 | 4, the Krt8 transitional programme | `docs/DOUBLETS_AND_SCRUBLET.md`; the displaced Krt8 trajectory |
| M5b | Strunz et al. 2020, *Nature Communications*: the Krt8+ transitional state | [10.1038/s41467-020-17358-3](https://doi.org/10.1038/s41467-020-17358-3) | 32678092 | 4, the Krt8 transitional programme | same |
| M6 | Tsukui et al. 2020, *Nature Communications*: collagen-producing lung cell atlas | [10.1038/s41467-020-15647-5](https://doi.org/10.1038/s41467-020-15647-5) | 32317643 | 4, fibroblast states | trials E4, C9, C10 |
| M7 | Vaughan et al. 2015, *Nature*: lineage-negative progenitors after major injury | [10.1038/nature14112](https://doi.org/10.1038/nature14112) | 25533958 | 4, the KRT5 dysplastic programme | not used |
| M8 | van den Brink et al. 2017, *Nature Methods*: dissociation-induced gene expression in tissue subpopulations | [10.1038/nmeth.4437](https://doi.org/10.1038/nmeth.4437) | 28960196 | 2, the dissociation-stress list any stress-like programme reading needs | named by Choi 2020 attack A3 (trial D7), not attempted because the list is in the paper's supplement and not on disk |


Outside the roadmap but already in the repository: Kadur Lakshminarasimha
Murthy et al. 2022, *Nature* (DOI
[10.1038/s41586-022-04541-3](https://doi.org/10.1038/s41586-022-04541-3),
PMID 35355018), the source paper for GSE178360, analysed in
[`Thesis/ungated_murthy_2022/GSE178360/`](ungated_murthy_2022/GSE178360/README.md). The HLCA lists
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
- **Re-ranking of 2026-09-15.** The owner re-ranked the reading order.
  Consequences for this roadmap: papers 3 and 6 (Nabhan 2018 and
  2023) move to the front of the queue, with proposal Nb1 as their trial;
  Gate 3A and 3B are paused rather than closed, since their opening
  conditions stand and are exactly what proposals S1 and D1 test; the
  analysis contract and the KRT8 pilot of the analysis gate are still
  unwritten and remain the stop rule's precondition. Paper 5 was completed
  out of order before this re-ranking and is not affected.
- **Gate 2 branches (2026-09-15).** Gate 2 is split into 2C (the Choi
  axis: paper 5 done, papers 12 and 13 added), 2N (Nabhan: papers 6 and
  14) and 2W (Wagner: papers 15 and 16, added because the method behind
  proposal W1 had no paper in the original order). Folder names keep the
  plain gate label (`gate2_`); the branch letter lives in the table and in
  `ROADMAP.json`. Papers in 2C and 2N are read while interpreting the pilot
  figures, as before; 2W is read alongside proposal W1.

## Folder contract

```
Thesis/
  README.md                       this index
  ROADMAP.json                    the same table, machine-readable
  gateG_NN_firstauthor_year/       G is 1, 2 (any Gate 2 branch), 3A or 3B
    README.md                     study note: five questions, plus the extracts the roadmap names for that paper
    *.json                        parameters and decision criteria extracted from the paper, reviewable
    PIPELINE_FRAMING.md           (when a paper is a methods reference) what this repository adopts, adapts, declines
    ANALYSIS_TRIAL_PLAN.md        pre-registered trial(s): rule, threshold, dataset, success and failure criteria, outcome
    trials/                       trial scripts and their logged artefacts (small tables and JSON are tracked)
    <accession>/                  the deposit itself, when this paper produced it (paper 1 holds GSE262927)
  ungated_firstauthor_year/       a source paper outside the roadmap, made only to hold its deposit beside it
    README.md                     a pointer note, not a study note; the study note waits for the reading (Murthy 2022 holds GSE178360)
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
