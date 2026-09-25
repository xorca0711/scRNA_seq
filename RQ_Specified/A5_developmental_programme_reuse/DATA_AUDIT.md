# A5 data audit: a replicated adult injury time course

25 September 2026. **Metadata only; no count was read and no score exists.** The
owner chose this audit at stage 4 of the
[shared component contract](../A5_A11_shared_component_contract/README.md)
(DEVELOPMENT decision 37). The tables come from
`scripts/01_audit_strunz_units.py` and its run record, `tables/audit_run.json`.

## Why adult injury animals are enough now

A5 asks whether adult alveolar repair reuses part of a developmental epithelial
programme. The earlier design compared neonatal and adult wells, so it needed
independent neonatal animals, and no local dataset has any. The contract changed
that. The developmental list now comes from an outside developmental study, so it
no longer has to be measured in neonatal animals here. Reuse can be tested in adult
injury alone: are the developmental genes higher in transitional cells than in type
2 cells of the same injured mouse?

Local adult injury data cannot carry that test. The viral-injury cohort reaches at
most five animals, spread over five days, and only by pooling a mixed label that
also appears in uninjured lung.

## The cohort: Strunz et al. 2020, high-resolution epithelial time course

GEO series GSE141259, from
[10.1038/s41467-020-17358-3](https://doi.org/10.1038/s41467-020-17358-3). The
sample design was already stored in the A1 metadata extract, so no new metadata
query was needed. Two small files were fetched to the ignored `cache/` folder.

| File | Bytes | Role |
|---|--:|---|
| `GSE141259_HighResolution_cellinfo.csv.gz` | 1,424,472 | author label, sample and day for every cell |
| `GSE141259_HighResolution_genes.txt.gz` | 70,000 | the 24,051 deposited gene names |

Their hashes are in the audit script and its run record. The 76 MB count matrix was
not fetched; it is scoring input.

**Each sample is one mouse.** The methods state that EpCAM-sorted epithelium was
sampled daily to day 13 and at later points to day 54, at 18 time points with two
replicate mice each, 36 mice in all. This comes from the published text, not from
the sample names. The deposit holds 32 of those samples.

## Units

The contrast needs at least 30 transitional cells and 30 type 2 cells in the same
mouse, by the authors' labels, and at least three such mice.

| Result | Value |
|---|--:|
| Deposited samples | 32 |
| Control samples | 2 |
| Injured mice meeting both floors | 26 |
| Days those mice cover | 14, from day 2 to day 21 |

The floors are counted before any depth filter. The Drop-seq libraries are shallow,
so the scoring plan must state its depth handling, and the eligible count may fall.
It starts far above the floor of three.

## Coverage

Every gene of every A5 module is present in the deposited gene list, including the
old RIKEN symbols that cost coverage elsewhere.

| Module | Genes | Present |
|---|--:|--:|
| Development-specific | 94 | 94 |
| Development-specific without identity genes | 51 | 51 |
| Shared remodelling | 12 | 12 |
| Development and injury pair | 5 | 5 |

## Provenance correction after review, before scoring

The original audit called the development-specific module non-circular and left
its Strunz experiment provenance open. This is superseded. The cached author
supplement description (`41467_2020_17358_MOESM3_ESM.pdf`) explicitly assigns
Supplementary Data 3 to the high-resolution epithelial dataset. The contract uses
that workbook's `cell_types_2` sheet for injury and identity markers.

The five shared genes were positively selected in this cohort. The 94-gene module
and 51-gene identity-excluded variant were negatively filtered using the same
cohort. Negative filtering can be conservative but is not independent selection.
The revised [plan](PLAN.md) uses Guo's external 99-gene module, Guo-only identity
exclusions (57 genes), and external control exclusions (53 genes). Old modules
remain descriptive. Strunz's discussion already reports poor correspondence with
Guo's developmental signature; this prior analysis now informs the test.

## Identity is the main rival

Forty-three of the 94 development-specific genes are type 1 or type 2 identity
genes. Transitional cells lean toward type 1, so a higher score in them could
simply reflect type 1-directed identity, which is the rival the contract already
found in the development-injury pair. The 51-gene version without identity genes
is fully covered and should be the primary module or at least a declared secondary.

## Historical design choices, now fixed in PLAN.md

1. **The reference population.** The authors split type 2 cells into resting and
   activated. Activated type 2 cells respond to injury themselves, so the choice
   changes the contrast.
2. **Identity handling.** Either make the 51-gene module primary or adjust for type
   1 identity in a declared way.
3. **Depth.** State a count budget suited to shallow Drop-seq libraries.
4. **Time.** Use one window of injury days or all eligible days, declared in
   advance, not chosen after seeing which days show an effect.

## A metadata conflict, recorded

The two control samples are listed as day 0 in the GEO sample characteristics but
as day 14, with sample names beginning NC, in the cell table. They are probably
PBS controls harvested at day 14. The within-mouse contrast does not use them, so
the conflict does not affect A5. It is recorded rather than resolved.

## Revised verdict

A feasible source-defined signature test, with 26 mice for resting AT2 and 24
for activated AT2 before depth filtering. The primary now uses activated AT2 and
independently sourced gene selection. The [plan](PLAN.md) fixes depth, time,
reference, identity/control handling and inferential limits before new counts.
