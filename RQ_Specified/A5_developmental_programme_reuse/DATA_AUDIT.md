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

## Circularity

The adult injury list in the contract came from this same study. That makes the
choice of module decisive.

- **The development-specific genes are not circular here.** They come from the Guo
  list and were never selected in this study's data.
- **That test is conservative.** Building the development-specific module removed
  the developmental genes that were among this study's top 400 transitional markers,
  so the module has already lost the genes most likely to separate these cells.
- **The development-injury pair may be circular here.** Its five genes were selected
  as transitional markers by this study. Whether that selection used this
  high-resolution series or the separate whole-lung experiment is not yet
  established; the published list is named only as a supplementary sheet. If it
  used the whole-lung mice, these 36 mice are independent of it.

## Identity is the main rival

Forty-three of the 94 development-specific genes are type 1 or type 2 identity
genes. Transitional cells lean toward type 1, so a higher score in them could
simply reflect type 1-directed identity, which is the rival the contract already
found in the development-injury pair. The 51-gene version without identity genes
is fully covered and should be the primary module or at least a declared secondary.

## Design choices the A5 plan must fix before scoring

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

## Verdict

The Strunz high-resolution time course is a feasible A5 cohort. It has 26
independent injured mice meeting the cell floors, full gene coverage, and a
development-specific module that is not circular in this data. Before scoring, A5
needs its own pre-registration fixing the four choices above, and one provenance
question remains open: which Strunz dataset produced the transitional marker list.
