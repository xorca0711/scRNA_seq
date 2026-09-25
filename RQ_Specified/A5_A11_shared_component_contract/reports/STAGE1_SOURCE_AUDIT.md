# Stage 1: developmental source audit

25 September 2026. Stage 1 of the [shared component contract](../PLAN.md) required
an independently sourced developmental signature. The existing specificity module
recorded that none was available locally. This audit searched externally and
records every candidate examined, including the six not used.

**Outcome: source found.** The selected list is the author-defined signature of a
mixed type 1 and type 2 population in normal mouse lung at postnatal day 1, from
Guo et al. 2019. Nothing in this audit scored any expression data.

## Selection criteria, fixed before searching

1. Normal development, with no injury model contributing to the list.
2. A laboratory that did not define the injury or lesion lists.
3. A complete author-defined gene list, not genes picked from a figure.
4. Retrievable through an open archive, so the file can be hashed.
5. Biologically the right object: a developmental state of mixed or intermediate
   alveolar identity, the analogue of the adult transitional state.

## Candidates

| Study | DOI | Outcome | Reason |
|---|---|---|---|
| Treutlein et al. 2014 | [10.1038/nature13173](https://doi.org/10.1038/nature13173) | not used | Its two supplementary tables are genes expressed in every cell and a qPCR primer panel. Neither is a state marker list. |
| Negretti et al. 2021 | [10.1242/dev.199512](https://doi.org/10.1242/dev.199512) | not retrievable | Not open access in Europe PMC, and the publisher blocks full-text XML. |
| Frank et al. 2019 | [10.1073/pnas.1813952116](https://doi.org/10.1073/pnas.1813952116) | not retrievable | Not open access in Europe PMC, and the publisher blocks full-text XML. |
| Zepp et al. 2021 | [10.1126/science.abc3172](https://doi.org/10.1126/science.abc3172) | not used | The author manuscript supplement contains figures and legends only. |
| Penkala et al. 2021 | [10.1016/j.stem.2021.04.026](https://doi.org/10.1016/j.stem.2021.04.026) | not a list source | A neonatal injury study, so it fails criterion 1. Recorded below as a biological rival to A5. |
| Hurskainen et al. 2021 | [10.1038/s41467-021-21865-2](https://doi.org/10.1038/s41467-021-21865-2) | rejected on inspection | Its five epithelial clusters are type 2, type 1, a Lyz1-high type 2, ciliated and club. None is transitional, and the markers pool normal and hyperoxia-injured lungs. |
| **Guo et al. 2019** | [10.1038/s41467-018-07770-1](https://doi.org/10.1038/s41467-018-07770-1) | **selected** | Normal postnatal day 1 lung, an independent laboratory, a complete author list per subtype, and an open archive. |

Hurskainen was rejected after its marker file was opened, not on the expectation
that it would fail. Pooling injured cells into a developmental list would build the
injury programme into the development side, which is the circularity A5 must avoid.

## The selected list

The authors list 100 signature genes per subtype, ranked by a binomial test with a
false discovery rate below 0.1. The mixed population has 460 cells.

A population with mixed identity in droplet data can be doublets. A doublet of a
type 1 and a type 2 cell would give a signature close to the union of both identity
lists, and would pull type 2 genes in heavily. The list was checked against the
authors' own subtype lists only, before any comparison with injury or lesion lists.

| Check | Genes of 100 |
|---|--:|
| Unique to this population across all 20 author subtypes | 54 |
| Shared with the author type 1 list | 37 |
| Shared with the author type 2 list | 6 |
| Operational genes present | 1 |

The operational gene is Cldn4, which the freeze removes. The lopsided overlap and
the majority of unique genes argue against a doublet union. They do not prove the
population is a single state, and the check uses the authors' top 100 lists, so
unique means absent from other lists rather than unexpressed elsewhere.

The unique genes include Areg, Egfr, Erbb3, Itgb6, Shh and Fzd5. They also include
Rprm, a p53 target, which is why the p53 control axis must be scored separately.
These names are recorded to show the list is not trivial. They are not a finding.

## Limitations of this selection

- **Selection by accessibility.** Negretti 2021 covers embryonic day 12 to
  postnatal day 14 with over 100,000 cells and would likely be the stronger source.
  It was excluded because it is not openly retrievable, not on scientific grounds.
  Retrieving its tables manually would make a useful sensitivity source.
- **One timepoint.** The list describes postnatal day 1, around birth. It does not
  describe the alveolar stage that follows. The neonatal test wells are at postnatal
  day 9, a different stage, which lowers circularity but also means the list may
  miss later maturation.
- **Perinatal stress.** The source paper reports unfolded protein response
  activation at birth. Part of any sharing with injury could be generic stress,
  which the Hallmark control axes exist to absorb.
- **Old symbols.** The list uses 2019 annotation, including RIKEN names. Symbols
  are kept exactly, so any that no longer match lower coverage, as they should.

## Biological rival recorded for A5

Penkala et al. 2021 report that after neonatal injury, type 1 cells reprogram into
type 2 cells, and that type 2 cells regenerate type 1 cells only in the mature lung.
Plasticity there is a trait acquired with age. If that holds generally, the neonatal
and adult lung repair in opposite lineage directions. A5 concerns a shared
transcriptional component, not lineage direction, so this does not refute it. It
does mean a shared component would not imply a shared route, and A5 must say so.

## Retrieval record

| Archive | Bytes | SHA-256 |
|---|--:|---|
| Guo 2019 supplementary archive | 12,979,470 | `5942c900289b15b35ee1405067263b35037ee815eb5678fb5258985fad6afd51` |
| Hurskainen 2021 supplementary archive | 66,277,402 | `396074f75839006a1f37244a5f7a5aa7af838026047c767ed35d93cc107907da` |

Both archives came from the Europe PMC supplementary files service, the same route
the repository used for the Strunz lists. Only the Guo member file is cached for the
freeze, under the untracked `sources/` directory; its hash is tracked in the
specification. Article metadata came from PubMed.
