# Trial M0: what the two multiome deposits contain

Gate 0 of the DATP epigenetics branch. Read before anything is fitted.

## Why this branch is not on the Choi deposit

Choi 2020 supports its Il1r1-positive AT2 claim with ATAC-seq, but GSE144598 deposits two bigwig coverage tracks with no peaks, no reads and one pooled sample per group. Trial D0 recorded it as unusable. These two deposits are the nearest public data that can carry a chromatin question about the transitional state at single-cell resolution.

## Files

|  | accession | file | cells | genes | peaks | libraries_in_file | annotation_rows | annotation_rows_unmatched | genes_with_a_promoter_peak |
|---|---|---|---|---|---|---|---|---|---|
| 0 | GSE310539 | totalaggr | 39849 | 32287 | 178178 | 4 | 215628 | 0 | 19086 |
| 1 | GSE247130 | P9 | 23458 | 32286 | 169112 | 2 | 204390 | 0 | 18821 |
| 2 | GSE247130 | 7wk | 18427 | 32286 | 148364 | 2 | 179857 | 0 | 18217 |
| 3 | GSE247130 | SeV | 22409 | 32286 | 208941 | 2 | 252141 | 0 | 19471 |

## Libraries, with the assumed suffix map

The barcode suffix to library map is the GEO sample order. It is the
cellranger-arc default and it is NOT deposited, so it is an assumption.
Trial M1 corroborates it against a marker the data can settle and
reports the outcome; where corroboration fails the libraries are named
L1 to L4 and no condition name is attached.

|  | accession | file | barcode_suffix | assumed_library | genotype | treatment | cells | stage |
|---|---|---|---|---|---|---|---|---|
| 0 | GSE310539 | totalaggr | 1 | wildtype_PBS | wildtype | PBS | 7340 |  |
| 1 | GSE310539 | totalaggr | 2 | wildtype_SeV | wildtype | SeV | 8093 |  |
| 2 | GSE310539 | totalaggr | 3 | AP1mut_PBS | AP-1 mutant | PBS | 13622 |  |
| 3 | GSE310539 | totalaggr | 4 | AP1mut_SeV | AP-1 mutant | SeV | 10794 |  |
| 4 | GSE247130 | P9 | 1 | P9_mutant | Cebpa mutant |  | 12186 | P9 |
| 5 | GSE247130 | P9 | 2 | P9_control | control |  | 11272 | P9 |
| 6 | GSE247130 | 7wk | 1 | 7wk_mutant | Cebpa mutant |  | 7589 | 7 weeks |
| 7 | GSE247130 | 7wk | 2 | 7wk_control | control |  | 10838 | 7 weeks |
| 8 | GSE247130 | SeV | 1 | SeV_mutant | Cebpa mutant |  | 11773 | SeV infected |
| 9 | GSE247130 | SeV | 2 | SeV_control | control |  | 10636 | SeV infected |

## Peak annotation classes, per file

Peaks are called per aggregate. A peak identifier is meaningful only
inside its own file and peak identity is never compared across files.

|  | file | promoter | distal | intergenic |
|---|---|---|---|---|
| 0 | totalaggr | 22234 | 193394 | 0 |
| 1 | P9 | 21778 | 182612 | 0 |
| 2 | 7wk | 20998 | 158859 | 0 |
| 3 | SeV | 22853 | 229288 | 0 |

## Replication, under the rule inherited from trial D0

|  | accession | contrast | held_fixed | libraries_per_group | within_group_replication |
|---|---|---|---|---|---|
| 0 | GSE310539 | treatment | genotype | 1 | False |
| 1 | GSE310539 | genotype | treatment | 1 | False |
| 2 | GSE247130 | genotype | stage | 1 | False |
| 3 | GSE247130 | stage | genotype | 1 | False |
| 4 | both | cell state inside one library | library, animal, batch, peak set, sequencing run | not applicable | not required; the comparison is inside one library |

## The reading

Ten libraries across four files and 104,143 cells, and
**not one contrast between conditions carries within-group replication**,
in either deposit. This is the fifth and sixth deposit in this project
with that ceiling. Nothing about genotype, treatment or stage is testable
from these files, and trial M1 does not attempt it.

What is admissible is a comparison between cell states inside one
library, where the animal, the batch, the peak set and the sequencing
run are held fixed by construction. The precedent in this repository is
trial D3, an ordering computed inside a library. The reading is then
required to hold in every library separately rather than pooled, and the
number of libraries is the number of independent chances the claim had
to fail, not a sample size.

The peak to gene join is clean: 0 annotation rows
across all four files fail to match a peak in their own matrix.
