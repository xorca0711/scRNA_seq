# Analysis trial plan: Choi, Lee et al. 2020 (roadmap paper 2, Gate 1)

Every trial in this folder writes its frozen rules to a run record before it
reads any data, and completes the record afterwards with inputs, package
versions and outputs. A rule is never moved after the result is seen; where a
rule turns out to be wrong, the first outcome stays in the record and a
corrected pass sits beside it, which is the precedent set in the Cardoso folder
(trials C1 and C1b, C2 and C2b).

Trial identifiers in this folder are the D series, for the DATP state this
paper defines. The study note is [`README.md`](README.md) and the
machine-readable extract is
[`choi_2020_extracts.json`](choi_2020_extracts.json).

| Trial | Gate | Question | State |
|---|---|---|---|
| D0 | 0 | What does the deposit contain, before anything is fitted | complete |
| D1 | 1 | Are the five states and the AT2 to DATP to AT1 ordering recoverable with the labels held out | not started |

---

## D0. Data reality check (Gate 0)

Script [`trials/d0_data_reality_check.py`](trials/d0_data_reality_check.py);
artefacts in [`trials/d0_data_reality_check/`](trials/d0_data_reality_check/).

### Pre-registration

Rules in full in the script docstring and in the run record. In short: scope is
the two single-cell accessions; design is parsed from the deposited GEO SOFT
family files rather than transcribed; a feature whose identifier is not an
Ensembl gene ID is not a gene; a matrix whose column count equals the 10x v2
whitelist of 737,280 is the raw droplet matrix rather than called cells, and
Gate 0 reports that rather than calling cells itself; gene spaces are compared
as exact ordered lists of identifier and symbol pairs; marker presence is
checked for every set in the extract including the negative conditions; and a
contrast has within-group replication only if at least two libraries share
every experimental variable except the one contrasted.

### Outcome

**One gene space across all eight libraries**, 27,998 features, identical
ordered identifier and symbol lists. Nothing in this folder needs a gene-space
intersection, which is a better position than the Cardoso deposit, where three
distinct spaces forced one.

**No non-gene features anywhere.** This is worth more than it looks. The
lineage-tracing design turns on a tdTomato reporter, and tdTomato is not
counted as a feature in the deposit. The Tomato-positive and Tomato-negative
split is therefore a fact about the sort and cannot be verified from the
matrix. A reanalysis has to trust the library labels for lineage, which is the
opposite of the Cardoso deposit, where the BSD selection marker was counted and
could be used as an independent sort check (trial C0).

**Six of the eight matrices are raw 10x barcode whitelists, not called cells.**
All six lineage-traced libraries deposit the full 737,280-barcode matrix; the
two organoid libraries are already filtered, at 2,101 and 3,066 cells. So cell
calling for the in vivo half is this repository's job, and the cell count will
not match the paper's unless Cell Ranger 2.0.2's caller is reproduced, which it
will not be. Gate 0 deliberately did not call cells: that needs a frozen
threshold and belongs to Gate 1.

To size the problem without pre-empting it, the inventory reports how many
barcodes clear a plain count floor. At 500 counts the six libraries hold
between 5,336 and 12,592 barcodes.

| Library | Barcodes over 500 counts | Median genes per such barcode |
|---|--:|--:|
| PBS Tomato | 6,877 | 1,061 |
| PBS nonTomato | 12,592 | 1,197 |
| Day 14 Tomato | 5,336 | 1,825 |
| Day 14 nonTomato | 9,213 | 891 |
| Day 28 Tomato | 5,693 | 1,689 |
| Day 28 nonTomato | 11,177 | 968 |

Two things in that table matter for Gate 1. The Tomato-negative libraries hold
roughly twice as many barcodes as their Tomato-positive partners at every time
point, so any composition comparison across the sort is comparing very
different cell numbers. And the median genes per barcode ranges from 891 to
1,825, a twofold spread across libraries that are going to be compared to each
other. Detection fractions are depth-sensitive, so a cross-library detection
comparison in this deposit carries a depth term, and the depth control that
trial E6 of the Cardoso folder found to be load-bearing should be carried here
by default.

**Every marker gene in the extract is present in every library**, including
`Sprr1a` and `AW112010`, which are the two DATP markers most likely to be
missing from an older annotation. The negative-condition genes (`Pdpn`,
`Hopx`, `Cav1`) are present too, so the DATP definition can be scored in full
rather than from its positive half.

**No contrast in this deposit carries within-group replication.** The frozen
replicate rule was applied to three contrasts and none passed: time point with
sort held constant, sort with time point held constant, and organoid treatment.
There is one library per condition throughout.

| Accession | Contrast | Held constant | Within-group replication |
|---|---|---|---|
| GSE145031 | time point | sort | No |
| GSE145031 | sort | time point | No |
| GSE144468 | treatment | none | No |

**GSE144598 was not parsed, and that is the result.** The ATAC-seq
supplementary files are bigwig coverage tracks, with no peaks and no reads, so
the epigenetic half of the Il1r1-subset claim cannot be re-derived from the
deposit at all. Raw reads would have to come from SRA. GSE144553 is the
SuperSeries and holds no files of its own.

### A reporting defect in the first run, disclosed

The first run's summary opened with a single boolean, "matrices are raw 10x
barcode whitelists rather than called cells: False", because the check was an
all-or-nothing test across the eight libraries. That reads as though none of
them are raw, when in fact six of them are and only the two organoid libraries
are filtered. The underlying table was correct and no number changed. The
summary now names both groups and their counts, and the run record reports two
lists instead of one boolean. The defect was in how the result was stated, not
in what was measured, and it is recorded here because a Gate 0 headline that
misleads is worse than one that is merely incomplete.

### Claims from D0

| Claim | Class |
|---|---|
| The deposit has one gene space across all eight libraries, 27,998 features | Descriptive only |
| No contrast in this deposit carries within-group replication | Descriptive only, and it is the ceiling for the whole folder |
| Six of eight matrices are raw barcode whitelists, so cell calling is this repository's job | Descriptive only |
| The tdTomato reporter is not a counted feature, so the lineage split cannot be verified from the matrix | Descriptive only, and a limit on every lineage claim a reanalysis could make here |
| Every marker gene in the extract, including the negative conditions, is present in every library | Descriptive only |
| The epigenetic claim can be re-derived from the deposit | Not establishable; the ATAC deposit is coverage tracks only |
| Median genes per barcode spans 891 to 1,825 across libraries that will be compared | Descriptive only, and the reason Gate 1 must carry a depth control |

---

## D1. The states and the ordering (Gate 1), not started

The question Gate 1 will ask, stated now so the rules can be frozen before the
data are opened: with the paper's labels held out of every fitting step, does a
blind pipeline recover the five states the paper reports (homeostatic AT2,
cycling AT2, primed AT2, DATP, AT1), and does the AT2 to primed AT2 to DATP to
AT1 ordering appear in a trajectory rooted in AT2 cells?

Three things about the design are already fixed by D0. Cell calling needs a
frozen threshold, because six libraries are raw droplet matrices. The primed
AT2 state is defined by loss of `Etv5`, `Abca3` and `Cebpa` rather than by a
positive marker, so it cannot be scored the way the other four are, and the
scoring rule has to say what it does about that before it is run. And the
ordering question is the one part of this paper a deposit without replication
can genuinely address, because an ordering inside a library needs no
between-group comparison, whereas any statement about how much of a state each
time point holds is a difference between two libraries.
