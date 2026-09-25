# Trials index: the Axin2 and Il1r1 branch

The branch note is [`../README.md`](../README.md) and the Route B outcome is
[`../ROUTE_B_OUTCOME.md`](../ROUTE_B_OUTCOME.md). The full pre-registration for
any trial is the docstring at the top of its script.

**This branch mostly assessed rather than ran.** Its conclusion is that the
definitive comparison cannot be made from public data, and the trials below
exist to test the routes that might have substituted for it. Six of the seven
returned NOT COMPUTABLE or closed a route.

## Measurements, which compare no states and set no thresholds

| Script | Question | What it returned |
|---|---|---|
| [`feasibility_of_the_joint_question.py`](feasibility_of_the_joint_question.py) | Is locus-level chromatin less sparse than the transcript for these two genes | **Yes for Axin2**, 4.1 and 10.4 fold, with three controls ordering correctly in both deposits |
| [`module_detection_check.py`](module_detection_check.py) | Does a Wnt-target module clear the per-cell detection floor Axin2 alone does not | **Yes**, 14.6 to 28.4 per cent against 11.1 to 18.0 for Axin2 alone |

## Route B, locus co-accessibility: three passes, closed

| Trial | Instrument | What it returned |
|---|---|---|
| [A1](a1_locus_co_accessibility.py) | binary: is at least one of the locus's peaks detected | **Unreadable.** Its physical control Krt8 against Krt18 did not clear, and was the wrong control: nearest-gene annotation splits the intergenic peaks rather than sharing them |
| [A1b](a1b_split_locus_control.py) | the same, gated on a split-locus control | **Not computable**, the split-locus control failed at z +0.51 to +2.05 |
| [A1c](a1c_graded_statistic.py) | graded: the fraction of the locus's peaks detected, correlated across cells | **Not computable** at z +0.24 to +1.62, and Route B closed under A1c's own pre-registered rule. The graded form does clear on Etv5 against Abca3 at z +3.55 and Krt8 against Krt18 at z +3.59, so the instrument has sensitivity its gate could not certify |

## Route A, England et al.: the data does not exist

| Trial | Question | What it returned |
|---|---|---|
| [A2](a2_does_il1r1_deletion_move_wnt.py) | Does deleting Il1r1 move the Wnt programme in Kras-mutant clones | **Not computable.** First, the England Figure S2 data this route was built on was never deposited. Then the replacement refused at its own label rule: **Il1r1 transcript does not distinguish the homozygous deletion from the heterozygous one**, so the genotype labels cannot be verified from the matrices |

## Shared code

| File | Role |
|---|---|
| [`write_route_b_summary.py`](write_route_b_summary.py) | regenerates the Route B outcome page from the three trials' tracked CSVs |

The heavy machinery (the streaming reader, the peak-to-gene join, downsampling)
is imported from the sibling branch at
[`../../datp_epigenetics/trials/multiome_utils.py`](../../datp_epigenetics/trials/multiome_utils.py)
rather than duplicated.

## Conventions

The same as the sibling branch: frozen rules in the docstring before any data is
read, a run record per trial, thresholds never moved after a result is seen, an
unreadable-if clause on every threshold, the library as the statistical unit,
and every tracked artefact a table, a figure or a run record. Two further rules
this branch added at its own cost:

- **A positive control must be demonstrated, not asserted.** A1 assumed one was
  an arithmetic identity and it was not; A1b asserted its replacement could not
  fail and it did.
- **Check that the data exists before proposing a route on it.** Route A was
  proposed on a deposit that does not contain the figure it was built around.
