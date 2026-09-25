# A1 robustness batch: source influence and promoter dependence

25 September 2026. The owner authorized the next executable batch after the
Notion summary was shortened. This is a source-informed descriptive extension
of the [second batch](SECOND_BATCH_REPORT.md), delivered after PR #70 was merged.
No primary RNA model, historical numerical result or claim grade was replaced.

## What changed in the interpretation

The CDKN1A acetylation lead is **promoter-dependent**: at a second annotated
transcript start, the iATCs-versus-iAT2 H3K27ac/H3 effect reverses in CUT1.
H3K27me3/H3 remains lower at both starts in both preparations. HPCS composition
is particularly source-sensitive in the Slc4a11 14wk group. Its chase comparison
also has a design limitation that recovering mouse identities alone cannot fix:
chase is completely represented by source-library indicators.

These results narrow the candidate claims. They do not establish a common
epigenetic state, a transition rate or regulatory causality.

## Design identifiability

The [source-level design audit](../tables/robustness_2026-09-25/hpcs/design_rank_audit.tsv)
uses the existing 22 source aliases and six libraries. Every library contains
only one chase interval. Adding chase to unrestricted library indicators adds
zero rank: 6 → 6 overall, 5 → 5 within Slc4a11, and 2 → 2 within Hopx.
Thus the proposed library-adjusted chase coefficient is not separately estimable.
A model omitting or constraining library effects would require additional
assumptions; none was fitted. This check does not assert that the magnitude of
a library effect is known.

Driver adds one algebraic dimension overall (6 → 7), supported by the mixed
IGO17402 library. That is not population validation: biological-unit identities
are still unresolved, and a within-library contrast does not establish
transportability across ages or experiments. The complete
[library/group support table](../tables/robustness_2026-09-25/hpcs/library_group_support.tsv)
retains these restrictions.

## HPCS source influence

The original 5,333 cells, 22 sources and all eight biological RNA-state labels
were retained. Each source was omitted once within its group. Across three
annotation fields this gives **22 distinct source omissions and 704
source/category omission rows**, with 192 group/category summaries. No cells
were resampled as independent animals, and no p-values were calculated.

| Deposited group | Equal-source HPCS % | Cell-pooled HPCS % | Equal-source omission range % | Cell-pooled omission range % |
|---|---:|---:|---:|---:|
| Slc4a11 6wk_3d | 93.08 | 92.83 | 92.24–93.98 | 92.50–93.33 |
| Slc4a11 8wk | 24.00 | 23.99 | 22.74–25.72 | 23.22–25.17 |
| Slc4a11 12wk_3d | 67.86 | 66.82 | 65.78–72.12 | 62.54–73.42 |
| Slc4a11 14wk | 43.73 | 31.79 | 38.54–51.84 | 29.36–51.11 |
| Hopx_12wk3d | 4.77 | 4.81 | 3.98–5.56 | 3.98–5.56 |
| Hopx_12wk14d | 4.32 | 3.84 | 2.00–6.65 | 2.00–6.65 |

These are sensitivity ranges, not confidence intervals or temporal effect
estimates. Omitting `BF1303_B0301_GFP+`, which supplies 811/991 cells in the
14wk group, raises the pooled HPCS fraction from 31.79% to 51.11% among the
remaining 180 cells. The equal-source fraction changes from 43.73% to 51.84%.
Neither weighting scheme is promoted to a mouse-level estimand while pool
membership remains unknown. Hopx groups have only two source aliases, so each
omission leaves a single source.

[All source counts](../tables/robustness_2026-09-25/hpcs/source_annotation_counts.tsv),
[group summaries](../tables/robustness_2026-09-25/hpcs/group_influence_summary.tsv),
and [every omission](../tables/robustness_2026-09-25/hpcs/leave_one_source_out.tsv)
include zero categories and explicit denominators.

## Annotation robustness and its boundary

The pinned author definition notebooks are 69,091,594 and 49,691,834 bytes,
above the existing 20-MiB document ceiling. The acquisition preflight stopped
before downloading either. The [source inventory](robustness_source_inventory.json)
records this hold; it does not claim that mappings are absent from those files.

The [eligibility refinement](../config/robustness_annotation_scope.md), frozen
before numerical analysis, therefore retains raw deposited categories. It
compares partitions using label-invariant adjusted Rand index (ARI), complete
source/group contingencies, and separately labelled literal cluster-code changes.
No majority-vote biological mapping was learned from these same cells.
The stringent field's `other` category remains explicit; it is not discarded.

| Partition pair, all 5,333 cells | ARI |
|---|---:|
| `cell type` / `clusterK12` | 0.6331 |
| `cell type` / `clusterK12_stringent` | 0.5723 |
| `clusterK12` / `clusterK12_stringent` | 0.7148 |

The two cluster columns have different literal codes in 1,282/5,333 cells
(24.04%). Literal disagreement is not a validated biological misclassification
rate. Different granularity and related source annotations also affect ARI;
these metrics are not accuracy, fate evidence or independent validation.
Semantic recoding of cluster categories into an alternative HPCS fraction stays
held. [All 87 partition comparisons](../tables/robustness_2026-09-25/hpcs/partition_agreement.tsv)
and [joint count tables](../tables/robustness_2026-09-25/hpcs/annotation_contingencies.tsv)
allow the unresolved correspondence to be inspected.

## Alternative transcript starts

The exact cached native-CHM13 RefSeq Liftoff v5.3 annotation supplies 150 eligible
parent-linked transcript features at 46 distinct gene/TSS positions for the 23
frozen loci. Eleven genes have alternative positions; the other twelve have
only the baseline in this annotation, so they do not provide an alternative-TSS
test. All transcript IDs, zero-based starts, strands and parent checks are in
the [TSS manifest](../tables/robustness_2026-09-25/histone/tss_manifest.tsv) and
[transcript eligibility table](../tables/robustness_2026-09-25/histone/transcript_eligibility.tsv).
No MANE preference or dominant isoform was inferred, and no coordinates were
transferred from another genome assembly.

The 24 cached tracks were quantified at both ±1-kb and ±5-kb windows, producing
2,208 signal rows and 1,104 contrasts. The original matched-H3, positive-signal
and ≥80% coverage rules remain unchanged. Every baseline value reproduces the
original tables; all raw mark, H3 and adjusted effects remain available.

| Mark | Direction changes / eligible alternative-TSS comparisons, ±1 kb | ±5 kb |
|---|---:|---:|
| H3K27ac | 13 / 87 | 18 / 92 |
| H3K4me3 | 7 / 75 | 10 / 84 |
| H3K27me3 | 18 / 90 | 20 / 92 |

Each comparison is an alternative TSS against its baseline for a particular
gene, preparation and state contrast. These overlapping windows and repeated
preparations are correlated, not independent observations. Some sign changes
are near zero; no sign-change threshold becomes a significance test.

For CDKN1A, the two plus-strand starts are chr6:36,497,107 and chr6:36,499,358
(zero-based CHM13 coordinates, 2,251 bp apart). At ±1 kb, iATCs minus iAT2:

| Mark | CUT1 baseline → alternative | CUT2 baseline → alternative |
|---|---:|---:|
| H3K27ac/H3 | +0.837 → −0.654 | +1.273 → +0.860 |
| H3K4me3/H3 | +0.019 → −0.376 | +0.795 → +0.724 |
| H3K27me3/H3 | −1.383 → −2.566 | −0.888 → −1.783 |

The broad-window acetylation effects already disagree between preparations and
retain that disagreement at the alternative start. The defensible follow-up is
therefore a promoter-resolved candidate, not gene-wide CDKN1A activation.
KRT8 and other loci also have reversals; they are retained in the
[complete contrasts](../tables/robustness_2026-09-25/histone/contrasts.tsv) and
[sensitivity summary](../tables/robustness_2026-09-25/histone/sensitivity_summary.tsv).

## Verification and next decisions

Four hand-calculated tests cover unequal-source weighting, invalid/zero
denominators, label-invariant agreement and parent/strand-aware TSS coordinates.
The [independent verification](robustness_verification.json) reconstructed
source counts and all 704 omission rows, checked all 87 ARIs against scikit-learn,
reproduced all 1,104 original histone windows, checked every new signal/control
and contrast, and compared 264 alternative windows with base-resolution reader
values. Prior output hashes remain unchanged. Both new PNGs were visually
inspected; see the [gallery](../figures/README.md).
The pinned BigWig API emitted a deprecation notice for its `missing` argument;
verification still completed successfully.

The [repository validation record](robustness_validation_commands.json) records
Python compilation, 28 lightweight tests (one expected scientific-runtime skip),
all 15 A1 tests in the scientific environment without skips, 18 numeric claim
bindings, Nb1 provenance checks and repository validation. These software and
arithmetic checks do not add biological replication.

The [executed-source archive](execution_sources/robustness_2026-09-25/manifest.json)
preserves the code, helper, frozen scope and annotation eligibility refinement.
The new [run records](hpcs_robustness_run.json) and
[histone run](histone_tss_robustness_run.json) identify input/output hashes.

The next biological comparisons require explicit identity and design gates;
repeated robustness calculations will not supply independent animals or fate
measurements. Update priorities are in
[REMAINING_ANALYSIS_OPTIONS.md](REMAINING_ANALYSIS_OPTIONS.md). The Notion page
intentionally contains only a concise reasoning summary, without this analysis
menu or a critique section.
