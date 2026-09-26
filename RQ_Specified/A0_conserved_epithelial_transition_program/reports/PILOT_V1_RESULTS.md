# A0 scientific pilot: lung discovery passes; intestinal transfer does not

Completed 26 September 2026. **The bounded scientific pilot is complete with a
negative primary transfer decision.** A 50-gene programme satisfies the fixed
repair/development discovery rules, but does not distinguish the intestinal
intermediate from both endpoints. This is not a finding that all conserved
transition processes are absent, and it does not test causal modulation of fate.

## Why this is the relevant test

An intermediate can differ from its starting population simply because it has
begun acquiring a destination identity. A general intermediate-associated
programme should therefore exceed **both** reference states, first within each
lung context and then with unchanged membership in another tissue. This prevents
a positive stem-to-progenitor difference from being mistaken for a conserved
transient peak.

The existing A5 result supports recruitment of an external developmental
signature in repair; A5/A11 overlap is pairwise, not a universal core. A1's
regulatory/outcome evidence and A10's measured growth do not provide this
cross-tissue test or a common fate endpoint. Their reports were checked before
launch. The [recovery report](SOURCE_RECOVERY.md) connects these precedents to
the source decisions and biological rationale.

## Cohorts and prospective decisions

| Role | Operational biological comparison | Independent units | Primary cells |
|---|---|---:|---:|
| D1: [Strunz, GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259) | Injured AT2 → Krt8+ ADI → AT1; days 10–15, source normal controls excluded | 9 mice | 3,463 |
| D2: [Sountoulidis 2023](https://doi.org/10.1038/s41556-022-01064-x) | Distal epi_cl2 → intermediate epi_cl1 → proximal-secretory epi_cl0 | 11 donors | 3,560 |
| V1: [Haber 2017](https://doi.org/10.1038/nature24489) | Stem → late enterocyte progenitor → mature proximal enterocyte | 3 mice | 1,013 |

Arrows name source-supported developmental comparisons; they do not trace each
sequenced cell. D2 is an early human airway axis, not mouse alveolar maturation.
Spatial evidence supports intermediate position; inferred trajectories do not
establish individual ancestry. Haber late progenitors are not asserted to be
exclusively proximal. These limitations constrain the tested prediction.

The original Negretti/immature-enterocyte proposal remained below its eligibility
requirements. The replacements were chosen from source biology and coverage
**before programme effects**, retaining the 30-cell/state and three-unit floors.
The original negative feasibility tables are unchanged. D1 also retains
387 activated-AT2 cells as a possible control, outside the primary three-state count.
They were scored, but no activated-AT2 specificity contrast was evaluated.

The source-QC populations were retained at positive integer library totals. This
pilot does not independently redo ambient-RNA/doublet calling. The D2 export was
checked against 255,340 raw GEO values at ten fixed matched barcodes; all agreed.
Every retained cell's persisted count total was independently recomputed. D2
library totals exactly match deposited metadata for all 3,560 selected cells.
A 39-capture crosswalk maps to 17 source donors; captures do not count as donors.

- `b544cbc`: exact cohorts, exclusions, ortholog rules and thresholds before effects.
- `d079150` / `e5ae0b6`: discovery code and compact output before discovery;
  `5e1576f` / `04261f5`: storage/parsing improvements only.
- `f1a1df7`: the selected 50-gene instrument committed before V1 scores.
- `e3e12a0`: transfer code and conditional specificity stop committed before V1 scores.

The [fixed configuration](../config/pilot_v1.json),
[transfer execution details](../TRANSFER_EXECUTION.md) and input hashes are the
prospective record. The term prospective applies to this continuation; previously
inspected Strunz/A5 evidence is not untouched validation.

## Discovery

The intersection contains **11,590 conservative one-to-one mouse–human
orthologs**. Library totals use all deposited genes in each assay. For each unit
and state, sum counts, calculate log2(CPM + 1), then compare the intermediate
with each endpoint. Each of the four contrasts must have a median at least 0.1
and a positive direction in at least two-thirds of units; intermediate CPM must
be at least 1 in at least two-thirds of units in both discovery contexts.

After the fixed annotation-marker and symbol exclusions, **134 genes qualify**.
Rank by the smallest of the four median effects, then human symbol, and retain
at most 50: the [frozen programme](../tables/pilot_v1/frozen_programme.json)
contains 50 genes. The 20-gene minimum and all thresholds were unchanged.
The complete effect table includes every candidate and unit, not only successes.
Discovery-context score differences are selected estimates, not independent
confirmation of conservation or evidence of fate control.

## Frozen transfer

Scores are the mean within-cell percentile ranks of the fixed genes over the
same ortholog universe, followed by unit/state means. The table reports paired
intermediate-minus-endpoint medians in **percentile-rank score points** (100 times
the stored 0–1 score difference), not percentage changes in expression. Effects
are interpreted within their source; cross-assay magnitude comparisons have no
calibrated biological scale.

| Context | Endpoint reference | n | Median difference | Positive units | Smallest leave-one-out median | Fixed criterion |
|---|---|---:|---:|---:|---:|---|
| D1 | start | 9 | +3.7170 | 9/9 | +3.6460 | Pass |
| D1 | destination | 9 | +2.7349 | 9/9 | +2.6250 | Pass |
| D2 | start | 11 | +3.3167 | 10/11 | +3.2554 | Pass |
| D2 | destination | 11 | +3.7912 | 11/11 | +3.7213 | Pass |
| V1 | start | 3 | +0.2694 | 2/3 | +0.0786 | Pass |
| V1 | destination | 3 | -0.0914 | 1/3 | -0.2929 | Fail |

V1 must have a positive median, at least two-thirds positive mice and a positive
median after **every** single-mouse omission for **both** endpoints. It passes
against stem cells but fails against mature enterocytes. The three
intermediate-minus-mature differences are +0.4065, −0.0914 and −0.4945 score
points for Control-Mouse1, Control-Mouse3 and Control-Mouse4, respectively.
The primary decision is therefore **STOP_PRIMARY_TRANSFER_NOT_SUPPORTED**.
With three mice, these are descriptive directional criteria, not a confirmatory
p-value, an equivalence test or strong evidence of biological absence.

D1's equal-day means remain positive (+3.6794 vs AT2, +2.6700 vs AT1), and every
leave-one-day-out median is positive. Those diagnostics do not repair the failed
new-tissue endpoint. The [unit table](../tables/pilot_v1/transfer_unit_differences.tsv)
and [leave-one-out table](../tables/pilot_v1/transfer_leave_one_out.tsv) retain
all estimates.

## Adaptive stop and interpretation

P0/P1 recovery and freezing, P2 discovery, and P3 frozen transfer are executed.
P4 stress/cycle exclusions, matched random modules, low-cycle strata and depth
rescoring are **pruned**, under the conditional rule declared before V1 scoring.
There is no positive transferred association whose specificity needs establishing.
No alternate gene set, state boundary or score was tried after the failure.
The negative decision is complete for this pilot; P4 is not an unfinished queued job.

What is established is narrower than the motivating universal hypothesis: these
lung settings yield a shared selected RNA signature, but its frozen primary
score does not peak at the specified intestinal intermediate relative to both
references. Generic activation, source labels, dissociation and other explanations
for the lung selection remain unresolved because specificity was not tested.
The negative transfer cannot distinguish tissue specificity from differences in
branch, species, platform, source QC or operational state definitions. Excluding
known markers does not eliminate expression-label circularity. Repeated libraries,
cells and genes do not add independent biological replication.

The broader conservation/fate hypothesis remains open. Testing a smaller module,
a different branch, regulatory activity rather than RNA membership, or causal
fate modulation would require a separately justified design and suitable new
evidence. Those are new investigations, not missing calculations to complete
this fixed pilot.

## Verification and artifacts

- Discovery: 197,037 saved-summary checks and 3,320 independent raw-count effect
  comparisons across 83 probe genes; maximum error 4.99e-12.
- Transfer: 8,423 cell scores reaggregated, 302 arithmetic/hash checks and 78
  independent tied-rank calculations directly from counts; maximum score error
  6.67e-16. No scoring helper was reused for the independent rank calculation.
- All 47 original feasibility artifacts match their recorded hashes, directly
  or through the four-document preservation archive.
- Four figures were visually inspected; source, output and review hashes are
  recorded. Portable repository tests check delivered evidence without raw caches.

Open the [figure gallery](../figures/pilot_v1/README.md),
[reproduction guide](../REPRODUCING.md),
[discovery verification](../tables/pilot_v1/verification.json),
[transfer verification](../tables/pilot_v1/transfer_verification.json), and
[current machine-readable status](../tables/pilot_v1/pilot_status.json).
Raw public sources and working count arrays remain ignored. Two interrupted
preparation attempts are preserved separately; neither calculated programme
effects. Numerical outcomes, selection thresholds and historical claim grades
were not rewritten to improve the result.
