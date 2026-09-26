# A0 exploratory expression pilot — staged amendment

Authorized 2026-09-26: proceed with the agreed exploratory analysis and reassess
whether the next stage is worthwhile after each stage. This supersedes the
earlier requirement to qualify all three contexts before any expression analysis.
This run does not itself fulfill the original three-context design; the separate
pilot_v1 uses different developmental and intestinal cohorts. Previously reported feasibility results are historical.

## E1 — Establish a reproducible repair-associated program

Use the nine GSE141259 enriched epithelial source mice with complete AT2,
Krt8+ ADI and AT1 triplets during days 10–15. Exclude NC-labeled controls.
Retain author QC and labels; check count/annotation alignment, nonnegative integer
counts, nonzero libraries and mitochondrial/UMI distributions. No new cell filters
are chosen from program outcomes. Author annotations are not blind and can retain
circularity even after removal of known label genes.

Sum raw counts by mouse/state, then calculate log2(CPM+1) intermediate-minus-AT2
and intermediate-minus-AT1 effects. Every mouse has equal weight. Eligible genes
must have ≥10% detection in ADI in ≥two-thirds of discovery mice; ≥two-thirds of
paired effects positive for EACH endpoint; and median effects ≥0.25 for EACH
endpoint. Exclude the frozen label-gene list and mitochondrial/ribosomal symbols.
The exclusion is implemented conservatively as mitochondrial and Rpl/Rps symbol
prefixes. Rank by the smaller of the two median effects, then the smaller positive fraction,
then gene symbol. Select at most 50 genes; require at least 20. This discovers a
repair-associated candidate, not a cross-context consensus program.

Score each gene against 1,000 fixed reference genes sampled with seed 0 from the
common feature-identity intersection of repair, development and intestine after
label/mitochondrial/Rpl/Rps-prefix exclusions. This small metadata-only check avoids
selecting reference genes absent from a target assay. Reference sampling uses
feature identity only, before expression outcomes. For each cell,
each program gene gets the fraction of reference genes below it plus half of ties;
the module score is their mean. Omit program/reference overlap from the reference
set for each module. This is a within-cell rank statistic, not a significance test.
Transfer uses the same reference symbols; report feature/detection coverage.

Leave one mouse out: reapply the entire gene filter/ranking on the other eight,
then score the held-out mouse. Report both endpoint differences and differences
divided by the pooled endpoint-cell score SD (descriptive magnitude, not biological
uncertainty). Reused author labels are a limitation of this internal validation.

**Checkpoint E1:** proceed to specificity checks if a ≥20-gene candidate exists,
at least six of nine held-out mice have positive effects against each endpoint,
and both median standardized effects are ≥0.25. Otherwise report a weak or unstable
source signal and stop downstream transfer of this candidate. These are exploratory
investment criteria, not calibrated proof thresholds or effect-size power claims.

## E2 — Determine what the source signal measures

Compare source effects with the existing, provenance-recorded ES1 modules: published
ADI/AT2/AT1 lists, p53, hypoxia, inflammatory/TNF signaling, E2F/G2M and AP1 panels
where available. These are controls/benchmarks, not new discoveries. Score-size
differences prohibit interpreting raw differences across unrelated modules as
relative biological strengths.

Use 200 fixed expression-bin- and size-matched random gene sets, within-mouse depth
matching, source score–depth/mitochondrial correlations, and endpoint synthetic
mixtures as descriptive technical checks. Mixtures and random sets are not animals.
Show a candidate variant excluding genes in the prespecified generic control union
if ≥10 genes remain, without replacing the frozen primary program. Do not regress
away biological stress indiscriminately. Compare activated AT2 when available,
explicitly acknowledging that it is not a proven stress-only negative control.

**Checkpoint E2:** proceed to descriptive transfer if the signal retains both
endpoint directions after depth matching and is not accounted for by endpoint
mixtures. Stress overlap alone does not block transfer: the normal-tissue check
can distinguish an injury-associated response from broader reuse. If technical
effects account for the signal, stop; if stress dominates, narrow the interpretation
and carry that limitation into transfer. Record all equivocal controls.

## E3 — Frozen descriptive transfer

Freeze the full-repair module, reference genes, score definition, and control
variants before any developmental/intestinal expression outcomes are examined.

- Development: author Negretti epithelial annotations; postnatal AT2 → Transitional
  → AT1 branch only. Public normalized SCT values may support a descriptive rank
  check if raw counts cannot be aligned economically. Disclose this limitation.
  Show every capture group and its state counts; the one complete group does not
  provide replication. Do not reinterpret barcode groups as mice.
- Intestine: Haber author Stem → Enterocyte.Immature.Proximal →
  Enterocyte.Mature.Proximal. Use the four source-mapped mice, distinguishing the
  two passing the original 30-cell floor. Any other batch-level display is separate
  and explicitly non-independent. No post-outcome branch merging or gene retuning.

At least 80% program feature coverage is required for interpretable scoring;
report reference coverage separately (require ≥80%). Zero expression is different
from an absent feature. Underfloor groups may have descriptive point estimates,
with counts shown and no cell-based inferential claims. They do not become
qualified biological replicates. Do not pool ages to manufacture replication.

**Checkpoint E3:** consistent positive effects against both endpoints motivate
better replicated validation. Endpoint-specific, absent or inconsistent effects
motivate narrowing or stopping this shared-RNA-program direction, conditional on
assay coverage and normalization limits. Two mice or one developmental group
cannot establish broad conservation; no causal/universal conclusion is available.

## E4 — Closeout

Generate gene/effect/control tables, sample-level figures, a report, hashes and a
machine-readable stage decision log. Keep observations, limitations and investment
recommendations separate. Do not force all stages to run after a failed checkpoint.
No unsupervised atlas integration, new trajectory reconstruction, GRN fitting or
raw-read reprocessing is needed for this pass. No biological p-values will treat
cells as independent replicates. Defaults are fixed in exploratory_config.json.

## Dated post-transfer addition — E3a (2026-09-26)

After E3, the near-zero primary intestinal score and positive pre-transfer
31-gene variant justified a targeted depth sensitivity check before deciding
whether expansion was warranted. Reuse existing frozen primary, variant and
published ADI scores; match equal state counts within pooled log-UMI quintiles
inside each mapped mouse, without replacement. Use seed 0 and 200 draws to show
matching sensitivity. Report matched cell counts and 5–95% draw ranges as
technical sensitivity, never biological confidence intervals. No gene or branch
retuning is allowed. This addition was selected after transfer outcomes and is
not represented as part of the pre-transfer frozen design.

The completed check supports narrowing and closeout; details are in the
[current report](reports/EXPLORATORY_PILOT_REPORT.md) and
[stage decisions](stage_decisions.json).

## Publication context

This is a separate repair-first exploratory run, executed from the earlier
checkout and published beside the already merged [pilot_v1](reports/PILOT_V1_RESULTS.md)
from PR #81. The source data overlap; these are not independent replications.
Unlike pilot_v1's two-context discovery and three-mouse intestinal branch, this
run discovers only in repair, describes mouse developmental SCT scores, and tests
the proximal intestinal branch with two qualifying mice. Its program, scoring,
controls and stage decisions are separate; it does not reopen pilot_v1's pruned P4.
The earlier replication gaps describe this run's chosen cohorts, not all A0 work.
