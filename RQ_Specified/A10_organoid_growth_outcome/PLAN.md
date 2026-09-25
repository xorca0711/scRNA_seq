# A10 analysis plan: epithelial and niche programmes against measured organoid growth

25 September 2026. **Stage 1 only is authorized to run. Nothing is fitted.** The
biology and the argument are in [RATIONALE.md](RATIONALE.md); the machine-readable
decisions in [`config/a10_outcome_contract.json`](config/a10_outcome_contract.json).
This plan follows the register's [A10 card](../../RESEARCH_QUESTIONS.md#a10) and the
[dataset gate](../../docs/NEXT_DATASET_GATE.md), which set the joins and the
prohibitions this plan inherits.

## Structure

Four stages. Each one can stop the work, and stopping is a result.

| Stage | Work | Reads | Can it stop the work |
|---|---|---|---|
| 1 | Identity and join audit | design, imaging, species QC; no counts | Yes. Unresolved units restrict everything after it to description |
| 2 | Freeze the model and programme set | nothing new | No, but it fixes what stage 3 may do |
| 3 | Nested fits with whole-preparation holdout | counts, 303 MB | Yes, on a failed assay control |
| 4 | Report | own outputs | No |

## Stage 1: the identity and join audit

The bounded first deliverable. It answers whether the design can carry a model at
all, and it is the only stage authorized now.

Questions it must answer, each from the files rather than from label structure:

1. **Do the three tables join on one well key?** Design is keyed by plate and well,
   imaging by plate, replicate, well and day, species QC by library, plate,
   replicate and well. Report unmatched rows in every direction, duplicate keys and
   wells present in one table but not another.
2. **What is the biological unit?** The replicate index is the candidate. Establish
   from evidence whether those indices are independent preparations, that is separate
   type 2 isolations and fibroblast lots, or repeats of one preparation. A label like
   `1-1` does not settle it. If the deposit cannot settle it, say so.
3. **Which species is which.** The species-assignment fields distinguish two
   genomes. Establish which corresponds to the mouse epithelium and which to the
   human fibroblasts from the files, not from convention, because the whole
   epithelial and niche split depends on it.
4. **How much read assignment is ambiguous.** Report the per-library fractions of
   ambiguous, both and neither, since these bound how cleanly the two compartments
   separate.
5. **Per-target replication.** Report libraries per target, and flag the uneven ones.
   Two targets are already known to be extreme: one appears in 87 libraries and one
   in 30, against four for most. Their role must be established, not assumed.
6. **The known asymmetry.** One plate replicate has imaging but no RNA libraries.
   Report it and any other such gap.

**Stop rule for stage 1.** If wells cannot be joined, the work stops. If wells join
but preparations cannot be established, stages 2 and 3 may still run, and every
result is labelled a within-screen descriptive association with no claim of
independent replication. That restriction is the register's, not a new one.

## Stage 2: what gets frozen, and why

**Outcome.** Primary is mean organoid area at day 14, conditional on mean organoid
area at day 7. Organoid count and area proportion are declared secondary. This is
the register's own choice and is not revisited after seeing results.

**Baseline model.** Growth has mundane determinants, and they go in first:

- day-7 mean area, because a bigger start gives a bigger finish;
- plate and plate replicate, for batch;
- **well composition**, from the species read fractions. This one matters. A well
  with relatively more fibroblast material yields different bulk expression on both
  sides, so without it a programme score could simply be reporting how much of each
  cell type was in the well rather than a programme state.

**Programme set.** The frozen modules from the
[shared component contract](../A5_A11_shared_component_contract/README.md), plus the
identity and control axes it already uses. They are reused rather than invented for
two reasons: they are already frozen, and they were defined without any reference to
this screen, so they cannot be tuned to it. The epithelial side is mouse and uses
them directly; the fibroblast side is human and maps through the same frozen strict
one-to-one ortholog table, under the same assayed-fraction gate of 0.7.

**Tests, in a fixed order.** Two block comparisons, not one test per module, so the
multiplicity stays small:

1. **Epithelial increment.** Baseline versus baseline plus the epithelial programme
   block.
2. **Niche increment.** The winner of step 1 versus that plus the fibroblast
   programme block.

Per-module associations are exploratory and reported as such, with their own
multiplicity correction and no decision attached.

## Stage 3: fitting, and how it is evaluated

**Held-out evaluation by whole preparation.** Every fold holds out all wells of one
preparation. Sibling wells share isolation, fibroblast lot, plate handling and
often target, so a random well split would leak and would inflate any increment.
This is the single most important guard in the plan.

**Metric and margin.** Held-out predictive error on the primary outcome, compared
between nested models, with the margin for a useful improvement declared in the
contract before fitting.

**Assay controls.** Targets with published effects in this assay serve as positive
controls. If they show no effect, the assay or the join is wrong and the fits do not
get interpreted. Their recovery is validation and is reported as validation, never
as a finding of this repository.

**Decision rules.**

- **Supported** if the increment exceeds the declared margin and reproduces across
  held-out preparations.
- **Precisely absent** if the confidence bound excludes the declared margin. This
  retires the specified association for this assay, nothing broader.
- **Inconclusive** otherwise, including when preparations are too few.
- **Descriptive only** if stage 1 left the biological unit unresolved.

## What this plan refuses to do

- No claim of forecasting. Expression and outcome are both from day 14.
- No treating wells, guides or organoids as biological replicates.
- No reading organoid area as mature type 1 fate or as repair in a lung.
- No new programme defined from this screen's own expression and then tested on it.
- No use of this screen as a test of A5 or A11. The cross-links in the rationale are
  opportunities for those questions' own plans.
- No study note on the source paper, whose reading by the owner is pending.

## Order of work

1. This plan, its rationale and its contract are committed before stage 1 runs.
2. Stage 1 runs on the three cached metadata files and its tables are committed.
3. The owner retains or revises the plan in light of what stage 1 found.
4. Only then is the 303 MB count table fetched and stages 2 and 3 specified in
   detail. Fetching it earlier would prejudge the gate.
