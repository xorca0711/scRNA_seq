# A2 analysis plan: does removing the epithelial ligand move the fibroblast

**Superseded in substance, 26 September 2026.** The staging below stands, but
[`config/a2_stage2_freeze_v2.json`](config/a2_stage2_freeze_v2.json) and
[reports/STAGE2_WITHDRAWN.md](reports/STAGE2_WITHDRAWN.md) are the authority on the
statistic, the reference set, the control wells, the endpoint order and the
discriminating set. In particular: there is no p-value and no alpha, the primary
endpoint is the five-gene activation score rather than the Hallmark set, the control
set is TIGIT plus TDTOMATO, Erbb4 is dropped, and leg 2 is exploratory. The text below
is preserved as declared.

Declared 26 September 2026. **Stage 1 only is authorized at declaration.** The
biology is in [RATIONALE.md](RATIONALE.md), the machine-readable decisions in
[`config/a2_delivery_contract.json`](config/a2_delivery_contract.json), and the
question in the [A2 card](../../RESEARCH_QUESTIONS.md#a2). This plan inherits the
prohibitions the [dataset gate](../../docs/NEXT_DATASET_GATE.md) and
[A10 stage 1](../A10_organoid_growth_outcome/reports/STAGE1_IDENTITY_AUDIT.md)
placed on GSE307112, and the measurement contracts
[MC2 to MC4](../../docs/RQ_MEASUREMENT_CONTRACTS.md#mc2).

## Structure

Six stages. Each can stop the work, and stopping is a result.

| Stage | Work | Reads | Can it stop the work |
|---|---|---|---|
| 0 | Frame the hypothesis, rewrite the card, freeze this contract | nothing | Yes. The owner may reject the reframing |
| 1 | Instrument, eligibility and power audit | screen metadata, existing A10 tables, mouse transcript rows for two unvalidated receptors | Yes, on coverage, layout or validation failure |
| 2 | Freeze the test | nothing new | No, but it fixes what stage 3 may do |
| 3 | Leg 1, the within-plate knockout contrast | human endpoint genes from the count workbook | Yes, on a failed totals check |
| 4 | Leg 2, the recipient-licensing correlation | GSE136831 counts and metadata | Yes, by its own depth rule |
| 5 | Report, and propose claim wording without grading it | own outputs | No |

Legs 1 and 2 test different predictions in different data and are reported
separately. Neither validates the other, and a positive leg 1 with a refused
leg 2 is a legitimate outcome that must be reported as such.

## Stage 0: the framing, which is what this commit contains

The card is rewritten around delivery rather than abundance, and the closed
abundance record is stated in the rationale with its claim identifiers. Nothing
is computed. The owner retains or rejects the reframing before stage 1 runs.
Rejection is recorded in DEVELOPMENT.md as a rejection, not as a revision.

## Stage 1: instrument, eligibility and power, before any endpoint exists

The bounded first deliverable. It answers whether the design can carry the test,
and it must not compute either endpoint on any well. Each item is established
from files, not from label structure or from this plan's own assertions.

1. **Plate-3 layout, recomputed.** Per plate-replicate unit, the wells for Areg,
   Egfr, Erbb2, Erbb3, Erbb4, Itgb6 and the control, and the unit total. The
   rationale states one well per axis target per unit, two control wells per unit
   and 60 wells per unit; stage 1 confirms or contradicts that from the deposit.
2. **Knockout validation for the axis.** Report the recorded mouse-transcript
   reduction for Areg, Egfr, Erbb3, Itgb6 and Hbegf from the existing A10 table.
   Erbb2 and Erbb4 have no recorded row, because A10 tested only genes that were
   also module members, so stage 1 extracts those two gene rows from the mouse
   sheet and reports them on the same footing. A target whose transcript does not
   fall is treated as unvalidated, following the recorded case of a target whose
   transcript rose.
3. **Endpoint coverage in the human index.** Verify both endpoint gene sets
   against the workbook's human gene index. Declared floors: at least 0.90 of the
   Hallmark TGF-beta signalling members, and all five members of the repository's
   fibroblast activation score. Report the missing symbols by name.
4. **Covariate audit.** For plate-3 wells, the distribution of fibroblast total
   counts and of epithelial read fraction, and where each axis well sits in it.
   This is the C51 precaution applied before the endpoint exists: a well that is
   extreme in depth cannot carry a programme reading.
5. **Precedent check.** Confirm from the repository that no per-target fibroblast
   endpoint contrast has been computed or reported in any prior analysis, and
   record the files inspected with their hashes. If such a contrast already
   exists, the test is not blind and stage 2 must be re-specified around what has
   already been seen.
6. **Power, computed before any endpoint value is read.** For the declared rank
   statistic with four units, the smallest attainable one-sided p, and the
   within-unit rank pattern needed to reach the declared threshold. This is
   reported whatever it shows, so that a later reader cannot mistake a weak design
   for a strong one.
7. **The cross-species assumption.** Attempt a sequence-level comparison of the
   mouse and human amphiregulin EGF-like domains. If it cannot be performed from
   local resources, record the assumption as carried and unverified rather than
   silently assumed.

**Stop rules.** Four conditions can stop the work, and the seven items above are
reported whatever they show. If the plate-3 layout differs from the record, stop and
re-specify. If endpoint coverage fails its floor, the endpoint is not replaced by
an unfrozen alternative; either a pre-declared alternative applies or the work
stops. If the Areg transcript reduction is absent or positive, leg 1 is not run,
because a perturbation that did not perturb cannot test necessity. If the
precedent check finds a prior per-target contrast, stage 2 is re-specified in the
open. If items 1 to 6 pass, every result downstream is still labelled a
within-screen descriptive association while preparation independence is
unresolved.

## Stage 2: what gets frozen, and why

Nothing in this stage reads data. It fixes what stage 3 is allowed to do.

**Unit.** The plate-replicate unit. A well is not a biological replicate, a guide
is not a biological replicate and an organoid is not a biological replicate. Four
units on plate 3 carry leg 1.

**Endpoints, both defined outside this screen and before it.**

- **Primary:** the human Hallmark TGF-beta signalling set, scored as mean log2
  CPM over its members present in the human index, with the well's own human
  total as the denominator. It is mechanism-matched: the recorded action of
  amphiregulin on mesenchymal cells is TGF-beta activation.
- **Co-primary confirmatory:** the repository's own five-gene fibroblast
  activation score, COL1A1, ACTA2, POSTN, CTHRC1 and TNC, scored the same way.
  It was frozen for trial E6 in a different dataset, it is the score behind C50,
  and it carries the C51 depth precedent with it.

Neither endpoint may be modified, re-weighted or swapped after any value is seen.
A block of Hallmark sets is not used, because one declared endpoint per mechanism
keeps the multiplicity family at two.

**Adjustment, declared because C51 made it load-bearing.** Within each unit,
regress the endpoint on log fibroblast total counts and on epithelial read
fraction across that unit's wells, and use the residual. The unadjusted score is
reported as a sensitivity, never as the primary.

**Statistic.** Within each unit, the rank of the Areg well's residual among all
wells of that unit. The four ranks are combined by their sum, with an exact
one-sided null from the uniform rank distribution. The other perturbed wells in
the unit supply the empirical reference, which is conservative when they act on
the same pathway, and the two control wells anchor the unperturbed location.

**Direction, declared in advance.** Lower. Removing delivered ligand lowers the
recipient's TGF-beta response. A shift in the opposite direction is reported and
does not support the hypothesis.

**Discriminating contrasts, reported with no decision attached.** The same
statistic for the Egfr, Erbb2, Erbb3, Erbb4 and Itgb6 wells, and for the Hbegf
well on plate 4. The Hbegf comparison crosses plates and is labelled weaker for
that reason. These contrasts are how the rivals are read, and they are
interpretive because none of them has its own pre-declared threshold.

**Decision rules.**

- **Supported:** the primary endpoint clears its declared threshold, the
  co-primary agrees in direction, and the receptor contrasts do not show the same
  shift. The reading is a within-screen descriptive association, not causation in
  a lung.
- **Inconclusive:** any other outcome, including a primary that clears while the
  co-primary contradicts it, and including a shift that the covariate sensitivity
  removes.
- **Precise absence:** **unavailable at this design.** Four units and a partial
  knockout cannot bound the effect, and no bound will be reported as absence.
- **Descriptive only:** in every case, while preparation independence is
  unresolved.

## Stage 3: leg 1, the within-plate contrast

1. Extract the endpoint genes from the human sheet of the count workbook,
   streaming it once, and verify the per-library totals against the A10
   extraction record. The script refuses to run if the totals disagree, which is
   the same fail-closed check A10 used.
2. Compute both endpoints per plate-3 well, residualize within unit, rank, and
   combine.
3. Write tables, a run record with the input hashes, and a report. The script
   refuses to overwrite an existing output; a superseded attempt is preserved in
   place, following the A5 and A11 precedent.
4. Report the covariate sensitivity, the unadjusted statistic and the per-unit
   ranks alongside the combined result, so that a single number cannot stand
   alone.

## Stage 4: leg 2, recipient licensing in an existing cohort

**Instrument reused unchanged from trial E6:** GSE136831, donors restricted to
IPF and Control, fibroblasts as Fibroblast and Myofibroblast with mural cells
excluded, a 50-cell floor per compartment per donor, a 10-donor minimum, median
deposited genes per cell as the depth measure, and the same five-gene activation
score. The donor is the unit.

**New estimand.** Does the fibroblast activation score track the recipient's own
TGF-beta activation machinery, measured in the same fibroblasts as the detection
fractions of ITGAV, ITGB1, ITGB8, LTBP1, TGFB1 and THBS1, more closely than it
tracks epithelial AREG detection, which C50 already found it does not track.

**Controls that are not optional.**

- The C51 depth rule applies unchanged: if both members of a pair correlate with
  their own compartment depth at absolute rho of 0.4 or more, the pair is
  reported as confounded and is not read, whatever its own correlation was.
- A frozen negative-control gene set in the same compartment, declared in the
  contract, guards against the artefact that any two gene sets measured in the
  same cells correlate through shared detection.
- The IPF-only stratum is reported next to the pooled result, because a
  correlation present only when controls are pooled in is consistent with both
  variables tracking disease.

**Declared reading.** If recipient machinery predicts activation where ligand
abundance did not, and the depth rule does not refuse the pair, the hypothesis
gains observational support with no direction and no proximity. If neither
predicts, the donor level cannot see the effect at this sample size, and the
report states the smallest effect it could have detected. If both predict
equally, the abundance rival is not excluded. If the depth rule refuses the pair,
that refusal is the result and is reported as one.

## Stage 5: report and register

The report states what is established, what is not, and which rival each contrast
excluded. Claim wording is proposed for the register; grading is the owner's
decision and no row is added by the assistant. If leg 1 is inconclusive and leg 2
is refused by its own control, that pair of outcomes is the finding, and it
belongs in NEGATIVE_RESULTS.md through the usual generation path.

## What this plan refuses to do

- No claim about proximity, distance or gradient from a single-well co-culture.
- No secreted-protein or receptor-engagement claim; C36 stays Not established
  until protein or perturbation evidence on the fibroblast side exists.
- No reading of a receptor-knockout shift as fibroblast causation.
- No wells, guides or organoids as biological replicates.
- No endpoint defined from this screen's own expression and then tested on it.
- No re-grading of C37, C49 or C50, and no new claim rows.
- No study note on the screen's source paper, whose reading is pending.
- No reinterpretation of A10's fits, and no use of A10's outcome as evidence here.
- No substitution of an alternative depth budget or endpoint for the declared one
  as though it were an independent confirmation.

## Order of work

1. This plan, its rationale, the contract and the rewritten card are committed
   before stage 1 runs.
2. Stage 1 runs on the cached screen metadata, the existing A10 tables and the two
   receptor gene rows, and its tables are committed.
3. The owner retains or revises the plan in light of what stage 1 found.
4. Only then is stage 2 committed and stage 3 run. Leg 2 may run in parallel with
   stage 3 once stage 2 is frozen, because it reads a different dataset and
   answers a different prediction.
