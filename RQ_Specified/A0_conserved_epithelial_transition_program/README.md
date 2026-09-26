# A0 — Conserved epithelial transition programme

**Scientific pilot complete, 26 September 2026: primary intestinal transfer not
supported.** Discovery in 9 repair mice and 11 developmental donors produced a
frozen 50-gene programme. In 3 intestinal mice, the intermediate exceeded stem
cells but failed the mature-enterocyte comparison. Specificity analyses were
pruned by the rule declared before transfer scoring.

Read the [scientific result](reports/PILOT_V1_RESULTS.md) and
[figure gallery](figures/pilot_v1/README.md). This completed pilot does not resolve
universality or demonstrate causal modulation of fate.

## Biological question and sequence

Could epithelia reuse some of the cellular work needed to move between identities,
and could that shared process help determine their fate? The observable RNA test
comes first: define a biologically supported intermediate, compare it with both
endpoints within each mouse/donor, discover one common lung signature, then freeze
it before testing another tissue. Only a positive transfer warrants the planned
specificity tests. Causal fate modulation would require additional linked
perturbation/fate measurements.

The [A0 question card](../../RESEARCH_QUESTIONS.md#a0) places this pilot alongside
A1, A5, A11 and A10. Those analyses motivate it without supplying its transfer or
causal endpoint. No independent cells or technical captures are counted as animals.

| Evidence | Purpose |
|---|---|
| [Pilot result](reports/PILOT_V1_RESULTS.md) | Biological rationale, all primary effects, adaptive stop and limits |
| [Source recovery](reports/SOURCE_RECOVERY.md) | Prior analyses, rejected candidates, donor identities and selected branches |
| [Frozen configuration](config/pilot_v1.json) | Fixed cohorts, thresholds, gene exclusions and planned controls |
| [Transfer execution](TRANSFER_EXECUTION.md) | Score definition and conditional stop declared before V1 effects |
| [Exact 50-gene programme](tables/pilot_v1/frozen_programme.json) | Membership committed before transfer |
| [Figures](figures/pilot_v1/README.md) | Coverage, discovery gates, gene effects and frozen transfer |
| [Current status](tables/pilot_v1/pilot_status.json) | Executed stages and explicitly pruned P4 |
| [Reproduction](REPRODUCING.md) | Inputs, scripts, safeguards and verification |

## Historical feasibility audit

The [original report](reports/PILOT_REPORT.md), `readiness.json`, `validation.json`
and `execution_record.json` describe the 25 September audit, when no programme
had been learned. Their blocked status is historical. Original tables are
unchanged; the [documentation archive](history/feasibility_documentation.zip)
preserves the original bytes of updated entry points. The
[continuation plan](CONTINUATION_PLAN.md) records the subsequent input recovery.
The original [plan](PLAN.md) and [candidate list](DATASET_CANDIDATES.md) remain
available with repaired links to the canonical `Research Article/` directory.

## Separate repair-first exploratory analysis

This is a separate repair-first exploratory run, executed from the earlier
checkout and published beside the already merged [pilot_v1](reports/PILOT_V1_RESULTS.md)
from PR #81. The source data overlap; these are not independent replications.
Unlike pilot_v1's two-context discovery and three-mouse intestinal branch, this
run discovers only in repair, describes mouse developmental SCT scores, and tests
the proximal intestinal branch with two qualifying mice. Its program, scoring,
controls and stage decisions are separate; it does not reopen pilot_v1's pruned P4.
The earlier replication gaps describe this run's chosen cohorts, not all A0 work.

Read its [report](reports/EXPLORATORY_PILOT_REPORT.md), [stage decisions](stage_decisions.json), [status](exploratory_readiness.json) and [reproduction instructions](EXPLORATORY_REPRODUCING.md). The full 50-gene program remains near zero in the selected intestinal comparison; an initially positive 31-gene control variant weakens after depth matching. The decision is to narrow the interpretation and stop broader expansion of this signature.
