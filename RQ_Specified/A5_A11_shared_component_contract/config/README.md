# Specification files

[`shared_component.json`](shared_component.json) is the reviewable contract. It
fixes the decisions that must not move once results are visible. The partition
rule inside it was committed before the freeze ran, and the freeze record proves
the freeze read those exact bytes.

## What the specification fixes

- The three module identifiers and which question owns each test.
- The label genes excluded from every module, and why.
- The Hallmark control axes, which are scored separately and never merged in.
- Two eligibility gates, reported per module per species arm per target dataset
  before any score is computed: an assayed source fraction of at least 0.7, where
  genes are lost both at ortholog mapping and at assay presence, and a floor of
  three complete biological units per paired contrast.
- The frozen strict one-to-one ortholog table to reuse, rather than a new mapping.
- The species rule: freeze once from one source definition, then map per arm.
- The unit rules, including the prohibition on transferring units across arms.
- Separate multiplicity families per owning question.

## Where membership lives

Module membership fields are `null` in this file by design. This file holds the
decisions; the freeze script wrote the membership to `../tables/frozen_modules.json`
and refuses to overwrite it. Keeping the two apart means a result cannot quietly
rewrite the rule that produced it.

The developmental source that was missing when this file was drafted is now
recorded under `source_lists`, with its archive and member hashes. The
[source audit](../reports/STAGE1_SOURCE_AUDIT.md) explains the choice.

## Reading order

Read [PLAN.md](../PLAN.md) first for the argument, then this specification for
the exact values. Where the two disagree, the specification is authoritative for
values and the plan is authoritative for scope.
