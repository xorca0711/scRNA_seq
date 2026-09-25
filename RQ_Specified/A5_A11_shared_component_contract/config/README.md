# Specification files

[`shared_component.json`](shared_component.json) is the reviewable contract. It
fixes the decisions that must not move once results are visible, and it leaves
open the parts that require a sourcing job the owner has not yet authorized.

## What the specification fixes now

- The three module identifiers and which question owns each test.
- The label genes excluded from every module, and why.
- The Hallmark control axes, which are scored separately and never merged in.
- The ortholog coverage gate at 0.7, reported per module per species arm, before
  any score is computed.
- The species rule: freeze once from one source definition, then map per arm.
- The unit rules, including the prohibition on transferring units across arms.
- Separate multiplicity families per owning question.

## What stays open

Module membership is `null` in all three modules. A null list is a missing input,
never an empty biological result. Membership is written only by the stage 1 and
stage 2 sourcing and freezing jobs, and only if the owner retains this draft.

The development-specific module carries `membership_state:
blocked_missing_source`. Its blocking input is an independently sourced
developmental maturation signature, which the existing specificity module records
as unavailable in its current pass.

## Reading order

Read [PLAN.md](../PLAN.md) first for the argument, then this specification for
the exact values. Where the two disagree, the specification is authoritative for
values and the plan is authoritative for scope.
