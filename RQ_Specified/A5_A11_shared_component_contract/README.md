# Shared epithelial component contract (A5 and A11)

**Status: draft, nothing executed.** Read the
[prospective plan](PLAN.md) and the
[machine-readable specification](config/shared_component.json). No gene list is
frozen, no score has been computed, no figure exists and no claim grade changes.
This folder contains a specification awaiting the owner's retain or reject.

## What this is

Enabling work owned jointly by [A5](../../RESEARCH_QUESTIONS.md#a5) and
[A11](../../RESEARCH_QUESTIONS.md#a11). Both hypotheses need the same object, a
frozen shared epithelial remodelling component. Defined separately, their two
results describe different baselines and neither constrains the other. Defined
once, a three-way partition becomes statable: shared, development-specific and
lesion-specific.

## What this is not

- **Not a new research question.** It takes no new register identifier. Its
  status matches the enabling source-identity question A12-S1.
- **Not a merge of A5 and A11.** The two keep separate species arms, biological
  units, contexts, multiplicity families and decisions. The register's merge rule
  requires a shared test and a shared decision, which these do not have.
- **Not a result.** Passing this contract licenses two later tests. It
  establishes no mechanism, cell identity, ancestry or fate.

## Current gate

The development-specific module is **blocked on a missing input**. The existing
[specificity module](../../Research%20Article/epithelial_state_specificity/README.md)
records that an independently sourced developmental maturation signature is not
available in its current pass. Sourcing one is stage 1, and it must come from a
study that did not define the injury or lesion lists. Otherwise sharing is built
in by construction and the A5 test is circular.

## Layout

| Path | Contents |
|---|---|
| [PLAN.md](PLAN.md) | Prospective plan, stages, species and unit rules, prohibited readings |
| [config/shared_component.json](config/shared_component.json) | Reviewable specification: modules, coverage gate, unit rules |
| `config/README.md` | What the specification fixes and what stays open |

Scripts, tables, reports and a gallery are absent because nothing has run. They
follow the question-specific layout in the
[structure contract](../../docs/REPOSITORY_STRUCTURE.md) if the owner retains
this draft.

## Evidence this draft rests on

Every number below is read from existing tracked evidence, not recomputed here.

| Fact | Value | Source |
|---|---|---|
| Ortholog coverage gate | 0.7 | Completed human specificity trial |
| Lists already failing that gate | 0.641 and 0.667 | Same trial |
| Lists already passing it | 0.780 and 0.791 | Same trial |
| A11 paired patients, lesion contrast | 23 | Same trial |
| A5 animals passing both group floors | 1 of 25 | Specificity module |
| Local human developmental lung deposits | none found | Scan of local GEO family records |
