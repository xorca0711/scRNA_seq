# Shared epithelial component contract (A5 and A11)

**Status: stage 1 complete, stage 2 rule pre-registered.** The owner authorized
execution. Read the [prospective plan](PLAN.md), the
[machine-readable specification](config/shared_component.json) and the
[stage 1 source audit](reports/STAGE1_SOURCE_AUDIT.md). No expression score has
been computed and no claim grade changes.

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

## Stage 1 outcome

The missing developmental input is now sourced. The existing
[specificity module](../../Research%20Article/epithelial_state_specificity/README.md)
recorded that no independent developmental signature was available locally. Seven
external candidates were examined and one selected: the author-defined signature
of a mixed type 1 and type 2 population in normal mouse lung at postnatal day 1,
from Guo et al. 2019. The audit records why the other six were not used, and that
the selection is partly by accessibility.

## Layout

| Path | Contents |
|---|---|
| [PLAN.md](PLAN.md) | Prospective plan, stages, species and unit rules, prohibited readings |
| [config/shared_component.json](config/shared_component.json) | Reviewable specification: sources, exclusions, partition rule, gates |
| `config/README.md` | What the specification fixes and what stays open |
| [reports/STAGE1_SOURCE_AUDIT.md](reports/STAGE1_SOURCE_AUDIT.md) | Every candidate examined, the doublet check and the limitations |
| `sources/README.md` | Identity of the untracked cached source spreadsheet |

Scripts and tables arrive with stage 2 and stage 3 and follow the
question-specific layout in the
[structure contract](../../docs/REPOSITORY_STRUCTURE.md).

## Evidence this draft rests on

Every number below is read from existing tracked evidence, not recomputed here.

| Fact | Value | Source |
|---|---|---|
| Assayed source fraction gate | 0.7 | Completed human specificity trial |
| Complete-unit floor per paired contrast | 3 | Same trial |
| Lists already failing the fraction gate | 0.641 and 0.667 | Same trial |
| Lists already passing it | 0.780 and 0.791 | Same trial |
| ADI holdout ortholog fraction, which a mapping-only gate would pass | 0.859 | Same trial |
| A11 paired patients, lesion contrast | 23 | Same trial |
| A5 animals passing both group floors | 1 of 25 | Specificity module |
| Local human developmental lung deposits | none found | Scan of local GEO family records |
