# Shared epithelial component contract (A5 and A11)

**Status: complete. The owner retained the frozen modules on 25 September 2026**
(DEVELOPMENT decision 37). They are now the fixed inputs for A5 and A11, whose own
work continues in [A5's folder](../A5_developmental_programme_reuse/README.md) and
[A11's folder](../A11_lesion_programme_addition/README.md). Start with the
[stages 2 and 3 report](reports/STAGE2_3_REPORT.md), then the
[stage 1 source audit](reports/STAGE1_SOURCE_AUDIT.md). No expression score was
computed here and no claim grade changes.

Start with [the biological logic and revised analysis sequence](BIOLOGICAL_LOGIC.md).

- **A11:** eight eligible Kim patient pairs; the original 23 are discovery.
  Its amended plan separates lesion association from relative activation.
- **A5:** an external Guo gene set avoids Strunz-derived filtering for the primary
  test. The old 94/51-gene variants remain descriptive in Strunz.
- **Shared:** 5 development–injury and 7 injury–lesion genes, with no three-way
  intersection. Biological unrelatedness is not established by limited overlap.

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
| [reports/STAGE2_3_REPORT.md](reports/STAGE2_3_REPORT.md) | Frozen modules, eligibility per question, limitations and stage 4 decisions |
| [scripts/README.md](scripts/README.md) | How to rerun the freeze and the coverage report |
| `tables/` | Frozen modules, membership, overlap, coverage, unit floors and run records |
| [tables/stage3_attempt1_refused/](tables/stage3_attempt1_refused/README.md) | The first coverage attempt, refused by its own precedent check |
| `sources/README.md` | Identity of the untracked cached source spreadsheet |

The layout follows the question-specific pattern in the
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
