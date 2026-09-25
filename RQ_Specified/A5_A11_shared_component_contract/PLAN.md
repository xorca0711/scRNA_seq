# Shared epithelial component contract: prospective plan

25 September 2026. **Draft for the owner's retain or reject decision. Nothing
has executed. No gene list is frozen, no score has been computed, and no claim
grade changes.** The machine-readable specification is
[`config/shared_component.json`](config/shared_component.json).

This contract is enabling work owned jointly by [A5](../../RESEARCH_QUESTIONS.md#a5)
and [A11](../../RESEARCH_QUESTIONS.md#a11). It is **not** a new research question
and takes no new identifier in the register, following the rule that A1 states
for itself: link to the existing questions rather than creating another global
register. Its status matches the enabling source-identity question A12-S1.

## Why a shared contract rather than a merged question

A5 and A11 have the same logical shape, a shared component plus a
context-specific addition. The register's merge rule permits merging only when
two hypotheses share a discriminating test and a decision, not merely a topic.
A5 and A11 fail that rule on four axes.

| Axis | A5 | A11 |
|---|---|---|
| Species arm | mouse | human |
| Biological unit | animal | patient |
| Contexts compared | neonatal development, adult injury | lesion, normal, fibrosis, repair |
| Decision | whether development contributes to a shared component | whether lesions add a component beyond it |

Either hypothesis can hold without the other, so they keep separate decisions.

What they genuinely share is the **object both tests need**, a frozen shared
remodelling component. Defined twice, independently, the two results describe
different baselines and neither constrains the other. Defined once, a three-way
partition becomes statable: shared, development-specific, lesion-specific. That
partition is the scientific gain from joining, and it is the whole purpose here.

## What the contract decides, and what it must not

It decides module eligibility only:

- whether each module clears the ortholog coverage gate in each species arm;
- whether the three modules are disjoint after label-gene and overlap exclusion;
- whether an independent developmental maturation signature exists at all.

It decides neither hypothesis. It establishes no mechanism, no cell identity, no
ancestry and no fate. Passing this contract licenses two separate later tests;
it is not itself a result.

## Stages

### Stage 1: source the missing developmental signature

The gating input. The existing specificity module records that an independently
sourced developmental maturation signature is **not available** in its current
pass. Until one exists, the development-specific module cannot be frozen and A5
has no addition to test.

The signature must come from a developmental study that did not define the
injury or lesion lists. If the same source defines both, sharing is built in by
construction and the A5 test becomes circular. Record the source, locator,
retrieval date and hash. If no suitable independent source is found, stop and
report that, rather than reusing an injury list as a developmental proxy.

### Stage 2: freeze three disjoint modules

Extend the existing module freeze rather than writing a new one. Reuse
`freeze_modules.py`, the tracked `modules.json` provenance pattern and the
existing overlap table, whose stated purpose is preventing overlapping
signatures from counting as independent corroboration.

Exclude the group-defining transcripts consistently across genotypes, as the
current specificity pass already does. A component must not contain the genes
that selected its cells. Score the three Hallmark control axes separately; they
are rivals to be survived, never members of a module.

### Stage 3: report coverage before scoring anything

Report retention after strict one-to-one ortholog mapping, per module and per
species arm, against a 0.7 gate. Publish the table before any score is computed,
so the gate cannot move after results are visible.

The precedent from the completed human specificity trial shows this gate bites:

| Source list | Coverage | Verdict |
|---|--:|---|
| ADI holdout | 0.641 | ineligible |
| Six-gene DATP/PATS holdout | 0.667 | ineligible |
| Complete HPCS source list | 0.780 | eligible |
| Overlap-reduced HPCS list | 0.791 | eligible |

A module that fails in one arm is recorded as ineligible for that arm. Do not
substitute a mouse gene for a missing human ortholog, do not impute, and do not
quietly retreat to the cross-species intersection. Missing source genes are not
zeros.

### Stage 4: stop for the owner's decision

No fitting in this contract. A5 and A11 then run their own tests against the
frozen baseline, in their own units, with their own multiplicity families.

## Species, units and multiplicity

**Species.** A scan of every local GEO family record found developmental terms
only in the mouse multiome deposit. There is no local human fetal or
developmental lung. So the developmental arm is mouse and the lesion arm is
human, and each module is frozen once from one source definition then mapped per
arm. Intersecting across species shrinks both arms rather than one.

**Units.** A5 has one of 25 external-study animals passing both group floors,
which permits component definition and forbids confirmation. A11 has 23 paired
patients for the lesion contrast, with the same patient on both sides. That
pairing is what answers the donor, batch and much of the composition rival by
construction, and it is the reason A11 can carry a test that A5 currently cannot.

Sharing a component definition does not transfer units. The human patient count
never lends power to the mouse arm. A joined presentation invites exactly that
error, so it is prohibited explicitly.

**Multiplicity.** Separate families per question. One joint family spanning two
species and four contexts inflates the correction and could suppress both arms.

## Prohibited readings

- Do not order development, repair, fibrosis and neoplasia as one progression.
  The register forbids an assumed repair-to-cancer trajectory, and a contract
  spanning both questions makes that misreading easier, not harder.
- A similar module score across contexts is not common ancestry or transformation.
- A frozen module is not a cell-type or malignancy classifier.
- A shared component is not evidence that no context-specific programme exists.
- The overlap-reduced signature is a scored module only. The A1 robustness batch
  found pooled fractions moving from 31.79 to 51.11 percent when the dominant
  source library is omitted, which qualifies its provenance as an identity.

## Next-layer proposal

Chromatin is the natural next layer for separating a shared component from a
context-specific addition, and two multiome deposits with peak annotations are
already cached locally.

**Feasibility verdict: descriptive only, not a confirmation arm.** Each
condition is one pooled library, so the accessibility data cannot supply
biological units. Accessibility also does not measure histone marks or
methylation. Treat it as a motivating panel whose units are disclosed.
