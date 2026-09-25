# Shared epithelial component contract: prospective plan

Current continuation: [biological logic and amended child tests](BIOLOGICAL_LOGIC.md). Unit-readiness statements below describe the original contract stage.

25 September 2026. **Complete; the owner retained the frozen modules at stage 4.**
The owner authorized execution of this plan. No expression score has been computed and
no claim grade changes. The machine-readable specification is
[`config/shared_component.json`](config/shared_component.json). Results are in the
[source audit](reports/STAGE1_SOURCE_AUDIT.md) and the
[stages 2 and 3 report](reports/STAGE2_3_REPORT.md). The plan text below is kept as
written before execution, with dated notes where a result changed its reading.

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

**Note after stage 2.** At the level of the three source lists there is no single
shared object. Development shares 5 genes with adult injury, adult injury shares 7
with lesions, and development shares none with lesions. A5 needs the first pair and
A11 the second, and the two are disjoint. The partition still gives both questions
one frozen basis, and the per-gene flags route each to its pair. The phrase "one
shared baseline" overstates what the lists support; see the
[stages 2 and 3 report](reports/STAGE2_3_REPORT.md#consequence-for-the-joint-framing).

## What the contract decides, and what it must not

It decides module eligibility only:

- whether each module clears the assayed source fraction gate and the
  complete-unit floor in each species arm;
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

**Outcome.** Seven candidates were examined; one was selected. The list is the
author-defined signature of a mixed type 1 and type 2 population in normal mouse
lung at postnatal day 1, from Guo et al. 2019. A doublet check against the
authors' own subtype lists found 54 of its 100 genes unique to it. Negretti 2021
would likely be stronger but is not openly retrievable, so the selection is partly
by accessibility. The [source audit](reports/STAGE1_SOURCE_AUDIT.md) records all
seven, the limitations, and a biological rival to A5 found along the way.

### Stage 2: freeze three disjoint modules

**Pre-registered rule, committed before any intersection was computed.** After
removing the operational genes from all three source lists identically:

| Module | Membership |
|---|---|
| Shared remodelling | genes in at least two of the developmental, injury and lesion lists |
| Development-specific | genes in the developmental list only |
| Lesion-specific | genes in the lesion list only |
| Injury residual | genes in the injury list only; a reference, not a contract module |

Shared means recurring in more than one independently defined context. The rule is
symmetric and privileges no pair, and it was fixed before its module sizes were
known. Every gene carries flags for all three lists, so A5 can inspect the part
shared between development and injury as a declared secondary. Genes in all three
lists are tagged, not split out. The type 1 and type 2 identity lists are not
partition inputs; each module's overlap with them is reported instead, which tests
the rival that a shared component is mere type 1 identity.

Two changes from the draft, both made before any intersection was computed. The
lesion module is source-defined rather than discovered in a patient split, because
a patient split is part of A11's own validation and running it here would begin
A11's fit. And Slc4a11, the high-plasticity selection marker, joins the operational
exclusions to match the precedent run.

**Reuse means reading frozen output.** The existing freeze script refuses to run
once its output exists, by design. This contract therefore reads the frozen
`modules.json` and the frozen high-plasticity specification rather than rerunning
or editing the script that produced them.

Exclude the group-defining transcripts consistently across genotypes, as the
current specificity pass already does. A component must not contain the genes
that selected its cells. Score the three Hallmark control axes separately; they
are rivals to be survived, never members of a module.

### Stage 3: report coverage before scoring anything

**Correction recorded during execution.** The first draft of this stage described
a single gate on ortholog mapping. Inspecting the completed human specificity run
showed the 0.7 gate applies to the **assayed source gene fraction**, which loses
genes at two points, and that a separate unit floor also applies. The precedent
values quoted in the draft were already assayed fractions; only the description
of what they measured was wrong. The specification now records both gates.

**Gate 1, assayed source fraction.** Count source genes that both map to a strict
one-to-one ortholog and appear in the target dataset's gene index. Divide by the
original source list length, never by the mapped subset. Threshold 0.7. This is
dataset-specific, so stage 3 reads each target dataset's gene index but not its
counts, and still computes no score.

**Gate 2, complete-unit floor.** At least three complete biological units must
contribute both sides of a paired contrast. A module can clear gate 1 and fail
here. The floor is an eligibility rule, not a power justification.

The distinction between the two loss stages is not cosmetic. Precedent from the
completed run:

| Source list | Ortholog fraction | Assayed fraction | Verdict |
|---|--:|--:|---|
| ADI published holdout | 0.859 | 0.641 | ineligible |
| Six-gene DATP/PATS holdout | 0.667 | 0.667 | ineligible |
| HPCS author top 100 | 0.790 | 0.780 | eligible |
| HPCS without ADI or operational markers | 0.802 | 0.791 | eligible |
| AT2 published holdout | 0.889 | 0.864 | eligible |
| AT1 published holdout | 0.882 | 0.842 | eligible |

The ADI holdout would pass a mapping-only gate and is correctly rejected once
assay presence is applied. A mapping-only report would have admitted it.

Reuse the frozen strict one-to-one ortholog table from the completed run rather
than building a second mapping. A module that fails in one arm is recorded as
ineligible for that arm. Do not substitute a mouse gene for a missing human
ortholog, do not impute, and do not quietly retreat to the cross-species
intersection. Missing source genes are not zeros.

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
