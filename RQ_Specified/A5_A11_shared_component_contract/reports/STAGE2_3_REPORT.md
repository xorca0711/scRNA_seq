# Stages 2 and 3: frozen modules and eligibility

25 September 2026. Executed under the owner's authorization of the
[contract plan](../PLAN.md). Every number below is read from the tracked tables in
[`tables/`](../tables/) and their run records. **No expression score was computed.**
Stage 3 read gene names only. No claim grade changes.

## Summary

- **Modules frozen under the pre-registered rule.** The rule was committed before
  any intersection was computed. The freeze asserted that the modules are disjoint
  and cover the union of the three source lists.
- **The shared module is two unrelated overlaps, not one programme.** Development
  and adult injury share 5 genes. Adult injury and lesion share 7. Development and
  lesion share none, and no gene appears in all three lists.
- **A11 is eligible for its test. A5 is not eligible for confirmation.** Every
  module each question needs clears the gene coverage gate. A11 also clears the
  independent-unit floor. A5 fails it, because no local dataset holds independent
  neonatal animals.
- **The coverage method reproduces both completed runs exactly.** Its first attempt
  failed three of 33 checks and refused to report. That attempt is
  [preserved](../tables/stage3_attempt1_refused/README.md) with its cause.

## Stage 2: the frozen modules

After removing the operational genes, the three source lists hold 99
developmental, 398 injury and 98 lesion genes.

| Module | Genes | Role |
|---|--:|---|
| Shared remodelling | 12 | contract module, both questions |
| Development-specific | 94 | contract module, A5 |
| Lesion-specific | 91 | contract module, A11 |
| Injury residual | 386 | reference only |

The shared module splits cleanly by origin:

| Origin | Genes | Members |
|---|--:|---|
| Development and injury | 5 | Areg, Itgb6, Lamc2, Ndnf, S100a14 |
| Injury and lesion | 7 | Actn1, Bax, Cxcl16, Gdf15, Klf6, Ly6a, Nars |
| Development and lesion | 0 | none |
| All three | 0 | none |

### What the two overlaps contain

**Development and injury share mostly type 1 identity.** Four of the five genes
appear in the published type 1 identity list; Areg alone does not. At the level of
these lists, what the developmental and adult transitional states share looks like
type 1-directed differentiation. That is the identity rival the plan declared in
advance. It argues against reading this overlap as a distinct plasticity programme.

**Injury and lesion share a stress signature.** Bax and Gdf15 are canonical p53
targets. Klf6 sits in the Hallmark hypoxia and inflammatory sets, Cxcl16 is a
chemokine, and Ly6a is interferon-inducible. Only Bax and Klf6 fall inside the
frozen Hallmark control sets, so those sets undercount stress genes. The overlap is
more plausibly generic stress than a shared lesion programme.

### How much weight the list overlaps can bear

Little, and in one direction only. Capped marker lists from different studies,
platforms and contrasts share few genes even when the biology overlaps, so low
overlap does not show the biology differs. As a rough chance reference, lists of
these sizes drawn from 15,000 to 20,000 expressed genes would share:

| Pair | Expected by chance | Observed |
|---|--:|--:|
| Development and injury | 2.0 to 2.6 | 5 |
| Injury and lesion | 2.0 to 2.6 | 7 |
| Development and lesion | 0.5 to 0.6 | 0 |

The two observed overlaps sit about two to three and a half times above chance. Zero shared between
development and lesion is what chance alone predicts, so it is not evidence of
dissimilarity. This reference is descriptive, not a test, and depends on the
assumed gene universe.

### Consequence for the joint framing

The contract's premise was one shared object serving both questions. At the level
of these lists, the sharing each question cares about is a different gene set: A5
needs the development-injury pair and A11 needs the injury-lesion pair, and the two
are disjoint. The adult injury programme sits between development and neoplasia,
sharing different genes with each.

The contract still earns its place. It gives both questions one partition, one set
of exclusions and one coverage method, so their results stay comparable. The
pre-registered per-gene flags already route each question to its own pair as a
declared secondary. What should be retired is the phrase "a single shared
baseline". The accurate statement is one partition with two relevant pairwise
components.

## Stage 3: eligibility

### Gate 1, assayed source fraction

Threshold 0.7. The mouse arm needs no ortholog mapping. The human arm maps through
the frozen strict one-to-one table, then checks presence in the cohort's gene index.

| Module | GSE247130 multiome, neonatal and adult | GSE310539 multiome, adult | GSE262927 external injury | GSE308103 human lesions |
|---|--:|--:|--:|--:|
| Shared remodelling | 1.000 | 1.000 | 1.000 | 0.833 |
| Development-specific | 0.979 | 0.979 | 0.957 | 0.894 |
| Lesion-specific | 0.945 | 0.945 | 0.967 | 0.791 |
| Development and injury pair | 1.000 | 1.000 | 1.000 | 1.000 |
| Injury and lesion pair | 1.000 | 1.000 | 1.000 | 0.714 |
| Injury residual, reference | 0.979 | 0.979 | 0.974 | 0.635 |

Every module a question needs passes in the datasets that question uses. The
injury-lesion pair passes the human arm at the margin: five of seven genes survive.

The injury residual fails the human arm, as the full adult injury list did in the
precedent run. Both lose genes at assay presence rather than at ortholog mapping.
The reference module is not scored by either question.

The losses have identifiable causes:

- **Ly6a** has no human one-to-one ortholog; it is the mouse-specific Sca-1 gene.
- **Nars** fails through symbol drift. The gene is now Nars1, so the published
  symbol does not map. The frozen symbol policy forbids alias repair, so the loss
  stands and is recorded here.
- **Most other absences** are obsolete RIKEN and predicted-gene names, as the source
  audit anticipated.

### Gate 2, complete independent units

Floor: three complete biological units per contrast. Counts are cited from tracked
records, not recomputed.

| Question | Contrast | Units | Gate |
|---|---|--:|---|
| A11 | Lesion minus normal | 23 patients | pass |
| A11 | In situ minus normal | 12 patients | pass |
| A11 | Atypical hyperplasia minus normal | 8 patients | pass |
| A11 | Minimally invasive minus normal | 4 patients | pass |
| A5 | Adult injury, external animals meeting both floors | 1 animal | fail |
| A5 | Neonatal development, independent animals | 0 animals | fail |

Human patient counts come from the primary configuration of the completed human
run. The three lesion-versus-precursor contrasts also pass and reuse the same
patients. The neonatal count is zero because the external cohort has no neonatal
arm and every neonatal multiome unit is a pooled library.

### Verdict per question

- **A11 is eligible for its test.** The shared and lesion-specific modules clear
  both gates in the human cohort, with 23 paired patients for the main contrast.
- **A5 is not eligible for confirmation.** Its modules are well covered in every
  mouse dataset, so coverage is not the obstacle. The obstacle is units. A5 can be
  scored descriptively inside pooled libraries, but no local data can confirm it.

## Provenance chain

The partition rule was committed in `6775754` before any intersection was computed.
The freeze record stores the hash of the configuration it read, and that hash equals
the hash of the configuration in that commit, so the freeze used the committed rule
byte for byte. The configuration was edited afterwards only to record status. The
coverage report refuses to run unless the frozen modules still match the freeze
record.

## Instrument check and the refused first attempt

Before reporting any new module, the script must reproduce the coverage values of
both completed runs. The corrected run matched all 33: 27 mouse values across three
cohorts and 6 human values.

The first attempt failed three mouse checks, all Hallmark sets on the external
cohort, all by a few genes. It had read the raw deposit's gene index. The precedent
read a post-QC object whose gene set had been filtered. The corrected run reads that
object and verifies its bytes against the precedent's recorded input hash. This
matters because any A5 scoring of that cohort would reuse the same object.

## Limitations

- **List overlap is a weak measure of shared biology.** It is conservative, so the
  small shared module does not show that the programmes differ. Expression tests
  decide A5 and A11.
- **The developmental source was partly chosen by accessibility.** The
  [source audit](STAGE1_SOURCE_AUDIT.md) names the likely stronger source that could
  not be retrieved openly.
- **Two small modules.** The two pairwise components hold 5 and 7 genes. Scores
  built on so few genes are fragile, and one module passes gate 1 by one gene.
- **Symbol policy costs genes.** Refusing alias repair keeps coverage honest but
  loses genes such as Nars that a curated update would recover.
- **Control sets undercount stress.** The frozen Hallmark sets flag only two of the
  seven injury-lesion genes, although most are stress-responsive.

## Stage 4: decisions for the owner

**Outcome, 25 September 2026.** The owner chose the recommended option in each
question: retain the frozen modules; pre-register A11's test in the Kim 2020 cohort
and run only its eligibility gates first; audit a replicated adult injury time
course for A5; update both register cards. DEVELOPMENT decision 37 records this.
The list below is kept as it was proposed.

The contract stops here. These are proposals, and none has run.

1. **Retain or reject the frozen modules.** If retained, they become the fixed
   inputs for both questions.
2. **Specify the A11 test in A11's own plan.** It needs a frozen estimand, the
   patient as unit, the lesion-versus-normal contrast as primary, a multiplicity
   family, a meaningful effect margin and a validation split by whole patients.
3. **Decide A5's path.** Either accept a descriptive analysis inside pooled
   libraries, or treat independent neonatal animals as a data-acquisition gate. A
   developmental atlas with several animals per timepoint could supply those units,
   but it must not also be the source of the developmental list, or the test becomes
   circular.
4. **Consider register updates.** The A5 and A11 cards could record the sourced
   developmental list, the disjoint pairwise finding and each question's
   eligibility. They are not edited here, because the register is the owner's.

## Next-layer proposals

| Layer | Question it would answer | Feasibility verdict |
|---|---|---|
| Chromatin accessibility | Are the development-injury genes such as Areg and Itgb6 accessible in both neonatal and adult transitional cells? | Descriptive only. The multiome deposits are local, but each condition is one pooled library. |
| Protein and tissue position | Does the integrin beta 6 or GDF15 protein track the lesion contrast in tissue? | Unverified. Depends on whether either target is in the processed spatial panel, which has not been checked. |
| Lineage tracing | Do neonatal and adult transitional cells take the same route? | Not feasible with public data here. The recorded rival suggests the direction differs by age, and testing that needs new tracing experiments. |
