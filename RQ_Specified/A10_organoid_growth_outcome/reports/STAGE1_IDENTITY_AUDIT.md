# A10 stage 1: identity and join audit

25 September 2026. **Metadata only; no count was read.** Produced by
`scripts/01_audit_identities.py` with its run record in `tables/audit_run.json`.
The plan and the biology are in [PLAN.md](../PLAN.md) and
[RATIONALE.md](../RATIONALE.md).

## Verdict

**The joins hold. The biological unit does not, and the deposit cannot fix it.**

Under the register's own rule, that restricts any analysis of this screen to
within-screen descriptive association, with no claim of independent replication.
The restriction is not a failure of the audit; it is the audit working.

## What joined, and what needed normalizing

Three key formats disagreed between tables, and all three were normalized rather
than worked around:

| Issue | Detail |
|---|---|
| Well labels | The design table is unpadded, for example B2; imaging and species quality control are padded, for example B02 |
| Replicate labels | Imaging uses `plate3-1`; species quality control uses `1-1` |
| No replicate column in the design | Correct by design: the design is the plate layout and applies to every replicate of that plate, so it joins on plate and well only |

After normalizing, the joins are good:

| Check | Result |
|---|--:|
| RNA libraries | 886 |
| Libraries with both day-7 and day-14 imaging | 885 |
| Libraries matched to a design row | 856 |
| Duplicate keys in any table | 0 |
| Wells in the design with no RNA library | 0 |

The 30 unmatched libraries are exactly the 30 of one target that is absent from the
design table, which is consistent with a control added to plates rather than laid out
as a candidate. Eight plate and well combinations appear in imaging and species
quality control but not in the design. One plate replicate has imaging but no RNA.
One well has a single imaging day rather than two.

## Why the biological unit is unresolved

The candidate unit is the replicate index: four plates with up to four replicates,
16 combinations, 15 of them holding RNA. Nothing in the three cached tables says
whether a replicate index is a separate type 2 isolation with its own fibroblast lot,
or a repeat of one preparation.

A single sample record was then probed to see whether the deposit's sample-level
metadata would settle it. It does not. The only design field the record carries is:

```
!Sample_characteristics_ch1 = batch: plate1-1
!Sample_characteristics_ch1 = genotype: Met CRISPR/Cas9 knockout
!Sample_source_name_ch1 = lung mouse alveolar type 2 cells and primary human fibroblasts
```

The deposit calls the label a batch and provides no isolation, donor or lot
identifier, and no field separating a biological repeat from a technical one.

**So fetching all 886 sample records would not resolve this, and is not worth doing
for that purpose.** The probe was one request and it settled the question. The source
name does confirm the two-species design from the deposit itself rather than from the
paper.

What would resolve it: the number of independent type 2 isolations and fibroblast
lots, which is a methods fact, or an author request. Reading the source paper's
methods for that one design fact is an owner decision, because the owner's reading of
that paper is recorded as not started and a previous session's note on an unread
paper was rejected. It is not taken unilaterally here.

## Two further items left open, deliberately

**Which assigned genome is which species.** The two assignment fields are named by
role, not species, and no cached file maps either one to mouse epithelium or human
fibroblast. Guessing would silently swap the epithelial and niche sides, which is the
core of the hypothesis. A stage 3 check settles it without guessing: the perturbed
genes are mouse, so knocking one out should reduce that gene in whichever compartment
is mouse. Until then the two sides are not labelled.

**Which targets are controls.** Not established from these files. The obvious
candidate is the target absent from the design and present in all 15 units, but that
is an inference from replication structure, not a statement in the deposit.

## Read assignment is not clean, and it bounds the design

Median fractions of input reads, across 886 libraries:

| Assignment | Median fraction |
|---|--:|
| Ambiguous | 0.034 |
| Both | 0.112 |
| Neither | 0.069 |

About a fifth of reads are not assigned to one species. That is a real bound on how
cleanly the epithelial and niche transcriptomes separate, and it is a reason the
species read fractions belong in the baseline model as composition covariates rather
than being treated as quality-control trivia.

## One target dominates the replication

Libraries per target are four for most, but two targets are extreme: one appears in
87 libraries across all 15 units, another in 30 across all 15. A target present in
every unit at high multiplicity behaves like a repeated reference rather than a
candidate. Its role must be established before any screen-wide summary, because a
target with 87 libraries would otherwise dominate a pooled fit.

That target is also a lesion-module member, and the same name carries an unresolved
hold in A1 over sample identity. The coincidence is worth noticing and is not
evidence of anything.

## What stage 2 must now assume

- Results are descriptive within-screen associations unless the unit question is
  answered.
- The 30 control libraries and the extreme-replication target are handled explicitly,
  not pooled silently.
- The epithelial and niche sides stay unlabelled until the species check in stage 3.
- Composition covariates are mandatory, not optional.

## What was not done

No count was read. The 303 MB count table was not fetched. No model was fitted, no
programme was scored, and no claim grade changed. One sample record was retrieved as
a bounded probe, and no author contact was made.
