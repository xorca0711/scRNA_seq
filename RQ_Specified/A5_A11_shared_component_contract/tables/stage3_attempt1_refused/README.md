# Stage 3, first attempt: refused by its own precedent check

These three files are the complete output of the first stage 3 run. They are kept
because a failed check is a result, not a mistake to hide.

The script must reproduce the coverage values of both completed runs before it may
report any new module. The first attempt ran 33 checks. Thirty matched exactly,
including every human check and every mouse check on the two multiome cohorts.
Three failed, and all three had the same pattern.

| Dataset | Module | Precedent | First attempt |
|---|---|--:|--:|
| External viral-injury cohort | Hallmark hypoxia | 196 | 199 |
| External viral-injury cohort | Hallmark inflammatory response | 195 | 197 |
| External viral-injury cohort | Hallmark p53 pathway | 197 | 198 |

**Cause.** For this cohort the first attempt read gene names from the raw
per-sample deposit files. The precedent read them from the post-QC object in the
viral-injury paper folder, whose gene set was filtered during quality control. A
few weakly expressed stress and immune genes are absent there. Epithelial marker
lists matched exactly on the same cohort, which fits that explanation, and reading
the precedent's code confirmed it.

**Consequence.** The script reported no coverage, as designed. The corrected run
reads the post-QC object and verifies that its bytes match the precedent's recorded
input hash. The corrected output sits one level up.

**Why it matters.** Any A5 scoring of the external cohort reuses that post-QC
object. Coverage computed from the raw deposit would have described a gene set the
actual test never sees.
