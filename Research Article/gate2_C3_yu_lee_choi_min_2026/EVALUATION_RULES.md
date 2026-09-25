# How to read eligibility and surviving evidence

The owner requested completion of the remaining planned analyses with the
declared criteria and sensitivity checks, followed by evaluation of what
survives. No numerical criterion was changed in response to the observed
number of discoveries.

## Different reasons for an absent headline result

| Status | Meaning | What remains useful |
|---|---|---|
| Tested, no q < 0.05 discovery | The design and assay supported the test; the declared multiplicity correction did not yield a discovery | Effect estimates, individual units and sensitivity results; this is not proof of no effect |
| Insufficient biological replication | Too few eligible mice, donors or paired patients for the declared inferential step | Individual-unit measurements and explicitly descriptive comparisons |
| Insufficient cell or gene coverage | A specific subtype, signature or contrast does not meet the declared coverage rule | Coverage tables and predeclared lower/higher-threshold sensitivities |
| Missing identifying information | Treatment history, source-state labels or independent pathology regions cannot be established | Analyses of broader, identified compartments; lowering a numeric cutoff cannot reconstruct the missing information |
| Conditional target analysis ineligible | A receiver does not supply ten mapped significant targets in the declared direction | Its completed DE/pathway results; no target-threshold relaxation after seeing results |
| Assay does not measure the claim | RNA compatibility does not measure secretion, receptor activation, cell contact, flux or ancestry | A hypothesis for a functional or spatial experiment |

These are contrast- and endpoint-level statuses, not whole-dataset quality
verdicts. A dataset may support one question and leave another unanswered.

## Which standards are universal and which were chosen here?

Biological replication, identifiable comparisons and faithful assay
interpretation are fundamental. Accounting for between-replicate variation
is supported by [Squair et al., 2021](https://doi.org/10.1038/s41467-021-25960-2).
The exact thresholds of three units, 50 niche cells, 100 alveolar cells for
the fraction endpoint, and ten mapped targets are project decisions rather
than universal field requirements. Passing a minimum does not guarantee
power or reliable mechanistic interpretation.

The primary use of estimated inter-gene correlation in CAMERA is also a
project-level modeling choice. The fixed-0.01 setting is a declared
sensitivity. Both were evaluated through the official implementation; this
does not make either setting a universally required standard. The completed
[correlation comparison](trials/u6_completion/pathway_correlation_sensitivity.csv)
shows substantial dependence on this assumption. In particular, the human
family has 0/279 primary discoveries versus 148/279 in the fixed-correlation
sensitivity. That contrast warrants methodological interpretation, not a
quiet replacement of the primary result or a claim that no biology exists.

Primary analyses retain their original settings. Existing 30/50/100 niche
cell floors, fraction-denominator sensitivities, expression thresholds,
resource alternatives, sampling seeds, prior counts and annotation-confidence
checks expose dependence on those choices. Their results remain labelled as
sensitivities rather than replacing a primary result because they look more
favourable. For the mouse fraction question, the 50-cell sensitivity retains
all seven early animals, but does not supply missing KAC annotations.

## Final synthesis

Evaluate separately whether a result is reproducible within its technical
checks, supported by independent biological units, consistent across cohorts,
specific to the proposed state/context, and supported by intervention or
function. Do not collapse these dimensions into one pass/fail score. Shared
programs, observational RNA edges, spatial mixtures and intervention outcomes
answer different questions. Related assays and reused specimens do not add
independent replication.

Preserve the full tables, negative and ineligible cases, source definitions,
and original frozen specifications. Any later exploratory revision should be
documented alongside the original analysis with its scientific reason.
