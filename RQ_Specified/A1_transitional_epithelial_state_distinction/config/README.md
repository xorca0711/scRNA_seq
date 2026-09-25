# Prospective contracts

`samples.json` preserves a **candidate** title-derived pairing for the first
processed-count pilots. `candidate_pair_id` must never substitute for a verified
`biological_unit_id`. Leave unknown pool membership and source identity null.
Resolve identities using paper supplements and deposited sample records;
record the exact evidence before `identity_verified` becomes true.

For each count library add an exact `counts_column`. Confirm disjoint source
members across independent units, complete positive/negative pairs, tissue,
genotype, harvest/sort, batch and assay QC. Technical-library aggregation requires
a documented adapter and its input hashes, not simply a renamed column.

`contrasts.json` contains three **proposed**, unfrozen contrasts: paired TIGIT
bulk ATAC, and paired CD44-sorted RNA separately within WT and mutant mice.
Their case is positive and reference is negative. A minimum of three genuine
independent pairs is an eligibility floor, not a power calculation. The mutant
and WT tests are separate estimands; a genotype-by-sort interaction is a later
model, not the difference between their significance labels.

The initial count-model defaults are TMM, edgeR quasi-likelihood with source-unit
blocking, robust dispersion estimation, `filterByExpr` with min.count 10 and
min.total.count 15, and BH within each predefined contrast. They must be reviewed
against assay QC, library composition and the chosen feature universe before
freezing. For ATAC, audit common peaks, blacklist/CNV effects, background/depth
and normalization; raw count format alone is insufficient. Any combined claim
across contrasts needs a declared joint testing family rather than choosing
whichever per-contrast result is significant.

Only after plan review and verification should execution be authorized and a
contrast marked frozen. The current files intentionally prevent model fitting.
