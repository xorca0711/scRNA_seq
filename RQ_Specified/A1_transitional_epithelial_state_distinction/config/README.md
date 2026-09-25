# Prospective contracts

Latest: `closure_analysis_contract.json` freezes HPCS mapping/abstention checks
and the recovered CD44 paired/interaction analysis before new numerical results.
`evidence_closure_scope.md` bounds acquisition and adaptive stopping decisions.
The older `samples.json`/`contrasts.json` below remain historical candidate
contracts. CD44 uses its new verified manifest and contract; TIGIT stays held.

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

The owner authorized revised execution on 25 September 2026. The three original
paired contrasts still prevent model fitting because their identity/QC gates
remain unresolved. `first_batch.json` freezes descriptive PCA, technical peak
audit and the source-lineage denominator rules. `ire1_kira8.json` separately
freezes the verified unpaired ten-mouse RNA comparison, marker panel and pathway
eligibility before expression fitting. It does not bypass the paired contracts.
See the [completed batch](../reports/FIRST_BATCH_REPORT.md) for actual outcomes.
