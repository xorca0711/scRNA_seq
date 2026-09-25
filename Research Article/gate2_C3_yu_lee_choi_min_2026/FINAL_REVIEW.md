# Adversarial review and execution release

24 September 2026. The owner authorized launch after self-review and fixes.
This supersedes the earlier confirmation hold. Scientific release gates remain.

## Findings and changes before treatment effects

| Finding | Consequence | Change |
|---|---|---|
| Small-n exact tests were specified without checking attainable P values | Three versus four units yield only 35 allocations; four versus four yield 70. For a two-sided absolute mean-difference test the best possible raw P is at least 1/35 and 2/70 respectively. Holm across the two endpoints cannot reject either at 0.05 even in ideal separation | Retain declared endpoints/correction; explicitly report this limitation. Do not change to one-sided tests or parametric tests to obtain significance. Mouse fraction results emphasize effect sizes and robustness; no equivalence claim from non-significance |
| Subtype coverage can change under treatment | Complete-case within-state estimates condition on surviving/captured states and may be selected by treatment | Report attrition and eligible-unit counts by arm; keep whole-compartment and subtype estimates distinct; no missing-state imputation. Label within-state effects conditional on coverage |
| A custom LR statistic can look like validated communication inference | Min-subunit logCPM compatibility is an analysis choice, sensitive to low counts and complex bottlenecks | Label it exploratory RNA compatibility; keep ligand/subunit effects visible, show prior.count 0.5/1/2 sensitivity, and do not claim a validated communication probability |
| Gene reuse can create edge-response association | Receptor/ligand and classification genes may contribute to both axes | Remove edge components and annotation genes from response-score sensitivity; specify the retained genes and coverage. NicheNet from the same RNA is not independent validation |
| Author code depends on external annotated objects | Raw public matrices do not supply a verified KAC label or all subtype metadata | Acquire public GEO counts and perform QC; block KAC inference until annotation provenance is sufficient. Do not label KRT8-high cells KAC by fiat |
| The Zenodo entry was treated too optimistically | Public metadata exist but the files are restricted; an older author README DOI returns 404 | Mark access as restricted/unavailable, not a large public download pending. No restricted-file retrieval attempted |
| The core is enriched for epithelium | Macrophage/fibroblast yield may be inadequate and proportions are not whole-lung abundance | Gate each sender/receiver pair after coverage audit; retain the separate IPF niche arm, without substituting it for blockade evidence |

The exact-test bounds follow the combinatorial design, not an observed effect.
For LR BH testing the discrete P grid is reported with family size and attainable
thresholds; discovery is not promised from the small mouse series. CAMERA has
different model assumptions and is not a workaround for a failed permutation test.

## Source resolution at launch

The [primary methods](https://doi.org/10.1016/j.ccell.2025.10.004) identify four
mice per treatment/time strategy, one QC exclusion leaving 31 libraries, and
epithelial enrichment from lung single-cell suspensions. This supports one
retained library per animal and a fresh-cell interpretation despite the GEO
FFPE characteristic. Original animal identifiers and the exclusion's detailed
reason remain unavailable. GSM IDs are traceable analysis-unit aliases, not
invented animal IDs. No replicate ordinal is treated as a matched pair.

The methods mention immediate and delayed treatment strategies. The GEO
start field assigns end-of-NNK to all libraries; the later endpoint needs
its strategy crosswalk resolved. Therefore the initial acquisition/QC batch
is the seven 3-month IgG and anti-IL-1beta libraries. The 7-month contrast
remains gated until its schedule is established.

[Author code](https://github.com/FuduanPeng/LungPCA_Code) is frozen at tree
`fba87dbee8dbae497b15b81a4b0fa79edccd17b2`. Figure 7 plotting code loads an
external annotated object and is not an end-to-end annotation pipeline.
Its Figure 7F denominator is epithelial cells; our alveolar denominator is a
new estimand and must not be presented as exact figure reproduction. A source-
denominator sensitivity is required once valid labels exist.
The [Zenodo record](https://zenodo.org/records/17172149) marks files restricted.
Source retrievals/hashes and unsuccessful endpoints are in
[U2 provenance](trials/u2_source_audit/run_record.json). An HTTP-200 error body
from the public BioC endpoint is rejected as unavailable full text.

## Released work and timing

The authorized run starts with a checkpointed public-count acquisition and
QC/compartment-feasibility stage. It cannot silently pass an annotation,
sample-design or resource-freeze gate into treatment inference. QC candidate
marker gates are provisional and never promoted to KAC or definitive subtypes.
Runtime estimates use measured transfer and per-library processing rates.
The complete multi-cohort interpretation has no defensible fixed end time
while annotation access and later-treatment mapping remain unresolved.
