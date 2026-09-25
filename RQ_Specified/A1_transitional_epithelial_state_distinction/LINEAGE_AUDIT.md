# Lineage and functional evidence: critique and revised analysis

25 September 2026. The owner challenged stages 3–4 and authorized the revised
main analysis. This audit changes the scientific design; it does not lower the
replication standard. Source-paper findings remain distinct from results computed
in this repository.

## Why the previous stages were inadequate

Stage 3 led with trajectory tools before defining a measured lineage endpoint.
It omitted the most direct DATP and PATS tracing experiments from its operational
design. Stage 4 grouped protein, location and function as if they could all
validate the same distinction. Neither tissue proximity nor regional protein
similarity establishes ancestry, future fate or epigenetic memory. A collection
of attractive plots would not answer A1.

## Established experiments and what their labels mean

| Study | Actual experimental anchor | What A1 can use; boundary |
|---|---|---|
| [Choi 2020, DATPs](https://pmc.ncbi.nlm.nih.gov/articles/PMC7487779/) | Pre-existing AT2 lineage plus injury-associated Ndrg1-CreERT2 and Krt8-CreERT2 tracing; day-9 transitional populations and day-28 differentiated descendants; organoid reformation and IL-1β manipulation | Separate AT2 origin from tracing a transitional population. AT1 contribution and labelled AT2 descendants are distinct endpoints; organoid capacity is not an in-vivo single-cell reversal measurement. Exact induction schedules are figure-specific and must be transcribed before quantitative comparison. |
| [Kobayashi 2020, PATS](https://www.nature.com/articles/s41556-020-0542-8) | Sftpc lineage establishes origin; Krt19-labelled cells are assessed for KRT8 and AGER during repair. CTGF-positive AT2-lineage cells provide chromatin measurements | Strong same-study bridge between a traced regenerative population and measured histone/TP53 profiles, although sequencing and tracing are not the same cells. Extended Data 4 provides numerical source data and reports three mice. Krt19 is not globally specific to alveolar PATS. |
| [Strunz 2020, Krt8 ADI](https://www.nature.com/articles/s41467-020-17358-3) | Sftpc-CreERT2 and Sox2-CreERT2 are induced before injury with washout; labelled origins are assessed among Krt8-high alveolar cells. Figure 7 uses two mice per lineage and multiple sampled regions | Origin convergence is measured. The counts/velocity atlas does not genetically trace every ADI cell to its later destination. Regions cannot become biological replicates; two separate driver experiments do not form a jointly calibrated mixture fraction. |
| [Kathiriya 2022, alveolar–basal intermediates](https://www.nature.com/articles/s41556-021-00809-4) | Purified human AT2 organoids and transplantation into injured mice; mesenchymal context changes alveolar versus basal phenotypes; ABI1/ABI2 precede basal output in culture | Experimental conversion capacity and niche dependence, not endogenous human IPF genetic lineage tracing. Human nuclear labelling tracks graft origin; cross-sectional human tissue similarity does not prove its ancestry. |
| [Auyeung 2022, IRE1α](https://pmc.ncbi.nlm.nih.gov/articles/PMC8957349/) | Krt19-CreERT2 pulse labelling on injury days 3–4, day-14 AGER endpoint under KIRA8; epithelial IRE1α perturbation and integrin/TGF-β measurements; GSE190821 | A direct transition-exit and profibrotic-function anchor for stage 4. Preserve mouse-level summaries of sampled fields. Reduced state abundance and enhanced labelled AT1 output are distinct effects; neither alone identifies a histone mechanism. |
| [Chan 2026, HPCS](https://www.nature.com/articles/s41586-025-09985-x) | Slc4a11-directed inducible tracing in established KP tumours; early baseline versus later descendants; Hopx AT1-like state comparator; growth reporters and selective ablation | Separate reporter-positive current state, permanently traced history, descendant composition and expansion. Some scRNA conditions pool two mice. Cell diversity/phenotypic volume alone is not independent-mouse evidence, and ablation tests population dependency rather than a chromatin mechanism. |
| [Krt8 perturbation, 2023](https://www.jci.org/articles/view/165612) | Genetic perturbation and injury responses accompanying transitional-state accumulation | A functional dependency anchor. A changed state frequency can reflect altered entry, proliferation, survival or exit; it is not automatically a measured lineage-transition rate. |
| [CD44 intermediate-state study, 2025](https://www.nature.com/articles/s41467-025-63735-1) | Protein-sorted AT2 populations and fibroblast-response experiments in defined models | Tests phenotype and profibrotic capability; CD44 sorting is not lineage tracing. Expression of mediators is not the conditioned-medium response itself. |

Do not force DATP, PATS, ADI, ABI and HPCS into a common categorical classifier.
These experiments support different origins, histories, contexts and endpoints.
Their shared RNA programmes motivate testing; they are not ground truth for a
universal state taxonomy.

## Revised stage 3: test observed lineage endpoints first

Construct one row per experiment: driver, permanently labelled compartment,
current-state reporter, induction window, injury/tumour age, harvest/chase,
mouse/pool identity, anatomical restriction, denominator and measured endpoint.
Keep three quantities separate:

1. Origin contribution: labelled cells among the transitional population.
2. Descendant composition: endpoint phenotype among permanently labelled cells.
3. Tissue contribution or expansion: labelled descendants among all endpoint
   cells, or change in lineage abundance. This has a different denominator.

The first measured tracing check is the downloadable Kobayashi Extended Data 4
source table. Recover its panel/time/unit mapping before any inference. Summarize
per-animal values if the table identifies them; otherwise preserve the source's
reported unit and keep the result descriptive. Do not invent animal IDs from
spreadsheet row positions or treat terminal harvests as repeated measurements.

For HPCS, only compare traces when GFP history, mScarlet current state, chase and
source pool are explicit in the deposited metadata. Verify early label purity,
growth/death selection and sorting gates before interpreting descendant spread.
Do not turn a two-mouse pooled sample into two observations or re-use the same
barcode pool as an independent validation set. The existing HPCS RNA projection
is not lineage validation.

Trajectory analysis follows only if it adds a test: does inferred ordering
agree with a held-out measured time/lineage endpoint? Freeze roots from
experimental origin, assess alternatives by sample exclusion and keep observed
destinations separate from model probabilities. No de-novo velocity pipeline or
omnibus cross-study trajectory is required for the first batch.

## Revised stage 4: test an endpoint, then corroborate its phenotype

Prioritize the same-model or same-study link: traced regenerative state with
histone profiles (PATS), sorted epithelial phenotype with fibroblast response
(CD44), and traced tumour state with descendant/ablation outcomes (HPCS). Separate
published functional results from accessible numerical endpoint reanalysis.

For CD44, first estimate within-source positive-minus-negative RNA differences,
stratified by verified genotype. A difference in significance between strata is
not an interaction. If the complete design supports it, freeze an explicit
sort-by-genotype interaction as a later comparison. Do not use CD44 alone to
call a cell pathological: healthy CD44-positive populations are a necessary
counterexample. Gene sets defining the gate cannot independently validate it.

Include the IRE1α perturbation study as a higher-priority functional anchor than
untargeted regional proteomics: it combines a state-labelled population with an
intervention and differentiation endpoint. Audit its expression files separately
from microscopy endpoints; they are not automatically the same experimental units.

IMC and regional proteomics become **optional phenotype/context extensions**.
They proceed only when a frozen biological comparison and donor/ROI crosswalk
exist. They cannot validate ancestry or rescue an unsupported epigenetic claim.
Do not download the complete IMC archive or begin segmentation simply because
those methods were listed in the original plan.

## First execution batch and claim ceiling

- Processed counts: verify column/sample mapping, integer scale, feature universe,
  library size and pairing. Fit only fully verified contrasts. Otherwise produce
  labelled descriptive sample profiles without p-values.
- Direct histone evidence: the downloaded GSE141635 headers reveal different
  caller settings (homeostasis `size 1000 / minDist 2000`; CTGF-positive
  `size 4000 / minDist 4000`). Audit interval geometry as a technical result.
  Remove raw peak-count/overlap differences from biological state evidence until
  common-region quantification or consistently recalled peaks are available.
- Measured tracing: reproduce the numerical PATS endpoint summaries if the source
  workbook maps to the stated populations; preserve the source unit hierarchy.
- Figures: sample PCA, source-block profiles, predefined perturbation markers
  and measured lineage endpoints. The incompatible histone intervals receive a
  technical table, not a biological overlap figure. No inferred arrow is shown
  as a measured fate transition.

This batch initiates real analysis while keeping broader raw-read, single-cell,
spatial and perturbation work contingent on its specific data requirements.

The recovered tracing workbook identifies BleoD12 and individual mice, but its
control rows have numerator = denominator = 0 and stored percentages = 0. These
are undefined lineage fractions. Preserve the published fields for provenance,
recompute percentages with explicit missingness, and do not manufacture a
control-versus-injury fate contrast from them.

## Execution outcome

The [first-batch report](reports/FIRST_BATCH_REPORT.md) records the completed
PATS endpoint reconstruction, ten-mouse IRE1α RNA analysis, two source-block
PCAs and histone-caller audit. These are source reproductions and within-study
analyses, not independent validation of a universal state taxonomy.

The IRE1α paper's Figure 6 establishes the day-3/4 pulse and day-14 AGER endpoint;
its Figure 1K and GEO records identify the distinct day-7 epithelial RiboTag
assay. See the [primary figure legends](https://pubmed.ncbi.nlm.nih.gov/35170357/).
The downloaded Chan workbook contains growth/ablation source panels from
Figures 3, 4 and Extended Data 10, not the Figure 2 lineage-composition values.
It is inventoried but has not been fitted as a lineage-transition dataset.

## Second-batch resolution, 25 September 2026

The [published PATS text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7461628/)
specifies Krt19-CreER induction seven days after injury and day-12 harvest:
a nominal five-day interval. The ED4 reconstruction retains three mice per
marker, separate marker denominators and undefined 0/0 controls. Histone ChIP
uses CTGF-positive cells at day 12; TP53 ChIP uses a separate Sftpc-lineage/
CTGF-positive day-8 preparation. They are distinct experiments.

The HPCS hold is now narrower: author notebooks and bounded recovery of the
GEO observation metadata enabled a [source-composition reconstruction](tables/hpcs_source_composition/source_state_counts.tsv)
with 5,333 traced cells and 22 source labels. It reproduces author group counts
and selected saved fractions without reading expression or refitting annotation.
The [source manifest](tables/hpcs_source_composition/source_manifest.tsv) preserves
driver/chase and conflicting age fields. Mouse/pool independence is still
unverified; current mScarlet is absent for traced rows. The measured tracing
assignment and current RNA-state labels do not supply a transition rate or
same-cell current-reporter measurement. See the [second-batch report](reports/SECOND_BATCH_REPORT.md).
