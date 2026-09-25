# A1 second batch: frozen scope and decisions

25 September 2026. Owner authorized proceeding after the remaining-work review.
This is a source-informed sensitivity/direct-mark follow-up; first-batch results
have already been inspected. It is not an unseen-data preregistration.

## IRE1α stability

Keep the completed ten-mouse primary model and every original output unchanged.
The expression universe is the first run's 14,811 retained genes. Re-estimate TMM,
dispersion and robust edgeR quasi-likelihood fits for each of ten leave-one-mouse-out
subsets using `~ batch + sex + group`. Check rank and residual degrees of freedom;
record inestimable fits rather than choosing a favourable alternative. Do not
drop a mouse from the primary analysis based on these results.

Fit the S061 within-batch sensitivity (three mice per arm) with `~ sex + group`.
S135 has two mice per arm, below the existing three-unit inferential floor:
report descriptive mean normalized-CPM log2 ratios with a frozen 0.5-CPM offset,
without p-values or FDR. Do not interpret nonsignificance across unequal batches
as a treatment-by-batch interaction.

Retain the eight frozen markers, plus all four original FDR discoveries explicitly
labelled post-selection. Report direction retention, logFC range and maximum
change from the primary estimate across omissions, alongside whole-universe
effect correlations. These ranges are sensitivity ranges, not confidence
intervals. Keep all tested genes available, with BH within each fit; no union of
sensitivity discoveries is promoted as a new discovery family.

For the two originally eligible pathways, retain the original memberships and
estimated-correlation CAMERA method. Keep terminal UPR held by its original
coverage rule. Report all omission/eligible within-batch tests; neither the best
p-value nor an alternative correlation model replaces the primary result.

## Direct marks and identity audits

Retrieve bounded processed inputs and source metadata for PATS GSE141635,
Tsutsui GSE289683/GSE291333 and the GSE150527 differentiation reference. Record
assembly, normalization, mark/control, state, preparation and replicate labels
before quantification. T2T-CHM13, mm10 and hg19 coordinates are separate; do not
join them by an unverified coordinate conversion. Inspect matching H3/input and
normalization evidence. Do not compare incompatible called-peak abundance.

Where compatible signal is available, freeze locus definitions before extracting
values; report direct-mark profiles with their actual replication limit. Two
replicate labels do not pass the three-unit inferential floor. No global
acetylation conclusion from CPM tracks without the required spike-in evidence.
Normal-differentiation methylation is a reference, not a purified transitional
population or a replicated DMR test.

Audit TIGIT pool independence and CD44 matrix/GSM/genotype mappings using public
deposits, methods and source code. Update verified identities only with explicit
evidence; title order and expected marker expression cannot resolve mappings.
If unresolved, preserve the hold and record the evidence needed. Conditional
paired models require a separate frozen assay/QC contract once identity passes.

Reconcile PATS pulse/chase and search for exact HPCS descendant tables. Additional
trajectory and protein/spatial analysis proceeds only when a specific measured
endpoint and independent biological units support its test. Archive-wide raw
sequencing or IMC acquisition is outside this bounded follow-up until resource
and sample mapping are verified.

Each package writes fresh tables, figures and a run record, with input, contract,
source and output hashes. Report completed work separately from held tests.
