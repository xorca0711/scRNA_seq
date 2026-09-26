# Next dataset gate: connect expression to measured outcomes

**Status update, 26 September 2026:** the first candidate's joins and two A10
model specifications are complete. Start with the
[revised result](../RQ_Specified/A10_organoid_growth_outcome/reports/STAGE4_REVISED_REPORT.md)
and [current handoff](../PROGRESS.md). Biological preparation identities remain
unresolved. The dated gate below is retained as provenance, not an unexecuted
instruction or confirmation of independent replication.

Decision date: 22 September 2026. This is a prospective analysis gate following
the corrections, not a declaration of independent confirmation. Public GEO
metadata were retrieved and saved as [versioned extracts](remediation/2026-09-22/dataset_metadata/).
The [candidate manifest](../analysis/claims/dataset_candidates.json) records
eligibility, missing design information and permissible interpretations.

## First candidate: GSE307112 organoid perturbation screen

The [GSE307351 superseries](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307351)
contains 886 RNA libraries in [GSE307112](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307112)
and four spatial samples in [GSE307128](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307128).
The RNA samples combine mouse AT2 cells and human fibroblasts; processing assigns
reads by species. The deposited file descriptions connect well-level RNA to
day-7/day-14 organoid counts, area and coverage. This supports an outcome-linked
question, whereas RNA state scores alone do not establish repair.

**Decision: prioritize a bounded design-and-QC pilot.** The RNA unit is a culture
well, not a single cell or automatically an independent animal. There are 15
distinct deposited RNA plate labels; their relationship to biological isolation
batches must be established. Spatial titles indicate two WT and two Nkx2.1-KO
blocks, without resolving whether these represent four independent animals.
The spatial subseries repeats the organoid design text despite spatial sample
processing: use sample-level assay metadata and resolve the discrepancy.

Before fitting an outcome model:

1. Join `plate_design`, `imaging_outputs`, `xenome_stats` and gene-count columns
   by unique well identifiers. Report missing wells, duplicate keys, controls,
   species ambiguity and plate coverage before inspecting target effects.
2. Recover biological preparation, donor/isolation and technical-repeat IDs.
   Separate target, guide and plate. Do not count wells or guides as animals.
3. Freeze a small programme set and an imaging endpoint. Prefer day-14 organoid
   area conditional on day-7 area, with number and coverage as declared secondary
   endpoints. Area measures growth/morphology, not mature AT1 fate.
4. Model plate/block structure only where estimable. Hold out entire biological
   batches for prediction; random well splits may leak shared preparations.
   Keep known Nkx2.1 results as positive controls, not new discoveries.
5. If biological identities remain unresolved, publish descriptive within-screen
   associations and technical robustness only. Do not promote a repair mechanism.

This starts order 5. It does not download sequencing reads or launch an
unbounded reanalysis. The imaging linkage and species separation make this a
stronger next feasibility target than another ligand ranking on the same atlas.

## Second candidate: GSE242510 myeloid–fibroblast coculture

[GSE242510](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242510) contains six
PIPseq libraries: two each for WT fibroblasts alone, WT fibroblast–macrophage
coculture, and P2rx4-KO fibroblast–WT macrophage coculture. Deposited labels do
not identify animal pooling or establish that the two repetitions are biological.
The perturbation is fibroblast P2rx4, not an Arg1 deletion in these RNA samples.

**Decision: descriptive secondary pilot; not a confirmation cohort.** First map
library preparation to animals and identify comparable fibroblast/macrophage
states. Any programme contrast must distinguish coculture exposure, genotype
and composition. Even two independent preparations per group would remain
below the current three-unit minimum, which itself is not a power justification.
RNA cannot substitute for ornithine flux, protein activity or fibrosis outcomes.

## Existing-data limits and portfolio sequence

ES1's GSE262927 check yields only one animal passing both frozen state floors;
it cannot close the confirmation gap. Bulk ATAC peak presence cannot establish
same-cell reopening. Additional analysis of these deposits should answer a
specific sensitivity question rather than create nominal replication.

Next: complete the first candidate's identifier/replication audit; run the
bounded outcome pilot if its design supports it; then seek independent biological
confirmation. The [portfolio summary](PORTFOLIO_SUMMARY.md) and four-figure guide
already present the completed corrections. Scientific selection is documented
here; private outreach rankings and contact plans remain outside the repository.
