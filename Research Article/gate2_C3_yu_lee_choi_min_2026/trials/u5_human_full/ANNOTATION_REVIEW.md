# Full-cohort reference-label review

Reviewed on 24 September 2026 after all 75 libraries were processed, before
interpreting paired pathway or compatibility outputs. The frozen mapping
retains 555,480 QC nuclei. Primary finest-label uncertainty ≤0.2 retains
308,646 nuclei (55.6%); the ≤0.3 sensitivity retains 362,208 (65.2%). Each
library measures 1,731 of the 2,000 reference-model genes (86.55%).

The reviewed broad populations have coherent multi-marker support in all
five histology groups:

- AT2-like: EPCAM detection 48.6–62.0%, SFTPC 93.2–99.8%, ABCA3 87.2–92.9%.
- Fibroblasts: COL1A1 80.9–89.3%, DCN 82.9–95.0%, LUM 73.9–88.0%,
  PDGFRA 47.4–58.2%.
- Macrophages: C1QA 92.7–97.0%, C1QB 92.5–96.8%, CD68 93.7–97.4%,
  TYROBP 90.1–94.5%, CSF1R 74.2–83.5%.

These are pooled within-population marker diagnostics, not biological
replicates or significance tests. The corresponding per-patient counts and
retention table remain available. Finest subtype names are reference
candidates; broad populations provide the main interpretive view.

Substantial cross-lineage marker RNA remains. For example, SFTPC is detected
in 60.6–78.2% of mapped fibroblasts and 71.8–92.6% of mapped macrophages.
This could reflect assay background or mixed profiles; the marker table
alone cannot identify its cause. Broad identity is supported by the full
marker pattern and atlas mapping, not SFTPC positivity alone. No unvalidated
background subtraction or post-result reassignment was introduced.

Confidence retention differs by histology: median patient-level retention is
67.2% in normal, 67.9% in AAH, 64.5% in AIS, 66.0% in MIA and 50.6% in LUAD.
The all-QC source analysis separately shows that unassigned cells carry
median IL1B count fractions of 52–72% across histologies. Confident-subtype
results therefore concern a selected population and cannot establish the
dominant source across all recovered cells. The ≤0.3 sensitivity measures
dependence on confidence filtering; it does not resolve biological identity
in every unassigned cell.

**Decision:** release descriptive niche and patient-level pathway results
for these reference-compatible candidate populations with the limitations
above. Do not claim independently validated fine subtypes, normal versus
malignant identity, KAC/HPCS assignment, tissue abundance or an unbiased
representation of every lesional state. The unavailable KAC endpoint remains
unreleased. No mapping weights, thresholds or labels changed in this review.

[Full marker values](full_cohort_marker_review.csv) ·
[Broad marker check](broad_lineage_marker_review.csv) ·
[Patient retention](patient_histology_annotation_retention.csv) ·
[All-QC source report](../u5_human_sources/REPORT.md) ·
[Review record](annotation_review_summary.json).
