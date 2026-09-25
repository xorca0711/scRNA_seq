# Human paired niche analysis

All 75 libraries from 23 deposited patient labels were processed. 555,480 nuclei passed full-assay QC; 308,646 met the primary reference-label confidence criterion. Repeated libraries were pooled within patient and histology. AAH, AIS, MIA and LUAD comparisons remain separate.

## Analysis and interpretation

Pathways use raw-count pseudobulks, filterByExpr, TMM, voom, a patient-blocked design, and official CAMERA. Global BH covers all eligible set/subtype/broad-view/histology contrasts within each declared analysis family. At least three complete patients are required. Broad and subtype views overlap and are not independent replications.

The primary estimated-correlation family contains 279 tests, with 0 q < 0.05. The fixed-0.01 sensitivity has 148 q < 0.05; it does not replace the primary result.

The descriptive RNA-compatibility arm contains 1080 primary resource/view contrasts, with paired patient values and omission ranges independently reconstructed. Native LIANA yielded 123,950 score rows across resources and sensitivities. Pair coverage, native expression eligibility and missing edges are retained separately; absent edges are never scored as zero.

Confidence 0.3, 30/100-cell floors, largest-library selection, prior counts 0.5/2, both LR resources, expression thresholds and cell caps remain labelled sensitivities. The full raw assay supplies library totals and normalization; the computational gene panel never replaces the assay denominator.

## Identity and assay limitations

These are reference-compatible candidate compartments, including AT2-like cells. The atlas does not establish nonmalignant status and does not supply a KAC or HPCS classifier. Low-confidence/unassigned cells remain in QC records and are excluded from subtype comparisons. Confidence retention by histology and full-cohort multi-marker support are reported separately in the full annotation review. Unassigned cells carry median IL1B count fractions of 52-72% across histologies; these subtype results cannot establish the dominant source across all recovered cells. Surfactant RNA is widespread outside epithelial assignments, so no cell type is assigned from SFTPC alone. Gene expression and native LR magnitude do not establish secretion, activation, physical contact or a feedback mechanism. Cross-sectional paired lesions do not prove a progression trajectory.

The planned KAC/NF-kB joint association remains gated by the unavailable source-defined KAC measure. Human HPCS/ADI/ISR program comparisons are a separate descriptive specificity extension and are not substituted for that endpoint. Conditional ligand-target refits use their own target and expression gates.

![Human paired pathways](../../figures/human_paired_recipient_pathways.png)

![Human niche contrasts](../../figures/human_paired_niche_RNA_contrasts.png)

[Pathways](camera.csv), [paired eligibility](paired_eligibility.csv), [primary RNA contrasts](primary_compatibility_contrasts.csv), [individual patients](primary_compatibility_patient_values.csv), [native-call coverage](liana_call_coverage.csv), [sensitivity summary](liana_sensitivity_summary.csv), [annotation review](../u5_human_full/ANNOTATION_REVIEW.md), [checks](validation.json).
