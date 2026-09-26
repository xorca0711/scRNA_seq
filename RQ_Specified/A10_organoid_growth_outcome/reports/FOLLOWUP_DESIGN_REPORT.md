# A10 follow-up design diagnostic

26 September 2026. Phase A was declared in commit `53f5c7a`. Script 07 reused the
hashed metadata and tracked joins; no model was fitted. See
[diagnostic run](../tables/followup_v1/diagnostic_run.json) and
[plan](../FOLLOWUP_PLAN.md).

## Findings that change the next analysis

- **885 eligible wells remain.** RNA/imaging joins reproduce the stored area
  columns exactly, including the one missing day-14 observation. The 886 GEO
  library identifiers match the RNA libraries one to one via explicit description
  fields, not by interpreting human-readable titles.
- **Target allocation is mostly plate-specific.** Of 203 targets, only MECOM,
  RNF43, TIGIT and TDTOMATO appear on multiple plates. After the original
  TIGIT/TDTOMATO exclusion, only MECOM bridges plates 1/2 and RNF43 bridges 1/4;
  plate 3 has no target overlap with any other plate. A target-adjusted four-plate
  contrast is not identifiable in this restricted design. A whole-plate holdout
  tests simultaneous plate and target-distribution shift.
- **The inspected imaging label does not explain the plate split.** All 1,883
  imaging rows use `SAM1`. This does not prove identical acquisition/calibration.
  GEO describes mean area, count and covered fraction but provides no physical
  area units or transformation definition. `area * count / coverage` is not a
  fixed quantity, so no conversion to total cell mass or physical area is inferred.
  Continue only with the deposited numeric area statistic, retaining the original
  transformation and adding a fixed untransformed-statistic sensitivity.
- **Observed differences are descriptive.** Median deposited day-14 area is
  8.345/8.342/8.763/8.787 across plates 1–4. Median epithelial RNA fraction is
  0.746/0.699/0.789/0.762. These differences do not by themselves explain why the
  model fits plates 2/4 poorly; target mix and unmeasured preparation remain rivals.
- **Preparation identity is still unresolved after the complete GEO sample audit.**
  All 886 records expose cell type, genotype and batch, but no structured
  animal/isolation/donor/lot/preparation identifiers. The common processing text
  describes RNA preparation and species assignment, not a sample-to-preparation
  crosswalk. This is a bounded statement about the inspected public metadata,
  not proof that source-paper methods or authors cannot resolve it.
- **Guide sequences do not measure editing success.** The layout has three guide
  columns. It does not supply well-specific editing efficiency or establish that
  these are independently replicated guide experiments. The previous transcript
  diagnostic is not promoted to functional validation.

Sources: [all sample fields](../tables/followup_v1/diagnostic_geo_sample_fields.tsv),
[unique metadata catalogue](../tables/followup_v1/diagnostic_geo_field_catalogue.tsv),
[target overlap](../tables/followup_v1/diagnostic_target_overlap.tsv),
[group summaries](../tables/followup_v1/diagnostic_group_summary.tsv),
[segmentation allocation](../tables/followup_v1/diagnostic_segmentation.tsv) and
[guide layout](../tables/followup_v1/diagnostic_guide_layout.tsv). The public
[GEO family metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307112/soft/GSE307112_family.soft.gz)
was retrieved on 26 September; its hash is in the run record.

## Gate decisions

Proceed with descriptive programme decomposition and a joint plate/target shift
stress test. Do not call the latter independent preparation validation or use it
to separate biological from technical plate effects. Freeze the amended metric,
programme blocks and interpretation rules before those fits.

Prune a target-adjusted plate-effect model after removing TIGIT/TDTOMATO, target
functional conclusions without editing-efficiency evidence, and expansion into
additional pathways merely to rescue transfer. The historical unperformed BH
procedure remains unperformed; the new amendment explicitly uses descriptive
effect sizes without p-values, confidence intervals or biological replicate claims.

Technical record: the first metadata-identity assertion stopped before any output
because GEO titles are descriptive sentences. Parsing the exact `Library name:`
description field resolved it without changing the cohort or scientific criteria.
