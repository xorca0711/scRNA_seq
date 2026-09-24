# Human IL-1 source and context profiles

All QC cells from 75 libraries are represented, including cells whose finest reference-label uncertainty exceeds the primary 0.2 cutoff. Repeated libraries are pooled within patient and histology. Raw counts agree exactly with the retained source panel in all 150 library/cutoff checks. Both 0.2 and 0.3 confidence views are preserved.

The figure shows the fraction of observed IL1B counts assigned to each broad population. This is a recovered-cell and library-dependent RNA allocation, not a tissue-composition correction or a secretion measurement. Patient-level detection fractions and full-library-normalized mean expression accompany raw count fractions. Missing or unassayed panel genes are recorded rather than filled with biological zeros.

A substantial fraction of IL1B RNA belongs to cells without confident finest-level labels. The medians below are patient-level count fractions, not pooled cohort fractions. This limits attribution to specific cell types and makes the confidently mapped niche comparisons conditional on label retention; it does not invalidate the measured RNA in the unassigned population.

| Histology | Patients | Median IL1B count fraction in unassigned cells |
|---|---:|---:|
| AAH | 8 | 51.7% |
| AIS | 12 | 68.3% |
| LUAD | 23 | 72.1% |
| MIA | 4 | 66.1% |
| normal | 23 | 63.5% |

IL1RN, IL1R2 and SIGIRR describe antagonist/decoy context. NLRP3, PYCARD, CASP1 and GSDMD are RNA measurements of processing-related components; they cannot establish inflammasome assembly, cleavage, release or mature IL-1β protein. Receptor and antagonist expression should be read alongside ligand RNA before proposing an activation model.

![IL1B source allocation](../../figures/human_IL1B_sources.png)

[Patient source fractions](IL1B_source_fractions.csv), [broad processing/receptor/antagonist context](broad_IL1_context.csv), [source summary](IL1B_source_summary.csv), [assay coverage](source_panel_coverage.csv), [count checks](all_cell_count_parity.csv).
