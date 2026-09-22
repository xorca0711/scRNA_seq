# ES1 results: epithelial-state specificity

Descriptive RNA measurements at 2,000 UMI per retained cell. No fate, reversibility, chromatin or causal result is inferred.

Frozen definitions: 16 panels/modules, including the complete 400-gene published ADI, AT2 and AT1 marker lists and label-free variants. Complete DATP/PATS and developmental maturation signatures remain unavailable; short panels, published marker lists and full Hallmark modules are distinguished.

## Complete published ADI list, excluding labeling genes

The original list contains 400 genes. Removing Cldn4 and Krt8 leaves 398; 390 are present in each multiome feature universe. All reported effects below are labelled-minus-reference detection points, not percent change in expression.

| Unit | Labelled cells, seeds 17 / 29 | Difference, seeds 17 / 29 | Both seeds evaluable |
|---|---:|---:|---|
| P9_control | 182 / 209 | +1.04 / +1.21 | True |
| P9_Cebpa_mutant | 412 / 446 | +2.62 / +2.88 | True |
| 7wk_control | 30 / 26 | +3.05 / +3.85 | False |
| 7wk_Cebpa_mutant | 78 / 61 | +4.04 / +3.81 | True |
| SeV_control | 32 / 36 | +7.61 / +6.53 | True |
| SeV_Cebpa_mutant | 222 / 214 | +5.94 / +7.07 | True |
| wildtype_PBS | 8 / 14 | +5.35 / +2.94 | False |
| wildtype_SeV | 58 / 73 | +5.90 / +5.73 | True |
| AP1mut_PBS | 15 / 11 | +3.75 / +4.49 | False |
| AP1mut_SeV | 46 / 44 | +9.52 / +9.61 | True |

In control-genotype wells the ADI contrast is 1.04 / 1.21 points at P9 and 7.61 / 6.53 in injured adults. This is descriptive evidence that the same two-transcript call has different panel specificity across contexts; it does not identify an injury-specific program or a developmental mechanism.

The seven-week control crosses the cell floor between seeds (30 versus 26 labelled cells). Its contrast is not robustly evaluable. Wildtype PBS and AP1-mutant PBS do not reach the floor in either seed.

External sample check: 25 animals scored; 1 meet both 30-cell group floors for this panel.
This is an external study consistency check on previously analyzed data, not independent confirmation of developmental-program reuse. No neonatal samples are present.
EEM-scRNA-167 at 11.0 dpi has 40/193 labelled/reference cells in seed 17 and 35/197 in seed 29, with ADI differences +12.64/+11.58 points. A single eligible animal does not provide replicated confirmation.

## Interpretation limits

- Full DATP/PATS and developmental maturation signatures and a development-by-injury factorial cohort are still missing. The ADI source caps its published list at 400 markers; no marker list is an exhaustive program.
- Each multiome condition is one pooled library. Technical seeds are sensitivity checks, never replicate animals.
- The full five-gene DATP marker panel overlaps the label and is only a circular control. Primary specificity readings use label-free panels.
- The reference requires Sftpc detection and thus operationally enriches AT2 cells; these are conditional group contrasts, not unbiased state prevalence.
- A score within one panel is comparable across wells at this budget; absolute scores across different panels are not comparable measures of program strength.
- Low cell coverage means insufficient measurement, not absence of a biological program.

## Files

`module_scores.csv`, `within_unit_effects.csv`, `between_well_contrasts.csv`, `retention.csv`, `module_overlap.csv`, `seed_sensitivity.csv`, `external_population_inventory.csv`, `verification.json`, `run_record.json` and `es1_specificity.png` are compact tracked outputs. The larger `gene_detection.csv` is a regenerable local artifact; source XLSX and plot PDF are also untracked. Source factual gene lists and hashes are preserved in `modules.json`.
