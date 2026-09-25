# Full-cell IPF LR robustness and source context

Both cohorts were evaluated with native LIANA CellChat-like magnitude. All-cell inputs retain each cell’s original full-library total and full-assay maximum, preserving normalization and native magnitude scaling. Two-donor/two-resource comparisons against full-gene inputs had zero score differences. The computational maximum column is not a gene and never appears in an interaction.

| Cohort | Annotated IPF/control donors | Annotated cells | Primary LR donors | Native calls | Empty calls |
|---|---:|---:|---:|---:|---:|
| GSE136831 | 60 | 243472 | 25 | 350 | 0 |
| GSE135893 | 22 | 89326 | 17 | 236 | 13 |

## Interpretation

The full-cell, two-seed 500-cell cap, 30/50/100-cell floor, 5/10/20% detection and two-resource outputs are all retained. Sensitivities vary one choice at a time. Availability overlap and rank correlation are different checks; correlation is calculated only on common edges and cannot conceal absent edges. No missing edge is imputed as zero. The cap seeds use stable donor/subtype seeds and are not an exact replay of the original pilot sampling.

All annotated compartments were scanned for source, receptor, inhibitor and processing-gene RNA. The tables distinguish within-label expression/detection from proportions of recovered cells. Neither quantity measures mature cytokine secretion, receptor activation or absolute tissue abundance. WNT-family resource edges are explicitly exploratory.

Native LIANA 1.10.0 raises an indexing error when no edge passes expression eligibility. In the second cohort these cases were independently checked across every resource subunit and fixed pair, then recorded as empty calls. Completed cached calls were retained; the original failed run and code are preserved under `.history`.

Complete ≥50-cell triads, distinct from pair coverage: GSE135893 IPF: 2, GSE135893 control: 1, GSE136831 IPF: 2, GSE136831 control: 4. Joint associations require at least ten complete patients, so this coverage must be checked for each proposed association.

![Full-cell robustness and sources](../../figures/ipf_full_cell_robustness_and_sources.png)

Tables: [sensitivity summary](sensitivity_summary.csv), [donor sensitivity](donor_sensitivity.csv), [source/context summary](all_label_context_summary.csv), [triad coverage](complete_triad_coverage.csv), [independent checks](validation.json). Per-cohort folders contain donor-level IL-1 context, cell composition, native score caches and exact resource definitions. This package does not complete the human-lesion, mouse-perturbation or ligand-target arms.
