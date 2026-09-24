# Initial analysis batch

Generated 2026-09-24T08:27:47.570941+00:00. The owner authorized staged
execution after final review. This report covers completed modules; it is not
the final multi-cohort study report.

## Completed mouse input and annotation diagnostics

Seven early IgG/anti-IL-1beta libraries were acquired and processed. Of
28,748 deposited cells, 28,243
pass the frozen gene-count/mitochondrial filters. No duplicate matrix
coordinates were found. The acquisition/QC stage took 107.58 seconds; the
treatment-blind clustering stage took 243.27 seconds and yielded 33 review
clusters. Source data and preparation inputs were preserved.

The marker review supports distinct stromal, immune and epithelial populations,
but several myeloid clusters need more specific labels and the alveolar group
does not yet establish a defensible KAC class. The provisional marker gates
are not used as a substitute for validated KAC or macrophage subtype identity.
No mouse treatment-effect conclusion has been released.

## Completed IPF pathway analysis

| Cohort | Primary tests at 50-cell floor | Primary global q < 0.05 | Fixed-correlation sensitivity q < 0.05 |
|---|---:|---:|---:|
| GSE136831 | 49 | 0 | 14 |
| GSE135893 | 36 | 0 | 12 |

The declared primary analysis estimates residual inter-gene correlation.
Neither cohort has a pathway passing its global 0.05 FDR threshold. Fixing
correlation at 0.01 changes the results substantially. These sensitivity
findings do not replace the primary analysis or establish a robust niche
mechanism. Failure to pass FDR does not establish biological absence.

Fibroblast inflammatory-set directions differ between these cohorts in the
broad-compartment comparisons. Their mixture and small control samples need
to remain visible. At the primary floor, broad fibroblast IPF/control n is
18/4 in GSE136831 and 7/3 in GSE135893. Most fine fibroblast subtypes lack
three controls, so they are explicitly ineligible. Broad and subtype views
are not independent biological replications; AT2 broad/sole-subtype rows
are duplicate views and were retained in the declared correction family.
The floor-30 sensitivity is conditional on the cache's original broad-
compartment >=50-cell donor selection, not recovery of all 30-49-cell donors.

![IPF pathway primary and sensitivity results](figures/ipf_pathway_primary_and_sensitivity.png)

The table-integrity check independently recalculated global BH families and
verified sample/assay eligibility. These checks do not validate a biological
claim. An initial output merge failed because CAMERA returns different
columns for the two correlation settings; its schema was corrected and both
cohorts completed with the same frozen statistical choices. The failed run
record is preserved.

## Completed initial descriptive ligand-receptor run

The GSE136831 ligand-receptor run completed at 17:43 KST on 24 September 2026,
in 1,187 seconds (about 20 minutes). It extracted 25 eligible donors
(15 IPF, 10 controls) from full raw counts with deposited labels. It includes
macrophages omitted from the old epithelial-stromal cache. The first
descriptive LIANA pass caps each donor/
subtype at 500 seeded cells, preserves fixed source-target pairs, and compares
consensus with CellChatDB. It does not provide biological-replicate P values,
an all-cell sensitivity result or a causal communication claim.

All 25 donors completed both resources: 3,449 consensus rows and 884
CellChatDB rows across 50 nonempty tables. These are donor/edge observations,
not counts of independent or significant interactions. Full-source count
parity passed. [Output checks](trials/u5_ipf_liana/validation.json) verified
file completeness, fixed source/target and resource membership, unique rows,
score bounds, detection eligibility and disabled cell-permutation P values.
These checks do not establish biological validity or interpret missing edges.

## Continuing and gated work

The initial compute job has ended. Custom RNA-compatibility contrasts,
full-cell and other sensitivities, ligand-target interpretation and new gallery
rendering are not yet completed.

KAC identity, later-treatment strategy mapping, human lesion crosswalks,
spatial access and required context-specific coverage remain release gates.
Zenodo's linked annotated files are restricted. These are scientific/access
conditions, not an outstanding request for user approval.

Sources and execution evidence: [final review](FINAL_REVIEW.md),
[mouse QC](trials/u3_acquire_qc/run_record.json),
[cluster diagnostics](trials/u3_cluster_review/run_record.json),
[pathway specification](trials/u5_ipf_spec/specification.json),
[pathway run](trials/u5_ipf_pathways/run_record.json),
[table checks](trials/u5_ipf_pathways/validation.json), and
[ligand-receptor run](trials/u5_ipf_liana/run_record.json).
