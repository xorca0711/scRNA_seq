# IPF ligand–target prioritization

Completed for eligible receivers using the checksum-verified NicheNet v2 human prior (Zenodo record 7074291). The implemented statistic is the official Pearson correlation between prior target potential and a binary DE-target vector; this is a base-R implementation of that statistic, not a claim that the complete nichenetr package or its AUPR metrics ran. Full-data DE fits reproduce the preceding pathway run, and the Pearson calculation was independently checked against centered vectors.

## Findings and limits

GSE136831 AT2 and GSE135893 broad fibroblasts fail the ten-mapped-target gate. No ranking is inferred for these receivers. IL1A/IL1B qualify for the GSE136831 fibroblast receiver; neither is the leading fit to its upregulated targets. Their better ranks for downregulated fibroblast genes do not establish inhibition or repression: the prior is unsigned. The two cohorts therefore do not provide a replicated IL-1-specific recipient programme.

Both focused triad and source-agnostic rankings are restricted to the planned ligand families with available source-panel measurements. They are not a genome-wide screen of every possible ligand. The full prior scores are retained separately and are not promoted into expression-supported findings. A prior receptor missing from the retained panel is marked unevaluable, not biologically absent. IL1R2/SIGIRR are excluded from activating candidate gates; canonical IL-1 and TGF-beta candidates require their essential receptor partners.

## Donor-omission stability

Every omission reruns receiver filtering, TMM, voom, DE/BH and target selection. Source/receiver expression eligibility is checked again. Omission fits with fewer than three donors per arm or fewer than ten mapped targets remain missing. The stability table reports planned versus eligible omissions; ranks can depend on the number of eligible candidates and are not confidence intervals. Broad and subtype views overlap and are not independent evidence.

| Cohort | Receiver | Direction | Ligand | Primary rank | Eligible/planned omissions | Rank range |
|---|---|---|---|---:|---:|---|
| GSE136831 | Fibroblast | down | IL1A | 1 | 11/14 | 1–4 |
| GSE136831 | Fibroblast | down | IL1B | 2 | 11/14 | 2–5 |
| GSE136831 | Fibroblast | up | IL1A | 8 | 10/14 | 2–8 |
| GSE136831 | Fibroblast | up | IL1B | 6 | 10/14 | 2–8 |
| GSE136831 | __broad__ | down | IL1A | 2 | 22/22 | 1–3 |
| GSE136831 | __broad__ | down | IL1B | 1 | 22/22 | 1–3 |
| GSE136831 | __broad__ | up | IL1A | 8 | 22/22 | 7–8 |
| GSE136831 | __broad__ | up | IL1B | 7 | 22/22 | 4–7 |

![Ligand-target eligibility and fit](../../figures/ipf_ligand_target_eligibility_and_fit.png)

Tables: [primary candidate ranks](primary_candidate_rankings.csv), [target/expression eligibility](combined_eligibility.csv), [fixed source/receiver support](fixed_source_receiver_support.csv), [ranking stability](candidate_ranking_stability.csv), [prior receptor-panel coverage](prior_edge_panel_coverage.csv), [plotted values](figure_core_target_values.csv), [checks](validation.json). Per-cohort records preserve target counts, background coverage, individual omission status and full-data DE parity.
