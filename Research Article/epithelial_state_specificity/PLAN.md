# Epithelial-state specificity: ES1

Frozen before ES1 scoring on 2026-09-22. This is a descriptive follow-up to
the repository audit, not a preregistered confirmatory experiment. Earlier
two-marker results and the source papers are already known.

## Question and boundaries

Do complete **available marker panels** beyond Cldn4/Krt8 show the same
within-library association with the two-transcript group in neonatal control,
adult control, injured control and Cebpa-mutant epithelium? Compare AT2 identity,
early and late AT1 markers, and full published pathway modules separately.
Neither a marker-panel score nor reduced AT2 identity identifies cell fate,
chromatin closure, productive repair, arrest, or reversibility.

The current repository has Choi's short DATP marker panels and the M branch's
combined DATP/PATS marker panel. These are **not full DATP/PATS/ADI signatures**.
Before any ES1 score was computed, the Strunz 2020 Supplementary Data 3 table
was recovered through Europe PMC's supplementary archive. Its `cell_types_2`
sheet supplies all 400 reported positive markers each for Krt8 ADI, AT2 and
AT1. ES1 freezes each complete published list and a version excluding Cldn4,
Krt8, Sftpc and Cebpa. The first three are labeling genes; the last is removed
consistently across genotypes. The source itself caps each list at 400 genes:
these are complete **reported marker lists**, not exhaustive biological
programs. The initial local-only freeze is preserved in
`modules_initial_no_state_signatures.json`; it was not scored.

Full DATP/PATS signatures and a developmental maturation signature remain
missing. `freeze_modules.py` records this source gate explicitly; the complete
MSigDB Hallmark controls and the adult AT2 marker list do not fill it.

## Locked implementation

1. Use the four cached multiome count files (10 wells); use the corrected
   GSE247130 suffix map retained in C116. Every experimental condition has one
   pooled library. Two source deposits share a laboratory and are consistency
   checks, not independent replication.
2. Keep nuclei with at least 2,000 RNA UMIs. Remove nuclei whose summed
   Scgb1a1/Scgb3a2/Foxj1/Krt5 exceeds Sftpc/Sftpb/Lamp3 (the earlier relative
   airway screen). This operational compartment may retain contaminants;
   publish retention and gene coverage rather than calling it pure AT2.
3. Jointly downsample the selected gene counts without replacement to a
   2,000-UMI total budget using sequential hypergeometric draws. Unselected
   genes remain the residual category. Seed 17 is primary; seed 29 is a
   technical sensitivity check, never a biological replicate. RNA features
   with duplicate gene symbols are summed before scoring; none is chosen
   arbitrarily. This input-format correction preceded any score output.
4. Label the two-transcript group by sampled Cldn4 > 0 and Krt8 > 0.
   Reference: sampled Sftpc > 0 and not two-transcript positive. Score each
   module as the mean of binary gene detections. Report a label-free variant
   of the DATP panel and exclude Cebpa from the common AT2 instrument.
   The raw DATP panel shares its two defining genes with the label, so is a
   circular positive control only. No classifier accuracy is claimed.
5. Publish scores by library/group, every gene's detection, group sizes,
   coverage and module overlap. An individual group needs at least 30 cells
   for an interpretable score; smaller groups stay visible with an explicit
   flag. Require at least 80% gene coverage. Scores and differences are
   descriptive, with no cell-as-replicate tests, null-derived confidence
   intervals or 'not significant means absent' conclusion.
6. The primary quantity is labelled-minus-reference detection points within
   each well. Between-well differences are descriptive contrasts: P9 versus
   seven-week **within control genotype**; adult injury versus adult control
   **within genotype**; mutant versus control **within stage**. They do not
   estimate independent injury, age or genotype effects. No neonatal-injury
   arm exists, so no age-by-injury interaction is identifiable.
7. Freeze modules before scoring and hash their inputs. Never select genes
   from this run. Generate summaries from output tables. Preserve historical
   M-series files unchanged.

## Independent confirmation gate

GSE262927 is a different laboratory/study from the two multiome sources and
has identifiable mice. Check its cached epithelial object and metadata for
adult alveolar cells, raw integer counts, sample identities, and evaluable
groups. Use only the deposited annotated cohort and its AT2, AT1, AT1_AT2 and
Alveolar_transitional labels (exact metadata vocabulary checked before scoring)
to audit population coverage; do not reopen its displaced trajectory.
If counts support it, apply exactly the frozen panels, depth and two-marker
rule and report per-animal within-sample effects and eligible-animal counts.
The cohort was analyzed previously by this repository, so it is external to
the discovery studies but **not an untouched held-out cohort**. Its different
injury and assay limit transportability. It has no neonatal arm and cannot
independently confirm developmental-program reuse, cellular reversibility or
repair-versus-fibrosis discrimination. If an arm fails the 30-cell group floor,
report insufficient coverage, not biological absence.

## Planned figure

A heatmap of within-library differences in the independent marker panels;
technical-seed agreement; and any evaluable external per-animal estimates.
No embedding, trajectory or fate labels inferred from the same module score.
