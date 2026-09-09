# Pipeline framing from the HLCA: which criteria this repository adopts, adapts, or declines

Status: **proposal, pending owner retain/reject. Nothing in this file has
been run.** It reads the HLCA's decision criteria
([`README.md`](README.md), [`integration_benchmark.json`](integration_benchmark.json))
against the pipeline that actually ran here
([`../../docs/PIPELINE_AS_RUN.md`](../../docs/PIPELINE_AS_RUN.md),
[`../../docs/ANALYSIS_RATIONALE.md`](../../docs/ANALYSIS_RATIONALE.md)) and
says, criterion by criterion, what a future version of this pipeline should
do differently and why. The first concrete trial is in
[`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md).

## 1. Where the two designs agree and where they part

| Design question | HLCA (Sikkema 2023) | This repository today | Consequence for the future pipeline |
|---|---|---|---|
| What is a "batch" | the *dataset*, chosen by a PCexpl test against shuffled labels; never the donor, precisely so that inter-individual variation survives integration | mouse: sample = biological group, so no correction; human: *donor* is the batch, Harmony by logged override | The human choice is what the HLCA would call over-correction. It remains defensible only because the three donors are healthy replicates with no condition variable, and that reasoning must stay written next to the decision. A future run should apply the PCexpl split test before deciding what "batch" means, and log it. |
| How the integration method is chosen | a scored benchmark on a subset (12 scIB metrics, 0.4 batch / 0.6 bio) plus rare-cell recall and precision | a design argument plus one replicate-mixing statistic (within-group enrichment 1.374) plus a donor-driven-clustering check | Adopt a pre-registered metric panel so the choice is measured, not argued (section 2.1). The design argument stays as the reason "none" is always a candidate. |
| Labels during integration | inside the model (scANVI is semi-supervised on harmonised labels) | held out of every fitting step; used only afterwards as an answer key | Keep held-out. Do **not** adopt scANVI with deposited labels: the median cluster purity of 0.947 would become circular. Admissible variants are scVI (unsupervised) or scANVI seeded only with this repository's own blind annotations. |
| Clustering | nested Leiden: r 0.01 at level 1, r 0.2 below, k 30/30/15/10/10 | Leiden 0.3 once on the atlas, then subset-and-recluster for the focused analyses | Same motif. Adopt the explicit level naming (1.2 = third child of cluster 1) and record k and r per level in decisions.json. |
| Cluster quality screen | label entropy > 0.56 flags disagreement or doublets; donor entropy < 0.43 flags donor-private clusters | the count of clusters with > 75% of cells from one donor; purity against deposited labels | Adopt both entropies (section 2.2). The single-donor count is a coarse version of donor entropy with an implicit threshold; make the threshold explicit and dataset-specific. |
| Marker genes | per-sample, per-cell-type pseudo-bulk t-test with cross-sample consistency filters (80% of pseudo-bulks, 50% of cells, 20% out-group) | per-cell Wilcoxon, documented as a ranking device rather than a test | Adopt the pseudo-bulk route. The statistical unit becomes the animal or donor, which is the owner's standing rule. |
| Covariate or condition effects on genes | pseudo-bulk mixed model with dataset as random effect; at least 40 samples per cell type; VIF > 5 excludes | binned programme summaries along pseudotime; no formal DE | By the HLCA's own standard a 25-sample time course cannot attribute gene-level variance to covariates; any such claim stays "Exploratory". Time-point contrasts are still possible with the animal as the unit, but with the sample count stated next to every P value. |
| Annotating new data | map to a reference (scArches) and transfer labels with an uncertainty score; abstain above 0.2 | blind marker panels, then post-hoc grading against deposited labels | Adopt uncertainty as a first-class output for the human series; "unknown" is a legitimate annotation and is where novel states hide. |
| Ambient RNA | not corrected in the core (one dataset arrived pre-corrected) | not corrected (raw droplet matrices absent for the mouse series) | No change; the caveat travels with both. |
| Doublets | removed post hoc as high-entropy clusters co-expressing lineage markers | Scrublet per capture before merging | Keep Scrublet; add the entropy screen as an independent cross-check on what Scrublet missed. |
| Gene space across datasets | harmonised to HGNC symbols via an Ensembl release graph (79 to 107); duplicates summed; 2,000 consensus HVGs; missing HVGs zero-filled in queries | Ensembl-ID intersection for the human series; symbols for the mouse series | Zero-filling missing HVGs is what makes mapping to the HLCA possible at all; note that it also silently lowers the query's information content. |

## 2. Criteria to adopt, pre-registered

Freeze these before the next run touches data. Each item names the rule, the
value, the datasets it applies to, and the key under which it must be logged.

### 2.1 Integration-selection panel

- **Candidates.** `none`, Harmony, scVI. scANVI only with blind labels.
  `none` is always a candidate because in the mouse design the batch and the
  condition are the same variable.
- **Metrics.** The four batch metrics (PCR batch, batch ASW, graph iLISI,
  graph connectivity) and the eight bio metrics (NMI, ARI, cell-type ASW,
  isolated-label F1 and silhouette, graph cLISI, cell-cycle conservation, HVG
  conservation). Deposited author labels are used for the label-dependent
  metrics as **evaluation only**; they never enter fitting.
- **Score.** Category means, overall 0.4 batch + 0.6 bio, exactly as the
  HLCA, so the numbers are comparable with the published benchmark.
- **Adoption rule.** A method is adopted only if (a) it beats `none` on the
  overall score and (b) it does not lower recall for the rare populations
  named in advance: ionocyte, tuft and neuroendocrine for the human series;
  Krt8-high transitional epithelium and the injury-associated capillary state
  for the mouse series.
- **Log.** `decisions.json: integration_benchmark.{candidates, metrics,
  weights, scores_per_candidate, rare_cell_recall, chosen, rule_applied}`.

### 2.2 Cluster entropy screen

- **Label entropy.** Shannon entropy (natural log) of deposited-label
  fractions per cluster, labelled cells only, NA when fewer than 20% of the
  cluster is labelled. High above 0.56. Every flagged cluster must be either
  re-annotated at a coarser level or tested as a doublet cluster
  (co-expression of two lineage panels). Cross-check against the three
  clusters already marked CONTRADICTED in the mouse annotation table.
- **Donor entropy.** Same formula over animals or donors. The HLCA threshold
  (0.43) encodes a 95/5 split among 107 donors and does **not** transfer; the
  threshold is recomputed per dataset as the entropy of a cluster with 95% of
  cells from one donor and 5% spread evenly over the remaining n minus 1. For
  the human series (n = 3) that is about 0.23; for the mouse it depends on the
  stratum (next point).
- **Mouse stratification.** In the time course a cluster private to two
  animals may simply be a state that exists at one time point (the injury
  capillary state peaks at 25 dpi). Donor entropy for the mouse is therefore
  computed *within* time point, so that it flags animal-private artefacts and
  not time-point biology. This is the single most important adaptation.
- **Log.** `decisions.json: cluster_entropy.{label_threshold, donor_threshold,
  donor_threshold_formula, stratification, flagged_clusters, actions}`.

### 2.3 Pseudo-bulk markers

Per animal or donor, per cluster; at least 10 cells (3 for populations under
100 cells); t-test versus the rest; keep genes expressed in at least 80% of
pseudo-bulks, in at least 50% of cells on average, and in at most 20% of
out-group pseudo-bulks; relax in the documented rounds only if no marker
survives. Cluster marker tables then carry a "supported by N of M animals"
column instead of a cell-level P value.

### 2.4 Uncertainty-bearing annotation for the human series

Map GSE178360 to the HLCA core with scArches; transfer labels with k = 50;
label cells with u > 0.2 as unknown; compare transferred labels with this
repository's blind annotations and with the deposited author labels; report
the confusion at HLCA level 3 and 4. The AT0 population is the named
question: HLCA carries AT0 (n = 1,440) and pre-TB secretory (n = 4,393) as
explicit classes, so label transfer tests the AT0 call by a route independent
of marker gates and of Scrublet.

### 2.5 Sample-count gate for covariate and condition claims

Any gene-level claim attributed to a covariate or condition carries the
number of samples behind it; below 40 samples per cell type the claim is
labelled Exploratory in the README claims table, following the HLCA's own
rule. This applies to every mouse time-point contrast.

## 3. Criteria declined or adapted, with the reason

| Criterion | Decision | Reason |
|---|---|---|
| Supervised scANVI on deposited labels | declined | circular with the held-out validation that gives this repository its credibility |
| Donor as batch for the human series | kept, relabelled | it is over-correction by the HLCA's definition; acceptable for three healthy replicates with no condition variable, and only while that is true |
| Mapping the mouse series to the HLCA | declined | human reference, adult only; the authors say cross-species mapping needs method development |
| 6,000-gene FULL feature set | not needed | the winning configuration used 2,000 HVGs |
| CCF anatomical score | not applicable | GSE178360 is distal parenchyma throughout (CCF about 0.97); no gradient was sampled |
| sc-LDSC, CIBERSORTx | out of scope for now | need GWAS summary statistics or bulk cohorts; neither is part of the portfolio question |
| Fixed donor-entropy threshold 0.43 | adapted | depends on donor count; recomputed per dataset and per stratum |

## 4. The future pipeline skeleton with the HLCA gates inserted

```
deposited counts
  -> per-sample MAD QC (existing) -> Scrublet per capture (existing) -> merge
  -> [G1] batch-unit test: PCexpl of candidate technical covariates vs 10 shuffles; split if > 1.5 SD
  -> consensus HVG ranking across samples (cell_ranger flavour per sample) -> 2,000
  -> [G2] integration panel on {none, Harmony, scVI}: 12 scIB metrics, 0.4/0.6, rare-cell recall; labels for evaluation only
  -> nested Leiden with k and r recorded per level
  -> [G3] entropy screen: label entropy > 0.56 => re-annotate or doublet test; donor entropy (within stratum) < threshold(n) => flag
  -> blind annotation (existing) -> post-hoc grading against deposited labels (existing)
  -> [G4] pseudo-bulk markers with per-animal support counts
  -> [G5] condition or covariate claims only from sample-level models; < 40 samples => Exploratory
  -> human series only: [G6] scArches to the HLCA core; u > 0.2 => unknown; confusion vs blind and deposited labels; AT0 check
```

Every gate writes its inputs, thresholds and outcome to decisions.json so
that the per-dataset reports and the pipeline record can be regenerated from
artefacts, as now.

## 5. Environment feasibility: established on 2026-09-09 by trial S2

PyTorch 2.14.0 (CPU) and scvi-tools 1.5.0.post1 were installed into the
emulated interpreter with the owner's authorisation, all as prebuilt
wheels; the HLCA model converted and the mapping ran at about 40 seconds
per epoch for 28k cells. The `scarches` package itself does not import
with anndata 0.13 and is not needed. The facts below were written before
that run and are kept for the record:

- scVI, scANVI and scArches need PyTorch. The scib-metrics package needs JAX.
  Neither is in `analysis/requirements.txt`.
- This machine is Windows 11 on ARM64 with no C compiler; JAX publishes no
  Windows ARM64 wheels. The emulated x86-64 interpreter at
  `.venv-x64/Scripts/python.exe` is the only route, and PyTorch CPU wheels for
  win_amd64 exist. Training speed under emulation is unknown. Any install
  must be attempted with uv against that interpreter and refused if it would
  require compilation.
- The original scib package computes graph LISI through a compiled helper;
  scib-metrics reimplements PCR batch, batch silhouette, kNN iLISI and cLISI,
  graph connectivity, NMI, ARI and isolated-label metrics in Python.
  Cell-cycle conservation and HVG overlap are short scanpy computations. The
  full 12-metric panel is therefore reachable without R.
- The HLCA reference model on Zenodo (10.5281/zenodo.7599104) was built with
  scArches 0.3.5 and scvi-tools 0.8.1; loading it with current versions is a
  compatibility risk that must be tested before the mapping trial is planned
  in detail. The HLCA mapping repository provides a tutorial for this.

## 6. Open decisions for the owner

1. Accept "deposited labels for evaluation only" inside the integration
   panel, or keep integration selection label-free and use labels only for
   the final grade?
2. Accept the per-dataset, per-stratum donor-entropy threshold as written in
   section 2.2?
3. (Answered 2026-09-09: the environment work was authorised and trial S2
   ran.) Decide how the human AT0 headline is re-worded in the light of
   S2, and whether cluster 22 is corrected to mast cells.
4. If the panel prefers `none` for the human series, does the biological
   argument for Harmony (one cell type fragmenting into donor-private
   clusters) still override it, and on what recorded evidence?
