# A1 study and assay map

Latest closure audit: [report](reports/EVIDENCE_CLOSURE_REPORT.md) and
[established analysis references](reports/ANALYSIS_REFERENCE_MAP.md). CD44's
count/sample identities are recovered from original SRA filenames; paired
contrasts and genotype interaction are complete. HPCS's biological labels and
confidence-abstention definition are verified. These updates supersede the
older identity/annotation caveats in the catalog below.

New primary follow-up: [Auyeung 2025](https://doi.org/10.1172/JCI184522),
[GSE243124](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243124), provides
day-10 epithelial scRNA after saline, bleomycin, epithelial IRE1 knockout or
KIRA8. Each group is pooled into one GEM library (2/2/3/3 contributing mice).
Replicated treatment inference requires demultiplexing or additional libraries.
The linked [GSE243129](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243129)
is neonatal Tgfbr2/hyperoxia. SAGE's current text still promises future
deposition; no new perturbation accession was recovered.

Catalog audit: 25 September 2026. Accession links lead to primary repositories;
paper links establish study context. GSM counts below are **library records,
not independent animals or donors**. Public metadata and listed files do not
establish that a usable matrix, sample crosswalk or independent replication has
passed inspection. The [metadata audit](metadata/README.md) records those limits.

## Direct regulatory measurements

| Study/resource | Deposited evidence | Role and limitation for A1 |
|---|---|---|
| [Kobayashi PATS study, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7461628/); [GSE141635](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141635) | 20 GSM: H3K4me3, H3K27ac, H3K36me3 and H3 in injured CTGF-positive epithelium versus homeostatic AT2; TP53/input in the transitional population; two replicate labels | Direct histone/TF evidence. Injured state, injury exposure and sorting are confounded; pooled tracks alone do not support replicated differential occupancy. |
| [Tsutsui et al., 2026](https://www.nature.com/articles/s41467-026-68909-z); [GSE289683](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289683), [GSE291333](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291333) | 12 + 12 GSM: iATC/iAT2/iAT1, two replicate labels, H3K27ac/H3 or H3K4me3/H3K27me3 | Direct modification comparisons in an iPSC model, not evidence for all human donors. Keep culture preparations and cell-line independence distinct. |
| Same study: [GSE290014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE290014), [GSE289846](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289846) | Three scATAC libraries; six scRNA libraries | Separate assays, not paired multiome. Baseline RNA culture day 7 differs from ATAC day 14; state/culture alignment must be explicit. |
| Same study: [GSE289676](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289676), [GSE289678](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289678), [GSE289679](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289679), [GSE289682](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289682) | Epithelial and fibroblast bulk RNA; 9/9/6/6 GSM. BLM ± p300/CBP inhibitors; TRF2 epithelial DOX ± CBP30; fibroblast PBS versus DOX coculture | Response/compartment comparisons after preparation IDs are resolved. GSE289682 is not an inhibitor contrast. |
| Same study: [PRJDB37980](https://www.ebi.ac.uk/ena/browser/view/PRJDB37980), [PRJDB37982](https://www.ebi.ac.uk/ena/browser/view/PRJDB37982), [PRJDB37983](https://www.ebi.ac.uk/ena/browser/view/PRJDB37983) | ENA lists 26/6/32 runs; RNA/CUT&Tag-related sequencing | Run existence verified; mark/treatment/pair mapping is incomplete. Sequence runs are not biological replicates. Do not infer a mapping from the order in the article's data statement. |
| [Zhou et al., 2021](https://link.springer.com/article/10.1186/s12864-021-08152-6); [GSE150527](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150527) | 48 GSM: histone marks, FAIRE, CTCF, RNA and WGBS across human AT2-to-AT1 culture | Normal differentiation comparator. RNA has three donor labels, selected histones two, WGBS only one donor at D0/D4/D6. No population DMR inference. D4 culture is not a purified DATP/PATS sample. |
| [Marjanovic HPCS study, 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7745838/); [GSE154966](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154966) | Eight bulk ATAC GSM and deposited counts/merged peaks; four apparent TIGIT-positive/negative source blocks | First processed-count candidate. Resolve the pooled `106621_106642` source and independence of all blocks. Cancer context and Kras/Trp53 genotype constrain interpretation. |
| Same study: [GSE154965](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154965), [GSE154978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154978) | Two scATAC GSM with fragments; six sorted/transplant scRNA GSM | scATAC titles and file suffixes conflict; label-dependent analysis is on hold. Transplant samples require source/outcome pairing, not assumed matching. |
| Choi DATP-associated ATAC: [GSE144598](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144598) | Four ATAC GSM; two pooled processed BigWigs; SRA relations listed | IL1R1-positive AT2 versus control AT2, not a purified DATP comparison. Raw-read reprocessing is potentially possible; the processed-track limitation is not absence of raw data. |
| Existing local multiome: [GSE310539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310539), [GSE247130](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247130) | Eight/twelve RNA+ATAC GSM correspond to four/six wells | Paired RNA/accessibility and depth controls. Preserve the repository's corrected well/genotype map; pooled animals and one well per condition do not support animal-level effects. Same laboratory is not independent replication. |
| [GSE309751](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309751) | 14 bulk ATAC GSM; post-viral time points and Kras perturbation; processed peak calls and SRA relations | Two 49-day control GSM have contradictory treatment fields. Peak-call presence alone is not quantitative accessibility or same-cell reopening. |
| [GSE327686](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327686), [GSE327565](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327565) | Eight RNA/ATAC GSM in four multiome conditions; separate 14-GSM organoid bulk-RNA series | GSE327686 series design describes bulk cultured AT2 despite whole-lung RNA/ATAC GSM; resolve before use. GSE327565 must not be labelled multiome. |

## Lineage, dynamics, protein and function

| Study/resource | Distinguishing measurement | Planned use and limitation |
|---|---|---|
| [Strunz et al., 2020](https://www.nature.com/articles/s41467-020-17358-3); [GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259) | Injury chronology, epithelial enrichment and published lineage/velocity analyses; 60 GSM across preparation types | Reconstruct time-associated state paths with sample-aware uncertainty. Resolve preparation overlap. Ordinary count matrices cannot supply splicing layers or measured lineage histories. [Author code](https://github.com/theislab/2019_Strunz). |
| [Auyeung et al., 2022](https://pubmed.ncbi.nlm.nih.gov/35170357/); [GSE190821](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190821) | Krt19 pulse/chase and AGER endpoint under IRE1α inhibition; separate day-7 RiboTag RNA; 48 GSM representing epithelial and whole-lung libraries from 24 mice | First RNA contrast completed: five KIRA8 versus five vehicle mice with batch/sex adjustment. Axum8 is a separate antibody control. RiboTag expression does not identify traced DATPs, fate transitions or chromatin changes. |
| [Chan et al., 2026](https://www.nature.com/articles/s41586-025-09985-x); [GSE277777](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277777) | HPCS tracing/ablation, GEX and HTO, processed objects; 32 GSM | Distinguish lineage origin, chase and perturbation from a shared RNA programme. Some animals are pooled. HTO is sample hashing, not a protein-phenotyping panel. Existing A11 projection is not a reanalysis of lineage outcomes. |
| [Krt8 perturbation study](https://www.jci.org/articles/view/165612); [GSE223302](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223302) | 24 bulk-RNA GSM, WT/KO time course | Time-by-genotype response after replicate verification. Links molecular state to perturbation; bulk expression alone cannot prove irreversible fate. |
| [Aberrant intermediate epithelial cells, 2025](https://www.nature.com/articles/s41467-025-63735-1); [GSE273123](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273123) | 16 CD44-positive/negative sorted AT2 RNA GSM, WT and SP-C mutant; published fibroblast-response assays | Protein-enriched populations provide an orthogonal selection axis. Four candidate pairs per genotype require animal-ID confirmation. This is not mass-spectrometry proteomics. |
| Same study: [GSE273122](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE273122), [GSE272862](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272862) | Four mouse scRNA GSM; two iAT2 genotype/corrected GSM | Mouse titles and genotype characteristics are reversed relative to each other: hold genotype tests. One library per human genotype limits inference. |
| [Ke developmental transition study](https://doi.org/10.1016/j.devcel.2024.11.017); [GSE254356](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE254356) | Three embryonic sorted epithelial scRNA GSM and published lineage evidence | Developmental reference coordinated with A5; two-versus-one preparation design, not an adult repair validation cohort. |
| [Konkimalla et al., 2023](https://pmc.ncbi.nlm.nih.gov/articles/PMC10762634/); [GSE218665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218665), [GSE218666](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE218666), [GSE235212](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235212) | BHT injury, AT1 ablation and AT2 Mdm2 deletion; 2/4/1 GSM | Tissue topology, mechanics and perturbation anchors. Deposits alone do not encode every imaging or fate endpoint; one sample per time is descriptive. Shared external controls are not new independent cohorts. |
| [IPF regenerating-niche atlas, 2025](https://www.nature.com/articles/s41467-025-61880-1); [Zenodo 10930946](https://zenodo.org/records/10930946) | 33-channel imaging mass cytometry; 22 MCD files (~15.55 GB) | Cell proteins, morphology and neighbourhoods. Prefer [processed cell data](https://mdv.molbiol.ox.ac.uk/projects/mdv_project/7430) after donor/ROI/mask/channel audit; 22 files are not 22 donors. Pathology stages are cross-sectional. [Code](https://github.com/LingPeiHo/Ho_Taylor_Byrne-SpOOx-2.0). |
| [MUC5B regional proteomics study, 2025](https://insight.jci.org/articles/view/189636); [PXD058626](https://www.ebi.ac.uk/pride/archive/projects/PXD058626) | Laser-capture microdissection mass spectrometry of defined lung regions, including epithelium over fibroblastic foci | Regional protein corroboration with donor/genotype blocking. Project metadata verified; processed abundance matrix and donor-region map remain unverified. Mixtures are not a purified transitional-cell proteome. |
| [SAGE Perturb-seq preprint, 2026](https://www.biorxiv.org/content/10.64898/2026.03.13.711474v1) | Time-resolved TF perturbation and repair/pathological epithelial branches | Literature anchor. Version 1 promises deposition; no new public accession was verified. Human reference accessions cited there are not its new perturbation data. |

No matched 3D-contact assay or directly relevant epithelial CITE-seq deposit
passed this bounded audit. These remain coverage gaps, not claims that no such
study exists. Re-query them when the first executable contrasts are frozen.

The [lineage audit](LINEAGE_AUDIT.md) distinguishes origin tracing, state pulse/
chase, experimental conversion and perturbation across eight anchor studies.
The [first-batch report](reports/FIRST_BATCH_REPORT.md) supersedes planning-only
status for the executed source tables without promoting the other catalog entries.

The [second-batch report](reports/SECOND_BATCH_REPORT.md) adds verified native
CHM13 histone profiles from GSE289683/GSE291333 (two preparations, one B2-3 line),
GSE150527 hg19 domain-overlap context (one donor), and IRE1α omission/batch
sensitivities. The 23-locus histone table includes H3 and window sensitivities;
these are descriptive comparisons, not independently replicated state classes.
PATS tracing timing is resolved, while deposited histone scaling remains held.
The expanded GEO/ENA audits do not resolve TIGIT pools or CD44 matrix identities.
GSE277777 observation metadata now supports a descriptive 22-source descendant
composition table; mouse/pool independence, current mScarlet and one age conflict
remain unresolved. No paired multimodal assay is implied by these separate studies.
