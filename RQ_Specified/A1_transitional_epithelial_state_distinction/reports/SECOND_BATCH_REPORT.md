# A1 second batch: regulatory profiles and perturbation stability

25 September 2026. Continuation authorized by the owner's request to resume A1.
This is a source-informed follow-up to the [first batch](FIRST_BATCH_REPORT.md),
not unseen-data confirmation. The primary IRE1α analysis, first render, executed
source archive, frozen loci and thresholds remain unchanged. A1 remains open:
these studies do not compare DATP, PATS, ADI, ABI and HPCS side by side.

## What this batch establishes

The induced alveolar-cell experiment provides direct, native-assembly histone
profiles at the 23 frozen loci. Several locus effects agree between its two
preparations, but some directions depend on H3 normalization or promoter width.
The IRE1α RNA response has relatively stable marker directions and unstable
significance. Normal-differentiation methylation-domain overlap supplies a
separate genomic reference. Neither these results nor the PATS tracing endpoint
establishes a universal epigenetic state taxonomy or chromatin causality.

## J1 — IRE1α sensitivity verified

All ten distinct omissions retain nine mice (4 versus 5), a full-rank four-column
design and five residual degrees of freedom. S061 has three mice per arm, rank
three and three residual degrees of freedom. Each fit retains exactly the
primary 14,811 genes. Output/input hashes, original counts/sample MD5 values,
normalization sample order, full-table/focus joins, all-gene BH values and
eligible pathway memberships were checked. No mouse was removed from the primary
analysis. See the [verification record](second_batch_verification.json) and
[fit diagnostics](../tables/ire1_stability_2026-09-25/fit_diagnostics.tsv).

| Frozen marker | Primary log2 FC | Omission minimum–maximum | Same direction / 10 |
|---|---:|---:|---:|
| Itgb6 | −0.811 | −1.016 to −0.683 | 10 |
| Krt8 | −0.566 | −0.720 to −0.436 | 10 |
| Krt19 | −0.604 | −0.775 to −0.430 | 10 |
| Cldn4 | −0.981 | −1.187 to −0.806 | 10 |
| Ager | +0.255 | −0.121 to +0.447 | 9 |
| Sftpc | +0.659 | +0.225 to +0.831 | 10 |
| Cdkn1a | −0.903 | −1.055 to −0.806 | 10 |
| Krt7 | −0.526 | −0.592 to −0.497 | 10 |

Direction is KIRA8 minus vehicle. Ranges are omission sensitivity ranges, not
confidence intervals. All eight markers failed primary all-gene FDR 0.05.
The four primary FDR hits remain a separate **post-selection** group in the
[complete focus table](../tables/second_batch_verified/ire1_focus_stability.tsv).
They retain their primary direction in all omissions, but pass FDR in 9, 10, 4
and 6 omissions respectively in the table's order. This is not a new discovery
family or independent replication.

Whole-family significant-gene counts range from 1 to 264 across omissions,
versus 4 in the primary model; omission Pearson effect correlations are
0.928–0.973. S061 yields 112 discoveries but is a smaller, different sample subset.
S135 has only two mice per arm and receives descriptive normalized-CPM ratios
with the specified 0.5 offset, without p-values. Krt19 and Sftpc directions in
S135 differ from the primary estimates. Unequal significance across batches is
not a treatment-by-batch interaction.

TGF-β points downward throughout but passes FDR 0.05 only when omitting mouse
148 (q=0.04490; primary q=0.12463). The original 47 TGF-β and 94 UPR genes and
estimated-correlation CAMERA method are retained; terminal UPR remains below
its original coverage rule. All tests are in the
[pathway comparison](../tables/second_batch_verified/ire1_pathway_comparison.tsv).
The strongest permissible interpretation is directional follow-up, not a
confirmed pathway or a replacement primary analysis. Day-7 epithelial RiboTag
RNA can reflect cell composition and translation-associated changes; it is
not sorted DATP RNA or the separate day-14 traced-AGER endpoint.

## J2 — direct histone profiles verified

GSE289683/GSE291333 contribute 24 CPM tracks: three induced states, two CUT
labels, and H3K27ac, H3K4me3, H3K27me3 and H3. The source cell line is B2-3.
The 24 headers agree, and every represented chromosome length matches the
[official native CHM13 size table](https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.chrom.sizes).
No hg19/mm10 coordinates enter this comparison. Annotation uses exact symbols,
primary chromosomes and strand-aware gene-boundary TSSs; it does not identify
the dominant transcript's TSS.

Hand-calculated tests cover clipping, half-open boundaries, gaps, zero signal,
overlaps, nested/adjacent domains, negative-strand TSSs and ambiguous/missing
symbols. Seventy-two saved promoter measurements (all tracks at SFTPC, AGER
and CDKN1A) agree with independent base-resolution reader values. All 1,104
signal rows, matching state/CUT H3 controls, ratio eligibility and contrast
arithmetic were checked. A ratio needs ≥80% coverage in both tracks and positive
mark/H3 means; uncovered bases contribute zero to the window mean, while
coverage remains separately visible. No pseudocount rescues a held ratio.

| Mark | Concordant preparations / eligible locus–contrast pairs at ±1 kb | Eligible preparation–contrasts at ±1 kb / 92 | Same direction at ±1 and ±5 kb / eligible at both | H3 reverses raw-mark direction at ±1 kb |
|---|---:|---:|---:|---:|
| H3K27ac | 36 / 45 | 91 | 78 / 91 | 26 / 91 |
| H3K4me3 | 36 / 42 | 86 | 75 / 86 | 15 / 86 |
| H3K27me3 | 38 / 42 | 87 | 78 / 87 | 18 / 87 |

There are 46 possible gene–contrast pairs per mark, each with two preparations.
These are correlated comparisons, **not independent samples**. Eligible ±5-kb
preparation–contrasts number 92, 88 and 92 respectively. Full values for **all
23 loci, both contrasts, both windows and both preparations**, including held
ratios, are in the [window/H3 sensitivity table](../tables/second_batch_verification/histone_window_and_H3_sensitivity.tsv).
The [held-ratio table](../tables/identity_audit/held_histone_ratios.tsv) lists
every ineligible comparison; [signal tables](../tables/direct_marks_2026-09-25/histone_signals.tsv)
retain the underlying means and coverage. No locus was omitted for discordance.

The predefined CDKN1A promoter has higher H3K27ac/H3 in iATCs versus iAT2
(+0.84, +1.27 log2) and lower H3K27me3/H3 (−1.38, −0.89). Raw H3K27ac also rises
(+0.82, +0.49), so the acetylation direction is not created by H3 adjustment,
although its magnitude changes. SFTPC H3K27ac/H3 is lower than iAT2 in both
preparations (−1.51, −2.40); KRT8 is higher (+0.40, +0.39).

Counterexamples matter. CLDN4 raw H3K27ac is higher than iAT2 (+0.36, +0.41),
whereas its H3-adjusted values are lower (−0.94, −0.89). KRT8 also changes from
higher raw signal to lower H3-adjusted signal relative to iAT1. CLDN4
H3K27me3/H3 versus iAT2 disagrees between preparations (+0.85, −1.27), and
CDKN1A H3K4me3/H3 versus iAT1 disagrees (−0.69, +0.34). AGER H3K27me3 ratios
are held for CUT2 at the narrow promoter. These observations argue against
treating all transitional markers as one uniform regulatory programme.

The result is a locus-specific descriptive lead in one iPSC line and two
preparations per state. Neither CPM nor mark/H3 is absolute or spike-in-scaled
occupancy. No global acetylation conclusion, population p-value, single-cell
bivalency, causal mechanism or independent state classifier is claimed.

## J3 — normal-differentiation methylation reference verified

All 414 domain intersections were independently reconstructed with per-base
boolean unions. The source is one donor at D0/D4/D6, separate hg19 annotation
and ±1-kb promoters. Both coordinate interpretations remain in the
[source table](../tables/direct_marks_2026-09-25/methylation_domain_overlap.tsv).
Author starts at 1 suggest one-based inclusive intervals; metadata do not
conclusively establish this convention. The alternative BED interpretation
changes an overlap fraction by at most 0.0005 (one base in a 2,000-base window).

UMR/LMR/PMD overlap is genomic context, not promoter-wide CpG methylation.
The domain classes are displayed separately and need not form a mutually
exclusive partition of a promoter. Segment methylation summaries exist in
some source tables but were not used to estimate promoter methylation here.
No replicated DMR test, purified transitional-state assignment or injury/tumour
fate conclusion follows from this normal-culture reference.

## J4 — PATS timing resolved; deposited histone scaling held

The [published PATS article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7461628/)
specifies a Krt19-CreER pulse seven days after injury and a day-12 harvest:
a nominal five-day interval, not an estimated transition rate. The original
ED4 reconstruction remains three mice per marker, with separate denominators
and undefined zero-denominator controls. The numerical first batch is unchanged.

The histone experiment sorts CTGF-positive cells on day 12; TP53 ChIP uses a
distinct day-8 sort with Sftpc lineage labelling. The methods describe
H3-normalized peak calling. The cached generic MintChIP pipeline calls HOMER
`makeUCSCfile` for bedGraph output, but neither that generic code nor the peak
method establishes each deposited track's scaling factors and paired controls.
Quantitative between-condition histone comparisons remain held pending that
contract. Four bedGraphs were acquired; the two larger CTGF-positive mark files
remain outside the frozen 256-MiB/file ceiling. Additional download cannot
resolve the normalization question by itself.

## J5–J6 — identity audit complete, inference held

The [GEO/ENA audit table](../tables/identity_audit/GEO_ENA_sample_audit.tsv)
contains all 24 relevant libraries, linked sample/experiment identifiers and
source URLs. The fresh audit inspected 48 ENA XML documents.

**TIGIT ATAC.** The published article describes four mouse replicates and a
paired `~ Mouse + Tigit_status` analysis. The deposited count columns map to
four named source blocks, but sample/experiment XML repeats `106621_106642`
without an explicit pool-membership/non-overlap statement. A compound alias
does not prove a pool, nor establish independence from YY1181, YY1916 and
106623. Required input: each positive/negative column's actual animal or
disjoint pool membership and pairing. No paired inference was fitted.

**CD44 RNA.** The 16 matrix columns use R26_8, R26_7, R26_2, R26_1, OG7, OG6,
OG5 and OG3, each in positive/negative sorts. GEO/ENA identifies WT1–4 and
mutant1–4, but ENA experiment/library aliases contain GSM labels only.
All 16 submitted-filename fields are empty in the
[remaining raw-file audit](../tables/identity_audit/CD44_submitted_file_audit.tsv).
The source workbook and published code likewise provide no exact crosswalk.
Required input: matrix column → GSM → genotype → mouse → CD44 sort.
The separate GSE273122 inconsistency was not imported into GSE273123.
No genotype was guessed and no significance difference was called an interaction.

## J7 — descendant source composition recovered; biological-unit inference held

The [pinned author code](https://github.com/dbetel/HPCS_LUAD/tree/b52d53c984e21d3bb3a163041fdb3f56b54c19c0)
identifies the Figure-2 input as `combined_data.h5ad`. A bounded metadata-only
read of the corresponding GEO object recovered 28,402 observation rows using
7,405,568 transferred bytes in 113 byte-range requests (158.61 seconds), from
an 8,045,265,618-byte object. No expression matrix or complete sequencing archive
was read. The [recovery record](hpcs_metadata_recovery.json) contains the remote
ETag, modification time, every fetched range/hash and the cached metadata hash.

Under a [separate source-informed contract](../config/hpcs_descendant_reconstruction.json),
the reconstruction retains 5,333 author-designated traced cells, 22 source
labels and all eight RNA-state categories. Each source label is verified against
the corresponding notebook's stored hash-feature output and its GEO library.
Group totals, 11 early HPCS fractions and four Hopx means reproduce saved
notebook outputs. The [run record](hpcs_source_composition_run.json),
[source manifest](../tables/hpcs_source_composition/source_manifest.tsv) and
[all numerators/denominators](../tables/hpcs_source_composition/source_state_counts.tsv)
make the reconstruction reviewable.

| Driver / deposited group | Chase | Source labels | Retained cells | Mean source HPCS fraction | Cell-pooled HPCS fraction |
|---|---:|---:|---:|---:|---:|
| Slc4a11 / 6wk_3d | 3 days | 5 | 586 | 93.08% | 92.83% |
| Slc4a11 / 8wk | 14 days | 4 | 742 | 24.00% | 23.99% |
| Slc4a11 / 12wk_3d | 3 days | 6 | 877 | 67.86% | 66.82% |
| Slc4a11 / 14wk | 14 days | 3 | 991 | 43.73% | 31.79% |
| Hopx / Hopx_12wk3d | 3 days | 2 | 1,226 | 4.77% | 4.81% |
| Hopx / Hopx_12wk14d | 14 days | 2 | 911 | 4.32% | 3.84% |

The 14wk Slc4a11 group illustrates why denominators matter: one source supplies
811 of 991 cells, and the equal-source mean differs from the cell-pooled result.
Both are descriptive. The [complete group table](../tables/hpcs_source_composition/group_descriptive_summary.tsv)
also retains all other categories, source ranges and zero counts. These are
descendant **RNA-label compositions among retained trace-sorted cells**, not
new measurements of cell production, clone size or transition rates.

The original broad data-availability hold is narrowed, but biological inference
still needs an explicit mouse/pool-membership crosswalk. Do not promote the 22
source labels to 22 independent animals. Current mScarlet status is absent from
the traced observation rows; GFP history uses the author's trace-sorted library
assignment, not an EGFP RNA threshold. The separate reporter-sorted Hopx controls
are excluded. No same-cell current-state reporter conclusion is available.

Age fields require care. The 8wk/14wk Slc4a11 harvest labels correspond to
6wk/12wk induction labels in the figure code. The generic `time` field for
Hopx_12wk3d inherits `6wk_3d` from its mixed sequencing lane; explicit driver/group
assignments identify the 12wk Hopx group. For IGO17543, the README says harvest
at 12 weeks whereas GEO says 14 weeks; the code's `12wk14d` label does not settle
that discrepancy. These fields are preserved and flagged rather than silently
harmonized. Independent temporal or driver comparisons remain held pending
pool/source verification and a clarified induction/harvest schedule. No cell-level
p-values or author cell-level tests were reused.

## J8 — synthesis and next discriminating evidence

The present evidence separates four measurements: induced-cell histone profiles,
epithelial perturbation RNA, measured PATS tracing, and HPCS trace-linked RNA-state
composition. The normal-differentiation domain reference adds genomic context.
They come from different experiments and cannot establish same-cell multimodal
linkage, regulatory causality or equivalence of named transitional populations.

The CDKN1A mark pattern is a useful locus-specific follow-up, but denominator
sensitivity and marker counterexamples prevent a generic transitional-epigenome
claim. IRE1α direction stability motivates follow-up while discovery instability
limits confirmation. The recovered HPCS descendants provide an explicit endpoint
table for later source verification, not independent fate inference today.

Next inputs are precise: PATS deposited-track scaling and H3 association; TIGIT
animal/pool memberships; CD44 matrix-to-GSM identities; and HPCS source/pool and
age reconciliation, with current reporter data if that endpoint is intended.
An author query may be drafted from these requirements; none has been sent.
Optional trajectory, spatial and raw-read extensions remain separate work.

## Evidence and reproducibility

- [Scientific verification](second_batch_verification.json): immutable inputs,
  fits, interval/control checks and domain reconstruction.
- [Delivery evidence audit](delivery_evidence_audit.json): official chromosome
  sizes, fresh identity evidence and original/corrected presentation checks.
- [Corrected render record](second_batch_verified_run.json): fresh PNG/SVG files
  and exactly unchanged summary-table hashes and numerical inputs.
- [Original render record](second_batch_summary_run.json) and
  [executed sources](execution_sources/second_batch_2026-09-25/manifest.json)
  remain at their original paths. The source archive is evidence, not a rerun route.
- [Gallery](../figures/README.md) contains the checked figures and assay-specific
  captions. All five new/corrected PNGs were visually inspected; SVG XML and text labels were
  checked. No R model was refit to correct the legend.
- [Continuation source archive](execution_sources/verified_batch_2026-09-25/manifest.json)
  preserves the exact executed verification/render/recovery/reconstruction scripts.
- [Validation record](verification_commands.json): Python compilation; repository
  tests (24 run, one expected scientific-runtime skip); all 11 A1 tests in the
  scientific environment without skips; 18 claim bindings; Nb1 output verification;
  repository validation; and `git diff --check` passed. Final repository
  validation counted 2,325 checks. These checks do not imply biological replication.

Assay payloads and retrieved documents remain ignored caches. Selected output
tables and reports are deliverable evidence; historical claim grades are unchanged.
