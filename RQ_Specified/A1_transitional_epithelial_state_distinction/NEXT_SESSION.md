# A1 — current handoff after the robustness batch

Updated 25 September 2026. The owner requested a concise Notion page and then
execution of the next plans. The Notion change and feasible robustness batch
are complete. Start with [ROBUSTNESS_REPORT.md](reports/ROBUSTNESS_REPORT.md).
PR #70 is merged; new branch: `codex/a1-robustness-analysis`.

Scripts 22–25 have finished; no process needs resuming. Source omissions and
annotation partition checks use the same 5,333 cells / 22 source aliases.
Native transcript-start sensitivity uses 46 positions over 23 loci. Old output
hashes remain valid and 1,104 baseline windows reproduce; 264 alternative windows
passed independent base-resolution checks. Both new figures were inspected.

The main new limit is substantive: CDKN1A acetylation reverses in CUT1 at the
alternative promoter. HPCS chase cannot be separated from unrestricted library
effects (rank 6 → 6 overall), and the 14wk composition is source-sensitive.
Keep the revised interpretation; do not promote the earlier locus lead.

Biological annotation recoding is held: the two pinned definition notebooks
exceed 20 MiB. Raw categories, including `other`, and label-invariant partition
metrics were retained. Do not infer a cluster-to-HPCS map from these same cells.

Next packages are conditional: verified units and an estimable design for TIGIT,
CD44 and HPCS inference; PATS track scaling/controls; independent state-specific
regulatory/fate data. The [remaining-analysis assessment](reports/REMAINING_ANALYSIS_OPTIONS.md)
separates completed robustness from these missing inputs. No author contact was
sent. Do not add the analysis menu or a critique to the concise Notion page.

## Historical second-batch handoff (superseded by the status above)

Updated 25 September 2026 after the owner resumed A1. **The feasible second
batch is now verified and reported.** Start with
[SECOND_BATCH_REPORT.md](reports/SECOND_BATCH_REPORT.md), [JOBS.md](JOBS.md)
and the [checked gallery](figures/README.md). No background process needs resuming.
Delivery branch: `codex/a1-regulatory-followup`; the owner has requested a PR
covering these accumulated results. See the
[remaining-analysis assessment](reports/REMAINING_ANALYSIS_OPTIONS.md) for
proposed next packages, which have not run.

J1–J3 are verified. Script 16 now requires `--run-id` and refuses overwrites;
the corrected presentation is `figures/second_batch_verified/`. Original figures,
records and numerical inputs are preserved. Scripts 17 and 20 record the
scientific and delivery audits. PATS pulse day 7 / harvest day 12 is documented.

J7 advanced beyond the original handoff: script 19 recovered only `/obs` from
the 8.05-GB GEO object (7.41 MB transferred), and script 21 reconstructed
5,333 traced cells across 22 source labels. Tables and figure are under
`hpcs_source_composition/`. This is descriptive source reproduction; current
mScarlet is unavailable and independent animals/pools are not verified.

Remaining work needs new evidence, not repeated fits:

- PATS: exact bedGraph scaling and H3 association.
- TIGIT: explicit animals/pool members and non-overlap for all four source blocks.
- CD44: R26/OG count columns → GSM → genotype → mouse → sort.
- HPCS: source-alias-to-mouse/pool crosswalk, IGO17543 induction/harvest age
  reconciliation, and current mScarlet labels if a reporter-persistence endpoint
  is intended. Generic `time` inherits the wrong lane label for Hopx_12wk3d;
  the explicit author driver/group assignment is preserved in the source manifest.

Do not rerun completed scripts to obtain the same evidence. No author contact
was sent and no claim grade changed. Publication of the delivery branch was
requested after the analysis continuation.
The report records the completed validation and remaining limitations.

## Historical pre-verification handoff (superseded by the status above)

25 September 2026. The owner requested **a plan for a separate session**, after
asking to divide the remaining A1 work into individual jobs. Active analysis
stopped at this handoff. All launched numerical processes have returned
successfully; no background job needs resuming. Final scientific verification,
interpretation and repository delivery are still pending.

## Start here

Use the **existing local checkout** `X:/GitHub/scRNA_seq`, branch
`codex/a1-regulatory-followup`, based on merged PR #69 (`bce5660`). This follow-up
is currently **uncommitted**. A fresh worktree or remote checkout will not contain
these outputs, untracked scripts, or ignored assay caches. Do not switch/reset
branches, clean the checkout, or blindly rerun a numerical entrypoint.

Read this document, [JOBS.md](JOBS.md), [PLAN.md](PLAN.md),
[the frozen second-batch scope](config/second_batch.md), and
[the locus contract](config/direct_mark_loci.json). Then inspect Git status and
the three completed run records below. Account for unrelated `.claude/` files
without staging or deleting them; the sandbox may not read the global ignore.

The scientific target remains a biological regulatory distinction among
RNA-similar transitional epithelial populations. DATP, PATS, ADI, ABI and HPCS
are study-defined populations, not interchangeable classes. Direct histone
evidence, chromatin accessibility, functional perturbation and measured lineage
each answer different parts of A1. Completing this batch does not settle a
universal state taxonomy or demonstrate chromatin causality.

## Completed computation to reuse

| Work | Actual output | Execution evidence |
|---|---|---|
| IRE1α robustness | Ten leave-one-mouse-out fits, S061 within-batch fit, S135 descriptive estimates; frozen 14,811-gene universe | [Run record](reports/ire1_stability_run.json), [tables](tables/ire1_stability_2026-09-25/fit_diagnostics.tsv) |
| Direct histone profiles | 24 CPM BigWigs; H3K27ac, H3K4me3, H3K27me3 and matched H3; 23 loci, two windows, 1,104 measurements | [Run record](reports/direct_mark_run.json), [manifest](tables/direct_marks_2026-09-25/sample_manifest.tsv) |
| Methylation-domain reference | One-donor D0/D4/D6 UMR/LMR/PMD overlaps; 414 intersections including two coordinate conventions | [Table](tables/direct_marks_2026-09-25/methylation_domain_overlap.tsv); same direct-mark run |
| Summary and figures | IRE1 stability, histone locus tracks, promoter heatmap and methylation-domain heatmap; PNG/SVG | [Summary run](reports/second_batch_summary_run.json), [summary tables](tables/second_batch_summary/ire1_focus_stability.tsv), `figures/second_batch/` |

The first-batch output hashes and original counts/sample-manifest MD5 values
were verified by the summary run. The primary model/results were not replaced.
Exact executed source bytes for scripts 11, 12, 15 and 16 are preserved in
[the execution archive](reports/execution_sources/second_batch_2026-09-25/manifest.json).

**Presentation issue already identified:** the first IRE1 figure's legend
overlaps two marker rows. Script 16 has a one-line legend-position correction
that has **not been rendered**. Its current source hash therefore differs from
the completed render record; the archived script matches that record exactly.
Preserve the initial render record/figures before producing a corrected render,
and verify that numerical inputs are unchanged. Three figures were visually
inspected; the methylation panel still needs inspection. Full repository checks
and direct interval/annotation tests have not run for this batch.

## Jobs, in execution order

### J1 — finish IRE1α sensitivity verification and interpretation

1. Verify the run's output hashes, ten distinct omissions, 9 mice per omission,
   balanced eligible arms, design rank/residual df, and unchanged primary inputs.
   Confirm the retained gene universe and pathway membership match the original.
2. Check summary joins against the complete effect tables. Keep all four
   original FDR hits labelled **post-selection**; do not add them to the frozen
   eight-marker family or promote a union of sensitivity discoveries.
3. Correct/rerender the overlapping legend with a fresh presentation record;
   inspect PNG and SVG readability. Do not refit R for a figure-layout change.
4. Write a package report with effect direction, omission ranges and batch
   limitations. These ranges are not confidence intervals.

Observed, pending final report: seven of eight frozen markers keep their primary
direction in all ten omissions; Ager does so in nine. Significant-gene counts
range from 1 to 264 across omissions versus 4 in the primary model. TGF-β points
downward throughout, but passes FDR 0.05 in only the omission of mouse 148
(`q=0.04490`; primary `q=0.12463`). S061 has 112 whole-family FDR hits; S135 has
only two mice per arm and remains descriptive. These observations support
directional follow-up, not a replacement primary analysis or confirmed pathway.

**Done when:** effect/pathway summaries and corrected figures are checked,
primary provenance remains intact, and the report states these limits.

### J2 — finish direct histone verification and biological interpretation

1. Test interval clipping, half-open coordinates, gaps, non-overlap and
   negative-strand TSS handling with hand-calculated cases. Independently compare
   a few BigWig interval means to base-resolution reader values.
2. Verify all 24 track headers share native CHM13 chromosome sizes; check the
   exact-symbol gene annotation, sample/mark/CUT pairing, zero handling and
   matched H3 ratios. No hg19/mm10 coordinates may enter this comparison.
3. Compare ±1-kb promoter results with the already computed ±5-kb sensitivity;
   distinguish raw mark effects from changes introduced by the H3 denominator.
   Report all 23 loci, held ratios and both preparations; do not select only
   concordant loci. The gene-boundary TSS is not necessarily the dominant isoform.
4. Write the report and captions; keep two-preparation/one-iPSC-line results
   descriptive. CPM and mark/H3 ratios are not absolute or spike-in-scaled
   occupancy, and no population p-values are licensed.

Observed promoter direction concordance: H3K27ac 36/45 eligible comparisons,
H3K4me3 36/42, H3K27me3 38/42 (46 possible per mark = 23 genes × 2 contrasts).
These are comparisons, **not independent samples**. For the predefined CDKN1A
locus, iATCs versus iAT2 show higher H3K27ac/H3 (+0.84, +1.27 log2) and lower
H3K27me3/H3 (−1.38, −0.89). This is a locus-specific descriptive lead, not proof
that every named transitional state has the same epigenetic programme.

**Done when:** numerical coordinate/control checks pass and biological leads,
counterexamples and denominator sensitivity are reported together.

### J3 — finish the normal-differentiation methylation reference

1. Verify domain intersections and inspect the fourth figure.
2. Preserve both coordinate conventions: author tables start at 1, suggesting
   R/GRanges output; the alternative BED interpretation changes an overlap by at
   most 0.0005 in the present table. State the assumption rather than claiming
   metadata verification that has not happened.
3. Explain UMR/LMR/PMD overlap as genomic context. Author UMR/LMR tables also
   contain segment methylation summaries, but this run did not estimate
   promoter-wide CpG methylation from them. No replicated DMR inference.

**Done when:** the reference panel is verified and described separately from
injury/tumour transitional-state evidence. This job has no sample-identity hold.

### J4 — finalize PATS tracing; isolate the histone hold

Resolved: the open-access published article (`PMC7461628`) explicitly specifies
a Krt19-CreER pulse **7 days after injury** and a **day-12 harvest**. Update
lineage records accordingly; this is a nominal five-day interval, not an
estimated transition rate. The ED4 source reconstruction remains 3 mice per
marker, separate marker denominators, with undefined zero-denominator controls.

The histone preparation is a separate experiment: CTGF-positive cells on day 12;
TP53 ChIP uses a distinct day-8 sort. Do not conflate them. The methods describe
H3-normalized peak calling, but that does not establish the scaling of every
deposited bedGraph. The generic MintChIP code is cached for inspection.

Four bedGraphs downloaded; CTGF-positive H3K4me3/H3K36me3 exceed the frozen
256-MiB per-file ceiling. Establish exact track normalization/control association
before any quantitative between-condition histone comparison. Do not increase
downloads simply to bypass the normalization hold. If a verified route requires
larger data, record measured size/runtime and freeze that extension first.

**Done when:** timing is documented and the histone comparison either has a
verified common quantification contract or an explicit, narrowly stated hold.

### J5 — resolve the TIGIT ATAC biological units

Inspect cached `hpcs_bioc.xml` (Marjanovic 2020, `PMC7745838`), GEO/ENA metadata
and any linked author sample tables. The article describes four mouse
replicates and a paired `~Mouse + Tigit_status` analysis. That is useful evidence,
but does not yet reconcile the deposited alias `106621_106642` with YY1181,
YY1916 and 106623. No explicit pool-member/non-overlap crosswalk was recovered.

Required: map each positive/negative count column to a verified animal or
disjoint pool, with source citation. If resolved, update the candidate manifest,
freeze ATAC-specific filtering/background and CNV/depth checks, then fit the
eligible paired contrast. Three verified independent units per arm/pair is the
existing minimum, not evidence of sufficient power. Do not infer identity from
title order, filename resemblance or a favourable PCA.

**If still blocked:** save the exact missing crosswalk and completed audit;
continue J6. An author request can be drafted, but do not send it without an
explicit instruction to contact the authors.

### J6 — resolve CD44 counts versus genotype and mouse

GSE273123 provides 16 libraries/8 deposited source aliases. The source workbook
and the published KatzenLab code were downloaded and inspected. Figure-3 source
data contain panels 3b/3d/3h, not the RNA matrix's column crosswalk; the code does
not resolve R26/OG identifiers to the WT1–4/mutant1–4 GSM labels. Do not propagate
the separate GSE273122 title/characteristic contradiction into this accession.

Required: exact matrix column → GSM → genotype → biological unit → CD44 sort
map. Audit remaining raw/processed headers or author supplements only where
they could provide that mapping. Then freeze within-genotype paired CD44
contrasts. A difference between two significance calls is not an interaction.

**If still blocked:** retain descriptive PCA and a precise crosswalk request;
continue J7. Do not fit guessed genotype labels.

### J7 — recover HPCS descendant data independently

Chan 2026's downloaded source workbook contains growth/ablation panels from
Figures 3, 4 and Extended Data 10, not Figure-2 descendant composition. The
article links [HPCS_LUAD](https://github.com/dbetel/HPCS_LUAD),
[Zenodo 17662770](https://doi.org/10.5281/zenodo.17662770) and GSE277777.
The repository tree is already cached as `cache/followup_sources/hpcs_code_tree.json`
but **has not yet been inspected**. Start there, downloading only relevant
notebooks or compact processed metadata, not the complete sequencing archive.

Recover GFP history, mScarlet current state, driver, 3-/14-day chase, tumour age,
source mouse/pool and denominator before analysis. Some conditions pool two
mice. Do not treat cells as biological replicates or growth as descendant
composition. If compact verified labels exist, freeze a source-level endpoint
reconstruction; otherwise record the specific missing metadata/data object.

**Done when:** a defensible measured-endpoint table is analysed, or its distinct
external-data hold is documented without holding J1–J6 hostage.

### J8 — synthesize and deliver the batch

Create `reports/SECOND_BATCH_REPORT.md`, update A1 README/PLAN/STUDY_MAP and the
question-specific gallery, then update the canonical A1 card in
`RESEARCH_QUESTIONS.md`, PROGRESS and AI_CONTEXT. Separate:

- regulatory profiles actually measured;
- perturbation effects and their stability;
- measured lineage evidence from different experiments;
- unresolved tests and the exact input needed for each.

Do not invent same-cell multimodal linkage or upgrade historical claim grades.
Run the repository-required checks below, plus scientific checks relevant to
the new code. Stage only A1-related work and intended documentation. Review
large regenerated tables before Git inclusion. The current session created no
commit or PR for this batch; decide delivery in the resumed session's scope.

## Environment and commands

```powershell
Set-Location 'X:/GitHub/scRNA_seq'
$a1Base = 'RQ_Specified/A1_transitional_epithelial_state_distinction'
$a1Python = 'C:/Users/dream/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$a1Rscript = 'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'
git status --short
git branch --show-current
```

Python 3.12.14. Scientific packages are supplied by
`analysis/scripts/run_with_environment.py --site-packages .venv-x64/Lib/site-packages`.
The isolated `cache/pybigtools/` contains pybigtools 0.3.0 (Windows wheel).
R uses the existing `.tools/R-library` with edgeR 4.10.5, limma 3.68.5 and statmod.
Assay payloads, annotations and source workbooks remain ignored under `cache/`;
selected counts remain under `processed/`. Do not stage these large inputs.

**Do not rerun completed numerical scripts 11 or 15.** They refuse output
overwrites. Script 16 is a renderer/summarizer but currently rewrites its output
folder and record: archive or add a fresh-output option before its next use.
The script archive is evidence, not an alternative execution entrypoint.

Repository checks after the actual edits:

```powershell
& $a1Python -m compileall -q -x '[/\\]\.tools[/\\]' analysis 'Research Article' RQ_Specified
& $a1Python -m unittest discover -s analysis/tests -q
& $a1Python analysis/scripts/claim_contract.py --check
& $a1Python 'Research Article/gate1_03_nabhan_2018/nb1/verify_outputs.py'
& $a1Python analysis/scripts/validate_repository.py
```

Run meaningful scientific-runtime tests separately; lightweight CI skips some
scientific integration checks. Record failures and resolve them without rerunning
unaffected expensive models. Source inventories retain URLs/hashes and download
failures. The UCSC reference TLS issue was resolved using Windows curl's trusted
certificate store; certificate validation was never disabled.

## Prompt to paste into the next session

> Continue A1 in the existing local checkout `X:/GitHub/scRNA_seq` on
> `codex/a1-regulatory-followup`. Read
> `RQ_Specified/A1_transitional_epithelial_state_distinction/NEXT_SESSION.md`
> and `JOBS.md`. Execute J1–J8 step by step, report each job as it finishes,
> and keep external-data holds separate. Reuse the completed analyses; preserve
> all original evidence and thresholds. Finish verification, reports and figures
> before proposing additional computation.
