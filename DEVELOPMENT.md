# Development and scientific ownership

This project was developed through an AI-assisted research workflow. I
directed the scientific design, dataset selection, analytical constraints,
validation standards, and biological interpretation; Claude Code assisted
extensively with implementation, execution, debugging, and documentation. AI
assistance is intentionally visible in the git history — the
`Co-Authored-By` trailers are disclosure, and nothing has been rewritten to
hide them.

This page is for human readers. The machine-oriented counterpart — repository
rules, environment constraints, and pitfalls a future AI session must not
violate — is [`AI_CONTEXT.md`](AI_CONTEXT.md).

## Who decided what

| Responsibility | Primary authority |
|---|---|
| Scientific question and dataset selection | Me |
| Experimental-design interpretation (cohorts, confounds, sort ratios) | Me |
| Analytical constraints and validation criteria | Me |
| Biological interpretation of results | Me |
| Implementation and refactoring | AI-assisted |
| Execution of computational runs | AI-assisted |
| Figures, tables, and report generation | AI-assisted |
| Review of intermediate and final outputs | Me |
| Retain / revise / reject decision on every result | Me |
| Literature extraction into study notes and reviewable configs (`Thesis/`) | AI-assisted; my review pending |
| Focused reproductions of the source paper's phase and myeloid claims (`phase_timecourse/`, `myeloid_focus/`) | AI-assisted, rules frozen before each run; my review pending |
| Repository framing, and what is displaced as established elsewhere | Me (instruction 2026-09-10); AI-assisted execution |

AI execution never meant automatic acceptance. Results were reviewed between
sessions, and several were sent back: one finding was refuted and rewritten
(decision 6 below), one annotation approach was demoted after it failed
validation (decision 5), and documentation that overstated what had run was
re-scoped rather than left to stand.

## The decisions that shaped the analysis

Abstract claims of supervision are cheap; these are the concrete decisions,
each verifiable in the repository's artefacts.

**1 · Raw data determines the workflow.** The governing brief forbade
assuming anything about the deposited data — format, species, metadata,
gene-space compatibility, reporter features, and QC thresholds all had to be
detected, not presumed. This constraint is why the pipeline found things a
template would have missed: the `SiteA`/`SiteB` lineage-reporter contigs
hiding in the mouse gene space, the human sample quantified against a
different annotation build (forcing an Ensembl-ID intersection), and
non-unique human gene symbols.

**2 · Mouse: no batch correction, because the design forbids it.** Every
mouse sample belongs to exactly one experimental group, so sample identity
and the influenza time course are the same variable — "correcting" on sample
would delete the biology under study. I had the question reframed to one that
is purely technical: do replicate animals *within* a group fail to mix?
Measured within-group replicate enrichment was 1.374 (1.0 = perfect mixing):
no correction. Reading the paper afterwards confirmed the authors integrated
nothing either.

**3 · Human: Harmony, against the pipeline's own default.** The three human
donors are healthy biological replicates, so donor separation *is* technical
— but with no condition metadata, the pipeline's automated rule declined to
recommend integration (`integration_recommended: false`). The deciding
evidence was biological: single proposed cell types were fragmenting into
donor-private clusters in the uncorrected embedding — one cell type is not
several cell types in several donors. Harmony was applied as an explicit,
logged override (`--integration harmony`), the unintegrated embedding was
retained for comparison, and the donor-driven-clustering check now reports
false. Together with decision 2 this is the point: the answer is not "batch
correction good" or "bad" — **the experimental design decides.**

**4 · Author labels held out.** The deposited annotations were excluded from
every clustering and trajectory step and used only afterwards, as an answer
key. That is what makes the median cluster purity of 0.947 and the pseudotime
ordering of the labels *validation* rather than circular confirmation.

**5 · Contradicted annotations retained, not overwritten.** The blind
marker-panel annotator disagrees with the deposited labels on 3 of 29 mouse
clusters; those proposals are flagged `[CONTRADICTED]` in the tables and on
the UMAP rather than silently corrected. The instructive failure is kept on
display: cluster 0's panel score said "transitional epithelium", but it is
89.2% CAP1 endothelium expressing an interferon program — a *state* that
fooled a *type* classifier.

**6 · A convenient finding was refuted and the refutation kept.** An early
audit suggested Scrublet was over-removing the human AT0 population at up to
2× background — a finding that would have flattered the project's critical
posture. I had the test itself examined: a co-expression gate cannot audit a
co-expression detector. A stricter gate (SFTPC⁺ SCGB3A2⁺ EPCAM⁺,
lineage-negative) reversed the conclusion — AT0 flagged at 3.9% vs a 6.3%
baseline, below background — and the full sequence, including the wrong first
pass, is documented in [`docs/DOUBLETS_AND_SCRUBLET.md`](docs/DOUBLETS_AND_SCRUBLET.md).

**7 · The trajectory required a commissioned re-analysis.** The whole-atlas
embedding does not answer the regeneration question — at atlas resolution the
alveolar states are clusters, not an ordered process. I commissioned a
focused analysis restricted to the 25-sample annotated cohort's alveolar
epithelium (5,694 cells), with PAGA and diffusion pseudotime rooted in AT2,
which recovered the AT2 → Krt8⁺ transitional → AT1 ordering with the labels
held out.

**8 · The capillary injury state required compartment-specific reclustering.**
The source study reports that the persistent injury state (iCAP) is not
resolvable at top-level clustering, so it was never expected to appear as an
atlas cluster. The focused capillary analysis (43,359 cells, reclustered
within-compartment) recovered it, including its defining behaviour: it
emerges after infection and does not resolve by one year.

**9 · The CAP2 tracing result is uninformative, not negative.** The
CAP2-specific Cre lines label only 2–8% of endothelium, so their near-zero
trace rates in the injury state cannot be distinguished from insufficient
labelling. The report says "uninformative" where "negative" would have been
the stronger — and unsupportable — claim. The Kit line, which labels CAP1,
traces the injury state at 33–53% per animal, and that is the claim actually
made.

**10 · The "not done" boundary is explicit.** No ambient-RNA correction (the
required raw droplet matrices are absent for the mouse series); no formal
trajectory-DE model; no whole-lung composition claims from MACS-enriched
material (the cell-type ratio is a sort ratio); and the tool-reference pages
in `docs/` describe the published method, not what ran — the generated
[`docs/PIPELINE_AS_RUN.md`](docs/PIPELINE_AS_RUN.md) is the authoritative
used/not-used record.

**11 · Reference criteria are adopted by pre-registration, not by copying.**
(2026-09-09, under review.) The HLCA paper (Sikkema et al. 2023) was read as
the roadmap's "reference framework". Its integration benchmark, entropy
thresholds, marker filters, sample-count rule and label-transfer uncertainty
cutoff were extracted into a reviewable JSON file, and the AI session then
proposed, criterion by criterion, what this pipeline should adopt, adapt or
decline (`Thesis/gate1_04_sikkema_2023_hlca/PIPELINE_FRAMING.md`). Two
proposals are already fixed by earlier decisions: supervised scANVI on the
deposited labels is declined because it would make decision 4 circular, and
the HLCA donor-entropy threshold is not copied because it encodes 107 donors.
A first trial applied the label-entropy and donor-entropy rules to the
tracked cluster tables with the thresholds frozen from the paper before the
tables were opened. My retain/reject decisions on the extracted material,
and any rejected AI output, will be recorded here.

**12 · Portfolio framing follows the target labs, not the source paper's
injury narrative.** (Owner instruction, 2026-09-09.) For the UC Berkeley
outreach portfolio the interferon and influenza context is not the point;
results are framed by cell state, repair and niche biology, macrophage and
monocyte states, annotation robustness and curation hygiene, matched to the
target labs. The data do not change; the write-ups do. Trials S3 to S5 and
their plan follow this rule.

**13 · Rule revisions are disclosed, not silently applied.** Trial S3's
first run used a two-marker minimum and a compartment set that the source
sheet does not support; eight identities were lost and a spurious
compartment appeared. The rules were corrected, the first-run outcome was
kept in the run record, and a hierarchical assignment was added as a
post hoc sensitivity rather than swapped in as the primary. The two schemes
disagree on AT0, and that disagreement is the result, not a nuisance to be
resolved by picking the scheme that flatters the earlier AT0 candidate.

**14 · A headline that the reference route does not support is flagged,
not softened.** (2026-09-09.) With my authorisation the AI session
installed PyTorch and scvi-tools into the emulated interpreter and mapped
the human series onto the HLCA core by scArches surgery (trial S2). The
mapping reproduces the HLCA authors' own transfer of the same cells at
99.2% (level 3), so the environment is sound, and it says the "AT0
candidate analogue" subcluster is mostly AT2 or uncertain, with AT0 a
minority of the series. That contradicts the wording of the human headline
in `FINDINGS.md` and the portfolio PDF. The contradiction is recorded in
`PROGRESS.md` with a proposed re-wording; the decision to retain, re-word
or reject is mine and is pending. The `scarches` package could not be
imported with the pinned anndata and was removed; the surgery uses the
scvi-tools implementation and the deviation is disclosed in the trial.

**15 · The paper's phase structure is reproduced from the trace, not from
the classifier.** (2026-09-10, under review.) The source paper's central
descriptive claim, that proliferation after injury runs in phases (immune,
then epithelium and mesenchyme, then endothelium), was reproduced from the
tracked metadata alone, with the expected peak windows frozen in the run
record before the table was opened. The Ki67-trace peak falls in the
paper's window for four of five lineages; lymphoid cells peak one harvest
later. The deposited cell-cycle call was carried as a cross-check and
disagrees for four of five lineages because it calls most lymphocytes
cycling; it is reported, not used. The myeloid compartment was then taken
from the blind atlas clustering (clusters 5, 17, 24), re-embedded with the
labels held out, and graded afterwards, at the resolution fixed from trial
S5; its per-animal composition reproduces the paper's Figure 3 (aMAC loss
and iMON expansion at 6 dpi, reconstitution by 19 to 42 dpi). Both analyses
use the deposited labels descriptively and say so, and neither computes a
P value, because the active-repair days carry two animals each. Retain,
re-word or reject is mine and pending.

**16 · Batch is tested on a key that crosses time, never on the animal, and
a rule that selected the wrong object is disclosed, not swapped.**
(2026-09-10, under review.) The reviewer's question about the myeloid
result is whether the 6 dpi inflammatory-monocyte state is a one-day batch
island. Correcting on sample cannot answer it (one animal per sample and
day). The paper's Table S3 shows that on every active-repair day and at 42
dpi the two replicate animals came from different infection rounds, so
round is a technical key that crosses time; it also carries Ki67-Cre
dosage, which makes correcting on it the conservative direction. The
compartment was embedded with and without Harmony on round, on the days
that carry both rounds. The pre-registered survival rule defined "the iMON
state" as the subcluster with the most iMON-labelled cells, and the first
run showed that this picks the 11 to 19 dpi monocyte state in both
embeddings (3.9% and 22.0% of its cells from 6 dpi), not the 6 dpi state
the question is about. As in decision 13, the first-run outcome stays in
the run record, the definition anchored on the 6 dpi cells is added as a
labelled post hoc reading, and both are reported. The same session read the
Ki67 trace by tamoxifen window for the rebuilt alveolar macrophage pool and
let two of its three pre-registered checks fail closed on a 30-cell floor
rather than lower the floor after seeing the tables. Retain or reject is
mine and pending.

**17 · The repository is an analysis log, and established material is
displaced rather than deleted.** (Owner instruction, 2026-09-10.) Two
kinds of content had accumulated that are not part of the ongoing analysis:
portfolio-curation material (a thesis-aware PDF and its generator) and the
Krt8-high transitional work (the alveolar trajectory, its time course, and
the human KRT8 reference-aligned panels), which is established elsewhere in
my G-SURF submission. Both were moved out of the main narrative into
`archive/` with a note on what moved, when and why; their artefacts and
scripts stay under `analysis/` because the validator checks their numbers
and the scripts regenerate them. The README was rewritten to open with the
claims table, to state the working question rather than a recovery
exercise, and to name H1N1 once as the injury model of one series, because
the roadmap ahead is not an influenza project. Claims were softened where
they overreached ("never resolves" became "has not resolved by 366 dpi").
The checks were renamed from portfolio to repository checks. Nothing was
deleted and no number changed.

## How outputs were reviewed

Every run writes its decisions to machine logs (`decisions.json`,
`analysis_log.txt`), and the per-dataset reports and pipeline record are
*generated* from those artefacts — numbers in the documentation cannot drift
from what was computed, and a claim I couldn't trace to an artefact was
treated as unverified and removed. Work was reviewed between sessions against
[`PROGRESS.md`](PROGRESS.md), landed through pull requests, and known issues
were carried forward in writing rather than dropped.

## Where the honest failures live

- The refuted Scrublet/AT0 over-removal claim, with both rounds of the test:
  [`docs/DOUBLETS_AND_SCRUBLET.md`](docs/DOUBLETS_AND_SCRUBLET.md)
- The three contradicted cluster annotations, flagged in
  `analysis/GSE262927/tables/cluster_annotation_proposals.csv`
- Thresholds that never bound, doublet calls that are a ranking rather than a
  detection, and every other caveat: [`FINDINGS.md § 6`](FINDINGS.md#5--negative-results-and-self-audits)
  and the per-dataset reports
