# Development and scientific ownership

This project was developed through an AI-assisted research workflow. I
directed the scientific design, dataset selection, analytical constraints,
validation standards, and biological interpretation; Claude Code assisted
extensively with implementation, execution, debugging, and documentation. AI
assistance is intentionally visible in the git history, the
`Co-Authored-By` trailers are disclosure, and nothing has been rewritten to
hide them.

This page is for human readers. The machine-oriented counterpart, repository
rules, environment constraints, and pitfalls a future AI session must not
violate, is [`AI_CONTEXT.md`](AI_CONTEXT.md).

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
| Reading the Cardoso 2026 deposit in gates, and stopping Gate 1 when it failed (`Thesis/gate2_05_cardoso_2026/`) | Me (instruction and gate design 2026-09-12); AI-assisted execution, rules frozen before each run; my review pending |
| Relaxing the displaced-material rule for the reference-aligned primary-marker panels (PR #12) | Me (decision 2026-09-13, against the agent's recommendation to hold them back); the agent had preserved them on a branch rather than reverting them, and flagged both the displaced-material rule and the non-validated palette; the palette is still unaddressed |
| Pre-registering the post hoc Fst and Runx2 lead and running it (trial C9) | Me (instruction 2026-09-13); AI-assisted execution. The agent named the untestable half before running, reported that the effect size its own trial computed is not usable, and disclosed that T1 and T5 were not composable; my review pending on rows C65 to C72 |
| Entering roadmap paper 2 (Choi 2020) on 2026-09-13: study note, extract, trial D0 | **Agent-proposed, approved by a one-word "proceed" on an agent-written list, not a specific instruction of mine** (corrected 2026-09-15). Later that day I said to leave the paper until I had read it. On 2026-09-15 I REJECTED the AI-written study note (withdrawn; kept in git history at PR #19), relocated trial D0 and its extract under the Cardoso folder as an extension, and reserved the paper-2 folder for after my own reading. Rows C58 to C64 stand, review pending |
| Re-entering roadmap paper 2 (Choi 2020) on 2026-09-15: study note, extract, trials D1 to D7 and the corrected passes D2b and D5b | Me (instruction 2026-09-15, after reading the paper: list the possible analyses, structure the pipeline, reproduce the paper's result first, then attack its claims, proceed until nothing remains; the backbone check trial by trial; the pull request and merge); AI-assisted execution with every rule frozen before its object was opened. The agent disclosed the Sftpc-clause defect in its own rule twice and left the thresholds where they were; my review pending on rows C85 to C104 |
| Revising the reading order (Gate 2 branches 2C, 2N, 2W for the placement targets; papers 12 to 16 and methods references M1 to M7 added) | Me (instruction 2026-09-15); AI-assisted execution with every identifier verified against PubMed |
| Promoting the generic deposit readers out of the Cardoso helper module into the shared one | Agent proposal, my approval implied by the instruction to proceed; a gate1 folder importing from a gate2 folder was the wrong dependency direction. Verified by importing all 19 Cardoso trial modules and re-running two trials to identical results; cardoso_utils re-exports so nothing written against it changed |
| Choosing the three follow-up questions off the agent's ranked list, and accepting three negative answers plus two disclosed rule defects (trials E6, C7, C8) | Me (selection 2026-09-13); AI-assisted execution with readings frozen before the data were opened; the agent reported the margin that refuses its own C7 answer and the zero-inflation that voids its own C8 mixture test; my review pending on rows C49 to C57 |
| Extending the Hbegf lead into four public datasets once the deposit was exhausted, and accepting a frozen rule's refutation of my own agent's best lead (trials E1 to E4) | Me (instruction 2026-09-13, including the condition that the extensions run only if the first did not refute); AI-assisted execution, readings frozen before the matrices were opened; my review pending on rows C37 to C48 and on the weakening of C29 |
| Opening two branches of paper 2 on other laboratories' multiome deposits (2026-09-20): the transitional state in chromatin, and the paper's own closing Axin2 and Il1r1 question | Me (instruction 2026-09-20); AI-assisted execution. The agent refused five times before producing a statistic and retracted the one reading it produced; I corrected Route C as out of focus and directed Route A to run. Rows C116 to C150, review pending |
| Auditing rows C116 to C150 with eight adversaries and independent verifiers, and accepting all 45 confirmed findings | Me (instruction 2026-09-20: "launch a multi-agent attack workflow to testify all the claims made"); AI-assisted execution. The agent's own register prose was the thing found wrong; corrections read from tables, trial M4 added, rows C151 to C154; review pending |
| Rewriting the root README around a ledger drawn from the register, listing all fifteen opened deposits, and writing RESEARCH_QUESTIONS.md as the question-first entry point for outreach | Me (instruction 2026-09-20, including the choice to drop the Leiden UMAP and to organise by question rather than by trial); AI-assisted execution. A first draft cited seven wrong register rows, caught by printing each cited row; a second adversarial workflow checked the landed text. My retain or reject on the outreach framing pending |

AI execution never meant automatic acceptance. Results were reviewed between
sessions, and several were sent back: one finding was refuted and rewritten
(decision 6 below), one annotation approach was demoted after it failed
validation (decision 5), and documentation that overstated what had run was
re-scoped rather than left to stand.

## The decisions that shaped the analysis

Abstract claims of supervision are cheap; these are the concrete decisions,
each verifiable in the repository's artefacts.

**1 · Raw data determines the workflow.** The governing brief forbade
assuming anything about the deposited data, format, species, metadata,
gene-space compatibility, reporter features, and QC thresholds all had to be
detected, not presumed. This constraint is why the pipeline found things a
template would have missed: the `SiteA`/`SiteB` lineage-reporter contigs
hiding in the mouse gene space, the human sample quantified against a
different annotation build (forcing an Ensembl-ID intersection), and
non-unique human gene symbols.

**2 · Mouse: no batch correction, because the design forbids it.** Every
mouse sample belongs to exactly one experimental group, so sample identity
and the influenza time course are the same variable, "correcting" on sample
would delete the biology under study. I had the question reframed to one that
is purely technical: do replicate animals *within* a group fail to mix?
Measured within-group replicate enrichment was 1.374 (1.0 = perfect mixing):
no correction. Reading the paper afterwards confirmed the authors integrated
nothing either.

**3 · Human: Harmony, against the pipeline's own default.** The three human
donors are healthy biological replicates, so donor separation *is* technical,
but with no condition metadata, the pipeline's automated rule declined to
recommend integration (`integration_recommended: false`). The deciding
evidence was biological: single proposed cell types were fragmenting into
donor-private clusters in the uncorrected embedding, one cell type is not
several cell types in several donors. Harmony was applied as an explicit,
logged override (`--integration harmony`), the unintegrated embedding was
retained for comparison, and the donor-driven-clustering check now reports
false. Together with decision 2 this is the point: the answer is not "batch
correction good" or "bad", **the experimental design decides.**

**4 · Author labels held out.** The deposited annotations were excluded from
every clustering and trajectory step and used only afterwards, as an answer
key. That is what makes the median cluster purity of 0.947 and the pseudotime
ordering of the labels *validation* rather than circular confirmation.

**5 · Contradicted annotations retained, not overwritten.** The blind
marker-panel annotator disagrees with the deposited labels on 3 of 29 mouse
clusters; those proposals are flagged `[CONTRADICTED]` in the tables and on
the UMAP rather than silently corrected. The instructive failure is kept on
display: cluster 0's panel score said "transitional epithelium", but it is
89.2% CAP1 endothelium expressing an interferon program, a *state* that
fooled a *type* classifier.

**6 · A convenient finding was refuted and the refutation kept.** An early
audit suggested Scrublet was over-removing the human AT0 population at up to
2× background, a finding that would have flattered the project's critical
posture. I had the test itself examined: a co-expression gate cannot audit a
co-expression detector. A stricter gate (SFTPC⁺ SCGB3A2⁺ EPCAM⁺,
lineage-negative) reversed the conclusion, AT0 flagged at 3.9% vs a 6.3%
baseline, below background, and the full sequence, including the wrong first
pass, is documented in [`docs/DOUBLETS_AND_SCRUBLET.md`](docs/DOUBLETS_AND_SCRUBLET.md).

**7 · The trajectory required a commissioned re-analysis.** The whole-atlas
embedding does not answer the regeneration question, at atlas resolution the
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
the stronger, and unsupportable, claim. The Kit line, which labels CAP1,
traces the injury state at 33–53% per animal, and that is the claim actually
made.

**10 · The "not done" boundary is explicit.** No ambient-RNA correction (the
required raw droplet matrices are absent for the mouse series); no formal
trajectory-DE model; no whole-lung composition claims from MACS-enriched
material (the cell-type ratio is a sort ratio); and the tool-reference pages
in `docs/` describe the published method, not what ran, the generated
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
scripts stay in place (then under `analysis/`; the artefacts moved with their
series to `Thesis/` on 2026-09-21, decision 25) because the validator checks
their numbers and the scripts regenerate them. The README was rewritten to open with the
claims table, to state the working question rather than a recovery
exercise, and to name H1N1 once as the injury model of one series, because
the roadmap ahead is not an influenza project. Claims were softened where
they overreached ("never resolves" became "has not resolved by 366 dpi").
The checks were renamed from portfolio to repository checks. Nothing was
deleted and no number changed.

**18 · A gate that fails is a result, and a rule that fails is disclosed
three times rather than repaired once.** (2026-09-12, under review.) I
instructed that the Cardoso 2026 paper be entered out of roadmap order and
analysed in gates, with the standing condition that an unexpected result be
reported and acted on rather than finished around. Gate 0 read the deposit
before anything was fitted and found the constraint that governs everything
after it: every deposited mouse library pools three mice and each genotype
contributes one library per sort, so no genotype contrast in that deposit
carries within-group replication. The same trial closed one of my own
proposed branches: the mesenchymal and immune libraries are a single time
point, so the paper's fibroblast-before-macrophage ordering cannot be tested
transcriptomically, and that was reported as the answer rather than
approximated with something weaker.

Gate 1 then returned "not recovered" against its own pre-registered rule, and
the rule turned out to be at fault in four separable ways: it selected
fibroblasts by a confidence floor that excluded the candidates; two genes of
the paper's marker set are mural markers, so the score peaked on smooth
muscle; a third is a myeloid transcription factor, so it also flagged a sort
contaminant; and the population it was hunting missed the Tnc floor by one
thousandth. None of these was fixed in place. The first-run outcome stands in
the record, the corrected passes sit beside it as C1b to C1d, and the
threshold was not moved after the fact, because a rule moved to fit the
result it just failed is not a rule. The same defect recurred in Gate 2b and
was handled the same way, in C2b.

Two findings came out of that stopping rather than out of the original plan.
The mesenchymal sort carries 6.5% off-target cells, among them 184 mutant
epithelial cells that are 88% Areg-positive, almost entirely in the tumour
arm; the paper is unaffected because it took its epithelium from a separate
series, but a reanalysis computing signalling inside that one library would
have been reading a contaminant as the source. And on the Areg-deletion arm,
four of the six genes in the paper's fibrotic set fall with the ligand while
Pdgfrb and Runx1 do not, which suggests the state has separable parts. Both
are mine to retain or reject, and both are recorded as Exploratory until I do.

**19 · Leave the deposit to get a testable unit, and let a rule frozen in
advance refute your own best lead.** (2026-09-13, under review.) Decision 18
recorded that the Cardoso deposit cannot test anything, and the agent's leads
from it were accordingly all Descriptive only or Exploratory. I directed that
the most interesting of them, the Hbegf lead, be extended into public data
rather than written up from one deposit, and that the extensions run only if
the first one did not refute the working picture. Four extensions were
pre-registered. Three things about how they turned out are worth recording.

First, moving datasets is what made a test possible. GSE131907 has eleven
donors with paired tumour and normal lung, so the unit becomes the donor
instead of the library, and the first admissible test in this whole section
returned AREG detection higher in epithelium than in myeloid cells within
donor, p = 0.0020. Every earlier row in the section is descriptive because of
a design choice in the deposit, not because of anything the analysis did, and
this is the demonstration.

Second, and this is the part I want kept prominently, a rule frozen before the
data were opened refuted the agent's own best lead. The three-tier fibrotic
response of C29 had been the most promising thing the Areg arm produced, and
its proposed reading was a second, tumour-specific signal driving the retained
Runx1 and Pdgfrb tier. Trial E4 asked what bleomycin alone does to the same
genes in sorted mesenchyme with no oncogene present, under a rule requiring
both injured animals to exceed both controls and with the interpretation of
each outcome written down in advance. Both genes cleared it. The
pre-specified reading therefore applied with no discretion left: the retained
tier is what an activated lung fibroblast does after any injury. The numbers of
C29 stand and its interpretation is gone, which is the correct outcome and the
reason for freezing readings rather than only thresholds.

Third, the agent reported two results that cut against its own earlier
statements without being asked to look for them. The compartment gate that E1
called "neutrophil" was 83 per cent deposited myeloid cells in a deposit with
no neutrophils at all, which changes which pre-named outcome that trial hit;
and the AREG source that E1's test established as epithelial-over-myeloid is
qualified at subtype resolution, where dendritic cells sit above both tumour
epithelial states. Both are in the trial log next to the original wording
rather than replacing it.

What the extensions did not deliver is also on the record. E2's only testable
comparison failed to detect the difference it was built on, seven donors and
p = 0.297, and the Zhao prediction could not be tested in either fibrosis
cohort because no comparison cleared the five-donor floor and the two cohorts
disagree in direction. Those are reported as Not established rather than as
trends. The one thing that replicated across three datasets, a
dendritic-cell and monocyte ligand source, is explicitly recorded as
established immunology rather than as a finding of this repository, and is
kept because it constrains a reading, not because it is new.

Mine to retain or reject: the weakening of C29, rows C37 to C48, and whether
the myeloid ligand source is worth a paragraph in the eventual writeup.

**20 · Three chosen questions, three negative answers, and two rules that
failed instead of the biology.** (2026-09-13, under review.) I picked three
jobs off the agent's ranked list: a donor-level test of the paper's axis, the
identity of the mesenchymal-sort contaminant, and whether the Areg-independent
tier is a population or a gradient. None returned a positive result, which is
the ordinary outcome of asking precise questions, and two of them failed
because the pre-registered rule was badly chosen rather than because the data
were silent. That distinction is the reason this entry exists.

**What the agent got right about its own rules.** In trial C7 the rule said to
name the best-correlating reference state and to report the margin. It did
both. The margin was 0.0056 while the same measure separated epithelium from
fibroblasts by 0.415, so the measure has compartment resolution and no state
resolution, and the named winner is meaningless. The agent reported the winner
as the rule demanded, then refused to read it, and left the threshold alone. In
trial C8 the mixture test preferred two components, which the frozen reading
calls a subpopulation, and the agent showed why that is arithmetic: 31.5 per
cent of the cells detect neither of the two genes in the score, so a spike at
zero guarantees the second component. It noted that trial C1 had used the same
test on a centred score where the problem does not arise, so the defect is the
reuse and not the original. Neither rule was edited after the fact.

**The one thing that was edited, and why that was correct.** Trial E6's depth
rule said any pair whose two members both correlate with sequencing depth is
not read. The first implementation applied it only to the primary tests, so a
significant control pair came through unmarked. The agent changed the code to
match the rule text, not the rule text to match the result, and re-ran from
cached values. I am satisfied that is a bug fix rather than a threshold move.

**The result I want kept visible.** E6's only significant correlation was a
control: epithelial TGFA against fibroblast activation, p = 0.045. Both of its
variables track sequencing depth, so the rule threw it out. Had that number
landed on AREG instead of TGFA it would have read as confirmation of the
paper's axis, and nothing but the pre-registered control would have caught it.
The primary test itself is a null, rho 0.348 at p = 0.112 over 22 donors, and
the agent reported the smallest effect the test could have seen, about rho
0.43, rather than implying the axis is absent. It also said plainly that
pooling all epithelium dilutes the state the paper's claim is about, so the
null does not contradict the paper.

**What the three jobs bought.** List A of the next-step document is now empty:
every question answerable from data on disk has been answered, five of the six
negatively or with a correction. The surviving lead is post hoc and labelled as
such: the co-organisation in the fibrotic set sits on the tier that falls, Fst
with Runx2 at a ratio of 1.735 in the control arm and gone after deletion,
while the retained genes are independently distributed in both arms. It rests
on 24 double-positive cells in a single library, and the agent said so in the
same sentence as the finding.

Mine to retain or reject: rows C49 to C57, and whether the post hoc Fst and
Runx2 lead is worth a pre-registered trial of its own.

**21 · Paper 2 was entered on a "proceed", not on a direction, and its note
is withdrawn.** (2026-09-13; corrected 2026-09-15.) The first version of this
entry said I directed a return to roadmap order starting with Choi 2020. The
session transcript says otherwise, and the record has to match it. After the
E6, C7 and C8 jobs the agent listed what remained and put "returning to
roadmap order with Choi 2020" at the end of that list as its own
recommendation. I asked whether anything more was to go, the agent repeated
the list, and I answered "proceed". The agent read that as approval of every
item and entered paper 2, study note, extract and trial D0, all before I had
read the paper. Later the same day I said I would open a session on Choi 2020
after actually reading it and that it should be left; the agent honoured that
from then on, and D1 never ran.

On 2026-09-15 I settled it. The AI-written study note is REJECTED: a note on a
paper I have not read is not mine to retain, whatever its quality, and the
folder for paper 2 will be added when I have read the paper and decided what
analysis it deserves. The note stays in git history (PR #19) and is not on
main. Trial D0, its artefacts and the marker extract are kept and relocated
under the Cardoso folder as an extension, because the deposit facts it
established (one library per condition, six raw whitelists, no counted
reporter, coverage-only ATAC) bear directly on the DATP state the Cardoso rows
lean on. The trial keeps its identifier so its run record and rows C58 to C64
stay true; nothing was re-run.

What the trial found stands as recorded: the same replication ceiling as the
Cardoso deposit, in a second key paper in a row, which is worth keeping as a
pattern rather than a surprise. These are deposits from labs whose conclusions
rest on genetics and imaging, where the transcriptome is a map rather than the
evidence, and the deposit reflects that. The reading note also stands: the
agent found that the PMC web rendering strips italicised gene symbols and
switched to the Europe PMC XML rather than filling marker sets from memory. And
the reporting defect disclosed in the first run, a single boolean that read as
"none are raw" when six of eight were, is unchanged in the record.

The rule this adds: an agent may run a deposit reality check when I say so, but
the study note for a roadmap paper is written or directed by me after I have
read the paper. "Proceed" against an agent-written list is approval of the
list, and the agent should say which items it is treating as approved before it
starts the ones that commit a reading on my behalf.

Mine to retain or reject: rows C58 to C64 only. Gate 1 on this paper is not on
the table until my own folder exists.

**22 · A frozen rule has to say what would make its own answer unreadable.**
(2026-09-13, under review.) I asked for the post hoc Fst and Runx2 lead from
trial C8 to be given a pre-registration of its own. It was, and the trial did
three things I want on the record.

It named what it could not do before it ran. Whether Areg deletion depletes
that population is not testable with anything on disk, one library per genotype
and no other Areg-flox fibroblast dataset, and the trial said so in its own
docstring rather than producing a number that looked like an answer.

It separated existence from magnitude, and only one of them survived. Above
chance co-occurrence replicates in both bleomycin animals against a null that
holds sequencing depth fixed by construction, p = 0.005 and 0.010, in a dataset
with no oncogene. The size of the effect does not: about thirty double-positive
cells per library cannot support a stable ratio, and the agent said that rather
than quoting the ratios it had computed.

And for the third time in this folder, writing the numbers out exposed a defect
in a rule rather than in a result. T1's criterion rested on a ratio that T5
forbade reading, so the two rules were not composable; the implementation took
the conservative branch, its verdict stands, and the agent flagged that the
pre-registered reading attached to that verdict is not what the data did. The
same shape appeared in C7, where the profile test named a winner without
requiring a margin, and in C8, where a mixture test was applied to a
zero-inflated score. The lesson I am taking from the three together is the one
the agent stated: freezing a threshold is not enough, and a rule also has to
say in advance what would make its own answer unreadable.

The magnitude criterion in T4 is a fourth instance of the same thing, caught
the same way. Nine of twenty-seven markers cleared "at least half the reference
magnitude", and the median marker sits at 0.535 against that 0.5 line while the
comparison libraries carry half the genes per cell. The criterion was measuring
depth. What does survive is stronger and simpler: twenty-six of twenty-seven
markers point the same way in both animals.

Mine to retain or reject: rows C65 to C72, and whether the matrix and mechanics
programme those cells carry (Piezo2, Ltbp2, P4ha3, Sdc1, Prrx2) is worth a
pre-registration of its own, given that it points away from the Areg axis
rather than into it.

**23 · Paper 2 entered at my direction after reading; the same rule defect
twice, and a state that does not separate.** (2026-09-15, under review.)
Having read Choi 2020, I asked for the possible gene analyses to be listed, a
separate folder with a structured pipeline, the paper's results reproduced
first and its claims attacked after, until nothing remained; and, part way
through, for the core logic backbone to be fetched from Notion and the
repository and checked against the plan trial by trial. That backbone check
is a table in the plan. Trials D1 to D7 ran the same day, every rule frozen in
a run record before its object was opened, and I said to proceed to the pull
request and merge when the write-up was done.

Two things belong here. First, the annotation rule I approved froze two
clauses on Sftpc detection below 0.5 (one for contaminants, one for AT1), and
in a lineage-sorted AT2 library every cluster detects Sftpc in every cell, so
the clauses could never fire. The first pass called a 337-cell AT1 cluster
hAT2 and let a ciliated cluster wear the primed-AT2 label; the organoid pass
kept a 586-cell stromal cluster, the size of the one the paper removed, and
missed a 481-cell AT1 cluster, which made the paper's Figure 7 reading look
not computable when it was. Both outcomes stand in their records, no
threshold moved, and corrected passes D2b and D5b sit beside them with one
added condition: the Sftpc clauses apply only where Sftpc separates clusters.
This is the shape of decision 22 again, the fifth and sixth instances and the
first outside the Cardoso folder: a threshold frozen on a variable that does
not vary in the deposit. The rule I take from it: before a clause is frozen,
check that the variable it tests spans the threshold somewhere in the data
the rule will see. Trial D0's marker table had Sftpc present in every library
and nobody read it that way.

Second, what reproduces and what does not. Four of the paper's five states,
the DATP time course, the hAT2 to DATP to AT1 ordering inside one library,
DATP's programmes in vivo and IL-1beta's shift of the organoid epithelium all
come back from the deposit as descriptions, and the organoid cell counts land
within 8 percent of the paper's. The primed AT2 state does not: at three
resolutions, in a sub-clustering of DATP and at the cell level, no group of
cells loses Etv5, Abca3 and Cebpa while gaining the inflammatory genes, and
the organoid cluster the paper calls 77 percent primed carries the DATP
markers with its identity genes intact. Whether that is a graded state the
paper's cluster averages made discrete, or an identity ratio my rule set too
strictly, cannot be settled without a new pre-registration; the folder
proposes one (E6) and does not run it.

Mine to retain or reject: rows C85 to C104; whether to pre-register E6;
whether to download the dissociation list (M8) and run attack A3.

**24 · The register was audited by adversaries against its own artefacts,
and the entry point for outreach was written question-first.** (2026-09-20,
under review.) I asked for the two branches of paper 2 to be run, then for a
multi-agent attack on every claim they made, then for the repository to be
scanned and cleaned with the README updated first and a new document of core
questions and remarkable phenotypes that could stand behind a cold email to a
principal investigator, rather than the register's trial-by-trial alignment.

The audit is the decision that matters. Eight adversaries and forty-seven
verifiers found that the branches' trials were mostly honest and the register
written from them was not: four numbers had been quoted that no run produced,
including the Sox9 corroboration of the branch's one validated result, which
was never in the frozen panel and existed only as a typed string; three
magnitudes were wrong against their own tables; one row named the wrong
laboratory; one was refuted by the branch's own artefacts. I accepted every
confirmed finding. The corrections were made with figures read from the CSVs,
a new trial (M4) logs the corroboration that had been asserted, frozen
docstrings and run records were left untouched with a banner beside them, and
the audit itself is rows C151 to C154. The rule I take from it is mechanical
and now binding: no number enters a summary or the register except by
formatting from the table that holds it.

Two smaller decisions of mine on the same day. I withdrew Route C of the
Axin2 branch as out of focus, because a bulk array cannot answer a
co-occurrence question at any level of replication (row C146), and I directed
that Route A run rather than wait for my reading of England 2025; it then
turned out to rest on data never deposited (row C148). For the README I chose
to replace the Leiden UMAP with a figure of the register's own shape, drawn
from `CLAIMS.md` by a script, because that shape is what this repository is.
`RESEARCH_QUESTIONS.md` went through a second adversarial fact-check before it
landed (119 findings, 17 independently verified before the session limit, the
rest checked by hand); its confirmed defects, chiefly register statuses quoted
one grade too strong and a range typed from the wrong pass, and their fifty-eight
fixes are listed in PROGRESS item 36. My retain or reject on the outreach
framing is pending.

**25 · The deposits moved beside their papers, the README was left open, and
a cost rule for multi-agent work was set.** (2026-09-21.) I asked for the
analysis folders to move into `Thesis/`, one folder per dataset, or into the
paper folder when that paper produced the deposit. GSE262927 therefore sits
under paper 1. GSE178360 sits under a folder for Murthy 2022, which is outside
the roadmap, and that folder holds a pointer note only, because decision 23
stands: no study note before I have read the paper. I asked for the README to
stay open-ended, because there is more to go, and it now says what is not
done rather than reading as a finished portfolio. And I set a rule for the
assistant, recorded in `AI_CONTEXT.md` and in the cross-project harness file:
before launching a multi-agent workflow, weigh its cost; do the deterministic
part with a script; deploy agents only where a script cannot judge, and then
with at most three lenses. The fact-check of the outreach documents the day
before had used six checkers and forty-seven verifiers and hit the usage
limit for what a script catches faster. The move itself was done and checked
by script, with the validator, the compiler and a path-resolution check as
the gates.

## How outputs were reviewed

Every run writes its decisions to machine logs (`decisions.json`,
`analysis_log.txt`), and the per-dataset reports and pipeline record are
*generated* from those artefacts, numbers in the documentation cannot drift
from what was computed, and a claim I couldn't trace to an artefact was
treated as unverified and removed. Work was reviewed between sessions against
[`PROGRESS.md`](PROGRESS.md), landed through pull requests, and known issues
were carried forward in writing rather than dropped.

## Where the honest failures live

- The refuted Scrublet/AT0 over-removal claim, with both rounds of the test:
  [`docs/DOUBLETS_AND_SCRUBLET.md`](docs/DOUBLETS_AND_SCRUBLET.md)
- The three contradicted cluster annotations, flagged in
  `Thesis/gate1_01_niethamer_2025/GSE262927/tables/cluster_annotation_proposals.csv`
- Thresholds that never bound, doublet calls that are a ranking rather than a
  detection, and every other caveat: [`FINDINGS.md § 6`](FINDINGS.md#5--negative-results-and-self-audits)
  and the per-dataset reports
