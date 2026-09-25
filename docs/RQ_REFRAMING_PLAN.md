# Plan to restructure the research questions

**Implementation update, 25 September 2026:** the owner subsequently requested
the actual rewrite. See the [current register](../RESEARCH_QUESTIONS.md),
[measurement contracts](RQ_MEASUREMENT_CONTRACTS.md) and
[implementation record](migrations/2026-09-25-rq-reframing/README.md). All A1–A14
remain, without a four-question cap. The dated text below is the historical
planning proposal; its grouping recommendations are not the current register.

25 September 2026. **Review draft; no new analysis launched and no question or
claim status changed.** Based on the owner's `RQ_FRAMING_PROPOSAL.md`, the current
question/claim registers, source tables and figure-generator code. PRs #66 and
#67 are merged; this plan starts from that state on `codex/rq-framing-plan`.
The attached note is a proposal to evaluate, not an instruction to adopt its
assertions or decision menu.

Repository baseline: merged commit `f8c3de7`. The reviewed owner-provided note
has SHA-256 `5aaaeF621f3cb3af4021a871292e2d08b591759be02cd2ed81e4dc9264ffc51d`
(hexadecimal is case-insensitive); its original file remains unchanged outside
the repository. Documentation/link and existing evidence-binding validation
passed 2,109 checks for this planning update.

## Recommendation

First inventory the project's own biological observations and the distinct
hypotheses they support or motivate. Separate reusable measurement checks, then
choose the front-page grouping. **There is no four-question survival limit.**
The A1, A5, A10 and A14 cards below are worked examples of possible priorities,
not an exhaustive list or a rejection of the other biological propositions.
Preserve every existing A identifier and its evidence in a migration crosswalk.
Merge hypotheses only when they have the same discriminating test and decision,
not just a shared topic.

This is a substantive prioritization, not a declarative rewrite of all fifteen
headlines. A hypothesis can be valuable while untested. A technical result can
be useful without becoming a biological mechanism. The organising question
remains productive epithelial repair versus persistent remodelling; cancer is
a distinct comparison context, not an assumed later point on that trajectory.

The immediate analytical priority after the restructuring should be the bounded
A10 outcome-data join, alongside identification of a compatible direct-mark
comparison for A1. Do not begin another unrestricted score, depth or ligand-rank
sweep merely because it is computationally available.

## Make the project's own biological interpretation explicit

The attachment's ownership argument is central to this revision. The project
should state what its analyses reveal about biological systems and use those
observations to formulate its own hypotheses. A list of failures of labels,
references or databases is not an adequate account of the scientific contribution.
The first draft underemphasized this by leading with consolidation and caveats.

Use two distinct sentences in each entry: **what we observed**, then **the
biological hypothesis that observation supports or motivates**. State the
strength of that connection. Put the discriminating experiment next; attach
measurement rivals without making them the entire identity of the question.

| This project's observation | Biological interpretation to foreground | Test of the stronger hypothesis |
|---|---|---|
| Transitional-marker expression occurs in neonatal controls; ES1 also examines label-excluded programme enrichment | Part of the epithelial transition programme may be reused during alveolar maturation and adult repair | Shared components beyond generic stress/cycling in independent developmental and injury comparisons |
| The overlap-reduced HPCS score increases in non-neoplastic repair/developmental libraries and IPF pairs as well as lesion comparisons | Neoplasia-associated transcript programmes contain features also expressed during non-neoplastic epithelial remodelling; this motivates a shared-plasticity hypothesis | Identify reproducible shared components, then test context-specific additions separately |
| Both reference identity scores and labelled-reference contrasts change in Cebpa-mutant wells (C167) | Cebpa loss may attenuate identity across the AT2 population rather than selectively preserving intermediates | Replicated genotype-by-state comparison with a state definition independent of the scored genes |
| Broad ADI/AT1 scores have shared definitions, while the small late-marker component is unstable (C168) | A maturation model should distinguish a shared transition component from an additional endpoint-associated component | Test independent information about measured maturation; list overlap alone does not prove biological continuity |

The first three examples contain biological distributions that can motivate
positive hypotheses; the fourth begins with a definitional dependency and
provides weaker biological support. They should not all be called instrument
failures or assigned the same evidential strength. The observations remain the
project's results even when a proposed mechanism is unconfirmed. Reframing
changes scientific emphasis without changing numerical evidence or granting a
new hypothesis an old claim's status.

## Biological candidates retained beyond the four worked cards

This inventory is for selection, not a declaration that each hypothesis has
passed an inferential test. Data availability controls execution readiness;
it does not determine whether a biological question is worth retaining.

| Existing IDs | Biological proposition retained | Place in the revision |
|---|---|---|
| A1 | Regulatory programmes distinguish RNA-overlapping intermediates and relate to different measured outcomes | Candidate main question; direct epigenetic measurements remain central |
| A5 | A developmental transition component is reused during adult repair | Its own hypothesis, independently testable from the next row |
| A11 | Lesions add a reproducible epithelial component beyond the shared transition programme | Its own hypothesis; possible thematic grouping with A5 |
| A7 | Cebpa loss produces a population-wide identity shift versus a state-selective change | Biological modifier hypothesis; nesting under A1 is optional |
| A8 | A maturation-specific component adds endpoint information beyond a shared transition component | Biological hypothesis; overlap auditing is a method prerequisite |
| A3 | Late injury-associated macrophage programmes differ from age-matched uninjured programmes | Temporal hypothesis; persistence of the same cells is a separate tracing claim |
| A6 | IPF changes a programme within a comparable macrophage state beyond changing subtype proportions | Cell-intrinsic versus composition hypothesis, not eliminated by moving harmonization checks |
| A12 | Recipient receptor/inhibitor context explains response differences beyond ligand RNA alone | Recipient-context hypothesis; association and causal reception need different tests |
| A13 | Fibroblast programmes add information about epithelial state beyond macrophage IL1B | Distinct joint-model hypothesis, conditional on complete triads and adequate precision |
| A10 | Molecular responses add information about measured organoid outcomes | Predictive hypothesis with a bounded public-data pilot |
| A14 | Exposure duration affects recovery after withdrawal; fibroblast reception modifies recovery | Two hypotheses with separate decisions, potentially tested in one factorial experiment |
| A4 | Wnt-associated maintenance and IL-1-responsive transition occur sequentially within a traced population | Conditional temporal hypothesis requiring activity/history measurements |
| A2/A9 | Source contributions and receiver competence vary with biological context | Optional functional hypotheses if exposure/receiver endpoints can be specified; resource coverage alone remains method work |

A12-S1 remains an enabling source-identity question until a biological state can
be defined independently of reference-label failure. None of these candidates
must be discarded to reach a predetermined number of headings. The final RQ
count should follow nonredundant predictions and scientific priority.

## What to retain and correct in the attached proposal

Retain its proposition-first presentation, compulsory rivals, explicit
discriminating observations, retirement criteria and separation of method work.
Retain the prohibition on upgrading existing evidence through new wording.

Three parts need correction:

1. **Falsifiability is not restricted to prediction headlines.** Current A6,
   A7 and A8 already specify tests that could weaken an explanation. Descriptive
   and methodological questions can also close under bounded designs. The
   problem is their prominence, scope and connection to the organising biology.
2. **A new proposition does not inherit an old claim's status.** Refuting
   injury-exclusive Krt8/Cldn4 expression motivates developmental reuse; it does
   not validate that stronger hypothesis. Existing C119 remains refuted as
   written. New hypotheses start explicitly untested or provisionally motivated,
   while each linked evidence row retains its own status.
3. **Several replacements go beyond reframing.** Temporal chromatin ordering,
   lineage persistence, receptor dimers and absence of a cancer-specific
   programme require different observations. They cannot be presented as the
   biological interpretation already established by the same numbers.

The attachment also predates the merged expansion of A1: its current remit
already includes direct histone marks, lineage and functional evidence. Its
reference to PROGRESS item 41 is historical. The new first batch adds measured
PATS endpoint reconstruction and IRE1α RiboTag results, not evidence for a
general regulatory taxonomy. See the [current A1 report](../RQ_Specified/A1_transitional_epithelial_state_distinction/reports/FIRST_BATCH_REPORT.md).

## Four worked hypothesis cards, not a survival shortlist

These are retrospective hypothesis formulations informed by inspected data.
Future tests require dated specifications and validation units not used to
choose the hypothesis. Directions, clinically/biologically meaningful margins
and power cannot be supplied by changing the headline.

### A1 — Regulatory states and differentiation potential

**Question:** Do overlapping epithelial transition RNA programmes conceal
distinct regulatory states with different differentiation potential?

**Working hypothesis:** Within a defined injury/model context, epithelial
intermediates with overlapping RNA markers contain reproducible chromatin
programmes associated with different independently measured maturation or
persistence outcomes. Regulatory measurements add information beyond the shared
RNA programme. Accessible regions and each measured histone mark are separate
features; a universal DATP/PATS/HPCS classifier is not presumed.

**Why it is live:** Existing multiome results reveal limitations of marker-only
definitions; they have not demonstrated the hypothesis. C167 motivates a
genotype-by-state alternative, and C168 identifies dependence between gene sets.
The PATS source reconstruction demonstrates an accessible measured endpoint.
Neither it nor the separate IRE1α RNA response establishes chromatin–fate coupling.

**Rivals:** A shared stress/maturation continuum; generic cycling or identity
loss; source selection, depth or genotype changing the reference; different
assays and cohorts masquerading as state differences.

**Discriminating test:** First quantify a compatible, controlled direct-mark or
accessibility contrast within a model using independent preparations. Then test
whether frozen regulatory features improve association/prediction of measured
AT1 protein/morphology or labelled AT1 contribution beyond RNA and design
covariates in an independent cohort or linked experimental units. A cross-study
juxtaposition cannot substitute for the second test. No verified joint
regulatory-and-fate dataset currently supports claiming this full test ready.

**What would weaken or retire it:** A sufficiently precise validation shows no
meaningful regulatory distinction or no incremental outcome information within
the declared assay/context. Those are separate possible failures. A compatible
assay positive control must work, and confidence bounds must exclude a
prespecified useful effect; sparse pooled data or failed coverage are inconclusive.

**Practical scope:** Compatible histone quantification is the next molecular
task. Keep AT2 silencing-before-closure as one optional temporal mechanism to
test if suitable longitudinal assays are found, not the entire A1 proposition.
The current promoter/depth results do not establish distal permissiveness or
temporal order. Cebpa interaction and AT1-unique endpoint information become
subquestions here rather than separate front-page RQs.

### A5 — Developmental reuse and context-specific additions

**Question:** Is part of adult epithelial repair a redeployment of a developmental
transition programme, and what additional programmes distinguish persistent or
neoplastic contexts?

**Working hypotheses, tested separately:** A5a: a label-independent component of
the neonatal alveolar transition programme recurs in adult injury beyond generic
stress and proliferation. A5b, retaining the A11 cross-reference: lesion-associated
epithelium contains a reproducible component beyond that shared programme.
Support for one does not establish the other.

**Why it is live:** C119 shows that the two-transcript call is not injury
exclusive, while ES1 examines fuller label-excluded programmes. The reduced-HPCS
score is also elevated across several contexts. These motivate programme reuse
and context-specificity tests; they do not establish common identity or absence
of cancer-specific biology. Technical seeds and pooled libraries do not add
independent biological replication.

**Rivals:** Generic stress/cycling; overlap among gene definitions; mixtures of
different cell states; genotype or preparation confounding with age/context;
species mapping and cohort-specific measurement.

**Discriminating test:** Freeze label-excluded developmental, injury and general
stress/cycle modules, preserving source definitions and mapping loss. Test the
shared component in replicated wild-type neonatal and adult-injury comparisons.
Evaluate any proposed lesion-specific addition on independent patients/cohorts,
using a shared-programme baseline and keeping selection and evaluation separate.
Keep developmental and lesion comparisons as distinct tests rather than one
cross-species trajectory or pooled classifier.

**What would weaken or retire it:** A5a weakens if the declared shared component
does not recur with adequate precision after the stated competing programmes
are accounted for. A5b weakens if a frozen addition fails independent evaluation
with precision sufficient to exclude useful discrimination. Neither a shared
score nor an underpowered null establishes that no other programme exists.

**Practical scope:** Reuse existing effects as hypothesis-generating displays.
Seek independent, adequately sampled comparisons before extending the classifier.
The existing one-library conditions cannot establish a population interaction.

### A10 — Molecular responses linked to an observed functional phenotype

**Question:** Do epithelial and fibroblast molecular responses explain
perturbation-associated organoid growth beyond baseline imaging and experimental
batch?

**Working hypothesis:** A small frozen programme set carries incremental
information about day-14 organoid area conditional on day-7 area, reproducible
across held-out biological preparations. Epithelial-only and epithelial-plus-
fibroblast models answer different, prespecified comparisons.

**Why it is live:** GSE307112 has deposited well-linked RNA and imaging resources;
the [existing dataset gate](NEXT_DATASET_GATE.md) identifies the necessary joins.
This is a route to an endpoint beyond another expression score. The relationship
has not been tested in this repository.

**Rivals:** Baseline size, plate/preparation effects, guide/target identity,
viability or cell composition, species-assignment errors and leakage among wells
from the same preparation.

**Discriminating test:** Join plate design, imaging, species-assignment QC and
counts; resolve guide, plate and biological-preparation IDs. Verify RNA harvest
timing: contemporaneous RNA can support held-out association/prediction, not a
claim of forecasting future growth. Compare with baseline-imaging/plate models
using complete preparation holdouts and an appropriately small feature set.
Preserve known perturbation effects as positive controls, not new discoveries.

**What would weaken or retire it:** Incremental held-out performance is bounded
below a prespecified practically useful improvement, with adequate independent
preparations and working controls. Failure to establish preparation identity
stops the inferential pilot. Wells cannot repair that missing replication.

**Practical scope:** Highest-priority public-data feasibility pilot after the
RQ revision. Its first deliverable is a verified join and eligibility decision,
not a fit. Organoid growth/morphology is the endpoint; mature AT1 fate, in vivo
repair and IL-1-specific causality require other measurements.

### A14 — Signal duration, recipient context and recovery

**Question:** Does sustained IL-1 exposure impair epithelial recovery after
withdrawal, and does fibroblast reception contribute independently of epithelial
reception?

**Working hypotheses, tested separately:** Longer exposure decreases recovery
of a specified mature epithelial endpoint after withdrawal; fibroblast-specific
IL1R1 perturbation modifies that epithelial recovery while epithelial exposure/
reception is controlled. Current Wnt activity can be an additional mechanism only
if a specific perturbation/temporal comparison is included.

**Why it is live:** The review-motivated RNA analyses show context dependence and
motivate recipient-centred testing. A12/A13 supply observational components and
alternative models, not demonstrated signalling or feedback. Published tracing
and intervention studies inform controls; the local datasets do not measure
the requested withdrawal outcome.

**Rivals:** Direct epithelial signalling alone; altered cell survival or entry
into the transitional state; incomplete washout, residual ligand or inadequate
target engagement; stromal changes unrelated to IL-1 reception.

**Discriminating test:** Verified transient/sustained exposure and washout,
time-matched controls, recipient-specific perturbations and independent cultures/
animals. Measure viability, mature-cell yield/function and labelled descendants,
with direct evidence of exposure and target engagement. Use a factorial or
otherwise identifiable design to separate epithelial and fibroblast effects.

**What would weaken or retire it:** With adequate precision and positive
controls, exposure duration has no meaningful effect on recovery, or fibroblast
perturbation has no incremental effect under controlled epithelial reception.
These retire the corresponding component, not both automatically. A failed
perturbation or insufficient follow-up is inconclusive.

**Practical scope:** Mechanistic follow-up requiring suitable external or new
experimental data. Do not launch a local RNA proxy analysis and label it the
withdrawal test. Cancer conversion is outside this endpoint.

## Complete migration crosswalk

Keep stable A IDs and anchors; do not rename study folders, historical scripts,
claim IDs or frozen runs as part of this documentation change.

| Current entry | Proposed destination | Biological content retained | Content moved out of the main question |
|---|---|---|---|
| A1 | Core A1 | Regulatory distinctions and independently measured outcomes | Depth/reference diagnostics linked as supporting contracts |
| A2 | Conditional source hypothesis supporting A14/niche work | Context-dependent AREG source contributions, if paired with a defined receiver endpoint | Resource ranks, expression eligibility and molecule-budget checks; no asserted epithelial dominance |
| A3 | Conditional myeloid branch with A6 | State abundance and within-state remodelling over measured time | Confound/mixture audit; individual-cell persistence reserved for actual tracing |
| A4 | Temporal/signalling subquestion of A14 | Current Wnt activity, IL-1 response and transition order in a defined model | Transcript co-detection and reporter-proxy eligibility |
| A5 | Core A5 | Developmental reuse, with label-free tests | Marker exclusivity and sampling checks become supporting evidence |
| A6 | Conditional myeloid branch with A3 | A cell-intrinsic programme shift beyond subtype fractions | Harmonization, composition and correlation-model diagnostics |
| A7 | Context/modifier subquestion of A1 | Genotype-by-state effects on identity | Reference selection and one-well limitations |
| A8 | Outcome-specific subquestion of A1 | Whether AT1-unique features add mature-endpoint information | List overlap and score dependence |
| A9 | Supporting receiver contract; optional receptor-function experiment | Receiver competence only with an appropriate measured response | Receptor representation, resource coverage and RNA detection; no dimer inference |
| A10 | Core A10 | Independently evaluated molecular–organoid outcome association | Data joins and sample-ID gate remain linked prerequisites |
| A11 | Core A5, distinct context-specific subhypothesis | Additional lesion-associated programmes beyond a shared component | Overlap/annotation safeguards; no universal absence claim |
| A12 | Conditional recipient model supporting A14 | Receptor/inhibitor context versus ligand-output explanations | RNA compatibility scores, enrichment model sensitivity and ligand eligibility |
| A12-S1 | Source-annotation measurement contract | Unassigned source candidates can motivate a later phenotype study | No biological state defined by failure of a chosen annotation threshold |
| A13 | Conditional niche model supporting A14; relevant comparator in A10 | Fibroblast programmes adding information beyond macrophage IL1B | Triad completeness and estimable joint-model requirements |
| A14 | Core A14 | Withdrawal, recipient-specific effects and measured recovery | Experimental design stays explicitly prospective |

The A3/A6 branch is retained; its front-page priority has not been decided.
A persistent population and persistence of the same cells are different targets.
For a within-state test, avoid defining “noncycling” solely through the very
cell-cycle score being tested. Age-matched controls, harmonized populations and
independent biological units are prerequisites for the corresponding inference.

AREG and Wnt biology are not discarded. They are candidate source/recipient or
temporal mechanisms with their own tests, rather than automatically expanded
into a separate main question by the availability of another rank plot.

## Specific proposed assertions that should not become headlines

| Attachment assertion | Assessment against current evidence | Safer use |
|---|---|---|
| AT2 RNA loss occurs over permissive distal chromatin, with closure later | C131/C133 do not establish distal state or time order; the old genotype comparator is not a clean temporal control | Optional A1 temporal hypothesis with explicit longitudinal requirements |
| Epithelium sets the AREG available to fibroblasts | C37 is measurement-sensitive; C45 supports competing myeloid sources; RNA is not secreted availability | Context-specific source/receiver experiment, not an existing supported conclusion |
| Late interstitial macrophages are the same persisting population | Cross-sectional proportions do not identify ancestry or replacement | Population-state persistence and cell persistence become separate estimands |
| AT1 identity is largely an extension of transition because 119 genes overlap | C168 quantifies shared definition, not developmental continuity or fate | Test added endpoint information from disjoint/shared components |
| Fibroblast EGFR is predominantly homodimeric | C114 is refuted as a resource-absence claim; receptor coexpression cannot measure dimer composition | RNA competence screen; direct complex/activation assay for the mechanistic claim |
| Neoplasia adds no separable epithelial programme | Shared-HPCS elevation does not exclude an additional programme | Test a bounded candidate addition on independent validation data; report an inconclusive null honestly |
| Unassigned IL1B-positive cells form a distinct source architecture | Annotation failure is partly determined by reference and threshold; ambient/mixed profiles remain alternatives | Resolve source identity before formulating a biological state claim |

These objections do not prevent proposing the mechanisms. They prevent treating
the proposal's motivating measurements as direct tests of those mechanisms.
Novelty remains to be checked against the literature for whichever narrow
hypotheses are retained; reframing alone creates neither novelty nor evidence.

## Proposed document and figure structure

`RESEARCH_QUESTIONS.md` remains the single biological-question index. Its opening
table should show the selected biological hypotheses, endpoint, evidence state,
execution readiness and next discriminating task. Follow with compact hypothesis
cards, then conditional branches and an A-ID crosswalk. Implementation detail
stays under `RQ_Specified/`; paper-specific evidence stays under `Research Article/`.

Create `docs/RQ_MEASUREMENT_CONTRACTS.md` during implementation as a concise
index of reusable checks, linking existing authoritative methods and reports
rather than duplicating their protocols: biological units/confounding; depth/
reference definition; source annotation; ligand/resource representation;
gene-set dependence/inference; figure/provenance ownership. Each core hypothesis
links only the contracts that could alter its decision. Completed checks have
bounded applicability and a stopping rule; they are not perpetual analyses.

Use the attachment's five fields, augmented with a short execution footer:
**hypothesis; motivating evidence with unchanged status; rivals; discriminating
observation/model and unit; retirement/inconclusive criteria; data readiness and
next action.** Keep hypothesis status separate from claim status. Avoid describing
methodological uncertainty as the sole biological motivation.

Figure selection follows the test:

| Core question | Primary scientific figure when data support it | Current material and limitation |
|---|---|---|
| A1 | Controlled direct-mark/accessibility effects plus measured maturation/lineage endpoint; incremental outcome information if linked data exist | New PATS and IRE1α panels are separate assay/study observations; historical depth panels are supporting diagnostics |
| A5 | Independent-unit shared versus context-specific programme effects, then held-out evaluation | Existing developmental/lesion figures motivate tests; their libraries, patients and genes retain their actual dependence |
| A10 | Held-out observed-versus-predicted outcome, effect uncertainty and comparison with the baseline model | Current figure is a design schematic; no invented growth result |
| A14 | Exposure/withdrawal design, replicate-level recovery time course and recipient-specific outcome contrasts | Current figure is a prospective schematic; outcome panels need new data |

UMAP localizes heterogeneity and PCA diagnoses sample structure; neither is a
mandatory main figure or evidence of the mechanism. Resource-coverage and depth
panels belong beside the measurement contracts or source galleries. They remain
available, with their warnings, when readers assess a result.

### Repair caption ownership before any redraw

The attachment's defect is confirmed in `16_research_question_figures.py`:
`write_captions()` still inserts missing A1–A5 blocks by heading; the active
register has no A2 marker and uses script 18's newer A2 figure. A rerun could
reintroduce the old A2 panel and overwrite current caption wording. Script 16's
A3 title also says “after repair” without a measured recovery endpoint.

Recommended ownership: the root question narrative is authored; generators own
figure files and generated caption data/fragments in the appropriate gallery.
Remove implicit edits to the root register from a default render. Preserve
numeric bindings to saved tables and identify one active generator per figure.
Do not solve the drift by abandoning provenance or manually changing numeric
captions without bindings.

First synchronize active caption content and ownership rules, retire script
16's A2 output from the current gallery, and add explicit figure-only selection.
Then redraw only figures whose titles/captions need correction from verified
cached inputs or saved tables. Rename A3's title to describe sampling after
infection, without asserting recovery or cell persistence. Record new figure
hashes while preserving original run records and scientific-table bytes. A1's
historical title can remain a diagnostic title if its scope is explicit; no
cosmetic full-data rerun is necessary. No generator is run during this planning
turn.

## Implementation sequence and acceptance criteria

| Step | Concrete deliverable | Acceptance / stopping rule |
|---|---|---|
| 1. Freeze the crosswalk | Snapshot current A headings, evidence links, figure ownership and C-row statuses; classify every entry using the table above | All A1–A14 plus A12-S1 accounted for; no deleted evidence or assumed human acceptance |
| 2. Rewrite the front page | Selected hypothesis cards and conditional branches, with linked measurement contracts; no fixed question quota | Each has one bounded primary test, biological unit, comparator, failure criterion and explicit readiness; nested subhypotheses have separate decisions |
| 3. Synchronize plans and handoffs | Update A1 workspace scope, A10 pilot pointer, PROGRESS, AI_CONTEXT, README navigation and structure documentation | Core/public-data questions distinguish exploratory evidence from pending endpoint tests; thresholds and finished runs unchanged |
| 4. Repair figure ownership | Safe generator/caption ownership, current A2 retained, corrected A3 title and selected displays | Rerender cannot restore retired blocks; scientific table hashes identical; figure-specific provenance and visual review pass |
| 5. Review before analysis resumes | Compact change report: crosswalk, selected hypotheses, data gaps and next approved tests | No fit launched simply to populate a new headline; hypotheses formulated from inspected data labelled retrospective |
| 6. Resume targeted work | A10 join/replicate gate; compatible A1 regulatory input audit; later A5 validation and A14 endpoint acquisition | Each task answers a discriminating question or reaches a documented data stop; source/receiver models require their existing coverage floors |

Implementation needs document/link checks and targeted tests for generator
ownership and unchanged numerical inputs. It does not require rerunning the
completed biological analyses. No new threshold, equivalence margin or effect
size is chosen in this framing plan. Existing replication and joint-model floors
remain, with power/precision justified for each new test.

The present deliverable is this plan. The canonical RQ register, original
proposal, figures, analysis code, scientific tables and claim statuses remain
unchanged. Remaining analysis is paused while the owner considers the revised
research framing.
