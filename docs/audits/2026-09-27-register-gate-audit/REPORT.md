# Register data-gate audit, 27 September 2026

Every remaining question in the register states a gate, and for most of them the gate is
data the repository does not hold. This audit asks once, for each of them, whether a
public deposit meets that gate. It adds no claim, grades nothing and fits no model.

**Method.** Two declared queries per question against NCBI GEO DataSets, restricted to
series, with the top six summaries of each kept. The queries, hit counts and every
candidate are in [`search_run.json`](search_run.json) and
[`candidates.tsv`](candidates.tsv), written by
[`run_gate_search.py`](run_gate_search.py), which retrieves and judges nothing. The
verdicts below are written by hand against that table and each names the accessions it
rests on.

**A15 is excluded on purpose.** Its gate was audited by the session that owns it, which
searched GEO, PRIDE, the Image Data Resource and the BioImage Archive and recorded the
result in `RQ_Specified/A15_epithelial_integrin_tgfb_activation/reports/PUBLIC_DATA_SEARCH.md`.
Repeating it here would waste effort and risk a contradictory verdict.

**Titles are stored verbatim.** Three candidate titles in the table carry an en dash
because the archive wrote them that way. They are retrieved evidence, so they are not
edited; the repository's own prose in this report uses none.

**Three limits on every verdict.** A hit is a candidate, not an eligible cohort: eligibility
needs per-sample metadata and cell counts, which this pass did not download. An absence
bounds these queries, not the archive. And GEO indexes sequencing deposits, so a gate
asking for protein, imaging or physiology is one GEO is the wrong place to answer.

## Verdicts

| Question | Verdict | Rests on |
|---|---|---|
| A3 | **Candidate meets the gate as stated** | GSE303646 |
| A5 | Candidate found, eligibility unverified | GSE264278, GSE201698 |
| A11 | Candidate found, eligibility unverified | GSE198864 |
| A13 | Candidate found, eligibility unverified | GSE233844, GSE122960 |
| A1 | Candidate found, but it does not close the gate | GSE215824, GSE286986 |
| A4 | Still blocked, constraint nameable | nothing new |
| A7 | Still blocked, constraint nameable | nothing in lung |
| A8 | Still blocked, wrong archive | GSE27964 only |
| A9 | Still blocked for protein | nothing |
| A12 | Still blocked, constraint nameable | nothing |
| A14 | Still blocked, nearest is underpowered | GSE295566 |

## The one gate a deposit meets

**A3 asks whether injury leaves a macrophage programme beyond normal aging, and its card
needs age-matched controls and comparable sampling.** GSE303646 is 56 mouse samples titled
as single-cell RNA sequencing of lung regeneration in young and aged mice after injury.
That is the design the gate names: both ages, injury, one study, so sampling is comparable
by construction, and 56 samples leaves room for animal-level units.

**One warning belongs with it, and it is the same warning that retired A6.** A3 as worded
asks about a programme beyond aging, and the boundary between a changed cell composition
and a changed within-state programme is set by how finely the macrophage states are
annotated. A deposit does not fix that. A3 becomes worth running only if its endpoint is a
within-state programme with states defined independently of the programme being tested, and
with the injury-by-age interaction as the estimand rather than a main effect. Without that,
the answer will move with the annotation convention, which is why A6 was rejected.

## Three candidates that need an eligibility check before anything else

**A5 needs another eligible adult injury cohort for generality.** GSE264278 is 18 mouse
samples with a time course of both bleomycin and lipopolysaccharide injury, and GSE201698
is a six-sample bleomycin time course. The check is whether either holds at least three
animals per time point with type 2 capture, which is what A5's frozen signature needs. Note
that GSE141259, also returned, is the cohort A5 already used.

**A11 needs a non-neoplastic injury comparator with paired patients and epithelium.**
GSE198864 is 91 human lung samples on limited permissiveness for SARS-CoV-2 with
virus-induced expansion. It is human non-neoplastic injury with epithelium and it is not in
this repository, so it does not carry the discovery contamination that the two IPF cohorts
on disk do. The check is paired structure and epithelial cells per donor. The queries found
nothing better, and the second A11 query returned no hits at all, which is worth recording
as a bound on the search rather than on the archive.

**A13 needs at least ten patients with complete macrophage, fibroblast and epithelial
triads.** GSE233844 is 38 human samples on progressive idiopathic pulmonary fibrosis, and
GSE122960 is 17. Both are candidates for triad counting at the 50-cell floor. GSE136831,
also returned, is already used. The check is arithmetic on per-patient cell counts and it
can fail cleanly, which is the right shape for a first step.

## The gate a candidate does not close

**A1 needs matched replicated regulation-to-fate linkage.** GSE215824, eight mouse samples
on alveolar epithelial progenitors driving regeneration through chromatin topology, and
GSE286986, a 22-sample lung multiome, both supply a regulatory measurement. Neither links
that measurement to a fate outcome in the same replicated animals, which is the part A1 has
been missing all along. So the honest verdict is that the regulatory half is available and
the linkage half is not, and A1 stays where it is.

## The gates that remain blocked, with their constraints named

- **A4** needs measured activity history and a lineage-linked response. Every candidate
  returned is either already in this repository or another cross-sectional snapshot, and a
  snapshot cannot carry history. The constraint is a recorder or a serial-sampling design,
  not a deposit.
- **A7** needs replicated age-matched animals with a Cebpa genotype and a state contrast.
  The Cebpa queries returned human activating-RNA and leukaemia studies and nothing in
  lung. The current data is one well per condition, so the interaction A7 asks about is not
  estimable, and no public deposit changes that.
- **A8** needs independent mature endpoints. GSE27964 carries clonogenicity and
  myofibroblastic differentiation, which is functional, but it is an old bulk series. The
  real constraint is that a mature endpoint is morphometry, physiology or protein, and GEO
  is the wrong archive to settle it.
- **A9** needs receptor protein or activation data on fibroblasts. No lung fibroblast
  surface-protein deposit surfaced. The card's own reading stands: the RNA screen is
  possible and the protein evidence is not available.
- **A12** needs measured interleukin-1 pathway activation. Every candidate measures
  transcript after a perturbation, which is what the card already says is insufficient.
- **A14** needs signal withdrawal with a recipient-specific intervention. GSE295566 is the
  nearest, on persistence of a fibroblast population during delayed resolution, but it is
  four mouse samples and has no recipient-specific arm.

## What this audit changes

Before it, eleven questions were gated on unspecified data. After it, one has a deposit
that matches its gate as written, three have named candidates whose eligibility is an
arithmetic check, one has half its gate met, and six are blocked for a reason that can be
stated in a sentence. That is the whole value of the pass: a blocked question with a
nameable constraint can be retired or re-specified, while a blocked question with a vague
gate sits in the register forever.

The recommended order, on biological content rather than on convenience: settle A13's triad
arithmetic first, because it is a counting exercise that can fail in an afternoon; then
A11's comparator eligibility, because A11 already has a positive replicated result waiting
on specificity; then A5's cohort eligibility. A3 should not be run until its endpoint is
re-specified as a within-state interaction, or it will inherit A6's defect. Nothing here is
a decision; the owner holds question selection.
