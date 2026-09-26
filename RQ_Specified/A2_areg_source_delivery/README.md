# A2: does the fibroblast response to AREG depend on delivery or on abundance

**Status, 26 September 2026: framing declared, stage 1 authorized, nothing
scored.** No endpoint has been computed in either dataset named below, and no
claim row has changed.

This analysis replaces A2's abundance question with a delivery question. The
register card is [A2](../../RESEARCH_QUESTIONS.md#a2); the biology and the closed
abundance record are in [RATIONALE.md](RATIONALE.md); the stages and their stop
rules are in [PLAN.md](PLAN.md); the frozen decisions are in
[`config/a2_delivery_contract.json`](config/a2_delivery_contract.json).

## The hypothesis in one sentence

AREG's contribution to a fibroblast response is set by where the ligand is
released relative to a competent recipient, not by which compartment transcribes
the most of it.

## Why the previous question is closed

Seven register rows answer the abundance version, and none supports a source
hierarchy: C37, C39, C45, C40, C48, C49 and C50. C51 records why, in this
repository, a donor-level correlation on these variables is hard to read at all:
the one significant pair in that trial tracked sequencing depth and a frozen rule
refused it. The rationale lists each row with its status. None of them is
re-graded here.

## What makes the question testable now

The mechanism is short-range and recipient-licensed. Amphiregulin activates
integrin alphaV on mesenchymal stromal cells and releases bioactive TGF-beta from
latent complexes, driving myofibroblast differentiation
([Minutti 2019](https://doi.org/10.1016/j.immuni.2019.01.008)), and the
fibroblast arm of TGF-beta signalling needs amphiregulin
([Zhou 2012](https://doi.org/10.1074/jbc.M112.356824)). A ligand that converts a
store the recipient already holds predicts that tissue-level ligand abundance
will not track the response, which is what C49 and C50 found.

Two legs follow, each able to fail alone:

| Leg | Prediction | Data | Unit |
|---|---|---|---|
| 1 | Removing the epithelial ligand lowers a frozen fibroblast TGF-beta programme; removing epithelial receptors does not | GSE307112 organoid knockout screen | 4 plate-replicate units on one plate |
| 2 | The fibroblast response tracks the recipient's own TGF-beta activation machinery better than it tracks epithelial ligand | GSE136831, the trial E6 instrument reused | donor |

## Why the screen can separate the rivals

The screen perturbs the mouse epithelium only and leaves the human fibroblasts
unedited, with reads assigned by species. Removing the ligand removes what the
fibroblast can receive. Removing epithelial EGFR, ERBB2, ERBB3 or ERBB4 removes
only epithelial reception. Removing epithelial ITGB6 removes epithelial TGF-beta
activation. Those three contrasts distinguish the delivery hypothesis from the
autocrine rival and from the epithelial-activation rival, and all six targets sit
on one plate with one well per target in each of its four replicate units.

## The gate in one sentence

Stage 1 asks whether the plate layout, the knockout validation, the endpoint gene
coverage, the depth profile and the attainable power support the test at all, and
it computes no endpoint value while asking.

## Layout

| Path | Contents |
|---|---|
| `RATIONALE.md` | the biological argument, the closed abundance record, the other-layer verdicts |
| `PLAN.md` | six stages, their stop rules and the order of work |
| `config/a2_delivery_contract.json` | endpoints, adjustment, statistic, thresholds, prohibitions |
| `tables/` | stage outputs, written with run records; created when stage 1 runs |
| `reports/` | stage reports; created when stage 1 runs |

## Three things a later session must not do

1. **Do not compute either endpoint by target before stage 2 is committed.** The
   pre-registration is the only thing that makes four wells per target readable,
   and stage 1 includes a precedent check that records whether any such contrast
   already exists.
2. **Do not read a null as absence.** The Areg knockout lowers the mouse
   transcript by 1.042 log2 CPM and leaves it at 6.226. The contract declares
   precise absence unavailable at this design, before any test.
3. **Do not upgrade either leg past its evidence.** Leg 1 is a within-screen
   descriptive association while preparation independence is unresolved, and leg
   2 is correlational with no direction and no proximity, and may be refused by
   its own depth control.
