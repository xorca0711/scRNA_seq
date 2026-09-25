# A10 rationale: why link epithelial programmes to measured organoid growth

25 September 2026. Written before any model is fitted. This document carries the
biology and the argument. The analysis structure is in [PLAN.md](PLAN.md), and the
machine-readable decisions in
[`config/a10_outcome_contract.json`](config/a10_outcome_contract.json).

## The biology this question sits in

The alveolus is where gas exchange happens, across a barrier one cell thick. Two
epithelial cells build it. Alveolar type 1 cells are enormous and thin and cover
almost the whole surface. Alveolar type 2 cells are compact, make surfactant, and
act as the facultative progenitor: they self-renew and give rise to type 1 cells
after injury. A competent progenitor subset is distinguishable inside the type 2
population rather than every type 2 cell being equivalent
([10.1038/nature25786](https://doi.org/10.1038/nature25786)).

**The decisive point for this question is that type 2 cells are not autonomous.**
Their stemness is held by signals from neighbouring mesenchyme. Single Wnt-producing
fibroblast niche cells maintain type 2 stemness, and leaving that niche is what
permits type 1 differentiation ([10.1126/science.aam6603](https://doi.org/10.1126/science.aam6603)).
The same Wnt axis drives alveolar formation in development
([10.1016/j.celrep.2016.11.001](https://doi.org/10.1016/j.celrep.2016.11.001)), and
targeting Frizzled receptors can push alveolar regeneration
([10.1016/j.cell.2023.05.022](https://doi.org/10.1016/j.cell.2023.05.022)).
Inflammatory signals from other neighbours also change what type 2 cells do: IL-1
and TNF contribute to a niche that enhances regeneration
([10.1016/j.stemcr.2019.02.013](https://doi.org/10.1016/j.stemcr.2019.02.013)), and
IL-1 drives type 2 cells into a transitional state on the way to type 1
([10.1016/j.stem.2020.06.020](https://doi.org/10.1016/j.stem.2020.06.020)).

So alveolar repair is a conversation, not a solo performance. A gene knocked out in
a type 2 cell can change repair in two different ways: by changing what that cell
can do, or by changing what it asks of its niche and how the niche answers.

## Why an organoid screen can separate those two

The alveolosphere assay is that conversation in a dish. Type 2 cells only grow into
organoids when cultured with fibroblasts, so the assay has the dependency built in.

This deposit adds two features that make it unusually informative:

- **Two species in one well.** Mouse type 2 cells are cultured with human
  fibroblasts. Reads can therefore be assigned to a species, so one well yields the
  epithelial transcriptome and the fibroblast transcriptome separately. The
  perturbation is made in the mouse cells only, before co-culture. So the epithelial
  side reads the cell-autonomous consequence, and the fibroblast side reads how the
  niche responded to a changed neighbour. Very few designs let you read both halves
  of a dialogue from the same sample.
- **A measured outcome.** Imaging gives organoid count, mean area and area
  proportion at day 7 and at day 14.

That second feature is what makes A10 different from every other question in this
register. A1, A5, A6, A11 and A12 all compare RNA states with other RNA states. Here
a programme can be set against something the experiment actually measured about
tissue building. The deposit holds 886 RNA libraries and 203 perturbation targets.

## The hypothesis, and what would make it interesting

**Hypothesis.** The epithelial programme state in a well carries information about
how much organoid tissue that well built, beyond what baseline size and plate batch
already explain. Separately, the fibroblast response programme may carry a further
increment.

The second half is the more interesting one. If fibroblast response adds information
after the epithelial state is accounted for, that is evidence that the niche's
reaction relates to the outcome, not merely the epithelial cell's condition. Given
that type 2 stemness depends on fibroblast Wnt, that is the expected shape of the
biology, and it has not been tested here.

## Logical flow, step by step

1. Type 2 stemness depends on the fibroblast niche, so repair capacity is a joint
   property of two cell types.
2. A perturbation in the epithelium can therefore act on growth through the
   epithelial cell or through the niche's response to it.
3. Neither route can be distinguished from RNA alone, because RNA has no outcome in
   it. An assay with a measured growth endpoint is required.
4. This screen has one, and its species separation gives both sides of the dialogue
   from the same well.
5. **But a well is not an animal, and not necessarily an independent preparation.**
   Before any model, the joins between design, imaging, species QC and counts must
   hold, and the biological unit must be established from evidence.
6. Only then is the question answerable: does the epithelial programme add
   information about growth beyond baseline size and batch, and does the fibroblast
   response add more?
7. Evaluation must hold out whole preparations. Sibling wells from one preparation
   share everything that matters, so splitting them would leak.

Step 5 is the gate, and it is where this work starts. The first deliverable is a
verified join and an eligibility decision, not a fit.

## What the screen cannot be made to say

- **RNA is harvested at the outcome time.** Both come from day 14, so any
  association is concurrent. This cannot become a claim about forecasting growth.
- **Organoid area measures growth and morphology.** It is not mature type 1 fate,
  not barrier function, and not repair in a living lung.
- **A well is not a mouse.** The deposit has 886 libraries; it does not have 886
  animals. Counting wells or guides as replicates is the error this design most
  invites.
- **Known results are controls, not discoveries.** The screen includes targets whose
  effects are already published. Recovering them validates the assay and must be
  reported as validation.
- **A null for one target may be a failed knockout.** Without editing efficiency per
  well, absence of an effect is not absence of a role.
- **Fibroblast response is not fibroblast causation.** Reading a changed fibroblast
  transcriptome shows the niche responded, not that the response caused the growth
  difference. That needs a fibroblast-side perturbation.

## Two connections to work already in this repository

Both were found by checking the target list against frozen modules, and both are
opportunities rather than results.

**The screen perturbs genes that A5 and A11 only score.** The contract's frozen
modules and this screen's targets overlap:

| Frozen module | Genes | Perturbed in this screen |
|---|--:|---|
| Shared, development and injury | 5 | Areg, Itgb6 |
| Development-specific | 94 | 11, including Egfr, Erbb3, Fzd5, Tfcp2l1 |
| Lesion-specific | 91 | Tgif1, Tigit |
| Shared, injury and lesion | 7 | none |

This is a weak and indirect cross-check, and it must be labelled as such. Whether
knocking out a gene changes organoid growth says nothing directly about whether a
developmental programme is reused in adult repair. But if a module is biologically
load-bearing rather than a list artefact, some of its perturbable members should
move a growth outcome. That is worth looking at after A10's own question is settled,
never as a substitute for A5's or A11's own tests.

**The screen partly reaches questions that were blocked for lack of perturbation
data.** Amphiregulin and its receptors are targets here, which touches the source
and receiver questions A2 and A9; amphiregulin's role in repair extends well beyond
epithelium ([10.1016/j.immuni.2015.01.020](https://doi.org/10.1016/j.immuni.2015.01.020)).
An interleukin-1 ligand is a target, though the receptor is not, so A12's recipient
question is only partly touched. These are leads for those questions' own plans.

## One governance matter, stated plainly

The source paper for this deposit is roadmap paper 14
([10.1073/pnas.2606113123](https://doi.org/10.1073/pnas.2606113123)), and the
roadmap records the owner's reading of it as not started. The stop rule says later
papers are read when a pilot makes their branch relevant.

Therefore this work uses **the deposit, not the paper's claims**. No study note is
written, no finding of that paper is repeated here as established, and nothing in
this analysis depends on its conclusions. The precedent is explicit: when an earlier
session wrote a study note for a paper the owner had not read, the owner rejected it
and reserved the paper folder until after his own reading (DEVELOPMENT decision 21).

One consequence for interpretation: because the paper is unread here, which targets
the authors treated as controls and which as candidates is not established from the
publication. It must come from the deposited files or be left open.

## What a finished A10 would and would not license

A reproducible increment, holding out whole preparations, would support one specific
statement: in this assay, epithelial programme state carries information about
organoid growth beyond baseline size and batch, and possibly fibroblast response
carries more. That is a statement about a co-culture assay.

It would not establish that the same programmes govern repair in a living lung, and
it would not identify a mechanism. Its value is that it is the first question in this
register whose answer is anchored to something the experiment measured about tissue
building, rather than to another transcriptional score.
