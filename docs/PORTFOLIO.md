# Portfolio reading guide

For a short entry point, use the [one-page summary and four figures](PORTFOLIO_SUMMARY.md).

This project demonstrates independent computational work on public lung data:
atlas recovery, measurement correction, sample-aware analysis and an explicit
record of conclusions that did not survive scrutiny. Its contribution is the
reproducible analysis and reasoning; the public datasets belong to their source
studies. [Development and scientific ownership](../DEVELOPMENT.md) discloses
AI assistance and which decisions were delegated.

## Case study 1: recover an atlas and expose its annotation errors

**Question:** can an independent pipeline recover the deposited mouse cell
types without fitting to their labels?

**Work:** sample QC, doublet handling, clustering, marker annotation, comparison
with deposited labels, per-animal summaries and lineage-tracing checks.

**Evidence:** median cluster purity is 0.947 across labelled cells; three of 29
marker annotations disagree with the deposited types. Those disagreements
identify where state programmes confuse a type classifier. They remain in the
record rather than being silently renamed.

Read [Findings](../FINDINGS.md), the
[pipeline record](PIPELINE_AS_RUN.md) and claims C1/C8.

## Case study 2: correct an attractive ligand-source claim

**Question:** which expression-derived source and ligand-ranking observations
survive different annotations, molecule-depth budgets and resource definitions?

**Work:** donor-level source comparisons, explicit epithelial senders and
fibroblast targets, actual cell-count floors, receptor-complex-aware resource
inspection and saved ranking provenance.

**Why it matters:** a permissive “myeloid” gate included many T/NK cells, while
an epithelial-to-fibroblast scan admitted mural senders. A strong-looking
result can therefore be a population-definition problem. The corrected analysis
separates this from resource coverage and from unmeasured functional signalling.

Read the [corrected ligand analysis](../analysis/corrections/ligand/README.md)
for numerical results and figures; historical claims C37 and C80–C115 remain
traceable in the register.

## Case study 3: distinguish repeated programmes from validated mechanisms

**Question:** do pathway directions survive independent samples, proper
gene-set inference and plausible composition or design alternatives?

**Work:** pseudobulk reconstruction, official reference CAMERA, normalization
sensitivity, leave-one-sample-out influence and subtype composition checks.
W1's no-hit result survives the reference implementation. The previous G1
DNA-replication lead remains limited by phase/age/round confounding and by its
reassessment with sample-level uncertainty.

Read the [statistical correction](../analysis/corrections/statistics/README.md).
The connected [epithelial specificity project](../Thesis/epithelial_state_specificity/README.md)
adds explicit injury, development and genotype contrasts, source-defined gene
sets and an external-study eligibility assessment. It does not infer cellular
fate from a chromatin proxy.

## How to assess the work

Start from a question, follow its figure to the saved data, then inspect the
sample definition, code and limitations. The
[structured claim index](../analysis/claims/manifest.json) exposes which claims
have numeric bindings and which are only checked for narrative/path consistency.
[Reproduction instructions](../REPRODUCIBILITY.md) separate lightweight review
from scientific reruns requiring large public inputs.

The portfolio is strongest as evidence of research judgment and reproducible
implementation. No claim is strengthened by the intended lab or application.
Mechanistic relevance and new dataset selection must satisfy the design gates
in [Research questions](../RESEARCH_QUESTIONS.md).
