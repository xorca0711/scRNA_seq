# Branch note: the DATP state in another genomic layer

A branch of roadmap paper 2 (Choi, Lee et al. 2020), opened 2026-09-20 at the
owner's direction. It asks what the transitional state looks like in a layer
other than messenger RNA, and it is on other people's data for a reason that is
itself a finding.

**Current scope after the 2026-09-22 audit.** The next analysis is the
[epithelial-state specificity project](../../epithelial_state_specificity/README.md),
which joins the neonatal two-marker caveat to a source-backed RNA module
comparison. Reduced AT2 accessibility alone cannot distinguish productive AT1
differentiation from arrest. The historical trials below retain their outputs;
the earlier closed-versus-reversible fate interpretation is not retained.

**Why this branch is not on the Choi deposit.** Choi 2020 supports the
epigenetic half of its Il1r1-positive AT2 claim with ATAC-seq (Figures 5 and 6),
but GSE144598 deposits two bigwig coverage tracks and nothing else: no peaks, no
reads, one pooled sample per group. Trial D0 recorded that accession as
unusable on 2026-09-13 and it remains unusable. The epigenetic claim of the
parent paper cannot be re-derived from the parent paper.

**What is used instead.** Two 10x multiome deposits, in which the same nucleus
gives both a transcriptome and a chromatin accessibility profile.

| Deposit | Paper | Content |
|---|---|---|
| [GSE310539](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310539) | Lynch, Noun, Yang, Zhou, Chen, Evans, Kadara, Chen. *AP-1 mediated chromatin changes govern alveolar type 2 cell transition in lung injury-repair.* Am J Respir Cell Mol Biol 2026. DOI [10.1093/ajrcmb/aanag157](https://doi.org/10.1093/ajrcmb/aanag157), preprint [10.1101/2025.10.25.684549](https://doi.org/10.1101/2025.10.25.684549) | 4 libraries, 39,849 cells, 32,287 genes and 178,178 peaks. Sorted CDH1-positive epithelium, 14 days after Sendai virus or PBS, wildtype and Fos/Fosb/Junb triple mutant |
| [GSE247130](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247130) | Hassan, Chen. *CEBPA restricts alveolar type 2 cell plasticity during development and injury-repair.* Nat Commun 2024;15:4148. DOI [10.1038/s41467-024-48632-3](https://doi.org/10.1038/s41467-024-48632-3), PMID [38755149](https://pubmed.ncbi.nlm.nih.gov/38755149/) | 6 libraries in 3 files, 64,294 cells. Sun1GFP AT2-lineage sorted, Cebpa mutant against control, at P9, at 7 weeks, and 14 days after Sendai virus |

---

## What each source paper does and does not already answer

This matters, because a reanalysis that re-derives a published result and calls
it confirmation is worth nothing. **Both papers are chromatin papers.** The gap
is narrow, specific, and different in the two.

**Lynch et al. 2026** is titled for its chromatin work. Bulk ATAC-seq at
baseline, 14 and 49 days after infection (a separate accession, GSE309751);
chromVAR motif accessibility per substate; AP-1 motif accessibility spiking in
the CLDN4-positive T2 substate; three AP-1-motif peaks at the Cldn4 locus;
and a chromatin module that collapses in the AP-1 mutant. What they do NOT do
is anchor accessibility to genes: they note in passing that their transitional
substates lose AT2 chromatin features *less dramatically* than they lose AT2
RNA, never quantify it, and never state it for the CLDN4-positive substate.
**Anything this branch finds on GSE310539 is registered as a quantification of
their observation, not as independent confirmation.**

**Hassan and Chen 2024** is also a chromatin paper: a five-timepoint
developmental scATAC series (GSE264098), 30,410 differential peaks in four motif
classes, 10,621 differential peaks on neonatal Cebpa deletion, and CEBPA and
NKX2-1 ChIP-seq (GSE247271). But **none of it touches their transitional
cells.** Their Sendai figure carries no ATAC panel, and consistently their two
Sendai ATAC samples deposit no peak files where their other four deposit two
each. For that deposit the question is genuinely unoccupied.

**The two deposits are not independent, and an earlier draft of this note said
they were.** Jichao Chen is a contributor on both GEO series and the contact
laboratory for GSE247130, so these are two first authors in one laboratory, not
two laboratories. They share the SftpcCreER and RosaSun1GFP lineage tools, the
Sendai virus model, E-cadherin-positive sorting, the 10x Multiome kit and
cellranger-arc on mm10; and on this side they share the analysis code, the
seeds, the vendor peak annotation and the budget rule. What genuinely differs is
the first author, the mouse cohort, the institution, the year and the genotype.
Agreement between them is a consistency check, not a replication, and nothing in
this folder may be described as replicated across independent data.

**Neither paper calls its state DATP on its own authority.** Lynch uses the word
once, in the introduction, listing DATP, PATS and ADI as prior names with
"overlapping yet distinct" signatures, and computes no DATP score. Hassan and
Chen overlay Choi's 89-gene DATP score on their cluster once, in a figure
legend, without argument. Treating either state as the DATP state is an
assumption, and this folder does not make it.

---

## The ceiling both deposits share, and it is the fifth and sixth time

Trial M0 applied the replicate rule this repository inherited from trial D0: a
contrast has within-group replication only if at least two libraries share every
experimental variable except the one being contrasted.

**Not one contrast between conditions in either deposit carries within-group
replication.** Ten libraries, 104,143 cells, one library per condition, and each
library pools two mice before loading so the mouse is not recoverable from
barcodes. Nothing about Sendai infection, AP-1 deletion, Cebpa deletion or
developmental stage is testable here, and nothing in this folder attempts it.

What is admissible is a comparison between cell states **inside one library**,
where animal, batch, peak set and sequencing run are held fixed by construction.
The precedent is trial D3. Everything in this folder is registered Descriptive
only.

---

## Claims, organised by what they license

Grouped by the decision each one governs rather than by the trial that produced
it, because this table is for planning the next analysis. Every number sits in
the trial summary the row links to, and the repository-wide row identifier is
given so the register and this page cannot drift apart.

**Read this table as four questions.** What do I have to accept about the data
(A)? What is the state, at the level the data can speak (B)? What does the
chromatin say (C)? And what may I no longer do (D)?

---

### A. Constraints any future use of these deposits inherits

| # | Established | Status | What follows for planning |
|---|---|---|---|
| A1 | GSE247130's deposited barcode order is inverted: suffix 1 is the control, suffix 2 the Cebpa mutant ([C116](../../../CLAIMS.md)) | **Validated** | Invert the map before touching that deposit. Generally: never attach a condition name to a deposited sample order without corroborating it on a knockout's own target gene first. This cost nothing to check and would have inverted every genotype downstream |
| A2 | Neither deposit carries within-group replication for any between-condition contrast ([C117](../../../CLAIMS.md)) | **Not establishable** | No genotype, treatment or stage claim can ever come from these files. Only within-library state comparisons are admissible and everything is Descriptive only. This is the fifth and sixth deposit in the project with that ceiling, so assume it of the next one until shown otherwise |
| A3 | Both source papers are chromatin papers. The untouched gap is accessibility of the transitional state in GSE247130, and gene-anchored rather than motif-level accessibility in GSE310539 | **Descriptive only** | Any chromatin result on GSE310539 must be registered as a re-derivation of Lynch et al., not as confirmation. Only GSE247130 can yield something new in this layer |
| A4 | The two deposits share a laboratory, lineage tools, injury model, kit, aligner and this analysis code | **Descriptive only** | Agreement between them is consistency, not replication. Nothing in this branch may be called replicated, and a genuinely independent deposit is worth more than a third well here |

### B. What the transitional state is, at the level transcripts can speak

| # | Established | Status | What follows for planning |
|---|---|---|---|
| B1 | The CLDN4-positive KRT8-positive group loses the AT2 identity programme in RNA: 12.6 to 17.1 detection points, five or more genes, with the AT1 arm flat in the same cells ([C118](../../../CLAIMS.md)) | **Descriptive only** | The state's defining transcriptional move is **loss of AT2 identity, not gain of AT1**. Design the next trial around the AT2 programme as the moving part; an AT1-anchored design would have measured nothing |
| B2 | That transcript call cannot separate neonatal developmental immaturity from injury-induced transition ([C119](../../../CLAIMS.md)) | **Descriptive only** | **This answers the owner's Q1 and the answer is no.** Any future state definition must either add a marker that development does not use, or hold developmental stage fixed by design. Do not treat Krt8 and Cldn4 as damage-associated on their own |
| B3 | The source papers define this state by immunostaining (KRT8-high by antibody); three-prime single-nucleus counting cannot reproduce that call ([C124](../../../CLAIMS.md)) | **Descriptive only** | Stop trying to rebuild an antibody-defined state from 3' counts. Name the transcript-level group as its own object and say what it is not. The same mismatch is the reason the Axin2 proposal needs care: Nabhan defines Axin2-positive by lineage reporter, not transcript |

### C. What the chromatin says, and what it does not

| # | Established | Status | What follows for planning |
|---|---|---|---|
| C1 | The background-centered AT2 distal contrast is negative in all three wells reaching the statistic: percentile 0.000, 0.007 and 0.000 of 300 matched gene sets, at -37.8, -18.3 and -15.5 per cent of reference accessibility ([C131](../../../CLAIMS.md)); raw reductions are -31.8, -16.0 and -8.5 per cent | **Exploratory** | No well passes both gates. The uninjured comparator is Cebpa mutant, which can show genuine genetically induced identity change. Accessibility loss alone does not favour arrest over productive differentiation. The reference-cell sham interval is not an animal-level confidence interval or a biological power estimate |
| C2 | Whether the AT2 identity programme's chromatin actually closes ([C121](../../../CLAIMS.md), [C133](../../../CLAIMS.md)) | **Not established** | No well satisfies both gates: the well whose AT2 arm clears has a positive control that does not fire, and vice versa. Needs more labelled cells, or peaks re-called on the labelled cells from the fragments files, before it is worth re-asking |
| C3 | "Silenced but not closed" ([C120](../../../CLAIMS.md)) | **Retracted, superseded** | Never write this sentence again. The positive control that licensed it was unstable, and its direction was an artefact of one gene carrying three peaks |
| C4 | AP-1 dependence of anything above ([C128](../../../CLAIMS.md)) | **Not established** | The infected AP-1 mutant well held 44 labelled cells against a floor of 50. Six cells short of the branch's most interesting genotype, so a design that recovers a few hundred more labelled cells buys that arm outright |

### D. Methodological constraints that now bind every trial in this repository

These are the transferable results. Each was learned here at the cost of a
refused trial, and each forbids a design that would otherwise look reasonable.

| # | Refuted | What follows for planning |
|---|---|---|
| D1 | That "do transitional cells differ from AT2 cells in accessibility" is a question ([C127](../../../CLAIMS.md)) | It cannot be false: labels and chromatin come from one nucleus. It returned p = 0/200 on a depth-matched contrast containing no transitional cell. **Before freezing any design, construct the contrast that contains no biology and check the design returns nothing on it** |
| D2 | That a cell-permutation sham is a sufficient null for a gene-set claim ([C130](../../../CLAIMS.md)) | Permuting cells asks whether the SPLIT is special; it cannot ask whether the GENES are. Every gene-set arm needs a second null over gene sets matched on peak or transcript count |
| D3 | That a positive control tested once licenses reading a null ([C126](../../../CLAIMS.md)) | Test it under resampling and leave-one-out first. A sensitivity control that is itself unstable licenses nothing |
| D4 | That presence of an abundant secreted transcript is identity ([C123](../../../CLAIMS.md)) | In single-nucleus data it is ambient, and a presence rule built on one is a depth filter under another name. Use within-cell ratios |
| D5 | That three sham standard deviations is a bound the data support ([C130](../../../CLAIMS.md), and see M2) | It is the 50-per-cent-power detection floor and it errs unsafely. Report a confidence interval on a relative scale the reader can interpret |
| D6 | That two separately written rules are composable ([C125](../../../CLAIMS.md)) | Check the arithmetic of every rule pair before freezing. A budget retaining 60 per cent and a gate forbidding 20 per cent loss cannot both hold |
| D7 | That one instrument was being used across deposits ([C129](../../../CLAIMS.md)) | A gene dropped in one arm and kept in another makes two arms that were never comparable, and here that one gene carried the sign. Freeze the gene list once, for every well |

---

## The two questions the owner asked this branch

**"Does the transition state show a character specific to regeneration or
disease?"** At the transcript level, for the markers that define it, no. The two
neonatal P9 wells carry more CLDN4-positive KRT8-positive cells (3.69 and 8.07
per cent) than any injured well (at most 5.09 per cent). Krt8 and Cldn4 are
expressed across immature postnatal alveolar epithelium, so the marker set is
shared with normal development and is not by itself damage-associated. Hassan
and Chen's own argument predicts this and never states it as a limitation of
the marker.

**"Is there an epigenetic character to the transition state?"** Not established
from these matrices. The binding constraints are measured, not guessed: the
labelled groups hold 64 and 320 cells; the smallest effect the design would
have called real is 0.004 to 0.012 in peak detection fraction, and the effects
present are of that same size; and the sensitivity control that would have
licensed reading a null is itself unstable across resampling seeds.

The one bounded statement the data support is in the infected Cebpa-mutant well
alone: AT2 distal accessibility there differs from reference by about -4 per
cent of reference accessibility, with a 95 per cent interval of roughly -10 to
+2 per cent. In the infected wildtype well the same interval is about -15 to
+27 per cent, which excludes nothing worth excluding. Put on the bar the trial
actually cleared, it detected no change at a level corresponding to a 45 per
cent loss of distal accessibility at the AT2 loci in GSE310539 and a 12 per cent
loss in GSE247130.

Trial M3 then rebuilt the instrument and the null, at an adversarial review's
direction, and the direction stopped being a lean. With Cebpa dropped from the
AT2 arm in every well so that both deposits run one eight-gene instrument, and
with the offset taken from 300 gene sets matched on distal peak count rather
than from a median that was zero by construction, **the AT2 arm is negative in
all three computable wells and sits at percentile 0.000, 0.007 and 0.000 of
those 300 matched sets.** It is still Not established, because no well
satisfies both gates: the well whose AT2 arm clears both nulls has a positive
control that does not fire, and the well whose positive control fires has an
AT2 arm short of the floor. But the question now has a direction and a bound.

Where the data lean at all, they lean the opposite way from the withdrawn
reading. Under every offset estimator other than the degenerate median, the AT2
arm in the Cebpa-mutant well moves further toward CLOSING, reaching z = -3.55
under a multiplicative correction, which would have read as closed. Under that
same correction the positive control falls below its floor and the well would
have been refused. So the well reads either not computable or closed, never
retained. **Nothing in this folder supports the idea that AT2 chromatin stays
open in the transitional state.**

---

## What would settle it, and it is not another threshold

Both deposits carry `atac_fragments.tsv.gz` (2.6 GB for GSE310539) which this
branch did not download. Fragments permit peaks to be re-called on the labelled
cells themselves rather than inherited from a peak set ascertained on the
majority population, plus per-cell ATAC quality control, TSS enrichment and
footprinting. That is the next instrument.

See [`ANALYSIS_TRIAL_PLAN.md`](ANALYSIS_TRIAL_PLAN.md) for the trial-by-trial
narrative and for the proposals in other genomic layers.
