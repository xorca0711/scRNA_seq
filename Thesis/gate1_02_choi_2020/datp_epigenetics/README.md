# Branch note: the DATP state in another genomic layer

A branch of roadmap paper 2 (Choi, Lee et al. 2020), opened 2026-09-20 at the
owner's direction. It asks what the transitional state looks like in a layer
other than messenger RNA, and it is on other people's data for a reason that is
itself a finding.

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

## Claims table

| Claim | Status | Where |
|---|---|---|
| GSE247130's barcode-suffix map is inverted relative to the GEO sample order: suffix 1 is the control and suffix 2 is the Cebpa mutant | **Validated** | [M1](trials/m1_closed_or_merely_silenced/m1_summary.md) |
| The CLDN4-positive KRT8-positive alveolar group loses the AT2 identity programme in RNA by 12.6 to 17.1 detection points, across five or more genes, with the AT1 programme flat in the same cells, in both deposits and at four resampling seeds. The deposits share a laboratory, so this is consistency, not replication | **Descriptive only** | [M1e](trials/m1e_per_well_budget/m1e_summary.md), [M2](trials/m2_robustness_of_the_m1e_reading/m2_summary.md) |
| A CLDN4-positive KRT8-positive transcript call cannot separate neonatal developmental immaturity from injury-induced transition | **Descriptive only** | [M1c](trials/m1c_label_at_the_depth_available/m1c_summary.md) |
| In the CLDN4-positive group the AT2 identity programme is silenced in RNA but its chromatin stays open | **Retracted, superseded** by M2 | [M1e](trials/m1e_per_well_budget/m1e_summary.md) |
| The AT2 identity programme's distal chromatin is lower in the labelled group than in 300 matched random gene sets, in all three computable wells, at percentile 0.000, 0.007 and 0.000, and by -37.8, -18.3 and -15.5 per cent of reference accessibility | **Exploratory** | [M3](trials/m3_one_instrument_and_the_right_null/m3_summary.md) |
| Whether the AT2 identity programme's chromatin closes in the transitional state. M3 gives the question a direction and a bound; no well satisfies both gates | **Not established** | [M2](trials/m2_robustness_of_the_m1e_reading/m2_summary.md), [M3](trials/m3_one_instrument_and_the_right_null/m3_summary.md) |
| Anything about AP-1 dependence of the above | **Not established**; the infected AP-1 mutant well held 44 labelled cells, below the floor | [M1e](trials/m1e_per_well_budget/m1e_summary.md) |

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
