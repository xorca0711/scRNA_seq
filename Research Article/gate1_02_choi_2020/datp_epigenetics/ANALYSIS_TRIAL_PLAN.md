# Analysis trial plan: the DATP state in another genomic layer

The narrative for the trials in [`trials/`](trials/README.md), written after
they ran, and the proposals that follow from what they returned. The claims
table is in [`README.md`](README.md).

---

## 1. What the branch set out to ask, and why the obvious question was dropped

The first design asked whether transitional cells differ from AT2 cells in
chromatin accessibility. It was written, attacked from five directions, and
abandoned before it ran, because **that question cannot be false.** The labels
come from the transcriptome of the same nucleus whose chromatin is being read,
both modalities share that nucleus's quality, and both source papers defined the
state partly from its chromatin. A rule whose null is false by construction has
no state of the world in which it reads "no", which is the defect shape this
repository has now disclosed eight times.

The check that settled it was empirical, not rhetorical. A pseudo-contrast built
from reference cells only, depth-matched, containing no transitional cell at
all, returned 3,393 peaks over a 0.10 detection difference against a permuted
null mean of 2,617 at p = 0/200 in GSE310539, and the same in two further files,
with 100 per cent of crossing peaks pointing at the deeper group. A statistic
that returns certainty on data containing no biology is not a statistic.

The replacement question has two outcomes and both were possible:

> In the transitional state, is the AT2 identity programme **closed** at the
> chromatin level, or merely **silenced** at the RNA level?

Closed favours a committed, arrested state. Silenced-but-open favours a
transient, reversible one. That is not a manufactured dichotomy: Choi 2020 names
the state transient by construction, and Hassan and Chen write of their own data
that it is "compatible with two parallel states with only the AT1-like cells
transitioning to AT1 cells and KRT8/CLDN4+ cells being arrested".

---

## 2. Five refusals, each from a different guard

The trial refused five times before producing a number, and each refusal was a
different defect. They are listed because the repository's rule is that a
threshold is never moved after its result is seen, so every one of them is still
in the record with its first outcome.

| Trial | Refused because | Shape of the defect |
|---|---|---|
| M1 | The airway filter dropped any cell with a non-zero Scgb1a1, Scgb3a2, Foxj1 or Krt5 count. Scgb1a1 is ambient in 73 to 100 per cent of nuclei, so it removed 99.7 per cent of cells, and its severity tracked depth (Scgb3a2 detection 0.950 in the deepest well, 0.028 in the shallowest) | A threshold that measures library depth |
| M1b | Cldn4 detection AND Krt8 above the 99th percentile cannot exceed a few tenths of a per cent by arithmetic, against a reported 12 per cent | A threshold whose own arithmetic forbids the answer |
| M1c | The neonatal P9 wells held more of the labelled group (3.69, 8.07 per cent) than any injured well (5.09 per cent), so the whole-trial negative control fired | The marker set is not specific to the condition |
| M1d | The budget rule retained 60 per cent of every well while the drop gate forbade losing more than 20 per cent | Two rules that are not composable |
| M2 | The sensitivity control that licensed reading a null was itself unstable | A positive control tested once |

M1's refusal produced the branch's one validated result, below. M1c's refusal
produced one of the two answers the owner asked for.

---

## 3. What stands

### 3.1 GSE247130's deposited barcode order is wrong (Validated)

Rule R2 required the suffix-to-condition map, which is the GEO sample order and
is not deposited, to be checked on an independent axis before any name was used.
The check is that a knockout carries less of its own target gene.

Cebpa is higher in the suffix the GEO order calls the Cebpa knockout in all
three files, by 14.5x at P9, 10.9x at seven weeks and
4.5x in the infected file. A conditional knockout cannot carry more
of its own target than its control, so suffix 1 is the control and suffix 2 is
the mutant, the reverse of the deposited order. Two further axes agree: Cldn4 is
4.2-fold higher in suffix 2 of the infected file, which is the
expansion of transitional cells the source paper reports for the mutant, and
Sox9 is 7.1x, 2.8x and 2.5x higher in
suffix 2 across the three files, much the strongest in the neonatal one, which
is the stage-specificity that paper reports.

**Those figures were wrong when this branch first wrote them, and the Sox9 axis
was not measured at all.** Trial M1's frozen corroboration panel was Cldn4, Fos
and Cebpa; its write loop iterated over a fourth gene behind a guard that
dropped it silently, and the Sox9 figure existed only as a hard-coded string.
The register said Cebpa was 11 to 14 times higher in all three files when the
infected file is 4.5 times, and it said Cldn4 was 4.5-fold higher in that file
when the logged table gives 4.2. All four are recomputed and logged in
[trial M4](trials/m4_the_corroboration_the_register_claimed/m4_summary.md). The
conclusion about the suffix map is unaffected; the account of how well it was
shown was not. Rows C116 and C151.

GSE310539 is corroborated as deposited: Fos is 3.7-fold
higher in the wildtype wells than in the mutant ones, and Cldn4 rises
10.5-fold with infection in the wildtype pair and 5.3-fold in
the mutant pair. The observation that Fosb and Junb do not discriminate the
genotypes, so that a check built on those two would have passed the wrong
answer, is recorded in trial M1's own output and is not recomputed here.

### 3.2 The AT2 identity programme is silenced in RNA (Descriptive only)

In the two injured wells that reached the statistic, the CLDN4-positive
KRT8-positive alveolar group loses the AT2 identity programme by **12.6 to 17.1
detection points** against a sham band of twenty random splits of reference
cells, at every one of four downsampling seeds. The AT1 arm does not move in any
seed of any well. Per gene the loss is carried by at least five genes: Etv5
falls 0.227 and 0.239, Napsa 0.272 and 0.167, Slc34a2 0.271 and 0.189, Abca3
0.182 and 0.212, Lamp3 0.160 and 0.115.

Descriptive only, and not negotiable: each well is one library pooling two mice.

### 3.3 The marker set is not specific to injury (Descriptive only)

This answers the owner's first question, and the answer is no. Krt8 and Cldn4
are expressed across immature postnatal alveolar epithelium, so a
CLDN4-positive KRT8-positive transcript call cannot separate neonatal
developmental immaturity from injury-induced transition. Hassan and Chen's own
argument predicts this, since their case is that neonatal AT2 cells are plastic
and that Cebpa deletion returns mature cells toward that state, but they never
state it as a limitation of the marker.

---

## 4. What was withdrawn, and why it is recorded rather than deleted

M1e read **silenced but not closed** in both test wells, on the strength of its
rule R9: the AT2 chromatin null may be read as retention only behind a
transitional-marker arm that fired. M2 shows that licence was not earned.

* Across four seeds, the transitional chromatin arm clears in **five of eight**
  seed-and-well combinations, not eight of eight.
* It fails leave-one-out on two of its six genes in each injured well, and per
  gene it is carried by Sfn and Ndrg1. In wildtype_SeV, Sfn alone returns 0.0357
  over seven peaks against an arm mean of 0.0097.
* The offset correction was not inert, it was **absent**. M1e subtracted the
  genome-wide **median** distal difference, which came out exactly 0.0 in every
  well because more than half of all distal peaks are detected in neither group
  at this budget. That is not the same as there being nothing to remove: the
  mean global distal difference is +0.00080, +0.00034 and +0.00201 in the three
  wells, positive in all three, across two deposits, two peak atlases and two
  budgets. A consistently signed nuisance added to two different true values is
  exactly what produces arm values that disagree in sign, so the sign
  disagreement between wells damns the analysis rather than rescuing it.
  Recomputed against the mean, the AT2 arm in SeV_Cebpa_mutant moves from z =
  -1.88 and not clearing, to z = -3.23 and clearing **in the closing direction**.
* Three sham standard deviations, the smallest value the rule would have called
  a clearance, is 0.0043 to 0.0120 in detection fraction. The effects present are
  0.004 to 0.012. The design works at its own floor.

An adversarial review reached the same place by a different route and sharpened
it. The defence that the two arms have comparable sham standard deviations, and
therefore comparable power, is a category error: the AT2 arm carries 129 distal
peaks over 9 genes and the transitional arm 50 over 6, so equal variance implies
the AT2 peaks are roughly four times more accessible at baseline. The two arms
have equal variance for opposite reasons.

A second reviewer found the deeper version of the same problem: **the AT2 arm is
bounded and the positive control is not.** The reference-side baseline of the
AT2 arm is about 0.0185 per peak in GSE310539 and 0.0380 in GSE247130, so the
arm cannot fall below those values even if every distal element at all nine AT2
loci shut completely, whereas the transitional arm can be driven arbitrarily
high by de novo opening. The defence that an effect the size of the transitional
arm's would have cleared therefore describes a counterfactual requiring 58 per
cent of all distal accessibility at those loci to vanish in fourteen days. An
unbounded statistic was being used to certify the power of a bounded one.

The honest statement of what was measured is that **no change was detected at a
bar corresponding to a 45 per cent loss of distal accessibility at the AT2 loci
in GSE310539, and a 12 per cent loss in GSE247130.** The second is a bound worth
having and the first excludes almost nothing. Converted to the scale the biology
is on, only one well supports a bound worth having: in the infected Cebpa mutant,
AT2 distal accessibility differs from reference by about -4 per cent of
reference accessibility with a 95 per cent interval of roughly -10 to +2 per
cent, whereas in the infected wildtype the interval is about -15 to +27 per cent
and excludes nothing.

A third reviewer pushed the offset further and found the direction. Under a
background matched on baseline accessibility the AT2 arm in the Cebpa-mutant
well reaches z = -2.89, and under a multiplicative correction z = -3.55, which
would clear as **closed**; under that same correction the transitional positive
control falls to z = +2.53, below the frozen floor, so rule R9 would have
refused the well. Across estimators that well reads either not computable or
closed, and never what M1e reported. Whatever is in these data, it is not
evidence that AT2 chromatin stays open.

A fourth checked independence and found an error of mine. Jichao Chen is a
contributor on both GEO series and the contact laboratory for GSE247130, so the
two deposits are one laboratory with two first authors, sharing the lineage
tools, the injury model, the sorting, the kit, the aligner, and on this side the
code, seeds, annotation and budget rule. Agreement between them is consistency,
not replication, and the branch documents have been corrected.

**Sentences this folder may not contain**, recorded so they cannot creep back
in: that AT2 chromatin "does not move", is "unchanged", "retained", "preserved"
or "stable"; that the programme is "silenced in RNA but not closed in
chromatin"; that any result here "replicated in two independent deposits", since the deposits share a laboratory;
or that the AT2 arm "had the power to see a transitional-sized effect and did
not".

---

## 4b. Trial M3, and the direction that will not go away

A five-lens adversarial review of M2 found three defects in M2 itself. This
repository did not find them, and the record says so.

**The two deposits were never running the same instrument.** Rule R12 removed
Cebpa from the AT2 arm where it is genetically deleted and kept it everywhere
else. Cebpa carries three distal peaks but the largest positive per-gene value
in wildtype_SeV, and it is the only reason that well's AT2 arm was positive at
all. M2's own leave-one-out had already said so and nobody read it that way.

**The sham band is the wrong null for a gene-set arm.** It permutes cells, so it
asks whether this split is special; it cannot ask whether these genes are
special. Three hundred gene sets matched gene by gene on distal peak count and
evaluated on the real split have a mean of +0.00069 and +0.00205, not zero.
That is also why the genome-wide median offset came out exactly 0.0: a fifth to
a third of arm peaks are detected in no cell of either group.

**The equivalence statistic erred unsafely.** Three sham standard deviations is
the 50-per-cent-power detection floor, not a bound the data support.

With one eight-gene instrument, a gene-set null, five hundred shams and
confidence intervals on the relative scale, **the AT2 identity arm is negative
in all three computable wells and sits at percentile 0.000, 0.007 and 0.000 of
300 matched random gene sets**, at -37.8, -18.3 and -15.5 per cent of reference
accessibility.

**And it is still Not established, because the two wells fail different gates.**
wildtype_SeV has a positive control that clears (z 3.46 and 3.05) and an AT2 arm
that does not (z -2.16 and -2.73). SeV_Cebpa_mutant has an AT2 arm that clears
both nulls at z -3.41 and -3.51 in the closing direction, behind a positive
control that does not clear (z 1.96 and 1.68), so rule R9 refuses it. Each well
holds half of what a reading needs.

The reason the mutant well's positive control now fails is itself a result:
re-centred on the gene-set null, its transitional arm falls from clearing to
z = 1.96, so **most of what M1e counted as signal there was the uncorrected
global shift.** Only wildtype_SeV has a transitional signal that survives proper
centring, at +61.8 per cent of reference.

The review's predictions were pre-registered under rule R22 and came out mixed:
the gene-set null means were predicted almost exactly (+0.0009 and +0.0021
against +0.00069 and +0.00205), the AT2 direction and the percentile-zero result
were confirmed and stronger than predicted, and the transitional arm came out
weaker than predicted, which is why a well the review expected to pass is now
refused.

---

## 5. What to do next, in this layer

**Proposal A, and it is the only one here with a testable contrast.**
GSE309751, the bulk ATAC-seq accompanying Lynch et al., carries **genuine
biological replicates**: 14-day mock PBS n = 2, 14 dpi n = 3, 49-day mock n = 2,
49 dpi n = 2, Kras control n = 2, Kras induced n = 3, all on sorted
lineage-labelled AT2 cells. This is the first chromatin contrast anywhere in
this project with within-group replication at the animal level, so the mouse
can be the unit and a real test is possible.

Two questions it can carry that nothing else here can:

1. Do AT2 identity loci lose accessibility 14 days after Sendai virus, tested
   across mice with a threshold frozen from the mock arm?
2. **Do they recover by 49 days?** That is the reversibility question, and
   reversibility is what "transient" in damage-associated *transient* progenitor
   actually asserts. Lynch report that AP-1 peaks open and reclose; whether the
   AT2 programme's own loci recover is not reported.

Two limits to state in the pre-registration rather than discover: bulk on sorted
AT2 cells cannot resolve the transitional substate, so a state-specific effect
is diluted roughly eightfold by the cells that did not transition; and the
deposit gives per-sample MACS2 broadPeak calls rather than a count matrix, so
peak presence across mice is testable while a properly normalised differential
accessibility test needs the reads from SRA. n = 2 per group at 49 days is thin
and the write-up must say so.

**Proposal B.** Download the `atac_fragments.tsv.gz` files this branch skipped
(2.6 GB for GSE310539, plus GSE247130's three). Fragments allow peaks to be
re-called on the labelled cells themselves instead of inherited from a peak set
ascertained on the majority population, which is the artefact M2 could bound but
not remove, plus per-cell ATAC quality control, TSS enrichment and footprinting.

**Proposal C.** GSE247271, the CEBPA and NKX2-1 ChIP-seq from Hassan and Chen,
is a different genomic layer again: protein-DNA occupancy. It can ask whether
the AT2 identity loci examined here are CEBPA-bound, which would connect Choi's
definition of primed AT2 by loss of Etv5, Abca3 and Cebpa to a physical
mechanism. Two replicates per antibody, intersected; descriptive.

---

**Proposal D, and it is a lead rather than a proposal because this branch has
no artefact behind it.** An adversarial recomputation reports that the labelled
cells detect 12 to 13 per cent more distinct distal peaks than reference cells
at identical in-peak fragment depth. If that holds it is a global statement
about the transitional state's chromatin, broader accessibility rather than
different accessibility, and it would be more interesting than anything this
branch tried to measure gene by gene. It is recorded here as Not established,
because no trial in this folder computed it, and the first job of a trial M3
would be to compute it with a sham band and a negative control.

## 6. Other genomic layers, with feasibility stated plainly

The owner asked for proposals in further layers after each analysis. These are
sorted by whether they can actually be done, not by how interesting they sound.

### 6.1 Transcript isoforms: NOT FEASIBLE on any data in this branch

Every matrix in this folder is 10x three-prime single-nucleus counting. It reads
the last few hundred bases of a transcript and cannot distinguish isoforms,
alternative first exons, retained introns or alternative polyadenylation beyond
the three-prime end. This is a limitation to state, not a plan to execute, and
it is the same limitation already recorded for Tnc variable FNIII domain usage
in the Cardoso branch.

It would take new data of a different kind: long-read single-cell (PacBio Kinnex
or ONT) on injured mouse lung. No such dataset is held here and none was found
for the transitional state. If one is wanted, the search is for long-read
single-cell lung injury data in GEO and SRA, and the honest prior is that it
probably does not exist yet for this model.

### 6.2 Protein: NOT FEASIBLE here, and it is the layer the biology most needs

Both source papers define their states by immunostaining, not by transcript,
and this branch has already measured the gap that opens as a result: at the
depth available, a KRT8-high call of the kind Lynch makes by antibody cannot be
made from counts, because Krt8 has a median count of 0 or 1 per nucleus.

Two protein questions are already outstanding in this repository and neither can
be answered by any deposit held here. Trial C13 refuted the transcriptional
explanation for the sort contaminant and left surface EpCAM dimming as a
post-translational hypothesis, which needs protein. And ADAM17, the shared
sheddase in that account, is controlled post-translationally, so its mRNA is a
weak proxy for its activity.

What would serve: CITE-seq or Abseq with a surface panel on injured mouse lung,
or targeted flow cytometry for surface EpCAM on sorted transitional cells. Both
are experiments, not reanalyses.

### 6.3 Spatial: FEASIBLE IN PRINCIPLE, and it is the highest-value missing layer

Lynch make a spatial claim this branch cannot touch: the CLDN4-positive substate
forms small clusters restricted to "AT2-less" regions undergoing in-situ repair,
while KRT8-high cells in de novo compensatory-growth regions are
CLDN4-negative. If that holds, then the transcriptional label used throughout
this branch is pooling two spatially distinct populations, which would explain
part of the heterogeneity M2 measured.

Spatial data also carries a standing obligation from elsewhere in this project:
a ligand-receptor hit without proximity evidence is a hypothesis, not an
interaction, and the Cardoso branch has accumulated several. A spatial dataset
of injured mouse lung would let the AREG-EGFR axis be asked in the one way that
matters, which is whether the cells are near each other.

The concrete next step is a GEO survey for Visium, Xenium, MERFISH or CosMx on
injured or fibrotic mouse lung, run the way the GSE310539 survey was run, with
the usual question asked first: does any of it carry more than one animal per
condition.

### 6.4 DNA methylation and histone marks: NOT AVAILABLE, and the authors say so

Hassan and Chen name this gap themselves, writing that the duration of chromatin
closure might lead to less reversible changes in histone modifications, DNA
methylation or higher-order chromatin structure. They present no such data, and
none exists for this state in either deposit. Recorded as Not established rather
than proposed, because there is nothing to reanalyse.
