# Branch note: Axin2-positive and Il1r1-positive AT2 cells, and whether they are the same cells

A branch of roadmap paper 2 (Choi, Lee et al. 2020), opened 2026-09-20 at the
owner's direction. It assesses a question the owner proposed and does not run a
trial on it, because the assessment concludes that the definitive version cannot
be run on any data that exists.

**The question is Choi 2020's own, and this is its closing Discussion sentence
verbatim:**

> Recently, Axin2+AT2 cells have been identified as a distinct subset of AT2
> cells (Nabhan et al., 2018; Zacharias et al., 2018). Related to the potential
> role of interconnectivity between IL-1beta and Wnt signaling in fate decision
> of AT2 cells, comparison of Il1r1+AT2 and Axin2+AT2 cells will be helpful to
> understand their relationships during alveolar regeneration.

Six years on, nobody has made that comparison. The gene symbols above were
recovered from the Europe PMC full-text XML for PMC7487779; both the PMC web
rendering and the PMC full-text API strip the italic tags and lose them.

---

## 1. The headline

| Question | Status |
|---|---|
| Are Axin2-positive AT2 cells a distinct subset? | **Not established.** The field's own data lean toward a state |
| Are Il1r1-positive AT2 cells a distinct subset? | **Not established**, for the same reasons, and the authors say so |
| Do the two populations overlap? | **Not established.** One hedged supplementary panel says they do |
| Can the definitive comparison be made from public data? | **No.** No deposit carries both readouts in the same cells |
| Can anything useful be done on data already here? | **Yes**, and section 5 costs them. One of the three was withdrawn as out of focus and replaced |

---

## 2. Why the premise is shakier than the question assumes

Both populations are defined by a tamoxifen-inducible lineage reporter. Neither
is a transcript call, an antibody stain, or a surface marker, and that governs
everything below.

| | Axin2-positive | Il1r1-positive |
|---|---|---|
| Allele | Axin2-CreERT2 knock-in, crossed to a Cre-dependent reporter | Il1r1-P2A-eGFP-IRES-CreERT2 knock-in, crossed to Rosa26-tdTomato |
| Fraction at homeostasis | **1 per cent** (Nabhan) against **20 per cent** (Zacharias), and about 29 per cent quoted since | 15 per cent |
| Fraction after injury | 73 per cent (diphtheria toxin) | 60 per cent (bleomycin, day 14) |
| Labels other cell types | mesenchyme in the Wnt-responsive work | airway ciliated cells and some mesenchyme, so AT2 identity needs a co-marker |

**The Axin2 fraction is unsettled by a factor of twenty to thirty**, between two
different knock-in alleles at different tamoxifen doses. Han's JCI review states
it plainly: what controls the range between distinct Cre lines remains unknown.
Eight years on, nobody has explained it, so the fraction is currently a property
of the allele and the dose rather than of the lung. Note also that Axin2 is a
negative regulator of Wnt, so an Axin2 knock-in allele is itself a heterozygous
Wnt sensitiser.

**Nabhan's stability argument cannot settle subset against state, structurally.**
Once a cell recombines it is labelled permanently, so a cell that is Axin2-positive
for a week, recombines, then switches Axin2 off is still scored as a permanent
member of a stable subpopulation. The experiment tests whether repeated tamoxifen
recruits more cells; it cannot test whether cells enter and leave the state.

**And the same paper carries three lines of evidence for a state.** After
diphtheria-toxin ablation 73 per cent of AT2 cells express Axin2, so bulk AT2
cells turn the marker on. AT2 cells within 15 micrometres of a Wnt5a-expressing
fibroblast express Axin2 far more often than those further away, so marker
status tracks position relative to a ligand source. Adding recombinant Wnt3a or
Wnt5a switches Axin2 on in ordinary cultured AT2 cells and Dkk3 switches it off.
The Morrisey laboratory's own developmental data (Frank 2016) asked this
directly with a three-week chase and concluded that some cells lose
Wnt-responsiveness while others gain it.

**Il1r1 is in exactly the same shape**, and Choi's Limitations section says so:
the rise from 15 to 60 per cent is compatible with expansion of a pre-existing
pool or with fresh Cre firing in bulk AT2 during repair, and washout periods
longer than 16 days would be needed to separate them. The lineage-tracing
quantifications also plot one dot per histological section rather than per
mouse, so by this repository's rules those figures describe sections.

**Neither marker is established as a discrete subset. Both behave as
injury-inducible signalling states layered on the AT2 compartment.**

---

## 3. What is actually known about the overlap

**There is exactly one direct measurement in the primary literature**, and it is
in a paper already on this roadmap: England et al. 2025 (roadmap paper 12,
doi:10.1016/j.stem.2025.01.011), whose deposit GSE247505 is already held here.
Figure S2L to S2N reports that Axin2 is specifically upregulated in the
Il1r1-lineage ZsGreen-positive fraction, and the authors write that this
suggests Axin2-positive cells are enriched in the Il1r1-positive subset.

So the only direct evidence says the two are **positively correlated, not
opposed**, which is the opposite of what phenotype-level reasoning predicts. It
is also a supplementary-figure expression measurement on a bulk-sorted fraction,
at homeostasis only, hedged with "suggesting", with no contingency table: what
fraction of Axin2-positive cells are Il1r1-positive, the converse, and whether
it survives injury are all unreported.

**The apparent opposition is an inference, not a measurement.** Wnt maintains
AT2 identity and blocks conversion to AT1 (Nabhan 2018); IL-1beta drives AT2
cells into the DATP state and on to AT1 (Choi 2020). Read naively that places
them at opposite poles of one axis. Nobody has measured Wnt activity and IL1R1
status in the same cell, and three findings complicate the simple picture:

* **Katsura 2019** (doi:10.1016/j.stemcr.2019.02.013) found IL-1 and TNFalpha
  *enhance* AT2 proliferation while preserving differentiation capacity. Acute
  IL-1 looks stemness-promoting and chronic IL-1 looks differentiation-blocking,
  so neither is cleanly the anti-Wnt pole.
* **Aumiller 2013** (doi:10.1165/rcmb.2012-0524OC) found that canonical Wnt in
  primary AT2 cells *induces IL-1beta*, at about a two-log fold change after
  WNT3a, confirmed by qPCR, ELISA and in vivo reporter. **This is the strongest
  mechanistic bridge between the two axes that exists, it runs Wnt to IL-1
  rather than the reverse, and the alveolar regeneration literature almost never
  cites it.**
* **IL-1 to Wnt in AT2 is not established.** Choi's Discussion proposes it and
  does not test it, and no subsequent paper does.

The niche geography is also unreconciled: the Wnt source is a single
Pdgfra-positive fibroblast in a juxtacrine niche, the IL-1beta source is an
interstitial macrophage, and nobody has imaged whether one AT2 cell can occupy
both niches. That is the spatial version of the same question.

---

## 4. Why public data cannot settle it

* The GEO DataSets intersection of Axin2 and Il1r1 returns **zero records across
  all organisms**. The two literatures are disjoint by construction.
* **No Il1r1 reporter mouse has ever been deposited in GEO.** Every search for
  an Il1r1 lineage, CreERT2 or tdTomato reporter returns nothing.
* Axin2 appears in GEO only as a reporter sort, always at homeostasis or in
  development. Il1r1 appears only as an antibody sort (GSE144598, coverage
  tracks only) or as a genetic deletion (GSE247503 and GSE247504).
* A PubMed title and abstract search for Axin2 co-occurring with IL-1 or Il1r1
  returns exactly one paper in all of PubMed, and it is about oligodendrocytes.

**And the transcripts cannot stand in for the reporters.** Measured on Choi's
own deposit (GSE145031), Axin2 is detected in 4.6, 4.3 and 5.6 per cent of
AT2-lineage cells at PBS, day 14 and day 28, at about one UMI per positive cell;
Il1r1 in 5.8, 7.1 and 7.1 per cent at about 1.1 UMI. At one molecule per
positive cell a per-cell call is a Poisson coin flip rather than a phenotype.
Nabhan's own supplement concedes the related point for the reporter: single-cell
sequencing of GFP-negative AT2 cells confirmed the Axin2-CreERT2 reporter does
not label all Axin2-positive AT2 cells, so the reporter is a
recombination-limited floor rather than a measurement of the fraction. The paper
also never publishes a discriminating transcriptomic signature for the labelled
cells, so there is no classifier to port.

---

## 5. What can still be done, with costs

**Route A. Re-ask England's Figure S2 properly, on data already on this machine.**
GSE247505 is in `raw_data/`. The only direct evidence in the field is one hedged
supplementary panel at homeostasis; this repository could ask whether Axin2 is
higher in the Il1r1-lineage fraction with the **mouse** as the unit, and whether
it survives the Kras and injury arms rather than homeostasis alone.
*Cost: low, data local. Ceiling: two libraries per group, which is the ceiling on
the Il1r1 side everywhere.*

**Route B. Ask it in chromatin rather than in transcript, and this is the one
this branch measured.** Accessibility aggregated over a gene's linked peaks is
much less sparse than its transcript, which is exactly how this escapes the
floor that kills Route A's RNA version. Measured in
[`feasibility_of_the_joint_question/`](feasibility_of_the_joint_question/feasibility_summary.md):

| Gene | transcript detection | any linked peak | gain |
|---|---|---|---|
| **Axin2** | 5.3 and 3.8 per cent | **22.1 and 39.4 per cent** (15 to 16 peaks) | **4.1 and 10.4 fold** |
| Il1r1 | 24.9 and 18.7 per cent | 29.9 and 51.6 per cent (39 to 40 peaks) | 1.2 and 2.8 fold |
| Sftpc (control) | 100 per cent | 14.5 and 27.4 per cent (2 peaks) | 0.1 and 0.3 fold |

The controls order correctly in both deposits, which is what makes the
measurement trustworthy: chromatin helps exactly where the transcript is sparse
and the gene carries many linked peaks, and it hurts where the transcript is
already universal. The multiome deposits this repository already holds
(GSE310539 and GSE247130) carry both loci in the same nuclei.

**RUN, AND CLOSED. See [`ROUTE_B_OUTCOME.md`](ROUTE_B_OUTCOME.md).** Three
passes with three instruments (A1, A1b, A1c) all returned NOT COMPUTABLE,
because the positive control that gates the reading never cleared. The gate is
now known to have been mis-specified, so the closure is a decision taken under a
rule frozen before the trial ran rather than a demonstration that the approach
cannot work. What the passes did establish is that the graded form of the
statistic has real sensitivity: Etv5 against Abca3 clears at z = +3.55 and Krt8
against Krt18 at z = +3.59. The outcome page states what a correct gate would
be. *Register status: Not established, never Refuted.*

**Route C. WITHDRAWN, and the reason is worth more than the route was.** An
earlier version of this page proposed GSE150957, a bulk array of Wnt-high,
Wnt-low and Wnt-negative distal lung epithelium sorted from the same animal
across six independent experiments, on the strength of its being the only paired
design anywhere with the mouse as the unit.

The owner pointed out that it is out of focus, and it is, decisively:
**bulk cannot answer a co-occurrence question at any level of replication.** The
question is whether two markers sit in the same cell; a sorted fraction reports
an average over cells and is silent on that. Three further mismatches sit on top
of the fatal one: it sorts on TCF/Lef:H2B-GFP, a Wnt activity reporter, which is
not the Axin2 population Nabhan and Zacharias defined; it is whole distal
epithelium rather than AT2; and it is elastase and emphysema, a fifth injury
model in a literature already confounded by four.

The route was ranked here partly because it had the mouse as its unit, in a
project where almost nothing does. That is a virtue of its statistics and not of
its relevance, and **replication of a measurement that cannot address the
question is worth very little.** Recorded as row C146 rather than deleted.

**Route C prime, the replacement, and it is same-cell.** The assessment refused a
transcript route because Axin2 is detected in 4 to 5 per cent of AT2-enriched
cells at about one molecule each. A Wnt-target MODULE is not Axin2. Measured in
[`module_detection_check/`](module_detection_check/module_summary.md) on
GSE262927, the one deposit in this project with genuine per-animal replication:

| | detection across six samples |
|---|---|
| Axin2 alone | 11.1 to 18.0 per cent |
| **two or more of six Wnt target genes** | **14.6 to 28.4 per cent** |
| **Il1r1 or Il1rap** | **19.9 to 40.2 per cent** |

Both sides clear the floor, per cell, in an atlas with replicate libraries per
timepoint. Three cautions, all of them earned:

1. **Clearing the floor is necessary and not sufficient.** Route B cleared its
   own by a factor of 4.1 to 10.4 and still closed after three passes, because
   the spread of the matched null swamped an effect its own control showed was
   real.
2. **It is a different question, and must be asked in its own words.** A
   Wnt-target module reports current target expression; the reporter reports
   transcription during a tamoxifen window plus everything descended from it.
   This asks whether cells with active Wnt target expression also carry IL-1
   receptor components. It does not ask whether Choi's Il1r1-lineage cells and
   Nabhan's Axin2-lineage cells are the same cells.
3. **Any trial needs a positive control that is demonstrated rather than
   asserted**, which is what trials A1 and A1b paid for: one control was assumed
   to be an arithmetic identity and was not, and its replacement was asserted to
   be unable to fail and did.

*Cost: low, data local, per-animal replication available. Not yet run.*

**Route D. The experiment that would settle it.** Two CreERT2 drivers cannot be
independently co-induced in one mouse, so a joint trace is not available. It
needs one CreERT2 lineage plus a non-Cre transcriptional reporter (Axin2-lacZ or
Axin2-GFP), or an intersectional Dre and Cre design, scoring four quadrants
within SPC-positive AT2 at homeostasis and after injury, with the mouse as the
unit and thresholds frozen from uninjured controls. The reporter-artefact
control Choi's own Limitations specifies, a tamoxifen washout longer than 16
days, would have to run alongside, and would also adjudicate the Nabhan against
Zacharias discrepancy.
*Cost: new mice. This is the honest answer to "how would you actually know".*

---

## 6. The caveat Route B carries, stated before anyone runs it

**An accessible Axin2 locus is not a Wnt-responsive cell.** Accessibility reports
that a locus is in a configuration permitting expression, not that the cell is
signalling now, and it is slower and more permissive than transcription. Route B
is a weaker proxy for the lineage reporter than the transcript would be if the
transcript worked. It is on the table because the transcript does not work, and
any trial built on it must say so in its own pre-registration rather than in a
footnote.

The same caution that governs the neighbouring
[`datp_epigenetics/`](../datp_epigenetics/README.md) branch applies here, and it
is row B3 of that branch's claims table: the source papers define their
populations by a method transcript counting cannot reproduce, and a reanalysis
that pretends otherwise is measuring a different object under a borrowed name.

---

## 7. One thing worth knowing before the analysis contract is written

Roadmap paper 13, Yu, Lee, Choi and Choi 2026 in *Seminars in Immunology*
(doi:10.1016/j.smim.2026.102050), is the local PI's current framing of IL-1beta
and is flagged in the roadmap as read-before-the-contract. Its 162-item
reference list cites Choi 2020 and contains **no Nabhan 2018, no Zacharias 2018,
no Zepp 2017 and no Frank 2016**, and neither Wnt, Axin2 nor beta-catenin
appears in its title, abstract, author keywords or MeSH terms. On that evidence
the question the laboratory posed in 2020 has dropped out of its current framing
rather than being advanced.

**Stated limit:** the review body is closed access with no PMC record and no
local copy, so this rests on the resolved reference list and the indexed
metadata, not on reading it. A passing textual mention carrying no citation
cannot be excluded, so the claim is specifically that the AEP and Axin2 primary
literature is absent from the bibliography.
