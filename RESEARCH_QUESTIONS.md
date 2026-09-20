# Research questions and the phenotypes behind them

This repository asks which epithelial and immune-state programmes distinguish
productive lung repair from persistent remodelling after injury, by re-analysing
fifteen public single-cell and multiome deposits from eleven studies of lung
injury, repair, fibrosis and tumour initiation, in mouse and human, with the
authors' annotations held out of every fitting step. Five questions below are
ones it can now ask well. Each states what the data say, what they do not, and
a piece of work sized to a semester that would move it.

This page is written for a principal investigator reading cold, organised by
question rather than by trial. Every number is read from a tracked table and
carries its row in the claims register, [`CLAIMS.md`](CLAIMS.md); everything
refuted or unestablished is in [`NEGATIVE_RESULTS.md`](NEGATIVE_RESULTS.md);
the trial that produced each number is one link away.

**How to read the statuses.** Rules are pre-registered and frozen before data
are opened; a refuted claim stays on display; the statistical unit is the animal
or donor, never the cell. Of 154 register rows, 8 are validated, 65 descriptive
or exploratory, 38 not established, 42 refuted or retracted and 1 displaced
(established outside this repository). That ratio is the point: the product is
calibrated evidence, and the negative results are load-bearing.

---

## Part A. Five questions the repository can now ask well

The pieces of work named under A1 and A3 are Part D's D1 and D2; the A2 spatial
survey and the A4 mouse design are not semester deliverables and are labelled as
such; A5 names its own.

### A1. When an AT2 cell enters the transitional state, is its AT2 identity programme closed at the chromatin level, or only silenced?

*Why it matters.* The transitional (Krt8-high, Cldn4-positive) state has been
named three times (DATP, PATS, ADI). Whether it is a reversible waypoint or an
arrested dead end is the difference between a regeneration intermediate and a
fibrosis seed, and neither of the two source papers that hold the right data
tests it (the branch note quotes what each says).

*What the data say.* Across two 10x multiome deposits (GSE310539, GSE247130),
the Cldn4-positive Krt8-positive alveolar group loses the AT2 identity
programme in RNA by 12.3 to 16.7 detection points at every one of four resampling
seeds (12.6 and 17.1 at the single M1e budget), carried by five or more genes
(Etv5, Abca3, Slc34a2, Napsa, Lamp3), with the AT1 programme flat in the same
cells (C118). On one common instrument and
a null built over 300 peak-count-matched random gene sets, the same group's AT2
distal chromatin sits at percentile 0.000, 0.007 and 0.000 in every well that
reached the statistic, at 15 to 38 per cent below reference accessibility
(C131). Where the statistic can be read at all it leans toward closing, but the
uninjured control well leans the same way at a larger effect (18 against 15 per
cent) than one of the two test wells, so the lean is Exploratory and not yet
label-specific.

*What they do not say.* The chromatin question is Not established (C121): no
well clears both its positive-control gate and its AT2 gate at once (C133), one of the three wells is an uninjured
negative control showing the same direction, and the two deposits share a
laboratory, so agreement is consistency rather than replication. A first
reading of "silenced but not closed" was retracted when its positive control
turned out to fire in five of eight resampling combinations (C120, C126).

*The semester piece (D1).* GSE309751, the bulk ATAC-seq accompanying one of
those papers (assessed and not opened; see REFERENCES.md), deposits per-sample
peak calls for two to three mice per group at 14 and 49 days after injury. It
is the only chromatin contrast in this project with animal-level replication.
Whether any AT2 identity loci that are less accessible at 14 days than in mock
are so again at 49 is a direct test of reversibility, and the first test of
closing itself with the mouse as the unit. Two limits are stated now: bulk
sorted AT2 dilutes a state-specific effect roughly eightfold, and the 49-day
groups hold two mice. A normalised differential test needs the reads from SRA
and is out of scope; peak presence at frozen loci is not.

### A2. Where does amphiregulin come from in the human lung, and how much of a ligand-receptor ranking is the database rather than the tissue?

*Why it matters.* Amphiregulin-to-EGFR is the axis a 2026 *Nature* paper places
at the centre of early fibrotic niches in tumour initiation, read from an
epithelium-centric sort. Ligand-receptor inference is the most reused and least
audited step in single-cell work.

*What the data say.* In human lung adenocarcinoma (GSE131907), AREG detection
is higher in epithelial than in myeloid cells within the same donor, 0.336
against 0.215 over 11 donors, paired Wilcoxon p = 0.0020 (C37, the one tested
claim in that branch). But dendritic cells and monocytes carry AREG and HBEGF
at or above the epithelial states in three human datasets across two diseases,
a source that no sort of labelled epithelium can see (C45). Holding cells,
donors and scoring fixed and varying only the curated resource, a
whole-database scan surfaces matrix-to-CD44 pairs at the head in one resource
of five (abundance share 0.80, then 0.40, 0.27, 0.20, 0.20), so that
domination is CellChatDB's, not the tissue's (C111); CellChatDB shares zero of
its top fifteen pairs with italk (C113). AREG is first among the EGFR ligands
in every resource that contains them (C112).

*What they do not say.* Enrichment of AREG or HBEGF in tumour over matched
normal lung is refuted at compartment level in adenocarcinoma (C40), and whether
either rises in fibrotic over control lung is Not established (C48). No
donor-level coupling of the epithelium-to-fibroblast axis was detected across
22 donors, 18 of them fibrotic (rho 0.348, p = 0.112 pooled; rho 0.276,
p = 0.268 within the 18), a null that bounds nothing weaker than rho 0.43, so
the axis is Not established rather than absent (C49). A communication inference
without proximity evidence is a hypothesis.

*Not a semester deliverable.* What would move this is spatial transcriptomics
of injured or fibrotic mouse lung with more than one animal per condition, to
ask the one question this axis has never been asked: whether the
ligand-producing and receptor-bearing cells are near each other. The survey for
such a dataset is a week; the analysis depends on what it finds.

### A3. Which injury-associated states are still present a year after repair, and where did they come from?

*Why it matters.* Persistence after apparent recovery is what separates repair
from remodelling. The question needs a long time course with lineage tracing,
which is rare.

*What the data say.* In a 25-animal mouse time course to 366 days (GSE262927;
33 libraries including an eight-sample lineage-driver cohort), an
injury-associated capillary endothelial state is nearly absent at baseline
(median per-animal 2.0 per cent), peaks at 37.5 per cent at 25 dpi, and is still
21.7 per cent at 366 dpi (C3, Validated on artefact-checked per-animal numbers). The Kit lineage
traces that state at 33 to 53 per cent per animal, supporting a CAP1 origin; the
CAP2-specific drivers label 2 to 8 per cent of endothelium and are uninformative
rather than negative (C4). In the myeloid compartment, alveolar macrophages
fall from 30.8 to 4.6 per cent of myeloid cells at 6 dpi and rebuild to 49.0 per
cent by 42 dpi while inflammatory monocytes rise from 1.9 to 56.0 per cent and
return; the 2 to 3 dpi tamoxifen window labels most of the rebuilt pool (median
79.7 per cent against 29 to 60 for the three later windows) (C12, C15); and
interstitial macrophages rise through 90 dpi (2.8, 10.3, 14.0 per cent) before
falling back to 7.7 at 366, an Exploratory reading on three to four animals
(C13).

*What they do not say.* The proliferation phases and the macrophage rebuild are
Descriptive only, with two animals per active-repair day and no P value. The
marrow-inheritance check on the rebuilt pool is Not established under a frozen
30-cell floor, with references evaluable in 3 of 8 animals (C16). The injury model is a respiratory virus; the analysis reads
cell states, niches and kinetics, not the infection.

*The semester piece (D2).* Sample-aware pseudobulk of the myeloid compartment
per animal and phase (proposal W1), defining the reconstituting macrophage
programme and comparing it with the published myeloid-to-mesenchymal ARG1 and
ornithine circuit of lung fibrosis.

### A4. Are the Wnt-responsive and the IL-1-responsive AT2 cells the same cells?

*Why it matters.* One literature says a Wnt-responsive (Axin2-positive) AT2
subset holds stemness and is restrained from AT1 conversion; another says an
IL-1-responsive (Il1r1-positive) subset is what enters the transitional state.
If both mark the facultative progenitor they pull opposite ways on one fate
decision. The paper that defined DATP closes by proposing exactly this
comparison. The one measurement since is a bulk qPCR panel on sorted cells at
homeostasis (England et al. 2025, C136); nobody has made it with the mouse as
the unit or after injury.

*What the data say.* Neither population is established as a distinct subset:
both are tamoxifen-inducible lineage reporters, the Axin2 fraction differs 20
to 30 fold between two knock-in alleles and nobody has explained the gap in
eight years, and both papers' own data lean toward an injury-inducible state
(C134, C135). The only direct evidence for overlap is a three-gene bulk qPCR
panel; it was never a sequencing experiment, so it cannot be re-examined (C136,
C148). No deposit carries the two reporter-defined populations in the same
cells; the registered form of that claim was refuted by this branch's own
artefacts, which measure both transcripts in the same cells, and the narrower
statement is what survives (C137). Both transcripts sit at four to seven per
cent detection at about one molecule per positive cell, which is a coin flip, not a
phenotype (C138). One mechanistic bridge exists and is almost never cited:
canonical Wnt induces IL-1beta in primary AT2 cells (Aumiller 2013). The
reverse direction, which the DATP paper hypothesised, is untested anywhere.

*What they do not say.* A locus-accessibility route was opened (Axin2's linked
peaks are detected 4 to 10 times more often than its transcript, C139) and then
closed after three passes, each under a rule frozen before that pass ran,
because no pass's positive control cleared (C142). The Il1r1 conditional deletion is
not visible in the deposited count matrices of the one deposit that has it, and
why is Not established, so that deposit's genotype labels rest on metadata
alone (C149).

*Not a semester deliverable.* This one is honest about needing mice: one CreERT2
lineage plus a non-Cre Axin2 reporter, four quadrants within SPC-positive AT2
at homeostasis and after injury, with the tamoxifen washout the DATP paper's
own limitations section specifies. It fits an 18-week window only as a
computational design and power calculation, which is offered as such.

### A5. Is the transitional-state marker set damage-associated, or a developmental programme re-used?

*Why it matters.* "Damage-associated" is in the state's name. If its defining
markers are on in normal postnatal lung, the name is a hypothesis.

*What the data say.* A Cldn4-positive Krt8-positive transcript call labels 3.69
and 8.07 per cent of cells in the two uninjured neonatal (P9) wells. The larger
of the two exceeds every injured adult well (at most 5.09 per cent) and the
smaller exceeds three of the four, so a whole-trial negative control fired and
the trial refused (C119). Krt8 and Cldn4 are expressed across
immature alveolar epithelium; the marker set is shared with normal development.
The paper this is read against argues that neonatal AT2 cells are plastic and
never states it as a limitation of the marker.

*What they do not say.* At transcript level, with these two genes, the state
has no character specific to regeneration. Whether a fuller module separates
neonatal from injury-induced cells is untested.

*The semester piece.* Score the published DATP, PATS and ADI gene sets, plus a
frozen developmental AT2 maturation set, across the P9, seven-week and infected
wells of GSE247130, asking which genes are injury-specific once development is
held fixed. It stops if the developmental set does not separate the P9 well
from the seven-week control under a threshold frozen from the seven-week well.

---

## Part B. Remarkable bioinformatic phenotypes

Things the data did that a computational biologist would want to know about,
independent of the biology they were found in. Each is a general lesson with a
number behind it.

1. **A response state masquerading as a cell type.** Blind clustering recovered
   the deposited mouse cell types at median purity 0.947 over 107,626 labelled
   cells, and was contradicted for 3 of 29 clusters. The instructive one is a
   cluster the marker panel called transitional epithelium that is 89 per cent
   capillary endothelium whose top data-derived markers are interferon-stimulated
   genes (Iigp1, Gbp4, Isg15): a response state, not a lineage, reported here as
   what the data show. Kept on display as a failure mode of panel annotation
   (C1, C8).

2. **Sequencing depth is the pervasive confound, and it hides in every
   detection fraction.** Three separate refutations in this repository were
   depth in disguise: a magnitude criterion that turned out to measure library
   depth (C69), a presence rule on an ambient transcript that removed a median
   of 99.9 per cent of nuclei across ten wells with a severity that tracked
   depth (C123), and a peak-counting statistic refused before it ran, because
   the labels and the chromatin come from one nucleus and a per-peak detection
   fraction is monotone in depth, so its null was false by construction (C127).
   The correction adopted was to downsample every cell to one budget and to
   build the null over matched gene sets rather than over cells (C130); under
   it the chromatin question is still Not established (C121).

3. **Curated ligand-receptor databases disagree at the head.** Same cells, same
   donors, same scoring: CellChatDB shares zero of its top fifteen pairs with
   italk and three with cellphonedb; three of the four non-CellChat resources
   agree at Jaccard 0.43 to 0.58 and cellphonedb shares none of its top fifteen
   with any of them (C113). A top-hits list from such a scan is at
   least as much a property of the resource as of the tissue.

4. **Deposit errors found by corroborating the metadata against the data.**
   A multiome deposit's barcode-suffix order is inverted relative to its GEO
   sample order: the conditional knockout carries 4.5 to 14.5 times more of its
   own target gene than its control in all three files; Sox9 agrees in all
   three and Cldn4 in the infected one (C116). A lineage-traced deposit's six
   matrices are raw barcode whitelists and its reporter is not a counted feature
   (C60, C61). A conditional Il1r1 deletion is not visible in the deposited
   counts, and why is Not established, so the genotype labels cannot be checked
   (C149). An ATAC deposit is coverage tracks only. Each of
   these would have inverted or voided a downstream result.

5. **One library per condition, six times.** Six of the deposits analysed
   (eight of the fifteen accessions) carry no within-group replication for any
   between-condition contrast, so nothing in them about genotype, treatment or
   stage is testable and every between-condition result is at most Descriptive
   only (C20, C59, C117). The repository's response is
   the within-library state comparison, which needs no between-group
   replication, and a register that never lets a descriptive row wear a
   validated label.

6. **A marker signature splits into tiers under perturbation.** Under deletion
   of the ligand a paper names as the driver, its six-gene fibrotic signature
   does not move as one unit: Runx1 and Pdgfrb are retained, Fst and Runx2 lose
   most of their detection (C29). The tiers stand; the reading of them as a
   second signal was refuted when Runx1 and Pdgfrb rose with bleomycin injury
   alone in both replicates (C42, trial E4).

7. **A co-expression gate cannot audit a co-expression detector.** A crude gate
   suggested doublet removal was deleting a novel human epithelial population at
   1.26 times the background rate, and a basal population at twice it; a stricter
   gate that could actually separate the first from doublets overturned it (3.9
   against 6.3 per cent), and the "lost" cells carried about twice the UMIs of
   the cells kept (C7).

8. **The transitional state is defined by immunostaining and cannot be rebuilt
   from three-prime counts.** Krt8-high is an antibody call; a transcript label
   requiring Cldn4 detection and Krt8 above the uninjured 99th percentile (3 and
   7 raw counts in the two deposits) returned 0.03 to 0.66 per cent of cells
   across ten wells at a common depth budget, so the published fraction cannot
   be rebuilt from counts. Four trials were spent learning that the
   transcript-level group must be named as its own object (C124).

9. **A register audited by adversaries against its own artefacts.** Eight
   agents, each finding sent to an independent verifier that defaulted to
   refuting it: 45 defects were confirmed or partly confirmed across 14 of the
   35 rows audited, because summary prose had been typed while the dataframe
   holding the true values sat in scope (C151). The trials were mostly honest; the register written from
   them was not, and the fix is mechanical.

---

## Part C. What a visiting student brings to this

The methods are the ones the questions above needed, demonstrated on the
deposits named:

- **Pre-registration as code.** Rules frozen in a run record before data are
  opened, with a clause on every threshold stating what would make its own
  answer unreadable; a corrected pass sits beside a wrong one, never over it.
- **Deposit auditing.** Structural reality checks before any fitting (C0, D0,
  M0): feature spaces, barcode whitelists, reporter contigs, metadata against
  data.
- **Blind pipeline graded against held-out labels**, with per-dataset batch
  decisions made from measured replicate mixing rather than by default.
- **Reference mapping** with scArches to the Human Lung Cell Atlas, reading
  uncertainty per donor rather than transferring labels blindly.
- **Donor- and animal-level statistics** with depth controls, sham bands and
  matched-gene-set nulls; no P value where a group holds two animals.
- **Ligand-receptor inference** run across five curated resources with a guard
  that refuses an unreadable ranking.
- **Single-nucleus multiome**: streaming readers for matrices too large to
  load, locus-level accessibility, depth-equalised co-detection.
- Outside this repository: **imaging quantification** (multiplex
  immunofluorescence, animal-level counts) from prior wet-lab work, which is where several of the questions
  above will have to be settled.

---

## Part D. Three deliverables that fit inside eighteen weeks

Each is on data already public, states its unit, and has a stop condition.

| | Deliverable | Data | Unit | Stops if |
|---|---|---|---|---|
| D1 | Reversibility of AT2 identity chromatin, 14 against 49 days after injury (answers A1) | GSE309751, per-sample bulk ATAC peak calls, 2 to 3 mice per group | mouse | the mock arms do not separate from injured at 14 days under a threshold frozen from mock |
| D2 | The reconstituting alveolar macrophage programme by animal and phase, against the published ARG1 and ornithine circuit (answers A3) | GSE262927, 25 animals | animal | fewer than 3 animals in any tested phase; the two-animal active-repair days enter only as a ranking, never as a tested arm |
| D3 | Ligand-receptor candidates between AT2 cells and the fibroblast niche at 11 dpi, only after pseudobulk (the prerequisite for A2 on mouse data) | GSE262927 | animal | candidate pairs not expressed in both partners in at least 2 animals |

---

## Part E. Limits stated first

- The ceiling on between-condition inference is set by each deposit's design.
  The mouse deposits of the roadmap papers (GSE145031, GSE144468, GSE316241,
  GSE316243, GSE316244, GSE310539, GSE247130) carry one library per condition
  with animals pooled, so nothing between their conditions is tested and the
  register does not let it be (C20, C59, C117). Tested rows come only from
  where the unit is replicated: the human deposits with the donor as the unit,
  GSE132771 and GSE262927 with the animal, and GSE247505 with two libraries per
  group. GSE309751 has animal-level replication and has not been opened.
- The mouse time course is a respiratory-virus injury model. The analysis
  reads repair kinetics, cell states and niches; the infection is not the
  subject.
- Isoform-level questions cannot be asked of any data held here (all
  three-prime short-read). Protein-level questions (surface EpCAM, ADAM17
  activity) are open and need a different assay.
- The Cardoso 2026 rows (C19 to C57, C65 to C84) and every row since await the
  owner's retain-or-reject review, recorded in `DEVELOPMENT.md` when made.

*Generated context: the ledger figure in the README is drawn from the register
by `analysis/scripts/15_claims_ledger_figure.py`; nothing on this page is a
number that does not also appear in `CLAIMS.md`, or in `REFERENCES.md` for a
deposit assessed and not opened.*
