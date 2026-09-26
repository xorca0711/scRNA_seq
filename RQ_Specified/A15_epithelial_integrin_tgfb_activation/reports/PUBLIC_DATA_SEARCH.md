# A15 public data search: can any public dataset test the proposal

27 September 2026. Searched before the question was written, because a question whose test
cannot run must say so on its card rather than in a footnote. The requirement set by A2's
own synthesis is specific: **more than one well or unit per target, positions that vary,
and a measure of activated TGF-beta rather than a transcript proxy.** Repositories
searched: the GEO DataSets index, PRIDE, the Image Data Resource and the BioImage Archive.
Literature searched on PubMed, to establish which assay types would satisfy the readout
requirement and to recover the clinical record.

## Conclusion

**No public dataset can test the proposal, and the binding gap is the readout.** No deposit
in any repository searched pairs an epithelial integrin perturbation with a measurement of
activated TGF-beta. Activated TGF-beta is measured at the bench, by reporter bioassay,
by immunoassay on conditioned medium, or by receptor-proximal signalling in the recipient,
and none of those measurements is deposited as reusable data in the repositories searched.
Transcriptomic deposits can only supply another transcript proxy, which is the thing the
founding observation already had.

**One deposit comes closest and is worth naming precisely.** GSE190821 blocks integrin
beta6 in vivo with a replicated design, and it cannot test the hypothesis: the fibroblast
compartment is not separated and there is no activation readout. It can bound one rival,
and that use is recorded as layer proposal 1 in [RATIONALE.md](../RATIONALE.md), not as
authorized work.

## What was searched, with counts

GEO DataSets, through the Entrez `gds` index. Counts are total hits at the time of the
search; where a query is broad, the number inspected is stated.

| # | Query | Hits | Inspected |
|---|---|--:|---|
| Q1 | `itgb6` | 55 | top 15 |
| Q2 | `itgb6 AND lung` | 9 | all |
| Q3 | `"integrin beta 6" AND (fibroblast OR fibrosis)` | 0 | none to inspect |
| Q4 | `avb6 OR "alpha v beta 6" OR alphavbeta6` | 4 | all |
| Q5 | `(Itgb6 OR ITGB6) AND (knockout OR deficient OR null OR siRNA OR shRNA OR CRISPR)` | 19 | top 15 |
| Q6 | `(Itgb6) AND (single cell OR scRNA OR conditional OR floxed OR Cre)` | 6 | all |
| Q7 | `3G9 OR "anti-integrin beta6" OR "anti-alphavbeta6"` | 68 | top 8 |
| Q8 | `bexotegrast OR PLN-74809 OR BG00011 OR STX-100` | 2 | all |
| Q9 | `"latent TGF-beta" AND activation` | 1 | all |
| Q10 | `(phospho-SMAD OR pSMAD OR "SMAD2 phosphorylation") AND (screen OR perturbation)` | 3 | all |
| Q11 | `(CRISPR screen) AND (TGF-beta reporter OR SMAD reporter OR "TGF-beta signaling reporter")` | 5 | all |
| Q12 | `(co-culture OR coculture) AND (CRISPR screen OR perturbation screen) AND (fibroblast)` | 4 | all |
| Q13 | `(alveolar organoid OR alveolosphere) AND (fibroblast)` | 1,447 | top 12 |
| Q14 | `(alveolosphere OR "alveolar organoid") AND (TGF-beta OR TGFB OR integrin)` | 45 | top 12 |
| Q15 | `(epithelial AND fibroblast AND co-culture) AND (TGF-beta OR TGFB1)` | 82 | top 12 |

Other repositories:

| Repository | Queries | Result |
|---|---|---|
| PRIDE, proteomics | `ITGB6`, `alphavbeta6`, `integrin beta 6` | 3, 0 and 0 projects. The three are integrin-adhesion-complex proteomics in HER2 breast cancer lines and mouse skeletal muscle; none perturbs an epithelial integrin or measures TGF-beta activation |
| Image Data Resource | `TGF-beta`, `SMAD`, `integrin` | Annotation hits only, from existing image screens. No screen pairing an integrin perturbation with a TGF-beta activation or SMAD translocation readout |
| BioImage Archive, through BioStudies | `TGF-beta reporter OR SMAD OR integrin beta 6` | 1,622 hits, top 8 inspected. Structural biology of integrin complexes and unrelated imaging studies |

The Q8 result is worth stating plainly: the two hits are a glioma methylation cohort and
one of its biopsies, matched on an unrelated token. **Neither clinical programme against
this integrin has deposited participant-level data**, so the trial evidence in
[RATIONALE.md](../RATIONALE.md) is published summary results only.

## The six nearest candidates, and why each falls short

| Deposit | Size | What it is | Why it cannot test A15 |
|---|--:|---|---|
| **GSE190821** | 48 libraries, 24 mice | Lung epithelial translatome and whole lung under bleomycin, with the 3G9 anti-integrin-beta6 antibody (4 mice) against the inert antibody Axum8 (7 mice), plus an IRE1a arm; RiboTag with ShhCre, day 7 | The nearest. It has replication and the mouse as the unit, but **no activation readout**, and its only fibroblast-containing sample is whole-lung homogenate, which confounds the fibroblast compartment with fibrosis extent. No ligand arm exists, so the partition cannot be formed |
| GSE298207 | 6 | Reconstructed human epidermis; inflammation driving TGF-beta1 activation through the alphaVbeta6 integrin and mechanotransduction | Skin keratinocytes, no fibroblast recipient in the deposit, and the deposited measurement is RNA |
| GSE224401 | 12 | ITGB6 induced through a TET-On system in breast myoepithelial cells | Gain of function on the epithelial side only, with the epithelial transcriptome as the readout and no recipient cell |
| GSE2255, indexed as GDS1874 | 15 | Alveolar macrophages from the Itgb6-null mouse emphysema model, on arrays | An Itgb6-null lung, but the readout is the macrophage, and the phenotype in that model is Mmp12-driven emphysema rather than fibroblast activation |
| GSE253494 | 4 | Itgb6 interference in renal proximal tubular epithelial cells under hypoxia and reoxygenation | Wrong organ, epithelial readout only, no recipient, and four samples |
| GSE307128 | 4 | A sibling series of the screen A2 used, in the same superseries GSE307351, mouse, type recorded as Other | Four samples, and its content is not a replicated integrin perturbation with an activation readout. Its series summary states the source paper's conclusions; that text is deliberately not used here, since the owner has not read roadmap paper 14 |

## Why the readout is the binding constraint, not the design

The design requirement is not hard to satisfy in principle. Perturbation screens with
several guides or wells per gene are common, and GSE190821 shows that a replicated in vivo
integrin blockade exists. What does not exist in a reusable deposit is the measurement.
Activated TGF-beta has to be distinguished from the latent pool, and the assays that do
that are a reporter cell bioassay, an immunoassay specific to the active dimer, or a
receptor-proximal signalling measurement in the recipient. Those appear in figures, not in
sequence or image archives. PRIDE holds no phosphoproteomics of this axis under an integrin
perturbation, and the imaging archives hold no SMAD translocation screen under one.

This is a statement about what was searched, not a proof of absence. The searches above
are keyword searches over indexed metadata, so a deposit that describes the same assay in
different words could have been missed.

## What this means for the question

1. **Readiness is blocked, and the blocking constraint is nameable.** That is better than
   a vague "needs more data": the card can say which measurement is missing and which
   assay would supply it.
2. **One rival can be bounded today.** GSE190821's epithelial arm can ask whether blocking
   the integrin moves the epithelium's own programme, which is the "changed epithelium,
   not the integrin" rival. It cannot address the hypothesis, and the report says so in
   both places.
3. **The spatial layer is shared with A2 and must be counted once.** Both questions
   propose human spatial transcriptomics, and A2 proposed it first as its own decisive
   layer. Whichever question runs it, it is one cohort audit and one frozen proximity
   definition, not two.

## What was not done

No dataset was downloaded. No author was contacted. The screen's source paper, roadmap
paper 14 (DOI 10.1073/pnas.2606113123), remains unread by the owner and was not consulted;
the sibling deposit's series summary was read while searching GEO and is not relied on.
One structural fact was taken from the deposit's own sample titles and is recorded in
[RATIONALE.md](../RATIONALE.md): the four ITGB6 libraries of GSE307112 are plate3-1 through
plate3-4 at well C02, which is the fixed-position confound stated in A2's withdrawn freeze,
now confirmed for this target by name.
