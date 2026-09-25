# Prepared requests for evidence that remains unavailable

25 September 2026. These are reviewable drafts, **not sent messages**. Public
primary records have already resolved the HPCS mouse crosswalk and Hopx harvest
time; neither is requested again. Sending author messages requires a separate
explicit instruction from the owner.

## TIGIT / HPCS 2020 data team

We are reanalyzing GSE154966 using the published mouse-paired TIGIT contrast.
Could you provide the animal IDs contributing to each ATAC source block
(`YY1181`, `YY1916`, `106623`, and `106621_106642`), whether each block is one
animal or a pool, and whether any animals occur in more than one block? Please
also confirm that each block's positive and negative gates came from the same
biological source. The source labels and paired model alone do not resolve
membership of the compound block.

## HPCS 2026 data team

Supplementary Table 4 now resolves all 22 retained traced-cell source aliases
to mouse tags and confirms IGO17543 harvest at 14 weeks. Is current mScarlet
fluorescence available for the trace-sorted cells in the combined scRNA object,
with a cell barcode or index-sort crosswalk? We need to distinguish permanent
lineage inheritance from ongoing reporter-positive state. If recorded, please
also clarify why BO1534, BR1311, BL1241 and BH1719 have no retained traced alias
in the published combined traced subset. An absent alias is not being treated
as a biological zero. New chase inference would additionally require a design
with sequencing libraries spanning multiple chase groups.

## PATS / MintChIP data team

For the deposited GSE141635 histone and H3 bedGraphs, could you provide the exact
per-library normalization procedure and scale factors, matched mark-to-H3
associations, treatment of uncovered bins and any spike-in scaling? We are
keeping the day-12 histone/sort experiment separate from the day-8 TP53/input
experiment and from lineage-tracing endpoint panels.

## Tsutsui p300/CBP data team

Full SRA sample records resolve the 26 PRJDB37980 and 32 PRJDB37983 CUT&Tag
libraries. Could you provide per-preparation peak counts or normalized tracks,
the sample-to-H3 mapping and E. coli scaling factors used for Figs. 8/9?
Their captions specify spike-in normalization, whereas baseline deposited
tracks use CPM. We would also appreciate preparation IDs linking chromatin,
RNA and phenotypic panels, if such experiments were actually paired. Separate
experiments will remain separate in our analysis.

## TP53 / Mdm2 preprint data team

The manuscript states that GSE335749/335750 are public, but GEO currently shows
both as private until 1 June 2027. Could the full processed count matrices and
sample-to-animal/origin/perturbation manifest be released? The supplementary
3,984/686/868-gene sets include decreases as well as increases, and 109 of the
shared significant genes have opposite signs. Please clarify the “upregulated”
wording and provide unfiltered DE results if available. Full inputs would allow
a direct origin-by-perturbation interaction, avoiding a comparison of separate
significance calls. Mouse-level lineage/imaging counts with nested-field IDs
would permit independent verification of the measured outcome.
