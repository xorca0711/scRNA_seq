# Nabhan 2018: reproduction from deposited adult mesenchymal FPKM

**Status: descriptive source-data reproduction, completed 2026-09-22.**
The deposited matrix recovers frequent Wnt5a expression alongside fibroblast
markers. The exact published statement that 74% of Wnt5a-positive cells express
Pdgfra is **not reproduced under the four declared thresholds**: the conditional
fraction ranges from 72.2% to 90.3%. An unstated positivity cutoff and incompletely
specified original filtering prevent an exact numerical replication. This is
neither independent validation of the paper nor a refutation of its experimental
niche mechanism.

## Cohort and source

[GSE109444](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109444)
deposits 47 adult mouse lung mesenchymal single-cell transcriptomes and one
approximately 200-cell bulk control. The bulk sample, **GSM2943054 / PTC1_BTL28**,
is excluded from every cell fraction and stored separately. The deposited
metadata describe CD45/Pecam/EpCAM depletion followed by Fluidigm C1 capture,
SMARTer amplification and Cufflinks 2.0.2 gene-level FPKM against mm10. These are
FPKM estimates, not UMI counts. No independent-animal mapping is supplied in the
sample metadata; the 47 cells are not treated as 47 biological replicates.

The reference is Nabhan et al., *Science* (2018),
[doi:10.1126/science.aam6603](https://doi.org/10.1126/science.aam6603), particularly
Figure 2D and the single-cell methods in the
[author-hosted paper and supplement](https://desailab.stanford.edu/sites/g/files/sbiybj24296/files/media/file/sc1.pdf)
(combined PDF pp. 3 and 15–16). The paper's 74% denominator is **Wnt5a-positive
cells**, not all cells or Pdgfra-positive cells. Its ligand/marker positivity
threshold is not specified. GEO currently links an additional later publication;
the downloaded metadata version and exact input hashes are preserved here.

## Fixed analysis choices

[specification.json](specification.json) was written after inspecting the source
methods and one input file, before the full-cohort summaries. It is a retrospective
reproduction specification, not a preregistration.

- All 47 deposited cells remain in the primary analysis. The 13-column tracking
  files are **headerless**; the parser preserves their first gene record and
  requires exactly one record for every target symbol in every sample.
- The Figure 2D panel contains Col1a2, Pdgfra, all 19 Wnt genes, Axin2, Ppia, Actb
  and Ubc. Ppia is the actual symbol; the figure caption's apparent “Ppla” is not
  used as a gene. Ten additional genes provide fibroblast, Wnt-processing and
  epithelial/immune/endothelial expression context.
- Detection means FPKM **>0**, **>=0.1**, **>=1**, or **>=5**, with the same cutoff
  for both genes in each conditional co-detection. FPKM >=1 is a transparent
  display convention anchored to the housekeeping QC threshold; it is not
  asserted to be the paper's Wnt positivity rule. No cutoff was selected to
  approximate 74%.
- The heatmap uses **log2(1+FPKM)**, with cells ordered by decreasing Wnt5a, then
  Pdgfra, then accession. The paper labels its transform “log2 FPKM” without
  specifying zero handling. Its clustering and original cell order are not
  reconstructed. No PCA or new clustering is needed to answer this bounded
  expression question.
- The paper reports removing cells failing detection of at least three of four
  housekeeping genes (Actb, Gapdh, Ubc, Ppia; FPKM >=1), or falling more than
  three standard deviations below the mean. The latter scale/implementation and
  rejected-cell universe cannot be recovered. All 47 deposited cells pass the
  recoverable three-of-four diagnostic; no additional cell exclusion is applied.

## Results

| FPKM positivity rule, both genes | Wnt5a-positive cells | Wnt5a/Pdgfra co-detected | Pdgfra-positive among Wnt5a-positive |
|---|---:|---:|---:|
| >0 | 31 | 28 | 90.3% |
| >=0.1 | 31 | 28 | 90.3% |
| >=1 | 26 | 21 | 80.8% |
| >=5 | 18 | 13 | 72.2% |

At FPKM >=1, Wnt5a is detected in 26/47 cells, Wnt2b in 9, Wnt11 in 7,
Wnt2 and Wnt9a in 3 each, and Wnt4 in 2. Other Wnts are not uniformly absent:
Wnt5b, Wnt6 and Wnt10b each occur in one cell. The full 19-gene panel avoids
turning the paper's selected examples into an exhaustive ligand list.

Col1a2, Col1a1 and Mgp are detected in all 47 cells; Pdgfra in 33. Axin2 is
detected in 6, including 5 of the 35 cells with any Wnt detected. Porcn occurs
in 5 and Wls in 36; among the 35 any-Wnt-positive cells the respective counts
are 3 and 29. These are transcript co-detection observations, not functional
evidence that a cell secretes an active Wnt or signals to a nearby recipient.

The additional context panel exposes an unresolved source-data caveat:
**Sftpc is detected in 46/47 cells** (median 11.20 FPKM), despite strong fibroblast
markers, whereas Epcam is detected in 4 and Pecam1, Ptprc and Ager in none at
FPKM >=1. This matrix cannot distinguish epithelial transcript carryover,
mixed-cell material, biological expression or other technical explanations.
We do not silently discard these cells, assign them AT2 identity, or claim a
doublet/ambient-RNA correction that the processed data cannot justify.

![Source expression panel](figures/source_expression_panel.png)

![Threshold sensitivity](figures/detection_sensitivity.png)

## Interpretation and limits

The defensible retained observation is: **Wnt5a is the most frequently detected
Wnt in this deposited adult mesenchymal cell panel, and a majority of
Wnt5a-positive cells also express Pdgfra across the four declared cutoffs.**
The exact 74% is not recovered as an executable result. Its discrepancy can
reflect a different threshold, filtering, matrix version or other unpublished
implementation choice; this analysis cannot select among these explanations.
The qualitative observation should not be promoted to modern fibroblast subtype
identity, proximity to an AT2 cell, independent-animal replication, or proof of
a Wnt-dependent niche. Those require the paper's separate spatial and perturbation
experiments or a new study. FPKM magnitudes and zero rates should not be compared
directly with UMI-based external cohorts without an explicit measurement model.

## Reproduction and evidence

Download the three original files listed with URLs and SHA-256 hashes in
[provenance.json](provenance.json) into `raw_data/GSE109444/`. The raw directory
and source PDF remain gitignored. The analysis reads the tar directly without
extracting or modifying members. Run with Python, NumPy, pandas and Matplotlib:

```powershell
$env:PYTHONPATH='X:/GitHub/scRNA_seq/.venv-x64/Lib/site-packages'
$python='C:/Users/dream/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $python "Research Article/gate1_03_nabhan_2018/source_reproduction/reproduce.py"
& $python -m unittest discover -s "Research Article/gate1_03_nabhan_2018/source_reproduction" -p test_reproduction.py
& $python "Research Article/gate1_03_nabhan_2018/source_reproduction/verify_outputs.py"
```

The first command regenerates derived artifacts. The verifier is read-only.
Four focused tests cover headerless input, bulk exclusion, threshold boundaries
and conditional denominators. Saved cohort/count summaries and all input,
code and output hashes are checked by the verifier. Both PNG figures were
visually inspected at delivery; SVGs are provided for export.

| Artifact | Purpose |
|---|---|
| [summary.json](summary.json) | Machine-readable result and limitations |
| [provenance.json](provenance.json) | URLs, exact input/member hashes, versions and output hashes |
| [single_cell_fpkm.csv](tables/single_cell_fpkm.csv) | 47-cell by 35-gene matrix |
| [bulk_control_fpkm.csv](tables/bulk_control_fpkm.csv) | Separate bulk control; never a cell |
| [panel_gene_records.csv](tables/panel_gene_records.csv) | Exact symbol, locus, FPKM bounds and status for all panel observations |
| [sample_metadata.csv](tables/sample_metadata.csv) | Accession/title, source type, missing animal mapping |
| [detection_by_threshold.csv](tables/detection_by_threshold.csv) | All 35 genes at four fixed thresholds |
| [conditional_codetection.csv](tables/conditional_codetection.csv) | Explicit joint counts and conditioning denominators |
| [housekeeping_qc_diagnostic.csv](tables/housekeeping_qc_diagnostic.csv) | Recoverable partial QC, not a reconstructed full filter |
| [parse_validation.csv](tables/parse_validation.csv) | Preserved first gene, row counts and target-status audit |
| [plot_cell_order.csv](tables/plot_cell_order.csv) | Accession map for every plotted column |
| [Expression SVG](figures/source_expression_panel.svg), [sensitivity SVG](figures/detection_sensitivity.svg) | Editable vector figures |

No root claim, roadmap, other analysis branch or original input was changed by
this reproduction.
