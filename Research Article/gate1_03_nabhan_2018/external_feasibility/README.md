# Nb1 independent-cohort feasibility gate

Screen frozen 2026-09-22. **Acquired: GSE129605, four saline and four bleomycin day-11 samples. Expression inference remains HOLD until animal provenance, annotation, and compartment coverage pass.** This is a candidate for generalizing Wnt-source/response *expression* across injury models. It cannot establish an autocrine loop, Wnt secretion, spatial proximity, or necessity/sufficiency, and day 11 cannot resolve the acute 1–3-day switch proposed from Nabhan's DT/hyperoxia experiments.

The minimum is three independent animals per arm, each with at least 50 retained AT2 cells and 50 retained fibroblasts. These are feasibility floors, not evidence of adequate statistical power. Cells do not replace animals; count thresholds will not be lowered after looking at Wnt effects.

## Frozen inclusion and acquisition

[Primary GEO GSE129605](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129605) describes whole-lung single-cell isolation on day 11 from male C57BL/6 mice aged 12–13 weeks, after bleomycin or saline. Eight separately indexed samples enter this contrast; three nintedanib-treated samples are excluded. The primary paper is [Peyser et al., 2019](https://doi.org/10.1165/rcmb.2018-0313OC). The deposited numeric identifiers distinguish samples; independent mouse identity still requires a final provenance check, rather than assuming that every GSM is an animal.

| Arm | Deposited IDs | GEO samples |
|---|---|---|
| Saline | 947170, 947172, 955736, 955737 | GSM3716976–GSM3716979 |
| Bleomycin | 947173, 947174, 947176, 955738 | GSM3716980–GSM3716983 |

Exact file URLs and sizes are frozen in `GSE129605_download_manifest.csv`; all 11 records, with inclusion flags, are in `GSE129605_sample_map.csv`. The selected 24 files total **83,692,684 bytes** (79.8 MiB). The complete 11-sample archive is 118,149,120 bytes; it is unnecessary. The primary GEO file list contains counts, genes and barcodes, but no per-cell annotation file. A numerical sample identifier is not treated as a documented batch label.

`fetch_selected.py` fetches only the selected files into `raw_data/GSE129605/`, checks the exact deposited byte sizes, and records SHA-256 checksums. It performs no biological analysis. Run with the documented repository Python launcher or any working Python 3. No package installation is required. It is separate from `build_feasibility.py`, which rebuilds metadata evidence without network access.

All 24 selected files were acquired successfully, containing **13,673 deposited cells across eight samples**. The source-byte totals and SHA-256 hashes are in `GSE129605_download_receipt.json`. `audit_download.py` passed checks of matrix dimensions, positive integer count values, entry counts, and unique barcodes within each sample. `GSE129605_acquisition_audit.csv` records deposited cell counts and explicitly leaves AT2/fibroblast counts **UNESTABLISHED**. The raw matrix cache is local and ignored by Git; compact manifests, receipts, audits, and scripts preserve the reproducible handoff.

**Exact unresolved gate:** the GEO deposit has no per-cell annotation file, and this screen has not established a one-to-one independent mouse mapping beyond separately indexed samples. Thus the intended ≥3 animals per arm × ≥50 AT2 and ≥50 fibroblasts criterion has not passed. No substitute annotation model or marker-based labels were fitted in this bounded task.

## Early-time alternative fails the measured coverage gate

[GSE141259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141259), the whole-lung component of [Strunz et al., 2020](https://doi.org/10.1038/s41467-020-17358-3), initially appears suitable: three day-3 bleomycin samples and seven PBS samples. The small deposited cell-annotation table gives the following coverage. Activated AT2 and myofibroblasts are included in this *generous* screen; their later analytic grouping would require freezing separately.

| Day-3 sample | AT2 including activated | Fibroblasts including myofibroblasts | Both ≥50 |
|---|---:|---:|---|
| muc4631 / GSM4200072 | 20 | 4 | No |
| muc4643 / GSM4200084 | 153 | 16 | No |
| muc4653 / GSM4200091 | 76 | 3 | No |

All seven PBS samples have fewer than 50 fibroblasts, even including myofibroblasts. Only one PBS sample is explicitly day 3 (71 AT2, 2 fibroblasts). Pooling other control days does not rescue the per-animal compartment floor. **HOLD**; no count-matrix download is justified for this contrast. `GSE141259_compartment_coverage.csv` preserves all 28 whole-lung samples, summing to 29,297 cells. Author labels were used directly; no Wnt-based labeling was introduced. Source metadata is only 383,214 bytes; the available whole-lung count matrix is approximately 39 MB but was not fetched.

## Other finite exclusions

| Cohort | Primary evidence | Decision for this injury-versus-control task |
|---|---|---|
| GSE202325 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202325): three young and three aged samples at each of days 3 and 9, all PR8 infected | HOLD: no uninfected arm. Age contrasts answer a different question. |
| GSE292515 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292515): young and old day-7 PR8, each age with three infected and two PBS samples | HOLD: control replication below three within age; combining ages does not make the intended age-matched design. |
| GSE184854 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE184854): two records, WT and CCR2-null, each mixing days 3/7/21, with a separate hashtag spreadsheet | HOLD: demultiplexing and independent control units not established in this bounded screen; two GEO records are not two animals. |

The repository's existing GSE262927 gives two animals at each early time and inadequate AT2 coverage at the inherited 50-cell floor. Existing GSE145031/GSE310539/GSE247130 condition pools cannot supply independent animals; GSE307112 organoid wells are not automatically independent biological preparations. These datasets remain useful for descriptive work but do not fill the requested validation design.

## Stop rules before any Nb1 inference

1. Establish one-to-one mouse/sample mapping and any known preparation batches for the eight selected GSE129605 samples. Do not silently equate barcode libraries with independent animals.
2. Obtain deposited author annotations if available. Otherwise record an annotation requirement and freeze a separately validated annotation procedure. Wnt ligands and Wnt-response genes must not determine the AT2/fibroblast labels used to test the hypothesis.
3. Report the post-QC animal-by-compartment table. Fewer than three eligible animals per arm means HOLD for this comparison. Do not borrow animals from the drug arm or merge cells across animals.
4. If eligible, freeze Wnt-source and response estimands, gene coverage, normalization and uncertainty at animal level before examining direction. A day-11 result supports or challenges expression generalization only. It does not validate Nabhan's acute timing or functional autocrine mechanism.

`feasibility.json` records decisions, remaining gates, source hashes, and exact transfer size. No Nb1 expression hypothesis has been tested by this screen.
