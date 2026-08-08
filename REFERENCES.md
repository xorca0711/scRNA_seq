# References

The primary paper, the second-dataset paper, and the five method papers that make up the non-Seurat components of the primary paper's pipeline. Metadata verified against PubMed.

---

## Target paper

**Niethamer TK, Planer JD, Morley MP, Babu A, Zhao G, Basil MC, Cantu E, Frank DB, Diamond JM, Nottingham AN, Li S, Sharma A, Hallquist H, Levin LI, Zhou S, Vaughan AE, Morrisey EE.**
*Longitudinal single-cell profiles of lung regeneration after viral infection reveal persistent injury-associated cell states.*
**Cell Stem Cell** 2025;32(2):302–321.e6.
DOI: [10.1016/j.stem.2024.12.002](https://doi.org/10.1016/j.stem.2024.12.002) · PMID [39818203](https://pubmed.ncbi.nlm.nih.gov/39818203/) · PMC [PMC11805657](https://pmc.ncbi.nlm.nih.gov/articles/PMC11805657/)

Preprint: *A longitudinal atlas of post-viral lung regeneration reveals persistent injury-associated cell states.* bioRxiv [2024.05.24.595801](https://www.biorxiv.org/content/10.1101/2024.05.24.595801v1)

**Data:** GEO [GSE262927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE262927)
**Mouse lines:** MMRRC 72089 (Car4<sup>CreERT2</sup>), 72090 (Ednrb<sup>CreERT2</sup>)

---

## Second dataset

**Kadur Lakshminarasimha Murthy P, Sontake V, Tata A, et al.**
*Human distal lung maps and lineage hierarchies reveal a bipotent progenitor.*
**Nature** 2022;604(7904):111–119.
DOI: [10.1038/s41586-022-04541-3](https://doi.org/10.1038/s41586-022-04541-3) · PMID [35355018](https://pubmed.ncbi.nlm.nih.gov/35355018/)

**Data:** GEO [GSE178360](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178360) — the three healthy-donor distal-lung samples analysed in [`analysis/GSE178360/`](analysis/GSE178360/README.md)

---

## Method papers, in pipeline order

### 1 · Ambient RNA removal — SoupX

**Young MD, Behjati S.**
*SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data.*
**GigaScience** 2020;9(12):giaa151.
DOI: [10.1093/gigascience/giaa151](https://doi.org/10.1093/gigascience/giaa151) · PMID [33367645](https://pubmed.ncbi.nlm.nih.gov/33367645/) · PMC [PMC7763177](https://pmc.ncbi.nlm.nih.gov/articles/PMC7763177/)
Software: [github.com/constantAmateur/SoupX](https://github.com/constantAmateur/SoupX) · Licence: CC BY 4.0

### 2 · Doublet detection — Scrublet

**Wolock SL, Lopez R, Klein AM.**
*Scrublet: computational identification of cell doublets in single-cell transcriptomic data.*
**Cell Systems** 2019;8(4):281–291.e9.
DOI: [10.1016/j.cels.2018.11.005](https://doi.org/10.1016/j.cels.2018.11.005) · PMID [30954476](https://pubmed.ncbi.nlm.nih.gov/30954476/) · PMC [PMC6625319](https://pmc.ncbi.nlm.nih.gov/articles/PMC6625319/)
Software: [github.com/AllonKleinLab/scrublet](https://github.com/AllonKleinLab/scrublet)

### 3 · Doublet detection — scds

**Bais AS, Kostka D.**
*scds: computational annotation of doublets in single-cell RNA sequencing data.*
**Bioinformatics** 2020;36(4):1150–1158.
DOI: [10.1093/bioinformatics/btz698](https://doi.org/10.1093/bioinformatics/btz698) · PMID [31501871](https://pubmed.ncbi.nlm.nih.gov/31501871/) · PMC [PMC7703774](https://pmc.ncbi.nlm.nih.gov/articles/PMC7703774/)
Software: Bioconductor [scds](https://bioconductor.org/packages/scds/) (doi:10.18129/B9.bioc.scds)

### 4 · Trajectory inference — Slingshot

**Street K, Risso D, Fletcher RB, Das D, Ngai J, Yosef N, Purdom E, Dudoit S.**
*Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics.*
**BMC Genomics** 2018;19(1):477.
DOI: [10.1186/s12864-018-4772-0](https://doi.org/10.1186/s12864-018-4772-0) · PMID [29914354](https://pubmed.ncbi.nlm.nih.gov/29914354/) · PMC [PMC6007078](https://pmc.ncbi.nlm.nih.gov/articles/PMC6007078/)
Software: Bioconductor [slingshot](https://bioconductor.org/packages/slingshot/) · Licence: CC BY 4.0

### 5 · Trajectory-based DE — tradeSeq

**Van den Berge K, Roux de Bézieux H, Street K, Saelens W, Cannoodt R, Saeys Y, Dudoit S, Clement L.**
*Trajectory-based differential expression analysis for single-cell sequencing data.*
**Nature Communications** 2020;11(1):1201.
DOI: [10.1038/s41467-020-14766-3](https://doi.org/10.1038/s41467-020-14766-3) · PMID [32139671](https://pubmed.ncbi.nlm.nih.gov/32139671/) · PMC [PMC7058077](https://pmc.ncbi.nlm.nih.gov/articles/PMC7058077/)
Software: Bioconductor [tradeSeq](https://bioconductor.org/packages/tradeSeq/) · Licence: CC BY 4.0

---

## Other software used by the target paper

| Tool | Role | Version in paper |
|---|---|---|
| [STAR / STARsolo](https://github.com/alexdobin/STAR) | alignment, UMI counting | 2.7.9a |
| [Seurat](https://satijalab.org/seurat/) | analysis framework | 4.9 (pre-release) |
| R | language | 4.3 |
| [CellChat](https://github.com/jinworks/CellChat) | ligand–receptor analysis | 1.6 |
| [clusterProfiler](https://bioconductor.org/packages/clusterProfiler/) · enrichplot · DOSE | GSEA / GO | not stated |
| [ComplexHeatmap](https://bioconductor.org/packages/ComplexHeatmap/) · viridis · [MetBrewer](https://github.com/BlakeRMills/MetBrewer) | figures | not stated |
| [deepTools](https://deeptools.readthedocs.io/) | *Ntrk2* coverage tracks | not stated |
| [AnimalTFDB](https://guolab.wchscu.cn/AnimalTFDB4/) | transcription factor list | not stated |
| [LungDamage](https://github.com/WALIII/LungDamage) | MATLAB histology damage scoring | — |

Reference genome: mm39 / GRCm39.

---

## A note on the PDFs

The PDFs of these papers are **not** included in this repository. SoupX, Slingshot and tradeSeq are CC BY 4.0 and freely redistributable; Scrublet, scds and the Cell Stem Cell paper are not. All six are open to read at the PMC links above.

---

*Bibliographic metadata retrieved from PubMed.*
