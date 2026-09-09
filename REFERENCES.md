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

## Roadmap papers (`Thesis/`)

The ordered paper roadmap followed in [`Thesis/README.md`](Thesis/README.md);
DOIs and PMIDs verified against PubMed on 2026-09-09. Only papers with a
study note in the repository get a full citation block; the rest are listed
so the order is visible here too.

### Gate 1, paper 4: the Human Lung Cell Atlas

**Sikkema L, Ramirez-Suastegui C, Strobl DC, Gillett TE, Zappia L, Madissoon E, Markov NS, Zaragosi L-E, et al.; Lung Biological Network Consortium; Luecken MD, Theis FJ.**
*An integrated cell atlas of the lung in health and disease.*
**Nature Medicine** 2023;29(6):1563-1577.
DOI: [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) · PMID [37291214](https://pubmed.ncbi.nlm.nih.gov/37291214/) · Licence: CC BY 4.0
Data: cellxgene collection `6f6d381a-7701-4781-935c-db10d30de293`; reference model Zenodo [10.5281/zenodo.7599104](https://doi.org/10.5281/zenodo.7599104)
Code: [github.com/LungCellAtlas/HLCA](https://github.com/LungCellAtlas/HLCA) · [HLCA_reproducibility](https://github.com/LungCellAtlas/HLCA_reproducibility) · [mapping_data_to_the_HLCA](https://github.com/LungCellAtlas/mapping_data_to_the_HLCA)
Study note: [`Thesis/gate1_04_sikkema_2023_hlca/README.md`](Thesis/gate1_04_sikkema_2023_hlca/README.md)

Benchmark framework used by the HLCA:
**Luecken MD, Buttner M, Chaichoompu K, et al.** *Benchmarking atlas-level data integration in single-cell genomics.*
**Nature Methods** 2022;19(1):41-50. DOI: [10.1038/s41592-021-01336-8](https://doi.org/10.1038/s41592-021-01336-8) · PMID [34949812](https://pubmed.ncbi.nlm.nih.gov/34949812/) · PMC [PMC8748196](https://pmc.ncbi.nlm.nih.gov/articles/PMC8748196/)

### The remaining roadmap papers (no study note yet)

| Order | Gate | Paper | DOI | PMID |
|--:|---|---|---|---|
| 2 | 1 | Choi J, et al. Inflammatory signals induce AT2 cell-derived damage-associated transient progenitors that mediate alveolar regeneration. *Cell Stem Cell* 2020 | [10.1016/j.stem.2020.06.020](https://doi.org/10.1016/j.stem.2020.06.020) | 32750316 |
| 3 | 1 | Nabhan AN, et al. Single-cell Wnt signaling niches maintain stemness of alveolar type 2 cells. *Science* 2018 | [10.1126/science.aam6603](https://doi.org/10.1126/science.aam6603) | 29420258 |
| 5 | 2 | Cardoso, Lee, et al. Early fibrotic niches establish tumour-permissive microenvironments. *Nature* 2026 | [10.1038/s41586-026-10399-6](https://doi.org/10.1038/s41586-026-10399-6) | 42020743 |
| 6 | 2 | Nabhan AN, et al. Targeted alveolar regeneration with Frizzled-specific agonists. *Cell* 2023 | [10.1016/j.cell.2023.05.022](https://doi.org/10.1016/j.cell.2023.05.022) | 37321220 |
| 7 | 3A | Saxton RA, et al. Structure-based decoupling of the pro- and anti-inflammatory functions of interleukin-10. *Science* 2021 | [10.1126/science.abc8433](https://doi.org/10.1126/science.abc8433) | 33737461 |
| 8 | 3A | Saxton RA, et al. The tissue protective functions of interleukin-22 can be decoupled from pro-inflammatory actions through structure-based design. *Immunity* 2021 | [10.1016/j.immuni.2021.03.008](https://doi.org/10.1016/j.immuni.2021.03.008) | 33852830 |
| 9 | 3B | DuPage M, et al. The chromatin-modifying enzyme Ezh2 is critical for the maintenance of regulatory T cell identity after activation. *Immunity* 2015 | [10.1016/j.immuni.2015.01.007](https://doi.org/10.1016/j.immuni.2015.01.007) | 25680271 |
| 10 | 3B | Wang D, et al. Targeting EZH2 reprograms intratumoral regulatory T cells to enhance cancer immunity. *Cell Reports* 2018 | [10.1016/j.celrep.2018.05.050](https://doi.org/10.1016/j.celrep.2018.05.050) | 29898397 |
| 11 | 3B | Zhang, et al. Intratumoral Treg cell ablation elicits NK cell-mediated control of CD8 T cell-resistant tumors. *Science Immunology* 2026 | [10.1126/sciimmunol.adx4411](https://doi.org/10.1126/sciimmunol.adx4411) | 41961946 |

---

## A note on the PDFs

The PDFs of these papers are **not** included in this repository. SoupX, Slingshot and tradeSeq are CC BY 4.0 and freely redistributable; Scrublet, scds and the Cell Stem Cell paper are not. All six are open to read at the PMC links above.

---

*Bibliographic metadata retrieved from PubMed.*
