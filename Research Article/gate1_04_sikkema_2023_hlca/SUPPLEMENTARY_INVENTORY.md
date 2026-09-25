# Supplementary material inventory for Sikkema et al. 2023

The files live outside the repository (PDF, XLSX and ZIP are gitignored; the
CSV is 11 MB and is not committed either). Sheet-to-table numbering below is
partly inferred from the gaps in the workbook (sheets 8, 9, 10, 12 and 17 are
absent from the XLSX and correspond in size and content to the ZIP and CSV
items); verify against the paper's Supplementary Information index before
citing a table number.

| Local file | Size | Paper item | Contents | Reusable for |
|---|---:|---|---|---|
| `Sikkema-2023-An-integrated-cell-atlas-of-the-lun.pdf` | 21 MB, 47 pages | main text, Methods, Extended Data Figs 1 to 10 | text extracted with pdftotext for this note | citation, parameters |
| `41591_2023_2327_MOESM1_ESM.pdf` | 5.6 MB, 12 pages | Supplementary Figures 1 to 10 | SF1 benchmark ranking (35 runs, 12 metrics); SF2 mislabel Sankeys per compartment; SF3 rare-cell markers and recall/precision for scANVI vs Harmony vs RPCA; SF4 covariate correlations and the 40-sample rule; SF5 BMI gene sets; SF6 endothelial clusters after cancer mapping; SF7 GWAS controls; SF8 deconvolution validation; SF9 spatial cell types; SF10 uncertainty ROC and per-feature uncertainty | the benchmark ranking and the calibration curve |
| `41591_2023_2327_MOESM2_ESM.pdf` | 4.4 MB | Nature Portfolio Reporting Summary | image-only; pdftotext returns nothing | none without OCR |
| `Dataset info/41591_2023_2327_MOESM3_ESM.xlsx` | 0.8 MB | Supplementary Tables 1 to 7, 11, 13 to 16 | sheet 1 dataset overview (37 studies, core/extension flag, condition, subjects, samples, cells, platform, genome build, GEO accession); sheet 2 sample overview (745 samples, 26 columns incl. donor, age, sex, condition, sampling, dissociation, harmonised ethnicity); sheet 3 CCF conversion; sheet 4 label harmonisation (288 rows, 5 levels x 22 datasets); sheet 5 manual annotations with n per level; sheet 6 marker genes (61 types, marker / marker_for / reference); sheet 7 batch-associated genes (28,527 genes, fraction of datasets affected); 11a to 11c deconvolution tests; 13a and 13b IPF signatures; 14 IPF alveolar-fibroblast cluster genes; 15a to 15g MDM cluster genes; 16 marker-selection parameters per round | dataset manifest (Tata_unpubl = GSE178360), the hierarchy, the marker sets, the batch-gene blacklist, the SPP1 MDM and IPF fibroblast signatures |
| `41591_2023_2327_MOESM4_ESM.zip` | 73 MB (170 MB unpacked) | Supplementary Table 8 | covariate-model gene results, one CSV per cell type (30 files) | sex, age, BMI, smoking, CCF coefficients per gene |
| `41591_2023_2327_MOESM5_ESM.zip` | 5.3 MB | Supplementary Table 9 | GO BP enrichment per cell type x covariate (TSV) | programme-level covariate effects |
| `41591_2023_2327_MOESM6_ESM.zip` | 0.5 MB | Supplementary Table 10 | CIBERSORTx signature matrices: airway, nose, parenchyma | bulk deconvolution |
| `41591_2023_2327_MOESM7_ESM.zip` | 1.0 MB | Supplementary Table 12 | full differential expression, high- vs low-uncertainty IPF alveolar macrophages and alveolar fibroblasts | disease-signature derivation |
| `Dataset info/41591_2023_2327_MOESM8_ESM.csv` | 11 MB | Supplementary Table 17 | gene-name harmonisation: rows = target (Ensembl 107 HGNC) names, columns = original names per dataset (37 columns) | mapping any GRCh38 build onto the HLCA gene space |

What is **not** in the local folder: the HLCA counts and embedding
(cellxgene), the reference model (Zenodo), and the code (GitHub). The
mapping trial needs the Zenodo model and the HLCA mapping tutorial.
