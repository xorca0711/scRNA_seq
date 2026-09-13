# AI context

A self-contained, machine-oriented operating file for AI assistants working
in this repository: status, environment, rules, and pitfalls that must not be
violated. Human-readable counterparts: [`README.md`](README.md) (overview),
[`FINDINGS.md`](FINDINGS.md) (results), [`DEVELOPMENT.md`](DEVELOPMENT.md)
(AI-assisted development disclosure and scientific ownership),
[`PROGRESS.md`](PROGRESS.md) (session state).

## Structured summary

```yaml
project:
  name: scRNA_seq
  type: analysis log of public single-cell RNA-seq reanalysis (not a portfolio)
  question: >
    Which epithelial and immune-state programmes distinguish productive lung
    repair from persistent remodelling after injury? Stage 0 asked whether the
    published biology of one injury series could be recovered from raw counts;
    Stage 1 follows the target-lab modules (Thesis/gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md).
  repository: https://github.com/xorca0711/scRNA_seq
  status: >
    Stage 0 complete (PR #1 to #4); Stage 1 follow-ups and trial S2 merged
    (PR #6 to #8), owner review pending; repository hygiene 2026-09-10
    (archive/ holds displaced material; README reframed as an analysis log)

datasets:
  - accession: GSE262927
    role: primary
    species: mouse
    design: lung injury time course (H1N1 as the injury model), uninjured -> 366 dpi
    samples: 33            # 25-sample annotated atlas + 8-sample lineage-tracing cohort
    cells_analysed: 162175 # 212,701 barcodes -> 169,807 post-QC -> 162,175 post-Scrublet
    source_paper: "Niethamer et al., Cell Stem Cell 2025, doi:10.1016/j.stem.2024.12.002"
  - accession: GSE178360
    role: secondary
    species: human
    design: healthy distal lung, 3 donors
    cells_analysed: 27729  # 36,464 -> 29,605 post-QC -> 27,729 post-Scrublet
    integration: Harmony (primary embedding)
    source_paper: "Kadur Lakshminarasimha Murthy et al., Nature 2022, doi:10.1038/s41586-022-04541-3"
  # Added 2026-09-12 for the Gate 2 paper. Downloaded, inventoried by trial C0,
  # analysed under Thesis/gate2_05_cardoso_2026/, NOT under analysis/.
  - accession: GSE316241
    role: Cardoso 2026 mesenchyme (Gate 1 target)
    species: mouse
    design: Confetti vs Red2Kras, 2 weeks post-induction, 2 libraries, 3 mice pooled each
    barcodes: 13226
  - accession: GSE316243
    role: Cardoso 2026 niche (immune plus stroma)
    species: mouse
    design: Confetti vs Red2Kras, 2 libraries, 3 mice pooled each
    barcodes: 7836
  - accession: GSE316244
    role: Cardoso 2026 Areg-flox arm (niche and RFP+ epithelium)
    species: mouse
    design: Areg-flox/+ vs Areg-flox/flox, 4 libraries, 3 mice pooled each
    barcodes: 33756
  - accession: GSE310335
    role: Cardoso 2026 human KRASG12D alveolar organoids
    species: human
    design: control vs KRASG12D, 2 libraries
    barcodes: 9408
  - accession: GSE131907
    role: "Kim et al. 2020 human LUAD; trials E1 and E1b. The first dataset in the Cardoso work whose unit permits a test"
    species: human
    design: 11 donors with paired tumour and normal lung in scope; deposited log2 TPM matrix
    barcodes: 88144
    why_it_matters: >
      eleven donors make a paired within-donor test admissible, which no Cardoso library can support;
      note the deposited annotation is tissue-exclusive, so AT2 exists only in normal lung and the
      tumour states tS1 to tS3 only in tumour lung, and a within-tumour contrast has zero pairs
  - accession: GSE136831
    role: "Adams et al. 2020 human IPF atlas; trial E2, aberrant basaloid as the transitional state"
    species: human
    design: 60 donors in scope, IPF and control, deposited cell-type grouping
    barcodes: 312928
  - accession: GSE135893
    role: "Habermann et al. 2020 human pulmonary fibrosis; trial E3, KRT5-/KRT17+ and Transitional AT2"
    species: human
    design: 22 donors in scope, IPF and control, deposited cell-type grouping
    barcodes: 114396
  - accession: GSE132771
    role: "Tsukui et al. 2020 mouse bleomycin, Col1a1-GFP sorted mesenchyme; trial E4, the injury control"
    species: mouse
    design: 2 bleomycin and 2 untreated GFP-positive libraries, one animal each
    why_it_matters: >
      the injury comparison the Cardoso paper itself used; it is what refuted the second-signal
      reading of claim C29, because Runx1 and Pdgfrb rise with injury alone
  - accession: GSE247505
    role: "England et al. 2025 (companion paper, ref 7): RFP+ and YFP+ lineage-labelled epithelium"
    species: mouse
    design: 20 libraries; 4 days, 2 weeks, 12 weeks; two replicate libraries per arm
    barcodes: 59581
    why_it_matters: >
      the only biological replication and the only time course in the reusable set, and the
      epithelial half of the paper's CellChat object; the Cardoso data-availability statement
      does not name it

stack:
  language: Python 3.12 only   # no R, no Seurat, R unavailable on this machine
  core: [scanpy, Scrublet, harmonypy==0.0.10, PAGA, diffusion pseudotime]
  environment_constraint: >
    Windows ARM64 host; numba/llvmlite/leidenalg have no ARM64 wheels and no C
    compiler exists, so everything runs on an emulated x86-64 CPython at
    .venv-x64/Scripts/python.exe. harmonypy pinned to 0.0.10 (pure Python).
  seeds: 0 everywhere (NumPy, PCA, UMAP, Leiden, Scrublet)
  lockfile: analysis/requirements.txt

method_principles:
  - raw data determines the workflow (format/species/thresholds detected, not assumed)
  - per-sample MAD-derived QC thresholds, each with a recorded rationale string
  - Scrublet per capture before merging; non-bimodal histograms fall back to
    the 10x expected-rate quantile (fallback logged per sample)
  - batch correction treated as a hypothesis test, decided per dataset
  - deposited author labels HELD OUT of all clustering/trajectory; used only
    afterwards as an answer key
  - every decision machine-logged (decisions.json); reports and
    docs/PIPELINE_AS_RUN.md are GENERATED from artefacts, never hand-written

key_results:
  cluster_purity_vs_author_labels: 0.947   # median, 107,626 labelled mouse cells
  mouse_clusters: 29                        # Leiden 0.3, no batch correction
  trajectory:
    method: PAGA + diffusion pseudotime rooted in AT2 (25-sample cohort, 5,694 alveolar cells)
    label_ordering_heldout: {AT2: 0.013, transitional: 0.179, AT1_AT2: 0.237, AT1: 0.327}  # median dpt per author label
    transitional_abundance_median_per_animal: {peak: "27.4% at 11 dpi", resolved: "0.3% at 366 dpi"}
    icap_persistence_median_per_animal: {homeostasis: "2.0%", peak: "37.5% at 25 dpi", one_year: "21.7%"}
  lineage_tracing:
    trace_call_reproduction: "100.0000% agreement with author labels over 107,626 cells"
    kit_line_icap_traced: "33-53% per animal -> supports CAP1 origin"
    cap2_lines: "uninformative, not negative (label only 2-8% of endothelium)"
  human_integration: >
    Harmony primary by explicit, logged override (--integration harmony); the
    automated rule did not recommend it. After: 4/31 clusters >75% one donor,
    donor_driven_clustering_check false. The uncorrected baseline count from
    an earlier run (20/41) is NOT in current artefacts, do not cite it.
  negative_results:
    - "Scrublet AT0 over-removal claim REFUTED by a stricter gate (3.9% vs 6.3% baseline)"
    - "marker-panel annotation contradicted by deposited labels in 3/29 clusters (0, 22, 25)"
    - "MAD upper bounds never bind; zero cells removed for excess counts/genes"

not_done:
  - ambient-RNA correction (SoupX; raw droplet matrices absent for GSE262927)
  - formal trajectory DE (tradeSeq-style)
  - cell-cycle regression
  - reading the GSE178360 .RDS objects (need R; they hold author annotations)

thesis_roadmap:
  directory: Thesis/
  index: Thesis/README.md        # order from the owner's Notion PI Target Map; also Thesis/ROADMAP.json
  rule: one folder per paper, added one at a time in roadmap order; study note + extracted JSON + pre-registered trial
  tracked: notes, JSON, small trial tables; PDFs and XLSX in Thesis/ are gitignored
  local_pdfs: "C:/Users/dream/Documents/AC_document/External Thesis/SAP_Thesis study/Gate_1-2_Universal/ (per gate); older ones directly under Thesis/"
  done:
    - gate1_01_niethamer_2025 (pointer to docs/ and analysis/GSE262927)
    - gate1_04_sikkema_2023_hlca (note, integration_benchmark.json, PIPELINE_FRAMING.md, trials S1 to S5 with run records; owner review pending; S2 result contradicts the human AT0 headline, see PROGRESS item 15)
    - gate2_05_cardoso_2026 (note, cardoso_2026_extracts.json, trials C0 to C3 with run records; entered out of order on the owner's instruction 2026-09-12; Gate 1 returned "not recovered" and stopped to characterise, see PROGRESS item 23)
  next: gate1_02_choi_2020, gate1_03_nabhan_2018; owner decisions on PROGRESS items 12 to 23
  s2_environment: .venv-x64 also holds torch 2.14.0 (CPU) and scvi-tools 1.5.0.post1 (frozen in trials/s2_reference_mapping/requirements_s2_env.txt); scarches package removed (incompatible with anndata 0.13); HLCA reference files under trials/s2_reference_mapping/reference/ are gitignored (embedding 2.37 GB, MD5 4aa9167707141dd884ff0202b3ab1205)

pitfalls_for_ai_assistants:
  - "raw_data/ is read-only, 7.8 GB, gitignored. NEVER modify or commit it. Never grep/walk it recursively."
  - "docs/PIPELINE_AS_RUN.md and both analysis/*/README.md are GENERATED. Edit the generators (analysis/scripts/05_write_pipeline_as_run.py, 03_write_report.py) and re-run; never hand-edit."
  - "use .venv-x64/Scripts/python.exe for all Python; the native interpreter cannot import scanpy."
  - "the whole-atlas mouse object merges two incompatible experiments; condition/trace claims must come from the focused analyses (regeneration_focus/, lineage_tracing_cohort/)."
  - "mouse composition reflects a MACS sort ratio (85-90% CD45-), not the lung; compare only within a compartment."
  - "docs/ tool pages (SoupX, scds, Slingshot, tradeSeq) describe the PUBLISHED method, not this pipeline; only Scrublet was used."
  - "where marker-panel annotation and deposited labels conflict, trust the deposited label."
  - "*.h5ad, processed/, sample_shards/ are gitignored and regenerable; figures and small tables are tracked."
  - "Thesis/ is tracked and link-checked; never commit a PDF or XLSX there. Add papers one at a time in Thesis/README.md order."
  - "Thresholds taken from a paper are frozen in the trial plan BEFORE the trial reads any table; the HLCA donor-entropy threshold (0.43) must be recomputed per dataset and, for the mouse, within time point."
  - "scvi-tools and torch are installed in the emulated .venv-x64 only (owner-authorised 2026-09-09); never on the native ARM64 interpreter. The scarches package does not import with anndata 0.13; use scvi.model.SCANVI.load_query_data for surgery. Train in the background: about 40 s per epoch for 28k cells."
  - "HLCA label transfer confidently mislabels neutrophils as classical monocytes; uncertainty does not flag absent identities that resemble present ones."
  - "Framing (owner instruction 2026-09-09, tightened 2026-09-10): no influenza or interferon narrative; H1N1 is the injury model of one series, not the subject. Describe results by cell state, niche, macrophage and monocyte states, annotation robustness, curation hygiene. The repository is an analysis log, not a portfolio."
  - "Displaced material (archive/DISPLACED.md, 2026-09-10): the Krt8-high transitional trajectory, the human KRT8 reference-aligned panels and the portfolio PDF are established elsewhere (the owner's G-SURF submission). Do not extend them here; their artefacts and scripts stay in place and are validated. ONE EXCEPTION, by owner decision on 2026-09-13 (PR #12): the primary-marker dotplot, violin and per-cluster table for KRT8, CLDN4, KRT17 and SFN are on main under analysis/GSE178360/epithelial_subanalysis/figures/reference_aligned/. Treat that as a one-off the owner authorised, not as a general relaxation, and note those figures do not use the validated palette."
  - "Repository checks: analysis/scripts/validate_repository.py (renamed from validate_portfolio.py on 2026-09-10) and .github/workflows/repository-checks.yml. The paper's own workflow document moved to docs/WORKFLOW_Niethamer2025.md."
  - "Trial scripts under Thesis/**/trials read the processed .h5ad objects row-wise (trial_utils.read_csr_rows); never load the 2.1 GB mouse object fully. Use absolute paths; the shell cwd can change between calls."
  - "Cardoso 2026 (Thesis/gate2_05_cardoso_2026/): EVERY deposited mouse library pools three mice and each genotype contributes one library per sort, so no genotype contrast in that deposit has within-group replication. Describe directions; never compute a P value on a genotype contrast there. The exception is GSE247505 (England 2025), which has two replicate libraries per arm."
  - "Cardoso 2026: the mesenchymal and immune libraries are all a single time point (2 weeks). The paper's claim that fibroblast reprogramming precedes macrophage remodelling is imaging-only and CANNOT be tested from the deposit. Do not attempt a transcriptomic ordering there."
  - "Cardoso 2026 deposits carry three distinct gene spaces (CellRanger 3.0.2 to 8.0.0 against GRCm38, plus GRCh38 for the organoids). GSE316241 and GSE316244 share one; GSE316243 and GSE247505 share another; integration across them needs an explicit Ensembl-ID intersection (30,406 genes)."
  - "GSE316244 carries one non-gene feature, BSD (the reporter construct's selection marker). The C0 rule removes any feature whose ID is not an Ensembl gene ID and carries it per cell; BSD is detected in about 39% of RFP-sorted cells and 3 to 5% of niche cells, so it is a sort check, never a cell-type call."
  - "CellChat cannot be run in this repository (R-only, no R on this machine). Do not describe any Python ligand-receptor computation as a CellChat rerun; trial C3 re-derives the expression fact the claim rests on instead, and says what it is not."
  - "A compartment marker gate named after a cell type is NOT that cell type. Trial E1's 'neutrophil' gate turned out to be 83% deposited myeloid cells in a deposit that annotates no neutrophils (dissociation loses them), which changed which pre-named outcome the trial hit. Always crosstab a gate against deposited labels before naming a result after it."
  - "Compartment-level gating dilutes rare states. In GSE131907, AREG looked no higher in tumour than in normal epithelium at compartment level, because all epithelium pools malignant states with normal AT2, club and ciliated cells. Regroup by deposited subtype (trial E1b) before concluding anything about a transitional state."
  - "The Hbegf lead does NOT transfer across species. In mouse (trial C6) the Areg-independent Hbegf share is endothelial and mesenchymal; in three human datasets (E1, E1b, E2, E3) the ligand is myeloid-dominant. EGFR is mesenchymal in mouse and in human fibrosis but epithelial in human adenocarcinoma. Do not write a cross-species Hbegf argument from this repository's results."
  - "Claim C29's three-tier fibrotic response keeps its numbers and has lost its interpretation. Trial E4 showed Runx1, Pdgfrb, Tnc, Fst, Runx2, Hbegf and Egfr are all injury-generic in bleomycin mesenchyme with no oncogene. Never describe the Areg-independent tier as a tumour-specific second signal."
  - "Dendritic cells and monocytes carry AREG and HBEGF at or above the epithelial states in every human lung dataset here. This is established immunology (Zaiss et al. 2015, doi:10.1016/j.immuni.2015.01.020), not a finding of this repository. Cite it as a constraint on epithelium-centric readings, never as a new result."
  - "The human atlases deposit one merged MatrixMarket file each (1.0 to 2.0 GB gzipped), too large to load whole on 15.6 GB. Use Thesis/gate2_05_cardoso_2026/trials/mtx_stream.py to extract selected gene rows in one streaming pass; memory stays in the hundreds of megabytes. Cache extracted vectors as .npz so later trials cost seconds."
  - "When a trial declares a reading in its frozen rules and does not compute it, supply it from that trial's own tracked tables in a separate trial (C2b for C2, E2b for E2). Do not edit the original trial's rules after the fact."
  - "Rank correlation of whole mean expression profiles has COMPARTMENT resolution and NO state resolution. Trial C7 separated an epithelial query from a fibroblast control by 0.415 of Spearman rho and separated epithelial states from each other by 0.0056, with the top four matches spanning three different state calls. Use marker scoring to name a state; use profile correlation only to confirm a compartment, and always report the margin to second place."
  - "Never fit a Gaussian mixture to a raw two-gene mean to test for a subpopulation. Trial C8's score had 31.5% of cells detecting neither gene, so a spike at zero guarantees a two-component win regardless of biology. Trial C1's use of the same test was sound because score_genes output is centred and continuous. If a state-versus-gradient test is needed on raw values, use co-detection against independence (C1c and C8 T1) with a depth-split control."
  - "Any donor-level or animal-level correlation in this repository must carry a sequencing-depth control on BOTH variables, and the rule must be applied to control pairs as well as to the primary test. In trial E6 the only significant correlation was a control pair (epithelial TGFA against fibroblast activation, p = 0.045) whose two members both track depth at 0.507 and 0.412; the frozen rule refused to read it. The same number on AREG would have read as confirmation of the paper's axis."
  - "The paper's epithelium-to-fibroblast axis shows NO donor-level coupling in human fibrosis that 22 donors could detect (trial E6: rho 0.348, p = 0.112; rho 0.43 is the smallest the test could call). Do not describe this as evidence against the paper: pooling all epithelium dilutes the transitional state the claim is about, and the state-resolved version is blocked because that state clears the 50-cell floor in only 7 donors."
  - "The mesenchymal-sort contaminant of GSE316241 (cluster 11, claim C26) is a MIXTURE, not a state: 41.8% DATP-like and 38% AT2 of 184 cells, below the 50% modal floor. Say 'about four in ten score as the paper's signalling population', never 'the contaminant is the DATP-like state'."
  - "Claim C29's tiers are amplitudes, not populations (trial C8: Runx1 with Pdgfrb co-detected at ratio 1.099, inside the 0.80 to 1.25 independence band). The co-organisation is on the FALLING tier instead: Fst with Runx2 at 1.735 in Areg-flox/+ and 0.816 after deletion. That last observation is post hoc and rests on 24 double-positive cells in one library; label it as such."
  - "The tracked analysis/raw_data_inventory.* files describe the Stage 0 downloads (51 files) and the validator checks that count. The Cardoso downloads were added to raw_data/ afterwards and are inventoried by trial C0, not by that file. Do not re-run 01_scan_raw_data.py without also updating the validator."

reproduce:
  - python analysis/scripts/01_scan_raw_data.py
  - python analysis/scripts/run_scrna_analysis.py --dataset GSE262927
  - python analysis/scripts/run_scrna_analysis.py --dataset GSE178360 --integration harmony
  - python analysis/scripts/06_regeneration_focus.py
  - python analysis/scripts/07_lineage_tracing_cohort.py
  - python analysis/scripts/10_phase_timecourse.py   # tracked metadata only; no scanpy needed
  - python analysis/scripts/11_myeloid_focus.py      # reads the 2.1 GB object row-wise
  - python analysis/scripts/12_amac_trace_by_window.py   # tracked myeloid metadata only
  - python analysis/scripts/13_myeloid_batch_sensitivity.py   # row-wise object + Table S3 (or the tracked round table)
  - python analysis/scripts/03_write_report.py --dataset GSE262927
  - python analysis/scripts/05_write_pipeline_as_run.py
```

## Authorship

Authorship, contribution, and review-process disclosure live in
[`DEVELOPMENT.md`](DEVELOPMENT.md), that file is the human-facing record.
For an AI assistant, the operative rules are: this repository's scientific
direction, validation standards, and retain/reject decisions belong to the
human project owner; AI sessions implement under those constraints. Do not
rewrite git history or remove `Co-Authored-By` trailers, AI assistance is
deliberately visible.
