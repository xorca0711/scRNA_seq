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
  type: independent reanalysis of public single-cell RNA-seq data (portfolio)
  question: >
    Can the published biology of lung regeneration after influenza injury be
    recovered from the raw deposited count matrices by an independently built
    pipeline — and where it cannot, why not?
  repository: https://github.com/xorca0711/scRNA_seq
  status: complete; all work merged to main (PR #1, PR #2)

datasets:
  - accession: GSE262927
    role: primary
    species: mouse
    design: H1N1 influenza injury time course, uninjured -> 366 dpi
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

stack:
  language: Python 3.12 only   # no R, no Seurat — R unavailable on this machine
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
    an earlier run (20/41) is NOT in current artefacts — do not cite it.
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
  next: gate1_02_choi_2020, gate1_03_nabhan_2018; owner decisions on PROGRESS items 12 to 21
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
  - "Portfolio framing for the UC Berkeley PI targets: no interferon or influenza narrative (owner instruction 2026-09-09); describe results by cell state, niche, macrophage and monocyte states, annotation robustness, curation hygiene."
  - "Trial scripts under Thesis/**/trials read the processed .h5ad objects row-wise (trial_utils.read_csr_rows); never load the 2.1 GB mouse object fully. Use absolute paths; the shell cwd can change between calls."

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
[`DEVELOPMENT.md`](DEVELOPMENT.md) — that file is the human-facing record.
For an AI assistant, the operative rules are: this repository's scientific
direction, validation standards, and retain/reject decisions belong to the
human project owner; AI sessions implement under those constraints. Do not
rewrite git history or remove `Co-Authored-By` trailers — AI assistance is
deliberately visible.
