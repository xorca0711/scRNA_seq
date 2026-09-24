# AI context

A self-contained, machine-oriented operating file for AI assistants working
in this repository: status, environment, rules, and pitfalls that must not be
violated. Human-readable counterparts: [`README.md`](README.md) (overview),
[`FINDINGS.md`](FINDINGS.md) (results), [`DEVELOPMENT.md`](DEVELOPMENT.md)
(AI-assisted development disclosure and scientific ownership),
[`PROGRESS.md`](PROGRESS.md) (session state).

## Current Nabhan analysis

The owner completed Nabhan 2018 and authorized the analysis sequence on
2026-09-22. [The Nabhan analysis](Thesis/gate1_03_nabhan_2018/README.md) contains
the source FPKM reproduction, local Nb1 raw-count analysis and independent
cohort eligibility screen. Use its frozen protocol and run records; older
Stage 2 references saying Nb1 is unrun are historical. Deposited fibroblast
labels replace the old proposed cluster list. Local sampling cannot support
the acute source switch or a replicated baseline/day-11 AT2 contrast.

## Current correction authority and runtime

On 2026-09-24 the owner opened the IL-1beta review branch under
`Thesis/gate2_C3_yu_lee_choi_min_2026/` (2C item 3, stable paper 13), then
requested final self-attack, improvements and launch. This supersedes the
earlier confirmation hold; both contracts now authorize staged execution.
Read FINAL_REVIEW.md and the stage run records before continuing. Seven early
GSE300288 IgG/anti-IL-1beta libraries passed acquisition/QC (28,243 retained
cells); treatment-blind cluster diagnostics are complete. Donor-level IPF
macrophage/fibroblast/AT2 pathway fits completed in two existing cohorts:
neither has primary global q < 0.05; fixed-correlation sensitivity findings
are not robust replacements. The initial GSE136831 per-donor LIANA run
completed at 17:43 KST (25 donors: 15 IPF/10 control; 500-cell/subtype cap).
Both resources produced 25 tables each; output-integrity checks passed.
This initial compute job has ended. Differential contrasts, sensitivities,
figure rendering and interpretation remain pending. Read the paper's
INITIAL_RUN_REPORT.md and current run records before reporting status.
Scientific gates remain: validated KAC/subtype annotation, resource freeze,
coverage, and later treatment-start crosswalk. Zenodo annotated files are
restricted; no restricted-file retrieval or author contact is authorized.
The mouse two-endpoint exact permutation/Holm family cannot reject at 0.05
at the deposited n; report effect sizes and limits without changing tests.
Keep the root README as navigation and figures in paper-level galleries.

On 2026-09-22 the owner explicitly delegated orders 1–4 of the audit improvement
sequence, claim reclassification, relevant analyses and replots for portfolio
use. Decision 32 and the remediation record supersede conflicting current-state
summaries below; historical owner decisions remain preserved. New scientific
runs live under `analysis/corrections/` and `Thesis/epithelial_state_specificity/`.
Use `analysis/scripts/run_with_environment.py` with a compatible working
Python and the existing x64 site-packages when the old venv launchers fail.
See REPRODUCIBILITY.md. Do not confuse a successful structural check with
full scientific reproduction, or a post-audit specification with unseen-data
preregistration. Main and subagent costs should remain tied to decisions.

## Structured summary

```yaml
project:
  name: scRNA_seq
  type: analysis log
  purpose: >
    Hypothesis generation from an integrative reanalysis of public lung
    single-cell and multiome data, applying frameworks newer than the source
    papers to surface phenotypes and data distributions (owner restatement,
    2026-09-22).
  question: >
    Current question: which epithelial and macrophage programme changes repeat
    across independent samples after accounting for composition, genotype and
    measurement quality? Repair versus pathological remodelling requires
    independently measured outcomes and is a follow-up question. Stage 0
    asked whether the published biology of one injury series could be
    recovered from raw counts; Stage 1 follows the paper's phase and myeloid
    claims (Thesis/gate1_01_niethamer_2025/ANALYSIS_TRIAL_PLAN.md).
  repository: https://github.com/xorca0711/scRNA_seq
  status: >
    Stage 0 complete (PR #1 to #4); Stage 1 follow-ups and trials S1 to S5
    merged (PR #5 to #8); Cardoso 2026 (Gate 2, branch 2C) and the Choi 2020 deposit check merged
    (PR #10 to #24); roadmap re-ranked on
    2026-09-15 (Nabhan papers next, Gate 3 paused); Choi 2020 re-entered and
    run as trials D0 to D7 (PR #25 to #33); two branches of Choi 2020 opened
    on 2026-09-20 on multiome deposits (PR #35 to #41), their register rows
    audited by eight adversaries and corrected (PR #42); README, references
    and this file brought current on 2026-09-20 with RESEARCH_QUESTIONS.md
    added as the question-first entry point; the two original series moved
    into Thesis/ beside their source papers on 2026-09-21. Owner retain/reject
    review pending on most register rows (PROGRESS items 12 to 37)

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
  # Added 2026-09-12 for the Gate 2 paper (reading-order branch 2C since 2026-09-15). Downloaded, inventoried by trial C0,
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
  - accession: GSE145031
    role: "Choi et al. 2020 (roadmap paper 2), scRNA-seq of AT2 lineage-traced epithelium; trial D0"
    species: mouse
    design: 6 libraries, three time points (PBS, day 14, day 28 after bleomycin) by two sorts (Tomato-positive, Tomato-negative); ONE library per condition
    why_it_matters: >
      the paper that defines the DATP state the Cardoso rows lean on; its deposit has the same
      replication ceiling, and six of its matrices are raw 10x barcode whitelists rather than called cells
  - accession: GSE144468
    role: "Choi et al. 2020, scRNA-seq of AT2 organoids, control against IL-1beta; trial D0"
    species: mouse
    design: 2 libraries, already filtered to 2,101 and 3,066 called cells
  - accession: GSE144598
    role: "Choi et al. 2020, ATAC-seq of AT2 subsets"
    species: mouse
    design: bigwig coverage tracks only; NOT reusable without going to SRA for raw reads
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
  # Added 2026-09-20 for the two branches of Choi 2020. Analysed under
  # Thesis/gate1_02_choi_2020/datp_epigenetics/ and axin2_il1r1/, NOT under analysis/.
  - accession: GSE310539
    role: "Lynch et al. 2026 multiome; the transitional state in chromatin (trials M0 to M4, A1 to A1c)"
    species: mouse
    design: 10x multiome, one aggregate of 4 GEM wells (wildtype and AP-1 mutant by PBS and Sendai), 39,849 nuclei; ONE well per condition, two mice pooled
    peaks: 178178
    note: barcode-suffix map corroborated as deposited on Fos and Cldn4
  - accession: GSE247130
    role: "Hassan and Chen 2024 multiome; the same branches"
    species: mouse
    design: 10x multiome, three aggregates of 2 GEM wells (Cebpa mutant and control at P9, 7 weeks, Sendai), 64,294 nuclei; ONE well per condition
    warning: "the deposited suffix order is INVERTED: suffix 1 is the control, suffix 2 the mutant (claim C116, trial M4); multiome_utils.py applies the inversion"

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
  index: Thesis/README.md        # the owner's reading order; also Thesis/ROADMAP.json
  rule: one folder per paper, added one at a time in roadmap order; study note + extracted JSON + pre-registered trial
  tracked: notes, JSON, small trial tables; PDFs and XLSX in Thesis/ are gitignored
  local_pdfs: "C:/Users/dream/Documents/AC_document/External Thesis/SAP_Thesis study/Gate_1-2_Universal/ (per gate); older ones directly under Thesis/"
  done:
    - gate1_01_niethamer_2025 (pointer to docs/ and Thesis/gate1_01_niethamer_2025/GSE262927; follow-ups N1 to N4 run 2026-09-10; Stage 2 proposals W1, S1, Nb1, D1, V1 written, not run)
    - gate1_02_choi_2020 (re-entered 2026-09-15 at the owner's direction after reading: study note, choi_2020_extracts.json, trials D0 to D7 and the corrected passes D2b and D5b with run records; the AI-written note of 2026-09-13 was withdrawn first, decision 21, and the re-entry is decision 23; owner review pending on rows C58 to C64 and C85 to C104)
    - gate1_04_sikkema_2023_hlca (note, integration_benchmark.json, PIPELINE_FRAMING.md, trials S1 to S5 with run records; owner review pending; S2 result contradicts the human AT0 headline, see PROGRESS item 15)
    - gate2_05_cardoso_2026 (note, cardoso_2026_extracts.json, trials C0 to C12 and E1 to E6 with run records; entered out of order on the owner's instruction 2026-09-12; Gate 1 returned "not recovered"; list A exhausted; owner review pending, see PROGRESS items 23 to 29)
  next: gate1_03_nabhan_2018, then gate2_06_nabhan_2023, with proposal Nb1 as their trial (re-ranking of 2026-09-15, recorded in Thesis/README.md gate rules and ROADMAP.json); Gate 3A and 3B paused; Gate 2 carries branches 2C (Choi axis), 2N (Nabhan) and 2W (Wagner, papers 15 and 16 added 2026-09-15); methods references M1 to M8 are listed in Thesis/README.md and read at the step that uses them; owner decisions on PROGRESS items 12 to 31
  s2_environment: .venv-x64 also holds torch 2.14.0 (CPU) and scvi-tools 1.5.0.post1 (frozen in trials/s2_reference_mapping/requirements_s2_env.txt); scarches package removed (incompatible with anndata 0.13); HLCA reference files under trials/s2_reference_mapping/reference/ are gitignored (embedding 2.37 GB, MD5 4aa9167707141dd884ff0202b3ab1205)

pitfalls_for_ai_assistants:
  - "raw_data/ is read-only, 7.8 GB, gitignored. NEVER modify or commit it. Never grep/walk it recursively."
  - "docs/PIPELINE_AS_RUN.md and both series README.md files (Thesis/gate1_01_niethamer_2025/GSE262927/README.md, Thesis/ungated_murthy_2022/GSE178360/README.md) are GENERATED. Edit the generators (analysis/scripts/05_write_pipeline_as_run.py, 03_write_report.py) and re-run; never hand-edit."
  - "Layout since 2026-09-21: each deposit lives beside its source paper under Thesis/ (GSE262927 under gate1_01_niethamer_2025/, GSE178360 under ungated_murthy_2022/, whose README is a pointer note because the paper is outside the roadmap and unread). analysis/ holds only the shared pipeline, config, repository-level figures and the raw-data inventory. pipeline_utils.SERIES_DIRS is the one place the two locations are written. Run records and logs written before the move keep their analysis/GSE... paths: they are artefacts, not pointers, and are not edited."
  - "The five figures in RESEARCH_QUESTIONS.md and their caption blocks (between <!-- rq-figure:A# --> markers) are GENERATED by analysis/scripts/16_research_question_figures.py from the analysed objects; every caption number is formatted from a CSV beside the figure. They are visual aids for register rows, not evidence: nothing in them is tested, the multiome embeddings are RNA-only per well with no batch correction, and the owner asked for paper-style panels (embeddings, feature and violin panels, dotplots), not summary bars. Edit the script and re-run; never hand-edit the blocks. A1 panels g to i are the per-nucleus promoter reading that depth dominates, kept beside the detection-at-budget heatmap by owner decision 29; never cite panel g alone, because on its own it shows the retracted direction (C120)."
  - "Owner instruction 2026-09-21 on multi-agent workflows: weigh the cost before launching one. Do the deterministic part (id checks, number checks, path rewrites) with a script; deploy agents only where a script cannot judge (wording, framing, a fresh-eyes critique), at most three lenses, never one verifier per finding. Spend tokens where they change the answer."
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
  - "Framing (owner instruction 2026-09-09, tightened 2026-09-10): no influenza or interferon narrative; H1N1 is the injury model of one series, not the subject. Describe results by cell state, niche, macrophage and monocyte states, annotation robustness, curation hygiene. The repository is an analysis log for hypothesis generation."
  - "Trial W1 (2026-09-22, rows C162 to C164): in GSE262927 any contrast of 90 against 366 dpi is confounded with about nine months of age (no aged uninjured animal exists), with harvest date even inside one infection round, and inside round 2022-12-06 with sex. Never read a difference there as persistence after injury. Three animals are Ki67Cre/Cre (EEM-scRNA-125, -127, -236); W1 excluded them and removed Mki67 from every set. Proliferation gene sets co-vary across animals (inter-gene correlation 0.30 to 0.45), so a gene-permutation null overstates them by a factor of tens; use a correlation-respecting set test (CAMERA) or animal-label permutation, and say which."
  - "Purpose (owner restatement 2026-09-22, DEVELOPMENT decision 27): repository text describes one purpose only, hypothesis generation from an integrative reanalysis of public lung single-cell and multiome data. The owner's personal planning, including how the resulting questions are used, lives in private notes and never in a tracked file. Reading-order branches (2C, 2N, 2W, 3A, 3B) are named by theme, not by laboratory."
  - "Displaced material (archive/DISPLACED.md, 2026-09-10): the Krt8-high transitional trajectory, the human KRT8 reference-aligned panels and the summary PDF under archive/portfolio_2026-08/ are established outside this repository. Do not extend them here; their artefacts and scripts stay in place and are validated. ONE EXCEPTION, by owner decision on 2026-09-13 (PR #12): the primary-marker dotplot, violin and per-cluster table for KRT8, CLDN4, KRT17 and SFN are on main under Thesis/ungated_murthy_2022/GSE178360/epithelial_subanalysis/figures/reference_aligned/. Treat that as a one-off the owner authorised, not as a general relaxation. Those panels were redrawn on the validated palette on the same day and their \"AT0 candidate\" label now reads \"SFTPC+SCGB3A2+ (mostly AT2)\", matching claim C6."
  - "liana 1.10.0 is installed in .venv-x64 for trial C12 (CellChat's resource and scoring logic in Python). Installing it DOWNGRADED pandas from 3.0.5 to 2.3.3. That was checked rather than assumed: the validator and three trials were re-run and reproduced identical results, with only the recorded version string changing. Older run records therefore name pandas 3.0.5, which is correct history, not drift."
  - "CellChat itself has never been run in this repository and cannot be: it is R-only and there is no R here. Trial C3 re-derives the expression fact the paper's communication claim rests on, and trial C12 scores CellChat's own resource through liana with permutations switched off. Neither is a CellChat rerun and both say so. Never describe a Python ligand-receptor computation in this repository as CellChat."
  - "The validated figure palette lives in analysis/config/palette.json and is the ONE source of truth; Thesis/gate2_05_cardoso_2026/trials/viz_style.py reads it. Never hard-code figure colours and never reach for viridis or another default ramp. Categorical slots 1 to 4 pass the dataviz validator on the light surface, with one contrast warning that obliges visible labels or a table view; the sequential ramp is for magnitude only."
  - "Choi 2020 (Thesis/gate1_02_choi_2020/): ONE library per condition in both single-cell accessions, so no contrast there carries within-group replication either. Describe directions; never compute a P value on a between-condition comparison. The trajectory question is the exception, because an ordering inside one library needs no between-group replication."
  - "Six of the eight Choi-2020 matrices are RAW 10x barcode whitelists (737,280 columns), not called cells. Cell calling is this repository's job and needs a frozen threshold; the cell count will not match the paper's, which used Cell Ranger 2.0.2. The two organoid libraries are already filtered."
  - "The tdTomato reporter is NOT a counted feature in the Choi-2020 deposit, so the Tomato-positive and Tomato-negative split cannot be verified from the matrix; a reanalysis must trust the library labels. This is the opposite of the Cardoso deposit, where the BSD selection marker gave an independent sort check."
  - "The Choi-2020 DATP definition includes a NEGATIVE condition (low Pdpn, Hopx, Cav1). A score built from the positive genes alone will not separate DATPs from mature AT1 cells. The primed AT2 state is defined by LOSS of Etv5, Abca3 and Cebpa rather than by a positive marker, so it cannot be scored the same way as the other states. Ndrg1 sits in both the DATP marker set and the hypoxia programme, so it is not independent evidence for both."
  - "The PMC web rendering strips italicised gene symbols, which silently empties every marker set in a paper. Read full text from the Europe PMC XML instead: https://www.ebi.ac.uk/europepmc/webservices/rest/PMCID/fullTextXML. GEO accessions are often stripped too; recover them with eutils elink from the PMID to the gds database, since the GEO web pages return reCAPTCHA."
  - "Generic deposit readers (read_mtx_triplet, parse_soft, qc_metrics, ENSEMBL_ID_RE) live in Thesis/gate1_04_sikkema_2023_hlca/trials/trial_utils.py beside RunRecord, as of 2026-09-13. Thesis/gate2_05_cardoso_2026/trials/cardoso_utils.py re-exports them for the trials already written against it. Do not add a second implementation, and do not import across paper folders."
  - "NEGATIVE_RESULTS.md is GENERATED from CLAIMS.md by analysis/scripts/14_write_negative_results.py. Never edit it by hand; change the register row and re-run the script. It collects every refuted, not-established, not-establishable and retracted row, and it deliberately gives no count of how many refute this repository's own claims, because that is a judgement the register does not encode."
  - "The .venv-x64 environment was rebuilt on 2026-09-20 after its base interpreter, which uv had placed under a tmp directory, was deleted; every scanpy, leidenalg, numba, torch and liana trial stopped running and nothing recorded what had been installed. Its exact contents are now tracked in analysis/config/requirements-x64.txt with the rebuild recipe in its header. Python 3.12.13, win-amd64, 128 distributions as of 2026-09-22 (pydeseq2 0.5.4 and its six dependencies added for trial W1, owner-authorised; gseapy 1.3.1, installed for G1 and G2 on 2026-09-21, pinned late). The interpreter now lives under uv's durable default, AppData/Roaming/uv/python, not under any tmp path."
  - "To check a venv is genuinely x86-64, use sysconfig.get_platform() and expect 'win-amd64'. Do NOT use platform.machine(): on Windows it reads PROCESSOR_ARCHITECTURE and reports ARM64 even from an x86-64 interpreter running under emulation, so it will tell you the rebuild failed when it succeeded."
  - "Branch lifecycle: land work on a Claude/<topic> branch, open a pull request, merge it, then DELETE the branch both locally and on the remote. Twenty-two merged branches had accumulated on the remote by 2026-09-13 because nothing wrote this down; they were deleted on 2026-09-14 after checking that main contained every commit. A merged pull request keeps its diff on GitHub after its branch is gone, and the branch can be restored from the pull request page, so deleting is safe and reversible. Before deleting any branch, verify it with git rev-list --count main..<branch> and expect 0."
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
  - "The paper's epithelium-to-fibroblast axis did not meet the donor-level correlation significance criterion in 22 human fibrosis/control donors (trial E6: rho 0.348, p = 0.112; the approximately 0.43 value is a significance threshold, not an exclusion bound or 80%-power calculation). Do not describe this as evidence against the paper: pooling all epithelium dilutes the transitional state the claim is about, and the state-resolved version is blocked because that state clears the 50-cell floor in only 7 donors."
  - "The mesenchymal-sort contaminant of GSE316241 (cluster 11, claim C26) is a MIXTURE, not a state: 41.8% DATP-like and 38% AT2 of 184 cells, below the 50% modal floor. Say 'about four in ten score as the paper's signalling population', never 'the contaminant is the DATP-like state'."
  - "Claim C29's tiers are amplitudes, not populations (trial C8: Runx1 with Pdgfrb co-detected at ratio 1.099, inside the 0.80 to 1.25 independence band). The co-organisation is on the FALLING tier instead: Fst with Runx2 at 1.735 in Areg-flox/+ and 0.816 after deletion. That last observation is post hoc and rests on 24 double-positive cells in one library; label it as such."
  - "The tracked analysis/raw_data_inventory.* files describe the Stage 0 downloads (51 files) and the validator checks that count. The Cardoso downloads were added to raw_data/ afterwards and are inventoried by trial C0, not by that file. Do not re-run 01_scan_raw_data.py without also updating the validator."
  - "Do not write a study note for a roadmap paper the owner has not read. The owner writes or directs the note after reading (DEVELOPMENT decision 21, 2026-09-15). A deposit reality check may run when instructed, but the paper's folder is the owner's. When the owner answers 'proceed' to an agent-written list, name which items are being treated as approved before starting any that commit a reading on the owner's behalf."
  - "Two gate namespaces. Reading-order gates (1, 2C, 2N, 2W, 3A, 3B) say when a paper is read and live in Thesis/README.md and ROADMAP.json; the trial gates inside Thesis/gate2_05_cardoso_2026 (0, 1, 2a, 2b, 2d) are the owner's analysis gates for that paper alone. Never conflate them, and never rename a trial gate to match a reading-order branch."
  - "Choi 2020 trials: regenerable objects live under raw_data/GSE145031/choi_trials/ and raw_data/GSE144468/choi_trials/ (gitignored); the annotation rule and its thresholds live in Thesis/gate1_02_choi_2020/trials/choi_utils.py and were frozen before D2 ran; never move a threshold after a result, add a corrected pass beside it (D2b is that pass: every cluster of a lineage-sorted library detects Sftpc, so D2's Sftpc-low clauses never fired and an AT1 cluster was called hAT2; annotate_clusters_b applies those clauses only where Sftpc separates clusters, and D3 to D7 read its object tomato_annotated_d2b.h5ad; D5b is the same pass on the organoids, read by D6). The cell filter is the paper's own (500 to 7,000 genes, 2,000 UMI) applied to the raw whitelist; Cell Ranger 2.0.2's caller is not reproduced, so cell counts are comparable in size to the paper's, not identical."

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
  - python analysis/scripts/14_write_negative_results.py   # regenerates NEGATIVE_RESULTS.md from CLAIMS.md
  - python analysis/scripts/16_research_question_figures.py [--replot]   # the five RESEARCH_QUESTIONS.md figures and their caption blocks; needs raw_data/ and final_clustered.h5ad; --replot reuses analysis/figures/rq/processed/
  - python Thesis/gate2_05_cardoso_2026/trials/c12_cellchatdb_full_resource_scan.py   # needs liana; streams GSE136831 once and caches under raw_data/
  - python Thesis/gate1_02_choi_2020/trials/d1_cells_and_qc.py   # then d2, d2b, d3, d4, d5, d5b, d6, d7 in order; x64 venv; about 15 minutes for D1, under 5 for each of the rest
```

## Authorship

Authorship, contribution, and review-process disclosure live in
[`DEVELOPMENT.md`](DEVELOPMENT.md), that file is the human-facing record.
For an AI assistant, the operative rules are: this repository's scientific
direction, validation standards, and retain/reject decisions belong to the
human project owner; AI sessions implement under those constraints. Do not
rewrite git history or remove `Co-Authored-By` trailers, AI assistance is
deliberately visible.
