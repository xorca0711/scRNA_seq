# Ligand correction, September 2026

This is a separately specified post-audit correction of C12/C14 and a post hoc
annotation/depth sensitivity for C37. It preserves every historical trial and
raw input. The original outcomes were known when `specification.json` was
written; this is neither independent preregistration nor a new validation cohort.

## What is corrected

- **C12/C14:** only deposited epithelial types can send to Fibroblast or
  Myofibroblast. A donor must supply at least 50 cells in each actual compartment;
  subtype groups still need ten cells. Every scored source/target combination
  and the winning combination for each donor/molecular pair are retained.
- **C37:** use deposited major types and the original E1 AREG cache, then stream
  the existing raw UMI matrix once for molecule totals and AREG counts. At the
  primary sensitivity budget of 1,000 molecules, a cell contributes the exact
  hypergeometric probability of retaining at least one AREG molecule. Group means
  are compared within donor; cells below the budget are excluded and each group
  must still contain at least 50 cells. Budgets 500 and 2,000 assess sensitivity.

The resource scan retains the original LIANA CellChat-like scoring and the
maximum subtype-pair score per donor. Ties now have explicit deterministic
pair/source/target ordering; the ordinal-rank rule is unchanged. Median ranks
are retained for pairs present in at least half the eligible donors. Resource
size changes the rank universe. The historical abundance guard is reported as a
composition flag only: it cannot establish or rule out actual communication.
Exact EGFR pairs and EGFR-containing complexes are reported separately.

Both analyses remain descriptive or exploratory. Source and target RNA do not
measure proximity, ligand processing, protein signalling or functional effect.
The corrected resource analysis pools IPF and control donors as the historical
analysis did. UMI-depth matching can also condition on biological RNA-content
differences; it changes the eligible cell population and does not eliminate
subtype composition, annotation uncertainty, ambient RNA or dissociation bias.

## Run

Use the repository's scientific Python environment, with numpy, pandas, scipy,
anndata, scanpy, liana and matplotlib. The run records report installed versions.
In this Windows session the bundled CPython 3.12 x64 interpreter can reuse the
existing `.venv-x64/Lib/site-packages`, despite the old venv launcher being broken:

```powershell
$env:PYTHONPATH = (Resolve-Path .venv-x64/Lib/site-packages).Path
$py = 'C:/Users/dream/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
& $py analysis/corrections/ligand/test_corrections.py
& $py analysis/corrections/ligand/c37_annotation_depth.py
& $py analysis/corrections/ligand/c12_c14_corrected.py
& $py analysis/corrections/ligand/verify_outputs.py
& $py analysis/corrections/ligand/plot_corrections.py
```

No download, cloud service or original-trial rerun is required. The resource scan
reads one donor's existing H5AD matrix at a time and checkpoints each resource
and donor. Checkpoints are reused only when code, specification, data,
environment fingerprint and scored-file hash match. Raw-derived per-cell depth
and restart checkpoints are ignored by Git; compact donor tables, source/target
score tables, run records and figures are the evidence artifacts.

The five targeted tests exercise the sender allowlist, the actual target floor,
the exact molecule-sampling expectation against scipy's hypergeometric
distribution, invalid-count rejection, and barcode-safe raw-depth extraction.
Completed-output checks independently verify all sources/targets, floors,
maximum scores, ranks and reconstructed donor-level summaries. They do not
establish the biology.

## Outputs

- `results/lr/donor_eligibility.csv`: source/target counts for every donor.
- `results/lr/*_all_source_target_scores.csv.gz`: all allowed scored combinations.
- `results/lr/*_donor_pair_winners.csv`: scores, ranks and winning source/target.
- `results/lr/resource_summary.csv` and `historical_rank_comparison.csv`:
  corrected resource findings and comparison to historical C14.
- `results/lr/egfr_component_pairs.csv`: explicit complex-aware EGFR scope.
- `results/c37/per_donor.csv` and `summary.json`: annotation and UMI sensitivities.
- `results/figures/`: PNG and SVG comparison figures.

Results and exact claim proposals are recorded in `RESULTS.md` after completion.
