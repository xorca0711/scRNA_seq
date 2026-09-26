# Reproduce the separate repair-first exploratory analysis

This run is distinct from [pilot_v1](REPRODUCING.md). It reuses overlapping public
source data, so its agreement with that pilot is not independent replication.
Its result is in [the exploratory report](reports/EXPLORATORY_PILOT_REPORT.md).
Its historical complete numerical validation is in [exploratory_validation.json](exploratory_validation.json).
The [publication record](exploratory_publication_record.json) distinguishes original
execution from later packaging and namespace/path adaptations.

## Inputs and reproduction

Use Python 3.12 with NumPy, pandas, SciPy and matplotlib. Exact executed versions
and hashes are in [exploratory_execution_record.json](exploratory_execution_record.json).
Public downloads and processed matrices are ignored. Obtain the exact files in
[expression_source_manifest.json](expression_source_manifest.json), and metadata
in [source_manifest.json](source_manifest.json), verifying saved hashes. The
helper `scripts/fetch_expression.py NAME URL --max-bytes LIMIT` records provenance;
pass the recorded Accept header when required. Developmental API queries are in
[developmental_expression_manifest.json](developmental_expression_manifest.json).
Their response is the author's normalized SCT data export, not raw counts.

Retain the frozen configuration, primary and control-variant files unchanged.
The configuration's pre-migration `Thesis/` paths and Windows separators are
resolved by `resolve_source_path` to canonical `Research Article/` paths when
needed; source hashes remain unchanged. The control definitions are also frozen
in the configuration. The referenced MSigDB GMT stays read-only under `raw_data/`.
Do not regenerate the frozen configuration after examining outcomes.

Use a separate reproduction checkout and preserve delivered exploratory outputs
as reference copies before running: the scripts write their own outputs. They
write `exploratory_readiness.json` and `exploratory_decisions.json`, and do not
change pilot_v1 or historical feasibility status. Do not run the unnumbered
historical feasibility generators over either completed analysis.

From the repository root with scientific Python active:

```powershell
$a0Scripts = 'RQ_Specified/A0_conserved_epithelial_transition_program/scripts'
python "$a0Scripts/exploratory_repair.py"
# Review the E1 gate in EXPLORATORY_PLAN.md before proceeding.
python "$a0Scripts/exploratory_specificity.py"
# Review the E2 gate before proceeding.
python "$a0Scripts/fetch_developmental_transfer.py"
python "$a0Scripts/exploratory_transfer.py"
# The recorded E3 decision justified this explicit post-transfer sensitivity.
python "$a0Scripts/exploratory_transfer_check.py"
python "$a0Scripts/plot_exploratory.py"
python "$a0Scripts/write_exploratory_report.py"
# Visually review changed PNGs and update exploratory_figure_review.json.
python "$a0Scripts/verify_exploratory.py"
```

The complete original numerical run checked expression/count alignment, frozen
definitions, independent rank probes and saved effects. Publishing the unchanged
results does not repeat that expensive analysis. Portable checks independently
check saved arithmetic, sample floors, matching summaries and preservation of
every original execution artifact (including archived pre-publication bytes):

```powershell
python -m unittest discover -s analysis/tests -p test_a0_exploratory_contract.py -q
```

No cell-based biological p-values were computed. Random-gene sets, mixtures and
matching draws do not increase the number of biological units.
