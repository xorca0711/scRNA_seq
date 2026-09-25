# A5: developmental-gene recruitment in adult repair

Read the [biological rationale](../A5_A11_shared_component_contract/BIOLOGICAL_LOGIC.md)
and [prospective plan](PLAN.md). **The revised test is complete.** All 24 primary mice have positive paired
changes: mean +0.735 detection percentage points (95% CI 0.579–0.892).
The external identity- and identity/control-excluded modules remain positive
under the declared Holm family. This supports partial signature recruitment,
not shared lineage or repair function. Read the
[full results](../A5_A11_shared_component_contract/reports/REVISED_TEST_RESULTS.md).

The primary compares an external Guo 99-gene signature in transitional versus
activated AT2 cells of the same injured mouse. External identity exclusions leave
57 genes, then stress/cycling exclusions leave 53. These address distinct rivals.
The original Strunz-filtered 94/51-gene variants remain descriptive.

All 24 primary-reference mice and 26 resting-reference mice retain eligibility
after the fixed 500-UMI depth filter. The test uses days 2–21 and at least 30 cells per arm.
See [the corrected cohort audit](DATA_AUDIT.md).

- `config/strunz_test_contract.json`: fixed design.
- `tables/external_test_modules.json`: additional source-defined modules.
- `scripts/02_freeze_external_test.py`: source freeze, no counts.
- `tables/external_freeze_run.json`: provenance and hashes.

Original audit outputs remain in `tables/`; completed scores, gates and inference
are in `tables/test_v1/`.

## Reproduce the completed test

Use the repository Python launcher with a compatible scientific environment.
Run `scripts/02_freeze_external_test.py --source-root SOURCE --data-root DATA`
only in a clean output directory; the tracked frozen definitions already exist.
`scripts/03_score_external_test.py --source-root SOURCE` reads metadata from
SOURCE and the downloaded matrix/barcodes from this folder's ignored `cache/`.
Then run `Rscript scripts/04_inference.R REPO_ROOT`.
The scripts refuse to overwrite results. Source URLs/hashes are recorded in
`tables/audit_run.json`, `tables/external_freeze_run.json` and `tables/count_retrieval.json`.
