# A5: developmental-gene recruitment in adult repair

Read the [biological rationale](../A5_A11_shared_component_contract/BIOLOGICAL_LOGIC.md)
and [prospective plan](PLAN.md). The owner authorized revisions and execution.

The primary compares an external Guo 99-gene signature in transitional versus
activated AT2 cells of the same injured mouse. External identity exclusions leave
57 genes, then stress/cycling exclusions leave 53. These address distinct rivals.
The original Strunz-filtered 94/51-gene variants remain descriptive.

Metadata gives 24 primary-reference candidate mice (26 for resting AT2), before
500-UMI depth filtering. The test uses days 2–21 and at least 30 cells per arm.
See [the corrected cohort audit](DATA_AUDIT.md).

- `config/strunz_test_contract.json`: fixed design.
- `tables/external_test_modules.json`: additional source-defined modules.
- `scripts/02_freeze_external_test.py`: source freeze, no counts.
- `tables/external_freeze_run.json`: provenance and hashes.

Original audit outputs remain in `tables/`; new results use `tables/test_v1/`.
