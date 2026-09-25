# Scripts

Both scripts refuse to overwrite their outputs. To rerun one, move its outputs
aside first and keep them; a superseded run is evidence.

| Script | Stage | Reads | Writes |
|---|---|---|---|
| `01_freeze_modules.py` | 2 | three source lists, identity and control axes; no expression data | `frozen_modules.json`, `module_membership.tsv`, `module_overlap.tsv`, `freeze_run.json` |
| `02_report_coverage.py` | 3 | gene names of four target datasets; never counts | `instrument_check.tsv`, `gene_universes.tsv`, `coverage_gate.tsv`, `coverage_gene_detail.tsv`, `unit_floor.tsv`, `coverage_run.json` |

The freeze checks the hash of the cached developmental spreadsheet and refuses a
different file. The coverage report checks that the frozen modules still match the
freeze record, then reproduces 33 coverage values from the two completed runs. It
reports nothing new unless all 33 match.

## Environment

Both run in the x86-64 scientific environment through the repository wrapper.

```powershell
$py = 'C:/Users/dream/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/python.exe'
$sp = 'X:/GitHub/scRNA_seq/.venv-x64/Lib/site-packages'
$c  = 'RQ_Specified/A5_A11_shared_component_contract/scripts'
& $py analysis/scripts/run_with_environment.py --site-packages $sp "$c/01_freeze_modules.py"
& $py analysis/scripts/run_with_environment.py --site-packages $sp "$c/02_report_coverage.py" --data-root 'X:/GitHub/scRNA_seq'
```

`--data-root` names the checkout that holds `raw_data/` and the ignored paper
caches. It defaults to the repository root, which is correct once this branch runs
in the main checkout. It was set explicitly because the work ran in a separate
worktree without those ignored directories.
