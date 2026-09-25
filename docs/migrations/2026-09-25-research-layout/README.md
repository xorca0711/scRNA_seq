# Research article and question workspace migration

This change renames `Thesis/` to `Research Article/` and establishes
`RQ_Specified/A1_transitional_epithelial_state_distinction/`. The detailed A1
specification is a subsequent change. No scientific analysis was run here.

The [manifest](manifest.json) records each relocated tracked file and each
updated existing file. `original_sha256` is the hash of its pre-migration bytes.
Unchanged files retain exactly those bytes, including scientific tables, figures,
historical run records, frozen specifications and archived scripts. Modified
active code and documentation also have an `original_bytes` archive and a
`migrated_sha256`. These are migration-time hashes, not permanent locks on future
development. Git retains the full path history.

Use `resolve_repo_path` from `analysis/lib/repository_paths.py` for current inputs
whose stored paths predate this rename. `recorded_file` additionally verifies an
expected SHA-256 and can select the explicitly archived original bytes when
checking a past run. It never treats different bytes as the same identity.
An archived script is evidence of what ran; it is not a substitute executable
for a new analysis. New runs must record the current code and current paths.

The native directory move includes ignored local caches; ignore rules retain
their exclusions. Raw data, environments, credentials and plan drafts are not
part of this PR. Existing shared question figures stay in `analysis/figures/rq/`.

Validation uses compilation, provenance unit tests, claim bindings, saved Nb1
evidence, local links, the migration manifest, and read-only saved IL-1 figure
checks. These checks verify relocation and existing evidence, not full analysis
reproduction.
