# A5: developmental programme reuse in adult repair

Question-specific work for [A5](../../RESEARCH_QUESTIONS.md#a5). The canonical
hypothesis stays in the register.

**Status: cohort audited; no A5 plan or score yet.** Read the
[data audit](DATA_AUDIT.md).

## Where A5 stands

- **The developmental list exists.** The
  [shared component contract](../A5_A11_shared_component_contract/README.md) sourced
  it and froze a 94-gene development-specific module, which the owner retained.
- **Neonatal animals are no longer required.** With an outside developmental list,
  reuse can be tested inside adult injury, comparing transitional and type 2 cells
  of the same mouse.
- **A replicated cohort exists.** The Strunz 2020 time course has 26 independent
  injured mice meeting the cell floors, with full gene coverage, and the
  development-specific genes were never selected in its data.
- **The main rival is identity.** Forty-three of the 94 genes are type 1 or type 2
  identity genes, so the plan must separate reuse from type 1-directed identity.

## Next

Write A5's pre-registration, fixing the reference population, identity handling,
depth budget and time window, then fetch the 76 MB count matrix and score. None of
that has started.

## Layout

| Path | Contents |
|---|---|
| [DATA_AUDIT.md](DATA_AUDIT.md) | Cohort audit: units, coverage, circularity, rivals and open items |
| `scripts/01_audit_strunz_units.py` | The audit, from cell labels and gene names only |
| `tables/` | Per-mouse cell counts, module coverage and the run record |
| `cache/` | The two fetched metadata files; ignored, hashes recorded in the script |
