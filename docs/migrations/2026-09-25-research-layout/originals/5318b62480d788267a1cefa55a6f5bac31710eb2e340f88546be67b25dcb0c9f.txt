# Preparation report: U0 and U1

24 September 2026. These results establish input/design feasibility, not a
biological result. The owner subsequently requested a complete planning
checkpoint before the full run. No further biological execution is released.

## U0: public metadata audit completed

- Ten series, **371 sample records**, checked against each series' sample-ref
  inventory. Source URLs, retrieval times and SHA-256 hashes are in the
  per-accession JSON files and [run record](u0_geo_design_audit/run_record.json).
- [GSE300288 group counts](u0_geo_design_audit/GSE300288_group_counts.csv): early
  anti-IL-1beta has **3**, early IgG has **4**, and every other group has **4**
  deposited libraries. Biological identity/pooling remains unresolved.
- [Human pairing audit](u0_geo_design_audit/human_pairing_coverage.csv):
  GSE308103 has **23** patient labels and **22** candidate normal/precursor/LUAD
  sets. GSE307534 has **25** patient labels and **24** precursor/LUAD pairs.
  These are title-derived links, not a verified specimen crosswalk.
- GSE267228 contains anti-CD8 and IgG samples, not anti-IL-1beta samples.
- All 31 GSE300288 records contain an FFPE characteristic alongside single-cell
  library/Cell Ranger processing metadata. The contradiction is preserved in
  [the design table](u0_geo_design_audit/GSE300288_design.csv).

Initial GEO SOFT requests included large expression payloads and were stopped
at the 20 MB response cap. MINiML XML supplied compact metadata successfully.
Those failed SOFT attempts produced no retained scientific results. A separate
attempt to retrieve the Peng full-text XML from Europe PMC returned HTTP 500;
the paper's published data availability and GEO records are verified, but
the detailed animal/capture crosswalk remains a gate rather than an assumed fact.

## U1: two technical count-input checks completed

The [specification](u1_input_spec.json) was written before opening the counts.
Selection was one control and one anti-IL-1beta example at the early endpoint,
both labelled replicate 2. This **does not imply paired animals**.

| Deposited sample | Features | Barcodes | Stored entries | Total counts |
|---|---:|---:|---:|---:|
| GSM9057712, IgG example | 31,053 | 3,565 | 5,219,251 | 14,306,254 |
| GSM9057713, anti-IL-1beta example | 31,053 | 4,247 | 6,475,466 | 17,544,592 |

Both examples passed dimension, unique feature/barcode, coordinate-bound,
nonnegative-integer-count and declared-entry-count checks, with no zero-count
barcodes. Sparse coordinate duplicates were not tested. These observations
do not establish biological QC, assay provenance, valid animal replication,
cell-state identity or treatment efficacy.

Generated [input table](u1_count_input_preflight/input_qc.csv) and
[hashed run record](u1_count_input_preflight/run_record.json) retain the evidence.
The six compressed input files remain under ignored `cache/count_inputs/`.
Their combined compressed size is 44,834,880 bytes (about 42.8 MiB); this is
an observed two-library input cost, not an estimate of the complete project.

## Reproducibility of preparation only

Run with a working Python 3 interpreter; these utilities use the standard
library. Commands below are documentation, not a scheduled or automatic run.

```powershell
python Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u0_geo_design_audit.py GSE267226 GSE267228 GSE141259 GSE277777 GSE222901 GSE300288 GSE300293 GSE307534 GSE308103 GSE307529
python Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u0_design_tables.py
python Thesis/gate2_C3_yu_lee_choi_min_2026/trials/u1_count_input_preflight.py --download
```

For local verification, U0 accepts `--offline`; U1 without `--download` reads
only cached inputs. Neither utility computes a gene-expression effect or
biological significance test. The [full plan](../ANALYSIS_TRIAL_PLAN.md) names
the next stages and the owner's confirmation requirement.

## Verification at the planning checkpoint

Repository validation passed 1,284 checks; the unchanged claim contract passed
18 numeric bindings; all seven repository unit tests passed. Independent
checks matched the ten saved metadata hashes and XML sample sets to the
371-row inventory. All 115 tracked/new repository Python source files compiled.
A broad directory compile also encountered a pre-existing Python-2 example
inside the ignored portable R/Tcl tools; that bundled dependency was excluded
from the source-only compilation and was not edited. Private annotation and
large count caches are ignored by Git. These are structural/preparation checks,
not scientific reproduction of the proposed analyses.

## Version-2 planning revision

After the owner identified missing macrophage/fibroblast LR and enrichment
analyses, the Body page was re-fetched and its context distinctions mapped
into the required niche specification. This revision changed documentation,
dataset-role gates and proposed contracts only. Repository validation passed
1,297 checks; version/eligibility/authorization consistency assertions and
`git diff --check` passed. The earlier input checks were not rerun, and no
biological comparison or full acquisition was launched. Owner review remains
the boundary before execution.
