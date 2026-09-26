# Reproduce the A0 scientific pilot

The current run is specified by [pilot_v1.json](config/pilot_v1.json), executed in
[tables/pilot_v1](tables/pilot_v1), and interpreted in the
[pilot results](reports/PILOT_V1_RESULTS.md). The original feasibility audit is
historical; its output tables and JSON records have not been rewritten.
Original documentation is preserved in
[the archive](history/feasibility_documentation.zip), with a
[migration record](history/documentation_migration.json).

## Inputs and environment

Use Python 3.12 with NumPy, pandas, SciPy, h5py and matplotlib. The executed
versions are recorded in `tables/pilot_v1/environment.json`. Count preparation streams the
large source export and writes integer count arrays; it does not load the full
163,236-cell developmental atlas into RAM. Allow roughly 2 GB of working memory
and several GB of disk space. Runtime depends strongly on text parsing and disk
speed; downloaded caches should be reused.

The [continuation manifest](continuation_source_manifest.json) records public
URLs, byte counts and SHA-256 hashes. Its entries include inspected but rejected
candidates, which are not analysis inputs. The definitive input subset and exact
hashes are in [preparation.json](tables/pilot_v1/preparation.json), plus the
[GEO/author identity check](tables/pilot_v1/source_identity.json).

| Input | Expected location in a reproduction checkout |
|---|---|
| Strunz cell metadata | A0 `cache/sources/GSE141259_HighResolution_cellinfo.csv.gz` |
| Strunz genes, barcodes, raw count matrix | A5 `cache/GSE141259_HighResolution_{genes.txt,barcodes.txt,rawcounts.mtx}.gz` |
| Sountoulidis author metadata, integer matrix, independent GEO raw H5 and GEO SOFT | A0 `cache/continuation_v1/` using the filenames in the run record |
| Haber full UMI matrix | A0 `cache/continuation_v1/Haber_counts.txt.gz` |
| MGI one-to-one homology source | A0 `cache/continuation_v1/MGI_MouseHuman.rpt` |
| Haber mouse mapping | Tracked `tables/v1_haber_verified_mouse_mapping.csv` |

The Strunz files are public GSE141259 HighResolution supplements. The original
[source manifest](source_manifest.json) records the metadata URL; the existing
A5 source audit/scoring scripts identify the other supplements. Obtain the
recorded bytes and verify SHA-256 before reproducing: a current download from a
mutable URL, especially MGI, may differ. Do not silently substitute updated maps.
The helper `scripts/continuation_fetch.py NAME URL --max-bytes LIMIT` downloads
bounded public sources into the continuation cache without overwriting them.

## Run in a separate reproduction checkout

The delivered tables are execution evidence. Scripts refuse to overwrite their
own output files. In a separate checkout, preserve the delivered
`tables/pilot_v1/` and `figures/pilot_v1/` directories as reference copies outside
their normal paths before running. Start with no `processed/pilot_v1/` directory.
No primary-checkout files need to be deleted. Keep the historical mouse mapping.

From the reproduction repository root, with its scientific Python active:

```powershell
$a0Scripts = 'RQ_Specified/A0_conserved_epithelial_transition_program/scripts'
python "$a0Scripts/10_verify_source_identity.py"
python "$a0Scripts/11_prepare_pilot.py" --legacy-root . --a5-source-root .
python "$a0Scripts/12_discover_pilot.py"
```

The execution used reusable inputs in other local checkouts through the two
explicit root arguments; these paths are provenance, not portable requirements.
All biological units, state choices, count floors, ortholog rules and thresholds
are fixed in the configuration before discovery. The full deposited gene set
supplies each pseudobulk library total; the shared ortholog set supplies candidates.

Read `tables/pilot_v1/frozen_programme.json` before continuing. A
`STOP_NO_QUALIFYING_COMMON_MODULE` decision prunes transfer and specificity
scoring. It is not permission to choose a smaller module, change the branches or
relax the thresholds. A valid module must be committed before any V1 programme
scores are examined. The delivered run passed discovery, froze the 50-gene module,
and ran `13_transfer_pilot.py` once. Transfer failed the mature-endpoint criterion,
so P4 was pruned under [the pre-V1 rule](TRANSFER_EXECUTION.md). Any additional
analysis needs a separately justified plan.

```powershell
python "$a0Scripts/14_verify_pilot.py"
# Commit the valid frozen_programme.json before this step.
python "$a0Scripts/13_transfer_pilot.py"
python "$a0Scripts/17_verify_transfer.py"
python "$a0Scripts/15_plot_pilot.py"
```

The numerical verifier independently recomputes all saved discovery decisions
and probes raw counts with scalar arithmetic. Portable tests run without raw
caches or scientific packages: they check saved arithmetic, biological-unit
eligibility, hashes and the preservation of the original feasibility evidence.
CI does not redownload or rerun the scientific analysis. Inspect the rendered
figures and record their current hashes after visual review, using the delivered
`figure_review.json` as the record format. Then run the portable checks:

```powershell
python -m unittest discover -s analysis/tests -p test_a0_pilot_contract.py -q
```

Two interrupted preparation attempts are documented under
`tables/preparation_attempt1_interrupted/` and
`tables/preparation_attempt2_interrupted/`; their partial counts remain ignored.
The changes were limited to sequential storage and equivalent integer parsing.
No programme effects were computed during those attempts.

## Historical audit

The old `readiness.json`, `validation.json`, `execution_record.json` and
`reports/PILOT_REPORT.md` describe the 25 September feasibility audit. They are
not current execution status. Its unnumbered audit/report scripts regenerate
those historical outputs and some reference the pre-migration `Thesis/` layout;
do not run them over this completed pilot. Commit `4ed38ff` and the documentation
archive preserve that original state. Current source-paper analyses live in
`Research Article/`; canonical links in the old plan/candidate list were repaired
without changing their scientific decisions.
