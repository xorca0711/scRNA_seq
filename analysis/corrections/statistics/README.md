# Statistical corrections for C158, C161 and W1

This folder preserves a retrospective correction separately from frozen trials.
The [protocol](PROTOCOL.md), [machine-readable summary](summary.json),
[run record](run_record.json) and [comparison figure](statistical_corrections.png)
describe the analysis and its limits. `cache/` contains regenerable counts and
complete gene-set tables; `.tools/` contains a local, ignored R installation.

## Results that change the claims

- **C158:** neither DNA-replication set passes the full-family correction after
  adjustment for infection round, sex and genotype. With estimated correlation,
  the two adjusted q-values are both 0.999; with fixed correlation 0.01 they are
  0.239 and 0.627. The fixed-correlation analysis supports the unadjusted contrast
  only. The heterozygote and same-round sensitivities also fail for both settings.
  Keep the observed direction as descriptive; the inferred programme difference
  is **Not established**. Adding sacrifice day produces a rank-deficient design:
  age/injury time cannot be identified separately from phase and round.
- **C161:** zero of the 254 frozen candidates passes corrected discovery and
  validation under estimated-correlation CAMERA. This does not depend solely on
  stricter multiplicity: discovery has zero hits even with the original separate
  compartment/collection families. However, CAMERA's default fixed correlation
  0.01 supports **30** of the frozen candidates in both cohorts: one AT2, one
  fibroblast and 28 macrophage sets. **28** of these overlap the historical 62.
  The appropriate current classification is **Exploratory, method-sensitive**,
  not robust validation and not evidence of absence. The secondary setting was
  requested after inspecting the primary correction; it is not selected as a
  replacement primary analysis.
- **C163/W1:** the no-hit result survives the actual limma implementation under
  both the historical transform and voom/TMM. Estimated residual correlations
  match the historical Python calculation to about 1.2e-15, confirming input and
  design alignment. The reference analysis uses 7 and 6 residual degrees of
  freedom in A and B. No conclusion of equivalence follows.

Gene-level C162 and C164 were not rerun here. Their descriptive status and
age/sex/processing caveats remain. The corrected EMT fibroblast positive control
is not significant under the primary model (nominal p=0.0538 in discovery and
0.216 in validation); the historical positive-control observation is preserved,
but should not be presented as a newly passing corrected gate.

## Stability and composition

[LODO summaries](tables/GSE136831_lodo_summary.csv) and the
[validation LODO summary](tables/GSE135893_lodo_summary.csv) refit the primary
model once per omitted donor. Omissions leaving fewer than three donors in an arm
are explicitly flagged; they are diagnostics, not new admissible inference.
Discovery directions survive every omission for all 254 candidates. Validation
has a direction reversal for at least one omission in 54 of 244 eligible sets;
only one of the historical 62 replicated sets, macrophage HALLMARK_DNA_REPAIR,
has such a reversal. Stable direction is useful evidence but does not resolve
the correlation-sensitive significance or establish a causal state change.
The [combined candidate table](tables/g2_corrected_replication.csv) joins the
historical, primary, fixed-correlation and LODO results by compartment and set.

[Subtype eligibility](tables/GSE135893_subtype_eligibility.csv) exposes important
limits: validation proliferating macrophages retain three IPF donors and only
one control, and no validation fibroblast subtype has three controls at the
50-cell floor. Discovery subtype analyses give exploratory candidate-family hits
within both macrophage labels, but validation's sufficiently populated macrophage
label has none. Different annotations and power prevent a composition-only
conclusion. The [cell fractions](tables/subtype_cell_fractions.csv) and
[transcript contributions](tables/subtype_transcript_contributions.csv) are
descriptive decompositions, not causal mediation estimates.

The count preparation recovers the original donor/animal membership, cell totals
and compartment library totals exactly. The same ten frozen validation candidates
that were historically untested remain ineligible under the expressed-gene intersection; all 254 remain
listed, with missing statistics and a reason for the ten untested entries. The
validation multiplicity correction still includes the entire 254-candidate
family. Full discovery tests cover 10,597 eligible hypotheses, avoiding a second
discovery FDR calculation on selected hits alone.

## Reproduction

Use the repository Python launcher, or the working bundled Python 3.12 x64 with
the existing `.venv-x64/Lib/site-packages` stack. The preparation scripts support
the latter without using the broken relocated virtual-environment executable.
Install/extract R with `setup_r.ps1`; it verifies published archive checksums and
extracts locally with innoextract, without running the Windows installer. The
initial conventional install attempt rolled back because the sandbox disallowed
an uninstall registry key; portable extraction then succeeded. Package binaries
were obtained from official CRAN/Bioconductor repositories and MD5-checked by R.

From the repository root:

```text
python analysis/corrections/statistics/prepare_counts.py mouse
python analysis/corrections/statistics/prepare_counts.py GSE136831
python analysis/corrections/statistics/prepare_counts.py GSE135893
Rscript --vanilla analysis/corrections/statistics/camera_corrected.R mouse
Rscript --vanilla analysis/corrections/statistics/camera_corrected.R GSE136831
Rscript --vanilla analysis/corrections/statistics/camera_corrected.R GSE135893
Rscript --vanilla analysis/corrections/statistics/camera_corrected.R lodo
Rscript --vanilla analysis/corrections/statistics/camera_corrected.R subtypes
Rscript --vanilla analysis/corrections/statistics/camera_corrected.R fixed001
python analysis/corrections/statistics/summarize_statistics.py
python analysis/corrections/statistics/check_statistics.py
```

Here `Rscript` means `.tools/R-portable/app/bin/Rscript.exe` under this folder.
The reference versions are **R 4.6.1, limma 3.68.5 and edgeR 4.10.5**, recorded
in [R-session.txt](R-session.txt). No new remote compute or paid service is used.
Downloading current package binaries in the future can change versions; compare
the session record before interpreting a rerun as exact reproduction.

The LODO loop completed all 139 donor omissions and wrote both cohorts' outputs.
Extending the open R source with the secondary sensitivity while that loop ran
caused a trailing parser error after those writes. The final saved source was
parsed separately, session metadata refreshed, and complete donor/candidate
coverage checked. The run record records this recovery; numerical refits were
not discarded or silently described as an uninterrupted successful process.

Focused checks verify historical count identity, independent W1 residual
correlation agreement, full-discovery multiplicity, explicit ineligible
candidates, age identifiability, one omission per donor, transcript-share
conservation and provenance hashes. These are biological-unit and analysis
bookkeeping checks, not a second implementation of CAMERA.
