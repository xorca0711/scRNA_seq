# Biological research-question rewrite

25 September 2026. Implemented after the owner requested the actual rewrite of
the plan merged in PR #68. Baseline: `c3bcc695b974c1653fe7ce557796ee939f9aacf2`.

The [root register](../../../RESEARCH_QUESTIONS.md) now foregrounds biological
hypotheses, motivating observations, rivals, tests, retirement/inconclusive
criteria and readiness. All A1–A14 remain; A12-S1 is enabling, and A14 has two
separable hypotheses. There is no four-question limit. The
[measurement contracts and crosswalk](../../RQ_MEASUREMENT_CONTRACTS.md) and
[figure gallery](../../../analysis/figures/rq/README.md) hold the supporting detail.
Historical claim grades and numerical results do not change with the wording.

## Preserved sources

- [Previous register](RESEARCH_QUESTIONS.before.md.txt): exact original bytes;
  historical text, not a second current register. Relative links inside this
  snapshot refer to its original repository-root location.
- [Previous generator](16_research_question_figures.before.py.txt): exact source,
  kept as text rather than an executable way to bypass the current guard.
- [Previous A3 image](rq_a3_persistence.before.png): original title/labels retained
  as historical provenance, not the current interpretation.
- [Preservation manifest](preservation_manifest.json): SHA-256 values for these
  snapshots and 1,353 unchanged evidence files observed before editing. These are
  local byte-preservation checks for this migration, not a rule against later
  legitimate scientific updates. Do not rewrite the manifest to mask changes.

## Caption ownership and A3 correction

The root register and gallery captions are authored. Script 16 generates only
selected A1/A3/A4/A5 assets, plotted tables and panel facts in a fresh run folder.
It has no root-caption writer, missing-block insertion or active A2 renderer.
Script 18 owns the current A2. Cache-only mode fails if caches are absent/invalid;
explicit embedding rebuilds use fresh caches. Original run records remain intact.

Command used with the local scientific environment:

```powershell
python analysis/scripts/16_research_question_figures.py --figures A3 --replot
```

The [render record](../../../analysis/figures/rq/renders/20260925T064155204852Z/run_record.json)
identifies the executed source, palette, existing embedding, source tables,
versions, panel facts and output hashes. Only the A3 image was promoted to its
established path. Display labels now describe the sampling intervals and days
after infection instead of asserted repair/resolution/homeostasis. Historical
table keys are unchanged so old numerical records remain comparable.

Both newly rendered A3 tables match their originals byte-for-byte. The image
was visually checked for readable labels, complete panels and unchanged plotted
measurements. No embedding or scientific model was refit. All other figure bytes,
scientific tables, claim grades and original run records are preserved.

## Verification

Regression tests exercise the actual render entrypoints with substituted plotting
calls: an A3-only render preserves authored documents and old outputs, retired A2
cannot be selected, and invalid/missing caches never trigger deletion/rebuild.
All 14 former heading anchors remain valid alongside stable short IDs. The link
validator now recognizes explicit anchors outside code fences. Verification
passed 2,249 repository checks, 18 claim bindings, Nb1 provenance checks and
20 unit tests (one scientific integration module skipped in the lightweight
environment). Initial anchor failures and their resolution are recorded in
[verification.json](verification.json). The same record binds the promoted image
and the preserved evidence; it does not promote any biological hypothesis.
