# A1: useful follow-up analysis and what it can resolve

25 September 2026. Prepared for delivery of the verified second batch. These
are proposed next packages, not executed results or relaxed eligibility rules.
The completed work and exact remaining holds are in
[SECOND_BATCH_REPORT.md](SECOND_BATCH_REPORT.md) and [JOBS.md](../JOBS.md).

Further analysis can reduce uncertainty in the measured patterns. It cannot
recover an unrecorded animal identity, measure an absent reporter, add donors
to a one-line experiment or establish chromatin causality from separate cohorts.

| Priority / package | Available evidence and proposed work | Gap it can reduce | Boundary / stopping rule |
|---|---|---|---|
| 1. HPCS robustness using cached metadata | Compare equal-source versus cell-pooled composition; omit one source label at a time; describe agreement among the deposited `cell type`, `clusterK12` and `clusterK12_stringent` labels after verifying their author-defined correspondence. Report every group, including low-cell sources. | Whether a few large sources or annotation choices dominate the descendant pattern. The 14wk Slc4a11 group already has 811/991 cells from one source. | Source omissions are sensitivity analyses, not independent-mouse validation. Related annotations are not independent classifiers. No cell-level p-values, forced label equivalence or inferred pool identities. |
| 2. Histone promoter-definition robustness | Retain the existing raw-mark/H3 decomposition and ±1/±5-kb results. Audit transcript/MANE annotations and, where unambiguous, quantify prespecified alternate TSS windows from the cached tracks. Report all 23 loci and a descriptive ranking of sensitivity. | Whether the CDKN1A lead and counterexamples depend on a gene-boundary TSS or a particular promoter window. | A post hoc robustness ranking is exploratory. No new FDR family, dominant-isoform claim, global occupancy conclusion or population inference from two preparations in one line. |
| 3. Targeted identity recovery | Inspect any remaining explicit sample sheets, NCBI SRA run/submission records or original file-name metadata not exposed in the completed GEO/ENA audit. Match identifiers only through explicit source records. For HPCS, reconcile source/pool membership and the IGO17543 age discrepancy against experimental records. | Could unlock TIGIT pairing, genotype-stratified CD44 RNA or biological-unit HPCS summaries if a genuine crosswalk is recovered. | The completed audit found no crosswalk; recovery is uncertain. Do not repeat broad searches or infer identity from column order, expression, PCA or filenames that merely resemble one another. If absent, request a short author-supplied crosswalk. |
| 4. PATS common-quantification route | First recover the executed normalization settings and mark-to-H3 association. If deposited tracks remain uninterpretable but raw libraries and controls are identifiable, estimate their actual storage/runtime and freeze a shared alignment, filtering and quantification procedure before processing. | Can replace incompatible deposited peak calls with a common within-study measurement contract. | A larger download alone does not resolve normalization. Raw processing is a new resource-bounded package; it does not remove injury/sort confounding or create independent biological replication. |
| 5. Independent regulatory replication | Screen a narrowly chosen additional cohort for matched state definitions, direct marks/accessibility, compatible controls and verified donors/preparations before acquiring payloads. Freeze the transported loci/programme and the comparison. | Whether a candidate distinction extends beyond the B2-3 induced-cell experiment. | An unrelated public assay is not automatic validation; no screened cohort is declared eligible here. Population replication requires genuinely independent units. Accessibility cannot substitute for a specific histone modification. |

## What requires missing information or new measurements

- **TIGIT:** positive/negative library → animal or explicit disjoint pool,
  with non-overlap across the four source blocks. Removing the compound alias
  does not automatically prove that the remaining aliases are independent mice.
- **CD44:** R26/OG matrix column → GSM → genotype → mouse → sort. A pooled
  eight-source description would not answer the intended genotype-specific test.
- **HPCS:** mouse/pool members and induction/harvest ages. mScarlet RNA, if
  retrieved later, would measure transcript abundance rather than the missing
  current-state protein gate; it cannot reconstruct that gate without validation.
- **PATS:** deposited-track scale/control metadata, or a verified raw-library
  route under a new resource contract.
- **Methylation and causal state distinction:** independent donors and matched
  state/perturbation/fate measurements. More overlap calculations cannot supply
  biological replication or a same-cell multimodal causal test.

The recommended next small batch is packages 1 and 2 plus one bounded pass at
the specific uninspected identity records in package 3. Freeze those scopes
before numerical work. Author contact is a separate action and has not been sent.
Whole-archive reprocessing, another omnibus embedding, and repeated testing of
the same primary RNA model do not address the recorded missing inputs.

Evidence: [HPCS source manifest](../tables/hpcs_source_composition/source_manifest.tsv),
[source/group counts](../tables/hpcs_source_composition/group_descriptive_summary.tsv),
[histone sensitivities](../tables/second_batch_verification/histone_window_and_H3_sensitivity.tsv),
[GEO/ENA audit](../tables/identity_audit/GEO_ENA_sample_audit.tsv), and the
[author HPCS code](https://github.com/dbetel/HPCS_LUAD/tree/b52d53c984e21d3bb3a163041fdb3f56b54c19c0).
