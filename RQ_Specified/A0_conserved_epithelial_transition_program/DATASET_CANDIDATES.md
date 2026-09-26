# A0 dataset candidates and selection requirements

Updated: 2026-09-25. **Extended feasibility audit complete; no primary A0 cohort frozen or scored.**
The tables below preserve the original shortlist. Current eligibility is in
[dataset_audit.csv](tables/dataset_audit.csv) and the
[current report](reports/PILOT_REPORT.md); source records and coverage have
now been checked for the named candidates and bounded alternatives.

## Current audit findings

- D1 GSE262927: only one complete author-defined AT2/Alveolar_transitional/AT1
  triplet at the 30-cell floor. GSE141259 enriched epithelium offers 10 complete
  libraries across several days. The paper supports separate mice; a proposed
  days 10–15 interval has nine complete triplets. P1 remains unfrozen.
- D2 Ke et al. resolves to GSE254356, with two enriched and one broad epithelial
  library. Pool identities and comparable replicated contrasts are unresolved.
  The GSE165063/GSE160876 author viewer has 10,918 epithelial cells and only one
  complete barcode-suffix group; mouse/hash IDs are absent. GSE149563 lacks
  recovered matching intermediate-state annotations. Neither is selected.
- V1 Haber et al. resolves to GSE92332. Its 7,216-cell atlas header has 10 batches;
  four control mice map explicitly, and two pass the proposed proximal-branch
  coverage floor. Remaining mouse mappings and transition evidence are unresolved.
- V1 skin fallback GSE67602 has 1,422 cells with author state labels in 34 capture
  batches, while the paper reports 19 mice. The mouse-to-batch map was not recovered.

## Required roles

| Role | First candidate/search direction | Evidence already available | Unresolved before selection |
|---|---|---|---|
| D1: lung repair | GSE262927, Niethamer et al., influenza injury time course | Local processed data and alveolar analysis; animal identifiers documented | Coverage of all three independently annotated states within suitable time strata; stressed controls; audit without reusing the ES1 two-marker definition as biological truth |
| D2: normal lung development | Data underlying Ke et al., *Morphogenesis and regeneration share a conserved core transition cell state program that controls lung epithelial cell fate* | Published lineage and functional evidence for a developmental transition program | Accession, downloadable matrix, biological replication, stage/batch structure and an explicit branch with both endpoints |
| V1: another epithelium | Normal intestinal differentiation; normal epidermal differentiation as an alternative | Candidate biological settings only | Identify a primary study with independently supported intermediate states, processed scRNA-seq, biological IDs and adequate paired state coverage; no accession selected |

Prefer mouse datasets in all roles because the first local repair candidate is
mouse. Exact V1 selection must precede the transfer test and use design and
coverage, not whichever tissue produces the strongest program score.

## Existing resources with known limitations

| Resource | Useful contribution | Why it cannot automatically fill a primary role |
|---|---|---|
| GSE262927 / existing ES1 external check | Local counts, animal metadata and repair context | Previously analyzed; only one animal passed ES1's specific two-marker/depth/group-floor rule. Coverage under independently supported A0 states is unknown; previous exposure must be disclosed |
| GSE145031, Choi 2020 | AT2 lineage-traced repair context and DATP reference | One library per condition, with pooled mice; sorting arms do not supply independent biological replication |
| GSE247130 and GSE310539 multiome data | Development/injury benchmarks and technical sensitivities | One pooled library per condition; no replicated age-by-injury design. The P9 material is not an independently validated developmental trajectory by itself |
| ES1 `modules.json` and source marker tables | Traceable ADI/AT2/AT1 lists, short DATP/PATS panels, pathway controls | Published marker overlap and prior inspection preclude treating these as independent discoveries or untouched validation |

The ES1 sample floor failure is a measurement warning. It neither proves that
GSE262927 is unusable under every definition nor authorizes choosing labels to
maximize A0's signal.

## Metadata to capture for each candidate

1. Primary-paper URL/DOI; accession and matrix/annotation URLs; file hashes.
2. Species, tissue, assay, disease/injury, genotype, age/stage and sampling time.
3. Animal/donor/isolation IDs, pooling, capture/library IDs and technical repeats.
4. Counts versus normalized expression, gene identifiers and available annotations.
5. Exact starting/intermediate/destination groups and source transition evidence.
6. Biological-unit counts and cells per state within the proposed comparison strata.
7. Batch/stage/condition confounding; available stressed-cell and cycling controls.
8. Prior local or investigator exposure to the data, labels and candidate effects.
9. Proposed role, eligibility decision, exclusions and unresolved questions.

Use the requirements in [PLAN.md](PLAN.md#3-dataset-eligibility). Missing
biological IDs or insufficient replication yields a descriptive resource, not
an eligible primary transfer test. If no dataset fits a role, retain that gap
and report it before expanding the search or changing the scientific question.

## Relevant local and primary sources

- [Niethamer project entry](../../Research%20Article/gate1_01_niethamer_2025/README.md)
- [Choi study and deposit notes](../../Research%20Article/gate1_02_choi_2020/README.md)
- [ES1 results and external sample coverage](../../Research%20Article/epithelial_state_specificity/results/SUMMARY.md)
- [ES1 module definitions](../../Research%20Article/epithelial_state_specificity/modules.json)
- [Ke et al. primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11945641/)
