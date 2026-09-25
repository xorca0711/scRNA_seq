# Metadata evidence and unresolved identities

The 30 `GSE*.json` extracts preserve the catalog response URL, retrieval time,
response SHA-256 and sample fields. Full GEO responses are cached locally under
ignored `tmp/a1_epigenetic_metadata/`. These are public design metadata, not
downloaded sequencing matrices or a verified biological-unit manifest. No
contact details are retained. PRIDE and Zenodo extracts have their own source
and scope fields. `ENA_projects.json` records run inventories and the unresolved
crosswalk, rather than implying that runs are independent samples.

GSE190821 was added during the lineage-stage revision. Its explicit mouse IDs,
compartment and treatment fields resolve the first five-versus-five KIRA8 RNA
contrast; see the [manifest](../tables/ire1/sample_manifest.tsv). The broader
catalog is not thereby validated. Figshare 17293883 supplies published gene
signatures; exact Ensembl lookups provide stable gene-ID labels, not coordinate
liftover to GRCm38. Mapping failures and expression coverage are retained in
the [gene-set audit](../tables/ire1/gene_set_mapping_tested.tsv).

Use [the study map](../STUDY_MAP.md) for assay roles and [the plan](../PLAN.md)
for interpretation. `config/samples.json` initially contains candidate pair IDs
from titles, with biological identity explicitly unverified. Titles are useful
for finding a crosswalk; they do not certify paired animals.

## Holds that must be resolved before the affected contrast

| Resource | Observed inconsistency or missing evidence | Required resolution |
|---|---|---|
| GSE154965 | GSM4685041 is titled TIGIT-negative but has `T_plus` files; GSM4685042 is titled TIGIT-positive but has `T_minus` files | Independent public sample/file mapping; do not swap labels based on expected expression |
| GSE273122 | WT titles/filenames have SP-C Mutant characteristics; C121G titles/filenames have SP-C WT characteristics | Author sample table or corrected repository crosswalk; do not transfer this uncertainty into unrelated GSE273123 without checking |
| GSE309751 | GSM9278491 and GSM9278492 have 49-day mock/PBS titles and peak filenames but Sendai-infected treatment characteristics | Confirm controls from source metadata before any time/treatment comparison |
| GSE327686 | Series-level design describes cultured sorted AT2 bulk RNA, while GSM/library labels and matrices indicate whole-lung RNA/ATAC | Resolve tissue, time and assay mapping; do not use the copied design as authoritative |
| GSE289846 / GSE290014 | Baseline RNA and ATAC were collected at different culture days | Preserve the time mismatch; no invented same-cell or matched-time comparison |
| GSE154966 / GSE273123 | Titles supply apparent sort pairs; pool membership, disjoint biological sources and exact pairing have not been verified | Recover biological-unit evidence, especially the HPCS two-ID pooled source |
| PRJDB37980/37982/37983 | Public run listings available; exact mark/treatment/preparation mapping incomplete | Join accession-level sample metadata to supplementary design, without using article list order as a mapping |
| IMC / PXD058626 | Raw/project catalogs exist, but complete processed-cell/protein and donor-region crosswalks have not been inspected | Obtain processed matrix and biological-unit identifiers before inference |

Counts of GSMs, files, ROIs, cells or sequencing runs must not be reported as
independent sample counts. Unknowns remain null; a hold is not evidence of an
absent biological effect. Datasets below the inferential replication floor can
still produce explicitly descriptive profiles.

Run `scripts/02_preflight.py` from any directory to rebuild the local catalog
and contrast-readiness report without network access or scientific analysis.
