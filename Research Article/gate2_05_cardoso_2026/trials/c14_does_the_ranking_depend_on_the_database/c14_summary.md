# Trial C14: does the ranking depend on the database

Method, object, donors and compartments held constant; only the resource varies.

Reading: the domination is partly a curation property: the guard fires in 1 of 5 resources, and the ones that stay readable are the finding. AREG is first among the EGFR ligands in every resource that contains them, so that result did not depend on CellChatDB.

## Per resource

| resource | pairs | abundance_share_top15 | guard_fires | n_egfr_ligands_present | areg_rank_among_egfr_ligands | egfr_ligand_order | top_pair |
|---|---|---|---|---|---|---|---|
| cellchatdb | 313 | 0.8 | True | 5 | 1.0 | AREG, HBEGF, TGFA, EREG, BTC | COL4A2 to CD44 |
| cellphonedb | 216 | 0.4 | False | 3 | 1.0 | AREG, EREG, BTC | APP to CD74 |
| consensus | 790 | 0.267 | False | 0 |  | None | TIMP1 to CD63 |
| connectomedb2020 | 406 | 0.2 | False | 5 | 1.0 | AREG, HBEGF, TGFA, EREG, BTC | TIMP1 to CD63 |
| italk | 539 | 0.2 | False | 5 | 1.0 | AREG, HBEGF, TGFA, EREG, BTC | TIMP1 to CD63 |

## Overlap of the top 15 lists

| resource_a | resource_b | shared_of_top15 | jaccard |
|---|---|---|---|
| cellchatdb | cellphonedb | 3 | 0.111 |
| cellchatdb | connectomedb2020 | 1 | 0.034 |
| cellchatdb | consensus | 1 | 0.034 |
| cellchatdb | italk | 0 | 0.0 |
| cellphonedb | connectomedb2020 | 0 | 0.0 |
| cellphonedb | consensus | 0 | 0.0 |
| cellphonedb | italk | 0 | 0.0 |
| connectomedb2020 | consensus | 9 | 0.429 |
| connectomedb2020 | italk | 9 | 0.429 |
| consensus | italk | 11 | 0.579 |

**Not tested:** whether any pair is a real interaction. Co-expression carries no proximity, in any resource.
