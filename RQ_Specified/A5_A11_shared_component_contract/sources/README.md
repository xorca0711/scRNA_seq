# Cached source files

The spreadsheet here is not tracked, by the repository-wide ignore rule for
spreadsheets. Its identity is tracked instead, in the
[specification](../config/shared_component.json).

| File | Origin | SHA-256 |
|---|---|---|
| `guo_2019_supplementary_data_2.xlsx` | Guo et al. 2019, Supplementary Data 2, archive member `41467_2018_7770_MOESM5_ESM.xlsx` | `81e6916773f91fa78669b02186a3078c24ff2ab418697ab9415ad2003e0702f8` |

To rebuild it, download the Europe PMC supplementary archive for PMC6318311 and
extract that member. The archive itself hashes to
`5942c900289b15b35ee1405067263b35037ee815eb5678fb5258985fad6afd51`. The freeze
script checks the member hash and refuses to run on a different file.
