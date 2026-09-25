# HPCS bounded metadata recovery

25 September 2026, source-informed continuation of J7. The pinned author
notebooks identify `write/combined_data.h5ad` as the Figure-2 input. GEO lists
`GSE277777_combined_data.h5ad` at approximately 7.5 GB. The full object is outside
this follow-up. Inspect only its HDF5 observation metadata through verified HTTP
byte ranges, with a hard 64-MiB cumulative response ceiling and a 15-minute
runtime ceiling. Stop if ranges are unsupported or the cap is reached. Never
read the expression matrix, fit a model, or execute the author notebooks.

The metadata inspection is an eligibility audit, not a frozen endpoint test.
Recover source labels, group/chase, state labels and trace labels if present;
record absent fields and crosswalk conflicts explicitly. Numerical reconstruction
requires a separate endpoint contract after biological units are verified.
Supplementary author notebooks may be fetched individually below 20 MiB each
to resolve hash-to-source labels; the complete Zenodo archive is excluded.
