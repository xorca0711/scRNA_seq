# Trial E6: donor-level coupling of the AREG to EGFR axis in human fibrosis

Donors with both compartments above the 50-cell floor: 22 (IPF 18, Control 4).
Smallest |rho| a two-sided test could call here: pooled 0.428, IPF only 0.475.

## Pre-registered readings

| test | stratum | rho | p_value | depth_confounded | reading |
|---|---|---|---|---|---|
| T8 primary | pooled | 0.3484 | 0.11207 | False | no donor-level coupling detectable at this sample size; report the detectable effect size |
| T9 | pooled | -0.1496 | 0.506287 | False | no donor-level coupling detectable at this sample size; report the detectable effect size |
| T10 ligand control | pooled | 0.2434 | 0.275112 | True | both variables correlate with depth at or above the frozen line; not read |
| T10 ligand control | pooled | 0.432 | 0.044692 | True | both variables correlate with depth at or above the frozen line; not read |

## All correlations, including the three control families

| comparison | x | y | n_donors | note | rho | p_value | significant |
|---|---|---|---|---|---|---|---|
| T8 primary | epi_AREG | fib_EGFR | 22 |  | 0.3484 | 0.11207 | False |
| T9 | epi_AREG | fib_activation | 22 |  | -0.1496 | 0.506287 | False |
| T10 ligand control | epi_TGFA | fib_EGFR | 22 |  | 0.2434 | 0.275112 | False |
| T10 ligand control | epi_TGFA | fib_activation | 22 |  | 0.432 | 0.044692 | True |
| T11 depth | epi_AREG | epi_depth | 22 | epithelial depth | 0.1237 | 0.583509 | False |
| T11 depth | epi_TGFA | epi_depth | 22 | epithelial depth | 0.5065 | 0.016154 | True |
| T11 depth | fib_EGFR | fib_depth | 22 | fibroblast depth | 0.4918 | 0.020082 | True |
| T11 depth | fib_activation | fib_depth | 22 | fibroblast depth | 0.4116 | 0.056989 | False |
| T12 IPF only | epi_AREG | fib_EGFR | 18 |  | 0.2755 | 0.268425 | False |
| T12 IPF only | epi_AREG | fib_activation | 18 |  | -0.0134 | 0.957864 | False |
