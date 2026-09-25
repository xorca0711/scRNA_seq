# Superseded power calculation

These two files are the power calculation committed with the pre-registration in
`e6b9dc9`. They are kept because a corrected input is a record, not something to
hide.

They set the effect size from the discovery run's narrower type 2 label. The
discovery run's own primary is the broad type 2 compartment label: its validation
record counts its primary contrasts as broad, and its report and the register quote
the broad figure. The corrected calculation one level up uses the broad label as
primary and keeps this narrower label as a sensitivity row.

The correction was made before any Kim score existed. Only eligibility had run,
which reads gene names and cell labels. The two labels give nearly the same
effect, so the power figures barely change:

| Discovery label | Standardized effect | Exact Wilcoxon power, 8 patients |
|---|--:|--:|
| Narrower type 2, these files | 0.840 | 0.464 |
| Narrower type 2, corrected run's sensitivity row | 0.840 | 0.466 |
| Broad compartment, corrected primary | 0.850 | 0.474 |

The narrower label's two values differ only by simulation noise. The corrected
script simulates the broad label first from the same seed, so the random draws for
the narrower label change. With 20,000 simulations the standard error of a power
near 0.47 is about 0.0035.
