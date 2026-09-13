# Trial C8: subpopulation or gradient

T1, flox/flox Runx1 with Pdgfrb: observed 0.2246 against 0.2044 expected, ratio 1.0986, consistent with independence.
T3, mixture BIC: one component 5173.42, two components 4442.64, favouring 2.
T4, depth shift of that ratio: 0.0426 against a limit of 0.25, passes.

Reading: unresolved; both numbers reported

## Co-detection against independence

| pair | n_cells | detection_a | detection_b | observed_both | expected_if_independent | ratio | genotype | tier | verdict |
|---|---|---|---|---|---|---|---|---|---|
| Runx1 with Pdgfrb | 2751 | 0.5823 | 0.4787 | 0.3333 | 0.2788 | 1.1957 | Areg-flox/+ | retained | consistent with independence |
| Fst with Runx2 | 2751 | 0.3635 | 0.2185 | 0.1378 | 0.0794 | 1.7348 | Areg-flox/+ | falling | same cells more than chance |
| Runx1 with Fst | 2751 | 0.5823 | 0.3635 | 0.2443 | 0.2117 | 1.154 | Areg-flox/+ | mixed | consistent with independence |
| Runx1 with Pdgfrb | 3700 | 0.5059 | 0.4041 | 0.2246 | 0.2044 | 1.0986 | Areg-flox/flox | retained | consistent with independence |
| Fst with Runx2 | 3700 | 0.1324 | 0.06 | 0.0065 | 0.0079 | 0.8163 | Areg-flox/flox | falling | consistent with independence |
| Runx1 with Fst | 3700 | 0.5059 | 0.1324 | 0.0549 | 0.067 | 0.8188 | Areg-flox/flox | mixed | consistent with independence |

## Depth control, flox/flox split at the median genes per cell

| pair | n_cells | detection_a | detection_b | observed_both | expected_if_independent | ratio | half | tier | median_genes_per_cell |
|---|---|---|---|---|---|---|---|---|---|
| Runx1 with Pdgfrb | 1851 | 0.3652 | 0.3425 | 0.1356 | 0.1251 | 1.084 | shallow | retained | 3035.0 |
| Fst with Runx2 | 1851 | 0.1264 | 0.0319 | 0.0027 | 0.004 | 0.6704 | shallow | falling | 3035.0 |
| Runx1 with Fst | 1851 | 0.3652 | 0.1264 | 0.0378 | 0.0462 | 0.8191 | shallow | mixed | 3035.0 |
| Runx1 with Pdgfrb | 1849 | 0.6468 | 0.4657 | 0.3137 | 0.3012 | 1.0414 | deep | retained | 4134.0 |
| Fst with Runx2 | 1849 | 0.1385 | 0.0882 | 0.0103 | 0.0122 | 0.8419 | deep | falling | 4134.0 |
| Runx1 with Fst | 1849 | 0.6468 | 0.1385 | 0.0719 | 0.0896 | 0.8032 | deep | mixed | 4134.0 |
