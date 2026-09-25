# Trial C9: the co-expressing Fst and Runx2 fibroblast population

T1, existence in injury with both controls: FAILS.
T2, injury-induced by the both-replicates rule: True.
T4, markers from the reference arm that replicate in both bleomycin animals: 9 of 30.

Reading: the C8 observation does not replicate in an independent injury dataset and the lead closes, whatever T3 shows the population is injury-induced rather than resident Depth control failed for Bleo1_GFPp, Bleo2_GFPp, whose ratio is not read.

**Not tested here:** whether Areg deletion depletes this population. One library per genotype, and no other Areg-flox fibroblast dataset exists.

## Per library

| series | library | group | n_gated | median_genes_per_cell | det_Fst | det_Runx2 | double_positive_share | expected_if_independent | ratio | permutation_p | n_double_positive | depth_shift | depth_control_passes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| GSE316244 | Expt3_Het_niche | Areg-flox/+ | 2751 | 3637.0 | 0.3635 | 0.2185 | 0.1378 | 0.0794 | 1.7348 | 0.005 | 379 | 0.0605 | True |
| GSE316244 | Expt3_Hom_niche | Areg-flox/flox | 3700 | 3612.0 | 0.1324 | 0.06 | 0.0065 | 0.0079 | 0.8163 | 0.9453 | 24 | 0.1715 | True |
| GSE132771 | Bleo1_GFPp | bleomycin | 2443 | 1799.0 | 0.1375 | 0.0336 | 0.0131 | 0.0046 | 2.8374 | 0.005 | 32 | 1.1009 | False |
| GSE132771 | Bleo2_GFPp | bleomycin | 2056 | 1938.0 | 0.1717 | 0.0452 | 0.018 | 0.0078 | 2.3172 | 0.01 | 37 | 0.4521 | False |
| GSE132771 | UT1_GFPp | untreated | 2355 | 1563.0 | 0.0318 | 0.0025 | 0.0 | 0.0001 | 0.0 | 1.0 | 0 | 0.0 | True |
| GSE132771 | UT2_GFPp | untreated | 1913 | 1653.0 | 0.0314 | 0.0031 | 0.0 | 0.0001 | 0.0 | 1.0 | 0 | 0.0 | True |

## What else marks them in the Areg-flox/+ arm (descriptive, one library)

| gene | det_double_positive | det_double_negative | difference |
|---|---|---|---|
| Piezo2 | 0.8786 | 0.1844 | 0.6942 |
| Kif26b | 0.7916 | 0.0974 | 0.6941 |
| Ltbp2 | 0.876 | 0.2112 | 0.6647 |
| Chst11 | 0.7731 | 0.1642 | 0.6089 |
| Rnf149 | 0.6913 | 0.0994 | 0.5919 |
| Basp1 | 0.7757 | 0.1844 | 0.5913 |
| Lhfpl2 | 0.7573 | 0.1668 | 0.5905 |
| Dock5 | 0.7784 | 0.1988 | 0.5795 |
| Bmper | 0.8443 | 0.2649 | 0.5794 |
| Rftn1 | 0.847 | 0.2734 | 0.5736 |
| Elk3 | 0.81 | 0.2368 | 0.5733 |
| Sdc1 | 0.6438 | 0.0746 | 0.5692 |
| Ctsc | 0.905 | 0.3427 | 0.5623 |
| Gtf2ird1 | 0.752 | 0.1975 | 0.5545 |
| Rab31 | 0.8364 | 0.2838 | 0.5526 |
| Rai14 | 0.6939 | 0.1419 | 0.552 |
| Cemip2 | 0.723 | 0.1799 | 0.5431 |
| P4ha3 | 0.6332 | 0.0929 | 0.5404 |
| Prrx2 | 0.686 | 0.1465 | 0.5395 |
| Ror2 | 0.8153 | 0.278 | 0.5373 |

## Whether those markers replicate across the two bleomycin animals

| gene | reference_difference | replicates | diff_Bleo1_GFPp | diff_Bleo2_GFPp |
|---|---|---|---|---|
| Piezo2 | 0.6942 | True | 0.4501 | 0.464 |
| Kif26b | 0.6941 | False | 0.509 | 0.307 |
| Ltbp2 | 0.6647 | True | 0.5041 | 0.4295 |
| Chst11 | 0.6089 | False | 0.0978 | 0.0249 |
| Rnf149 | 0.5919 | True | 0.4426 | 0.4871 |
| Basp1 | 0.5913 | True | 0.7051 | 0.4152 |
| Lhfpl2 | 0.5905 | False | 0.223 | 0.1342 |
| Dock5 | 0.5795 | False | 0.2369 | 0.01 |
| Bmper | 0.5794 | False | 0.2115 | 0.4082 |
| Rftn1 | 0.5736 | False | 0.1583 | 0.25 |
| Elk3 | 0.5733 | False | 0.3594 | 0.2703 |
| Sdc1 | 0.5692 | True | 0.538 | 0.3395 |
| Ctsc | 0.5623 | False | 0.4846 | 0.2285 |
| Gtf2ird1 | 0.5545 | False | 0.221 | 0.2402 |
| Rab31 | 0.5526 | False | 0.2978 | 0.1248 |
| Rai14 | 0.552 | False | 0.0983 | 0.0556 |
| Cemip2 | 0.5431 | False |  |  |
| P4ha3 | 0.5404 | True | 0.5269 | 0.4549 |
| Prrx2 | 0.5395 | True | 0.364 | 0.4513 |
| Ror2 | 0.5373 | False | 0.2077 | 0.1118 |
