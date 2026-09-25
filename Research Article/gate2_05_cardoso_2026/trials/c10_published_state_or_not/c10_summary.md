# Trial C10: the published pathological fibroblast, or not

T1, pathological score separating double-positive from double-negative: PASSES the 0.5 floor.
T4, largest effect among the four sets in the reference library: pathological.
T3, distinct population by the frozen cluster rule: False.

Reading: these are the published pathological fibroblast; trial C9's signature is that state reached by an unusual route, and the lead closes to a method note The depth control is not evaluable for Bleo1_GFPp, Bleo2_GFPp, where the shallow half holds fewer than ten double-positive cells; that is a limit of the control at this cell count, not evidence against the effect. They are a graded state rather than a distinct cluster.

**Not tested:** anything mechanical, anything spatial, and whether Areg deletion depletes the population.

## Standardised differences between the groups

| library | group | n_gated | n_double_positive | n_double_negative | d_pathological | d_alveolar | d_adventitial | d_smooth_muscle | d_pathological_shallow | d_pathological_deep | depth_control | depth_control_passes | n_double_positive_shallow | n_double_positive_deep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Expt3_Het_niche | Areg-flox/+ | 2751 | 379 | 1529 | 1.8127 | -0.756 | -0.5937 | 0.062 | 1.3451 | 1.8188 | passes | True | 37 | 342 |
| Bleo1_GFPp | bleomycin | 2443 | 32 | 2057 | 1.4926 | -1.0349 | -0.1707 | 0.2772 |  | 1.5169 | not evaluable | False | 1 | 31 |
| Bleo2_GFPp | bleomycin | 2056 | 37 | 1647 | 1.4613 | -1.2634 | -0.225 | 0.0462 |  | 1.5199 | not evaluable | False | 4 | 33 |

## Saturation audit of the pathological set

| library | gene | set | det_double_positive | det_double_negative | difference | unrankable |
|---|---|---|---|---|---|---|
| Expt3_Het_niche | Col3a1 | pathological | 1.0 | 0.9869 | 0.0131 | True |
| Expt3_Het_niche | Cthrc1 | pathological | 0.4485 | 0.0203 | 0.4283 | False |
| Expt3_Het_niche | Fn1 | pathological | 0.9763 | 0.7305 | 0.2457 | False |
| Expt3_Het_niche | Postn | pathological | 0.3061 | 0.1805 | 0.1256 | False |
| Expt3_Het_niche | Spp1 | pathological | 0.6623 | 0.191 | 0.4713 | False |
| Expt3_Het_niche | Tnc | pathological | 0.6702 | 0.2191 | 0.4511 | False |
| Bleo1_GFPp | Col3a1 | pathological | 1.0 | 0.9703 | 0.0297 | True |
| Bleo1_GFPp | Cthrc1 | pathological | 0.6875 | 0.1327 | 0.5548 | False |
| Bleo1_GFPp | Fn1 | pathological | 0.9688 | 0.7428 | 0.2259 | False |
| Bleo1_GFPp | Postn | pathological | 0.6562 | 0.1721 | 0.4842 | False |
| Bleo1_GFPp | Spp1 | pathological | 0.4688 | 0.1736 | 0.2952 | False |
| Bleo1_GFPp | Tnc | pathological | 0.75 | 0.227 | 0.523 | False |
| Bleo2_GFPp | Col3a1 | pathological | 0.973 | 0.9672 | 0.0058 | True |
| Bleo2_GFPp | Cthrc1 | pathological | 0.6216 | 0.0577 | 0.5639 | False |
| Bleo2_GFPp | Fn1 | pathological | 0.9459 | 0.7377 | 0.2082 | False |
| Bleo2_GFPp | Postn | pathological | 0.5946 | 0.1724 | 0.4222 | False |
| Bleo2_GFPp | Spp1 | pathological | 0.6757 | 0.3388 | 0.3369 | False |
| Bleo2_GFPp | Tnc | pathological | 0.7027 | 0.2617 | 0.441 | False |

## Where the double-positives sit after blind clustering of the reference library

| cluster | n_cells | n_double_positive | fraction_of_cluster | fraction_of_all_double_positive | enrichment_over_library | mean_score_pathological |
|---|---|---|---|---|---|---|
| 4 | 196 | 159 | 0.8112 | 0.4195 | 5.888 | 1.6707 |
| 9 | 256 | 80 | 0.3125 | 0.2111 | 2.268 | 0.8177 |
| 7 | 281 | 86 | 0.306 | 0.2269 | 2.221 | 0.866 |
| 10 | 126 | 11 | 0.0873 | 0.029 | 0.634 | 0.4306 |
| 0 | 407 | 24 | 0.059 | 0.0633 | 0.428 | 0.4335 |
| 1 | 262 | 8 | 0.0305 | 0.0211 | 0.222 | 0.6437 |
| 2 | 231 | 5 | 0.0216 | 0.0132 | 0.157 | 0.408 |
| 3 | 226 | 3 | 0.0133 | 0.0079 | 0.096 | 0.1773 |
| 6 | 531 | 3 | 0.0056 | 0.0079 | 0.041 | 0.3805 |
| 5 | 69 | 0 | 0.0 | 0.0 | 0.0 | -0.0989 |
| 8 | 166 | 0 | 0.0 | 0.0 | 0.0 | 0.9857 |
