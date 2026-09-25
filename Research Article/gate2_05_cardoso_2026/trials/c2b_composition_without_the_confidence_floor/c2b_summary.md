# Trial C2b output: the same composition, without the confidence floor

Trial C2 scored two of its five pre-registered directions as 0.0% against 0.0%, because the
reprogrammed fibroblast and Cd177-positive calls never reached the 50% mode floor. That is
not a measurement of absence. Removing the floor, and carrying the weakest mode fraction
beside every number so the reader can see how weak each call is:

## Fibroblasts

| compartment | genotype | call | n_cells | pct_of_compartment | n_clusters | weakest_mode_fraction |
|---|---|---|---|---|---|---|
| fibroblasts | Areg-flox/+ | adventitial fibroblast | 1210 | 38.36 | 3 | 0.65 |
| fibroblasts | Areg-flox/+ | alveolar fibroblast | 1087 | 34.46 | 3 | 0.614 |
| fibroblasts | Areg-flox/+ | reprogrammed fibroblast | 857 | 27.17 | 2 | 0.375 |
| fibroblasts | Areg-flox/flox | adventitial fibroblast | 1615 | 35.35 | 3 | 0.65 |
| fibroblasts | Areg-flox/flox | alveolar fibroblast | 2034 | 44.53 | 3 | 0.614 |
| fibroblasts | Areg-flox/flox | reprogrammed fibroblast | 919 | 20.12 | 2 | 0.375 |

## RFP+ epithelium

| compartment | genotype | call | n_cells | pct_of_compartment | n_clusters | weakest_mode_fraction |
|---|---|---|---|---|---|---|
| RFP+ epithelium | Areg-flox/+ | AT1_like | 1019 | 16.61 | 2 | 0.44 |
| RFP+ epithelium | Areg-flox/+ | AT2 | 1798 | 29.31 | 3 | 0.417 |
| RFP+ epithelium | Areg-flox/+ | DATP_like | 2896 | 47.2 | 4 | 0.463 |
| RFP+ epithelium | Areg-flox/+ | cycling | 422 | 6.88 | 2 | 0.491 |
| RFP+ epithelium | Areg-flox/flox | AT1_like | 487 | 9.34 | 2 | 0.44 |
| RFP+ epithelium | Areg-flox/flox | AT2 | 3101 | 59.46 | 3 | 0.417 |
| RFP+ epithelium | Areg-flox/flox | DATP_like | 1178 | 22.59 | 4 | 0.463 |
| RFP+ epithelium | Areg-flox/flox | cycling | 449 | 8.61 | 2 | 0.491 |

## The directions C2 could not score

| direction | flox_plus_pct | flox_flox_pct | direction_met | scorable |
|---|---|---|---|---|
| reprogrammed_fibroblast_share_falls_in_hom | 27.17 | 20.12 | True | True |
| Cd177_positive_share_falls_in_hom | 0.0 | 0.0 | None | False |

## Cluster-level depletion, the stronger statement

A share computed over weakly called clusters is weaker than one cluster's own count.
Clusters most depleted in the Areg-flox/flox library:

| compartment | cluster | call | mode_fraction | n_flox_plus | n_flox_flox | ratio_flox_plus_over_flox_flox |
|---|---|---|---|---|---|---|
| RFP+ epithelium | 1 | DATP_like | 0.936 | 1712 | 111 | 15.42 |
| niche | 14 | mesothelium | 0.429 | 353 | 34 | 10.38 |
| niche | 19 | alveolar macrophage | 0.86 | 310 | 33 | 9.39 |
| niche | 13 | smooth muscle | 0.305 | 99 | 19 | 5.21 |
| RFP+ epithelium | 8 | cycling | 0.769 | 54 | 11 | 4.91 |
| RFP+ epithelium | 7 | DATP_like | 0.695 | 675 | 139 | 4.86 |
| niche | 17 | reprogrammed fibroblast | 0.384 | 210 | 45 | 4.67 |
| niche | 6 | monocyte/interstitial macrophage | 0.448 | 428 | 183 | 2.34 |

One library per genotype per sort: these are two libraries, not two groups of animals.

