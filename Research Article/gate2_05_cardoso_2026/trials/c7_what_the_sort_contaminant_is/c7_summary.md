# Trial C7: what the mesenchymal-sort contaminant is

T1 best match: AT1_like (rho 0.8778, margin to second 0.0056).
T2 direct scoring modal call: DATP_like {'AT2': 70, 'DATP_like': 77, 'AT1_like': 13, 'Cd177_positive': 4, 'cycling': 20}.
T3 control: cluster 11 best 0.8778 against cluster 14 best 0.4628, passes the 0.1 margin.

Reading: the leak is ordinary alveolar epithelium; the caution is milder The two routes disagree: profile correlation says AT1_like and direct scoring says DATP_like, so neither is read as settled.

## Reference clusters of the RFP-sorted epithelium

| cluster | n_cells | call | mode_fraction | confident | top_score |
|---|---|---|---|---|---|
| 0 | 2164 | DATP_like | 0.873 | True | 0.831 |
| 1 | 1615 | DATP_like | 0.972 | True | 0.8478 |
| 2 | 499 | cycling | 0.673 | True | 1.2848 |
| 3 | 2260 | AT2 | 0.881 | True | 1.5114 |
| 4 | 2305 | AT2 | 0.999 | True | 2.1828 |
| 5 | 973 | AT1_like | 0.593 | True | 0.8128 |
| 6 | 119 | cycling | 0.538 | True | -0.0018 |
| 7 | 819 | DATP_like | 0.673 | True | 0.6079 |
| 8 | 25 | AT2 | 0.8 | True | 0.7919 |
| 9 | 456 | AT2 | 1.0 | True | 2.2176 |
| 10 | 115 | AT2 | 0.93 | True | 1.534 |

## Profile correlations

| query | n_query_cells | reference_cluster | reference_call | n_reference_cells | genes_compared | rho |
|---|---|---|---|---|---|---|
| control cluster 14 | 937 | 2 | cycling | 499 | 11075 | 0.4628 |
| control cluster 14 | 937 | 6 | cycling | 119 | 8908 | 0.4567 |
| control cluster 14 | 937 | 3 | AT2 | 2260 | 10189 | 0.4499 |
| control cluster 14 | 937 | 7 | DATP_like | 819 | 9939 | 0.4474 |
| control cluster 14 | 937 | 5 | AT1_like | 973 | 9760 | 0.44 |
| control cluster 14 | 937 | 1 | DATP_like | 1615 | 9952 | 0.4385 |
| control cluster 14 | 937 | 0 | DATP_like | 2164 | 9963 | 0.4259 |
| control cluster 14 | 937 | 8 | AT2 | 25 | 10408 | 0.4116 |
| control cluster 14 | 937 | 4 | AT2 | 2305 | 9222 | 0.3766 |
| control cluster 14 | 937 | 9 | AT2 | 456 | 9217 | 0.3735 |
| control cluster 14 | 937 | 10 | AT2 | 115 | 7481 | 0.1043 |
| query cluster 11 | 184 | 5 | AT1_like | 973 | 9751 | 0.8778 |
| query cluster 11 | 184 | 7 | DATP_like | 819 | 9908 | 0.8722 |
| query cluster 11 | 184 | 0 | DATP_like | 2164 | 9891 | 0.8657 |
| query cluster 11 | 184 | 8 | AT2 | 25 | 10268 | 0.8631 |
| query cluster 11 | 184 | 1 | DATP_like | 1615 | 9973 | 0.8331 |
| query cluster 11 | 184 | 2 | cycling | 499 | 10782 | 0.8321 |
| query cluster 11 | 184 | 3 | AT2 | 2260 | 10156 | 0.827 |
| query cluster 11 | 184 | 4 | AT2 | 2305 | 9935 | 0.731 |
| query cluster 11 | 184 | 9 | AT2 | 456 | 9958 | 0.7266 |
| query cluster 11 | 184 | 6 | cycling | 119 | 10307 | 0.5642 |
| query cluster 11 | 184 | 10 | AT2 | 115 | 9587 | 0.5151 |
