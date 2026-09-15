# Trial D6: programmes and the responder gene

## invivo

Readings: {'P1': {'p53': True, 'arrest': True, 'hypoxia': True, 'ifng_response': True}, 'P1_all': True, 'P2': {'glycolysis': True, 'Hif1a_detection': False}, 'P2_all': False, 'P3': True, 'P5': {'P1': {'p53': True, 'arrest': True, 'hypoxia': True, 'ifng_response': True}, 'P1_all': True, 'P2': {'glycolysis': True, 'Hif1a_detection': False}, 'P2_all': False, 'P3': False}, 'P4_Il1r1_by_state': {'hAT2': 0.0618, 'cAT2': 0.1182, 'DATP': 0.0614, 'AT1': 0.1039}, 'P4_G2M_claim': 'not attempted (needs an external phase gene list)'}

| state | n | score_p53 | score_arrest | score_hypoxia | score_hypoxia_without_Ndrg1 | score_ifng_response | score_glycolysis | det_Hif1a | det_Il1r1 | det_Ndrg1 | det_Pgk1 | det_Pkm | det_Slc16a3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hAT2 | 9546 | -0.1437 | -0.0855 | -0.0431 | -0.0703 | -0.0104 | -0.0181 | 0.2631 | 0.0618 | 0.0263 | 0.1266 | 0.707 | 0.001 |
| cAT2 | 313 | -0.0611 | -0.0848 | -0.0567 | -0.1125 | -0.0062 | 0.0283 | 0.476 | 0.1182 | 0.0256 | 0.3962 | 0.901 | 0.0 |
| DATP | 880 | 0.1438 | 0.846 | 0.0115 | -0.0656 | 0.2953 | 0.2526 | 0.3864 | 0.0614 | 0.1614 | 0.2693 | 0.9386 | 0.0148 |
| AT1 | 337 | -0.2633 | 0.0261 | -0.1261 | -0.182 | 0.1213 | 0.0943 | 0.1573 | 0.1039 | 0.0119 | 0.1098 | 0.7537 | 0.003 |

## organoid

Readings: {'P1': {'p53': False, 'arrest': False, 'hypoxia': True, 'ifng_response': False}, 'P1_all': False, 'P2': {'glycolysis': False, 'Hif1a_detection': False}, 'P2_all': False, 'P3': True, 'P4_Il1r1_by_state': {'hAT2': 0.0267, 'cAT2': 0.1013, 'DATP': 0.033, 'AT1': 0.1954}, 'P4_G2M_claim': 'not attempted (needs an external phase gene list)'}

| state | n | score_p53 | score_arrest | score_hypoxia | score_hypoxia_without_Ndrg1 | score_ifng_response | score_glycolysis | det_Hif1a | det_Il1r1 | det_Ndrg1 | det_Pgk1 | det_Pkm | det_Slc16a3 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| hAT2 | 935 | -0.0675 | -0.3393 | -0.0276 | 0.0147 | -0.0661 | -0.1763 | 0.6235 | 0.0267 | 0.0257 | 0.1807 | 0.8086 | 0.0032 |
| cAT2 | 79 | 0.008 | -0.2776 | 0.1398 | 0.1344 | 0.1261 | 0.005 | 0.7722 | 0.1013 | 0.2658 | 0.3418 | 0.9747 | 0.0633 |
| DATP | 2121 | -0.1479 | -0.2635 | 0.155 | 0.2089 | 0.0126 | 0.1131 | 0.7261 | 0.033 | 0.2018 | 0.2126 | 0.9486 | 0.166 |
| AT1 | 481 | -0.1581 | 0.0323 | -0.0449 | -0.0582 | 0.2206 | 0.1777 | 0.6008 | 0.1954 | 0.0769 | 0.2495 | 0.9771 | 0.1705 |

