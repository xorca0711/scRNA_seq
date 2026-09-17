# Trial C13: Epcam transcript in the transitional state

T4 sanity, Areg higher in DATP-like everywhere: True.
T1, Epcam lower in DATP-like in every library: False; median drop -0.0663 against a required 0.1. Met: False.

Reading: the transcript does not explain the escape, so any surface dimming would have to be post-transcriptional; consistent with the shedding account and NOT evidence for it The depth control failed in at least one library, whose T1 is not read.

**Not measured:** surface EpCAM protein, shedding, ADAM17 activity. Those need flow cytometry.

## Epcam and context, per library

| series | library | arm | gene | det_DATP_like | det_AT2 | difference |
|---|---|---|---|---|---|---|
| GSE247505 | Expt1_2wRFPr1 | KrasG12D 2w | Epcam | 0.9656 | 0.9265 | 0.0391 |
| GSE247505 | Expt1_2wRFPr2 | KrasG12D 2w | Epcam | 0.982 | 0.8943 | 0.0877 |
| GSE247505 | Expt1_4dRFPr2 | KrasG12D 4d | Epcam | 0.6674 | 0.6456 | 0.0218 |
| GSE316244 | Expt3_Het_RFP | Areg-flox/+ | Epcam | 0.9912 | 0.9249 | 0.0663 |
| GSE316244 | Expt3_Hom_RFP | Areg-flox/flox | Epcam | 0.9896 | 0.9108 | 0.0788 |
| GSE247505 | Expt1_2wRFPr1 | KrasG12D 2w | Areg | 0.9703 | 0.7132 | 0.2571 |
| GSE247505 | Expt1_2wRFPr2 | KrasG12D 2w | Areg | 0.9741 | 0.5432 | 0.4309 |
| GSE247505 | Expt1_4dRFPr2 | KrasG12D 4d | Areg | 0.893 | 0.3143 | 0.5787 |
| GSE316244 | Expt3_Het_RFP | Areg-flox/+ | Areg | 0.9887 | 0.6945 | 0.2942 |
| GSE316244 | Expt3_Hom_RFP | Areg-flox/flox | Areg | 0.9523 | 0.3066 | 0.6457 |

## Exploratory sheddase context (mRNA is a poor proxy for ADAM17 activity)

| series | library | arm | gene | det_DATP_like | det_AT2 | difference |
|---|---|---|---|---|---|---|
| GSE247505 | Expt1_2wRFPr1 | KrasG12D 2w | Adam17 | 0.3828 | 0.3864 | -0.0036 |
| GSE247505 | Expt1_2wRFPr2 | KrasG12D 2w | Adam17 | 0.2076 | 0.2973 | -0.0897 |
| GSE247505 | Expt1_4dRFPr2 | KrasG12D 4d | Adam17 | 0.1651 | 0.1315 | 0.0336 |
| GSE316244 | Expt3_Het_RFP | Areg-flox/+ | Adam17 | 0.4363 | 0.5749 | -0.1386 |
| GSE316244 | Expt3_Hom_RFP | Areg-flox/flox | Adam17 | 0.4803 | 0.5621 | -0.0818 |
| GSE247505 | Expt1_2wRFPr1 | KrasG12D 2w | Rhbdf1 | 0.3391 | 0.2089 | 0.1302 |
| GSE247505 | Expt1_2wRFPr2 | KrasG12D 2w | Rhbdf1 | 0.2116 | 0.1672 | 0.0444 |
| GSE247505 | Expt1_4dRFPr2 | KrasG12D 4d | Rhbdf1 | 0.1442 | 0.0793 | 0.0649 |
| GSE316244 | Expt3_Het_RFP | Areg-flox/+ | Rhbdf1 | 0.2544 | 0.2112 | 0.0432 |
| GSE316244 | Expt3_Hom_RFP | Areg-flox/flox | Rhbdf1 | 0.2949 | 0.173 | 0.1219 |
| GSE247505 | Expt1_2wRFPr1 | KrasG12D 2w | Rhbdf2 | 0.4531 | 0.2997 | 0.1534 |
| GSE247505 | Expt1_2wRFPr2 | KrasG12D 2w | Rhbdf2 | 0.2355 | 0.2317 | 0.0038 |
| GSE247505 | Expt1_4dRFPr2 | KrasG12D 4d | Rhbdf2 | 0.2744 | 0.1315 | 0.1429 |
| GSE316244 | Expt3_Het_RFP | Areg-flox/+ | Rhbdf2 | 0.5026 | 0.4749 | 0.0277 |
| GSE316244 | Expt3_Hom_RFP | Areg-flox/flox | Rhbdf2 | 0.6061 | 0.4427 | 0.1634 |
| GSE247505 | Expt1_2wRFPr1 | KrasG12D 2w | Timp3 | 0.0953 | 0.4954 | -0.4001 |
| GSE247505 | Expt1_2wRFPr2 | KrasG12D 2w | Timp3 | 0.0539 | 0.3227 | -0.2688 |
| GSE247505 | Expt1_4dRFPr2 | KrasG12D 4d | Timp3 | 0.0047 | 0.1195 | -0.1148 |
| GSE316244 | Expt3_Het_RFP | Areg-flox/+ | Timp3 | 0.1482 | 0.5825 | -0.4343 |
| GSE316244 | Expt3_Hom_RFP | Areg-flox/flox | Timp3 | 0.1184 | 0.5888 | -0.4704 |
