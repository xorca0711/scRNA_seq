#!/usr/bin/env python
"""Regenerate the Route B outcome page from the three trials' tracked artefacts."""

from __future__ import annotations

import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "datp_epigenetics" / "trials"))
from multiome_utils import df_to_markdown  # noqa: E402

A1 = pd.read_csv(HERE / "a1_locus_co_accessibility" / "a1_pairs.csv")
A1B = pd.read_csv(HERE / "a1b_split_locus_control" / "a1b_pairs.csv")
A1C = pd.read_csv(HERE / "a1c_graded_statistic" / "a1c_pairs.csv")
R1C = pd.read_csv(HERE / "a1c_graded_statistic" / "a1c_readings.csv")

at2c = A1C[(A1C.cell_set == "AT2")]
ctrl = at2c[at2c.pair == "split-locus median"]
bio = at2c[at2c.pair == "Etv5|Abca3"]
phys = at2c[at2c.pair == "Krt8|Krt18"]
prim = at2c[at2c.pair == "Axin2|Il1r1"]
modu = at2c[at2c.pair == "WNT_MODULE|Il1r1"]
tr = A1C[(A1C.cell_set == "transitional")]

lines = [
    "# Route B: locus co-accessibility. Closed on these deposits, under a rule that was frozen first",
    "",
    "**Outcome: NOT COMPUTABLE in every well of every pass. Route B is closed on",
    "GSE310539 and GSE247130 under the decision rule trial A1c froze before it ran,",
    "and there is no fourth pass.**",
    "",
    "Route B was the one route this branch opened rather than closed. The assessment",
    "measured that at least one of Axin2's linked peaks is detected in 22.1 and 39.4",
    "per cent of cells against a transcript detected in 5.3 and 3.8 per cent, a gain",
    "of 4.1 and 10.4 fold, so the chromatin readout escapes the floor that kills the",
    "transcript readout. That was true, and it was not enough.",
    "",
    "## Three passes, three instruments, one gate",
    "",
    "| Trial | Instrument | Its positive control | Outcome |",
    "|---|---|---|---|",
    "| A1 | binary: is at least one of the locus's peaks detected | Krt8 against Krt18, 24 kb apart | **Unreadable.** The control did not clear. It was the wrong control: the vendor annotation assigns each peak to its NEAREST gene, so the two loci largely do not share peaks and the control was a biological question rather than an arithmetic identity |",
    "| A1b | the same, with a better control | a split-locus control: two random halves of ONE gene's peak set, 20 genes, median | **Not computable.** The split-locus control did not clear either, at z +0.51 to +2.05 |",
    "| A1c | graded: the FRACTION of the locus's peaks detected, correlated across cells | the same split-locus control | **Not computable**, at z +0.24 to +1.62, and Route B closed under A1c's own pre-registered rule |",
    "",
    "## What the graded statistic did show, and it is the useful part",
    "",
    "The binary form could not detect co-accessibility between anything. The graded",
    "form can:",
    "",
    df_to_markdown(pd.concat([phys, bio])[["well", "pair", "correlation", "null_mean",
                                           "null_sd", "z", "clears"]].round(5)),
    "",
    "**Etv5 against Abca3 clears at z = +3.55 in one well and Krt8 against Krt18 at",
    "z = +3.59 in another.** Those are two AT2 identity genes on different",
    "chromosomes, and two keratins 24 kilobases apart. So the graded instrument has",
    "real sensitivity to co-accessibility, both the kind that comes from shared",
    "regulation and the kind that comes from physical proximity. It is not a dead",
    "instrument; it is an instrument whose gate was set by a control that could not",
    "pass it.",
    "",
    "## The defect in the gate, disclosed rather than repaired",
    "",
    "The split-locus control was introduced in A1b as one that could not fail for",
    "biological reasons. That claim was too strong, in two ways this folder did not",
    "see until A1c had run:",
    "",
    "1. **It is judged against the wrong null.** The matched-pair null is built for",
    "   a pair with Axin2's peak count on one side and Il1r1's on the other, roughly",
    "   15 and 38. A split-locus control on a gene with 30 peaks is 15 against 15.",
    "   Its two scores are each estimated from half as many peaks as the null's",
    "   larger side, so its correlation is attenuated relative to the distribution",
    "   it is being compared with.",
    "2. **Peaks within one gene's annotation are not necessarily one regulatory",
    "   domain.** A promoter peak and a distal peak 100 kilobases away are both",
    "   linked to the same gene by the vendor annotation and need not be",
    "   co-accessible at all.",
    "",
    "So the closure of Route B is **a decision taken under a frozen rule, not a",
    "demonstration that the approach cannot work.** The rule was frozen before the",
    "trial ran and it is honoured here rather than revised, because revising a gate",
    "after it refuses is how a result gets manufactured. But the register status is",
    "Not established, never Refuted, and the next section says what a correct gate",
    "would be.",
    "",
    "## What the numbers were, with the caveat that none of them is readable",
    "",
    "Recorded because a reader is entitled to see what the trial saw, and because a",
    "later pass with a correct gate must not be able to present these as new.",
    "",
    df_to_markdown(pd.concat([prim, modu])[["well", "pair", "correlation", "null_mean",
                                            "null_sd", "z", "percentile"]].round(5)),
    "",
    "In the AT2 sets the Axin2 against Il1r1 correlation is **negative in every well**",
    "and the Wnt module against Il1r1 is near zero. In the one transitional set large",
    "enough to compute, both flip positive: Axin2 against Il1r1 at r = +0.0697",
    "(z = +1.61) and the Wnt module at r = +0.0917 (z = +2.18), on 773 cells.",
    "Neither clears, both are unreadable, and the direction of a difference between",
    "cell states is exactly the kind of thing an underpowered design invents.",
    "",
    "## What would make Route B computable, for whoever picks it up",
    "",
    "1. **A gate matched to the null it is judged against.** Build the split-locus",
    "   control from genes with about 53 peaks, split 15 against 38 rather than in",
    "   half, so the control has the same peak-count profile as the pair under test.",
    "   Restrict the split to peaks within 50 kilobases of each other so the halves",
    "   really are one regulatory domain.",
    "2. **Use the fragments files.** Both deposits carry `atac_fragments.tsv.gz`,",
    "   2.6 GB for GSE310539, which this branch has not downloaded. They allow peaks",
    "   to be re-called on the analysis cells, per-cell quality control and TSS",
    "   enrichment, none of which the filtered matrix supports.",
    "3. **Do not raise the depth budget to buy significance.** The budget was held at",
    "   the twentieth percentile through all three passes precisely so that it could",
    "   not be tuned, and it should stay frozen for a fourth.",
]

(HERE.parent / "ROUTE_B_OUTCOME.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
print("wrote ROUTE_B_OUTCOME.md")
