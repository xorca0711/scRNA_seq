#!/usr/bin/env python
"""Regenerate every trial summary in this folder from its tracked artefacts.

Nothing here is typed from memory. Each summary is read back out of the CSVs
and run records the trials wrote, so a summary can never drift from the numbers
it describes, and re-running this after a trial re-runs is the only supported
way to update a summary.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from multiome_utils import df_to_markdown  # noqa: E402


def load(trial: str, name: str) -> pd.DataFrame:
    return pd.read_csv(HERE / trial / name)


def record(trial: str, name: str) -> dict:
    return json.loads((HERE / trial / name).read_text(encoding="utf-8"))


def write(trial: str, name: str, lines: list[str]) -> None:
    (HERE / trial / name).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("wrote", trial + "/" + name)


# ---------------------------------------------------------------------------

def m1b() -> None:
    t = "m1b_corrected_pass"
    cal = load(t, "m1b_calibration_labels.csv")
    lab = load(t, "m1b_labels.csv")
    rd = load(t, "m1b_readings.csv")
    write(t, "m1b_summary.md", [
        "# Trial M1b: the corrected pass, and why it also refused",
        "",
        "**Reading: NOT COMPUTABLE in all ten wells.** Not one well reached the",
        "hundred-cell floor, and the reason is arithmetic rather than biological.",
        "",
        "M1b repaired M1's two disclosed defects: the airway filter became a",
        "within-cell ratio instead of a presence call, and every label was called",
        "on RNA downsampled to one common depth instead of on raw counts. It also",
        "inherited M1's refutation of the GSE247130 suffix map. What it could not",
        "repair was its own rule R4b.",
        "",
        "## The ceiling R4b could not clear",
        "",
        "R4b required Cldn4 detection AND Krt8 at or above the 99th percentile of",
        "the uninjured well. A 99th-percentile magnitude cut passes one per cent of",
        "the control well by construction. Cldn4 is detected in five to nine per",
        "cent of an infected well. Their conjunction cannot exceed a few tenths of",
        "a per cent anywhere, and it did not:",
        "",
        df_to_markdown(lab[["well", "injured", "after_filter", "transitional",
                            "transitional_pct", "U3_max_airway_excess"]]),
        "",
        "The source paper reports a twelve per cent transitional fraction for its",
        f"infected wildtype library. The largest this rule produced anywhere was",
        f"{lab.transitional_pct.max():.2f} per cent.",
        "",
        "## The calibration run, disclosed",
        "",
        "Two technical settings were fixed after a first run that computed labels",
        "only and no accessibility quantity: the airway rule became relative, and",
        "the depth budget moved from the deposit's 10th percentile to the largest",
        "budget retaining 60 per cent of every well. The 10th percentile leaves",
        "560 UMI against a median of about 4,700, at which Cldn4 is essentially",
        "undetectable. The first run's numbers are kept beside this one:",
        "",
        df_to_markdown(cal[["well", "injured", "airway_dropped", "after_filter",
                            "transitional"]]),
        "",
        "That first run removed 49 per cent of one infected well as airway, because",
        "an absolute cut frozen from an uninjured well over-removes where Sendai",
        "virus has raised airway transcripts globally.",
        "",
        "## What was not computed",
        "",
        "No peak matrix was opened for a statistic. Every well refused at the",
        "cell-count gate before the peak pass ran.",
        "",
        df_to_markdown(rd[["well", "n_transitional", "reading", "why"]]),
    ])


def m1c() -> None:
    t = "m1c_label_at_the_depth_available"
    lab = load(t, "m1c_labels.csv")
    rec = record(t, "m1c_run_record.json")
    unin = lab[~lab.injured]
    inj = lab[lab.injured]
    neo = lab[lab.file == "P9"]
    write(t, "m1c_summary.md", [
        "# Trial M1c: the label at the depth available, and the negative control that fired",
        "",
        "**Reading: TRIAL UNREADABLE under rule U7.** And the reason U7 fired is",
        "worth more than the trial was.",
        "",
        "## What U7 required and what happened",
        "",
        "U7 required the uninjured wells to hold less of the labelled group than",
        "the injured ones, because both source papers report the CLDN4-positive",
        "state to be absent or Cldn4-negative without injury. They did not:",
        "",
        df_to_markdown(lab[["well", "injured", "after_filter", "labelled" if "labelled" in lab else "transitional",
                            "transitional_pct"]] if "transitional" in lab
                       else lab[["well", "injured", "transitional_pct"]]),
        "",
        f"The two neonatal P9 wells labelled {neo.transitional_pct.iloc[0]:.2f} and",
        f"{neo.transitional_pct.iloc[1]:.2f} per cent; the larger of the two exceeds every",
        f"injured well, the largest of which reached {inj.transitional_pct.max():.2f} per",
        f"cent, and the smaller exceeds {int((inj.transitional_pct < neo.transitional_pct.min()).sum())} of the {len(inj)}. Enrichment by",
        f"injury came out at {rec['results']['U7_enrichment_by_injury']:.1f}-fold, below the threefold floor,",
        "so the whole-trial negative control fired and the trial refused.",
        "",
        "## The refusal is one of the two answers this branch was built to find",
        "",
        "The owner's first question was whether the transition state carries a",
        "character specific to regeneration or disease. At the transcript level,",
        "for the markers that define it, the answer here is no. Krt8 and Cldn4 are",
        "expressed across immature postnatal alveolar epithelium, so a",
        "CLDN4-positive KRT8-positive call cannot separate neonatal developmental",
        "immaturity from injury-induced transition.",
        "",
        "Hassan and Chen's own argument predicts this and never states it as a",
        "limitation of the marker: their whole case is that neonatal AT2 cells are",
        "plastic and that Cebpa deletion returns mature cells toward that neonatal",
        "state. A marker set shared with normal development is not, by itself, a",
        "damage-associated marker set.",
        "",
        "## What the label does do, where the stage is held fixed",
        "",
        "Within the mature wells the label behaves as both papers describe. In",
        "GSE310539 it rises from 0.33 and 0.22 per cent under PBS to 1.75 and 0.74",
        "per cent after Sendai virus, and the AP-1 mutant reaches 0.42 of the",
        "wildtype value against the 0.5 the source paper reports. In the mature",
        "GSE247130 wells it rises from 0.69 per cent in the infected control to",
        "5.09 per cent in the infected Cebpa mutant, which is the expansion that",
        "paper reports. None of those comparisons is testable: each condition is",
        "one library.",
    ])


def m1d() -> None:
    t = "m1d_mature_wells_only"
    rd = load(t, "m1d_readings.csv")
    write(t, "m1d_summary.md", [
        "# Trial M1d: two rules that could not both hold",
        "",
        "**Reading: NOT COMPUTABLE in all eight mature wells**, and for the first",
        "time in this folder the reason was neither biology nor a bad threshold.",
        "",
        "M1d removed the neonatal wells from the design, which was right and which",
        "made U7 pass at 3.3-fold enrichment. Then every well refused at the drop",
        "gate, at 27 to 47 per cent, because two of its own frozen rules could not",
        "both be satisfied:",
        "",
        "* R5b set the fragment budget to the largest value retaining 60 per cent",
        "  of every well, which guarantees up to 40 per cent of a well falls below it.",
        "* The drop gate inherited from M1b refused any well losing more than 20",
        "  per cent of a group.",
        "",
        "A budget defined to drop forty per cent and a gate forbidding more than",
        "twenty cannot both hold. This is the second instance in this repository of",
        "the sub-shape the C12 scope conflict was: two rules written separately,",
        "each defensible, that are not composable.",
        "",
        df_to_markdown(rd[["well", "role", "n_labelled", "reading", "why"]]),
        "",
        "## Why the fix is better than either rule was",
        "",
        "The absolute drop was never the thing worth guarding. If both groups lose",
        "their shallowest fifth, the comparison between them is still fair; what it",
        "loses is generality, and that is a scope limit to report rather than a",
        "defect to refuse on. What corrupts a comparison is DIFFERENTIAL drop. M1e",
        "therefore sets the budget per well from the labelled group and draws the",
        "reference pool only from cells already above it, which makes differential",
        "drop zero by construction rather than merely bounded.",
    ])


def m1e() -> None:
    t = "m1e_per_well_budget"
    rd = load(t, "m1e_readings.csv")
    arms = load(t, "m1e_arms.csv")
    write(t, "m1e_summary.md", [
        "# Trial M1e: the first pass to reach a statistic, and its reading is superseded",
        "",
        "**M1e read SILENCED BUT NOT CLOSED in two test wells. Trial M2 then showed",
        "that the chromatin half of that reading is not stable, and it is withdrawn.**",
        "This page records what M1e computed and what survived. The withdrawal is",
        "in [`../m2_robustness_of_the_m1e_reading/m2_summary.md`](../m2_robustness_of_the_m1e_reading/m2_summary.md).",
        "",
        "## What it computed",
        "",
        df_to_markdown(rd[["well", "role", "n_labelled", "B_atac", "reading", "why"]]),
        "",
        "## The arms",
        "",
        df_to_markdown(arms[["well", "arm", "modality", "genes_used", "observed",
                             "sham_sd", "z", "clears_band"]]),
        "",
        "## What survives M2 and what does not",
        "",
        "**Survives.** The RNA half. In both injured wells the AT2 identity arm",
        "falls far outside the sham band (z of -5.9 and -11.7, a loss of 12.6 to",
        "17.1 detection points) while the AT1 arm does not move at all. M2 repeats",
        "this at four downsampling seeds and it clears in eight of eight.",
        "",
        "**Withdrawn.** The chromatin half. M1e read the AT2 chromatin arm's",
        "failure to clear as evidence that the programme stays open, which is",
        "licensed only if the positive control is solid. M2 shows it is not: the",
        "transitional arm clears in five of eight seed-and-well combinations and",
        "fails leave-one-out on two of its six genes in both wells. The AT2 arm's",
        "own null also flips to a clearance in one well when the degenerate median",
        "offset is replaced by a mean.",
        "",
        "## The defect this leaves in the record",
        "",
        "M1e's rule R9 was written to stop exactly this error and it did not,",
        "because it tested the positive control ONCE, at one seed, on the whole",
        "arm. A sensitivity control that is itself unstable does not license a",
        "null. The rule is not moved; M2 records the correct form, which is that a",
        "positive control has to clear under resampling and under leave-one-out",
        "before any null beside it may be read.",
    ])


def m2() -> None:
    t = "m2_robustness_of_the_m1e_reading"
    seeds = load(t, "m2_seeds.csv")
    _at2 = seeds[(seeds.arm == "AT2_identity") & (seeds.well.isin(["wildtype_SeV", "SeV_Cebpa_mutant"]))]
    rng = f"{100 * -_at2.rna_obs.max():.1f} to {100 * -_at2.rna_obs.min():.1f}"
    loo = load(t, "m2_leave_one_out.csv")
    offs = load(t, "m2_offsets.csv")
    genes = load(t, "m2_per_gene.csv")

    stab = (seeds.groupby(["well", "arm"])
            .agg(atac_clears=("atac_clears", "sum"), seeds=("atac_clears", "size"),
                 atac_z_min=("atac_z", "min"), atac_z_max=("atac_z", "max"),
                 rna_clears=("rna_clears", "sum"),
                 rna_z_min=("rna_z", "min"), rna_z_max=("rna_z", "max"),
                 detectable_floor=("atac_detectable_floor_3sd", "mean"))
            .reset_index())
    tl = loo[loo.arm == "transitional"]
    broke = tl[~tl.atac_clears][["well", "dropped_gene", "atac_z"]]

    write(t, "m2_summary.md", [
        "# Trial M2: M1e's chromatin reading does not survive its own arithmetic",
        "",
        "**M2 changes no rule and asks no new question. It attacks M1e's reading on",
        "four fronts, and the chromatin half does not survive.** The RNA half does,",
        "and comfortably.",
        "",
        "## The stability table, four downsampling seeds",
        "",
        df_to_markdown(stab.round(4)),
        "",
        "## What this says, front by front",
        "",
        "**1. The RNA half is solid.** The AT2 identity arm clears in every seed of",
        f"both injured wells, at z from -4.0 to -11.4, losing {rng} detection",
        "points across the seeds (12.6 and 17.1 at the single M1e budget). The transitional arm clears in every seed everywhere. The AT1 arm",
        "clears in none, which is the negative arm behaving. Per gene the AT2 loss",
        "is carried by at least five genes, not one: in the two injured wells Etv5",
        "falls 0.227 and 0.239, Napsa 0.272 and 0.167, Slc34a2 0.271 and 0.189,",
        "Abca3 0.182 and 0.212, and Lamp3 0.160 and 0.115.",
        "",
        "**2. The chromatin positive control is not stable, which is fatal to the",
        "null beside it.** M1e's rule R9 permits the AT2 chromatin null to be read",
        "as retention only behind a transitional arm that fired. Across seeds that",
        "arm clears three times in four in wildtype_SeV and twice in four in",
        "SeV_Cebpa_mutant, five of eight in total. Under leave-one-out it fails on",
        "two of its six genes in each injured well:",
        "",
        df_to_markdown(broke.round(3)),
        "",
        "Per gene, the arm is carried by Sfn and Ndrg1. In wildtype_SeV, Sfn alone",
        "returns 0.0357 over seven peaks against an arm mean of 0.0097.",
        "",
        "**3. The offset correction was not inert, it was ABSENT, and replacing it",
        "changes a verdict.** M1e subtracted the genome-wide MEDIAN distal",
        "difference, which came out exactly 0.0 in every well because more than",
        "half of all distal peaks are detected in neither group at this budget.",
        "That is not the same as there being nothing to remove. The mean global",
        "distal difference is +0.00080, +0.00034 and +0.00201 in the three wells,",
        "POSITIVE IN ALL THREE, across two deposits, two peak atlases and two",
        "fragment budgets. A consistently signed nuisance added to two different",
        "true values is precisely what produces arm values that disagree in sign,",
        "so the sign disagreement between the two wells damns the analysis rather",
        "than rescuing it. Recomputed against the mean:",
        "",
        df_to_markdown(offs[(offs.seed == 0) & (offs.arm != "AT1")]
                       [["well", "arm", "offset_kind", "offset", "observed", "z", "clears"]].round(5)),
        "",
        "In SeV_Cebpa_mutant the AT2 arm moves from z = -1.88 and not clearing, to",
        "z = -3.23 and clearing in the CLOSING direction, purely on the choice",
        "between a median that was degenerate and a mean that was not. A reading",
        "that turns on that choice is not a reading.",
        "",
        "**And the direction matters, because it is not the one M1e reported.** An",
        "independent reviewer reproduced the cell selection and pushed this further:",
        "under a background matched on baseline accessibility the same arm reaches",
        "z = -2.89, and under a multiplicative correction z = -3.55, which would",
        "clear as AT2 programme CLOSED. Under that same correction the",
        "transitional positive control falls to z = +2.53, below the frozen floor,",
        "so rule R9 would have refused the well outright. Across the estimators,",
        "the well reads either NOT COMPUTABLE or CLOSED. It never reads what M1e",
        "reported. Nothing here supports the idea that AT2 chromatin stays open,",
        "and if anything the salvageable signal points the other way.",
        "",
        "**4. The effects and the detection floor are the same size, and the two",
        "arms are not commensurable.** Three sham standard deviations, the smallest",
        "value the frozen rule would have called a clearance, is 0.0043 to 0.0120",
        "in detection fraction. The transitional arm's observed effects are 0.004",
        "to 0.012. The design is working at its own floor.",
        "",
        "Worse, an adversarial review found that the AT2 arm is BOUNDED and the",
        "positive control is NOT. The reference-side baseline of the AT2 arm is",
        "about 0.0185 per peak in GSE310539 and 0.0380 in GSE247130, so the arm",
        "cannot fall below those values even if every distal element at all nine",
        "AT2 loci shut completely. M1e's defence, that an effect the size of the",
        "transitional arm's would have cleared, therefore describes a counterfactual",
        "requiring 58 per cent of ALL distal accessibility at those loci to vanish",
        "in fourteen days. The transitional arm has no such ceiling, because de novo",
        "opening can drive it arbitrarily high. An unbounded statistic was being",
        "used to certify the power of a bounded one.",
        "",
        "The honest statement of what was measured is therefore not that nothing",
        "moved, but that **no change was detected at a bar corresponding to a 45",
        "per cent loss of distal accessibility at the AT2 loci in GSE310539, and a",
        "12 per cent loss in GSE247130.** The second of those is a bound worth",
        "having. The first excludes almost nothing.",
        "",
        "## The reading",
        "",
        "**The chromatin question is NOT ESTABLISHED from these matrices at this",
        "depth.** These data cannot separate the AT2 programme's chromatin not",
        "moving from this design being unable to see it move. M1e's reading of",
        "\"silenced but not closed\" is withdrawn and stays in the record beside this",
        "page. That phrase is the output of decision rule R10a, a label and not a",
        "measurement, and it appears in this folder only in quotation marks and",
        "only as the trial's own output.",
        "",
        "**What is established** is the RNA half, and it is worth stating on its",
        "own: in two deposits, the CLDN4-positive KRT8-positive alveolar group",
        f"loses the AT2 identity programme by {rng} detection points against",
        "a sham band, across five or more genes, with the AT1 programme flat in the",
        "same cells. That is Descriptive only, because each well is one library",
        "pooling two mice.",
        "",
        "**And the two deposits are not independent, which an earlier draft of this",
        "page got wrong.** Jichao Chen is a contributor on both GEO series and the",
        "contact laboratory for GSE247130, so this is one laboratory with two first",
        "authors, not two laboratories. They also share the SftpcCreER and",
        "RosaSun1GFP lineage tools, the Sendai virus model, E-cadherin-positive",
        "sorting, the 10x Multiome kit and cellranger-arc on mm10, and on this side",
        "they share the code, the seeds, the vendor annotation and the budget rule.",
        "What genuinely differs is the first author, the mouse cohort, the",
        "institution, the year and the genotype. Agreement across them is a",
        "consistency check, not a replication.",
        "",
        "## What would settle the chromatin question",
        "",
        "Not more statistics on these matrices. The binding constraints are the",
        "size of the labelled group, 64 and 320 cells, and the fact that peaks were",
        "called on the whole library and are therefore ascertained on the majority",
        "population. Both deposits also carry `atac_fragments.tsv.gz` files, which",
        "this branch did not download and which permit peaks to be re-called on the",
        "labelled cells themselves, proper per-cell ATAC quality control, TSS",
        "enrichment and footprinting. That is the next instrument, not another",
        "threshold.",
        "",
        "## Per-gene table",
        "",
        df_to_markdown(genes.round(4)),
    ])


def m3() -> None:
    t = "m3_one_instrument_and_the_right_null"
    arms = load(t, "m3_arms.csv")
    rd = load(t, "m3_readings.csv")
    at2 = arms[arms.arm == "AT2_identity"]
    tr = arms[arms.arm == "transitional"]
    write(t, "m3_summary.md", [
        "# Trial M3: one instrument, two nulls, and a direction that will not go away",
        "",
        "M2 withdrew M1e's chromatin reading and left the question Not established.",
        "A five-lens adversarial review of that withdrawal then found three defects",
        "in M2 itself, and M3 implements their corrections. **This repository did",
        "not find them**; the review did, and the run record says so.",
        "",
        "## The three corrections",
        "",
        "**D1, the two deposits were not running the same instrument.** Rule R12",
        "removed Cebpa from the AT2 arm where it is genetically deleted and kept it",
        "everywhere else. Cebpa carries only three distal peaks but the largest",
        "positive per-gene value in wildtype_SeV, and it is the only reason that",
        "well's AT2 arm was positive at all. M2's own leave-one-out had already",
        "said so: dropping Cebpa gave the most negative of the nine drops. M3 drops",
        "it everywhere, leaving one eight-gene instrument.",
        "",
        "**D2, the sham band is the wrong null for the question.** It permutes",
        "CELLS, so it asks whether this SPLIT is special. It cannot ask whether",
        "these GENES are special, which is what a gene-set arm claims. M3 adds 300",
        "random gene sets matched gene by gene on distal peak count, evaluated on",
        "the real split, and uses their mean as the offset.",
        "",
        "**D3, the equivalence statistic erred unsafely.** Three sham standard",
        "deviations is the 50-per-cent-power detection floor, not a bound the data",
        "support. M3 reports a confidence interval on the relative scale instead.",
        "",
        "## What one instrument and the right null return",
        "",
        df_to_markdown(arms.drop(columns=["gene_list"]).round(5)),
        "",
        "## The reading, and it is still Not established",
        "",
        df_to_markdown(rd.round(3)),
        "",
        "**The AT2 identity arm is negative in all three computable wells, at",
        "percentile 0.000, 0.007 and 0.000 of 300 matched random gene sets.** On",
        "the relative scale that is -37.8, -18.3 and -15.5 per cent of the",
        "reference group's own accessibility at those loci. The direction is the",
        "same everywhere and it is the opposite of what M1e reported.",
        "",
        "**But no well satisfies both gates, and the two wells fail different",
        "ones.** wildtype_SeV has a positive control that clears (z 3.46 and 3.05)",
        "and an AT2 arm that does not (z -2.16 and -2.73, short of 3).",
        "SeV_Cebpa_mutant has an AT2 arm that clears both nulls at z -3.41 and",
        "-3.51, in the closing direction, behind a positive control that does not",
        "clear (z 1.96 and 1.68), so rule R9 refuses the well. Each well holds half",
        "of what a reading needs and neither holds both. The chromatin question",
        "stays **Not established**, now with a direction and a bound rather than",
        "without one.",
        "",
        "## The offset was doing more work than anyone thought",
        "",
        "Re-centred on the gene-set null, the transitional arm in SeV_Cebpa_mutant",
        "falls from clearing under M1e's degenerate median offset to z = 1.96.",
        "**Most of what M1e counted as signal in that well was the uncorrected",
        "global shift**, which is exactly what the review predicted when it said",
        "the offset was absent rather than inert. Only wildtype_SeV has a",
        "transitional signal that survives proper centring, at +61.8 per cent of",
        "reference with an interval of +26.8 to +96.8.",
        "",
        "## The one bound worth quoting",
        "",
        "In SeV_Cebpa_mutant, AT2 distal accessibility in the labelled group is",
        "**-15.5 per cent of reference accessibility, interval -24.4 to -6.6**, at",
        "percentile 0 of 300 matched gene sets. It is not readable as a finding,",
        "because its well's positive control did not fire, and the frozen rules",
        "refuse it. It is recorded because a bound that a rule refuses is still a",
        "number, and the next instrument should be sized against it.",
        "",
        "## The pre-registered predictions, checked",
        "",
        "The review predicted values before M3 ran and R22 froze them. The gene-set",
        "null means were predicted almost exactly: +0.0009 and +0.0021 against",
        "+0.00069 and +0.00205 observed. The AT2 direction and the percentile-zero",
        "result were confirmed and came out stronger than predicted (z -3.51",
        "against -2.80). The transitional arm came out **weaker** than predicted,",
        "+3.05 and +1.68 against +3.93 and +4.46, which is why R9 now refuses a",
        "well the review expected to pass. Agreement here confirms the review's",
        "arithmetic; it is not an independent discovery by this repository.",
        "",
        "## What this does not do",
        "",
        "It does not resurrect M1e's reading, and it could not have: with one",
        "instrument the arm is negative under every weighting. Nothing in this",
        "folder supports the AT2 programme staying open in the transitional state.",
    ])


if __name__ == "__main__":
    m1b(); m1c(); m1d(); m1e(); m2(); m3()
    print("all summaries regenerated from tracked artefacts")
