"""Generate the A0 feasibility closeout, figures, and explicit downstream gate."""
from pathlib import Path
from string import Template
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]


def md_table(frame):
    rows = ["| " + " | ".join(map(str, frame.columns)) + " |",
            "| " + " | ".join("---" for _ in frame.columns) + " |"]
    rows += ["| " + " | ".join(map(str, r)) + " |" for r in frame.itertuples(index=False, name=None)]
    return "\n".join(rows)


def main():
    audit = json.loads((BASE / "extended_audit_record.json").read_text())
    developmental = pd.read_csv(BASE / "tables/d2_negretti_encoded_group_coverage.csv")
    gut = pd.read_csv(BASE / "tables/v1_haber_verified_mouse_coverage.csv")
    repair = pd.read_csv(BASE / "tables/d1_strunz_proposed_window.csv")
    planning = pd.read_csv(BASE / "tables/developmental_capture_planning.csv")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(14, 6.5), gridspec_kw={"width_ratios": [1, 1.6, 0.85]})
    r = repair[repair.in_proposed_day10_15_window].sort_values(["day", "sample_id"])
    panels = [(r, ["AT2", "Krt8+ ADI", "AT1"], r.sample_id, "Repair: proposed days 10–15\nSource mouse libraries"),
        (developmental, ["AT2", "Transitional", "AT1"], developmental.encoded_library_group,
         "Development: encoded library groups\nAnimal/hash mapping unavailable"),
        (gut, ["Stem", "Enterocyte.Immature.Proximal", "Enterocyte.Mature.Proximal"], gut.mouse_id,
         "Intestine: four mapped mice\nProximal enterocyte branch")]
    for ax, (frame, states, labels, title) in zip(axes, panels):
        y = np.arange(len(frame))
        for i, (state, color, label) in enumerate(zip(states, ["#3978a8", "#c77634", "#52864f"],
                                                    ["Starting state", "Intermediate", "Destination"])):
            ax.scatter(frame[state], y + (i - 1) * 0.2, s=27, color=color, label=label)
        ax.axvline(30, color="#8c3434", linestyle="--", linewidth=1)
        ax.set(xscale="symlog", xlim=(-0.2, float(frame[states].max().max()) * 1.4), xlabel="Cells per group (symlog scale)",
               yticks=y, yticklabels=labels, title=title)
        ax.invert_yaxis()
        ax.grid(axis="x", alpha=0.15)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.935), ncol=3, frameon=False)
    fig.suptitle("A0 feasibility: both endpoint comparisons require 30 cells in every state", fontsize=14)
    fig.text(0.5, 0.015, "Dots are annotation counts, not expression effects. Developmental groups are not verified biological replicates.",
             ha="center", fontsize=10)
    fig.tight_layout(rect=(0, 0.05, 1, 0.86))
    fig.savefig(BASE / "figures/extended_feasibility_coverage.png", dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.2))
    x = np.arange(len(planning))
    a = planning.retained_epithelial_cells_per_independent_unit.to_numpy()
    b = planning.target_if_intermediate_fraction_halves.to_numpy()
    ax.bar(x - .18, a, .36, color="#3978a8", label="Observed pooled cell fractions")
    ax.bar(x + .18, b, .36, color="#adbfd0", label="Intermediate fraction halved")
    for position, value in zip(x - .18, a):
        ax.text(position, value + 250, f"{value:,}", ha="center", fontsize=9)
    ax.set(xticks=x, xticklabels=planning.timepoint, ylim=(0, 22000),
           ylabel="Retained epithelial cells per new independent sample",
           title="Illustrative capture requirement: ≥30 cells in each state\n95% lower bound under a fixed multinomial sampling model")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(axis="y", alpha=.15)
    fig.text(.5, .018, "Capture planning only. Fractions depend on stage, sorting and filtering; these are not gene-effect power estimates.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, .045, 1, 1))
    fig.savefig(BASE / "figures/developmental_capture_planning.png", dpi=160)
    plt.close(fig)

    readiness = {
        "current_status": "FEASIBILITY_COMPLETE_PRIMARY_ANALYSIS_BLOCKED",
        "biological_decision": "UNRESOLVED",
        "full_pilot_complete": False,
        "primary_analysis_ready": False,
        "expression_program_learned": False,
        "transfer_test_performed": False,
        "completed": ["Initial and extended public metadata audit", "Published-state coverage counts",
                      "Repair-window feasibility", "Developmental capture planning", "Report and reproducibility checks"],
        "roles": {
            "D1": {"status": "feasible_candidate_not_frozen", "accession": "GSE141259",
                   "complete_triplets_proposed_window": int(repair.coverage_pass_in_proposed_window.sum())},
            "D2": {"status": "blocked_coverage_and_independent_unit_metadata",
                   "complete_developmental_encoded_groups": int(developmental.complete_at_floor.sum())},
            "V1": {"status": "blocked_replication_and_mouse_mapping",
                   "complete_verified_intestinal_mice": int(gut.triplet_eligible.sum())}},
        "blocked_stages": {
            "P1": "Need qualifying D2 and V1 cohorts, unit IDs, raw counts and independent state evidence before freezing all cohorts/settings.",
            "P2": "Cannot learn a two-context program without eligible D2; one-context selection would answer a different question.",
            "P3": "Requires a frozen discovery program and an eligible held-out tissue cohort.",
            "P4": "Requires the program and evaluable cohorts; no specificity evidence has been generated."},
        "unblocking_requirements": [
            "D2: cell-to-animal or independent-pool map and source-supported start/intermediate/destination states; at least 3 independent units, each with at least 30 cells/state.",
            "V1: verified batch-to-mouse map that yields at least 3 qualifying units, or an eligible independent non-lung epithelial cohort.",
            "Both: public raw count matrices matching those annotations and documented transition evidence beyond an unvalidated pseudotime label."],
        "not_claimed": ["All public datasets exhausted", "Universal process established", "Conservation disproved", "P1–P4 completed"]}
    (BASE / "readiness.json").write_text(json.dumps(readiness, indent=2) + "\n", encoding="utf-8")

    cohort_audit = pd.read_csv(BASE / "tables/dataset_audit.csv", keep_default_na=False)
    d1 = (cohort_audit.accession == "GSE141259") & (cohort_audit.role == "D1_alternative")
    cohort_audit.loc[d1, "replication_evidence"] = "Primary paper describes two independent mice per sampled time; available subset contains 32 libraries versus 36 mice reported"
    cohort_audit.loc[d1, "status"] = "feasible_candidate_not_frozen"
    cohort_audit.loc[d1, "remaining_requirement"] = "Use paired within-mouse effects across a prespecified repair interval; proposed days 10-15 yield 9 triplets. Full P1 waits for D2/V1."
    d2 = cohort_audit.accession == "GSE165063+GSE160876"
    cohort_audit.loc[d2, "coverage_evidence"] = "10918 author-annotated epithelial cells; 14 barcode-suffix groups; only 1 group meets all three 30-cell floors"
    cohort_audit.loc[d2, "replication_evidence"] = "Public viewer omits mouse/hash IDs; suffixes indicate capture groups, not verified independent mice; age pooling is not replication"
    cohort_audit.loc[d2, "status"] = "not_eligible_current_public_annotations"
    cohort_audit.loc[d2, "remaining_requirement"] = "Recover original unit/hash identities and >=3 qualifying independent units; matching raw counts also needed. Do not merge pools to cross floors."
    skin = {"role": "V1_fallback", "accession": "GSE67602", "setting": "normal mouse epidermal differentiation",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5052454/", "expression_downloaded_for_A0": False,
        "coverage_evidence": "1422 cells, 34 capture batches, published state labels; no batch-to-mouse map recovered",
        "replication_evidence": "Paper reports 19 mice; 34 batches cannot be counted as 34 animals",
        "status": "unit_mapping_unresolved", "remaining_requirement": "Recover mouse-to-batch map, then audit one source-supported IFE trajectory; no primary state merging or gene scoring performed"}
    if not (cohort_audit.accession == "GSE67602").any():
        cohort_audit = pd.concat([cohort_audit, pd.DataFrame([skin])], ignore_index=True)
    cohort_audit.to_csv(BASE / "tables/dataset_audit.csv", index=False)

    age_table = pd.read_csv(BASE / "tables/d2_negretti_age_coverage.csv")
    dtable = developmental[["encoded_library_group", "AT2", "Transitional", "AT1", "complete_at_floor"]].copy()
    dtable.columns = ["Encoded group (not mouse ID)", "AT2", "Intermediate", "AT1", "All ≥30"]
    ptable = planning[["timepoint", "intermediate_fraction", "retained_epithelial_cells_per_independent_unit",
                      "target_if_intermediate_fraction_halves"]].copy()
    ptable.intermediate_fraction = ptable.intermediate_fraction.map(lambda x: f"{100*x:.2f}%")
    ptable.columns = ["Age", "Observed intermediate fraction", "Cells/sample: point scenario", "Cells/sample: fraction halved"]
    values = {"repair_n": int(repair.coverage_pass_in_proposed_window.sum()),
        "developmental_table": md_table(dtable), "capture_table": md_table(ptable),
        "gut_n": int(gut.triplet_eligible.sum()), "devo_cells": f"{len(developmental):,}",
        "annotation_cells": f"{audit['negretti']['cells']:,}"}
    output = Template((BASE / "scripts/feasibility_report.md.in").read_text(encoding="utf-8")).substitute(values)
    (BASE / "reports/PILOT_REPORT.md").write_text(output, encoding="utf-8")
    print("Generated feasibility report, 2 figures, updated candidate audit and readiness gate; primary analysis remains blocked.")


if __name__ == "__main__":
    main()
