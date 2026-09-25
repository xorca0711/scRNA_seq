"""Plot RQ evidence from saved tables; no refitting or new biological inference."""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis/figures/rq"
LIGAND = "analysis/corrections/ligand/results/"
STATS = "analysis/corrections/statistics/tables/"
ES = "Thesis/epithelial_state_specificity/results/"
TEAL, ORANGE, BLUE = "#087e8b", "#cf663b", "#4967a4"
INPUTS, OUTPUTS = {}, []
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 12, "axes.titleweight": "bold", "svg.fonttype": "none",
    "svg.hashsalt": "lung-rq-evidence", "savefig.facecolor": "white"})


def read(relative):
    path = ROOT / relative
    INPUTS[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path)


def save(fig, name):
    for suffix in ("png", "svg"):
        path = OUT / f"{name}.{suffix}"
        fig.savefig(path, dpi=180, bbox_inches="tight",
                    **({"metadata": {"Date": None}} if suffix == "svg" else {}))
        if suffix == "svg":
            # Matplotlib emits spaces at path-line ends; normalize display XML.
            path.write_text("\n".join(x.rstrip() for x in path.read_text(encoding="utf-8").splitlines()) + "\n", encoding="utf-8")
        OUTPUTS.append(str(path.relative_to(ROOT)).replace("\\", "/"))
    plt.close(fig)


def a2():
    d = read(LIGAND + "c37/per_donor.csv")
    path = ROOT / LIGAND / "c37/summary.json"
    INPUTS[str(path.relative_to(ROOT)).replace("\\", "/")] = hashlib.sha256(path.read_bytes()).hexdigest()
    included = json.loads(path.read_text())["detection_at_1000"]["donors"]
    d = d[d.donor.isin(included)]
    rank = read(LIGAND + "lr/resource_summary.csv")
    assert len(d) == 20 and d.groupby("donor").size().eq(2).all()
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.9), layout="constrained")
    for ax, field, title in zip(axes[:2], ["raw_detection", "detection_at_1000"],
                                ["a  Observed detection", "b  Detection at 1,000 UMIs"]):
        pairs = d.pivot(index="donor", columns="compartment", values=field)
        for _, r in pairs.iterrows():
            ax.plot([0, 1], 100*r[["myeloid", "epithelial"]].to_numpy(),
                    "-o", color=TEAL, alpha=.55, lw=1, ms=4)
        ax.set(xticks=[0, 1], xticklabels=["Myeloid", "Epithelial"],
               ylabel="AREG detection (%)", ylim=(0, 100), xlim=(-.3, 1.3), title=title)
        ax.text(.04, .94, "10 paired donors", transform=ax.transAxes, va="top")
    r = rank[rank.areg_median_rank.notna()].copy()
    axes[2].scatter(r.areg_median_rank, np.arange(len(r)), s=48, color=ORANGE)
    for i, row in enumerate(r.itertuples()):
        axes[2].text(row.areg_median_rank+1, i, f"{row.areg_median_rank:g} / {row.pairs}", va="center", fontsize=9)
    axes[2].set(yticks=np.arange(len(r)), yticklabels=r.resource,
                xlabel="AREG–EGFR median donor rank", xlim=(0, 65),
                title="c  Resource-dependent rank")
    axes[2].invert_yaxis()
    fig.suptitle("AREG source detection and ligand rank", fontsize=15)
    save(fig, "rq_a2_source_rank")


def a6():
    fractions = read(STATS + "subtype_cell_fractions.csv")
    shares = read(STATS + "subtype_transcript_contributions.csv")
    f = fractions[fractions.compartment.eq("macrophages")]
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4), layout="constrained")
    for ax, cohort, title in zip(axes[:2], ["GSE136831", "GSE135893"],
                                ["a  Discovery composition", "b  Validation composition"]):
        t = f[f.cohort.eq(cohort)]
        assert len(t) == 2
        for j, group in enumerate(["control", "IPF"]):
            vals = 100*t[f"mean_cell_fraction_{group}"].to_numpy()
            ax.bar(np.arange(2)+(j-.5)*.32, vals, width=.30,
                   color=[TEAL, ORANGE][j], label=group)
            for x, v in zip(np.arange(2)+(j-.5)*.32, vals):
                ax.text(x, v+1.5, f"{v:.1f}", ha="center", fontsize=9)
        labels = [s.replace("_", "\n").replace("Proliferating ", "Proliferating\n") for s in t.label]
        ax.set(xticks=[0, 1], xticklabels=labels, ylim=(0, 109),
               ylabel="Mean donor cell fraction (%)", title=title)
        ax.set_xlabel(f"{cohort}: {int(t.iloc[0].n_IPF)} IPF / {int(t.iloc[0].n_control)} control donors", fontsize=9)
    p = f[(f.cohort.eq("GSE135893")) & f.label.eq("Proliferating Macrophages")].iloc[0]
    t = shares[(shares.cohort.eq("GSE135893")) & shares.label.eq("Proliferating Macrophages")].set_index("set")
    for j, group in enumerate(["control", "IPF"]):
        vals = [p[f"mean_cell_fraction_{group}"], t.loc["HALLMARK_E2F_TARGETS", f"mean_set_transcript_share_{group}"],
                t.loc["HALLMARK_G2M_CHECKPOINT", f"mean_set_transcript_share_{group}"]]
        axes[2].plot(np.arange(3), np.array(vals)*100, "o-", color=[TEAL, ORANGE][j], label=group)
    axes[2].set(xticks=[0, 1, 2], xticklabels=["Cells", "E2F\ntranscripts", "G2M\ntranscripts"],
                ylabel="Mean donor share (%)", ylim=(0, 19), title="c  Proliferating-state contribution")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=2, frameon=False)
    fig.suptitle("Macrophage composition and programme contributions", fontsize=15)
    save(fig, "rq_a6_composition")


def a7():
    scores = read(ES + "module_scores.csv")
    s = scores[(scores.module.eq("AT2_published_holdout")) & scores.interpretable]
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.2), sharey=True, layout="constrained")
    for ax, prefix, title in zip(axes, ["P9", "SeV"], ["a  P9, uninjured", "b  Adult, SeV injury"]):
        for genotype, color in [("control", TEAL), ("Cebpa_mutant", ORANGE)]:
            t = s[s.unit.eq(f"{prefix}_{genotype}")].pivot(index="seed", columns="group", values="score")
            assert t.shape == (2, 2) and set(t.index) == {17, 29}
            y = 100*t[["reference", "two_marker"]]
            ax.plot([0, 1], y.mean(), color=color, lw=1.5)
            for seed, marker, offset in [(17, "o", -.025), (29, "D", .025)]:
                ax.scatter(np.array([0, 1])+offset, y.loc[seed], color=color, marker=marker, s=42, zorder=3)
        ax.set(xticks=[0, 1], xticklabels=["Reference AT2", "Cldn4/Krt8-labelled"],
               xlim=(-.25, 1.25), title=title)
    axes[0].set_ylabel("AT2 holdout mean gene detection (%)")
    handles = [Line2D([], [], color=TEAL, label="Control"), Line2D([], [], color=ORANGE, label="Cebpa mutant"),
               Line2D([], [], color=".4", marker="o", ls="", label="Seed 17"), Line2D([], [], color=".4", marker="D", ls="", label="Seed 29")]
    fig.legend(handles=handles, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("AT2 programme contrasts across genotype and cell state", fontsize=14)
    save(fig, "rq_a7_genotype_reference")


def a8():
    overlap = read(ES + "module_overlap.csv")
    effects = read(ES + "within_unit_effects.csv")
    scores = read(ES + "module_scores.csv")
    modules = ["ADI_published_holdout", "AT1_published_holdout", "AT2_published_holdout"]
    matrix = np.zeros((3, 3), dtype=int)
    for i, a in enumerate(modules):
        matrix[i, i] = int(scores[scores.module.eq(a)].total_genes.iloc[0])
        for j, b in enumerate(modules):
            if i != j:
                row = overlap[(overlap.module_a.eq(a) & overlap.module_b.eq(b)) | (overlap.module_a.eq(b) & overlap.module_b.eq(a))]
                assert len(row) == 1
                matrix[i, j] = int(row.intersection.iloc[0])
    assert matrix[0, 1] == 119
    fig, axes = plt.subplots(1, 2, figsize=(10.3, 4.2), layout="constrained", gridspec_kw={"width_ratios": [1, 1.45]})
    axes[0].imshow(matrix, cmap="Blues", vmin=0, vmax=400)
    for (i, j), value in np.ndenumerate(matrix):
        axes[0].text(j, i, str(value), ha="center", va="center", color="white" if value>220 else "black")
    axes[0].set(xticks=range(3), xticklabels=["ADI", "AT1", "AT2"], yticks=range(3),
                yticklabels=["ADI", "AT1", "AT2"], title="a  Shared holdout genes")
    units = ["P9_control", "P9_Cebpa_mutant", "SeV_control", "SeV_Cebpa_mutant", "EEM-scRNA-167"]
    mods = ["ADI_published_holdout", "AT1_published_holdout", "AT1_late_panel"]
    v, flip = np.zeros((5, 3)), np.zeros((5, 3), dtype=bool)
    for i, unit in enumerate(units):
        for j, module in enumerate(mods):
            t = effects[effects.unit.eq(unit) & effects.module.eq(module) & effects.interpretable]
            assert len(t) == 2 and set(t.seed) == {17, 29}
            v[i, j] = t.difference.mean()*100
            flip[i, j] = t.difference.min()<0<t.difference.max()
    im = axes[1].imshow(v, cmap="RdBu_r", vmin=-15, vmax=15, aspect="auto")
    for (i, j), value in np.ndenumerate(v):
        axes[1].text(j, i, f"{value:+.1f}"+("*" if flip[i, j] else ""), ha="center", va="center",
                     color="white" if abs(value)>10 else "black")
    axes[1].set(xticks=range(3), xticklabels=["ADI holdout", "AT1 holdout", "Late AT1\n4-gene panel"],
                yticks=range(5), yticklabels=["P9 control", "P9 Cebpa mutant", "SeV control", "SeV Cebpa mutant", "External mouse 167"],
                title="b  Labelled minus reference detection")
    fig.colorbar(im, ax=axes[1], label="Detection difference (percentage points)", shrink=.8)
    fig.suptitle("Programme overlap and epithelial-state contrasts", fontsize=15)
    save(fig, "rq_a8_signature_specificity")


def a9():
    d = read(LIGAND + "lr/canonical_egfr_receptor_coverage.csv")
    assert not d.duplicated(["resource", "receptor", "ligand"]).any()
    resources = ["cellchatdb", "cellphonedb", "connectomedb2020", "italk", "consensus"]
    ligands = ["AREG", "HBEGF", "EREG", "TGFA", "BTC", "EGF", "EPGN"]
    pairs = [(r, receptor) for r in resources for receptor in sorted(d[d.resource.eq(r)].receptor.unique())]
    m = d.pivot(index=["resource", "receptor"], columns="ligand", values="n_donors_scored").reindex(index=pd.MultiIndex.from_tuples(pairs), columns=ligands)
    fig, ax = plt.subplots(figsize=(10.4, 5.6), layout="constrained")
    cmap = plt.colormaps["Blues"].copy(); cmap.set_bad("#eeeeee")
    im = ax.imshow(m.to_numpy(), cmap=cmap, vmin=0, vmax=22, aspect="auto")
    for i, (resource, receptor) in enumerate(pairs):
        for j, ligand in enumerate(ligands):
            row = d[d.resource.eq(resource) & d.receptor.eq(receptor) & d.ligand.eq(ligand)]
            if row.empty:
                label = "—"
            else:
                r = row.iloc[0]
                label = f"{int(r.n_donors_scored)}/22"
            ax.text(j, i, label, ha="center", va="center", color="white" if m.iloc[i,j]>14 else "black", fontsize=10)
    ax.set(xticks=range(7), xticklabels=ligands, yticks=range(len(pairs)), yticklabels=[f"{r}\n{receptor}" for r, receptor in pairs],
           title="Canonical EGFR ligand coverage across resources")
    fig.colorbar(im, ax=ax, label="Donors with a score (retention requires ≥11)", shrink=.85)
    save(fig, "rq_a9_receptor_coverage")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    for plot in (a2, a6, a7, a8, a9):
        plot()
    record = {"scope": "Descriptive figures from saved tables; no new hypothesis tests or embeddings.",
              "inputs_sha256": INPUTS,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "outputs_sha256": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in OUTPUTS},
              "checks": ["10 donor pairs selected by primary-budget eligibility", "Two deposited macrophage states per cohort", "Two technical seeds per eligible ES1 comparison", "119 ADI/AT1 shared holdout genes", "Unique resource/receptor/ligand coverage rows"]}
    (OUT / "rq_evidence_figures.json").write_text(json.dumps(record, indent=2)+"\n", encoding="utf-8")
    print("Saved five empirical figures (PNG and SVG), with input/output hashes.")
