#!/usr/bin/env python
"""Does the 6 dpi inflammatory-monocyte state survive batch correction on the
strongest technical key available?

11_myeloid_focus.py embedded the myeloid compartment without batch correction
(root decision 2: every sample is one animal at one day, so correcting on
sample would remove the time course). The reviewer's question is whether the
dense 6 dpi iMON state is biology or a one-day batch island. Correcting on
sample cannot answer it. The paper's Table S3 records the date of influenza
infection, and for every active-repair day and the 42 dpi harvest the two
replicate animals came from two different infection rounds (October 2021 and
July 2022), so infection round is a technical key that crosses time on those
days. It also carries the Ki67-Cre dosage (Cre/Cre in most October 2021
animals, Cre/+ in all July 2022 animals), so it absorbs date and genotype
together. Days with a single round (0, 90, 366 dpi) cannot be tested and are
excluded.

The subset is re-embedded twice with identical rules, uncorrected and after
Harmony on infection round, and the two embeddings are compared on the
pipeline's own kNN mixing metric, on the grading against the held-out labels,
and on a pre-registered rule for the survival of the 6 dpi iMON state. Rules
are frozen in the run record before the object is opened. The unit is the
animal; no P values are computed.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import sys  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from focus_utils import REPO, RunRecord, df_to_markdown, read_csr_rows, shannon  # noqa: E402
from pipeline_utils import ANALYSIS, RANDOM_SEED, log, save_fig, setup_matplotlib  # noqa: E402

warnings.filterwarnings("ignore")
plt = setup_matplotlib()
import anndata as ad  # noqa: E402
import scanpy as sc  # noqa: E402
from anndata.io import read_elem  # noqa: E402

sc.settings.verbosity = 0
np.random.seed(RANDOM_SEED)

SERIES = ANALYSIS / "GSE262927"
OBJ = SERIES / "processed" / "final_clustered.h5ad"
MYELOID = SERIES / "myeloid_focus"
OUT = MYELOID / "batch_sensitivity"
FIG = OUT / "figures"
TAB = OUT / "tables"
TABLE_S3 = REPO / "Thesis" / "Primary" / "Longitudinal single-cell profiles of lung regeneration" / "xlsx raw" / "mmc4.xlsx"
ROUND_CSV = TAB / "sample_infection_round.csv"

CLUSTERS = ["5", "17", "24"]
DAYS_TESTED = [6, 11, 19, 25, 42]
MYELOID_LABELS = ["aMAC", "iMAC", "iMON", "cMON", "pMON", "cDC1", "cDC2", "maDC", "Neutrophil"]
LABEL_COLOURS = {"aMAC": "#1f77b4", "iMAC": "#17becf", "iMON": "#2ca02c", "cMON": "#98df8a",
                 "pMON": "#bcbd22", "cDC1": "#d62728", "cDC2": "#ff9896", "maDC": "#e377c2",
                 "Neutrophil": "#9467bd", "other lineage": "#7f7f7f"}
ROUND_COLOURS = {"2021-10-21": "#e6550d", "2022-07-06": "#3182bd"}
PRIMARY_RES = 0.5
K = 30
REPS = {"X_pca": "uncorrected", "X_pca_harmony": "Harmony on infection round"}

RULES = {
    "object": "analysis/GSE262927/processed/final_clustered.h5ad, read row-wise; never loaded fully",
    "cells": "atlas Leiden clusters 5, 17, 24 within the 25-sample Ki67 atlas (as 11_myeloid_focus.py), restricted "
             f"to days where both infection rounds are present: {DAYS_TESTED} dpi. 0, 90 and 366 dpi are "
             "single-round and excluded",
    "batch_key": "infection round = Date of Influenza Infection in the paper's Table S3 (mmc4.xlsx, local and "
                 "gitignored; the sample-to-round map is written to tables/sample_infection_round.csv and tracked). "
                 "The two rounds also differ in Ki67-Cre dosage (Cre/Cre in most 2021-10-21 animals, Cre/+ in all "
                 "2022-07-06 animals), so the key absorbs date and genotype together",
    "why_not_sample": "sample = one animal at one day; correcting on it removes the time course (root decision 2)",
    "embedding": {"hvg": {"flavor": "seurat_v3", "n_top_genes": 2000, "layer": "counts", "batch_key": None},
                  "scale_max": 10, "n_pcs": 30, "knn_k": K,
                  "leiden": {"resolution": PRIMARY_RES, "flavor": "igraph", "n_iterations": 2, "random_state": 0},
                  "umap_random_state": 0,
                  "note": "identical to 11_myeloid_focus.py; computed on X_pca (uncorrected) and on X_pca_harmony"},
    "harmony": {"package": "harmonypy 0.0.10 (pinned)", "vars_use": ["round"], "random_state": 0,
                "input": "the 30 uncorrected PCs"},
    "mixing_metric": "the pipeline's within-group replicate rule (run_scrna_analysis.py inspect_batch) with group = "
                     f"day and replicate = round: for each cell, among its {K} nearest neighbours that lie in the same "
                     "day, the fraction sharing its round; averaged over the day's cells and divided by the expected "
                     "value under mixing (sum of squared round frequencies within the day). 1.0 = mixed; the pipeline "
                     "treats > 2.0 as failure to mix",
    "imon_state_rule": "in each embedding, the Leiden subcluster holding the most iMON-labelled cells is the iMON "
                       "state. It survives if its iMON purity among labelled cells is >= 0.75 AND it contains >= 20 "
                       "percent of the iMON-labelled cells of EACH 6 dpi animal. Reported for both embeddings; the "
                       "question is whether the Harmony embedding still passes",
    "imon_state_rule_revision": "DISCLOSED after the first run (2026-09-10): the frozen definition 'subcluster "
                                "holding the most iMON-labelled cells' selects the 11 to 19 dpi monocyte state in both "
                                "embeddings (3.9% and 22.0% of its cells from 6 dpi), not the 6 dpi state the question "
                                "is about. The first-run outcome is kept below (first_run_outcome). A second, post hoc "
                                "definition anchored on the 6 dpi cells is reported alongside and labelled post hoc: the "
                                "subcluster holding the most iMON-labelled cells FROM 6 dpi animals, judged by the same "
                                "purity and per-animal criteria. Both definitions are reported for both embeddings",
    "first_run_outcome": {"frozen_rule_passes": {"uncorrected": False, "Harmony on infection round": True},
                          "selected_subcluster_pct_cells_from_6dpi": {"uncorrected": 3.9, "Harmony on infection round": 22.0},
                          "adjusted_rand_index": 0.837, "n_subclusters": {"uncorrected": 14, "Harmony on infection round": 15},
                          "same_round_knn_enrichment_overall": {"uncorrected": 1.311, "Harmony on infection round": 1.072}},
    "agreement": "adjusted Rand index between the two Leiden partitions; grading (purity, label entropy, threshold "
                 "0.56) against the held-out deposited labels in both embeddings",
    "unit": "animal; descriptive; no P values",
    "seed": RANDOM_SEED,
}


def sample_round_table(rec: RunRecord) -> pd.DataFrame:
    if TABLE_S3.exists():
        rec.add_input(TABLE_S3)
        x = pd.read_excel(TABLE_S3, sheet_name=0, header=0)
        x.columns = [str(c) for c in x.columns]
        col = {k: next(c for c in x.columns if k in c) for k in
               ("Sequencing ID", "GEO Sample Name", "Date of Influenza Infection", "Euthanasia (days", "Genotype", "Sex")}
        t = pd.DataFrame({
            "sample_id": x[col["Sequencing ID"]].astype(str).str.strip().str.replace("_", "-", regex=False),
            "geo_sample_name": x[col["GEO Sample Name"]].astype(str).str.strip(),
            "infection_date": pd.to_datetime(x[col["Date of Influenza Infection"]], errors="coerce").dt.strftime("%Y-%m-%d"),
            "euthanasia": x[col["Euthanasia (days"]].astype(str).str.strip(),
            "genotype": x[col["Genotype"]].astype(str).str.strip(),
            "sex": x[col["Sex"]].astype(str).str.strip(),
        })
        t["round"] = t["infection_date"].fillna("uninfected")
        t = t.drop(columns=["infection_date"])
        TAB.mkdir(parents=True, exist_ok=True)
        t.to_csv(ROUND_CSV, index=False)
        rec.add_output(ROUND_CSV)
        rec.set("sample_round_source", "Table S3 (mmc4.xlsx), rewritten to tables/sample_infection_round.csv")
    elif ROUND_CSV.exists():
        rec.add_input(ROUND_CSV)
        t = pd.read_csv(ROUND_CSV, dtype=str)
        rec.set("sample_round_source", "tracked tables/sample_infection_round.csv (Table S3 not present locally)")
    else:
        raise SystemExit("neither Table S3 (mmc4.xlsx) nor tables/sample_infection_round.csv is available")
    return t


def load_subset(rec: RunRecord, rounds: pd.DataFrame) -> ad.AnnData:
    if not OBJ.exists():
        raise SystemExit(f"source object not found: {OBJ}")
    rec.add_input(OBJ)
    log(f"opening {OBJ.name} row-wise")
    with h5py.File(OBJ, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        n_cols = int(f["X"].attrs["shape"][1])
        day = pd.to_numeric(obs["sacrifice_day"], errors="coerce")
        mask = (obs["leiden_cluster"].astype(str).isin(CLUSTERS).to_numpy()
                & obs["has_author_metadata"].astype(str).eq("True").to_numpy()
                & day.isin(DAYS_TESTED).to_numpy())
        idx = np.where(mask)[0]
        X = read_csr_rows(f["X"], idx, n_cols)
        C = read_csr_rows(f["layers"]["counts"], idx, n_cols)
    sub = ad.AnnData(X=X, obs=obs.iloc[idx].copy(), var=var.copy(), layers={"counts": C})
    sub.obs_names = obs.index[idx]
    sub.obs["sample_id"] = sub.obs["sample_id"].astype(str)
    sub.obs["day"] = pd.to_numeric(sub.obs["sacrifice_day"], errors="coerce").astype(int)
    rmap = dict(zip(rounds["sample_id"], rounds["round"]))
    sub.obs["round"] = sub.obs["sample_id"].map(rmap)
    if sub.obs["round"].isna().any():
        missing = sorted(sub.obs.loc[sub.obs["round"].isna(), "sample_id"].unique())
        raise SystemExit(f"no infection round for samples {missing}; refusing to guess")
    per_day = sub.obs.groupby("day")["round"].nunique()
    if (per_day < 2).any():
        raise SystemExit(f"a tested day has a single round: {per_day.to_dict()}")
    lab = sub.obs["author_celltype"].astype(str)
    lin = sub.obs["author_lineage"].astype(str)
    label = np.where(lin == "Myeloid", lab, "other lineage")
    cats = [k for k in MYELOID_LABELS + ["other lineage"] if k in set(label)]
    sub.obs["label"] = pd.Categorical(label, categories=cats)
    sub.uns["label_colors"] = [LABEL_COLOURS[c] for c in cats]
    rcats = sorted(sub.obs["round"].unique())
    sub.obs["round"] = pd.Categorical(sub.obs["round"], categories=rcats)
    sub.uns["round_colors"] = [ROUND_COLOURS.get(c, "#7f7f7f") for c in rcats]
    rec.set("n_cells", int(sub.n_obs))
    rec.set("n_animals", int(sub.obs["sample_id"].nunique()))
    rec.set("animals_by_day_and_round", {str(d): g.groupby("round", observed=True)["sample_id"].unique().apply(list).to_dict()
                                         for d, g in sub.obs.groupby("day")})
    rec.set("cells_by_round", sub.obs["round"].value_counts().to_dict())
    log(f"subset: {sub.n_obs:,} cells, {sub.obs['sample_id'].nunique()} animals, rounds {rcats}")
    return sub


def embed_both(sub: ad.AnnData, rec: RunRecord) -> None:
    sc.pp.highly_variable_genes(sub, flavor="seurat_v3", n_top_genes=2000, layer="counts")
    emb = sub[:, sub.var["highly_variable"]].copy()
    sc.pp.scale(emb, max_value=10)
    sc.tl.pca(emb, n_comps=30, random_state=RANDOM_SEED)
    sub.obsm["X_pca"] = emb.obsm["X_pca"]
    del emb
    import harmonypy  # noqa: E402
    ho = harmonypy.run_harmony(sub.obsm["X_pca"], sub.obs, ["round"], random_state=RANDOM_SEED)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != sub.n_obs:
        Z = Z.T
    sub.obsm["X_pca_harmony"] = np.ascontiguousarray(Z, dtype=np.float32)
    for rep in REPS:
        nb = f"nb_{rep}"
        sc.pp.neighbors(sub, n_neighbors=K, use_rep=rep, key_added=nb, random_state=RANDOM_SEED)
        sc.tl.leiden(sub, resolution=PRIMARY_RES, key_added=f"leiden_{rep}", neighbors_key=nb,
                     flavor="igraph", n_iterations=2, random_state=RANDOM_SEED, directed=False)
        sc.tl.umap(sub, neighbors_key=nb, random_state=RANDOM_SEED)
        sub.obsm[f"X_umap_{rep}"] = sub.obsm["X_umap"].copy()
        log(f"{REPS[rep]}: {sub.obs[f'leiden_{rep}'].nunique()} subclusters")
    rec.set("n_subclusters", {rep: int(sub.obs[f"leiden_{rep}"].nunique()) for rep in REPS})


def round_mixing(X: np.ndarray, days: np.ndarray, rounds: np.ndarray) -> dict:
    from sklearn.neighbors import NearestNeighbors
    nn = NearestNeighbors(n_neighbors=K + 1).fit(X)
    _, ind = nn.kneighbors(X)
    ind = ind[:, 1:]
    out = {}
    for d in DAYS_TESTED:
        m = days == d
        in_day = days[ind][m] == d
        same = rounds[ind][m] == rounds[m][:, None]
        denom = in_day.sum(axis=1)
        ok = denom > 0
        frac = (same & in_day)[ok].sum(axis=1) / denom[ok]
        fr = pd.Series(rounds[m]).value_counts(normalize=True)
        exp = float((fr ** 2).sum())
        out[d] = {"n_cells": int(m.sum()), "cells_with_same_day_neighbours": int(ok.sum()),
                  "observed_same_round_fraction": round(float(frac.mean()), 4), "expected_if_mixed": round(exp, 4),
                  "enrichment": round(float(frac.mean()) / exp, 3)}
    # overall same-round enrichment across all tested cells, as the pipeline's knn_same_sample_enrichment
    same_all = (rounds[ind] == rounds[:, None]).mean(axis=1)
    freq = pd.Series(rounds).value_counts(normalize=True)
    expected = freq.reindex(rounds).to_numpy()
    out["overall"] = {"knn_same_round_enrichment": round(float(np.mean(same_all / expected)), 3)}
    return out


def grade(sub: ad.AnnData, key: str) -> tuple[pd.DataFrame, dict]:
    lab = sub.obs["label"].astype(str)
    labelled = lab != "other lineage"
    tab = pd.crosstab(sub.obs.loc[labelled, key].astype(str), lab[labelled])
    rows = []
    for scl in sorted(sub.obs[key].astype(str).unique(), key=int):
        mask = sub.obs[key].astype(str) == scl
        n_all = int(mask.sum())
        if scl in tab.index:
            counts = tab.loc[scl]
            n_lab = int(counts.sum())
            purity = float(counts.max() / n_lab)
            ent = shannon((counts / n_lab).to_numpy())
            top = str(counts.idxmax())
        else:
            n_lab, purity, ent, top = 0, np.nan, np.nan, ""
        rd = sub.obs.loc[mask, "round"].astype(str).value_counts(normalize=True)
        rows.append({"subcluster": scl, "n_cells": n_all, "n_labelled": n_lab, "top_label": top,
                     "purity": round(purity, 3) if n_lab else None,
                     "label_entropy": round(ent, 3) if n_lab else None,
                     "entropy_flag": ("HIGH" if ent > 0.56 else "low") if n_lab else "NA",
                     "n_animals": int(sub.obs.loc[mask, "sample_id"].nunique()),
                     "max_round_fraction": round(float(rd.iloc[0]), 3), "dominant_round": rd.index[0],
                     "pct_cells_from_6dpi": round(100 * float((sub.obs.loc[mask, "day"] == 6).mean()), 1)})
    df = pd.DataFrame(rows)
    n_lab_total = int(labelled.sum())
    in_pure = int(df.loc[df["purity"].fillna(0) >= 0.75, "n_labelled"].sum())
    summary = {"fraction_labelled_cells_in_pure_subclusters": round(in_pure / n_lab_total, 4),
               "cell_weighted_mean_entropy": round(float(np.nansum(df["label_entropy"].fillna(0) * df["n_labelled"]) / n_lab_total), 4),
               "n_subclusters": int(len(df))}
    return df, summary


DEFINITIONS = {"frozen_rule_most_iMON_cells": "frozen rule: subcluster with the most iMON-labelled cells",
               "post_hoc_most_6dpi_iMON_cells": "post hoc (disclosed): subcluster with the most iMON-labelled cells from 6 dpi animals"}


def imon_state(sub: ad.AnnData, key: str) -> dict:
    """Both definitions of 'the iMON state', judged by the same criteria."""
    lab = sub.obs["label"].astype(str).to_numpy()
    cl = sub.obs[key].astype(str).to_numpy()
    imon = lab == "iMON"
    labelled = lab != "other lineage"
    six = (sub.obs["day"] == 6).to_numpy()
    samples = sub.obs["sample_id"].astype(str).to_numpy()
    rounds = sub.obs["round"].astype(str).to_numpy()
    picks = {"frozen_rule_most_iMON_cells": pd.Series(cl[imon]).value_counts().idxmax(),
             "post_hoc_most_6dpi_iMON_cells": pd.Series(cl[imon & six]).value_counts().idxmax()}
    out = {}
    for name, top in picks.items():
        in_top = cl == top
        purity = float((lab[in_top & labelled] == "iMON").mean())
        per_animal = {}
        for a in sorted(set(samples[six])):
            m = six & (samples == a) & imon
            per_animal[a] = round(float(in_top[m].mean()), 3) if m.sum() else np.nan
        survives = bool(purity >= 0.75 and all(v >= 0.20 for v in per_animal.values()))
        rd = pd.Series(rounds[in_top]).value_counts(normalize=True).round(3).to_dict()
        out[name] = {"definition": DEFINITIONS[name], "subcluster": top, "n_cells": int(in_top.sum()),
                     "iMON_purity_among_labelled": round(purity, 3),
                     "fraction_of_each_6dpi_animals_iMON_cells_in_subcluster": per_animal,
                     "pct_cells_from_6dpi": round(100 * float(six[in_top].mean()), 1),
                     "round_composition": rd, "survives": survives}
    return out


def main() -> int:
    for d in (OUT, FIG, TAB):
        d.mkdir(parents=True, exist_ok=True)
    rec = RunRecord(OUT / "run_record.json", "myeloid compartment: sensitivity of the 6 dpi iMON state to Harmony "
                    "on infection round", RULES,
                    notes="rules frozen before Table S3 and the object were opened; labels held out of both embeddings")
    rec.add_output(OUT / "run_record.json")
    rounds = sample_round_table(rec)
    sub = load_subset(rec, rounds)
    embed_both(sub, rec)

    days = sub.obs["day"].to_numpy()
    rnd = sub.obs["round"].astype(str).to_numpy()
    mixing = {rep: round_mixing(sub.obsm[rep], days, rnd) for rep in REPS}
    rows = []
    for d in DAYS_TESTED:
        row = {"day": d}
        for rep, name in REPS.items():
            m = mixing[rep][d]
            row[f"enrichment_{name}"] = m["enrichment"]
            row[f"observed_same_round_fraction_{name}"] = m["observed_same_round_fraction"]
        row["expected_if_mixed"] = mixing["X_pca"][d]["expected_if_mixed"]
        row["n_cells"] = mixing["X_pca"][d]["n_cells"]
        rows.append(row)
    mix_df = pd.DataFrame(rows)
    p = TAB / "round_mixing_by_day.csv"
    mix_df.to_csv(p, index=False)
    rec.add_output(p)
    rec.set("round_mixing", {rep: {str(k): v for k, v in m.items()} for rep, m in mixing.items()})

    grades, checks = {}, {}
    for rep in REPS:
        key = f"leiden_{rep}"
        df, summary = grade(sub, key)
        grades[rep] = {"table": df, **summary}
        p = TAB / f"grading_{rep}.csv"
        df.to_csv(p, index=False)
        rec.add_output(p)
        checks[rep] = imon_state(sub, key)
    from sklearn.metrics import adjusted_rand_score
    ari = float(adjusted_rand_score(sub.obs["leiden_X_pca"].astype(str), sub.obs["leiden_X_pca_harmony"].astype(str)))
    rec.set("adjusted_rand_index_uncorrected_vs_harmony", round(ari, 4))
    rec.set("grading", {rep: {k: v for k, v in g.items() if k != "table"} for rep, g in grades.items()})
    rec.set("imon_state_check", checks)
    frozen_pass = {REPS[rep]: checks[rep]["frozen_rule_most_iMON_cells"]["survives"] for rep in REPS}
    posthoc_pass = {REPS[rep]: checks[rep]["post_hoc_most_6dpi_iMON_cells"]["survives"] for rep in REPS}
    verdict_frozen = (f"frozen rule (selects the 11 to 19 dpi monocyte state, see rule revision): "
                      f"uncorrected {frozen_pass['uncorrected']}, Harmony {frozen_pass['Harmony on infection round']}")
    verdict = ("the 6 dpi iMON state survives Harmony on infection round (post hoc definition, disclosed)"
               if posthoc_pass["Harmony on infection round"] else
               "the 6 dpi iMON state does NOT survive Harmony on infection round (post hoc definition, disclosed)")
    rec.set("verdict_frozen_rule", verdict_frozen)
    rec.set("verdict", verdict)
    pd.DataFrame([{"embedding": REPS[rep], "definition": name,
                   **{k: (v if not isinstance(v, dict) else str(v)) for k, v in c.items() if k != "definition"}}
                  for rep, defs in checks.items() for name, c in defs.items()]).to_csv(TAB / "imon_state_check.csv", index=False)
    rec.add_output(TAB / "imon_state_check.csv")
    cells = sub.obs[["sample_id", "day", "round", "label", "leiden_X_pca", "leiden_X_pca_harmony"]].copy()
    cells.insert(0, "cell_id", sub.obs_names)
    cells.to_csv(TAB / "batch_sensitivity_cell_metadata.csv", index=False)
    rec.add_output(TAB / "batch_sensitivity_cell_metadata.csv")
    ct = pd.crosstab(sub.obs["leiden_X_pca"].astype(str), sub.obs["leiden_X_pca_harmony"].astype(str))
    p = TAB / "subcluster_crosstab_uncorrected_vs_harmony.csv"
    ct.to_csv(p)
    rec.add_output(p)

    # figures
    fig, axes = plt.subplots(2, 3, figsize=(19, 11))
    for i, rep in enumerate(REPS):
        sub.obsm["X_umap"] = sub.obsm[f"X_umap_{rep}"]
        sc.pl.umap(sub, color="round", ax=axes[i, 0], show=False, size=10,
                   title=f"{REPS[rep]}: infection round")
        sc.pl.umap(sub, color="label", ax=axes[i, 1], show=False, size=10, legend_fontsize=7,
                   title=f"{REPS[rep]}: deposited label (held out)")
        sc.pl.umap(sub, color=f"leiden_{rep}", ax=axes[i, 2], show=False, size=10, legend_loc="on data",
                   title=f"{REPS[rep]}: Leiden {PRIMARY_RES}; 6 dpi iMON state = sub "
                         f"{checks[rep]['post_hoc_most_6dpi_iMON_cells']['subcluster']}")
    fig.suptitle(f"Myeloid compartment, {sub.n_obs:,} cells from {DAYS_TESTED} dpi; "
                 f"same-round kNN enrichment overall {mixing['X_pca']['overall']['knn_same_round_enrichment']} "
                 f"(uncorrected) vs {mixing['X_pca_harmony']['overall']['knn_same_round_enrichment']} (Harmony); "
                 f"ARI between partitions {ari:.3f}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    for p in save_fig(fig, FIG, "UMAP_uncorrected_vs_harmony_by_round"):
        rec.add_output(p)

    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(DAYS_TESTED))
    for j, (rep, name) in enumerate(REPS.items()):
        ax.bar(x + (j - 0.5) * 0.36, [mixing[rep][d]["enrichment"] for d in DAYS_TESTED], width=0.34, label=name,
               color=["#999999", "#3182bd"][j])
    ax.axhline(1.0, color="black", lw=0.8, ls="--")
    ax.axhline(2.0, color="#d62728", lw=0.8, ls=":")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d} dpi" for d in DAYS_TESTED])
    ax.set_ylabel("same-round kNN enrichment within day\n(1 = mixed; dotted = pipeline threshold 2)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Do the two infection rounds mix within each day?", fontsize=10)
    fig.tight_layout()
    for p in save_fig(fig, FIG, "round_mixing_by_day"):
        rec.add_output(p)

    # README
    chk_rows = []
    for rep, defs in checks.items():
        for name, c in defs.items():
            chk_rows.append({"embedding": REPS[rep], "definition": c["definition"], "subcluster": c["subcluster"],
                             "n_cells": c["n_cells"], "iMON purity": c["iMON_purity_among_labelled"],
                             "share of each 6 dpi animal's iMON cells": ", ".join(f"{a}: {v}" for a, v in c["fraction_of_each_6dpi_animals_iMON_cells_in_subcluster"].items()),
                             "% cells from 6 dpi": c["pct_cells_from_6dpi"],
                             "round composition": ", ".join(f"{k}: {v}" for k, v in c["round_composition"].items()),
                             "passes criteria": c["survives"]})
    gr_rows = [{"embedding": REPS[rep], **{k: v for k, v in g.items() if k != "table"}} for rep, g in grades.items()]
    lines = [
        "# Batch sensitivity of the myeloid embedding: Harmony on infection round (generated)",
        "",
        "Generated by `analysis/scripts/13_myeloid_batch_sensitivity.py`; frozen rules, inputs and results in",
        "[`run_record.json`](run_record.json). Owner retain/reject review pending.",
        "",
        "**Question.** Is the dense 6 dpi inflammatory-monocyte (iMON) state of `../README.md` biology or a",
        "one-day batch island? Correcting on sample cannot answer it (one animal per sample and day). The",
        "paper's Table S3 gives the infection round, and on every tested day the two replicate animals came",
        "from different rounds, so the round is a technical key that crosses time. It also carries Ki67-Cre",
        "dosage (Cre/Cre in most 2021 animals, Cre/+ in all 2022 animals). The compartment was embedded twice",
        "with identical rules, uncorrected and after Harmony on round, and compared.",
        "",
        f"Cells: {rec.record['results']['n_cells']:,} from {rec.record['results']['n_animals']} animals at {DAYS_TESTED} dpi "
        f"(0, 90 and 366 dpi are single-round and excluded). Sample-to-round map: `tables/sample_infection_round.csv`.",
        "",
        "## 1. Do the rounds mix within each day",
        "",
        "![Round mixing by day](figures/round_mixing_by_day.png)",
        "",
        df_to_markdown(mix_df, index=False),
        "",
        f"Overall same-round kNN enrichment: {mixing['X_pca']['overall']['knn_same_round_enrichment']} uncorrected, "
        f"{mixing['X_pca_harmony']['overall']['knn_same_round_enrichment']} after Harmony (1 = mixed; the pipeline's "
        "failure threshold is 2).",
        "",
        "## 2. Does the 6 dpi iMON state survive",
        "",
        "![UMAP uncorrected vs Harmony](figures/UMAP_uncorrected_vs_harmony_by_round.png)",
        "",
        df_to_markdown(pd.DataFrame(chk_rows), index=False),
        "",
        "**Rule revision, disclosed.** The frozen definition (subcluster with the most iMON-labelled cells)",
        "selects the 11 to 19 dpi monocyte state in both embeddings, not the 6 dpi state the question is",
        "about; its outcome is kept in the run record (`first_run_outcome`, `verdict_frozen_rule`). The",
        "post hoc definition anchored on the 6 dpi cells is reported alongside and judged by the same",
        "purity and per-animal criteria.",
        "",
        f"**Verdict: {verdict}.** {verdict_frozen}. Adjusted Rand index between the two Leiden "
        f"partitions: {ari:.3f}.",
        "",
        "## 3. Grading against the held-out labels in both embeddings",
        "",
        df_to_markdown(pd.DataFrame(gr_rows), index=False),
        "",
        "Per-subcluster tables: `tables/grading_X_pca.csv`, `tables/grading_X_pca_harmony.csv`; the",
        "partition crosstab: `tables/subcluster_crosstab_uncorrected_vs_harmony.csv`; per-cell subcluster",
        "assignments in both embeddings: `tables/batch_sensitivity_cell_metadata.csv`.",
        "",
        "## Caveats",
        "",
        "- Infection round and Ki67-Cre dosage are confounded; the correction removes both, which is the",
        "  conservative direction for this question.",
        "- Days with one round (0, 90, 366 dpi) are not tested; the persistence claims at 90 and 366 dpi",
        "  rest on the uncorrected embedding and on within-day replicate agreement only.",
        "- Two animals per active-repair day; descriptive, no P values.",
        "",
    ]
    p = OUT / "README.md"
    p.write_text("\n".join(lines), encoding="utf-8")
    rec.add_output(p)
    rec.finish()
    print(mix_df.to_string(index=False))
    print(pd.DataFrame(chk_rows).to_string(index=False))
    print(f"\n{verdict}; {verdict_frozen}; ARI {ari:.3f}")
    print(f"wrote {OUT.relative_to(ANALYSIS.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
