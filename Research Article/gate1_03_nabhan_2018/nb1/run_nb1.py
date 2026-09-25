"""Animal-level descriptive Wnt analysis; see PROTOCOL.md before interpreting."""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy import sparse
from scipy.special import gammaln
from anndata._io.specs import read_elem
from adjustText import adjust_text
from importlib.metadata import version

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "Research Article/gate1_01_niethamer_2025/GSE262927"
MATRIX = DATA / "processed/final_clustered.h5ad"
METADATA = DATA / "tables/cell_metadata.csv"
ANIMALS = DATA / "myeloid_focus/batch_sensitivity/tables/sample_infection_round.csv"
MODULES = ROOT / "Research Article/epithelial_state_specificity/modules.json"
TYPES = ["AT2", "AT1_AT2", "AT1", "AF1", "AF2", "Adventitial_fibroblast", "Peribronchial_fibroblast"]
WNTS = ["Wnt1", "Wnt2", "Wnt2b", "Wnt3", "Wnt3a", "Wnt4", "Wnt5a", "Wnt5b", "Wnt6", "Wnt7a", "Wnt7b", "Wnt8a", "Wnt8b", "Wnt9a", "Wnt9b", "Wnt10a", "Wnt10b", "Wnt11", "Wnt16"]
PANEL = WNTS + ["Porcn", "Wls", "Axin2", "Lef1", "Mki67", "Top2a"]
BUDGET = 1000


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def expected_detection(total, target, budget=BUDGET):
    """Exact E[detected] after sampling budget UMIs without replacement."""
    n, k = np.broadcast_arrays(np.asarray(total, float), np.asarray(target, float))
    if np.any(k < 0) or np.any(k > n) or np.any(n < budget):
        raise ValueError("Invalid counts or library below sampling budget")
    out = np.zeros(n.shape)
    sure = (k > 0) & (n - k < budget)
    out[sure] = 1
    valid = (k > 0) & ~sure
    nv, kv = n[valid], k[valid]
    lp = gammaln(nv-kv+1) + gammaln(nv-budget+1) - gammaln(nv-kv-budget+1) - gammaln(nv+1)
    out[valid] = np.clip(-np.expm1(lp), 0, 1)
    return out


def savefig(fig, name):
    for ext in ["png", "svg"]:
        fig.savefig(HERE / "figures" / f"{name}.{ext}", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plots(units, expr, scores):
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False,
                         "axes.spines.right": False, "svg.hashsalt": "nb1-20260922"})
    order = units[["sample_id", "day"]].drop_duplicates().sort_values(["day", "sample_id"])
    counts = units.pivot(index="sample_id", columns="cell_type", values="cells").reindex(index=order.sample_id, columns=TYPES).fillna(0)
    fig, ax = plt.subplots(figsize=(10, 9))
    im = ax.imshow(np.log10(counts.to_numpy()+1), aspect="auto", cmap="Blues")
    ax.set(xticks=np.arange(len(TYPES)), xticklabels=[t.replace("_fibroblast", "\nfibroblast") for t in TYPES],
           yticks=np.arange(len(order)), yticklabels=[f"{r.sample_id} | day {r.day:g}" for r in order.itertuples()],
           title="Nb1 | Animals and compartment coverage\nBold counts pass the primary 50-cell floor")
    for i in range(len(counts)):
        for j in range(len(TYPES)):
            n = int(counts.iloc[i,j])
            ax.text(j, i, str(n), ha="center", va="center", fontsize=7,
                    fontweight="bold" if n >= 50 else "normal", color="white" if n > 1000 else "black")
    fig.colorbar(im, ax=ax, label="log10(cells + 1)", shrink=.6)
    savefig(fig, "01_animal_coverage")

    genes = ["Wnt5a", "Wnt2", "Wnt2b", "Wnt4", "Wnt9a", "Wnt7b", "Porcn", "Wls", "Axin2", "Lef1"]
    eligible = expr.loc[expr.eligible_50 & expr.gene.isin(genes)].copy()
    eligible["log2_cpm"] = np.log2(eligible.cpm + 1)
    summary = eligible.groupby(["cell_type", "day", "gene"], observed=True).agg(
        mean_log2_cpm=("log2_cpm", "mean"), animals=("sample_id", "nunique"),
        mean_detection_1000=("expected_detection_1000", "mean")).reset_index()
    summary.to_csv(HERE / "tables/figure_source_summary.csv", index=False)
    row_index = pd.MultiIndex.from_tuples([(t,d) for t in TYPES for d in sorted(summary.day.unique())
                                          if ((summary.cell_type==t)&(summary.day==d)).any()], names=["cell_type","day"])
    fig, axs = plt.subplots(1,2,figsize=(13, max(9, len(row_index)*.22)), sharey=True)
    for ax, value, title, cmap, vmax in zip(axs, ["mean_log2_cpm","mean_detection_1000"],
            ["Mean across animals: log2(CPM + 1)", "Mean across animals: expected detection at 1,000 UMIs"],
            ["viridis", "magma"], [None, 1]):
        a = summary.pivot(index=["cell_type","day"],columns="gene",values=value).reindex(index=row_index,columns=genes)
        im = ax.imshow(a.to_numpy(), aspect="auto", cmap=cmap, vmin=0, vmax=vmax)
        ax.set(xticks=np.arange(len(genes)),xticklabels=genes,title=title)
        ax.tick_params(axis="x",rotation=70)
        fig.colorbar(im,ax=ax,shrink=.35)
    ns = summary.groupby(["cell_type","day"]).animals.first()
    axs[0].set(yticks=np.arange(len(row_index)),yticklabels=[f"{t} | d{d:g} | n={ns.loc[t,d]}" for t,d in row_index])
    fig.suptitle("Nb1 | Wnt source and response by deposited compartment\nOnly animal × compartment units with ≥50 cells; descriptive averages; unequal sampling across days", y=1.01)
    savefig(fig, "02_source_response")

    at2 = expr.loc[(expr.cell_type=="AT2") & expr.gene.isin(["Wnt7b","Porcn","Wls","Axin2","Lef1","Mki67"])]
    days = sorted(at2.day.unique())
    daypos = {d:i for i,d in enumerate(days)}
    jitter = {s:(i-(len(g)-1)/2)*.065 for _,g in order.groupby("day") for i,s in enumerate(g.sample_id)}
    fig, axs = plt.subplots(2,3,figsize=(12,7),sharex=True)
    for ax, gene in zip(axs.flat,["Wnt7b","Porcn","Wls","Axin2","Lef1","Mki67"]):
        part = at2.loc[at2.gene==gene]
        for r in part.itertuples():
            ax.scatter(daypos[r.day]+jitter[r.sample_id], np.log2(r.cpm+1),
                       facecolors="#146b8a" if r.eligible_50 else "none", edgecolors="#146b8a",s=35)
        ax.set(title=gene,ylabel="log2(CPM + 1)",xticks=range(len(days)),xticklabels=[f"{d:g}" for d in days])
    for ax in axs[-1]:
        ax.set_xlabel("Days after infection (categorical spacing)")
    fig.suptitle("Nb1 | AT2 expression, one point per animal\nFilled: ≥50 cells; open: below floor; no time-course test or inferred autocrine signaling")
    fig.tight_layout()
    savefig(fig,"03_at2_per_animal")

    part = scores.loc[scores.cell_type=="AT2"]
    fig, axs = plt.subplots(1,2,figsize=(11,4.8))
    norm = matplotlib.colors.Normalize(vmin=0,vmax=366)
    for ax, x, label in zip(axs,["at2_identity_mean_log2_cpm","proliferation_mean_log2_cpm"],
            ["AT2 published holdout: mean log2(CPM + 1)","Mki67 / Top2a: mean log2(CPM + 1)"]):
        labels=[]
        for r in part.itertuples():
            color=plt.cm.viridis(norm(r.day))
            ax.scatter(getattr(r,x),r.response_mean_log2_cpm,facecolors=[color] if r.eligible_50 else "none",edgecolors=[color],s=45)
            labels.append(ax.text(getattr(r,x),r.response_mean_log2_cpm,r.sample_id.rsplit("-",1)[-1],fontsize=6))
        ax.set(xlabel=label,ylabel="Axin2 / Lef1: mean log2(CPM + 1)")
        np.random.seed(20260922)
        adjust_text(labels,ax=ax,iter_lim=500,arrowprops={"arrowstyle":"-","color":".5","lw":.4})
    fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap="viridis"),ax=axs,label="Day",shrink=.7)
    fig.suptitle("Nb1 | AT2 identity, proliferation and Wnt-response transcripts\nAnimal labels shown; open points below 50 cells; confounding prevents causal interpretation")
    savefig(fig,"04_at2_dimensions")

    # Every paper-defined fibroblast ligand is shown; no outcome-selected subset.
    fibro=expr.loc[expr.eligible_50 & expr.cell_type.isin(["AF1","AF2"])]
    fig,axs=plt.subplots(1,5,figsize=(13,4),sharey=True)
    pairs=[]
    for ax,gene in zip(axs,["Wnt5a","Wnt2","Wnt2b","Wnt4","Wnt9a"]):
        one=fibro.loc[fibro.gene==gene]
        for endpoint in ["cpm","expected_detection_1000"]:
            wide=one.pivot(index="sample_id",columns="cell_type",values=endpoint).dropna()
            for sample,row in wide.iterrows():
                pairs.append({"sample_id":sample,"gene":gene,"endpoint":endpoint,"AF1":row.AF1,"AF2":row.AF2,
                              "direction":"AF1>AF2" if row.AF1>row.AF2 else "AF2>AF1" if row.AF2>row.AF1 else "tie"})
                if endpoint=="cpm":
                    ax.plot([0,1],np.log2(row.to_numpy()+1),color="#407a93",alpha=.35,lw=.7,marker="o",markersize=2)
        ax.set(title=gene,xticks=[0,1],xticklabels=["AF1","AF2"],xlim=(-.3,1.3))
    axs[0].set_ylabel("log2(CPM + 1)")
    fig.suptitle("Nb1 | Paper-defined Wnt ligands in paired alveolar fibroblast compartments\nEach line joins the same animal; both compartments ≥50 cells; descriptive comparisons")
    fig.tight_layout()
    savefig(fig,"05_fibroblast_pairs")
    pd.DataFrame(pairs).to_csv(HERE/"tables/paired_fibroblast_ligands.csv",index=False)


def main():
    for d in ["tables","figures"]:
        (HERE/d).mkdir(exist_ok=True)
    module = next(m for m in json.loads(MODULES.read_text())["modules"] if m["name"]=="AT2_published_holdout")
    identity = module["genes"]
    genes = list(dict.fromkeys(PANEL + identity))
    with h5py.File(MATRIX) as f:
        obs, var = read_elem(f["obs"]), read_elem(f["var"])
    meta = pd.read_csv(METADATA,low_memory=False).set_index("cell_id")
    assert meta.index.is_unique and obs.index.is_unique
    assert set(meta.index)==set(obs.index), "Matrix and metadata have different cells"
    meta = meta.loc[obs.index].copy()
    for col in ["sample_id","author_celltype"]:
        # The H5AD uses literal "NA" for absent author labels; CSV parsing uses NaN.
        left = meta[col].astype(object).fillna("").astype(str).replace("NA", "")
        right = obs[col].astype(object).fillna("").astype(str).replace("NA", "")
        assert np.array_equal(left,right), f"Metadata mismatch: {col}"
    meta["day"] = pd.to_numeric(meta.sacrifice_day,errors="coerce")
    selected = meta.author_celltype.isin(TYPES) & meta.day.notna() & meta.sample_id.notna()
    assert meta.loc[selected,"has_author_metadata"].astype(str).str.lower().isin(["true","1"]).all()
    cells = meta.loc[selected].copy()
    cells["cell_type"] = cells.author_celltype
    units = cells.groupby(["sample_id","day","cell_type"],observed=True).size().rename("cells").reset_index().sort_values(["day","sample_id","cell_type"]).reset_index(drop=True)
    assert not units.duplicated(["sample_id","cell_type"]).any(), "Animal has multiple days"
    animals = pd.read_csv(ANIMALS)
    assert animals.sample_id.is_unique
    units = units.merge(animals[["sample_id","genotype","sex","round"]],on="sample_id",how="left",validate="many_to_one")
    assert units[["genotype","sex","round"]].notna().all().all(), "Missing animal covariates"
    units["unit_id"] = [f"{s}|{t}" for s,t in zip(units.sample_id,units.cell_type)]
    units["heterozygous"] = units.genotype.str.contains("Ki67Cre/+",regex=False)
    lookup = {u:i for i,u in enumerate(units.unit_id)}
    group = np.full(len(meta),-1,int)
    group[selected] = [lookup[f"{s}|{t}"] for s,t in zip(cells.sample_id,cells.cell_type)]
    symbols = var.gene_symbol.astype(str).to_numpy()
    gene_lookup = {g:i for i,g in enumerate(genes)}
    feature_rows = np.array([i for i,g in enumerate(symbols) if g in gene_lookup])
    feature_cols = np.array([gene_lookup[symbols[i]] for i in feature_rows])
    mapping = sparse.csr_matrix((np.ones(len(feature_rows),dtype=np.int64),(feature_rows,feature_cols)),shape=(len(var),len(genes)))
    present = np.array([g in set(symbols) for g in genes])
    n_units = len(units)
    pb = np.zeros((n_units,len(genes)),dtype=np.int64)
    positive = np.zeros((n_units,len(PANEL)),dtype=np.int64)
    expected = np.zeros_like(positive,dtype=float)
    library = np.zeros(n_units,dtype=np.int64)
    metadata_library = np.zeros(n_units,dtype=np.int64)
    budget_cells = np.zeros(n_units,dtype=np.int64)
    seen = np.zeros(n_units,dtype=np.int64)
    count_match = 0
    with h5py.File(MATRIX) as f:
        layer = f["layers/counts"]
        assert layer.attrs["encoding-type"]=="csr_matrix"
        assert tuple(layer.attrs["shape"])==(len(meta),len(var))
        ptr = layer["indptr"][:]
        for start in range(0,len(meta),4000):
            stop = min(start+4000,len(meta))
            take = np.flatnonzero(group[start:stop]>=0)
            if not len(take):
                continue
            lo,hi = int(ptr[start]),int(ptr[stop])
            data = layer["data"][lo:hi]
            assert np.isfinite(data).all() and (data>=0).all() and (data==np.floor(data)).all(), "Noninteger raw counts"
            block = sparse.csr_matrix((data.astype(np.int64),layer["indices"][lo:hi],ptr[start:stop+1]-lo),shape=(stop-start,len(var)))[take]
            total = np.asarray(block.sum(axis=1)).ravel()
            assert (total>0).all()
            target = (block@mapping).toarray()
            assert (target.sum(axis=1)<=total).all()
            groups = group[start:stop][take]
            np.add.at(pb,groups,target)
            np.add.at(library,groups,total)
            np.add.at(seen,groups,1)
            np.add.at(positive,groups,target[:,:len(PANEL)]>0)
            enough = total>=BUDGET
            np.add.at(budget_cells,groups[enough],1)
            np.add.at(expected,groups[enough],expected_detection(total[enough,None],target[enough,:len(PANEL)]))
            original = pd.to_numeric(meta.iloc[start:stop].iloc[take].total_counts).to_numpy()
            assert np.isfinite(original).all() and (original==np.floor(original)).all()
            np.add.at(metadata_library,groups,original.astype(np.int64))
            count_match += int((total==original).sum())
    assert np.array_equal(seen,units.cells.to_numpy())
    assert (positive<=seen[:,None]).all()
    assert (expected<=budget_cells[:,None]+1e-6).all()
    units["library_umis"] = library
    units["metadata_pre_gene_filter_umis"] = metadata_library
    units["library_minus_metadata_umis"] = library-metadata_library
    units["cells_ge_1000_umis"] = budget_cells
    for floor in [20,50,100]:
        units[f"eligible_{floor}"] = units.cells>=floor
    units.to_csv(HERE/"tables/animal_compartment_coverage.csv",index=False)
    coverage=[]
    for subset, frame in [("all_genotypes",units),("heterozygous_only",units.loc[units.heterozygous])]:
        for floor in [20,50,100]:
            for (day,typ), frame2 in frame.groupby(["day","cell_type"],observed=True):
                coverage.append(dict(subset=subset,floor=floor,day=day,cell_type=typ,animals_present=len(frame2),animals_eligible=int((frame2.cells>=floor).sum())))
    pd.DataFrame(coverage).to_csv(HERE/"tables/eligibility_sensitivity.csv",index=False)
    cpm = pb/library[:,None]*1e6
    records=[]
    for i,r in enumerate(units.to_dict("records")):
        for j,gene in enumerate(PANEL):
            records.append({**r,"gene":gene,"feature_present":bool(present[j]),
                "raw_count":int(pb[i,j]) if present[j] else np.nan,
                "cpm":cpm[i,j] if present[j] else np.nan,
                "positive_cells":int(positive[i,j]) if present[j] else np.nan,
                "positive_fraction":positive[i,j]/seen[i] if present[j] else np.nan,
                "expected_detection_1000":expected[i,j]/budget_cells[i] if present[j] and budget_cells[i] else np.nan})
    expr = pd.DataFrame(records)
    expr.to_csv(HERE/"tables/gene_expression_by_animal.csv",index=False)
    scores = units.copy()
    for name,panel in [("at2_identity",identity),("response",["Axin2","Lef1"]),("proliferation",["Mki67","Top2a"])]:
        columns = [gene_lookup[g] for g in panel if present[gene_lookup[g]]]
        assert columns, f"No measured genes for {name}"
        scores[f"{name}_genes_present"] = len(columns)
        scores[f"{name}_genes_requested"] = len(panel)
        scores[f"{name}_mean_log2_cpm"] = np.log2(cpm[:,columns]+1).mean(axis=1)
    scores.to_csv(HERE/"tables/animal_marker_summaries.csv",index=False)
    pd.DataFrame({"gene":genes,"present":present,"input_feature_count":[int((symbols==g).sum()) for g in genes]}).to_csv(HERE/"tables/panel_coverage.csv",index=False)
    plots(units,expr,scores)
    inputs=[MATRIX,METADATA,ANIMALS,MODULES,HERE/"PROTOCOL.md",Path(__file__)]
    outputs=sorted((HERE/"tables").glob("*.csv"))+sorted((HERE/"figures").glob("*"))
    record={"analysis":"Nb1","date":"2026-09-22","status":"descriptive_only",
        "inputs":[{"path":p.relative_to(ROOT).as_posix(),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in inputs],
        "outputs":[{"path":p.relative_to(HERE).as_posix(),"sha256":sha256(p)} for p in outputs],
        "dimensions":{"input_cells":len(meta),"input_features":len(var),"selected_cells":len(cells),"units":len(units),"animals":units.sample_id.nunique()},
        "checks":{"metadata_alignment":"exact cell IDs and sample/author labels (NA sentinel normalized)","raw_count_integrity":"finite nonnegative integers", "unit_cell_totals":"exact", "detection_bounds":"passed", "selected_cells_matching_metadata_total_counts_exactly":count_match,
                  "current_layer_minus_metadata_umis":int((library-metadata_library).sum()),"max_unit_relative_loss":float(((metadata_library-library)/metadata_library).max())},
        "identity_panel_provenance":module,"environment":{"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__,"h5py":h5py.__version__,"matplotlib":matplotlib.__version__,"adjustText":version("adjustText")}}
    (HERE/"run_record.json").write_text(json.dumps(record,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"complete","dimensions":record["dimensions"],"checks":record["checks"]}))


if __name__=="__main__":
    main()
