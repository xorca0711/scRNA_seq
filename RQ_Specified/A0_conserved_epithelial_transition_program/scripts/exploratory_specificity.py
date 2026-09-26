"""E2: descriptive specificity and technical challenges to the frozen repair program."""
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from exploratory_common import BASE, TABLES, PROC, save_json, load_repair, rank_scores, paired_effects, excluded


def main():
    matrix, all_obs, genes, config = load_repair()
    program = json.loads((BASE / "frozen_repair_program.json").read_text())["genes"]
    rng = np.random.default_rng(config["seed"])
    mask = all_obs.state.isin(["AT2","Krt8+ ADI","AT1"]).to_numpy()
    x, obs = matrix[mask], all_obs.loc[mask].reset_index(drop=True)
    score, _ = rank_scores(x,genes,program,config["reference_genes"])
    effects = paired_effects(obs,score)
    medians = effects.groupby("endpoint").difference.median()
    core_scores = {"candidate":score}
    modules = config["control_modules"]
    generic = set().union(*(set(modules[n]) for n in config["generic_control_names"]))
    reduced = [g for g in program if g not in generic]
    save_json(BASE / "frozen_control_variant.json", {"genes":reduced,
        "definition":"Frozen repair program minus prespecified generic-control union",
        "parent_program":"frozen_repair_program.json", "not_a_replacement_primary_program":True,
        "transfer_outcomes_seen":False, "generic_overlap_genes":sorted(set(program)&generic)})
    benchmark_effects, coverage = [], []
    for name,module in modules.items():
        s,c = rank_scores(x,genes,module,config["reference_genes"])
        core_scores[name] = s
        e = paired_effects(obs,s);e["module"] = name
        benchmark_effects.append(e)
        coverage.append({"module":name,**c})
    if len(reduced) >= 10:
        s,c = rank_scores(x,genes,reduced,config["reference_genes"])
        core_scores["candidate_without_generic_controls"] = s
        e=paired_effects(obs,s);e["module"]="candidate_without_generic_controls"
        benchmark_effects.append(e);coverage.append({"module":"candidate_without_generic_controls",**c})
    pd.concat(benchmark_effects,ignore_index=True).to_csv(TABLES/"repair_control_effects.csv",index=False)
    pd.DataFrame(coverage).to_csv(TABLES/"repair_control_coverage.csv",index=False)
    overlap = pd.DataFrame({"gene":program,"in_generic_union":[g in generic for g in program]})
    for n in config["generic_control_names"]:
        overlap[n] = overlap.gene.isin(modules[n])
    overlap.to_csv(TABLES/"program_control_overlap.csv",index=False)
    print("Control modules scored; beginning technical challenges.",flush=True)

    # Match without replacement within mouse and pooled UMI quintile.
    matched_indexes=[]
    for unit, f in obs.groupby("unit"):
        edges=np.unique(np.quantile(np.log1p(f.measured_umi),np.linspace(0,1,6)))
        bins=pd.cut(np.log1p(f.measured_umi),edges,include_lowest=True,labels=False)
        for b in sorted(bins.dropna().unique()):
            within=f.loc[bins==b]
            groups=[within.index[within.state==state].to_numpy() for state in ["AT2","Krt8+ ADI","AT1"]]
            n=min(map(len,groups))
            if n:
                for group in groups:matched_indexes.extend(rng.choice(group,n,replace=False).tolist())
    matched_indexes=np.sort(matched_indexes)
    depth=paired_effects(obs.loc[matched_indexes],score[matched_indexes],dataset="repair_depth_matched")
    depth.to_csv(TABLES/"repair_depth_matched_effects.csv",index=False)
    obs.loc[matched_indexes,["cell_barcode","unit","state","measured_umi"]].to_csv(PROC/"depth_matched_cells.csv",index=False)
    correlations=[]
    for (unit,state),f in obs.groupby(["unit","state"]):
        for metric in ["measured_umi","measured_genes","measured_mito_fraction"]:
            correlations.append({"unit":unit,"state":state,"metric":metric,"n_cells":len(f),
                "spearman_rho":float(spearmanr(score[f.index],f[metric]).statistic)})
    pd.DataFrame(correlations).to_csv(TABLES/"repair_score_quality_correlations.csv",index=False)

    cycle=(core_scores["HALLMARK_E2F_TARGETS"]+core_scores["HALLMARK_G2M_CHECKPOINT"])/2
    low_cycle=np.zeros(len(obs),bool)
    for unit,f in obs.groupby("unit"):
        low_cycle[f.index]=cycle[f.index] <= np.median(cycle[f.index])
    cycle_effects=paired_effects(obs.loc[low_cycle],score[low_cycle],dataset="repair_low_cycle")
    cycle_effects.to_csv(TABLES/"repair_low_cycle_effects.csv",index=False)

    # Synthetic normalized endpoint mixtures are measurement controls, never animals.
    mix_rows=[]
    for unit,f in obs.groupby("unit"):
        ai=rng.choice(f.index[f.state=="AT2"],200,replace=True)
        di=rng.choice(f.index[f.state=="AT1"],200,replace=True)
        adi_mean=float(score[f.index[f.state=="Krt8+ ADI"]].mean())
        for weight in [.25,.5,.75]:
            synthetic=(x[ai].multiply((weight/obs.measured_umi.to_numpy()[ai])[:,None])
                + x[di].multiply(((1-weight)/obs.measured_umi.to_numpy()[di])[:,None])).tocsr()
            s,_=rank_scores(synthetic,genes,program,config["reference_genes"])
            mix_rows.append({"unit":unit,"AT2_weight":weight,"synthetic_pairs":200,
                "mixture_mean":float(s.mean()),"intermediate_mean":adi_mean,"intermediate_minus_mixture":adi_mean-float(s.mean())})
    mixture=pd.DataFrame(mix_rows)
    mixture.to_csv(TABLES/"repair_endpoint_mixture_checks.csv",index=False)

    # Activated AT2 is a contextual comparator, not a validated stress-only state.
    activated_mask=all_obs.state.isin(["AT2 activated","Krt8+ ADI","AT1"]).to_numpy()
    activated_score,_=rank_scores(matrix[activated_mask],genes,program,config["reference_genes"])
    activated=paired_effects(all_obs.loc[activated_mask],activated_score,start="AT2 activated",dataset="repair_activated_AT2_comparator")
    activated.to_csv(TABLES/"repair_activated_AT2_effects.csv",index=False)
    print("Depth, cycle, mixture and activated-AT2 checks complete; beginning 200 matched gene sets.",flush=True)

    pb=np.load(PROC/"repair_pseudobulk.npz")["logcpm"].mean(axis=0)
    bins=pd.qcut(pd.Series(pb),10,labels=False,duplicates="drop").to_numpy()
    index={g:i for i,g in enumerate(genes)}
    allowed=[g for g in genes if g not in program and not excluded(g,config)]
    pools={b:[g for g in allowed if bins[index[g]]==b] for b in np.unique(bins)}
    random_rows=[];random_members=[]
    for iteration in range(config["random_gene_sets"]):
        matched=[]
        for g in program:
            pool=[v for v in pools[bins[index[g]]] if v not in matched]
            assert pool
            matched.append(str(rng.choice(pool)))
        s,_=rank_scores(x,genes,matched,config["reference_genes"])
        e=paired_effects(obs,s)
        for endpoint,frame in e.groupby("endpoint"):
            random_rows.append({"random_set":iteration,"endpoint":endpoint,"median_difference":float(frame.difference.median())})
        random_members += [{"random_set":iteration,"gene":g} for g in matched]
        if (iteration+1)%50==0:print(f"Matched sets scored: {iteration+1}/200",flush=True)
    random=pd.DataFrame(random_rows)
    random.to_csv(TABLES/"repair_matched_random_set_effects.csv",index=False)
    pd.DataFrame(random_members).to_csv(PROC/"matched_random_gene_membership.csv",index=False)
    random_summary={e:{"candidate_median_difference":float(medians[e]),
        "random_median_difference":float(f.median_difference.median()),
        "fraction_random_at_least_candidate":float((f.median_difference >= medians[e]).mean())} for e,f in random.groupby("endpoint")}
    depth_summary={e:{"mice_with_overlap":len(f),"positive_mice":int((f.difference>0).sum()),
        "median_difference":float(f.difference.median()),"qualified_matched_mice":int(f.all_states_at_least_30.sum()),
        "minimum_matched_cells_per_state":int(f.n_intermediate.min())} for e,f in depth.groupby("endpoint")}
    mixture_worst=mixture.groupby("unit").intermediate_minus_mixture.min()
    pass_technical=all(v["median_difference"]>0 and v["positive_mice"]>=6 for v in depth_summary.values()) and float(mixture_worst.median())>0
    result={"stage":"E2","generic_overlap_genes":len(set(program)&generic),"non_generic_candidate_genes":len(reduced),
        "depth_matching":depth_summary,"matched_random_sets":random_summary,
        "mixture_check":{"positive_mice_against_all_three_mixture_weights":int((mixture_worst>0).sum()),
            "median_worst_intermediate_minus_mixture":float(mixture_worst.median())},
        "checkpoint_passed":bool(pass_technical),"next_stage_recommended":"E3 descriptive transfer" if pass_technical else "stop_or_narrow_technical_signal",
        "limitations":["Matched random sets use training-context data; tail fractions are not inferential p-values",
            "Control membership overlap is incomplete knowledge, not a clean separation of stress biology",
            "Depth matching and low-cycle restriction may reduce cell counts below the original floor",
            "Synthetic endpoint mixtures do not exhaust real doublets or ambient RNA",
            "No independent cell lineage validation or overparameterized donor-adjusted model"]}
    save_json(BASE/"stage_E2_result.json",result)
    np.savez_compressed(PROC/"repair_cell_scores.npz",**core_scores)
    obs.to_csv(PROC/"repair_scored_obs.csv",index=False)
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":
    main()
