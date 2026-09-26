"""Post-transfer, targeted intestinal depth check; never retunes the candidate."""
import json
import numpy as np
import pandas as pd
from exploratory_common import BASE, PROC, TABLES, save_json


def main():
    config=json.loads((BASE/'exploratory_config.json').read_text())
    states=config['transfer']['intestinal_states']
    obs=pd.read_csv(PROC/'intestine_scored_obs.csv')
    scores=np.load(PROC/'intestine_transfer_scores.npz')
    modules=['candidate','candidate_without_generic_controls','ADI_published_holdout']
    rng=np.random.default_rng(0)
    rows=[];quality=[]
    for unit,f in obs[obs.state.isin(states)].groupby('unit'):
        edges=np.unique(np.quantile(np.log1p(f.total_umi),np.linspace(0,1,6)))
        bins=pd.cut(np.log1p(f.total_umi),edges,include_lowest=True,labels=False)
        strata=[]
        for b in sorted(bins.dropna().unique()):
            within=f.loc[bins==b]
            groups=[within.index[within.state==s].to_numpy() for s in states]
            n=min(map(len,groups))
            if n:strata.append((groups,n))
        nmatched=sum(n for _,n in strata)
        assert nmatched>0
        for iteration in range(200):
            chosen=[[] for _ in states]
            for groups,n in strata:
                for i,group in enumerate(groups):chosen[i].extend(rng.choice(group,n,replace=False).tolist())
            for i,s in enumerate(states):
                quality.append({'unit':unit,'iteration':iteration,'state':s,'matched_cells':len(chosen[i]),
                    'mean_log1p_umi':float(np.log1p(obs.loc[chosen[i],'total_umi']).mean())})
            for module in modules:
                for endpoint,index in [(states[0],0),(states[2],2)]:
                    delta=scores[module][chosen[1]].mean()-scores[module][chosen[index]].mean()
                    rows.append({'unit':unit,'iteration':iteration,'module':module,'endpoint':endpoint,
                        'matched_cells_per_state':nmatched,'all_states_at_least_30':nmatched>=30,'difference':float(delta)})
    draws=pd.DataFrame(rows)
    draws.to_csv(PROC/'intestine_depth_matching_draws.csv',index=False)
    quality=pd.DataFrame(quality)
    quality.to_csv(PROC/'intestine_depth_matching_quality_draws.csv',index=False)
    summary=draws.groupby(['unit','module','endpoint']).agg(matched_cells_per_state=('matched_cells_per_state','first'),
        all_states_at_least_30=('all_states_at_least_30','first'),median_difference=('difference','median'),
        draw_q05=('difference',lambda x:x.quantile(.05)),draw_q95=('difference',lambda x:x.quantile(.95)),
        fraction_draws_positive=('difference',lambda x:(x>0).mean())).reset_index()
    summary.to_csv(TABLES/'intestine_depth_matched_sensitivity.csv',index=False)
    quality.groupby(['unit','state']).agg(matched_cells=('matched_cells','first'),
        median_mean_log1p_umi=('mean_log1p_umi','median')).reset_index().to_csv(TABLES/'intestine_depth_matching_quality.csv',index=False)
    result={'stage':'E3a','status':'targeted_post_transfer_sensitivity_completed',
        'reason':'Primary intestinal transfer is near zero while the pre-transfer generic-excluded variant is positive; check whether depth explains the discrepancy before investment decision.',
        'program_retuned':False,'new_gene_selection':False,'sampling_iterations':200,
        'interpretation':'Sampling quantiles describe dependence on depth-matching draws, not biological confidence intervals; mice remain the biological units.',
        'matched_cells_per_mouse':summary.drop_duplicates('unit')[['unit','matched_cells_per_state']].to_dict('records'),
        'qualified_after_matching':int(summary.drop_duplicates('unit').all_states_at_least_30.sum()),
        'summary_table':'tables/exploratory/intestine_depth_matched_sensitivity.csv'}
    save_json(BASE/'stage_E3a_result.json',result)
    print(summary.to_string(index=False),flush=True)


if __name__=='__main__':main()
