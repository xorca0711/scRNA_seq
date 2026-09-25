"""Expression-supported planned-ligand rankings and donor-omission stability.

Focused triad and source-agnostic views use the same predeclared ligand family
panel, not the entire ligandomics space. Unmeasured genes are not zero-imputed.
"""
from pathlib import Path
import sys,json
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]
from u5_liana_robustness import LABELS,DIRECTIONS

def main():
    import numpy as np,pandas as pd
    out=PAPER/'trials/u5_ligand_targets';support=[];rankings=[];eligibility=[];stability=[];edge_coverage=[];expected_omissions={}
    candidates=['IL1A','IL1B','AREG','HBEGF','TGFB1','CCL2','CXCL12','CXCL1','CXCL2','FGF7','FGF10','VEGFA']
    for cohort in ['GSE136831','GSE135893']:
        here=out/cohort;state=json.loads((here/'run_record.json').read_text());assert state['status']=='completed_refits_candidate_filtering_pending'
        scores=pd.read_csv(here/'all_prior_pearson_scores.csv.gz');e=pd.read_csv(here/'target_eligibility.csv');units=pd.read_csv(here/'receiver_units.csv')
        for (comp,label),u in units.groupby(['compartment','label']):expected_omissions[(cohort,comp,label)]=u.donor.nunique()
        prior_lr=pd.read_csv(here/'prior_ligand_receptor.csv')[['from','to']].drop_duplicates()
        p=pd.read_csv(PAPER/'trials/u5_liana_robustness'/cohort/'all_label_source_recipient_profiles.csv.gz')
        p=p[~p.label.str.contains('multiplet|doublet|unknown|unassigned',case=False,regex=True)].copy()
        ligands=sorted((set(candidates)|{x for x in p.gene.unique() if x.startswith('WNT')})&set(scores.ligand))
        measured=set(p.gene);labels=LABELS[cohort]
        for edge in prior_lr[prior_lr['from'].isin(ligands)].itertuples(index=False,name=None):
            genes=set(edge[0].split('_')+edge[1].split('_'));missing=genes-measured
            edge_coverage.append(dict(cohort=cohort,ligand=edge[0],receptor=edge[1],retained_in_targeted_panel=not bool(missing),not_retained_genes=';'.join(sorted(missing)),note='not retained is not biological absence; targeted panel follows frozen LIANA resources'))
        comp_map={'AT2':'epithelial','fibroblasts':'fibroblast','macrophages':'myeloid'}
        # Keep donor-level detection and full-library denominators, not group means.
        sources={label:frame.set_index(['donor','gene']) for label,frame in p[p.disease=='IPF'].groupby('label')}
        for (comp,label),ee in e.groupby(['compartment','label']):
            primary=ee[ee.omitted_donor=='__none__']
            if not primary.eligible.any():
                for row in primary.itertuples():eligibility.append(dict(cohort=cohort,compartment=comp,label=label,direction=row.direction,target_eligible=False,mapped_targets=row.mapped_targets,scope='both',eligible_ligands=0,reason=row.reason))
                continue
            receiver_labels=[k for k,v in labels.items() if v==comp_map[comp]] if label=='__broad__' else [label]
            receiver=p[(p.disease=='IPF')&p.label.isin(receiver_labels)].groupby(['donor','gene'],observed=True).agg(cells=('cells','sum'),detected_cells=('detected_cells','sum'))
            receiver['detection_fraction']=receiver.detected_cells/receiver.cells
            pair_support={}
            for ligand in ligands:
                receptor_entries=prior_lr.loc[prior_lr['from']==ligand,'to'].unique().tolist()
                receptor_entries=[r for r in receptor_entries if r not in ['IL1R2','SIGIRR']]
                # Require essential signalling partners for the canonical IL-1/TGF-beta families.
                if ligand in ['IL1A','IL1B']:receptor_entries=['IL1R1_IL1RAP'] if any(r in ['IL1R1','IL1RAP'] for r in receptor_entries) else []
                if ligand=='TGFB1':receptor_entries=['TGFBR1_TGFBR2'] if any(r in ['TGFBR1','TGFBR2'] for r in receptor_entries) else []
                for source_label,source in sources.items():
                    if ligand not in source.index.get_level_values('gene'):continue
                    possible_donors=[];matching=[]
                    for donor in sorted(set(source.index.get_level_values('donor'))&set(receiver.index.get_level_values('donor'))):
                        s=source.loc[(donor,ligand)]
                        if s.cells<50 or s.detection_fraction<.1:continue
                        good=[]
                        for receptor in receptor_entries:
                            subunits=receptor.split('_')
                            if not set(subunits)<=measured:continue
                            if all((donor,g) in receiver.index and receiver.loc[(donor,g),'cells']>=50 and receiver.loc[(donor,g),'detection_fraction']>=.1 for g in subunits):good.append(receptor)
                        if good:possible_donors.append(donor);matching.extend(good)
                    focused=(labels.get(source_label),comp_map[comp]) in DIRECTIONS
                    pair_support[(ligand,source_label)]=set(possible_donors)
                    support.append(dict(cohort=cohort,compartment=comp,label=label,ligand=ligand,source_label=source_label,focused_direction=focused,matching_receptors=';'.join(sorted(set(matching))),eligible_IPF_donors=len(possible_donors),donors=';'.join(possible_donors),minimum_required=3))
            for scope in ['focused_triad','source_agnostic_planned_panel']:
                for row in primary.itertuples():
                    subset=scores[(scores.compartment==comp)&(scores.label==label)&(scores.direction==row.direction)]
                    allowed_full=[]
                    for ligand in ligands:
                        donor_sets=[v for (l,s),v in pair_support.items() if l==ligand and (scope!='focused_triad' or (labels.get(s),comp_map[comp]) in DIRECTIONS)]
                        # No best subtype per donor: source labels remain distinct;
                        # a ligand qualifies if a fixed source label has ≥3 donors.
                        if any(len(v)>=3 for v in donor_sets):allowed_full.append(ligand)
                    eligibility.append(dict(cohort=cohort,compartment=comp,label=label,direction=row.direction,target_eligible=bool(row.eligible),mapped_targets=row.mapped_targets,scope=scope,eligible_ligands=len(allowed_full) if row.eligible else 0,reason=row.reason))
                    if not row.eligible:continue
                    for omitted,z in subset.groupby('omitted_donor'):
                        allowed=[]
                        for ligand in allowed_full:
                            donor_sets=[v for (l,s),v in pair_support.items() if l==ligand and (scope!='focused_triad' or (labels.get(s),comp_map[comp]) in DIRECTIONS)]
                            if any(len(v-{omitted})>=3 for v in donor_sets):allowed.append(ligand)
                        z=z[z.ligand.isin(allowed)&np.isfinite(z.pearson)].copy();z['cohort']=cohort;z['scope']=scope
                        z['rank']=z.pearson.rank(ascending=False,method='min');z['ranked_ligands']=len(z)
                        rankings.append(z)
    r=pd.concat(rankings,ignore_index=True)
    assert not r.duplicated(['cohort','compartment','label','direction','scope','omitted_donor','ligand']).any()
    full=r[r.omitted_donor=='__none__'].copy();full.to_csv(out/'primary_candidate_rankings.csv',index=False)
    r.to_csv(out/'candidate_rankings_with_omissions.csv',index=False)
    for keys,g in r.groupby(['cohort','compartment','label','direction','scope','ligand']):
        f=g[g.omitted_donor=='__none__'];loo=g[g.omitted_donor!='__none__']
        if f.empty:continue
        stability.append(dict(zip(['cohort','compartment','label','direction','scope','ligand'],keys))|dict(primary_pearson=float(f.pearson.iloc[0]),primary_rank=float(f['rank'].iloc[0]),planned_omissions=expected_omissions[keys[:3]],eligible_omissions=len(loo),complete_omission_coverage=len(loo)==expected_omissions[keys[:3]],min_rank=float(loo['rank'].min()) if len(loo) else None,max_rank=float(loo['rank'].max()) if len(loo) else None,minimum_pearson=float(loo.pearson.min()) if len(loo) else None,maximum_pearson=float(loo.pearson.max()) if len(loo) else None))
    pd.DataFrame(stability).to_csv(out/'candidate_ranking_stability.csv',index=False)
    pd.DataFrame(support).to_csv(out/'fixed_source_receiver_support.csv',index=False)
    pd.DataFrame(eligibility).to_csv(out/'combined_eligibility.csv',index=False)
    pd.DataFrame(edge_coverage).to_csv(out/'prior_edge_panel_coverage.csv',index=False)
    qa=dict(status='passed',checked_utc=datetime.now(timezone.utc).isoformat(),primary_candidate_rows=len(full),rows_with_omissions=len(r),source_pairs=len(support),checks=['unique ranking keys','same-donor expression support >=50 cells and >=10% detection','fixed source label >=3 IPF donors','at least ten mapped DE targets','full DE parity and centered-vector Pearson checks in R'])
    (out/'validation.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa))

if __name__=='__main__':main()
