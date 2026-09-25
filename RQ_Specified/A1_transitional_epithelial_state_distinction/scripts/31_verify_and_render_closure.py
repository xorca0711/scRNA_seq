"""Independent table checks and compact views of the completed closure analyses."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import gzip
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a1_robustness import sha

BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[1]


def bh(p):
    p=np.asarray(p,float);order=np.argsort(p,kind='stable');values=np.minimum.accumulate((p[order]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]
    result=np.empty(len(p));result[order]=np.minimum(values,1);return result


def main():
    out=BASE/'tables/closure_verification';figdir=BASE/'figures/evidence_closure';rp=BASE/'reports/closure_verification_run.json'
    if out.exists() or figdir.exists() or rp.exists():raise SystemExit('Refusing to overwrite closure verification')
    checked=0
    for name in ['closure_identity_run.json','cd44_closure_run.json']:
        record=json.loads((BASE/'reports'/name).read_text())
        for p,h in record['output_sha256'].items():assert sha(BASE/p)==h,p;checked+=1
    baseline=json.loads((BASE/'tables/evidence_closure/prior_numerical_sha256.json').read_text())
    # The gallery README was intentionally updated; every prior numerical artifact stays fixed.
    documentation_path=(BASE/'figures/README.md').relative_to(ROOT).as_posix()
    baseline.pop(documentation_path)
    for p,h in baseline.items():assert sha(ROOT/p)==h,p
    table=BASE/'tables/cd44_closure';m=pd.read_csv(table/'sample_manifest.tsv',sep='\t')
    source=pd.read_csv(BASE/'cache/inputs/GSE273123/GSE273123_Count_Matrix.txt.gz',sep='\t',index_col=0)
    used=pd.read_csv(table/'counts_input.tsv.gz',sep='\t',index_col=0)
    pd.testing.assert_frame_equal(source,used,check_names=False)
    assert set(m.mouse)=={'WT1','WT2','WT3','WT4','Mutant1','Mutant2','Mutant3','Mutant4'}
    assert list(source.columns)==list(m.sample_id)
    design=pd.read_csv(table/'design.tsv',sep='\t',index_col=0)
    assert np.linalg.matrix_rank(design.to_numpy())==10 and len(design)-design.shape[1]==6
    # Independent synthetic cell means show that the direct contrast is a difference of paired effects.
    # Use a common baseline within each mouse.
    synthetic=np.array([sorted(m.mouse.unique()).index(mouse)+({'WT':2,'Mutant':5}[g] if s=='positive' else 0)
                        for mouse,g,s in zip(m.mouse,m.genotype,m['sort'])])
    coef=np.linalg.lstsq(design.to_numpy(),synthetic,rcond=None)[0]
    ix={n:i for i,n in enumerate(design.columns)}
    assert abs(coef[ix['positiveWT']]-2)<1e-10 and abs(coef[ix['positiveMutant']]-5)<1e-10
    effects=pd.read_csv(table/'all_gene_effects.tsv',sep='\t',dtype={'feature_id':str})
    for _,d in effects.groupby('contrast'):assert np.allclose(bh(d.PValue),d.FDR,rtol=1e-10,atol=1e-15)
    assert np.allclose(bh(effects.PValue),effects.FDR_all_three,rtol=1e-10,atol=1e-15)
    wide=effects.pivot(index='feature_id',columns='contrast',values='logFC')
    assert np.allclose(wide.Mutant_positive_minus_negative-wide.WT_positive_minus_negative,wide.Mutant_minus_WT_CD44_interaction,rtol=1e-10,atol=1e-12)
    retained=pd.read_csv(table/'feature_filter.tsv',sep='\t',dtype={'feature_id':str})
    assert set(retained.loc[retained.retained,'feature_id'])==set(wide.index)
    focus=pd.read_csv(table/'focus_effects.tsv',sep='\t');loo=pd.read_csv(table/'focus_leave_one_pair_out.tsv',sep='\t')
    assert len(loo)==8*8*3 and set(loo.omitted_mouse)==set(m.mouse) and (loo.residual_df==5).all()
    rows=[]
    for (symbol,contrast),d in loo.groupby(['symbol','contrast']):
        p=focus[(focus.symbol==symbol)&(focus.contrast==contrast)].iloc[0]
        rows.append(dict(symbol=symbol,contrast=contrast,primary_logFC=p.logFC,joint_FDR=p.FDR_all_three,
                         omissions=len(d),same_direction=int((np.sign(d.logFC)==np.sign(p.logFC)).sum()),
                         min_logFC=d.logFC.min(),max_logFC=d.logFC.max()))
    # Independently count annotation properties directly from CSV; do not reuse audit helper.
    with gzip.open(BASE/'cache/followup_sources/hpcs_combined_obs.csv.gz','rt') as h:obs=list(csv.DictReader(h))
    cfg=json.loads((BASE/'config/closure_analysis_contract.json').read_text())
    tr=[r for r in obs if r['batch'] in ['traced','Hopx-MACD_traced']]
    assert len(tr)==5333 and all(cfg['hpcs']['mapping'][r['newleiden']]==r['cell type'] for r in obs)
    assert sum(r['clusterK12_stringent']=='other' for r in tr)==1282
    assert all(r['clusterK12_stringent'] in [r['clusterK12'],'other'] for r in obs)
    out.mkdir();figdir.mkdir()
    summary=pd.DataFrame(rows);summary.to_csv(out/'cd44_focus_stability.tsv',sep='\t',index=False)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    symbols=cfg['cd44']['focus_symbols'];contrasts=cfg['cd44']['contrasts']
    fig,axes=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[1.15,1]})
    vals=focus.pivot(index='symbol',columns='contrast',values='logFC').loc[symbols,contrasts]
    im=axes[0].imshow(vals,cmap='RdBu_r',vmin=-3,vmax=3,aspect='auto')
    q=focus.pivot(index='symbol',columns='contrast',values='FDR_all_three').loc[symbols,contrasts]
    for i in range(8):
        for j in range(3):axes[0].text(j,i,f'{vals.iloc[i,j]:+.2f}'+(' *' if q.iloc[i,j]<.05 else ''),ha='center',va='center',color='white' if abs(vals.iloc[i,j])>1.5 else '#202020')
    axes[0].set(yticks=range(8),yticklabels=symbols,xticks=range(3),xticklabels=['WT\nCD44+ − CD44−','Mutant\nCD44+ − CD44−','Interaction\nMutant − WT'],title='A  Shared markers are not disease-specific')
    fig.colorbar(im,ax=axes[0],shrink=.75,label='Model log2 fold change',pad=.03)
    paired=pd.read_csv(table/'marker_paired_differences.tsv',sep='\t')
    for k,(group,color,shift) in enumerate([('WT','#356B9B',-.15),('Mutant','#B25E32',.15)]):
        for i,symbol in enumerate(symbols):
            d=paired[(paired.genotype==group)&(paired.symbol==symbol)]
            axes[1].scatter(d.difference,np.repeat(i+shift,len(d)),color=color,s=20,alpha=.8,label=group if i==0 else None)
    axes[1].axvline(0,c='#555555',lw=.7)
    axes[1].set(yticks=range(8),yticklabels=symbols,ylim=(7.7,-.7),xlabel='Within-mouse CD44+ − CD44− log2CPM',title='B  Each dot is one paired mouse')
    axes[1].legend(frameon=False);axes[1].grid(axis='x',alpha=.15)
    fig.suptitle('CD44-sorted AT2 RNA: four paired mice per genotype',fontsize=15)
    fig.text(.5,.02,'* BH FDR < 0.05 over all 46,548 gene–contrast tests. Panel B is descriptive TMM log2CPM (prior count 2).\nSame published cohort; association with a sorting gate does not measure regulatory mechanism or lineage fate.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.1,1,.95))
    for ext in ['png','svg']:fig.savefig(figdir/f'a1_cd44_context.{ext}',dpi=200)
    plt.close(fig)
    abst=pd.read_csv(BASE/'tables/evidence_closure/hpcs_abstention_summary.tsv',sep='\t')
    d=abst[abst.level=='state'].sort_values('abstention_fraction')
    fig,ax=plt.subplots(figsize=(9,5.7));ax.barh(d.label,d.abstention_fraction*100,color='#547D91')
    for i,r in enumerate(d.itertuples()):ax.text(100*r.abstention_fraction+.8,i,f'{r.stringent_other}/{r.cells}',va='center',fontsize=9)
    ax.set(xlim=(0,80),xlabel='Cells assigned “other” by the stringent K12 classifier (%)',title='HPCS labels: confidence abstention depends on the RNA state')
    ax.invert_yaxis();fig.text(.5,.025,'All 5,333 traced cells retained in denominators. 1,282 abstentions; zero switches between retained K12 codes.\nState names come from newleiden; this is not independent annotation validation or cell-level inference.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.1,1,1))
    for ext in ['png','svg']:fig.savefig(figdir/f'a1_hpcs_abstention.{ext}',dpi=200)
    plt.close(fig)
    record=dict(status='completed',utc=datetime.now(timezone.utc).isoformat(),prior_numerical_files_unchanged=len(baseline),
                output_hashes_verified=checked,BH_tests_verified=len(effects),interaction_identities_verified=len(wide),
                HPCS_rows_verified=len(obs),source_sha256=sha(Path(__file__)),
                output_sha256={p.relative_to(BASE).as_posix():sha(p) for folder in [out,figdir] for p in folder.iterdir()})
    rp.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':main()
