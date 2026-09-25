"""Summarize completed second-batch outputs and render figures without refitting."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE=Path(__file__).resolve().parents[1]
S=BASE/'tables/ire1_stability_2026-09-25'
D=BASE/'tables/direct_marks_2026-09-25'
OUT=BASE/'tables/second_batch_summary'
FIG=BASE/'figures/second_batch'

def sha(p):
    with p.open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def read(p):return pd.read_csv(p,sep='\t')
def save(fig,name):
    for ext in ['png','svg']:fig.savefig(FIG/f'{name}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

def main():
    global OUT,FIG
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id',required=True,help='Fresh presentation directory; existing evidence is never overwritten')
    args=parser.parse_args()
    if not re.fullmatch(r'[a-zA-Z0-9_-]+',args.run_id):parser.error('Use letters, digits, underscores or hyphens')
    OUT=BASE/'tables'/args.run_id;FIG=BASE/'figures'/args.run_id
    record_path=BASE/'reports'/f'{args.run_id}_run.json'
    if any(p.exists() for p in [OUT,FIG,record_path]):
        raise SystemExit('Refusing to overwrite an existing presentation run')
    OUT.mkdir(parents=True);FIG.mkdir(parents=True)
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    inputs={}
    for run in ['ire1_run.json','ire1_stability_run.json','direct_mark_run.json']:
        p=BASE/'reports'/run;r=json.loads(p.read_text());assert r['status']=='completed'
        inputs[p.relative_to(BASE).as_posix()]=sha(p)
        for name,digest in r['output_sha256'].items():assert sha(BASE/name)==digest,name;inputs[name]=digest
    # Verify the processed counts and sample order against the first run's own
    # recorded MD5 contract, independently of the wrapper's newer SHA256 record.
    contract=read(BASE/'tables/ire1/input_contract.tsv').set_index('key').value
    md5_checks={}
    for key,relative in [('counts','processed/ire1/counts.tsv'),('samples','tables/ire1/sample_manifest.tsv')]:
        with (BASE/relative).open('rb') as h:digest=hashlib.file_digest(h,'md5').hexdigest()
        assert digest==contract[key+'_md5'],key
        md5_checks[key]=digest
    cfg=json.loads((BASE/'config/direct_mark_loci.json').read_text())
    inputs['config/direct_mark_loci.json']=sha(BASE/'config/direct_mark_loci.json')
    primary=read(BASE/'tables/ire1/gene_effects.tsv').set_index('feature_id')
    markers=read(BASE/'tables/ire1/marker_mapping.tsv')
    effects=read(S/'focus_effects.tsv');small=read(S/'S135_descriptive.tsv').set_index('feature_id')
    names=dict(zip(markers.feature_id,markers.symbol))
    focus=list(markers.feature_id)+[g for g in primary.index[primary.FDR<.05] if g not in names]
    rows=[]
    for gene in focus:
        p=primary.loc[gene];loo=effects[(effects.feature_id==gene)&(effects.kind=='leave_one_mouse_out')]
        within=effects[(effects.feature_id==gene)&(effects.run=='S061')].iloc[0]
        rows.append(dict(feature_id=gene,label=names.get(gene,gene),selection='frozen_marker' if gene in names else 'post_selection_primary_FDR_hit',
            primary_logFC=p.logFC,primary_FDR=p.FDR,LOMO_min_logFC=loo.logFC.min(),LOMO_max_logFC=loo.logFC.max(),
            LOMO_same_direction=int((np.sign(loo.logFC)==np.sign(p.logFC)).sum()),LOMO_n=len(loo),
            LOMO_FDR05=int((loo.FDR<.05).sum()),LOMO_min_FDR=loo.FDR.min(),LOMO_max_FDR=loo.FDR.max(),
            LOMO_max_abs_change=float(abs(loo.logFC-p.logFC).max()),S061_logFC=within.logFC,S061_FDR=within.FDR,
            S135_descriptive_log2_ratio=small.loc[gene,'log2_ratio_offset_0.5']))
    summary=pd.DataFrame(rows);summary.to_csv(OUT/'ire1_focus_stability.tsv',sep='\t',index=False)
    pathway=read(S/'pathway_sensitivity.tsv');primary_pathway=read(BASE/'tables/ire1/gene_set_camera.tsv')
    primary_pathway['run']='primary';pathway=pd.concat([primary_pathway,pathway],ignore_index=True)
    pathway.to_csv(OUT/'ire1_pathway_comparison.tsv',sep='\t',index=False)
    diag=read(S/'fit_diagnostics.tsv')
    contrasts=read(D/'histone_contrasts.tsv');hist=[]
    for key,g in contrasts.groupby(['gene','gene_group','mark','window','contrast'],sort=False):
        g=g.set_index('replicate');one=g.loc['CUT1'];two=g.loc['CUT2'];valid=bool(one.eligible and two.eligible)
        hist.append(dict(zip(['gene','gene_group','mark','window','contrast'],key))|dict(
            CUT1_log2_ratio_difference=one.log2_ratio_difference,CUT2_log2_ratio_difference=two.log2_ratio_difference,
            both_eligible=valid,direction_concordant=bool(np.sign(one.log2_ratio_difference)==np.sign(two.log2_ratio_difference)) if valid else None,
            CUT1_raw_mark_log2_ratio=one.raw_mark_log2_ratio,CUT2_raw_mark_log2_ratio=two.raw_mark_log2_ratio,
            CUT1_H3_log2_ratio=one.H3_log2_ratio,CUT2_H3_log2_ratio=two.H3_log2_ratio))
    hist=pd.DataFrame(hist);hist.to_csv(OUT/'histone_replicate_concordance.tsv',sep='\t',index=False)
    # Eight frozen markers remain the main effect panel; post-selection hits are
    # fully available in the summary, not silently added to the marker family.
    fig,axes=plt.subplots(1,3,figsize=(14,5.5),gridspec_kw={'width_ratios':[1.3,1.1,1]})
    ax=axes[0]
    for y,row in enumerate(summary[summary.selection=='frozen_marker'].itertuples()):
        vals=effects[(effects.feature_id==row.feature_id)&(effects.kind=='leave_one_mouse_out')].logFC
        ax.scatter(vals,np.full(len(vals),y),s=22,color='#377E96',alpha=.65)
        ax.scatter(row.primary_logFC,y,s=55,color='black',marker='D',zorder=4)
        ax.scatter(row.S061_logFC,y,s=42,color='#CA8C2C',marker='s',zorder=3)
    ax.axvline(0,color='.65',lw=.8);ax.set_yticks(range(len(markers)),markers.symbol)
    ax.invert_yaxis();ax.set(xlabel='KIRA8 − vehicle (log2 fold change)',title='A  Predefined epithelial markers')
    ax.legend(handles=[Line2D([],[],marker='D',color='black',ls='',label='Primary: 10 mice'),Line2D([],[],marker='o',color='#377E96',ls='',label='Each omission: 9 mice'),Line2D([],[],marker='s',color='#CA8C2C',ls='',label='S061: 3 vs 3 mice')],frameon=False,loc='upper center',bbox_to_anchor=(.5,-.16),fontsize=8)
    runs=['primary']+list(diag.run)
    for gene_set,color,label in [('TGFbeta','#8661A6','TGF-β'),('Gene Ontology 0006986','#36917E','UPR')]:
        g=pathway[pathway.gene_set==gene_set].set_index('run').loc[runs]
        axes[1].plot(-np.log10(g.FDR),range(len(runs)),marker='o',ls='',color=color,label=label)
    axes[1].axvline(-np.log10(.05),color='.5',ls='--',lw=.8)
    axes[1].set_yticks(range(len(runs)),[r.replace('omit_','Omit ') for r in runs]);axes[1].invert_yaxis()
    axes[1].set(xlabel='−log10 pathway FDR',title='B  Eligible pathways');axes[1].legend(frameon=False,fontsize=8)
    axes[2].barh(range(len(diag)),diag.genes_FDR05,color=['#377E96']*10+['#CA8C2C'])
    axes[2].axvline(4,color='black',ls='--',lw=1,label='Primary = 4')
    axes[2].set_yticks(range(len(diag)),[r.replace('omit_','Omit ') for r in diag.run]);axes[2].invert_yaxis()
    axes[2].set(xlabel='Genes with whole-family FDR < 0.05',title='C  Discovery-count sensitivity');axes[2].legend(frameon=False,fontsize=8)
    fig.suptitle('IRE1α inhibition: effect directions and sample sensitivity',fontsize=15)
    fig.text(.5,.005,'GSE190821 · robust edgeR QL and estimated-correlation CAMERA · each omission is a sensitivity fit, not a new cohort\nDots are point estimates, not confidence intervals. S135 (2 vs 2) is descriptive only in the source table.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.16,1,.95));save(fig,'a1_ire1_stability')
    # Histone CPM tracks retain every state/replicate, without merging replicates.
    tracks=read(D/'locus_tracks.tsv');colors={'iAT2':'#3388B0','iATCs':'#A54F85','iAT1':'#469C78'}
    marks=['H3K27ac','H3K4me3','H3K27me3','H3']
    fig,axes=plt.subplots(4,3,figsize=(12,10),sharex=True)
    for i,mark in enumerate(marks):
        for j,gene in enumerate(cfg['locus_panels']):
            ax=axes[i,j];block=tracks[(tracks.mark==mark)&(tracks.gene==gene)]
            for (state,rep),g in block.groupby(['state','replicate']):
                g=g.sort_values('offset_bp');ax.plot(g.offset_bp/1000,g.mean_CPM,c=colors[state],ls='-' if rep=='CUT1' else '--',lw=1,alpha=.85)
            ax.axvline(0,color='.75',lw=.7);ax.set_ylim(bottom=0)
            if i==0:ax.set_title(gene,fontsize=12)
            if j==0:ax.set_ylabel(mark+'\nCPM signal')
            if i==3:ax.set_xlabel('Position from gene-boundary TSS (kb)')
    handles=[Line2D([],[],color=c,label=s) for s,c in colors.items()]
    handles += [Line2D([],[],color='.35',ls=ls,label=r) for r,ls in [('CUT1','-'),('CUT2','--')]]
    fig.legend(handles=handles,ncol=5,loc='upper center',bbox_to_anchor=(.5,.955),frameon=False)
    fig.suptitle('Direct histone profiles across induced alveolar states',fontsize=15)
    fig.text(.5,.01,'GSE289683 / GSE291333 · T2T-CHM13v2.0 · 100-bp means · two preparations per state, one iPSC line\nAxes scale independently by locus/mark. CPM is not absolute or spike-in-scaled occupancy; H3 is shown as control.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.06,1,.915));save(fig,'a1_histone_loci')
    genes=[gene for group in cfg['gene_groups'].values() for gene in group]
    fig,axes=plt.subplots(1,4,figsize=(13,9.2),sharey=True)
    columns=[('iATCs_minus_'+s,r) for s in ['iAT2','iAT1'] for r in ['CUT1','CUT2']]
    for ax,mark in zip(axes,marks):
        if mark=='H3':
            sub=contrasts[(contrasts.mark=='H3K27ac')&(contrasts.window=='promoter')];value='H3_log2_ratio'
        else:
            sub=contrasts[(contrasts.mark==mark)&(contrasts.window=='promoter')];value='log2_ratio_difference'
        pivot=sub.pivot(index='gene',columns=['contrast','replicate'],values=value).reindex(index=genes,columns=pd.MultiIndex.from_tuples(columns))
        cmap=plt.get_cmap('RdBu_r').copy();cmap.set_bad('#DADADA')
        image=ax.imshow(pivot.to_numpy(),vmin=-3,vmax=3,cmap=cmap,aspect='auto')
        ax.set_title(mark+(' control' if mark=='H3' else ' / H3'))
        ax.set_xticks(range(4),['vs AT2\nCUT1','vs AT2\nCUT2','vs AT1\nCUT1','vs AT1\nCUT2'],rotation=40,ha='right')
        ax.set_yticks(range(len(genes)),genes)
        for boundary in [4.5,9.5,18.5]:ax.axhline(boundary,color='white',lw=1.5)
    fig.suptitle('iATCs relative to iAT2 and iAT1: promoter histone profiles',fontsize=15)
    fig.subplots_adjust(left=.09,right=.88,top=.92,bottom=.15,wspace=.15)
    cax=fig.add_axes([.91,.30,.015,.43]);fig.colorbar(image,cax=cax,label='Log2 ratio difference (display clipped at ±3)',extend='both')
    fig.text(.5,.025,'TSS ±1 kb · each CUT label shown separately · grey: the frozen signal-coverage/positive-control rule was not met\nFirst three panels: log2(mark/H3) state difference; fourth: log2(H3) state difference. No inferential tests.',ha='center',fontsize=9)
    save(fig,'a1_histone_promoters')
    domains=read(D/'methylation_domain_overlap.tsv');main=domains[domains.coordinate_assumption=='one_based_inclusive']
    fig,axes=plt.subplots(1,3,figsize=(9,8),sharey=True)
    for ax,kind in zip(axes,['UMR','LMR','PMD']):
        p=main[main.domain==kind].pivot(index='gene',columns='day',values='domain_overlap_fraction').reindex(genes)[['D0','D4','D6']]
        im=ax.imshow(p,aspect='auto',vmin=0,vmax=1,cmap='Blues');ax.set_title(kind)
        ax.set_xticks(range(3),p.columns);ax.set_yticks(range(len(genes)),genes)
    fig.suptitle('Normal differentiation: methylation-domain context',fontsize=14)
    fig.subplots_adjust(left=.13,right=.85,top=.92,bottom=.13,wspace=.12)
    cax=fig.add_axes([.89,.30,.02,.45]);fig.colorbar(im,cax=cax,label='Promoter fraction overlapping a called domain')
    fig.text(.5,.02,'GSE150527 · one donor, D0/D4/D6 · hg19 gene-boundary TSS ±1 kb\nDomain overlap is not CpG methylation percentage or a replicated DMR test.',ha='center',fontsize=9)
    save(fig,'a1_methylation_domains')
    record=dict(status='completed',utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(Path(__file__)),
        input_sha256=inputs,original_input_md5_verified=md5_checks,
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for folder in [OUT,FIG] for p in sorted(folder.iterdir())})
    record_path.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print(summary[['label','primary_logFC','LOMO_same_direction','LOMO_FDR05','S061_logFC']].to_string(index=False))
    print('Histone promoter pairs:',hist[hist.window=='promoter'].groupby('mark')[['both_eligible','direction_concordant']].sum().to_string())

if __name__=='__main__':main()
