"""Scientific figures for the executed exploratory pilot; no inferential error bars."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from exploratory_common import BASE, TABLES, FIG

COLORS = {'AT2':'#176B91','AT1':'#B64A39'}
plt.rcParams.update({'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False,
    'svg.fonttype':'none', 'savefig.facecolor':'white'})


def save(fig,name):
    for extension in ['png','svg']:
        fig.savefig(FIG/f'{name}.{extension}',dpi=170,bbox_inches='tight')
    plt.close(fig)


def discovery():
    f=pd.read_csv(TABLES/'repair_heldout_effects.csv')
    fig,axes=plt.subplots(1,2,figsize=(12,5.1),layout='constrained')
    units=sorted(f.unit.unique())
    for ax,value,title in zip(axes,['difference','standardized_difference'],
        ['A  Held-out intermediate − endpoint score','B  Same difference / pooled endpoint-cell SD']):
        for e,offset in [('AT2',-.13),('AT1',.13)]:
            d=f[f.endpoint==e].set_index('unit').loc[units]
            ax.scatter(d[value],np.arange(len(units))+offset,color=COLORS[e],label=f'vs {e}',s=36)
        ax.axvline(0,color='#555555',lw=.8)
        ax.set(yticks=np.arange(len(units)),yticklabels=units,xlabel='Rank-score difference' if value=='difference' else 'Descriptive standardized difference',title=title)
        ax.invert_yaxis();ax.grid(axis='x',alpha=.18);ax.legend(frameon=False,loc='lower right')
    fig.suptitle('Repair discovery: each point uses a program trained without that mouse',fontsize=13)
    fig.supxlabel('9 mice, days 10–15. Source labels reused; cell SD is a scale, not uncertainty across animals.',fontsize=9)
    save(fig,'repair_discovery')


def challenges():
    apparent=pd.read_csv(TABLES/'repair_apparent_effects.csv')
    depth=pd.read_csv(TABLES/'repair_depth_matched_effects.csv')
    low=pd.read_csv(TABLES/'repair_low_cycle_effects.csv')
    controls=pd.read_csv(TABLES/'repair_control_effects.csv')
    reduced=controls[controls.module=='candidate_without_generic_controls']
    mixture=pd.read_csv(TABLES/'repair_endpoint_mixture_checks.csv')
    random=pd.read_csv(TABLES/'repair_matched_random_set_effects.csv')
    fig,axes=plt.subplots(2,2,figsize=(13,8.2),layout='constrained')
    groups=[('Full source',apparent),('Depth matched',depth),('31-gene variant',reduced),('Low cycle',low)]
    ax=axes[0,0]
    for i,(label,data) in enumerate(groups):
        for endpoint,shift in [('AT2',-.13),('AT1',.13)]:
            d=data[data.endpoint==endpoint]
            q=d.all_states_at_least_30
            for good in [True,False]:
                values=d.loc[q==good,'difference']
                ax.scatter(np.full(len(values),i+shift),values,s=25,
                    facecolors=COLORS[endpoint] if good else 'none',edgecolors=COLORS[endpoint],alpha=.8)
            ax.plot([i+shift-.1,i+shift+.1],[d.difference.median()]*2,color=COLORS[endpoint],lw=2)
    ax.axhline(0,color='#555555',lw=.8)
    ax.set(xticks=range(4),xticklabels=['Full source\n9/9 qualify','Depth matched\n9/9 qualify','31-gene variant\n9/9 qualify','Low cycle\n1/9 qualifies'],
        ylabel='Intermediate − endpoint score',title='A  Source sensitivity checks')
    ax=axes[0,1]
    for i,(weight,d) in enumerate(mixture.groupby('AT2_weight')):
        ax.scatter(np.full(len(d),i),d.intermediate_minus_mixture,color='#387660',s=28)
        ax.plot([i-.15,i+.15],[d.intermediate_minus_mixture.median()]*2,color='#153E30',lw=2)
    ax.axhline(0,color='#555555',lw=.8)
    ax.set(xticks=range(3),xticklabels=['25% AT2','50% AT2','75% AT2'],ylabel='Intermediate − synthetic mixture score',
        title='B  Endpoint mixtures: 200 pairs per mouse / weight')
    ax=axes[1,0]
    for i,e in enumerate(['AT2','AT1']):
        d=random[random.endpoint==e]
        jitter=np.random.default_rng(i).uniform(-.15,.15,len(d))
        ax.scatter(i+jitter,d.median_difference,color=COLORS[e],s=7,alpha=.25)
        ax.scatter(i,apparent.loc[apparent.endpoint==e,'difference'].median(),marker='D',s=70,
            color=COLORS[e],edgecolors='black',zorder=3)
    ax.axhline(0,color='#555555',lw=.8)
    ax.set(xticks=[0,1],xticklabels=['vs AT2','vs AT1'],ylabel='Median paired score difference',
        title='C  200 matched gene sets (dots); candidate (diamonds)')
    ax=axes[1,1]
    config=json.loads((BASE/'exploratory_config.json').read_text())
    names=config['generic_control_names']
    labels={'HALLMARK_P53_PATHWAY':'p53','HALLMARK_HYPOXIA':'Hypoxia','HALLMARK_INFLAMMATORY_RESPONSE':'Inflammatory',
        'HALLMARK_TNFA_SIGNALING_VIA_NFKB':'TNF / NFκB','HALLMARK_E2F_TARGETS':'E2F','HALLMARK_G2M_CHECKPOINT':'G2M',
        'Immediate_early_AP1_panel':'AP1 / immediate early'}
    for e,offset in [('AT2',-.12),('AT1',.12)]:
        vals=[controls[(controls.module==n)&(controls.endpoint==e)].difference.median() for n in names]
        ax.scatter(vals,np.arange(len(names))+offset,color=COLORS[e],s=32)
    ax.axvline(0,color='#555555',lw=.8)
    ax.set(yticks=range(len(names)),yticklabels=[labels.get(n,n) for n in names],xlabel='Median intermediate − endpoint score',
        title='D  Generic-response panels: descriptive benchmarks')
    ax.invert_yaxis()
    legend=[Line2D([],[],marker='o',linestyle='none',color=COLORS[e],label=f'vs {e}') for e in ['AT2','AT1']]
    legend.append(Line2D([],[],marker='o',linestyle='none',color='#555555',markerfacecolor='none',label='Below 30 cells in ≥1 state'))
    fig.legend(handles=legend,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Specificity challenges: technical robustness does not establish transition-specific biology',fontsize=13)
    save(fig,'repair_specificity')


def transfer():
    config=json.loads((BASE/'exploratory_config.json').read_text())
    data={k:pd.read_csv(TABLES/f'{k}_transfer_effects.csv') for k in ['development','intestine']}
    fig,axes=plt.subplots(2,2,figsize=(14,10),layout='constrained')
    limits={0:[],1:[]}
    for index,ax in enumerate(axes.flat):
        row,col=divmod(index,2)
        dataset=['development','intestine'][col]
        module=['candidate','candidate_without_generic_controls'][row]
        f=data[dataset][data[dataset].module==module]
        states=config['transfer']['developmental_states' if dataset=='development' else 'intestinal_states']
        units=sorted(f.unit.unique(),key=lambda x:(int(x.split(':')[0][1:]),x) if dataset=='development' else (0,x))
        counts=f.drop_duplicates('unit').set_index('unit')
        for endpoint,offset,color in [(states[0],-.13,COLORS['AT2']),(states[2],.13,COLORS['AT1'])]:
            d=f[f.endpoint==endpoint].set_index('unit').loc[units]
            for i,(_,r) in enumerate(d.iterrows()):
                ax.scatter(r.difference,i+offset,s=44 if r.all_states_at_least_30 else 22,edgecolors=color,facecolors=color if r.all_states_at_least_30 else 'none',zorder=3)
            limits[index%2].extend(d.difference.tolist())
        labels=[f'{u}  ({int(counts.loc[u,"n_start"])}/{int(counts.loc[u,"n_intermediate"])}/{int(counts.loc[u,"n_end"])})' for u in units]
        ax.axvline(0,color='#555555',lw=.9)
        panel=['A','B','C','D'][index]
        context='Development: capture groups' if dataset=='development' else 'Intestine: verified mice'
        variant='50-gene primary' if row==0 else '31-gene control variant'
        ax.set(yticks=range(len(units)),yticklabels=labels,xlabel='Intermediate − endpoint rank score',
            title=f'{panel}  {context} | {variant}')
        ax.invert_yaxis();ax.grid(axis='x',alpha=.18)
    for col in range(2):
        pad=.01 if col==0 else .002
        for ax in axes[:,col]:ax.set_xlim(min(limits[col])-pad,max(limits[col])+pad)
    legend=[Line2D([],[],marker='o',linestyle='none',color=COLORS[e],label=label) for e,label in [('AT2','vs start (AT2 / Stem)'),('AT1','vs end (AT1 / Mature proximal)')]]
    legend.append(Line2D([],[],marker='o',linestyle='none',color='#555555',markerfacecolor='none',label='Below 30 cells in ≥1 state'))
    fig.legend(handles=legend,loc='outside lower center',ncol=3,frameon=False)
    fig.suptitle('Frozen repair program and pre-transfer variant: no gene or branch retuning\nLabels show start / intermediate / end cell counts. Development groups are not animals. Column scales differ.',fontsize=12)
    save(fig,'frozen_transfer')


if __name__=='__main__':
    discovery();challenges();transfer()
    print('Saved three exploratory figures as PNG and editable SVG.')
