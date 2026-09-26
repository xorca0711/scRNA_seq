"""Render the A10 follow-up from tracked summaries; no fits or intervals."""
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / 'tables/followup_v1'
FIG = HERE / 'figures'


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    FIG.mkdir(exist_ok=True)
    names = [FIG/'a10_followup.png',FIG/'a10_followup.svg']
    assert not any(p.exists() for p in names), 'Refusing to overwrite figure'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig, ax = plt.subplots(2,2,figsize=(12.4,8.4),layout='constrained')
    overlap = pd.read_csv(OUT/'diagnostic_target_overlap.tsv',sep='\t')
    labels = ['plate1','plate2','plate3','plate4']
    matrix = np.zeros((4,4))
    for _, row in overlap[overlap.scope=='drop_top_and_TDTOMATO'].iterrows():
        a,b=labels.index(row.plate_a),labels.index(row.plate_b)
        matrix[a,b]=matrix[b,a]=row.common_targets
    masked = np.ma.array(matrix,mask=np.eye(4,dtype=bool))
    ax[0,0].imshow(masked,cmap='Blues',vmin=0,vmax=2)
    for i in range(4):
        for j in range(4):
            ax[0,0].text(j,i,'—' if i==j else str(int(matrix[i,j])),ha='center',va='center',fontsize=14)
    ax[0,0].set(xticks=range(4),yticks=range(4),xticklabels=['P1','P2','P3','P4'],yticklabels=['P1','P2','P3','P4'])
    ax[0,0].set_title('A  Shared targets after TIGIT / TDTOMATO removal',loc='left',fontsize=11)
    ax[0,0].set_xlabel('Plate 3 has no overlap; plate holdout also changes target mix')

    summary = pd.read_csv(OUT/'model_comparisons.tsv',sep='\t')
    choices = ['proliferation_over_baseline','remaining_over_proliferation','full_growth_over_baseline']
    display = ['Proliferation\nover baseline','Remaining growth\nover proliferation','Full growth\nover baseline']
    x=np.arange(3)
    for offset,(evaluation,label,color) in zip([-.18,.18],[('within_group','Within observed groups','#3B769D'),('plate_shift','Whole-plate shift','#B86A42')]):
        frame=summary[(summary.setting=='primary')&(summary.scale=='inherited_log2')&(summary.evaluation==evaluation)].set_index('comparison')
        values=frame.loc[choices,'relative_error_reduction'].to_numpy()*100
        ax[0,1].bar(x+offset,values,.34,label=label,color=color)
    ax[0,1].axhline(2,color='gray',linestyle=':',linewidth=1)
    ax[0,1].set(xticks=x,xticklabels=display,ylabel='Reduction of reference squared error (%)')
    ax[0,1].legend(frameon=False,fontsize=9)
    ax[0,1].set_title('B  The six-score block adds beyond E2F / G2M',loc='left',fontsize=11)

    absolute=pd.read_csv(OUT/'model_absolute_performance.tsv',sep='\t')
    x=np.arange(4)
    for offset,(model,label,color) in zip([-.18,.18],[('baseline','Baseline','#B9BEC3'),('baseline+proliferation+remaining_growth','Full growth','#3B769D')]):
        frame=absolute[(absolute.setting=='primary')&(absolute.scale=='inherited_log2')&(absolute.model==model)].set_index('plate')
        values=frame.loc[labels,'r2_against_heldout_mean'].to_numpy()
        ax[1,0].bar(x+offset,values,.34,label=label,color=color)
    ax[1,0].axhline(0,color='black',linewidth=.8)
    ax[1,0].set(xticks=x,xticklabels=['P1','P2','P3','P4'],ylabel='Held-out R² against that plate’s observed mean')
    ax[1,0].legend(frameon=False,fontsize=9,loc='lower left')
    ax[1,0].set_title('C  Absolute performance remains poor on 3 / 4 plates',loc='left',fontsize=11)

    fc=pd.read_csv(OUT/'model_fold_comparisons.tsv',sep='\t')
    frame=fc[(fc.setting=='primary')&(fc.scale=='inherited_log2')&(fc.evaluation=='within_group')&(fc.comparison=='full_growth_over_baseline')].sort_values('fold')
    colors={'plate1':'#3B769D','plate2':'#B86A42','plate3':'#397F6D','plate4':'#9A659A'}
    positions=np.arange(len(frame))
    for position,(_,row) in zip(positions,frame.iterrows()):
        color=colors[row.fold.split('-rep')[0]]
        ax[1,1].plot([position,position],[0,100*row.relative_error_reduction],color=color,linewidth=1.5)
        ax[1,1].scatter(position,100*row.relative_error_reduction,color=color,s=30)
    ax[1,1].axhline(0,color='black',linewidth=.8)
    ax[1,1].set(xticks=positions,xticklabels=[v.replace('plate','P').replace('-rep','.') for v in frame.fold],ylabel='Full-growth error reduction within each group (%)')
    ax[1,1].tick_params(axis='x',rotation=60,labelsize=8)
    ax[1,1].set_title('D  Within-group gains remain heterogeneous',loc='left',fontsize=11)
    fig.suptitle('A10 follow-up: relative gains survive; independent repair prediction remains unestablished',fontsize=14)
    for path in names:
        fig.savefig(path,dpi=170,metadata={'Creator':'scRNA_seq A10 follow-up'})
    plt.close(fig)
    record={'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'inputs':{n:hashlib.sha256((OUT/n).read_bytes()).hexdigest() for n in ['diagnostic_target_overlap.tsv','model_comparisons.tsv','model_absolute_performance.tsv','model_fold_comparisons.tsv']},
            'outputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in names},
            'scope':'Recorded estimates only, no new fit or inferential intervals. Primary inherited-log2 scale except panel A target-removal diagnostic.'}
    (OUT/'figure_run.json').write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('Rendered a10_followup.png and a10_followup.svg')


if __name__ == '__main__':
    main()
