"""Report eligible NicheNet-prior rankings without an activation claim."""
from pathlib import Path
import json
PAPER=Path(__file__).resolve().parents[1];ROOT=PAPER.parents[1]

def main():
    import numpy as np,pandas as pd,matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    out=PAPER/'trials/u5_ligand_targets';assert json.loads((out/'validation.json').read_text())['status']=='passed'
    r=pd.read_csv(out/'primary_candidate_rankings.csv');stable=pd.read_csv(out/'candidate_ranking_stability.csv');eligible=pd.read_csv(out/'combined_eligibility.csv')
    palette=json.loads((ROOT/'analysis/config/palette.json').read_text())
    planned=[('GSE136831','AT2','ATII'),('GSE136831','fibroblasts','__broad__'),('GSE136831','macrophages','__broad__'),('GSE135893','AT2','AT2'),('GSE135893','fibroblasts','__broad__'),('GSE135893','macrophages','__broad__')]
    ligands=['IL1A','IL1B','AREG','HBEGF','TGFB1'];values=[];rows=[];labels=[]
    for cohort,comp,label in planned:
        for direction in ['up','down']:
            labels.append(f'{cohort} | {comp} | {direction}');line=[]
            for ligand in ligands:
                z=r[(r.cohort==cohort)&(r.compartment==comp)&(r.label==label)&(r.direction==direction)&(r.scope=='focused_triad')&(r.ligand==ligand)]
                e=eligible[(eligible.cohort==cohort)&(eligible.compartment==comp)&(eligible.label==label)&(eligible.direction==direction)]
                if len(z):
                    q=z.iloc[0];value=float(q.pearson);status='eligible';rank=float(q['rank']);n=int(q.ranked_ligands)
                else:
                    value=np.nan;rank=np.nan;n=0;status='target_ineligible' if e.empty or not e.target_eligible.any() else 'source_receiver_expression_ineligible'
                line.append(value);rows.append(dict(cohort=cohort,compartment=comp,label=label,direction=direction,ligand=ligand,pearson=value,status=status,rank=rank,ranked_ligands=n))
            values.append(line)
    arr=np.asarray(values);fig,ax=plt.subplots(figsize=(10,8));fig.set_facecolor(palette['surface']);ax.set_facecolor(palette['deemph'])
    cmap=LinearSegmentedColormap.from_list('signed_prior',[palette['categorical']['2'],palette['surface'],palette['categorical']['1']]);cmap.set_bad(palette['deemph'])
    maximum=max(.1,float(np.nanmax(abs(arr))));im=ax.imshow(arr,cmap=cmap,vmin=-maximum,vmax=maximum,aspect='auto')
    for i in range(len(labels)):
        for j in range(len(ligands)):
            q=rows[i*len(ligands)+j]
            text=f"{arr[i,j]:.3f}\nrank {int(q['rank'])}/{q['ranked_ligands']}" if np.isfinite(arr[i,j]) else ('Targets <10' if q['status']=='target_ineligible' else 'Expression\ngate')
            ax.text(j,i,text,ha='center',va='center',fontsize=8,color=palette['ink'])
    ax.set_xticks(range(len(ligands)),ligands);ax.set_yticks(range(len(labels)),labels);ax.tick_params(length=0,labelsize=9,colors=palette['ink'])
    for s in ax.spines.values():s.set_visible(False)
    bar=fig.colorbar(im,ax=ax,fraction=.035,pad=.025);bar.set_label('Pearson fit to the specified target set',fontsize=9)
    fig.suptitle('Ligand–target compatibility differs for up- and downregulated genes',x=.02,ha='left',fontsize=14,color=palette['ink'])
    fig.text(.02,.025,'NicheNet v2 human prior; focused fixed source–receiver candidates. Rank is among all eligible planned ligands, including exploratory families.\nA fit to downregulated targets does not mean that a ligand is inhibited or causes repression: this prior is unsigned.\nExpression gates require ≥3 IPF donors with the same source label, ≥50 cells per compartment and ≥10% subunit detection.\nTargets: donor-level DE q < 0.05, ≥10 prior-mapped genes. Observational prioritization; no independent or causal validation.',fontsize=8,color=palette['ink_2'])
    fig.subplots_adjust(left=.34,right=.92,top=.91,bottom=.18)
    for ext in ['png','svg']:fig.savefig(PAPER/'figures'/('ipf_ligand_target_eligibility_and_fit.'+ext),dpi=220,facecolor=fig.get_facecolor())
    pd.DataFrame(rows).to_csv(out/'figure_core_target_values.csv',index=False)
    relevant=stable[(stable.scope=='focused_triad')&stable.ligand.isin(['IL1A','IL1B'])]
    lines=['# IPF ligand–target prioritization','', 'Completed for eligible receivers using the checksum-verified NicheNet v2 human prior (Zenodo record 7074291). The implemented statistic is the official Pearson correlation between prior target potential and a binary DE-target vector; this is a base-R implementation of that statistic, not a claim that the complete nichenetr package or its AUPR metrics ran. Full-data DE fits reproduce the preceding pathway run, and the Pearson calculation was independently checked against centered vectors.','', '## Findings and limits','', 'GSE136831 AT2 and GSE135893 broad fibroblasts fail the ten-mapped-target gate. No ranking is inferred for these receivers. IL1A/IL1B qualify for the GSE136831 fibroblast receiver; neither is the leading fit to its upregulated targets. Their better ranks for downregulated fibroblast genes do not establish inhibition or repression: the prior is unsigned. The two cohorts therefore do not provide a replicated IL-1-specific recipient programme.','', 'Both focused triad and source-agnostic rankings are restricted to the planned ligand families with available source-panel measurements. They are not a genome-wide screen of every possible ligand. The full prior scores are retained separately and are not promoted into expression-supported findings. A prior receptor missing from the retained panel is marked unevaluable, not biologically absent. IL1R2/SIGIRR are excluded from activating candidate gates; canonical IL-1 and TGF-beta candidates require their essential receptor partners.','', '## Donor-omission stability','', 'Every omission reruns receiver filtering, TMM, voom, DE/BH and target selection. Source/receiver expression eligibility is checked again. Omission fits with fewer than three donors per arm or fewer than ten mapped targets remain missing. The stability table reports planned versus eligible omissions; ranks can depend on the number of eligible candidates and are not confidence intervals. Broad and subtype views overlap and are not independent evidence.','', '| Cohort | Receiver | Direction | Ligand | Primary rank | Eligible/planned omissions | Rank range |','|---|---|---|---|---:|---:|---|']
    for row in relevant.itertuples():lines.append(f'| {row.cohort} | {row.label} | {row.direction} | {row.ligand} | {row.primary_rank:g} | {row.eligible_omissions}/{row.planned_omissions} | {row.min_rank:g}–{row.max_rank:g} |')
    lines += ['', '![Ligand-target eligibility and fit](../../figures/ipf_ligand_target_eligibility_and_fit.png)','', 'Tables: [primary candidate ranks](primary_candidate_rankings.csv), [target/expression eligibility](combined_eligibility.csv), [fixed source/receiver support](fixed_source_receiver_support.csv), [ranking stability](candidate_ranking_stability.csv), [prior receptor-panel coverage](prior_edge_panel_coverage.csv), [plotted values](figure_core_target_values.csv), [checks](validation.json). Per-cohort records preserve target counts, background coverage, individual omission status and full-data DE parity.']
    (out/'REPORT.md').write_text('\n'.join(lines)+'\n');print('Ligand-target report and figure generated; primary rows',len(r))

if __name__=='__main__':main()
