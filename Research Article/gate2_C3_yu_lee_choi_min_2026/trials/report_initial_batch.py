"""Validate completed pathway tables and publish the initial batch report."""
from pathlib import Path
import json
import math
import sys
from datetime import datetime,timezone
PAPER=Path(__file__).resolve().parents[1]
ROOT=PAPER.parents[1]


def bh(values):
    order=sorted(range(len(values)),key=values.__getitem__)
    result=[0.0]*len(values);running=1.0
    for rank in range(len(values),0,-1):
        i=order[rank-1];running=min(running,values[i]*len(values)/rank);result[i]=running
    return result


def main():
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out=PAPER/'trials/u5_ipf_pathways'
    figures=PAPER/'figures';figures.mkdir(exist_ok=True)
    records=[];frames=[];checks=0
    for cohort in ['GSE136831','GSE135893']:
        frame=pd.read_csv(out/(cohort+'_camera.csv'))
        assert not frame.duplicated(['compartment','label','set','cell_floor','correlation_setting']).any();checks+=1
        assert frame.PValue.between(0,1).all() and frame.FDR_global.between(0,1).all();checks+=1
        assert (frame.n_IPF>=3).all() and (frame.n_control>=3).all();checks+=1
        assert (frame.NGenes>=10).all() and (frame.assayed_fraction>=.7-1e-12).all();checks+=1
        for _,part in frame.groupby(['cell_floor','correlation_setting']):
            expected=bh(part.PValue.tolist())
            assert np.allclose(expected,part.FDR_global,rtol=1e-10,atol=1e-12);checks+=1
        primary=frame[(frame.cell_floor==50)&(frame.correlation_setting=='estimated')]
        sensitivity=frame[(frame.cell_floor==50)&(frame.correlation_setting=='fixed001')]
        records.append({'cohort':cohort,'primary_tests':len(primary),'primary_q_lt_05':int((primary.FDR_global<.05).sum()),
                        'fixed001_q_lt_05':int((sensitivity.FDR_global<.05).sum())})
        frames.append(frame)
    combined=pd.concat(frames,ignore_index=True)
    qa={'checked_utc':datetime.now(timezone.utc).isoformat(),'checks_passed':checks,'cohorts':records,
        'scope':'output integrity, eligibility and global-BH recalculation; not independent scientific validation'}
    (out/'validation.json').write_text(json.dumps(qa,indent=2)+'\n')
    palette=json.loads((ROOT/'analysis/config/palette.json').read_text())
    broad=combined[(combined.cell_floor==50)&(combined.label=='__broad__')].copy()
    keys=list(dict.fromkeys(zip(broad.compartment,broad['set'])))
    fig,ax=plt.subplots(figsize=(10.8,10.3),facecolor=palette['surface'])
    ax.set_facecolor(palette['surface'])
    xlabels=[]
    for col,(cohort,setting) in enumerate([(c,s) for c in ['GSE136831','GSE135893'] for s in ['estimated','fixed001']]):
        xlabels.append(cohort+'\n'+('Primary: estimated r' if setting=='estimated' else 'Sensitivity: r = 0.01'))
        for row,(comp,name) in enumerate(keys):
            result=broad[(broad.cohort==cohort)&(broad.correlation_setting==setting)&(broad.compartment==comp)&(broad['set']==name)].iloc[0]
            q=float(result.FDR_global);direction=result.Direction
            color=palette['categorical']['1' if direction=='Up' else '2']
            size=35+65*min(-math.log10(max(q,1e-12)),4)
            ax.scatter(col,row,s=size,facecolors=color if q<.05 else palette['surface'],edgecolors=color,linewidths=1.2)
            ax.text(col+.14,row,f'{q:.2g}',fontsize=7.8,va='center',color=palette['ink'])
    rowlabels=[comp.replace('macrophages','Macrophage').replace('fibroblasts','Fibroblast')+' | '+name.replace('HALLMARK_','').replace('REACTOME_','').replace('_',' ').lower() for comp,name in keys]
    ax.set_yticks(range(len(keys)),rowlabels,fontsize=8)
    ax.set_xticks(range(4),xlabels,fontsize=8)
    ax.invert_yaxis();ax.set_xlim(-.45,3.6)
    ax.grid(axis='y',color=palette['grid'],linewidth=.5);ax.set_axisbelow(True)
    for spine in ax.spines.values():spine.set_visible(False)
    ax.tick_params(length=0,colors=palette['ink'])
    fig.suptitle('IPF pathway results depend on the correlation assumption',fontsize=15,x=.04,ha='left',color=palette['ink'])
    fig.text(.04,.938,'Donor pseudobulks; TMM / voom / CAMERA. Broad compartments shown; q corrected across all tested subtypes.',fontsize=8.7,color=palette['ink_2'])
    fig.text(.04,.048,'Numbers = global BH q. Filled = q < 0.05; open = not below threshold. Blue = Up; orange = Down.\nDirection is competitive gene-set enrichment, not measured pathway activation. Primary: no q < 0.05 in either cohort.\nPreviously inspected observational cohorts; fixed-correlation results do not replace the declared primary analysis.',fontsize=8.2,color=palette['ink_2'])
    fig.subplots_adjust(left=.45,right=.96,top=.88,bottom=.14)
    fig.savefig(figures/'ipf_pathway_primary_and_sensitivity.png',dpi=220,facecolor=fig.get_facecolor())
    plt.close(fig)
    qc=pd.read_csv(PAPER/'trials/u3_acquire_qc/sample_qc.csv')
    report=f'''# Initial analysis batch

Generated {datetime.now(timezone.utc).isoformat()}. The owner authorized staged
execution after final review. This report covers completed modules; it is not
the final multi-cohort study report.

## Completed mouse input and annotation diagnostics

Seven early IgG/anti-IL-1beta libraries were acquired and processed. Of
{int(qc.input_cells.sum()):,} deposited cells, {int(qc.qc_pass_cells.sum()):,}
pass the frozen gene-count/mitochondrial filters. No duplicate matrix
coordinates were found. The acquisition/QC stage took 107.58 seconds; the
treatment-blind clustering stage took 243.27 seconds and yielded 33 review
clusters. Source data and preparation inputs were preserved.

The marker review supports distinct stromal, immune and epithelial populations,
but several myeloid clusters need more specific labels and the alveolar group
does not yet establish a defensible KAC class. The provisional marker gates
are not used as a substitute for validated KAC or macrophage subtype identity.
No mouse treatment-effect conclusion has been released.

## Completed IPF pathway analysis

| Cohort | Primary tests at 50-cell floor | Primary global q < 0.05 | Fixed-correlation sensitivity q < 0.05 |
|---|---:|---:|---:|
'''
    for r in records:report+=f"| {r['cohort']} | {r['primary_tests']} | {r['primary_q_lt_05']} | {r['fixed001_q_lt_05']} |\n"
    report+='''
The declared primary analysis estimates residual inter-gene correlation.
Neither cohort has a pathway passing its global 0.05 FDR threshold. Fixing
correlation at 0.01 changes the results substantially. These sensitivity
findings do not replace the primary analysis or establish a robust niche
mechanism. Failure to pass FDR does not establish biological absence.

Fibroblast inflammatory-set directions differ between these cohorts in the
broad-compartment comparisons. Their mixture and small control samples need
to remain visible. At the primary floor, broad fibroblast IPF/control n is
18/4 in GSE136831 and 7/3 in GSE135893. Most fine fibroblast subtypes lack
three controls, so they are explicitly ineligible. Broad and subtype views
are not independent biological replications; AT2 broad/sole-subtype rows
are duplicate views and were retained in the declared correction family.
The floor-30 sensitivity is conditional on the cache's original broad-
compartment >=50-cell donor selection, not recovery of all 30-49-cell donors.

![IPF pathway primary and sensitivity results](figures/ipf_pathway_primary_and_sensitivity.png)

The table-integrity check independently recalculated global BH families and
verified sample/assay eligibility. These checks do not validate a biological
claim. An initial output merge failed because CAMERA returns different
columns for the two correlation settings; its schema was corrected and both
cohorts completed with the same frozen statistical choices. The failed run
record is preserved.

## Initial descriptive ligand-receptor run and remaining work

The GSE136831 ligand-receptor module uses 25 eligible donors from full raw
counts with deposited labels. It includes macrophages omitted from the old
epithelial-stromal cache. The first descriptive LIANA pass caps each donor/
subtype at 500 seeded cells, preserves fixed source-target pairs, and compares
consensus with CellChatDB. It does not provide biological-replicate P values,
an all-cell sensitivity result or a causal communication claim. Consult its
run record and output validation for completion status. Custom RNA-compatibility contrasts, further
sensitivities and ligand-target interpretation are not yet completed.

KAC identity, later-treatment strategy mapping, human lesion crosswalks,
spatial access and required context-specific coverage remain release gates.
Zenodo's linked annotated files are restricted. These are scientific/access
conditions, not an outstanding request for user approval.

Sources and execution evidence: [final review](FINAL_REVIEW.md),
[mouse QC](trials/u3_acquire_qc/run_record.json),
[cluster diagnostics](trials/u3_cluster_review/run_record.json),
[pathway specification](trials/u5_ipf_spec/specification.json),
[pathway run](trials/u5_ipf_pathways/run_record.json),
[table checks](trials/u5_ipf_pathways/validation.json), and
[ligand-receptor run](trials/u5_ipf_liana/run_record.json).
'''
    (PAPER/'INITIAL_RUN_REPORT.md').write_text(report,encoding='utf-8')
    print(json.dumps(qa),flush=True)


if __name__=='__main__':main()
