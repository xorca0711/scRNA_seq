"""All annotated native transcript starts, retaining the original TSS baseline."""
from pathlib import Path
from datetime import datetime, timezone
import importlib.util
import json
import math
import sys
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a1_robustness import sha, transcript_tss, sign

BASE = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('direct_marks', Path(__file__).with_name('15_quantify_direct_marks.py'))
direct = importlib.util.module_from_spec(spec)
spec.loader.exec_module(direct)


def main():
    out = BASE/'tables/robustness_2026-09-25/histone'
    figdir = BASE/'figures/robustness_2026-09-25'
    rp = BASE/'reports/histone_tss_robustness_run.json'
    if out.exists() or rp.exists():
        raise SystemExit('Refusing to overwrite histone TSS robustness')
    sys.path.insert(0, str(BASE/'cache/pybigtools'))
    import pybigtools
    cp=BASE/'config/robustness_batch.json'; cfg=json.loads(cp.read_text())['histone']
    inventory_path=BASE/'reports/direct_mark_input_inventory.json'
    inventory=json.loads(inventory_path.read_text())
    reference=next(r for r in inventory['files'] if r.get('accession')=='reference' and r.get('assembly')=='T2T-CHM13v2.0')
    annotation_path=BASE/reference['path']; assert sha(annotation_path)==reference['sha256']
    loci_path=BASE/'tables/direct_marks_2026-09-25/loci.tsv'
    loci=pd.read_csv(loci_path,sep='\t')
    loci=loci.query("assembly == 'T2T-CHM13v2.0' and status == 'eligible'").to_dict('records')
    assert len(loci)==23
    positions,audit=transcript_tss(annotation_path,loci)
    assert all(r['parent_gene_count']==1 for r in positions)
    assert {r['gene'] for r in audit if r['eligible']} == {r['gene'] for r in loci}
    sample_path=BASE/'tables/direct_marks_2026-09-25/sample_manifest.tsv'
    samples=pd.read_csv(sample_path,sep='\t')
    assert len(samples)==24
    original=json.loads((BASE/'reports/direct_mark_run.json').read_text())
    for name,digest in original['output_sha256'].items(): assert sha(BASE/name)==digest,name
    inputs={p.relative_to(BASE).as_posix():sha(p) for p in [cp,BASE/cfg['reuse_contract'],annotation_path,loci_path,sample_path,inventory_path,BASE/'reports/direct_mark_run.json']}
    signals=[]; chroms=None
    for source in samples.itertuples():
        path=BASE/source.source_path; assert sha(path)==source.sha256,source.source_path
        inputs[source.source_path]=source.sha256
        bw=pybigtools.open(str(path))
        if chroms is None: chroms=bw.chroms()
        assert bw.chroms()==chroms
        for r in positions:
            for width in cfg['half_widths_bp']:
                start=max(0,r['tss']-width);end=min(chroms[r['chrom']],r['tss']+width)
                value,covered=direct.interval_signal(bw.records(r['chrom'],start,end),start,end)
                signals.append(dict(**r,half_width_bp=width,start=start,end=end,state=source.state,
                    replicate=source.replicate,mark=source.mark,gsm=source.gsm,mean_CPM=value,covered_fraction=covered))
        bw.close()
    key=lambda r:(r['gene'],r['tss'],r['half_width_bp'],r['state'],r['replicate'])
    controls={key(r):r for r in signals if r['mark']=='H3'}
    for r in signals:
        h=controls[key(r)]
        eligible=r['mean_CPM']>0 and h['mean_CPM']>0 and r['covered_fraction']>=.8 and h['covered_fraction']>=.8
        r.update(H3_mean_CPM=h['mean_CPM'],H3_covered_fraction=h['covered_fraction'],ratio_eligible=eligible,
                 log2_mark_over_H3=math.log2(r['mean_CPM']/h['mean_CPM']) if eligible else None)
    index={(*key(r),r['mark']):r for r in signals}
    contrasts=[]
    for r in signals:
        if r['state']!='iATCs' or r['mark']=='H3':continue
        for comparator in ['iAT2','iAT1']:
            other=index[(r['gene'],r['tss'],r['half_width_bp'],comparator,r['replicate'],r['mark'])]
            eligible=r['ratio_eligible'] and other['ratio_eligible']
            contrasts.append(dict(gene=r['gene'],gene_group=r['gene_group'],chrom=r['chrom'],strand=r['strand'],
                tss=r['tss'],is_baseline=r['is_baseline'],half_width_bp=r['half_width_bp'],mark=r['mark'],
                replicate=r['replicate'],contrast='iATCs_minus_'+comparator,eligible=eligible,
                effect=r['log2_mark_over_H3']-other['log2_mark_over_H3'] if eligible else None,
                raw_mark_effect=math.log2(r['mean_CPM']/other['mean_CPM']) if r['mean_CPM']>0 and other['mean_CPM']>0 else None,
                H3_effect=math.log2(r['H3_mean_CPM']/other['H3_mean_CPM']) if r['H3_mean_CPM']>0 and other['H3_mean_CPM']>0 else None))
    contrasts=pd.DataFrame(contrasts)
    summaries=[]
    for keys,d in contrasts.groupby(['gene','gene_group','mark','contrast','half_width_bp','replicate'],sort=False):
        baseline=d[d.is_baseline].iloc[0]; valid=d[d.eligible]; alt=d[~d.is_baseline]; alt_valid=alt[alt.eligible]
        compare=alt_valid if baseline.eligible else alt_valid.iloc[:0]
        summaries.append(dict(zip(['gene','gene_group','mark','contrast','half_width_bp','replicate'],keys),
            total_tss=len(d),alternative_tss=len(alt),eligible_tss=len(valid),held_tss=len(d)-len(valid),
            eligible_alternatives=len(alt_valid),baseline_eligible=bool(baseline.eligible),baseline_effect=baseline.effect,
            minimum_effect=valid.effect.min(),maximum_effect=valid.effect.max(),
            effect_span=valid.effect.max()-valid.effect.min(),
            direction_changes=int(sum(sign(v)!=sign(baseline.effect) for v in compare.effect)),
            direction_comparisons=len(compare),max_abs_change=max([abs(v-baseline.effect) for v in compare.effect],default=None)))
    summary=pd.DataFrame(summaries)
    out.mkdir(parents=True);figdir.mkdir(parents=True,exist_ok=True)
    for name,frame in [('tss_manifest',pd.DataFrame(positions)),('transcript_eligibility',pd.DataFrame(audit)),
                       ('signals',pd.DataFrame(signals)),('contrasts',contrasts),('sensitivity_summary',summary)]:
        frame.to_csv(out/(name+'.tsv'),sep='\t',index=False)
    plt.rcParams.update({'font.size':9,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    genes=[r['gene'] for r in loci]; y=np.arange(len(genes))
    fig,axes=plt.subplots(1,3,figsize=(15,11),sharey=True)
    for ax,mark in zip(axes,['H3K27ac','H3K4me3','H3K27me3']):
        for offset,cut,color in [(-.13,'CUT1','#235789'),(.13,'CUT2','#c96a2b')]:
            d=summary[(summary.mark==mark)&(summary.contrast=='iATCs_minus_iAT2')&(summary.half_width_bp==1000)&(summary.replicate==cut)].set_index('gene').reindex(genes)
            ax.hlines(y+offset,d.minimum_effect,d.maximum_effect,color=color,lw=2)
            ax.scatter(d.baseline_effect,y+offset,c=color,s=16,label=cut,zorder=3)
        ax.axvline(0,color='gray',lw=.7,ls='--');ax.set_title(mark);ax.set_xlabel('log2 mark/H3 contrast')
        ax.set_yticks(y,genes);ax.grid(axis='y',alpha=.1)
    axes[0].invert_yaxis();axes[2].legend(frameon=False,loc='lower right')
    fig.suptitle('Alternative transcript starts: iATCs versus iAT2, ±1 kb promoters',fontsize=16)
    fig.text(.5,.025,'Dots: original gene-boundary TSS. Lines: range over all eligible annotated TSSs, not confidence intervals.\nOne iPSC line; two preparations. Held ratios remain missing. Both contrasts and window sizes are retained in the tables.',ha='center',fontsize=10)
    fig.tight_layout(rect=(0,.08,1,.95))
    paths=[]
    for ext in ['png','svg']:
        p=figdir/f'a1_histone_tss_robustness.{ext}';fig.savefig(p,dpi=180,bbox_inches='tight');paths.append(p)
    plt.close(fig)
    rp.write_text(json.dumps(dict(status='completed_descriptive_TSS_sensitivity',utc=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(Path(__file__)),helper_sha256=sha(Path(__file__).with_name('a1_robustness.py')),
        interval_source_sha256=sha(Path(__file__).with_name('15_quantify_direct_marks.py')),
        input_sha256=inputs,original_output_hashes_verified=len(original['output_sha256']),
        positions=len(positions),signals=len(signals),contrasts=len(contrasts),
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in [*out.iterdir(),*paths]}),indent=2)+'\n')
    print('Histone sensitivity completed:',len(positions),'unique gene/TSSs;',len(signals),'signals;',len(contrasts),'contrasts')


if __name__=='__main__':main()
