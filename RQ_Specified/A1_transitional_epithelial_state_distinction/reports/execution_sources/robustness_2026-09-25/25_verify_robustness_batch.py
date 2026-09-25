"""Verify source summaries, partition metrics and alternative-promoter arithmetic."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import json
import math
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score
from a1_robustness import sha

BASE = Path(__file__).resolve().parents[1]


def close(a,b):
    np.testing.assert_allclose(a,b,rtol=1e-9,atol=1e-10,equal_nan=True)


def main():
    rp=BASE/'reports/robustness_verification.json'
    out=BASE/'tables/robustness_2026-09-25/verification'
    if rp.exists() or out.exists():raise SystemExit('Refusing to overwrite verification')
    hashes={}
    for name in ['hpcs_source_composition_run','direct_mark_run','hpcs_robustness_run','histone_tss_robustness_run']:
        p=BASE/'reports'/f'{name}.json'; hashes[p.relative_to(BASE).as_posix()]=sha(p)
        record=json.loads(p.read_text())
        for path,digest in record['output_sha256'].items():
            assert sha(BASE/path)==digest,path
            hashes[path]=digest
    for name,script in [('hpcs_robustness_run','23_analyze_hpcs_robustness.py'),('histone_tss_robustness_run','24_analyze_histone_tss.py')]:
        record=json.loads((BASE/'reports'/f'{name}.json').read_text())
        assert sha(BASE/'scripts'/script)==record['code_sha256']
        assert sha(BASE/'scripts/a1_robustness.py')==record['helper_sha256']
        for path,digest in record['input_sha256'].items():
            # The complete BigWigs were hashed immediately before extraction by script 24.
            if not path.endswith('.bw'):assert sha(BASE/path)==digest,path
    h=BASE/'tables/robustness_2026-09-25/hpcs'
    read=lambda p:pd.read_csv(p,sep='\t')
    obs=pd.read_csv(BASE/'cache/followup_sources/hpcs_combined_obs.csv.gz',dtype=str,keep_default_na=False,index_col=0)
    data=obs[obs.batch.isin(['traced','Hopx-MACD_traced'])].copy()
    fields=['cell type','clusterK12','clusterK12_stringent']
    data[fields]=data[fields].replace('','UNASSIGNED')
    counts=read(h/'source_annotation_counts.tsv'); summaries=read(h/'group_influence_summary.tsv')
    omissions=read(h/'leave_one_source_out.tsv'); agreements=read(h/'partition_agreement.tsv')
    for row in counts.itertuples():
        d=data[data.Classification==row.source_label]
        assert row.denominator==len(d) and row.numerator==int((d[row.annotation]==str(row.category)).sum())
        close(row.fraction,row.numerator/row.denominator)
    for keys,d in counts.groupby(['annotation','source_label']):assert d.numerator.sum()==d.denominator.iloc[0]
    for row in omissions.itertuples():
        d=counts[(counts.annotation==row.annotation)&(counts.group==row.group)&(counts.category.astype(str)==str(row.category))&(counts.source_label!=row.omitted_source)]
        assert len(d)==row.remaining_sources and d.denominator.sum()==row.remaining_cells
        close([row.equal_source,row.cell_pooled],[d.fraction.mean(),d.numerator.sum()/d.denominator.sum()])
    for row in summaries.itertuples():
        d=counts[(counts.annotation==row.annotation)&(counts.group==row.group)&(counts.category.astype(str)==str(row.category))]
        loo=omissions[(omissions.annotation==row.annotation)&(omissions.group==row.group)&(omissions.category.astype(str)==str(row.category))]
        close([row.equal_source,row.cell_pooled,row.loo_equal_min,row.loo_equal_max,row.loo_pooled_min,row.loo_pooled_max],
              [d.fraction.mean(),d.numerator.sum()/d.denominator.sum(),loo.equal_source.min(),loo.equal_source.max(),loo.cell_pooled.min(),loo.cell_pooled.max()])
    for row in agreements.itertuples():
        d=data if row.level=='all' else data[data['Group' if row.level=='group' else 'Classification']==row.unit]
        close(row.adjusted_rand,adjusted_rand_score(d[row.left],d[row.right]))
    original=read(BASE/'tables/hpcs_source_composition/group_descriptive_summary.tsv')
    for row in original.itertuples():
        d=summaries[(summaries.annotation=='cell type')&(summaries.group==row.group)&(summaries.category==row.state)].iloc[0]
        close([row.source_mean_fraction,row.cell_pooled_fraction],[d.equal_source,d.cell_pooled])
    manifest=read(BASE/'tables/hpcs_source_composition/source_manifest.tsv')
    assert (manifest.groupby('source_library').chase_days.nunique()==1).all()
    ranks=read(h/'design_rank_audit.tsv')
    assert (ranks[ranks.effect=='chase_days'].added_rank==0).all()
    t=BASE/'tables/robustness_2026-09-25/histone'
    signals=read(t/'signals.tsv'); contrasts=read(t/'contrasts.tsv'); sensitivity=read(t/'sensitivity_summary.tsv')
    old=read(BASE/'tables/direct_marks_2026-09-25/histone_signals.tsv')
    old['half_width_bp']=old.window.map({'promoter':1000,'broad_promoter':5000})
    keys=['gene','gsm','half_width_bp']
    baseline=signals[signals.is_baseline].merge(old,on=keys,suffixes=('_new','_old'),validate='one_to_one')
    assert len(baseline)==len(old)==1104
    for col in ['tss','start','end','mean_CPM','covered_fraction','H3_mean_CPM','log2_mark_over_H3']:
        close(baseline[col+'_new'],baseline[col+'_old'])
    assert (baseline.ratio_eligible_new==baseline.ratio_eligible_old).all()
    control=signals[signals.mark=='H3'].set_index(['gene','tss','half_width_bp','state','replicate'])
    for row in signals.itertuples():
        other=control.loc[(row.gene,row.tss,row.half_width_bp,row.state,row.replicate)]
        eligible=min(row.covered_fraction,other.covered_fraction)>=.8 and min(row.mean_CPM,other.mean_CPM)>0
        assert row.ratio_eligible==eligible
        close(row.H3_mean_CPM,other.mean_CPM)
        close(row.log2_mark_over_H3,math.log2(row.mean_CPM/other.mean_CPM) if eligible else np.nan)
    for row in contrasts.itertuples():
        mask=(signals.gene==row.gene)&(signals.tss==row.tss)&(signals.half_width_bp==row.half_width_bp)&(signals.replicate==row.replicate)&(signals.mark==row.mark)
        d=signals[mask].set_index('state');a=d.loc['iATCs'];b=d.loc[row.contrast.removeprefix('iATCs_minus_')]
        assert row.eligible==bool(a.ratio_eligible and b.ratio_eligible)
        close(row.effect,a.log2_mark_over_H3-b.log2_mark_over_H3 if row.eligible else np.nan)
        if row.eligible:close(row.effect,row.raw_mark_effect-row.H3_effect)
    for row in sensitivity.itertuples():
        d=contrasts[(contrasts.gene==row.gene)&(contrasts.mark==row.mark)&(contrasts.contrast==row.contrast)&(contrasts.replicate==row.replicate)&(contrasts.half_width_bp==row.half_width_bp)]
        v=d[d.eligible];base=d[d.is_baseline].iloc[0];alt=d[(~d.is_baseline)&d.eligible]
        assert row.total_tss==len(d) and row.eligible_tss==len(v) and row.held_tss==len(d)-len(v)
        close([row.minimum_effect,row.maximum_effect],[v.effect.min(),v.effect.max()])
        assert row.direction_changes==(int((np.sign(alt.effect)!=np.sign(base.effect)).sum()) if base.eligible else 0)
    sys.path.insert(0,str(BASE/'cache/pybigtools'));import pybigtools
    positions=read(t/'tss_manifest.tsv');alt=positions[~positions.is_baseline].copy()
    alt['distance']=(alt.tss-alt.baseline_tss).abs()
    chosen=alt.sort_values(['gene','distance','tss']).groupby('gene').tail(1)
    checks=[]
    samples=read(BASE/'tables/direct_marks_2026-09-25/sample_manifest.tsv')
    for sample in samples.itertuples():
        bw=pybigtools.open(str(BASE/sample.source_path))
        for locus in chosen.itertuples():
            r=signals[(signals.gsm==sample.gsm)&(signals.gene==locus.gene)&(signals.tss==locus.tss)&(signals.half_width_bp==1000)].iloc[0]
            values=np.asarray(bw.values(r.chrom,int(r.start),int(r.end),missing=float('nan')))
            mean=np.nan_to_num(values,nan=0).mean(); coverage=np.isfinite(values).mean()
            close([mean,coverage],[r.mean_CPM,r.covered_fraction])
            checks.append(dict(gene=locus.gene,tss=locus.tss,gsm=sample.gsm,mean_CPM=mean,covered_fraction=coverage))
        bw.close()
    out.mkdir(parents=True)
    check_path=out/'alternative_tss_base_resolution.tsv';pd.DataFrame(checks).to_csv(check_path,sep='\t',index=False)
    rp.write_text(json.dumps(dict(status='passed',utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(Path(__file__)),
        verified_sha256=hashes,original_outputs_preserved=True,source_count_rows=len(counts),source_omissions=len(omissions),
        partition_metrics_checked_against_sklearn=len(agreements),original_histone_windows_reproduced=1104,
        histone_signal_rows=len(signals),histone_contrasts=len(contrasts),alternative_base_resolution_windows=len(checks),
        chase_library_alias_confirmed=True,
        output_sha256={check_path.relative_to(BASE).as_posix():sha(check_path)}),indent=2)+'\n')
    print('Verified:',len(omissions),'omissions;',len(agreements),'partition metrics;',len(signals),'signals;',len(checks),'alternate base-resolution windows')


if __name__=='__main__':main()
