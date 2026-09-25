"""Independently verify paired calculations, then render and report completed tests."""
import hashlib
import importlib.util
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
ROOT=HERE.parents[1]
A5=ROOT/'RQ_Specified/A5_developmental_programme_reuse'
A11=ROOT/'RQ_Specified/A11_lesion_programme_addition'


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def main():
    import numpy as np
    import pandas as pd
    from scipy import stats
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a5=A5/'tables/test_v1';a11=A11/'tables/test_v2'
    r5=pd.read_csv(a5/'inference.tsv',sep='\t')
    r11=pd.read_csv(a11/'inference.tsv',sep='\t')
    d5=pd.read_csv(a5/'paired_differences.tsv',sep='\t')
    d11=pd.read_csv(a11/'paired_differences.tsv',sep='\t')
    checks=[]
    def check(name, condition, **detail):
        assert bool(condition),name
        checks.append({'check':name,'pass':True,**detail})
    for name,path,rows,keys,score_key,case_key,ref_key in [
        ('A5',a5,d5,['sample_id','module'],'score','Krt8+ ADI',None),
        ('A11',a11,d11,['patient','module'],'score','lesion','normal')]:
        scores=pd.read_csv(path/('unit_scores.tsv' if name=='A5' else 'scores.tsv'),sep='\t')
        key='label' if name=='A5' else 'arm'
        pivot=scores.pivot(index=keys,columns=key,values=score_key)
        errors=[]
        for row in rows.to_dict('records'):
            pair=pivot.loc[tuple(row[k] for k in keys)]
            expected=pair[case_key]-pair[row['reference'] if ref_key is None else ref_key]
            got=row['difference_pp' if name=='A5' else 'difference']
            errors.append(abs(expected-got))
        check(name+' paired arithmetic',max(errors)<1e-10,n=len(errors),max_error=max(errors))
    for row in r5.to_dict('records'):
        d=d5.loc[(d5.module==row['module']) & (d5.reference==row['reference']),'difference_pp'].to_numpy()
        z=stats.ttest_1samp(d,0);ci=z.confidence_interval()
        check('A5 '+row['test']+' independent t inference',
              abs(z.pvalue-row['p_mean'])<1e-12 and abs(ci.low-row['mean_low'])<1e-10 and abs(ci.high-row['mean_high'])<1e-10)
    def dd(k):return d11[d11.module==k].set_index('patient').difference.sort_index()
    tests={'primary':dd('lesion_specific'),'beyond_shared':dd('lesion_specific')-dd('shared_remodelling'),
           'stress_excluded':dd('stress_excluded'),'injury_lesion_pair':dd('injury_lesion_pair')}
    for row in r11.to_dict('records'):
        d=tests[row['test']].to_numpy(); n=len(d)
        check('A11 '+row['test']+' no zero/tied absolute differences',len(np.unique(abs(d)))==n and (d!=0).all())
        ranks=stats.rankdata(abs(d)); observed=ranks[d>0].sum()
        # Enumerate every possible sign assignment: independent check of exact p.
        signs=((np.arange(2**n)[:,None] >> np.arange(n)) & 1)
        values=signs@ranks
        p=min(1.,2*min(np.mean(values<=observed),np.mean(values>=observed)))
        walsh=(d[:,None]+d[None,:])[np.triu_indices(n)]/2
        check('A11 '+row['test']+' exact sign enumeration and HL',
              abs(p-row['p_exact'])<1e-12 and abs(np.median(walsh)-row['HL'])<1e-8,
              permutations=2**n)
    def adjusted(ps,method):
        ps=np.array(ps);order=np.argsort(ps);sorted_p=ps[order];n=len(ps)
        vals=np.minimum.accumulate((sorted_p*n/np.arange(1,n+1))[::-1])[::-1] if method=='BH' else np.maximum.accumulate(sorted_p*np.arange(n,0,-1))
        ans=np.empty(n);ans[order]=np.minimum(vals,1);return ans
    check('A5 Holm family',np.allclose(adjusted(r5.p_mean.iloc[1:],'Holm'),r5.q_secondary.iloc[1:],atol=1e-12,rtol=0))
    check('A11 BH family',np.allclose(adjusted(r11.p_exact.iloc[1:],'BH'),r11.q_secondary.iloc[1:],atol=1e-12,rtol=0))
    spec=importlib.util.spec_from_file_location('a5_score',A5/'scripts/03_score_external_test.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    N=np.array([500,500,501,1000,1000,2000,10000,100000]);k=np.array([0,1,1,1,50,250,10,100])
    error=np.max(abs(mod.detection_probability(N,k,500)-(1-stats.hypergeom.pmf(0,N,k,500))))
    check('Exact-depth expectation vs hypergeometric PMF',error<1e-8,max_error=float(error))
    instrument=pd.read_csv(a11/'instrument_check.tsv',sep='\t')
    check('Original human instrument',len(instrument)==46 and instrument.abs_error.max()<=1e-6,max_error=float(instrument.abs_error.max()))
    for folder,record in [(a5,'score_run.json'),(a11,'prepare_run.json')]:
        run=json.loads((folder/record).read_text())
        for path,h in run['outputs'].items():
            p=folder/path if folder==a5 else A11/path
            check('Output hash '+str(p.relative_to(ROOT)),sha(p)==h)
    # Preserve original data-derived artifacts, including all shared frozen outputs.
    changed=subprocess.check_output(['git','diff','--name-only','origin/main'],cwd=ROOT,text=True).splitlines()
    protected=['RQ_Specified/A5_A11_shared_component_contract/tables/'+n for n in
               ['frozen_modules.json','module_membership.tsv','module_overlap.tsv','freeze_run.json','coverage_gate.tsv']]
    check('Original frozen artifacts unchanged',not(set(changed)&set(protected)))
    figures=HERE/'figures';figures.mkdir(exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(13.6,5.3),layout='constrained')
    labels5=['Full external signature','Without Guo identity genes','Without identity + controls','Full signature vs resting AT2']
    labels11=['Lesion-derived primary','Lesion change minus shared','Lesion, stress genes excluded','Injury–lesion pair']
    for ax,rs,labels,x,lo,hi,color,title in [
        (axes[0],r5,labels5,'mean','mean_low','mean_high','#176B70','A5 · developmental-gene recruitment'),
        (axes[1],r11,labels11,'HL','low','high','#8356A1','A11 · lesion-associated expression')]:
        y=np.arange(len(rs))[::-1]
        ax.errorbar(rs[x],y,xerr=np.vstack([rs[x]-rs[lo],rs[hi]-rs[x]]),fmt='o',color=color,capsize=4,markersize=7)
        ax.axvline(0,color='#7D8992',lw=1)
        ax.set_yticks(y,[f'{label}\n(n={n})' for label,n in zip(labels,rs.n)])
        ax.set_ylim(-.6,3.6);ax.set_title(title,loc='left',fontweight='bold',pad=16)
        ax.grid(axis='x',alpha=.15)
    axes[0].set_xlabel('Mean paired difference (detection percentage points)\n95% t intervals; 500-UMI expectation')
    axes[1].set_xlabel('Paired location shift (log2 CPM)\nHodges–Lehmann; exact 95% intervals')
    axes[1].axvline(.1,color='#C88938',ls=':',lw=1,label='Primary planning margin 0.10')
    axes[1].legend(loc='lower right',fontsize=8,frameon=False)
    fig.suptitle('Distinct questions and measurement scales; expression does not establish lineage or function',fontsize=12)
    for suffix in ['png','svg']:
        fig.savefig(figures/('a5_a11_results.'+suffix),dpi=180)
        if suffix=='svg':
            vector=figures/('a5_a11_results.'+suffix)
            vector.write_text('\n'.join(line.rstrip() for line in vector.read_text().splitlines())+'\n',encoding='utf-8',newline='\n')
    plt.close(fig)
    primary5=r5.iloc[0];primary11=r11.iloc[0]
    l5=pd.read_csv(a5/'omission_diagnostics.tsv',sep='\t')
    l11=pd.read_csv(a11/'omission_diagnostics.tsv',sep='\t')
    p5=d5[(d5.module=='Guo_AT1_AT2_external')&(d5.reference=='AT2 activated')]
    day=pd.read_csv(a5/'day_means.tsv',sep='\t')
    def fmt_table(rows,cols):
        lines=['| '+' | '.join(cols)+' |','|'+'---|'*len(cols)]
        for values in rows:lines.append('| '+' | '.join(map(str,values))+' |')
        return '\n'.join(lines)
    table5=fmt_table([[dict(zip(r5.test,labels5))[r.test],int(r.n),f'{r.mean:.3f}',f'[{r.mean_low:.3f}, {r.mean_high:.3f}]',f'{r.p_mean:.3g}',('—' if pd.isna(r.q_secondary) else f'{r.q_secondary:.3g}')] for r in r5.itertuples()],['Test','Mice','Mean pp','95% CI','p','Holm q'])
    table11=fmt_table([[dict(zip(r11.test,labels11))[r.test],int(r.n),f'{r.HL:.3f}',f'[{r.low:.3f}, {r.high:.3f}]',f'{r.p_exact:.4g}',('—' if pd.isna(r.q_secondary) else f'{r.q_secondary:.4g}')] for r in r11.itertuples()],['Test','Patients','HL log2 CPM','95% CI','Exact p','BH q'])
    shared_mean=float(dd('shared_remodelling').mean())
    direction_count=int((dd('lesion_specific')>0).sum())
    stronger=primary11.direction=='positive' and all((r11.set_index('test').loc[t,'HL']>0 and r11.set_index('test').loc[t,'q_secondary']<.05) for t in ['beyond_shared','stress_excluded'])
    magnitude_text={'meaningful_positive_supported':'positive shift exceeding the 0.10 planning margin supported', 'positive_margin_ruled_out':'positive shift of 0.10 ruled out', 'unresolved':'magnitude unresolved', 'unavailable':'magnitude inference unavailable'}[primary11.magnitude]
    report=f'''# A5/A11 revised tests: results and biological interpretation

The external developmental signature increases in adult transitional cells,
including after the planned external identity and control exclusions. A11's
primary direction is {primary11.direction}; the predefined stronger relative-activation
criterion is {'met' if stronger else 'not met'}. Interpret these as paired RNA
associations, with the distinct biological limits below.

25 September 2026. Plans and external A5 definitions were committed in `eb5317e`
before new expression scores. Implementation and the passed discovery check were
committed in `19fa116`. See [biological logic](../BIOLOGICAL_LOGIC.md),
[A5 plan](../../A5_developmental_programme_reuse/PLAN.md) and
[A11 plan](../../A11_lesion_programme_addition/PLAN.md).

## A5: partial developmental-signature recruitment

{table5}

Intervals are per-test; secondary decisions use the displayed Holm adjustment.

All {len(p5)} primary mice have positive paired differences. The equal-day mean is
{day.difference_pp.mean():.3f} pp; leave-one-mouse-out means range from
{l5.remaining_mean.min():.3f} to {l5.remaining_mean.max():.3f} pp. All 24 activated-AT2
and 26 resting-AT2 candidate mice still meet the fixed floors after raw depth >=500.
The primary and both external exclusion modules have complete source-gene coverage.

The source signature is recruited beyond activated AT2, including the 57-gene
identity-excluded and 53-gene identity/control-excluded versions. This weakens the
specific explanation that all signal comes from Guo-defined mature alveolar identity
or the six nominated Hallmarks. It does not remove unlisted stress programmes, ambient RNA or annotation effects,
validate clustering, establish lineage direction or measure repair function.
Effect size is a mean detection-probability difference at 500 UMI, not a fraction
of cells converted or a log-fold change. No meaningful-effect margin was invented.

Strunz's pre-established observation concerned poor overall correspondence of
developmental and injury signatures. A positive paired mean for an external
subset is compatible with that observation: some genes can be recruited without
the populations sharing global identity. The 94/51-gene Strunz-filtered modules
and pairwise overlaps are descriptive references, not independent confirmations.

## A11: report direction, magnitude and biological claim separately

{table11}

Intervals are per-test; secondary decisions use the displayed BH adjustment.

Primary direction: **{primary11.direction}**. Primary positive-margin classification:
**{magnitude_text}**. This is a pragmatic log2 CPM margin, not a clinical threshold.
The arithmetic mean sensitivity is {primary11['mean']:.3f} log2 CPM
(95% t interval [{primary11.mean_low:.3f}, {primary11.mean_high:.3f}]); it estimates
a mean rather than the primary location parameter. {direction_count}/8 patient
differences are positive. Omission means range from {l11.remaining_mean.min():.3f}
to {l11.remaining_mean.max():.3f}; omitting double-primary P0028 gives
{float(l11.loc[l11.omitted_patient=='P0028','remaining_mean'].iloc[0]):.3f}.
These are descriptive omissions of already-normalized scores.

The predefined stronger *narrow* criterion (positive primary plus positive
BH-significant beyond-shared and stress-excluded tests) is **{'met' if stronger else 'not met'}**.
The shared union's mean paired change is {shared_mean:.3f} log2 CPM; its own
direction must accompany any interpretation of the relative score contrast.
{'The result supports relative activation above this nominated shared score while surviving the specified stress-gene exclusions.' if stronger else 'A stronger claim of activation beyond the nominated shared response is not established by the predefined combined criterion.'}

This experiment compares pooled author tumour epithelium with normal AT2, a changed
population definition from discovery. It does not establish malignant identity or
specificity relative to non-neoplastic repair/fibrosis. A module's source-list
exclusivity and even a positive beyond-shared comparison do not establish a separate
mechanism. Three Hallmark exclusions do not exhaust stress/cycling biology.

## How the two results fit A1

A5 asks which externally defined developmental genes are recruited during repair.
A11 asks whether lesion-associated expression replicates and how it compares with
a nominated common response. Neither moves a cell along a presumed repair-to-cancer
trajectory. A1's observed opposing regional AP-1 directions and mixed-direction TP53 overlap
show why shared RNA must be separated from lineage, regulatory mediation and
successful repair. Existing evidence supplies that context without repeating pooled
chromatin analysis as if it were independent confirmation.

![Paired results](../figures/a5_a11_results.png)

## Verification and execution record

- All 46 original human broad/narrow discovery pairs reproduced within
  {instrument.abs_error.max():.3g}; the original all-histology TMM scope was retained.
- Paired arithmetic, A5 t intervals, A11 exact sign enumeration and HL estimates,
  multiplicity adjustment, depth expectation and output hashes were independently
  checked by this reporting script. Original frozen artifacts remain unchanged.
- Strunz's deposited matrix is cells-by-genes (32,559 by 24,051); the loader's first
  shape assertion stopped before scores. Explicit dimension-verified transposition
  corrected the orientation. A Python-to-R boolean parsing mismatch similarly
  stopped A5 inference before outputs; explicit parsing corrected it. Neither
  changed a population, gene set, threshold or test. Failed attempts produced no
  scientific result to replace.
- A5 source/raw retrieval and score records are under its `tables/`; Kim source,
  pseudobulk and normalization records are under `tables/test_v2/`. Raw inputs and
  large intermediates remain in ignored caches. R session files record versions.

## What remains scientifically unresolved

The requested transcriptional tests are complete. Independent replication in
another adult injury study would test generality of A5 beyond Strunz. A comparable
non-neoplastic human injury arm and supported epithelial identity would be needed
for A11 specificity. Neither new cohort acquisition nor a functional perturbation
is replaced by re-scoring existing pooled data. There is no result-dependent
subtype search, threshold change or automatic expansion into these new studies.
'''
    reports=HERE/'reports';reports.mkdir(exist_ok=True)
    (reports/'REVISED_TEST_RESULTS.md').write_text(report,encoding='utf-8',newline='\n')
    prepare=json.loads((a11/'prepare_run.json').read_text())
    source_root=Path(next(iter(prepare['inputs']))).parents[2]
    paper='Research Article/gate2_C3_yu_lee_choi_min_2026'
    instrument_inputs=[ROOT/paper/'trials/u5_human_niche/unc20_pooled_triad_units.csv',
        source_root/paper/'cache/u5_human_niche/unc20_pooled_triad_counts.csv.gz',
        ROOT/paper/'trials/u6_human_specificity/paired_program_values.csv',
        ROOT/paper/'trials/u6_human_specificity/human_module_genes.tsv',
        ROOT/paper/'trials/u5_human_paired_pathways.R']
    manifest={'verified_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,
              'instrument_inputs':{str(p):sha(p) for p in instrument_inputs},
              'script_sha256':sha(Path(__file__)),
              'analysis_code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in
                [A5/'scripts/03_score_external_test.py',A5/'scripts/04_inference.R',
                 A11/'scripts/02_reproduce_discovery.R',A11/'scripts/03_prepare_kim.py',A11/'scripts/04_score_kim.R',
                 HERE/'scripts/paired_inference.R']},
              'inference_outputs':{str(p.relative_to(ROOT)):sha(p) for folder in [a5,a11] for p in folder.iterdir() if p.is_file()}}
    (HERE/'tables/revised_tests_verification.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Verified',len(checks),'checks; report and two-panel figure written')


if __name__=='__main__':main()
