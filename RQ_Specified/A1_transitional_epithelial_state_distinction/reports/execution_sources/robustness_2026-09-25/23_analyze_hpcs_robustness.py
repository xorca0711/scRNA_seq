"""Descriptive source influence, partition concordance and design identifiability."""
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from a1_robustness import sha, weighted_summary, adjusted_rand

BASE = Path(__file__).resolve().parents[1]


def main():
    out = BASE/'tables/robustness_2026-09-25/hpcs'
    figdir = BASE/'figures/robustness_2026-09-25'
    rp = BASE/'reports/hpcs_robustness_run.json'
    if out.exists() or rp.exists():
        raise SystemExit('Refusing to overwrite HPCS robustness')
    cp = BASE/'config/robustness_batch.json'
    cfg = json.loads(cp.read_text())['hpcs']
    source_cfg = json.loads((BASE/cfg['reuse_contract']).read_text())
    original = json.loads((BASE/'reports/hpcs_source_composition_run.json').read_text())
    for name, digest in original['output_sha256'].items():
        assert sha(BASE/name) == digest, name
    obs_path = BASE/cfg['input']
    assert sha(obs_path) == original['input_sha256'][cfg['input']]
    obs = pd.read_csv(obs_path, index_col=0, dtype=str, keep_default_na=False)
    data = obs[obs.batch.isin(source_cfg['select_batches'])].copy()
    assert len(data) == 5333 and data.index.is_unique
    assert set(data.Sorting_Groups) == {'traced'} and (data.Reporter == '').all()
    manifest_path = BASE/'tables/hpcs_source_composition/source_manifest.tsv'
    manifest = pd.read_csv(manifest_path, sep='\t')
    assert manifest.source_label.is_unique and len(manifest) == 22
    assert set(manifest.source_label) == set(data.Classification)
    for row in manifest.itertuples():
        d = data[data.Classification == row.source_label]
        assert len(d) == row.retained_cells and set(d.Group) == {row.group}
    fields = cfg['annotations']
    data[fields] = data[fields].replace('', cfg['missing_label'])
    counts, summary, omissions = [], [], []
    for field in fields:
        categories = sorted(data[field].unique())
        for group in source_cfg['groups']:
            sources = sorted(manifest.loc[manifest.group == group, 'source_label'])
            for category in categories:
                values = []
                for source in sources:
                    d = data[data.Classification == source]
                    n, den = int((d[field] == category).sum()), len(d)
                    counts.append(dict(annotation=field, group=group, source_label=source,
                                       category=category, numerator=n, denominator=den, fraction=n/den))
                    values.append((n, den))
                equal, pooled = weighted_summary(values)
                local = []
                for i, source in enumerate(sources):
                    e, p = weighted_summary(values[:i] + values[i+1:])
                    row = dict(annotation=field, group=group, category=category, omitted_source=source,
                               remaining_sources=len(sources)-1, remaining_cells=sum(v[1] for j,v in enumerate(values) if i != j),
                               equal_source=e, cell_pooled=p, equal_change_pp=100*(e-equal), pooled_change_pp=100*(p-pooled))
                    omissions.append(row); local.append(row)
                summary.append(dict(annotation=field, group=group, category=category, sources=len(sources),
                    total_cells=sum(d for n,d in values), state_cells=sum(n for n,d in values),
                    equal_source=equal, cell_pooled=pooled, weighting_gap_pp=100*(equal-pooled),
                    loo_equal_min=min(r['equal_source'] for r in local), loo_equal_max=max(r['equal_source'] for r in local),
                    loo_pooled_min=min(r['cell_pooled'] for r in local), loo_pooled_max=max(r['cell_pooled'] for r in local),
                    max_abs_equal_change_pp=max(abs(r['equal_change_pp']) for r in local),
                    max_abs_pooled_change_pp=max(abs(r['pooled_change_pp']) for r in local)))
    partitions, joint = [], []
    scopes = [('all', 'all', data)]
    scopes += [('group', str(k), v) for k,v in data.groupby('Group')]
    scopes += [('source', str(k), v) for k,v in data.groupby('Classification')]
    for level, unit, d in scopes:
        for left, right in combinations(fields, 2):
            same_codes = left == 'clusterK12' and right == 'clusterK12_stringent'
            partitions.append(dict(level=level, unit=unit, left=left, right=right, cells=len(d),
                adjusted_rand=adjusted_rand(d[left].tolist(), d[right].tolist()),
                literal_code_disagreements=int((d[left] != d[right]).sum()) if same_codes else None,
                literal_code_disagreement_fraction=float((d[left] != d[right]).mean()) if same_codes else None,
                missing_left=int((d[left] == cfg['missing_label']).sum()), missing_right=int((d[right] == cfg['missing_label']).sum())))
            ct = pd.crosstab(d[left], d[right]).reindex(index=sorted(data[left].unique()), columns=sorted(data[right].unique()), fill_value=0)
            for a in ct.index:
                for b in ct.columns:
                    joint.append(dict(level=level, unit=unit, left=left, right=right,
                                      left_category=a, right_category=b, cells=int(ct.loc[a,b])))
    designs = []
    for driver in ['all', 'Slc4a11', 'Hopx']:
        d = manifest if driver == 'all' else manifest[manifest.driver == driver]
        base = pd.get_dummies(d.source_library, dtype=int).to_numpy()
        for effect in ['chase_days', 'driver']:
            extra = pd.get_dummies(d[effect].astype(str), dtype=int).to_numpy()
            rank = int(np.linalg.matrix_rank(base)); combined = int(np.linalg.matrix_rank(np.column_stack([base, extra])))
            designs.append(dict(scope=driver, effect=effect, source_labels=len(d), libraries=d.source_library.nunique(),
                levels=d[effect].nunique(), library_rank=rank, combined_rank=combined, added_rank=combined-rank,
                status='no_contrast' if d[effect].nunique()<2 else ('held_library_alias' if combined==rank else 'algebraically_estimable_only'),
                independent_biological_units_verified=False))
    out.mkdir(parents=True); figdir.mkdir(parents=True, exist_ok=True)
    frames = {'source_annotation_counts':pd.DataFrame(counts), 'group_influence_summary':pd.DataFrame(summary),
              'leave_one_source_out':pd.DataFrame(omissions), 'partition_agreement':pd.DataFrame(partitions),
              'annotation_contingencies':pd.DataFrame(joint), 'design_rank_audit':pd.DataFrame(designs)}
    frames['library_group_support'] = manifest.groupby(['source_library','group','driver','chase_days'],as_index=False).agg(source_labels=('source_label','nunique'),retained_cells=('retained_cells','sum'))
    for name, frame in frames.items():
        frame.to_csv(out/(name+'.tsv'), sep='\t', index=False)
    plt.rcParams.update({'font.size':10, 'svg.fonttype':'none', 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1,2,figsize=(14,5.6))
    d = frames['group_influence_summary'].query("annotation == 'cell type' and category == 'HPCS'").set_index('group').loc[list(source_cfg['groups'])]
    y = np.arange(len(d))
    for offset, name, color in [(-.1,'equal','#235789'),(.1,'pooled','#c96a2b')]:
        center = d['equal_source' if name=='equal' else 'cell_pooled']*100
        axes[0].hlines(y+offset,d[f'loo_{name}_min']*100,d[f'loo_{name}_max']*100,color=color,lw=2)
        axes[0].scatter(center,y+offset,c=color,label='Equal-source' if name=='equal' else 'Cell-pooled',s=35,zorder=3)
    axes[0].set_yticks(y,d.index); axes[0].invert_yaxis(); axes[0].set_xlabel('Author HPCS label (% of retained cells)')
    axes[0].set_title('Source omission sensitivity'); axes[0].legend(loc='lower right',frameon=False)
    d = frames['partition_agreement'].query("level == 'group'").copy()
    d['pair'] = d.left+' / '+d.right
    matrix=d.pivot(index='unit',columns='pair',values='adjusted_rand').reindex(list(source_cfg['groups']))
    im=axes[1].imshow(matrix.to_numpy(),vmin=-1,vmax=1,cmap='coolwarm',aspect='auto')
    axes[1].set_xticks(range(3),['State / K12','State / stringent','K12 / stringent'],rotation=20,ha='right')
    axes[1].set_yticks(y,matrix.index); axes[1].set_title('Partition agreement (ARI)')
    for i,row in enumerate(matrix.to_numpy()):
        for j,value in enumerate(row): axes[1].text(j,i,f'{value:.2f}',ha='center',va='center')
    fig.colorbar(im,ax=axes[1],fraction=.04,pad=.02)
    fig.suptitle('HPCS: robustness of deposited source and annotation summaries',fontsize=15)
    fig.text(.5,.025,'Bars are source-omission ranges, not confidence intervals. Related labels are not independent validation.\nSource aliases are not verified mice; chase and source library cannot be separated in this design.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.13,1,.93))
    figure_paths=[]
    for ext in ['png','svg']:
        p=figdir/f'a1_hpcs_robustness.{ext}';fig.savefig(p,dpi=180,bbox_inches='tight');figure_paths.append(p)
    plt.close(fig)
    inputs=[cp,BASE/cfg['reuse_contract'],BASE/'config/robustness_annotation_scope.md',obs_path,manifest_path,
            BASE/'reports/hpcs_source_composition_run.json',BASE/'reports/robustness_source_inventory.json']
    rp.write_text(json.dumps(dict(status='completed_descriptive_robustness',utc=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(Path(__file__)), helper_sha256=sha(Path(__file__).with_name('a1_robustness.py')),
        input_sha256={p.relative_to(BASE).as_posix():sha(p) for p in inputs},
        original_output_hashes_verified=len(original['output_sha256']), annotation_mapping='biological recoding held; label-invariant partition audit completed',
        selected_cells=5333,source_aliases=22,
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in [*out.iterdir(),*figure_paths]}),indent=2)+'\n')
    print('HPCS robustness completed:', len(summary), 'group/category summaries;', len(omissions), 'source omissions')


if __name__ == '__main__':
    main()
