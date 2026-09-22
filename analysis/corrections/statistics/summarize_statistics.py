"""Assemble correction tables, subtype contributions, provenance and one figure."""
from __future__ import annotations
import hashlib
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
site = ROOT / '.venv-x64/Lib/site-packages'
if site.is_dir():
    sys.path.insert(0, str(site))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import false_discovery_control

CACHE = HERE / 'cache'
TABLES = HERE / 'tables'
TRIALS = ROOT / 'Thesis/gate1_01_niethamer_2025/trials'


def read_sets():
    result = {}
    for name in ['h.all.v2024.1.Hs.symbols.gmt', 'c5.go.bp.v2024.1.Hs.symbols.gmt']:
        for line in (ROOT / 'raw_data/msigdb' / name).read_text().splitlines():
            fields = line.split('\t')
            result[fields[0]] = set(fields[2:])
    return result


def contributions(frozen):
    sets = read_sets()
    rows = []
    cellrows = []
    for cohort in ['GSE136831', 'GSE135893']:
        units = pd.read_csv(CACHE / f'{cohort}_subtypes_units.csv')
        whole = pd.read_csv(CACHE / f'{cohort}_units.csv')
        arr = np.load(CACHE / f'{cohort}_subtypes.npz')
        counts, genes = arr['counts'], arr['genes']
        for comp in whole.comp.unique():
            totalunits = whole[whole.comp == comp].reset_index(drop=True)
            labels = sorted(units.loc[units.comp == comp, 'label'].unique())
            for label in labels:
                su = units[(units.comp == comp) & (units.label == label)]
                positions = {r.donor: i for i, r in su.iterrows()}
                fractions = np.array([units.loc[positions[r.donor], 'cells'] / r.cells if r.donor in positions else 0.0 for r in totalunits.itertuples()])
                ipf = (totalunits.disease == 'IPF').to_numpy()
                cellrows.append({'cohort': cohort, 'compartment': comp, 'label': label,
                    'n_IPF': int(ipf.sum()), 'n_control': int((~ipf).sum()),
                    'mean_cell_fraction_IPF': float(fractions[ipf].mean()),
                    'mean_cell_fraction_control': float(fractions[~ipf].mean()),
                    'delta_cell_fraction': float(fractions[ipf].mean() - fractions[~ipf].mean())})
                for name in frozen.loc[frozen.compartment == comp, 'set']:
                    hit = np.isin(genes, list(sets[name]))
                    gene_counts = counts[:, hit].sum(1)
                    denominator = units.assign(set_counts=gene_counts).query('comp == @comp').groupby('donor').set_counts.sum()
                    shares = np.array([gene_counts[positions[r.donor]] / denominator[r.donor]
                        if r.donor in positions and denominator[r.donor] > 0 else 0.0 for r in totalunits.itertuples()])
                    rows.append({'cohort': cohort, 'compartment': comp, 'label': label, 'set': name,
                        'mean_set_transcript_share_IPF': float(shares[ipf].mean()),
                        'mean_set_transcript_share_control': float(shares[~ipf].mean()),
                        'delta_transcript_share': float(shares[ipf].mean() - shares[~ipf].mean())})
    pd.DataFrame(cellrows).to_csv(TABLES / 'subtype_cell_fractions.csv', index=False)
    pd.DataFrame(rows).to_csv(TABLES / 'subtype_transcript_contributions.csv', index=False)


def main():
    frozen = pd.read_csv(TRIALS / 'g2_gsea_ipf/g2_replication.csv')
    merged = frozen[['compartment', 'collection', 'set', 'discovery_nes', 'held_out_nes', 'discovery_fdr', 'held_out_fdr', 'replicates']].copy()
    multiplicity = []
    for cohort, short in [('GSE136831', 'discovery'), ('GSE135893', 'validation')]:
        frame = pd.read_csv(TABLES / f'{cohort}_camera_candidates.csv')
        allsets = pd.read_csv(CACHE / f'{cohort}_camera_allsets.csv.gz')
        allsets['collection'] = np.where(allsets['set'].str.startswith('HALLMARK_'), 'hallmark', 'gobp')
        for (comp, coll), group in allsets.groupby(['compartment', 'collection']):
            q = false_discovery_control(group.PValue.to_numpy())
            allsets.loc[group.index, 'FDR_original_family'] = q
            multiplicity.append({'cohort': cohort, 'compartment': comp, 'collection': coll,
                'tested_sets': len(group), 'min_FDR': float(q.min()), 'hits_FDR05': int((q < .05).sum())})
        frame = frame.merge(allsets[['compartment', 'set', 'FDR_original_family']], on=['compartment', 'set'], how='left', validate='one_to_one')
        keep = ['compartment', 'set', 'eligible', 'Direction', 'Correlation', 'PValue', 'FDR_global', 'FDR_original_family', 'NGenes', 'n_first', 'n_second', 'residual_df']
        frame = frame[keep].rename(columns={c: f'{short}_{c}' for c in keep if c not in ['compartment', 'set']})
        merged = merged.merge(frame, on=['compartment', 'set'], validate='one_to_one')
        fixed = pd.read_csv(TABLES / f'{cohort}_fixed001_camera_candidates.csv')
        fixed_keep = ['compartment','set','eligible','Direction','PValue','FDR_global']
        fixed = fixed[fixed_keep].rename(columns={c:f'{short}_fixed001_{c}' for c in fixed_keep if c not in ['compartment','set']})
        merged = merged.merge(fixed,on=['compartment','set'],validate='one_to_one')
        lodopath = TABLES / f'{cohort}_lodo_summary.csv'
        if lodopath.exists():
            lod = pd.read_csv(lodopath).drop(columns='cohort')
            lod = lod.rename(columns={c: f'{short}_lodo_{c}' for c in lod if c not in ['compartment', 'set']})
            merged = merged.merge(lod, on=['compartment', 'set'], how='left', validate='one_to_one')
    expected = np.where(merged.discovery_nes > 0, 'Up', 'Down')
    merged['corrected_replicated'] = ((merged.discovery_FDR_global < .05) & (merged.validation_FDR_global < .05)
        & (merged.discovery_Direction == expected) & (merged.validation_Direction == expected))
    merged['fixed001_replicated'] = ((merged.discovery_fixed001_FDR_global < .05) & (merged.validation_fixed001_FDR_global < .05)
        & (merged.discovery_fixed001_Direction == expected) & (merged.validation_fixed001_Direction == expected))
    merged.to_csv(TABLES / 'g2_corrected_replication.csv', index=False)
    pd.DataFrame(multiplicity).to_csv(TABLES / 'multiplicity_sensitivity.csv', index=False)
    contributions(frozen)
    w1 = pd.read_csv(TABLES / 'w1_reference_camera.csv')
    g1 = pd.read_csv(TABLES / 'g1_confounded_sensitivity.csv')
    overview = {
        'frozen_discovery_candidates': len(merged),
        'historical_replicated': int(merged.replicates.sum()),
        'corrected_replicated': int(merged.corrected_replicated.sum()),
        'fixed001_replicated': int(merged.fixed001_replicated.sum()),
        'fixed001_and_historical_replicated': int((merged.fixed001_replicated & merged.replicates).sum()),
        'discovery_candidates_eligible': int(merged.discovery_eligible.sum()),
        'validation_candidates_eligible': int(merged.validation_eligible.sum()),
        'discovery_candidates_global_FDR05': int((merged.discovery_FDR_global < .05).sum()),
        'validation_candidates_global_FDR05': int((merged.validation_FDR_global < .05).sum()),
        'discovery_hits_under_original_multiplicity_families': int(sum(r['hits_FDR05'] for r in multiplicity if r['cohort'] == 'GSE136831')),
        'by_compartment': merged.groupby('compartment')[['replicates', 'corrected_replicated', 'fixed001_replicated']].sum().astype(int).to_dict('index'),
        'w1': {f'{method}/{contrast}': {'sets_tested': len(f), 'hits': int((f.FDR_family < .05).sum()),
                'min_FDR': float(f.FDR_family.min()), 'max_FDR': float(f.FDR_family.max())}
            for (method, contrast), f in w1.groupby(['method', 'contrast'])},
        'g1_frozen_DNA_replication': g1[g1['set'].isin(['GOBP_MITOTIC_DNA_REPLICATION', 'GOBP_CELL_CYCLE_DNA_REPLICATION'])].to_dict('records')
    }
    overview['lodo'] = {}
    for short in ['discovery','validation']:
        col = f'{short}_lodo_direction_reversals'
        if col in merged:
            overview['lodo'][short] = {'sets_evaluated':int(merged[col].notna().sum()),
                'sets_with_any_direction_reversal':int((merged[col]>0).sum()),
                'historical_replicated_sets_with_reversal':int(((merged[col]>0)&merged.replicates).sum()),
                'sets_with_below_floor_omissions':int((merged[f'{short}_lodo_below_floor_omissions']>0).sum())}
    (HERE / 'summary.json').write_text(json.dumps(overview, indent=2) + '\n')

    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.6), gridspec_kw={'width_ratios': [1.1, 1, 1.25]})
    comps = ['AT2','fibroblasts','macrophages']
    counts = merged.groupby('compartment')[['replicates','corrected_replicated','fixed001_replicated']].sum().loc[comps]
    for offset, col, label, color in [(-.25,'replicates','Historical GSEA','#b1b6ba'),
        (0,'corrected_replicated','CAMERA estimated ρ','#007f8b'),(.25,'fixed001_replicated','CAMERA fixed ρ=.01','#b96d22')]:
        bars = axs[0].bar(np.arange(3)+offset,counts[col],width=.24,color=color,label=label)
        axs[0].bar_label(bars,padding=2,fontsize=7)
    axs[0].set_xticks(range(3),['AT2','Fibroblasts','Macrophages'])
    axs[0].set(ylabel='Sets passing in both cohorts',ylim=(0,max(counts.max())*1.28),
        title='G2: replicated inference depends on\nthe correlation assumption')
    axs[0].legend(frameon=False,fontsize=7,loc='upper left')
    for i, contrast in enumerate(['A', 'B']):
        subset = w1[(w1.method == 'voom_TMM') & (w1.contrast == contrast)]
        axs[1].scatter(np.full(len(subset), i), subset.FDR_family, s=34, color=['#007f8b', '#b96d22'][i], alpha=.7)
    axs[1].set_xticks([0,1], ['Repair vs\nresolution', 'Resolution vs\nlong-term'])
    axs[1].set(xlim=(-.5,1.5),ylim=(0,1),ylabel='Official CAMERA + voom FDR',title='W1: no set crosses the threshold')
    axs[1].axhline(.05,color='#b94f48',ls='--',lw=1)
    specs = ['all_unadjusted', 'all_adjusted', 'heterozygote_adjusted', 'december_heterozygote']
    g1fixed = pd.read_csv(TABLES / 'g1_confounded_sensitivity_fixed001.csv')
    for name, color in [('GOBP_CELL_CYCLE_DNA_REPLICATION','#007f8b'), ('GOBP_MITOTIC_DNA_REPLICATION','#7560a6')]:
        vals = g1[g1['set'] == name].set_index('specification').loc[specs]
        fixedvals = g1fixed[g1fixed['set'] == name].set_index('specification').loc[specs]
        shortname = 'Cell-cycle DNA' if 'CELL_CYCLE' in name else 'Mitotic DNA'
        axs[2].plot(range(4), vals.FDR_family, 'o:', color=color, label=shortname+'; estimated ρ')
        axs[2].plot(range(4), fixedvals.FDR_family, 's-', color=color, label=shortname+'; fixed ρ=.01')
    axs[2].set_xticks(range(4), ['Unadjusted','Batch/sex/\ngenotype','Heterozygote\nadjusted','Same-round\nheterozygotes'])
    axs[2].tick_params(axis='x',labelsize=8)
    axs[2].axhline(.05,color='#b94f48',ls='--',lw=1)
    axs[2].set(ylabel='Official CAMERA full-family FDR',ylim=(-.03,1.08),title='G1: neither correlation setting supports\nDNA-replication sets after adjustment')
    axs[2].legend(frameon=False,fontsize=6.5,loc='upper left',bbox_to_anchor=(.025,.80))
    fig.text(.015,.015,'Retrospective sensitivity; non-significance is not equivalence. G1 age and injury time remain inseparable.',fontsize=8,color='#555555')
    fig.tight_layout(rect=[0,.045,1,1])
    fig.savefig(HERE / 'statistical_corrections.png',dpi=180)
    plt.close(fig)

    originals = [TRIALS / p for p in ['gsea_utils.py','g1_gsea_by_phase.py','g2_gsea_ipf.py','w1_amac_pseudobulk_de.py','g2_gsea_ipf/g2_replication.csv']]
    code = list(HERE.glob('*.py')) + list(HERE.glob('*.R')) + list(HERE.glob('*.ps1')) + [HERE / 'PROTOCOL.md',HERE / 'README.md']
    outputs = list(TABLES.glob('*.csv')) + [HERE / 'summary.json',HERE / 'statistical_corrections.png',HERE / 'R-session.txt']
    def facts(p):
        return {'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
    record = {'completed_utc': datetime.now(timezone.utc).isoformat(), 'analysis_type':'retrospective correction and sensitivity',
        'python':{'executable':sys.executable,'version':platform.python_version(),
            'packages':{p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','anndata','h5py','matplotlib']}},
        'R_runtime':'official CRAN R 4.6.1, portable extraction; package versions in R-session.txt',
        'runtime_sources':{'R':'https://cran.r-project.org/bin/windows/base/R-4.6.1-win.exe',
            'R_official_md5':'7907f3a20ec8ec88cd0da279024b8e27',
            'innoextract':'https://constexpr.org/innoextract/files/innoextract-1.9-windows.zip',
            'innoextract_official_md5':'72d0d0dd874b6236eaa44411f4470ee1'},
        'method':'official limma CAMERA, estimated residual inter-gene correlation primary; fixed .01 sensitivity; voom/TMM primary; W1 historical-transform reference also retained',
        'code_snapshot_scope':'Final saved sources, including secondary fixed-correlation option added after primary runs; primary numerical algorithm unchanged by that addition.',
        'execution_recovery':[{'stage':'LODO','completed_refits':139,
            'issue':'After both cohort summaries had been written, the R interpreter encountered a trailing parse error because its open script was extended while the long loop ran.',
            'resolution':'Final saved R source separately parsed successfully; official package session metadata refreshed; output completeness checked against every donor and candidate. Numerical refits were complete and not rerun.'}],
        'prepared_inputs': [json.loads(p.read_text()) for p in sorted(HERE.glob('*_preparation.json'))],
        'gene_set_inputs':[facts(p) for p in sorted((ROOT / 'raw_data/msigdb').glob('*.gmt'))],
        'historical_inputs': [facts(p) for p in originals], 'code': [facts(p) for p in code], 'outputs':[facts(p) for p in outputs],
        'unresolved':['age vs injury time in mouse','unavailable human age/sex covariates','human subtype label harmonization','ambient RNA and doublets','selection from 50-cell floor','observational causality']}
    (HERE / 'run_record.json').write_text(json.dumps(record,indent=2) + '\n')
    print(json.dumps({k:v for k,v in overview.items() if k not in ['g1_frozen_DNA_replication','w1']},indent=2))


if __name__ == '__main__':
    main()
