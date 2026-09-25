"""Deposited source-block PCA and technical peak audit; no biological tests."""
import gzip
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parents[1]
ROOT = BASE.parents[1]


def main():
    if (BASE / 'reports/descriptive_run.json').exists():
        raise SystemExit('Existing descriptive run: archive its outputs before a new numerical run; use script 10 to render figures.')
    inventory = json.loads((BASE / 'reports/processed_input_inventory.json').read_text())['files']
    contract_path = BASE / 'config/first_batch.json'
    contract = json.loads(contract_path.read_text())
    assert contract['inference'].startswith('No p-values')
    tables = BASE / 'tables/descriptive'; tables.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5))
    plt.rcParams.update({'font.size': 10, 'svg.fonttype': 'none'})
    record = {'status': 'running', 'started_utc': datetime.now(timezone.utc).isoformat(),
              'contract_sha256': hashlib.sha256(contract_path.read_bytes()).hexdigest(),
              'source_files': [], 'count_datasets': [], 'inference': 'None; deposited aliases are not verified independent animals or genotypes'}
    for ax, accession in zip(axes, ['GSE154966', 'GSE273123']):
        entry = next(f for f in inventory if f['accession'] == accession and 'count' in f['path'].lower())
        path = ROOT / entry['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256']
        record['source_files'].append(entry)
        data = pd.read_csv(path, sep='\t', comment='#', index_col=0)
        if accession == 'GSE154966':
            x = data.iloc[:, 5:]
            assert list(data.columns[:5]) == ['Chr', 'Start', 'End', 'Strand', 'Length']
            labels = []
            for col in x:
                match = re.match(r'(.+)_Tigit_(Neg|Pos)_S\d+_', col)
                assert match, col
                labels.append({'sample_id': col, 'source_alias': match[1], 'sort_gate': match[2]})
        else:
            x = data
            labels = []
            for col in x:
                match = re.fullmatch(r'(.+)_CD44([+-])', col)
                assert match, col
                labels.append({'sample_id': col, 'source_alias': match[1], 'sort_gate': 'Pos' if match[2] == '+' else 'Neg'})
        a = x.to_numpy()
        assert not x.index.duplicated().any() and np.isfinite(a).all() and (a >= 0).all() and (a == np.floor(a)).all() and (a.sum(axis=0) > 0).all()
        meta = pd.DataFrame(labels)
        assert (meta.groupby(['source_alias', 'sort_gate']).size() == 1).all()
        assert (meta.groupby('source_alias').size() == 2).all()
        cpm = a / a.sum(axis=0) * 1e6
        keep = (cpm >= 1).sum(axis=1) >= (a.shape[1] + 1) // 2
        logcpm = np.log2(cpm[keep] + 1)
        order = np.argsort(-logcpm.var(axis=1, ddof=1), kind='stable')[:2000]
        features = x.index[keep][order]
        centered = logcpm[order].T - logcpm[order].T.mean(axis=0)
        u, s, vt = np.linalg.svd(centered, full_matrices=False)
        scores = u * s
        var = s**2 / (s**2).sum()
        meta['PC1'], meta['PC2'] = scores[:, 0], scores[:, 1]
        meta['raw_library_count'], meta['detected_features'] = a.sum(axis=0), (a > 0).sum(axis=0)
        meta.to_csv(tables / f'{accession}_source_PCA_QC.tsv', sep='\t', index=False)
        pd.DataFrame({'feature_id': features}).to_csv(tables / f'{accession}_PCA_features.tsv', sep='\t', index=False)
        pd.DataFrame({'component': np.arange(1, len(var) + 1), 'variance_fraction': var}).to_csv(tables / f'{accession}_PCA_variance.tsv', sep='\t', index=False)
        for _, block in meta.groupby('source_alias'):
            ax.plot(block.PC1, block.PC2, color='#CBCFD2', lw=.8, zorder=1)
        for gate, color, marker in [('Neg', '#64717C', 'o'), ('Pos', '#AA5B8A', '^')]:
            g = meta[meta.sort_gate == gate]
            ax.scatter(g.PC1, g.PC2, c=color, marker=marker, label=gate, s=45, zorder=2)
            for row in g.itertuples():
                ax.annotate(row.source_alias, (row.PC1, row.PC2), xytext=(3, 4), textcoords='offset points', fontsize=7)
        ax.set(xlabel=f'PC1 ({var[0]:.1%})', ylabel=f'PC2 ({var[1]:.1%})', title=f"{accession} · {'TIGIT ATAC' if accession == 'GSE154966' else 'CD44 RNA'}")
        ax.spines[['top', 'right']].set_visible(False)
        ax.legend(frameon=False, fontsize=9)
        record['count_datasets'].append({'accession': accession, 'features': len(x), 'columns': x.shape[1], 'PCA_features': len(features), 'source_blocks': meta.source_alias.nunique()})
    fig.suptitle('Descriptive sample structure from deposited count matrices', fontsize=14)
    fig.text(.5, .02, 'Lines connect deposited aliases, not certified independent animals. R26/OG genotype mapping unresolved.\nCPM + log2 transform; no batch correction or significance testing.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .95))
    for ext in ('png', 'svg'):
        fig.savefig(BASE / f'figures/a1_deposited_source_PCA.{ext}', dpi=220)
    plt.close(fig)
    peaks = []
    for entry in (f for f in inventory if f['accession'] == 'GSE141635'):
        path = ROOT / entry['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256']
        record['source_files'].append(entry)
        with gzip.open(path, 'rt') as f:
            headers = [line.strip() for line in f if line.startswith('#')]
        data = pd.read_csv(path, sep='\t', comment='#', header=None)
        assert data.shape[1] == 6 and data.iloc[:, 0].str.startswith('chr').all()
        widths = data.iloc[:, 2] - data.iloc[:, 1]
        assert (data.iloc[:, 1] >= 0).all() and (widths > 0).all()
        cmd = next(h for h in headers if 'cmd' in h)
        peaks.append({'file': path.name, 'called_intervals': len(data), 'median_width_bp': float(widths.median()),
                      'q25_width_bp': float(widths.quantile(.25)), 'q75_width_bp': float(widths.quantile(.75)),
                      'caller_command': cmd, 'biological_comparison': 'HELD: incompatible peak-caller settings'})
    pd.DataFrame(peaks).to_csv(tables / 'H3K4me3_technical_geometry.tsv', sep='\t', index=False)
    record.update(status='completed', finished_utc=datetime.now(timezone.utc).isoformat(),
                  histone_result='Different size/minDist/caller parameters preclude biological interpretation of interval count/overlap',
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (BASE / 'reports/descriptive_run.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record['count_datasets'], indent=2))


if __name__ == '__main__':
    main()
