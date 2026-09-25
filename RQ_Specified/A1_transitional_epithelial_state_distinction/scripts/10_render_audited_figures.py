"""Render existing results without rerunning fits or changing scientific tables."""
import hashlib
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

BASE = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inputs = {p.relative_to(BASE).as_posix(): sha(p) for p in (BASE / 'tables').rglob('*.tsv')}
    run = json.loads((BASE / 'reports/ire1_run.json').read_text())
    assert run['status'] == 'completed'
    assert all(sha(BASE / path) == digest for path, digest in run['output_sha256'].items())
    script = Path(__file__).with_name('07_ire1_epithelial_analysis.py')
    spec = importlib.util.spec_from_file_location('ire1', script)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    mod.draw_figures(BASE / 'tables/ire1', json.loads((BASE / 'config/ire1_kira8.json').read_text()),
                     json.loads((BASE / 'metadata/ensembl_marker_lookup.json').read_text()))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, accession, label in zip(axes, ['GSE154966', 'GSE273123'], ['TIGIT ATAC', 'CD44 RNA']):
        meta = pd.read_csv(BASE / f'tables/descriptive/{accession}_source_PCA_QC.tsv', sep='\t')
        var = pd.read_csv(BASE / f'tables/descriptive/{accession}_PCA_variance.tsv', sep='\t').variance_fraction
        for _, block in meta.groupby('source_alias'):
            ax.plot(block.PC1, block.PC2, color='#CBCFD2', lw=.8, zorder=1)
        for gate, color, marker in [('Neg', '#64717C', 'o'), ('Pos', '#AA5B8A', '^')]:
            g = meta[meta.sort_gate == gate]
            ax.scatter(g.PC1, g.PC2, c=color, marker=marker, label=gate, s=45, zorder=2)
        ax.set(xlabel=f'PC1 ({var.iloc[0]:.1%})', ylabel=f'PC2 ({var.iloc[1]:.1%})', title=f'{accession} · {label}')
        ax.legend(frameon=False, fontsize=9)
    fig.suptitle('Descriptive sample structure from deposited counts', fontsize=14)
    fig.text(.5, .025, 'Lines join deposited aliases; independence and genotype mapping remain unresolved.\nNo significance tests. Exact source labels and coordinates are in the accompanying tables.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .10, 1, .93))
    for ext in ('png', 'svg'):
        fig.savefig(BASE / f'figures/a1_deposited_source_PCA.{ext}', dpi=220)
    plt.close(fig)
    # Rewriting the marker summary in the plotting helper must preserve its values/bytes.
    assert all(sha(BASE / path) == digest for path, digest in inputs.items())
    record = {'utc': datetime.now(timezone.utc).isoformat(), 'status': 'completed',
              'scope': 'Presentation-only rerender; remove overlapping PCA point labels, keep identities in source tables',
              'input_table_sha256': inputs, 'script_sha256': {p.name: sha(p) for p in (Path(__file__), script)},
              'figures': {p.name: sha(p) for p in (BASE / 'figures').glob('a1_*.*')}}
    (BASE / 'reports/figure_render.json').write_text(json.dumps(record, indent=2) + '\n')
    print('Figures rendered; all scientific-table hashes unchanged.')


if __name__ == '__main__':
    main()
