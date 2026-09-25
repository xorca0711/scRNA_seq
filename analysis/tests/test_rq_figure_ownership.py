"""Exercise the render boundary without importing the scientific plotting stack."""
import argparse
import ast
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time
import statistics
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/16_research_question_figures.py'


def entrypoints(root):
    # Execute the actual boundary functions; scientific renderers are substituted
    # at their call boundary so these regression tests run in lightweight CI.
    tree = ast.parse(SCRIPT.read_text(encoding='utf-8'))
    names = {'main', 'cached', 'sha256', 'ensure_promoter_score'}
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    nodes += [n for n in tree.body if isinstance(n, ast.Assign) and
              any(isinstance(t, ast.Name) and t.id == 'ACTIVE_FIGURES' for t in n.targets)]
    env = dict(Path=Path, argparse=argparse, datetime=datetime, timezone=timezone,
               hashlib=hashlib, json=json, time=time, REPO=root, REBUILD=False,
               CACHE=root/'analysis/figures/rq/processed', OUT=root/'analysis/figures/rq',
               RECORD={}, __file__=str(SCRIPT), __doc__='test renderer',
               vs=SimpleNamespace(apply=Mock()), plt=Mock(), log=Mock())
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(SCRIPT), 'exec'), env)
    return env


class FigureOwnershipTests(unittest.TestCase):
    def test_curated_numeric_captions_match_their_saved_tables(self):
        folder = SCRIPT.parents[1] / 'figures/rq'
        gallery = (folder/'README.md').read_text(encoding='utf-8')
        def rows(name):
            with (folder/name).open(encoding='utf-8', newline='') as handle:
                return list(csv.DictReader(handle))
        def section(key, next_key):
            return gallery.split(f'<a id="{key}"></a>', 1)[1].split(f'<a id="{next_key}"></a>', 1)[0]
        a1 = section('a1', 'a2')
        values = rows('rq_a1_detection_at_budget.csv')
        for field in ['rna_detection_pct', 'atac_detection_pct']:
            means = [statistics.mean(float(x[field]) for x in values if x['group'] == group)
                     for group in ['PBS reference', 'SeV transitional']]
            self.assertIn(f'{means[0]:.0f}%', a1)
            self.assertIn(f'{means[1]:.0f}%', a1)
        a3 = section('a3', 'a4')
        for row in rows('rq_a3_icap_by_day.csv'):
            if float(row['day']) in [0, 25, 366]:
                self.assertIn(f"{float(row['median_icap_pct']):.1f}%", a3)
        a4 = section('a4', 'a5')
        values = rows('rq_a4_detection.csv')
        for gene in ['Axin2', 'Il1r1']:
            fractions = [float(x['pct_detected']) for x in values if x['gene'] == gene]
            self.assertIn(f'{min(fractions):.1f}–{max(fractions):.1f}%', a4)
        a5 = section('a5', 'a6')
        for row in rows('rq_a5_wells.csv'):
            self.assertIn(f"{int(row['n_nuclei']):,}", a5)
            self.assertIn(f"{float(row['transitional_pct']):.2f}%", a5)

    def test_selected_a3_render_preserves_authored_text_and_previous_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            env = entrypoints(root)
            protected = {
                'RESEARCH_QUESTIONS.md': b'Authored A2 without an old marker.\n',
                'analysis/figures/rq/README.md': b'Curated captions.\n',
                'analysis/figures/rq/run_record.json': b'{"original": true}\n',
                'analysis/figures/rq/rq_a2_source_rank.png': b'current A2',
                'analysis/figures/rq/rq_a3_persistence.png': b'previous A3',
                'analysis/config/palette.json': b'{}',
            }
            for name, data in protected.items():
                p = root/name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_bytes(data)
            mouse = root/'mouse'
            for name in ['myeloid_focus/tables/myeloid_cell_metadata.csv',
                         'regeneration_focus/tables/icap_abundance_per_sample.csv']:
                p = mouse/name
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text('source\n1\n')
            cap = SimpleNamespace(isbacked=True, file=SimpleNamespace(close=Mock()))
            env.update(MOUSE=mouse, cached=Mock(return_value=cap), build_capillary=Mock())
            def render(_):
                (env['OUT']/'rq_a3_persistence.png').write_bytes(b'new A3')
                env['RECORD']['panels'] = {'A3': {'myeloid_cells': 9997}}
            env['figure_a3'] = Mock(side_effect=render)
            for name in ['figure_a1', 'figure_a4', 'figure_a5', 'figure_a2']:
                env[name] = Mock(side_effect=AssertionError('Unselected figure executed'))
            with patch.object(sys, 'argv', ['render', '--figures', 'A3', '--replot']), \
                 patch('importlib.metadata.version', return_value='test'):
                self.assertEqual(env['main'](), 0)
            for name, data in protected.items():
                self.assertEqual((root/name).read_bytes(), data, name)
            self.assertNotEqual(env['OUT'], root/'analysis/figures/rq')
            record = json.loads((env['OUT']/'run_record.json').read_text())
            self.assertEqual(record['selected_figures'], ['A3'])
            self.assertEqual(record['render_mode'], 'existing_embeddings')
            self.assertEqual(set(record['output_sha256']), {'rq_a3_persistence.png'})
            env['cached'].assert_called_once_with('capillary', env['build_capillary'], backed='r')
            cap.file.close.assert_called_once()

    def test_retired_a2_cannot_be_selected(self):
        with tempfile.TemporaryDirectory() as directory:
            env = entrypoints(Path(directory))
            with patch.object(sys, 'argv', ['render', '--figures', 'A2', '--replot']), \
                 patch('sys.stderr'), self.assertRaises(SystemExit) as raised:
                env['main']()
            self.assertEqual(raised.exception.code, 2)
            self.assertFalse((Path(directory)/'analysis').exists())

    def test_missing_or_invalid_cache_does_not_rebuild_or_delete(self):
        with tempfile.TemporaryDirectory() as directory:
            env = entrypoints(Path(directory))
            builder = Mock(side_effect=AssertionError('Unexpected rebuild'))
            reader = Mock(side_effect=ValueError('Corrupt cache'))
            fake_ad = SimpleNamespace(read_h5ad=reader)
            p = env['CACHE']/'capillary.h5ad'
            with patch.dict(sys.modules, {'anndata': fake_ad}):
                with self.assertRaises(FileNotFoundError):
                    env['cached']('capillary', builder)
                p.parent.mkdir(parents=True)
                p.write_bytes(b'invalid but preserved cache')
                with self.assertRaises(ValueError):
                    env['cached']('capillary', builder)
                self.assertEqual(p.read_bytes(), b'invalid but preserved cache')
            builder.assert_not_called()

    def test_missing_promoter_data_cannot_rewrite_a_replot_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            env = entrypoints(Path(directory))
            adata = SimpleNamespace(obsm={}, obs={}, write_h5ad=Mock())
            with self.assertRaisesRegex(ValueError, 'never rewrites'):
                env['ensure_promoter_score'](adata)
            adata.write_h5ad.assert_not_called()


if __name__ == '__main__':
    unittest.main()
