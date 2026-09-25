"""Relocation must preserve identities and reject missing/corrupted evidence."""
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from analysis.lib.repository_paths import MANIFEST, recorded_file, resolve_repo_path


class RepositoryPathTests(unittest.TestCase):
    def test_old_windows_and_posix_paths_resolve(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            expected = root / 'Research Article/paper/results.csv'
            for name in ['Thesis/paper/results.csv', r'Thesis\paper\results.csv',
                         'Research Article/paper/results.csv']:
                self.assertEqual(resolve_repo_path(root, name), expected.resolve())
            with self.assertRaises(ValueError):
                resolve_repo_path(root, '../outside.csv')

    def test_verified_archive_is_historical_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / 'Research Article/paper/run.py'
            source.parent.mkdir(parents=True)
            source.write_bytes(b'new code')
            archive = root / 'old.py.txt'
            archive.write_bytes(b'old code')
            expected = hashlib.sha256(archive.read_bytes()).hexdigest()
            manifest = root / MANIFEST
            manifest.parent.mkdir(parents=True)
            manifest.write_text(json.dumps({'files': [dict(new_path='Research Article/paper/run.py',
                original_sha256=expected, original_bytes='old.py.txt')]}))
            self.assertEqual(recorded_file(root, 'Thesis/paper/run.py', expected), archive)
            archive.write_bytes(b'corrupted archive')
            with self.assertRaises(ValueError):
                recorded_file(root, 'Thesis/paper/run.py', expected)
            with self.assertRaises(ValueError):
                recorded_file(root, 'Thesis/paper/unrelated.py', expected)

    def test_changed_scientific_output_is_not_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / 'result.csv'
            output.write_bytes(b'value\n1\n')
            expected = hashlib.sha256(output.read_bytes()).hexdigest()
            self.assertEqual(recorded_file(root, 'result.csv', expected), output)
            output.write_bytes(b'value\n2\n')
            with self.assertRaises(ValueError):
                recorded_file(root, 'result.csv', expected)
