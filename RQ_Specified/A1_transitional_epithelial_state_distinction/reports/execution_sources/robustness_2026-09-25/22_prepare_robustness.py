"""Retrieve bounded, pinned annotation-definition evidence; never execute it."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import re
import urllib.request

BASE = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def main():
    report = BASE / 'reports/robustness_source_inventory.json'
    if report.exists():
        raise SystemExit('Refusing to overwrite source inventory')
    tree_path = BASE / 'cache/followup_sources/hpcs_code_tree.json'
    tree = json.loads(tree_path.read_text())
    wanted = ['tracing_depletion_analysis/01_Concatenate_cells.ipynb',
              'tracing_depletion_analysis/02_Analyze_control_DT_depletion.ipynb']
    entries = []
    for name in wanted:
        item = next(row for row in tree['tree'] if row['path'] == name)
        path = BASE / 'cache/robustness_sources' / name
        url = 'https://raw.githubusercontent.com/dbetel/HPCS_LUAD/' + tree['sha'] + '/' + name
        if item['size'] > 20 * 1024**2:
            # A held document is evidence of a resource limit, not an empty source.
            entries.append(dict(url=url, bytes=item['size'], git_blob_sha=item['sha'],
                                status='held_full_notebook_exceeds_20_MiB'))
            continue
        if not path.exists():
            with urllib.request.urlopen(url, timeout=45) as response:
                data = response.read(20 * 1024**2 + 1)
            assert len(data) == item['size']
            assert hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest() == item['sha']
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        data = path.read_bytes()
        assert len(data) == item['size']
        assert hashlib.sha1(b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest() == item['sha']
        entries.append(dict(path=path.relative_to(BASE).as_posix(), url=url,
                            bytes=len(data), sha256=sha(path), git_blob_sha=item['sha']))
    report.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(Path(__file__)), contract_sha256=sha(BASE/'config/robustness_batch.json'),
        tree_sha256=sha(tree_path), author_commit=tree['sha'], files=entries), indent=2)+'\n')
    print('Source inventory:', [(r.get('status', 'verified'), r['bytes']) for r in entries])


if __name__ == '__main__':
    main()
