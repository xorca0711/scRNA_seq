"""Bounded public-file retrieval for the A0 continuation; immutable cache records."""
from __future__ import annotations
import argparse
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('name')
    p.add_argument('url')
    p.add_argument('--max-bytes', type=int, default=10_000_000)
    a = p.parse_args()
    if Path(a.name).name != a.name:
        p.error('name must be a basename')
    folder = BASE / 'cache/continuation_v1'
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / a.name
    record = target.with_name(target.name + '.provenance.json')
    if target.exists():
        r = json.loads(record.read_text())
        assert r['requested_url'] == a.url
        assert hashlib.sha256(target.read_bytes()).hexdigest() == r['sha256']
        print(json.dumps({'file': a.name, 'status': 'verified_cache', 'bytes': target.stat().st_size}))
        return
    req = urllib.request.Request(a.url, headers={'User-Agent': 'A0-public-research/1.0'})
    temporary = target.with_name(target.name + '.part')
    if temporary.exists():
        raise FileExistsError(temporary)
    digest, size = hashlib.sha256(), 0
    try:
        with urllib.request.urlopen(req, timeout=45) as response:
            if int(response.headers.get('Content-Length', 0)) > a.max_bytes:
                raise ValueError('Response exceeds declared byte limit')
            with temporary.open('xb') as output:
                while chunk := response.read(4 * 1024 * 1024):
                    size += len(chunk)
                    if size > a.max_bytes:
                        raise ValueError('Response exceeds declared byte limit')
                    output.write(chunk)
                    digest.update(chunk)
            r = {'requested_url': a.url, 'resolved_url': response.url,
                 'retrieved_at_utc': datetime.now(timezone.utc).isoformat(),
                 'bytes': size, 'sha256': digest.hexdigest(), 'cache_file': a.name}
        temporary.rename(target)
    except Exception:
        if temporary.exists():
            temporary.unlink()  # Only the temporary file created by this invocation.
        raise
    record.write_text(json.dumps(r, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({'file': a.name, 'status': 'downloaded', 'bytes': size}))

if __name__ == '__main__':
    main()
