"""Bounded, append-only primary-source acquisition; never execute author code."""
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import urllib.request

BASE = Path(__file__).resolve().parents[1]
CACHE = BASE / 'cache/regulatory_fate'
REPORT = BASE / 'reports/regulatory_fate_sources.json'
CAP = 20 * 1024**2
TOTAL_CAP = 150 * 1024**2


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('name')
    ap.add_argument('url')
    args = ap.parse_args()
    target = (CACHE / args.name).resolve()
    if not target.is_relative_to(CACHE.resolve()):
        raise SystemExit('Target must remain in the acquisition cache')
    record = json.loads(REPORT.read_text()) if REPORT.exists() else {'files': []}
    previous = [r for r in record['files'] if r['name'] == args.name and r['status'] == 'retrieved']
    if previous:
        row = previous[-1]
        assert row['url'] == args.url and sha(target.read_bytes()) == row['sha256']
        print(args.name, 'verified existing', row['bytes'])
        return
    used = sum(r.get('bytes', 0) for r in record['files'] if r['status'] == 'retrieved')
    row = dict(name=args.name, url=args.url, utc=datetime.now(timezone.utc).isoformat(),
               fetcher_sha256=sha(Path(__file__).read_bytes()),
               scope_sha256=sha((BASE / 'config/regulatory_fate_scope.md').read_bytes()))
    try:
        assert not target.exists(), 'Unrecorded file exists; inspect before reuse'
        req = urllib.request.Request(args.url, headers={'User-Agent': 'A1 scientific evidence audit'})
        with urllib.request.urlopen(req, timeout=60) as response:
            data = response.read(min(CAP, TOTAL_CAP - used) + 1)
        assert 0 < len(data) <= CAP and used + len(data) <= TOTAL_CAP, 'Acquisition ceiling exceeded'
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        row.update(status='retrieved', bytes=len(data), sha256=sha(data), path=target.relative_to(BASE).as_posix())
    except Exception as exc:
        row.update(status='retrieval_failed', error=str(exc))
    record['files'].append(row)
    REPORT.write_text(json.dumps(record, indent=2) + '\n')
    print(args.name, row['status'], row.get('bytes', row.get('error')))


if __name__ == '__main__':
    main()
