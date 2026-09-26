"""Stream one public research input into A0 cache, retaining source integrity."""
import argparse
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    parser.add_argument("url")
    parser.add_argument("--max-bytes", type=int, default=1_000_000_000)
    parser.add_argument("--accept", default="*/*")
    args = parser.parse_args()
    assert Path(args.name).name == args.name
    folder = BASE / "cache/expression"
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / args.name
    provenance = target.with_name(target.name + ".provenance.json")
    if target.exists():
        assert provenance.exists(), "Existing input lacks provenance"
        expected = json.loads(provenance.read_text())
        digest = hashlib.file_digest(target.open("rb"), "sha256").hexdigest()
        assert digest == expected["sha256"] and target.stat().st_size == expected["bytes"]
        print(json.dumps({"status": "cached_verified", "file": args.name}), flush=True)
        return
    part = target.with_name(target.name + ".part")
    digest = hashlib.sha256()
    size = 0
    request = urllib.request.Request(args.url, headers={"Accept": args.accept, "User-Agent": "A0-exploratory-scRNA/1.0"})
    with urllib.request.urlopen(request, timeout=60) as response, part.open("wb") as stream:
        for chunk in iter(lambda: response.read(1024 * 1024), b""):
            size += len(chunk)
            if size > args.max_bytes:
                raise ValueError("Input exceeds configured byte limit")
            stream.write(chunk)
            digest.update(chunk)
        resolved = response.url
    record = {"requested_url": args.url, "resolved_url": resolved, "bytes": size,
              "sha256": digest.hexdigest(), "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
              "request_accept": args.accept, "retrieval_mode": "complete_file", "hash_scope": "saved_bytes"}
    part.replace(target)
    provenance.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "downloaded", "file": args.name, "bytes": size, "sha256": digest.hexdigest()}), flush=True)


if __name__ == "__main__":
    main()
