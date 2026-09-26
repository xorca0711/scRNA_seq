"""Cache a small public source, with byte limit and provenance; no expression download."""
from __future__ import annotations

import argparse
import hashlib
import gzip
import json
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("name")
    parser.add_argument("url")
    parser.add_argument("--max-bytes", type=int, default=5_000_000)
    parser.add_argument("--accept", default="*/*")
    parser.add_argument("--gzip-first-line", action="store_true",
                        help="Save only a decompressed header; never the complete matrix")
    args = parser.parse_args()
    if Path(args.name).name != args.name:
        parser.error("name must be a basename")
    folder = BASE / "cache/sources"
    folder.mkdir(parents=True, exist_ok=True)
    target = folder / args.name
    if target.exists():
        print(json.dumps({"status": "cached", "file": str(target.relative_to(BASE))}))
        return
    request = urllib.request.Request(args.url, headers={"User-Agent": "A0-scRNA-research-audit/1.0", "Accept": args.accept})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            if args.gzip_first_line:
                with gzip.GzipFile(fileobj=response) as decoded:
                    payload = decoded.readline(args.max_bytes + 1)
            else:
                payload = response.read(args.max_bytes + 1)
            resolved = response.url
    except urllib.error.URLError as error:
        raise SystemExit(f"Metadata retrieval failed: {error.reason}") from None
    if len(payload) > args.max_bytes:
        raise ValueError("Source exceeds the metadata byte budget")
    if args.gzip_first_line and not payload.endswith(b"\n"):
        raise ValueError("No complete header within metadata budget")
    target.write_bytes(payload)
    record = {"requested_url": args.url, "resolved_url": resolved,
              "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
              "bytes": len(payload), "sha256": hashlib.sha256(payload).hexdigest(),
              "retrieval_mode": "decompressed_first_line_only" if args.gzip_first_line else "complete_file",
              "hash_scope": "saved_bytes"}
    target.with_suffix(target.suffix + ".provenance.json").write_text(
        json.dumps(record, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "downloaded", "file": str(target.relative_to(BASE)), **record}))


if __name__ == "__main__":
    main()
