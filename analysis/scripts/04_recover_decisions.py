#!/usr/bin/env python
"""
Rebuild logs/decisions.json from the pipeline's own run logs.

`run_scrna_analysis.py` now persists its decisions after every `record()` and
carries them forward when a run resumes at a later stage, so a fresh run has no
need for this.  It exists because this dataset was produced by several resumed
runs made *before* that persistence was added, and the decisions from the
earlier stages survive only in the run logs.

Each `record()` call emits exactly one line of the form

    [  606.0s]   DECISION <key> = <value>

so the logs are a faithful transcript.  Logs are applied in the order given,
later ones overriding earlier ones, which reproduces the order the stages
actually ran in.

Usage:
    python 04_recover_decisions.py --out ../logs/decisions.json LOG [LOG ...]
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

PATTERN = re.compile(r"^\[\s*[\d.]+s\]\s+DECISION (\S+) = (.*)$")


def parse(log_path: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    if not log_path.exists():
        return found
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        m = PATTERN.match(line)
        if m:
            found[m.group(1)] = m.group(2).strip()
    return found


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--out", required=True)
    ap.add_argument("--merge-existing", action="store_true",
                    help="keep keys already present in --out unless a log "
                         "supplies a newer value")
    args = ap.parse_args()

    out = Path(args.out)
    merged: dict[str, str] = {}
    if args.merge_existing and out.exists():
        merged.update(json.loads(out.read_text(encoding="utf-8")))
        print(f"starting from {len(merged)} existing key(s)")

    for lp in args.logs:
        got = parse(Path(lp))
        print(f"{Path(lp).name}: {len(got)} decision(s)")
        merged.update(got)

    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(merged, indent=2), encoding="utf-8")
    print(f"wrote {out} with {len(merged)} decision(s)")
    missing = [k for k in ["species", "n_hvg", "hvg_method", "n_pcs",
                           "batch_correction", "leiden_resolution", "n_clusters",
                           "normalization", "marker_test", "metadata_table"]
               if k not in merged]
    if missing:
        print(f"WARNING still missing: {missing}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
