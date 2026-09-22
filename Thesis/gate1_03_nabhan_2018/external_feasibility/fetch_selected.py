"""Fetch the frozen eight-sample GSE129605 contrast, bounded to 150 MB.

Uses only Python standard library. Existing files must match GEO byte size;
failed downloads retain a .part suffix. No expression analysis or cell labels.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import csv
import hashlib
import json
from urllib.request import urlopen

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
DEST = REPO / "raw_data" / "GSE129605"


def main():
    with (HERE / "GSE129605_download_manifest.csv").open() as f:
        manifest = list(csv.DictReader(f))
    assert len(manifest) == 24
    assert sum(int(r["expected_bytes"]) for r in manifest) < 150_000_000
    DEST.mkdir(parents=True, exist_ok=True)

    def fetch(row):
        dest = DEST / row["name"]
        expected = int(row["expected_bytes"])
        if not dest.exists() or dest.stat().st_size != expected:
            partial = dest.with_suffix(dest.suffix + ".part")
            size = 0
            with urlopen(row["url"], timeout=90) as response, partial.open("wb") as out:
                while chunk := response.read(1024 * 1024):
                    size += len(chunk)
                    if size > expected:
                        raise ValueError(f"Unexpected file size: {row['name']}")
                    out.write(chunk)
            if size != expected:
                raise ValueError(f"Incomplete download: {row['name']} {size}/{expected}")
            partial.replace(dest)
        return {**row, "observed_bytes": dest.stat().st_size,
                "sha256": hashlib.sha256(dest.read_bytes()).hexdigest()}

    with ThreadPoolExecutor(max_workers=3) as pool:
        outputs = list(pool.map(fetch, manifest))
    (HERE / "GSE129605_download_receipt.json").write_text(json.dumps(outputs, indent=2) + "\n")
    print(json.dumps({"files": len(outputs), "bytes": sum(r["observed_bytes"] for r in outputs), "destination": str(DEST)}))


if __name__ == "__main__":
    main()
