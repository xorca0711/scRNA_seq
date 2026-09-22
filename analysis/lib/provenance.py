"""Standard-library provenance helpers; no analysis-package imports required."""
from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def archive_existing_record(path: Path) -> Path | None:
    """Keep the exact previous bytes, addressed by content, before a new run."""
    path = Path(path)
    if not path.exists():
        return None
    data = path.read_bytes()
    digest = hashlib.sha256(data).hexdigest()
    destination = path.parent / ".history" / path.stem / f"{digest}.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists():
        if destination.read_bytes() != data:
            raise ValueError(f"Run-record archive mismatch: {destination}")
    else:
        with destination.open("xb") as handle:
            handle.write(data)
    return destination


def code_identity(root: Path, entrypoint: str | None = None) -> dict:
    """Record commit, tracked diff and entrypoint; explicitly expose dirty runs."""
    identity: dict = {}
    try:
        identity["commit"] = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=root, stderr=subprocess.DEVNULL,
            text=True).strip()
        diff = subprocess.check_output(
            ["git", "diff", "HEAD", "--binary"], cwd=root, stderr=subprocess.DEVNULL)
        identity["tracked_diff_sha256"] = hashlib.sha256(diff).hexdigest()
        identity["tracked_changes_present"] = bool(diff)
        identity["untracked_files_present"] = bool(subprocess.check_output(
            ["git", "ls-files", "--others", "--exclude-standard"], cwd=root,
            stderr=subprocess.DEVNULL).strip())
    except (OSError, subprocess.CalledProcessError):
        identity["commit"] = None
    if entrypoint and Path(entrypoint).is_file():
        path = Path(entrypoint).resolve()
        identity["entrypoint"] = str(path.relative_to(root)) if path.is_relative_to(root) else str(path)
        identity["entrypoint_sha256"] = sha256_file(path)
    return identity


def write_json_atomic(path: Path, value: dict, default=None) -> None:
    """Avoid leaving a partially written current run record after interruption."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, default=default) + "\n", encoding="utf-8")
    temporary.replace(path)
