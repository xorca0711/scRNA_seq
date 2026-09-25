"""Resolve pre-migration paths without changing historical run identities.

Use resolve_repo_path for live inputs. Use recorded_file only to verify a past
run: it may return the archived original script, never a substitute live input.
"""
from pathlib import Path
import hashlib
import json

MANIFEST = 'docs/migrations/2026-09-25-research-layout/manifest.json'


def resolve_repo_path(root: Path, relative: str) -> Path:
    root = Path(root).resolve()
    name = str(relative).replace('\\', '/')
    if name.startswith('Thesis/'):
        name = 'Research Article/' + name[len('Thesis/'):]
    path = (root / name).resolve()
    if not path.is_relative_to(root):
        raise ValueError(f'Path outside repository: {relative}')
    return path


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def recorded_file(root: Path, relative: str, expected_sha256: str) -> Path:
    """Return bytes matching a recorded identity or fail; not a rerun check."""
    path = resolve_repo_path(root, relative)
    if path.is_file() and digest(path) == expected_sha256:
        return path
    manifest = root / MANIFEST
    if manifest.exists():
        for entry in json.loads(manifest.read_text(encoding='utf-8'))['files']:
            if (resolve_repo_path(root, entry['new_path']) == path
                    and entry['original_sha256'] == expected_sha256
                    and entry.get('original_bytes')):
                archive = resolve_repo_path(root, entry['original_bytes'])
                if archive.is_file() and digest(archive) == expected_sha256:
                    return archive
    raise ValueError(f'Recorded file missing or changed: {relative}')
