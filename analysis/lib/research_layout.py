"""Lightweight checks for canonical RQ placement and migration integrity."""
from pathlib import Path
import hashlib
import json
import re
from analysis.lib.repository_paths import recorded_file


def check_research_layout(root: Path, result) -> None:
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    paper = root / 'Research Article/gate2_C3_yu_lee_choi_min_2026'
    package = root / 'analysis/figures/rq/il1b_context'
    questions = (root / 'RESEARCH_QUESTIONS.md').read_text(encoding='utf-8')
    ids = re.findall(r'^### (A\d+)\.', questions, flags=re.M)
    # A15 is proposed and pending the owner's retain or reject; if it is rejected,
    # drop the card and return this range to 15.
    result.equal(ids, [f'A{i}' for i in range(16)], 'canonical question sequence')
    result.require(not (paper / 'DERIVED_RESEARCH_QUESTIONS.md').exists(),
                   'duplicate current paper-local RQ register')
    roadmap = json.loads((root / 'Research Article/ROADMAP.json').read_text(encoding='utf-8'))
    for item in roadmap['papers']:
        folder = item.get('folder')
        if folder:
            result.require(not folder.startswith('Research Article/') and not folder.endswith('/'),
                           f'paper {item["order"]}: folder must be relative to Research Article')
            result.require((root / 'Research Article' / folder).is_dir(), f'missing paper folder: {folder}')
    manifest = json.loads((package / 'relocation_manifest.json').read_text())
    by_path = {x['new_path']: x for x in manifest['files']}
    for name in manifest['preserved_scientific_files']:
        result.equal(digest(root / name), by_path[name]['original_sha256'], f'migrated scientific bytes: {name}')
    prep = json.loads((package / 'preparation_run_record.json').read_text())
    result.equal(digest(package / '.history/layout_2026-09-25/u7_prepare_proposal_figures.py.txt'),
                 prep['code']['entrypoint_sha256'], 'original preparation code retained')
    render = json.loads((package / 'render_run_record.json').read_text())
    entrypoint = render['code']['entrypoint'].replace('\\', '/')
    result.equal(digest(recorded_file(root, entrypoint, render['code']['entrypoint_sha256'])),
                 render['code']['entrypoint_sha256'], 'historical rendering code retained')
    review = json.loads((package / 'visual_review.json').read_text())
    pngs = [x for x in render['outputs'] if x['path'].endswith('.png')]
    result.equal(len(pngs), 6, 'shared IL-1 figure count')
    for output in render['outputs']:
        result.equal(digest(root / output['path']), output['sha256'], f'rendered figure: {output["path"]}')
    for output in pngs:
        record = review['figures'].get(Path(output['path']).name, {})
        result.equal(record.get('status'), 'visually_reviewed', 'current PNG review status')
        result.equal(record.get('png_sha256'), output['sha256'], 'current PNG review hash')
