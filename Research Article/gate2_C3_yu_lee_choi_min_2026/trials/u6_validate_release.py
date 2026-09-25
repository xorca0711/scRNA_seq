"""Validate completed stages, local links and the exact paper-local release files."""
from pathlib import Path
from datetime import datetime, timezone
import json
import hashlib
import re
import subprocess
import sys
from urllib.parse import unquote

PAPER = Path(__file__).resolve().parents[1]
ROOT = PAPER.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.provenance import sha256_file, code_identity, write_json_atomic
from analysis.lib.repository_paths import recorded_file


def main():
    import pandas as pd
    out = PAPER / 'trials/u6_completion'
    validations = ['u5_ipf_pathways', 'u5_ipf_liana', 'u5_ipf_compatibility', 'u5_liana_robustness',
                   'u5_ligand_targets', 'u4_mouse_niche', 'u6_specificity', 'u6_isr_extension',
                   'u6_ipf_specificity', 'u5_spatial_context', 'u5_human_niche',
                   'u5_human_sources', 'u5_human_ligand_targets', 'u6_human_specificity']
    checked = []
    for name in validations:
        path = PAPER / 'trials' / name / 'validation.json'
        result = json.loads(path.read_text())
        if name == 'u5_ipf_pathways':
            assert result['checks_passed'] == 20 and len(result['cohorts']) == 2, result
        else:
            assert result['status'] == 'passed', (name, result)
        checked.append(dict(stage=name, validation_sha256=sha256_file(path)))
    pipeline = json.loads((PAPER / 'trials/u5_human_niche/paired_pipeline_run_record.json').read_text())
    assert pipeline['status'] == 'completed_computations_require_reports'
    spec = json.loads((PAPER / 'trials/u5_human_niche/paired_pipeline_specification.json').read_text())
    for name, expected in spec['scripts'].items():
        relative = (PAPER / 'trials' / name).relative_to(ROOT).as_posix()
        assert sha256_file(recorded_file(ROOT, relative, expected)) == expected, f'Historical frozen script missing: {name}'
    annotation = json.loads((PAPER / 'trials/u5_human_full/annotation_review_summary.json').read_text())
    assert annotation['release_decision_status'] == 'reviewed_for_reference_compatible_broad_niche_interpretation'
    coverage = json.loads((out / 'pathway_coverage_validation.json').read_text())
    assert coverage['status'] == 'passed_identical_eligibility'
    assert json.loads((out / 'pathway_correlation_validation.json').read_text())['status'] == 'passed'
    listing = subprocess.run(['git', 'ls-files', '-c', '-o', '--exclude-standard', '--', PAPER.relative_to(ROOT).as_posix()],
                             cwd=ROOT, check=True, text=True, capture_output=True).stdout.splitlines()
    paths = sorted({ROOT / s for s in listing if (ROOT / s).is_file()})
    forbidden = [p for p in paths if '/cache/' in p.as_posix() or p.suffix.lower() in ['.h5', '.h5ad', '.npz', '.pdf', '.rds', '.gz']]
    assert not forbidden, forbidden
    broken = []
    for path in paths:
        if path.suffix != '.md':
            continue
        for target in re.findall(r'\]\(([^)]+)\)', path.read_text(encoding='utf-8')):
            target = target.strip('<>').split('#')[0]
            if not target or re.match(r'^[a-zA-Z][a-zA-Z0-9+.-]*:', target):
                continue
            actual = (path.parent / unquote(target)).resolve()
            generated_here = {out / 'release_validation.json', out / 'release_manifest.csv'}
            if not actual.exists() and actual not in generated_here:
                broken.append((str(path.relative_to(PAPER)), target))
    assert not broken, broken
    figures = sorted((PAPER / 'figures').glob('*.png'))
    reviews = json.loads((out / 'figure_review.json').read_text())
    legacy_png_only = []
    for path in figures:
        assert path.stat().st_size > 10000
        if not path.with_suffix('.svg').exists():
            # Preserve the already committed initial raster figure unchanged.
            assert path.name == 'ipf_pathway_primary_and_sensitivity.png', path
            previous = subprocess.run(['git', 'show', 'HEAD:' + path.relative_to(ROOT).as_posix()], cwd=ROOT,
                                      capture_output=True, check=True).stdout
            assert hashlib.sha256(previous).hexdigest() == sha256_file(path)
            legacy_png_only.append(path.name)
        reviewed = reviews['figures'][path.name]
        assert reviewed['status'] == 'visually_reviewed'
        assert reviewed['png_sha256'] == sha256_file(path), f'Figure changed after review: {path.name}'
    manifests = []
    for path in paths:
        if path.name in ['release_manifest.csv', 'release_validation.json']:
            continue
        manifests.append(dict(path=path.relative_to(PAPER).as_posix(), bytes=path.stat().st_size, sha256=sha256_file(path)))
    frame = pd.DataFrame(manifests)
    assert frame.path.is_unique
    frame.to_csv(out / 'release_manifest.csv', index=False)
    result = dict(status='passed', checked_utc=datetime.now(timezone.utc).isoformat(), code=code_identity(ROOT, __file__),
                  stages=checked, stage_validation_count=len(checked), frozen_human_scripts_verified=len(spec['scripts']),
                  local_markdown_links='all_exist', forbidden_large_or_private_inputs_tracked=False,
                  visually_reviewed_figures=len(figures), legacy_png_only_figures=legacy_png_only,
                  release_files=len(frame), release_bytes=int(frame.bytes.sum()),
                  manifest_sha256=sha256_file(out / 'release_manifest.csv'))
    write_json_atomic(out / 'release_validation.json', result)
    print(json.dumps({k: v for k, v in result.items() if k not in ['code', 'stages']}))


if __name__ == '__main__':
    main()
