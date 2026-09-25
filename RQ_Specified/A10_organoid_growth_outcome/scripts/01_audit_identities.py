"""Stage 1: identity and join audit for the alveolosphere screen.

Reads three small metadata tables only: perturbation design, imaging outcomes and
species-assignment quality control. Never reads the count table. Decides whether the
design can carry a model, or whether results are restricted to description.

Refuses to overwrite its outputs.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import platform
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parents[1]
CACHE = HERE / 'cache'
OUT = HERE / 'tables'
CONFIG = HERE / 'config/a10_outcome_contract.json'
OUTPUTS = ['join_summary.tsv', 'well_join.tsv', 'preparation_units.tsv',
           'target_replication.tsv', 'species_assignment.tsv', 'audit_run.json']
DESIGN = 'GSE307112_plate_design.csv.gz'
IMAGING = 'GSE307112_imaging_outputs.csv.gz'
SPECIES = 'GSE307112_xenome_stats.csv.gz'


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 24), b''):
            digest.update(block)
    return digest.hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT)).replace('\\', '/')


def pad_well(value: str) -> str:
    """Normalize a well label to one letter plus two digits. Design is unpadded; the other tables are padded."""
    m = re.fullmatch(r'([A-Za-z])0*(\d{1,2})', str(value).strip())
    return f'{m.group(1).upper()}{int(m.group(2)):02d}' if m else str(value).strip().upper()


def norm_replicate(plate: str, label: str) -> str:
    """Normalize a replicate label. Imaging uses 'plate1-1'; species QC uses '1-1'."""
    text = str(label).strip()
    if text.startswith('plate'):
        text = text[len('plate'):]
    return f'{plate}-rep{text.split("-")[-1]}'


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text('empty\n', encoding='utf-8')
        return
    header = list(rows[0].keys())
    lines = ['\t'.join(header)]
    for r in rows:
        lines.append('\t'.join('none' if r[h] is None or r[h] == '' else str(r[h]) for h in header))
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main() -> None:
    import pandas as pd

    OUT.mkdir(exist_ok=True)
    existing = [n for n in OUTPUTS if (OUT / n).exists()]
    if existing:
        raise SystemExit(f'Refusing to overwrite: {existing}')
    config = json.loads(CONFIG.read_text(encoding='utf-8'))
    for name, meta in config['cached_inputs']['files'].items():
        path = CACHE / name
        if not path.exists():
            raise SystemExit(f'Missing cached input {name}')
        if path.stat().st_size != meta['bytes'] or sha256(path) != meta['sha256']:
            raise SystemExit(f'{name} differs from the audited file; refusing to proceed')
    if (CACHE / 'GSE307112_gene_counts.xlsx').exists():
        raise SystemExit('The count table is present; this stage must not run with counts available')

    design = pd.read_csv(CACHE / DESIGN)
    imaging = pd.read_csv(CACHE / IMAGING)
    species = pd.read_csv(CACHE / SPECIES)

    design['well_key'] = design.well.map(pad_well)
    imaging['well_key'] = imaging.well.map(pad_well)
    species['well_key'] = species.well.map(pad_well)
    imaging['unit'] = [norm_replicate(p, r) for p, r in zip(imaging.plate, imaging.plate_replicate)]
    species['unit'] = [norm_replicate(p, r) for p, r in zip(species.plate, species.plate_rep)]

    findings = []

    def note(check, result, detail):
        findings.append({'check': check, 'result': result, 'detail': detail})

    note('design key format', 'normalized',
         f'design wells are unpadded (example {design.well.iloc[0]}), the other tables padded '
         f'(example {species.well.iloc[0]}); both normalized to one letter plus two digits')
    note('replicate label format', 'normalized',
         f'imaging uses {imaging.plate_replicate.iloc[0]!r} and species QC {species.plate_rep.iloc[0]!r}; '
         'both normalized to plate plus replicate index')
    note('design has no replicate column', 'by design',
         'the design table is the plate layout and applies to every replicate of that plate; '
         'it therefore joins on plate and well only')

    dup_design = int(design.duplicated(['plate', 'well_key']).sum())
    note('duplicate design keys', 'pass' if dup_design == 0 else 'fail', f'{dup_design} duplicated plate and well rows')
    dup_species = int(species.duplicated(['unit', 'well_key']).sum())
    note('duplicate library keys', 'pass' if dup_species == 0 else 'fail',
         f'{dup_species} duplicated unit and well rows in the species table')
    dup_imaging = int(imaging.duplicated(['unit', 'well_key', 'day']).sum())
    note('duplicate imaging keys', 'pass' if dup_imaging == 0 else 'fail',
         f'{dup_imaging} duplicated unit, well and day rows')

    d_keys = set(zip(design.plate, design.well_key))
    s_keys = set(zip(species.plate, species.well_key))
    i_keys = set(zip(imaging.plate, imaging.well_key))
    note('wells in species QC absent from design', 'report', f'{len(s_keys - d_keys)} plate and well combinations')
    note('wells in imaging absent from design', 'report', f'{len(i_keys - d_keys)} plate and well combinations')
    note('wells in design absent from species QC', 'report', f'{len(d_keys - s_keys)} plate and well combinations')

    units_imaging = set(imaging.unit)
    units_species = set(species.unit)
    note('replicate units with imaging but no RNA', 'report',
         ', '.join(sorted(units_imaging - units_species)) or 'none')
    note('replicate units with RNA but no imaging', 'report',
         ', '.join(sorted(units_species - units_imaging)) or 'none')

    day_counts = imaging.day.value_counts().to_dict()
    both_days = imaging.groupby(['unit', 'well_key']).day.nunique()
    note('imaging days per well', 'report',
         f'{day_counts}; wells with both days {int((both_days == 2).sum())}, with one day {int((both_days == 1).sum())}')

    rows = []
    for unit in sorted(units_imaging | units_species):
        iu = imaging[imaging.unit == unit]
        su = species[species.unit == unit]
        paired = iu.groupby('well_key').day.nunique()
        wells_both = set(paired[paired == 2].index)
        rows.append({'unit': unit, 'plate': unit.split('-rep')[0], 'replicate_index': unit.split('-rep')[-1],
                     'rna_libraries': int(len(su)), 'imaging_rows': int(len(iu)),
                     'wells_with_both_days': int(len(wells_both)),
                     'wells_with_rna_and_both_days': int(len(set(su.well_key) & wells_both)),
                     'distinct_targets': int(su.crispr_target.nunique()) if len(su) else 0})
    write_tsv(OUT / 'preparation_units.tsv', rows)
    complete = [r for r in rows if r['wells_with_rna_and_both_days'] > 0]
    note('candidate preparation units', 'report',
         f'{len(rows)} plate and replicate combinations; {len(complete)} have RNA wells with both imaging days')
    note('biological independence of replicate indices', 'UNRESOLVED',
         'the deposited metadata give plate and replicate labels only. Nothing in these three files states '
         'whether a replicate index is a separate type 2 isolation and fibroblast lot or a repeat of one '
         'preparation. Resolving it needs the series sample metadata or the authors. Until then, holdouts by '
         'replicate index cannot be called independent biological preparations.')

    joined = species.merge(design.drop(columns=['well', 'plate.1'], errors='ignore'),
                           on=['plate', 'well_key'], how='left', validate='many_to_one')
    wide = imaging.pivot_table(index=['unit', 'well_key'], columns='day',
                               values=['organoids_count', 'organoids_area_mean', 'organoids_area_prop'])
    wide.columns = [f'{a}_{b}' for a, b in wide.columns]
    wide = wide.reset_index()
    full = joined.merge(wide, on=['unit', 'well_key'], how='left')
    area_cols = [c for c in full.columns if c.startswith('organoids_area_mean')]
    usable = full.dropna(subset=area_cols) if len(area_cols) == 2 else full.iloc[0:0]
    note('libraries joined to design', 'report',
         f'{int(full.gene_name.notna().sum())} of {len(full)} libraries matched a design row by plate and well')
    note('libraries with both imaging days', 'report',
         f'{len(usable)} of {len(full)} libraries have day-7 and day-14 mean area')
    write_tsv(OUT / 'well_join.tsv', full[['library name', 'plate', 'unit', 'well_key', 'crispr_target', 'gene_name',
                                           'gene_id'] + sorted(c for c in wide.columns if c not in ('unit', 'well_key'))
                                          ].to_dict('records'))

    vc = species.crispr_target.value_counts()
    trows = []
    for target, n in vc.items():
        gene = design.loc[design.gene_name.astype(str).str.upper() == str(target).upper(), 'gene_name']
        trows.append({'crispr_target': target, 'libraries': int(n),
                      'design_gene_name': gene.iloc[0] if len(gene) else None,
                      'in_design': bool(len(gene)),
                      'units_present': int(species[species.crispr_target == target].unit.nunique()),
                      'flag': 'extreme_replication' if n > 12 else ('below_4_libraries' if n < 4 else '')})
    write_tsv(OUT / 'target_replication.tsv', trows)
    extreme = [r for r in trows if r['flag'] == 'extreme_replication']
    note('targets with extreme replication', 'report',
         '; '.join(f"{r['crispr_target']} in {r['libraries']} libraries across {r['units_present']} units" for r in extreme)
         or 'none')
    note('targets absent from the design table', 'report',
         ', '.join(r['crispr_target'] for r in trows if not r['in_design']) or 'none')
    note('control identity', 'UNRESOLVED',
         'which targets the authors treated as controls is not established from these files, and the source '
         'paper is unread by the owner, so it is not taken from there either.')

    read_cols = [c for c in species.columns if c.startswith('xenome_numReads')]
    total = species['xenome_numReadsInput']
    srows = []
    for col in read_cols:
        if col == 'xenome_numReadsInput':
            continue
        frac = species[col] / total
        srows.append({'field': col, 'median_fraction_of_input': round(float(frac.median()), 4),
                      'min': round(float(frac.min()), 4), 'max': round(float(frac.max()), 4)})
    write_tsv(OUT / 'species_assignment.tsv', srows)
    amb = [r for r in srows if r['field'].endswith(('Ambiguous', 'Both', 'Neither'))]
    note('unassigned read fraction', 'report',
         '; '.join(f"{r['field'].replace('xenome_numReads','')} median {r['median_fraction_of_input']}" for r in amb))
    note('which assigned genome is which species', 'UNRESOLVED',
         'the two assignment fields are named by role, not by species, and these files do not map either to '
         'mouse epithelium or human fibroblast. A stage 3 check settles it without guessing: the perturbed '
         'genes are mouse, so knocking one out should reduce that gene in whichever compartment is mouse. '
         'Until settled, the epithelial and niche sides must not be labelled.')

    write_tsv(OUT / 'join_summary.tsv', findings)

    unresolved = [f['check'] for f in findings if f['result'] == 'UNRESOLVED']
    failed = [f['check'] for f in findings if f['result'] == 'fail']
    verdict = ('stop_joins_failed' if failed else
               'descriptive_only_unit_unresolved' if 'biological independence of replicate indices' in unresolved else
               'model_permitted')
    record = {
        'stage': 1, 'purpose': 'identity and join audit; metadata only, no counts read',
        'completed_utc': dt.datetime.now(dt.timezone.utc).isoformat(),
        'script': rel(Path(__file__).resolve()), 'script_sha256': sha256(Path(__file__).resolve()),
        'config_sha256': sha256(CONFIG),
        'python': sys.version.split()[0], 'platform': platform.platform(), 'pandas': pd.__version__,
        'inputs': {rel(CACHE / n): {'bytes': (CACHE / n).stat().st_size, 'sha256': sha256(CACHE / n)}
                   for n in [DESIGN, IMAGING, SPECIES]},
        'counts_read': False,
        'libraries': int(len(species)), 'targets': int(species.crispr_target.nunique()),
        'candidate_units': len(rows), 'libraries_with_both_imaging_days': int(len(usable)),
        'unresolved': unresolved, 'failed_checks': failed, 'verdict': verdict,
        'outputs': {n: sha256(OUT / n) for n in OUTPUTS if n != 'audit_run.json'},
    }
    (OUT / 'audit_run.json').write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: record[k] for k in ['libraries', 'targets', 'candidate_units',
                                             'libraries_with_both_imaging_days', 'unresolved', 'verdict']}, indent=2))


if __name__ == '__main__':
    main()
