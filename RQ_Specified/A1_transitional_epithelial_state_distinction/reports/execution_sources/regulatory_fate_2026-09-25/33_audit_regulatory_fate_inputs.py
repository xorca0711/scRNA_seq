"""Resolve primary identities and sequencing designs without changing old data."""
import csv
import json
import re
import xml.etree.ElementTree as ET
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
import openpyxl
from a1_regulatory_fate import BASE, describe, sha, write_tsv

OUT = BASE / 'tables/regulatory_fate'
CACHE = BASE / 'cache/regulatory_fate'


def main():
    report = BASE / 'reports/regulatory_fate_input_audit.json'
    assert not report.exists(), 'Refusing to overwrite audit'
    inputs = [CACHE / 'hpcs_resources.xlsx',
              BASE / 'tables/hpcs_source_composition/source_manifest.tsv',
              BASE / 'tables/hpcs_source_composition/source_state_counts.tsv',
              BASE / 'metadata/ENA_projects.json']
    with inputs[1].open() as f:
        manifest = list(csv.DictReader(f, delimiter='\t'))
    with inputs[2].open() as f:
        counts = list(csv.DictReader(f, delimiter='\t'))
    wb = openpyxl.load_workbook(inputs[0], read_only=True, data_only=True)
    mice, section = [], ''
    for rownum, row in enumerate(wb['Mice'].values, 1):
        if row[0] and str(row[0]).startswith('IGO'):
            section = str(row[0])
        elif row[0] and not row[1]:
            section = str(row[0])
        if len(row) >= 13 and row[1] and str(row[1]) != 'Tag':
            mice.append(dict(tag=str(row[1]), section=section, row=rownum,
                             sex=str(row[2]), harvest_week=row[4],
                             slc4a11=str(row[8]), rosa26=str(row[10]), hopx=str(row[11])))
    resolved = []
    for source in manifest:
        tag = source['source_label'].split('_')[0]
        hit = [m for m in mice if m['tag'] == tag and source['source_library'] in m['section']]
        assert len(hit) == 1, (tag, hit)
        m = hit[0]
        driver = 'Hopx' if m['hopx'] != '+/+' else 'Slc4a11'
        assert driver == source['driver']
        resolved.append(dict(source_label=source['source_label'], mouse_tag=tag,
            group=source['group'], library=source['source_library'], gsm=source['gsm'],
            driver=driver, chase_days=int(source['chase_days']),
            harvest_week_as_listed=m['harvest_week'], sex=m['sex'],
            slc4a11_genotype=m['slc4a11'], rosa26_genotype=m['rosa26'], hopx_genotype=m['hopx'],
            retained_cells=int(source['retained_cells']), source_sheet='Mice',
            source_range=f'A{m["row"]}:M{m["row"]}',
            evidence='exact animal tag plus library and driver agreement',
            current_reporter='not available for these traced cells'))
    assert len(resolved) == len({r['mouse_tag'] for r in resolved}) == 22
    assert sum(r['retained_cells'] for r in resolved) == 5333
    write_tsv(OUT / 'hpcs_mouse_crosswalk.tsv', resolved)
    libraries = {r['library'] for r in resolved}
    used = {r['mouse_tag'] for r in resolved}
    absent = [dict(mouse_tag=m['tag'], source_section=m['section'], source_row=m['row'],
                   status='listed in source experiment; no retained traced alias in the recovered 5333-cell subset; reason not inferred')
              for m in mice if any(lib in m['section'] for lib in libraries) and m['tag'] not in used]
    assert len(absent) == 4
    write_tsv(OUT / 'hpcs_listed_without_retained_trace.tsv', absent)
    groups = defaultdict(list)
    for r in counts:
        groups[(r['group'], r['state'])].append(float(r['fraction']))
    summaries = []
    for (group, state), vals in groups.items():
        members = [r for r in resolved if r['group'] == group]
        assert len(vals) == len(members)
        summaries.append(dict(group=group, state=state, **describe(vals),
            male=sum(r['sex'] == 'M' for r in members), female=sum(r['sex'] == 'F' for r in members),
            libraries=','.join(sorted({r['library'] for r in members})),
            estimand='equal named-mouse fraction; no chase-effect inference'))
    write_tsv(OUT / 'hpcs_named_mouse_summary.tsv', summaries)

    catalog = json.loads(inputs[3].read_text())
    sequence_rows = []
    for project in catalog['projects']:
        p = CACHE / (project['project'] + '_sra.xml')
        inputs.append(p)
        packages = ET.parse(p).findall('.//EXPERIMENT_PACKAGE')
        expected = {r['experiment_accession']: r for r in project['runs']}
        assert {r.find('EXPERIMENT').get('accession') for r in packages} == set(expected)
        for pkg in packages:
            experiment = pkg.find('EXPERIMENT').get('accession')
            sample = pkg.find('SAMPLE')
            attrs = {x.findtext('TAG'): x.findtext('VALUE') for x in sample.findall('.//SAMPLE_ATTRIBUTE')}
            title = sample.findtext('TITLE')
            sequence_rows.append(dict(project=project['project'], experiment=experiment,
                run=expected[experiment]['run_accession'], biosample=expected[experiment]['sample_accession'],
                sample_title=title, antibody=attrs.get('antibody', ''),
                biological_replicate=attrs.get('biological_replicate', ''),
                treatment=attrs.get('treatment', ''), sample_name=attrs.get('sample_name', ''),
                listed_fastq_bytes=sum(int(x) for x in expected[experiment]['fastq_bytes'].split(';')),
                source=p.relative_to(BASE).as_posix(),
                raw_reads_downloaded=False))
    assert len(sequence_rows) == 64
    write_tsv(OUT / 'tsutsui_sequencing_crosswalk.tsv', sequence_rows)
    availability = []
    for accession in ['GSE335749', 'GSE335750']:
        p = CACHE / ('tp53_' + accession + '.txt')
        inputs.append(p)
        text = p.read_text()
        hit = re.search(r'Accession .*?is currently private and is scheduled to be released on ([^<]+)', text)
        assert hit, 'Reassess availability instead of carrying forward an old hold'
        availability.append(dict(accession=accession, observed_status='private',
                                 scheduled_release=hit.group(1).rstrip('.'),
                                 manuscript_claim='publicly available',
                                 checked_utc=datetime.now(timezone.utc).isoformat(),
                                 source_sha256=sha(p)))
    write_tsv(OUT / 'tp53_accessibility.tsv', availability)
    report.write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),
        code_sha256=sha(Path(__file__)), helper_sha256=sha(Path(__file__).with_name('a1_regulatory_fate.py')),
        contract_sha256=sha(BASE / 'config/regulatory_fate_analysis.json'),
        inputs={p.relative_to(BASE).as_posix(): sha(p) for p in inputs},
        resolved_mice=22, listed_mice_without_retained_trace=4, sequencing_records=64,
        outputs={p.relative_to(BASE).as_posix(): sha(p) for p in OUT.glob('*.tsv')}), indent=2) + '\n')
    print('Resolved 22 mice, recorded 4 other listed mice, mapped 64 sequencing records; no raw reads downloaded.')


if __name__ == '__main__':
    main()
