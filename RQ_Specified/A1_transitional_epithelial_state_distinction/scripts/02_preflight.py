"""Audit saved catalogs and candidate designs; no downloads or model fitting."""
import json
from a1_contract import STUDY,readiness


def main():
    catalog=[]
    for path in sorted((STUDY/'metadata').glob('GSE*.json')):
        data=json.loads(path.read_text(encoding='utf-8'))
        ids=[s['accession'] for s in data['samples']]
        assert len(ids)==len(set(ids))==data['sample_count'],path
        assert set(ids)==set(data['series']['sample_id']),path
        catalog.append(dict(accession=data['accession'],GSM_records=len(ids),title=data['series']['title'][0],
                            source_url=data['source_url'],response_sha256=data['response_sha256']))
    samples=json.loads((STUDY/'config/samples.json').read_text())['samples']
    contracts=json.loads((STUDY/'config/contrasts.json').read_text())
    contrasts=[]
    for c in contracts['contrasts']:
        holds,selected=readiness(c,samples)
        contrasts.append(dict(id=c['id'],design_ready=not holds,holds=holds,sample_records=len(selected)))
    report=dict(status='metadata_audit_complete',scientific_analysis_executed=False,
                execution_authorized=contracts['execution_authorized'],catalog=catalog,contrasts=contrasts)
    directory=STUDY/'reports';directory.mkdir(exist_ok=True)
    (directory/'catalog_audit.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    lines=['# A1 catalog and design preflight','',
           'Metadata-only report. No assay matrix has been analysed and no biological effect estimated.','',
           f'{len(catalog)} GEO series; {sum(x["GSM_records"] for x in catalog)} GSM records (not a biological-unit count).','',
           '| Prospective contrast | Sample records | Status |','|---|---:|---|']
    for c in contrasts:lines.append(f'| {c["id"]} | {c["sample_records"]} | '+('Design ready; execution separate' if c['design_ready'] else 'Held: identity/processed-input verification pending')+' |')
    lines.extend(['','Exact holds and catalog sources: [catalog_audit.json](catalog_audit.json).',
                  'See [the study map](../STUDY_MAP.md) for descriptive assays and other work packages.',''])
    (directory/'PREFLIGHT.md').write_text('\n'.join(lines),encoding='utf-8')
    print(f'Metadata audit: {len(catalog)} series; {len(contrasts)} proposed contrasts; '+
          f'{sum(c["design_ready"] for c in contrasts)} designs ready. Scientific analysis not run.')


if __name__=='__main__':main()
