"""Retrieve only frozen scoring features from the author's public cellxgene API."""
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
import urllib.parse
import urllib.request
from exploratory_common import BASE, TABLES, save_json, sha
import pandas as pd


def main():
    result=json.loads((BASE/'stage_E2_result.json').read_text())
    assert result['checkpoint_passed'], 'E2 must be reviewed before transfer retrieval'
    config=json.loads((BASE/'exploratory_config.json').read_text())
    frozen=json.loads((BASE/'frozen_repair_program.json').read_text())
    variants=json.loads((BASE/'frozen_control_variant.json').read_text())
    feature=pd.read_csv(TABLES/'developmental_feature_index.csv')
    modules={'candidate':frozen['genes'],'candidate_without_generic_controls':variants['genes'],
        **{name:config['control_modules'][name] for name in config['generic_control_names']},
        'ADI_published_holdout':config['control_modules']['ADI_published_holdout']}
    requested=set(config['reference_genes']).union(*(set(x) for x in modules.values()))
    genes=[g for g in feature.gene if g in requested]
    folder=BASE/'cache/expression/negretti_chunks';folder.mkdir(exist_ok=True)
    endpoint='https://lungcells.app.vumc.org/public/sucre/mouse_development_epithelium_revision/api/v0.2/data/var'
    batches=[genes[i:i+32] for i in range(0,len(genes),32)]

    def fetch(batch):
        tag=hashlib.sha256(json.dumps(batch).encode()).hexdigest()[:16]
        path=folder/f'{tag}.bin';provenance=folder/f'{tag}.json'
        if path.exists() and provenance.exists():
            record=json.loads(provenance.read_text());assert record['sha256']==sha(path)
            return record
        url=endpoint+'?'+urllib.parse.urlencode([('var:gene',gene) for gene in batch])
        request=urllib.request.Request(url,headers={'Accept':'application/octet-stream','User-Agent':'A0-exploratory-scRNA/1.0'})
        with urllib.request.urlopen(request,timeout=60) as response:
            payload=response.read(5_000_001)
        assert len(payload)<=5_000_000
        path.write_bytes(payload)
        record={'file':str(path.relative_to(BASE)).replace('\\','/'),'requested_genes':batch,
            'requested_url':url,'bytes':len(payload),'sha256':hashlib.sha256(payload).hexdigest(),
            'retrieved_at_utc':datetime.now(timezone.utc).isoformat(),'expression_type':'author SCT data-slot export'}
        save_json(provenance,record)
        return record

    records=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        for i,record in enumerate(pool.map(fetch,batches)):
            records.append(record)
            if (i+1)%10==0:print(f'Developmental feature batches retrieved: {i+1}/{len(batches)}',flush=True)
    save_json(BASE/'developmental_expression_manifest.json',{'frozen_program_sha256':sha(BASE/'frozen_repair_program.json'),
        'n_requested_features':len(genes),'chunks':records})
    print(json.dumps({'features':len(genes),'chunks':len(records),'bytes':sum(r['bytes'] for r in records)}),flush=True)


if __name__=='__main__':
    main()
