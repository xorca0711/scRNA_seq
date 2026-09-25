"""Bounded source-document retrieval for identity and lineage audit (no raw reads)."""
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import hashlib, json, urllib.request

BASE = Path(__file__).resolve().parents[1]
SOURCES = {
    'cd44_source.xlsx': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-025-63735-1/MediaObjects/41467_2025_63735_MOESM4_ESM.xlsx',
    'tsutsui_source.xlsx': 'https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-68909-z/MediaObjects/41467_2026_68909_MOESM9_ESM.xlsx',
    'KatzenLab_repos.json': 'https://api.github.com/users/KatzenLab/repos?per_page=100',
    'katzen_tree.json': 'https://api.github.com/repos/KatzenLab/SFTPC_Mutant_Mice-Code/git/trees/main?recursive=1',
    'katzen_Manuscript_code': 'https://raw.githubusercontent.com/KatzenLab/SFTPC_Mutant_Mice-Code/0be3a649b6c8c72a5d64787a9ac6f3023ddba52d/Manuscript_code',
    'hpcs_code_tree.json': 'https://api.github.com/repos/dbetel/HPCS_LUAD/git/trees/main?recursive=1',
    'mintchip_tree.json': 'https://api.github.com/repos/jianhong/MintChIP/git/trees/master?recursive=1',
    'mintchip_main.nf': 'https://raw.githubusercontent.com/jianhong/MintChIP/b0b7eb7f8c0c4ee6f9c3c90f17b1b1207708f71a/main.nf',
    'pats_bioc.xml': 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC7461628/unicode',
    'hpcs_bioc.xml': 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC7745838/unicode',
}

def fetch(item):
    name, url = item
    p = BASE/'cache/followup_sources'/name
    p.parent.mkdir(parents=True, exist_ok=True)
    try:
        cached = p.exists()
        if not cached:
            req = urllib.request.Request(url, headers={'User-Agent':'scRNA_seq public research'})
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read(60*1024**2+1)
            if len(data)>60*1024**2: raise ValueError('Source exceeds 60 MiB ceiling')
            p.write_bytes(data)
        return dict(name=name, url=url, path=p.relative_to(BASE).as_posix(),
                    bytes=p.stat().st_size, sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                    cache_reused=cached, checked_utc=datetime.now(timezone.utc).isoformat())
    except Exception as exc:
        return dict(name=name,url=url,error=str(exc))

if __name__=='__main__':
    with ThreadPoolExecutor(max_workers=3) as pool:
        records=list(pool.map(fetch,SOURCES.items()))
    report=BASE/'reports/followup_source_inventory.json'
    report.write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
    for r in records: print(r['name'],r.get('bytes',r.get('error')))
