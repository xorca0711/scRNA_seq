"""Verify scored estimands against per-gene evidence, source lists and floors."""
from pathlib import Path
import json
import hashlib
import sys
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0, str(ROOT))
from analysis.lib.repository_paths import recorded_file

def main():
    modules=json.loads((HERE/'modules.json').read_text(encoding='utf-8'))['modules']
    defs={m['name']:m for m in modules}
    s=pd.read_csv(HERE/'results/module_scores.csv')
    d=pd.read_csv(HERE/'results/gene_detection.csv',low_memory=False)
    e=pd.read_csv(HERE/'results/within_unit_effects.csv')
    assert s.score.dropna().between(0,1).all()
    assert s.coverage.between(0,1).all()
    assert not s.duplicated(['cohort','unit','seed','group','module']).any()
    assert (s.interpretable==((s.n_cells>=30)&(s.coverage>=.8))).all()
    lookup={k:g.set_index('gene') for k,g in d.groupby(['cohort','unit','group'])}
    checked=0
    for r in s[s.seed==17].itertuples():
        g=lookup[(r.cohort,r.unit,r.group)]
        genes=[x for x in defs[r.module]['genes'] if x in g.index and g.loc[x,'present']]
        expected=g.loc[genes,'detection'].mean()
        assert np.isclose(r.score,expected,equal_nan=True,atol=1e-12)
        assert len(genes)==r.n_genes
        checked+=1
    assert np.allclose(e.difference,e.score_two_marker-e.score_reference,equal_nan=True)
    expected=(e.n_two_marker>=30)&(e.n_reference>=30)&(e.coverage_two_marker>=.8)&(e.coverage_reference>=.8)
    assert (e.interpretable==expected).all()
    for m in modules:
        p=recorded_file(ROOT,m['provenance']['path'],m['provenance']['sha256'])
        assert hashlib.sha256(p.read_bytes()).hexdigest()==m['provenance']['sha256']
    main=e[e.seed==17].merge(e[e.seed==29],on=['cohort','unit','module'],suffixes=('_17','_29'))
    main['both_seeds_evaluable']=main.interpretable_17 & main.interpretable_29
    main['same_direction']=np.sign(main.difference_17)==np.sign(main.difference_29)
    main[['cohort','unit','module','difference_17','difference_29','both_seeds_evaluable','same_direction']].to_csv(HERE/'results/seed_sensitivity.csv',index=False)
    result={'status':'passed','primary_scores_reconstructed':checked,'source_modules_checked':len(modules),
            'two_seed_evaluable':int(main.both_seeds_evaluable.sum()),
            'same_direction_among_two_seed_evaluable':int((main.same_direction&main.both_seeds_evaluable).sum()),
            'checks':['scores reconstructed from gene detections','source file hashes','gene coverage',
                      'within-unit effect arithmetic','30-cell group floor','no duplicate unit/group/seed/module keys']}
    (HERE/'results/verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':main()
