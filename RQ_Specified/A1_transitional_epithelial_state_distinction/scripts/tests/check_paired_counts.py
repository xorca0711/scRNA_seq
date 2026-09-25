"""Synthetic software check only; never included among scientific outputs."""
import csv,hashlib,json,subprocess
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from a1_contract import ROOT as root
out=root/'tmp/a1_synthetic_qll';out.mkdir(exist_ok=True)
samples=out/'samples.tsv';counts=out/'counts.tsv'
with samples.open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['sample_id','biological_unit_id','group'])
 for i in range(8):w.writerow([f's{i}',f'unit{i//2}','positive' if i%2 else 'negative'])
with counts.open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['feature_id',*[f's{i}' for i in range(8)]])
 for g in range(500):
  values=[25+g%40+(g*(i+1)%19)+i//2*3 for i in range(8)]
  if g<20:values=[v*4 if i%2 else v for i,v in enumerate(values)]
  w.writerow([f'gene{g}',*values])
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
spec=dict(contrast=dict(status='frozen',input_scale='raw_integer_counts',reference='negative',case='positive',minimum_independent_pairs=3,min_count=10,min_total_count=15),counts_path=str(counts.resolve()),samples_path=str(samples.resolve()),counts_sha256=sha(counts),samples_sha256=sha(samples))
(out/'input_contract.json').write_text(json.dumps(spec))
parameters=dict(spec['contrast'],counts_path=spec['counts_path'],samples_path=spec['samples_path'],
                counts_md5=hashlib.md5(counts.read_bytes()).hexdigest(),samples_md5=hashlib.md5(samples.read_bytes()).hexdigest())
with (out/'input_contract.tsv').open('w',newline='') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['key','value']);w.writerows(parameters.items())
r=root/'analysis/corrections/statistics/.tools/R-portable/app/bin/Rscript.exe'
script=root/'RQ_Specified/A1_transitional_epithelial_state_distinction/scripts/04_fit_paired_counts.R'
result=subprocess.run([str(r),str(script),str(out)],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
if result.returncode:print(result.stdout);raise SystemExit(result.returncode)
rows=list(csv.DictReader((out/'paired_effects.tsv').open(),delimiter='\t'))
assert len(rows)==500
assert all(float(row['logFC'])>1 for row in rows if int(row['feature_id'][4:])<20)
summary=next(csv.DictReader((out/'fit_summary.tsv').open(),delimiter='\t'));assert summary['coefficient_tested']=='grouppositive'
assert int(summary['independent_pairs'])==4
print('Synthetic software check passed: 500 features, 4 pairs, intended case-minus-reference direction. No real study data analysed.')
