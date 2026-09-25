"""Native-assembly, exact-interval histone profiles and methylation-domain context.

Run through analysis/scripts/run_with_environment.py with the scientific site-packages.
Requires the isolated cache/pybigtools installation (version 0.3.0).
"""
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict
import gzip, hashlib, json, math, re, sys
from urllib.parse import unquote

BASE=Path(__file__).resolve().parents[1]
OUT=BASE/'tables/direct_marks_2026-09-25'

def sha(p):
    with p.open('rb') as h: return hashlib.file_digest(h,'sha256').hexdigest()

def interval_signal(records,start,end):
    """Integrate sorted non-overlapping half-open records, clipping both edges."""
    if end<=start: raise ValueError('Empty window')
    integral=0.; covered=0; previous=start
    for left,right,value in records:
        left=max(start,int(left));right=min(end,int(right))
        if right<=left: continue
        if left<previous: raise ValueError('Overlapping/unsorted BigWig records')
        if not math.isfinite(value) or value<0: raise ValueError('Invalid CPM value')
        integral+=(right-left)*value;covered+=right-left;previous=right
    return integral/(end-start),covered/(end-start)

def union_length(intervals,start,end):
    total=0;last=start
    for left,right in sorted(intervals):
        left=max(left,start,last);right=min(right,end)
        if right>left:total+=right-left;last=right
    return total

def annotations(path,assembly,groups):
    wanted={g:group for group,genes in groups.items() for g in genes}
    found=defaultdict(list)
    with gzip.open(path,'rt') as h:
        for line in h:
            if line.startswith('#'):continue
            f=line.rstrip().split('\t')
            if assembly=='T2T-CHM13v2.0':
                if len(f)!=9 or f[2]!='gene':continue
                a=dict(v.split('=',1) for v in f[8].split(';') if '=' in v)
                symbol=unquote(a.get('gene_name',a.get('gene','')))
                chrom,left,right,strand=f[0],int(f[3])-1,int(f[4]),f[6]
            else:
                if len(f)<13:continue
                symbol=f[12];chrom,strand,left,right=f[2],f[3],int(f[4]),int(f[5])
            if symbol in wanted and re.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)',chrom):
                found[symbol].append((chrom,left,right,strand))
    rows=[]
    for symbol,group in wanted.items():
        loci=list(set(found[symbol]))
        if assembly!='T2T-CHM13v2.0' and len({(v[0],v[3]) for v in loci})==1:
            loci=[(loci[0][0],min(v[1] for v in loci),max(v[2] for v in loci),loci[0][3])]
        row=dict(assembly=assembly,gene=symbol,gene_group=group,annotation_records=len(found[symbol]))
        if len(loci)!=1:row.update(status='held_missing_or_multiple_loci')
        else:
            c,l,r,s=loci[0]
            row.update(chrom=c,gene_start=l,gene_end=r,strand=s,tss=l if s=='+' else r-1,status='eligible')
        rows.append(row)
    return rows

def main():
    import pandas as pd
    sys.path.insert(0,str(BASE/'cache/pybigtools'))
    import pybigtools
    if OUT.exists():raise SystemExit('Refusing to overwrite direct-mark output')
    cp=BASE/'config/direct_mark_loci.json';cfg=json.loads(cp.read_text())
    ip=BASE/'reports/direct_mark_input_inventory.json';inventory=json.loads(ip.read_text())
    inputs=[r for r in inventory['files'] if 'path' in r and r['accession']!='GSE141635' and r.get('assembly')!='mm10']
    for row in inputs:
        assert sha(BASE/row['path'])==row['sha256'],row['path']
    record=dict(status='running',started_utc=datetime.now(timezone.utc).isoformat(),
                code_sha256=sha(Path(__file__)),contract_sha256=sha(cp),inventory_sha256=sha(ip),
                input_sha256={r['path']:r['sha256'] for r in inputs},pybigtools_version='0.3.0')
    rp=BASE/'reports/direct_mark_run.json';rp.write_text(json.dumps(record,indent=2)+'\n')
    OUT.mkdir(parents=True)
    def write(rows,name):pd.DataFrame(rows).to_csv(OUT/name,sep='\t',index=False)
    refs={r['assembly']:BASE/r['path'] for r in inputs if r['accession']=='reference'}
    loci=[]
    for assembly,p in refs.items():loci.extend(annotations(p,assembly,cfg['gene_groups']))
    write(loci,'loci.tsv')
    human=[r for r in loci if r['assembly']=='T2T-CHM13v2.0' and r['status']=='eligible']
    signals=[];tracks=[];sample_rows=[];chrom_reference=None
    for source in inputs:
        if not source['path'].endswith('.bw'):continue
        state,cut,mark=source['title'].split('_')
        sample_rows.append(dict(gsm=source['gsm'],state=state,replicate=cut,mark=mark,
                                assembly=source['assembly'],normalization='CPM; not spike-in scaled',
                                source_path=source['path'],sha256=source['sha256']))
        bw=pybigtools.open(str(BASE/source['path']));chroms=bw.chroms()
        if chrom_reference is None:chrom_reference=chroms
        assert chroms==chrom_reference,'Track chromosome sizes differ'
        for locus in human:
            common={k:locus[k] for k in ['assembly','gene','gene_group','chrom','strand','tss']}
            common.update(gsm=source['gsm'],state=state,replicate=cut,mark=mark)
            for window,width in [('promoter',cfg['promoter_half_width_bp']),('broad_promoter',cfg['broad_promoter_half_width_bp'])]:
                start=max(0,locus['tss']-width);end=min(chroms[locus['chrom']],locus['tss']+width)
                value,coverage=interval_signal(bw.records(locus['chrom'],start,end),start,end)
                signals.append({**common,'window':window,'start':start,'end':end,'mean_CPM':value,'covered_fraction':coverage})
            if locus['gene'] in cfg['locus_panels']:
                width=cfg['track_half_width_bp'];step=cfg['track_bin_bp']
                for offset in range(-width,width,step):
                    start=locus['tss']+offset;end=start+step
                    assert 0<=start<end<=chroms[locus['chrom']]
                    value,coverage=interval_signal(bw.records(locus['chrom'],start,end),start,end)
                    midpoint=offset+step/2
                    tracks.append({**common,'offset_bp':midpoint if locus['strand']=='+' else -midpoint,
                                   'start':start,'end':end,'mean_CPM':value,'covered_fraction':coverage})
        bw.close()
        print(source['gsm'],source['title'],'quantified',flush=True)
    write(sample_rows,'sample_manifest.tsv');write(tracks,'locus_tracks.tsv')
    h3={(r['gene'],r['state'],r['replicate'],r['window']):r for r in signals if r['mark']=='H3'}
    for row in signals:
        control=h3[row['gene'],row['state'],row['replicate'],row['window']]
        row['H3_mean_CPM']=control['mean_CPM'];row['H3_covered_fraction']=control['covered_fraction']
        usable=row['covered_fraction']>=.8 and control['covered_fraction']>=.8 and row['mean_CPM']>0 and control['mean_CPM']>0
        row['ratio_eligible']=usable
        row['log2_mark_over_H3']=math.log2(row['mean_CPM']/control['mean_CPM']) if usable else None
    write(signals,'histone_signals.tsv')
    contrasts=[]
    index={(r['gene'],r['state'],r['replicate'],r['window'],r['mark']):r for r in signals}
    for row in signals:
        if row['state']!='iATCs' or row['mark']=='H3':continue
        for comparator in ['iAT2','iAT1']:
            other=index[row['gene'],comparator,row['replicate'],row['window'],row['mark']]
            eligible=row['ratio_eligible'] and other['ratio_eligible']
            contrasts.append(dict(gene=row['gene'],gene_group=row['gene_group'],mark=row['mark'],window=row['window'],
                replicate=row['replicate'],contrast='iATCs_minus_'+comparator,eligible=eligible,
                log2_ratio_difference=row['log2_mark_over_H3']-other['log2_mark_over_H3'] if eligible else None,
                raw_mark_log2_ratio=math.log2(row['mean_CPM']/other['mean_CPM']) if row['mean_CPM']>0 and other['mean_CPM']>0 else None,
                H3_log2_ratio=math.log2(row['H3_mean_CPM']/other['H3_mean_CPM']) if row['H3_mean_CPM']>0 and other['H3_mean_CPM']>0 else None))
    write(contrasts,'histone_contrasts.tsv')
    # Author tabular domains start at 1 (R/GRanges convention). Also expose the
    # alternative BED interpretation; at most one base changes per boundary.
    domains=[];hg19=[r for r in loci if r['assembly']=='hg19' and r['status']=='eligible']
    for source in inputs:
        if source['accession']!='GSE150527':continue
        day=re.search(r'D[046]',source['path']).group()
        d=pd.read_csv(BASE/source['path'],sep='\t')
        if 'type' not in d:d['type']='PMD'
        for locus in hg19:
            start=max(0,locus['tss']-cfg['promoter_half_width_bp']);end=locus['tss']+cfg['promoter_half_width_bp']
            for kind,sub in d[d.chr==locus['chrom']].groupby('type'):
                for origin,shift in [('one_based_inclusive',1),('zero_based_half_open_sensitivity',0)]:
                    intervals=[(int(r.start)-shift,int(r.end)) for r in sub.itertuples() if r.end>start and r.start-shift<end]
                    covered=union_length(intervals,start,end)
                    domains.append(dict(gene=locus['gene'],gene_group=locus['gene_group'],day=day,domain=kind,
                        chrom=locus['chrom'],start=start,end=end,coordinate_assumption=origin,
                        domain_overlap_fraction=covered/(end-start),source_path=source['path']))
    write(domains,'methylation_domain_overlap.tsv')
    record.update(status='completed',finished_utc=datetime.now(timezone.utc).isoformat(),
        output_sha256={p.relative_to(BASE).as_posix():sha(p) for p in sorted(OUT.iterdir())})
    rp.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    print('completed',len(signals),'histone windows;',len(domains),'domain intersections')

if __name__=='__main__':main()
