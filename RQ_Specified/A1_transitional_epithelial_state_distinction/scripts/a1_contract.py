"""Shared validation for prospective paired bulk-count contrasts."""
from collections import Counter,defaultdict
import hashlib
from pathlib import Path

STUDY=Path(__file__).resolve().parents[1]
ROOT=STUDY.parents[1]


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()


def readiness(contrast,samples,root=ROOT):
    """Require complete verified pairs; never turn title IDs into replicates."""
    reasons=[]
    wanted=contrast['sample_ids']
    selected=[s for s in samples if s['sample_id'] in wanted]
    if len(set(wanted))!=len(wanted) or Counter(s['sample_id'] for s in selected)!=Counter(wanted):
        reasons.append('Missing or duplicate sample IDs')
    units=defaultdict(list)
    members={}
    for sample in selected:
        sid=sample['sample_id']
        if sample.get('accession')!=contrast['accession'] or sample.get('assay')!=contrast['assay']:
            reasons.append(f'{sid}: accession/assay mismatch')
        if (sample.get('identity_verified') is not True or not sample.get('identity_evidence')
                or not sample.get('biological_unit_id') or not sample.get('unit_kind')
                or not sample.get('member_ids')):
            reasons.append(f'{sid}: biological identity/pool membership unverified')
        if not sample.get('counts_column'):reasons.append(f'{sid}: count-column crosswalk missing')
        unit=sample.get('biological_unit_id')
        if unit:
            units[unit].append(sample.get('group'))
            composition=set(sample.get('member_ids') or [])
            if unit in members and members[unit]!=composition:
                reasons.append(f'{unit}: inconsistent pool membership')
            members[unit]=composition
    if len(units)<max(3,contrast['minimum_independent_pairs']):reasons.append('Too few verified independent pairs')
    for unit,groups in units.items():
        if Counter(groups)!=Counter([contrast['reference'],contrast['case']]):
            reasons.append(f'{unit}: requires one aggregated library per sort; incomplete/repeated pair')
    names=list(members)
    for i,a in enumerate(names):
        for b in names[i+1:]:
            if members[a]&members[b]:reasons.append(f'{a}/{b}: overlapping biological members')
    path=contrast.get('counts_path')
    if not path or not contrast.get('counts_sha256'):reasons.append('Processed count payload/path/hash not verified')
    else:
        target=(root/path).resolve()
        if not target.is_relative_to(root.resolve()):reasons.append('Count path is outside repository')
        elif not target.is_file() or digest(target)!=contrast['counts_sha256']:reasons.append('Count file missing or hash mismatch')
    if contrast.get('input_scale')!='raw_integer_counts':reasons.append('Raw integer-count scale unverified')
    if contrast.get('assay_qc_verified') is not True or not contrast.get('assay_qc_evidence'):
        reasons.append('Assay QC, feature universe and normalization suitability not verified')
    if len({s.get('genotype') for s in selected})>1:
        reasons.append('Paired pilot requires one genotype stratum')
    columns=[s.get('counts_column') for s in selected if s.get('counts_column')]
    if len(set(columns))!=len(columns):reasons.append('Count column reused')
    return sorted(set(reasons)),selected
