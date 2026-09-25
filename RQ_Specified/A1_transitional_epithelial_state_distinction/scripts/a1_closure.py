"""Identity and annotation operations for the adaptive closure batch."""
import re


def cd44_key(value, kind):
    pattern = r'(.+)_CD44([+-])' if kind == 'column' else r'HTS053_\d+_(.+)_CD44(minus|plus)_S\d+_R[12]_001.fastq.gz'
    match = re.fullmatch(pattern, value)
    if not match:
        raise ValueError(f'Unrecognized CD44 {kind}: {value}')
    source, gate = match.groups()
    source = re.sub(r'^R26_([1278])$', r'R26\1', source)
    return source, 'positive' if gate in ('+', 'plus') else 'negative'


def annotation_audit(rows, mapping):
    total = len(rows)
    if not total:
        raise ValueError('Empty annotation scope')
    return dict(cells=total,
                unmapped_newleiden=sum(r['newleiden'] not in mapping for r in rows),
                cell_type_mismatches=sum(mapping.get(r['newleiden']) != r['cell type'] for r in rows),
                stringent_other=sum(r['clusterK12_stringent'] == 'other' for r in rows),
                literal_changes=sum(r['clusterK12'] != r['clusterK12_stringent'] for r in rows),
                retained_code_changes=sum(r['clusterK12_stringent'] != 'other' and
                                          r['clusterK12'] != r['clusterK12_stringent'] for r in rows))
