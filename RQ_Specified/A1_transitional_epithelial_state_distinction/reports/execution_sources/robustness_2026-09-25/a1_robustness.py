"""Small deterministic operations shared by the A1 robustness analyses."""
from collections import Counter, defaultdict
from urllib.parse import unquote
import gzip
import hashlib
import math


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def weighted_summary(counts):
    """Input: (numerator, denominator) for each source, including zero counts."""
    if not counts or any(d <= 0 or not 0 <= n <= d for n, d in counts):
        raise ValueError('Invalid or empty source counts')
    return sum(n / d for n, d in counts) / len(counts), sum(n for n, d in counts) / sum(d for n, d in counts)


def adjusted_rand(left, right):
    """Partition agreement; invariant to arbitrary category names, no p-value."""
    if len(left) != len(right) or len(left) < 2:
        raise ValueError('Need aligned labels for at least two cells')
    pairs = lambda v: sum(n * (n - 1) // 2 for n in v)
    total = len(left) * (len(left) - 1) // 2
    joint = pairs(Counter(zip(left, right)).values())
    a, b = pairs(Counter(left).values()), pairs(Counter(right).values())
    expected, maximum = a * b / total, (a + b) / 2
    return 1. if maximum == expected else (joint - expected) / (maximum - expected)


def transcript_tss(path, loci):
    """Retain explicit parent-linked transcripts in the frozen native gene loci."""
    wanted = {r['gene']: r for r in loci}
    parents = defaultdict(set)
    transcripts = []
    with gzip.open(path, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) != 9 or f[2] not in ('gene', 'transcript', 'mRNA', 'ncRNA', 'lnc_RNA'):
                continue
            a = {k: unquote(v) for k, v in (x.split('=', 1) for x in f[8].split(';') if '=' in x)}
            symbol = a.get('gene_name', a.get('gene', ''))
            if symbol not in wanted:
                continue
            r = wanted[symbol]
            left, right = int(f[3]) - 1, int(f[4])
            if (f[0], f[6]) != (r['chrom'], r['strand']):
                continue
            if f[2] == 'gene':
                if (left, right) == (int(r['gene_start']), int(r['gene_end'])):
                    parents[symbol].add(a['ID'])
            else:
                transcripts.append((symbol, left, right, a.get('Parent', '').split(','), a.get('ID', ''), f[2]))
    found = defaultdict(set)
    audit = []
    for symbol, left, right, parent_ids, ident, kind in transcripts:
        r = wanted[symbol]
        eligible = (len(parents[symbol]) == 1 and bool(set(parent_ids) & parents[symbol])
                    and int(r['gene_start']) <= left < right <= int(r['gene_end']) and bool(ident))
        tss = left if r['strand'] == '+' else right - 1
        audit.append(dict(gene=symbol, transcript_id=ident, feature=kind, start=left, end=right,
                          tss=tss, parent_ids=';'.join(parent_ids), eligible=eligible))
        if eligible:
            found[symbol, tss].add(ident)
    rows = []
    for symbol, r in wanted.items():
        positions = {int(r['tss'])} | {t for g, t in found if g == symbol}
        for tss in sorted(positions):
            ids = sorted(found.get((symbol, tss), set()))
            rows.append(dict(gene=symbol, gene_group=r['gene_group'], chrom=r['chrom'], strand=r['strand'],
                             tss=tss, baseline_tss=int(r['tss']), is_baseline=tss == int(r['tss']),
                             transcript_count=len(ids), transcript_ids=';'.join(ids),
                             parent_gene_count=len(parents[symbol])))
    return rows, audit


def sign(value):
    return int(value > 0) - int(value < 0)
