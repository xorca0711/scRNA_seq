"""Small, explicit calculations for nested source measurements."""
import csv
import hashlib
import itertools
import math
import statistics
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_tsv(path, rows):
    rows = list(rows)
    assert rows, 'Empty output requires an explicit schema'
    path = Path(path)
    assert not path.exists(), f'Refusing to overwrite {path}'
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def describe(values):
    assert values and all(math.isfinite(x) for x in values)
    return dict(n=len(values), mean=statistics.mean(values),
                sd=statistics.stdev(values) if len(values) > 1 else None,
                minimum=min(values), maximum=max(values))


def pooled_fraction(positives, totals):
    assert len(positives) == len(totals) and positives
    assert all(0 <= k <= n and n > 0 for k, n in zip(positives, totals))
    return sum(positives) / sum(totals)


def region_interaction(wt, mutant):
    """Each observation is one mouse's de-novo minus in-situ percentage."""
    return statistics.mean(mutant) - statistics.mean(wt)


def exact_permutation(wt, mutant):
    """Two-sided difference of means; all group allocations, no plus-one fudge."""
    all_values = list(wt) + list(mutant)
    observed = region_interaction(wt, mutant)
    null = []
    for idx in itertools.combinations(range(len(all_values)), len(wt)):
        chosen = set(idx)
        a = [x for i, x in enumerate(all_values) if i in chosen]
        b = [x for i, x in enumerate(all_values) if i not in chosen]
        null.append(region_interaction(a, b))
    extreme = sum(abs(x) >= abs(observed) - 1e-12 for x in null)
    return dict(effect_pp=observed, permutations=len(null), extreme=extreme,
                p_two_sided=extreme / len(null))


def sign_class(a, b):
    if a > 0 and b > 0:
        return 'both_up'
    if a < 0 and b < 0:
        return 'both_down'
    if a * b < 0:
        return 'opposite'
    return 'includes_zero'
