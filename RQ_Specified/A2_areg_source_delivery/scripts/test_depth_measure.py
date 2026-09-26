"""Check the depth-standardised expectation against scipy's hypergeometric distribution.

The C37 correction carries the same test, for the same reason: the measure in
07_leg2_depth_standardised.py is a closed-form expression built from log gamma functions,
and a sign or an off-by-one in it would change every number it produces without failing
loudly. This exercises it against an independent implementation on small cases where the
answer can also be checked by hand, and on the edge cases the formula branches on.

Run with the bundled x64 interpreter and the venv on PYTHONPATH.
"""
from __future__ import annotations

import sys

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
from scipy.special import gammaln
from scipy.stats import hypergeom


def expectation(total: float, count: float, budget: int) -> float:
    """The expression 07_leg2_depth_standardised.py uses, isolated for testing."""
    rest = total - count
    if rest < budget:
        return 1.0
    p_none = np.exp(
        gammaln(rest + 1) - gammaln(rest - budget + 1)
        + gammaln(total - budget + 1) - gammaln(total + 1)
    )
    return float(1.0 - p_none)


def main() -> None:
    failures = []

    def check(label: str, got: float, want: float, tol: float = 1e-9) -> None:
        ok = abs(got - want) <= tol
        print("%s %-58s got %.12f want %.12f" % ("PASS" if ok else "FAIL", label, got, want))
        if not ok:
            failures.append(label)

    # 1. Against scipy over a grid. P(detect at least one) = 1 - P(X = 0) where X is
    #    hypergeometric with population total, successes count, draws budget.
    grid = 0
    for total in [10, 25, 100, 1000, 5000]:
        for count in [0, 1, 2, 5, 17]:
            for budget in [1, 2, 5, 10, 50]:
                if budget > total or count > total:
                    continue
                want = 1.0 - float(hypergeom.pmf(0, total, count, budget))
                got = expectation(float(total), float(count), budget)
                grid += 1
                if abs(got - want) > 1e-9:
                    failures.append("grid total=%d count=%d budget=%d" % (total, count, budget))
    print("PASS %-58s %d combinations agree with scipy" % ("grid against scipy.stats.hypergeom", grid)
          if not failures else "FAIL grid against scipy")

    # 2. Hand-checkable cases.
    check("no counts gives zero detection", expectation(1000, 0, 100), 0.0)
    check("drawing the whole library always detects", expectation(100, 1, 100), 1.0)
    # one success in 10, drawing 1: probability 1/10
    check("one of ten, single draw", expectation(10, 1, 1), 0.1)
    # two successes in 10, drawing 1: 2/10
    check("two of ten, single draw", expectation(10, 2, 1), 0.2)
    # one success in 10, drawing 2: 1 - C(9,2)/C(10,2) = 1 - 36/45 = 0.2
    check("one of ten, two draws", expectation(10, 1, 2), 0.2)

    # 3. The branch the formula takes when the non-success pool is smaller than the budget:
    #    every draw must include a success, so detection is certain.
    check("non-success pool below the budget is certain", expectation(100, 95, 10), 1.0)
    check("boundary, non-success pool exactly the budget", expectation(100, 90, 10),
          1.0 - float(hypergeom.pmf(0, 100, 90, 10)))

    # 4. Monotonicity, which is what makes the measure a depth correction at all: at a fixed
    #    budget more counts cannot lower detection, and at fixed counts a deeper library
    #    cannot raise it.
    rising = [expectation(5000, c, 1000) for c in [0, 1, 5, 20, 100]]
    print("PASS %-58s %s" % ("detection rises with counts at a fixed budget",
                             all(a <= b for a, b in zip(rising, rising[1:]))))
    if not all(a <= b for a, b in zip(rising, rising[1:])):
        failures.append("monotone in counts")
    falling = [expectation(t, 10, 1000) for t in [1000, 2000, 5000, 20000]]
    print("PASS %-58s %s" % ("detection falls as the library deepens at fixed counts",
                             all(a >= b for a, b in zip(falling, falling[1:]))))
    if not all(a >= b for a, b in zip(falling, falling[1:])):
        failures.append("monotone in depth")

    # 5. Rate invariance, and the regime in which it actually holds. A gene at a fixed rate
    #    should give the same standardised value at every library size. That is true only
    #    when the budget is small relative to the smallest library, and it degenerates as
    #    the budget approaches the library size, because a cell sampled in its entirety
    #    returns its raw detection. This deposit is floored at 1,000 molecules, so the
    #    1,000-molecule budget the specification inherited from C37 sits exactly at that
    #    degenerate edge. The numbers below are a property of the formula, computed with no
    #    cohort value, and they are reported as a limitation rather than used to pick a
    #    budget after the fact.
    rate = 0.002
    sizes = [1000, 1500, 2000, 3185, 5723, 9961, 20000, 94444]
    spreads = {}
    for budget in [1000, 500, 100]:
        vals = [expectation(t, rate * t, budget) for t in sizes]
        spreads[budget] = max(vals) - min(vals)
        print("     budget %4d: spread across library sizes %.4f, binomial limit %.4f"
              % (budget, spreads[budget], 1 - (1 - rate) ** budget))
    print("PASS %-58s %.4f" % ("a small budget is rate invariant", spreads[100])
          if spreads[100] < 0.01 else "FAIL a small budget is rate invariant")
    if spreads[100] >= 0.01:
        failures.append("rate invariance at a small budget")
    print("PASS %-58s %.4f" % ("the inherited 1000 budget is not, and this is recorded",
                               spreads[1000])
          if spreads[1000] > spreads[100] else "FAIL degeneration at the large budget")
    if not spreads[1000] > spreads[100]:
        failures.append("degeneration ordering")

    print()
    if failures:
        print("%d failures: %s" % (len(failures), failures[:6]))
        raise SystemExit(1)
    print("all checks pass; the measure matches scipy and behaves as a depth correction")


if __name__ == "__main__":
    main()
