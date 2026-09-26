#!/usr/bin/env python3
r"""Differential check of the REMARK in TworowD4Kernel/SnpLatticePoints.lean:

    the clause `x >= 0` of Q is redundant given the rest.

Claim: if x in Z^l satisfies x([l]) = |lhat| and x(S) <= Lambda_{|S|} for every
S subseteq [l], then x >= 0 automatically.  Reason on paper: for each i,
x_i = |lhat| - x([l]\{i}) >= Lambda_l - Lambda_{l-1} = lhat_l >= 0.

This script does NOT assume that reason; it enumerates integer vectors with
negative entries allowed and reports any counterexample.  A count of 0 violations
together with a nonzero count of ENUMERATED feasible points is the evidence
(an empty feasible set would make the check vacuous).
"""
from itertools import combinations, product

def partitions(n, maxlen, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0:
        yield ()
        return
    if maxlen == 0:
        return
    for first in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - first, maxlen - 1, first):
            yield (first,) + rest

feasible = 0
violations = []
checked_pairs = 0
for l in range(1, 5):
    for d in range(0, 8):
        for lhat in partitions(d, l):
            lh = list(lhat) + [0] * (l - len(lhat))
            Lam = [sum(lh[:r]) for r in range(l + 1)]           # Lambda_0..Lambda_l
            subsets = [S for r in range(l + 1) for S in combinations(range(l), r)]
            checked_pairs += 1
            lo = -(d + 2)
            for x in product(range(lo, d + 1), repeat=l):
                if sum(x) != d:
                    continue
                if any(sum(x[i] for i in S) > Lam[len(S)] for S in subsets):
                    continue
                feasible += 1
                if any(xi < 0 for xi in x):
                    violations.append((l, lhat, x))

print(f"(lhat, l) pairs enumerated : {checked_pairs}")
print(f"feasible integer points    : {feasible}   (nonzero => check not vacuous)")
print(f"violations of x >= 0       : {len(violations)}")
for v in violations[:5]:
    print("   ", v)
