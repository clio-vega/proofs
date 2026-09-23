"""Faithfulness of the Lean transport (2026-09-23 LEAN session).

`lean/tworow_d4_kernel/TworowD4Kernel/AffineAdditive.lean` does NOT build the affine
symmetric group.  It *defines* `Additive S T` to be the criterion of `lem:add` applied
to the block model `lem:blocks`.  That transport is only as good as the claim that the
resulting predicate IS `l(u_S u_T) = |S| + |T|`.

This script re-implements the Lean definitions verbatim (fuel-bounded `firstGapAux`,
`firstGap`, `uS`, `Additive` with the per-element run condition over lifts 0..n-1) and
compares them against `affine.py`'s length, computed by Shi's inversion formula on the
window of the composed affine permutation -- an independent route that never touches
the block model.

Result 2026-09-23: 0 disagreements over all 87376 pairs with 2 <= n <= 8.
"""
import sys; sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from affine import *
from itertools import combinations


def firstGapAux(S, n, j, f):
    k = 0
    while f > 0 and (j % n) in S:
        j += 1; k += 1; f -= 1
    return k


def firstGap(S, n, j):
    return j + firstGapAux(S, n, j, n)


def uS(S, n, j):
    return j - 1 if ((j - 1) % n) in S else firstGap(S, n, j)


def Additive(S, T, n):
    """Verbatim `AffineAdditive.Additive`."""
    if len(S) >= n or len(T) >= n:
        return False
    return all(uS(S, n, k) < uS(S, n, firstGap(T, n, k + 1))
               for k in range(n) if (k % n) in T)


def shi(S, T, n):
    """The group-theoretic statement, via Shi's inversion formula."""
    if len(S) >= n or len(T) >= n:
        return False
    return length(window(compose(u_S(S, n), u_S(T, n)), n), n) == len(S) + len(T)


def subs(n):
    return [frozenset(c) for r in range(n + 1) for c in combinations(range(n), r)]


if __name__ == '__main__':
    grand = graddis = 0
    for n in range(2, 9):
        tot = dis = 0
        for S in subs(n):
            for T in subs(n):
                tot += 1
                if Additive(S, T, n) != shi(S, T, n):
                    dis += 1
                    print('  DISAGREE', n, sorted(S), sorted(T))
        print(f'n={n}: {tot} pairs, {dis} disagreements')
        grand += tot; graddis += dis
    print(f'total {grand} pairs, {graddis} disagreements')
