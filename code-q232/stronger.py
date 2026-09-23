"""Q232: test the STRONGER statement my Case-analysis appears to prove.

Q232 (bracketed):  (S,T) additive, i in S, i+1 in T  ==>  (S-i, T-(i+1)) additive.
STRONGER:          (S,T) additive, i in S  (NO hypothesis on i+1)
                                          ==>  (S-i, T\{i+1}) additive.
When i+1 not in T this reads: (S-i, T) is additive.  That is the new content.
"""
import sys
sys.path.insert(0, '/home/clio/projects/proofs/code-q220')
from affine import *
from itertools import combinations


def is_additive(S, T, n):
    if len(S) >= n or len(T) >= n:
        return False
    return length(window(compose(u_S(S, n), u_S(T, n)), n), n) == len(S) + len(T)


def subs(n):
    return [frozenset(c) for r in range(n) for c in combinations(range(n), r)]


print("STRONGER claim: (S,T) additive, i in S  ==>  (S-{i}, T-{i+1}) additive")
print(f"{'n':>2} {'add.pairs':>10} {'unbracketed (i+1 notin T)':>26} {'NON-additive':>13}")
tot = bad = 0
wit = []
for n in range(3, 9):
    sl = subs(n)
    npairs = ncase = nbad = 0
    for S in sl:
        for T in sl:
            if not is_additive(S, T, n):
                continue
            npairs += 1
            for i in range(n):
                if i not in S:
                    continue
                if (i + 1) % n in T:        # that is the already-tested bracketed case
                    continue
                ncase += 1
                if not is_additive(frozenset(S) - {i}, frozenset(T), n):
                    nbad += 1
                    if len(wit) < 8:
                        wit.append((n, sorted(S), sorted(T), i))
    print(f"{n:>2} {npairs:>10} {ncase:>26} {nbad:>13}")
    tot += ncase
    bad += nbad
print(f"total: {tot} UNbracketed instances n=3..8, {bad} non-additive reductions")
for w in wit:
    print("  witness n=%d S=%s T=%s i=%d" % w)
