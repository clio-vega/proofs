"""Two independent refutations of the literal (C), plus the replacement count."""
from itertools import combinations
from beads import *
import collections

def all_beadsets(n):
    return [frozenset(c) for k in range(1, n) for c in combinations(range(n), k)]
def hops(A, n):
    return sum(1 for p in A if (p + 1) % n not in A)

print("(K1) STRUCTURAL: is #hops(nu) even an invariant of v = u_S u_T?")
first = None; bad = 0; tot = 0
for n in range(3, 8):
    for S, T, win in additive_pairs(n):
        hs = set()
        for A in all_beadsets(n):
            nu = act_word(A, word_cd(T, n), n)
            if nu is None: continue
            if act_word(nu, word_cd(S, n), n) is None: continue
            hs.add(hops(nu, n))
        if not hs: continue
        tot += 1
        if len(hs) > 1:
            bad += 1
            if first is None:
                first = (n, sorted(S), sorted(T), sorted(hs))
    print(f"   n<={n}: realisable pairs {tot}, of which #hops(nu) is NOT constant: {bad}")
print("   smallest witness:", first)

print()
print("(K2) COUNT: #Thm4.1 moves vs #hops(nu), and the replacement count")
first2 = None
c = collections.Counter()
for n in range(3, 8):
    for S, T, win in additive_pairs(n):
        E = clio_letters(S, T, n)
        repl = sum(1 for r in runs(S, n) if (r[0] + 1) % n not in T) if len(S) < n else 0
        for A in all_beadsets(n):
            nu = act_word(A, word_cd(T, n), n)
            if nu is None: continue
            lam = act_word(nu, word_cd(S, n), n)
            if lam is None: continue
            c['triples'] += 1
            c['eq_hops_nu'] += (len(E) == hops(nu, n))
            c['eq_hops_lam'] += (len(E) == hops(lam, n))
            c['eq_hops_A'] += (len(E) == hops(A, n))
            c['eq_repl'] += (len(E) == repl)
            if len(E) != hops(nu, n) and first2 is None:
                first2 = (n, sorted(S), sorted(T), sorted(A), sorted(nu), sorted(lam),
                          len(E), hops(nu, n))
print("  ", dict(c))
print("   smallest (C)-failure  (n,S,T,A,nu,lam,#moves,#hops(nu)):", first2)
