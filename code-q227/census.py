"""Q227 census.  For every length-additive (S,T) in S~_n, n=3..7:

  chain   A --u_T--> nu --u_S--> lambda     (order verified in order.py)

Realisable (S,T) := exists a bead set A with the whole chain nonzero.
Claims tested:
  (R1)  realisable  <=>  v = u_S u_T is 321-avoiding
  (C1)  on the realisable locus, every run of S is usable and e = m
  (C2)  #Thm4.1 moves  =  #runs of S           (on the realisable locus)
  (C3)  #Thm4.1 moves  =  #available forward hops of nu   [the literal (C)]
  (C4)  excess moves (Thm4.1 minus union_x etilde_1) occur only off the locus
"""
from itertools import combinations
from beads import *

def all_beadsets(n):
    out = []
    for k in range(1, n):
        for c in combinations(range(n), k):
            out.append(frozenset(c))
    return out

def hops(A, n):
    return sum(1 for p in A if (p + 1) % n not in A)

def is321(win, n):
    """w has no x<y<z in Z with w(x)>w(y)>w(z).  Window w(1..n); w(i+kn)=w(i)+kn.
    Enough to search x in 1..n and y,z in a bounded range."""
    def W(i):
        q, r = divmod(i - 1, n)
        return win[r] + q * n
    lo, hi = 1 - 2 * n, 3 * n
    for x in range(1, n + 1):
        for y in range(x + 1, hi):
            if W(y) >= W(x):
                continue
            for z in range(y + 1, hi):
                if W(z) < W(y):
                    return False
    return True

tot = {}
for n in range(3, 8):
    pairs = list(additive_pairs(n))
    beads = all_beadsets(n)
    stats = dict(pairs=len(pairs), real=0, r_and_321=0, c1_fail=0,
                 c2_fail=0, c3_fail=0, c3_ok=0,
                 exc_pairs=0, exc_real=0, e_gt_m_offlocus=0, offlocus=0)
    for S, T, win in pairs:
        wS, wT = word_cd(S, n), word_cd(T, n)
        realA = []
        for A in beads:
            nu = act_word(A, wT, n)
            if nu is None:
                continue
            lam = act_word(nu, wS, n)
            if lam is not None:
                realA.append((A, nu, lam))
        real = bool(realA)
        a321 = is321(win, n)
        stats['real'] += real
        stats['r_and_321'] += (real == a321)
        nruns = len(runs(S, n)) if len(S) < n else 0
        letters = clio_letters(S, T, n)
        mv = clio_moves(S, T, n)
        # e = m for every usable run?
        e_is_m = all(r[0] in letters for r in runs(S, n) if
                     any(x in letters for x in r)) if len(S) < n else True
        allusable = (len(mv) == nruns)
        if real:
            if not (e_is_m and allusable):
                stats['c1_fail'] += 1
            if len(mv) != nruns:
                stats['c2_fail'] += 1
            for (A, nu, lam) in realA:
                if len(mv) == hops(nu, n):
                    stats['c3_ok'] += 1
                else:
                    stats['c3_fail'] += 1
        else:
            stats['offlocus'] += 1
            if not e_is_m:
                stats['e_gt_m_offlocus'] += 1
        # excess vs Morse-Schilling
        ims = set()
        for x in X_set(S, T, n):
            r = ms_etilde(S, T, x, n)
            if r:
                ims.add((r[1], r[2]))
        exc = mv - ims
        if exc:
            stats['exc_pairs'] += 1
            if real:
                stats['exc_real'] += 1
    tot[n] = stats
    print(n, stats)
print()
print("total additive pairs 3..7:", sum(t['pairs'] for t in tot.values()))
print("excess pairs by n:", [tot[n]['exc_pairs'] for n in range(3, 8)])
