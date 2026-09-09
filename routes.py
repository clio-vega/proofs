"""
Route-level mechanism.  Two DISJOINT bead moves b->b+e and c->c+f.
CORRECTED: reversing the order shifts the e-height by -k and the f-height by +k,
but each order attributes them to different variables, so the WEIGHT RATIO is (ts)^k.
Hence the weight ratio of the two orders is exactly (ts)^k.
"""
from itertools import combinations

def ht(M, x, y):           # beads strictly between x and y
    return len([z for z in M if x < z < y])

bad = 0; tot = 0; kvals = {}
W = range(-6, 12)
import itertools, random
random.seed(11)
for trial in range(60000):
    M = frozenset(random.sample(list(W), random.randint(4, 9)))
    e = random.randint(1,5); f = random.randint(1,5)
    cand = [(b,c) for b in M for c in M if b != c]
    if not cand: continue
    b,c = random.choice(cand)
    # need both moves legal and independent (4 distinct sites)
    if len({b, c, b+e, c+f}) != 4: continue
    if (b+e) in M or (c+f) in M: continue
    # order 1: do e-move first, then f-move
    M1 = (M - {b}) | {b+e}
    h_e_first = ht(M, b, b+e)
    h_f_after = ht(M1, c, c+f)
    # order 2: do f-move first, then e-move
    M2 = (M - {c}) | {c+f}
    h_f_first = ht(M, c, c+f)
    h_e_after = ht(M2, b, b+e)
    k = (1 if b+e > c and b+e < c+f else 0) - (1 if b > c and b < c+f else 0)
    tot += 1
    d_e = h_e_after - h_e_first      # change in e-height when f goes first
    d_f = h_f_after - h_f_first      # change in f-height when e goes first
    kvals.setdefault((d_e, d_f, k), 0)
    kvals[(d_e, d_f, k)] += 1
    # order-1 total weight t^{h_e_first} s^{h_f_after};
    # order-2 total weight t^{h_e_after} s^{h_f_first}
    # claim: (h_e_after - h_e_first) == (h_f_after - h_f_first) == k
    if not (d_e == k and d_f == k):
        bad += 1
        if bad < 5:
            print(f"  MISMATCH M={sorted(M)} b={b} e={e} c={c} f={f}: "
                  f"d_e={d_e} d_f={d_f} k={k}")
print(f"tested {tot} disjoint move-pairs;  failures of (d_e == d_f == k): {bad}")
print(f"observed (d_e, d_f, k) triples: {sorted(kvals)}")
