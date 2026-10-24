"""T1. Verify the PAIR FORM of B(T):

    B(T) = N(T) / D(T),
    N(T) = prod over unordered pairs {a,b}, a<=b, with {a,b} cap T != empty, of (w_a w_b - 1)
    D(T) = e^r * (prod_{a in T} w_a) * prod_{a in T, t notin T} (w_a - w_t),    e = prod_all w.

This is the identity on which Theorem A (structural locus) rests.  Checked against B()
from the 2026-10-22 alternant expansion, exact Fraction arithmetic only.
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import B

def N_pair(w, T):
    n = len(w); Ts = set(T); out = F(1)
    for a in range(n):
        for b in range(a, n):
            if a in Ts or b in Ts:          # {a,b} meets T  (a==b: a in T)
                out *= (w[a]*w[b] - 1)
    return out

def D_pair(w, T):
    n = len(w); r = len(T); Ts = set(T)
    e = F(1)
    for x in w: e *= x
    out = e**r
    for a in T: out *= w[a]
    for a in T:
        for t in range(n):
            if t not in Ts: out *= (w[a] - w[t])
    return out

random.seed(20261024)
def rand_w(n, bound=9):
    """distinct, nonzero, NO hygiene exclusions -- we WANT the special loci reachable."""
    while True:
        w = [F(random.randint(-bound, bound), random.randint(1, 4)) for _ in range(n)]
        if any(x == 0 for x in w): continue
        if len(set(w)) != n: continue
        return w

ok = bad = 0
fails = []
for n in range(1, 7):
    for r in range(1, min(n, 4)+1):
        for trial in range(40):
            w = rand_w(n)
            for T in itertools.combinations(range(n), r):
                try:
                    lhs = B(w, T)
                except ZeroDivisionError:
                    continue
                d = D_pair(w, T)
                if d == 0: continue
                rhs = N_pair(w, T)/d
                if lhs == rhs: ok += 1
                else:
                    bad += 1
                    if len(fails) < 4: fails.append((n, r, T, w, lhs, rhs))
print(f"T1 pair form of B(T):  {ok}/{ok+bad} agree   (n<=6, r<=4, 40 random M per cell)")
for f in fails: print("  FAIL", f)
