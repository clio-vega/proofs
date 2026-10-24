"""T10.  Is the P1 falsifier a LAW or an accident?

In M={-7,2,3,7,-2} (lambda=2 negation pairs, 1 free element, r=3=n-lambda) the four
surviving terms cancelled in two pairs, matched by the involution sigma that swaps
each negation pair and fixes the free elements.  Test whether this is a law:
  lambda negation pairs, f free elements, n = 2*lambda+f, r = lambda+f, m_j = 2(j-1).
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt

random.seed(4242)
POOL = [F(2),F(3),F(5),F(7),F(11),F(13),F(4),F(9),F(8),F(5,2),F(7,3),F(9,4),F(6),F(10)]

print(f"{'lam':>3} {'f':>2} {'n':>2} {'r':>2} {'draws':>6} {'Psi=0':>6}   verdict")
for lam in range(1, 4):
    for f in range(0, 3):
        n = 2*lam + f; r = lam + f
        if r < 1 or r > n: continue
        theta = list(range(r)); m = [theta[j]+j for j in range(r)]
        if len(set(m)) != r: continue
        zer = tot = 0
        for _ in range(25):
            vals = random.sample(POOL, lam + f)
            w = []
            for v in vals[:lam]: w += [v, -v]
            w += vals[lam:]
            if len(set(w)) != n or any(x == 0 for x in w): continue
            # keep it reciprocal-free so no structural vanishing can be confused for this
            if any(w[i]*w[j] == 1 for i in range(n) for j in range(i, n)): continue
            try:
                a, b = psi_headsum(w, theta), psi_alt(w, theta)
            except ZeroDivisionError: continue
            assert a == b, (w, a, b)
            tot += 1; zer += (a == 0)
        if tot == 0: continue
        verdict = "ALWAYS ZERO -- law" if zer == tot else ("never" if zer == 0 else "MIXED")
        print(f"{lam:>3} {f:>2} {n:>2} {r:>2} {tot:>6} {zer:>6}   {verdict}")

print()
print("refine: fix f=1, sweep lambda and r around the threshold r = n - lambda")
print(f"{'lam':>3} {'f':>2} {'n':>2} {'r':>2} {'n-lam':>5} {'draws':>6} {'Psi=0':>6}")
for lam in range(1, 4):
    f = 1; n = 2*lam + f
    for r in range(1, n+1):
        theta = list(range(r)); m = [theta[j]+j for j in range(r)]
        if len(set(m)) != r: continue
        zer = tot = 0
        for _ in range(20):
            vals = random.sample(POOL, lam + f)
            w = []
            for v in vals[:lam]: w += [v, -v]
            w += vals[lam:]
            if len(set(w)) != n or any(x == 0 for x in w): continue
            if any(w[i]*w[j] == 1 for i in range(n) for j in range(i, n)): continue
            try: a = psi_headsum(w, theta)
            except ZeroDivisionError: continue
            tot += 1; zer += (a == 0)
        if tot: print(f"{lam:>3} {f:>2} {n:>2} {r:>2} {n-lam:>5} {tot:>6} {zer:>6}")
