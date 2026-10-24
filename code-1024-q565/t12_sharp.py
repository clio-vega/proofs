"""T12.  The SHARP negation law, read off the t11 grid.

Grid: M2 fills r >= n-lambda+1 (Prop:neg).  M3 fills exactly f+2 <= r <= f+lambda.
Since n - lambda = lambda + f, BOTH are the single condition

        r >= f + 2          (f = number of UNPAIRED elements)

Reading: a term survives only if T takes at most one member of each negation
pair, so T must touch at least r-f pairs.  r >= f+2 says every surviving T
touches TWO OR MORE pairs -- and then the terms cancel.  r = f+1 means some T
touches only one pair, and no vanishing is seen.
"""
import sys, itertools, random
sys.path.insert(0, '/home/clio/projects/proofs/code-1023-q559')
from fractions import Fraction as F
from psi import psi_headsum, psi_alt, B, alternant_ratio

random.seed(90210)
POOL = [F(2),F(3),F(5),F(7),F(11),F(13),F(4),F(9),F(8),F(5,2),F(7,3),F(9,4),F(6),F(10),
        F(17),F(19),F(23),F(12,5),F(15,4),F(21,2)]

def draw(lam, f):
    vals = random.sample(POOL, lam + f)
    w = []
    for v in vals[:lam]: w += [v, -v]
    w += vals[lam:]
    n = len(w)
    if len(set(w)) != n or any(x == 0 for x in w): return None
    if any(w[i]*w[j] == 1 for i in range(n) for j in range(i, n)): return None
    return w

print("SHARP LAW:  r >= f+2  =>  Psi = 0   (theta = (0,1,...,r-1))")
print(f"{'lam':>3} {'f':>2} {'n':>2} {'r':>2} {'r>=f+2':>7} {'draws':>6} {'Psi=0':>6} {'M2':>3} {'M3':>3}  verdict")
ok = bad = 0
for lam in range(1, 5):
    for f in range(0, 4):
        n = 2*lam + f
        if n > 7: continue
        for r in [f+1, f+2, f+3]:
            if r < 1 or r > n: continue
            theta = list(range(r)); m = [theta[j]+j for j in range(r)]
            if len(set(m)) != r: continue
            z = t = m2 = m3 = 0
            for _ in range(10):
                w = draw(lam, f)
                if w is None: continue
                try:
                    v = psi_headsum(w, theta); v2 = psi_alt(w, theta)
                except ZeroDivisionError: continue
                assert v == v2, (w, theta, v, v2)
                t += 1
                if v != 0: continue
                z += 1
                live = [T for T in itertools.combinations(range(n), r)
                        if B(w, T)*alternant_ratio([w[i] for i in T], m) != 0]
                if live: m3 += 1
                else: m2 += 1
            if t == 0: continue
            pred = (r >= f + 2)
            agree = (z == t) if pred else (z == 0)
            ok += agree; bad += (not agree)
            print(f"{lam:>3} {f:>2} {n:>2} {r:>2} {str(pred):>7} {t:>6} {z:>6} {m2:>3} {m3:>3}"
                  f"  {'ok' if agree else '<<< LAW FAILS'}")
print()
print(f"sharp law  r >= f+2  <=>  Psi = 0 :  {ok}/{ok+bad} cells")
print("(lambda=1..4, f=0..3, n<=8, all r, 10 draws per cell, two implementations)")
