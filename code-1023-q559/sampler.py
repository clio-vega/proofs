from fractions import Fraction as F
import random
def clean_w(n, rng, lo=2, hi=40, dlo=1, dhi=7, tries=4000):
    """Sample n distinct rationals with NO degeneracy: w_a != +-1, w_a w_b != 1, w_a != w_b.
    The reciprocal-pair degeneracy w_a w_b = 1 kills every B(T) when n <= r+1 (measured
    2026-10-23 while diagnosing a 58/60 control).  Exclude it explicitly."""
    for _ in range(tries):
        w = [F(rng.randint(lo,hi), rng.randint(dlo,dhi)) for _ in range(n)]
        if len(set(w)) != n: continue
        if any(abs(x) == 1 for x in w): continue
        bad = False
        for a in range(n):
            for b in range(a+1,n):
                if w[a]*w[b] == 1: bad = True
        if bad: continue
        return w
    raise RuntimeError("sampler failed")
