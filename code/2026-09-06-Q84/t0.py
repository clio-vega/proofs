"""t0.py -- the one specialization the hook witness cannot see.

Theorem 1 gives  (-1)^d t^{e-2}(1+t),  which vanishes at t=-1 (all e) and at
t=0 (e>=3).  So the hook family is BLIND at t=0 for e>=3 --
[[a-function-can-have-a-LADDER-of-degenerate-slices]].  R_e(0) = "add e boxes
in a single row".  Is IT of finite order?  Sweep other lam and other a-words.
"""
import itertools
from math import comb
from fast import op_p, ribbons

def op_R0(state, e):
    out = {}
    for lam, w in state.items():
        for (mu, h) in ribbons(lam, e):
            if h == 0: out[mu] = out.get(mu, 0) + w
    return {k: v for k, v in out.items() if v}

def nested(X, e, aa, lam):
    d = len(aa); out = {}
    for r in range(d + 1):
        for S in itertools.combinations(range(d), r):
            st = {lam: 1}
            for i in range(d):
                if i not in S: st = op_p(st, aa[i])
            st = X(st, e)
            for i in S: st = op_p(st, aa[i])
            for k, v in st.items(): out[k] = out.get(k, 0) + (-1)**r * v
    return {k: v for k, v in out.items() if v}

def parts(n, cap=None):
    if cap is None: cap = n
    if n == 0: return [()]
    return [(p,) + r for p in range(min(n, cap), 0, -1) for r in parts(n - p, p)]

lams = [l for n in range(0, 6) for l in parts(n)]
for e in (3, 4):
    print('R_%d(0)  ("add %d boxes in one row")' % (e, e))
    for d in range(1, 5):
        best = None; nz = 0; tot = 0
        for aa in itertools.product((1, 2), repeat=d):
            for lam in lams:
                tot += 1
                r = nested(op_R0, e, aa, lam)
                if r:
                    nz += 1
                    if best is None: best = (aa, lam, dict(sorted(r.items())[:2]))
        print('   C^(%d): %d/%d nonzero' % (d, nz, tot),
              ('  witness a=%s lam=%s -> %s' % best) if best else '  <== VANISHES')
