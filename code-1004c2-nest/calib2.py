"""CALIBRATION 2: the y/width route must reproduce gen.slice_sum on real slices,
   and the recorded identities (eq:w, eq:Lam, k=(||y||-sigma)/2) must hold."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gen
from ycoord import *
from collections import Counter

st = Counter()
for m in range(2, 6):
    for n in range(m, m+5):
        for (mu, lam) in gen.pairs(n, m, 7):
            g = gaps(lam, n, m); G = n - m
            assert sum(g) == G, (g, G)
            st['pairs'] += 1
            for b in range(0, 2*n+1):
                nus = gen.slice_nus(mu, lam, n, m, b)
                if not nus: continue
                off, co = gen.slice_sum(mu, lam, n, m, b)
                if not co: continue
                st['slices'] += 1
                # y-route
                ws = []
                sigs = set()
                for nu in nus:
                    y = to_y(nu, lam, n, m)
                    w = wmap(y, g)
                    # eq:w check against the direct L_i,R_i widths
                    lr = gen.LR(nu, lam, n, m)
                    wdir = tuple(R-L+1 for (L, R) in lr)
                    st['eqw_ok' if w == wdir else 'eqw_BAD'] += 1
                    if min(w) <= 0: continue
                    sigs.add(sum(y))
                    # eq:Lam
                    Lam2 = G - sum(abs(t) for t in y)
                    st['eqLam_ok' if Lam2 == sum(w)-m else 'eqLam_BAD'] += 1
                    # k = (||y||-sigma)/2
                    st['eqk_ok' if 2*defect(y) == sum(abs(t) for t in y) - sum(y) else 'eqk_BAD'] += 1
                    ws.append(w)
                st['one_sigma' if len(sigs) <= 1 else 'MULTI_SIGMA'] += 1
                cs = centred_sum(ws)
                # compare with gen's coefficient list (trim zeros both sides)
                def trim(v):
                    nz = [i for i, c in enumerate(v) if c]
                    return tuple(v[nz[0]:nz[-1]+1]) if nz else ()
                st['sum_ok' if trim(cs) == trim(co) else 'sum_BAD'] += 1
print(dict(st))
